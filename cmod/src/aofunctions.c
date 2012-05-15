/* Numpy extension to bin up a Numpy image/array  */
//squash an image from dimsIn pixels into array of size dimsOut pixels.  Note, the input array must be an integer number of times larger than the output array.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw.h>




int binimg(void* inaddr,int ind1,int ind2,int ins1,int ins2,char intype,
	   void* outaddr,int outd1,int outd2,int outs1,int outs2,char outtype){
  //Requires input array address, dimensions, strides and data type, as well as
  //same for the output array.
  //You could call this function using ctypes, and for a numpy array would use eg:
  //binimg(arrin.ctypes.get_data(),arrin.shape[0],arrin.shape[1],arrin.strides[0],arrin.strides[1],arrin.dtype.char,
  //       arrout.ctypes.get_data(),arrout.shape[0],...)
  int nbin1,nbin2,i,j,k,l,kk,ll;
  double pxlout;
  nbin1=ind1/outd1;
  nbin2=ind2/outd2;
  if(intype!='d' || outtype!='d'){
    printf("binimg not yet working for non-double data types\n");
    return -1;
  }
  for(i=0; i<outd1; i++){
    for(j=0; j<outd2; j++){
      pxlout=0.;
      k=i*nbin1;
      for(kk=0; kk<nbin1; kk++){
	l=j*nbin2;
	for(ll=0; ll<nbin2; ll++){
	  pxlout+=*(double *)(inaddr+k*ins1+l*ins2);
	  l++;
	}
	k++;
      }
      *(double*)(outaddr+i*outs1+j*outs2)=pxlout;
    }
  }
  return 0;
}
typedef struct{
  void* fftarr;
  void* fft_plan;
  int dim1;
  int dim2;
  int size;
} mkimg_plan;

void* createmkimgplan(int imgoutd1,int imgoutd2){
  //set up the plan and the FFT data store... (avoids mallocing every time).
  //possibly change this to ACML - or change only if ACML exists?
  mkimg_plan* p;
  if((p=(mkimg_plan*)malloc(sizeof(mkimg_plan)))!=NULL){
    p->fftarr=malloc(sizeof(fftw_complex)*imgoutd1*imgoutd2);
    p->fftw_plan=fftw2d_create_plan(imgoutd1,imgoutd2,FFTW_FORWARD,FFTW_MEASURE | FFTW_IN_PLACE);//is this a pointer?
    p->dim1=imgoutd1;
    p->dim2=imgoutd2;
    p->size=imgoutd1*imgoutd2*sizeof(fftw_complex);
  }
  return (void*)p;
}
int freemkimgplan(void* p){
  //de-malloc the plan
  mkimg_plan *plan;
  plan=(mkimg_plan*)p;
  if(plan!=NULL){
    fftw2d_destroy_plan(plan->fftw_plan);//does this function exist?
    free(plan->fftarr);
    free(plan);
  }
  p=NULL;//doesn't do anything since only a reference...
  return 0;
}

int mkimg(void* p,void* phsinaddr,int phsind1,int phsind2,int phsins1,int phsins2,char phsintype,
	  void* imgoutaddr,int imgoutd1,int imgoutd2,int imgouts1,int imgouts2,char imgouttype,
	  void* pupaddr,int pupd1,int pupd2,int pups1,int pups2,char puptype){
  //inputs are a plan, the input phase address, dimensions, strides and type, the output image address dimensions strides and type and the pupil address dimensions strides and type.
  //note fft is always in float64 - have to include sfftw.h if want float32...
  mkimg_plan *plan;
  int i,j,ii;
  fftw_real pupil,phs,img;
  plan=(mkimg_plan*)p;
  //clear the temporary storage fft array...
  memset(plan->fftarr,0,plan->size);
  //Note, place the typecode if statements outside the loops... more efficient.
  if(phsintype=="d" && puptype=="d"){
    for(i=0;i<phsind1; i++){//prepare cos/sine of phase...
      for(j=0; j<phsind2; j++){
	pupil=*(double*)(pupaddr+i*pups1+j*pups2);
	phs=*(double*)(phsinaddr+i*phsins1+j*phsins2);
	ii=i*imgoutd1+j;
	plan->fftarr[ii].re=cos(phs)*pupil;
	plan->fftarr[ii].re=sin(phs)*pupil;
      }
    }
  }else if(phsintype=="f" && puptype=="f"){
    for(i=0;i<phsind1; i++){//prepare cos/sine of phase...
      for(j=0; j<phsind2; j++){
	pupil=(double)*(float*)(pupaddr+i*pups1+j*pups2);
	phs=(double)*(float*)(phsinaddr+i*phsins1+j*phsins2);
	ii=i*imgoutd1+j;
	plan->fftarr[ii].re=cos(phs)*pupil;
	plan->fftarr[ii].re=sin(phs)*pupil;
      }
    }
  }else if(phsintype=="f" && puptype=="c"){
    for(i=0;i<phsind1; i++){//prepare cos/sine of phase...
      for(j=0; j<phsind2; j++){
	pupil=(double)*(char*)(pupaddr+i*pups1+j*pups2);
	phs=(double)*(float*)(phsinaddr+i*phsins1+j*phsins2);
	ii=i*imgoutd1+j;
	plan->fftarr[ii].re=cos(phs)*pupil;
	plan->fftarr[ii].re=sin(phs)*pupil;
      }
    }
  }else{
    printf("Illegal types for arrays - phsin must be double or float, and pupil must be the same or char in which case phsin must be float.\n");
    return -1;
  }
  //do the fft
  fftwnd_one(plan->fftw_plan,plan->fftarr,NULL);
  if(imgouttype=="d"){
    for(i=0;i<imgoutd1;++i){//and prepare the image...
      for(j=0;j<imgoutd2;++j){
	ii=i*imgoutd1+j;
	img=(plan->fftarr[ii].re*plan->fftarr[ii].re)+(plan->fftarr[ii].im*plan->fftarr[ii].im);
	*(double *)(imgoutaddr+i*imgouts1+j*imgouts2) = (double)img;
      }
    }
  }else if(imgouttype=="f"){
    for(i=0;i<imgoutd1;++i){//and prepare the image...
      for(j=0;j<imgoutd2;++j){
	ii=i*imgoutd1+j;
	img=(plan->fftarr[ii].re*plan->fftarr[ii].re)+(plan->fftarr[ii].im*plan->fftarr[ii].im);
	*(float *)(imgoutaddr+i*imgouts1+j*imgouts2) = (float)img;
      }
    }
  }else{
    printf("Illegal type for imgout array - should be float or double\n");
    return -1;
  }
  return 0;
}
