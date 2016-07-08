#include <Python.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <complex.h>
#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#include <string.h>//memcpy etc.
#include "numpy/arrayobject.h"
#include "qsort.h"
static PyObject *CentError;

//#define simmalloc malloc
/*
  module to convert phase to centroids...
  optionally, multithreaded (for multicore processors).
  
  The idea here is that there would be an initialise call, which initialises stuff, and returns a python long (a pointer to the data structure).  You can then do a run call, which uses this and does the computation...
  
*/

#define CALSOURCE 1
#define SIG 2
#define ADDPOISSON 3
#define NCEN 4
#define READNOISE 5
#define READBG 6
#define NOISEFLOOR 7
#define SKYBRIGHTNESS 8
#define PXLPOWER 9
#define SPOTPSF 10
#define NTHREADS 11
#define OPTICALBINNING 12
#define CENTWEIGHT 13
#define CORRELATIONCENTROIDING 14
#define CORRTHRESH 15
#define CORRPATTERN 16
#define CALDATA 17
#define REFCENTS 18
#define CALCOEFF 19
#define USEBRIGHTEST 20
#define INTEGSTEPS 21
typedef struct{
  int n;
float *corr;
  void *centstr;
} centrunstruct;
typedef struct{
  int nthreads;
  int fftsize;
  int psfsize;//this is the size of the spot psf, or equal to clipsize otherwise
int fsize;//size of fftArrays allocated.
  int clipsize;//size after clipping FFT... img is then created by binning this down to nimg.
  int nimg;
  int phasesize;
  int preBinningFactor;
  int preBinningFactorOrig;
  int paddedsize;//either fftsize or if not binning before convolving, can be equal to psfsize.
  int nsubaps;
  //int nsubapsLatest;//how many subaps in the latest DMA.
  //int subapCnt;//how many subaps we've completed in total.
  int ncen;
  float readNoise;
  float readBg;
  int addPoisson;
  float noiseFloor;
  float sig;
  float skybrightness;
  float *skybrightnessArr;
  int calsource;
  int nintegrations;
  int maxIntegrations;
  int spotpsfDim;
  int centWeightDim;
  float *centWeight;
  int centWeightSize;
  float pxlPower;//not yet used.
  //int showCCDImg;
  //int allCents;
  //int centTotCnt;
  //uint64 phsAddr;
  //uint64 pupAddr;
  //uint64 spotpsfAddr;
  //uint64 centsAddr;
  //uint64 ccdImgAddr;
  unsigned long seed;
  //int nimg_v;
  int fftpxls;
  //int phasesize_v;
  int phasepxls;
  //int phasepxls_v;
  int imgpxls;
  //int imgpxls_v;
  //int nsubapInBuffer;
  //int dmatag;
  //dmagetstruct *phsDmaInfo;
  //dmagetstruct *pupDmaInfo;
  //dmagetstruct *psfDmaInfo;
int hllSize;
int fftArraySize;
int planSize;
  float **hll;//store the high light level image
  float *phs;//store the phase
  float *pupfn;//store the pupil
  float *spotpsf;//array 2 or 4D holding the PSFs - used for LGS spot elongation.
  float *sigArr;//optional array holding subap sig values if varying across ap - used for LGS.
  float *cents;
  //int ccdImgSize;
  //float *ccdImg;
  float *tiltfn;
  float *bimg;
  float *fracSubArea;
  int *subflag;
  pthread_t* threadid;
  centrunstruct **runinfo;
  int *subapStart;
  complex float **fftArrays;
  fftwf_plan fftplan;
  fftwf_plan rcfftplan;
  fftwf_plan crfftplan;
  complex float *fftPsf;
  float *binthispxl;
  float *binnextpxl;
  int *binindx;
  int opticalBinning;
  gsl_rng **gslRand;
  //float *imgmask;
  //float *ximgmask;
  //float *yimgmask;
  //int speid;
  pthread_attr_t *pthread_attr;
  float corrThresh;
  float *corrPattern;
  int correlationCentroiding;
  float *corrimg;//the output - ie what we'd be calculating centroids from.
  int corrsize;
  int corrPatternSize;
  fftwf_plan corrPlan;
  fftwf_plan invCorrPlan;
  float *refCents;//reference centroids
  int *calBounds;//for calibration of centroids.
  float *calData;
  float *calSteps;
  int calNSteps;
  int threshType;
  int imageOnly;
  float *calCoeff;
  int calNCoeff;
  int useBrightest;
  int *useBrightestArr;
  float **sortarr;
  int phaseStep;
  float fitMx[30];
  int fitMxParam;
  int gaussianFit;
  int parabolicFit;
  float gaussianMinVal;
  float gaussianReplaceVal;
} centstruct;



void flipQuad(float *hll, int clipsize, int h, int w){
  //flip a single quadrant... hll should be pointer to the start of this quad.
  int h2=h/2,w2=w/2;
  int i,j;
  float tmp;
  for(i=0; i<h2; i++){
    for(j=0; j<w; j++){
      tmp=hll[i*clipsize+j];
      hll[i*clipsize+j]=hll[(h-1-i)*clipsize+w-1-j];
      hll[(h-1-i)*clipsize+w-1-j]=tmp;
    }
  }
  if(h%2!=0){//do the last row
    i=h2;
    for(j=0; j<w2; j++){
      tmp=hll[i*clipsize+j];
      hll[i*clipsize+j]=hll[(h-1-i)*clipsize+w-1-j];
      hll[(h-1-i)*clipsize+w-1-j]=tmp;
    }
  }
}
int flipArray(int clipsize,float *hll,int mirror,float *tmparr){
  //reorganises the FFT'd array...
  //if mirror, will mirror in x and y (same as util.fliparray2).
  //If mirror==0, flips the quadrants but keeps them in the same place.  
  //Clipsize can be odd if mirror==0
  //If mirror==1, keeps quadrants unflipped, but swaps them round.
  //If mirror==1, clipsize must be even 
  //tmparr can be NULL or an float array with clipsize/2 elements.  Only used if mirror==1.
  int i,rtval=0,dofree=0;
  int f2;
  //float tmp;
  f2=clipsize/2;
  if(mirror==0){
    int h1,h2,h3,h4,w1,w2,w3,w4;
    w1=w2=clipsize/2;
    w3=w4=(clipsize+1)/2;
    h1=h3=w1;
    h2=h4=w3;
    flipQuad(hll,clipsize,h1,w1);
    flipQuad(&hll[h1*clipsize],clipsize,h2,w2);
    flipQuad(&hll[w1],clipsize,h3,w3);
    flipQuad(&hll[h1*clipsize+w1],clipsize,h4,w4);

      /*for(i=0; i<f2; i++){
      for(j=0; j<f4; j++){
	//first quadrant
	tmp=hll[i*clipsize+j];//store it temporarily.
	hll[i*clipsize+j]=hll[(f2-1-i)*clipsize+f2-j-1];//swap the values.
	hll[(f2-1-i)*clipsize+f2-j-1]=tmp;
	//second quadrant
	tmp=hll[(i+f2)*clipsize+j];//store it temporarily.
	hll[(i+f2)*clipsize+j]=hll[(clipsize-1-i)*clipsize+f2-j-1];//swap the values.
	hll[(clipsize-1-i)*clipsize+f2-j-1]=tmp;
	//3rd quad
	tmp=hll[i*clipsize+j+f2];//store it temporarily.
	hll[i*clipsize+j+f2]=hll[(f2-1-i)*clipsize+clipsize-j-1];//swap the values.
	hll[(f2-1-i)*clipsize+clipsize-j-1]=tmp;
	//4th quad
	tmp=hll[(i+f2)*clipsize+j+f2];//store it temporarily.
	hll[(i+f2)*clipsize+j+f2]=hll[(clipsize-1-i)*clipsize+clipsize-j-1];//swap the values.
	hll[(clipsize-1-i)*clipsize+clipsize-j-1]=tmp;
      }
      }*/
  }else{//swap the quadrants over.
    if(clipsize%2==0){
      if(tmparr==NULL){
	if((tmparr=malloc(sizeof(float)*f2))==NULL){//temporary storage allocation failed
	  printf("fliparray - temporary storage allocation failed, results will be wrong.\n");
	  rtval=-1;
	}else{
	  dofree=1;
	}
      }
      if(tmparr!=NULL){//temporary storage allocation successful
	for(i=0; i<f2; i++){
	  //1st and 4th quads
	  memcpy(tmparr,&hll[i*clipsize],sizeof(float)*f2);//copy to temporary storage
	  memcpy(&hll[i*clipsize],&hll[(i+f2)*clipsize+f2],sizeof(float)*f2);
	  memcpy(&hll[(i+f2)*clipsize+f2],tmparr,sizeof(float)*f2);
	  //2nd and 3rd quads
	  memcpy(tmparr,&hll[i*clipsize+f2],sizeof(float)*f2);//copy to temporary storage
	  memcpy(&hll[i*clipsize+f2],&hll[(i+f2)*clipsize],sizeof(float)*f2);
	  memcpy(&hll[(i+f2)*clipsize],tmparr,sizeof(float)*f2);
	}
	if(dofree==1)
	  free(tmparr);
      }    
    }else{//slightly more complex way...
      //first do a quadrand flip, and then flip this.  This is the
      //same result, but is correct for odd sized arrays...
      flipArray(clipsize,hll,0,NULL);
      flipQuad(hll,clipsize,clipsize,clipsize);
    }
  }
  return rtval;
}
int addPowerSpec(int fftsize,float *cdata, int psfsize,float *hll,int additive){
  //takes real and imaginary parts and sums them squared into hll.
  //The cdata is size fftsize x fftsize while hll is size psfsize x psfsize.  Note, that due to flip being called after this, we take the corners, not the middle.
  //cdata should be an array of size sizeof(float)*2*fftsize*fftsize, ie containing both real and imag parts interleaved.
  //cdata should be a complex float* array cast to float*.
  //If psfsize>fftsize, data is zeropadded around the edges.
  //If psfsize<fftsize (e.g. when not using a psf, and when clipping), then is clipped.
  int i,j,indx;
  //int ff=fftsize*fftsize;
  int c2;//=clipsize/2;
  int cIsOdd;//=clipsize%2;
  if(psfsize>fftsize){
    c2=fftsize/2;
    cIsOdd=fftsize%2;
  }else{
    c2=psfsize/2;
    cIsOdd=psfsize%2;
  }
  if(additive==0){//zero the first iteration...
    for(i=0; i<c2; i++){
      for(j=0; j<c2; j++){
	indx=(i*fftsize+j)*2;
	hll[i*psfsize+j]=cdata[indx]*cdata[indx]+cdata[indx+1]*cdata[indx+1];
	indx=(i*fftsize+j+fftsize-c2-cIsOdd)*2;
	hll[i*psfsize+j+psfsize-c2-cIsOdd]=cdata[indx]*cdata[indx]+cdata[indx+1]*cdata[indx+1];
	indx=((i+fftsize-c2-cIsOdd)*fftsize+j)*2;
	hll[(i+psfsize-c2-cIsOdd)*psfsize+j]=cdata[indx]*cdata[indx]+cdata[indx+1]*cdata[indx+1];
	indx=((i+fftsize-c2-cIsOdd)*fftsize+j+fftsize-c2-cIsOdd)*2;
	hll[(i+psfsize-c2-cIsOdd)*psfsize+j+psfsize-c2-cIsOdd]=cdata[indx]*cdata[indx]+cdata[indx+1]*cdata[indx+1];
      }
    }
    if(cIsOdd){
      for(i=0; i<c2; i++){
	indx=(i*fftsize+fftsize-1)*2;
	hll[i*psfsize+psfsize-1]=cdata[indx]*cdata[indx]+cdata[indx+1]*cdata[indx+1];
	indx=((fftsize-1)*fftsize+i)*2;
	hll[(psfsize-1)*psfsize+i]=cdata[indx]*cdata[indx]+cdata[indx+1]*cdata[indx+1];
	indx=((i+fftsize-c2-cIsOdd)*fftsize+fftsize-1)*2;
	hll[(i+psfsize-c2-cIsOdd)*psfsize+psfsize-1]=cdata[indx]*cdata[indx]+cdata[indx+1]*cdata[indx+1];
	indx=((fftsize-1)*fftsize+i+fftsize-c2-cIsOdd)*2;
	hll[(psfsize-1)*psfsize+i+psfsize-c2-cIsOdd]=cdata[indx]*cdata[indx]+cdata[indx+1]*cdata[indx+1];
      }
      indx=(fftsize*fftsize-1)*2;
      hll[psfsize*psfsize-1]=cdata[indx]*cdata[indx]+cdata[indx+1]*cdata[indx+1];
      
    }
    //for(i=0; i<ff; i++){
    //  hll[i]=cdata[i*2]*cdata[i*2]+cdata[i*2+1]*cdata[i*2+1];
    // }
  }else{
    for(i=0; i<c2; i++){
      for(j=0; j<c2; j++){
	indx=(i*fftsize+j)*2;
	hll[i*psfsize+j]+=cdata[indx]*cdata[indx]+cdata[indx+1]*cdata[indx+1];
	indx=(i*fftsize+j+fftsize-c2-cIsOdd)*2;
	hll[i*psfsize+j+psfsize-c2-cIsOdd]+=cdata[indx]*cdata[indx]+cdata[indx+1]*cdata[indx+1];
	indx=((i+fftsize-c2-cIsOdd)*fftsize+j)*2;
	hll[(i+psfsize-c2-cIsOdd)*psfsize+j]+=cdata[indx]*cdata[indx]+cdata[indx+1]*cdata[indx+1];
	indx=((i+fftsize-c2-cIsOdd)*fftsize+j+fftsize-c2-cIsOdd)*2;
	hll[(i+psfsize-c2-cIsOdd)*psfsize+j+psfsize-c2-cIsOdd]+=cdata[indx]*cdata[indx]+cdata[indx+1]*cdata[indx+1];
      }
    }
    if(cIsOdd){
      for(i=0; i<c2; i++){
	indx=(i*fftsize+fftsize-1)*2;
	hll[i*psfsize+psfsize-1]+=cdata[indx]*cdata[indx]+cdata[indx+1]*cdata[indx+1];
	indx=((fftsize-1)*fftsize+i)*2;
	hll[(psfsize-1)*psfsize+i]+=cdata[indx]*cdata[indx]+cdata[indx+1]*cdata[indx+1];
	indx=((i+fftsize-c2-cIsOdd)*fftsize+fftsize-1)*2;
	hll[(i+psfsize-c2-cIsOdd)*psfsize+psfsize-1]+=cdata[indx]*cdata[indx]+cdata[indx+1]*cdata[indx+1];
	indx=((fftsize-1)*fftsize+i+fftsize-c2-cIsOdd)*2;
	hll[(psfsize-1)*psfsize+i+psfsize-c2-cIsOdd]+=cdata[indx]*cdata[indx]+cdata[indx+1]*cdata[indx+1];
      }
      indx=(fftsize*fftsize-1)*2;
      hll[psfsize*psfsize-1]+=cdata[indx]*cdata[indx]+cdata[indx+1]*cdata[indx+1];
    }
  }
  return 0;
}

int addPowerSpecCentral(int fftsize,float *cdata, int psfsize,float *hll,int additive){
  //NEW, for when using tiltfn rather than fliparray.  This version takes the central part, rather then the corners.
  
  //takes real and imaginary parts and sums them squared into hll.
  //The cdata is size fftsize x fftsize while hll is size psfsize x psfsize.
  
  //cdata should be an array of size sizeof(float)*2*fftsize*fftsize, ie containing both real and imag parts interleaved.
  //cdata should be a complex float* array cast to float*.
  //If psfsize>fftsize, data is zeropadded around the edges.
  //If psfsize<fftsize (e.g. when not using a psf, and when clipping), then is clipped.
  int i,j,indx;
  //int ff=fftsize*fftsize;
  int c;
  int rowstart=0;
  int rowoff=0;
  if(psfsize>fftsize){
    rowstart=(psfsize-fftsize)/2;//colstart same since square.
    rowoff=0;//coloff same since square.
    c=fftsize;
  }else{
    rowstart=0;
    rowoff=(fftsize-psfsize)/2;
    c=psfsize;
  }
  if(additive==0){
    for(i=0;i<c;i++){
      for(j=0;j<c;j++){
	indx=((rowoff+i)*fftsize+rowoff+j)*2;
	hll[(i+rowstart)*psfsize+rowstart+j]=cdata[indx]*cdata[indx]+cdata[indx+1]*cdata[indx+1];
      }
    }
  }else{
    for(i=0;i<c;i++){
      for(j=0;j<c;j++){
	indx=((rowoff+i)*fftsize+rowoff+j)*2;
	hll[(i+rowstart)*psfsize+rowstart+j]+=cdata[indx]*cdata[indx]+cdata[indx+1]*cdata[indx+1];
      }
    }
  }
  return 0;
}

/*
int doFFT(int fftsize,int phasesize,int forward,float *re,float *im,workstruct *workbuf){
  //phasesize gives info about zero padding, and allows us to not do some of the ffts.
  //make use of phasesize - ie the first for loop should be up to p4 only 
  //(and memset parts of wr/wi to zero).
  //Does a 2D FFT.

  vector float *ar,*ai,*br,*bi,*wr,*wi;
  int f4,p4,i,j;
  ar=(vector float *)re;
  ai=(vector float *)im;
  br=(vector float *)workbuf->ffttmpre;//temporary array fftsize * fftsize in size.
  bi=(vector float *)workbuf->ffttmpim;
  wr=(vector float *)workbuf->singlefftre;//size fftsize (float)
  wi=(vector float *)workbuf->singlefftim;

  //printf("FFT Pointers: %p==%p, %p==%p %p==%p %p==%p %p==%p %p==%p\n",ar,re,ai,im,br,workbuf->ffttmpre,bi,workbuf->ffttmpim,wr,workbuf->singlefftre,wi,workbuf->singlefftim);
  //forward=1-forward;
  f4=fftsize/4;//fftsize always a multiple of 4.
  p4=(phasesize+3)/4;//divided by 4, rounded up.
  memset(br,0,fftsize*fftsize*sizeof(float));
  memset(bi,0,fftsize*fftsize*sizeof(float));
  for(i=0; i<p4; i++){
    fft_2d(&ar[fftsize*i],&ai[fftsize*i],wr,wi,forward);
    for(j=0; j<fftsize; j++){
      br[i+f4*j]=wr[j];
      bi[i+f4*j]=wi[j];
    }
  }
  for(i=0; i<f4; i++){
    fft_2d(&br[fftsize*i],&bi[fftsize*i],wr,wi,forward);
    for(j=0; j<fftsize; j++){
      ar[i+f4*j]=wr[j];
      ai[i+f4*j]=wi[j];
    }
  }
  return 0;
}
*/
int computeHll(float *phs,int phasesize,int niters,int fftsize,int paddedsize,float *pup,float *tiltfn,complex float *fftarr,float *hll,centstruct *c){
  int iter,i,j,rtval=0,indx;
  float tmp;
  memset(hll,0,sizeof(float)*paddedsize*paddedsize);
  for(iter=0; iter<niters; iter++){
    //for each simulation iteration to be binned...
    memset(fftarr,0,fftsize*fftsize*sizeof(complex float));
    if(c->phaseStep==1){
      for(i=0; i<phasesize; i++){
	for(j=0; j<phasesize; j++){
	  tmp=phs[(iter*phasesize*phasesize+i*phasesize+j)]-tiltfn[i*phasesize+j];
	  fftarr[i*fftsize+j]=pup[i*phasesize+j]*(cos(tmp)+I*sin(tmp));
	  //workbuf->phsIm[i*fftsize+j]=pup[i*phasesize+j]*sin(tmp);
	}
      }
    }else{
      for(i=0; i<phasesize; i++){
	for(j=0; j<phasesize; j++){
	  indx=(iter*phasesize*phasesize+i*phasesize+j)*c->phaseStep;
	  tmp=phs[indx]-tiltfn[i*phasesize+j];
	  fftarr[i*fftsize+j]=pup[i*phasesize+j]*phs[indx+1]*(cos(tmp)+I*sin(tmp));
	  //workbuf->phsIm[i*fftsize+j]=pup[i*phasesize+j]*sin(tmp);
	}
      }
    }
    fftwf_execute_dft(c->fftplan,fftarr,fftarr);
    //rtval|=doFFT(fftsize,phasesize,1,workbuf->phsRe,workbuf->phsIm,workbuf);//do the 2d fft
    rtval|=addPowerSpecCentral(fftsize,(float*)fftarr,paddedsize,hll,iter);//changed from addPowerSpec when extra tilt added to tiltfn. compute the power spectrum (ie high light level img) 
  }
  //printf("%d %d\n",fftsize,paddedsize);
  //rtval|=flipArray(paddedsize,hll,0,NULL);//removed when extra tilt added to tiltfn.
  return rtval;
}
void multArrArr(int datasize,float *data1,float *data2){
  int i;
  for(i=0; i<datasize; i++){
    data1[i]*=data2[i];
  }
}
int doConvolution(int psfsize,float *img,complex float *tmpbuf,complex float *fftPsf,centstruct *c,int threadno,int clipsize){
  //do a convolution of img with psf, and return the result in img.
  //Note, this destroys psf array.
  int rtval=0;
  int i,off;
  //Note, the flip and forward FFT of the spot pattern has already been taken during initialisation of the module.  It is this FFT'd version that is passed here, not the real thing.
  //First, FFT the hll img.
  //printf("Starting convolution %d\n",threadno);
  fftwf_execute_dft_r2c(c->rcfftplan,img,tmpbuf);
  //memset(workbuf->fftbuf,0,sizeof(float)*clipsize*clipsize);
  //rtval |= doFFT(clipsize,clipsize,1,img,workbuf->fftbuf,workbuf);
  //multiply the arrays together.
  //multArrArr(2*clipsize*(clipsize/2+1),(float*)tmpbuf,(float*)fftPsf);
  for(i=0; i<psfsize*(psfsize/2+1); i++){
    tmpbuf[i]*=fftPsf[i];
  }
  fftwf_execute_dft_c2r(c->crfftplan,tmpbuf,img);
  //And now need to clip img from psfsize to clipsize.
  if(clipsize<psfsize){
    off=(psfsize-clipsize)/2;
    for(i=0; i<clipsize; i++){
      memmove(&img[i*clipsize],&img[(i+off)*psfsize+off],clipsize*sizeof(float));
    }
  }


  /*if(psf!=NULL){//if psf==NULL, if means that the psf doesn't change between subaps, and that the inverse has already been computed, and stored in workbuf.

    rtval|=flipArray(clipsize,psf,1,workbuf->fliparr);
    memset(workbuf->invPsfIm,0,sizeof(float)*clipsize*clipsize);
    rtval|=doFFT(clipsize,clipsize,1,psf,workbuf->invPsfIm,workbuf);    //result is stored back in psf.
    //now multiply the arrays together on an element by element basis...
    multArrArr_v(clipsize*clipsize,img,psf);
    multArrArr_v(clipsize*clipsize,workbuf->fftbuf,workbuf->invPsfIm);
  }else{//multiply the arrays
    multArrArr_v(clipsize*clipsize,img,workbuf->invPsfRe);
    multArrArr_v(clipsize*clipsize,workbuf->fftbuf,workbuf->invPsfIm);
  }
  rtval|=doFFT(clipsize,clipsize,0,img,workbuf->fftbuf,workbuf);*/
  return rtval;
}
int prepareBinImage(int clipsize,  int nimg, centstruct *c){
  //set up some arrays used when binning the image.
  float *thispxl;//[clipsize];
  float *nextpxl;//[clipsize];
  int *indx;//[clipsize];
  int cnt=0,i,excess,rtval=0;
  //compute the helper arrays.
  if(c->binthispxl!=NULL)
    free(c->binthispxl);
  if(c->binnextpxl!=NULL)
    free(c->binnextpxl);
  if(c->binindx!=NULL)
    free(c->binindx);
  c->binindx=NULL;
  c->binnextpxl=NULL;
  c->binthispxl=NULL;
  if(clipsize>0){
    thispxl=malloc(sizeof(float)*clipsize);
    nextpxl=malloc(sizeof(float)*clipsize);
    indx=malloc(sizeof(int)*clipsize);
    c->binthispxl=thispxl;
    c->binnextpxl=nextpxl;
    c->binindx=indx;
    if(thispxl==NULL || nextpxl==NULL || indx==NULL){
      rtval=-1;
    }else{
      for(i=0; i<clipsize; i++){
	cnt+=nimg;
	if(cnt>=clipsize){
	  excess=cnt-clipsize;
	  thispxl[i]=(nimg-excess)/(float)nimg;
	  nextpxl[i]=excess/(float)nimg;
	  indx[i]=1;
	  cnt=excess;
	}else{
	  thispxl[i]=1.;//tmp[i]=(nimg&0x3f)
	  nextpxl[i]=0.;
	  indx[i]=0;//possibly better to implement this using bits, ands and ror or rol.
	}
      }
    }
  }
  return rtval;
}
int binImage(int clipsize, float *imgin,int nimg,float *imgout,centstruct *c){
  //bin the imgin into ingout.
  int i,j,ii,jj;
  float *thispxl,*nextpxl;
  int *indx;
  ii=0;
  if(clipsize==nimg){//no binning needed
    if(imgout!=imgin)
      memcpy(imgout,imgin,sizeof(float)*clipsize*clipsize);
  }else{
    memset(imgout,0,sizeof(float)*nimg*nimg);
    thispxl=c->binthispxl;
    nextpxl=c->binnextpxl;
    indx=c->binindx;
    for(i=0; i<clipsize-1; i++){
      jj=0;
      for(j=0; j<clipsize-1; j++){
	imgout[ii*nimg+jj]+=thispxl[i]*thispxl[j]*imgin[i*clipsize+j];
	if(jj+1<nimg)
	  imgout[ii*nimg+jj+1]+=thispxl[i]*nextpxl[j]*imgin[i*clipsize+j+1];
	if(ii+1<nimg){
	  imgout[(ii+1)*nimg+jj]+=nextpxl[i]*thispxl[j]*imgin[(i+1)*clipsize+j];
	  if(jj+1<nimg)
	    imgout[(ii+1)*nimg+jj+1]+=nextpxl[i]*nextpxl[j]*imgin[(i+1)*clipsize+j+1];
	}
	jj+=indx[j];
      }
      imgout[ii*nimg+jj]+=thispxl[i]*thispxl[j]*imgin[i*clipsize+j];
      if(ii+1<nimg)
	imgout[(ii+1)*nimg+jj]+=nextpxl[i]*thispxl[j]*imgin[(i+1)*clipsize+j];
      ii+=indx[i];
    }
    jj=0;
    for(j=0; j<clipsize-1; j++){
      imgout[ii*nimg+jj]+=thispxl[i]*thispxl[j]*imgin[i*clipsize+j];
      if(jj+1<nimg)
	imgout[ii*nimg+jj+1]+=thispxl[i]*nextpxl[j]*imgin[i*clipsize+j+1];
      jj+=indx[j];
    }
    imgout[ii*nimg+jj]+=thispxl[i]*thispxl[j]*imgin[i*clipsize+j];
  }
  return 0;
}

float getSum(int datasize,float *data){
  //note, if data is an image, datasize should be dimensions squared, eg nimg*nimg.
  int i;
  float sum=0.;
  for(i=0; i<datasize; i++){
    sum+=data[i];
  }
  return sum;
}
void multArr(int datasize,float *data,float multiplier){
  int i;
  for(i=0; i<datasize; i++){
    data[i]*=multiplier;
  }
}
void addArr(int datasize,float *data,float adder){
  int i;
  for(i=0; i<datasize; i++){
    data[i]+=adder;
  }
}
void addArrArr(int datasize,float *data,float *adder){
  int i;
  for(i=0; i<datasize; i++){
    data[i]+=adder[i];
  }
}


void poissonise(int datasize,float *data,const gsl_rng *r){
  int i;
  for(i=0; i<datasize; i++)
    data[i]=(float)gsl_ran_poisson(r,(double)data[i]);
}


void readCCD(int datasize,float *data,float noiseMean,float readNoise,float floor,const gsl_rng *r,int threshType){
  //adds read noise, and floors the data.
  int i;
  if(readNoise==0){
    if(threshType==0){
      for(i=0; i<datasize; i++){
	data[i]+=noiseMean;
	data[i]=data[i]<floor?0.:data[i]-floor;
      }
    }else if(threshType==1){
      for(i=0; i<datasize; i++){
	data[i]+=noiseMean;
	data[i]=data[i]<floor?0.:data[i];
      }
    }else{
      printf("Unknown threshold type - not adding read noise\n");
    }
  }else{
    if(threshType==0){//this one is probably best.  With floor being at least twice readNoise.
      for(i=0; i<datasize; i++){
	data[i]+=(float)gsl_ran_gaussian(r,(double)readNoise)+noiseMean;
	data[i]=data[i]<floor?0.:data[i]-floor;
      }
    }else if(threshType==1){
      for(i=0; i<datasize; i++){
	data[i]+=(float)gsl_ran_gaussian(r,(double)readNoise)+noiseMean;
	data[i]=data[i]<floor?0.:data[i];
      }
    }else{
      printf("Unknown threshold type - not adding read noise\n");
    }
  }
}

void opticallyBin(int nimg,int ncen,float *img){
  float v0=0;
  int s,e,i,j;
  //first bin the first row of pixels into v0.
  //Then bin all the columns into the first row.
  //Then bin all rows but the first, and store in the second row.
  //Then copy v0 into the first element of second row.
  s=(nimg-ncen)/2;//starting point
  e=nimg-s;//end point.
  
  for(i=s; i<e; i++){
    v0+=img[nimg*s+i];
  }
  //now sum columns into the first row... - ie the x part.
  if(s!=0)
    memset(img,0,sizeof(float)*ncen);
  for(i=s; i<e; i++){
    for(j=(s==0?1:s); j<e; j++){
      img[i-s]+=img[j*nimg+i];
    }
  }
  //and now sum the rows (have already done the first one).
  for(i=s+1; i<e; i++){
    if(s!=0)
      img[i-s+ncen]=0;
    for(j=(s==0?1:s); j<e; j++){
      img[i-s+ncen]+=img[i*nimg+j];
    }
  }
  img[ncen]=v0;
  //finally clear all the other pixels - makes better image.
  memset(&(img[2*ncen]),0,sizeof(float)*(nimg*nimg-2*ncen));
}

void computeOBCents(int ncen,float *img,float *cx,float *cy){
  int i;
  float sx=0.,sy=0.;
  *cx=0.;
  *cy=0.;
  for(i=0; i<ncen; i++){
    sx+=img[i];
    *cx+=img[i]*i;
    sy+=img[i+ncen];
    *cy+=img[i+ncen]*i;
  }
  if(sx>0){
    *cx/=sx;
    *cx-=(ncen-1)/2.;
  }
  if(sy>0){
    *cy/=sy;
    *cy-=(ncen-1)/2.;
  }
}
void computeCoG(int nimg,int ncen,float const * const img,float *centx,float *centy,float *weighting,int centreOnBrightest){
  //If centreOnBrightest set (i.e. if correlation centroiding), the ncenxncen array will be centred around the brightest pixel.
  int i,j,sx,sy,ex,ey,maxpos;
  float tmp;
  float sum=0.;
  float mx;
  float cx=0,cy=0;
  //*centx=0.;
  //*centy=0.;
  if(centreOnBrightest==0){
    sx=(nimg-ncen)/2;//starting point
    ex=nimg-sx;//end point
    sy=sx;
    ey=ex;
  }else{
    //find brightest pixel:
    mx=0.;
    maxpos=0;
    for(i=0;i<nimg*nimg;i++){
      if(img[i]>mx){
	mx=img[i];
	maxpos=i;
      }
    }
    sx=(int)(maxpos%nimg-ncen/2.+0.5);
    sy=(int)(maxpos/nimg-ncen/2.+0.5);
    if(sx<0)sx=0;
    if(sy<0)sy=0;
    if(sx>nimg-ncen)sx=nimg-ncen;
    if(sy>nimg-ncen)sy=nimg-ncen;
    ex=sx+ncen;
    ey=sy+ncen;
    
  }
  if(weighting==NULL){
    for(i=sy; i<ey; i++){
      for(j=sx; j<ex; j++){
	tmp=img[i*nimg+j];
	sum+=tmp;
	cx+=tmp*j;
	cy+=tmp*i;
      }
    }
  }else{
    for(i=sy; i<ey; i++){
      for(j=sx; j<ex; j++){
	tmp=img[i*nimg+j]*weighting[i*nimg+j];
	sum+=tmp;
	cx+=tmp*j;
	cy+=tmp*i;
      }
    }
  }
  if(sum>0){
    cx/=sum;
    cy/=sum;
    if(centreOnBrightest!=0){
      cx+=sx;
      cy+=sy;
    }
    cx-=(nimg-1)/2.;
    cy-=(nimg-1)/2.;
    *centx=cx;
    *centy=cy;
  }else{
    *centx=0;
    *centy=0;
  }
}

inline void makeFitVector(float *vec,float *subap,int nimg,int ncen,int startx,int starty){
  //vec should have size [6].
  int x,y;
  float val;
  for(y=0;y<6;y++)//probably faster than memset.
    vec[y]=0;
  for(y=0;y<ncen;y++){
    for(x=0;x<ncen;x++){
      val=subap[(y+starty)*nimg + (x+startx)];
      vec[0]+=val*x*x;
      vec[1]+=val*x*y;
      vec[2]+=val*y*y;
      vec[3]+=val*x;
      vec[4]+=val*y;
      vec[5]+=val;
    }
  }
}

void makeFitMx(float *fitMx,int ncen){//fitMx is 5x6.
  int y,x;
  double mx[36],inva[36];
  memset(mx,0,sizeof(double)*36);
  for(y=0;y<ncen;y++){
    for(x=0;x<ncen;x++){
      mx[0*6+0]+=x*x*x*x;
      mx[0*6+1]+=x*x*x*y;
      mx[0*6+2]+=x*x*y*y;
      mx[0*6+3]+=x*x*x;
      mx[0*6+4]+=x*x*y;
      mx[0*6+5]+=x*x;
      mx[1*6+0]+=x*x*x*y;
      mx[1*6+1]+=x*x*y*y;
      mx[1*6+2]+=x*y*y*y;
      mx[1*6+3]+=x*x*y;
      mx[1*6+4]+=x*y*y;
      mx[1*6+5]+=x*y;
      mx[2*6+0]+=x*x*y*x;
      mx[2*6+1]+=x*y*y*y;
      mx[2*6+2]+=y*y*y*y;
      mx[2*6+3]+=x*y*y;
      mx[2*6+4]+=y*y*y;
      mx[2*6+5]+=y*y;
      mx[3*6+0]+=x*x*x;
      mx[3*6+1]+=x*x*y;
      mx[3*6+2]+=x*y*y;
      mx[3*6+3]+=x*x;
      mx[3*6+4]+=x*y;
      mx[3*6+5]+=x;
      mx[4*6+0]+=x*x*y;
      mx[4*6+1]+=x*y*y;
      mx[4*6+2]+=y*y*y;
      mx[4*6+3]+=x*y;
      mx[4*6+4]+=y*y;
      mx[4*6+5]+=y;
      mx[5*6+0]+=x*x;
      mx[5*6+1]+=x*y;
      mx[5*6+2]+=y*y;
      mx[5*6+3]+=x;
      mx[5*6+4]+=y;
      mx[5*6+5]+=1;
      
    }
  }
  //now invert.
  int s;
  gsl_matrix_view m = gsl_matrix_view_array(mx, 6, 6);
  gsl_matrix_view inv = gsl_matrix_view_array(inva,6,6);
  gsl_permutation * p = gsl_permutation_alloc (6);
  gsl_linalg_LU_decomp (&m.matrix, p, &s);    
  gsl_linalg_LU_invert (&m.matrix, p, &inv.matrix);
  gsl_permutation_free (p);
  for(y=0;y<30;y++){//don't need the last row.
    fitMx[y]=(float)inva[y];
  }
}
void sgemv(int m,int n,float *mx, float *vec, float *res){
  int i,j;
  float tmp;
  for(i=0;i<m;i++){
    tmp=0;
    for(j=0;j<n;j++){
      tmp+=mx[i*n+j]*vec[j];
    }
    res[i]=tmp;
  }
}

void computeFit(const int nimg,const int ncen,float *img,float *centx,float *centy, float *fitMx,int *fitMxParam,int doGaussian,float gaussianMinVal,float gaussianReplaceVal){
  int i,mxpos=0;
  float sum=0;
  float mx=0;
  int startx,starty;
  float minflux=0;
  //Do a parabolic/gaussian fit...
  if(doGaussian){
    //take the log of the data
    for(i=0;i<nimg*nimg;i++){
      sum+=img[i];
      if(img[i]>mx){
	mxpos=i;
	mx=img[i];
      }
      if(img[i]<gaussianMinVal)
	img[i]=gaussianReplaceVal;
      else
	img[i]=logf(img[i]);
    }
  }else{
    for(i=0;i<nimg*nimg;i++){
      sum+=img[i];
      if(img[i]>mx){
        mxpos=i;
	mx=img[i];
      }
    }
  }
  if(sum>=minflux){
    float vec[6];
    float res[5];
    //Should this be centred around brightest pixel, or CoG?????
    //Since usually used with correlation images, brightest pixel probably ok (and maybe best - use it for now!)
    startx=(int)(mxpos%nimg-ncen/2.+0.5);
    starty=(int)(mxpos/nimg-ncen/2.+0.5);
    if(startx<0)startx=0;
    if(starty<0)starty=0;
    if(startx>nimg-ncen)startx=nimg-ncen;
    if(starty>nimg-ncen)starty=nimg-ncen;
    makeFitVector(vec,img,nimg,ncen,startx,starty);
    //dot the matrix with the vector.
    if(*fitMxParam!=ncen){
      printf("centmodule: Creating fit matrix size %d\n",ncen);
      makeFitMx(fitMx,ncen);//shape of matrix doens't change, but contents does.
      *fitMxParam=ncen;
    }
    sgemv(5,6,fitMx,vec,res);
    if(doGaussian){
      *centx=-res[3]/(2*res[0]);
      *centy=-res[4]/(2*res[2]);
    }else{
      *centx=(res[1]*res[4]/(2*res[2])-res[3])/(2.*res[0]-res[1]*res[1]/(2.*res[2]));
      *centy=-(res[4]+res[1]*(*centx))/(2.*res[2]);
    }
    (*centy)-=ncen/2.-0.5;
    (*centx)-=ncen/2.-0.5;
    (*centx)+=startx+ncen/2.-nimg/2.;
    (*centy)+=starty+ncen/2.-nimg/2.;
  }else{
    *centx=0;
    *centy=0;
  }
  //printf("parabolicFit: %d %d %d %g %g %g %g\n",startx,starty,mxpos,mxpos%nimg-nimg/2.-0.5,mxpos/nimg-nimg/2.-0.5,*centx,*centy);
}


/*
int testNan(int n,float* arr){
  int i;
  int got=0;
  for(i=0; i<n; i++){
    if(isnan(arr[i])){
      printf("testNan %d ",i);
      got=1;
    }
  }
  return got;
  }*/

#define B(y,x) corrPattern[(y)*corrsize+x]

int calcCorrelation(int nimg,int corrsize,float *corrPattern,float *bimg,float *corrimg,float *tmp,fftwf_plan corrPlan,  fftwf_plan invCorrPlan){
  /*This code has been copied from the DARC RTC code, for maximum efficiency.
   */
  int i=0,j,n,m,neven,meven;
  float *a;
  float r1,r2,r3,r4,r5,r6,r7,r8;
  //This is how the plans should be created (elsewhere).  Will need a different plan for each different sized subap (see subapLocation).  
  //c->corrPlan=fftwf_plan_r2r_2d(nimg,nimg,tmpsubap,tmpsubap2,FFTW_R2HC, FFTW_R2HC, FFTW_ESTIMATE);
  //c->invCorrPlan=fftwf_plan_r2r_2d(nimg,nimg,tmpsubap2,tmpsubap2,FFTW_HC2R, FFTW_HC2R, FFTW_ESTIMATE);
  
  //FFT the SH image.
  if(nimg==corrsize)
    memcpy(tmp,bimg,sizeof(float)*nimg*nimg);
  else{
    //copy with padding...
    memset(tmp,0,sizeof(float)*corrsize*corrsize);
    for(i=0;i<nimg;i++){
      memcpy(&tmp[((corrsize-nimg)/2)*corrsize+(corrsize-nimg)/2+i*corrsize],&bimg[i*nimg],sizeof(float)*nimg);
    }
  }
  fftwf_execute_r2r(corrPlan,tmp,tmp);
  //memcpy(bimg,tmp,sizeof(float)*nimg*nimg);
  
  //Now multiply by the reference...
  //This is fairly complicated due to the half complex format.  If you need to edit this, make sure you know what you're doing.
  a=tmp;
  n=corrsize;
  m=corrsize;
  a[0]*=B(0,0);
  neven=(n%2==0);
  meven=(m%2==0);
  if(neven){
    a[n/2]*=B(0,n/2);
  }
  if(meven){
    a[n*m/2]*=B(m/2,0);
    if(neven){
      a[n*m/2+n/2]*=B(m/2,n/2);
    }
  }
  for(i=1; i<(n+1)/2; i++){
    r1=a[i]*B(0,i)-a[n-i]*B(0,n-i);
    r2=a[i]*B(0,n-i)+a[n-i]*B(0,i);
    a[i]=r1;
    a[n-i]=r2;
    if(meven){
      r3=a[m/2*n+i]*B(m/2,i)-a[m/2*n+n-i]*B(m/2,n-i);
      r4=a[m/2*n+i]*B(m/2,n-i)+a[m/2*n+n-i]*B(m/2,i);
      a[m/2*n+i]=r3;
      a[m/2*n+n-i]=r4;
    }
  }
  
  for(i=1; i<(m+1)/2; i++){
    //do the 4 rows/cols that only require 2 values...
    r5=a[i*n]*B(i,0)-a[(m-i)*n]*B(m-i,0);
    r6=a[i*n]*B(m-i,0)+a[(m-i)*n]*B(i,0);
    a[i*n]=r5;
    a[(m-i)*n]=r6;
    if(neven){
      r7=a[i*n+n/2]*B(i,n/2)-a[(m-i)*n+n/2]*B(m-i,n/2);
      r8=a[i*n+n/2]*B(m-i,n/2)+a[(m-i)*n+n/2]*B(i,n/2);
      a[i*n+n/2]=r7;
      a[(m-i)*n+n/2]=r8; // IS THIS A BUG?? SHOULD IT BE "=r8" INSTEAD?? (UB, 2012Aug08)  Yes, I think so - changed to r8 on 120830 by agb.
    }
    
    for(j=1; j<(n+1)/2; j++){
      //and now loop over the rest.
      r1=a[i*n+j]*B(i,j)+a[(m-i)*n+n-j]*B(m-i,n-j)-a[i*n+n-j]*B(i,n-j)-a[(m-i)*n+j]*B(m-i,j);
      r2=a[i*n+j]*B(m-i,n-j)+a[(m-i)*n+n-j]*B(i,j)+a[(m-i)*n+j]*B(i,n-j)+a[i*n+n-j]*B(m-i,j);
      r3=a[i*n+j]*B(i,n-j)-a[(m-i)*n+n-j]*B(m-i,j)+a[i*n+n-j]*B(i,j)-a[(m-i)*n+j]*B(m-i,n-j);
      r4=a[i*n+j]*B(m-i,j)-a[(m-i)*n+n-j]*B(i,n-j)+a[(m-i)*n+j]*B(i,j)-a[i*n+n-j]*B(m-i,n-j);
      a[i*n+j]=r1;
      a[(m-i)*n+n-j]=r2;
      a[i*n+n-j]=r3;
      a[(m-i)*n+j]=r4;
    }
  }
  //and now do the inverse fft...
  fftwf_execute_r2r(invCorrPlan,tmp,tmp);
  if(corrimg!=NULL)
    memcpy(corrimg,tmp,sizeof(float)*corrsize*corrsize);
  return 0;

}
#undef B
int thresholdCorrelation(int npxl,float thresh,float *img){
  int i;
  float mx;
  if(thresh<=1. && thresh>0.){
    //threshold as a fraction of maximum.
    //So - find the maximum.
    mx=0;
    for(i=0; i<npxl; i++){
      if(img[i]>mx)
	mx=img[i];
    }
    thresh=mx*thresh;
  }
  for(i=0; i<npxl; i++){
    img[i]=img[i]<thresh?0:img[i]-thresh;
  }
  return 0;
}
float interp(float x,int n,float *xp,float *fp){
  //n is number of entries of xp minus 1.
  int i;
  if(x<=xp[0])
    return fp[0];
  if(x>=xp[n])
    return fp[n];
  i=1;
  while(x>xp[i])
    i++;
  //now do the linear interp
  return (x-xp[i-1])*(fp[i]-fp[i-1])/(xp[i]-xp[i-1])+fp[i-1];
}

int applyCentroidCalibrationInterpolation(float *cx,int n,float *coeff){
  int i;
  float res=0;//coeff[0];
  float xc=1.;
  for(i=0;i<n;i++){
    res+=xc*coeff[i];
    xc*=*cx;   //what a cool line of code!!!
  } 
  *cx=res;//replace with the calibrated centroid.
  return 0;
}

int applyCentroidCalibration(float *cx,int n, float *calData,float *calSteps){
  //This is a c version of applyCalibrationUnique
  //What should calData/bound/step be?
  //n is the number of entries in calSteps etc minus 1.  So, it is the index of the last entry.

  //if(*cx>calData[calBound[0,i,j,1],0,i,j] || *cx<calData[calBound[0,i,j,0],0,i,j])
  if(*cx>calData[n] || *cx<calData[0])
    printf("Warning: centroid with value %g is outside calibrated bounds\n",*cx);
  *cx=interp(*cx,n,calData,calSteps);
  return 0;
}

int selectBrightestPixels(int npxl,float *img,int useBrightest,float *sort){
  int subtract=0;
  float thr,sub;
  int i;
  if(useBrightest<0){
    useBrightest=-useBrightest;
    subtract=1;
  }
  if(useBrightest>0 && useBrightest<npxl){
    memcpy(sort,img,sizeof(float)*npxl);
#define cmp(a,b) ((*a)<(*b))
    QSORT(float,sort,npxl,cmp);
#undef cmp
    //The threshold to use is:
    thr=sort[npxl-useBrightest];
    if(subtract){//want to subtract the next brightest pixel
      subtract=npxl-useBrightest-1;
      while(subtract>=0 && sort[subtract]==thr)
	subtract--;
      if(subtract>=0)
	sub=sort[subtract];
      else
	sub=0;
      for(i=0; i<npxl; i++){
	if(img[i]<thr)
	  img[i]=0;
	else
	  img[i]-=sub;
      }
    }else{
      for(i=0; i<npxl; i++){
	if(img[i]<thr)
	  img[i]=0;
      }
    }
  }
  return 0;
}

int prebinImage(int paddedsize,float *hll,int preBinningFactor,int psfsize){
  //Bin from paddedsize downto binnedsize, and then clip to psfsize.
  //Can do the clipping by binning the correct parts...

  int binnedsize;
  int scaled=paddedsize/preBinningFactor;
  float tmp;
  int offset=0;
  int i,j,ii,jj;
  binnedsize=scaled;
  if(scaled>psfsize){//need to clip.
    offset=(scaled-psfsize)/2;
    binnedsize=psfsize;
  }

  for(i=0;i<binnedsize;i++){
    for(j=0;j<binnedsize;j++){
      tmp=0;
      for(ii=0;ii<preBinningFactor;ii++){
	for(jj=0;jj<preBinningFactor;jj++){
	  tmp+=hll[(i*preBinningFactor+ii+offset)*paddedsize+j*preBinningFactor+jj+offset];
	}
      }
      hll[i*binnedsize+j]=tmp;
    }
  }
  if(scaled<psfsize){//need to zero-pad: move image to centre of array.
    offset=(psfsize-scaled)/2;
    for(i=scaled-1;i>=0;i--){
      memcpy(&hll[(i+offset)*psfsize+offset],&hll[i*scaled],sizeof(float)*scaled);
    }
    //now clear round the edge.
    for(i=0;i<offset;i++)//below
      memset(&hll[i*psfsize],0,sizeof(float)*psfsize);
    for(i=offset+scaled;i<psfsize;i++)//above
      memset(&hll[i*psfsize],0,sizeof(float)*psfsize);
    for(i=offset;i<offset+scaled;i++){
      memset(&hll[i*psfsize],0,sizeof(float)*offset);//left
      memset(&hll[i*psfsize+offset+scaled],0,sizeof(float)*(psfsize-offset-scaled));//right
    }
  }

  return 0;
}

int centroidsFromPhase(centrunstruct *runinfo){
  centstruct *c;
  int threadno;
  int i;
  int npxl,nexposed;
  int error=0;
  float nphspxl;
  float totsig;
  float *bimg,*cweight;
  int indx,n;
  int imageOnly;
  c=(centstruct*)runinfo->centstr;
  imageOnly=c->imageOnly;
  threadno=runinfo->n;
  //printf("threadno %d doing subaps %d to %d (c %p)\n",threadno,c->subapStart[threadno],c->subapStart[threadno+1],c);
  //compute the requested centroids.
  for(i=c->subapStart[threadno]; i<c->subapStart[threadno+1]; i++){
    //memset(&(c->bimg[i*npxl]),0,npxl*sizeof(float));
    npxl=c->imgpxls;
    if(c->subflag[i]==1){//subap is full enough to use.  ie subflag[i]==1...

      nphspxl=c->fracSubArea[i];//getSum(c->phasepxls,&c->pup[i*c->phasepxls])/(float)c->phasepxls;//fraction of active pixels.
      //printf("computehll %d %d\n",threadno,i);
      //phaseStep==1 for phaseOnly, 2 for phaseamp.
      error|=computeHll(&(c->phs[i*c->phasepxls*c->maxIntegrations*c->phaseStep]),c->phasesize,c->nintegrations,c->fftsize,c->paddedsize,&(c->pupfn[i*c->phasepxls]),c->tiltfn,c->fftArrays[threadno],c->hll[threadno],c);//paddedsize was psfsize.
      //Optional - do a binning here, before convolution.  Why?  To reduce the size necessary for the convolution and thus improve performance.  Note, a binning after convolution is also usually necessary to give correct results.  This option is useful when trying to get larger pixel scales (i.e. larger phase size) without needing huge psfs.
      
      if(c->preBinningFactor>1)//The binned image becomes paddedsize/prebinfact, and is then further clipped (or, usually padded) to psfsize.
	error|=prebinImage(c->paddedsize,c->hll[threadno],c->preBinningFactor,c->psfsize);


      if(c->spotpsfDim==2){
	error|=doConvolution(c->psfsize,c->hll[threadno],c->fftArrays[threadno],c->fftPsf,c,threadno,c->clipsize);
	//error|=doConvolution(c->fftsize,c->hll,NULL,workbuf);//inverse already computed and stored in workbuf.
      }else if(c->spotpsfDim==4){//get the right psf for this subap.
	error|=doConvolution(c->psfsize,c->hll[threadno],c->fftArrays[threadno],&(c->fftPsf[i*c->psfsize*(c->psfsize/2+1)]),c,threadno,c->clipsize);
	//error|=doConvolution(c->fftsize,c->hll,&(c->spotpsf[i*c->fftpxls]),workbuf);
      }
      //if(testNan(c->fftsize*c->fftsize,c->hll[threadno])){
      //printf("Got NaN %d\n",i);
      //}
      error|=binImage(c->clipsize,c->hll[threadno],c->nimg,&(c->bimg[i*npxl]),c);
      //if(testNan(c->nimg*c->nimg,&(c->bimg[i*npxl]))){
      //printf("Got NaN bimg %d\n",i);
      //}
      //need to set the unneeded parts of bimg to zero for the vector sum.
      totsig=getSum(npxl,&(c->bimg[i*npxl]));
      if(totsig>0){//normalise the intensity in each subap
	if(c->sigArr==NULL){
	  multArr(npxl,&(c->bimg[i*npxl]),c->sig*nphspxl/totsig);
	}else{
	  multArr(npxl,&(c->bimg[i*npxl]),c->sigArr[i]/totsig);
	}
      }
      //if(testNan(c->nimg*c->nimg,&(c->bimg[i*npxl]))){
      //printf("Got NaN bimg2 %d\n",i);
      //}
      //if(c->calsource==0){//not a calibration source, so add noise...
      //As of 100312, have taken calsource out, and now add background, but not noise (shot or read) if a calibration source.  The floor is also included regardless
      if(c->skybrightness*nphspxl>0)
	addArr(npxl,&(c->bimg[i*npxl]),c->skybrightness*nphspxl);
      if(c->skybrightnessArr!=NULL)
	addArrArr(npxl,&(c->bimg[i*npxl]),&(c->skybrightnessArr[i*npxl]));
      if(c->opticalBinning==1){
	opticallyBin(c->nimg,c->ncen,&(c->bimg[i*npxl]));
	nexposed=c->ncen*2;
      }else{
	nexposed=npxl;
      }
      if((c->addPoisson==1 && c->calsource==0) && (totsig>0 || c->skybrightness*nphspxl>0 || c->skybrightnessArr!=NULL))
	poissonise(nexposed,&(c->bimg[i*npxl]),c->gslRand[threadno]);
      readCCD(nexposed,&(c->bimg[i*npxl]),c->readBg,c->readNoise*(c->calsource==0),c->noiseFloor,c->gslRand[threadno],c->threshType);
      if(c->calsource==0)
	selectBrightestPixels(nexposed,&c->bimg[i*npxl],c->useBrightestArr==NULL?c->useBrightest:c->useBrightestArr[i],c->sortarr[threadno]);
 	//}else{//a calibration source.
	//if(c->opticalBinning==1){//optically bin the calibration source...
      // opticallyBin(c->nimg,c->ncen,&(c->bimg[i*npxl]));
	//	}
	//}
      //if(testNan(c->nimg*c->nimg,&(c->bimg[i*npxl]))){
      //printf("Got NaN bimg3 %d\n",i);
      //}
      if(imageOnly==0){
	if(c->correlationCentroiding==1 && c->corrPattern!=NULL){//compute the correlation
	  //This shouldn't be used with optical binning...
	  //The correlation image may be bigger than nimg.  So, zeropad bimg, then do the correlation.  After that, use ncen.
	  int npxlcorr=c->corrsize*c->corrsize;
	  bimg=&(c->corrimg[i*npxlcorr]);
	  calcCorrelation(c->nimg,c->corrsize,&(c->corrPattern[i*npxlcorr]),&(c->bimg[i*npxl]),bimg,runinfo->corr,c->corrPlan,c->invCorrPlan);
	  npxl=npxlcorr;//The subaps may have grown!!!
	  thresholdCorrelation(npxl,c->corrThresh,bimg);
	}else{
	  bimg=&(c->bimg[i*npxl]);
	}
	//now apply the centroid mask and compute the centroids.
	if(c->opticalBinning==1){
	  computeOBCents(c->ncen,&(c->bimg[i*npxl]),&(c->cents[i*2]),&(c->cents[i*2+1]));
	}else if(c->parabolicFit==1 || c->gaussianFit==1){//do a parabolic/gaussian fit
	  computeFit(c->corrsize,c->ncen,bimg,&c->cents[i*2],&c->cents[i*2+1],c->fitMx,&c->fitMxParam,c->gaussianFit,c->gaussianMinVal,c->gaussianReplaceVal);

	}else{//do a CoG
	  if(c->centWeightDim==2){//weighted CoG
	    cweight=c->centWeight;
	  }else if(c->centWeightDim==4 && c->centWeight!=NULL){//weighted CoG
	    cweight=&c->centWeight[i*npxl];
	  }else{//cog
	    cweight=NULL;
	  }
	  //Note: corrsize==nimg unless convolving with something larger.  
	  computeCoG(c->corrsize,c->ncen,bimg,&(c->cents[i*2]),&(c->cents[i*2+1]),cweight,c->correlationCentroiding);
	}
	//and now apply the calibration... (linearisation)
	if(c->calData!=NULL){
	  indx=c->calBounds[i*2];//c->calBounds[0,i,j,0];
	  n=c->calBounds[i*2+1]-indx;
	  applyCentroidCalibration(&(c->cents[i*2]),n,&c->calData[i*c->calNSteps+indx],&c->calSteps[indx]);//calData[0,i,j,indx]
	  indx=c->calBounds[(c->nsubaps+i)*2];//[1,i,j,0];
	  n=c->calBounds[(c->nsubaps+i)*2+1]-indx;
	  applyCentroidCalibration(&(c->cents[i*2+1]),n,&c->calData[(c->nsubaps+i)*c->calNSteps+indx],&c->calSteps[indx]);
	}else if(c->calCoeff!=NULL){//calibrate using a polynomial...
	  
	  applyCentroidCalibrationInterpolation(&(c->cents[i*2]),c->calNCoeff,&c->calCoeff[i*2*c->calNCoeff]);
	  applyCentroidCalibrationInterpolation(&(c->cents[i*2+1]),c->calNCoeff,&c->calCoeff[(i*2+1)*c->calNCoeff]);
	}
	//and now subtract reference centroids...
	if(c->refCents!=NULL){
	  c->cents[i*2]-=c->refCents[i*2];
	  c->cents[i*2+1]-=c->refCents[i*2+1];
	  //printf("%g %g - after refsub\n",c->refCents[i*2],c->refCents[i*2+1]);
	}
      }
    }else{
      //return ignored centroids too.
      if(c->imageOnly==0){
	c->cents[i*2]=0.;
	c->cents[i*2+1]=0.;
      }else{//unused subaps still have noise...
	//add random pixel noise...
	nexposed=npxl;
	memset(&(c->bimg[i*npxl]),0,nexposed*sizeof(float));
	readCCD(nexposed,&(c->bimg[i*npxl]),c->readBg,c->readNoise*(c->calsource==0),c->noiseFloor,c->gslRand[threadno],c->threshType);
	
      }
      //if(testNan(c->nimg*c->nimg,&(c->bimg[i*npxl]))){
      //printf("Got NaN bimg4 %d\n",i);
      //}
    }
  }
  return error;
}


int checkFloatContigArr(PyArrayObject *arr){
  //checks whether an array is contiguous or not, and whether it is float..
  int rtval=0;//return 0 on success.
  if(arr->descr->type_num!=NPY_FLOAT)
    rtval=1;
  if(PyArray_ISCONTIGUOUS(arr)==0)//contiguous?
    rtval=1;
  return rtval;
}
int checkIntContigArr(PyArrayObject *arr){
  //checks whether an array is contiguous or not, and whether it is float..
  int rtval=0;//return 0 on success.
  if(arr->descr->kind!='i' || arr->descr->elsize!=sizeof(int))
    rtval=1;
  if(PyArray_ISCONTIGUOUS(arr)==0)//contiguous?
    rtval=1;
  return rtval;
}
int setCalData(centstruct *c,PyObject *dataTuple){
  PyArrayObject *calStep,*calBounds,*calData;
  if(!PyTuple_Check(dataTuple)){
    if(dataTuple==Py_None){
      c->calSteps=NULL;
      c->calData=NULL;
      c->calBounds=NULL;
      c->calNSteps=0;
      return 0;
    }else{
      printf("Error: centmodule - calibration stuff not a tuple\n");
      return -1;
    }
  }
  if((calData=(PyArrayObject*)PyTuple_GetItem(dataTuple,0))==NULL){
    printf("Error: centmodule - didn't get calibrateData\n");
    return -1;
  }
  if((calBounds=(PyArrayObject*)PyTuple_GetItem(dataTuple,1))==NULL){
    printf("Error: centmodule - didn't get calibrateBounds\n");
    return -1;
  }
  if((calStep=(PyArrayObject*)PyTuple_GetItem(dataTuple,2))==NULL){
    printf("Error: centmodule - didn't get calibrateStep\n");
    return -1;
  }
  if(checkFloatContigArr(calData)!=0){
    printf("Error: centmodule - calData should be float32 contiguous\n");
    return -1;
  }
  if(checkFloatContigArr(calStep)!=0){
    printf("Error: centmodule - calStep should be float32 contiguous\n");
    return -1;
  }
  if(checkIntContigArr(calBounds)!=0){
    printf("Error: centmodule - calBounds should be int32 contiguous\n");
    return -1;
  }
  if(calStep->nd!=1){
    printf("Error: centmodule - calStep should be 1d\n");
    return -1;
  }
  c->calNSteps=calStep->dimensions[0];
  if(calData->nd!=4 || calData->dimensions[0]!=2 || calData->dimensions[1]*calData->dimensions[2]!=c->nsubaps || calData->dimensions[3]!=c->calNSteps){
    printf("Error: centmodule - calData should be 4d, (2,nsubx,nsubx,steps)\n");
    return -1;
  }
  if(calBounds->nd!=4 || calBounds->dimensions[0]!=2 || calBounds->dimensions[1]*calBounds->dimensions[2]!=c->nsubaps || calBounds->dimensions[3]!=2){
    printf("Error: centmodule - calBounds should be 4d (2,nsubx,nsubx,2)\n");
    return -1;
  }
  c->calSteps=(float*)calStep->data;
  c->calBounds=(int*)calBounds->data;
  c->calData=(float*)calData->data;
  return 0;
}

int setCalCoeff(centstruct *c,PyArrayObject *calCoeff){
  if(calCoeff==NULL){
    c->calCoeff=NULL;
    c->calNCoeff=0;
    return 0;
  }
  if(checkFloatContigArr(calCoeff)!=0){
    printf("Error: centmodule - calCoeff should be float32 contiguous\n");
    return -1;
  }
  if(calCoeff->nd!=4 || calCoeff->dimensions[2]!=2 || calCoeff->dimensions[1]*calCoeff->dimensions[0]!=c->nsubaps){
    printf("Error: centmodule - calCoeff should be 4d (nsubx,nsubx,2,ncoeff)\n");
    return -1;
  }
  c->calNCoeff=(int)calCoeff->dimensions[3];
  c->calCoeff=(float*)calCoeff->data;
  return 0;
}


int setCentWeight(centstruct *c,PyObject *centWeightObj){
  PyArrayObject *centWeightArr;
  c->centWeightDim=0;
  c->centWeight=NULL;
  if(PyArray_Check(centWeightObj)==1){
    centWeightArr=(PyArrayObject*)centWeightObj;
    if(centWeightArr->nd==2){
      c->centWeightDim=2;
    }else if(centWeightArr->nd==4){
      c->centWeightDim=4;
    }else{
      printf("Error: centmodule - centWeight not 2d or 4d\n");
      return -1;
    }
    if(checkFloatContigArr(centWeightArr)!=0){
      printf("Error: centmodule - centWeight array not float or contig\n");
      return -1;
    }
    c->centWeightSize=centWeightArr->dimensions[c->centWeightDim-2];
    if(centWeightArr->dimensions[c->centWeightDim-2]!=c->corrsize || centWeightArr->dimensions[c->centWeightDim-1]!=c->corrsize){
      printf("Error: centmodule - centweight array not correct dimensions (need to be nsubx,nsubx,corrsize,corrsize\n");
      return -1;
    }
    c->centWeight=(float*)centWeightArr->data;
  }
  return 0;
}

int setSpotPsf(centstruct *c,PyObject *spotpsfObj){
  PyArrayObject *spotpsfArr;
  c->spotpsfDim=0;
  c->spotpsf=NULL;
  c->psfsize=c->clipsize;
  c->preBinningFactor=1;//no pre binning if no psf.
  c->paddedsize=c->clipsize;//After initial fft, psf gets resized (padded or clipped) to paddedsize (which should be a multiple of preBinningFactor)
  if(PyArray_Check(spotpsfObj)==1){
    spotpsfArr=(PyArrayObject*)spotpsfObj;
    if(spotpsfArr->nd==2){
      c->preBinningFactor=c->preBinningFactorOrig;
      c->spotpsfDim=2;
      c->psfsize=spotpsfArr->dimensions[0];
    }else if(spotpsfArr->nd==4){
      c->preBinningFactor=c->preBinningFactorOrig;
      c->spotpsfDim=4;
      c->psfsize=spotpsfArr->dimensions[2];
      //printf("setSpotPsf got psfsize %d\n",c->psfsize);
    }else{
      printf("Error: centmodule - spotpsf not 2d or 4d\n");
      return -1;
    }
    if(c->preBinningFactor<=1){//no prebinning
      c->paddedsize=c->psfsize;
    }else{
      c->paddedsize=c->fftsize;
    }
    if(c->paddedsize%c->preBinningFactor!=0){
      printf("Error: paddedsize(%d) modulo preBinningFactor(%d) should be zero\n",c->paddedsize,c->preBinningFactor);
      return -1;
    }
    if(c->preBinningFactor>1 && c->psfsize>c->paddedsize){
      printf("Error: centmodule - spotpsf dimensions should be equal to or less than paddedsize (%d)\n",c->paddedsize);
      return -1;
    }
    if(c->psfsize<c->clipsize){
      printf("Error: centmodule - spotpsf dimensions should be at least clipsize\n");
      return -1;
    }
    c->spotpsf=(float*)spotpsfArr->data;
    if(checkFloatContigArr(spotpsfArr)!=0){
      printf("Error: centmodule - spotpsf array not float or contiguous\n");
      return -1;
    }
  }
  return 0;
}

int prepareSpotPsf(centstruct *c){
  //int clipsize=c->clipsize;
  int psfsize=c->psfsize;
  float *tmparr;
  int i;
  if(c->fftPsf!=NULL)
    fftwf_free(c->fftPsf);
  c->fftPsf=NULL;
  //printf("psfsize, fsize %d %d, %p %d %d\n",c->psfsize,c->fsize,c->fftArrays,c->spotpsfDim,c->nthreads);
  if(c->psfsize>c->fsize){
    //need to reallocate fftArrays...
    c->fsize=c->fftsize>c->psfsize?c->fftsize:c->psfsize;
    //c->fsize=c->psfsize;
    if(c->fftArraySize!=c->fsize){
      for(i=0; i<c->nthreads; i++){
	if(c->fftArrays!=NULL){
	  if(c->fftArrays[i]!=NULL)
	    fftwf_free(c->fftArrays[i]);
	  c->fftArrays[i]=fftwf_malloc(sizeof(complex float)*c->fsize*c->fsize);
	  if(c->fftArrays[i]==NULL)
	    printf("Oops - malloc failed in prepareSpotPsf - this could be catastrophic...\n");
	}
      }
      c->fftArraySize=c->fsize;
    }
  }
  if(c->hllSize!=c->paddedsize && c->hll!=NULL){
    if(c->paddedsize<c->psfsize)
      printf("Oops - this could be bad... in centmodule.c - paddedsize<psfsize\n");
    for(i=0;i<c->nthreads;i++){
      if(c->hll[i]!=NULL)
	free(c->hll[i]);
      c->hll[i]=fftwf_malloc(sizeof(float)*c->paddedsize*c->paddedsize);
      if(c->hll[i]==NULL)
	printf("Oops - malloc failed for hll - this could be catastrophic... (centmodule.c\n");
    }
    c->hllSize=c->paddedsize;
  }
  if(c->planSize!=c->psfsize){
    fftwf_destroy_plan(c->rcfftplan);
    fftwf_destroy_plan(c->crfftplan);
    c->rcfftplan=fftwf_plan_dft_r2c_2d(c->psfsize,c->psfsize,c->hll[0],c->fftArrays[0],FFTW_MEASURE);
    c->crfftplan=fftwf_plan_dft_c2r_2d(c->psfsize,c->psfsize,c->fftArrays[0],c->hll[0],FFTW_MEASURE);
    c->planSize=c->psfsize;
  }
  if(c->spotpsfDim==2){
    c->fftPsf=fftwf_malloc(sizeof(complex float)*psfsize*(psfsize/2+1));
    if(c->fftPsf==NULL){
      printf("Failed to allocate memory for fftPsf\n");
      return -1;
    }
    memcpy(c->hll[0],c->spotpsf,sizeof(float)*psfsize*psfsize);
    flipArray(psfsize,c->hll[0],1,NULL);
    fftwf_execute_dft_r2c(c->rcfftplan,c->hll[0],c->fftPsf);
  }else if(c->spotpsfDim==4){
    c->fftPsf=fftwf_malloc(sizeof(complex float)*psfsize*(psfsize/2+1)*c->nsubaps);
    if(c->fftPsf==NULL){
      printf("Failed to allocate memory for fftPsf\n");
      return -1;
    }
    if((tmparr=malloc(sizeof(float)*(psfsize/2+1)))==NULL){
      printf("Failed to allocate temporary memory for tmparr\n");
      return -1;
    }
    for(i=0; i<c->nsubaps; i++){
      memcpy(c->hll[0],&(c->spotpsf[i*psfsize*psfsize]),sizeof(float)*psfsize*psfsize);
      flipArray(psfsize,c->hll[0],1,tmparr);
      fftwf_execute_dft_r2c(c->rcfftplan,c->hll[0],&(c->fftPsf[i*psfsize*(psfsize/2+1)]));
    }
    free(tmparr);
  }
  return 0;
}

int setSigArr(centstruct *c,PyObject *sigObj){
  PyArrayObject *sigArr;
  if(PyArray_Check(sigObj)==1){
    if(checkFloatContigArr((PyArrayObject*)sigObj)!=0){
      printf("Error: centmodule - sig array must be float and contiguous\n");
      return -1;
    }
    sigArr=(PyArrayObject*)sigObj;
    if(sigArr->nd==1 && sigArr->dimensions[0]==c->nsubaps){
      if(c->sigArr==NULL){
	if((c->sigArr=malloc(sizeof(float)*c->nsubaps))==NULL){
	  printf("Failed to alloc sigarr\n");
	  return -1;
	}
      }
      memcpy(c->sigArr,sigArr->data,sizeof(float)*c->nsubaps);
    }else{
      printf("sigArr should be 1D, size %d\n",c->nsubaps);
      return -1;
    }
  }else{
    return 1;
  }
  return 0;
}


int setupThreads(centstruct *c,int nthreads){
  int i,subapCnt,subapsLeft;
  int fftsize=c->fftsize;
  //First free any previous memory...
  //printf("setting up threads\n");
  for(i=0; i<c->nthreads; i++){
    if(c->runinfo!=NULL && c->runinfo[i]!=NULL){
      if(c->runinfo[i]->corr!=NULL){
	fftwf_free(c->runinfo[i]->corr);
      }
      free(c->runinfo[i]);
    }
    if(c->fftArrays!=NULL && c->fftArrays[i]!=NULL)
      fftwf_free(c->fftArrays[i]);
    if(c->hll!=NULL && c->hll[i]!=NULL)
      fftwf_free(c->hll[i]);
    if(c->gslRand!=NULL && c->gslRand[i]!=NULL)
      gsl_rng_free(c->gslRand[i]);
    if(c->sortarr!=NULL && c->sortarr[i]!=NULL)
      free(c->sortarr[i]);
  }
  if(c->hll!=NULL)
    free(c->hll);
  if(c->runinfo!=NULL)
    free(c->runinfo);
  if(c->threadid!=NULL)
    free(c->threadid);
  if(c->subapStart!=NULL)
    free(c->subapStart);
  if(c->fftArrays!=NULL)
    free(c->fftArrays);
  if(c->gslRand!=NULL)
    free(c->gslRand);
  if(c->sortarr!=NULL)
    free(c->sortarr);
  c->hll=NULL;
  c->runinfo=NULL;
  c->threadid=NULL;
  c->subapStart=NULL;
  c->fftArrays=NULL;
  c->sortarr=NULL;
  //and now set up...
  c->nthreads=nthreads;
  if(nthreads>0){
    c->hll=malloc(sizeof(float*)*nthreads);
    c->runinfo=malloc(sizeof(centrunstruct*)*nthreads);
    c->threadid=malloc(sizeof(pthread_t)*nthreads);
    c->subapStart=malloc(sizeof(int)*(nthreads+1));
    c->fftArrays=malloc(sizeof(complex float*)*nthreads);
    c->sortarr=malloc(sizeof(float*)*nthreads);

    //c->fftplan=malloc(sizeof(fftw_plan)*nthreads);
    subapCnt=0;
    subapsLeft=c->nsubaps;
    for(i=0; i<nthreads; i++){
      c->sortarr[i]=malloc(sizeof(float)*c->imgpxls);
      c->runinfo[i]=malloc(sizeof(centrunstruct));
      c->runinfo[i]->centstr=c;
      c->runinfo[i]->n=i;
      c->runinfo[i]->corr=fftwf_malloc(sizeof(float)*c->corrsize*c->corrsize);
      c->subapStart[i]=subapCnt;
      subapCnt+=subapsLeft/(nthreads-i);
      subapsLeft-=subapsLeft/(nthreads-i);
      c->fsize=fftsize>c->psfsize?fftsize:c->psfsize;
      c->fftArraySize=c->fsize;
      c->fftArrays[i]=fftwf_malloc(sizeof(complex float)*c->fsize*c->fsize);
      c->hll[i]=fftwf_malloc(sizeof(float)*c->paddedsize*c->paddedsize);
      c->hllSize=c->paddedsize;
    }
    c->subapStart[nthreads]=c->nsubaps;
    //setup random number generator...
    c->gslRand=malloc(sizeof(gsl_rng*)*nthreads);
    for(i=0; i<nthreads; i++){
      c->gslRand[i]=gsl_rng_alloc(gsl_rng_mt19937);
      if(c->seed==0){
	//printf("centmodule: Initialising random seed with time+%d\n",i);
	gsl_rng_set(c->gslRand[i],(unsigned long)time(0)+i);
      }else{
	gsl_rng_set(c->gslRand[i],c->seed+i);
      }
    }


  }
  //printf("done setting up threads\n");
  return 0;
}

PyObject *py_update(PyObject *self,PyObject *args){
  //update a value in the struct.
  int code;
  long lval;
  //int i;
  double dval;
  centstruct *c;
  PyObject *obj;
  PyArrayObject *aobj;
  int rt;
  if(!PyArg_ParseTuple(args,"liO",&c,&code,&obj)){
    printf("centmodule: parsing parameters failed.\nUsage: centstruct object, code for value to be changed, new value\n");
    return NULL;
  }
  PyErr_Clear();
  switch(code){
    case CALSOURCE:
      lval=PyInt_AsLong(obj);
      if((lval==-1) && PyErr_Occurred()){
	printf("centmodule: Error extracting integer value for calsource\n");
	return NULL;
      }
      c->calsource=lval;
      break;
    case SIG:
      if((rt=setSigArr(c,obj))==1){//failed to set as array - try single value
	dval=PyFloat_AsDouble(obj);
	if(PyErr_Occurred()){
	  printf("centmodule: Error extracting float value for sig\n");
	  return NULL;
	}
	c->sig=dval;
      }else if(rt==-1)
	return NULL;//failed
      break;
    case ADDPOISSON:
      lval=PyInt_AsLong(obj);
      if((lval==-1) && PyErr_Occurred()){
	printf("centmodule: Error extracting integer value for addpoisson\n");
	return NULL;
      }
      c->addPoisson=lval;
      break;
    case NCEN:
      lval=PyInt_AsLong(obj);
      if((lval==-1) && PyErr_Occurred()){
	printf("centmodule: Error extracting integer value for ncen\n");
	return NULL;
      }
      c->ncen=lval;
      break;
    case READNOISE:
      dval=PyFloat_AsDouble(obj);
      if(PyErr_Occurred()){
	printf("centmodule: Error extracting float value for readnoise\n");
	return NULL;
      }
      c->readNoise=dval;
      break;
    case READBG:
      dval=PyFloat_AsDouble(obj);
      if(PyErr_Occurred()){
	printf("centmodule: Error extracting float value for readbg\n");
	return NULL;
      }
      c->readBg=dval;
      break;
    case NOISEFLOOR:
      dval=PyFloat_AsDouble(obj);
      if(PyErr_Occurred()){
	printf("centmodule: Error extracting float value for noisefloor\n");
	return NULL;
      }
      c->noiseFloor=dval;
      break;
    case SKYBRIGHTNESS:
      if(PyArray_Check(obj)){
	aobj=(PyArrayObject*)obj;
	if(checkFloatContigArr(aobj)!=0){
	  printf("Error - centmodule - changing skybrightness must be float contiguous array or float\n");
	  return NULL;
	}else{
	  if(aobj->nd!=3 || aobj->dimensions[0]!=c->nsubaps || aobj->dimensions[1]!=c->nimg || aobj->dimensions[2]!=c->nimg){
	    printf("Error - centmodule - skybrightness should be 3D array shape nsubaps,nimg,nimg\n");
	    return NULL;
	  }else{
	    c->skybrightness=0;
	    c->skybrightnessArr=(float*)(aobj->data);
	  }
	}
      }else{
	dval=PyFloat_AsDouble(obj);
	if(PyErr_Occurred()){
	  printf("centmodule: Error extracting float value for skybrightness\n");
	  return NULL;
	}
	c->skybrightness=(float)dval;
	c->skybrightnessArr=NULL;
      }
      break;
    case PXLPOWER:
      dval=PyFloat_AsDouble(obj);
      if(PyErr_Occurred()){
	printf("centmodule: Error extracting float value for pxlpower\n");
	return NULL;
      }
      c->pxlPower=dval;
      break;
    case SPOTPSF:
      if(setSpotPsf(c,obj)==-1)
	return NULL;
      if(prepareSpotPsf(c)==-1)
	return NULL;
      break;
    case NTHREADS:
      lval=PyInt_AsLong(obj);
      if((lval==-1) && PyErr_Occurred()){
	printf("centmodule: Error extracting int value for nthreads\n");
	return NULL;
      }
      setupThreads(c,(int)lval);
      break;
    case OPTICALBINNING:
      lval=PyInt_AsLong(obj);
      if(lval==-1 && PyErr_Occurred()){
	printf("centmodule: Error extracting int value for opticalbinning\n");
	return NULL;
      }
      c->opticalBinning=(int)lval;
      break;
    case CENTWEIGHT:
      if(setCentWeight(c,obj)==-1){
	printf("centmodule: Error in setCentWeight\n");
	return NULL;
      }
      break;
    case CORRELATIONCENTROIDING:
      lval=PyInt_AsLong(obj);
      if((lval==-1) && PyErr_Occurred()){
	printf("centmodule: Error extracting int value for correlationCentroiding\n");
	return NULL;
      }
      c->correlationCentroiding=(int)lval;
      if(c->correlationCentroiding)
	c->corrsize=c->corrPatternSize;
      else
	c->corrsize=c->nimg;
      if(c->centWeight!=NULL && c->centWeightSize!=c->corrsize){
	printf("WARNING: centWeightSize INCORRECT: Ignoring centWeight until updated.\n");
	c->centWeight=NULL;
      }
      break;
    case CORRTHRESH:
      dval=PyFloat_AsDouble(obj);
      if(PyErr_Occurred()){
	printf("centmodule: Error extracting float value for corrThresh\n");
	return NULL;
      }
      c->corrThresh=(float)dval;
      break;
    case CORRPATTERN:
      if(PyArray_Check(obj)){
	aobj=(PyArrayObject*)obj;
	if(checkFloatContigArr(aobj)!=0){
	  printf("Error - centmodule - changing corrPattern must be float contiguous array\n");
	  return NULL;
	}
	c->corrPattern=(float*)aobj->data;
	if(c->corrPatternSize!=aobj->dimensions[2]){
	  //Here, insist that the pattern doesn't change size - simply because then the corrimg doesn't change also.
	  printf("Error - centmodule - corrPattern has changed size - please don't do that (or recode)!\n");
	  return NULL;
	  /*
	  for(i=0;i<c->nthreads;i++){
	    if(c->runinfo[i]->corr!=NULL)
	      fftwf_free(c->runinfo[i]->corr);
	    c->runinfo[i]->corr=fftwf_malloc(sizeof(float)*c->corrPatternSize*c->corrPatternSize);
	  }
	  fftwf_plan_free(c->corrPlan);
	  fftwf_plan_free(c->invCorrPlan);
	  c->corrPlan=fftwf_plan_r2r_2d(c->corrPatternSize,c->corrPatternSize,c->runinfo[0]->corr,c->runinfo[0]->corr,FFTW_R2HC, FFTW_R2HC, FFTW_ESTIMATE);
	  c->invCorrPlan=fftwf_plan_r2r_2d(c->corrPatternSize,c->corrPatternSize,c->runinfo[0]->corr,c->runinfo[0]->corr,FFTW_HC2R, FFTW_HC2R, FFTW_ESTIMATE);
	  c->corrPatternSize=aobj->dimensions[2];
	  if(c->correlationCentroiding)
	    c->corrsize=c->corrPatternSize;
	  if(c->centWeight!=NULL && c->centWeightSize!=c->corrsize){
	    printf("WARNING: centWeightSize INCORRECT: Ignoring centWeight until updated.\n");
	    c->centWeight=NULL;
	  }
	  */
	}
      }else{
	printf("Warning - setting corrPattern to NULL\n");
	c->corrPattern=NULL;
      }
      break;
    case CALDATA:
      if(setCalData(c,obj)==-1){
	printf("centmodule: Error in setCalData\n");
	return NULL;
      }
      break;
    case REFCENTS:
      if(PyArray_Check(obj)){
	aobj=(PyArrayObject*)obj;
	if(checkFloatContigArr(aobj)!=0){
	  printf("Error - centmodule - refCents must be float contig array\n");
	  return NULL;
	}
	c->refCents=(float*)aobj->data;
      }else{
	printf("Warning - setting refCents to NULL\n");
	c->refCents=NULL;
      }
      break;
    case CALCOEFF:
      if(PyArray_Check(obj)){
	aobj=(PyArrayObject*)obj;
	if(setCalCoeff(c,aobj)==-1){
	  printf("centmodule: Error in setCalCoeff\n");
	  return NULL;
	}
      }else if(obj==Py_None){
	setCalCoeff(c,NULL);
      }else{
	printf("centmodule: Error in setCalCoeff\n");
	return NULL;
      }
      break;
    case USEBRIGHTEST:
      if(PyInt_Check(obj)){
	c->useBrightest=(int)PyInt_AsLong(obj);
	c->useBrightestArr=NULL;
      }else if(PyArray_Check(obj)){
	if(checkIntContigArr((PyArrayObject*)obj)!=0){
	  printf("Error: centmodule - useBrightest must be int or int contiguous array\n");
	  return NULL;
	}
	c->useBrightest=0;
	c->useBrightestArr=(int*)(((PyArrayObject*)obj)->data);
	if(((PyArrayObject*)obj)->nd!=1 || ((PyArrayObject*)obj)->dimensions[0]!=c->nsubaps){
	  printf("Error: centmodule - useBrightest array should be 1D of shape nsubaps (nsubx^2)\n");
	  return NULL;
	}
      }else{
	printf("centmodule: Error in set useBrightest\n");
	return NULL;
      }
      break;
    case INTEGSTEPS://should always be smaller than the value this is initialised with.
      if(PyInt_Check(obj)){
	c->nintegrations=(int)PyInt_AsLong(obj);
	if(c->nintegrations>c->maxIntegrations){
	  printf("centmodule: Error - nintegrations>%d (max integrations)\n",c->maxIntegrations);
	  return NULL;
	}
      }else{
	printf("centmodule: Error in set integsteps\n");
	return NULL;
      }
      break;
    default:
      printf("centmodule: Code %d not recognised, nothing updated\n",code);
      return NULL;
      break;
  }
  Py_INCREF(Py_None);
  return Py_None;
}

PyObject *py_initialise(PyObject *self,PyObject *args){
  int nthreads,nsubaps, ncen,fftsize,clipsize,nimg,phasesize,addPoisson,calsource;
  int nintegrations,correlationCentroiding;
  int i,j;
  size_t stacksize;
  int scope;
  unsigned long seed;
  float readnoise,readbg,noiseFloor,sig=0.,pxlPower,tmp,corrThresh;
  //float *skybrightnessArr=NULL;
  double dval;
  int opticalBinning,imageOnly,preBinningFactor;
  int gaussianFit,parabolicFit;
  float gaussianMinVal,gaussianReplaceVal;
  int threshType;
  float *sigArr=NULL;
  PyObject *spotpsfObj,*sigObj,*skybrightnessObj,*centWeightObj,*corrPatternObj,*corrimgObj,*useBrightestObj;
  PyArrayObject *phs,*pupfn,*cents,*fracSubArea,*subflag,*bimg,*sigArrObj,*aobj,*corrPattern=NULL,*corrimg=NULL;
  centstruct *c;
  if(!PyArg_ParseTuple(args,"iiiiiiiffifOOifilO!O!OO!O!O!O!iOifOOiiOiiiff",&nthreads,
		       &nsubaps,&ncen,&fftsize,&clipsize,
		       &nimg,&phasesize,&readnoise,&readbg,&addPoisson,&noiseFloor,
		       &sigObj,&skybrightnessObj,&calsource,
		       &pxlPower,&nintegrations,&seed,
		       &PyArray_Type,&phs,&PyArray_Type,
		       &pupfn,&spotpsfObj,&PyArray_Type,&cents,
		       &PyArray_Type,&subflag,&PyArray_Type,&bimg,&PyArray_Type,&fracSubArea,&opticalBinning,&centWeightObj,&correlationCentroiding,&corrThresh,&corrPatternObj,&corrimgObj,&threshType,&imageOnly,&useBrightestObj,&preBinningFactor,&parabolicFit,&gaussianFit,&gaussianMinVal,&gaussianReplaceVal)){
    printf("Usage: nthreads,nsubaps,ncen,fftsize,clipsize,nimg,phasesize,readnoise,readbg,addpoisson,noisefloor,sig,skybrightness,calsource,pxlpower,nintegrations,seed,phs,pupfn,spotpsf,cents,subflag,bimg,fracsubarea,opticalbinning,centWeight,correlationCentroiding,corrThresh,corrPattern,corrimg,threshType,imageOnly,useBrightest,preBinningFactor,parabolicFit,gaussianFit,gaussianMinVal,gaussianReplaceVal\n");
    return NULL;
  }
  if((c=malloc(sizeof(centstruct)))==NULL){
    printf("Failed to malloc centstruct\n");
    return NULL;
  }
  memset(c,0,sizeof(centstruct));
  c->clipsize=clipsize;
  c->nthreads=nthreads;
  c->nsubaps=nsubaps;
  c->ncen=ncen;
  c->fftsize=fftsize;
  c->nimg=nimg;
  c->phasesize=phasesize;
  c->preBinningFactor=preBinningFactor;
  c->preBinningFactorOrig=preBinningFactor;
  c->parabolicFit=parabolicFit;
  c->gaussianFit=gaussianFit;
  c->gaussianMinVal=gaussianMinVal;
  c->gaussianReplaceVal=gaussianReplaceVal;
    
  //Now check that the arrays are okay...
  if(setSpotPsf(c,spotpsfObj)==-1){
    return NULL;
  }
  if(checkFloatContigArr(cents)!=0){
    printf("Error: centmodule - cents must be float and contiguous\n");
    return NULL;
  }
  if(checkFloatContigArr(phs)!=0){
    printf("Error: centmodule - phs must be float and contiguous\n");
    return NULL;
  }
  j=1;
  for(i=0;i<phs->nd;i++)
    j*=phs->dimensions[i];
  if(j==phasesize*phasesize*nsubaps*nintegrations)
    c->phaseStep=1;
  else if(j==phasesize*phasesize*nsubaps*nintegrations*2)
    c->phaseStep=2;//phaseamp mode.
  else{
    printf("Error: centmodule - phs size is wrong\n");
    return NULL;
  }
    
  if(checkFloatContigArr(pupfn)!=0){
    printf("Error: centmodule - pupfn must be float and contiguous\n");
    return NULL;
  }
  if(checkIntContigArr(subflag)!=0){
    printf("Error: centmodule - subflag must be int32 and contiguous\n");
    return NULL;
  }
  if(checkFloatContigArr(bimg)!=0){
    printf("Error: centmodule - bimg must be float and contiguous\n");
    return NULL;
  }
  if(checkFloatContigArr(fracSubArea)!=0){
    printf("Error: centmodule - fracSubArea must be float and contiguous\n");
    return NULL;
  }

  if(PyFloat_Check(sigObj)){
    sig=(float)PyFloat_AsDouble(sigObj);
  }else if(PyArray_Check(sigObj)){
    if(checkFloatContigArr((PyArrayObject*)sigObj)!=0){
      printf("Error: centmodule - sig array must be float and contiguous\n");
      return NULL;
    }
    sigArrObj=(PyArrayObject*)sigObj;
    sigArr=(float*)(sigArrObj->data);
    if(sigArrObj->nd!=1 || sigArrObj->dimensions[0]!=nsubaps){
      printf("Error: centmodule - sig arr dimensions must be 1D (nsubx^2)\n");
      printf("sigArrObj nd=%d (nsubaps=%d)\n",sigArrObj->nd,nsubaps);
      for(i=0; i<sigArrObj->nd; i++){
	printf("%d %ld\n",i,(long)sigArrObj->dimensions[i]);
      }
      return NULL;
    }
  }else{
    printf("Error: centmodule - sig object must be float or float array\n");
    return NULL;
  }
  if(PyArray_Check(skybrightnessObj)){
    aobj=(PyArrayObject*)skybrightnessObj;
    if(checkFloatContigArr(aobj)!=0){
      printf("Error - centmodule - skybrightness must be float contiguous array or float\n");
      return NULL;
    }else{
      if(aobj->nd!=3 || aobj->dimensions[0]!=nsubaps || aobj->dimensions[1]!=nimg || aobj->dimensions[2]!=nimg){
	printf("Error - centmodule - skybrightness should be 3D array shape nsubaps,nimg,nimg\n");
	return NULL;
      }else{
	c->skybrightness=0;
	c->skybrightnessArr=(float*)(aobj->data);
      }
    }
  }else{
    dval=PyFloat_AsDouble(skybrightnessObj);
    if(PyErr_Occurred()){
      printf("centmodule: Error extracting float value for skybrightness\n");
      return NULL;
    }
    c->skybrightness=(float)dval;
    c->skybrightnessArr=NULL;
  }
  //printf("checking corr arrays\n");
  if(PyArray_Check(corrPatternObj)){
    //printf("corrPattern is array\n");
    corrPattern=(PyArrayObject*)corrPatternObj;
    c->corrPattern=(float*)corrPattern->data;
    //printf("got corrPattern array\n");
    c->corrsize=corrPattern->dimensions[2];
    c->corrPatternSize=corrPattern->dimensions[2];//save it incase switch in and out of correlation mode.
    if(checkFloatContigArr(corrPattern)!=0){
      printf("Error: centmodule - corrPattern must be float and contiguous\n");
      return NULL;
    }
  }
  if(PyArray_Check(corrimgObj)){
    //printf("corrimg is array\n");
    corrimg=(PyArrayObject*)corrimgObj;
    c->corrimg=(float*)corrimg->data;
    if(checkFloatContigArr(corrimg)!=0){
      printf("Error: centmodule - corrimg must be float and contiguous\n");
      return NULL;
    }
  }
  if(PyInt_Check(useBrightestObj)){
    c->useBrightest=(int)PyInt_AsLong(useBrightestObj);
    c->useBrightestArr=NULL;
  }else if(PyArray_Check(useBrightestObj)){
    if(checkIntContigArr((PyArrayObject*)useBrightestObj)!=0){
      printf("Error: centmodule - useBrightest must be int or int contiguous array\n");
      return NULL;
    }
    c->useBrightest=0;
    c->useBrightestArr=(int*)(((PyArrayObject*)useBrightestObj)->data);
    if(((PyArrayObject*)useBrightestObj)->nd!=1 || ((PyArrayObject*)useBrightestObj)->dimensions[0]!=nsubaps){
      printf("Error: centmodule - useBrightest array should be 1D of shape nsubaps (nsubx^2)\n");
      return NULL;
    }
  }else{
    printf("Error: centmodule - useBrightest must be int or int array\n");
    return NULL;
  }
    

  c->readNoise=readnoise;
  c->readBg=readbg;
  c->threshType=threshType;
  c->addPoisson=addPoisson;
  c->noiseFloor=noiseFloor;
  c->sig=sig;
  c->sigArr=sigArr;
  //c->skybrightness=skybrightness;
  c->calsource=calsource;
  c->pxlPower=pxlPower;
  c->nintegrations=nintegrations;
  c->maxIntegrations=nintegrations;
  c->seed=seed;
  c->phs=(float*)phs->data;//can include amplitude data too (interleaved, phase first).
  c->pupfn=(float*)pupfn->data;
  c->cents=(float*)cents->data;
  c->subflag=(int*)subflag->data;
  c->bimg=(float*)bimg->data;
  c->fracSubArea=(float*)fracSubArea->data;
  c->fftpxls=fftsize*fftsize;
  c->phasepxls=phasesize*phasesize;
  c->imgpxls=nimg*nimg;
  c->opticalBinning=opticalBinning;
  c->imageOnly=imageOnly;
  c->corrThresh=corrThresh;
  //ns=c->nsubaps/c->nthreads;//number of subaps per thread (roughly).
  //nx=c->nsubaps-ns*c->nthreads;//the number that didn't get divided correctly...

  c->tiltfn=malloc(sizeof(float)*phasesize*phasesize);
  memset(c->tiltfn,0,sizeof(float)*phasesize*phasesize);
  tmp=M_PI/fftsize*(fftsize+1.);//was M_PI/fftsize and a fliparray.
  for(i=0; i<phasesize; i++){
    for(j=0; j<phasesize; j++){
      c->tiltfn[i*phasesize+j]=tmp*(i+j+1-phasesize);
    }
  }
  c->correlationCentroiding=correlationCentroiding;
  if(correlationCentroiding){
    c->corrsize=c->corrPatternSize;
  }else{
    c->corrsize=c->nimg;
  }
  if(setCentWeight(c,centWeightObj)==-1)
    return NULL;

  setupThreads(c,nthreads);
  //setup plans for correlation...
  c->corrPlan=fftwf_plan_r2r_2d(c->corrsize,c->corrsize,c->runinfo[0]->corr,c->runinfo[0]->corr,FFTW_R2HC, FFTW_R2HC, FFTW_ESTIMATE);
  c->invCorrPlan=fftwf_plan_r2r_2d(c->corrsize,c->corrsize,c->runinfo[0]->corr,c->runinfo[0]->corr,FFTW_HC2R, FFTW_HC2R, FFTW_ESTIMATE);
  c->pthread_attr=malloc(sizeof(pthread_attr_t));
  pthread_attr_init(c->pthread_attr);
  pthread_attr_getstacksize(c->pthread_attr,&stacksize);
  pthread_attr_setscope(c->pthread_attr,PTHREAD_SCOPE_SYSTEM);
  pthread_attr_getscope(c->pthread_attr,&scope);
  //printf("centmodle - Stacksize: %ld %d %d\n",stacksize,scope,PTHREAD_SCOPE_SYSTEM);
  stacksize/=8;
  pthread_attr_setstacksize(c->pthread_attr,stacksize);
  //we can use the same plan for each thread, since they are thread safe.  And use 
  //fftw_execute_dft for on different arrays for thread safty.
  c->fftplan=fftwf_plan_dft_2d(fftsize,fftsize,c->fftArrays[0], c->fftArrays[0],
			      FFTW_FORWARD,FFTW_MEASURE);//use fftw_execute_dft(fftplan,in,out) to thread exec
  //c->rcfftplan=fftwf_plan_dft_r2c_2d(c->psfsize,c->psfsize,c->hll[0],c->fftArrays[0],FFTW_MEASURE);
  //c->crfftplan=fftwf_plan_dft_c2r_2d(c->psfsize,c->psfsize,c->fftArrays[0],c->hll[0],FFTW_MEASURE);

  //now perform ffts for the spot PSFs if necessary...
  if(prepareSpotPsf(c)==-1)
    return NULL;
  c->binthispxl=NULL;
  c->binnextpxl=NULL;
  c->binindx=NULL;
  prepareBinImage(clipsize,nimg,c);
  //printf("centmodule: finished initialise\n");
  return Py_BuildValue("l",(long)c);//return a pointer to centstruct.
}

PyObject *py_free(PyObject *self,PyObject *args){
  centstruct *c;
  if(!PyArg_ParseTuple(args,"l",&c)){
    printf("Usage: centstruct\n");
    return NULL;
  }
  

  free(c->tiltfn);
  fftwf_destroy_plan(c->corrPlan);
  fftwf_destroy_plan(c->invCorrPlan);
  fftwf_destroy_plan(c->fftplan);
  fftwf_destroy_plan(c->rcfftplan);
  fftwf_destroy_plan(c->crfftplan);
  setupThreads(c,0);
  c->spotpsfDim=0;
  prepareSpotPsf(c);
  prepareBinImage(0,0,c);
  if(c->fftPsf!=NULL)
    fftwf_free(c->fftPsf);
  free(c);
  Py_INCREF(Py_None);
  return Py_None;
}

PyObject *py_run(PyObject *self,PyObject *args){
  centstruct *c;
  centrunstruct **runinfo;
  pthread_t *thread;
  int nthreads,i;
  float dtime=0.;
  struct timeval t1,t2;
  //pid_t pid;
  //pid_t *pidlist;
  //pidlist=malloc(sizeof(pid_t)*nthreads);
  if(!PyArg_ParseTuple(args,"l",&c)){
    printf("Usage: centstruct\n");
    return NULL;
  }
  //memset(c->bimg,0,sizeof(float)*c->nsubaps*c->imgpxls);
  gettimeofday(&t1,NULL);
  nthreads=c->nthreads;
  runinfo=c->runinfo;//malloc(sizeof(centrunstruct)*nthreads);
  thread=c->threadid;//malloc(sizeof(pthread_t)*nthreads);
  //first create the threads and set them running, and then run ourself.
  for(i=1; i<nthreads; i++){
    pthread_create(&thread[i],c->pthread_attr,(void*)centroidsFromPhase,runinfo[i]);
  }
  //runinfo[0]->centstr=c;
  //runinfo[0]->n=0;
  centroidsFromPhase(runinfo[0]);//run ourself...
  //wait for thread completion.
  for(i=1; i<nthreads; i++){
    pthread_join(thread[i],NULL);
  }
  gettimeofday(&t2,NULL);
  dtime=(t2.tv_sec*1000000+t2.tv_usec-t1.tv_sec*1000000-t1.tv_usec)/1e6;
  //free(runinfo);
  //free(thread);
  return Py_BuildValue("f",dtime);
}

PyObject *py_testcorrelation(PyObject *self,PyObject *args){
  fftwf_plan corrPlan;
  fftwf_plan invCorrPlan;
  float *tmp;
  PyArrayObject *img,*corr;
  int nimg;
  //Note - corr must already be shifted and fft'd in half complex format.
  //Note, corr gets replaced by the correlated image.
  if(!PyArg_ParseTuple(args,"O!O!",&PyArray_Type,&img,&PyArray_Type,&corr)){
    printf("usage: img, corr\n");
    return NULL;
  }
  if(checkFloatContigArr(img)!=0){
    printf("Error: centmodule - img must be float and contiguous\n");
    return NULL;
  }
  if(checkFloatContigArr(corr)!=0){
    printf("Error: centmodule - corr must be float and contiguous\n");
    return NULL;
  }
  if(img->nd!=2 || corr->nd!=2 || img->dimensions[0]!=img->dimensions[1] || corr->dimensions[0]!=img->dimensions[0] || corr->dimensions[1]!=img->dimensions[0]){
    printf("Error: centmodule - corr/img wrong shape\n");
    return NULL;
  }
  nimg=img->dimensions[0];

  tmp=fftwf_malloc(sizeof(float)*nimg*nimg);
  corrPlan=fftwf_plan_r2r_2d(nimg,nimg,tmp,tmp,FFTW_R2HC, FFTW_R2HC, FFTW_ESTIMATE);
  invCorrPlan=fftwf_plan_r2r_2d(nimg,nimg,tmp,tmp,FFTW_HC2R, FFTW_HC2R, FFTW_ESTIMATE);
  calcCorrelation(nimg,nimg,(float*)corr->data,(float*)img->data,(float*)corr->data,tmp,corrPlan,invCorrPlan);
  fftwf_free(tmp);
  fftwf_destroy_plan(corrPlan);
  fftwf_destroy_plan(invCorrPlan);
  return Py_BuildValue("");
}



static PyMethodDef centMethods[] = {
  {"run",  py_run, METH_VARARGS,"Run the centroid algorithm."},
  {"initialise",py_initialise,METH_VARARGS,"Initialise the cent module"},
  {"update",py_update,METH_VARARGS,"Update values used for centroiding"},
  {"free",py_free,METH_VARARGS,"Free the arrays allocated (invalidates the centstruct pointer)"},
  {"testcorrelation",py_testcorrelation,METH_VARARGS,"Test correlation"},
  {NULL, NULL, 0, NULL}        /* Sentinel */
};
//PyMODINIT_FUNC 
void initcent(void)
{
  PyObject *m;
  PyImport_AddModule("cent");
  m=Py_InitModule("cent", centMethods);
  import_array();
  CentError = PyErr_NewException("cent.error", NULL, NULL);
  Py_INCREF(CentError);
  PyModule_AddObject(m, "error", CentError);
}


int
main(int argc, char *argv[])
{
    // Pass argv[0] to the Python interpreter 
  //printf("main\n");
  Py_SetProgramName(argv[0]);
  
  // Initialize the Python interpreter.  Required. 
  Py_Initialize();
  
  // Add a static module 
  initcent();
  return 0;
}

