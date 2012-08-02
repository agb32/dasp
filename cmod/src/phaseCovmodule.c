/* 
 *   Python extension module to calculate subaperture tilt covariances 
 *   for non-Kolmogorov turbulence (either Boreman-Dainty or Von Karman)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>

#include "Python.h"
#include "numpy/arrayobject.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_gamma.h>

#include "nr.h"
//#include "nrutil.h"


/* Kolmogorov structure function
   (r, r_0) */
double Kolstrucfunc(double x, double r0)
{
  return 6.88*pow(x/r0,5./3.);
}



/* Boreman Dainty structure function (needs multiplying by gamma_beta elsewhere) 
   (r, rho_0, beta) */
double BDstrucfunc(double x, double rho, double beta)
{
  return pow(x/rho,beta-2.);
}



/* function to evaluate gamma_beta (Rao et al. 2000) */
double gammabeta(double beta)
{
  double g0,g2,g4,gz;

  g0=gsl_sf_gamma(beta/2.);
  g2=gsl_sf_gamma((beta+2.)/2.);
  g4=gsl_sf_gamma((beta+4.)/2.);
  gz=gsl_sf_gamma(beta+1.);

  return pow(2.,beta-1.)*g2*g2*g4/(g0*gz);
}



/* Von Karman structure function (Jenkins 1998) 
   (r, r0, mu) */
double VKstrucfunc(double x, double r0, double mu)
{
  double A,C,K,brack;

  if (x < 1e-6)
    return 0.;
  else {
    A=2.*pow(M_PI*x/mu,5./6.)/gsl_sf_gamma(5./6.);
    //B=pow(x/mu,5./6.);
    K=gsl_sf_bessel_Knu(5./6.,2.*M_PI*x/mu);
    brack=1.-A*K;
    C=pow(mu/r0,5./3.);
    return 0.17253*C*brack;
  }
}
/*
static PyObject *kolmogorovArray(self,args)
     PyObject *self, *args;
{//compute array with kolmogorov values stored in 't.
  PyArrayObject	*phaseCov,*modes,*output;

  int      g,h,i,j,k,l,nmodes,nx,ny,odi,odj,pcdi,pcdj,mdi,mdj,mdk;
  float partint;

  //inputs needed:
  //phiphi, the phase covariance matrix, shape (y,x)
  //the mirror modes, a 3d numeric array with shape (nmodes, y,x)


  if (!PyArg_ParseTuple(args, "O!ddO!", &PyArray_Type ,&pupilfn,&scale,&r0,&PyArray_Type,&output)){
      printf("Inputs should be phase covariance (2d array) and modes (3d array) and output (2d array)\n");
      return NULL;
  }
}
*/
/*********** THEORY *****************
let a_i=integ( phs(r) Z_i(r) dr)
where Z_i is zernike mode i and phs is the phase screen.
So a_i is the zernike coefficients for the phase.
Then, the Noll matrix is M_ij = < a_i a_j > where <> means average over time.

For more general modes, M_ij = < a_i a_j > - <a_i><a_j>

Order of integrals can be changed.

So (letting x,y be dummy variables instead of r:
 int_t(int_x(p(x,t)Z_i(x)dx) int_y(p(y,t)Z_j(y)dy) dt)
=int_y(int_x(int_t(p(x,t)p(y,t)dt)Z_i(x)dx)Z_j(y)dy)

Let f(x,y)=int_t(p(x,t)p(y,t)dt)
Then, since shift invariant, let:
f(x-y)=int_t(p(x-y,t),p(0,t)dt)
and let r=x-y, so dr=dx when y is constant

So, in the covariance function, f is phaseCov, and Z is modes.
Then, when doing the integral, we do:
int_y(int_r(f(r)Z_i(r-y)dr) Z_j(y)dy)
which is what the function below does.

Scaling:  The integrals are averages not sums.  If using zernikes, they should be normalised to one.
The phaseCov (f) array should also be averaged.


************************************/
typedef struct{
  PyArrayObject *phaseCov;
  PyArrayObject *modes;
  PyArrayObject *modeCoords;
  PyArrayObject *vig;
  PyArrayObject *output;
  int start;
  int end;
  int msg;//0 or equal to end-start if this thread should print...
  int nact;
  int ny;
  int nx;
  float *coords;
  float *xin;
  int di;
  float *pc;
  int pcdi;
  int pcdj;
  float *out;
  int odi;
  int odj;
  float *mode;
}runStruct;

int covWorker(runStruct *runInfo){
  PyArrayObject *phaseCov,*modes,*output;
  int start,end,msg;
  int nmodes,ny,nx,odi,odj,pcdi,pcdj,mdi,mdj,mdk,g,h,i,k,j,l;
  float partint;
  phaseCov=runInfo->phaseCov;
  modes=runInfo->modes;
  output=runInfo->output;
  start=runInfo->start;
  end=runInfo->end;
  msg=runInfo->msg;
  nmodes=modes->dimensions[0];
  ny=modes->dimensions[1];
  nx=modes->dimensions[2];
  odi=output->strides[0];
  odj=output->strides[1];
  pcdi=phaseCov->strides[0];
  pcdj=phaseCov->strides[1];
  mdi=modes->strides[0];
  mdj=modes->strides[1];
  mdk=modes->strides[2];
  for(g=start; g<end; g++){
    //sym h=g. 
    for(h=g; h<nmodes; h++){
      if(msg){
	printf("Creating phase covariance: %d %d / %d %d \r",g,h,msg,nmodes);
	fflush(NULL);
      }
      (*((float*)(output->data+g*odi+h*odj)))=0.;
      for(i=0; i<ny; i++){
	for(j=0; j<nx; j++){
	  partint=0.;
	  for(k=0; k<ny; k++){
	    for(l=0; l<nx; l++){
	      partint+=(*((float*)(phaseCov->data+abs(k-i)*pcdi+abs(l-j)*pcdj))) * (*((float*)(modes->data+h*mdi+k*mdj+l*mdk)));
	    }
	  }
	  (*((float*)(output->data+g*odi+h*odj)))+=partint*(*((float*)(modes->data+g*mdi+i*mdj+j*mdk)));
	}
      }
      //symmetric, so copy the transpose part.
      (*((float*)(output->data+h*odi+g*odj)))=(*((float*)(output->data+g*odi+h*odj)));
      //and sym about other diagnal too...
      //(*((float*)(output->data+(nmodes-h-1)*odi+(nmodes-g-1)*odj)))=(*((float*)(output->data+g*odi+h*odj)));
      //(*((float*)(output->data+(nmodes-g-1)*odi+(nmodes-h-1)*odj)))=(*((float*)(output->data+g*odi+h*odj)));
    }
  }
  return 0;
}



static PyObject *covariance(PyObject *self,PyObject *args)
{
  PyArrayObject	*phaseCov,*modes,*output;

  int nmodes,nx,ny,odi,odj,pcdi,pcdj,mdi,mdj,mdk;
  int t;
  //float partint;
  int nthreads=1;
  pthread_t *thread;
  runStruct *threadInfo;
  int NperThread2;
  int s,e;
  //float nPerThread;
  //inputs needed:
  //phiphi, the phase covariance matrix, shape (y,x)
  //the mirror modes, a 3d numeric array with shape (nmodes, y,x)

  //If you change this function, then please also increment the version 
  //number in aosim/util/phaseCoveriance.py

  if (!PyArg_ParseTuple(args, "O!O!O!|i", &PyArray_Type ,&phaseCov,&PyArray_Type ,&modes,&PyArray_Type,&output,&nthreads)){
      printf("Inputs should be phase covariance (2d array) and modes (3d array) and output (2d array), nthreads (optional)\n");
      return NULL;
  }
  if(phaseCov->nd!=2 || modes->nd!=3 || output->nd!=2){
      printf("Dimensions for inputs should be 2,3,2\n");
      return NULL;
  }
  if(phaseCov->descr->type_num!=NPY_FLOAT ||
     modes->descr->type_num!=NPY_FLOAT ||
     output->descr->type_num!=NPY_FLOAT){
      printf("Arrays should be 32 bit floating point\n");
      return NULL;
  }
  nmodes=modes->dimensions[0];
  ny=modes->dimensions[1];
  nx=modes->dimensions[2];
  odi=output->strides[0];
  odj=output->strides[1];
  pcdi=phaseCov->strides[0];
  pcdj=phaseCov->strides[1];
  mdi=modes->strides[0];
  mdj=modes->strides[1];
  mdk=modes->strides[2];
  if(phaseCov->dimensions[0]!=ny || phaseCov->dimensions[1]!=nx || output->dimensions[0]!=nmodes || output->dimensions[1]!=nmodes){
      printf("Shape of phaseCov or output is wrong\n");
      return NULL;
  }
  Py_BEGIN_ALLOW_THREADS;
  thread=(pthread_t*)malloc(sizeof(pthread_t)*nthreads);
  threadInfo=(runStruct*)malloc(sizeof(runStruct)*nthreads);
  //nleft=nmodes;
  NperThread2=(nmodes*nmodes)/nthreads;//twice the number of points evaluated per thread.
  //nPerThread=(nmodes*nmodes)/4./nthreads;//no of points evaluated per thread.
  threadInfo[0].phaseCov=phaseCov;
  threadInfo[0].modes=modes;
  threadInfo[0].output=output;
  threadInfo[0].start=0;
  s=threadInfo[0].start;
  //threadInfo[0].end=nleft/nthreads;
  //threadInfo[0].end=(int)(nmodes+0.5-sqrt(nmodes*nmodes-NperThread2+s*s-2*nmodes*s));//nleft/nthreads;
  //e=(int)(0.5+(nmodes-sqrt(nmodes*nmodes-4*(nPerThread+s*nmodes-s*s)))/2);
  e=(int)(nmodes+0.5-sqrt(nmodes*nmodes-NperThread2+s*s-2*nmodes*s));
  if(nthreads==1)
    //e=(nmodes+1)/2;//just in case of floating point round errors.
    e=nmodes;
  threadInfo[0].end=e;
  //printf("phasecov %d (%d elements)\n",e,(e-s)*(nmodes-e-s));
  printf("phasecov %d (%g elements)\n",e,(2*nmodes-s-e)/2.*(e-s));
  threadInfo[0].msg=e-s;
  //nleft=nmodes-threadInfo[0].end;
  for(t=1; t<nthreads; t++){
    threadInfo[t].phaseCov=phaseCov;
    threadInfo[t].modes=modes;
    threadInfo[t].output=output;
    threadInfo[t].start=e;//threadInfo[t-1].end;
    //threadInfo[t].end=threadInfo[t].start+nleft/(nthreads-t);
    s=e;//threadInfo[t].start;
    //e=(int)(0.5+(nmodes-sqrt(nmodes*nmodes-4*(nPerThread+s*nmodes-s*s)))/2);
    e=(int)(nmodes+0.5-sqrt(nmodes*nmodes-NperThread2+s*s-2*nmodes*s));
    if(t==nthreads-1)
      //e=(nmodes+1)/2;
      e=nmodes;
    threadInfo[t].end=e;//(int)(nmodes+0.5-sqrt(nmodes*nmodes-NperThread2+s*s-2*nmodes*s));
    //nleft=nmodes-threadInfo[t].end;
    threadInfo[t].msg=0;
    printf("phasecov %d (%g elements)\n",e,(2*nmodes-s-e)/2.*(e-s));
    
    //printf("phasecov %d (%g elements)\n",threadInfo[t].end,(2*nmodes-threadInfo[t].start-threadInfo[t].end)/2.*(threadInfo[t].end-threadInfo[t].start));
    pthread_create(&thread[t],NULL,(void*)covWorker,&threadInfo[t]);
  }

  covWorker(&threadInfo[0]);
  for(t=1; t<nthreads; t++){
    pthread_join(thread[t],NULL);
  }
  //printf("freeing\n");
  free(thread);
  free(threadInfo);
  //printf("freed\n");

/*  for(g=0; g<nmodes; g++){
    for(h=0; h<nmodes; h++){
      printf("Creating phase covariance: %d %d / %d  \r",g,h,nmodes);
      fflush(NULL);
      (*((float*)(output->data+g*odi+h*odj)))=0.;
      for(i=0; i<ny; i++){
	for(j=0; j<nx; j++){
	  partint=0.;
	  for(k=0; k<ny; k++){
	    for(l=0; l<nx; l++){
	      partint+=(*((float*)(phaseCov->data+abs(k-i)*pcdi+abs(l-j)*pcdj))) * (*((float*)(modes->data+h*mdi+k*mdj+l*mdk)));
	    }
	  }
	  (*((float*)(output->data+g*odi+h*odj)))+=partint*(*((float*)(modes->data+g*mdi+i*mdj+j*mdk)));
	}
      }
    }
    }*/
  printf("Creating phase covariance: Done               \n");
  Py_END_ALLOW_THREADS;
  return Py_BuildValue("");
}


int covWorkerLocal(runStruct *runInfo){
  PyArrayObject *phaseCov,*modes,*modeCoords,*output,*vig;
  int start,end,msg;
  int nmodes,ny,nx,odi,odj,pcdi,pcdj,mdi,mdj,mdk,g,h,i,k,j,l,vdi;
  int goffsetx,goffsety,hoffsetx,hoffsety;
  int mcdi,mcdj;
  float partint;
  phaseCov=runInfo->phaseCov;
  modes=runInfo->modes;
  vig=runInfo->vig;
  if(vig==NULL)
    vdi=0;
  else
    vdi=vig->strides[0];
  modeCoords=runInfo->modeCoords;
  mcdi=modeCoords->strides[0];
  mcdj=modeCoords->strides[1];
  output=runInfo->output;
  start=runInfo->start;
  end=runInfo->end;
  msg=runInfo->msg;
  nmodes=modes->dimensions[0];
  ny=modes->dimensions[1];
  nx=modes->dimensions[2];
  odi=output->strides[0];
  odj=output->strides[1];
  pcdi=phaseCov->strides[0];
  pcdj=phaseCov->strides[1];
  mdi=modes->strides[0];
  mdj=modes->strides[1];
  mdk=modes->strides[2];
  for(g=start; g<end; g++){
    //sym h=g. 
    goffsetx=*((int*)(modeCoords->data+g*mcdi));
    goffsety=*((int*)(modeCoords->data+g*mcdi+mcdj));
    for(h=g; h<nmodes; h++){
      hoffsetx=*((int*)(modeCoords->data+h*mcdi));
      hoffsety=*((int*)(modeCoords->data+h*mcdi+mcdj));
      if(msg){
	printf("Creating phase covariance: %d %d / %d %d \r",g,h,msg,nmodes);
	fflush(NULL);
      }
      if(vig==NULL || *((int*)(vig->data+g*vdi))==1 || *((int*)(vig->data+h*vdi))==1){
	//a vignetted mode... or not caring about vignetting
	(*((float*)(output->data+g*odi+h*odj)))=0.;
	for(i=0; i<ny; i++){
	  for(j=0; j<nx; j++){
	    partint=0.;
	    for(k=0; k<ny; k++){
	      for(l=0; l<nx; l++){
		partint+=(*((float*)(phaseCov->data+abs(k+hoffsety-i-goffsety)*pcdi+abs(l+hoffsetx-j-goffsetx)*pcdj))) * (*((float*)(modes->data+h*mdi+k*mdj+l*mdk)));
	      }
	    }
	    (*((float*)(output->data+g*odi+h*odj)))+=partint*(*((float*)(modes->data+g*mdi+i*mdj+j*mdk)));
	  }
	}
	//symmetric, so copy the transpose part.
	(*((float*)(output->data+h*odi+g*odj)))=(*((float*)(output->data+g*odi+h*odj)));
      }
    }
  }
  return 0;
}

static PyObject *covarianceLocal(PyObject *self,PyObject *args)
{
  //A version of above that uses modes that are localised... scales much better for large pupils.
  PyArrayObject	*phaseCov,*modes,*output,*modeCoords,*vig;
  PyObject *vigObj;
  int nmodes,nx,ny;
  //int odi,odj,pcdi,pcdj,mdi,mdj,mdk;
  int t;
  //float partint;
  int nthreads=1;
  pthread_t *thread;
  runStruct *threadInfo;
  int NperThread2;
  int s,e;
  //float nPerThread;
  //inputs needed:
  //phiphi, the phase covariance matrix, shape (y,x)
  //the mirror modes, a 3d numeric array with shape (nmodes, y,x)

  //If you change this function, then please also increment the version 
  //number in aosim/util/phaseCoveriance.py

  if (!PyArg_ParseTuple(args, "O!O!O!OO!|i", &PyArray_Type ,&phaseCov,&PyArray_Type ,&modes,&PyArray_Type,&modeCoords,&vigObj,&PyArray_Type,&output,&nthreads)){
      printf("Inputs should be phase covariance (2d array) and modes (3d array), modeCoords (2d array), vignette flag (None or 1d array) and output (2d array), nthreads (optional)\n");
      return NULL;
  }
  if(phaseCov->nd!=2 || modes->nd!=3 || output->nd!=2 || modeCoords->nd!=2){
      printf("Dimensions for inputs should be 2,3,2\n");
      return NULL;
  }
  if(phaseCov->descr->type_num!=NPY_FLOAT ||
     modes->descr->type_num!=NPY_FLOAT ||
     output->descr->type_num!=NPY_FLOAT ||
     (modeCoords->descr->kind!='i' || modeCoords->descr->elsize!=sizeof(int))
    ){
      printf("Arrays should be 32 bit floating point (32 bit int for modecoords)\n");
      return NULL;
  }

  nmodes=modes->dimensions[0];
  ny=modes->dimensions[1];
  nx=modes->dimensions[2];
  if(phaseCov->dimensions[0]!=phaseCov->dimensions[1] || output->dimensions[0]!=nmodes || output->dimensions[1]!=nmodes || modeCoords->dimensions[0]!=nmodes || modeCoords->dimensions[1]!=2){
      printf("Shape of phaseCov or output or modeCoords is wrong\n");
      return NULL;
  }
  if(PyArray_Check(vigObj)){
    vig=(PyArrayObject*)vigObj;
    if(vig->nd!=1 || vig->descr->kind!='i' || vig->descr->elsize!=sizeof(int) || vig->dimensions[0]!=nmodes){
      printf("vig array should be 1D int32 size nmodes\n");
      return NULL;
    }
  }else
    vig=NULL;

  Py_BEGIN_ALLOW_THREADS;
  thread=(pthread_t*)malloc(sizeof(pthread_t)*nthreads);
  threadInfo=(runStruct*)malloc(sizeof(runStruct)*nthreads);
  //nleft=nmodes;
  NperThread2=(nmodes*nmodes)/nthreads;//twice the number of points evaluated per thread.
  //nPerThread=(nmodes*nmodes)/4./nthreads;//no of points evaluated per thread.
  threadInfo[0].phaseCov=phaseCov;
  threadInfo[0].modes=modes;
  threadInfo[0].vig=vig;
  threadInfo[0].modeCoords=modeCoords;
  threadInfo[0].output=output;
  threadInfo[0].start=0;
  s=threadInfo[0].start;
  //threadInfo[0].end=nleft/nthreads;
  //threadInfo[0].end=(int)(nmodes+0.5-sqrt(nmodes*nmodes-NperThread2+s*s-2*nmodes*s));//nleft/nthreads;
  //e=(int)(0.5+(nmodes-sqrt(nmodes*nmodes-4*(nPerThread+s*nmodes-s*s)))/2);
  e=(int)(nmodes+0.5-sqrt(nmodes*nmodes-NperThread2+s*s-2*nmodes*s));
  if(nthreads==1)
    //e=(nmodes+1)/2;//just in case of floating point round errors.
    e=nmodes;
  threadInfo[0].end=e;
  //printf("phasecov %d (%d elements)\n",e,(e-s)*(nmodes-e-s));
  printf("phasecov %d (%g elements)\n",e,(2*nmodes-s-e)/2.*(e-s));
  threadInfo[0].msg=e-s;
  //nleft=nmodes-threadInfo[0].end;
  for(t=1; t<nthreads; t++){
    threadInfo[t].phaseCov=phaseCov;
    threadInfo[t].modes=modes;
    threadInfo[t].vig=vig;
    threadInfo[t].modeCoords=modeCoords;
    threadInfo[t].output=output;
    threadInfo[t].start=e;//threadInfo[t-1].end;
    //threadInfo[t].end=threadInfo[t].start+nleft/(nthreads-t);
    s=e;//threadInfo[t].start;
    //e=(int)(0.5+(nmodes-sqrt(nmodes*nmodes-4*(nPerThread+s*nmodes-s*s)))/2);
    e=(int)(nmodes+0.5-sqrt(nmodes*nmodes-NperThread2+s*s-2*nmodes*s));
    if(t==nthreads-1)
      //e=(nmodes+1)/2;
      e=nmodes;
    threadInfo[t].end=e;//(int)(nmodes+0.5-sqrt(nmodes*nmodes-NperThread2+s*s-2*nmodes*s));
    //nleft=nmodes-threadInfo[t].end;
    threadInfo[t].msg=0;
    printf("phasecov %d (%g elements)\n",e,(2*nmodes-s-e)/2.*(e-s));
    
    //printf("phasecov %d (%g elements)\n",threadInfo[t].end,(2*nmodes-threadInfo[t].start-threadInfo[t].end)/2.*(threadInfo[t].end-threadInfo[t].start));
    pthread_create(&thread[t],NULL,(void*)covWorkerLocal,&threadInfo[t]);
  }

  covWorkerLocal(&threadInfo[0]);
  for(t=1; t<nthreads; t++){
    pthread_join(thread[t],NULL);
  }
  //printf("freeing\n");
  free(thread);
  free(threadInfo);
  //printf("freed\n");


  printf("Creating phase covariance: Done               \n");
  Py_END_ALLOW_THREADS;
  return Py_BuildValue("");
}

#define NRANSI
#include "nrutil.h"

void splin4(float *x1a, float *x2a, float *ya,int di, float *y2a, int m, int n, int mout, int nout,float *x1, float *x2, float *mxout,int ddi,int ddj)
{
  int i,j;
  float *ytmp,*yytmp;
  float y=0;
  void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);
  void splint(float xa[], float ya[], float y2a[], int n, float x, float *y);
  
  ytmp=vector(1,m);
  yytmp=vector(1,m);
  for (i=0;i<nout;i++){
    for (j=0;j<m;j++)
      splint(x2a-1,&ya[j*di-1],&y2a[j*n-1],n,x2[i],&yytmp[j+1]);
    spline(x1a-1,yytmp,m,1.0e30,1.0e30,ytmp);
    for (j=0;j<mout;j++){
      splint(x1a-1,yytmp,ytmp,m,x1[j],&y);
      mxout[i*ddj+ddi*j]=y; 
    }
  }
  free_vector(yytmp,1,m);
  free_vector(ytmp,1,m);
}
#undef NRANSI


int covWorkerQuick(runStruct *runInfo){
  int g,h,i,j,k,l,nact,m,n;
  float yshift,xshift;
  float *coords,*pc,*xin,*yin,*xout,*yout,*deriv,*tmparr,*mode,*mode2,*out;
  float partint;
  int di,pcdi,pcdj,odi,odj,msg;
  out=runInfo->out;
  odi=runInfo->odi;
  odj=runInfo->odj;
  mode=runInfo->mode;
  nact=runInfo->nact;
  coords=runInfo->coords;
  m=runInfo->ny;
  n=runInfo->nx;
  yin=xin=runInfo->xin;
  di=runInfo->di;
  pc=runInfo->pc;
  pcdi=runInfo->pcdi;
  pcdj=runInfo->pcdj;
  msg=runInfo->msg;
  yout=malloc(sizeof(float)*m);
  xout=malloc(sizeof(float)*n);
  tmparr=malloc(sizeof(float)*m*n);
  deriv=malloc(sizeof(float)*m*n);

  for(g=runInfo->start; g<runInfo->end; g++){
    for(h=0; h<nact; h++){
      yshift=coords[g]-(int)coords[g];
      xshift=coords[nact+h]-(int)coords[nact+h];
      if(msg){
	printf("%d, %d shift %g %g coords %g %g     \r",g,h,yshift,xshift,coords[g],coords[h]);
	fflush(NULL);
      }
      if(yshift!=0 || xshift!=0){//do an interpolation of the mode...
	//printf("\ndoing interpolation\n");
	for(i=0; i<m; i++){
	  yout[i]=yin[i]-yshift;
	  xout[i]=xin[i]-xshift;
	}
	mode2=tmparr;
	for(i=0; i<m; i++)
	  spline(&xin[-1],&mode[i*n-1],n,1e30,1e30,&deriv[i*n-1]);
	splin4(yin,xin,mode,m,deriv,m,n,m,n,yout,xout,mode2,n,1);
      }else{
	mode2=mode;
      }
      for(i=0; i<m; i++){
	for(j=0; j<n; j++){
	  partint=0.;
	  for(k=0; k<m; k++){
	    for(l=0; l<n; l++){
	      partint+=pc[abs(k-i+(int)coords[g])*pcdi+abs(l-j+(int)coords[nact+h])*pcdj]*mode2[k*n+l];
	    }
	  }
	  out[g*odi+h*odj]+=partint*mode[i*n+j];
	}
      }
    }
  }
  free(yout);
  free(xout);
  free(tmparr);
  free(deriv);
  return 0;
}

static PyObject *covarianceQuick(PyObject *self,PyObject *args)
{
  //A version of above that uses a single generic mode to do the bulk of the computation.
  PyArrayObject	*phaseCov,*pymode,*output,*modeCoords,*outPhs;
  PyObject *phaseCovObj=NULL;
  //PyObject *vigObj;
  //int nmodes,nx,ny;
  //int odi,odj,pcdi,pcdj,mdi,mdj,mdk;
  //int t;
  //float partint;
  //int nthreads=1;
  pthread_t *thread;
  runStruct *threadInfo;
  float *coords;
  //int NperThread2;
  int s,e,t,i,j,m,n;
  float *pc,*xin,*yin,*mode,*out,*outbig=NULL;
  int pcdi,pcdj,di,dj,odi,odj;
  int nthreads=1;
  int nact;
  /*
  int g,h,j,k,l;
  float yshift,xshift;
  float partint;
  float *tmparr,*mode2,*deriv,*yout,*xout;
  */
  //float nPerThread;
  //inputs needed:
  //phiphi, the phase covariance matrix, shape (y,x)
  //the mirror modes, a 3d numeric array with shape (nmodes, y,x)

  //If you change this function, then please also increment the version 
  //number in aosim/util/phaseCoveriance.py

  if (!PyArg_ParseTuple(args, "O!O!O!O!|Oi", &PyArray_Type ,&phaseCov,&PyArray_Type ,&pymode,&PyArray_Type,&modeCoords,&PyArray_Type,&output,&phaseCovObj,&nthreads)){
      printf("Inputs should be phase covariance (2d array) and mode (2d array), modeCoords (2d array), and output (2d array), phaseCov output (optional, array or None), nthreads (optional)\n");
      return NULL;
  }
  if(phaseCov->nd!=2 || pymode->nd!=2 || output->nd!=2 || modeCoords->nd!=2){
      printf("Dimensions for inputs should be 2,2,2,2\n");
      return NULL;
  }
  if(phaseCov->descr->type_num!=NPY_FLOAT ||
     pymode->descr->type_num!=NPY_FLOAT ||
     output->descr->type_num!=NPY_FLOAT ||
     modeCoords->descr->type_num!=NPY_FLOAT
     ){
      printf("Arrays should be 32 bit floating point\n");
      return NULL;
  }
  nact=output->dimensions[0];
  coords=(float*)(modeCoords->data);
  if(modeCoords->dimensions[0]!=2 || modeCoords->dimensions[1]!=nact || modeCoords->strides[1]!=sizeof(float) || modeCoords->strides[0]!=nact*sizeof(float)){
    printf("modeCoords should have shape (2,nact) and be contiguous\n");
    return NULL;
  }
  if(pymode->dimensions[0]!=pymode->dimensions[1]){
    printf("mode should be square\n");
    return NULL;
  }
  if(pymode->strides[1]!=sizeof(float) || pymode->strides[0]!=sizeof(float)*pymode->dimensions[1]){
    printf("Mode should be contig\n");
    return NULL;
  }
  if(phaseCovObj!=NULL && PyArray_Check(phaseCovObj)){
    outPhs=(PyArrayObject*)phaseCovObj;
    if(outPhs->descr->type_num!=NPY_FLOAT){
      printf("phaseCov output must be float32\n");
      return NULL;
    }
    j=1;
    for(i=0; i<outPhs->nd; i++){
      j*=outPhs->dimensions[i];
    }
    if(j!=nact*nact*nact*nact){
      printf("outPhs wrong size\n");
      return NULL;
    }
    if(!PyArray_ISCONTIGUOUS(outPhs)){
      printf("outPhs not contiguous\n");
      return NULL;
    }
    outbig=(float*)outPhs->data;
  }

  pc=(float*)phaseCov->data;
  pcdi=phaseCov->strides[0]/sizeof(float);
  pcdj=phaseCov->strides[1]/sizeof(float);
  m=pymode->dimensions[0];
  n=pymode->dimensions[1];
  di=pymode->strides[0];
  dj=pymode->strides[1];
  mode=(float*)pymode->data;
  odi=output->strides[0]/sizeof(float);
  odj=output->strides[1]/sizeof(float);
  out=(float*)output->data;
  xin=malloc(sizeof(float)*n);
  for(i=0; i<n; i++){
    xin[i]=i;
  }
  yin=xin;
  Py_BEGIN_ALLOW_THREADS;
  thread=(pthread_t*)malloc(sizeof(pthread_t)*nthreads);
  threadInfo=(runStruct*)malloc(sizeof(runStruct)*nthreads);
  threadInfo[0].pc=pc;
  threadInfo[0].pcdi=pcdi;
  threadInfo[0].pcdj=pcdj;
  threadInfo[0].out=out;
  threadInfo[0].odi=odi;
  threadInfo[0].odj=odj;
  threadInfo[0].mode=mode;
  threadInfo[0].nact=nact;
  threadInfo[0].coords=coords;
  threadInfo[0].ny=m;
  threadInfo[0].nx=n;
  threadInfo[0].xin=xin;
  threadInfo[0].start=s=0;
  e=nact/nthreads;
  if(nthreads==1)
    e=nact;
  threadInfo[0].end=e;
  threadInfo[0].msg=e-s;
  for(t=1; t<nthreads; t++){
    memcpy(&threadInfo[t],&threadInfo[0],sizeof(runStruct));
    threadInfo[t].start=e;//threadInfo[t-1].end;
    s=e;//threadInfo[t].start;
    e=(nact-s)/(nthreads-t);
    if(t==nthreads-1)
      e=nact;
    threadInfo[t].end=e;//(int)(nmodes+0.5-sqrt(nmodes*nmodes-NperThread2+s*s-2*nmodes*s));
    threadInfo[t].msg=0;
    pthread_create(&thread[t],NULL,(void*)covWorkerQuick,&threadInfo[t]);
  }
  covWorkerQuick(&threadInfo[0]);
  for(t=1; t<nthreads; t++){
    pthread_join(thread[t],NULL);
  }
  //printf("freeing\n");
  free(thread);
  free(threadInfo);
  printf("\nCreating phase covariance: Done               \n");
  Py_END_ALLOW_THREADS;




  
  /*
  deriv=malloc(sizeof(float)*m*n);

  yout=malloc(sizeof(float)*m);
  xout=malloc(sizeof(float)*n);
  tmparr=malloc(sizeof(float)*m*n);
  printf("doing covarianceQuick\n");
  for(g=0; g<nact; g++){
    for(h=0; h<nact; h++){
      yshift=coords[g]-(int)coords[g];
      xshift=coords[nact+h]-(int)coords[nact+h];
      printf("%d, %d shift %g %g coords %g %g     \n",g,h,yshift,xshift,coords[g],coords[h]);
      fflush(NULL);
      if(yshift!=0 || xshift!=0){//do an interpolation of the mode...
	printf("\ndoing interpolation\n");
	for(i=0; i<m; i++){
	  yout[i]=yin[i]-yshift;
	  xout[i]=xin[i]-xshift;
	}
	mode2=tmparr;
	for(i=0; i<m; i++)
	  spline(&xin[-1],&mode[i*n-1],n,1e30,1e30,&deriv[i*n-1]);
	splin4(yin,xin,mode,m,deriv,m,n,m,n,yout,xout,mode2,n,1);
      }else{
	mode2=mode;
      }
      for(i=0; i<m; i++){
	for(j=0; j<n; j++){
	  partint=0.;
	  for(k=0; k<m; k++){
	    for(l=0; l<n; l++){
	      partint+=pc[abs(k-i+(int)coords[g])*pcdi+abs(l-j+(int)coords[nact+h])*pcdj]*mode2[k*n+l];
	    }
	  }
	  out[g*odi+h*odj]+=partint*mode[i*n+j];
	}
      }
    }
  }
  printf("\ndone covarianceQuick\n");
  free(tmparr);
  free(deriv);
  free(yout);
  free(xout);*/
  free(xin);
  
  if(outbig!=NULL){
    int g,h,i,j,k,l;
    g=nact*nact;
    h=g*nact;
    for(i=0; i<nact; i++){
      for(j=0; j<nact; j++){
	for(k=0; k<nact; k++){
	  for(l=0; l<nact; l++){
	    outbig[i*h+j*g+k*nact+l]=out[abs(i-k)*odi+abs(j-l)*odj];
	  }
	}
      }
    }

  }

  return Py_BuildValue("");
}


/* Calculate subaperture tilt covariances for Kolmogorov power spectrum
   (npup, nsubx, d, r_0, lambda) */
static PyObject *tiltcov_Kolmogorov(PyObject *self,PyObject *args)
{
  long     npup,nsubx;                /* inputs */
  double   d,r0,lamda;

  PyArrayObject   *raw;               /* output */
  npy_intp             dims[2] = {0,0};

  int      i, j, ia, ja, ib, jb, n2, n3, n4;
  double   scaling, x, dbl_intgrl, phiphi, tiltcov;
  double   *tilt, *rxy, *ra_intgrl, *rb_intgrl, *D_phi;

  if (!PyArg_ParseTuple(args, "llddd", &npup, &nsubx, &d, &r0, &lamda)) {
    return NULL;
  }

  dims[0] = nsubx;
  dims[1] = nsubx;
  //raw = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE);    /* create output array */
  raw = (PyArrayObject *)PyArray_New(&PyArray_Type,2, dims, NPY_DOUBLE,NULL,NULL,0,0,NULL);    /* create output array */

  /* allocate arrays needed for calculation */
  n2=npup*npup;
  n3=n2*npup;
  n4=n3*npup;
  tilt=calloc(npup, sizeof(double));
  rxy=calloc(npup, sizeof(double));
  ra_intgrl=calloc(n2, sizeof(double));
  rb_intgrl=calloc(n2, sizeof(double));
  D_phi=calloc(n4, sizeof(double));

  scaling = 3.*pow(lamda/(M_PI*d),2.);

  for(i = 0; i < npup; ++i) {
    rxy[i] = ((double)(i-(npup/2)) + 0.5)/((double) npup);
    tilt[i] = 2.*sqrt(3.)*rxy[i];
  }

  for (i = 0; i < nsubx; ++i) {
    for (j = 0; j < nsubx; ++j) {

      /* printf("%d %d\n",i,j); */

      dbl_intgrl=0.;
      for (ia = 0; ia < npup; ++ia) {
	for (ja = 0; ja < npup; ++ja) {
	  for (ib = 0; ib < npup; ++ib) {
	    for (jb = 0; jb < npup; ++jb) {
	      x = sqrt(pow(((double) i)-rxy[ia]+rxy[ib],2.)+pow(((double) j)-rxy[ja]+rxy[jb],2.));
	      D_phi[n3*ia + n2*ja + npup*ib + jb] = Kolstrucfunc(x,r0/d);
	      ra_intgrl[npup*ib + jb] += D_phi[n3*ia + n2*ja + npup*ib + jb];
	      rb_intgrl[npup*ia + ja] += D_phi[n3*ia + n2*ja + npup*ib + jb];
	      dbl_intgrl += D_phi[n3*ia + n2*ja + npup*ib + jb];
	    }
	  }
	}
      }
      tiltcov=0.;
      for (ia = 0; ia < npup; ++ia) {
	for (ja = 0; ja < npup; ++ja) {
	  for (ib = 0; ib < npup; ++ib) {
	    for (jb = 0; jb < npup; ++jb) {
	      phiphi = 0.5*((ra_intgrl[npup*ib + jb] + rb_intgrl[npup*ia + ja])/((double) n2));
	      phiphi -= 0.5*D_phi[n3*ia + n2*ja + npup*ib + jb];
	      phiphi -= 0.5*dbl_intgrl/pow((double) npup,4.);
	      tiltcov += phiphi*tilt[ia]*tilt[ib];
	    }
	  }
	}
      }

      *(double *)(raw->data+i*raw->strides[0]+j*raw->strides[1]) = scaling*tiltcov/((double) n4);
    }
  }

  /* deallocate arrays */
  free(tilt);
  free(rxy);
  free(ra_intgrl);
  free(rb_intgrl);
  free(D_phi);

  return PyArray_Return(raw);
}



/* Calculate subaperture tilt covariances for Boreman-Dainty generalised power spectrum
   (npup, nsubx, d, rho_0, beta, lambda) */
static PyObject *tiltcov_BoremanDainty(PyObject *self,PyObject *args)
{
  long  npup,nsubx;    /* inputs */
  double  d,r0,beta,lamda;
 
  PyArrayObject  *raw;    /* output */
  npy_intp  dims[2] = {0,0};

  int  i, j, ia, ja, ib, jb, n2, n3, n4;
  double  scaling, gambet, x, dbl_intgrl, phiphi, tiltcov;
  double  *tilt, *rxy, *ra_intgrl, *rb_intgrl, *D_phi;

  if (!PyArg_ParseTuple(args, "lldddd", &npup, &nsubx, &d, &r0, &beta, &lamda)) {
    return NULL;
  }

  dims[0] = nsubx;
  dims[1] = nsubx;
  raw = (PyArrayObject *)PyArray_SimpleNew(2, dims, PyArray_DOUBLE);    /* create output array */

  /* allocate arrays needed for calculation */
  n2=npup*npup;
  n3=n2*npup;
  n4=n3*npup;
  tilt=calloc(npup, sizeof(double));
  rxy=calloc(npup, sizeof(double));
  ra_intgrl=calloc(n2, sizeof(double));
  rb_intgrl=calloc(n2, sizeof(double));
  D_phi=calloc(n4, sizeof(double));

  scaling = 3.*pow(lamda/(M_PI*d),2.);

  for(i = 0; i < npup; ++i) {
    rxy[i] = ((double)(i-(npup/2)) + 0.5)/((double) npup);
    tilt[i] = 2.*sqrt(3.)*rxy[i];
  }

  gambet=gammabeta(beta);
  for (i = 0; i < nsubx; ++i) {
    for (j = 0; j < nsubx; ++j) {

      /* printf("%d %d\n",i,j); */

      dbl_intgrl=0.;
      for (ia = 0; ia < npup; ++ia) {
	for (ja = 0; ja < npup; ++ja) {
	  for (ib = 0; ib < npup; ++ib) {
	    for (jb = 0; jb < npup; ++jb) {
	      x = sqrt(pow(((double) i)-rxy[ia]+rxy[ib],2.)+pow(((double) j)-rxy[ja]+rxy[jb],2.));
	      D_phi[n3*ia + n2*ja + npup*ib + jb] = gambet*BDstrucfunc(x,r0/d,beta);
	      ra_intgrl[npup*ib + jb] += D_phi[n3*ia + n2*ja + npup*ib + jb];
	      rb_intgrl[npup*ia + ja] += D_phi[n3*ia + n2*ja + npup*ib + jb];
	      dbl_intgrl += D_phi[n3*ia + n2*ja + npup*ib + jb];
	    }
	  }
	}
      }
      tiltcov=0.;
      for (ia = 0; ia < npup; ++ia) {
	for (ja = 0; ja < npup; ++ja) {
	  for (ib = 0; ib < npup; ++ib) {
	    for (jb = 0; jb < npup; ++jb) {
	      phiphi = 0.5*((ra_intgrl[npup*ib + jb] + rb_intgrl[npup*ia + ja])/((double) n2));
	      phiphi -= 0.5*D_phi[n3*ia + n2*ja + npup*ib + jb];
	      phiphi -= 0.5*dbl_intgrl/pow((double) npup,4.);
	      tiltcov += phiphi*tilt[ia]*tilt[ib];
	    }
	  }
	}
      }

      *(double *)(raw->data+i*raw->strides[0]+j*raw->strides[1]) = scaling*tiltcov/((double) n4);
    }
  }

  /* deallocate arrays */
  free(tilt);
  free(rxy);
  free(ra_intgrl);
  free(rb_intgrl);
  free(D_phi);

  return PyArray_Return(raw);
}



/* Calculate subaperture tilt covariances for Von Karman power spectrum
   (npup, nsubx, d, r_0, L_0, lambda) */
static PyObject *tiltcov_VonKarman(PyObject *self,PyObject *args)
     
{
  long     npup,nsubx;                /* inputs */
  double   d,r0,L0,lamda;

  PyArrayObject   *raw;               /* output */
  npy_intp             dims[2] = {0,0};

  int      i, j, ia, ja, ib, jb, n2, n3, n4;
  double   scaling, x, dbl_intgrl, phiphi, tiltcov;
  double   *tilt, *rxy, *ra_intgrl, *rb_intgrl, *D_phi;

  if (!PyArg_ParseTuple(args, "lldddd", &npup, &nsubx, &d, &r0, &L0, &lamda)) {
    return NULL;
  }

  dims[0] = nsubx;
  dims[1] = nsubx;
  raw = (PyArrayObject *)PyArray_SimpleNew(2, dims, PyArray_DOUBLE);    /* create output array */

  /* allocate arrays needed for calculation */
  n2=npup*npup;
  n3=n2*npup;
  n4=n3*npup;
  tilt=calloc(npup, sizeof(double));
  rxy=calloc(npup, sizeof(double));
  ra_intgrl=calloc(n2, sizeof(double));
  rb_intgrl=calloc(n2, sizeof(double));
  D_phi=calloc(n4, sizeof(double));

  scaling = 3.*pow(lamda/(M_PI*d),2.);

  for(i = 0; i < npup; ++i) {
    rxy[i] = ((double)(i-(npup/2)) + 0.5)/((double) npup);
    tilt[i] = 2.*sqrt(3.)*rxy[i];
  }

  for (i = 0; i < nsubx; ++i) {
    for (j = 0; j < nsubx; ++j) {

      /* printf("%d %d\n",i,j); */

      dbl_intgrl=0.;
      for (ia = 0; ia < npup; ++ia) {
	for (ja = 0; ja < npup; ++ja) {
	  for (ib = 0; ib < npup; ++ib) {
	    for (jb = 0; jb < npup; ++jb) {
	      x = sqrt(pow(((double) i)-rxy[ia]+rxy[ib],2.)+pow(((double) j)-rxy[ja]+rxy[jb],2.));
	      D_phi[n3*ia + n2*ja + npup*ib + jb] = VKstrucfunc(x,r0/d,L0/d);
	      ra_intgrl[npup*ib + jb] += D_phi[n3*ia + n2*ja + npup*ib + jb];
	      rb_intgrl[npup*ia + ja] += D_phi[n3*ia + n2*ja + npup*ib + jb];
	      dbl_intgrl += D_phi[n3*ia + n2*ja + npup*ib + jb];
	    }
	  }
	}
      }
      tiltcov=0.;
      for (ia = 0; ia < npup; ++ia) {
	for (ja = 0; ja < npup; ++ja) {
	  for (ib = 0; ib < npup; ++ib) {
	    for (jb = 0; jb < npup; ++jb) {
	      phiphi = 0.5*((ra_intgrl[npup*ib + jb] + rb_intgrl[npup*ia + ja])/((double) n2));
	      phiphi -= 0.5*D_phi[n3*ia + n2*ja + npup*ib + jb];
	      phiphi -= 0.5*dbl_intgrl/pow((double) npup,4.);
	      tiltcov += phiphi*tilt[ia]*tilt[ib];
	    }
	  }
	}
      }

      *(double *)(raw->data+i*raw->strides[0]+j*raw->strides[1]) = scaling*tiltcov/((double) n4);
    }
  }

  /* deallocate arrays */
  free(tilt);
  free(rxy);
  free(ra_intgrl);
  free(rb_intgrl);
  free(D_phi);

  return PyArray_Return(raw);
}



/* define a methods table for the module */
static PyMethodDef tiltcov_methods[] = {
                                        {"covariance", covariance, METH_VARARGS},
                                        {"covarianceLocal", covarianceLocal, METH_VARARGS},
                                        {"covarianceQuick", covarianceQuick, METH_VARARGS},

					{"Kolmogorov", tiltcov_Kolmogorov, METH_VARARGS},
                                        {"BoremanDainty", tiltcov_BoremanDainty, METH_VARARGS},
                                        {"VonKarman", tiltcov_VonKarman, METH_VARARGS},
                                        {NULL, NULL} };



/* initialisation - register the methods with the Python interpreter */
void initphaseCov(void)
{
        (void) Py_InitModule("phaseCov", tiltcov_methods);
        import_array();
}

