
/* Numpy extension to do matrix interpolation  using Num. Rec. Cubic spline interp */

#include <stdio.h>
#include <math.h>
#include "Python.h"

#include "numpy/arrayobject.h"//lib/python2.5/site-packages/numpy/core/include/numpy/arrayobject.h
/*
#define NUMERIC
#ifdef NUMERIC
#include "Numeric/arrayobject.h"
#else
#include "numarray/arrayobject.h"
#endif
*/
#include "nr.h"
#include "nrutil.h"
#include <time.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

/* =============================================================================*/

#define NRANSI
#include "nrutil.h"


/* Linear interpolation in 2D (requires 64-bit arrays) */
static PyObject *gslLinIntrp(PyObject *self,PyObject *args){
	PyArrayObject	*pyin,*pyout;
	npy_intp		*pyinDims,*pyoutDims;
	double			*pyinData,*pyoutData;

	double		*x1,*x2,*y1,*ytmp;
	int			i,j,n1,n2;

	gsl_interp 	*interpObj;
	gsl_interp_accel *interpAcc ;
	size_t  ninterp;


	/* handle input arrays */
	if (!PyArg_ParseTuple (args, "O!O!", &PyArray_Type, &pyin, &PyArray_Type, &pyout))
		return NULL;
	
	pyinDims=PyArray_DIMS(pyin);
	pyoutDims=PyArray_DIMS(pyout);
	n1=(int) pyinDims[0];
	n2=(int) pyoutDims[0];
	pyinData = (double *)PyArray_DATA(pyin);
	pyoutData = (double *)PyArray_DATA(pyout);


	/* allocate populate working arrays */
	x1=calloc(n1, sizeof(double));
	x2=calloc(n2, sizeof(double));
	y1=calloc(n1, sizeof(double));
	ytmp=calloc(n1*n2, sizeof(double));

	for (i=0; i<n1; ++i) {
		x1[i] = (double)i / (double)(n1-1);
	}
	for (i=0; i<n2; ++i) {
		x2[i] = (0.5 + (double)i) / (double)n2;
	}


	/* allocate interpolation stuff */
	ninterp=(size_t)n1;
	interpObj=gsl_interp_alloc(gsl_interp_linear, ninterp);
	interpAcc=gsl_interp_accel_alloc();

	
	/* do interpolation pass in first dimension */
	for (j=0; j<n1; ++j) {
		for (i=0; i<n1; ++i) {
			y1[i]=pyinData[i*n1+j];
		}
		gsl_interp_init(interpObj, x1, y1, ninterp);

		for (i=0; i<n2; ++i) {
			ytmp[j*n2+i]=gsl_interp_eval(interpObj, x1, y1, x2[i], interpAcc);
		}	
	}
	

	/* do interpolation pass in second dimension and put results in output array */
	for (j=0; j<n2; ++j) {
		for (i=0; i<n1; ++i) {
			y1[i]=ytmp[i*n2+j];
		}
		gsl_interp_init(interpObj, x1, y1, ninterp);
		for (i=0; i<n2; ++i) {
			pyoutData[j*n2+i]=gsl_interp_eval(interpObj, x1, y1, x2[i], interpAcc);
		}	
	}	
	
	
	/* tidy up */
	//printf("intrp.interp2d: Tidying up\n");
	gsl_interp_accel_free(interpAcc);
	gsl_interp_free(interpObj);
	free(x1);
	free(x2);
	free(y1);
	free(ytmp);
	
	return Py_BuildValue("");   /* return None */
}



/* Cubic spline interpolation in 2D (requires 64-bit arrays) */
static PyObject *gslCubSplineIntrp(PyObject *self,PyObject *args){
  PyArrayObject	*pyin,*pyout;
  PyObject *pyyin,*pyxin,*pyyout,*pyxout;
  npy_intp		*pyinDims,*pyoutDims;
  double			*pyinData,*pyoutData;
  float *pyinDataf,*pyoutDataf;
  
  //double *yin=NULL,*xin=NULL,*yout=NULL,*xout=NULL;

  int		zpad=10;   /* no of pixels to zero-pad input by */
  
  double		*x1,*x2,*x3,*x4,*y1,*y2,*ytmp;
  double *x1free=NULL,*x2free=NULL;
  double		dx1,dx2,dx3,dx4;
  int			i,j,n1x,n1y,n2x,n2y,n1xp,n1yp;
  double *x1usr=NULL,*x2usr=NULL,*x3usr=NULL,*x4usr=NULL;
  gsl_spline 	*interpObj;
  //gsl_interp 	*interpObj;
  gsl_interp_accel *interpAcc ;
  //size_t  ninterp;
  int s1y,s1x,s2y,s2x,insize,outsize;
  
  /* handle input arrays */
  if (!PyArg_ParseTuple (args, "O!OOOOO!|i", &PyArray_Type, &pyin, &pyyin,&pyxin,&pyyout,&pyxout,&PyArray_Type, &pyout,&zpad)){
    printf("Usage: mxin, yin,xin,yout,xout,mxout,zpad (optional, default=10)\nxin/yin should be eg x1=(numpy.arange(n)/(n-1)).astype('d') where n is mxin.shape[0]+npad*2, and xout/yout should be (numpy.arange(nn)+0.5)*(mxin.shape[0]-1)/nn/(n-1)+x1[zpad]\nAny of yin/xin/yout/xout can be None.");
    return NULL;
  }
  if(PyArray_NDIM(pyin)!=2 || PyArray_NDIM(pyout)!=2){
    printf("in and out must be 2d\n");
    return NULL;
  }
  pyinDims=PyArray_DIMS(pyin);
  pyoutDims=PyArray_DIMS(pyout);
  n1y=(int) pyinDims[0];
  n1x=(int) pyinDims[1];
  n2y=(int) pyoutDims[0];
  n2x=(int) pyoutDims[1];
  pyinData = (double *)PyArray_DATA(pyin);
  pyoutData = (double *)PyArray_DATA(pyout);
  pyinDataf = (float *)PyArray_DATA(pyin);
  pyoutDataf = (float *)PyArray_DATA(pyout);
  switch(PyArray_TYPE(pyin)){
  case NPY_FLOAT:
    insize=4;
    break;
  case NPY_DOUBLE:
    insize=8;
    break;
  default:
    insize=0;
    break;
  }
  switch(PyArray_TYPE(pyout)){
  case NPY_FLOAT:
    outsize=4;
    break;
  case NPY_DOUBLE:
    outsize=8;
    break;
  default:
    outsize=0;
    break;
  }
  if(insize==0 || outsize==0){
    printf("in and out must be float32 or float64\n");
    return NULL;
  }
  s1y=PyArray_STRIDE(pyin,0)/insize;
  s1x=PyArray_STRIDE(pyin,1)/insize;
  s2y=PyArray_STRIDE(pyout,0)/outsize;
  s2x=PyArray_STRIDE(pyout,1)/outsize;
  n1xp=n1x+2*zpad;
  n1yp=n1y+2*zpad;

  //Now check the optional input arrays...
  if(PyArray_Check(pyyin)){
    if(!PyArray_ISCONTIGUOUS((PyArrayObject*)pyyin) || PyArray_TYPE((PyArrayObject*)pyyin)!=NPY_DOUBLE){
      printf("yin must be contiguous, float64\n");
      return NULL;
    }
    if(PyArray_NDIM((PyArrayObject*)pyyin)==1 && PyArray_DIM((PyArrayObject*)pyyin,0)==n1yp){//correct shape...
      x1usr=(double*)PyArray_DATA((PyArrayObject*)pyyin);
    }else{
      printf("yin is wrong shape, should be %d\n",n1yp);
      return NULL;
    }
  }
  if(PyArray_Check(pyxin)){
    if(!PyArray_ISCONTIGUOUS((PyArrayObject*)pyxin) || PyArray_TYPE((PyArrayObject*)pyxin)!=NPY_DOUBLE){
      printf("xin must be contiguous, float64\n");
      return NULL;
    }
    if(PyArray_NDIM((PyArrayObject*)pyxin)==1 && PyArray_DIM((PyArrayObject*)pyxin,0)==n1xp){//correct shape...
      x3usr=(double*)PyArray_DATA((PyArrayObject*)pyxin);
    }else{
      printf("xin is wrong shape, should be %d\n",n1xp);
      return NULL;
    }
  }
  if(PyArray_Check(pyyout)){
    if(!PyArray_ISCONTIGUOUS((PyArrayObject*)pyyout) || PyArray_TYPE((PyArrayObject*)pyyout)!=NPY_DOUBLE){
      printf("yout must be contiguous, float64\n");
      return NULL;
    }
    if(PyArray_NDIM((PyArrayObject*)pyyout)==1 && PyArray_DIM((PyArrayObject*)pyyout,0)==n2y){//correct shape...
      x2usr=(double*)PyArray_DATA((PyArrayObject*)pyyout);
    }else{
      printf("yout is wrong shape, should be %d\n",n2y);
      return NULL;
    }
  }
  if(PyArray_Check(pyxout)){
    if(!PyArray_ISCONTIGUOUS((PyArrayObject*)pyxout) || PyArray_TYPE((PyArrayObject*)pyxout)!=NPY_DOUBLE){
      printf("yin must be contiguous, float64\n");
      return NULL;
    }
    if(PyArray_NDIM((PyArrayObject*)pyxout)==1 && PyArray_DIM((PyArrayObject*)pyxout,0)==n2x){//correct shape...
      x4usr=(double*)PyArray_DATA((PyArrayObject*)pyxout);
    }else{
      printf("yin is wrong shape, should be %d\n",n2x);
      return NULL;
    }
  }


  /* allocate and (partially) populate working arrays */
  //printf("intrp.interp2d: n1 = %d, n2 = %d, n1p = %d\n",n1,n2,n1p);
  if(x1usr==NULL || x3usr==NULL){
    x1=malloc((n1yp>n1xp?n1yp:n1xp)*sizeof(double));
    x3=x1;
    x1free=x1;
  }
  if(x2usr==NULL || x4usr==NULL){
    x2=malloc((n2y>n2x?n2y:n2x)*sizeof(double));
    x4=x2;
    x2free=x2;
  }
  y1=malloc((n1yp>n1xp?n1yp:n1xp)*sizeof(double));//size n1yp needed
  y2=y1;//size n1xp
  ytmp=malloc(n1x*n2y*sizeof(double));
  
  dx1=1./(n1yp-1.);
  dx2=dx1 * (n1y-1.) / n2y;
  dx3=1./(n1xp-1.);
  dx4=dx3*(n1x-1.)/n2x;

  if(x1usr==NULL){
    for (i=0; i<n1yp; ++i)
      x1[i] = i * dx1;
  }else
    x1=x1usr;
  if(x2usr==NULL){
    for (i=0; i<n2y; ++i)
      x2[i] = x1[zpad] + (0.5+i) * dx2;
  }else
    x2=x2usr;
  for (i=0; i<zpad; ++i) {
    y1[i] = 0.;
    y1[zpad+n1y+i]=0.;
  }
  interpObj=gsl_spline_alloc(gsl_interp_cspline_periodic, (size_t)n1yp);
  interpAcc=gsl_interp_accel_alloc();
  
  // do interpolation pass in first dimension 
  for (j=0; j<n1x; ++j) {
    if(insize==8){
      for (i=0; i<n1y; ++i) {
	y1[zpad+i]=pyinData[i*s1y+j*s1x];
	//y1[zpad+i]=pyinData[i*n1x+j];
      }
    }else{//its a float... cast to double
      for (i=0; i<n1y; ++i) 
	y1[zpad+i]=(double)pyinDataf[i*s1y+j*s1x];
    }      
    gsl_spline_init(interpObj, x1, y1,(size_t)n1yp);
    //gsl_interp_init(interpObj, x1, y1, ninterp);
    
    for (i=0; i<n2y; ++i) {
      //ytmp[j*n2+i]=gsl_interp_eval(interpObj, x1, y1, x2[i], interpAcc);
      ytmp[j+i*n1x]=gsl_spline_eval(interpObj,x2[i], interpAcc);
    }
  }
  if(x3usr==NULL){
    for(i=0; i<n1xp; i++)
      x3[i]=i*dx3;
  }else
    x3=x3usr;
  if(x4usr==NULL){
    for(i=0; i<n2x; i++)
      x4[i]=x3[zpad]+(0.5+i)*dx4;
  }else
    x4=x4usr;
  for (i=0; i<zpad; ++i) {
    y2[i] = 0.;
    y2[zpad+n1x+i]=0.;
  }
  
  // do interpolation pass in second dimension and put results in output array
  gsl_spline_free(interpObj);
  interpObj=gsl_spline_alloc(gsl_interp_cspline_periodic, (size_t)n1xp);
  for (j=0; j<n2y; ++j) {
    memcpy(&y2[zpad],&ytmp[j*n1x],n1x*sizeof(double));
    gsl_spline_init(interpObj, x3, y2,(size_t)n1xp);
    //gsl_interp_init(interpObj, x1, y1, ninterp);
    if(outsize==8){
      for (i=0; i<n2x; ++i) {
	//pyoutData[j*n2+i]=gsl_interp_eval(interpObj, x1, y1, x2[i], interpAcc);
	pyoutData[j*s2y+i*s2x]=gsl_spline_eval(interpObj,x4[i], interpAcc);
	//pyoutData[j*n2x+i]=gsl_spline_eval(interpObj,x4[i], interpAcc);
      }	
    }else{//output is float32
      for (i=0; i<n2x; ++i)
	pyoutDataf[j*s2y+i*s2x]=(float)gsl_spline_eval(interpObj,x4[i], interpAcc);
    }

  }	
  
  /* tidy up */
  gsl_interp_accel_free(interpAcc);
  gsl_spline_free(interpObj);
  //gsl_interp_free(interpObj);
  if(x1free!=NULL)
    free(x1free);
  if(x2free!=NULL)
    free(x2);
  free(y1);
  free(ytmp);
  
  return Py_BuildValue("");   /* return None */
}


static PyObject *cubIntrpOrig(self,args)
	PyObject *self, *args;
{
	PyArrayObject	*pyin,*pyout;
	npy_intp		*pyinDims,*pyoutDims;
	double			*pyinData,*pyoutData;

	int		zpad=10;   /* no of pixels to zero-pad input by */

	double		*x1,*x2,*y1,*ytmp;
	double		dx1,dx2;
	int			i,j,n1,n2,n1p;

	gsl_interp 	*interpObj;
	gsl_interp_accel *interpAcc ;
	size_t  ninterp;


	/* handle input arrays */
	if (!PyArg_ParseTuple (args, "O!O!", &PyArray_Type, &pyin, &PyArray_Type, &pyout))
		return NULL;
	
	pyinDims=PyArray_DIMS(pyin);
	pyoutDims=PyArray_DIMS(pyout);
	n1=(int) pyinDims[0];
	n2=(int) pyoutDims[0];
	pyinData = (double *)PyArray_DATA(pyin);
	pyoutData = (double *)PyArray_DATA(pyout);


	/* allocate and (partially) populate working arrays */
	n1p=n1+2*zpad;
	//printf("intrp.interp2d: n1 = %d, n2 = %d, n1p = %d\n",n1,n2,n1p);
	x1=calloc(n1p, sizeof(double));
	x2=calloc(n2, sizeof(double));
	y1=calloc(n1p, sizeof(double));
	ytmp=calloc(n1*n2, sizeof(double));

	dx1=1./(double)(n1p-1);
	dx2=dx1 * (double)(n1-1) / (double)(n2);

	for (i=0; i<n1p; ++i) {
		x1[i] = (double)i * dx1;
	}
	for (i=0; i<n2; ++i) {
		x2[i] = x1[zpad] + (0.5+(double)i) * dx2;
	}
	for (i=0; i<zpad; ++i) {
		y1[i] = 0.;
		y1[zpad+n1+i]=0.;
	}
	

	/* allocate interpolation stuff (use periodic cspline to make sure influence functions are identical) */
	ninterp=(size_t)n1p;
	interpObj=gsl_interp_alloc(gsl_interp_cspline_periodic, ninterp);
	interpAcc=gsl_interp_accel_alloc();

	
	/* do interpolation pass in first dimension */
	for (j=0; j<n1; ++j) {
		for (i=0; i<n1; ++i) {
			y1[zpad+i]=pyinData[i*n1+j];
		}
		gsl_interp_init(interpObj, x1, y1, ninterp);

		for (i=0; i<n2; ++i) {
			ytmp[j*n2+i]=gsl_interp_eval(interpObj, x1, y1, x2[i], interpAcc);
		}	
	}
	
	
	/* do interpolation pass in second dimension and put results in output array */
	for (j=0; j<n2; ++j) {
		for (i=0; i<n1; ++i) {
			y1[zpad+i]=ytmp[i*n2+j];
		}
		gsl_interp_init(interpObj, x1, y1, ninterp);
		for (i=0; i<n2; ++i) {
			pyoutData[j*n2+i]=gsl_interp_eval(interpObj, x1, y1, x2[i], interpAcc);
		}	
	}	
	
	
	/* tidy up */
	gsl_interp_accel_free(interpAcc);
	gsl_interp_free(interpObj);
	free(x1);
	free(x2);
	free(y1);
	free(ytmp);
	
	return Py_BuildValue("");   /* return None */
}
//splint - agb version to avoid exit if h==0.
int splintagb(float xa[], float ya[], float y2a[], int n, float x, float *y)
{
	void nrerror(char error_text[]);
	int klo,khi,k;
	float h,b,a;

	klo=1;
	khi=n;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if (h == 0.0){
	  printf("Bad xa input to routine splint in interpmodule\n");
	  printf("SIMULATION NEEDS TO BE HALTED\n");
	  return 1;
	}
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
	return 0;
}
/* (C) Copr. 1986-92 Numerical Recipes Software 21%. */


void splin3(float x1a[], float x2a[], float **ya, float **y2a, int m, int n, int mout, int nout,
	float x1[], float x2[], float **mxout)
{
	int i,j;
	float *ytmp,*yytmp,*y;
	void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);
	void splint(float xa[], float ya[], float y2a[], int n, float x, float *y);

	ytmp=vector(1,m);
	yytmp=vector(1,m);
	y=calloc(1,sizeof(float)); 

	for (i=1;i<=nout;i++){

	  for (j=1;j<=m;j++)
	    //if(splintagb(x2a,ya[j],y2a[j],n,x2[i],&yytmp[j]))
	    // return 1;
	    splint(x2a,ya[j],y2a[j],n,x2[i],&yytmp[j]);
	  spline(x1a,yytmp,m,1.0e30,1.0e30,ytmp);

	  for (j=1;j<=mout;j++){
	    //if(splintagb(x1a,yytmp,ytmp,m,x1[j],y))
	    //return 1;
	    splint(x1a,yytmp,ytmp,m,x1[j],y);
	    mxout[i][j]=*y; 
	  }
	}

	free(y);
	free_vector(yytmp,1,m);
	free_vector(ytmp,1,m);
	//return 0;
}
//splin4 - agb version to avoid excessive data copying
void splin4(float x1a[], float x2a[], float *ya,int di, float **y2a, int m, int n, int mout, int nout,float x1[], float x2[], float *mxout,int ddi,int ddj)
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
      //if(splintagb(x2a-1,&ya[j*di-1],y2a[j+1],n,x2[i],&yytmp[j+1]))
      //return 1;
      splint(x2a-1,&ya[j*di-1],y2a[j+1],n,x2[i],&yytmp[j+1]);

    spline(x1a-1,yytmp,m,1.0e30,1.0e30,ytmp);
    for (j=0;j<mout;j++){
      //if(splintagb(x1a-1,yytmp,ytmp,m,x1[j],&y))
      //return 1;
      splint(x1a-1,yytmp,ytmp,m,x1[j],&y);
      mxout[i*ddj+ddi*j]=y; 
    }
  }
  free_vector(yytmp,1,m);
  free_vector(ytmp,1,m);
  //return 0;
}

#undef NRANSI


//static PyObject *interp_fastmxinterp(PyObject *self,PyObject *args){
static PyObject *fastmxinterp(PyArrayObject *pymxin,PyArrayObject *pyyin,PyArrayObject *pyxin,PyArrayObject *pyyout,PyArrayObject *pyxout,PyArrayObject *pymxout){
  //does the matrix interpolation without copying.  Much more strict on input types.  Should return same result as mxinterp.
  //Actually not really much faster, though should be very slightly.  Better use of memory too.
  //PyArrayObject *pymxin,*pyyin,*pyxin,*pyyout,*pyxout,*pymxout;
  int m,n,di,ddi,nyout,nxout,j;
  float **deriv;
  //if (!PyArg_ParseTuple (args, "O!O!O!O!O!O!", &PyArray_Type ,&pymxin, &PyArray_Type ,&pyyin, &PyArray_Type ,&pyxin, &PyArray_Type ,&pyyout, &PyArray_Type ,&pyxout, &PyArray_Type ,&pymxout)) {
  //return NULL;
  //}
  if(pymxin->descr->type_num==NPY_FLOAT && 
     pyyin->descr->type_num==NPY_FLOAT && 
     pyxin->descr->type_num==NPY_FLOAT && 
     pyyout->descr->type_num==NPY_FLOAT && 
     pyxout->descr->type_num==NPY_FLOAT && 
     pymxout->descr->type_num==NPY_FLOAT){
    //all float so okay.
    //Check dimensionality
    if(pyyin->nd==1 && pyxin->nd==1 && pyyout->nd==1 && pyxout->nd==1 && pymxin->nd==2 && pymxout->nd==2){
      //okay
      //Check that all are contiguous.
      if(pyyin->strides[0]==sizeof(float) &&
	 pyxin->strides[0]==sizeof(float) &&
	 pyyout->strides[0]==sizeof(float) &&
	 pyxout->strides[0]==sizeof(float) &&
	 pymxin->strides[1]==sizeof(float)){
	//okay
	//Now do the calc.
	//First, a handwritten of splie2.
	//splie2(yin,xin,mxin,ny,nx,deriv);
	//void splie2(float x1a[], float x2a[], float **ya, int m, int n, float **y2a)
	//x1a==yin, x2a==xin, ya==mxin, m=ny, n=nx, y2a==deriv
	m=pymxin->dimensions[0];
	n=pymxin->dimensions[1];
	di=pymxin->strides[0];
	ddi=pymxout->strides[0];
	nyout=pyyout->dimensions[0];
	nxout=pyxout->dimensions[0];
	deriv=matrix(1,m,1,n);

	for (j=0;j<m;j++)
	  spline((float*)(pyxin->data-sizeof(float)),(float*)(pymxin->data+j*di-sizeof(float)),n,1.0e30,1.0e30,deriv[j+1]);
  	splin4((float*)(pyyin->data),(float*)(pyxin->data),(float*)(pymxin->data),di/sizeof(float),deriv,m,n,nyout,nxout,(float*)(pyyout->data),(float*)(pyxout->data),(float*)(pymxout->data),pymxout->strides[0]/sizeof(float),pymxout->strides[1]/sizeof(float));
	  //return NULL;
	free_matrix(deriv,1,m,1,n);

      }else{
	//printf("fastmxinterp: contiguous condition not met\n");
	return NULL;
      }

    }else{
      //printf("fastmxinterp: Array dimensionality wrong\n");
      return NULL;
    }
       

  }else{
    //printf("fastmxinterp inputs must be float32\n");
    return NULL;
  }
  
  return Py_BuildValue(""); 


}
/* ============================================================================= */



/* Matrix interpolation function. Takes Python array and vectors of x,y coords for interpolated points */
/* Returns Python array of interp vales */
//This is a bicubic spline interpolation.
static PyObject *interp_mxinterp(self,args)
	PyObject *self, *args;

{
	PyArrayObject	*pymxin,*pymxout,*pyxin,*pyyin,*pyxout,*pyyout;
	PyObject *tmp;
	int		i,j,nd,di,dj,nx,ny,nxout,nyout,dims[2];
	//clock_t         t1,t2;
	/* float		*z; */
	float		*xin,*yin,*xout,*yout;
	float		**mxin,**mxout,**deriv;


	/* z=calloc(1,sizeof(float)); */

	if (!PyArg_ParseTuple (args, "O!O!O!O!O!O!", &PyArray_Type ,&pymxin, &PyArray_Type ,&pyyin, &PyArray_Type ,&pyxin, &PyArray_Type ,&pyyout, &PyArray_Type ,&pyxout, &PyArray_Type ,&pymxout)) {
		return NULL;
	}
	//First see if we can do the non-copy version (needs float inputs).
	if((tmp=fastmxinterp(pymxin,pyyin,pyxin,pyyout,pyxout,pymxout))!=NULL){
	  return tmp;
	}


/* get input array dimensions */


	nd=pymxin->nd;
	ny=pymxin->dimensions[0];
	nx=pymxin->dimensions[1];
	di=pymxin->strides[0];
	dj=pymxin->strides[1];

	/* printf("mxin %d %d %d %d \n",	pyyin->dimensions[0],pyxin->dimensions[0],pyyout->dimensions[0],pyxout->dimensions[0]); */


/* copy to a C array */

	mxin=matrix(1,ny,1,nx);
	deriv=matrix(1,ny,1,nx);
	if(pymxin->descr->type_num==NPY_DOUBLE){
	    /* mxin=alloc2d_float(ny,nx); */
	    for(i=0;i<ny;++i){
		for(j=0;j<nx;++j){
		    mxin[i+1][j+1] = (float)(*(double *)(pymxin->data + i*di + j*dj));
		    
		}
	    }
	}else if(pymxin->descr->type_num==NPY_FLOAT){
	    for(i=0;i<ny;++i){
		for(j=0;j<nx;++j){
		    mxin[i+1][j+1] = (float)(*(float *)(pymxin->data + i*di + j*dj));
		}
	    } 
	}else{
	    printf("ERROR in interpmodule - input array not double or float\n");
	    return NULL;
	}
/* copy python x,y coord
   vectors to C vectors */
	if(pyyin->descr->type_num!=NPY_DOUBLE || pyxin->descr->type_num!=NPY_DOUBLE || pyyout->descr->type_num!=NPY_DOUBLE || pyxout->descr->type_num!=NPY_DOUBLE){
	  printf("mxinterp: pyyin/xin/etc should be float64 unless trying to do the non-copy method (when they should be float32) - if you are trying to do the non-copy method, then some other of your inputs are currently wrong\n");
	  return NULL;
	}
	yin=vector(1,ny);
	di=pyyin->strides[0];
	for(i=0;i<ny;++i) yin[i+1] = (float)(*(double *)(pyyin->data + i*di));

	xin=vector(1,nx);
	dj=pyxin->strides[0];
	for(j=0;j<nx;++j) xin[j+1] = (float)(*(double *)(pyxin->data + j*dj));



	nyout=pyyout->dimensions[0];
	di=pyyout->strides[0];
	yout=vector(1,nyout);
	for(i=0;i<nyout;++i) yout[i+1] = (float)(*(double *)(pyyout->data + i*di));

	nxout=pyxout->dimensions[0];
	dj=pyxout->strides[0];
	xout=vector(1,nxout);
	for(j=0;j<nxout;++j) xout[j+1] = (float)(*(double *)(pyxout->data + j*dj));



/* do the interpolation */


/* 	t1=clock(); */
/* 	printf("time#1 = %ld\n",t1); */

	splie2(yin,xin,mxin,ny,nx,deriv);

/* 	t1=clock(); */
/* 	printf("time#2 = %ld\n",t1); */

	mxout=matrix(1,nyout,1,nxout);

 /* quick version of splin2 for interpd points regular grid */
  	splin3(yin,xin,mxin,deriv,ny,nx,nyout,nxout,yout,xout,mxout);
	  //return NULL;

/*  	for(i=1;i<=nyout;++i){ */
/*  	 for(j=1;j<=nxout;++j){ */
/*  		splin2(yin,xin,mxin,deriv,ny,nx,yout[i],xout[j],z); */
/*  		mxout[i][j]=*z; */
/*  	 } */
/*  	} */

/* 	t1=clock(); */
/* 	printf("time#3 = %ld\n",t1); */

/* create new python array for output */

	dims[0]=nxout;
	dims[1]=nyout;
/* 	pymxout=(PyArrayObject *)PyArray_FromDims(2,dims,PyArray_DOUBLE); */


/* populate output Numpy array with the interpolated vals */

/* 	t1=clock(); */



	dj=pymxout->strides[0];
	di=pymxout->strides[1];

	if(pymxout->descr->type_num==NPY_DOUBLE){
	    for(i=0;i<nyout;++i){
		for(j=0;j<nxout;++j){
		    *(double *)(pymxout->data+i*di+j*dj) = (double)mxout[i+1][j+1]; 	/* !!!! this fn transposes the input mx - why ? !!! */
		}
	    }
	}else if(pymxout->descr->type_num==NPY_FLOAT){
	    for(i=0;i<nyout;++i){
		for(j=0;j<nxout;++j){
		    *(float *)(pymxout->data+i*di+j*dj) = (float)mxout[i+1][j+1]; 	/* !!!! this fn transposes the input mx - why ? !!! */
		}
	    }
	}else{
	    printf("ERROR in interpmodule - output array should be double or float\n");
	    return NULL;
	}



	free_vector(yin,1,ny);
	free_vector(xin,1,nx);
	free_vector(yout,1,nyout);
	free_vector(xout,1,nxout);

	free_matrix(mxin,1,ny,1,nx);
	free_matrix(mxout,1,nyout,1,nxout);
	free_matrix(deriv,1,ny,1,nx);



/* Return the Numpy array */

	return Py_BuildValue(""); 
}


static PyObject *interp_bicubicinterp(self,args)
	PyObject *self, *args;

{
    void bcucof(float y[], float y1[], float y2[], float y12[], 
		float d1, float d2,float **c);
	PyArrayObject	*pymxin,*pymxout,*pyxin,*pyyin,*pyxyin;
	int		i,j,k,l,m,nd,di,dj,nx,ny,nxout,nyout,dims[2];
	int diy,djy,dix,djx,dixy,djxy,diout,djout;
	//clock_t         t1,t2;
	/* float		*z; */
	float		*zin,*dy,*dx,*dxy;
	float		**c;
	int nyinterp,nxinterp;
	int iilast=-1,jjlast=-1,ii=-1,jj=-1,iin,jjn;
	float iif,jjf,xstep,ystep;
	float yinterppos,xinterppos,ansy;
	void bcucof(float y[], float y1[], float y2[], float y12[], 
		    float d1, float d2,float **c);
	if (!PyArg_ParseTuple (args, "O!O!O!O!O!", &PyArray_Type ,&pymxin, 
			       &PyArray_Type ,&pyyin, //y gradient of function
			       &PyArray_Type ,&pyxin, //x gradient of function
			       &PyArray_Type ,&pyxyin, //xy gradient of func.
			       &PyArray_Type ,&pymxout)) {
		return NULL;
	}

        // get input array dimensions 
	nd=pymxin->nd;
	ny=pymxin->dimensions[0];
	nx=pymxin->dimensions[1];
	di=pymxin->strides[0];
	dj=pymxin->strides[1];

	diy=pyyin->strides[0];
	djy=pyyin->strides[1];
	dix=pyxin->strides[0];
	djx=pyxin->strides[1];
	dixy=pyxyin->strides[0];
	djxy=pyxyin->strides[1];

	nyout=pymxout->dimensions[0];
	nxout=pymxout->dimensions[1];
	diout=pymxout->strides[0];
	djout=pymxout->strides[1];

	/* printf("mxin %d %d %d %d \n",	pyyin->dimensions[0],pyxin->dimensions[0],pyyout->dimensions[0],pyxout->dimensions[0]); */
        /* copy to a C array */
	nyinterp=nyout/(ny-1);
	nxinterp=nxout/(nx-1);
	
	if(pymxin->descr->type_num!=NPY_FLOAT ||
	   pyyin->descr->type_num!=NPY_FLOAT ||
	   pyxin->descr->type_num!=NPY_FLOAT ||
	   pyxyin->descr->type_num!=NPY_FLOAT ||
	   pymxout->descr->type_num!=NPY_FLOAT 
	    ){
	    printf("matrixes for bicubic interpolation must be float32.\n");
	    return NULL;
	}
	c=matrix(1,4,1,4);
	zin=vector(1,4);
	dy=vector(1,4);
	dx=vector(1,4);
	dxy=vector(1,4);

/*	for(i=0; i<ny-1; i++){
	    for(j=0; j<nx-1; j++){
		//iterate over grid points...
		zin[1]=(*(float*)(pymxin->data+i*di+j*dj));
		zin[2]=(*(float*)(pymxin->data+(i+1)*di+j*dj));
		zin[3]=(*(float*)(pymxin->data+(i+1)*di+(j+1)*dj));
		zin[4]=(*(float*)(pymxin->data+(i)*di+(j+1)*dj));
		dy[1]=(*(float*)(pyyin->data+i*diy+j*djy));
		dy[2]=(*(float*)(pyyin->data+(i+1)*diy+j*djy));
		dy[3]=(*(float*)(pyyin->data+(i+1)*diy+(j+1)*djy));
		dy[4]=(*(float*)(pyyin->data+i*diy+(j+1)*djy));
		dx[1]=(*(float*)(pyxin->data+i*dix+j*djx));
		dx[2]=(*(float*)(pyxin->data+(i+1)*dix+j*djx));
		dx[3]=(*(float*)(pyxin->data+(i+1)*dix+(j+1)*djx));
		dx[4]=(*(float*)(pyxin->data+i*dix+(j+1)*djx));
		dxy[1]=(*(float*)(pyxyin->data+i*dixy+j*djxy));
		dxy[2]=(*(float*)(pyxyin->data+(i+1)*dixy+j*djxy));
		dxy[3]=(*(float*)(pyxyin->data+(i+1)*dixy+(j+1)*djxy));
		dxy[4]=(*(float*)(pyxyin->data+i*dixy+(j+1)*djxy));
		//first get the coefficients...
		bcucof(zin,dy,dx,dxy,1.,1.,c);
		for(k=0; k<nyinterp; k++){
		    for(l=0; l<nxinterp; l++){
			yinterppos=(float)k/(float)(nyinterp-1);//t
			xinterppos=(float)l/(float)(nxinterp-1);//u
			ansy=0.;
			for(m=4; m>=1; m--){
			    ansy=yinterppos*ansy + ((c[m][4]*xinterppos+c[m][3])*xinterppos+c[m][2])*xinterppos+c[m][1];
			}
			//bcuint(zin,dy,dx,dxy,0,1,0,1,xinterppos,yinterppos,&zout,&dyout,&dxout);
			(*(float*)(pymxout->data+(i*nyinterp+k)*diout+(j*nxinterp+l)*djout))=ansy;
		    }
		}
	    }
	}
*/
	xstep=(float)(nx-1)/(float)(nxout-1);
	ystep=(float)(ny-1)/(float)(nyout-1);

	for(i=0; i<nyout; i++){
	    for(j=0; j<nxout; j++){
		//iterate over grid points...
		iilast=ii;
		jjlast=jj;
		iif=(float)i*ystep;
		jjf=(float)j*xstep;
		ii=(int)iif;
		jj=(int)jjf;
		iin=(int)ceil(iif);
		jjn=(int)ceil(jjf);
		if(ii==iin && ii!=ny-1)
		    iin+=1;
		if(jj==jjn && jj!=nx-1)
		    jjn+=1;
		if(iilast!=ii || jjlast!=jj){//recompute c...
		    zin[1]=(*(float*)(pymxin->data+ii*di+jj*dj));
		    zin[2]=(*(float*)(pymxin->data+(iin)*di+jj*dj));
		    zin[3]=(*(float*)(pymxin->data+(iin)*di+(jjn)*dj));
		    zin[4]=(*(float*)(pymxin->data+(ii)*di+(jjn)*dj));
		    dy[1]=(*(float*)(pyyin->data+ii*diy+jj*djy));
		    dy[2]=(*(float*)(pyyin->data+(iin)*diy+jj*djy));
		    dy[3]=(*(float*)(pyyin->data+(iin)*diy+(jjn)*djy));
		    dy[4]=(*(float*)(pyyin->data+ii*diy+(jjn)*djy));
		    dx[1]=(*(float*)(pyxin->data+ii*dix+jj*djx));
		    dx[2]=(*(float*)(pyxin->data+(iin)*dix+jj*djx));
		    dx[3]=(*(float*)(pyxin->data+(iin)*dix+(jjn)*djx));
		    dx[4]=(*(float*)(pyxin->data+ii*dix+(jjn)*djx));
		    dxy[1]=(*(float*)(pyxyin->data+ii*dixy+jj*djxy));
		    dxy[2]=(*(float*)(pyxyin->data+(iin)*dixy+jj*djxy));
		    dxy[3]=(*(float*)(pyxyin->data+(iin)*dixy+(jjn)*djxy));
		    dxy[4]=(*(float*)(pyxyin->data+ii*dixy+(jjn)*djxy));
		    //first get the coefficients...
		    bcucof(zin,dy,dx,dxy,1.,1.,c);
		    if(isnan(c[1][1])){
			printf("Got c as nan for ii iilast jj jjlast %d %d %d %d\nzin dy dx dxy\n",ii,iilast,jj,jjlast);
			for(k=1; k<=4; k++){
			    printf("%g %g %g %g\n",zin[k],dy[k],dx[k],dxy[k]);
			}

		    }
	        }
		yinterppos=iif-(float)ii;//(float)k/(float)(nyinterp-1);//t
		xinterppos=jjf-(float)jj;//(float)l/(float)(nxinterp-1);//u
		ansy=0.;
		for(m=4; m>=1; m--){
		    ansy=yinterppos*ansy + ((c[m][4]*xinterppos+c[m][3])*xinterppos+c[m][2])*xinterppos+c[m][1];
		}
		//bcuint(zin,dy,dx,dxy,0,1,0,1,xinterppos,yinterppos,&zout,&dyout,&dxout);
		if(isnan(ansy)){
		    printf("interpmodule: Got nan for %d %d\n",i,j);
		    printf("xinterppos,yinterppos %g %g\n",xinterppos,yinterppos);
		    printf("iif, jjf, ii, jj, iin, jjn %g %g %d %d %d %d\n",iif,jjf, ii,jj,iin,jjn);
		    printf("C:\n");
		    for(k=1; k<=4; k++){
			for(l=1; l<=4; l++){
			    printf("%g ",c[k][l]);
			}
			printf("\n");
		    }
		}
		(*(float*)(pymxout->data+i*diout+j*djout))=ansy;
	    }
	}


	free_matrix(c,1,4,1,4);
	free_vector(zin,1,4);
	free_vector(dy,1,4);
	free_vector(dx,1,4);
	free_vector(dxy,1,4);
	return Py_BuildValue(""); 
};




static PyObject *interp_linearshift(PyObject *self,PyObject *args){
  //This function takes an input array, a fractional x,y value (<1)
  //and shifts the array by this much into the output array (1 smaller
  //in each dimension).
  //approx 10x faster than python on 700MHz processor (AMD Athlon)
  PyArrayObject	*inarr,*outarr;
  float xshift,yshift,s1,s2,s3,s4;
  int i,j,ny,nx,di,dj,dii,djj;
  float fans;
  double dans;
  if (!PyArg_ParseTuple (args, "O!ffO!", &PyArray_Type, &inarr, &xshift,
			 &yshift, &PyArray_Type, &outarr)){
    printf("linearshift: input 2D array, x shift, y shift, output array\n");
    return NULL;
  }
  
  // get input array dimensions 
  if(inarr->nd!=2 || outarr->nd!=2){
    printf("input and output must be 2D arrays\n");
    return NULL;
  }
  if((inarr->descr->type_num!=NPY_FLOAT && inarr->descr->type_num!=NPY_DOUBLE) || inarr->descr->type_num!=outarr->descr->type_num){
    printf("input and output array types must be same and one of float or double\n");
    return NULL;
  }
  
  ny=outarr->dimensions[0];
  nx=outarr->dimensions[1];
  di=inarr->strides[0];
  dj=inarr->strides[1];
  dii=outarr->strides[0];
  djj=outarr->strides[1];
  if(inarr->dimensions[0]<ny+1 || inarr->dimensions[1]<nx+1){
    printf("WARNING (intermodule.c): input array must be at least 1 bigger than output array in each dimension\nContinuing, but output array won't be filled completely. (%d<%d+1, %d<%d+1)\n",inarr->dimensions[0],ny,inarr->dimensions[1],nx);
    ny=inarr->dimensions[0]-1;
    nx=inarr->dimensions[1]-1;
  }
  s1=(1.-xshift)*(1.-yshift);
  s2=(1.-xshift)*yshift;
  s3=xshift*(1.-yshift);
  s4=xshift*yshift;
  if(inarr->descr->type_num==NPY_FLOAT){
    for(i=0; i<ny; i++){
      for(j=0; j<nx; j++){
	fans=(*(float*)(inarr->data+i*di+j*dj))*s1;
	fans+=(*(float*)(inarr->data+(i+1)*di+j*dj))*s2;
	fans+=(*(float*)(inarr->data+i*di+(j+1)*dj))*s3;
	fans+=(*(float*)(inarr->data+(i+1)*di+(j+1)*dj))*s4;
	(*(float*)(outarr->data+i*dii+j*djj))=fans;
      }
    }
  }else{//PyArray_DOUBLE...
    for(i=0; i<ny; i++){
      for(j=0; j<nx; j++){
	dans=(*(double*)(inarr->data+i*di+j*dj))*(double)s1;
	dans+=(*(double*)(inarr->data+(i+1)*di+j*dj))*(double)s2;
	dans+=(*(double*)(inarr->data+i*di+(j+1)*dj))*(double)s3;
	dans+=(*(double*)(inarr->data+(i+1)*di+(j+1)*dj))*(double)s4;
	(*(double*)(outarr->data+i*dii+j*djj))=dans;
      }
    }
  }
  return Py_BuildValue("");
}

/* Matrix interpolation function. Takes Python array and vectors of x,y coords for interpolated points */
/* Returns Python array of interp vales */
//This is a linear interpolation.
static PyObject *interp_linearinterp(self,args)
	PyObject *self, *args;

{
	PyArrayObject	*pymxin,*pymxout,*pyxin,*pyyin,*pyxout,*pyyout;
	int		i,j;
	//clock_t         t1,t2;
	/* float		*z; */
	int xpos,xposa,xposb,ypos,yposa,yposb;
	float xa,xb,x,ya,yb,y,m,vala,valb;
	float *tmp;
	/* z=calloc(1,sizeof(float)); */

	if (!PyArg_ParseTuple (args, "O!O!O!O!O!O!", &PyArray_Type ,&pymxin, &PyArray_Type ,&pyyin, &PyArray_Type ,&pyxin, &PyArray_Type ,&pyyout, &PyArray_Type ,&pyxout, &PyArray_Type ,&pymxout)) {
		return NULL;
	}



	/* copy to a C array */
	if(pymxin->descr->type_num!=NPY_FLOAT || pymxout->descr->type_num!=NPY_FLOAT || pyxin->descr->type_num!=NPY_FLOAT || pyxout->descr->type_num!=NPY_FLOAT || pyyin->descr->type_num!=NPY_FLOAT || pyyout->descr->type_num!=NPY_FLOAT){
	    /* mxin=alloc2d_float(ny,nx); */
	  printf("ERROR in interpmodule - input arrays not float32\n");
	    return NULL;
	}
	/* copy python x,y coord vectors to C vectors */
	if(pymxin->nd!=2 || pymxout->nd!=2 || pyxin->nd!=1 || pyxout->nd!=1 || pyyin->nd!=1 || pyyout->nd!=1){
	  printf("linearinterp error - dimensions wrong\n");
	  return NULL;
	}
	if(pymxin->dimensions[0]!=pyyin->dimensions[0] || pymxin->dimensions[1]!=pyxin->dimensions[0] || pymxout->dimensions[0]!=pyyout->dimensions[0] || pymxout->dimensions[1]!=pyxout->dimensions[0]){
	  printf("linearinterp error - wrong shapes\n");
	  return NULL;
	}


	/* do the interpolation.
	   We first interpolate along one dimension, and then interpolate these
	   results along the other dimension.
	*/
	if((tmp=malloc(sizeof(float)*pyxout->dimensions[0]*pyyin->dimensions[0]))==NULL){
	  printf("Unable to allocate temporary memory\n");
	  return NULL;
	}
	xpos=0;//first position of points we have
	for(i=0; i<pyxout->dimensions[0]; i++){
	  x=*(float*)(pyxout->data+i*pyxout->strides[0]);//the requested x coord.
	  while(xpos<pyxin->dimensions[0] && *(float*)(pyxin->data+xpos*pyxin->strides[0])<x){
	    xpos++;
	  }
	  if(xpos==0){
	    xb=*(float*)(pyxin->data+0*pyxin->strides[0]);
	    xa=xb;
	    xposa=xposb=0;
	  }else if(xpos==pyxin->dimensions[0]){
	    xa=*(float*)(pyxin->data+(xpos-1)*pyxin->strides[0]);
	    xb=xa;
	    xposa=xposb=xpos-1;
	  }else{
	    xb=*(float*)(pyxin->data+xpos*pyxin->strides[0]);
	    xa=*(float*)(pyxin->data+(xpos-1)*pyxin->strides[0]);
	    xposa=xpos-1;
	    xposb=xpos;
	  }
	  //now do the interpolation along the first direction.
	  for(j=0;j<pyyin->dimensions[0]; j++){
	    vala=*(float*)(pymxin->data+xposa*pymxin->strides[1]+j*pymxin->strides[0]);
	    valb=*(float*)(pymxin->data+xposb*pymxin->strides[1]+j*pymxin->strides[0]);
	    if(xb==xa)
	      m=0;
	    else
	      m=(valb-vala)/(xb-xa);
	    tmp[j*pyxout->dimensions[0]+i]=m*(x-xa)+vala;
	  }
	}
	//now do the 2nd part of the interpolation...
	ypos=0;
	for(i=0; i<pyyout->dimensions[0]; i++){
	  //printf("i %d\n",i);
	  y=*(float*)(pyyout->data+i*pyyout->strides[0]);
	  while(ypos<pyyin->dimensions[0] && *(float*)(pyyin->data+ypos*pyyin->strides[0])<y){
	    ypos++;
	  }
	  if(ypos==0){
	    yb=((float*)pyyin->data)[0];
	    ya=yb;
	    yposa=yposb=0;
	  }else if(ypos==pyyin->dimensions[0]){
	    ya=*(float*)(pyyin->data+(ypos-1)*pyyin->strides[0]);
	    yb=ya;
	    yposb=yposa=ypos-1;
	  }else{
	    yb=*(float*)(pyyin->data+ypos*pyyin->strides[0]);
	    ya=*(float*)(pyyin->data+(ypos-1)*pyyin->strides[0]);
	    yposa=ypos-1;
	    yposb=ypos;
	  }
	  //now do interpolation along 2nd direction.
	  for(j=0; j<pyxout->dimensions[0]; j++){
	    vala=tmp[j+yposa*pyxout->dimensions[0]];
	    valb=tmp[j+yposb*pyxout->dimensions[0]];
	    if(yb==ya)
	      m=0;
	    else
	      m=(valb-vala)/(yb-ya);
	    *(float*)(pymxout->data+i*pymxout->strides[0]+j*pymxout->strides[1])=m*(y-ya)+vala;

	  }
	}
	free(tmp);
	return Py_BuildValue(""); 
}

/*
import util.dm
ms=util.dm.MirrorSurface("linear",64,9)
acts=numpy.random.random((9,9)).astype("f")
phs=ms.fit(acts)
plot(acts,0)
plot(phs,1)
*/


/* =============================================================================*/



/* define a methods table for this module */

static PyMethodDef interp_methods[] = 	{
					{"mxinterp", interp_mxinterp, METH_VARARGS}, 
					{"bicubicinterp", interp_bicubicinterp, METH_VARARGS}, 
					{"linearshift",interp_linearshift,METH_VARARGS},
					{"linearinterp", interp_linearinterp, METH_VARARGS}, 
					{"gslCubSplineInterp",gslCubSplineIntrp,METH_VARARGS},
					{"gslLinInterp",gslLinIntrp,METH_VARARGS},
					{"gslCubSplineInterpOrig",cubIntrpOrig,METH_VARARGS},

					{NULL, NULL} };


/* initialisation - register the methods with the Python interpreter */

void initinterp()
{
	(void) Py_InitModule("interp", interp_methods);
	import_array();
}

