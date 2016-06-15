
/* 
Numpy extension to do matrix interpolation using gsl Cubic spline interpolation.
UB, 2012 Aug 08: In bicubic_spline the Numerical Recipies function is still used
(could be cleaned up, but needs some effort).
*/

#include "Python.h"
#include <stdio.h>
#include <math.h>

#include "numpy/arrayobject.h"

//#include <time.h> // commented out UB, 30th Jul 2012
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include <pthread.h>     // for multi-threaded interpolation
#include "interpolate.h" // functions for multi-threaded interpolation

/* =============================================================================*/

/* Cubic spline interpolation in 2D (requires 64-bit arrays) */
// 
// Uses 
//     gsl_interp_cspline_periodic,
//                                  which requires that the y value 
// of the last point of the array is equal to the y value of the first point.
//
static PyObject *gslPeriodicCubSplineIntrp(PyObject *self,PyObject *args){
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
// END of gslPeriodicCubSplineIntrp


// 
// Urban Bitenc, 25 July 2012:
// This function was derived from gslPeriodicCubSplineIntrp. It replaced mxinterp,
// therefore it uses 
//                   gsl_interp_cspline 
//
// and does not have the zpad argument.
//
// It is MULTI-THREADED. For the default value of threads it always uses
// sysconf(_SC_NPROCESSORS_ONLN)/2 threads, regardless of the size of the input/output
// arrays.
//
// Upgraded: 1st Feb 2013 to handle cases where input and output arrays are the same.
//
static PyObject *gslCubSplineIntrp(PyObject *self,PyObject *args){

  PyArrayObject	*pyin,    *pyout;
  void          *pyinData,*pyoutData;  // void, so it can handle both float and double
  npy_intp      *pyinDims,*pyoutDims;
  PyObject *pyyin,*pyxin, *pyyout,*pyxout;

  double *x1,*x2,*x3,*x4,*ytmp;
  double *x1free=NULL,*x2free=NULL;
  double dx1,dx2,dx3,dx4;

  int i,j;
  int nInX, nInY, nOutX, nOutY;
  int sInY, sInX, sOutY, sOutX, insize, outsize;
  double *x1usr=NULL, *x2usr=NULL, *x3usr=NULL, *x4usr=NULL;

  int nThreads = 0;      // the number of threads (provided on input)
  int addToOutput=0;
  pthread_t*     thread; // vector of threads
  interp_data_t* params; // vector of parameters passed to the funct. executed in a thread
                         //  defined in interpolate.h

  int forEach;           // m/nThreads = number of lines to be processed by each thread
  int oneMore;           // m%nThreads = number of threads that have to process one line more
 
  // (1) PARSE THE INPUT arguments from python:
  if (!PyArg_ParseTuple (args, "O!OOOOO!i|i", 
			 &PyArray_Type, &pyin, 
			                &pyyin,
			                &pyxin,
			                &pyyout,
			                &pyxout,
			 &PyArray_Type, &pyout,
                 	                &addToOutput,
			                &nThreads)){
    printf("%s, line %d: parsing the arguments failed.\n", __FILE__, __LINE__);
    printf("A working example:\nxin,  yin  = numpy.arange(5).astype(TYPE),\n");
    printf("xout, yout = numpy.arange(9).astype(TYPE)/2.,\n");
    printf("mxin = numpy.zeros([5, 5], TYPE), mxout = numpy.zeros([9, 9], TYPE)\n");
    printf("TYPE must be either 'd' or 'f'; any of xin, yin, xout, yout can be 'None'.\n");
    return NULL;
  }

  // (2) EXTRACT SOME INFO AND CHECK THE INPUT

  // First check that input and output data are 2-D arrays:
  if(PyArray_NDIM(pyin)!=2 || PyArray_NDIM(pyout)!=2)
    {
      PyErr_SetString(PyExc_TypeError, "1st and 6th argument must be 2-D arrays.\n");
      return NULL;
    }

  // Extract data:
  pyinData  = (void *)PyArray_DATA(pyin);
  pyoutData = (void *)PyArray_DATA(pyout);
  pyinDims  = PyArray_DIMS(pyin);
  pyoutDims = PyArray_DIMS(pyout);
  nInY  = (int) pyinDims[0];
  nInX  = (int) pyinDims[1];
  nOutY = (int) pyoutDims[0];
  nOutX = (int) pyoutDims[1];

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

  // Check that input and output arrays are float or double:
  if(insize==0 || outsize==0){
    PyErr_SetString(PyExc_TypeError, "in and out arrays must be float32 or float64\n");
    return NULL;
  }

  // Get the strides
  sInY  = PyArray_STRIDE(pyin,0)/insize;
  sInX  = PyArray_STRIDE(pyin,1)/insize;
  sOutY = PyArray_STRIDE(pyout,0)/outsize;
  sOutX = PyArray_STRIDE(pyout,1)/outsize;

  // Further checks of the input arrays:
  //   (pyin and pyout have already been checked before up.)

  if(PyArray_Check(pyyin)){
    if(!PyArray_ISCONTIGUOUS((PyArrayObject*)pyyin) || PyArray_TYPE((PyArrayObject*)pyyin)!=NPY_DOUBLE){
      PyErr_SetString(PyExc_TypeError, "The second input parameter must be contiguous, float64\n");
      return NULL;
    }
    if(PyArray_NDIM((PyArrayObject*)pyyin)==1 && PyArray_DIM((PyArrayObject*)pyyin,0)>=nInY){//correct shape...
      x1usr=(double*)PyArray_DATA((PyArrayObject*)pyyin);
    }else{
      PyErr_SetString(PyExc_TypeError, "The second input parameter is wrong shape\n");
      return NULL;
    }
  }
  if(PyArray_Check(pyxin)){
    if(!PyArray_ISCONTIGUOUS((PyArrayObject*)pyxin) || PyArray_TYPE((PyArrayObject*)pyxin)!=NPY_DOUBLE){
      PyErr_SetString(PyExc_TypeError, "The third input parameter must be contiguous, float64\n");
      return NULL;
    }
    if(PyArray_NDIM((PyArrayObject*)pyxin)==1 && PyArray_DIM((PyArrayObject*)pyxin,0)>=nInX){//correct shape...
      x3usr=(double*)PyArray_DATA((PyArrayObject*)pyxin);
    }else{
      PyErr_SetString(PyExc_TypeError, "The third input parameter is wrong shape\n");
      return NULL;
    }
  }
  if(PyArray_Check(pyyout)){
    if(!PyArray_ISCONTIGUOUS((PyArrayObject*)pyyout) || PyArray_TYPE((PyArrayObject*)pyyout)!=NPY_DOUBLE){
      PyErr_SetString(PyExc_TypeError, "The 4th input parameter must be contiguous, float64\n");
      return NULL;
    }
    if(PyArray_NDIM((PyArrayObject*)pyyout)==1 && PyArray_DIM((PyArrayObject*)pyyout,0)>=nOutY){//correct shape...
      x2usr=(double*)PyArray_DATA((PyArrayObject*)pyyout);
    }else{
      PyErr_SetString(PyExc_TypeError, "The 4th input parameter is wrong shape\n");
      return NULL;
    }
  }
  if(PyArray_Check(pyxout)){
    if(!PyArray_ISCONTIGUOUS((PyArrayObject*)pyxout) || PyArray_TYPE((PyArrayObject*)pyxout)!=NPY_DOUBLE){
      PyErr_SetString(PyExc_TypeError, "The 5th input parameter must be contiguous, float64\n");
      return NULL;
    }
    if(PyArray_NDIM((PyArrayObject*)pyxout)==1 && PyArray_DIM((PyArrayObject*)pyxout,0)>=nOutX){//correct shape...
      x4usr=(double*)PyArray_DATA((PyArrayObject*)pyxout);
    }else{
      PyErr_SetString(PyExc_TypeError, "The 5th input parameter is wrong shape\n");
      return NULL;
    }
  }

  // (3) ALLOCATE AND (PARTIALLY) POPULATE WORKING ARRAYS
  // (This is Alastair's code, I don't touch.)
  //printf("intrp.interp2d: n1 = %d, n2 = %d, n1p = %d\n",n1,n2,n1p);
  if(x1usr==NULL || x3usr==NULL){
    x1=malloc((nInY > nInX ? nInY : nInX)*sizeof(double));
    x3=x1;
    x1free=x1;
  }
  if(x2usr==NULL || x4usr==NULL){
    x2=malloc((nOutY > nOutX ? nOutY : nOutX)*sizeof(double));
    x4=x2;
    x2free=x2;
  }
  ytmp  = malloc(nInX*nOutY*sizeof(double));
  
  dx1=1./(nInY-1.);
  dx2=dx1 * (nInY-1.) / nOutY;
  dx3=1./(nInX-1.);
  dx4=dx3*(nInX-1.)/nOutX;

  if(x1usr==NULL){
    for (i=0; i<nInY; ++i)
      x1[i] = i * dx1;
  }else
    x1=x1usr;
  if(x2usr==NULL){
    for (i=0; i<nOutY; ++i)
      x2[i] = x1[0] + (0.5+i) * dx2;
  }else
    x2=x2usr;

  if(x3usr==NULL){
    for(i=0; i<nInX; i++)
      x3[i]=i*dx3;
  }else
    x3=x3usr;
  if(x4usr==NULL){
    for(i=0; i<nOutX; i++)
      x4[i]=x3[0]+(0.5+i)*dx4;
  }else
    x4=x4usr;

  // (4) DETERMINE THE NUMBER OF THREADS:
  // If the input nThreads is < 1, use the default multi-threading.
  // Find out the number of cores (hyperthreads) and use half that number for the number of threads.
  // (On average this seems to make most sense -  see U.B.`s log book 1, p. 169)
  //
  if(nThreads < 1)
    {
      // Check the input matrix size. This was NOT investigated in details and optimised.
      // On MacBook the simulation for VLT (nInX = 16) runs 10% slower if multi-threaded.
      // On other computers (gpuserver, cpuserver, gpu2) it is only very little slower.
      // For ELT (nInX = 84) multi-threaded is faster.
      // The value of 40 used below lies just somewhere between 16 and 84.
      if( nInX < 40)
	{
	  nThreads = 1;
	}
      else
	{
	  // Find the number of cores and divide it by two:
	  nThreads = sysconf(_SC_NPROCESSORS_ONLN)/2;

	  // For the case that _SC_NPROCESSORS_ONLN returns 1, nThreads will at this point be 0.
	  // Or, if anything else funny happens and the value is even negative, correct it to 1
	  // and the process will run in a single thread. Print a warning so the user is aware of that.
	  if(nThreads < 1)
	    {
	      printf("WARNING: Multiple threading was requested for .gslCubSplineIntrp, but it\n");
	      printf("         will run in a single thread, because \"sysconf(_SC_NPROCESSORS_ONLN)\"");
	      printf("         returned %ld.", sysconf(_SC_NPROCESSORS_ONLN));
	      nThreads = 1;
	    }
	}
    }

  // At this point nThreads must be >= 1.
  // printf("nThreads = %d\n", nThreads);

  //  (5) THREADING:
  if(nThreads == 1)
    {
      // (5.1) NON-THREADED:

      // (a) Allocate the memory for the parameters that will be passed to the function:
      if( (params = malloc(sizeof(interp_data_t))) == NULL )
	{
	  PyErr_SetString(PyExc_MemoryError, "cmod.interpmodule.CubSplineIntrp: failed to malloc 'params'.");
	  if(x1free!=NULL) free(x1free);
	  if(x2free!=NULL) free(x2free);
	  free(ytmp);
	  return NULL;
	}
      // (b) Assign the values to the params for interpolation. 
      params->interpAcc = gsl_interp_accel_alloc(); // common for interpolation over both dimensions

      // Set the rest of parameters for interpolation in the Y dimension:
      params->inN     = nInY;
      params->outN    = nOutY;
      params->M       = nInX;
      params->inX     = x1;
      params->outX    = x2;
      params->inData  = pyinData;
      params->outData = (void*)ytmp;
      params->sx      = sInX;
      params->sy      = sInY;
      params->sytmp   = nInX;
      params->addToOutput=0;
      // (c) Do interpolation pass in Y:
      if(insize==8) interp_first_inDouble( params );
      else          interp_first_inFloat(  params );
  
      // (d) Set the params for interpolation in the X dimension:
      params->inN     = nInX;
      params->outN    = nOutX;
      params->M       = nOutY;
      params->inX     = x3;
      params->outX    = x4;
      params->inData  = (void*)ytmp;
      params->outData = pyoutData;
      params->sx      = sOutX;
      params->sy      = sOutY;
      params->addToOutput=addToOutput;

      // (e) Do interpolation pass in X and put results in output array:
      if(outsize==8) interp_second_outDouble( params );
      else           interp_second_outFloat(  params );

      // (f) Tidy up:
      gsl_interp_accel_free( params->interpAcc );
      free( params );
    }
  else
    {
      // (5.2) THREADED:

      // Important for Python:
      Py_BEGIN_ALLOW_THREADS;
 
      // (a) Allocate the memory for the 'params' vector (parameters of the 'interpolate' function):
      if((params = malloc(sizeof(interp_data_t)*nThreads))==NULL){
	PyErr_SetString(PyExc_MemoryError, "cmod.interp.gslCubSplineIntrp: Failed to malloc 'params'.\n");
	if(x1free!=NULL) free(x1free);
	if(x2free!=NULL) free(x2free);
	free(ytmp);
	return NULL;
      }

      // (b) Allocate the memory for 'thread' (the vector of threads):
      if((thread=malloc(sizeof(pthread_t)*nThreads))==NULL){
	PyErr_SetString(PyExc_MemoryError, "cmod.interp.gslCubSplineIntrp: Failed to malloc 'thread'.\n");
	free(params);
	if(x1free!=NULL) free(x1free);
	if(x2free!=NULL) free(x2free);
	free(ytmp);
	return NULL;
      }

      //     FOR THE FIRST DIMENSION:
      // (c) Determine the number of lines processed by each thread:
      //     the first 'oneMore' threads will process 'forEach+1' lines,
      //     the rest of the threads will process 'forEach' lines.
      forEach = nInX/nThreads;
      oneMore = nInX%nThreads;

      // (d) Loop over all the threads:
      for(i=0; i < nThreads; i++)
      	{
	  // (d.1) SET THE PARAMETERS to be passed to the thread function:
	  //     Allocate the accelerator - common for both dimensions:
      	  params[i].interpAcc = gsl_interp_accel_alloc();

	  //     The rest of the parameters for the interpolation in the first dimension:
      	  params[i].inN  = nInY;
      	  params[i].outN = nOutY;
      	  params[i].inX  = x1;
      	  params[i].outX = x2;
      	  params[i].sx   = sInX;
      	  params[i].sy   = sInY;
	  params[i].sytmp= nInX;
	  params[i].addToOutput=0;
      
      	  if(i < oneMore) params[i].M = forEach+1;  // the first oneMore threads
      	  else            params[i].M = forEach;    // the rest of the threads:
      	  if(i == 0){
      	    params[i].inData  = pyinData;
      	    params[i].outData = (void*)ytmp;
      	  }else{
      	    params[i].inData  = params[i-1].inData  + params[i-1].M*insize;
      	    params[i].outData = params[i-1].outData + params[i-1].M*sizeof(double);
      	  }
	  // Now all the parameters are assigned.

	  // (d.2) RUN in threads:
	  if(insize==8) pthread_create( &thread[i], NULL, (void*)interp_first_inDouble, &params[i] );
	  else          pthread_create( &thread[i], NULL, (void*)interp_first_inFloat,  &params[i] );
      	}

      // (e) Wait until all the threads are finished:
      for(i=0; i < nThreads; i++)
      pthread_join(thread[i], NULL);

      // to test the speed: call the function directly instead of in a thread
      // (only for nThreads = 1; comment out (d.2) and (e):
      //if(insize==8) interp_first_inDouble( params );
      //else          interp_first_inFloat(  params );


      //     FOR THE SECOND DIMENSION:
      // (f) Determine the number of lines processed by each thread in the second dimension:
      //     the first 'oneMore' threads will process 'forEach+1' lines,
      //     the rest of the threads will process 'forEach' lines.
      forEach = nOutY / nThreads;
      oneMore = nOutY % nThreads;

      // (g) Loop over all the threads:
      for(i=0; i < nThreads; i++)
      	{
      	  // (g.1) assign the values to the params (those that are different from the first dimension):
      	  params[i].inN     = nInX;
      	  params[i].outN    = nOutX;
      	  params[i].inX     = x3;
      	  params[i].outX    = x4;
      	  params[i].sx      = sOutX;
      	  params[i].sy      = sOutY;
	  params[i].addToOutput=addToOutput;
      	  if(i < oneMore) params[i].M = forEach+1;  // the first oneMore threads
      	  else            params[i].M = forEach;    // the rest of the threads:
      	  if(i == 0){
      	    params[i].inData  = (void*)ytmp;
      	    params[i].outData = pyoutData;
      	  }else{
      	    params[i].inData  = params[i-1].inData  + params[i-1].M*nInX*sizeof(double);
      	    params[i].outData = params[i-1].outData + params[i-1].M*nOutX*outsize;
      	  }
	  // Now all the parameters are assigned.

	  // (g.2) RUN in threads in second dimension and put results in output array 
      	  if(outsize==8) pthread_create( &thread[i], NULL, (void*)interp_second_outDouble, &params[i] );
       	  else           pthread_create( &thread[i], NULL, (void*)interp_second_outFloat,  &params[i] );
      	}

      // (h) Wait until all the threads are finished:
      for(i=0; i < nThreads; i++)
      pthread_join(thread[i], NULL);

      // to test the speed: call the function directly instead of in a thread
      // (only for nThreads = 1; comment out (g.2) and (h):
      //if(outsize==8) interp_second_outDouble( params );
      //else           interp_second_outFloat(  params );


      // (i) Tidy up:
      for(i=0; i < nThreads; i++){
	gsl_interp_accel_free( params[i].interpAcc );
	//	free( params[i].y1);
      }
      free( params );
      free( thread );

      // Important for Python:
      Py_END_ALLOW_THREADS;
    }
  // FINISHED (5)

  // (6) TIDY UP:
  if(x1free!=NULL) free(x1free);
  if(x2free!=NULL) free(x2free);
  free(ytmp);
  
  return Py_BuildValue("");   /* return None */
}
// END of gslCubSplineIntrp









static PyObject *gslCubSplineHex(PyObject *self,PyObject *args){
  /*To test this:
    xin=numpy.arange(9).astype(numpy.float64)*numpy.sqrt(3)/2.
    yin=numpy.arange(8).astype(numpy.float64)
    inarr=numpy.zeros((8,9),"f")
    shiftY=0.5
    xout=numpy.arange(64)/63.*7
    yout=numpy.arange(64)/63.*7
    outarr=numpy.zeros((64,64),"f")
    inarr[3,3]=1
    inarr[5,3]=1
    inarr[4,5]=1
    inarr[1]=1
    cmod.interp.gslCubSplineHex(inarr,yin,xin,yout,xout,shiftY,outarr,1)
    pylab.imshow(outarr,interpolation="nearest",cmap="gray")
    pylab.show()
    #If you want the shift in the other axis, just transpose the result!
   */
  
  PyArrayObject	*pyin,    *pyout;
  void          *pyinData,*pyoutData;  // void, so it can handle both float and double
  npy_intp      *pyinDims,*pyoutDims;
  PyObject *pyyin,*pyxin, *pyyout,*pyxout;

  double *x1,*x2,*x3,*x4,*ytmp;
  double *x1free=NULL,*x2free=NULL;
  double dx1,dx2,dx3,dx4;

  int i,j;
  int nInX, nInY, nOutX, nOutY;
  int sInY, sInX, sOutY, sOutX, insize, outsize;
  double *x1usr=NULL, *x2usr=NULL, *x3usr=NULL, *x4usr=NULL;
  
  int nThreads = 0;      // the number of threads (provided on input)
  pthread_t*     thread; // vector of threads
  interp_data_t* params; // vector of parameters passed to the funct. executed in a thread
                         //  defined in interpolate.h

  int forEach;           // m/nThreads = number of lines to be processed by each thread
  int oneMore;           // m%nThreads = number of threads that have to process one line more
  int last;//nunmber of lines processed by last thread.
  float shiftY; //the relative shift of alternating rows.  Set one of these to 0.5 to get hexagonal.
  // (1) PARSE THE INPUT arguments from python:
  if (!PyArg_ParseTuple (args, "O!OOOOfO!|i", 
			 &PyArray_Type, &pyin, 
			                &pyyin,
			                &pyxin,
			                &pyyout,
			                &pyxout,
			 &shiftY,
			 &PyArray_Type, &pyout, 
			                &nThreads)){
    printf("%s, line %d: parsing the arguments failed.\n", __FILE__, __LINE__);
    printf("A working example:\nxin,  yin  = numpy.arange(5).astype(TYPE),\n");
    printf("xout, yout = numpy.arange(9).astype(TYPE)/2.,\n");
    printf("mxin = numpy.zeros([5, 5], TYPE), mxout = numpy.zeros([9, 9], TYPE)\n");
    printf("TYPE must be either 'd' or 'f'; any of xin, yin, xout, yout can be 'None'.\n");
    return NULL;
  }

  // (2) EXTRACT SOME INFO AND CHECK THE INPUT

  // First check that input and output data are 2-D arrays:
  if(PyArray_NDIM(pyin)!=2 || PyArray_NDIM(pyout)!=2)
    {
      PyErr_SetString(PyExc_TypeError, "1st and 8th argument must be 2-D arrays.\n");
      return NULL;
    }

  // Extract data:
  pyinData  = (void *)PyArray_DATA(pyin);
  pyoutData = (void *)PyArray_DATA(pyout);
  pyinDims  = PyArray_DIMS(pyin);
  pyoutDims = PyArray_DIMS(pyout);
  nInY  = (int) pyinDims[0];
  nInX  = (int) pyinDims[1];
  nOutY = (int) pyoutDims[0];
  nOutX = (int) pyoutDims[1];

  switch(PyArray_TYPE(pyin)){
  case NPY_FLOAT:
    insize=4;
    break;
  default:
    insize=0;
    break;
  }
  switch(PyArray_TYPE(pyout)){
  case NPY_FLOAT:
    outsize=4;
    break;
  default:
    outsize=0;
    break;
  }

  // Check that input and output arrays are float or double:
  if(insize==0 || outsize==0){
    PyErr_SetString(PyExc_TypeError, "in and out arrays must be float32\n");
    return NULL;
  }

  // Get the strides
  sInY  = PyArray_STRIDE(pyin,0)/insize;
  sInX  = PyArray_STRIDE(pyin,1)/insize;
  sOutY = PyArray_STRIDE(pyout,0)/outsize;
  sOutX = PyArray_STRIDE(pyout,1)/outsize;

  // Further checks of the input arrays:
  //   (pyin and pyout have already been checked before up.)

  if(PyArray_Check(pyyin)){
    if(!PyArray_ISCONTIGUOUS((PyArrayObject*)pyyin) || PyArray_TYPE((PyArrayObject*)pyyin)!=NPY_DOUBLE){
      PyErr_SetString(PyExc_TypeError, "The second input parameter must be contiguous, float64\n");
      return NULL;
    }
    if(PyArray_NDIM((PyArrayObject*)pyyin)==1 && PyArray_DIM((PyArrayObject*)pyyin,0)>=nInY){//correct shape...
      x1usr=(double*)PyArray_DATA((PyArrayObject*)pyyin);
    }else{
      PyErr_SetString(PyExc_TypeError, "The second input parameter is wrong shape\n");
      return NULL;
    }
  }
  if(PyArray_Check(pyxin)){
    if(!PyArray_ISCONTIGUOUS((PyArrayObject*)pyxin) || PyArray_TYPE((PyArrayObject*)pyxin)!=NPY_DOUBLE){
      PyErr_SetString(PyExc_TypeError, "The third input parameter must be contiguous, float64\n");
      return NULL;
    }
    if(PyArray_NDIM((PyArrayObject*)pyxin)==1 && PyArray_DIM((PyArrayObject*)pyxin,0)>=nInX){//correct shape...
      x3usr=(double*)PyArray_DATA((PyArrayObject*)pyxin);
    }else{
      PyErr_SetString(PyExc_TypeError, "The third input parameter is wrong shape\n");
      return NULL;
    }
  }
  if(PyArray_Check(pyyout)){
    if(!PyArray_ISCONTIGUOUS((PyArrayObject*)pyyout) || PyArray_TYPE((PyArrayObject*)pyyout)!=NPY_DOUBLE){
      PyErr_SetString(PyExc_TypeError, "The 4th input parameter must be contiguous, float64\n");
      return NULL;
    }
    if(PyArray_NDIM((PyArrayObject*)pyyout)==1 && PyArray_DIM((PyArrayObject*)pyyout,0)>=nOutY){//correct shape...
      x2usr=(double*)PyArray_DATA((PyArrayObject*)pyyout);
    }else{
      PyErr_SetString(PyExc_TypeError, "The 4th input parameter is wrong shape\n");
      return NULL;
    }
  }
  if(PyArray_Check(pyxout)){
    if(!PyArray_ISCONTIGUOUS((PyArrayObject*)pyxout) || PyArray_TYPE((PyArrayObject*)pyxout)!=NPY_DOUBLE){
      PyErr_SetString(PyExc_TypeError, "The 5th input parameter must be contiguous, float64\n");
      return NULL;
    }
    if(PyArray_NDIM((PyArrayObject*)pyxout)==1 && PyArray_DIM((PyArrayObject*)pyxout,0)>=nOutX){//correct shape...
      x4usr=(double*)PyArray_DATA((PyArrayObject*)pyxout);
    }else{
      PyErr_SetString(PyExc_TypeError, "The 5th input parameter is wrong shape\n");
      return NULL;
    }
  }

  // (3) ALLOCATE AND (PARTIALLY) POPULATE WORKING ARRAYS
  // (This is Alastair's code, I don't touch.)
  //printf("intrp.interp2d: n1 = %d, n2 = %d, n1p = %d\n",n1,n2,n1p);
  if(x1usr==NULL || x3usr==NULL){
    x1=malloc((nInY > nInX ? nInY : nInX)*sizeof(double));
    x3=x1;
    x1free=x1;
  }
  if(x2usr==NULL || x4usr==NULL){
    x2=malloc((nOutY > nOutX ? nOutY : nOutX)*sizeof(double));
    x4=x2;
    x2free=x2;
  }
  ytmp  = malloc(nInX*nOutY*sizeof(double));
  
  dx1=1./(nInY-1.);
  dx2=dx1 * (nInY-1.) / nOutY;
  dx3=1./(nInX-1.);
  dx4=dx3*(nInX-1.)/nOutX;

  if(x1usr==NULL){
    for (i=0; i<nInY; ++i)
      x1[i] = i * dx1;
  }else
    x1=x1usr;
  if(x2usr==NULL){
    for (i=0; i<nOutY; ++i)
      x2[i] = x1[0] + (0.5+i) * dx2;
  }else
    x2=x2usr;

  if(x3usr==NULL){
    for(i=0; i<nInX; i++)
      x3[i]=i*dx3;
  }else
    x3=x3usr;
  if(x4usr==NULL){
    for(i=0; i<nOutX; i++)
      x4[i]=x3[0]+(0.5+i)*dx4;
  }else
    x4=x4usr;
  
  
  if(nThreads<1){//get a sensible default (half number of cores)
    if(nInX<40){
      nThreads = 1;
    }else{
      nThreads = sysconf(_SC_NPROCESSORS_ONLN)/2;
    }
    if(nThreads<1){
      printf("WARNING: Multiple threading was requested for .gslCubSplineHex, but it\n");
      printf("         will run in a single thread, because \"sysconf(_SC_NPROCESSORS_ONLN)\"");
      printf("         returned %ld.", sysconf(_SC_NPROCESSORS_ONLN));
      nThreads = 1;
    }
  }
 
  if(nThreads == 1){// NON-THREADED:
    if( (params = malloc(sizeof(interp_data_t))) == NULL ){
      PyErr_SetString(PyExc_MemoryError, "cmod.interpmodule.CubSplineIntrp: failed to malloc 'params'.");
      if(x1free!=NULL) free(x1free);
      if(x2free!=NULL) free(x2free);
      free(ytmp);
      return NULL;
    }
    // Assign the values to the params for interpolation. 
    params->interpAcc = gsl_interp_accel_alloc(); // common for interpolation over both dimensions
    
    // Set the rest of parameters for interpolation in the Y dimension:
    params->inN     = nInY;
    params->outN    = nOutY;
    params->M       = (nInX+1)/2;
    params->inX     = x1;
    params->outX    = x2;
    params->inData  = pyinData;
    params->outData = (void*)ytmp;
    params->sx      = sInX*2;
    params->sy      = sInY;
    params->sytmp   = nInX;
    params->sxtmp   = 2;
    // (c) Do interpolation pass in Y:
    interp_first_inFloatStep(  params );
    //and now the offset cols.
    params->M=nInX/2;
    for(i=0;i<nInY;i++)//offset the rows...
      x1[i]+=shiftY;
    params->inData+=sInX*sizeof(float);//offset the input data to the next col.
    params->outData+=sizeof(double);//offset the output data too.
    interp_first_inFloatStep(  params );
    //now do the 2nd dimension (interpolate along X).
    params->inN     = nInX;
    params->outN    = nOutX;
    params->M       = nOutY;
    params->inX     = x3;
    params->outX    = x4;
    params->inData  = (void*)ytmp;
    params->outData = pyoutData;
    params->sx      = sOutX;
    params->sy      = sOutY;
    interp_second_outFloat(  params );
    // (f) Tidy up:
    gsl_interp_accel_free( params->interpAcc );
    free( params );
  }else{// THREADED:
    Py_BEGIN_ALLOW_THREADS;
    if((params = malloc(sizeof(interp_data_t)*nThreads))==NULL){
      PyErr_SetString(PyExc_MemoryError, "cmod.interp.gslCubSplineHex: Failed to malloc 'params'.\n");
      if(x1free!=NULL) free(x1free);
      if(x2free!=NULL) free(x2free);
      free(ytmp);
      return NULL;
    }
    if((thread=malloc(sizeof(pthread_t)*nThreads))==NULL){
      PyErr_SetString(PyExc_MemoryError, "cmod.interp.gslCubSplineHex: Failed to malloc 'thread'.\n");
      free(params);
      if(x1free!=NULL) free(x1free);
      if(x2free!=NULL) free(x2free);
      free(ytmp);
      return NULL;
    }
    //     FOR THE FIRST DIMENSION:
    // (c) Determine the number of lines processed by each thread:
    //     the first 'oneMore' threads will process 'forEach+1' lines,
    //     the rest of the threads will process 'forEach' lines.
    //Threads should process an even number of lines (except for the last one).
    
    forEach = (nInX/nThreads+1)&(~1);//make it an even number
    last=nInX-(nThreads-1)*forEach;//last thread does this many rows
    //oneMore = nInX%nThreads;
    // (d) Loop over all the threads:
    for(i=0; i < nThreads; i++){
      // (d.1) SET THE PARAMETERS to be passed to the thread function:
      //     Allocate the accelerator - common for both dimensions:
      params[i].interpAcc = gsl_interp_accel_alloc();
      params[i].inN  = nInY;
      params[i].outN = nOutY;
      params[i].inX  = x1;
      params[i].outX = x2;
      params[i].sx   = sInX*2;
      params[i].sy   = sInY;
      params[i].sytmp= nInX;
      params[i].sxtmp=2;
      if(i<nThreads-1){
	params[i].M=forEach/2;
      }	else{
	params[i].M=(last+1)/2;
      }
      //if(i < oneMore){
      //	params[i].M = (forEach+1+1)/2;  // the first oneMore threads
      //}else{
      //params[i].M = (forEach+1)/2;    // the rest of the threads:
      //}
      if(i == 0){
	params[i].inData  = pyinData;
	params[i].outData = (void*)ytmp;
      }else{//note: insize=sizeof(float)
	params[i].inData  = params[i-1].inData  + forEach*insize;
	params[i].outData = params[i-1].outData + forEach*sizeof(double);
      }

      // Now all the parameters are assigned.
      pthread_create( &thread[i], NULL, (void*)interp_first_inFloatStep,  &params[i] );
    }
    for(i=0; i < nThreads; i++)
      pthread_join(thread[i], NULL);
    for(i=0;i<nInY;i++)
      x1[i]+=shiftY;
    for(i=0; i < nThreads; i++){
      if(i<nThreads-1)      params[i].M=forEach/2;
      else               params[i].M=last/2;
      params[i].inData+=sInX*sizeof(float);
      params[i].outData+=sizeof(double);
      pthread_create( &thread[i], NULL, (void*)interp_first_inFloatStep,  &params[i] );
    }
    for(i=0; i < nThreads; i++)
      pthread_join(thread[i], NULL);
    //     FOR THE SECOND DIMENSION:
    //     the first 'oneMore' threads will process 'forEach+1' lines,
    //     the rest of the threads will process 'forEach' lines.
    forEach = nOutY / nThreads;
    oneMore = nOutY % nThreads;
    for(i=0; i < nThreads; i++){
      params[i].inN     = nInX;
      params[i].outN    = nOutX;
      params[i].inX     = x3;
      params[i].outX    = x4;
      params[i].sx      = sOutX;
      params[i].sy      = sOutY;
      if(i < oneMore) params[i].M = forEach+1;  // the first oneMore threads
      else            params[i].M = forEach;    // the rest of the threads:
      if(i == 0){
	params[i].inData  = (void*)ytmp;
	params[i].outData = pyoutData;
      }else{
	params[i].inData  = params[i-1].inData  + params[i-1].M*nInX*sizeof(double);
	params[i].outData = params[i-1].outData + params[i-1].M*nOutX*outsize;
      }
      pthread_create( &thread[i], NULL, (void*)interp_second_outFloat,  &params[i] );
    }
    for(i=0; i < nThreads; i++)
      pthread_join(thread[i], NULL);
    // (i) Tidy up:
    for(i=0; i < nThreads; i++){
      gsl_interp_accel_free( params[i].interpAcc );
    }
    free( params );
    free( thread );
    Py_END_ALLOW_THREADS;
  }
  if(x1free!=NULL) free(x1free);
  if(x2free!=NULL) free(x2free);
  free(ytmp);
  return Py_BuildValue("");   /* return None */
}
// END of gslCubSplineHex






// Stolen from the Numerical Recipies:
void NRbcucof(float y[], float y1[], float y2[], float y12[], float d1, float d2,	float **c)
{
  static int wt[16][16]=
    {{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0},
     {-3,0,0,3,0,0,0,0,-2,0,0,-1,0,0,0,0},
     {2,0,0,-2,0,0,0,0,1,0,0,1,0,0,0,0},
     {0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0},
     {0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1},
     {0,0,0,0,2,0,0,-2,0,0,0,0,1,0,0,1},
     {-3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0},
     {9,-9,9,-9,6,3,-3,-6,6,-6,-3,3,4,2,1,2},
     {-6,6,-6,6,-4,-2,2,4,-3,3,3,-3,-2,-1,-1,-2},
     {2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0},
     {-6,6,-6,6,-3,-3,3,3,-4,4,2,-2,-2,-2,-1,-1},
     {4,-4,4,-4,2,2,-2,-2,2,-2,-2,2,1,1,1,1}};
  int l,k,j,i;
  float xx,d1d2,cl[16],x[16];

  d1d2=d1*d2;
  for (i=1;i<=4;i++) {
    x[i-1]=y[i];
    x[i+3]=y1[i]*d1;
    x[i+7]=y2[i]*d2;
    x[i+11]=y12[i]*d1d2;
  }
  for (i=0;i<=15;i++) {
    xx=0.0;
    for (k=0;k<=15;k++) xx += wt[i][k]*x[k];
    cl[i]=xx;
  }
  l=0;
  for (i=1;i<=4;i++)
    for (j=1;j<=4;j++) c[i][j]=cl[l++];
}
/* (C) Copr. 1986-92 Numerical Recipes Software 21%. */

// Stolen from Numerical Recipies:
float *NRvector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;

	v=(float *)malloc((size_t) ((nh-nl+1+1)*sizeof(float)));
	if (!v) printf("allocation failure in vector()");
	return v-nl+1;
}

// Stolen from Numerical Recipies:
float **NRmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((size_t)((nrow+1)*sizeof(float*)));
	if (!m) printf("allocation failure 1 in matrix()");
	m += 1;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(float *) malloc((size_t)((nrow*ncol+1)*sizeof(float)));
	if (!m[nrl]) printf("allocation failure 2 in matrix()");
	m[nrl] += 1;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

// Stolen from Numerical Recipies:
void NRfree_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
	free((char*) (v+nl-1));
}

// Stolen from Numerical Recipies:
void NRfree_matrix(float **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by matrix() */
{
	free((char*) (m[nrl]+ncl-1));
	free((char*) (m+nrl-1));
}


static PyObject *interp_bicubicinterp(PyObject* self, PyObject* args)
{
  PyArrayObject	*pymxin,*pymxout,*pyxin,*pyyin,*pyxyin;
  int		i,j,k,l,m,di,dj,nx,ny,nxout,nyout;
  int diy,djy,dix,djx,dixy,djxy,diout,djout;
  float		*zin,*dy,*dx,*dxy;
  float		**c;
  int iilast=-1,jjlast=-1,ii=-1,jj=-1,iin,jjn;
  float iif,jjf,xstep,ystep;
  float yinterppos,xinterppos,ansy;

  printf("WARNING-REMINDER: interp_bicubicinterp is still using bcucof from Numerical Recipies.\n    TODO: remove bcucoff from bicubicinterp!\n");

  if (!PyArg_ParseTuple (args, "O!O!O!O!O!", &PyArray_Type ,&pymxin, 
			 &PyArray_Type ,&pyyin, //y gradient of function
			 &PyArray_Type ,&pyxin, //x gradient of function
			 &PyArray_Type ,&pyxyin, //xy gradient of func.
			 &PyArray_Type ,&pymxout)) {
    return NULL;
  }

  // get input array dimensions 
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

  if(pymxin->descr->type_num!=NPY_FLOAT ||
     pyyin->descr->type_num!=NPY_FLOAT ||
     pyxin->descr->type_num!=NPY_FLOAT ||
     pyxyin->descr->type_num!=NPY_FLOAT ||
     pymxout->descr->type_num!=NPY_FLOAT 
     ){
    printf("matrixes for bicubic interpolation must be float32.\n");
    return NULL;
  }

  c=NRmatrix(1,4,1,4);
  zin=NRvector(1,4);
  dy=NRvector(1,4);
  dx=NRvector(1,4);
  dxy=NRvector(1,4);

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
	NRbcucof(zin,dy,dx,dxy,1.,1.,c);
	if(isnan(c[1][1])){
	  printf("Got c as nan for ii iilast jj jjlast %d %d %d %d\nzin dy dx dxy\n",ii,iilast,jj,jjlast);
	  for(k=1; k<=4; k++){
	    printf("%g %g %g %g\n",zin[k],dy[k],dx[k],dxy[k]);
	  }
	}
      }
      yinterppos=iif-(float)ii;
      xinterppos=jjf-(float)jj;
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
  NRfree_matrix(c,1,4,1,4);
  NRfree_vector(zin,1,4);
  NRfree_vector(dy,1,4);
  NRfree_vector(dx,1,4);
  NRfree_vector(dxy,1,4);
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
  // check the value of shift (shift must be between 0 and 1):
  if( xshift > 1 || yshift > 1 || xshift < 0 || yshift < 0 )
    {
      // set exception, print the faulty values and return
      PyErr_SetString(PyExc_ValueError, "xshift and yshift must be between 0 and 1. If you see this message, the function interp_linearshift must be extended for shifts outside the [0,1] interval.");
      printf("xshift = %f, yshift = %f\n", xshift, yshift);
      return NULL;
    }
  
  ny=outarr->dimensions[0];
  nx=outarr->dimensions[1];
  di=inarr->strides[0];
  dj=inarr->strides[1];
  dii=outarr->strides[0];
  djj=outarr->strides[1];
  if(inarr->dimensions[0]<ny+1 || inarr->dimensions[1]<nx+1){
    printf("WARNING (intermodule.c): input array must be at least 1 bigger than output array in each dimension\nContinuing, but output array won't be filled completely. (%d<%d+1, %d<%d+1)\n", (int)inarr->dimensions[0], ny, (int)inarr->dimensions[1], nx);
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
static PyObject *interp_linearinterp(PyObject* self, PyObject* args)
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
					{"bicubicinterp", interp_bicubicinterp, METH_VARARGS}, 
					{"linearshift",interp_linearshift,METH_VARARGS},
					{"linearinterp", interp_linearinterp, METH_VARARGS}, 
					{"gslPeriodicCubSplineInterp",gslPeriodicCubSplineIntrp,METH_VARARGS},
					{"gslCubSplineInterp",gslCubSplineIntrp,METH_VARARGS},
					{"gslCubSplineHex",gslCubSplineHex,METH_VARARGS},
					//
					// Not used anymore; commented out by UB on July 30th 2012:
					//
					//{"mxinterp", interp_mxinterp, METH_VARARGS}, 
					//{"gslLinInterp",gslLinIntrp,METH_VARARGS},
					//{"gslCubSplineInterpOrig",cubIntrpOrig,METH_VARARGS},
					//
					{NULL, NULL} };

/* initialisation - register the methods with the Python interpreter */

void initinterp(void)
{
	(void) Py_InitModule("interp", interp_methods);
	import_array();
}

























/* ==================================================================================== */
/* ==================================================================================== */
/*                                                                                      */
/*                OLD FUNCTIONS THAT ARE NOT USED ANY MORE                              */
/*                                                                                      */
/* ==================================================================================== */
/* ==================================================================================== */
// 
// - splintagb, splin3, splin4 and fastmxinterp are used only by interp_mxinterp
// - interp_mxinterp was replaced by gslCubSplineInterp
// - gslLinIntrp and cubIntrpOrig are not used anywhere in the simulation
// 
// Commented out and moved below here by UB on July 30th 2012
//


/* #define NRANSI */
/* #include "nrutil.h" */

/* static PyObject *cubIntrpOrig(PyObject* self, PyObject* args) */
/* { */
/* 	PyArrayObject	*pyin,*pyout; */
/* 	npy_intp		*pyinDims,*pyoutDims; */
/* 	double			*pyinData,*pyoutData; */

/* 	int		zpad=10;   /\* no of pixels to zero-pad input by *\/ */

/* 	double		*x1,*x2,*y1,*ytmp; */
/* 	double		dx1,dx2; */
/* 	int			i,j,n1,n2,n1p; */

/* 	gsl_interp 	*interpObj; */
/* 	gsl_interp_accel *interpAcc ; */
/* 	size_t  ninterp; */


/* 	/\* handle input arrays *\/ */
/* 	if (!PyArg_ParseTuple (args, "O!O!", &PyArray_Type, &pyin, &PyArray_Type, &pyout)) */
/* 		return NULL; */
	
/* 	pyinDims=PyArray_DIMS(pyin); */
/* 	pyoutDims=PyArray_DIMS(pyout); */
/* 	n1=(int) pyinDims[0]; */
/* 	n2=(int) pyoutDims[0]; */
/* 	pyinData = (double *)PyArray_DATA(pyin); */
/* 	pyoutData = (double *)PyArray_DATA(pyout); */


/* 	/\* allocate and (partially) populate working arrays *\/ */
/* 	n1p=n1+2*zpad; */
/* 	//printf("intrp.interp2d: n1 = %d, n2 = %d, n1p = %d\n",n1,n2,n1p); */
/* 	x1=calloc(n1p, sizeof(double)); */
/* 	x2=calloc(n2, sizeof(double)); */
/* 	y1=calloc(n1p, sizeof(double)); */
/* 	ytmp=calloc(n1*n2, sizeof(double)); */

/* 	dx1=1./(double)(n1p-1); */
/* 	dx2=dx1 * (double)(n1-1) / (double)(n2); */

/* 	for (i=0; i<n1p; ++i) { */
/* 		x1[i] = (double)i * dx1; */
/* 	} */
/* 	for (i=0; i<n2; ++i) { */
/* 		x2[i] = x1[zpad] + (0.5+(double)i) * dx2; */
/* 	} */
/* 	for (i=0; i<zpad; ++i) { */
/* 		y1[i] = 0.; */
/* 		y1[zpad+n1+i]=0.; */
/* 	} */
	

/* 	/\* allocate interpolation stuff (use periodic cspline to make sure influence functions are identical) *\/ */
/* 	ninterp=(size_t)n1p; */
/* 	interpObj=gsl_interp_alloc(gsl_interp_cspline_periodic, ninterp); */
/* 	interpAcc=gsl_interp_accel_alloc(); */

	
/* 	/\* do interpolation pass in first dimension *\/ */
/* 	for (j=0; j<n1; ++j) { */
/* 		for (i=0; i<n1; ++i) { */
/* 			y1[zpad+i]=pyinData[i*n1+j]; */
/* 		} */
/* 		gsl_interp_init(interpObj, x1, y1, ninterp); */

/* 		for (i=0; i<n2; ++i) { */
/* 			ytmp[j*n2+i]=gsl_interp_eval(interpObj, x1, y1, x2[i], interpAcc); */
/* 		}	 */
/* 	} */
	
	
/* 	/\* do interpolation pass in second dimension and put results in output array *\/ */
/* 	for (j=0; j<n2; ++j) { */
/* 		for (i=0; i<n1; ++i) { */
/* 			y1[zpad+i]=ytmp[i*n2+j]; */
/* 		} */
/* 		gsl_interp_init(interpObj, x1, y1, ninterp); */
/* 		for (i=0; i<n2; ++i) { */
/* 			pyoutData[j*n2+i]=gsl_interp_eval(interpObj, x1, y1, x2[i], interpAcc); */
/* 		}	 */
/* 	}	 */
	
	
/* 	/\* tidy up *\/ */
/* 	gsl_interp_accel_free(interpAcc); */
/* 	gsl_interp_free(interpObj); */
/* 	free(x1); */
/* 	free(x2); */
/* 	free(y1); */
/* 	free(ytmp); */
	
/* 	return Py_BuildValue("");   /\* return None *\/ */
/* } */
/* // END of cubIntrpOrig */

/* /\* Linear interpolation in 2D (requires 64-bit arrays) *\/ */
/* static PyObject *gslLinIntrp(PyObject *self,PyObject *args){ */
/* 	PyArrayObject	*pyin,*pyout; */
/* 	npy_intp *pyinDims,*pyoutDims; */
/* 	double *pyinData,*pyoutData; */

/* 	double *x1,*x2,*y1,*ytmp; */
/* 	int i,j,n1,n2; */

/* 	gsl_interp *interpObj; */
/* 	gsl_interp_accel *interpAcc ; */
/* 	size_t  ninterp; */


/* 	/\* handle input arrays *\/ */
/* 	if (!PyArg_ParseTuple (args, "O!O!", &PyArray_Type, &pyin, &PyArray_Type, &pyout)) */
/* 		return NULL; */
	
/* 	pyinDims=PyArray_DIMS(pyin); */
/* 	pyoutDims=PyArray_DIMS(pyout); */
/* 	n1=(int) pyinDims[0]; */
/* 	n2=(int) pyoutDims[0]; */
/* 	pyinData = (double *)PyArray_DATA(pyin); */
/* 	pyoutData = (double *)PyArray_DATA(pyout); */


/* 	/\* allocate populate working arrays *\/ */
/* 	x1=calloc(n1, sizeof(double)); */
/* 	x2=calloc(n2, sizeof(double)); */
/* 	y1=calloc(n1, sizeof(double)); */
/* 	ytmp=calloc(n1*n2, sizeof(double)); */

/* 	for (i=0; i<n1; ++i) { */
/* 		x1[i] = (double)i / (double)(n1-1); */
/* 	} */
/* 	for (i=0; i<n2; ++i) { */
/* 		x2[i] = (0.5 + (double)i) / (double)n2; */
/* 	} */


/* 	/\* allocate interpolation stuff *\/ */
/* 	ninterp=(size_t)n1; */
/* 	interpObj=gsl_interp_alloc(gsl_interp_linear, ninterp); */
/* 	interpAcc=gsl_interp_accel_alloc(); */

/* 	/\* do interpolation pass in first dimension *\/ */
/* 	for (j=0; j<n1; ++j) { */
/* 		for (i=0; i<n1; ++i) { */
/* 			y1[i]=pyinData[i*n1+j]; */
/* 		} */
/* 		gsl_interp_init(interpObj, x1, y1, ninterp); */

/* 		for (i=0; i<n2; ++i) { */
/* 			ytmp[j*n2+i]=gsl_interp_eval(interpObj, x1, y1, x2[i], interpAcc); */
/* 		}	 */
/* 	} */
	

/* 	/\* do interpolation pass in second dimension and put results in output array *\/ */
/* 	for (j=0; j<n2; ++j) { */
/* 		for (i=0; i<n1; ++i) { */
/* 			y1[i]=ytmp[i*n2+j]; */
/* 		} */
/* 		gsl_interp_init(interpObj, x1, y1, ninterp); */
/* 		for (i=0; i<n2; ++i) { */
/* 			pyoutData[j*n2+i]=gsl_interp_eval(interpObj, x1, y1, x2[i], interpAcc); */
/* 		}	 */
/* 	}	 */
	
	
/* 	/\* tidy up *\/ */
/* 	//printf("intrp.interp2d: Tidying up\n"); */
/* 	gsl_interp_accel_free(interpAcc); */
/* 	gsl_interp_free(interpObj); */
/* 	free(x1); */
/* 	free(x2); */
/* 	free(y1); */
/* 	free(ytmp); */
	
/* 	return Py_BuildValue("");   /\* return None *\/ */
/* } */
/* // END of gslLinIntrp */

/* //splint - agb version to avoid exit if h==0. */
/* int splintagb(float xa[], float ya[], float y2a[], int n, float x, float *y) */
/* { */
/* 	void nrerror(char error_text[]); */
/* 	int klo,khi,k; */
/* 	float h,b,a; */

/* 	klo=1; */
/* 	khi=n; */
/* 	while (khi-klo > 1) { */
/* 		k=(khi+klo) >> 1; */
/* 		if (xa[k] > x) khi=k; */
/* 		else klo=k; */
/* 	} */
/* 	h=xa[khi]-xa[klo]; */
/* 	if (h == 0.0){ */
/* 	  printf("Bad xa input to routine splint in interpmodule\n"); */
/* 	  printf("SIMULATION NEEDS TO BE HALTED\n"); */
/* 	  return 1; */
/* 	} */
/* 	a=(xa[khi]-x)/h; */
/* 	b=(x-xa[klo])/h; */
/* 	*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0; */
/* 	return 0; */
/* } */


/* /\* (C) Copr. 1986-92 Numerical Recipes Software 21%. *\/ */
/* void splin3(float x1a[], float x2a[], float **ya, float **y2a, int m, int n, int mout, int nout, */
/* 	float x1[], float x2[], float **mxout) */
/* { */
/* 	int i,j; */
/* 	float *ytmp,*yytmp,*y; */
/* 	void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]); */
/* 	void splint(float xa[], float ya[], float y2a[], int n, float x, float *y); */

/* 	ytmp=vector(1,m); */
/* 	yytmp=vector(1,m); */
/* 	y=calloc(1,sizeof(float));  */

/* 	for (i=1;i<=nout;i++){ */

/* 	  for (j=1;j<=m;j++) */
/* 	    //if(splintagb(x2a,ya[j],y2a[j],n,x2[i],&yytmp[j])) */
/* 	    // return 1; */
/* 	    splint(x2a,ya[j],y2a[j],n,x2[i],&yytmp[j]); */
/* 	  spline(x1a,yytmp,m,1.0e30,1.0e30,ytmp); */

/* 	  for (j=1;j<=mout;j++){ */
/* 	    //if(splintagb(x1a,yytmp,ytmp,m,x1[j],y)) */
/* 	    //return 1; */
/* 	    splint(x1a,yytmp,ytmp,m,x1[j],y); */
/* 	    mxout[i][j]=*y;  */
/* 	  } */
/* 	} */

/* 	free(y); */
/* 	free_vector(yytmp,1,m); */
/* 	free_vector(ytmp,1,m); */
/* 	//return 0; */
/* } */
/* //splin4 - agb version to avoid excessive data copying */
/* void splin4(float x1a[], float x2a[], float *ya,int di, float **y2a, int m, int n, int mout, int nout,float x1[], float x2[], float *mxout,int ddi,int ddj) */
/* { */
/*   int i,j; */
/*   float *ytmp,*yytmp; */
/*   float y=0; */
/*   void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]); */
/*   void splint(float xa[], float ya[], float y2a[], int n, float x, float *y); */
  
/*   ytmp=vector(1,m); */
/*   yytmp=vector(1,m); */
/*   for (i=0;i<nout;i++){ */
/*     for (j=0;j<m;j++) */
/*       //if(splintagb(x2a-1,&ya[j*di-1],y2a[j+1],n,x2[i],&yytmp[j+1])) */
/*       //return 1; */
/*       splint(x2a-1,&ya[j*di-1],y2a[j+1],n,x2[i],&yytmp[j+1]); */

/*     spline(x1a-1,yytmp,m,1.0e30,1.0e30,ytmp); */
/*     for (j=0;j<mout;j++){ */
/*       //if(splintagb(x1a-1,yytmp,ytmp,m,x1[j],&y)) */
/*       //return 1; */
/*       splint(x1a-1,yytmp,ytmp,m,x1[j],&y); */
/*       mxout[i*ddj+ddi*j]=y;  */
/*     } */
/*   } */
/*   free_vector(yytmp,1,m); */
/*   free_vector(ytmp,1,m); */
/*   //return 0; */
/* } */

/* #undef NRANSI */


/* //static PyObject *interp_fastmxinterp(PyObject *self,PyObject *args){ */
/* static PyObject *fastmxinterp(PyArrayObject *pymxin,PyArrayObject *pyyin,PyArrayObject *pyxin,PyArrayObject *pyyout,PyArrayObject *pyxout,PyArrayObject *pymxout){ */
/*   //does the matrix interpolation without copying.  Much more strict on input types.  Should return same result as mxinterp. */
/*   //Actually not really much faster, though should be very slightly.  Better use of memory too. */
/*   //PyArrayObject *pymxin,*pyyin,*pyxin,*pyyout,*pyxout,*pymxout; */
/*   int m,n,di,nyout,nxout,j; */
/*   float **deriv; */
/*   //if (!PyArg_ParseTuple (args, "O!O!O!O!O!O!", &PyArray_Type ,&pymxin, &PyArray_Type ,&pyyin, &PyArray_Type ,&pyxin, &PyArray_Type ,&pyyout, &PyArray_Type ,&pyxout, &PyArray_Type ,&pymxout)) { */
/*   //return NULL; */
/*   //} */
/*   if(pymxin->descr->type_num==NPY_FLOAT &&  */
/*      pyyin->descr->type_num==NPY_FLOAT &&  */
/*      pyxin->descr->type_num==NPY_FLOAT &&  */
/*      pyyout->descr->type_num==NPY_FLOAT &&  */
/*      pyxout->descr->type_num==NPY_FLOAT &&  */
/*      pymxout->descr->type_num==NPY_FLOAT){ */
/*     //all float so okay. */
/*     //Check dimensionality */
/*     if(pyyin->nd==1 && pyxin->nd==1 && pyyout->nd==1 && pyxout->nd==1 && pymxin->nd==2 && pymxout->nd==2){ */
/*       //okay */
/*       //Check that all are contiguous. */
/*       if(pyyin->strides[0]==sizeof(float) && */
/* 	 pyxin->strides[0]==sizeof(float) && */
/* 	 pyyout->strides[0]==sizeof(float) && */
/* 	 pyxout->strides[0]==sizeof(float) && */
/* 	 pymxin->strides[1]==sizeof(float)){ */
/* 	//okay */
/* 	//Now do the calc. */
/* 	//First, a handwritten of splie2. */
/* 	//splie2(yin,xin,mxin,ny,nx,deriv); */
/* 	//void splie2(float x1a[], float x2a[], float **ya, int m, int n, float **y2a) */
/* 	//x1a==yin, x2a==xin, ya==mxin, m=ny, n=nx, y2a==deriv */
/* 	m=pymxin->dimensions[0]; */
/* 	n=pymxin->dimensions[1]; */
/* 	di=pymxin->strides[0]; */
/* 	nyout=pyyout->dimensions[0]; */
/* 	nxout=pyxout->dimensions[0]; */
/* 	deriv=matrix(1,m,1,n); */

/* 	for (j=0;j<m;j++) */
/* 	  spline((float*)(pyxin->data-sizeof(float)),(float*)(pymxin->data+j*di-sizeof(float)),n,1.0e30,1.0e30,deriv[j+1]); */
/*   	splin4((float*)(pyyin->data),(float*)(pyxin->data),(float*)(pymxin->data),di/sizeof(float),deriv,m,n,nyout,nxout,(float*)(pyyout->data),(float*)(pyxout->data),(float*)(pymxout->data),pymxout->strides[0]/sizeof(float),pymxout->strides[1]/sizeof(float)); */
/* 	  //return NULL; */
/* 	free_matrix(deriv,1,m,1,n); */

/*       }else{ */
/* 	//printf("fastmxinterp: contiguous condition not met\n"); */
/* 	return NULL; */
/*       } */

/*     }else{ */
/*       //printf("fastmxinterp: Array dimensionality wrong\n"); */
/*       return NULL; */
/*     } */
       

/*   }else{ */
/*     //printf("fastmxinterp inputs must be float32\n"); */
/*     return NULL; */
/*   } */
  
/*   return Py_BuildValue("");  */


/* } */
/* /\* ============================================================================= *\/ */



/* /\* Matrix interpolation function. Takes Python array and vectors of x,y coords for interpolated points *\/ */
/* /\* Returns Python array of interp vales *\/ */
/* //This is a bicubic spline interpolation. */
/* static PyObject *interp_mxinterp(PyObject* self, PyObject* args) */
/* { */
/* 	PyArrayObject	*pymxin,*pymxout,*pyxin,*pyyin,*pyxout,*pyyout; */
/* 	PyObject *tmp; */
/* 	int		i,j,di,dj,nx,ny,nxout,nyout; */
/* 	float		*xin,*yin,*xout,*yout; */
/* 	float		**mxin,**mxout,**deriv; */


/* 	/\* z=calloc(1,sizeof(float)); *\/ */

/* 	if (!PyArg_ParseTuple (args, "O!O!O!O!O!O!", &PyArray_Type ,&pymxin, &PyArray_Type ,&pyyin, &PyArray_Type ,&pyxin, &PyArray_Type ,&pyyout, &PyArray_Type ,&pyxout, &PyArray_Type ,&pymxout)) { */
/* 		return NULL; */
/* 	} */
/* 	//First see if we can do the non-copy version (needs float inputs). */
/* 	if((tmp=fastmxinterp(pymxin,pyyin,pyxin,pyyout,pyxout,pymxout))!=NULL){ */
/* 	  return tmp; */
/* 	} */

/* /\* get input array dimensions *\/ */
/* 	ny=pymxin->dimensions[0]; */
/* 	nx=pymxin->dimensions[1]; */
/* 	di=pymxin->strides[0]; */
/* 	dj=pymxin->strides[1]; */

/* /\* copy to a C array *\/ */
/* 	mxin=matrix(1,ny,1,nx); */
/* 	deriv=matrix(1,ny,1,nx); */
/* 	if(pymxin->descr->type_num==NPY_DOUBLE){ */
/* 	    /\* mxin=alloc2d_float(ny,nx); *\/ */
/* 	    for(i=0;i<ny;++i){ */
/* 		for(j=0;j<nx;++j){ */
/* 		    mxin[i+1][j+1] = (float)(*(double *)(pymxin->data + i*di + j*dj)); */
		    
/* 		} */
/* 	    } */
/* 	}else if(pymxin->descr->type_num==NPY_FLOAT){ */
/* 	    for(i=0;i<ny;++i){ */
/* 		for(j=0;j<nx;++j){ */
/* 		    mxin[i+1][j+1] = (float)(*(float *)(pymxin->data + i*di + j*dj)); */
/* 		} */
/* 	    }  */
/* 	}else{ */
/* 	    printf("ERROR in interpmodule - input array not double or float\n"); */
/* 	    return NULL; */
/* 	} */
/* /\* copy python x,y coord */
/*    vectors to C vectors *\/ */
/* 	if(pyyin->descr->type_num!=NPY_DOUBLE || pyxin->descr->type_num!=NPY_DOUBLE || pyyout->descr->type_num!=NPY_DOUBLE || pyxout->descr->type_num!=NPY_DOUBLE){ */
/* 	  printf("mxinterp: pyyin/xin/etc should be float64 unless trying to do the non-copy method (when they should be float32) - if you are trying to do the non-copy method, then some other of your inputs are currently wrong\n"); */
/* 	  return NULL; */
/* 	} */
/* 	yin=vector(1,ny); */
/* 	di=pyyin->strides[0]; */
/* 	for(i=0;i<ny;++i) yin[i+1] = (float)(*(double *)(pyyin->data + i*di)); */

/* 	xin=vector(1,nx); */
/* 	dj=pyxin->strides[0]; */
/* 	for(j=0;j<nx;++j) xin[j+1] = (float)(*(double *)(pyxin->data + j*dj)); */



/* 	nyout=pyyout->dimensions[0]; */
/* 	di=pyyout->strides[0]; */
/* 	yout=vector(1,nyout); */
/* 	for(i=0;i<nyout;++i) yout[i+1] = (float)(*(double *)(pyyout->data + i*di)); */

/* 	nxout=pyxout->dimensions[0]; */
/* 	dj=pyxout->strides[0]; */
/* 	xout=vector(1,nxout); */
/* 	for(j=0;j<nxout;++j) xout[j+1] = (float)(*(double *)(pyxout->data + j*dj)); */

/* /\* do the interpolation *\/ */
/* 	splie2(yin,xin,mxin,ny,nx,deriv); */

/* 	mxout=matrix(1,nyout,1,nxout); */

/*  /\* quick version of splin2 for interpd points regular grid *\/ */
/*   	splin3(yin,xin,mxin,deriv,ny,nx,nyout,nxout,yout,xout,mxout); */
/* 	  //return NULL; */

/* /\* populate output Numpy array with the interpolated vals *\/ */

/* 	dj=pymxout->strides[0]; */
/* 	di=pymxout->strides[1]; */

/* 	if(pymxout->descr->type_num==NPY_DOUBLE){ */
/* 	    for(i=0;i<nyout;++i){ */
/* 		for(j=0;j<nxout;++j){ */
/* 		    *(double *)(pymxout->data+i*di+j*dj) = (double)mxout[i+1][j+1]; 	/\* !!!! this fn transposes the input mx - why ? !!! *\/ */
/* 		} */
/* 	    } */
/* 	}else if(pymxout->descr->type_num==NPY_FLOAT){ */
/* 	    for(i=0;i<nyout;++i){ */
/* 		for(j=0;j<nxout;++j){ */
/* 		    *(float *)(pymxout->data+i*di+j*dj) = (float)mxout[i+1][j+1]; 	/\* !!!! this fn transposes the input mx - why ? !!! *\/ */
/* 		} */
/* 	    } */
/* 	}else{ */
/* 	    printf("ERROR in interpmodule - output array should be double or float\n"); */
/* 	    return NULL; */
/* 	} */

/* 	free_vector(yin,1,ny); */
/* 	free_vector(xin,1,nx); */
/* 	free_vector(yout,1,nyout); */
/* 	free_vector(xout,1,nxout); */

/* 	free_matrix(mxin,1,ny,1,nx); */
/* 	free_matrix(mxout,1,nyout,1,nxout); */
/* 	free_matrix(deriv,1,ny,1,nx); */

/* /\* Return the Numpy array *\/ */

/* 	return Py_BuildValue("");  */
/* } */

