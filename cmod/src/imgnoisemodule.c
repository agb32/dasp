
/* Numpy extension to add Poisson noise to an image in a Numpy array */
/* Uses Numerical Recipes `Poidev' random number generator */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nr.h"
#include "nrutil.h"
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
static PyObject *imgnoise_shot(self,args)
	PyObject *self, *args;

{
	PyArrayObject	*imgin,*imgout;
	long		*idum;
	int		i,j,nd,di,dj,diout,djout,size,dims[2];
	float		val,mean;

 
	if (!PyArg_ParseTuple (args, "O!O!", &PyArray_Type, &imgin, &PyArray_Type, &imgout))
		return NULL;


	idum=calloc(1,sizeof(long));
	*idum=rand(); 
	/* printf("...%d %p %p\n",*idum,&idum,idum); */


/* get input array dimensions */

	nd=imgin->nd;
	dims[0]=imgin->dimensions[0];
	dims[1]=imgin->dimensions[1];
	di=imgin->strides[0];
	dj=imgin->strides[1];
	/* printf("%d %d %d %d %d \n",nd,dims[0],dims[1],di,dj); */
	diout=imgout->strides[0];
	djout=imgout->strides[1];



/* populate output image with random Poisson deviates with mean = input image pixel value */

	for(i=0;i<dims[0];++i){
	  for(j=0;j<dims[1];++j){
		mean = (float)(*(double *)(imgin->data + i*di + j*dj));
		*(double *)(imgout->data+i*diout+j*djout) = (double)poidev(mean,idum);

	  }
	}

	free(idum);

	return Py_BuildValue("");


}



static PyObject *imgnoise_shot1d(self,args)
	PyObject *self, *args;

{
	PyArrayObject	*imgin,*imgout;
	long		*idum;
	int		i,j,nd,di,dj,diout,djout,size,dims[1];
	float		val,mean;

 
	if (!PyArg_ParseTuple (args, "O!O!", &PyArray_Type, &imgin, &PyArray_Type, &imgout))
		return NULL;


	idum=calloc(1,sizeof(long));
	*idum=rand(); 
	/* printf("...%d %p %p\n",*idum,&idum,idum); */


/* get input array dimensions */

	nd=imgin->nd;
	dims[0]=imgin->dimensions[0];
	di=imgin->strides[0];
	diout=imgout->strides[0];


/* populate output image with random Poisson deviates with mean = input image pixel value */

	for(i=0;i<dims[0];++i){
		mean = (float)(*(double *)(imgin->data + i*di));
		*(double *)(imgout->data+i*diout) = (double)poidev(mean,idum);
	}

	free(idum);

	return Py_BuildValue("");


}
#ifndef Py_RETURN_NONE
#define Py_RETURN_NONE return Py_INCREF(Py_None), Py_None
#endif

static PyObject *imgnoise_seed(PyObject *self,PyObject *args){
    int theseed=0;
    if(!PyArg_ParseTuple(args,"i",&theseed))
	return NULL;
    if(theseed==0){
	printf("Warning, using a seed of zero will use current time\n");
	srand(time(0));
    }else
	srand(theseed);
    Py_RETURN_NONE;
}

/* define a methods table for the module */

static PyMethodDef imgnoise_methods[] = 	{
					{"shot", imgnoise_shot, METH_VARARGS}, 
					{"shot1d", imgnoise_shot1d, METH_VARARGS}, 
					{"seed",imgnoise_seed,METH_VARARGS},
					{NULL, NULL} };


/* initialisation - register the methods with the Python interpreter */

void initimgnoise()
{
	(void) Py_InitModule("imgnoise", imgnoise_methods);
	import_array();
}

