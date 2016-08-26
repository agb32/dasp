/*
#dasp, the Durham Adaptive optics Simulation Platform.
#Copyright (C) 2004-2016 Alastair Basden and Durham University.

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU Affero General Public License as
#published by the Free Software Foundation, either version 3 of the
#License, or (at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU Affero General Public License for more details.

#You should have received a copy of the GNU Affero General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/* Extension to add Poisson noise to an image in a Numpy array */
/* This function used Numerical Recipes `poidev' random number generator */
//
/* On 2012 Aug 2nd poidev was replaced by gsl_ran_poisson by UB. */

#include "Python.h"
#include <stdio.h>
#include "gsl/gsl_randist.h"
#include "gsl/gsl_rng.h"

#include "numpy/arrayobject.h"

/* I (UB) found out, that the Numerical Recipies version of this module
   always called poidev with the same seed, so the result for the same input
   was always the same.
   I implemented the same behaviour with gsl_ran_poisson, to there is no 
   possibility to provide a different seed.
   I added warnings to notify the user. 
   2013 Jun 06: Not to flood the user with warnings I (UB) print out this warning
   only 10 times and suppress further warnings.
*/
#define IMGNOISESHOTWARNINGS 10

#ifdef IMGNOISESHOTWARNINGS
int imgnoise_shot_warning = 0;
#endif

static PyObject *imgnoise_shot(PyObject *self, PyObject *args)
{
	PyArrayObject	*imgin,*imgout;
	int		i,j,di,dj,diout,djout,dims[2];
	float		mean;
	const gsl_rng_type * T;
	gsl_rng * r;

#ifdef IMGNOISESHOTWARNINGS
	if( imgnoise_shot_warning < IMGNOISESHOTWARNINGS ) // Print this warning only 10 times.
	  {
	    printf("WARNING: imgnoise_shot is called with the same seed every time.\n");
	    imgnoise_shot_warning++;
	  }
	else if (imgnoise_shot_warning == IMGNOISESHOTWARNINGS)
	  {
	    printf("WARNING: furhter warnings about the same seed for imgnoise_shot suppressed.");
	    imgnoise_shot_warning++;
	  }
#endif

	// parse the input parameters:
 	if (!PyArg_ParseTuple (args, "O!O!", &PyArray_Type, &imgin, &PyArray_Type, &imgout))
		return NULL;

       /* create a generator chosen by the environment variable GSL_RNG_TYPE */
       gsl_rng_env_setup();
       T = gsl_rng_default;
       r = gsl_rng_alloc (T);

/* get input array dimensions */

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
		*(double *)(imgout->data+i*diout+j*djout) = (double)gsl_ran_poisson (r, mean);
	  }
	}

	gsl_rng_free(r);
	return Py_BuildValue("");
}

static PyObject *imgnoise_shot1d(PyObject *self, PyObject *args)
{
	PyArrayObject	*imgin,*imgout;
	int		i,di,diout,dims[1];
	float		mean;
	const gsl_rng_type * T;
	gsl_rng * r;

	printf("WARNING: imgnoise_shot1d is called with the same seed every time.\n");

	// parse the input parameters:
	if (!PyArg_ParseTuple (args, "O!O!", &PyArray_Type, &imgin, &PyArray_Type, &imgout))
		return NULL;

       /* create a generator chosen by the environment variable GSL_RNG_TYPE */
       gsl_rng_env_setup();
       T = gsl_rng_default;
       r = gsl_rng_alloc (T);
     
/* get input array dimensions */

	dims[0]=imgin->dimensions[0];
	di=imgin->strides[0];
	diout=imgout->strides[0];

/* populate output image with random Poisson deviates with mean = input image pixel value */

	for(i=0;i<dims[0];++i){
		mean = (float)(*(double *)(imgin->data + i*di));
		*(double *)(imgout->data+i*diout) = (double)gsl_ran_poisson (r, mean);
	}

	gsl_rng_free(r);

	return Py_BuildValue("");
}

#ifndef Py_RETURN_NONE
#define Py_RETURN_NONE return Py_INCREF(Py_None), Py_None
#endif

static PyObject *imgnoise_seed(PyObject *self,PyObject *args){

    int theseed=0;

    // parse the input:
    if(!PyArg_ParseTuple(args,"i",&theseed))
	return NULL;

    printf("WARNING: function imgnoise_seed in imgnoisemodule is NOT IMPLEMENTED.\n");
    printf("         The same seed will be used every time when calling imgnoise_shot!\n");

    //printf("%d\n", time(0)); // this is time in seconds; if the function 
                               // is called two times consecutively within
                               // the same second, the seed will be the same!

    Py_RETURN_NONE;
}


/* define a methods table for the module */

static PyMethodDef imgnoise_methods[] = {
					{"shot", imgnoise_shot, METH_VARARGS}, 
					{"shot1d", imgnoise_shot1d, METH_VARARGS}, 
					{"seed",imgnoise_seed,METH_VARARGS},
					{NULL, NULL} };

/* initialisation - register the methods with the Python interpreter */

void initimgnoise(void)
{
	(void) Py_InitModule("imgnoise", imgnoise_methods);
	import_array();
}










