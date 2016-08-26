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

/* Numpy extension to make an image from a phase array */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef USE_FFTW_FLOAT
  #include <sfftw.h>
  #define MODULENAME "mkimgfloat"
  #define INITFUNCTION initmkimgfloat()
#else
  #include <fftw.h>
  #define MODULENAME "mkimg"
  #define INITFUNCTION initmkimg()
#endif
#include "Python.h"
//which if sfftw or fftw.h is used will determine whether fftw_real is float or double.

#include <numpy/arrayobject.h>
/*
#define NUMERIC
#ifdef NUMERIC
#include "Numeric/arrayobject.h"
#else
#include "numarray/arrayobject.h"
#endif
*/


static PyObject *mkimg_setup(self,args)
	PyObject *self, *args;

{
	PyArrayObject	 *imgout;
	int		diout,djout,size,nxout,nyout;
	fftwnd_plan p;

 
	if (!PyArg_ParseTuple (args, "O!",  &PyArray_Type, &imgout))
		return NULL;


/* get FFT array dimensions */

	nxout=imgout->dimensions[0];
	nyout=imgout->dimensions[1];


	p=fftw2d_create_plan(nxout,nyout,FFTW_FORWARD,FFTW_MEASURE | FFTW_IN_PLACE);
	//return the pointer...
	return Py_BuildValue("l",(long)p);


}



static PyObject *mkimg_mkimg(self,args)
	PyObject *self, *args;

{
	PyArrayObject	*phsin,*imgout,*pypupil;
	int		i,j,ii,nd,di,dj,diout,djout,size,nx,ny,nxout,nyout,ndim,isign;
	fftw_real img,phs,pupil;
	fftw_complex  *amp;

	fftwnd_plan p;
	//printf("cmod.mkimg depreciated\n");
 
	if (!PyArg_ParseTuple (args, "lO!O!O!",&p, &PyArray_Type, &phsin, &PyArray_Type, &imgout, &PyArray_Type, &pypupil))
		return NULL;
	//here, p is the pointer returned from mkimg_setup.

	//check types of the arrays.
#ifdef USE_FFTW_FLOAT
	if(phsin->descr->type_num!=NPY_FLOAT || imgout->descr->type_num!=NPY_FLOAT || pypupil->descr->type_num!=NPY_FLOAT){
	    printf("Inputs to mkimg must be of type float (current: %d %d %d)\n",phsin->descr->type_num,imgout->descr->type_num,pypupil->descr->type_num);
	    return NULL;
	}
#else
	if(phsin->descr->type_num!=NPY_DOUBLE || imgout->descr->type_num!=NPY_DOUBLE || pypupil->descr->type_num!=NPY_DOUBLE){
	    printf("Inputs to mkimg must be of type double (current: %d %d %d)\n",phsin->descr->type_num,imgout->descr->type_num,pypupil->descr->type_num);
	    return NULL;
	}
#endif

/* get input array dimensions, assume pupil is same size as phase array !!!! */

	nd=phsin->nd;
	nx=phsin->dimensions[0];
	ny=phsin->dimensions[1]; 
	di=phsin->strides[0];
	dj=phsin->strides[1];

	diout=imgout->strides[0];
	djout=imgout->strides[1];
	nxout=imgout->dimensions[0];
	nyout=imgout->dimensions[1];


	/* copy input phase to a C complex amp array  */
	
	amp = (fftw_complex *)calloc(nxout*nyout,sizeof(fftw_complex)); 
	
	for(i=0;i<ny;++i){
	    for(j=0;j<nx;++j){
		pupil = *(fftw_real *)(pypupil->data + i*di + j*dj);
		phs = *(fftw_real *)(phsin->data + i*di + j*dj);
		amp[i*nyout+j].re    = cos(phs)*pupil;
		amp[i*nyout+j].im    = sin(phs)*pupil;
	   }
	 }
	

	fftwnd_one(p,amp,NULL);


	/* Populate output Numpy array */


	ii=0;
	for(i=0;i<nyout;++i){
	  for(j=0;j<nxout;++j){

		ii=i*nyout+j;
		
		img = (amp[ii].re * amp[ii].re) + (amp[ii].im * amp[ii].im);

		*(fftw_real *)(imgout->data+i*diout+j*djout) = (fftw_real)img;
 		/* printf("%d %d %f  %f \n", i,j,img,*(fftw_real *)(imgout->data+i*diout+j*djout)); */


	  }
	}

	free(amp);

	return Py_BuildValue("");


}




/* define a methods table for the module */

static PyMethodDef mkimg_methods[] = 	{
					{"setup", mkimg_setup, METH_VARARGS},
					{"mkimg", mkimg_mkimg, METH_VARARGS}, 
					{NULL, NULL} };


/* initialisation - register the methods with the Python interpreter */
void INITFUNCTION
{
	(void) Py_InitModule(MODULENAME, mkimg_methods);
	import_array();
}


