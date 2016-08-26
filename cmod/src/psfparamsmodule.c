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

/* Calculate PSF azimuthal average etc. */

#include "Python.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "numpy/arrayobject.h"
#include "jose.h"

/* ==================  azav =================================================*/ 
/* find PSF azimuthal average, around peak pixel */


static PyObject *psfparams_azav(PyObject *self, PyObject *args)
{
	PyArrayObject	*pypsf,*pyazav,*pyencirc;
	int		i,j,ir,nx,ny,di,dj,imax=0,jmax=0;
	float		x,y,r,max,total,sum;
	float		**psf,*azav,*encirc,*npts;


	if (!PyArg_ParseTuple (args, "O!O!O!", &PyArray_Type ,&pypsf, &PyArray_Type ,&pyazav, &PyArray_Type ,&pyencirc))
		return NULL;


/* get input psf array dimensions */


	ny=pypsf->dimensions[0];
	nx=pypsf->dimensions[1];
	di=pypsf->strides[0];
	dj=pypsf->strides[1];
	//agb assume input type is a double (64 bit)?

/* copy to a C array */


	//printf("Performing copy of the input PSF\n");
	psf=alloc2d_float(ny,nx); 


	for(i=0;i<ny;++i){
	    for(j=0;j<nx;++j){
		psf[i][j] = (float)(*(double *)(pypsf->data + i*di + j*dj));
	   }
	 }


	//printf("Performing memory allocation of the arrays\n");

	azav=calloc(nx,sizeof(float));
	encirc=calloc(nx,sizeof(float));
	npts=calloc(nx,sizeof(float));




/* find peak pixel */


	max=0.;
	for(i=0;i<ny;++i){
	 for(j=0;j<nx;++j){
	  if(psf[i][j]>max){ 
		max=psf[i][j];
		imax=i;
		jmax=j;
	  }
	 }
	}




/* calculate azimuthal average  */


	//printf("Performing azimutal average\n");

	total=0;
	for(i=0;i<ny;++i){
	 y=(float)(i-imax);
	 for(j=0;j<nx;++j){
	  x=(float)(j-jmax);
	  r=sqrt(x*x + y*y);
	  ir=(int)r;
	  total+=psf[i][j];
	  if(ir<nx){//agb without this, we can overwrite end of array...
	      azav[ir] += psf[i][j];
	      npts[ir] += 1.;
	  }
	 }
	}

	sum=0;
	for(ir=0;ir<nx;++ir){
	 sum+=azav[ir];
	 encirc[ir]=sum;
	 if(npts[ir]>0)  azav[ir] /= npts[ir];
	}


/* encircled energy */

	//printf("Performing encircled energy\n");

	for(ir=0;ir<nx;++ir){
		encirc[ir] /= total;
	}




/* populate output Numpy arrays with the azav values */


	//printf("Copy to output arrays\n");

	di=pyazav->strides[0];
	for(i=0;i<nx;++i){
		*(double *)(pyazav->data+i*di)   = (double)azav[i]; 
		*(double *)(pyencirc->data+i*di) = (double)encirc[i]; 
	}

	//printf("Freeing azimutal average array\n");
	free(azav);
	//printf("Freeing Encircled energy array\n");
	free(encirc);
	//printf("Freeing npts array\n");
	free(npts);

	//printf("Freeing psf copy\n");

        for(i=0;i<ny;++i){
                free(psf[i]);
        }
	free(psf);


	return Py_BuildValue("");

	/* return Py_None; 			Don't Use this !!!! */


}



/* ===============================================================================*/ 


/* define a methods table for this module */

static PyMethodDef psfparams_methods[] = 	{
					{"azav", psfparams_azav, METH_VARARGS}, 
					{NULL, NULL} };


/* initialisation - register the methods with the Python interpreter */

void initpsfparams(void)
{
	(void) Py_InitModule("psfparams", psfparams_methods);
	import_array();
}

