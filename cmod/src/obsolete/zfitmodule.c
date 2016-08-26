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
// 
// This (zfitmodule.c) is apparently not used anywhere in aosim and will most
// probably not be
// used in the future. I (UB) moved it to the "obsolete" folder on 2012 Aug 02.
// Could probably be removed completely, but to stay on the safe side we keep
// it for a bit longer.
// 

/* find Zernike coefficients for Shack Hartmann centroid data */
/* uses NR SVD routine */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "nr.h"
#include "nrutil.h"
#include "Python.h"

#include "numpy/arrayobject.h"
#include "jose.h"

int		ndata,nz;
float		*w,*data,*coeff;
float		**u,**v;

#define	TOL	1.e-6



/* ================ Set up SDV matrices =========================================== */


static PyObject *zfit_setup(self,args)
	PyObject *self, *args;

{
	PyArrayObject	*py_zgrad;
	int		i,j,di,dj,idata,m,nd,status=0;
	float		thresh,wmax;
	float		**zgrad,**temp;

	if (!PyArg_ParseTuple (args, "O!", &PyArray_Type, &py_zgrad))
		return NULL;

/* get input data array dimensions (2*nsubaps*nzern)  */

	nd=py_zgrad->nd;
	ndata=py_zgrad->dimensions[0];
	nz=   py_zgrad->dimensions[1];
	/* printf("%d %d\n",ndata,nz); */
	di=py_zgrad->strides[0];
	dj=py_zgrad->strides[1];

/* allocate array space for NR SVD matrices */
	
	zgrad=matrix(1,ndata,1,nz);
	u=matrix(1,ndata,1,nz);
	v=matrix(1,nz,1,nz);
	w=vector(1,nz);
	data=vector(1,ndata);
	coeff=vector(1,nz);

	temp=alloc2d_float(ndata,nz);


/* convert input array to C array */

	for(i=0;i<ndata;++i){
		for(j=0;j<nz;++j){
		 zgrad[i+1][j+1] = (float)(*(double *)(py_zgrad->data + i*di + j*dj));
		}
	}


/* populate SVD matrix.  zgrad[][] starts at [1][1] */


/* do SVD, edit singular values	*/
//Decomposes zgrad into UWV^T, and puts U into zgrad, gives W as vector and V not V^T.
 	svdcmp(zgrad,ndata,nz,w,v); 

	wmax=0.;
	for(m=1;m<=nz;++m){
		/* printf("%d %f\n",m,w[m]); */
		 if(w[m]>wmax) wmax=w[m];
	}
	thresh=TOL*wmax;
	for(m=1;m<nz;++m){
		 if(w[m]<thresh) w[m]=0.;
	}
	for(idata=1;idata<=ndata;++idata){
	 for(m=1;m<=nz;++m){
		 u[idata][m] = zgrad[idata][m];
	 }
	}
        return Py_BuildValue("i",status);
}

 /* =============== Do SDV fit for this data. Return Zern coeffs =============================== */

static PyObject *zfit_fit(self,args)
	PyObject *self, *args;

{
	PyArrayObject	*py_data,*py_coeff;
	int	 	i,nd,di,diout,dims[1];
	float		*val;

	if (!PyArg_ParseTuple (args, "O!O!", &PyArray_Type, &py_data, &PyArray_Type, &py_coeff))
		return NULL;

/* get input data array dimensions */

	nd=py_data->nd;
	dims[0]=py_data->dimensions[0];		/* Check this is same as NZ !!! */
	if(dims[0]!=ndata) printf("Wrong number of data... %d %d\n",dims[0],ndata);
	di=py_data->strides[0];

	diout=py_coeff->strides[0];

/* convert input data to C array */

	for(i=0;i<ndata;++i){
		data[i+1] = (float)(*(double *)(py_data->data + i*di));
	}

/* do the SVD fit */

	svbksb(u,w,v,ndata,nz,data,coeff);

/* create new python array for output coeffs */

	nd=1;
	dims[0]=nz;
	/* py_coeff=(PyArrayObject *)PyArray_FromDims(nd,dims,PyArray_DOUBLE); */

 /* populate output array with Zernike coeff values */

	for(i=0;i<nz;++i){
		*(double *)(py_coeff->data+i*diout) = (double)coeff[i+1]; 
	}

/* Return the Numpy array */

	/* return PyArray_Return(py_coeff); */

	return Py_BuildValue("");

}

/* =============================================================================== */



/* define a methods table for this module */

static PyMethodDef zfit_methods[] = 	{
					{"setup", zfit_setup, METH_VARARGS}, 
					{"fit", zfit_fit, METH_VARARGS}, 
					{NULL, NULL} };


/* initialisation - register the methods with the Python interpreter */

void initzfit()
{
	(void) Py_InitModule("zfit", zfit_methods);
	import_array();
}
