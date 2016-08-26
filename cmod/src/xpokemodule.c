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

/* Do the pokin' for a Xintetics style DM  */

#include "Python.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "numpy/arrayobject.h"
#include "jose.h"



int              nx,nact,ndm;
int             *iy,*jx;
float          ***pokefn;

/*============================================================== */

/* Setup arrays for poke functions and actuator centers */
/* Inputs are actuator centers and Gaussian sigmas  */

static PyObject *xpoke_setup(PyObject *self, PyObject *args)
{
	PyArrayObject	*xain,*yain,*krin,*ktin;
	int		i,j,iact,di;
	float          x,y,xc,yc,rp,th,phi,rr,rt; // ra not used - commented out by UB, 2012Aug08
	float          *xa,*ya,*kr,*kt;


 
	if (!PyArg_ParseTuple (args, "iiffO!O!O!O!", &nx,&ndm,&xc,&yc,&PyArray_Type, &xain, &PyArray_Type, &yain, &PyArray_Type, &krin, &PyArray_Type, &ktin))
		return NULL;



/* get input array dimensions (number of actuators) */

	nact=xain->dimensions[0];
	di=xain->strides[0];

	printf("%d %d %d %f %f  \n",nact,nx,ndm,xc,yc); 

	iy=calloc(nact,sizeof(int));
	jx=calloc(nact,sizeof(int));
	xa=calloc(nact,sizeof(float));
	ya=calloc(nact,sizeof(float));
	kr=calloc(nact,sizeof(float));
	kt=calloc(nact,sizeof(float));



	pokefn=calloc(nact,sizeof(float**));
	for(iact=0;iact<nact;iact++){
	  pokefn[iact]=calloc(nx,sizeof(float*));
	  for(i=0;i<nx;i++){
		pokefn[iact][i]=calloc(nx,sizeof(float));
	  }
	}




/* Copy to C arrays */

	for(iact=0;iact<nact;++iact){
	  xa[iact]=(float)(*(double *)(xain->data + iact*di));
	  ya[iact]=(float)(*(double *)(yain->data + iact*di));
	  iy[iact]=(int)ya[iact];
	  jx[iact]=(int)xa[iact];
	  kr[iact]=(float)(*(double *)(krin->data + iact*di));
	  kt[iact]=(float)(*(double *)(ktin->data + iact*di));
	  /* printf("%d %d %d %f %f  %f %f \n",iact,iy[iact],jx[iact],xa[iact],ya[iact],kr[iact],kt[iact]); */
	}


/* Create poke functions in 3d array */

	for(iact=0;iact<nact;++iact){
	  for(i=0;i<nx;++i){
		for(j=0;j<nx;++j){
		  x=(float)(j-nx/2);
		  y=(float)(i-nx/2);
		  //ra = sqrt(  pow((xa[iact]-xc),2.)  +  pow((ya[iact]-yc),2.)  ); // ra not used - commented out by UB, 2012Aug08
		  rp = sqrt(  pow(x,2.)  +  pow(y,2.)  );
		  th   = atan((ya[iact]-yc)/(xa[iact]-xc));
		  phi = atan2(y,x);
		  rr = rp * cos(phi-th); 
		  rt = rp * sin(phi-th);

		  pokefn[iact][i][j] = exp( (-pow(rr,2.)) / (pow(kr[iact],2.)) ) * exp( (-pow(rt,2.)) / (pow(kt[iact],2.)) );

		}
	  }
	}


	return Py_BuildValue("");


}




/*============================================================== */

/* Do the pokin' and return the DM figure, for input poke amplitude vector */

static PyObject *xpoke_poke(PyObject *self, PyObject *args)
{
	PyArrayObject	*coeffin,*dmfig;
	int		i,j,ii,jj,iact,di,diout,djout; //,id=0;  id not used, commented out by UB, 2012 Aug 08
	float         coeff;



 
	if (!PyArg_ParseTuple (args, "O!O!", &PyArray_Type, &coeffin, &PyArray_Type, &dmfig))
		return NULL;

	di=coeffin->strides[0];
	diout=dmfig->strides[1];
	djout=dmfig->strides[0];

	for(iact=0;iact<nact;++iact){
	  coeff= *(double *)(coeffin->data + iact*di);

	  if(coeff>0.){
	  if( (iy[iact]-nx/2)>0 &&  (iy[iact]+nx/2)<ndm   &&   ( jx[iact]-nx/2)>0 &&  (jx[iact]+nx/2)<ndm   ){

		for(i=0;i<nx;i++){
		  for(j=0;j<nx;j++){
			ii=iy[iact]+i-nx/2;
			jj=jx[iact]+j-nx/2;
			*(double *)(dmfig->data+ii*diout+jj*djout) += coeff*pokefn[iact][i][j]; 
		  }
		}

	  }
	  else{

		for(i=0;i<nx;i++){
		  for(j=0;j<nx;j++){
			ii=iy[iact]+i-nx/2;
			jj=jx[iact]+j-nx/2;
			if( (ii>0) && (ii<ndm) && (jj>0) && (jj<ndm) ){ 
			  *(double *)(dmfig->data+ii*diout+jj*djout) += coeff*pokefn[iact][i][j]; 
			}
		  }
		}
		
	  }
	  }

	}
	//id=1; id not used, commented out by UB, 2012 Aug 08

	return Py_BuildValue("");


}

/*============================================================== */

/* define a methods table for the module */

static PyMethodDef xpoke_methods[] = 	{
					{"setup", xpoke_setup, METH_VARARGS}, 
					{"poke", xpoke_poke, METH_VARARGS},
					{NULL, NULL} };


/* initialisation - register the methods with the Python interpreter */

void initxpoke(void)
{
	(void) Py_InitModule("xpoke", xpoke_methods);
	import_array();
}

