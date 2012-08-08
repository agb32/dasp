
/* Numpy extension to create  Zernikes on a 2D Numpy array */

#include <Python.h>
#include <stdio.h>
#include <math.h>

#include "numpy/arrayobject.h"
#include "jose.h"


float fac(float n)
{
  int i;
  float f;

  f=n;

  if(n==0){
	f=1;
  }
  else{
	for(i=n-1;i>0;i--){
	  f = f*(float)i;
	}
  }

  return(f);

}

static PyObject *zernike_zern(PyObject *self,PyObject *args)
{
	PyArrayObject	*pyzern;
	int		i,j,di,dj,nx,ny,n,p1,q1,t1,s,status=0;
	float		x,y,r,theta,a,b,c,d,exp,coeff,sum;
	float		**zern;


	if (!PyArg_ParseTuple (args, "O!iii", &PyArray_Type ,&pyzern, &p1,&q1,&t1))
		return NULL;


/* get input array dimensions */

	ny=pyzern->dimensions[0];
	nx=pyzern->dimensions[1];
	di=pyzern->strides[0];
	dj=pyzern->strides[1];

	if(nx != ny){ 
		printf("Zernike array must be square !\n");
		status=1;
	        return Py_BuildValue("i",status);	
	}
	n=nx;


/* create the Zernike in a C array */

	zern=alloc2d_float(n,n);

	for(i=0;i<n;++i){ 
	   y=(float)(i-n/2)+0.5;
	   for(j=0;j<n;++j){
		x=(float)(j-n/2)+0.5;
		r=(sqrt(x*x + y*y))/((float)n/2);
		theta=atan2(y,x);

 		if(r<=1.){
		  sum=0.;
		  for(s=0;s<(p1-q1)/2+1;++s){ 
			a=pow(-1.,(float)s)*fac(p1-s);
			b=fac(s);
			c=fac((p1+q1)/2-s);
			d=fac((p1-q1)/2-s);
			coeff=a/(float)((b*c*d));
			exp=(float)(p1-2*s);
			sum=sum+coeff*pow(r,exp);
		  }

		  zern[i][j] = sqrt((float)(p1+1))*sum;

		  if(q1!=0) zern[i][j] = zern[i][j]*sqrt(2.);

		  if(t1==1){
			zern[i][j] = zern[i][j]*cos((float)(q1)*theta);
		  }
		  else{
			zern[i][j] = zern[i][j]*sin((float)(q1)*theta);
		  }
		}

	  }
	}


/* populate output Numpy array with the Zernike */

	for(i=0;i<ny;++i){
	  for(j=0;j<nx;++j){
		*(double *)(pyzern->data+i*di+j*dj) = (double)zern[i][j]; 
	  }
	}
	free2d_float(zern,n,n);

/* Return the Numpy array */

        return Py_BuildValue("i",status);
}





/* define a methods table for this module */

static PyMethodDef zernike_methods[] = 	{
					{"zern", zernike_zern, METH_VARARGS}, 
					{NULL, NULL} };


/* initialisation - register the methods with the Python interpreter */

void initzernike(void)
{
	(void) Py_InitModule("zernike", zernike_methods);
	import_array();
}

