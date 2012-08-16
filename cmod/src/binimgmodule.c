
/* Numpy extension to bin up a Numpy image/array  */
//squash an image from dimsIn pixels into array of size dimsOut pixels.
#include "Python.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "numpy/arrayobject.h"

static PyObject *binimg_binimg(PyObject *self, PyObject *args)
{
	PyArrayObject	*imgin,*imgout;
	int		i,j,i2,j2,id,jd,nbin,diin,djin,diout,djout,dimsin[2],dimsout[2];
	double		pixelin,pixelout;

 
	if (!PyArg_ParseTuple (args, "O!O!", &PyArray_Type, &imgin, &PyArray_Type, &imgout))
		return NULL;



/* get input and output array dimensions */

	dimsin[0]=imgin->dimensions[0];
	dimsin[1]=imgin->dimensions[1];
	diin=imgin->strides[0];
	djin=imgin->strides[1];

	dimsout[0]=imgout->dimensions[0];
	dimsout[1]=imgout->dimensions[1];
	diout=imgout->strides[0];
	djout=imgout->strides[1];

	nbin=dimsin[0]/dimsout[0];
	if(dimsin[0]%dimsout[0]!=0){
	  printf("cmod.binimg.binimg warning - not an integer bin factor\n");
	}

	/* printf("%d %d %d %d %d \n",dimsin[0],dimsin[1],dimsout[0],dimsout[1],nbin); */


/* do the binning */
	if(imgin->descr->type_num!=PyArray_DOUBLE && imgin->descr->type_num!=PyArray_FLOAT){
	    printf("binimg: input must be double or float\n");
	    return NULL;
	}
	if(imgout->descr->type_num!=PyArray_DOUBLE && imgout->descr->type_num!=PyArray_FLOAT){
	    printf("binimg: output must be double or float\n");
	    return NULL;
	}


	if(imgin->descr->type_num==PyArray_DOUBLE){
	    for(i=0;i<dimsout[0];++i){
		for(j=0;j<dimsout[1];++j){
		    pixelout=0.;
		    i2=i*nbin-1;
		    for(id=0;id<nbin;++id){
			++i2;
			j2=j*nbin-1;
			for(jd=0;jd<nbin;++jd){
			    ++j2;
			    pixelin=(double)(*(double *)(imgin->data + i2*diin + j2*djin));
			    pixelout+=pixelin;
			    /* printf("%d %d %d %d %d %d %f %f \n",i,j,id,jd,i2,j2,pixelin,pixelout); */
			}
		    }
		    if(imgout->descr->type_num==PyArray_DOUBLE)
			*(double *)(imgout->data+i*diout+j*djout) = pixelout; 
		    else if(imgout->descr->type_num==PyArray_FLOAT)
			*(float *)(imgout->data+i*diout+j*djout) = (float)pixelout; 
		}
	    }
	}else if(imgin->descr->type_num==PyArray_FLOAT){
	    for(i=0;i<dimsout[0];++i){
		for(j=0;j<dimsout[1];++j){
		    pixelout=0.;
		    i2=i*nbin-1;
		    for(id=0;id<nbin;++id){
			++i2;
			j2=j*nbin-1;
			for(jd=0;jd<nbin;++jd){
			    ++j2;
			    pixelin=(double)(*(float *)(imgin->data + i2*diin + j2*djin));
			    pixelout+=pixelin;
			    /* printf("%d %d %d %d %d %d %f %f \n",i,j,id,jd,i2,j2,pixelin,pixelout); */
			}
		    }
		    if(imgout->descr->type_num==PyArray_DOUBLE)
			*(double *)(imgout->data+i*diout+j*djout) = pixelout; 
		    else if(imgout->descr->type_num==PyArray_FLOAT)
			*(float *)(imgout->data+i*diout+j*djout) = (float)pixelout; 
		}
	    }
	}

	return Py_BuildValue("");


}




/* define a methods table for the module */

static PyMethodDef binimg_methods[] = 	{
					{"binimg", binimg_binimg, METH_VARARGS}, 
					{NULL, NULL} };


/* initialisation - register the methods with the Python interpreter */

void initbinimg(void)
{
	(void) Py_InitModule("binimg", binimg_methods);
	import_array();
}

