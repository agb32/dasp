
/**************************************************************************** 
*
*  Python extension to calculate reconstructed piston values using SOR algorithm 
*
****************************************************************************/ 


#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include "nr.h"
#include "nrutil.h"
#include "Python.h"


#define NUMERIC
#ifdef NUMERIC
#include "Numeric/arrayobject.h"
#else
#include "numarray/arrayobject.h"
#endif

/* Mirror Mask */


 int imirr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 8, 8, 8, 8, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 8, 8, 1, 1, 1, 1, 1, 1, 8, 8, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 2, 8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 8, 4, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 2, 8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 8, 4, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 0, 0, 0, 0,
		0, 0, 0, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 0, 0, 0,
		0, 0, 0, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7, 0, 0, 0,
		0, 0, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 0, 0,
		0, 0, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7, 0, 0,
		0, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 0,
		0, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 9, 9, 9, 9, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7, 0,
		0, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 0, 0, 0, 0, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7, 0,
		2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 0, 0, 0, 0, 0, 0, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4,
		6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7, 0, 0, 0, 0, 0, 0, 0, 0, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7,
		6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7, 0, 0, 0, 0, 0, 0, 0, 0, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7,
		6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7, 0, 0, 0, 0, 0, 0, 0, 0, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7,
		6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7, 0, 0, 0, 0, 0, 0, 0, 0, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7,
		3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 0, 0, 0, 0, 0, 0, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5,
		0, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 0, 0, 0, 0, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7, 0,
		0, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 8, 8, 8, 8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7, 0,
		0, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 0,
		0, 0, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7, 0, 0,
		0, 0, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 0, 0,
		0, 0, 0, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7, 0, 0, 0,
		0, 0, 0, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 0, 0, 0,
		0, 0, 0, 0, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 3, 9, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 9, 5, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 3, 9, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 9, 5, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 9, 9, 1, 1, 1, 1, 1, 1, 9, 9, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 9, 9, 9, 9, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };


#define  ISEGNUM 32
#define  ISEGSZ 8


float conv=0.01;


static PyObject *sor32_fit(self,args)
	PyObject *self, *args;

{
  PyArrayObject	*py_angx,*py_angy,*py_pist;
  int i,j,ix,iy;
  int index,ns;
  int iter;
  int ndx,dix,djx,dimsx[2],ndy,diy,djy,dimsy[2],nd,di,dj,dims[2];
  float w;
  float err, pl;
  float avpist;
  float diff;
  float h;
  float pist[ISEGNUM][ISEGNUM];
  float angx[ISEGNUM][ISEGNUM],angy[ISEGNUM][ISEGNUM];
  float p[ISEGNUM][ISEGNUM];
  float sx[ISEGNUM][ISEGNUM],sy[ISEGNUM][ISEGNUM];


  if (!PyArg_ParseTuple (args, "O!O!O!", &PyArray_Type, &py_angx, &PyArray_Type, &py_angy, &PyArray_Type, &py_pist))
	return NULL;



/* get input Python array dimensions */

	ndx=py_angx->nd;
	dimsx[0]=py_angx->dimensions[0];
	dimsx[1]=py_angx->dimensions[1];
	dix=py_angx->strides[0];
	djx=py_angx->strides[1];

	ndy=py_angy->nd;
	dimsy[0]=py_angy->dimensions[0];
	dimsy[1]=py_angy->dimensions[1];
	diy=py_angy->strides[0];
	djy=py_angy->strides[1];

	nd=py_pist->nd;
	dims[0]=py_pist->dimensions[0];
	dims[1]=py_pist->dimensions[1];
	di=py_pist->strides[0];
	dj=py_pist->strides[1];

	if(dimsx[0]!=ISEGNUM || dimsx[1]!=ISEGNUM || dimsy[0]!=ISEGNUM || dimsy[1]!=ISEGNUM){
		printf("Tilt array wrong size.");
	}
	if(dims[0]!=ISEGNUM || dims[1]!=ISEGNUM ){
		printf("Piston array wrong size.");
	}


/* convert input arrays to C array */

	for(i=0;i<ISEGNUM;++i){
		for(j=0;j<ISEGNUM;++j){
		    if(py_angx->descr->type_num==PyArray_DOUBLE)
			angx[j][i] = (float)(*(double *)(py_angx->data + i*dix + j*djx));
		    else
			angx[j][i] = (float)(*(float*)(py_angx->data + i*dix + j*djx));
		    if(py_angy->descr->type_num==PyArray_DOUBLE)
			angy[j][i] = (float)(*(double *)(py_angy->data + i*diy + j*djy));
		    else
			angy[j][i] = (float)(*(float *)(py_angy->data + i*diy + j*djy));
		    if(py_pist->descr->type_num==PyArray_DOUBLE)
			pist[j][i] = (float)(*(double *)(py_pist->data + i*di + j*dj));
		    else
			pist[j][i] = (float)(*(float *)(py_pist->data + i*di + j*dj));

		}
	}


  avpist=0;
  h=(ISEGSZ)*(2.0*3.14157);

  w=2.0/(1.0+sin(3.14157/(ISEGNUM+1.0)));

  ns=0;
    for (ix=0;ix<ISEGNUM;ix++) {

  for (iy=0;iy<ISEGNUM;iy++) {
      index=iy*ISEGNUM+ix;
      if(imirr[index]!=0) {
	p[ix][iy]=pist[ix][iy];
	sx[ix][iy]=0.5*angx[ix][iy];
	sy[ix][iy]=0.5*angy[ix][iy];
	ns++;
      } else {
	p[ix][iy]=0.0;
	sx[ix][iy]=0.0;
	sy[ix][iy]=0.0;
      }
    }
  } 
  
  iter=0;
  err=1.0;
  while(err>conv) {
	/*printf("... %d %f\n",iter,err);*/
    ns=0;   
    err=0;
    for (iy=0;iy<ISEGNUM;iy++) {
      for (ix=0;ix<ISEGNUM;ix++) {
	
	index=iy*ISEGNUM+ix;
	
	/*  corners  */
	
	if(imirr[index]==2) {
	  pl=p[ix][iy];
	  p[ix][iy]+=
	    w*((p[ix+1][iy]+p[ix][iy+1]
		+(sy[ix][iy+1]+sy[ix][iy]+sx[ix+1][iy]+sx[ix][iy])*h)/2
	       -p[ix][iy]);
	  diff=p[ix][iy]-pl;
	}
	
	if(imirr[index]==3) {
	  pl=p[ix][iy];
	  p[ix][iy]+=
	    w*((p[ix][iy-1]+p[ix+1][iy]
		+(-sy[ix][iy]-sy[ix][iy-1]+sx[ix+1][iy]+sx[ix][iy])*h)/2
	       -p[ix][iy]);
	  diff=p[ix][iy]-pl;
	}
	
	if(imirr[index]==4) {
	  pl=p[ix][iy];
	  p[ix][iy]+=
	    w*((p[ix][iy+1]+p[ix-1][iy]
		+(sy[ix][iy+1]+sy[ix][iy]-sx[ix][iy]-sx[ix-1][iy])*h)/2
	       -p[ix][iy]);
	  diff=p[ix][iy]-pl;
	}
	
	if(imirr[index]==5) {
	  pl=p[ix][iy];
	  p[ix][iy]+=
	    w*((p[ix][iy-1]+p[ix-1][iy]
		+(-sy[ix][iy]-sy[ix][iy-1]-sx[ix][iy]-sx[ix-1][iy])*h)/2
	       -p[ix][iy]);
	  diff=p[ix][iy]-pl;
	}
	
	
	/*     sides   */
	
	if(imirr[index]==6) {
	  pl=p[ix][iy];
	  p[ix][iy]+=
	    w*((p[ix][iy+1]+p[ix][iy-1]+p[ix+1][iy]
		+(sy[ix][iy+1]-sy[ix][iy-1]+sx[ix+1][iy]+sx[ix][iy])*h)/3
	       -p[ix][iy]);
	  diff=p[ix][iy]-pl;
	}
	
	if(imirr[index]==7) {
	  pl=p[ix][iy];
	  p[ix][iy]+=
	    w*((p[ix][iy+1]+p[ix][iy-1]+p[ix-1][iy]
		+(sy[ix][iy+1]-sy[ix][iy-1]-sx[ix][iy]-sx[ix-1][iy])*h)/3
	       -p[ix][iy]);
	  diff=p[ix][iy]-pl;
	}
	
	if(imirr[index]==8) {
	  pl=p[ix][iy];
	  p[ix][iy]+=
	    w*((p[ix][iy+1]+p[ix+1][iy]+p[ix-1][iy]
		+(sy[ix][iy+1]+sy[ix][iy]+sx[ix+1][iy]-sx[ix-1][iy])*h)/3
	       -p[ix][iy]);
	  diff=p[ix][iy]-pl;
	}
	
	if(imirr[index]==9) {
	  pl=p[ix][iy];
	  p[ix][iy]+=
	    w*((p[ix][iy-1]+p[ix+1][iy]+p[ix-1][iy]
		+(-sy[ix][iy]-sy[ix][iy-1]+sx[ix+1][iy]-sx[ix-1][iy])*h)/3
	       -p[ix][iy]);
	  diff=p[ix][iy]-pl;
	}

	
	/*     main body   */

	if(imirr[index]==1) {  
	  pl=p[ix][iy];
	  p[ix][iy]+=
	    w*((p[ix][iy+1]+p[ix][iy-1]+p[ix+1][iy]+p[ix-1][iy]
		+(sy[ix][iy+1]-sy[ix][iy-1]+sx[ix+1][iy]-sx[ix-1][iy])*h)/4
	       -p[ix][iy]);
	  diff=p[ix][iy]-pl;
	}
	
	if(imirr[index]!=0) {  
	  err+=fabs(diff); 
          ns++;

	}  
      }
    }

    err=err/(float)ns;
    
   iter++;
  }
 
  ns=0;
  for (iy=0;iy<ISEGNUM;iy++) {
    for (ix=0;ix<ISEGNUM;ix++) {
      index=iy*ISEGNUM+ix;
      pist[ix][iy]=p[ix][iy];
      /* avpist+=pist[ix][iy]*pist[ix][iy] */;     
	avpist+=pist[ix][iy];
    }
  } 
  avpist=avpist/736.0;
  /* printf("%f\n", avpist); */




 /* populate output Python array with piston values */

	for(i=0;i<ISEGNUM;++i){
	 for(j=0;j<ISEGNUM;++j){
	   /* *(double *)(py_pist->data+i*di+j*dj) = (double)((pist[j][i]-avpist)/(double)ISEGSZ); */
	   *(double *)(py_pist->data+i*di+j*dj) = (double)((pist[j][i]-avpist)/8.);
	 /* printf("%d %d %d %f %f %f\n",i,j,imirr[index],angx[i][j],angy[i][j],pist[i][j]); */

	  }
	}


 return Py_BuildValue("");



}


/* =============================================================================== */



/* define a methods table for this module */

static PyMethodDef sor32_methods[] = 	{
					{"fit", sor32_fit, METH_VARARGS}, 
					{NULL, NULL} };


/* initialisation - register the methods with the Python interpreter */

void initsor32()
{
	(void) Py_InitModule("sor32", sor32_methods);
	import_array();
}
