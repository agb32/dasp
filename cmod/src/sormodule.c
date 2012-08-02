
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
//#include "nr.h"
//#include "nrutil.h"
#include "Python.h"


#include "numpy/arrayobject.h"
//#define  ISEGNUM 32
#define  ISEGSZ 8


/*This can actually be done in python...
static PyObject *sor_createMask(PyObject *self,PyObject *args){
    PyObject *pupfno;
    PyArrayObject *pupfn;
    PyArrayObject *outarr;
    int *imirr;
    int size;
    if (!PyArg_ParseTuple (args, "O!", &PyArray_Type, &pupfno))
	return NULL;
    pupfn=(PyArrayObject*)pupfno;
    xx=pupfn->dimensions[1];
    yy=pupfn->dimensions[0];
    size=xx*yy;
    imirr=(int*)malloc(sizeof(int)*size);
    
    outarr=(PyArrayObject *)PyArray_FromDims(2,dims,PyArray_FLOAT);
*/
typedef struct{
  int dimirr;
  int isegnum;
  int isegnum2;
  PyArrayObject *imirr;
  double avpistDivisor;
  double conv;
  int maxiters;
  PyObject *tmparrobj;
  PyArrayObject *tmparr;
  float *p;
  float *ss;
  float w;
  float h;
} SorStruct;

static PyObject *sor_init(PyObject *self,PyObject *args){
  //initialise the SOR object.  Have found that passing too many
  //arguments to sor_fit causes a performance loss, and so it is
  //better to initialise a struct here.
  SorStruct *sorstr;
  double conv=0.01;
  PyObject *tmparrobj=NULL;
  PyArrayObject *imirr;
  double avpistDivisor;
  int maxiters=1000;

  if(!PyArg_ParseTuple(args,"O!d|diO",&PyArray_Type,&imirr,&avpistDivisor,&conv,&maxiters,&tmparrobj)){
    printf("sor_init Usage: sormask, avpistdivisor, convergence value (optional, default 0.01), max iters (optional, default 1000), tmparrobj (optional, will be malloced if not supplied)\n");
    return NULL;
  }
  sorstr=malloc(sizeof(SorStruct));
  sorstr->dimirr=imirr->strides[0];
  sorstr->isegnum=sqrt(imirr->dimensions[0]);
  sorstr->isegnum2=sorstr->isegnum*sorstr->isegnum;
  sorstr->imirr=imirr;
  sorstr->avpistDivisor=avpistDivisor;
  sorstr->tmparrobj=tmparrobj;
  sorstr->tmparr=NULL;
  sorstr->conv=conv;
  sorstr->maxiters=maxiters;
  

  if(imirr->dimensions[0]!=sorstr->isegnum2){
    printf("Mirror mask is wrong shape (should be 1D flat array size %d)\n",sorstr->isegnum2);
    free(sorstr);
    return NULL;
  }
  if(imirr->descr->kind!='i' || imirr->descr->elsize!=sizeof(int)){
    printf("Mirror mask should be int\n");
    free(sorstr);
    return NULL;
  }
  if(avpistDivisor>sorstr->isegnum2){
    printf("avpistDivisor should be double, current value %g should be less than %d\n",avpistDivisor,sorstr->isegnum2);
  }

  if(tmparrobj==NULL || !PyArray_Check(tmparrobj)){
    sorstr->p=(float*)malloc(sizeof(float*)*sorstr->isegnum2);
    sorstr->ss=(float*)malloc(sizeof(float*)*sorstr->isegnum2);
  }else{
    sorstr->tmparr=(PyArrayObject*)tmparrobj;
    if(sorstr->tmparr->nd==3 && sorstr->tmparr->descr->type_num==PyArray_FLOAT && sorstr->tmparr->dimensions[0]==2 && sorstr->tmparr->dimensions[1]==sorstr->isegnum && sorstr->tmparr->dimensions[2]==sorstr->isegnum){
      sorstr->p=(float*)(sorstr->tmparr->data);
      sorstr->ss=&(((float*)(sorstr->tmparr->data))[sorstr->isegnum2]);
    }else{
      printf("temporary array should have dimensions (2,%d,%d), type Float32\n",sorstr->isegnum,sorstr->isegnum);
      free(sorstr);
      return NULL;
    }
  }

  sorstr->h=(ISEGSZ)*(2.0*3.14157)/2;
  sorstr->w=2.0/(1.0+sin(3.14157/(sorstr->isegnum+1.0)));

  return Py_BuildValue("l",(long)sorstr);//return the address...
}
static PyObject *sor_setconv(PyObject *self,PyObject *args){
  SorStruct *sorstr;
  double conv;
  if(!PyArg_ParseTuple(args,"ld",&sorstr,&conv)){
    printf("setconv: usage sorstruct, convergence value\n");
    return NULL;
  }
  sorstr->conv=conv;
  return Py_BuildValue("i",0);
}

static PyObject *sor_free(PyObject *self,PyObject *args){
  SorStruct *sorstr;
  if(!PyArg_ParseTuple(args,"l",&sorstr)){
    printf("sor_free: usage - sorStruct (from init)\n");
    return NULL;
  }
  

  if(sorstr->tmparr==NULL){
    free(sorstr->p);
    free(sorstr->ss);
  }
  free(sorstr);
  return Py_BuildValue("i",0);
}


static PyObject *sor_fit(PyObject *self,PyObject *args){
    PyArrayObject *py_angx,*py_angy,*py_pist;
    int i,j,dimirr;
    //int ix,iy;
    int index,ns;
    int iter;
    int ndx,dix,djx,dimsx[2],ndy,diy,djy,dimsy[2],nd,di,dj,dims[2];
    float w;
    float err;
    float avpist;
    float diff;
    float h;
    float *p;
    float *ss;
    char *sormask;
    double conv=0.01;
    double avpistDivisor;
    int isegnum,isegnum2;
    int maxiters=1000;
    SorStruct *sorstr;
    if (!PyArg_ParseTuple (args, "O!O!O!l", &PyArray_Type, &py_angx, &PyArray_Type, &py_angy, &PyArray_Type, &py_pist,&sorstr)){
      printf("Need centx, centy, piston (output), sorStruct (from init)\n");//, imirr (mask), avpistDivisor, conv(optional),maxiters(optional),temporary array (optional)\n");
      return NULL;
    }
    //printf("%g\n",conv);
    
/* get input Python array dimensions */
    //dimirr=imirr->strides[0];
    dimirr=sorstr->dimirr;
    isegnum=sorstr->isegnum;    //isegnum=dimsx[0];
    isegnum2=sorstr->isegnum2;//    isegnum2=isegnum*isegnum;
    p=sorstr->p;
    ss=sorstr->ss;
    h=sorstr->h;
    w=sorstr->w;
    sormask=sorstr->imirr->data;
    avpistDivisor=sorstr->avpistDivisor;
    conv=sorstr->conv;
    maxiters=sorstr->maxiters;
    
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
    
    if(ndx!=2 || ndy!=2 || nd!=2){
      printf("Arrays are wrong size\n");
      return NULL;
    }
    if(dimsx[0]!=isegnum || dimsx[1]!=isegnum || dimsy[0]!=isegnum || dimsy[1]!=isegnum){
	printf("Tilt array wrong size.\n");
	return NULL;
    }
    if(dims[0]!=isegnum || dims[1]!=isegnum ){
	printf("Piston array wrong size.\n");
	return NULL;
    }
    /* convert input arrays to C array, appropriately scaled. */
    if(py_angx->descr->type_num==PyArray_DOUBLE && py_angy->descr->type_num==PyArray_DOUBLE && py_pist->descr->type_num==PyArray_DOUBLE){
      for(i=0;i<isegnum;++i){
	for(j=0;j<isegnum;++j){
	  index=i*isegnum+j;
	  switch(*(int*)(sormask+index*dimirr)){
	    case 2://(sy[index+isegnum]+sy[index]+sx[index+1]+sx[index])*h
	      ss[index]=(float)(h/2*(+(*(double *)(py_angy->data + (i+1)*diy + j*djy))  
				     +(*(double *)(py_angy->data + (i+0)*diy + j*djy)) 
				     +(*(double *)(py_angx->data + i*diy + (j+1)*djy)) 
				     +(*(double *)(py_angx->data + i*diy + (j+0)*djy)) ));
	      p[index] = (float)(*(double *)(py_pist->data + i*di + j*dj));
	      break;
	    case 3://(-sy[index]-sy[index-isegnum]+sx[index+1]+sx[index])*h
	      ss[index]=(float)(h/2*(-(*(double *)(py_angy->data + (i+0)*diy + j*djy))  
				     -(*(double *)(py_angy->data + (i-1)*diy + j*djy)) 
				     +(*(double *)(py_angx->data + i*diy + (j+1)*djy)) 
				     +(*(double *)(py_angx->data + i*diy + (j+0)*djy)) ));
	      p[index] = (float)(*(double *)(py_pist->data + i*di + j*dj));
	      break;
	    case 4://(sy[index+isegnum]+sy[index]-sx[index]-sx[index-1])*h
	      ss[index]=(float)(h/2*(+(*(double *)(py_angy->data + (i+1)*diy + j*djy))  
				     +(*(double *)(py_angy->data + (i+0)*diy + j*djy)) 
				     -(*(double *)(py_angx->data + i*diy + (j+0)*djy)) 
				     -(*(double *)(py_angx->data + i*diy + (j-1)*djy)) ));
	      p[index] = (float)(*(double *)(py_pist->data + i*di + j*dj));
	      break;
	    case 5://(-sy[index]-sy[index-isegnum]-sx[index]-sx[index-1])*h
	      ss[index]=(float)(h/2*(-(*(double *)(py_angy->data + (i+0)*diy + j*djy))  
				     -(*(double *)(py_angy->data + (i-1)*diy + j*djy)) 
				     -(*(double *)(py_angx->data + i*diy + (j+0)*djy)) 
				     -(*(double *)(py_angx->data + i*diy + (j-1)*djy)) ));
	      p[index] = (float)(*(double *)(py_pist->data + i*di + j*dj));
	      break;
	    case 6://(sy[index+isegnum]-sy[index-isegnum]+sx[index+1]+sx[index])*h
	      ss[index]=(float)(h/3*(+(*(double *)(py_angy->data + (i+1)*diy + j*djy))  
				     -(*(double *)(py_angy->data + (i-1)*diy + j*djy)) 
				     +(*(double *)(py_angx->data + i*diy + (j+1)*djy)) 
				     +(*(double *)(py_angx->data + i*diy + (j-0)*djy)) ));
	      p[index] = (float)(*(double *)(py_pist->data + i*di + j*dj));
	      break;
	    case 7://(sy[index+isegnum]-sy[index-isegnum]-sx[index]-sx[index-1])*h
	      ss[index]=(float)(h/3*(+(*(double *)(py_angy->data + (i+1)*diy + j*djy))  
				     -(*(double *)(py_angy->data + (i-1)*diy + j*djy)) 
				     -(*(double *)(py_angx->data + i*diy + (j+0)*djy)) 
				     -(*(double *)(py_angx->data + i*diy + (j-1)*djy)) ));
	      p[index] = (float)(*(double *)(py_pist->data + i*di + j*dj));
	      break;
	    case 8://(sy[index+isegnum]+sy[index]+sx[index+1]-sx[index-1])*h
	      ss[index]=(float)(h/3*(+(*(double *)(py_angy->data + (i+1)*diy + j*djy))  
				     +(*(double *)(py_angy->data + (i-0)*diy + j*djy)) 
				     +(*(double *)(py_angx->data + i*diy + (j+1)*djy)) 
				     -(*(double *)(py_angx->data + i*diy + (j-1)*djy)) ));
	      p[index] = (float)(*(double *)(py_pist->data + i*di + j*dj));
	      break;
	    case 9://(-sy[index]-sy[index-isegnum]+sx[index+1]-sx[index-1])*h
	      ss[index]=(float)(h/3*(-(*(double *)(py_angy->data + (i+0)*diy + j*djy))  
				     -(*(double *)(py_angy->data + (i-1)*diy + j*djy)) 
				     +(*(double *)(py_angx->data + i*diy + (j+1)*djy)) 
				     -(*(double *)(py_angx->data + i*diy + (j-1)*djy)) ));
	      p[index] = (float)(*(double *)(py_pist->data + i*di + j*dj));
	      break;
	    case 1://(sy[index+isegnum]-sy[index-isegnum]+sx[index+1]-sx[index-1])*h
	      ss[index]=(float)(h/4*(+(*(double *)(py_angy->data + (i+1)*diy + j*djy))  
				     -(*(double *)(py_angy->data + (i-1)*diy + j*djy)) 
				     +(*(double *)(py_angx->data + i*diy + (j+1)*djy)) 
				     -(*(double *)(py_angx->data + i*diy + (j-1)*djy)) ));
	      p[index] = (float)(*(double *)(py_pist->data + i*di + j*dj));
	      break;
	    default:
	      ss[index]=0.;
	      p[index]=0.;
	      break;
	  }
	}
      }
    }else if(py_angx->descr->type_num==PyArray_FLOAT && py_angy->descr->type_num==PyArray_FLOAT && py_pist->descr->type_num==PyArray_FLOAT){
      for(i=0;i<isegnum;++i){
	for(j=0;j<isegnum;++j){
	  index=i*isegnum+j;
	  switch(*(int*)(sormask+index*dimirr)){
	    case 2://(sy[index+isegnum]+sy[index]+sx[index+1]+sx[index])*h
	      ss[index]=(float)(h/2*(+(*(float *)(py_angy->data + (i+1)*diy + j*djy))  
				     +(*(float *)(py_angy->data + (i+0)*diy + j*djy)) 
				     +(*(float *)(py_angx->data + i*diy + (j+1)*djy)) 
				     +(*(float *)(py_angx->data + i*diy + (j+0)*djy)) ));
	      p[index] = (float)(*(float *)(py_pist->data + i*di + j*dj));
	      break;
	    case 3://(-sy[index]-sy[index-isegnum]+sx[index+1]+sx[index])*h
	      ss[index]=(float)(h/2*(-(*(float *)(py_angy->data + (i+0)*diy + j*djy))  
				     -(*(float *)(py_angy->data + (i-1)*diy + j*djy)) 
				     +(*(float *)(py_angx->data + i*diy + (j+1)*djy)) 
				     +(*(float *)(py_angx->data + i*diy + (j+0)*djy)) ));
	      p[index] = (float)(*(float *)(py_pist->data + i*di + j*dj));
	      break;
	    case 4://(sy[index+isegnum]+sy[index]-sx[index]-sx[index-1])*h
	      ss[index]=(float)(h/2*(+(*(float *)(py_angy->data + (i+1)*diy + j*djy))  
				     +(*(float *)(py_angy->data + (i+0)*diy + j*djy)) 
				     -(*(float *)(py_angx->data + i*diy + (j+0)*djy)) 
				     -(*(float *)(py_angx->data + i*diy + (j-1)*djy)) ));
	      p[index] = (float)(*(float *)(py_pist->data + i*di + j*dj));
	      break;
	    case 5://(-sy[index]-sy[index-isegnum]-sx[index]-sx[index-1])*h
	      ss[index]=(float)(h/2*(-(*(float *)(py_angy->data + (i+0)*diy + j*djy))  
				     -(*(float *)(py_angy->data + (i-1)*diy + j*djy)) 
				     -(*(float *)(py_angx->data + i*diy + (j+0)*djy)) 
				     -(*(float *)(py_angx->data + i*diy + (j-1)*djy)) ));
	      p[index] = (float)(*(float *)(py_pist->data + i*di + j*dj));
	      break;
	    case 6://(sy[index+isegnum]-sy[index-isegnum]+sx[index+1]+sx[index])*h
	      ss[index]=(float)(h/3*(+(*(float *)(py_angy->data + (i+1)*diy + j*djy))  
				     -(*(float *)(py_angy->data + (i-1)*diy + j*djy)) 
				     +(*(float *)(py_angx->data + i*diy + (j+1)*djy)) 
				     +(*(float *)(py_angx->data + i*diy + (j-0)*djy)) ));
	      p[index] = (float)(*(float *)(py_pist->data + i*di + j*dj));
	      break;
	    case 7://(sy[index+isegnum]-sy[index-isegnum]-sx[index]-sx[index-1])*h
	      ss[index]=(float)(h/3*(+(*(float *)(py_angy->data + (i+1)*diy + j*djy))  
				     -(*(float *)(py_angy->data + (i-1)*diy + j*djy)) 
				     -(*(float *)(py_angx->data + i*diy + (j+0)*djy)) 
				     -(*(float *)(py_angx->data + i*diy + (j-1)*djy)) ));
	      p[index] = (float)(*(float *)(py_pist->data + i*di + j*dj));
	      break;
	    case 8://(sy[index+isegnum]+sy[index]+sx[index+1]-sx[index-1])*h
	      ss[index]=(float)(h/3*(+(*(float *)(py_angy->data + (i+1)*diy + j*djy))  
				     +(*(float *)(py_angy->data + (i-0)*diy + j*djy)) 
				     +(*(float *)(py_angx->data + i*diy + (j+1)*djy)) 
				     -(*(float *)(py_angx->data + i*diy + (j-1)*djy)) ));
	      p[index] = (float)(*(float *)(py_pist->data + i*di + j*dj));
	      break;
	    case 9://(-sy[index]-sy[index-isegnum]+sx[index+1]-sx[index-1])*h
	      ss[index]=(float)(h/3*(-(*(float *)(py_angy->data + (i+0)*diy + j*djy))  
				     -(*(float *)(py_angy->data + (i-1)*diy + j*djy)) 
				     +(*(float *)(py_angx->data + i*diy + (j+1)*djy)) 
				     -(*(float *)(py_angx->data + i*diy + (j-1)*djy)) ));
	      p[index] = (float)(*(float *)(py_pist->data + i*di + j*dj));
	      break;
	    case 1://(sy[index+isegnum]-sy[index-isegnum]+sx[index+1]-sx[index-1])*h
	      ss[index]=(float)(h/4*(+(*(float *)(py_angy->data + (i+1)*diy + j*djy))  
				     -(*(float *)(py_angy->data + (i-1)*diy + j*djy)) 
				     +(*(float *)(py_angx->data + i*diy + (j+1)*djy)) 
				     -(*(float *)(py_angx->data + i*diy + (j-1)*djy)) ));
	      p[index] = (float)(*(float *)(py_pist->data + i*di + j*dj));
	      break;
	    default:
	      ss[index]=0.;
	      p[index]=0.;
	      break;
	  }
	}
      }
    }else{
      printf("sormodule - array inputs should be of same type\n");
      return NULL;
    }

    
 
    
    iter=0;
    err=1000.0;

    while(err>conv && iter<maxiters) {
      /*printf("... %d %f\n",iter,err);*/
      ns=0;   
      err=0;
      //for (iy=0;iy<isegnum;iy++) {
      //for (ix=0;ix<isegnum;ix++) {
      //  index=iy*isegnum+ix;
      for(index=0; index<isegnum2; index++){
	  /*  corners  */
	  switch((*(int*)(sormask+index*dimirr))){
	    case 2:
	      diff=w*((p[index+1]+p[index+isegnum])/2 
		      +ss[index]
		      -p[index]);
	      p[index]+=diff;
	      break;
	    case 3:
	      diff=w*((p[index-isegnum]+p[index+1])/2
		      +ss[index]
		      -p[index]);
	      p[index]+=diff;
	      break;
	    case 4:
	      diff=w*((p[index+isegnum]+p[index-1])/2
		      +ss[index]
		      -p[index]);
	      p[index]+=diff;
	      break;
	    case 5:
	      diff=w*((p[index-isegnum]+p[index-1])/2
		      +ss[index]
		      -p[index]);
	      p[index]+=diff;
	      break;
	      /*     sides   */
	    case 6:
	      diff=w*((p[index+isegnum]+p[index-isegnum]+p[index+1])/3
		      +ss[index]
		      -p[index]);
	      p[index]+=diff;
	      break;
	    case 7:
	      diff=w*((p[index+isegnum]+p[index-isegnum]+p[index-1])/3
		      +ss[index]
		      -p[index]);
	      p[index]+=diff;
	      break;
	    case 8:
	      diff=w*((p[index+isegnum]+p[index+1]+p[index-1])/3
		      +ss[index]
		      -p[index]);
	      p[index]+=diff;
	      break;
	    case 9:
	      diff=w*((p[index-isegnum]+p[index+1]+p[index-1])/3
		      +ss[index]
		      -p[index]);
	      p[index]+=diff;
	      break;
	      /*     main body   */
	    case 1:
	      diff=w*((p[index+1]+p[index-1]+p[index+isegnum]+p[index-isegnum])/4
		      +ss[index]
		      -p[index]);
	      p[index]+=diff;
	      break;
	    default:
	      diff=0.;
	      break;
	  }
	  if((*(int*)(sormask+index*dimirr))!=0) {  
	    err+=fabs(diff); 
	    ns++;
	  }  
	  //}
      }
      err=err/(float)ns;
      iter++;
    }
    avpist=0;
    //for (iy=0;iy<isegnum;iy++) {
    //for (ix=0;ix<isegnum;ix++) {
    //index=iy*isegnum+ix;
    for(index=0; index<isegnum2; index++){
      avpist+=p[index];
      //}
    } 
    avpist/=avpistDivisor;//736.0;
    /* populate output Python array with piston values */
    if(py_pist->descr->type_num==PyArray_DOUBLE){
      for(i=0;i<isegnum;++i){
	for(j=0;j<isegnum;++j){
	  *(double *)(py_pist->data+i*di+j*dj) = (double)((p[i*isegnum+j]-avpist)/8.);
	}
      }
    }else{
      for(i=0;i<isegnum;++i){
	for(j=0;j<isegnum;++j){
	  *(float *)(py_pist->data+i*di+j*dj)=((p[i*isegnum+j]-avpist)/8.);
	}
      }
    }
    return Py_BuildValue("if",iter,err);
}


static PyObject *sor_fitOld1(PyObject *self,PyObject *args){
    PyArrayObject *py_angx,*py_angy,*py_pist,*imirr;
    PyArrayObject *tmparr=NULL;
    int i,j,ix,iy,dimirr;
    int index,ns;
    int iter;
    int ndx,dix,djx,dimsx[2],ndy,diy,djy,dimsy[2],nd,di,dj,dims[2];
    float w;
    float err;
    float avpist;
    float diff;
    float h;
    float *p;
    float *sx;
    float *sy;
    double conv=0.01;
    double avpistDivisor;
    int isegnum;
    int maxiters=1000;
    if (!PyArg_ParseTuple (args, "O!O!O!O!d|diO!", &PyArray_Type, &py_angx, &PyArray_Type, &py_angy, &PyArray_Type, &py_pist, &PyArray_Type,&imirr,&avpistDivisor,&conv,&maxiters,&PyArray_Type,&tmparr)){
	printf("Need centx, centy, piston (output), imirr (mask), avpistDivisor, conv(optional),maxiters(optional),temporary array (optional)\n");
	return NULL;
    }
    //printf("%g\n",conv);
    
/* get input Python array dimensions */
    dimirr=imirr->strides[0];

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
    
    isegnum=dimsx[0];
    
    if(dimsx[0]!=isegnum || dimsx[1]!=isegnum || dimsy[0]!=isegnum || dimsy[1]!=isegnum){
	printf("Tilt array wrong size.\n");
	return NULL;
    }
    if(dims[0]!=isegnum || dims[1]!=isegnum ){
	printf("Piston array wrong size.\n");
	return NULL;
    }
    if(imirr->dimensions[0]!=isegnum*isegnum){
	printf("Mirror mask is wrong shape\n");
	return NULL;
    }
    if(imirr->descr->kind!='i' || imirr->descr->elsize!=sizeof(int)){
	printf("Mirror mask should be int\n");
	return NULL;
    }
    if(avpistDivisor>isegnum*isegnum){
	printf("avpistDivisor should be double, current value %g\n",avpistDivisor);
    }
    if(tmparr==NULL){
      p=(float*)malloc(sizeof(float*)*isegnum*isegnum);
      sx=(float*)malloc(sizeof(float*)*isegnum*isegnum);
      sy=(float*)malloc(sizeof(float*)*isegnum*isegnum);
    }else{
      if(tmparr->nd==3 && tmparr->descr->type_num==PyArray_FLOAT && tmparr->dimensions[0]==3 && tmparr->dimensions[1]==isegnum && tmparr->dimensions[2]==isegnum){
	p=(float*)tmparr->data;
	sx=&(((float*)tmparr->data)[isegnum*isegnum]);
	sy=&(((float*)tmparr->data)[2*isegnum*isegnum]);
      }else{
	printf("temporary array should have dimensions (3,%d,%d), type Float32\n",isegnum,isegnum);
	return NULL;
      }
    }
    /* convert input arrays to C array, appropriately scaled. */
    if(py_angx->descr->type_num==PyArray_DOUBLE && py_angy->descr->type_num==PyArray_DOUBLE && py_pist->descr->type_num==PyArray_DOUBLE){
      for(i=0;i<isegnum;++i){
	for(j=0;j<isegnum;++j){
	  index=i*isegnum+j;
	  if((*(int*)(imirr->data+index*dimirr))!=0){
	    sx[index] = ((float)(*(double *)(py_angx->data + i*dix + j*djx)))/2;
	    sy[index] = ((float)(*(double *)(py_angy->data + i*diy + j*djy)))/2;
	    p[index] = (float)(*(double *)(py_pist->data + i*di + j*dj));
	  }else{
	    sx[index]=0.;
	    sy[index]=0.;
	    p[index]=0.;
	  }
	}
      }
    }else if(py_angx->descr->type_num==PyArray_FLOAT && py_angy->descr->type_num==PyArray_FLOAT && py_pist->descr->type_num==PyArray_FLOAT){
      for(i=0;i<isegnum;++i){
	for(j=0;j<isegnum;++j){
	  index=i*isegnum+j;
	  if((*(int*)(imirr->data+index*dimirr))!=0){
	    sx[index] = (float)(*(float *)(py_angx->data + i*dix + j*djx))/2;
	    sy[index] = (float)(*(float *)(py_angy->data + i*diy + j*djy))/2;
	    p[index] = (float)(*(float *)(py_pist->data + i*di + j*dj));
	  }else{
	    sx[index]=0.;
	    sy[index]=0.;
	    p[index]=0.;
	  }
	}
      }
    }else{
      printf("sormodule - array inputs should be of same type\n");
      return NULL;
    }

    avpist=0;
    h=(ISEGSZ)*(2.0*3.14157);
    
    w=2.0/(1.0+sin(3.14157/(isegnum+1.0)));
 
    
    iter=0;
    err=1.0;
    while(err>conv && iter<maxiters) {
      /*printf("... %d %f\n",iter,err);*/
      ns=0;   
      err=0;
      for (iy=0;iy<isegnum;iy++) {
	for (ix=0;ix<isegnum;ix++) {
	  index=iy*isegnum+ix;
	  /*  corners  */
	  switch((*(int*)(imirr->data+index*dimirr))){
	    case 2:
	      //if((*(int*)(imirr->data+index*dimirr))==2) {
	      diff=w*((p[index+1]+p[index+isegnum]
		       +(sy[index+isegnum]+sy[index]+sx[index+1]+sx[index])*h)/2
		      -p[index]);
	      p[index]+=diff;
	      //}
	      break;
	    case 3:
	      //if((*(int*)(imirr->data+index*dimirr))==3) {
	      diff=w*((p[index-isegnum]+p[index+1]
		       +(-sy[index]-sy[index-isegnum]+sx[index+1]+sx[index])*h)/2
		      -p[index]);
	      p[index]+=diff;
	      //}
	      break;
	    case 4:
	      //if((*(int*)(imirr->data+index*dimirr))==4) {
	      diff=w*((p[index+isegnum]+p[index-1]
		       +(sy[index+isegnum]+sy[index]-sx[index]-sx[index-1])*h)/2
		      -p[index]);
	      p[index]+=diff;
	      //}
	      break;
	    case 5:
	      //if((*(int*)(imirr->data+index*dimirr))==5) {
	      diff=w*((p[index-isegnum]+p[index-1]
		       +(-sy[index]-sy[index-isegnum]-sx[index]-sx[index-1])*h)/2
		      -p[index]);
	      p[index]+=diff;
	      //}
	      break;
	      /*     sides   */
	    case 6:
	      //if((*(int*)(imirr->data+index*dimirr))==6) {
	      diff=w*((p[index+isegnum]+p[index-isegnum]+p[index+1]
		       +(sy[index+isegnum]-sy[index-isegnum]+sx[index+1]+sx[index])*h)/3
		      -p[index]);
	      p[index]+=diff;
	      //}
	      break;
	    case 7:
	      //if((*(int*)(imirr->data+index*dimirr))==7) {
	      diff=w*((p[index+isegnum]+p[index-isegnum]+p[index-1]
		       +(sy[index+isegnum]-sy[index-isegnum]-sx[index]-sx[index-1])*h)/3
		      -p[index]);
	      p[index]+=diff;
	      //}
	      break;
	    case 8:
	      //if((*(int*)(imirr->data+index*dimirr))==8) {
	      diff=w*((p[index+isegnum]+p[index+1]+p[index-1]
		       +(sy[index+isegnum]+sy[index]+sx[index+1]-sx[index-1])*h)/3
		      -p[index]);
	      p[index]+=diff;
	      //}
	      break;
	    case 9:
	      //if((*(int*)(imirr->data+index*dimirr))==9) {
	      diff=w*((p[index-isegnum]+p[index+1]+p[index-1]
		       +(-sy[index]-sy[index-isegnum]+sx[index+1]-sx[index-1])*h)/3
		      -p[index]);
	      p[index]+=diff;
	      //}
	      break;
	      /*     main body   */
	    case 1:
	      //if((*(int*)(imirr->data+index*dimirr))==1) {  
	      diff=w*((p[index+1]+p[index-1]+p[index+isegnum]+p[index-isegnum]
		       +(sy[index+isegnum]-sy[index-isegnum]+sx[index+1]-sx[index-1])*h)/4
		      -p[index]);
	      p[index]+=diff;
	      //}
	      break;
	    default:
	      diff=0.;
	      break;
	  }
	  if((*(int*)(imirr->data+index*dimirr))!=0) {  
	    err+=fabs(diff); 
	    ns++;
	  }  
	}
      }
      err=err/(float)ns;
      iter++;
    }
    
    ns=0;
    for (iy=0;iy<isegnum;iy++) {
      for (ix=0;ix<isegnum;ix++) {
	index=iy*isegnum+ix;
	avpist+=p[index];
      }
    } 
    avpist/=avpistDivisor;//736.0;
    /* populate output Python array with piston values */
    if(py_pist->descr->type_num==PyArray_DOUBLE){
      for(i=0;i<isegnum;++i){
	for(j=0;j<isegnum;++j){
	  *(double *)(py_pist->data+i*di+j*dj) = (double)((p[i*isegnum+j]-avpist)/8.);
	}
      }
    }else{
      for(i=0;i<isegnum;++i){
	for(j=0;j<isegnum;++j){
	  *(float *)(py_pist->data+i*di+j*dj)=((p[i*isegnum+j]-avpist)/8.);
	}
      }
    }
    if(tmparr==NULL){
      free(p);
      free(sx);
      free(sy);
    }
    return Py_BuildValue("if",iter,err);
}


static PyObject *sor_fitOld(PyObject *self,PyObject *args){
    PyArrayObject	*py_angx,*py_angy,*py_pist,*imirr;
    int i,j,ix,iy,dimirr;
    int index,ns;
    int iter;
    int ndx,dix,djx,dimsx[2],ndy,diy,djy,dimsy[2],nd,di,dj,dims[2];
    float w;
    float err, pl;
    float avpist;
    float diff;
    float h;
    //float pist[ISEGNUM][ISEGNUM];
    //float angx[ISEGNUM][ISEGNUM],angy[ISEGNUM][ISEGNUM];
    //float p[ISEGNUM][ISEGNUM];
    //float sx[ISEGNUM][ISEGNUM],sy[ISEGNUM][ISEGNUM];
    float **pist;
    float **angx;
    float **angy;
    float **p;
    float **sx;
    float **sy;
    double conv=0.01;
    double avpistDivisor;
    int isegnum;
    int maxiters=1000;
    if (!PyArg_ParseTuple (args, "O!O!O!O!d|di", &PyArray_Type, &py_angx, &PyArray_Type, &py_angy, &PyArray_Type, &py_pist, &PyArray_Type,&imirr,&avpistDivisor,&conv,&maxiters)){
	printf("Need centx, centy, piston (output), imirr (mask), avpistDivisor, conv(optional),maxiters(optional)\n");
	return NULL;
    }
    //printf("%g\n",conv);
    
/* get input Python array dimensions */
    dimirr=imirr->strides[0];

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
    
    isegnum=dimsx[0];
    
    if(dimsx[0]!=isegnum || dimsx[1]!=isegnum || dimsy[0]!=isegnum || dimsy[1]!=isegnum){
	printf("Tilt array wrong size.\n");
	return NULL;
    }
    if(dims[0]!=isegnum || dims[1]!=isegnum ){
	printf("Piston array wrong size.\n");
	return NULL;
    }
    if(imirr->dimensions[0]!=isegnum*isegnum){
	printf("Mirror mask is wrong shape\n");
	return NULL;
    }
    if(imirr->descr->kind!='i' || imirr->descr->elsize!=sizeof(int)){
	printf("Mirror mask should be int\n");
	return NULL;
    }
    if(avpistDivisor>isegnum*isegnum){
	printf("avpistDivisor should be double, current value %g\n",avpistDivisor);
    }

    pist=(float**)malloc(sizeof(float*)*isegnum);
    angx=(float**)malloc(sizeof(float*)*isegnum);
    angy=(float**)malloc(sizeof(float*)*isegnum);
    p=(float**)malloc(sizeof(float*)*isegnum);
    sx=(float**)malloc(sizeof(float*)*isegnum);
    sy=(float**)malloc(sizeof(float*)*isegnum);
    
    for(i=0; i<isegnum; i++){
	pist[i]=(float*)malloc(sizeof(float)*isegnum);
	angx[i]=(float*)malloc(sizeof(float)*isegnum);
	angy[i]=(float*)malloc(sizeof(float)*isegnum);
	p[i]=(float*)malloc(sizeof(float)*isegnum);
	sx[i]=(float*)malloc(sizeof(float)*isegnum);
	sy[i]=(float*)malloc(sizeof(float)*isegnum);
    }


/* convert input arrays to C array */

    for(i=0;i<isegnum;++i){
	for(j=0;j<isegnum;++j){
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
    
    w=2.0/(1.0+sin(3.14157/(isegnum+1.0)));
    
    ns=0;
    for (ix=0;ix<isegnum;ix++) {
	
	for (iy=0;iy<isegnum;iy++) {
	    index=iy*isegnum+ix;
	    if((*(int*)(imirr->data+index*dimirr))!=0) {
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
    while(err>conv && iter<maxiters) {
	/*printf("... %d %f\n",iter,err);*/
	ns=0;   
	err=0;
	for (iy=0;iy<isegnum;iy++) {
	    for (ix=0;ix<isegnum;ix++) {
		
		index=iy*isegnum+ix;
		
		/*  corners  */
		switch((*(int*)(imirr->data+index*dimirr))){
		  case 2:
			//if((*(int*)(imirr->data+index*dimirr))==2) {
		    pl=p[ix][iy];
		    p[ix][iy]+=
			w*((p[ix+1][iy]+p[ix][iy+1]
			    +(sy[ix][iy+1]+sy[ix][iy]+sx[ix+1][iy]+sx[ix][iy])*h)/2
			   -p[ix][iy]);
		    diff=p[ix][iy]-pl;
		    //}
		    break;
		    case 3:
			//if((*(int*)(imirr->data+index*dimirr))==3) {
		    pl=p[ix][iy];
		    p[ix][iy]+=
			w*((p[ix][iy-1]+p[ix+1][iy]
			    +(-sy[ix][iy]-sy[ix][iy-1]+sx[ix+1][iy]+sx[ix][iy])*h)/2
			   -p[ix][iy]);
		    diff=p[ix][iy]-pl;
		    //}
		    break;
		    case 4:
			//if((*(int*)(imirr->data+index*dimirr))==4) {
		    pl=p[ix][iy];
		    p[ix][iy]+=
			w*((p[ix][iy+1]+p[ix-1][iy]
			    +(sy[ix][iy+1]+sy[ix][iy]-sx[ix][iy]-sx[ix-1][iy])*h)/2
			   -p[ix][iy]);
		    diff=p[ix][iy]-pl;
		    //}
		    break;
		    case 5:
			//if((*(int*)(imirr->data+index*dimirr))==5) {
		    pl=p[ix][iy];
		    p[ix][iy]+=
			w*((p[ix][iy-1]+p[ix-1][iy]
			    +(-sy[ix][iy]-sy[ix][iy-1]-sx[ix][iy]-sx[ix-1][iy])*h)/2
			   -p[ix][iy]);
		    diff=p[ix][iy]-pl;
		    //}
		    break;
		    /*     sides   */
		    case 6:
			//if((*(int*)(imirr->data+index*dimirr))==6) {
		    pl=p[ix][iy];
		    p[ix][iy]+=
			w*((p[ix][iy+1]+p[ix][iy-1]+p[ix+1][iy]
			    +(sy[ix][iy+1]-sy[ix][iy-1]+sx[ix+1][iy]+sx[ix][iy])*h)/3
			   -p[ix][iy]);
		    diff=p[ix][iy]-pl;
		    //}
		    break;
		    case 7:
			//if((*(int*)(imirr->data+index*dimirr))==7) {
		    pl=p[ix][iy];
		    p[ix][iy]+=
			w*((p[ix][iy+1]+p[ix][iy-1]+p[ix-1][iy]
			    +(sy[ix][iy+1]-sy[ix][iy-1]-sx[ix][iy]-sx[ix-1][iy])*h)/3
			   -p[ix][iy]);
		    diff=p[ix][iy]-pl;
		    //}
		    break;
		    case 8:
			//if((*(int*)(imirr->data+index*dimirr))==8) {
		    pl=p[ix][iy];
		    p[ix][iy]+=
			w*((p[ix][iy+1]+p[ix+1][iy]+p[ix-1][iy]
			    +(sy[ix][iy+1]+sy[ix][iy]+sx[ix+1][iy]-sx[ix-1][iy])*h)/3
			   -p[ix][iy]);
		    diff=p[ix][iy]-pl;
		    //}
		    break;
		    case 9:
			//if((*(int*)(imirr->data+index*dimirr))==9) {
		    pl=p[ix][iy];
		    p[ix][iy]+=
			w*((p[ix][iy-1]+p[ix+1][iy]+p[ix-1][iy]
			    +(-sy[ix][iy]-sy[ix][iy-1]+sx[ix+1][iy]-sx[ix-1][iy])*h)/3
			   -p[ix][iy]);
		    diff=p[ix][iy]-pl;
		    //}
		    break;
		    /*     main body   */
		    case 1:
			//if((*(int*)(imirr->data+index*dimirr))==1) {  
		    pl=p[ix][iy];
		    p[ix][iy]+=
			w*((p[ix][iy+1]+p[ix][iy-1]+p[ix+1][iy]+p[ix-1][iy]
			    +(sy[ix][iy+1]-sy[ix][iy-1]+sx[ix+1][iy]-sx[ix-1][iy])*h)/4
			   -p[ix][iy]);
		    diff=p[ix][iy]-pl;
		    //}
		    break;
		    default:
		    break;
		}
		if((*(int*)(imirr->data+index*dimirr))!=0) {  
		    err+=fabs(diff); 
		    ns++;
		    
		}  
	    }
	}
	
	err=err/(float)ns;
	
	iter++;
    }
    
    ns=0;
    for (iy=0;iy<isegnum;iy++) {
	for (ix=0;ix<isegnum;ix++) {
	    index=iy*isegnum+ix;
	    pist[ix][iy]=p[ix][iy];
	    /* avpist+=pist[ix][iy]*pist[ix][iy] */;     
	    avpist+=pist[ix][iy];
	}
    } 
    avpist=avpist/avpistDivisor;//736.0;
    /* printf("%f\n", avpist); */
    
    
    
    
    /* populate output Python array with piston values */
    
    for(i=0;i<isegnum;++i){
	for(j=0;j<isegnum;++j){
	    /* *(double *)(py_pist->data+i*di+j*dj) = (double)((pist[j][i]-avpist)/(double)ISEGSZ); */
	    *(double *)(py_pist->data+i*di+j*dj) = (double)((pist[j][i]-avpist)/8.);
	    /* printf("%d %d %d %f %f %f\n",i,j,(*(int*)(imirr->data+index*dimirr)),angx[i][j],angy[i][j],pist[i][j]); */
	    
	}
    }
    

    for(i=0; i<isegnum; i++){
	free(pist[i]);
	free(angx[i]);
	free(angy[i]);
	free(p[i]);
	free(sx[i]);
	free(sy[i]);
    }

    free(pist);
    free(angx);
    free(angy);
    free(p);
    free(sx);
    free(sy);

    
    return Py_BuildValue("if",iter,err);
    
    
    
}




/* =============================================================================== */



/* define a methods table for this module */

static PyMethodDef sor_methods[] = 	{
					{"init", sor_init, METH_VARARGS}, 
					{"setconv", sor_setconv, METH_VARARGS}, 
					{"free", sor_free, METH_VARARGS}, 
					{"fit", sor_fit, METH_VARARGS}, 
					{"fitOld1", sor_fitOld1, METH_VARARGS}, 
					{"fitOld", sor_fitOld, METH_VARARGS}, 
					{NULL, NULL} };


/* initialisation - register the methods with the Python interpreter */

void initsor()
{
	(void) Py_InitModule("sor", sor_methods);
	import_array();
}
