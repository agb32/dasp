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
#include <Python.h>
#include <numpy/arrayobject.h>
#include <stringobject.h>
#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include <errno.h>
#include <stdlib.h>
#include <unistd.h>
#include <limits.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <pthread.h>
#include <assert.h>
#include <math.h>
#include <sys/mman.h>

#include "mvm.h"

/* Added by NAB 04/Apr/2013
 * reason: functions X-64 are not defined, replaced with synonyms to the normal
 * functions which are 64-bit aware.
 */
#ifdef __APPLE__
   /* Darwin/FreeBSD will have 64bit support built in so can define
    * macros here */
  #define _DARWIN_USE_64_BIT_INODE 
  #define open64 open
  #define truncate64 truncate
  #define stat64 stat
#endif
/* End of add
 */

static PyObject *UtilsError;

typedef struct{
  int dimensions[2];
  int strides[2];
  int ostrides[2];
  float *odata;
  float *idata;
  float c;//cos of angle.
  float s;//sin of angle
}rotateStruct;
typedef struct{
  int ystart;
  int yend;
  rotateStruct *rs;
}rotateStructThread;//per-thread information (start and end).


//for single threaded, ystart=0, yend=dimensions[0].
int rotateWorker(rotateStructThread *rst){
  rotateStruct *rs=rst->rs;
  int yy,xx,x1,x2,y1,y2;
  float x,y;
  float xold,yold;
  float xm,ym;
  float s=rs->s;
  float c=rs->c;
  int *dimensions=rs->dimensions;
  int *ostrides=rs->ostrides;
  int *strides=rs->strides;
  int yend=rst->yend;
  float val;
  float *idata=rs->idata;
  float *odata=rs->odata;
  for(yy=rst->ystart; yy<yend; yy++){
    y=yy-dimensions[0]/2.+0.5;
    for(xx=0; xx<dimensions[1]; xx++){
      x=xx-dimensions[1]/2.+0.5;
      yold=-s*x+c*y+dimensions[0]/2.-0.5;
      xold=c*x+s*y+dimensions[1]/2.-0.5;
      x1=(int)floor(xold);
      x2=x1+1;
      y1=(int)floor(yold);
      y2=y1+1;
      xm=xold-x1;
      ym=yold-y1;
      if(y2==dimensions[0]){
	y2=y1;
	ym=0;
      }
      if(x2==dimensions[1]){
	x2=x1;
	xm=0;
      }
      if(x1==-1){
	x1=0;
	xm=0;
      }
      if(y1==-1){
	y1=0;
	ym=0;
      }
      val=0.;
      if(y1>=0 && y1<dimensions[0]){
	if(x1>=0 && x1<dimensions[1])
	  val+=idata[y1*strides[0]+x1*strides[1]]*(1-xm)*(1-ym);
	if(x2<dimensions[1] && x2>=0)
	  val+=idata[y1*strides[0]+x2*strides[1]]*xm*(1-ym);
      }
      if(y2<dimensions[0] && y2>=0){
	if(x2<dimensions[1] && x2>=0)
	  val+=idata[y2*strides[0]+x2*strides[1]]*xm*ym;
	if(x1>=0 && x1<dimensions[1])
	  val+=idata[y2*strides[0]+x1*strides[1]]*(1-xm)*ym;
      }
      odata[yy*ostrides[0]+xx*ostrides[1]]=val;
    }
  }
  return 0;
}


//Rotates a 2D array.
static PyObject* rotateArray(PyObject *self,PyObject *args){
  PyArrayObject *inarr,*outarr=NULL;
  double angle;
  int i,x,y;
  int nthreads=0;
  float *idata,*odata;
  int idx,idy,odx,ody;
  rotateStruct rs;
  rotateStructThread *rst;
  pthread_t *thread;
  int ystart;
  if(!PyArg_ParseTuple(args,"O!d|O!i",&PyArray_Type,&inarr,&angle,&PyArray_Type,&outarr,&nthreads)){
    printf("Usage: numpy array, rotation angle, output array (optional - if not provided, will overwrite input).\n");
    return NULL;
  }
  if(inarr->nd!=2){
    printf("Input array must be 2D\n");
    return NULL;
  }
  if(inarr->descr->type_num!=NPY_FLOAT){
    printf("Input array must be float32\n");
    return NULL;
  }
  if(nthreads<=0){//get number of CPUs.
    nthreads=sysconf(_SC_NPROCESSORS_ONLN);
  }
  //printf("Using %d processors, angle %g\n",nthreads,angle);

  if(outarr!=NULL){
    if(outarr->nd!=2){
      printf("Output array must be 2D\n");
      return NULL;
    }
    if(outarr->descr->type_num!=NPY_FLOAT){
      printf("Output array must be float32\n");
      return NULL;
    }
    if(inarr->dimensions[0]!=outarr->dimensions[0] || inarr->dimensions[1]!=outarr->dimensions[1]){
      printf("Arrays must be same shape\n");
      return NULL;
    }
    odata=(float*)outarr->data;
    odx=outarr->strides[1]/sizeof(float);
    ody=outarr->strides[0]/sizeof(float);
  }else{
    if((odata=calloc(sizeof(float),inarr->dimensions[0]*inarr->dimensions[1]))==NULL){
      printf("Failed to allocate memory for temporary output\n");
      return NULL;
    }
    odx=1;
    ody=inarr->dimensions[1];
  }
  if((rst=malloc(sizeof(rotateStructThread)*nthreads))==NULL){
    printf("failed to malloc rst\n");
    if(outarr==NULL)
      free(odata);
    return NULL;
  }
  if((thread=malloc(sizeof(pthread_t)*nthreads))==NULL){
    printf("failed to malloc thread\n");
    free(rst);
    if(outarr==NULL)
      free(odata);
    return NULL;
  }


  idata=(float*)inarr->data;
  idx=inarr->strides[1]/sizeof(float);
  idy=inarr->strides[0]/sizeof(float);


  rs.dimensions[0]=inarr->dimensions[0];
  rs.dimensions[1]=inarr->dimensions[1];
  rs.strides[0]=idy;
  rs.strides[1]=idx;
  rs.ostrides[0]=ody;
  rs.ostrides[1]=odx;
  rs.odata=odata;
  rs.idata=idata;
  rs.c=(float)cos(angle*M_PI/180.);
  rs.s=(float)sin(angle*M_PI/180.);
  ystart=0;
  for(i=0; i<nthreads; i++){
    rst[i].rs=&rs;
    rst[i].ystart=ystart;
    rst[i].yend=ystart+(inarr->dimensions[0]-ystart)/(nthreads-i);
    ystart=rst[i].yend;
    pthread_create(&thread[i],NULL,(void*)rotateWorker,&rst[i]);
  }

  for(i=0; i<nthreads; i++){
    pthread_join(thread[i],NULL);
  }
    
  if(outarr==NULL){//copy from odata into idata.
    for(y=0; y<inarr->dimensions[0]; y++){
      for(x=0; x<inarr->dimensions[1]; x++){
	idata[y*idy+x*idx]=odata[y*ody+x*odx];
      }
    }
    free(odata);
  }
  free(rst);
  free(thread);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* CompressFloatArray(PyObject *self,PyObject *args){
  //Compress a float array by clipping the mantissa.  This is inplace, and the input array is destroyed.  Note, compression is lossy.  After this routine has completed, the input array should be viewed as type uint32.
  /*
import cmod.utils
import time
size=16384
steps=16
n=size/steps
if size%steps!=0:
 raise Exception("steps")

a=(numpy.random.random(size).astype("f")-0.5)*200
aa=a.copy()
au=numpy.zeros((n,),"f")
bits=16
cmod.utils.compressFloatArray(aa,bits)
t1=time.time()
for i in xrange(steps):
 cmod.utils.uncompressFloatArray(aa,au,bits,n*i)

t2=time.time()
au
a[-n:]
mm((a[-n:]-au)/a[-n:])
print t2-t1
  */
  PyArrayObject *arr;
  int bits;
  long i,pos;
  unsigned int mask;
  int offset,left;
  unsigned int *data;
  if(!PyArg_ParseTuple(args,"O!i",&PyArray_Type,&arr,&bits)){
    printf("Usage: numpy array, number of bits to compress each element too.\nThis destroys the input.\n");
    return NULL;
  }
  if(!PyArray_ISCONTIGUOUS(arr)){
    printf("Input array must be contiguous\n");
    return NULL;
  }
  if(arr->nd!=1){
    printf("Input array must be 1D\n");
    return NULL;
  }
  if(arr->descr->type_num!=NPY_FLOAT){
    printf("Array must be float32\n");
    return NULL;
  }
  if(bits<10){
    printf("Need at least 9 bits for sign and exponent\n");
    return NULL;
  }
  data=(unsigned int*)arr->data;
  mask=((1<<bits)-1)<<(32-bits);
  for(i=0; i<arr->dimensions[0]; i++){
    pos=(i*bits)/32;//current word to write to.
    offset=(i*bits)%32;//starting point in current word.
    left=bits-(32-offset);//no of bits needed to write into new word.
    if(offset==0)
      data[pos]=data[i]&mask;
    else
      data[pos]|=(data[i]&mask)>>offset;
    if(left>0)
      data[pos+1]=(data[i]&mask)<<(bits-left);
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *UncompressFloatArray(PyObject *self,PyObject *args){
  long off,i,pos;
  int bits,offset,left;
  unsigned int mask=0;
  PyArrayObject *arr,*out;
  unsigned int *data,*odata;
  if(!PyArg_ParseTuple(args,"O!O!il",&PyArray_Type,&arr,&PyArray_Type,&out,&bits,&off)){
    printf("Usage: numpy compressed array,numpy output array number of bits to compress each element too, offset index in compressed array.\n");
    return NULL;
  }
  if(!PyArray_ISCONTIGUOUS(arr)){
    printf("Input array must be contiguous\n");
    return NULL;
  }
  if(arr->nd!=1){
    printf("Input array must be 1D\n");
    return NULL;
  }
  if(arr->descr->elsize!=sizeof(int)){
    printf("Input array must be uint32 (or float, or anything 4 bytes)\n");
    return NULL;
  }
  if(!PyArray_ISCONTIGUOUS(out)){
    printf("Output array must be contiguous\n");
    return NULL;
  }
  if(out->nd!=1){
    printf("Output array must be 1D\n");
    return NULL;
  }
  if(out->descr->type_num!=NPY_FLOAT){
    printf("Output array must be float32\n");
    return NULL;
  }
  //for(i=0; i<bits; i++)
  //  mask|=1<<(31-i);
  mask=((1<<bits)-1)<<(32-bits);
  odata=(unsigned int*)out->data;
  data=(unsigned int*)arr->data;
  for(i=0; i<out->dimensions[0]; i++){
    pos=((i+off)*bits)/32;//word to read from
    offset=((i+off)*bits)%32;//offset in word to read from.
    left=bits-(32-offset);//no of bits needed from next word.
    odata[i]=(data[pos]<<offset);
    if(left>0)
      odata[i]|=(data[pos+1]>>(bits-left));
    odata[i]&=mask;
  }
  Py_INCREF(Py_None);
  return Py_None;
}

typedef struct{
  unsigned int *data;
  unsigned int *odata;
  int bits;
  long n;
  long off;
  int expMin;
  int expMax;
}UfwInfo;

int uncompressFloatWorker(UfwInfo *ufwInfo){//unsigned int *data, unsigned int *odata, int bits, long n, long off){
  //run as a worker thread...
  unsigned int mask;
  int offset,left;
  long pos,i;
  unsigned int *data=ufwInfo->data;
  unsigned int *odata=ufwInfo->odata;
  int bits=ufwInfo->bits;
  long n=ufwInfo->n;
  long off=ufwInfo->off;
  mask=((1<<bits)-1)<<(32-bits);
  for(i=0; i<n; i++){
    pos=((i+off)*bits)/32;//word to read from
    offset=((i+off)*bits)%32;//offset in word to read from.
    left=bits-(32-offset);//no of bits needed from next word.
    odata[i]=(data[pos]<<offset);
    if(left>0)
      odata[i]|=(data[pos+1]>>(bits-left));
    odata[i]&=mask;
  }
  return 0;
}

static PyObject *UncompressFloatArrayThreaded(PyObject *self,PyObject *args){
  long off,i;//,pos;
  int bits;//,offset,left;
  long ndone,nleft;
  //unsigned int mask=0;
  PyArrayObject *arr,*out;
  //unsigned int *data,*odata;
  int nthreads;
  UfwInfo *ufwInfo;
  pthread_t *thread;
  if(!PyArg_ParseTuple(args,"O!O!ili",&PyArray_Type,&arr,&PyArray_Type,&out,&bits,&off,&nthreads)){
    printf("Usage: numpy compressed array,numpy output array number of bits to compress each element too, offset index in compressed array, nthreads.\n");
    return NULL;
  }
  if((thread=malloc(sizeof(pthread_t)*nthreads))==NULL){
    printf("cmod.utils unable to malloc thread\n");
    return NULL;
  }
  if((ufwInfo=malloc(sizeof(UfwInfo)*nthreads))==NULL){
    printf("cmod.utils unable to malloc ufwInfo\n");
    return NULL;
  }
  if(!PyArray_ISCONTIGUOUS(arr)){
    printf("Input array must be contiguous\n");
    return NULL;
  }
  if(arr->nd!=1){
    printf("Input array must be 1D\n");
    return NULL;
  }
  if(arr->descr->elsize!=sizeof(int)){
    printf("Input array must be uint32 (or float, or anything 4 bytes)\n");
    return NULL;
  }
  if(!PyArray_ISCONTIGUOUS(out)){
    printf("Output array must be contiguous\n");
    return NULL;
  }
  if(out->nd!=1){
    printf("Output array must be 1D\n");
    return NULL;
  }
  if(out->descr->type_num!=NPY_FLOAT){
    printf("Output array must be float32\n");
    return NULL;
  }
  //for(i=0; i<bits; i++)
  //  mask|=1<<(31-i);
  nleft=out->dimensions[0];
  ndone=0;
  for(i=0; i<nthreads; i++){
    ufwInfo[i].data=(unsigned int*)arr->data;
    ufwInfo[i].bits=bits;
    ufwInfo[i].n=nleft/(nthreads-i);
    ufwInfo[i].off=off+ndone;
    ufwInfo[i].odata=&(((unsigned int*)out->data)[ndone]);

    nleft-=ufwInfo[i].n;
    ndone+=ufwInfo[i].n;
    pthread_create(&thread[i],NULL,(void*)uncompressFloatWorker,&ufwInfo[i]);


  }
  for(i=0; i<nthreads; i++){
    pthread_join(thread[i],NULL);
  }

  free(ufwInfo);
  free(thread);
  Py_INCREF(Py_None);
  return Py_None;
}


static PyObject* CompressFloatArrayAll(PyObject *self,PyObject *args){
  //Compress a float array by clipping the mantissa.  This is inplace, and the input array is destroyed.  Note, compression is lossy.  After this routine has completed, the input array should be viewed as type uint32.
  //This version also compresses the exponent, giving a further saving.  Exponent compression is not lossy.
  PyArrayObject *arr;
  int bits;
  long i,pos;
  unsigned int mask;
  unsigned int word;
  int offset,left,ebits,mbits,min,max;
  unsigned int *data;
  if(!PyArg_ParseTuple(args,"O!i",&PyArray_Type,&arr,&mbits)){
    printf("Usage: numpy array, number of bits to compress each mantissa element too (less than 23).\nThis destroys the input.\n");
    return NULL;
  }
  if(!PyArray_ISCONTIGUOUS(arr)){
    printf("Input array must be contiguous\n");
    return NULL;
  }
  if(arr->nd!=1){
    printf("Input array must be 1D\n");
    return NULL;
  }
  if(arr->descr->type_num!=NPY_FLOAT){
    printf("Array must be float32\n");
    return NULL;
  }
  if(mbits>22){
    printf("Specified mantissa bits>22\n");
    return NULL;
  }
  //First calculate how many bits will be required for exponent.
  min=0xff;
  max=0;
  data=(unsigned int*)arr->data;
  for(i=0; i<arr->dimensions[0]; i++){
    offset=(data[i]>>23)&0xff;
    if(offset>0 && offset<min)
      min=offset;
    if(offset>max)
      max=offset;
  }
  offset=max-min+1;//this is the number of different exponent values (the +1 is for the zero).
  ebits=1;
  while((offset>>ebits)>0){
    ebits++;
  }
  bits=mbits+ebits+1;
  printf("Got exponent ranging from %d to %d, will fit into %d bits.  Overall compression into %d bits\n",min,max,ebits,bits);

  //mask=((1<<bits)-1)<<(32-bits);
  mask=((1<<mbits)-1)<<(23-mbits);
  for(i=0; i<arr->dimensions[0]; i++){
    pos=(i*bits)/32;//current word to write to.
    offset=(i*bits)%32;//starting point in current word.
    left=bits-(32-offset);//no of bits needed to write into new word.
    //Compress data into the msb of word.
    if(((data[i]>>23)&0xff)==0)
      word=(data[i]&0x80000000)|((data[i]&mask)<<(8-ebits));
    else
      word=(data[i]&0x80000000)|((((data[i]>>23)&0xff)-min+1)<<(31-ebits))|((data[i]&mask)<<(8-ebits));
    if(offset==0)
      data[pos]=word;//data[i]&mask;
    else
      data[pos]|=word>>offset;
    if(left>0)
      data[pos+1]=(word)<<(bits-left);
  }
  return Py_BuildValue("ii",min,max);
}
static PyObject *UncompressFloatArrayAll(PyObject *self,PyObject *args){
  long off,i,pos;
  int bits,offset,left,mbits,min,max,ebits;
  unsigned int mask=0,word,emask;
  PyArrayObject *arr,*out;
  unsigned int *data,*odata;
  if(!PyArg_ParseTuple(args,"O!O!iiil",&PyArray_Type,&arr,&PyArray_Type,&out,&mbits,&min,&max,&off)){
    printf("Usage: numpy compressed array,numpy output array number of mantissa bits to compress each element too, exponent min, exponent max, offset index in compressed array.\n");
    return NULL;
  }
  if(!PyArray_ISCONTIGUOUS(arr)){
    printf("Input array must be contiguous\n");
    return NULL;
  }
  if(arr->nd!=1){
    printf("Input array must be 1D\n");
    return NULL;
  }
  if(arr->descr->elsize!=sizeof(int)){
    printf("Input array must be uint32 (or float, or anything 4 bytes)\n");
    return NULL;
  }
  if(!PyArray_ISCONTIGUOUS(out)){
    printf("Output array must be contiguous\n");
    return NULL;
  }
  if(out->nd!=1){
    printf("Output array must be 1D\n");
    return NULL;
  }
  if(out->descr->type_num!=NPY_FLOAT){
    printf("Output array must be float32\n");
    return NULL;
  }
  if(mbits>22){
    printf("Specified mantissa bits>22\n");
    return NULL;
  }
  offset=max-min+1;
  ebits=1;
  while((offset>>ebits)>0){
    ebits++;
  }

  bits=mbits+ebits+1;
  //for(i=0; i<bits; i++)
  //  mask|=1<<(31-i);
  //mask=((1<<bits)-1)<<(32-bits);
  mask=((1<<mbits)-1)<<(23-mbits);
  odata=(unsigned int*)out->data;
  data=(unsigned int*)arr->data;
  emask=((1<<ebits)-1);
  for(i=0; i<out->dimensions[0]; i++){
    pos=((i+off)*bits)/32;//word to read from
    offset=((i+off)*bits)%32;//offset in word to read from.
    left=bits-(32-offset);//no of bits needed from next word.
    word=(data[pos]<<offset);
    if(left>0)
      word|=data[pos+1]>>(bits-left);
    if(((word>>(31-ebits))&emask)==0)
      odata[i]=(word&0x80000000)|((word>>(8-ebits))&mask);
    else
      odata[i]=(word&0x80000000)|((((word>>(31-ebits))&emask)+min-1)<<(23))|((word>>(8-ebits))&mask);
    //odata[i]=(data[pos]<<offset);
    //if(left>0)
    //  odata[i]|=(data[pos+1]>>(bits-left));
    //odata[i]&=mask;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
/*
import util.FITS
rmx=util.FITS.Read("/var/ali/spmx2layer42mb104dense_rmxden1e-05.fits")[1]
rrmx=rmx.ravel().view(numpy.uint32)
rmxr=rmx.copy().ravel()
out=numpy.zeros(rmx.shape,"f")
import cmod.utils
mbits=20
mn,mx=cmod.utils.compressFloatArrayAll(rmxr,mbits)
cmod.utils.uncompressFloatArrayAll(rmxr,out.ravel(),mbits,mn,mx,0)
out
rmx
mm(out-rmx)
 */
int uncompressFloatAllWorker(UfwInfo *ufwInfo){//unsigned int *data, unsigned int *odata, int bits, long n, long off){
  //run as a worker thread...
  unsigned int mask,emask,word;
  int offset,left;
  long pos,i;
  unsigned int *data=ufwInfo->data;
  unsigned int *odata=ufwInfo->odata;
  int bits;//=ufwInfo->bits;
  int mbits=ufwInfo->bits;
  long n=ufwInfo->n;
  int ebits;
  int min=ufwInfo->expMin;
  long off=ufwInfo->off;
  offset=ufwInfo->expMax-ufwInfo->expMin+1;
  ebits=1;
  while((offset>>ebits)>0){
    ebits++;
  }

  bits=ufwInfo->bits+ebits+1;
  //for(i=0; i<bits; i++)
  //  mask|=1<<(31-i);
  //mask=((1<<bits)-1)<<(32-bits);
  mask=((1<<mbits)-1)<<(23-mbits);
  emask=((1<<ebits)-1);
  for(i=0; i<n; i++){
    pos=((i+off)*bits)/32;//word to read from
    offset=((i+off)*bits)%32;//offset in word to read from.
    left=bits-(32-offset);//no of bits needed from next word.
    word=(data[pos]<<offset);
    if(left>0)
      word|=data[pos+1]>>(bits-left);
    if(((word>>(31-ebits))&emask)==0)
      odata[i]=(word&0x80000000)|((word>>(8-ebits))&mask);
    else
      odata[i]=(word&0x80000000)|((((word>>(31-ebits))&emask)+min-1)<<(23))|((word>>(8-ebits))&mask);
  }
  return 0;
}

static PyObject *UncompressFloatArrayAllThreaded(PyObject *self,PyObject *args){
  long off,i;//,pos;
  int bits;//,offset,left;
  long ndone,nleft;
  //unsigned int mask=0;
  PyArrayObject *arr,*out;
  //unsigned int *data,*odata;
  int nthreads;
  UfwInfo *ufwInfo;
  pthread_t *thread;
  int expMin,expMax;
  if(!PyArg_ParseTuple(args,"O!O!iiili",&PyArray_Type,&arr,&PyArray_Type,&out,&bits,&expMin,&expMax,&off,&nthreads)){
    printf("Usage: numpy compressed array,numpy output array number of bits to compress each element too, offset index in compressed array, nthreads.\n");
    return NULL;
  }
  if((thread=malloc(sizeof(pthread_t)*nthreads))==NULL){
    printf("cmod.utils unable to malloc thread\n");
    return NULL;
  }
  if((ufwInfo=malloc(sizeof(UfwInfo)*nthreads))==NULL){
    printf("cmod.utils unable to malloc ufwInfo\n");
    return NULL;
  }
  if(!PyArray_ISCONTIGUOUS(arr)){
    printf("Input array must be contiguous\n");
    return NULL;
  }
  if(arr->nd!=1){
    printf("Input array must be 1D\n");
    return NULL;
  }
  if(arr->descr->elsize!=sizeof(int)){
    printf("Input array must be uint32 (or float, or anything 4 bytes)\n");
    return NULL;
  }
  if(!PyArray_ISCONTIGUOUS(out)){
    printf("Output array must be contiguous\n");
    return NULL;
  }
  if(out->nd!=1){
    printf("Output array must be 1D\n");
    return NULL;
  }
  if(out->descr->type_num!=NPY_FLOAT){
    printf("Output array must be float32\n");
    return NULL;
  }
  if(expMin<1 || expMax>255 || expMin>expMax){
    printf("ExpMin should be >0, expMax<256.\n");
    return NULL;
  }
  //for(i=0; i<bits; i++)
  //  mask|=1<<(31-i);
  nleft=out->dimensions[0];
  ndone=0;
  for(i=0; i<nthreads; i++){
    ufwInfo[i].data=(unsigned int*)arr->data;
    ufwInfo[i].bits=bits;
    ufwInfo[i].expMax=expMax;
    ufwInfo[i].expMin=expMin;
    ufwInfo[i].n=nleft/(nthreads-i);
    ufwInfo[i].off=off+ndone;
    ufwInfo[i].odata=&(((unsigned int*)out->data)[ndone]);

    nleft-=ufwInfo[i].n;
    ndone+=ufwInfo[i].n;
    pthread_create(&thread[i],NULL,(void*)uncompressFloatAllWorker,&ufwInfo[i]);


  }
  for(i=0; i<nthreads; i++){
    pthread_join(thread[i],NULL);
  }

  free(ufwInfo);
  free(thread);
  Py_INCREF(Py_None);
  return Py_None;
}



static PyObject* ArrayFromArray(PyObject *self,PyObject *args){
    //create an array from existing memory, and return new object.  The new array can use all or part of the old one, and can be different data type.
    char *type;
    PyObject* dim;
    npy_intp dims[5];
    int offset=0,i;
    int size,oldsize,itype;
    PyArrayObject *narr;
    PyArrayObject *newarr;
    PyObject *narray;
    dims[0]=dims[1]=dims[2]=dims[3]=dims[4]=-1;
    if(!PyArg_ParseTuple(args,"OOs|i",&narray,&dim,&type,&offset)){
	printf("Usage: Numeric array, new dimensions, new type, offset (bytes) into array to start new array(optional)\nCreates a new array from the memory used by the old array, sharing this memory.\n");
	return NULL;
    }
    narr=(PyArrayObject*)narray;
    if (!PyArg_ParseTuple(dim,"i|iiii",&dims[0],&dims[1],&dims[2],&dims[3],&dims[4])){
	printf("Usage2: numeric array,new dims(sizex,y,z,zz,zzz),new type\n");
	return NULL;
    }
    if(!PyArray_ISCONTIGUOUS(narr)){
	printf("Input array must be contiguous\n");
	return NULL;
    }
    size=dims[0]*(dims[1]>0?dims[1]:1)*(dims[2]>0?dims[2]:1)*(dims[3]>0?dims[3]:1)*(dims[4]>0?dims[4]:1);
  if(*type=='f'){
    size*=sizeof(float);
    itype=NPY_FLOAT;
  }else if(*type=='d'){
    size*=sizeof(double);
    itype=NPY_DOUBLE;
  }else if(*type=='D'){
    size*=sizeof(double)*2;
    itype=NPY_CDOUBLE;
  }else if(*type=='F'){
    size*=sizeof(float)*2;
    itype=NPY_CFLOAT;
  }else if(*type=='i'){
    size*=sizeof(int);
    itype=NPY_INT;
  }else if(*type=='s'){
    size*=sizeof(short);
    itype=NPY_SHORT;
  }else if(*type=='c'){
    size*=sizeof(char);
    itype=NPY_CHAR;
  }else if(*type=='u'){//changed from char to int on 060502.
    size*=sizeof(int);
    itype=NPY_UINT;
  }else if(*type=='b'){
    size*=sizeof(char);
    itype=NPY_UBYTE;
/*  }else if(*type=='1'){
    size*=sizeof(char);
    itype=NPY_SBYTE;*/
  }else if(*type=='l'){
    size*=sizeof(long);
    itype=NPY_LONG;
  }else{
    printf("cmod.utils.arrayFromArray: Type code %c not recognised - assuming int\n",*type);
    itype=NPY_INT;
    size*=sizeof(int);
  }  

  oldsize=1;
  for(i=0; i<narr->nd; i++){
      oldsize*=narr->dimensions[i];
  }
  switch(narr->descr->type_num){
      case NPY_CHAR:
      case NPY_UBYTE:
//      case NPY_SBYTE:
	  oldsize*=sizeof(char);
	  break;  
      case NPY_SHORT:
	  oldsize*=sizeof(short);
	  break;
      case NPY_INT:
	  oldsize*=sizeof(int);
	  break;
      case NPY_UINT:
	  oldsize*=sizeof(int);
	  break;
      case NPY_LONG:
	  oldsize*=sizeof(long);
	  break;
      case NPY_FLOAT:
	  oldsize*=sizeof(float);
	  break;
      case NPY_DOUBLE:
	  oldsize*=sizeof(double);
	  break;
      case NPY_CFLOAT:
	  oldsize*=sizeof(float)*2;
	  break;
      case NPY_CDOUBLE:
	  oldsize*=sizeof(double)*2;
	  break;
      case NPY_OBJECT:
	  oldsize*=sizeof(char*);
	  break;
      default:
	  printf("Datatype of old array not known, assuming int (%d)\n",narr->descr->type_num);
	  oldsize*=sizeof(int);
	  break;
  }
  if(size+offset>oldsize){
      printf("New array won't fit in old array %s so you should use a larger original array.\n",offset==0?"":"with specified offset");
      printf("%d\n",narr->nd);
      for(i=0; i<narr->nd; i++){
	  printf("%d\n",(int)(narr->dimensions[i]));
      }
      printf("%d %d %d %d\n",size+offset,oldsize,size,offset);
      for(i=0; i<5; i++){
	  printf("%d\n",(int)(dims[i]));
      }
      return NULL;
  }
  //newarr=(PyArrayObject*)PyArray_FromDimsAndData(dims[1]>0?(dims[2]>0?(dims[3]>0?(dims[4]>0?5:4):3):2):1,dims,itype,narr->data+offset);//this data probably shouldn't ever be freed?
  newarr=(PyArrayObject*)PyArray_SimpleNewFromData(dims[1]>0?(dims[2]>0?(dims[3]>0?(dims[4]>0?5:4):3):2):1,dims,itype,narr->data+offset);//this data probably shouldn't ever be freed?
  newarr->base=(PyObject*)narr;//so that decref can be called when newarr finished with.
  Py_INCREF((PyObject*)narr);//since it shares the data...
  return (PyObject*)newarr;
}

static PyObject* ArrayFromDiagonal(PyObject *self,PyObject *args){
    //create a vector from existing memory, the diagonal, and return new object.  The new array uses the diag of the old one.
/*Test with:
import cmod.utils,Numeric
a=Numeric.zeros((10,10))
d=cmod.utils.arrayFromDiagonal(a)
*/
  return NULL;//can be done using: a.ravel()[::a.shape[0]+1]
/*
    int dims[1];
    PyArrayObject *narr;
    PyArrayObject *newarr;
    PyObject *po;
    dims[0]=dims[1]=dims[2]=dims[3]=-1;
    if(!PyArg_ParseTuple(args,"O!",&PyArray_Type,&narr)){
	printf("Usage: Numeric array\nCreates a new array from the memory used by the diagonal of the old array, sharing this memory.\n");
	return NULL;
    }
    printf("this will segment don't use - instead use a.flat[::11]=xxx\n");
    if(narr->nd!=2){
	printf("Array must be 2D\n");
	return NULL;
    }
    printf("here2 %d %d\n",narr->strides[0],narr->strides[1]);
    dims[0]=narr->dimensions[0]>narr->dimensions[1]?narr->dimensions[1]:narr->dimensions[0];

    newarr=(PyArrayObject*)PyArray_FromDimsAndData(1,dims,narr->descr->type_num,narr->data);//this data probably shouldn't ever be freed?
    printf("here3 %d\n",newarr->strides[0]);
    newarr->strides[0]=narr->strides[0]+narr->strides[1];
    printf("here4\n");
    //newarr->base=(PyObject*)narr;//so that decref can be called when newarr finished with.
    printf("here5\n");
    //Py_INCREF((PyObject*)narr);//since it shares the data...
    printf("here6 %d %d %d %d %d\n",newarr->nd,newarr->dimensions[0],newarr->strides[0],narr->strides[0],narr->strides[1]);
    po=(PyObject*)newarr;
    printf("here7\n");
    return po;
*/
}

static PyObject *mmapArray(PyObject *self,PyObject *args){
  //create a Numeric array from mmap memory... this can allow huge arrays to be created (though access may be slow).
  npy_intp *dims=NULL;
  int ndim;
  struct stat64 st;
  size_t size=1;
  char *filename;
  void *data;
  PyArrayObject *newarr;
  PyObject *fileobj,*dimtuple;
  char *type;
  int itype,i,fd;
  if(!PyArg_ParseTuple(args,"OOs",&fileobj,&dimtuple,&type)){
    printf("Usage: Filename or None, dimensions, type\n");
    return NULL;
  }
  if(fileobj==Py_None){
    //make up a filename
    filename=tmpnam(NULL);
    if(filename)
      printf("Using temporary file %s for mmap (you should delete this when finished with)\n",filename);
    else{
      printf("Unable to get filename for mmap.\n");
      return NULL;
    }
  }else if(!(filename=PyString_AsString(fileobj))){
    printf("filename should be None or a string\n");
    return NULL;
  }
  if(PyTuple_Check(dimtuple)){
    ndim=PyTuple_Size(dimtuple);
    if(ndim<=0){
      printf("Dimensions invalid\n");
      return NULL;
    }
    dims=(npy_intp*)malloc(sizeof(npy_intp)*ndim);
    if(dims==NULL){
      printf("Unable to allocate dims array\n");
      return NULL;
    }
    for(i=0; i<ndim; i++){
      dims[i]=PyInt_AsLong(PyTuple_GetItem(dimtuple,i));
      if(dims[i]<=0){
	printf("Dims should be >0\n");
	free(dims);
	return NULL;
      }
      size*=dims[i];
    }
  }else{
    printf("Dims must be a tuple\n");
    return NULL;
  }

  switch(*type){
    case 'f':
      size*=sizeof(float);
      itype=NPY_FLOAT;
      break;
    case 'd':
      size*=sizeof(double);
      itype=NPY_DOUBLE;
      break;
    case 'D':
      size*=sizeof(double)*2;
      itype=NPY_CDOUBLE;
      break;
    case 'F':
      size*=sizeof(float)*2;
      itype=NPY_CFLOAT;
      break;
    case 'i':
      size*=sizeof(int);
      itype=NPY_INT;
      break;
    case 's':
      size*=sizeof(short);
      itype=NPY_SHORT;
      break;
    case 'c':
      size*=sizeof(char);
      itype=NPY_CHAR;
      break;
    case 'u'://changed from char to int on 060502.
      size*=sizeof(int);
      itype=NPY_UINT;
      break;
    case 'b':
      size*=sizeof(char);
      itype=NPY_UBYTE;
      break;
/*    case '1':
      size*=sizeof(char);
      itype=NPY_SBYTE;
      break;*/
    case 'l':
      size*=sizeof(long);
      itype=NPY_LONG;
      break;
    default:
      printf("cmod.utils.mmapArray: Type code %c not recognised - assuming int\n",*type);
      itype=NPY_INT;
      size*=sizeof(int);
      break;
  }
  //see if the file exists...
  if(stat64(filename,&st)==0){//file exists...
    //see if the file is big enough
    if(st.st_size<size){
      //append to file...
      printf("File not large enough for mmap - appending to it\n");
      if(truncate64(filename,size)<0){
	printf("Call (1) to truncate64 failed\n");
	perror(NULL);
      }
    }
  }else{//file doesn't exist..
    //create file and make it right size...
    printf("Creating file %s for mmap of size %ld bytes\n",filename,size);
    if((fd=open(filename,O_RDWR|O_CREAT,S_IRWXU|S_IRWXG|S_IRWXO))>=0){
      close(fd);
      if(truncate64(filename,size)<0){
	printf("Call (2) to truncate64 failed\n");
	perror(NULL);
      }
    }else{
      printf("Could not open file %s for creating and truncation\n",filename);
    }
  }
  if((fd=open64(filename,O_RDWR))>=0){
    data=mmap(0,size,PROT_READ|PROT_WRITE,MAP_SHARED,fd,0);
    close(fd);
  }else{
    data=MAP_FAILED;
    printf("Could not open file %s for mapping\n",filename);
  }
  if(data==MAP_FAILED){
    printf("mmap call failed (errno=%d)\n",errno);
    printf("Errno: %d %d %d %d %d %d %d\n",EBADF ,EACCES,EINVAL,ETXTBSY,EAGAIN ,ENOMEM,ENODEV);
    perror("Error was:\n");
    return NULL;
  }
  //newarr=(PyArrayObject*)PyArray_FromDimsAndData(ndim,dims,itype,data);//this data probably shouldn't ever be freed?
  newarr=(PyArrayObject*)PyArray_SimpleNewFromData(ndim,dims,itype,data);//this data probably shouldn't ever be freed?
  return PyArray_Return(newarr);
}

static PyObject *mmapArrayFree(PyObject *self,PyObject *args){//free a mmap Array.
  size_t len=1;
  int i;
  PyArrayObject *arr;
  if(!PyArg_ParseTuple(args,"O!",&PyArray_Type,&arr)){
    printf("Usage: Numeric array created by mmap\n");
    return NULL;
  }
  for(i=0; i<arr->nd; i++){
    len*=arr->dimensions[i];
  }
  if(munmap(arr->data,len)==-1){
    printf("cmod.utils: munmap call failed\n");
    return NULL;
  }
  //now, since we've removed the memory, we need to make sure the array remains safe.  So, create a single element array...
  switch(arr->descr->type_num){
      case NPY_CHAR:
      case NPY_UBYTE:
//      case NPY_SBYTE:
	  len=sizeof(char);
	  break;  
      case NPY_SHORT:
	  len=sizeof(short);
	  break;
      case NPY_INT:
	  len=sizeof(int);
	  break;
      case NPY_UINT:
	  len=sizeof(int);
	  break;
      case NPY_LONG:
	  len=sizeof(long);
	  break;
      case NPY_FLOAT:
	  len=sizeof(float);
	  break;
      case NPY_DOUBLE:
	  len=sizeof(double);
	  break;
      case NPY_CFLOAT:
	  len=sizeof(float)*2;
	  break;
      case NPY_CDOUBLE:
	  len=sizeof(double)*2;
	  break;
      case NPY_OBJECT:
	  len=sizeof(char*);
	  break;
      default:
	  printf("Datatype of old array not known, assuming int (%d)\n",arr->descr->type_num);
	  len=sizeof(int);
	  break;
  }
  for(i=0; i<arr->nd; i++){
    arr->dimensions[i]=1;
    arr->strides[i]=len;
  }
  arr->data=malloc(len);//this will never be freed, but since its only 1 byte, do we care?
  //Note, this isn't quite safe, the array probably shouldn't be accessed, since other parts of arr->descr have not been updated for the new typecode!
  Py_INCREF(Py_None);
  return Py_None;
}

void calcCentroid(PyArrayObject *narr,float* xcent,float* ycent){
    float tot;
    int i,j;
    int isint;
    isint=(narr->descr->kind=='i' && narr->descr->elsize!=sizeof(int));
    *xcent=*ycent=tot=0.;
    for(i=0; i<narr->dimensions[0]; i++){
	for(j=0; j<narr->dimensions[1]; j++){
	    if(isint){
		*xcent+=(float)(*(int*)(narr->data+i*narr->strides[0]+j*narr->strides[1]))*i;
		*ycent+=(float)(*(int*)(narr->data+i*narr->strides[0]+j*narr->strides[1]))*j;
		tot+=(float)(*(int*)(narr->data+i*narr->strides[0]+j*narr->strides[1]));
	    }else{
		*xcent+=(*(float*)(narr->data+i*narr->strides[0]+j*narr->strides[1]))*i;
		*ycent+=(*(float*)(narr->data+i*narr->strides[0]+j*narr->strides[1]))*j;
		tot+=(*(float*)(narr->data+i*narr->strides[0]+j*narr->strides[1]));
	    }
	}
    }
    if(tot==0.)
	tot=1.;
    *xcent/=tot;
    *ycent/=tot;
}

static PyObject* Centroid(PyObject *self,PyObject *args){
    //compute centroid values.
    PyObject *narray;
    PyArrayObject *narr;
    float xcent=0.,ycent=0.;
    if(!PyArg_ParseTuple(args,"O",&narray)){
	printf("Usage: Numeric array\n");
	return NULL;
    }
    narr=(PyArrayObject*)narray;
    if(narr->nd!=2){
	printf("Array must be 2D\n");
	return NULL;
    }
    if(narr->descr->type_num!=NPY_FLOAT && (narr->descr->kind!='i' || narr->descr->elsize!=sizeof(int))){
	printf("centroid: Array must be Float32 or INT32\n");
	return NULL;
    }
    calcCentroid(narr,&xcent,&ycent);
    return Py_BuildValue("ff",xcent,ycent);
}


static PyObject* RadialProfileEncircledEnergy(PyObject *self,PyObject *args){
    //compute the radial profile and encircled energy of an image.  If x,y values are given, will use these as the centre of the image.  If not, will use the maximum value as the centre of the image.
    npy_intp dims[2];
    int i,j;
    int xpos=-1,ypos=-1;
    int xmax,ymax,maxdist;
    float max,pval;
    PyArrayObject *narr;
    PyArrayObject *outarr;
    PyObject *narray;
    float *binarr;
    int *countarr;
    int npos;
    int indx;
    int maxtype=0;
    int isint;
    float halfMax;
    float fwhm,a,b,d50;
    if(!PyArg_ParseTuple(args,"O|iii",&narray,&xpos,&ypos,&maxtype)){
	printf("Usage: Numeric array, x and y centre position (optional), centre location type (0==max pxl value, 1==centroid)\n");
	return NULL;
    }
    narr=(PyArrayObject*)narray;
    if(narr->nd!=2){
	printf("radialProfileEncircledEnergy: Array must be 2D\n");
	return NULL;
    }
    if((narr->descr->type_num!=NPY_FLOAT) && (narr->descr->kind!='i' || narr->descr->elsize!=sizeof(int))){
	printf("radialProfileEncircledEnergy: Array must be Float32 or Int32\n");
	return NULL;
    }
    isint=(narr->descr->kind!='i' && narr->descr->elsize!=sizeof(int));
    if(xpos==-1 || ypos==-1){
	//first find maximum light position
	if(maxtype==0){//use max pxl value as centre.
	    max=0;
	    for(i=0; i<narr->dimensions[0]; i++){
		for(j=0; j<narr->dimensions[1]; j++){
		    if(isint)
			pval=(float)(*(int*)(narr->data+i*narr->strides[0]+j*narr->strides[1]));
		    else
			pval=(*(float*)(narr->data+i*narr->strides[0]+j*narr->strides[1]));
		    if(pval>max){
			max=pval;
			xpos=i;
			ypos=j;
		    }
		}
	    }
	}else{//use centroid location as centre...
	    float xcent,ycent;
	    calcCentroid(narr,&xcent,&ycent);
	    xpos=(int)(xcent+0.5);
	    ypos=(int)(ycent+0.5);
	}
	//printf("cmod.utils: Maximum found at (%d, %d)\n",xpos,ypos);
    }
    //now compute the maximum square distance from this point.
    xmax=xpos<narr->dimensions[0]-xpos-1 ? narr->dimensions[0]-xpos-1 : xpos;
    ymax=ypos<narr->dimensions[1]-ypos-1 ? narr->dimensions[1]-ypos-1 : ypos;
    maxdist=xmax*xmax+ymax*ymax;
    //printf("cmod.utils: Maximum distance is %d\n",maxdist);
    //create an array to hold the values.
    binarr=(float*)malloc(sizeof(float)*(maxdist+1));
    countarr=(int*)malloc(sizeof(int)*(maxdist+1));
    memset(binarr,0,sizeof(float)*(maxdist+1));
    memset(countarr,0,sizeof(float)*(maxdist+1));
    npos=0;//number of bins...
    for(i=0; i<narr->dimensions[0]; i++){
	for(j=0; j<narr->dimensions[1]; j++){
	    indx=(i-xpos)*(i-xpos)+(j-ypos)*(j-ypos);
	    if(isint)
		binarr[indx]+=(float)(*(int*)(narr->data+i*narr->strides[0]+j*narr->strides[1]));
	    else
		binarr[indx]+=(*(float*)(narr->data+i*narr->strides[0]+j*narr->strides[1]));
	    if(countarr[indx]==0)
		npos++;
	    countarr[indx]++;
	}
    }
    //printf("Have %d bins\n",npos);
    //now create Numeric array to hold the values...
    //for(i=0; i<maxdist+1; i++){
    //printf("%d %g\n",countarr[i],binarr[i]);
    //}
    dims[0]=3;
    dims[1]=npos;
    //outarr=(PyArrayObject *)PyArray_FromDims(2,dims,NPY_FLOAT);
    outarr=(PyArrayObject *)PyArray_SimpleNew(2,dims,NPY_FLOAT);
    indx=0;
    //put the values into the output array.
    for(i=0; i<=maxdist; i++){
	if(countarr[i]>0){
	    //store the profile... (mean)
	    (*(float*)(outarr->data+outarr->strides[0]*0+outarr->strides[1]*indx))=binarr[i]/countarr[i];
	    //store the enclosed energy (cumulative)
	    if(indx==0){
		(*(float*)(outarr->data+outarr->strides[0]*1+outarr->strides[1]*indx))=binarr[i];
	    }else{
		(*(float*)(outarr->data+outarr->strides[0]*1+outarr->strides[1]*indx))=binarr[i]+(*(float*)(outarr->data+outarr->strides[0]*1+outarr->strides[1]*(indx-1)));
	    }
	    //and store the distance from centre...
	    (*(float*)(outarr->data+outarr->strides[0]*2+outarr->strides[1]*indx))=sqrt((double)i);
	    
	    indx++;
	}
    }
    //printf("Values placed into output array (indx,npos=%d %d, same?)\n",indx,npos);
    //for(i=0; i<indx; i++){
//	printf("%g %g %g\n",(*(float*)(outarr->data+outarr->strides[0]*2+outarr->strides[1]*i)),(*(float*)(outarr->data+outarr->strides[0]*0+outarr->strides[1]*i)),(*(float*)(outarr->data+outarr->strides[0]*1+outarr->strides[1]*i)));
//    }
    free(binarr);
    free(countarr);
    //Now compute the fwhm...
    //Find the first point where the profile drops to half of the max.
    //Note, this is in pixel scale, so should be multiplied by a scaling factor for correct result.
    i=1;
    halfMax=(*(float*)(outarr->data))/2;
    while(i<npos){
	if((*(float*)(outarr->data+outarr->strides[1]*i))<halfMax)
	    break;
	i++;
    }
    if(i==npos)
	fwhm=1./0.;
    else{
	//do a linear interpolation to get the fwhm.
	//a=(y2-y1)/(x2-x1).
	a=((*(float*)(outarr->data+outarr->strides[1]*i))-(*(float*)(outarr->data+outarr->strides[1]*(i-1))))/
	    ((*(float*)(outarr->data+outarr->strides[0]*2+outarr->strides[1]*i))-(*(float*)(outarr->data+outarr->strides[0]*2+outarr->strides[1]*(i-1))));
	//b=y2-a*x2
	b=(*(float*)(outarr->data+outarr->strides[1]*i))-a*(*(float*)(outarr->data+outarr->strides[0]*2+outarr->strides[1]*i));
	fwhm=2*(halfMax-b)/a;
    }
    //printf("Got fwhm %g\n",fwhm);
    //now compute the diameter of the aperture with 50% of the energy.
    halfMax=(*(float*)(outarr->data+outarr->strides[0]+outarr->strides[1]*(npos-1)))/2;//half the total enclosed energy...
    i=0;
    while(i<npos){
	if((*(float*)(outarr->data+outarr->strides[0]+outarr->strides[1]*i))>halfMax)
	    break;
	i++;
    }
    //do linear interpolation to get d50...
    //a=(y2-y1)/(x2-x1).
    a=((*(float*)(outarr->data+outarr->strides[0]+outarr->strides[1]*i))-(*(float*)(outarr->data+outarr->strides[0]+outarr->strides[1]*(i-1))))/
	((*(float*)(outarr->data+outarr->strides[0]*2+outarr->strides[1]*i))-(*(float*)(outarr->data+outarr->strides[0]*2+outarr->strides[1]*(i-1))));
    //b=y2-a*x2
    b=(*(float*)(outarr->data+outarr->strides[0]+outarr->strides[1]*i))-a*(*(float*)(outarr->data+outarr->strides[0]*2+outarr->strides[1]*i));
    d50=(halfMax-b)/a;
    //printf("Got d50 %g\n",d50);
	
    return Py_BuildValue("Off",(PyObject*)outarr,fwhm,d50);

}

static PyObject *queryNumericArray(PyObject *self,PyObject *args){
  PyObject *Narray,*tup;
  PyArrayObject *NumericArray;
  //int i;
  int raiseError=1;
  if(!PyArg_ParseTuple(args, "O|i", &Narray,&raiseError)){
    printf("Usage: NumericArray (optional, raiseError=1 if non-contiguous)\n");
    return NULL;
  }
  NumericArray=(PyArrayObject*)Narray;
  if(raiseError==1 && (NumericArray->flags&1)==0){
    printf("Error - Array is non-contiguous\n");
    return NULL;
  }
  switch(NumericArray->nd){
      case 1:
	  tup=Py_BuildValue("(i)",NumericArray->strides[0]);
	  break;
      case 2:
	  tup=Py_BuildValue("(ii)",NumericArray->strides[0],NumericArray->strides[1]);
	  break;
      case 3:
	  tup=Py_BuildValue("(iii)",NumericArray->strides[0],NumericArray->strides[1],NumericArray->strides[2]);
	  break;
      case 4:
	  tup=Py_BuildValue("(iiii)",NumericArray->strides[0],NumericArray->strides[1],NumericArray->strides[2],NumericArray->strides[3]);
	  break;
      default:
	  tup=Py_BuildValue("()");
  }
  return Py_BuildValue("(lO)",(long)NumericArray->data,tup);
}

static PyObject *
arrayfrombuffer(PyObject *self, PyObject *args)
{
  return NULL;//probably no longer needed.
  /*
  PyArrayObject *arrayobj;
  PyObject *bufobj;
  char *type = "l";  // long int is default 
  void *data;
  int datalen;
  npy_intp nitems[1];
  PyArray_Descr *descr;

  if (!PyArg_ParseTuple(args, "O|s:arrayfrombuffer", &bufobj, &type)) {
    return NULL;
  }
  if (-1 == PyObject_AsWriteBuffer(bufobj, &data, &datalen)) {
    return NULL;
  }
  if (!(descr = PyArray_DescrFromType(*type))) {
    return NULL;
  }
  if (datalen % descr->elsize != 0) {
    PyErr_SetString(PyExc_ValueError, "buffer length must be a multiple of element size");
    return NULL;
  }
  nitems[0] = datalen / descr->elsize;
  // The Numeric docs say not to do this unless 'data' is eternal,
  // but they're wrong.  See longer comment below. 
  //if (!(arrayobj = PyArray_FromDimsAndData(1,&nitems,*type,data))) {
  if (!(arrayobj = PyArray_SimpleNewFromData(1,nitems,*type,data))) {
    return NULL;
  }
  // This prevents the bufobj from being inadvertently
  // deallocated before the array is, and arranges for the bufobj
  // to be deallocated when the array is if it's not otherwise
  // referenced.  Of course, you can still cause a core dump or
  // incorrect results by invalidating the bufobj's memory before
  // deleting the array --- for example, by calling close() on an
  // mmap object.  Them's the breaks.
  arrayobj->base = bufobj;
  Py_INCREF(bufobj);
  // do I need to incref arrayobj?! fromstring does... 
  return (PyObject*)arrayobj;*/
}

//
// Multi-threaded matrix-vector multiplication.
// Urban Bitenc, July 2nd 2012
//
// Possible usages:
// c = dot(a, b)
// c = dot(a, b, r)
// c = dot(a, b, N)
// c = dot(a, b, r, N)
// 
// a ... 2-D array - input matrix
// b ... 1-D array - input vector
// r ... 1-D array - result vector (if you provide it on input, the function does 
//                   not need to allocate memory; increase in speed of about 0.3%)
// N ... the number of threads:
//           N <= 0 : sysconf(_SC_NPROCESSORS_ONLN)/2 threads for big matrices,
//                    1 thread for small matrices (m*n < 500.000)
//           N >= 1 : N number of threads are used
//
//
static PyObject* dot(PyObject* self, PyObject* args)
{
  PyArrayObject *ArrayM; // First argument:  input 2-D array
  PyArrayObject *ArrayV; // Second argument: input 1-D array
  PyObject      *arg3;   // Third argument (optional). Can be the output array (PyArrayObject)
                              // or the number of threads (integer) or None. If it is the
                              // number of threads (int), the value in the fourth argument is ignored.
  int nThreads = 0;      // Fourth argument (optional): the number of threads (int).
                         // The default is 0, which means default multi-threaded.

  PyArrayObject *ArrayR; // The Results array.
  
  void* M;               // input Matrix data (from ArrayM)  - they are void so that
  void* V;               // input Vector data (from ArrayV)  - they can be used for
  void* R;               // output Vector data (to ArrayR)   - float and double
  int    n;              // size of V
  int    m;              // size of R
  int bytesPerElt;       // to store the number of bytes (4 or 8) per each element. This is necessary to
                         // make the same code work for both float and double.
  int inputFloat;        // 1 if input is float, 0 if input is double, undefined if it is neither.

  npy_intp* iSizes;      // sizes of the input array
  int i;                 // for-loop index

  int forEach;           // m/nThreads = number of lines to be processed by each thread
  int oneMore;           // m%nThreads = number of threads that have to process one line more

  int check_arg3 = 0;    // a flag to check the output array to make the code more readable

  pthread_t*  thread;    // for the vector of threads
  mvm_data_t* params;    // for the vector of parameters (mvm_data_t defined in mvm.h)


  // (1) PARSE THE INPUT arguments from python:
  if( !PyArg_ParseTuple(args, "O!O!|Oi",
  			&PyArray_Type, &ArrayM, // input 2-D matrix
  			&PyArray_Type, &ArrayV, // input 1-D vector
			&arg3,                  // optional: output 1-D vector OR number of threads
  			&nThreads) )            // optional: number of threads
    return NULL;

  // (2) CHECK THE INPUT:
  // check the first two arrays:
  if(ArrayM->nd != 2) // input  matrix must be 2-D
    {
      PyErr_SetString(PyExc_TypeError, "1st argument must be 2-D array.");
      return NULL;
    }
  if(ArrayV->nd != 1) // input  vector must be 1-D
    {
      PyErr_SetString(PyExc_TypeError, "2nd argument must be 1-D array.");
      return NULL;
    }
  if(!PyArray_ISCONTIGUOUS( ArrayM ) )
    {
      PyErr_SetString(PyExc_TypeError, "1st argument must be contiguous.");
      return NULL;
    }
  if(!PyArray_ISCONTIGUOUS( ArrayV ) )
    {
      PyErr_SetString(PyExc_TypeError, "2nd argument must be contiguous.");
      return NULL;
    }
  if( !(PyArray_TYPE( ArrayM )==NPY_FLOAT || PyArray_TYPE( ArrayM )==NPY_DOUBLE) )
    {
      PyErr_SetString(PyExc_TypeError, "1st argument must be an array of FLOATs or DOUBLEs.");
      return NULL;
    }
  if( PyArray_TYPE( ArrayM ) != PyArray_TYPE( ArrayV ) )
    {
      PyErr_SetString(PyExc_TypeError, "1st and 2nd argument must be of the same type (both float or both double).");
      return NULL;
    }
  if(ArrayM->dimensions[1] != ArrayV->dimensions[0]) // matrix and input vector sizes must agree
    {
      printf("N cols should equal length of second, but %ld != %ld\n",ArrayM->dimensions[1],ArrayV->dimensions[0]);
      PyErr_SetString(PyExc_ValueError, "Number of columns of 1st argument must agree with length of 2nd argument.");
      return NULL;
    }

  // Check, whether input is float or double and save this information.
  // If it is neither float or double, exit. If this happens, something is wrong with this code,
  // because if it is neither float nor double, it should have exited earlier already.
  if( PyArray_TYPE( ArrayM )==NPY_FLOAT  )
    {
      bytesPerElt = sizeof(float);
      inputFloat = 1;
    }
  else if( PyArray_TYPE( ArrayM )==NPY_DOUBLE )
    {
      bytesPerElt = sizeof(double);
      inputFloat = 0;
    }
  else{
    PyErr_SetString(PyExc_TypeError, "Input array must be float od double. At this point in the program it really should be either one or the other; this is a strange bug.");
    return NULL;
  }

  // Check the 3rd argument and set the outputArray accordingly:
  if( arg3->ob_type == &PyArray_Type ) // If it is an array
    {                  // Set the flag to check the output array (or generate it if not given)
      check_arg3 = 1;  // later. (You could check it straight away, but the code would be less readable.)
    }
  else if( PyInt_Check( arg3 ) ) // If it is an integer
    {
      // Take its value as the number of threads (and herewith ignore the value in argument 4):
      nThreads = PyInt_AsLong( arg3 );
    }
  else if( arg3 == Py_None ){} // If it is None, do nothing
  else // If it is not array, integer or none, raise exception:
    {
      PyErr_SetString(PyExc_TypeError, "3rd argument must be PyArrayObject or integer or None.");
      return NULL;
    }

  // Check the array for results, provided on input:
  if( check_arg3 )
    {
      if( ((PyArrayObject *)arg3)->nd != 1) // result vector must be 1-D
	{
	  PyErr_SetString(PyExc_TypeError, "3rd argument array must be 1-D.");
	  return NULL;
	}
      if(!PyArray_ISCONTIGUOUS( arg3 ) )
	{
	  PyErr_SetString(PyExc_TypeError, "3rd argument array must be contiguous.");
	  return NULL;
	}
      if( PyArray_TYPE( ArrayM ) != PyArray_TYPE( arg3 ) )
	{
	  PyErr_SetString(PyExc_TypeError, "The 3rd argument is an array and it must be of the same type (float or double) as the 1st argument.");
	  return NULL;
	}
      if(ArrayM->dimensions[0] != ((PyArrayObject *)arg3)->dimensions[0]) // matrix and output vector sizes must agree
	{
	  PyErr_SetString(PyExc_ValueError, "Number of lines of 1st argument must agree with length of 3rd argument array.");
	  return NULL;
	}
      // If all is fine, arg3 becomes the output array:
      ArrayR = (PyArrayObject *)arg3;
    }
  else // if arg 3 is not an array, you have to generate the output array:
    {
      // Generate the results array:
      // (At this point input really must be either float or double, this has been checked twice already.)
      if(inputFloat) ArrayR = (PyArrayObject *)PyArray_SimpleNew(1, &(ArrayM->dimensions[0]), NPY_FLOAT);
      else           ArrayR = (PyArrayObject *)PyArray_SimpleNew(1, &(ArrayM->dimensions[0]), NPY_DOUBLE);
    }

  // END OF "CHECKS OF THE INPUT"
  //
  // ArrayR is now the output array, which was either provided on input or
  // created above.
  //

  // (3) FOR EASE OF USE:
  // Extract the size of the input array:
  iSizes = ArrayM->dimensions; // sizes of array
  n = (int) iSizes[1]; // input vector size
  m = (int) iSizes[0]; // output vector size

  // Get the data from the python Array objects into c arrays
  M = ArrayM->data;
  V = ArrayV->data;
  R = ArrayR->data;


  // (4) DETERMINE THE NUMBER OF THREADS:
  // If the input nThreads is < 1, use the default multi-threading.
  // Find out the number of cores (hyperthreads) and use half that number for the number of threads.
  // (On average this seems to make most sense - see U.B.`s log book 1, p. 124.)
  if(nThreads < 1)
    {
      // Check the matrix size. If the matrix is too small, run in a single thread,
      // because this will be faster. The maximum matrix size to run single threaded
      // was determined empirically and is different on different computers - see UB's 
      // log book 1, p.127 (for MacBook und gpu2: 500.000; for cpuserver: 1.100.000, for
      // gpuserver: 340.000). Here we use one value - 500.000:
      if( m*n < 500000)
	{
	  nThreads = 1;
	}
      else
	{
	  // Find the number of cores and divide it by two:
	  nThreads = sysconf(_SC_NPROCESSORS_ONLN)/2;

	  // For the case that _SC_NPROCESSORS_ONLN returns 1, nThreads will at this point be 0.
	  // Or, if anything else funny happens and the value is even negative, correct it to 1
	  // and the process will run in a single thread. Print a warning so the user is aware of that.
	  if(nThreads < 1)
	    {
	      printf("WARNING: Multiple threading was requested for .dot, but the process will run in a single thread, because \"sysconf(_SC_NPROCESSORS_ONLN)\" returned %ld.", sysconf(_SC_NPROCESSORS_ONLN));
	      nThreads = 1;
	    }
	}
    }
  // printf("Number of CPUs: %d\n nThreads: %d\n", (int)sysconf(_SC_NPROCESSORS_ONLN), nThreads);
  //
  // nThreads is now determined (either from input or set to the default value).
  //

  // Important for Python:
  Py_BEGIN_ALLOW_THREADS;

  // (5) THREADING:
  if(nThreads == 1) // for nThreads == 1 you do not even create one thread, but just run it:
    {
      // The non-threaded case:

      // set the parameters
      params = malloc(sizeof(mvm_data_t));
      params->n = n; // input vector size
      params->m = m; // output vector size
      params->M = M;
      params->V = V;
      params->R = R;

      // multiply:
      if( inputFloat ) mvm_float( params );   // At this point input must be either float or double.
      else             mvm_double( params );  //    (This has been checked two times already.)

      // free the memory:
      free(params);
    }
  else
    {
      // Threaded:
      //printf(".dot: generating %d threads.", nThreads);
      // Allocate the memory for the 'params' vector
      // (parameters of the mvm function):
      if((params=malloc(sizeof(mvm_data_t)*nThreads))==NULL)
	{
	  PyErr_SetString(PyExc_MemoryError, "cmod.utils.dot: Failed to malloc 'params'.\n");
	  return NULL;
	}

      // Allocate the memory for 'thread' (the vector of threads):
      if((thread=malloc(sizeof(pthread_t)*nThreads))==NULL)
	{
	  PyErr_SetString(PyExc_MemoryError, "cmod.utils.dot: Failed to malloc 'thread'.\n");
	  free(params); // free the already allocated memory for 'params'
	  return NULL;
	}

      // To determine the number of lines processed by each thread:
      forEach = m/nThreads;
      oneMore = m%nThreads;

      // Run the jobs in threads:
      // The first 'oneMore' threads will process 'forEach+1' matrix lines,
      //     the rest of the threads will process 'forEach' lines.
      for(i=0; i < nThreads; i++)
	{
	  // Set the input parameters:
	  params[i].n = n;         // input vector size
	  params[i].V = V;         // input vector
	  if(i < oneMore) // the first oneMore threads
	      params[i].m = forEach+1; // output vector size
	  else  // the rest of the threads:
	      params[i].m = forEach;
	  if(i == 0)
	    {
	      params[i].M = M;
	      params[i].R = R;
	    }
	  else
	    {
	      params[i].M = params[i-1].M + params[i-1].m*n*bytesPerElt; // first element of the input matrix
	      params[i].R = params[i-1].R + params[i-1].m*bytesPerElt;   // output vector
	    }
	  // Run the job in a thread:
	  if( inputFloat) pthread_create( &thread[i], NULL, (void*)mvm_float,  &params[i] );
	  else            pthread_create( &thread[i], NULL, (void*)mvm_double, &params[i] );
	}

      // Wait untill all the threads are finished:
      for(i=0; i < nThreads; i++)
	{
	  pthread_join(thread[i], NULL);
	}

      // Free the memory allocated for 'thread' and 'params':
      free(params);
      free(thread);
    }
  // END OF THREADING

  // Important for Python:
  Py_END_ALLOW_THREADS;

  // Return the results vector:
  return Py_BuildValue("O", (PyObject*)ArrayR);
}


static PyObject* correlateDifferenceSquared(PyObject *self,PyObject *args){
  PyArrayObject *inarr,*refarr,*outarr=NULL;
  double angle;
  int i,j,x,y;
  int nthreads=0;
  float *idata,*odata;
  int idx,idy,odx,ody;
  rotateStruct rs;
  rotateStructThread *rst;
  pthread_t *thread;
  int ystart;
  int outsizex,outsizey;
  PyObject *outarrobj=NULL;
  npy_intp corrsize[2];
  float *outdata,*indata,*refdata;
  int outstridex,outstridey,instridex,instridey,refstridex,refstridey;
  int inx,iny,ncenx,nceny,mx,my,refx,refy;
  float s,d;
  int dx,dy,minx,miny,sinx,siny,srefx,srefy;

  if(!PyArg_ParseTuple(args,"O!O!|Oi",&PyArray_Type,&inarr,&PyArray_Type,&refarr,&outarrobj,&nthreads)){
    printf("Usage: input array, reference array, output array or correlation size\n");
    return NULL;
  }
  if(inarr->nd!=2){
    printf("Input array must be 2D\n");
    return NULL;
  }
  if(refarr->nd!=2){
    printf("Reference array must be 2D\n");
    return NULL;
  }
  if(inarr->descr->type_num!=NPY_FLOAT){
    printf("Input array must be float32\n");
    return NULL;
  }
  if(refarr->descr->type_num!=NPY_FLOAT){
    printf("Reference array must be float32\n");
    return NULL;
  }
  if(outarrobj==NULL){
    corrsize[0]=inarr->dimensions[0];
    corrsize[1]=inarr->dimensions[1];
    outarr=(PyArrayObject*)PyArray_ZEROS(2,corrsize,NPY_FLOAT,0);
  }else if(PyArray_Check(outarrobj)){//output array specified
    outarr=(PyArrayObject*)outarrobj;
    if(outarr->nd!=2){
      printf("output array must be 2d\n");
      return NULL;
    }
    if(outarr->descr->type_num!=NPY_FLOAT){
      printf("outarr must by float32\n");
      return NULL;
    }
    Py_INCREF(outarr);
    corrsize[0]=outarr->dimensions[0];
    corrsize[1]=outarr->dimensions[1];
  }else if(PyInt_Check(outarrobj)){//size of output specified
    corrsize[0]=corrsize[1]=(int)PyInt_AsLong(outarrobj);
    outarr=(PyArrayObject*)PyArray_ZEROS(2,corrsize,NPY_FLOAT,0);
  }else{
    printf("Output must be an array or int\n");
    return NULL;
  }
  inx=(int)inarr->dimensions[1];
  iny=(int)inarr->dimensions[0];
  nceny=corrsize[0];//the output size.
  ncenx=corrsize[1];
  instridex=inarr->strides[1]/sizeof(float);
  instridey=inarr->strides[0]/sizeof(float);
  refstridey=refarr->strides[0]/sizeof(float);
  refstridex=refarr->strides[1]/sizeof(float);
  outstridey=outarr->strides[0]/sizeof(float);
  outstridex=outarr->strides[1]/sizeof(float);
  refdata=(float*)refarr->data;
  outdata=(float*)outarr->data;
  indata=(float*)inarr->data;
  refy=refarr->dimensions[0];
  refx=refarr->dimensions[1];
  //The ref must be at least as big as the inarr.
  //Now do the correlation
  dx=(refx-inx)/2;
  dy=(refy-iny)/2;
  printf("dx, dy: %d %d\n",dx,dy);
  minx=refx<inx?refx:inx;
  miny=refy<iny?refy:iny;
  for(i=0;i<nceny;i++){
    my=miny;
    if(abs(i-nceny/2)>dy)//the overlap between input and ref.
      my-=abs(i-nceny/2)-dy;
    siny=0;
    if(nceny/2-i>dy)
      siny=nceny/2-i-dy;//starting point for input
    srefy=dy-(i-nceny/2);//starting point for ref.
    if(srefy<0)
      srefy=0;
    for(j=0;j<ncenx;j++){
      mx=minx;//the overlap between input and ref.
      if(abs(j-ncenx/2)>dx)
	mx-=abs(j-ncenx/2)-dx;
      sinx=0;
      if(ncenx/2-j>dx)
	sinx=ncenx/2-j-dx;//starting point for input
      srefx=dx-(j-ncenx/2);
      if(srefx<0)
	srefx=0;
      s=0;
      printf("%d %d: %d %d   %d %d  %d %d\n",i,j,my,mx,srefy,siny,srefx,sinx);
      for(y=0;y<my;y++){
	for(x=0;x<mx;x++){
	  d=refdata[(srefy+y)*refstridey+(srefx+x)*refstridex]-indata[(siny+y)*instridey+(sinx+x)*instridex];
	  s+=d*d;
	}
      }
      outdata[i*outstridey+j*outstridex]=s/(my*mx);
    }
  }
  /*
  for(i=0;i<nceny;i++){
    my=refy-abs(nceny/2-i);
    for(j=0;j<ncenx;j++){
      mx=refx-abs(ncenx/2-j);
      s=0;
      for(y=0;y<my;y++){
	for(x=0;x<mx;x++){
	  if(i<nceny/2){
	    if(j<nceny/2){
	      if(y>=refy || x>=refx || nceny/2-i+y>=iny || ncenx/2-j+x>=inx)
		printf("a  %d %d %d %d  %d %d %d %d\n",y,x,nceny/2-i+y,ncenx/2-j+x,i,j,y,x);
	      d=refdata[y*refstridey+x*refstridex]-indata[(nceny/2-i+y)*instridey+(ncenx/2-j+x)*instridex];
	    }else{
	      if(y>=refy || j-ncenx/2+x>=refx || nceny/2-i+y>=iny || x>=inx)
		printf("b  %d %d %d %d  %d %d %d %d\n",y,j-ncenx/2+x,nceny/2-i+y,x,i,j,y,x);
	      d=refdata[y*refstridey+(j-ncenx/2+x)*refstridex]-indata[(nceny/2-i+y)*instridey+x*instridex];
	    }
	  }else{
	    if(j<ncenx/2){
	      if(i-nceny/2+y>=refy ||x>=refx ||y>=iny ||ncenx/2-j+x>=inx)
		printf("c  %d %d %d %d  %d %d %d %d\n",i-nceny/2+y,x,y,ncenx/2-j+x,i,j,y,x);
	      d=refdata[(i-nceny/2+y)*refstridey+x*refstridex]-indata[(y)*instridey+(ncenx/2-j+x)*instridex];
	    }else{
	      if(i-nceny/2+y>=refy ||j-ncenx/2+x>=refx ||y>=iny||x>=inx)
		printf("d  %d %d %d %d  %d %d %d %d\n",i-nceny/2+y,j-ncenx/2+x,y,x,i,j,y,x);
	      d=refdata[(i-nceny/2+y)*refstridey+(j-ncenx/2+x)*refstridex]-indata[(y)*instridey+x*instridex];
	    }
	  }
	  s+=d*d;
	}
      }
      printf("%g %g %d %d\n",s/(mx*my),s,mx,my);
      outdata[i*outstridey+j*outstridex]=s/(my*mx);
    }
    }*/
    /*
    //compute number of pixels to correlate over:
    //Several cases to consider.
    //starting location is nceny/2-i from centre of ref.
    //So refy/2-(nceny/2-i).
    sy=refy/2-(nceny/2-i);//get the starting point in the ref array.
    my=iny;
    if(sy<0){
      my+=sy;//reduce my.
      sy=0;
    }
    //1. inarr is entirely within ref.
    
    if(iny<refy){
      if(refy/2>=nceny/2-i)
	
    
    my=corrsize[0]-abs(iny-refy)/2-i;
    
    
    my=corrsize[0]-abs(nceny/2-i);
    for(j=0;j<ncenx;j++){
      mx=corrsize[1]-abs(ncenx/2-j);
      s=0;
      for(y=0;y<my;y++){
	for(x=0;x<mx;x++){
	  if(i<ncen/2){
	    if(j<ncen/2){
	      d=refarr[y*corrsize+x]-inarr[(ncen/2-i+y)*nimg+ncen/2-j+x];
	    }else{
	      d=refarr[y*corrsize+j-ncen/2+x]-inarr[(ncen/2-i+y)*nimg+x];
	    }
	  }else{
	    if(j<ncen/2){
	      d=corrPattern[(i-ncen/2+y)*corrsize+x]-bimg[(y)*nimg+ncen/2-j+x];
	    }else{
	      d=corrPattern[(i-ncen/2+y)*corrsize+j-ncen/2+x]-bimg[(y)*nimg+x];
	    }
	  }
	  s+=d*d;
	}
      }
      output[(i+(nimg-ncen)/2)*corrsize+j+(nimg-ncen)/2]=s/(my*mx);
    }
    }*/
  return (PyObject*)outarr;
}
  
static PyMethodDef UtilsMethods[] = {
  {"rotateArray",rotateArray,METH_VARARGS,"Rotate a 2D array"},
  {"compressFloatArray",CompressFloatArray,METH_VARARGS,"Compress floating point"},
  {"uncompressFloatArray",UncompressFloatArray,METH_VARARGS,"uncompress floating point"},
  {"compressFloatArrayAll",CompressFloatArrayAll,METH_VARARGS,"Compress floating point"},
  {"uncompressFloatArrayAll",UncompressFloatArrayAll,METH_VARARGS,"uncompress floating point"},
  {"uncompressFloatArrayThreaded",UncompressFloatArrayThreaded,METH_VARARGS,"uncompress floating point"},
  {"uncompressFloatArrayAllThreaded",UncompressFloatArrayAllThreaded,METH_VARARGS,"uncompress floating point"},
  {"arrayFromArray",  ArrayFromArray, METH_VARARGS,
   "Create an array using the memory from an existing array."},
  {"arrayFromDiagonal",ArrayFromDiagonal,METH_VARARGS,
   "Return diagonal of array as a new array (shared data)"},
  {"mmapArray",mmapArray,METH_VARARGS,"Create a mmap'd array"},
  {"mmapArrayFree",mmapArrayFree,METH_VARARGS,"Free a mmap'd array"},
  {"radialProfileEncircledEnergy",RadialProfileEncircledEnergy,METH_VARARGS,
  "Compute the radial profile, encircled energy, FWHM value and d50 value.  Returns (a[3,x],fwhm,d50) where a[0] is radius, a[1] is profile, a[2] is enc energy"},
  {"centroid",Centroid,METH_VARARGS,"Compute centroid locations"},
  {"queryNumericArray",  queryNumericArray, METH_VARARGS,
   "Find out about a numeric array."},
  {"arrayfrombuffer", arrayfrombuffer, METH_VARARGS, "Create a Numeric array from an existing mmap.mmap object"},
  {"dot", dot, METH_VARARGS, "Matrix-Vector multiplication. Usage:\nc = dot(a, b, r, N)\n       a ... 2-D array - input matrix\n       b ... 1-D array - input vector\n       r ... 1-D array - result vector (optional)\n       N ... the number of threads (optional)"},
  {"correlateDifferenceSquared",correlateDifferenceSquared,METH_VARARGS,"Computed the difference squared correlation"},
  {NULL, NULL, 0, NULL}        /* Sentinel */
};
//PyMODINIT_FUNC 
void initutils(void)
{
  PyObject *m;
  PyImport_AddModule("utils");
  m=Py_InitModule("utils", UtilsMethods);
  import_array();
  UtilsError = PyErr_NewException("utils.error", NULL, NULL);
  Py_INCREF(UtilsError);
  PyModule_AddObject(m, "error", UtilsError);
}


int
main(int argc, char *argv[])
{
    // Pass argv[0] to the Python interpreter 
  Py_SetProgramName(argv[0]);
  
  // Initialize the Python interpreter.  Required. 
  Py_Initialize();
  
  // Add a static module 
  initutils();
  return 0;
}

