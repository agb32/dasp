#include <Python.h>
//#include </usr/local/include/python2.4/Numeric/arrayobject.h>
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
//taken from Numeric/Src/arrayobject.c - has comment that "obviously this needs some work".
//#define ISCONTIGUOUS(m) ((m)->flags & CONTIGUOUS)

static PyObject *UtilsError;

/*
static PyObject* StringFromArray(PyObject *self,PyObject *args){
  //This is dangerous to use unless you know what you're doing.  Python stirng objects should be immutable, but the one returned here can be mutated using the array.  This is typically only used when a tostring() would otherwise be used and the sting is instantly forgotton about.  Not sure what happens when the string is deleted... probably will destroy the array too...
  PyArrayObject *arr;
  PyStringObject *s;
  int i;
  if(!PyArg_ParseTuple(args,"O",&PyArray_Type,&arr)){
    printf("Usage: array\n");
    return NULL;
  }
  if(!PyArray_ISCONTIGUOUS(arr)){
    printf("must be contiguous\n");
    return NULL;
  }
  s=(PyStringObject*)PyString_FromString("");
  &(s->ob_sval[0])=arr->data;
  s->ob_shash=-1;
  s->ob_size=arr->descr->elsize;
  for(i=0; i<arr->nd; i++){
    s->ob_size*=arr->dimensions[i];
  }
  Py_INCREF(s);
  return (PyObject*)s;
  }*/

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
      truncate64(filename,size);
    }
  }else{//file doesn't exist..
    //create file and make it right size...
    printf("Creating file %s for mmap of size %ld bytes\n",filename,size);
    if((fd=open(filename,O_RDWR|O_CREAT,S_IRWXU|S_IRWXG|S_IRWXO))>=0){
      close(fd);
      if(truncate64(filename,size)<0){
	printf("Call to truncate64 failed\n");
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


static PyMethodDef UtilsMethods[] = {
  {"rotateArray",rotateArray,METH_VARARGS,"Rotate a 2D array"},
  {"compressFloatArray",CompressFloatArray,METH_VARARGS,"Compress floating point"},
  //{"stringFromArray",StringFromArray,METH_VARARGS,"Dangerous..."},
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

