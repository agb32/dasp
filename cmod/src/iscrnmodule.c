/*Module to create phase screen in c.
Used by science/iscrn.py
*/

#include "Python.h"
#include <stdio.h>
#include <math.h>

#include "numpy/arrayobject.h"

/*#include <gsl/gsl_cblas.h>*/
#include <cblas.h>
/* #include <atlas/cblas.h>
 * NAB 02/Apr/2013 */
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/* Older (<=1.5.1, perhaps <1.7.0) versions of numpy don't have IS_X_CONTIGUOUS
 * defined in "ndarraytypes.h", so it is added here.
 * Added NAB 09/Apr/2013
 */
#ifndef NPY_ARRAY_C_CONTIGUOUS
	#define NPY_ARRAY_C_CONTIGUOUS NPY_C_CONTIGUOUS
#endif
#ifndef NPY_ARRAY_F_CONTIGUOUS
	#define NPY_ARRAY_F_CONTIGUOUS NPY_F_CONTIGUOUS
#endif
#ifndef PyArray_IS_C_CONTIGUOUS(m)
	#define PyArray_IS_C_CONTIGUOUS(m) PyArray_CHKFLAGS(m, NPY_ARRAY_C_CONTIGUOUS)
#endif
#ifndef PyArray_IS_F_CONTIGUOUS(m)
	#define PyArray_IS_F_CONTIGUOUS(m) PyArray_CHKFLAGS(m, NPY_ARRAY_F_CONTIGUOUS)
#endif

#define nbCol 2

typedef enum CBLAS_ORDER CBLAS_ORDER;
typedef enum CBLAS_TRANSPOSE CBLAS_TRANSPOSE;

/*Needs:
r0 (which can change).
scrnYPxls, scrnXPxls
screen
nbCol
Ay
L0
A seed.
nbColToAdd==1
By
xstep
Bx
Ax
AStartX
AStartY

nrew,nadd,interppos (or newCols.next())

sendWholeScreen
maxColAdd
maxRowAdd

*/

static PyObject *iScrnError;

#define R0 1
//#define R0X 2
//#define R0Y 3
//#define XSTEP 4
#define YSTEP 5
//#define RANDN 6

typedef struct{
  double stepSize;
  int extraCols;
  double interpPosition;
  int nadd;
}NewRowsStruct;

typedef struct{
  float r0;
  //float r0x;
  //float r0y;
  float L0;
  int scrnXPxls;
  int scrnYPxls;
  //int nbCol;//2 always
  double *screen;
  //double *Ay;//covariance matrices stuff
  //double *By;
  //double *AStartY;
  double *Ax;
  double *Bx;
  double *AStartX;
  //double Xstep;//for adding a step in the phase.
  double Ystep;
  //double *XstepArr;
  double *YstepArr;
  double *randn;
  //double *randnpy;
  int sendWholeScreen;
  //int maxColAdd;
  //double colAdd;
  double rowAdd;
  int maxRowAdd;
  int nthreads;
  NewRowsStruct nrStruct;
  //int colsAdded;
  int rowsAdded;
  //double *colOutput;
  double *rowOutput;
  gsl_rng *rng;
  int insertPos;
}ScrnStruct;





int checkContigDoubleSize(PyArrayObject *a,int s){
  if(PyArray_SIZE(a)!=s)
    return 1;
  if(a->descr->type_num!=NPY_DOUBLE)
    return 1;
  if(!PyArray_IS_C_CONTIGUOUS(a))
    return 1;
  return 0;
}
int checkContigDouble(PyArrayObject *a){
  if(a->descr->type_num!=NPY_DOUBLE)
    return 1;
  if(!PyArray_IS_C_CONTIGUOUS(a))
    return 1;
  return 0;
}
int checkContigFloat(PyArrayObject *a){
  if(a->descr->type_num!=NPY_FLOAT)
    return 1;
  if(!PyArray_IS_C_CONTIGUOUS(a))
    return 1;
  return 0;
}
int checkFContigDoubleSize(PyArrayObject *a,int s){
  if(PyArray_SIZE(a)!=s)
    return 1;
  if(a->descr->type_num!=NPY_DOUBLE)
    return 1;
  if(!PyArray_IS_F_CONTIGUOUS(a))
    return 1;
  return 0;
}
int checkContigSize(PyArrayObject *a,int s){
  if(!PyArray_IS_C_CONTIGUOUS(a))
    return 1;
  if(PyArray_NBYTES(a)!=s)
    return 1;
  return 0;
}

int nextNewRows(NewRowsStruct *nr){
  double newpos=nr->interpPosition+nr->stepSize;
  int nremove=(int)newpos;
  nr->nadd=(int)(ceil(newpos)-nr->extraCols);
  nr->extraCols+=nr->nadd-nremove;
  nr->interpPosition=newpos-nremove;
  return nr->nadd;
}

inline int mvm(int m,int n,double alpha,double *A,int lda,int ldb,double *x,int incx,double beta,double *y,int incy){
  //performs:
  //y=beta*y + alpha* A dot x
  //todo();
  if(lda!=1)
    printf("Error I think...\n");
  cblas_dgemv(CblasRowMajor,CblasNoTrans,m,n,alpha,A,ldb,x,incx,beta,y,incy);
  return 0;
}
/*
int mvm2(int m,int n,double alpha,double *A,int lda,int ldb,double *x,int incx,double beta,double *y,int incy){
  //y=beta*y+alpha* A.x
  int i,j;
  double tmp;
  for(i=0;i<m;i++){
    y[i*incy]*=beta;
    tmp=0.;
    for(j=0;j<n;j++){
      tmp+=A[i*ldb+j*lda]*x[j*incx];
    }
    y[i*incy]+=alpha*tmp;
  }
  return 0;
  }*/


void addNewRows(ScrnStruct *ss,int nadd){
  float r0;
  int nX,nY;
  double *AZ;
  double *screen=ss->screen;
  int i,ip,indx;
  double coeffTurb;
  //printf("addNewRowOnEnd\n");
  r0=ss->r0;
  nY=ss->scrnYPxls;
  nX=ss->scrnXPxls;
  while(nadd>0){
    ip=ss->insertPos;
    //reset the working array (ie the row where we're adding new data).
    //memset(&ss->screen[ip*nX],0,sizeof(double)*nX);
    //No need - do it with beta in the mvm.
    for(i=0;i<nbCol;i++){//2
      indx=ip-nbCol+i;
      if(indx<0)//wrap
	indx+=nY;
      //Now dot Ax[:,pos:pos+nX] with oldPhi (screen[indx*nX], placing into screen[ip*nx].
      AZ=&screen[ip*nX];
      //printf("addNewRow %d %g %g %g ",ip,ss->Ax[i*nX],screen[indx*nX],AZ[ip*nX]);
      mvm(nX,nX,1.,&ss->Ax[i*nX],1,nX*nbCol,&screen[indx*nX],1,(double)(i!=0),AZ,1);
    }
    coeffTurb=powf(ss->L0/r0,5./6);
    //Now create a random array of standard normal distribution.
    for(i=0;i<nX;i++)
      ss->randn[i]=gsl_ran_ugaussian(ss->rng);//randn();
    //printf("l0 %g r0 %g cT %g %g %g %g %g ",ss->L0,r0,coeffTurb,ss->randn[0],ss->randn[1],ss->Bx[0],AZ[0]);
    mvm(nX,nX,coeffTurb,ss->Bx,1,nX,ss->randn,1,1.,AZ,1);//Dot By with random normal values, multiply by coeffTurb, and add result into AZ.
    //Add xstep...
    if(ss->YstepArr!=NULL){
      for(i=0;i<nX;i++)
	AZ[i]+=ss->YstepArr[i];
    }else if(ss->Ystep!=0){
      for(i=0;i<nX;i++)
	AZ[i]+=ss->Ystep;
    }
    nadd--;
    ss->insertPos++;
    if(ss->insertPos>=nY)//wrap
      ss->insertPos=0;
  }
  //printf("\n");
}



int addNewData(ScrnStruct *ss){
  //computeR0 already called in python.
  //Work out how many to add and remove.
  //now add the rows
  int nadd=nextNewRows(&ss->nrStruct);
  ss->rowsAdded=nadd;
  addNewRows(ss,nadd);
  return 0;
}

int prepareOutput(ScrnStruct *ss){
  //if not sending whole screen, copy the parts to be sent.
  int nX,nY,pos,i;
  //printf("prepareOutput %d %d\n",ss->maxColAdd,ss->colsAdded);
  if(ss->sendWholeScreen==0){
    nX=ss->scrnXPxls;
    nY=ss->scrnYPxls;
    for(i=0;i<ss->maxRowAdd;i++){
      pos=ss->insertPos-ss->maxRowAdd+i;
      if(pos<0)
	pos+=nY;//wrap
      memcpy(&ss->rowOutput[nX*i],&ss->screen[nX*pos],sizeof(double)*nX);
    }
  }
  //printf("Done\n");
  return 0;
}



PyObject *py_initialise(PyObject *self,PyObject *args){
  int nthreads;
  float r0,L0;
  int size,scrnXPxls,scrnYPxls,sendWholeScreen,maxRowAdd;
  double rowAdd;
  PyArrayObject *screen;
  PyArrayObject *Ax;
  PyArrayObject *Bx;
  PyArrayObject *AStartX;
  PyArrayObject *randn;
  PyObject *YstepObj;
  double *YstepArr=NULL;
  double Ystep=0.;
  int seed;
  double *rowOutput=NULL;
  ScrnStruct *ss;
  PyObject *rowOutputObj;
  npy_intp dimsize;
  if(!PyArg_ParseTuple(args,"iffiiiidiO!O!O!O!OO!O",&nthreads,&r0,&L0,&scrnXPxls,&scrnYPxls,&sendWholeScreen,&maxRowAdd,&rowAdd,&seed,&PyArray_Type,&screen,&PyArray_Type,&Ax,&PyArray_Type,&Bx,&PyArray_Type,&AStartX,&YstepObj,&PyArray_Type,&randn,&rowOutputObj)){
    printf("Args for scrnmodule.initialise should be:\n");
    printf("iffiiiiiffiO!O!O!O!O!O!O!OOO!OO,&nthreads,&r0,&L0,&scrnXPxls,&scrnYPxls,&sendWholeScreen,&maxRowAdd,&rowAdd,&seed,&PyArray_Type,&screen,&PyArray_Type,&Ax,&PyArray_Type,&Bx,&PyArray_Type,&AStartX,&YstepObj,&PyArray_Type,&randn,&rowOutputObj\n");
    return NULL;
  }

  size=scrnXPxls;
  if(checkContigDoubleSize(Ax,size*size*2)!=0){
    printf("Error, Ax must have shape ((scrnXPxls+1),(scrnXPxls+1)*2), be contiguous and float64\n");
    return NULL;
  }
  if(checkContigDoubleSize(AStartX,size*size*2)!=0){
    printf("Error, AStartX must have shape ((scrnXPxls+1),(scrnXPxls+1)*2), be contiguous and float64\n");
    return NULL;
  }
  if(checkContigDoubleSize(Bx,size*size)!=0){
    printf("Error, Bx must have shape ((scrnXPxls+1),(scrnXPxls+1)), be contiguous and float64\n");
    return NULL;
  }
  if(checkContigDoubleSize(screen,(scrnYPxls)*(scrnXPxls))!=0){
    printf("screen must have shape %d,%d be contiguous and float64\n",scrnYPxls,scrnXPxls);
    return NULL;
  }
  if(checkContigDoubleSize(randn,(scrnXPxls))!=0){
    printf("random array space must have shape %d becontiguous and float64\n",scrnXPxls+1);
    return NULL;
  }
  if(PyArray_Check(YstepObj)){
    if(checkContigDoubleSize((PyArrayObject*)YstepObj,scrnXPxls)!=0){
      printf("Ystep shoud be size scrnXPxls, contiguous, float64\n");
      return NULL;
    }
    YstepArr=(double*)(((PyArrayObject*)YstepObj)->data);
  }else{
    Ystep=(double)PyFloat_AsDouble(YstepObj);
  }
  if(PyArray_Check(rowOutputObj)){
    if(checkContigDoubleSize((PyArrayObject*)rowOutputObj,(scrnXPxls)*maxRowAdd)!=0){
      printf("rowOutput should be size (scrnXPxls+1)*maxRowAdd, contiguous, float64\n");
      return NULL;
    }
    rowOutput=(double*)(((PyArrayObject*)rowOutputObj)->data);
  }
  if(sendWholeScreen==0 && (rowOutput==NULL )){
    printf("If not sending whole screen, row output must be defined\n");
    return NULL;
  }

  if((ss=calloc(sizeof(ScrnStruct),1))==NULL){
    printf("Failed to calloc ScrnStruct\n");
    return NULL;
  }
  ss->rng=gsl_rng_alloc(gsl_rng_mt19937);
  if(seed==0)
    gsl_rng_set(ss->rng,(unsigned long)time(0));
  else
    gsl_rng_set(ss->rng,seed);
  ss->nthreads=nthreads;
  ss->r0=r0;
  ss->L0=L0;
  ss->scrnXPxls=scrnXPxls;
  ss->scrnYPxls=scrnYPxls;
  ss->sendWholeScreen=sendWholeScreen;
  ss->maxRowAdd=maxRowAdd;
  ss->rowAdd=rowAdd;
  ss->Ax=(double*)Ax->data;
  ss->AStartX=(double*)AStartX->data;
  ss->Bx=(double*)Bx->data;
  ss->screen=(double*)screen->data;
  ss->randn=(double*)randn->data;
  ss->YstepArr=YstepArr;
  ss->Ystep=Ystep;
  ss->nrStruct.stepSize=fabs(rowAdd);
  ss->rowOutput=rowOutput;
  dimsize=sizeof(ScrnStruct);
  return PyArray_SimpleNewFromData(1,&dimsize,NPY_UBYTE,ss);
}


PyObject *py_update(PyObject *self,PyObject *args){
  int code;
  ScrnStruct *ss;
  PyArrayObject *ssArr;
  PyObject *obj;
  if(!PyArg_ParseTuple(args,"O!iO",&PyArray_Type,&ssArr,&code,&obj)){
    printf("Usage: ScrnStruct object, code for value to be changed, new value\n");
    return NULL;
  }
  if(checkContigSize(ssArr,sizeof(ScrnStruct))!=0){
    printf("ScrnStruct should be initialised with the initialise method of scrn module\n");
    return NULL;
  }
  ss=(ScrnStruct*)(ssArr->data);
  PyErr_Clear();
  switch(code){
  case R0:
    ss->r0=PyFloat_AsDouble(obj);
    if(PyErr_Occurred()){
      printf("scrnmodule: Error extracting float value for r0\n");
      return NULL;
    }
    break;
  case YSTEP:
    if(PyArray_Check(obj)){
      if(checkContigDoubleSize((PyArrayObject*)obj,ss->scrnXPxls)!=0){
	printf("Ystep shoud be size scrnXPxls, contiguous, float64\n");
	return NULL;
      }
      ss->YstepArr=(double*)(((PyArrayObject*)obj)->data);
    }else{
      ss->Ystep=(double)PyFloat_AsDouble(obj);
      if(PyErr_Occurred()){
	printf("scrnmodule: Error extracting double value for ystep\n");
	return NULL;
      }
    }
    break;
    /*case RANDN:
    if(PyArray_Check(obj)){
      if(checkContigDoubleSize((PyArrayObject*)obj,ss->scrnXPxls>ss->scrnYPxls?ss->scrnXPxls+1:ss->scrnYPxls+1)!=0){
	printf("RANDNPY should be size max(scrnXPxls,scrnYPxls)+1, contiguous, float64\n");
	return NULL;
      }
      ss->randnpy=(double*)(((PyArrayObject*)obj)->data);
    }else{
      ss->randnpy=NULL;
    }
    break;*/
  default:
    printf("Unrecognised parameter %d in cmod.scrn.update\n",code);
    return NULL;
    break;
  }
  Py_INCREF(Py_None);
  return Py_None;
}  

PyObject *py_free(PyObject *self,PyObject *args){
  ScrnStruct *ss;
  PyArrayObject *ssArr;
  if(!PyArg_ParseTuple(args,"O!",&PyArray_Type,&ssArr)){
    printf("Usage: ScrnStruct\n");
    return NULL;
  }
  if(checkContigSize(ssArr,sizeof(ScrnStruct))!=0){
    printf("ScrnStruct should be initialised with the initialise method of scrn module\n");
    return NULL;
  }
  ss=(ScrnStruct*)ssArr->data;
  //Now free anything that needs freeing.  
  //is it worth changing the size/shape of sarr so that it can't be used again?
  gsl_rng_free(ss->rng);
  //Do we need to free ssArr?  I don't think so, but could be wrong...
  Py_INCREF(Py_None);
  return Py_None;
}

PyObject *py_run(PyObject *self,PyObject *args){
  ScrnStruct *ss;
  PyArrayObject *ssArr;
  if(!PyArg_ParseTuple(args,"O!",&PyArray_Type,&ssArr)){
    printf("Usage: ScrnStruct\n");
    return NULL;
  }
  if(checkContigSize(ssArr,sizeof(ScrnStruct))!=0){
    printf("ScrnStruct should be initialised with the initialise method of scrn module\n");
    return NULL;
  }
  ss=(ScrnStruct*)ssArr->data;
  Py_BEGIN_ALLOW_THREADS;
  addNewData(ss);
  prepareOutput(ss);
  Py_END_ALLOW_THREADS;
  //Py_INCREF(Py_None);
  return Py_BuildValue("i",ss->insertPos);//Py_None;
}
    
typedef struct{
  double *img;
  double *dimg;
  float *out;
  float r;
  float s;
  float c;
  int dim[2];
  int imgdim[2];
  int nthreads;
  int nblockx;
  int nblocky;
} interpStruct;

PyObject *py_initialiseInterp(PyObject *self,PyObject *args){
  npy_intp *dim;
  npy_intp *imgdim;
  float ang;
  int nthreads,nblockx,nblocky;
  float deg,r;
  PyArrayObject *imgObj, *dimgObj, *outObj;
  interpStruct *s;
  npy_intp dimsize;
  if(!PyArg_ParseTuple(args,"O!O!fO!fiii",&PyArray_Type,&imgObj,&PyArray_Type,&dimgObj,&deg,&PyArray_Type,&outObj,&r,&nthreads,&nblockx,&nblocky)){
    printf("Args for setupInterpolate should be:\n");
    printf("img,dimg,deg,out,r,nthread,nblockx,nblocky\n");
    return NULL;
  }
  if(PyArray_NDIM(imgObj)!=2 || PyArray_NDIM(outObj)!=2){
    printf("Arrays should be 2D\n");
    return NULL;
  }


  dim=PyArray_DIMS(outObj);
  imgdim=PyArray_DIMS(imgObj);

  if(checkContigDoubleSize(dimgObj,imgdim[0]*imgdim[1])!=0){
    printf("Error, dimg must have shape same as img (%ld, %ld), be contiguous and float64\n",imgdim[0],imgdim[1]);
    return NULL;
  }
  if(checkContigDouble(imgObj)!=0){
    printf("img should be contiguous and float64\n");
    return NULL;
  }
  if(checkContigFloat(outObj)!=0){
    printf("out should be contiguous and float32\n");
    return NULL;
  }
  if((s=calloc(sizeof(interpStruct),1))==NULL){
    printf("Unable to alloc interpStruct memory\n");
    return NULL;
  }
  dimsize=sizeof(interpStruct);
  s->img=(double*)PyArray_DATA(imgObj);
  s->dimg=(double*)PyArray_DATA(dimgObj);
  s->out=(float*)PyArray_DATA(outObj);
  ang=(float)(deg*M_PI/180.);
  s->r=r;
  s->s=(float)(r*sin(ang));
  s->c=(float)(r*cos(ang));
  s->dim[0]=(int)dim[0];
  s->dim[1]=(int)dim[1];
  s->imgdim[0]=(int)imgdim[0];
  s->imgdim[1]=(int)imgdim[1];
  s->nthreads=nthreads;
  s->nblockx=nblockx;
  s->nblocky=nblocky;
  return PyArray_SimpleNewFromData(1,&dimsize,NPY_UBYTE,s);  
}


PyObject *py_rotShiftWrapSplineImage(PyObject *self,PyObject *args){
  interpStruct *ss;
  int outofrange[4];
  float points[4];
  //float ang;
  float sx,sy;
  int wrappoint,nthreads,nblockx,nblocky;
  //PyArrayObject *imgObj, *dimgObj, *outObj;
  PyArrayObject *structObj;
  float s;
  float c;
  int *dim;
  int *imgdim;
  double *img;
  double *dimg;
  float *out;
  int i,yy,xx;
  float x,y,xold,yold;
  int x1;
  int y1,y2,y1x1,y2x1;
  float xm,ym;
  float k1,k2,Y1,Y2,a,b,val;
  float const0,const1,const2,const3;
  if(!PyArg_ParseTuple(args,"O!ffi",&PyArray_Type,&structObj,&sx,&sy,&wrappoint)){
    printf("Args for rotShiftWrapSplineImage should be:\n");
    printf("struct returned from initialiseInterp,shiftx,shifty,wrappoint\n");
    return NULL;
  }
  if(checkContigSize(structObj,sizeof(interpStruct))!=0){
    printf("interpStruct should be initialised with the initialiseInterp method of iscrn module\n");
    return NULL;
  }
  ss=(interpStruct*)(structObj->data);

  dim=ss->dim;
  imgdim=ss->imgdim;
  img=ss->img;
  dimg=ss->dimg;
  out=ss->out;

  if(wrappoint>imgdim[0] || wrappoint<0){
    printf("Illegal wrappoint in iscrnmodule %d\n",wrappoint);
    return NULL;
  }
  const0=-dim[0]/2.+0.5;
  const1=-dim[1]/2.+0.5;
  const2=imgdim[0]/2.-.5-sy+wrappoint;
  const3=imgdim[1]/2.-.5-sx;
  //ang=ss->ang;//(float)(deg*M_PI/180.);
  s=ss->s;//(float)(r*sin(ang));
  c=ss->c;//(float)(r*cos(ang));
  for(yy=0;yy<dim[0];yy++){
    y=yy+const0;//-dim[0]/2.+0.5;
    for(xx=0;xx<dim[1];xx++){
      x=xx+const1;//-dim[1]/2.+0.5;
      yold=-s*x+c*y+const2;//imgdim[0]/2.-.5-sy+wrappoint;
      xold=c*x+s*y+const3;//imgdim[1]/2.-.5-sx;
      x1=(int)floorf(xold);
      //First, we need to compute 4 splines in the y direction.  These are then used to compute in the x direction, giving the final value.
      y1=(int)floorf(yold);
      if(y1>=imgdim[0]){//wrap it
	y1-=imgdim[0];
	yold-=imgdim[0];
      }
      y2=y1+1;
      xm=xold-x1;
      ym=yold-y1;
      if(y2==imgdim[0])
	y2=0;
      x1--;
      y1x1=y1*imgdim[1]+x1;
      y2x1=y2*imgdim[1]+x1;
      //4 interpolations in Y direction using precomputed gradients.
      for(i=0;i<4;i++){//at x1-1, x1, x2, x2+1.
	if(x1+i>=0 && x1+i<imgdim[1]){
	  k1=(float)dimg[y1x1+i];
	  k2=(float)dimg[y2x1+i];
	  Y1=(float)img[y1x1+i];
	  Y2=(float)img[y2x1+i];
	  a=k1-(Y2-Y1);//k1*(X2-X1)-(Y2-Y1)
	  b=-k2+(Y2-Y1);//-k2*(X2-X1)+(Y2-Y1)
	  points[i]=((1-ym)*Y1+ym*Y2+ym*(1-ym)*(a*(1-ym)+b*ym));
	  outofrange[i]=0;
	}else{
	  outofrange[i]=1;
	}
      }
      //and now interpolate in X direction (using points).
      if(outofrange[0])
	k1=points[2]-points[1];
      else
	k1=(points[2]-points[0])*.5;
      if(outofrange[3])
	k2=points[2]-points[1];
      else
	k2=(points[3]-points[1])*.5;
      if(outofrange[1] || outofrange[2]){
	printf("Out of range y:%d x:%d %g %g %d %d,%d %d,%d %d\n",yy,xx,sx,sy,wrappoint,dim[0],dim[1],imgdim[0],imgdim[1],x1+1);
      }else{
	Y1=points[1];
	Y2=points[2];
	a=k1-(Y2-Y1);//k1*(X2-X1)-(Y2-Y1)
	b=-k2+(Y2-Y1);//-k2*(X2-X1)+(Y2-Y1)
	val=(1-xm)*Y1+xm*Y2+xm*(1-xm)*(a*(1-xm)+b*xm);
	out[yy*dim[1]+xx]+=val;
      }
    }
  }
  Py_INCREF(Py_None);
  return Py_None;
}




PyObject *py_rotShiftWrapSplineImageNoInit(PyObject *self,PyObject *args){
  int outofrange[4];
  float points[4];
  float ang,deg,r;
  float sx,sy;
  int wrappoint,nthreads,nblockx,nblocky;
  PyArrayObject *imgObj, *dimgObj, *outObj;
  float s;
  float c;
  npy_intp *dim;
  npy_intp *imgdim;
  double *img;
  double *dimg;
  float *out;
  int i,yy,xx;
  float x,y,xold,yold;
  int x1;
  int y1,y2,y1x1,y2x1;
  float xm,ym;
  float k1,k2,Y1,Y2,a,b,val;
  float const0,const1,const2,const3;

  if(!PyArg_ParseTuple(args,"O!O!fffiO!fiii",&PyArray_Type,&imgObj,&PyArray_Type,&dimgObj,&deg,&sx,&sy,&wrappoint,&PyArray_Type,&outObj,&r,&nthreads,&nblockx,&nblocky)){
    printf("Args for rotShiftWrapSplineImageNoInit should be:\n");
    printf("img,gradients,angle,shiftx,shifty,wrappoint,outputarray,nthreads,nblockx,nblocky\n");
    return NULL;
  }

  if(PyArray_NDIM(imgObj)!=2 || PyArray_NDIM(outObj)!=2){
    printf("Arrays should be 2D\n");
    return NULL;
  }


  dim=PyArray_DIMS(outObj);
  imgdim=PyArray_DIMS(imgObj);

  if(checkContigDoubleSize(dimgObj,imgdim[0]*imgdim[1])!=0){
    printf("Error, dimg must have shape same as img (%ld, %ld), be contiguous and float64\n",imgdim[0],imgdim[1]);
    return NULL;
  }
  if(checkContigDouble(imgObj)!=0){
    printf("img should be contiguous and float64\n");
    return NULL;
  }
  if(checkContigFloat(outObj)!=0){
    printf("out should be contiguous and float32\n");
    return NULL;
  }

  img=(double*)PyArray_DATA(imgObj);
  dimg=(double*)PyArray_DATA(dimgObj);
  out=(float*)PyArray_DATA(outObj);
  ang=(float)(deg*M_PI/180.);
  s=(float)(r*sin(ang));
  c=(float)(r*cos(ang));

  if(wrappoint>imgdim[0] || wrappoint<0){
    printf("Illegal wrappoint in iscrnmodule %d\n",wrappoint);
    return NULL;
  }
  const0=-dim[0]/2.+0.5;
  const1=-dim[1]/2.+0.5;
  const2=imgdim[0]/2.-.5-sy+wrappoint;
  const3=imgdim[1]/2.-.5-sx;
  for(yy=0;yy<dim[0];yy++){
    y=yy+const0;//-dim[0]/2.+0.5;
    for(xx=0;xx<dim[1];xx++){
      x=xx+const1;//-dim[1]/2.+0.5;
      yold=-s*x+c*y+const2;//imgdim[0]/2.-.5-sy+wrappoint;
      xold=c*x+s*y+const3;//imgdim[1]/2.-.5-sx;
      x1=(int)floorf(xold);
      //First, we need to compute 4 splines in the y direction.  These are then used to compute in the x direction, giving the final value.
      y1=(int)floorf(yold);
      if(y1>=imgdim[0]){//wrap it
	y1-=imgdim[0];
	yold-=imgdim[0];
      }
      y2=y1+1;
      xm=xold-x1;
      ym=yold-y1;
      if(y2==imgdim[0])
	y2=0;
      x1--;
      y1x1=y1*imgdim[1]+x1;
      y2x1=y2*imgdim[1]+x1;
      //4 interpolations in Y direction using precomputed gradients.
      for(i=0;i<4;i++){//at x1-1, x1, x2, x2+1.
	if(x1+i>=0 && x1+i<imgdim[1]){
	  k1=(float)dimg[y1x1+i];
	  k2=(float)dimg[y2x1+i];
	  Y1=(float)img[y1x1+i];
	  Y2=(float)img[y2x1+i];
	  a=k1-(Y2-Y1);//k1*(X2-X1)-(Y2-Y1)
	  b=-k2+(Y2-Y1);//-k2*(X2-X1)+(Y2-Y1)
	  points[i]=((1-ym)*Y1+ym*Y2+ym*(1-ym)*(a*(1-ym)+b*ym));
	  outofrange[i]=0;
	}else{
	  outofrange[i]=1;
	}
      }
      //and now interpolate in X direction (using points).
      if(outofrange[0])
	k1=points[2]-points[1];
      else
	k1=(points[2]-points[0])*.5;
      if(outofrange[3])
	k2=points[2]-points[1];
      else
	k2=(points[3]-points[1])*.5;
      if(outofrange[1] || outofrange[2]){
	printf("Out of range y:%d x:%d %g %g %d %d,%d %d,%d %d\n",yy,xx,sx,sy,wrappoint,(int)dim[0],(int)dim[1],(int)imgdim[0],(int)imgdim[1],x1+1);
      }else{
	Y1=points[1];
	Y2=points[2];
	a=k1-(Y2-Y1);//k1*(X2-X1)-(Y2-Y1)
	b=-k2+(Y2-Y1);//-k2*(X2-X1)+(Y2-Y1)
	val=(1-xm)*Y1+xm*Y2+xm*(1-xm)*(a*(1-xm)+b*xm);
	out[yy*dim[1]+xx]+=val;
      }
    }
  }
  Py_INCREF(Py_None);
  return Py_None;
}






/* define a methods table for this module */

static PyMethodDef iscrn_methods[] = 	{
  {"run", py_run, METH_VARARGS}, 
  {"initialise", py_initialise, METH_VARARGS}, 
  {"update", py_update, METH_VARARGS}, 
  {"free", py_free, METH_VARARGS}, 
  {"initialiseInterp",py_initialiseInterp,METH_VARARGS},
  {"rotShiftWrapSplineImage", py_rotShiftWrapSplineImage, METH_VARARGS},
  {"rotShiftWrapSplineImageNoInit", py_rotShiftWrapSplineImageNoInit, METH_VARARGS},
  {NULL, NULL} };


/* initialisation - register the methods with the Python interpreter */

void initiscrn(void){
  PyObject *m;
  PyImport_AddModule("iscrn");
  m=Py_InitModule("iscrn", iscrn_methods);
  import_array();
  iScrnError=PyErr_NewException("iscrn.error",NULL,NULL);
  Py_INCREF(iScrnError);
  PyModule_AddObject(m,"error",iScrnError);
}

int main(int argc,char **argv){
  Py_SetProgramName(argv[0]);
  Py_Initialize();
  initiscrn();
  return 0;
}
