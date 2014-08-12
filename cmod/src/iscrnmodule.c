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
      mvm(nX,nX,1.,&ss->Ax[i*nX],1,nX*nbCol,&screen[indx*nX],1,(double)(i!=0),AZ,1);
    }
    coeffTurb=powf(ss->L0/r0,5./6);
    //Now create a random array of standard normal distribution.
    for(i=0;i<nX;i++)
      ss->randn[i]=gsl_ran_ugaussian(ss->rng);//randn();
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
    
/* define a methods table for this module */

static PyMethodDef iscrn_methods[] = 	{
  {"run", py_run, METH_VARARGS}, 
  {"initialise", py_initialise, METH_VARARGS}, 
  {"update", py_update, METH_VARARGS}, 
  {"free", py_free, METH_VARARGS}, 
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
