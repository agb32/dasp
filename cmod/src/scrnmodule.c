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
/*Module to create phase screen in c.
Used by science/infScrn.py
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



//typedef enum CBLAS_ORDER CBLAS_ORDER;
//typedef enum CBLAS_TRANSPOSE CBLAS_TRANSPOSE;

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

static PyObject *ScrnError;

#define R0 1
#define R0X 2
#define R0Y 3
#define XSTEP 4
#define YSTEP 5
//#define RANDN 6

typedef struct{
  double stepSize;
  int extraCols;
  double interpPosition;
  int nadd;
}NewColsStruct;

typedef struct{
  float r0;
  float r0x;
  float r0y;
  float L0;
  int scrnXPxls;
  int scrnYPxls;
  //int nbCol;//2 always
  double *screen;
  double *Ay;//covariance matrices stuff
  double *By;
  double *AStartY;
  double *Ax;
  double *Bx;
  double *AStartX;
  double Xstep;//for adding a step in the phase.
  double Ystep;
  double *XstepArr;
  double *YstepArr;
  double *randn;
  //double *randnpy;
  int sendWholeScreen;
  int maxColAdd;
  double colAdd;
  double rowAdd;
  int maxRowAdd;
  int nthreads;
  NewColsStruct ncStruct[2];
  int colsAdded;
  int rowsAdded;
  double *colOutput;
  double *rowOutput;
  gsl_rng *rng;
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

int nextNewCols(NewColsStruct *nc){
  double newpos=nc->interpPosition+nc->stepSize;
  int nremove=(int)newpos;
  //int oldExtraCols;
  nc->nadd=(int)(ceil(newpos)-nc->extraCols);
  //oldExtraCols=nc->extraCols;
  nc->extraCols+=nc->nadd-nremove;
  nc->interpPosition=newpos-nremove;
  return nc->nadd;
}

static inline int mvm(int m,int n,double alpha,double *A,int lda,int ldb,double *x,int incx,double beta,double *y,int incy){
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

void addNewColumnsOnEnd(ScrnStruct *ss,int nadd){
  float r0;
  int nX,nY;
  int nToAdd;
  double *AZ;
  double *screen=ss->screen;
  int i;
  double coeffTurb;
  double *randn;
  //printf("addNewColOnEnd\n");
  if(ss->r0x!=0)
    r0=ss->r0x;
  else
    r0=ss->r0;
  nY=ss->scrnYPxls+1;
  nX=ss->scrnXPxls+1;
  while(nadd>0){
    if(nadd<ss->scrnXPxls-1){//do all...
      nToAdd=nadd;
      nadd=0;
    }else{//very fast turbulence!
      //do some.
      nToAdd=ss->scrnXPxls-1;
      nadd-=ss->scrnXPxls-1;
    }
    //First, shift the screen left by nToAdd pixels.
    //Because its contiguous, this works, and we're ignoring the last columns anyway.
    memmove(ss->screen,&ss->screen[nToAdd],sizeof(double)*(nY*nX-nToAdd));
    while(nToAdd>0){
      //printf("nToAdd %d\n",nToAdd);
      //now an mvm of Ay with last 2 columns of screen ravelled.
      //Best to do this in 2 parts.
      //Left half of Ay dot end but one column
      //Right half of Ay dot end column.
      //Then add the 2 vectors together.
      AZ=&screen[nX-nToAdd];//vector with lda of nX.  This is where result goes
      mvm(nY,nY,1.,ss->Ay,1,nY*2,&screen[nX-nToAdd-2],nX,0.,AZ,nX);//dot left half of Ay with end but one column of screen, put result in AZ
      mvm(nY,nY,1.,&ss->Ay[nY],1,nY*2,&screen[nX-nToAdd-1],nX,1.,AZ,nX);//dot right half of Ay with end column of screen, add result into AZ
      coeffTurb=powf(ss->L0/r0,5./6);
      randn=ss->randn;
      for(i=0;i<nY;i++){
	randn[i]=gsl_ran_ugaussian(ss->rng);//randn();
      }
      mvm(nY,nY,coeffTurb,ss->By,1,nY,randn,1,1.,AZ,nX);//Dot By with random normal values, multiply by coeffTurb, and add result into AZ.
      //Add xstep...
      if(ss->XstepArr!=NULL){
	for(i=0;i<nY;i++)
	  AZ[i*nX]+=ss->XstepArr[i];
      }else if(ss->Xstep!=0){
	for(i=0;i<nY;i++)
	  AZ[i*nX]+=ss->Xstep;
      }
      nToAdd--;
    }
  }
}

void addNewColumnsOnStart(ScrnStruct *ss,int nadd){
  float r0;
  int nX,nY;
  int nToAdd;
  double *AZ;
  double *screen=ss->screen;
  int i;
  double coeffTurb;
  //printf("addNewColOnStart\n");
  if(ss->r0x!=0)
    r0=ss->r0x;
  else
    r0=ss->r0;
  nY=ss->scrnYPxls+1;
  nX=ss->scrnXPxls+1;
  while(nadd>0){
    if(nadd<ss->scrnXPxls-1){//do all...
      nToAdd=nadd;
      nadd=0;
    }else{//very fast turbulence!
      //do some.
      nToAdd=ss->scrnXPxls-1;
      nadd-=ss->scrnXPxls-1;
    }
    //First, shift the screen right by nToAdd pixels.
    //Because its contiguous, this works, and we're ignoring the last columns anyway.
    memmove(&ss->screen[nToAdd],ss->screen,sizeof(double)*(nY*nX-nToAdd));
    while(nToAdd>0){

      //now an mvm of Ay with last 2 columns of screen ravelled.
      //Best to do this in 2 parts.
      //Left half of Ay dot end but one column
      //Right half of Ay dot end column.
      //Then add the 2 vectors together.
      AZ=&screen[nToAdd-1];//vector with lda of nX.  This is where result goes
      mvm(nY,nY,1.,ss->AStartY,1,nY*2,&screen[nToAdd],nX,0.,AZ,nX);//dot left half of Ay with end but one column of screen, put result in AZ
      mvm(nY,nY,1.,&ss->AStartY[nY],1,nY*2,&screen[nToAdd+1],nX,1.,AZ,nX);//dot right half of Ay with end column of screen, add result into AZ
      coeffTurb=powf(ss->L0/r0,5./6);
      //Now create a random array of standard normal distribution.
      for(i=0;i<nY;i++)
	ss->randn[i]=gsl_ran_ugaussian(ss->rng);//randn();
      mvm(nY,nY,coeffTurb,ss->By,1,nY,ss->randn,1,1.,AZ,nX);//Dot By with random normal values, multiply by coeffTurb, and add result into AZ.
      //Add xstep...
      if(ss->XstepArr!=NULL){
	for(i=0;i<nY;i++)
	  AZ[i*nX]+=ss->XstepArr[i];
      }else if(ss->Xstep!=0){
	for(i=0;i<nY;i++)
	  AZ[i*nX]+=ss->Xstep;
      }
      nToAdd--;
    }
  }
}


void addNewRowsOnEnd(ScrnStruct *ss,int nadd){
  float r0;
  int nX,nY;
  int nToAdd;
  double *AZ;
  double *screen=ss->screen;
  int i;
  double coeffTurb;
  //printf("addNewRowOnEnd\n");
  if(ss->r0y!=0)
    r0=ss->r0y;
  else
    r0=ss->r0;
  nY=ss->scrnYPxls+1;
  nX=ss->scrnXPxls+1;
  while(nadd>0){
    if(nadd<ss->scrnYPxls-1){//do all...
      nToAdd=nadd;
      nadd=0;
    }else{//very fast turbulence!
      //do some.
      nToAdd=ss->scrnYPxls-1;
      nadd-=ss->scrnYPxls-1;
    }
    //First, shift the screen up by nToAdd pixels.
    //Because its contiguous, this works, and we're ignoring the last columns anyway.
    memmove(ss->screen,&ss->screen[nToAdd*nX],sizeof(double)*(nX*(nY-nToAdd)));
    while(nToAdd>0){

      //now an mvm of Ax with last 2 rows of screen ravelled.
      //Since the rows are contiguous, can do this in one mvm.
      //Best to do this in 2 parts.
      AZ=&screen[nX*(nY-nToAdd)];//vector with lda of 1.  This is where result goes
      mvm(nX,nX*2,1.,ss->Ax,1,nX*2,&screen[nX*(nY-nToAdd-2)],1,0.,AZ,1);//dot of Ax with end but one and two rows of screen, put result in AZ
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
      nToAdd--;
    }
  }
}
void addNewRowsOnStart(ScrnStruct *ss,int nadd){
  float r0;
  int nX,nY;
  int nToAdd;
  double *AZ;
  double *screen=ss->screen;
  int i;
  double coeffTurb;
  //printf("addNewRowOnStart\n");
  if(ss->r0y!=0)
    r0=ss->r0y;
  else
    r0=ss->r0;
  nY=ss->scrnYPxls+1;
  nX=ss->scrnXPxls+1;
  while(nadd>0){
    if(nadd<ss->scrnYPxls-1){//do all...
      nToAdd=nadd;
      nadd=0;
    }else{//very fast turbulence!
      //do some.
      nToAdd=ss->scrnYPxls-1;
      nadd-=ss->scrnYPxls-1;
    }
    //First, shift the screen down by nToAdd pixels.
    //Because its contiguous, this works, and we're ignoring the last columns anyway.
    memmove(&ss->screen[nToAdd*nX],ss->screen,sizeof(double)*(nX*(nY-nToAdd)));
    while(nToAdd>0){

      //now an mvm of Ax with first 2 rows of screen ravelled.
      //Since the rows are contiguous, can do this in one mvm.
      //Best to do this in 2 parts.
      AZ=&screen[nX*(nToAdd-1)];//vector with lda of 1.  This is where result goes
      mvm(nX,nX*2,1.,ss->AStartX,1,nX*2,&screen[nX*nToAdd],1,0.,AZ,1);//dot of Ax with end but one and two rows of screen, put result in AZ
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
      nToAdd--;
    }
  }
}



int addNewData(ScrnStruct *ss){
  //computeR0 already called in python.
  //Work out how many to add and remove.
  int nadd=nextNewCols(&ss->ncStruct[0]);//,&nrem,&nadd,&interppos);
  //printf("addNewData %d %g\n",nadd,ss->ncStruct[0].stepSize);
  ss->colsAdded=nadd;
  if(ss->colAdd<0){
    addNewColumnsOnEnd(ss,nadd);
  }else{
    addNewColumnsOnStart(ss,nadd);
  }
  //now add the rows
  nadd=nextNewCols(&ss->ncStruct[1]);
  ss->rowsAdded=nadd;
  if(ss->rowAdd<0){
    addNewRowsOnEnd(ss,nadd);
  }else{
    addNewRowsOnStart(ss,nadd);
  }
  return 0;
}

int prepareOutput(ScrnStruct *ss){
  //if not sending whole screen, copy the parts to be sent.
  int nX,nY,pos,pos2,i,j;
  //printf("prepareOutput %d %d\n",ss->maxColAdd,ss->colsAdded);
  if(ss->sendWholeScreen==0){
    nX=ss->scrnXPxls+1;
    nY=ss->scrnYPxls+1;
    if(ss->colAdd<0){//added on end

      pos=(ss->maxColAdd-ss->colsAdded)*nY;
      //printf("pos %d\n",pos);
      for(i=0;i<ss->colsAdded;i++){
	pos2=nX-ss->colsAdded+i;
	//printf("pos2 %d\n",pos2);
	for(j=0;j<nY;j++){
	  //printf("pos,pos2 %d %d\n",pos,pos2);
	  ss->colOutput[pos]=ss->screen[pos2];
	  pos++;
	  pos2+=nX;
	}
      }
    }else{//added on start
      pos=0;
      for(i=0;i<ss->colsAdded;i++){
	pos2=i;
	for(j=0;j<nY;j++){
	  ss->colOutput[pos]=ss->screen[pos2];
	  pos++;
	  pos2+=nX;
	}
      }
    }
    //printf("preparing rows\n");
    if(ss->rowAdd<0){//added at bottom
      memcpy(&ss->rowOutput[nX*(ss->maxRowAdd-ss->rowsAdded)],&ss->screen[nX*(nY-ss->rowsAdded)],sizeof(double)*(ss->rowsAdded*nX));
    }else{//added at top
      memcpy(ss->rowOutput,ss->screen,sizeof(double)*ss->rowsAdded*nX);
    }
  }
  //printf("Done\n");
  return 0;
}



PyObject *py_initialise(PyObject *self,PyObject *args){
  int nthreads;
  float r0,L0;
  int size,scrnXPxls,scrnYPxls,sendWholeScreen,maxColAdd,maxRowAdd;
  double colAdd,rowAdd;
  PyArrayObject *screen;
  PyArrayObject *Ay;
  PyArrayObject *By;
  PyArrayObject *AStartY;
  PyArrayObject *Ax;
  PyArrayObject *Bx;
  PyArrayObject *AStartX;
  PyArrayObject *randn;
  PyObject *XstepObj;
  PyObject *YstepObj;
  double *XstepArr=NULL;
  double *YstepArr=NULL;
  double Xstep=0.;
  double Ystep=0.;
  int seed;
  double *rowOutput=NULL,*colOutput=NULL;
  ScrnStruct *ss;
  PyObject *colOutputObj;
  PyObject *rowOutputObj;
  npy_intp dimsize;
  if(!PyArg_ParseTuple(args,"iffiiiiiddiO!O!O!O!O!O!O!OOO!OO",&nthreads,&r0,&L0,&scrnXPxls,&scrnYPxls,&sendWholeScreen,&maxColAdd,&maxRowAdd,&colAdd,&rowAdd,&seed,&PyArray_Type,&screen,&PyArray_Type,&Ay,&PyArray_Type,&By,&PyArray_Type,&AStartY,&PyArray_Type,&Ax,&PyArray_Type,&Bx,&PyArray_Type,&AStartX,&XstepObj,&YstepObj,&PyArray_Type,&randn,&colOutputObj,&rowOutputObj)){
    printf("Args for scrnmodule.initialise should be:\n");
    printf("iffiiiiiffiO!O!O!O!O!O!O!OOO!OO,&nthreads,&r0,&L0,&scrnXPxls,&scrnYPxls,&sendWholeScreen,&maxColAdd,&maxRowAdd,&colAdd,&rowAdd,&seed,&PyArray_Type,&screen,&PyArray_Type,&Ay,&PyArray_Type,&By,&PyArray_Type,&AStartY,&PyArray_Type,&Ax,&PyArray_Type,&Bx,&PyArray_Type,&AStartX,&XstepObj,&YstepObj,&PyArray_Type,&randn,&colOutputObj,&rowOutputObj\n");
    return NULL;
  }

  size=scrnYPxls+1;
  if(checkContigDoubleSize(Ay,size*size*2)!=0){
    printf("Error, Ay must have shape ((scrnYPxls+1),(scrnYPxls+1)*2), be contiguous and float64\n");
    return NULL;
  }
  if(checkContigDoubleSize(AStartY,size*size*2)!=0){
    printf("Error, AStartY must have shape ((scrnYPxls+1),(scrnYPxls+1)*2), be contiguous and float64\n");
    return NULL;
  }
  if(checkContigDoubleSize(By,size*size)!=0){
    printf("Error, By must have shape ((scrnYPxls+1),(scrnYPxls+1)), be contiguous and float64\n");
    return NULL;
  }
  size=scrnXPxls+1;
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
  if(checkContigDoubleSize(screen,(scrnYPxls+1)*(scrnXPxls+1))!=0){
    printf("screen must have shape %d,%d be contiguous and float64\n",scrnYPxls+1,scrnXPxls+1);
    return NULL;
  }
  if(checkContigDoubleSize(randn,(scrnYPxls>scrnXPxls?scrnYPxls+1:scrnXPxls+1))!=0){
    printf("random array space must have shape %d becontiguous and float64\n",scrnYPxls>scrnXPxls?scrnYPxls+1:scrnXPxls+1);
    return NULL;
  }
  if(PyArray_Check(XstepObj)){
    if(checkContigDoubleSize((PyArrayObject*)XstepObj,scrnYPxls+1)!=0){
      printf("Xstep shoud be size scrnYPxls+1, contiguous, float64\n");
      return NULL;
    }
    XstepArr=(double*)(((PyArrayObject*)XstepObj)->data);
  }else{
    Xstep=(double)PyFloat_AsDouble(XstepObj);
  }
  if(PyArray_Check(YstepObj)){
    if(checkContigDoubleSize((PyArrayObject*)YstepObj,scrnXPxls+1)!=0){
      printf("Ystep shoud be size scrnXPxls+1, contiguous, float64\n");
      return NULL;
    }
    YstepArr=(double*)(((PyArrayObject*)YstepObj)->data);
  }else{
    Ystep=(double)PyFloat_AsDouble(YstepObj);
  }
  if(PyArray_Check(colOutputObj)){
    if(checkContigDoubleSize((PyArrayObject*)colOutputObj,(scrnYPxls+1)*maxColAdd)!=0){
      printf("colOutput should be size (scrnYPxls+1)*maxColAdd, contiguous, float64\n");
      return NULL;
    }
    colOutput=(double*)(((PyArrayObject*)colOutputObj)->data);
  }
  
  if(PyArray_Check(rowOutputObj)){
    if(checkContigDoubleSize((PyArrayObject*)rowOutputObj,(scrnXPxls+1)*maxRowAdd)!=0){
      printf("rowOutput should be size (scrnXPxls+1)*maxRowAdd, contiguous, float64\n");
      return NULL;
    }
    rowOutput=(double*)(((PyArrayObject*)colOutputObj)->data);
  }
  if(sendWholeScreen==0 && (rowOutput==NULL || colOutput==NULL)){
    printf("If not sending whole screen, row and col output must be defined\n");
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
  ss->maxColAdd=maxColAdd;
  ss->maxRowAdd=maxRowAdd;
  ss->colAdd=colAdd;
  ss->rowAdd=rowAdd;
  ss->Ay=(double*)Ay->data;
  ss->Ax=(double*)Ax->data;
  ss->AStartX=(double*)AStartX->data;
  ss->By=(double*)By->data;
  ss->Bx=(double*)Bx->data;
  ss->AStartY=(double*)AStartY->data;
  ss->screen=(double*)screen->data;
  ss->randn=(double*)randn->data;
  ss->XstepArr=XstepArr;
  ss->YstepArr=YstepArr;
  ss->Xstep=Xstep;
  ss->Ystep=Ystep;
  ss->ncStruct[0].stepSize=fabs(colAdd);
  ss->ncStruct[1].stepSize=fabs(rowAdd);
  ss->rowOutput=rowOutput;
  ss->colOutput=colOutput;
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
  case R0X:
    ss->r0x=PyFloat_AsDouble(obj);
    if(PyErr_Occurred()){
      printf("scrnmodule: Error extracting float value for r0x\n");
      return NULL;
    }
    break;
  case R0Y:
    ss->r0y=PyFloat_AsDouble(obj);
    if(PyErr_Occurred()){
      printf("scrnmodule: Error extracting float value for r0y\n");
      return NULL;
    }
    break;
  case XSTEP:
    if(PyArray_Check(obj)){
      if(checkContigDoubleSize((PyArrayObject*)obj,ss->scrnYPxls+1)!=0){
	printf("Xstep shoud be size scrnYPxls+1, contiguous, float64\n");
	return NULL;
      }
      ss->XstepArr=(double*)(((PyArrayObject*)obj)->data);
    }else{
      ss->Xstep=(double)PyFloat_AsDouble(obj);
      if(PyErr_Occurred()){
	printf("scrnmodule: Error extracting double value for xstep\n");
	return NULL;
      }

    }
    break;
  case YSTEP:
    if(PyArray_Check(obj)){
      if(checkContigDoubleSize((PyArrayObject*)obj,ss->scrnXPxls+1)!=0){
	printf("Ystep shoud be size scrnXPxls+1, contiguous, float64\n");
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
  Py_INCREF(Py_None);
  return Py_None;
}
    
/* define a methods table for this module */

static PyMethodDef scrn_methods[] = 	{
  {"run", py_run, METH_VARARGS}, 
  {"initialise", py_initialise, METH_VARARGS}, 
  {"update", py_update, METH_VARARGS}, 
  {"free", py_free, METH_VARARGS}, 
  {NULL, NULL} };


/* initialisation - register the methods with the Python interpreter */

void initscrn(void){
  PyObject *m;
  PyImport_AddModule("scrn");
  m=Py_InitModule("scrn", scrn_methods);
  import_array();
  ScrnError=PyErr_NewException("scrn.error",NULL,NULL);
  Py_INCREF(ScrnError);
  PyModule_AddObject(m,"error",ScrnError);
}

int main(int argc,char **argv){
  Py_SetProgramName(argv[0]);
  Py_Initialize();
  initscrn();
  return 0;
}
