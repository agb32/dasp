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
/*
Python to perform mmx using amd lapack.

Routine to perform SVD of a matrix using lapack.  Hopefully, similar to the scalapack version - I have a hunch that intel acml may be multi-threaded, so may be faster than scalapack on a single computer (eg 8 core mac).


In this version, reads an input file, and does svd of it.
The size of the input/output can vary, though this is still square.  Command line arguments are n

In later versions, will have Commandline arguments:
ncols, nrows, input data filename, output filename
The input and outputs are FITS files.

Perform the SVD, then write the results.
*/


#include <Python.h>
#include <numpy/arrayobject.h>
#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include <errno.h>
#include <stdlib.h>
#include <unistd.h>
#include <limits.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <assert.h>
#include <math.h>
#include <sys/mman.h>
#include <acml.h>
typedef int ACML_INT;

static PyObject *acmlError;

static PyObject* ludecomp(PyObject *self,PyObject *args){
  PyArrayObject *Aarray,*ipiv;
  ACML_INT info=0;
  ACML_INT n,m;
  if(!PyArg_ParseTuple(args,"O!O!",&PyArray_Type,&Aarray,&PyArray_Type,&ipiv)){
    printf("Usage: A, ipiv (dimension equal to smallest dimension of A).\n");
    return NULL;
  }
  if(Aarray->nd!=2){
    printf("A must be 2D\n");
    return NULL;
  }
  n=Aarray->dimensions[0];
  m=Aarray->dimensions[1];
  if(ipiv->nd!=1 || ipiv->dimensions[0]<(m<n?m:n)){
    printf("ipiv should be 1D with dimensions equal to min of dimension of A\n");
    return NULL;
  }
  if(!(Aarray->descr->type_num==NPY_FLOAT || Aarray->descr->type_num==NPY_DOUBLE)){
    printf("A must be float or double\n");
    return NULL;
  }
  if(sizeof(ACML_INT)!=ipiv->descr->elsize || ipiv->descr->kind!='i'){
    printf("ipiv should be integer type size %ld (has size %d, kind %c\n",sizeof(ACML_INT),ipiv->descr->elsize,ipiv->descr->kind);
    return NULL;
  }
  if(!PyArray_ISCONTIGUOUS(Aarray)){
    printf("A must be contiguous\n");
    return NULL;
  }
  if(!PyArray_ISCONTIGUOUS(ipiv)){
    printf("ipiv must be contiguous\n");
    return NULL;
  }
  info=0;
  Py_BEGIN_ALLOW_THREADS;
  if(Aarray->descr->type_num==NPY_FLOAT){
    sgetrf(m,n,(float*)Aarray->data,m,(ACML_INT*)ipiv->data,&info);
  }else{
    dgetrf(m,n,(double*)Aarray->data,m,(ACML_INT*)ipiv->data,&info);
  }
  Py_END_ALLOW_THREADS;
  if(info!=0){
    printf("Error in trf: info=%ld\n",(long int)info);
    return NULL;
  }
  return Py_BuildValue("l",(long)info);
}

static PyObject *luinv(PyObject *self,PyObject *args){
  PyArrayObject *Aarray,*ipiv,*workArray;
  PyObject *workObj;
  ACML_INT info=0;
  ACML_INT n,m,lwork;
  int checkWorkSize=0;
  char *workarr=NULL;
  if(!PyArg_ParseTuple(args,"O!O!O",&PyArray_Type,&Aarray,&PyArray_Type,&ipiv,&workObj)){
    printf("Usage: A containing results from ludecomp, ipiv (dimension equal to smallest dimension of A, type long, containing results from ludecomp), work array or None (in which case, this returns the ideal size for the work array).\n");
    return NULL;
  }
  if(Aarray->nd!=2){
    printf("A must be 2D\n");
    return NULL;
  }
  n=Aarray->dimensions[0];
  m=Aarray->dimensions[1];
  if(m!=n){
    printf("A must be square\n");
    return NULL;
  }
  if(ipiv->nd!=1 || ipiv->dimensions[0]<(m<n?m:n)){
    printf("ipiv should be 1D with dimensions equal to min of dimension of A\n");
    return NULL;
  }
  if(!(Aarray->descr->type_num==NPY_FLOAT || Aarray->descr->type_num==NPY_DOUBLE)){
    printf("A must be float or double\n");
    return NULL;
  }
  if(sizeof(ACML_INT)!=ipiv->descr->elsize || ipiv->descr->kind!='i'){
    printf("ipiv should be integer type size %ld\n",sizeof(ACML_INT));
    return NULL;
  }
  if(!PyArray_ISCONTIGUOUS(Aarray)){
    printf("A must be contiguous\n");
    return NULL;
  }
  if(!PyArray_ISCONTIGUOUS(ipiv)){
    printf("ipiv must be contiguous\n");
    return NULL;
  }

  if(Py_None==workObj){
    //must only be doing the work size check...
    checkWorkSize=1;
    lwork=-1;
    if((workarr=malloc(sizeof(double)))==NULL){
      printf("unable to allocate workarr\n");
      return NULL;
    }
  }else if(PyArray_Check(workObj)){
    workArray=(PyArrayObject*)workObj;
    workarr=workArray->data;
    if(workArray->nd!=1 || workArray->descr->type_num!=Aarray->descr->type_num || !PyArray_ISCONTIGUOUS(workArray)){
      printf("work array should be 1D, of same type as A and contiguous\n");
      return NULL;
    }
    lwork=workArray->dimensions[0];
  }else{
    printf("work must be array or None\n");
    return NULL;
  }



  info=0;
  Py_BEGIN_ALLOW_THREADS;
  if(Aarray->descr->type_num==NPY_FLOAT){
    sgetri(n,(float*)Aarray->data,n,(ACML_INT*)ipiv->data,&info);
  }else{
    dgetri(n,(double*)Aarray->data,n,(ACML_INT*)ipiv->data,&info);

  }
  Py_END_ALLOW_THREADS;

  if(checkWorkSize){
    if(Aarray->descr->type_num==NPY_FLOAT){
      lwork=(ACML_INT)(*((float*)workarr));
    }else{
      lwork=(ACML_INT)(*((double*)workarr));
    }
    free(workarr);
  }else{
    lwork=0;
  }
  if(info!=0){
    printf("Error in luinv: info=%ld\n",(long int)info);
    return NULL;
  }
  return Py_BuildValue("l",(long)lwork);


}

/*
static PyObject* svd(PyObject *self,PyObject *args){
  //Perform an SVD.
  ACML_INT info=0;
  ACML_INT n,lwork;
  PyArrayObject *Aarray,*Uarray,*VTarray,*workArray,*evalArray;
  PyObject *Uobj,*VTobj,*workObj,*evalObj;
  int usesdd=1;
  char jobu='N',jobvt='N',jobz='N';
  int doU=0,doVT=0,overwriteUA=0,overwriteVTA=0;
  ACML_INT *iwork=NULL;
  char *uarr=NULL,*vtarr=NULL,*workarr=NULL,*evalarr=NULL,*aarr=NULL;
  int checkWorkSize=0;

  if(!PyArg_ParseTuple(args,"O!OOOO|i",&PyArray_Type,&Aarray,&Uobj,&evalObj,&VTobj,&workObj,&usesdd)){
    printf("Usage: A, U or None, eval array, VT, work array or None, usesdd (optional, default=1)\nIf U==None, then A will be overwritten.  If work array==None, function does nothing except computes the optimal size for work and returns it.  usesdd can be set to choose between DGESVD and DGESDD (or the corresponding floating point versions, SGESVD and SGESDD).\n");
    printf("Sizeof(lwork) %ld\n",sizeof(lwork));
    return NULL;
  }
  if(Aarray->nd!=2 || Aarray->dimensions[0]!=Aarray->dimensions[1]){
    printf("A must be square\n");
    return NULL;
  }
  if(!(Aarray->descr->type_num==NPY_FLOAT || Aarray->descr->type_num==NPY_DOUBLE)){
    printf("A must be float or double\n");
    return NULL;
  }
  if(!PyArray_ISCONTIGUOUS(Aarray)){
    printf("A must be contiguous\n");
    return NULL;
  }
  aarr=Aarray->data;
  if(Py_None==Uobj){
    doU=0;
    Uarray=NULL;
  }else if(PyArray_Check(Uobj)){
    Uarray=(PyArrayObject*)Uobj;
    doU=1;
    if(Uarray->nd!=2 || Uarray->dimensions[0]!=Uarray->dimensions[1] || Uarray->dimensions[0]!=Aarray->dimensions[0] || Uarray->descr->type_num!=Aarray->descr->type_num || !PyArray_ISCONTIGUOUS(Uarray)){
      printf("U must be square, same type as A and contiguous\n");
      return NULL;
    }
    if(Uarray->data==Aarray->data){
      overwriteUA=1;//overwrite A.
    }
  }else{
    printf("U must be array or None\n");
    return NULL;
  }
  if(Py_None==VTobj){
    doVT=0;
    VTarray=NULL;
  }else if(PyArray_Check(VTobj)){
    VTarray=(PyArrayObject*)VTobj;
    doVT=1;
    if(VTarray->nd!=2 || VTarray->dimensions[0]!=VTarray->dimensions[1] || VTarray->dimensions[0]!=Aarray->dimensions[0] || VTarray->descr->type_num!=Aarray->descr->type_num || !PyArray_ISCONTIGUOUS(VTarray)){
      printf("VT must be square, same type as A and contiguous\n");
      return NULL;
    }
    if(VTarray->data==Aarray->data){
      overwriteVTA=1;//overwrite A.
      if(overwriteUA==1){
	printf("Error - cannot overwrite A with both U and VT\n");
	return NULL;
      }
    }
  }else{
    printf("VT must be array or None\n");
    return NULL;
  }
  if(Py_None==evalObj){
    //must only be doing the work size check...
    if(!(Py_None==workObj)){
      printf("Must supply eval array\n");
      return NULL;
    }
    evalArray=NULL;
    evalarr=NULL;
  }else if(PyArray_Check(evalObj)){
    evalArray=(PyArrayObject*)evalObj;
    evalarr=evalArray->data;
    if(evalArray->nd!=1 || evalArray->dimensions[0]!=Aarray->dimensions[0] || evalArray->descr->type_num!=Aarray->descr->type_num || !PyArray_ISCONTIGUOUS(evalArray)){
      printf("eval array should be 1D, size of A.shape[0], of same type as A and contiguous\n");
      return NULL;
    }
  }else{
    printf("eval must be array or None\n");
    return NULL;
  }
  if(Py_None==workObj){
    //must only be doing the work size check...
    checkWorkSize=1;
    lwork=-1;
    if((workarr=malloc(sizeof(double)))==NULL){
      printf("unable to allocate workarr\n");
      return NULL;
    }
  }else if(PyArray_Check(workObj)){
    workArray=(PyArrayObject*)workObj;
    workarr=workArray->data;
    if(workArray->nd!=1 || workArray->descr->type_num!=Aarray->descr->type_num || !PyArray_ISCONTIGUOUS(workArray)){
      printf("work array should be 1D, of same type as A and contiguous\n");
      return NULL;
    }
    lwork=workArray->dimensions[0];
  }else{
    printf("work must be array or None\n");
    return NULL;
  }
  
  if(doU){
    if(overwriteUA){
      uarr=NULL;
      if(usesdd){
	jobz='O';
      }else{
	jobu='O';
      }
    }else{//writing U into U.
      uarr=Uarray->data;
      if(usesdd){
	jobz='S';
      }else{
	jobu='S';
      }
    }
  }else{
    uarr=NULL;
    jobu='N';
    jobz='N';
  }
  if(doVT){
    if(overwriteVTA){
      if(usesdd){
	if(jobz=='S'){//have requested inplace for U
	  jobz='O';//here, actually, U will go onto A, and vt will go onto U.
	  vtarr=Uarray->data;
	  
	}else{
	  jobz='N';
	  vtarr=NULL;
	}
      }else{
	jobvt='O';
	vtarr=NULL;
      }
    }else{//not overwriting vt...
      if(usesdd){
	//not overwriting VT.  So, if jobz already S, okay, if O okay, if N okay.
	if(jobz=='N')
	  vtarr=NULL;
	else
	  vtarr=VTarray->data;
      }else{
	jobvt='S';
	vtarr=VTarray->data;
      }
    }
  }else{
    jobvt='N';
    jobz='N';
    vtarr=NULL;
  }
  n=Aarray->dimensions[0];
  info=0;
  
  if(checkWorkSize==0){
    if(usesdd){
      if((iwork=(ACML_INT*)malloc(sizeof(ACML_INT)*8*n))==NULL){
	printf("unable to allocate iwork\n");
	return NULL;
      }
    }
  }
  //printf("%c %c %c\n",jobu,jobvt,jobz);
  Py_BEGIN_ALLOW_THREADS;
  if(usesdd){
    printf("jobz: %c",jobz);
    if(Aarray->descr->type_num==NPY_FLOAT){
      sgesdd(&jobz,&n,&n,(float*)aarr,&n,(float*)evalarr,(float*)uarr,&n,(float*)vtarr,&n,(float*)workarr,&lwork,iwork,&info);
    }else{
      dgesdd(&jobz,&n,&n,(double*)aarr,&n,(double*)evalarr,(double*)uarr,&n,(double*)vtarr,&n,(double*)workarr,&lwork,iwork,&info);
    }
  }else{
    printf("jobu %c jobvt %c\n",jobu,jobvt);
    if(Aarray->descr->type_num==NPY_FLOAT){
      sgesvd(&jobu,&jobvt,&n,&n,(float*)aarr,&n,(float*)evalarr,(float*)uarr,&n,(float*)vtarr,&n,(float*)workarr,&lwork,&info);
    }else{
      dgesvd(&jobu,&jobvt,&n,&n,(double*)aarr,&n,(double*)evalarr,(double*)uarr,&n,(double*)vtarr,&n,(double*)workarr,&lwork,&info);
    }
  }
  Py_END_ALLOW_THREADS;
  
  if(checkWorkSize){
    if(Aarray->descr->type_num==NPY_FLOAT){
      lwork=(ACML_INT)(*((float*)workarr));
    }else{
      lwork=(ACML_INT)(*((double*)workarr));
    }
    free(workarr);
  }else{
    if(usesdd)
      free(iwork);
    lwork=0;
    
  }
  if(info!=0){
    printf("Error in svd: info=%ld\n",(long int)info);
    return NULL;
  }
  return Py_BuildValue("l",(long)lwork);

}
*/
/*
In FORTRAN - elements in a column are contiguous.
In C, elements in a row are contiguous.
Fortran description:
M
    (input) INTEGER
    The number of matrix rows to be operated on. 
N
    (input) INTEGER
    The number of matrix columns to be operated on. 
A
    (input/output) TYPE DEPENDS ON ROUTINE, array of dimension (LDA,N)
    A pointer to the beginning of the (sub)array to be sent. 
LDA
    (input) INTEGER
    The distance between two elements in matrix row, ==M usually.

C description:
M - number of columns
N - number of rows
LDA - usually == M
Basically, the data is transposed.
*/

static PyObject* gemm(PyObject *self,PyObject *args){
  //perform matrix multiply
  PyArrayObject *Aarray,*Barray,*Carray;
  double dalpha=1.,dbeta=0.;
  float falpha=1.,fbeta=0.;
  char transA='N',transB='N';
  ACML_INT lda,ldb,k,m,n,ldc;
  long wsize;

  if(!PyArg_ParseTuple(args,"O!O!O!|dd",&PyArray_Type,&Aarray,&PyArray_Type,&Barray,&PyArray_Type,&Carray,&dalpha,&dbeta)){
    printf("Usage: A, B, C, [alpha] [beta]\nResult is C = alpha * op(A) * op(B) + beta*C where op corresponds to a transpose if A or B are transposed, and alpha and beta are constants.\nIf given only A, B and C, then computes matrix multiply of A and B, putting result in C.");
    return NULL;
  }
  falpha=(float)dalpha;
  fbeta=(float)dbeta;
  if(Aarray->nd!=2 || Barray->nd!=2 || Carray->nd!=2){
    printf("arrays must be 2D\n");
    return NULL;
  }
  if(Aarray->descr->type_num!=NPY_FLOAT && Aarray->descr->type_num!=NPY_DOUBLE){
    printf("Must be float or double arrays\n");
    return NULL;
  }
  if(Barray->descr->type_num!=Aarray->descr->type_num || Carray->descr->type_num!=Aarray->descr->type_num){
    printf("Arrays must be of same type\n");
    return NULL;
  }
  /*if(!(Carray->flags&NPY_FORTRAN)){
    printf("Output array must be FORTRAN contiguous (use order='F' in call to numpy.zeros) or change the C.strides.\n");
    return NULL;
    }*/
  //now work out whether A, B are transposed or not, and check contiguous.
  wsize=Aarray->descr->type_num==NPY_FLOAT?4:8;
  if(Aarray->strides[1]==wsize)// && Aarray->strides[0]==wsize*Aarray->dimensions[1])
    transA='T';
  else if(Aarray->strides[0]==wsize)// && Aarray->strides[1]==wsize*Aarray->dimensions[0])
    transA='N';
  else{
    printf("A must be contiguous or transpose of contiguous along the stepping dimension\n");
    return NULL;
  }

  if(Barray->strides[1]==wsize)// && Barray->strides[0]==wsize*Barray->dimensions[1])
    transB='T';
  else if(Barray->strides[0]==wsize)// && Barray->strides[1]==wsize*Barray->dimensions[0])
    transB='N';
  else{
    printf("B must be contiguous or transpose of contiguous\n");
    return NULL;
  }
  if(Carray->strides[0]!=wsize){
    printf("C must be contiguous along first dimension - use order='F' in call to numpy.zeros, or change c.strides\n");
    return NULL;
  }

  //now check matrix sizes are okay.
  m=(ACML_INT)Aarray->dimensions[0];
  if(Carray->dimensions[0]!=m){
    printf("Rows of A and C not in agreement\n");
    return NULL;
  }
  n=(ACML_INT)Barray->dimensions[1];
  if(Carray->dimensions[1]!=n){
    printf("Cols of B and C not in agreement\n");
    return NULL;
  }
  k=(ACML_INT)Aarray->dimensions[1];
  if(Barray->dimensions[0]!=k){
    printf("Rows of B not equal to cols of A\n");
    return NULL;
  }
  if(transA=='N'){
    //lda=m;
    lda=Aarray->strides[1]/wsize;
  }else{
    //lda=k;
    lda=Aarray->strides[0]/wsize;
  }
  if(transB=='N'){
    //ldb=k;
    ldb=Barray->strides[1]/wsize;
  }else{
    //ldb=n;
    ldb=Barray->strides[0]/wsize;
  }
  ldc=(ACML_INT)Carray->strides[1]/wsize;
  printf("gemm values: %d %d %d %d %d %d\n",(int)m,(int)n,(int)k,(int)lda,(int)ldb,(int)ldc);
  if(Aarray->descr->type_num==NPY_FLOAT){
    sgemm(transA,transB,m,n,k,falpha,(float*)Aarray->data,lda,(float*)Barray->data,ldb,fbeta,(float*)Carray->data,ldc);
  }else{
    dgemm(transA,transB,m,n,k,dalpha,(double*)Aarray->data,lda,(double*)Barray->data,ldb,dbeta,(double*)Carray->data,ldc);
  }
  return Py_BuildValue("l",0);
}



static PyMethodDef acmlMethods[] = {
  //  {"svd",svd, METH_VARARGS,"Do SVD of a square array."},
  {"gemm",gemm, METH_VARARGS,"Do matrix-matrix multiply."},
  {"ludecomp",ludecomp,METH_VARARGS,"Do LU decomposition (destroys input)"},
  {"luinv",luinv,METH_VARARGS,"Use results from ludecomp to do inversion"},
  {NULL, NULL, 0, NULL}        /* Sentinel */
};
//PyMODINIT_FUNC 
void initacml(void)
{
  PyObject *m;
  PyImport_AddModule("acml");
  m=Py_InitModule("acml", acmlMethods);
  import_array();
  acmlError = PyErr_NewException("acml.error", NULL, NULL);
  Py_INCREF(acmlError);
  PyModule_AddObject(m, "error", acmlError);
}


int
main(int argc, char *argv[])
{
    // Pass argv[0] to the Python interpreter 
  Py_SetProgramName(argv[0]);
  
  // Initialize the Python interpreter.  Required. 
  Py_Initialize();
  
  // Add a static module 
  initacml();
  return 0;
}

