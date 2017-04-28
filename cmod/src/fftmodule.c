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
#include <stdio.h>
#include <fftw3.h>
#include <numpy/arrayobject.h>
//A module for fft related things.
//Note - fftw faster for out of place transforms, particularly when threading.
//Appears to be faster than acml.

int checkCFloatContigArr(PyArrayObject *arr){
  //checks whether an array is contiguous or not, and whether it is float..
  int rtval=0;//return 0 on success.
  if(arr->descr->type_num!=NPY_CFLOAT)
    rtval=1;
  if(PyArray_ISCONTIGUOUS(arr)==0)//contiguous?
    rtval=1;
  return rtval;
}


static PyObject* InitialiseThreading(PyObject *self,PyObject *args){
  //create an array from existing memory, and return new object.  The new array can use all or part of the old one, and can be different data type.
  static int initialised=0;
  int nthreads;
  if(!PyArg_ParseTuple(args,"i",&nthreads)){
    printf("Usage - nthreads\n");
    return NULL;
  }
  if(initialised==0){
    initialised=1;
    fftwf_init_threads();
    fftwf_plan_with_nthreads(nthreads);
  }else{
    //printf("fftmodule InitialiseThreading call ignored - already initialised\n");
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* Plan(PyObject *self,PyObject *args){
  fftwf_plan p;
  PyArrayObject *in,*out;
  if(!PyArg_ParseTuple(args,"O!O!",&PyArray_Type,&in,&PyArray_Type,&out)){
    printf("Usage - in array, out array\n");
    return NULL;
  }
  if(in->nd!=2 || out->nd!=2 || in->dimensions[0]!=out->dimensions[0] || in->dimensions[1]!=out->dimensions[1]){
    printf("fft2Plan - dimensions must be same for input and output, 2D array\n");
    return NULL;
  }
  if(checkCFloatContigArr(in)!=0 || checkCFloatContigArr(out)!=0){
    printf("input and output must be complex float contiguous\n");
    return NULL;
  }
  if(((long)(in->data))&0xf || ((long)(out->data))&0xf){
    printf("Warning - cmod.fftmodule - data may not be aligned for SIMD operation\n");
  }
  p=fftwf_plan_dft_2d((int)(in->dimensions[0]),(int)(in->dimensions[1]),(fftwf_complex*)(in->data),(fftwf_complex*)(out->data),FFTW_FORWARD,FFTW_DESTROY_INPUT | FFTW_MEASURE);//the fft is allowed to destoy the input.
  return Py_BuildValue("l",(long)p);
}

static PyObject* ExecutePlan(PyObject *self,PyObject *args){
  //execute an FFT from the plan.
  fftwf_plan p;
  if(!PyArg_ParseTuple(args,"l",&p)){
    printf("Usage - fft plan (returned from Plan())\n");
    return NULL;
  }
  fftwf_execute(p);
  Py_INCREF(Py_None);
  return Py_None;
}


static PyObject *FreePlan(PyObject *self,PyObject *args){
  fftwf_plan p;
  if(!PyArg_ParseTuple(args,"l",&p)){
    printf("Usage - fft plan (returned from Plan())\n");
    return NULL;
  }
  fftwf_destroy_plan(p);
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject *CleanUp(PyObject *self,PyObject *args){
  //should be called at the very end.
  static int initialised=0;
  if(initialised==0){
    initialised=1;
    fftwf_cleanup_threads();
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyMethodDef fftMethods[] = {
  {"InitialiseThreading",InitialiseThreading,METH_VARARGS,
   "Initialise the threading (arg=nthreads)"},
  {"Plan",Plan,METH_VARARGS,"Plan the FFT"},
  {"ExecutePlan",ExecutePlan,METH_VARARGS,"Do the FFT"},
  {"FreePlan",FreePlan,METH_VARARGS,"Free the plan"},
  {"CleanUp",CleanUp,METH_VARARGS,"Clean up at the end"},
  {NULL,NULL,0,NULL}
};
void initfft(void){
  //  PyObject *m; // commented out UB, 2012Aug08
  PyImport_AddModule("fft");
  (void)Py_InitModule("fft",fftMethods); // was: m = Py_... (UB, 2012Aug08); changed to remove warning "m set but not used"
  import_array();
}
int main(int argc,char *argv[]){
  printf("fftmodule main\n");
  Py_SetProgramName(argv[0]);
  
  // Initialize the Python interpreter.  Required. 
  Py_Initialize();
  
  // Add a static module 
  initfft();
  return 0;

}
