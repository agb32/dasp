/* Module numfftw3 : interface to functions in FFTW3 (single precision)
 * Wrapper between FFTW3 and NUMERIC (not NUMARRAY) arrays*/

%module numfftw3

%include typemaps.i

%{
#include <stdio.h>
#include <fftw3.h>
#include "Python.h"
#include "Numeric/arrayobject.h"
%}

%init %{
      import_array();
%}

enum fftw_r2r_kind_do_not_use_me {
     FFTW_R2HC=0, FFTW_HC2R=1, FFTW_DHT=2,
     FFTW_REDFT00=3, FFTW_REDFT01=4, FFTW_REDFT10=5, FFTW_REDFT11=6,
     FFTW_RODFT00=7, FFTW_RODFT01=8, FFTW_RODFT10=9, FFTW_RODFT11=10
};

#define FFTW_FORWARD (-1)
#define FFTW_BACKWARD (+1)

/* documented flags */
#define FFTW_MEASURE (0U)
#define FFTW_DESTROY_INPUT (1U << 0)
#define FFTW_UNALIGNED (1U << 1)
#define FFTW_CONSERVE_MEMORY (1U << 2)
#define FFTW_EXHAUSTIVE (1U << 3) /* NO_EXHAUSTIVE is default */
#define FFTW_PRESERVE_INPUT (1U << 4) /* cancels FFTW_DESTROY_INPUT */
#define FFTW_PATIENT (1U << 5) /* IMPATIENT is default */
#define FFTW_ESTIMATE (1U << 6)

/* undocumented beyond-guru flags */
#define FFTW_ESTIMATE_PATIENT (1U << 7)
#define FFTW_BELIEVE_PCOST (1U << 8)
#define FFTW_DFT_R2HC_ICKY (1U << 9)
#define FFTW_NONTHREADED_ICKY (1U << 10)
#define FFTW_NO_BUFFERING (1U << 11)
#define FFTW_NO_INDIRECT_OP (1U << 12)
#define FFTW_ALLOW_LARGE_GENERIC (1U << 13) /* NO_LARGE_GENERIC is default */
#define FFTW_NO_RANK_SPLITS (1U << 14)
#define FFTW_NO_VRANK_SPLITS (1U << 15)
#define FFTW_NO_VRECURSE (1U << 16)
#define FFTW_NO_SIMD (1U << 17)


/* Typemaps */
/* Conversion from Numeric Complex64 arrays to fftw_complex* data type */
%typemap(python,in) fftw_complex * {
    /* Used to convert Numeric array to C pointer*/
    PyArrayObject* arr;

    /* conversion from source to Numeric array */
    if (!PyArray_Check($input))
    {
        PyErr_SetString(PyExc_TypeError,"expected a Numeric array");
        return(NULL);
    }
   
    arr=(PyArrayObject*) PyArray_ContiguousFromObject($input,PyArray_CDOUBLE,0,0);
    if (arr==NULL)
    {
        PyErr_SetString(PyExc_MemoryError,"problem in PyArray_ContiguousFromObject (Complex64)");
        return NULL;
    }
    
    /* we adjust the pointers */
    $1=(fftw_complex *)arr->data;
    $input=(PyObject *)arr;
    Py_XDECREF(arr);    
}

/* Conversion from Numeric Float64 arrays to double* data type */
%typemap(python,in) double * {
    /* Used to convert Numeric array to C pointer*/
    PyArrayObject* arr;
    
    /* conversion from source to Numeric array */
    if (!PyArray_Check($input))
    {
        PyErr_SetString(PyExc_TypeError,"expected a Numeric array");
        return(NULL);
    }
    
    /* Make sure input is a Numeric array */
    arr=(PyArrayObject*) PyArray_ContiguousFromObject($input,PyArray_DOUBLE,0,0);
    if (arr==NULL)
    {
        PyErr_SetString(PyExc_MemoryError,"problem in PyArray_ContiguousFromObject (Float64)");
        return NULL;
    }
     
    /* we adjust the pointers */
    $1=(double *)arr->data;
    $input=(PyObject *)arr;
    Py_XDECREF(arr);
}

/* Conversion from numarray Complex32 arrays to fftwf_complex* data type */
%typemap(python,in) fftwf_complex * {
    /* Used to convert Numeric array to C pointer*/
    PyArrayObject* arr;
    
    /* conversion from source to Numeric array */
    if (!PyArray_Check($input))
    {
        PyErr_SetString(PyExc_TypeError,"expected a Numeric array");
        return(NULL);
    }

    arr=(PyArrayObject*) PyArray_ContiguousFromObject($input,PyArray_CFLOAT,0,0);
    if (arr==NULL)
    {
        PyErr_SetString(PyExc_MemoryError,"problem in PyArray_ContiguousFromObject (Complex32)");
        return NULL;
    }
    
    /* we adjust the pointers */
    $1=(fftwf_complex *)arr->data;
    $input=(PyObject *)arr;
    Py_XDECREF(arr);
}

/* Conversion from numarray Float32 arrays to float* data type */
%typemap(python,in) float * {
    /* Used to convert Numeric array to C pointer*/
    PyArrayObject* arr;
    
    /* conversion from source to Numeric array */
    if (!PyArray_Check($input))
    {
        PyErr_SetString(PyExc_TypeError,"expected a Numeric array");
        return(NULL);
    }

    arr=(PyArrayObject*) PyArray_ContiguousFromObject($input,PyArray_FLOAT,0,0);
    if (arr==NULL)
    {
        PyErr_SetString(PyExc_MemoryError,"problem in PyArray_ContiguousFromObject (Float32)");
        return NULL;
    }
    
    /* we adjust the pointers */
    $1=(float *)arr->data;
    $input=(PyObject *)arr;
    Py_XDECREF(arr);
}

/* Conversion from python tuple to int* data type */
%typemap(python,in) int * {
    /* Used to convert tuple to C pointer*/
    PyArrayObject* t;

    /* Check we have a tuple as input */
    if (!PyTuple_Check($input))
    {
        PyErr_SetString(PyExc_TypeError,"expected a tuple");
        return(NULL);
    }
    
    t=(PyArrayObject*) PyArray_ContiguousFromObject($input,PyArray_INT,0,0);
    if (t==NULL)
    {
        PyErr_SetString(PyExc_MemoryError,"problem in PyArray_ContiguousFromObject (tuple)");
        return NULL;
    }
    
    /* we adjust the pointers */
    $1=(int *)t->data;
    $input=(PyObject *)t;
    Py_XDECREF(t);   
}

/* Conversion from python string to FILE* data type */
%typemap(python,in) FILE * {
  PyObject* str;
  char *filename;
  FILE* f;

  /* check for string */
  if (!PyString_Check($input))
  {
      PyErr_SetString(PyExc_TypeError,"expected a string");
      return(NULL);
  }

  /* conversion from python object to string */
  filename=PyString_AsString($input);
  if (filename==NULL)
  {
      PyErr_SetString(PyExc_TypeError,"String copy problem");
      return(NULL);
  }

  /* opening of the file */
  f=fopen(filename,"a");
  if (f==NULL)
  {
      PyErr_SetString(PyExc_IOError,"File opening problem");
      return(NULL);
  }
  /* we adjust the pointers */
  $1=f;
  /*$input=(PyObject *)tuple;*/
}


%typemap(python,freearg) FILE * {
    if ($1!=NULL) fclose($1);
}


%include fftw64.i
