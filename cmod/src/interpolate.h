//
// Urban Bitenc, July 2012
//       improved Feb 1st 2013  (see comment below)
//
// Multi-threaded cubic interpolation
//
// This is used in aosim/cmod/src/interpmodule.c: gslCubSplineIntrp
//
// July 2012: only works if input and output matrices are NOT identical
//
// Feb 2013: works also if input and output matrices are the same, e.g.
//           cmod.interp.gslCubSplineInterp(mxout[0:5, 0:5], yin, xin, yout, xout, mxout, nThreads)
//           Main changes: added ...


#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <string.h>  // necessary for memcpy

// data structure necessary to pass the parameters to a function executed in a thread:
typedef struct {
  int M;             // how many input lines/columns are to be processed by this thread
  int inN;           // number of input points
  int outN;          // number of output points
  double* inX;       // vector of input X coordinates
  double* outX;      // vector of output X coordinates
  void* inData;      // input data - several lines of inY
  void* outData;     // output data - sevelar lines of outY

  int sx, sy;        // (i) strides of the input  array, needed for indexing in interp_first_*
                     // (ii)strides of the output array, needed for indexing in interp_second_*
  int sytmp;         // y stride of sytmp, needed in interp_first_*
  
                     // The following two variables are used only in the interp_first_inDouble etc.
                     // functions and could hence be made local. However, by adding them to params
                     // and initialising them in the main function enables you to re-use them which
                     // makes  the program about 1.5% faster
                     // (for a 100x100 input to 991x991 output array).

  gsl_interp_accel* interpAcc; // interpolation accelerator
  //  double*           y1;        // tmp vector to pass the data to the interpolator
} interp_data_t;

// function to perform the interpolation in a thread:
void interp_first_inDouble(  interp_data_t* arguments ); // interp. in Y if input  is Double
void interp_first_inFloat(   interp_data_t* arguments ); //  -''-      Y    -''-      Float
void interp_second_outDouble(interp_data_t* arguments ); //  -''-      X if output is Double
void interp_second_outFloat( interp_data_t* arguments ); //  -''-      X    -''-      Float
