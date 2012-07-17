//
// Urban Bitenc, July 2012
//
// Multi-threaded cubic interpolation
//
// This is used in aosim/cmod/src/interpmodule.c: gslCubSplineIntrp_UB 
//

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <string.h>  // necessary for memcpy

// data structure necessary to pass the parameters to a function executed in a thread:
typedef struct {
  int inN;           // number of input points
  int outN;          // number of output points
  int M;             // how many input lines are to be processed
  double* inX;       // vector of input X coordinates
  double* outX;      // vector of output X coordinates
  void* inData;      // input data - several lines of inY
  void* outData;     // output data - sevelar lines of outY
  int sx, sy;        // number of strides, needed for correct indexing
                     //        of the input or output matrix

  gsl_interp_accel *interpAcc; // interpolation accelerator
  double* y1;       // tmp vector to pass the data to the interpolator
} interp_data_t;

// function to perform the interpolation in a thread:
void interp_first_inDouble(  interp_data_t* arguments ); // interpolation in first dimension if the input  is Double
void interp_first_inFloat(   interp_data_t* arguments ); //        -''-      first    -''-   if     -''-      Float
void interp_second_outDouble(interp_data_t* arguments ); //        -''-      second   -''-   if the output is Double
void interp_second_outFloat( interp_data_t* arguments ); //        -''-      second   -''-   if     -''-      Float
