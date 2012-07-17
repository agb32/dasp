// 
// Urban Bitenc, 2012 July 13
// 
// Functions called from interpmodule.c :: gclCubSplineIntrp_UB
// 
// 
// 
// Some comments:
//  - params->y1 and params->interpAcc can easily be made local (declared and used only in the
//    functions below). However, by adding them to params and initialising them in the main function
//    the program is about 1.5% faster (for a 100x100 input to 991x991 output array), because in the
//    main program they are initialised only once. If you initialise them in the functions below, they
//    are initialised twice - once with interp_first... and once with interp_second...
//
//  - There are four functions used: first_inFloat, first_inDouble, second_outFloat, second_outDouble.
//       It would be easy to write on function with two if clauses instead of 4 functions with no
//       if clauses. I prefer 4 functions with no if clauses because of better readability.
//
//  - The only difference between interp_first_inDouble and interp_first_inFloat is this line:
//       params->y1[i] = (double)(((double*)(par...
//       params->y1[i] = (double)(((float* )(par...
//    I wrote comments for only one of them.
//
//  - The only difference between interp_second_inDouble and interp_second_inFloat is this line:
//    	 ((double*)(params->outData))[j*(params->sy)+i*(params->sx)]=(double)gsl_spline_e...
//       ((float* )(params->outData))[j*(params->sy)+i*(params->sx)]=(float )gsl_spline_e...
//
//
//

#include "interpolate.h"
#include <stdio.h>

// interpolation in first dimension if the input is Double
void interp_first_inDouble(  interp_data_t* params )
{
  int i, j;              // for looping
  gsl_spline *interpObj; // for interpolation

  // Here you define the interpolation method used:
  interpObj=gsl_spline_alloc(gsl_interp_cspline, (size_t)(params->inN));

  // For each line of array:
  for (j=0; j<(params->M); ++j){
    // Move the data into the temporary vector y1 to pass it to the interpolation:
    for (i=0; i<(params->inN); ++i)
      params->y1[i] = (double)(((double*)(params->inData))[ i*(params->sy) + j*(params->sx) ]);
    // Prepare the interpolation function:
    gsl_spline_init(interpObj, params->inX, params->y1,(size_t)(params->inN));
    // Evaluate the interpolation function and save the result in the intermediate array:
    for (i=0; i<(params->outN); ++i)
      ((double*)(params->outData))[j+i*(params->sy)] = gsl_spline_eval(interpObj, (params->outX)[i], params->interpAcc);
  }
  // Tidy up:
  gsl_spline_free(interpObj);
}

// interpolation in first dimension if the input is Float
void interp_first_inFloat(  interp_data_t* params )
{
  int i, j;
  gsl_spline *interpObj;

  interpObj=gsl_spline_alloc(gsl_interp_cspline, (size_t)(params->inN));

  for (j=0; j<(params->M); ++j){
    for (i=0; i<(params->inN); ++i)
      params->y1[i] = (double)(((float*)(params->inData))[ i*(params->sy) + j*(params->sx) ]);
    gsl_spline_init(interpObj, params->inX, params->y1,(size_t)(params->inN));
    for (i=0; i<(params->outN); ++i)
      ((double*)(params->outData))[j+i*(params->sy)] = gsl_spline_eval(interpObj, (params->outX)[i], params->interpAcc);

  }
  gsl_spline_free(interpObj);
}

// interpolation in second dimension if the input is Double
void interp_second_outDouble(  interp_data_t* params )
{
  int i, j;              // for looping
  gsl_spline *interpObj; // for interpolation

  // Here you define the interpolation method used:
  interpObj=gsl_spline_alloc(gsl_interp_cspline, (size_t)(params->inN));

  // For each line of array:
  for (j=0; j<(params->M); ++j){
    // Move the data into the temporary vector y1 to pass it to the interpolation:
    memcpy(params->y1, &((double*)(params->inData))[j*(params->inN)], (params->inN)*sizeof(double));
    // Prepare the interpolation function:
    gsl_spline_init(interpObj, params->inX, params->y1,(size_t)(params->inN));
    // Evaluate the interpolation function and save the result in the output array:
    for (i=0; i<(params->outN); ++i)
      ((double*)(params->outData))[j*(params->sy)+i*(params->sx)]=(double)gsl_spline_eval(interpObj,(params->outX)[i], params->interpAcc);
  }
  // Tidy up:
  gsl_spline_free(interpObj);
}

// interpolation in second dimension if the input is Float
void interp_second_outFloat(  interp_data_t* params )
{
  int i, j;
  gsl_spline *interpObj;

  interpObj=gsl_spline_alloc(gsl_interp_cspline, (size_t)(params->inN));

  for (j=0; j<(params->M); ++j){
    memcpy(params->y1, &((double*)(params->inData))[j*(params->inN)], (params->inN)*sizeof(double));
    gsl_spline_init(interpObj, params->inX, params->y1,(size_t)(params->inN));
    for (i=0; i<(params->outN); ++i)
      ((float*)(params->outData))[j*(params->sy)+i*(params->sx)]=(float)gsl_spline_eval(interpObj,(params->outX)[i], params->interpAcc);
  }
  gsl_spline_free(interpObj);
}

