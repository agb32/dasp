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
// 
// Urban Bitenc, 2012 July 13
//
// Upgrade 2013 Feb 01:
//        Also works also if input and output matrices are the same, e.g.
//        cmod.interp.gslCubSplineInterp(mxout[0:5, 0:5], yin, xin, yout, xout, mxout, nThreads)
// 
// Functions called from interpmodule.c :: gclCubSplineIntrp
// 
// 
// 
// Some comments:
//
//  - There are four functions used:
//       . interp_first_inFloat, 
//       . interp_first_inDouble, 
//       . interp_second_outFloat,
//       . interp_second_outDouble.
//     It would be easy to write one function with two if clauses instead of 4 functions with no
//     if clauses. I prefer 4 functions with no if clauses because of better readability.
//
//  - The only difference between interp_first_inDouble and interp_first_inFloat is this line:
//       params->y1[i] = (double)(((double*)(par...
//       params->y1[i] = (double)(((float* )(par...
//
//  - The only difference between interp_second_inDouble and interp_second_inFloat is this line:
//    	 ((double*)(params->outData))[j*(params->sy)+i*(params->sx)] = (double)gsl_spline_e...
//       ((float* )(params->outData))[j*(params->sy)+i*(params->sx)] = (float )gsl_spline_e...
//
//  - params->interpAcc can easily be made local (declared and used only in the
//    functions below). However, by adding it to params execution is a bit faster, probably 
//    about 0.8% faster for a 100x100 input and 991x991 output array.
//    Why: in the main program it is initialised only once and then used by both _first_ and
//    _second_, rather than being initialised twice.
//
//  - params->y1: "y1" used to be in "params" for the same reason as interpAcc (combined
//    improvement of 1.5%). In Jan 2013 I realised that _second_ actually does not need it,
//    so I made "y1" local to _first and deleted it from prams.
//
//

#include "interpolate.h"
#include <stdio.h>

//
///////////////// interpolation in first dimension if the input is Double
//
void interp_first_inDouble(  interp_data_t* params )
{
  int i, j;              // for looping: i for x-dim., j for y-dim.
  gsl_spline *interpObj; // for interpolation
  double* y1;            // to get the data from the input into a contiguous array
                         //        that is needed by gsl_spline_init
  double maxx=params->inX[params->inN-1];
  double minx=params->inX[0];
  double val;
  // Malloc y1 (for contigous input data):
  if( (y1 = malloc( params->inN * sizeof(double) )) == NULL )
    {
      printf("%s, line %d: failed to malloc %ld bytes for 'y1'.\n",
	     __FILE__, __LINE__, params->inN * sizeof(double));
      return;
    }

  // Define the interpolation method used:
  interpObj = gsl_spline_alloc(gsl_interp_cspline, (size_t)(params->inN));
  // For each column of the input:
  for (i=0; i<(params->M); ++i)
    {
      // Move the data into the temporary vector y1 to pass it to the interpolation:
      for (j=0; j<(params->inN); ++j)
	y1[j] = (double)(((double*)(params->inData))[ j*(params->sy) + i*(params->sx) ]);
      // Prepare the interpolation function:
      gsl_spline_init(interpObj, params->inX, y1,(size_t)(params->inN));
      // Evaluate the interpolation function and save the result in the intermediate array:
      for (j=0; j<(params->outN); ++j){
	val=params->outX[j];
	if(val>maxx)
	  val=maxx;
	if(val<minx)
	  val=minx;
	((double*)(params->outData))[i+j*(params->sytmp)] = 
	  gsl_spline_eval(interpObj, val, params->interpAcc);
      }
    }
  // Tidy up:
  gsl_spline_free( interpObj );
  free( y1 );
}


//
///////////////// interpolation in first dimension if the input is Float
//
void interp_first_inFloat(  interp_data_t* params )
{
  int i, j;              // for looping: i for x-dim., j for y-dim.
  gsl_spline *interpObj; // for interpolation
  double* y1;            // to get the data from the input into a contiguous array
                         //        that is needed by gsl_spline_init
  double maxx=params->inX[params->inN-1];
  double minx=params->inX[0];
  double val;
  // ? y1 must be double - required by gsl_spline_init ?

  // Malloc y1 (for contigous input data):
  if( (y1 = malloc( params->inN * sizeof(double) )) == NULL )
    {
      printf("%s, line %d: failed to malloc %ld bytes for 'y1'.\n",
	     __FILE__, __LINE__, params->inN * sizeof(double));
      return;
    }

  // Define the interpolation method used:
  interpObj = gsl_spline_alloc(gsl_interp_cspline, (size_t)(params->inN));

  // For each column of the input:
  for (i=0; i<(params->M); ++i)
    {
      // Move the data into the temporary vector y1 to pass it to the interpolation:
      for (j=0; j<(params->inN); ++j)
	y1[j] = (double)(((float*)(params->inData))[ j*(params->sy) + i*(params->sx) ]);
      // Prepare the interpolation function:
      gsl_spline_init(interpObj, params->inX, y1,(size_t)(params->inN));
      // Evaluate the interpolation function and save the result in the intermediate array:
      for (j=0; j<(params->outN); ++j){
	val=params->outX[j];
	if(val>maxx)
	  val=maxx;
	if(val<minx)
	  val=minx;
	((double*)(params->outData))[i+j*(params->sytmp)] = 
	  gsl_spline_eval(interpObj, val, params->interpAcc);
      }
    }
  // Tidy up:
  gsl_spline_free( interpObj );
  free( y1 );
}


//
////////////////// interpolation in second dimension if the input is Double
//
void interp_second_outDouble(  interp_data_t* params )
{
  int i, j;              // for looping
  gsl_spline *interpObj; // for interpolation
  double maxx=params->inX[params->inN-1];
  double minx=params->inX[0];
  double val;

  // Here you define the interpolation method used:
  interpObj=gsl_spline_alloc(gsl_interp_cspline, (size_t)(params->inN));

  // For each line of array:
  for (j=0; j<(params->M); ++j){
    // Prepare the interpolation function:
    // (You can input inData directly, no need to copy to y1, since it is contigous in X.)
    gsl_spline_init(interpObj, params->inX, 
		    &((double*)(params->inData))[j*(params->inN)], (size_t)(params->inN));
    // Evaluate the interpolation function and save the result in the output array:
    if(params->addToOutput){
      for (i=0; i<(params->outN); ++i){
	val=params->outX[i];
	if(val>maxx)
	  val=maxx;
	if(val<minx)
	  val=minx;
	((double*)(params->outData))[j*(params->sy)+i*(params->sx)] += 
	  (double)gsl_spline_eval(interpObj,val, params->interpAcc);
      }
    }else{
      for (i=0; i<(params->outN); ++i){
	val=params->outX[i];
	if(val>maxx)
	  val=maxx;
	if(val<minx)
	  val=minx;
	((double*)(params->outData))[j*(params->sy)+i*(params->sx)] = 
	  (double)gsl_spline_eval(interpObj,val, params->interpAcc);
      }
    }
  }
  // Tidy up:
  gsl_spline_free(interpObj);
}


//
///////////////// interpolation in second dimension if the input is Float
//
void interp_second_outFloat(  interp_data_t* params )
{
  int i, j;              // for looping
  gsl_spline *interpObj; // for interpolation
  double maxx=params->inX[params->inN-1];
  double minx=params->inX[0];
  double val;

  // Here you define the interpolation method used:
  interpObj=gsl_spline_alloc(gsl_interp_cspline, (size_t)(params->inN));

  // For each line of array:
  for (j=0; j<(params->M); ++j){
    // Prepare the interpolation function:
    // (You can input inData directly, no need to copy to y1, since it is contigous in X.)
    gsl_spline_init(interpObj, params->inX, 
		    &((double*)(params->inData))[j*(params->inN)], (size_t)(params->inN));
    // Evaluate the interpolation function and save the result in the output array:
    if(params->addToOutput){
      for (i=0; i<(params->outN); ++i){
	val=params->outX[i];
	if(val>maxx)
	  val=maxx;
	if(val<minx)
	  val=minx;
	((float*)(params->outData))[j*(params->sy)+i*(params->sx)] += 
	  (float)gsl_spline_eval(interpObj,val, params->interpAcc);
      }
    }else{
      for (i=0; i<(params->outN); ++i){
	val=params->outX[i];
	if(val>maxx)
	  val=maxx;
	if(val<minx)
	  val=minx;
	((float*)(params->outData))[j*(params->sy)+i*(params->sx)] = 
	  (float)gsl_spline_eval(interpObj,val, params->interpAcc);
      }
    }
  }
  // Tidy up:
  gsl_spline_free(interpObj);
}




//Used for hexagonal DMs
void interp_first_inFloatStep(  interp_data_t* params )
{
  int i, j;              // for looping: i for x-dim., j for y-dim.
  gsl_spline *interpObj; // for interpolation
  double* y1;            // to get the data from the input into a contiguous array
                         //        that is needed by gsl_spline_init
  double maxx=params->inX[params->inN-1];
  double minx=params->inX[0];
  double val;
  // ? y1 must be double - required by gsl_spline_init ?

  // Malloc y1 (for contigous input data):
  if( (y1 = malloc( params->inN * sizeof(double) )) == NULL )
    {
      printf("%s, line %d: failed to malloc %ld bytes for 'y1'.\n",
	     __FILE__, __LINE__, params->inN * sizeof(double));
      return;
    }

  // Define the interpolation method used:
  interpObj = gsl_spline_alloc(gsl_interp_cspline, (size_t)(params->inN));

  // For each column of the input:
  for (i=0; i<(params->M); ++i)
    {
      //printf("%d %p %p\n",i,params->inData,params->outData);
      // Move the data into the temporary vector y1 to pass it to the interpolation:
      for (j=0; j<(params->inN); ++j){
	//printf("%p\n",&(((float*)(params->inData))[j*(params->sy) + i*(params->sx)]) );
	y1[j] = (double)(((float*)(params->inData))[ j*(params->sy) + i*(params->sx) ]);
      }
      // Prepare the interpolation function:
      gsl_spline_init(interpObj, params->inX, y1,(size_t)(params->inN));
      // Evaluate the interpolation function and save the result in the intermediate array:
      for (j=0; j<(params->outN); ++j){
	val=params->outX[j];
	if(val>maxx)
	  val=maxx;
	if(val<minx)
	  val=minx;
	//printf("%p\n",&(((double*)(params->outData))[i*params->sxtmp+j*(params->sytmp)]));
	((double*)(params->outData))[i*params->sxtmp+j*(params->sytmp)] = 
	  gsl_spline_eval(interpObj, val, params->interpAcc);
      }
    }
  // Tidy up:
  gsl_spline_free( interpObj );
  free( y1 );
}


void interp_second_outFloatStep(  interp_data_t* params )
{
  int i, j;              // for looping
  gsl_spline *interpObj; // for interpolation
  double maxx=params->inX[params->inN-1];
  double minx=params->inX[0];
  double val;

  // Here you define the interpolation method used:
  interpObj=gsl_spline_alloc(gsl_interp_cspline, (size_t)(params->inN));
  printf("%p\n",params->outData);
  // For each line of array:
  for (j=0; j<(params->M); ++j){
    // Prepare the interpolation function:
    // (You can input inData directly, no need to copy to y1, since it is contigous in X.)
    gsl_spline_init(interpObj, params->inX, 
		    &((double*)(params->inData))[j*(params->inN)*2], (size_t)(params->inN));
    // Evaluate the interpolation function and save the result in the output array:
    //printf("%d %d %p\n",j,params->inN,&((double*)(params->inData))[j*(params->inN)*2]);
    for (i=0; i<(params->outN); ++i){
      val=params->outX[i];
      if(val>maxx)
	val=maxx;
      if(val<minx)
	val=minx;
      //printf("%p\n",&(((float*)(params->outData))[j*(params->sy)+i*(params->sx)]));
      ((float*)(params->outData))[j*(params->sy)+i*(params->sx)] = 
	(float)gsl_spline_eval(interpObj,val, params->interpAcc);
    }
  }
  // Tidy up:
  gsl_spline_free(interpObj);
}


