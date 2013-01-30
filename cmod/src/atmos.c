/*C code to do the atmos pupil phase.  Also does the shifting of new columns and rows for the new phase.

NOT YET COMPLETED.

This assumes the infAtmos data follows a circular buffer sort of this - i.e. the new cols aren't added at the end, but somewhere in the middle... to avoid lots of memcopies

*/

typedef struct{
  int nscrn;
  void **data;//the phase screens.
  int isFloat;
  int npup;
  interp_data_t* splineParams; // vector of parameters passed to the funct. executed i
}AtmosStruct;


atmosClose(){
  gsl_interp_accel_free( astr->splineParams);
}

atmosOpen(AtmosStruct *astr){
  astr->splineParams=gsl_interp_accel_alloc();
}

//These next few functions taken from GSL - but I want to make simplified floating point versions.

static inline void coeff_calc (const double c_array[], double dy, double dx, size_t index,  
            double * b, double * c, double * d)
{
  const double c_i = c_array[index];
  const double c_ip1 = c_array[index + 1];
  *b = (dy / dx) - dx * (c_ip1 + 2.0 * c_i) / 3.0;
  *c = c_i;
  *d = (c_ip1 - c_i) / (3.0 * dx);
}

static int cspline_eval (const void * vstate,
              const double x_array[], const double y_array[], size_t size,
              double x,
              gsl_interp_accel * a,
              double *y)
{
  const cspline_state_t *state = (const cspline_state_t *) vstate;

  double x_lo, x_hi;
  double dx;
  size_t index;
  
  if (a != 0)
    {
      index = gsl_interp_accel_find (a, x_array, size, x);
    }
  else
    {
      index = gsl_interp_bsearch (x_array, x, 0, size - 1);
    }
  
  /* evaluate */
  x_hi = x_array[index + 1];
  x_lo = x_array[index];
  dx = x_hi - x_lo;
  if (dx > 0.0)
    {
      const double y_lo = y_array[index];
      const double y_hi = y_array[index + 1];
      const double dy = y_hi - y_lo;
      double delx = x - x_lo;
      double b_i, c_i, d_i; 
      coeff_calc(state->c, dy, dx, index,  &b_i, &c_i, &d_i);
      *y = y_lo + delx * (b_i + delx * (c_i + delx * d_i));
      return GSL_SUCCESS;
    }
  else
    {
      *y = 0.0;
      return GSL_EINVAL;
    }
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
gsl_interp_init_float (gsl_interp * interp, const float x_array[], const float y_array[], size_t size)
{
  size_t i;

  if (size != interp->size)
    {
      GSL_ERROR ("data must match size of interpolation object", GSL_EINVAL);
    }

  for (i = 1; i < size; i++) 
    {
      if (!(x_array[i-1] < x_array[i])) 
        {
          GSL_ERROR ("x values must be monotonically increasing", GSL_EINVAL);
        }
    }

  interp->xmin = x_array[0];
  interp->xmax = x_array[size - 1];

  {
    int status = interp->type->init(interp->state, x_array, y_array, size);
    return status;
  }
}
void
gsl_interp_free (gsl_interp * interp)
{
  RETURN_IF_NULL (interp);

  if (interp->type->free)
    interp->type->free (interp->state);
  free (interp);
}

gsl_interp *gsl_interp_alloc (const gsl_interp_type * T, size_t size)
{
  gsl_interp * interp;

  if (size < T->min_size)
    {
      GSL_ERROR_NULL ("insufficient number of points for interpolation type",
                      GSL_EINVAL);
    }

  interp = (gsl_interp *) malloc (sizeof(gsl_interp));
  
  if (interp == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for interp struct", 
                      GSL_ENOMEM);
    }
  
  interp->type = T;
  interp->size = size;

  if (interp->type->alloc == NULL)
    {
      interp->state = NULL;
      return interp;
    }

  interp->state = interp->type->alloc(size);
  
  if (interp->state == NULL)
    {
      free (interp);          
      GSL_ERROR_NULL ("failed to allocate space for interp state", GSL_ENOMEM);
    };
    
  return interp;
}

gsl_spline *gsl_spline_alloc (const gsl_interp_type * T, size_t size)
{
  gsl_spline * spline = (gsl_spline *) malloc (sizeof(gsl_spline));
  
  if (spline == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for spline struct", 
                      GSL_ENOMEM);
    }
  
  spline->interp = gsl_interp_alloc (T, size);
  
  if (spline->interp == NULL)
    {
      free (spline);          
      GSL_ERROR_NULL ("failed to allocate space for interp", GSL_ENOMEM);
    };
    
  spline->x = (double *) malloc (size * sizeof(double));

  if (spline->x == NULL)
    {
      gsl_interp_free(spline->interp);
      free(spline);
      GSL_ERROR_NULL ("failed to allocate space for x", GSL_ENOMEM);
    }

  spline->y = (double *) malloc (size * sizeof(double));

  if (spline->y == NULL)
    {
      free(spline->x);
      gsl_interp_free(spline->interp);
      free(spline);
      GSL_ERROR_NULL ("failed to allocate space for y", GSL_ENOMEM);
    }
  
  spline->size = size;

  return spline;
}


gsl_spline_init_float (gsl_spline * spline, const float x_array[], const float y_array[], size_t size)
{
  if (size != spline->size)
    {
      GSL_ERROR ("data must match size of spline object", GSL_EINVAL);
    }
  
  memcpy (spline->x, x_array, size * sizeof(float));
  memcpy (spline->y, y_array, size * sizeof(float));

  {
    int status = gsl_interp_init_float (spline->interp, x_array, y_array, size);
    return status;
  }
}


cubicSplineInterp(SplineStruct *sstr,idata,odata){
  //done previously:  params->interpAcc=gsl_interp_accel_alloc();
      params->inN     = n1y;
      params->outN    = n2y;
      params->M       = n1x;
      params->inX     = x1;
      params->outX    = x2;
      params->inData  = pyinData;
      params->outData = (void*)ytmp;
      params->sx      = s1x;
      params->sy      = s1y;

      // (c) Do interpolation pass in first dimension:
      //interp_first_inFloat(params);
      interpObj=gsl_spline_alloc(gsl_interp_cspline, (size_t)n1y);
      
      for (j=0; j<n1x; ++j){
	for (i=0; i<n1y; ++i)
	  params->y1[i] = (double)(((float*)(params->inData))[ i*(params->sy) + j*(params->sx) ]);
	gsl_spline_init(interpObj, params->inX, params->y1,(size_t)(params->inN));
	for (i=0; i<(params->outN); ++i)
	  ((double*)(params->outData))[j+i*(params->sy)] = gsl_spline_eval(interpObj, (params->outX)[i], params->interpAcc);
	
      }
  gsl_spline_free(interpObj);
  
      // (d) Assign the values to the params for interpolation in the second dimension:
      params->inN     = n1x;
      params->outN    = n2x;
      params->M       = n2y;
      params->inX     = x3;
      params->outX    = x4;
      params->inData  = (void*)ytmp;
      params->outData = pyoutData;
      params->sx      = s2x;
      params->sy      = s2y;

      // (e) Do interpolation pass in second dimension and put results in output array:
      interp_second_outFloat(params);



}

int computePupilPhase(AtmosStruct *astr,int zeroData){
  //zeroData, if set means memset the output first.
  int nscrn=astr->nscrn;
  int i;
  if(zeroData)
    memset(outputData,0,sizeof(float)*npup*npup);
  for(iscrn=0;iscrn<nscrn;iscrn++){
    //first select the correct phase, then linear shift it.
    //Note, if an LGS, the output here may be smaller than npup.
    nx=astr->nx[iscrn];
    ny=astr->ny[iscrn];
    scrn=astr->data[iscrn];
    if(sourceAlt>0)
      odata=xxx;//temporary data to be interpolated for cone effect
    else
      odata=xxx;//pupil data.
    if(asFloat){
      for(i=0; i<ny; i++){
	//ii indexes into the infScrn (which is not shifted as new cols are added, so will probably wrap around
	ii=xxx;
	ii1=(ii+1)%something;
	ii*=scrnnx;
	ii1*=scrnny;
	io=i*nx;

	for(j=0; j<nx; j++){
	  jj=xxx;
	  jj1=(jj+1)%something;
	  fans=(float*)scrn[ii+jj]*s1;
	  fans+=(float*)scrn[ii1+jj]*s2;
	  fans+=(float*)scrn[ii+jj+1]*s3;
	  fans+=(float*)scrn[ii1+jj+1]*s4;
	  odata[io+j]=fans;
	}
      }
      if(computeUplinkTT){

      }
      if(sourceAlt>0){//a LGS - interpolate for cone effect.
	cubicSplineInterp(astr->splineStruct,odata,tmp);
	odata=tmp;
      }
      //now add odata to outputData.
      for(i=0;i<npup*npup;i++){
	outputData[i]+=odata[i];
      }
    }
    return 0;

  }
