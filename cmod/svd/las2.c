/*************************************************************************
                           (c) Copyright 2003
                              Douglas Rohde

                     adapted from SVDPACKC, which is

                           (c) Copyright 1993
                        University of Tennessee
                          All Rights Reserved                          
 *************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <fcntl.h>
#include <time.h>

#include "sparsemem.h"
#include "svdlib.h"
#include "svdutil.h"
#include "ritvec.h"
#define MAXLL 2

#define LMTNW   100000000 /* max. size of working area allowed  */

enum storeVals {STORQ = 1, RETRQ, STORP, RETRP};

static char *error_msg[] = {  /* error messages used by function    *
                               * check_parameters                   */
  NULL,
  "",
  "ENDL MUST BE LESS THAN ENDR",
  "REQUESTED DIMENSIONS CANNOT EXCEED NUM ITERATIONS",
  "ONE OF YOUR DIMENSIONS IS LESS THAN OR EQUAL TO ZERO",
  "NUM ITERATIONS (NUMBER OF LANCZOS STEPS) IS INVALID",
  "REQUESTED DIMENSIONS (NUMBER OF EIGENPAIRS DESIRED) IS INVALID",
  "6*N+4*ITERATIONS+1 + ITERATIONS*ITERATIONS CANNOT EXCEED NW",
  "6*N+4*ITERATIONS+1 CANNOT EXCEED NW", NULL};

//double **LanStore;
//float **LanStoreFloat;
double *OPBTemp;
double eps, eps1, reps, eps34;
long ierr;


/*
double rnm, anorm, tol;
FILE *fp_out1, *fp_out2;
*/

void   purge(long n, long ll, double *r, double *q, double *ra,  
             double *qa, double *wrk, double *eta, double *oldeta, long step, 
             double *rnmp, double tol,LanStoreStruct *lss);
void   ortbnd(double *alf, double *eta, double *oldeta, double *bet, long step,
              double rnm);
double startv(SMat A, double *wptr[], long step, long n,LanStoreStruct *lss,int nthreads);
void   store(LanStoreStruct *lss,long, long, long, double *);
void   storeFloat(LanStoreStruct *lss,long, long, long, float *);
//void   imtql2(long, long, double *, double *, double *,dSpMem *sp, float *zf, fSpMem *spf);
void   imtqlb(long n, double d[], double e[], double bnd[]);
void   write_header(long, long, double, double, long, double, long, long, 
                    long);
long   check_parameters(SMat A, long dimensions, long iterations, 
                        double endl, double endr, long vectors);
int    lanso(SMat A, long iterations, long dimensions, double endl,
             double endr, double *ritz, double *bnd, double *wptr[], 
             long *neigp, long n,LanStoreStruct *lss,int nthreads);
//long   ritvec(long n, SMat A, SVDRec R, double kappa, double *ritz, 
//              double *bnd, double *alf, double *bet, STORETYPE *w2, 
//              long steps, long neig,long size,double *genInv,long neigForGenInv,double minEig,double fracEig,int floatS);
long   lanczos_step(SMat A, long first, long last, double *wptr[],
                    double *alf, double *eta, double *oldeta,
                    double *bet, long *ll, long *enough, double *rnmp, 
                    double *tolp, long n,LanStoreStruct *lss,int nthreads);
void   stpone(SMat A, double *wrkptr[], double *rnmp, double *tolp, long n,LanStoreStruct *lss,int nthreads);
long   error_bound(long *, double, double, double *, double *, long step, 
                   double tol);
void   machar(long *ibeta, long *it, long *irnd, long *machep, long *negep);
float macharfloat(long *ibeta, long *it, long *irnd, long *machep, long *negep);

/***********************************************************************
 *                                                                     *
 *                        main()                                       *
 * Sparse SVD(A) via Eigensystem of A'A symmetric Matrix 	       *
 *                  (double precision)                                 *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   This sample program uses landr to compute singular triplets of A via
   the equivalent symmetric eigenvalue problem                         

   B x = lambda x, where x' = (u',v'), lambda = sigma**2,
   where sigma is a singular value of A,
                                                                     
   B = A'A , and A is m (nrow) by n (ncol) (nrow >> ncol),                
                                                                 
   so that {u,sqrt(lambda),v} is a singular triplet of A.        
   (A' = transpose of A)                                      
                                                            
   User supplied routines: svd_opa, opb, store, timer              
                                                        
   svd_opa(     x,y) takes an n-vector x and returns A*x in y.
   svd_opb(ncol,x,y) takes an n-vector x and returns B*x in y.
                                                                  
   Based on operation flag isw, store(n,isw,j,s) stores/retrieves 
   to/from storage a vector of length n in s.                   
                                                               
   User should edit timer() with an appropriate call to an intrinsic
   timing routine that returns elapsed user time.                      


   External parameters 
   -------------------

   Defined and documented in las2.h


   Local parameters 
   ----------------

  (input)
   endl     left end of interval containing unwanted eigenvalues of B
   endr     right end of interval containing unwanted eigenvalues of B
   kappa    relative accuracy of ritz values acceptable as eigenvalues
              of B
	      vectors is not equal to 1
   r        work array
   n	    dimension of the eigenproblem for matrix B (ncol)
   dimensions   upper limit of desired number of singular triplets of A
   iterations   upper limit of desired number of Lanczos steps
   nnzero   number of nonzeros in A
   vectors  1 indicates both singular values and singular vectors are 
	      wanted and they can be found in output file lav2;
	      0 indicates only singular values are wanted 
   		
  (output)
   ritz	    array of ritz values
   bnd      array of error bounds
   d        array of singular values
   memory   total memory allocated in bytes to solve the B-eigenproblem


   Functions used
   --------------

   BLAS		svd_daxpy, svd_dscal, svd_ddot
   USER		svd_opa, svd_opb, timer
   MISC		write_header, check_parameters
   LAS2		landr


   Precision
   ---------

   All floating-point calculations are done in double precision;
   variables are declared as long and double.


   LAS2 development
   ----------------

   LAS2 is a C translation of the Fortran-77 LAS2 from the SVDPACK
   library written by Michael W. Berry, University of Tennessee,
   Dept. of Computer Science, 107 Ayres Hall, Knoxville, TN, 37996-1301

   31 Jan 1992:  Date written 

   Theresa H. Do
   University of Tennessee
   Dept. of Computer Science
   107 Ayres Hall
   Knoxville, TN, 37996-1301
   internet: tdo@cs.utk.edu

 ***********************************************************************/

/***********************************************************************
 *								       *
 *		      check_parameters()			       *
 *								       *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------
   Function validates input parameters and returns error code (long)  

   Parameters 
   ----------
  (input)
   dimensions   upper limit of desired number of eigenpairs of B           
   iterations   upper limit of desired number of lanczos steps             
   n        dimension of the eigenproblem for matrix B               
   endl     left end of interval containing unwanted eigenvalues of B
   endr     right end of interval containing unwanted eigenvalues of B
   vectors  1 indicates both eigenvalues and eigenvectors are wanted 
            and they can be found in lav2; 0 indicates eigenvalues only
   nnzero   number of nonzero elements in input matrix (matrix A)      
                                                                      
 ***********************************************************************/

long check_parameters(SMat A, long dimensions, long iterations, 
		      double endl, double endr, long vectors) {
   long error_index;
   error_index = 0;

   if (endl >/*=*/ endr)  error_index = 2;
   else if (dimensions > iterations) error_index = 3;
   else if (A->cols <= 0 || A->rows <= 0) error_index = 4;
   /*else if (n > A->cols || n > A->rows) error_index = 1;*/
   else if (iterations <= 0 || iterations > A->cols || iterations > A->rows)
     error_index = 5;
   else if (dimensions <= 0 || dimensions > iterations) error_index = 6;
   if (error_index) 
     svd_error("svdLAS2 parameter error: %s\n", error_msg[error_index]);
   return(error_index);
}

/***********************************************************************
 *								       *
 *			  write_header()			       *
 *   Function writes out header of output file containing ritz values  *
 *								       *
 ***********************************************************************/

void write_header(long iterations, long dimensions, double endl, double endr, 
                  long vectors, double kappa, long nrow, long ncol, 
                  long vals) {
  printf("SOLVING THE [A^TA] EIGENPROBLEM\n");
  printf("NO. OF ROWS               = %6ld\n", nrow);
  printf("NO. OF COLUMNS            = %6ld\n", ncol);
  printf("NO. OF NON-ZERO VALUES    = %6ld\n", vals);
  printf("MATRIX DENSITY            = %6.2f%%\n", 
         ((float) vals / nrow) * 100 / ncol);
  /* printf("ORDER OF MATRIX A         = %5ld\n", n); */
  printf("MAX. NO. OF LANCZOS STEPS = %6ld\n", iterations);
  printf("MAX. NO. OF EIGENPAIRS    = %6ld\n", dimensions);
  printf("LEFT  END OF THE INTERVAL = %9.2E\n", endl);
  printf("RIGHT END OF THE INTERVAL = %9.2E\n", endr);
  printf("KAPPA                     = %9.2E\n", kappa);
  /* printf("WANT S-VECTORS?   [T/F]   =     %c\n", (vectors) ? 'T' : 'F'); */
  printf("\n");
  return;
}
void freeLanStore(LanStoreStruct *l){
  long i;
  if(!l)
    return;
  if (l->LanStore) {
    for (i = 0; i < l->iterations + l->maxll; i++)
      SAFE_FREE(l->LanStore[i]);
    SAFE_FREE(l->LanStore);
  }
  if (l->LanStoreFloat) {
    for (i = 0; i < l->iterations + l->maxll; i++)
      SAFE_FREE(l->LanStoreFloat[i]);
    SAFE_FREE(l->LanStoreFloat);
  }
}
/***********************************************************************
 *                                                                     *
 *				landr()				       *
 *        Lanczos algorithm with selective orthogonalization           *
 *                    Using Simon's Recurrence                         *
 *                       (double precision)                            *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   landr() is the LAS2 driver routine that, upon entry,
     (1)  checks for the validity of input parameters of the 
	  B-eigenproblem 
     (2)  determines several machine constants
     (3)  makes a Lanczos run
     (4)  calculates B-eigenvectors (singular vectors of A) if requested 
	  by user


   arguments
   ---------

   (input)
   n        dimension of the eigenproblem for A'A
   iterations   upper limit of desired number of Lanczos steps
   dimensions   upper limit of desired number of eigenpairs
   nnzero   number of nonzeros in matrix A
   endl     left end of interval containing unwanted eigenvalues of B
   endr     right end of interval containing unwanted eigenvalues of B
   vectors  1 indicates both eigenvalues and eigenvectors are wanted
              and they can be found in output file lav2; 
	    0 indicates only eigenvalues are wanted
   kappa    relative accuracy of ritz values acceptable as eigenvalues
	      of B (singular values of A)
   r        work array

   (output)
   j        number of Lanczos steps actually taken                     
   neig     number of ritz values stabilized                           
   ritz     array to hold the ritz values                              
   bnd      array to hold the error bounds


   External parameters
   -------------------

   Defined and documented in las2.h


   local parameters
   -------------------

   ibeta    radix for the floating-point representation
   it       number of base ibeta digits in the floating-point significand
   irnd     floating-point addition rounded or chopped
   machep   machine relative precision or round-off error
   negeps   largest negative integer
   wptr	    array of pointers each pointing to a work space


   Functions used
   --------------

   MISC         svd_dmax, machar, check_parameters
   LAS2         ritvec, lanso

 ***********************************************************************/

SVDRec svdLAS2A(SMat A, long dimensions) {
  double end[2] = {-1.0e-30, 1.0e-30};
  double kappa = 1e-6;
  //R can be NULL, or previously allocated... (useful if you want to use preallocated python arrays).
  //ritvecSize can be 0 or the maximum number of elements to use for this array.  If 0, will attempt to create the full array (which, for a large system, could mean you run out of memory).
  SVDRec R=NULL;
  long ritvecSize=0;
  int userNeig=0;
  int useSFloat=0;
  if (!A) {
    svd_error("svdLAS2A called with NULL array\n");
    return NULL;
  }
  return svdLAS2(A, dimensions, 0, end, kappa,R,ritvecSize,userNeig,NULL,0.,0.,0,useSFloat,0,0,1,0,1);
}


SVDRec svdLAS2(SMat A, long dimensions, long iterations, double end[2], 
               double kappa,SVDRec R,long ritvecSize,int userNeig,ArrUnion *genInv, float minEig, float fracEig,int useStoreFloat, int useSFloat,double minGIVal,int transposeGI,int prepareForGenInv,int considerSwap,int nthreads){
  char transpose = FALSE;
  long ibeta, it, irnd, machep, negep, n, i, steps, nsig, neig, m;
  double *wptr[10], *ritz, *bnd;
  int allocatedR=0;
  LanStoreStruct *lss;
  //SVDRec R = NULL;
  clock_t starttime;
  starttime=clock();
  svdResetCounters();
  printf("rows %u cols %u\n",A->rows,A->cols);
  m = svd_imin((long)A->rows, (long)A->cols);
  if (dimensions <= 0 || dimensions > m)
    dimensions = m;
  if (iterations <= 0 || iterations > m)
    iterations = m;
  if (iterations < dimensions) iterations = dimensions;
  printf("iterations %ld\n",iterations);
  /* Write output header */
  if (SVDVerbosity > 0)
    write_header(iterations, dimensions, end[0], end[1], TRUE, kappa, (long)A->rows, 
                 (long)A->cols, (long)A->vals);

  /* Check parameters */
  if (check_parameters(A, dimensions, iterations, end[0], end[1], TRUE))
    return NULL;

  /* If A is wide, the SVD is computed on its transpose for speed. */
  if (A->cols >= A->rows * 1.2) {
    if (SVDVerbosity > 0) printf("TRANSPOSING THE MATRIX FOR SPEED\n");
    transpose = TRUE;
    A = svdTransposeS(A);
  }
  
  n = (long)A->cols;
  /* Compute machine precision */ 
  if(useStoreFloat)
    eps=(double)macharfloat(&ibeta,&it,&irnd,&machep,&negep);
  else
    machar(&ibeta, &it, &irnd, &machep, &negep);
  printf("Using m=%ld, n=%ld, precision %g, starting lanso\n",m,n,eps);

  eps1 = eps * sqrt((double) n);
  reps = sqrt(eps);
  eps34 = reps * sqrt(reps);

  /* Allocate temporary space. */
  if (!(wptr[0] = svd_doubleArray(n, TRUE, "las2: wptr[0]"))) goto abort;
  if (!(wptr[1] = svd_doubleArray(n, FALSE, "las2: wptr[1]"))) goto abort;
  if (!(wptr[2] = svd_doubleArray(n, FALSE, "las2: wptr[2]"))) goto abort;
  if (!(wptr[3] = svd_doubleArray(n, FALSE, "las2: wptr[3]"))) goto abort;
  if (!(wptr[4] = svd_doubleArray(n, FALSE, "las2: wptr[4]"))) goto abort;
  if (!(wptr[5] = svd_doubleArray(n, FALSE, "las2: wptr[5]"))) goto abort;
  if (!(wptr[6] = svd_doubleArray(iterations, FALSE, "las2: wptr[6]"))) 
    goto abort;
  if (!(wptr[7] = svd_doubleArray(iterations, FALSE, "las2: wptr[7]"))) 
    goto abort;
  if (!(wptr[8] = svd_doubleArray(iterations, FALSE, "las2: wptr[8]"))) 
    goto abort;
  if (!(wptr[9] = svd_doubleArray(iterations + 1, FALSE, "las2: wptr[9]"))) 
    goto abort;
  /* Calloc may be unnecessary: */
  if (!(ritz    = svd_doubleArray(iterations + 1, TRUE, "las2: ritz"))) 
    goto abort;  
  /* Calloc may be unnecessary: */
  if (!(bnd     = svd_doubleArray(iterations + 1, TRUE, "las2: bnd"))) 
    goto abort;
  memset(bnd, 127, (iterations + 1) * sizeof(double));
  if(!(lss=(LanStoreStruct*)malloc(sizeof(LanStoreStruct))))
    goto abort;
  lss->maxll=MAXLL;
  lss->iterations=iterations;
  lss->LanStore=NULL;
  lss->LanStoreFloat=NULL;
  if(useStoreFloat){
    if (!(lss->LanStoreFloat = (float **) calloc(iterations + MAXLL, sizeof(float *))))
      goto abort;
  }else{
    if (!(lss->LanStore = (double **) calloc(iterations + MAXLL, sizeof(double *))))
      goto abort;
  }
  if (!(OPBTemp = svd_doubleArray((long)A->rows, FALSE, "las2: OPBTemp"))) 
    goto abort;

  /* Actually run the lanczos thing: */
  //A is input matrix.  iterations==m (typically), dimensions=m(typically), end[0,1] ar doubles, ritz is double array length iterations, as is bnd. wptr points to 10 arrays of varying lengths.  neig is an output, n is ncents (number of columns, ie matrix has shape m,n).
  steps = lanso(A, iterations, dimensions, end[0], end[1], ritz, bnd, wptr, 
                &neig, n,lss,nthreads);
  printf("Done lanso, time %gs\n",((double)(clock()-starttime)/CLOCKS_PER_SEC));
  /* Print some stuff. */
  if (SVDVerbosity > 0) {
    printf("NUMBER OF LANCZOS STEPS   = %6ld\n"
           "RITZ VALUES STABILIZED    = %6ld\n", steps + 1, neig);
  }
  if (SVDVerbosity > 2) {
    printf("\nCOMPUTED RITZ VALUES  (ERROR BNDS)\n");
    for (i = 0; i <= steps; i++)
      printf("%3ld  %22.14E  (%11.2E)\n", i + 1, ritz[i], bnd[i]);
  }

  SAFE_FREE(wptr[0]);
  SAFE_FREE(wptr[1]);
  SAFE_FREE(wptr[2]);
  SAFE_FREE(wptr[3]);
  SAFE_FREE(wptr[4]);
  SAFE_FREE(wptr[7]);
  SAFE_FREE(wptr[8]);

  /* Compute eigenvectors */
  kappa = svd_dmax(fabs(kappa), eps34);//find out if kappa or eps34 is larger
  //kappa=(fabs(kappa)>eps34?fabs(kappa):eps34);
  
  if(R==NULL){
    allocatedR=1;
    R = svdNewSVDRec();
    if (!R) {
      svd_error("svdLAS2: allocation of R failed");
      goto cleanup;
    }
    //here, allocate some potentially large arrays.  Need a way around this for eagle case...
    R->d  = /*svd_imin(nsig, dimensions)*/dimensions;
    R->Ut = svdNewDMat(R->d, (long)A->rows);
    R->S  = svd_floatArray(R->d, TRUE, "las2: R->s");
    R->Vt = svdNewDMat(R->d, (long)A->cols);
    if (!R->Ut || !R->S || !R->Vt) {
      svd_error("svdLAS2: allocation of R failed");
      goto cleanup;
    }
  }
  if(R->d!=dimensions){
    printf("R->d dimensions wrong, is %d, should be %ld\n",R->d,dimensions);
    goto cleanup;
  }
  if(R->Ut->rows!=R->d || R->Ut->cols!=A->rows){
    printf("R->Ut->rows is %ld, should be %d, cols is %ld should be %u\n",R->Ut->rows,R->d,R->Ut->cols,A->rows);
    goto cleanup;
  }
  if(R->Vt->rows!=R->d || R->Vt->cols!=A->cols){
    printf("R->Vt->rows is %ld, should be %d, cols is %ld should be %u\n",R->Vt->rows,R->d,R->Vt->cols,A->cols);
    goto cleanup;
  }
 
  //agb - here computes the U and V matricees...
  //n is ncents A is sparse matrix, R holds the results.  kappa is a
  //double, ritz is an array length iterations (typically m) as is
  //bnd, wptr are work arrays length m or n, steps is the number of
  //steps, and neig is the number of eigen values.
  //Note, for efficiency, you could set neig to less, eg the number of evals you wish to use...
  if(userNeig<=0 || userNeig>neig)
    userNeig=neig;

  //now decide which version of ritvec we want to use...
  if(lss->LanStore==NULL){//so using LanStoreFloat
    if(R->Vt->value==NULL){//float - ie valuef is assumed to be in use.
      if(useSFloat){//float
	nsig = ritvecfff(n, A, R, kappa, ritz, bnd, wptr[6], wptr[9], (float*)wptr[5],lss, steps, neig,ritvecSize,genInv,userNeig,minEig,fracEig,minGIVal,transposeGI,prepareForGenInv,considerSwap,nthreads,&ierr);
      }else{
	nsig = ritvecdff(n, A, R, kappa, ritz, bnd, wptr[6], wptr[9], (float*)wptr[5],lss, steps, neig,ritvecSize,genInv,userNeig,minEig,fracEig,minGIVal,transposeGI,prepareForGenInv,considerSwap,nthreads,&ierr);
      }
    }else{
      if(useSFloat){//float
	nsig = ritvecffd(n, A, R, kappa, ritz, bnd, wptr[6], wptr[9], (float*)wptr[5],lss, steps, neig,ritvecSize,genInv,userNeig,minEig,fracEig,minGIVal,transposeGI,prepareForGenInv,considerSwap,nthreads,&ierr);
      }else{
	nsig = ritvecdfd(n, A, R, kappa, ritz, bnd, wptr[6], wptr[9], (float*)wptr[5],lss, steps, neig,ritvecSize,genInv,userNeig,minEig,fracEig,minGIVal,transposeGI,prepareForGenInv,considerSwap,nthreads,&ierr);
      }
    }
  }else{//using LanStore (double).
    if(R->Vt->value==NULL){//float - ie valuef is assumed to be in use.
      if(useSFloat){//float
	nsig = ritvecfdf(n, A, R, kappa, ritz, bnd, wptr[6], wptr[9], wptr[5], lss,steps, neig,ritvecSize,genInv,userNeig,minEig,fracEig,minGIVal,transposeGI,prepareForGenInv,considerSwap,nthreads,&ierr);
      }else{
	nsig = ritvecddf(n, A, R, kappa, ritz, bnd, wptr[6], wptr[9], wptr[5], lss,steps, neig,ritvecSize,genInv,userNeig,minEig,fracEig,minGIVal,transposeGI,prepareForGenInv,considerSwap,nthreads,&ierr);
      }
    }else{
      if(useSFloat){//float
	nsig = ritvecfdd(n, A, R, kappa, ritz, bnd, wptr[6], wptr[9], wptr[5], lss,steps, neig,ritvecSize,genInv,userNeig,minEig,fracEig,minGIVal,transposeGI,prepareForGenInv,considerSwap,nthreads,&ierr);
      }else{
	nsig = ritvecddd(n, A, R, kappa, ritz, bnd, wptr[6], wptr[9], wptr[5], lss,steps, neig,ritvecSize,genInv,userNeig,minEig,fracEig,minGIVal,transposeGI,prepareForGenInv,considerSwap,nthreads,&ierr);
      }
    }
  }
  //and that is basically everything done.
  if (SVDVerbosity > 1) {
    printf("\nSINGULAR VALUES: ");
    svdWriteDenseArrayFloat(R->S, R->d, "-", FALSE);

    if (SVDVerbosity > 2) {
      printf("\nLEFT SINGULAR VECTORS (transpose of U): ");
      if(R->Ut->value!=NULL)
	svdWriteDenseMatrix(R->Ut, "-", SVD_F_DT);

      printf("\nRIGHT SINGULAR VECTORS (transpose of V): ");
      svdWriteDenseMatrix(R->Vt, "-", SVD_F_DT);
    }
  } else if (SVDVerbosity > 0)
    printf("SINGULAR VALUES FOUND     = %6d\n", R->d);

 cleanup:    
  for (i = 0; i <= 9; i++)
    SAFE_FREE(wptr[i]);
  SAFE_FREE(ritz);
  SAFE_FREE(bnd);
  freeLanStore(lss);
  SAFE_FREE(lss);
  SAFE_FREE(OPBTemp);

  /* This swaps and transposes the singular matrices if A was transposed. */
  if (R && transpose) {
    DMat T;
    svdFreeSMat(A);
    T = R->Ut;
    R->Ut = R->Vt;
    R->Vt = T;
  }
  R->err=ierr;
  return R;
abort:
  svd_error("svdLAS2: fatal error, aborting");
  return NULL;
}


/***********************************************************************
 *                                                                     *
 *                        ritvec()                                     *
 * 	    Function computes the singular vectors of matrix A	       *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   This function is invoked by landr() only if eigenvectors of the A'A
   eigenproblem are desired.  When called, ritvec() computes the 
   singular vectors of A and writes the result to an unformatted file.


   Parameters
   ----------

   (input)
   nrow       number of rows of A
   steps      number of Lanczos iterations performed
   fp_out2    pointer to unformatted output file
   n	      dimension of matrix A
   kappa      relative accuracy of ritz values acceptable as 
		eigenvalues of A'A
   ritz       array of ritz values
   bnd        array of error bounds
   alf        array of diagonal elements of the tridiagonal matrix T
   bet        array of off-diagonal elements of T
   w1, w2     work space
   shiftVt    int, whether to shift Vt.  Default==1 for normal use.
   (output)
   xv1        array of eigenvectors of A'A (right singular vectors of A)
   ierr	      error code
              0 for normal return from imtql2()
	      k if convergence did not occur for k-th eigenvalue in
	        imtql2()
   nsig       number of accepted ritz values based on kappa

   (local)
   s	      work array which is initialized to the identity matrix
	      of order (j + 1) upon calling imtql2().  After the call,
	      s contains the orthonormal eigenvectors of the symmetric 
	      tridiagonal matrix T

   Functions used
   --------------

   BLAS		svd_dscal, svd_dcopy, svd_daxpy
   USER		store
   		imtql2

 ***********************************************************************/

void rotateArray(double *a, int size, int x) {
  int i, j, n, start;
  double t1, t2;
  printf("Calling rotateArray with size %d, x %d\n",size,x);
  if (x == 0) return;
  j = start = 0;
  t1 = a[0];
  for (i = 0; i < size; i++) {
    n = (j >= x) ? j - x : j + size - x;
    t2 = a[n];
    a[n] = t1;
    t1 = t2;
    j = n;
    if (j == start) {
      start = ++j;
      t1 = a[j];
    }
  }
}
/*
long ritvec(long n, SMat A, SVDRec R, double kappa, double *ritz, double *bnd, 
            double *alf, double *bet, double *w2, long steps, long neig, long size, double *genInv,long neigForGenInv,double minEig,double fracEig,int floatS) {
  //size, added by agb is the maximum space that this array should take.  If the array needs to be larger than this, it will be implemented as a sparse memory matrix (sparsemem.h etc).
  //genInv, the dataspace for the generalised inverse added by agb.  Can be NULL.  If not NULL, V and U may be destroyed.  
  //neigForGenInv is the number of eigenvals to use in generalised inverse.
  //minEig is the value of minimum eigenval that will be accepted, or zero.
  //fracEig is the value of min eigenval as a fraction of max eigenval that will be accepted or zero.
  //If fracEig is set, this will overwrite minEig.
  //If minEig is set, this will overwrite neigForGenInv.
  //floatS, if set, will mean that float arrays are used for s or sp.  Otherwise, they will be double.
  long js, jsq, i, k,  id2, tmp, nsig, x,xoffset;
  double *s=NULL, *xv2, tmp0, tmp1, *w1;// = R->Vt->value[0];
  dSpMem *sp=NULL;
  fSpMem *spf=NULL;
  float *sf=NULL;
  int xx,yy;
  js = steps + 1;
  jsq = js * js;//agb - this could be large for EAGLE - approx nact**2. - ie about 8GB since double.
  printf("ritvec js=%ld, n=%ld\n",js,n);
  xoffset=n-neig;
  if(size==0 || jsq<size){
    printf("Allocating large array s, size %ld elements (%s)\n",jsq,floatS?"float":"double");
    if(floatS)
      sf=svd_floatArray(jsq,TRUE,"ritvec: sf");
    else
      s = svd_doubleArray(jsq, TRUE, "ritvec: s");//agb - large memory.
  }else{
    s=NULL;
    printf("Allocating sparse memory structure with %ld elements for a %ld x %ld array\n",size,js,js);
    //sp=smNewSparseMem(size,jsq,1);
    if(floatS)
      spf=smNewSparseMemFloat(size,js,js);
    else
      sp=smNewSparseMem(size,js,js);
  }
  xv2 = svd_doubleArray(n, FALSE, "ritvec: xv2");
  
  // initialize s to an identity matrix 
  if(s!=NULL){
    for (i = 0; i < jsq; i+= (js+1))
      s[i] = 1.0;
  }else if(sf!=NULL){
    for(i=0; i<jsq; i+=(js+1))
      sf[i]=1.;
  }else if(sp!=NULL){
    for(i=0; i<js; i++)
      smInsert(sp,i,i,1.0);
  }else{
    for(i=0; i<js; i++)
      smInsertFloat(spf,i,i,1.);
  }
  //svd_dcopy(js, alf, 1, w1, -1);//agb - w1=Vt->value[0].  copy from alf to w1 in reverse order (I think). alf is array of diagonal elements of the tridiagonal matrix T, w1 is the first column of Vt.  Agb replaced this with vecReverse.
  vecReverse(js,alf);//agb added this to use alf instead of w1 in imtql2.
  svd_dcopy(steps, &bet[1], 1, &w2[1], -1);
  
  // on return from imtql2(), w1 contains eigenvalues in ascending 
  // order and s contains the corresponding eigenvectors 
  printf("Starting imtql2\n");
  imtql2(js, js, alf, w2, s,sp,sf,spf);//w1=Vt->value[0] - agb - alf used to be w1.
  printf("Done imtql2 (ierr=%ld)\n",ierr);
  if (ierr) return 0;//agb - ierr is a global variable- so not thread safe.
  
  nsig = 0;
  x = 0;
  id2 = jsq - js;
  xx=0;//added by agb for 2D s array...
  for (k = 0; k < js; k++) {
    tmp = id2;
    yy=js-1;//added by agb for 2D s array...
    if (bnd[k] <= kappa * fabs(ritz[k]) && k > js-neig-1) {
      if (--x < 0)
	x = R->d - 1;
      //Note, there is no reason to have w1 as Vt.  It could be an allocated array, and the result at the end of each iteration saved to disk, in eg Vt_x.fits where x is the row number. (note in python Vt is Ut since transposed.
      if(x>=xoffset)
	w1 = R->Vt->value[x-xoffset];
      else{
	printf("Error: las2.c - x<xoffset (%ld<%ld)\n",x,xoffset);
	w1=R->Vt->value[0];
      }
      memset(w1,0,sizeof(double)*n);//w1=Vt->value[x]
      //for (i = 0; i < n; i++) 
      //w1[i] = 0.0;//w1=Vt->value[x]
      if(s!=NULL){
	for (i = 0; i < js; i++) {
	  store(n, RETRQ, i, w2);//retrieve w2 from the store.  This could be slow for a disk store.
	  svd_daxpy(n, s[tmp], w2, 1, w1, 1);//w1=s[tmp]*w2+w1 (w1=Vt->value[x],w2 are vectors) All of s gets accessed, but some values are much smaller than others...
	  tmp -= js;
	}
      }else if(sf!=NULL){
	for (i = 0; i < js; i++) {
	  store(n, RETRQ, i, w2);//retrieve w2 from the store.  This could be slow for a disk store.
	  svd_daxpy(n, (double)sf[tmp], w2, 1, w1, 1);//w1=s[tmp]*w2+w1 (w1=Vt->value[x],w2 are vectors) All of s gets accessed, but some values are much smaller than others...
	  tmp -= js;
	}
      }else if(sp!=NULL){//using sparse matrix for s...
	long e=sp->indptr[xx+1];
	for(i=sp->indptr[xx]; i<e; i++){
	  //row=sp->rowind[i];
	  //val=sp->data[i];
	  store(n,RETRQ,js-1-sp->rowind[i],w2);
	  svd_daxpy(n,sp->data[i],w2,1,w1,1);//W1=Vt->value[x]
	}
      }else{// if(spf!=NULL){
	long e=spf->indptr[xx+1];
	for(i=spf->indptr[xx]; i<e; i++){
	  //row=sp->rowind[i];
	  //val=sp->data[i];
	  store(n,RETRQ,js-1-spf->rowind[i],w2);
	  svd_daxpy(n,(double)spf->data[i],w2,1,w1,1);//W1=Vt->value[x]
	}
      }
      
      //agb - could we compute R->S (evals) here?
      svd_opb(A, w1, xv2, OPBTemp);//compute Vt->value[x].AT.A, output in xv2.
      tmp0 = svd_ddot(n, w1, 1, xv2, 1);//dot Vt with xv2 (result is a scalar).
      tmp0=sqrt(tmp0);
      if(x>=xoffset)
	R->S[x-xoffset] = tmp0;
      else
	printf("Error: las2.c - x<xoffset (%ld<%ld)\n",x,xoffset);
      //end R->S computation.
      nsig++;
    }
    id2++;
    xx++;//added by agb
  }
  SAFE_FREE(s);
  smFreeSparseMem(sp);
  SAFE_FREE(sf);
  smFreeSparseMemFloat(spf);
  //Note, from this stage onwards, Vt is not changed - it could simply by loaded from disk as needed.

  // agb - removed Rotate the singular vectors and values. 
  //rotateArray(R->Vt->value[0], R->Vt->rows * R->Vt->cols, x * R->Vt->cols);
  // x is now the location of the highest singular value. 
  xoffset=0;
  R->d = svd_imin(R->d, nsig);
  //printf("rotating at x=%ld (this point moved to zero.  R->d now %ld\n",x,R->d);//577 downto 27.
  for (x = 0; x < R->d; x++) {
    // multiply by matrix B first 
//    svd_opb(A, R->Vt->value[x+xoffset], xv2, OPBTemp);//compute Vt->value[x].AT.A, output in xv2.
//    tmp0 = svd_ddot(n, R->Vt->value[x+xoffset], 1, xv2, 1);//dot Vt with xv2 (result is a scalar).
//    //svd_daxpy(n, -tmp0, R->Vt->value[x+xoffset], 1, xv2, 1);//agb removed this.xv2=xv2+-tmp0*Vt->value[x]
//    tmp0 = sqrt(tmp0);
//    R->S[x] = tmp0;
    //xnorm = sqrt(svd_ddot(n, xv2, 1, xv2, 1));//agb removed this.
      
    // multiply by matrix A to get (scaled) left s-vector 
    //This is the first place that U is used.
    w1=R->Vt->value[x+xoffset];
    svd_opa(A, w1, R->Ut->value[x]);//A.Vt->value[x] result stored in Ut->value[x].
    tmp1=(R->S[x]==0.?0.:1./R->S[x]);//R->S used to be tmp0.
    svd_dscal(A->rows, tmp1, R->Ut->value[x], 1);
    if(genInv!=NULL){
      //perform the first part of the generalised inverse.  This is multiplying V by 1/S.
      //Note, this destroys V.
      svd_dscal(R->Vt->cols,tmp1,w1,1);
    }
    //xnorm *= tmp1;//agb removed this.
    //bnd[i] = xnorm;//agb removed this.
    //If wanted to compute the reconstruction matrix in place, here would be the place to do it.
    //In python, the reconstructor is V.1/w.UT.
    //  However, here V and U are swapped around (since A is transposed).  So, I would want to do:
    //  U.1/w.VT
    //  The 1/w.VT simply scales each row [i] of VT by 1/w[i].  And R->S[x] is the w.
    //  This will have been done in the previous stage.
    //  After this, can then dot the cols of VT with the rows of U, ie the cols of UT.
    //  So, I could do something like (assuming Vt can be destroyed):
    //  svd_dscal(neig,1./tmp0,R->Vt->value[x+xoffset],1);
    //  save(Vt->value[x+xoffset]);//store it back in memory.
    //  //We need to dot every col of Ut by every col of Vt to get the reconmx (or rather, 
    //  //the neig first values of each col).

  }
  if(fracEig>0.)
    minEig=R->S[0]*fracEig;
  if(minEig>0.){
    i=0;
    while(i<neig && R->S[i]>=minEig)
      i++;
    if(i>1)
      neigForGenInv=i-1;
  }
  //Initially, assuming all U and V stored in memory, for testing.  ie dot U with Vt.
  if(genInv!=NULL){
    printf("Computing generalised inverse: ncent (ncols)=%ld, nact (nrows)=%ld, neig=%ld, neigForGenInv=%ld\n",R->Ut->cols,R->Vt->rows,neig,neigForGenInv);
    printf("Ut->cols %ld Ut->rows %ld Vt->cols %ld Vt->rows %ld\n",R->Ut->cols,R->Ut->rows,R->Vt->cols,R->Vt->rows);
    for(i=0; i<R->Vt->cols; i++){//nact
      for(k=0; k<R->Ut->cols; k++){//ncent
	genInv[i*R->Ut->cols+k]=svd_ddot(neigForGenInv,&(R->Ut->value[0][k]),R->Ut->cols,&(R->Vt->value[0][i]),R->Vt->cols);
      }
    }
  }
  //for(i=0; i<neig; i++){//display Vt[:neig,9],Ut[:neig,9] - agrees with expected.
  // printf("%ld    %g    %g\n",i,R->Ut->value[0][9+i*R->Ut->cols],R->Vt->value[0][9+i*R->Vt->rows]);
  //}
  SAFE_FREE(xv2);
  return nsig;
}
*/
/***********************************************************************
 *                                                                     *
 *                          lanso()                                    *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   Function determines when the restart of the Lanczos algorithm should 
   occur and when it should terminate.

   Arguments 
   ---------

   (input)
   n         dimension of the eigenproblem for matrix B
   iterations    upper limit of desired number of lanczos steps           
   dimensions    upper limit of desired number of eigenpairs             
   endl      left end of interval containing unwanted eigenvalues
   endr      right end of interval containing unwanted eigenvalues
   ritz      array to hold the ritz values                       
   bnd       array to hold the error bounds                          
   wptr      array of pointers that point to work space:            
  	       wptr[0]-wptr[5]  six vectors of length n		
  	       wptr[6] array to hold diagonal of the tridiagonal matrix T
  	       wptr[9] array to hold off-diagonal of T	
  	       wptr[7] orthogonality estimate of Lanczos vectors at 
		 step j
 	       wptr[8] orthogonality estimate of Lanczos vectors at 
		 step j-1

   (output)
   j         number of Lanczos steps actually taken
   neig      number of ritz values stabilized
   ritz      array to hold the ritz values
   bnd       array to hold the error bounds
   ierr      (globally declared) error flag
	     ierr = 8192 if stpone() fails to find a starting vector
	     ierr = k if convergence did not occur for k-th eigenvalue
		    in imtqlb()
	     ierr = 0 otherwise


   Functions used
   --------------

   LAS		stpone, error_bound, lanczos_step
   MISC		svd_dsort2
   UTILITY	svd_imin, svd_imax

 ***********************************************************************/

//Typical inputs for this for reconstruction are: A is input matrix.
//iterations==m (typically), dimensions=m(typically), end[0,1] are
//doubles, ritz is double array length iterations, as is bnd. wptr
//points to 10 arrays of varying lengths.  neig is an output, n is
//ncents (number of columns, ie matrix has shape m,n).

int lanso(SMat A, long iterations, long dimensions, double endl,
          double endr, double *ritz, double *bnd, double *wptr[], 
          long *neigp, long n,LanStoreStruct *lss,int nthreads) {
  double *alf, *eta, *oldeta, *bet, *wrk, rnm, tol;
  long ll, first, last, ENOUGH, id2, id3, i, l, neig, j = 0, intro = 0;
  
  alf = wptr[6];
  eta = wptr[7];
  oldeta = wptr[8];
  bet = wptr[9];
  wrk = wptr[5];
  
  /* take the first step */
  stpone(A, wptr, &rnm, &tol, n,lss,nthreads);
  if (!rnm || ierr) return 0;
  eta[0] = eps1;
  oldeta[0] = eps1;
  ll = 0;
  first = 1;
  last = svd_imin(dimensions + svd_imax(8, dimensions), iterations);
  ENOUGH = FALSE;
  /*id1 = 0;*/
  while (/*id1 < dimensions && */!ENOUGH) {
    if (rnm <= tol) rnm = 0.0;
    
    /* the actual lanczos loop */
    j = lanczos_step(A, first, last, wptr, alf, eta, oldeta, bet, &ll,
                     &ENOUGH, &rnm, &tol, n,lss,nthreads);
    if (ENOUGH) j = j - 1;
    else j = last - 1;
    first = j + 1;
    bet[j+1] = rnm;
    
    /* analyze T */
    l = 0;
    for (id2 = 0; id2 < j; id2++) {
      if (l > j) break;
      for (i = l; i <= j; i++) if (!bet[i+1]) break;
      if (i > j) i = j;
      
      /* now i is at the end of an unreduced submatrix */
      svd_dcopy(i-l+1, &alf[l],   1, &ritz[l],  -1);
      svd_dcopy(i-l,   &bet[l+1], 1, &wrk[l+1], -1);
      
      imtqlb(i-l+1, &ritz[l], &wrk[l], &bnd[l]);
      
      if (ierr) {
        svd_error("svdLAS2: imtqlb failed to converge (ierr = %ld)\n", ierr);
        svd_error("  l = %ld  i = %ld\n", l, i);
        for (id3 = l; id3 <= i; id3++) 
          svd_error("  %ld  %lg  %lg  %lg\n", 
                    id3, ritz[id3], wrk[id3], bnd[id3]);
      }
      for (id3 = l; id3 <= i; id3++) 
        bnd[id3] = rnm * fabs(bnd[id3]);
      l = i + 1;
    }
    
    /* sort eigenvalues into increasing order */
    svd_dsort2((j+1) / 2, j + 1, ritz, bnd);

    /*    for (i = 0; i < iterations; i++)
      printf("%f ", ritz[i]);
      printf("\n"); */
    
    /* massage error bounds for very close ritz values */
    neig = error_bound(&ENOUGH, endl, endr, ritz, bnd, j, tol);
    *neigp = neig;
    
    /* should we stop? */
    if (neig < dimensions) {
      if (!neig) {
        last = first + 9;
        intro = first;
      } else last = first + svd_imax(3, 1 + ((j - intro) * (dimensions-neig)) /
                                     neig);
      last = svd_imin(last, iterations);
    } else ENOUGH = TRUE;
    ENOUGH = ENOUGH || first >= iterations;
    /* id1++; */
    /* printf("id1=%d dimen=%d first=%d\n", id1, dimensions, first); */
  }
  store(lss,n,STORQ,j,wptr[1]);
  //if(lss->LanStore==NULL){
  //  storeFloatUsingDouble(n,STORQ,j,wptr[1]);
  //}else{
  //  store(n, STORQ, j, wptr[1]);
  // }
  return j;
}


/***********************************************************************
 *                                                                     *
 *			lanczos_step()                                 *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   Function embodies a single Lanczos step

   Arguments 
   ---------

   (input)
   n        dimension of the eigenproblem for matrix B
   first    start of index through loop				      
   last     end of index through loop				     
   wptr	    array of pointers pointing to work space		    
   alf	    array to hold diagonal of the tridiagonal matrix T
   eta      orthogonality estimate of Lanczos vectors at step j   
   oldeta   orthogonality estimate of Lanczos vectors at step j-1
   bet      array to hold off-diagonal of T                     
   ll       number of intitial Lanczos vectors in local orthog. 
              (has value of 0, 1 or 2)			
   enough   stop flag			

   Functions used
   --------------

   BLAS		svd_ddot, svd_dscal, svd_daxpy, svd_datx, svd_dcopy
   USER		store
   LAS		purge, ortbnd, startv
   UTILITY	svd_imin, svd_imax

 ***********************************************************************/

long lanczos_step(SMat A, long first, long last, double *wptr[],
		  double *alf, double *eta, double *oldeta,
		  double *bet, long *ll, long *enough, double *rnmp, 
                  double *tolp, long n,LanStoreStruct *lss,int nthreads) {
   double t, *mid, rnm = *rnmp, tol = *tolp, anorm;
   long i, j;

   for (j=first; j<last; j++) {
      mid     = wptr[2];
      wptr[2] = wptr[1];
      wptr[1] = mid;
      mid     = wptr[3];
      wptr[3] = wptr[4];
      wptr[4] = mid;
      store(lss,n,STORQ,j-1,wptr[2]);
      if(j-1<MAXLL)
	store(lss,n,STORP,j-1,wptr[2]);
      //if(LanStore==NULL){
      //storeFloatUsingDouble(n, STORQ, j-1, wptr[2]);
      //if (j-1 < MAXLL)
      //  storeFloatUsingDouble(n, STORP, j-1, wptr[4]);
      //}else{
      //store(n, STORQ, j-1, wptr[2]);
      //if (j-1 < MAXLL)
      //  store(n, STORP, j-1, wptr[4]);
      //}
      bet[j] = rnm;

      /* restart if invariant subspace is found */
      if (!bet[j]) {
	rnm = startv(A, wptr, j, n,lss,nthreads);
	 if (ierr) return j;
	 if (!rnm) *enough = TRUE;
      }
      if (*enough) {
        /* added by Doug... */
        /* These lines fix a bug that occurs with low-rank matrices */
        mid     = wptr[2];
        wptr[2] = wptr[1];
        wptr[1] = mid;
        /* ...added by Doug */
        break;
      }

      /* take a lanczos step */
      t = 1.0 / rnm;
      svd_datx(n, t, wptr[0], 1, wptr[1], 1);
      svd_dscal(n, t, wptr[3], 1);
      svd_opb_t(A, wptr[3], wptr[0], OPBTemp,nthreads);
      svd_daxpy(n, -rnm, wptr[2], 1, wptr[0], 1);
      alf[j] = svd_ddot(n, wptr[0], 1, wptr[3], 1);
      svd_daxpy(n, -alf[j], wptr[1], 1, wptr[0], 1);

      /* orthogonalize against initial lanczos vectors */
      if (j <= MAXLL && (fabs(alf[j-1]) > 4.0 * fabs(alf[j])))
	 *ll = j;  
      for (i=0; i < svd_imin(*ll, j-1); i++) {
	store(lss,n, RETRP, i, wptr[5]);
	t = svd_ddot(n, wptr[5], 1, wptr[0], 1);
	store(lss,n, RETRQ, i, wptr[5]);
	svd_daxpy(n, -t, wptr[5], 1, wptr[0], 1);
	eta[i] = eps1;
	oldeta[i] = eps1;
      }

      /* extended local reorthogonalization */
      t = svd_ddot(n, wptr[0], 1, wptr[4], 1);
      svd_daxpy(n, -t, wptr[2], 1, wptr[0], 1);
      if (bet[j] > 0.0) bet[j] = bet[j] + t;
      t = svd_ddot(n, wptr[0], 1, wptr[3], 1);
      svd_daxpy(n, -t, wptr[1], 1, wptr[0], 1);
      alf[j] = alf[j] + t;
      svd_dcopy(n, wptr[0], 1, wptr[4], 1);
      rnm = sqrt(svd_ddot(n, wptr[0], 1, wptr[4], 1));
      anorm = bet[j] + fabs(alf[j]) + rnm;
      tol = reps * anorm;

      /* update the orthogonality bounds */
      ortbnd(alf, eta, oldeta, bet, j, rnm);

      /* restore the orthogonality state when needed */
      purge(n, *ll, wptr[0], wptr[1], wptr[4], wptr[3], wptr[5], eta, oldeta,
            j, &rnm, tol,lss);
      if (rnm <= tol) rnm = 0.0;
   }
   *rnmp = rnm;
   *tolp = tol;
   return j;
}

/***********************************************************************
 *                                                                     *
 *                          ortbnd()                                   *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   Funtion updates the eta recurrence

   Arguments 
   ---------

   (input)
   alf      array to hold diagonal of the tridiagonal matrix T         
   eta      orthogonality estimate of Lanczos vectors at step j        
   oldeta   orthogonality estimate of Lanczos vectors at step j-1     
   bet      array to hold off-diagonal of T                          
   n        dimension of the eigenproblem for matrix B		    
   j        dimension of T					  
   rnm	    norm of the next residual vector			 
   eps1	    roundoff estimate for dot product of two unit vectors

   (output)
   eta      orthogonality estimate of Lanczos vectors at step j+1     
   oldeta   orthogonality estimate of Lanczos vectors at step j        


   Functions used
   --------------

   BLAS		svd_dswap

 ***********************************************************************/

void ortbnd(double *alf, double *eta, double *oldeta, double *bet, long step,
            double rnm) {
   long i;
   if (step < 1) return;
   if (rnm) {
      if (step > 1) {
	 oldeta[0] = (bet[1] * eta[1] + (alf[0]-alf[step]) * eta[0] -
		      bet[step] * oldeta[0]) / rnm + eps1;
      }
      for (i=1; i<=step-2; i++) 
	 oldeta[i] = (bet[i+1] * eta[i+1] + (alf[i]-alf[step]) * eta[i] +
		      bet[i] * eta[i-1] - bet[step] * oldeta[i])/rnm + eps1;
   }
   oldeta[step-1] = eps1;
   svd_dswap(step, oldeta, 1, eta, 1);  
   eta[step] = eps1;
   return;
}

/***********************************************************************
 *                                                                     *
 *				purge()                                *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   Function examines the state of orthogonality between the new Lanczos
   vector and the previous ones to decide whether re-orthogonalization 
   should be performed


   Arguments 
   ---------

   (input)
   n        dimension of the eigenproblem for matrix B		       
   ll       number of intitial Lanczos vectors in local orthog.       
   r        residual vector to become next Lanczos vector            
   q        current Lanczos vector			           
   ra       previous Lanczos vector
   qa       previous Lanczos vector
   wrk      temporary vector to hold the previous Lanczos vector
   eta      state of orthogonality between r and prev. Lanczos vectors 
   oldeta   state of orthogonality between q and prev. Lanczos vectors
   j        current Lanczos step				     

   (output)
   r	    residual vector orthogonalized against previous Lanczos 
	      vectors
   q        current Lanczos vector orthogonalized against previous ones


   Functions used
   --------------

   BLAS		svd_daxpy,  svd_dcopy,  svd_idamax,  svd_ddot
   USER		store

 ***********************************************************************/

void purge(long n, long ll, double *r, double *q, double *ra,  
	   double *qa, double *wrk, double *eta, double *oldeta, long step, 
           double *rnmp, double tol,LanStoreStruct *lss) {
  double t, tq, tr, reps1, rnm = *rnmp;
  long k, iteration, flag, i;
  
  if (step < ll+2) return; 
  
  k = svd_idamax(step - (ll+1), &eta[ll], 1) + ll;
  if (fabs(eta[k]) > reps) {
    reps1 = eps1 / reps;
    iteration = 0;
    flag = TRUE;
    while (iteration < 2 && flag) {
      if (rnm > tol) {
        
        /* bring in a lanczos vector t and orthogonalize both 
         * r and q against it */
        tq = 0.0;
        tr = 0.0;
	for (i = ll; i < step; i++) {
	  store(lss,n,  RETRQ,  i,  wrk);
	  t   = -svd_ddot(n, qa, 1, wrk, 1);
	  tq += fabs(t);
	  svd_daxpy(n,  t,  wrk,  1,  q,  1);
	  t   = -svd_ddot(n, ra, 1, wrk, 1);
	  tr += fabs(t);
	  svd_daxpy(n, t, wrk, 1, r, 1);
	}
        svd_dcopy(n, q, 1, qa, 1);
        t   = -svd_ddot(n, r, 1, qa, 1);
        tr += fabs(t);
        svd_daxpy(n, t, q, 1, r, 1);
        svd_dcopy(n, r, 1, ra, 1);
        rnm = sqrt(svd_ddot(n, ra, 1, r, 1));
        if (tq <= reps1 && tr <= reps1 * rnm) flag = FALSE;
      }
      iteration++;
    }
    for (i = ll; i <= step; i++) { 
      eta[i] = eps1;
      oldeta[i] = eps1;
    }
  }
  *rnmp = rnm;
  return;
}


/***********************************************************************
 *                                                                     *
 *                         stpone()                                    *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   Function performs the first step of the Lanczos algorithm.  It also
   does a step of extended local re-orthogonalization.

   Arguments 
   ---------

   (input)
   n      dimension of the eigenproblem for matrix B

   (output)
   ierr   error flag
   wptr   array of pointers that point to work space that contains
	    wptr[0]             r[j]
	    wptr[1]             q[j]
	    wptr[2]             q[j-1]
	    wptr[3]             p
	    wptr[4]             p[j-1]
	    wptr[6]             diagonal elements of matrix T 


   Functions used
   --------------

   BLAS		svd_daxpy, svd_datx, svd_dcopy, svd_ddot, svd_dscal
   USER		store, opb
   LAS		startv

 ***********************************************************************/

void stpone(SMat A, double *wrkptr[], double *rnmp, double *tolp, long n,LanStoreStruct *lss,int nthreads) {
   double t, *alf, rnm, anorm;
   alf = wrkptr[6];

   /* get initial vector; default is random */
   rnm = startv(A, wrkptr, 0, n,lss,nthreads);
   if (rnm == 0.0 || ierr != 0) return;

   /* normalize starting vector */
   t = 1.0 / rnm;
   svd_datx(n, t, wrkptr[0], 1, wrkptr[1], 1);
   svd_dscal(n, t, wrkptr[3], 1);

   /* take the first step */
   svd_opb_t(A, wrkptr[3], wrkptr[0], OPBTemp,nthreads);
   alf[0] = svd_ddot(n, wrkptr[0], 1, wrkptr[3], 1);
   svd_daxpy(n, -alf[0], wrkptr[1], 1, wrkptr[0], 1);
   t = svd_ddot(n, wrkptr[0], 1, wrkptr[3], 1);
   svd_daxpy(n, -t, wrkptr[1], 1, wrkptr[0], 1);
   alf[0] += t;
   svd_dcopy(n, wrkptr[0], 1, wrkptr[4], 1);
   rnm = sqrt(svd_ddot(n, wrkptr[0], 1, wrkptr[4], 1));
   anorm = rnm + fabs(alf[0]);
   *rnmp = rnm;
   *tolp = reps * anorm;

   return;
}

/***********************************************************************
 *                                                                     *
 *                         startv()                                    *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   Function delivers a starting vector in r and returns |r|; it returns 
   zero if the range is spanned, and ierr is non-zero if no starting 
   vector within range of operator can be found.

   Parameters 
   ---------

   (input)
   n      dimension of the eigenproblem matrix B
   wptr   array of pointers that point to work space
   j      starting index for a Lanczos run
   eps    machine epsilon (relative precision)

   (output)
   wptr   array of pointers that point to work space that contains
	  r[j], q[j], q[j-1], p[j], p[j-1]
   ierr   error flag (nonzero if no starting vector can be found)

   Functions used
   --------------

   BLAS		svd_ddot, svd_dcopy, svd_daxpy
   USER		svd_opb, store
   MISC		random

 ***********************************************************************/

double startv(SMat A, double *wptr[], long step, long n,LanStoreStruct *lss,int nthreads) {
   double rnm2, *r, t;
   long irand;
   long id, i;

   /* get initial vector; default is random */
   rnm2 = svd_ddot(n, wptr[0], 1, wptr[0], 1);
   irand = 918273 + step;
   r = wptr[0];
   for (id = 0; id < 3; id++) {
      if (id > 0 || step > 0 || rnm2 == 0) 
	 for (i = 0; i < n; i++) r[i] = svd_random2(&irand);
      svd_dcopy(n, wptr[0], 1, wptr[3], 1);

      /* apply operator to put r in range (essential if m singular) */
      svd_opb_t(A, wptr[3], wptr[0], OPBTemp,nthreads);
      svd_dcopy(n, wptr[0], 1, wptr[3], 1);
      rnm2 = svd_ddot(n, wptr[0], 1, wptr[3], 1);
      if (rnm2 > 0.0) break;
   }

   /* fatal error */
   if (rnm2 <= 0.0) {
      ierr = 8192;
      return(-1);
   }
   if (step > 0) {
     for (i = 0; i < step; i++) {
       store(lss,n, RETRQ, i, wptr[5]);
       t = -svd_ddot(n, wptr[3], 1, wptr[5], 1);
       svd_daxpy(n, t, wptr[5], 1, wptr[0], 1);
     }
      /* make sure q[step] is orthogonal to q[step-1] */
      t = svd_ddot(n, wptr[4], 1, wptr[0], 1);
      svd_daxpy(n, -t, wptr[2], 1, wptr[0], 1);
      svd_dcopy(n, wptr[0], 1, wptr[3], 1);
      t = svd_ddot(n, wptr[3], 1, wptr[0], 1);
      if (t <= eps * rnm2) t = 0.0;
      rnm2 = t;
   }
   return(sqrt(rnm2));
}

/***********************************************************************
 *                                                                     *
 *			error_bound()                                  *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   Function massages error bounds for very close ritz values by placing 
   a gap between them.  The error bounds are then refined to reflect 
   this.


   Arguments 
   ---------

   (input)
   endl     left end of interval containing unwanted eigenvalues
   endr     right end of interval containing unwanted eigenvalues
   ritz     array to store the ritz values
   bnd      array to store the error bounds
   enough   stop flag


   Functions used
   --------------

   BLAS		svd_idamax
   UTILITY	svd_dmin

 ***********************************************************************/

long error_bound(long *enough, double endl, double endr, 
                 double *ritz, double *bnd, long step, double tol) {
  long mid, i, neig;
  double gapl, gap;
  
  /* massage error bounds for very close ritz values */
  mid = svd_idamax(step + 1, bnd, 1);

  for (i=((step+1) + (step-1)) / 2; i >= mid + 1; i -= 1)
    if (fabs(ritz[i-1] - ritz[i]) < eps34 * fabs(ritz[i])) 
      if (bnd[i] > tol && bnd[i-1] > tol) {
        bnd[i-1] = sqrt(bnd[i] * bnd[i] + bnd[i-1] * bnd[i-1]);
        bnd[i] = 0.0;
      }
  
  
  for (i=((step+1) - (step-1)) / 2; i <= mid - 1; i +=1 ) 
    if (fabs(ritz[i+1] - ritz[i]) < eps34 * fabs(ritz[i])) 
      if (bnd[i] > tol && bnd[i+1] > tol) {
        bnd[i+1] = sqrt(bnd[i] * bnd[i] + bnd[i+1] * bnd[i+1]);
        bnd[i] = 0.0;
      }
  
  /* refine the error bounds */
  neig = 0;
  gapl = ritz[step] - ritz[0];
  for (i = 0; i <= step; i++) {
    gap = gapl;
    if (i < step) gapl = ritz[i+1] - ritz[i];
    gap = svd_dmin(gap, gapl);
    if (gap > bnd[i]) bnd[i] = bnd[i] * (bnd[i] / gap);
    if (bnd[i] <= 16.0 * eps * fabs(ritz[i])) {
      neig++;
      if (!*enough) *enough = endl < ritz[i] && ritz[i] < endr;
    }
  }   
  return neig;
}

/***********************************************************************
 *                                                                     *
 *				imtqlb()			       *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   imtqlb() is a translation of a Fortran version of the Algol
   procedure IMTQL1, Num. Math. 12, 377-383(1968) by Martin and 
   Wilkinson, as modified in Num. Math. 15, 450(1970) by Dubrulle.  
   Handbook for Auto. Comp., vol.II-Linear Algebra, 241-248(1971).  
   See also B. T. Smith et al, Eispack Guide, Lecture Notes in 
   Computer Science, Springer-Verlag, (1976).

   The function finds the eigenvalues of a symmetric tridiagonal
   matrix by the implicit QL method.


   Arguments 
   ---------

   (input)
   n      order of the symmetric tridiagonal matrix                   
   d      contains the diagonal elements of the input matrix           
   e      contains the subdiagonal elements of the input matrix in its
          last n-1 positions.  e[0] is arbitrary	             

   (output)
   d      contains the eigenvalues in ascending order.  if an error
            exit is made, the eigenvalues are correct and ordered for
            indices 0,1,...ierr, but may not be the smallest eigenvalues.
   e      has been destroyed.					    
   ierr   set to zero for normal return, j if the j-th eigenvalue has
            not been determined after 30 iterations.		    

   Functions used
   --------------

   UTILITY	svd_fsign
   MISC		svd_pythag

 ***********************************************************************/

void imtqlb(long n, double d[], double e[], double bnd[])

{
   long last, l, m, i, iteration;

   /* various flags */
   long exchange, convergence, underflow;	

   double b, test, g, r, s, c, p, f;

   if (n == 1) return;
   ierr = 0;
   bnd[0] = 1.0;
   last = n - 1;
   for (i = 1; i < n; i++) {
      bnd[i] = 0.0;
      e[i-1] = e[i];
   }
   e[last] = 0.0;
   for (l = 0; l < n; l++) {
      iteration = 0;
      while (iteration <= 30) {
	 for (m = l; m < n; m++) {
	    convergence = FALSE;
	    if (m == last) break;
	    else {
	       test = fabs(d[m]) + fabs(d[m+1]);
	       if (test + fabs(e[m]) == test) convergence = TRUE;
	    }
	    if (convergence) break;
	 }
	    p = d[l]; 
	    f = bnd[l]; 
	 if (m != l) {
	    if (iteration == 30) {
	       ierr = l;
	       return;
	    }
	    iteration += 1;
	    /*........ form shift ........*/
	    g = (d[l+1] - p) / (2.0 * e[l]);
	    r = svd_pythag(g, 1.0);
	    g = d[m] - p + e[l] / (g + svd_fsign(r, g));
	    s = 1.0;
	    c = 1.0;
	    p = 0.0;
	    underflow = FALSE;
	    i = m - 1;
	    while (underflow == FALSE && i >= l) {
	       f = s * e[i];
	       b = c * e[i];
	       r = svd_pythag(f, g);
	       e[i+1] = r;
	       if (r == 0.0) underflow = TRUE;
	       else {
		  s = f / r;
		  c = g / r;
		  g = d[i+1] - p;
		  r = (d[i] - g) * s + 2.0 * c * b;
		  p = s * r;
		  d[i+1] = g + p;
		  g = c * r - b;
		  f = bnd[i+1];
		  bnd[i+1] = s * bnd[i] + c * f;
		  bnd[i] = c * bnd[i] - s * f;
		  i--;
	       }
	    }       /* end while (underflow != FALSE && i >= l) */
	    /*........ recover from underflow .........*/
	    if (underflow) {
	       d[i+1] -= p;
	       e[m] = 0.0;
	    }
	    else {
	       d[l] -= p;
	       e[l] = g;
	       e[m] = 0.0;
	    }
	 } 		       		   /* end if (m != l) */
	 else {

            /* order the eigenvalues */
	    exchange = TRUE;
	    if (l != 0) {
	       i = l;
	       while (i >= 1 && exchange == TRUE) {
	          if (p < d[i-1]) {
		     d[i] = d[i-1];
		     bnd[i] = bnd[i-1];
	             i--;
	          }
	          else exchange = FALSE;
	       }
	    }
	    if (exchange) i = 0;
	    d[i] = p;
	    bnd[i] = f; 
	    iteration = 31;
	 }
      }			       /* end while (iteration <= 30) */
   }				   /* end for (l=0; l<n; l++) */
   return;
}						  /* end main */

/***********************************************************************
 *                                                                     *
 *				imtql2()			       *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   imtql2() is a translation of a Fortran version of the Algol
   procedure IMTQL2, Num. Math. 12, 377-383(1968) by Martin and 
   Wilkinson, as modified in Num. Math. 15, 450(1970) by Dubrulle.  
   Handbook for Auto. Comp., vol.II-Linear Algebra, 241-248(1971).  
   See also B. T. Smith et al, Eispack Guide, Lecture Notes in 
   Computer Science, Springer-Verlag, (1976).

   This function finds the eigenvalues and eigenvectors of a symmetric
   tridiagonal matrix by the implicit QL method.


   Arguments
   ---------

   (input)                                                             
   nm     row dimension of the symmetric tridiagonal matrix           
   n      order of the matrix                                        
   d      contains the diagonal elements of the input matrix        
   e      contains the subdiagonal elements of the input matrix in its
            last n-1 positions.  e[0] is arbitrary	             
   z      contains the identity matrix				    
                                                                   
   (output)                                                       
   d      contains the eigenvalues in ascending order.  if an error
            exit is made, the eigenvalues are correct but unordered for
            for indices 0,1,...,ierr.				   
   e      has been destroyed.					  
   z      contains orthonormal eigenvectors of the symmetric   
            tridiagonal (or full) matrix.  if an error exit is made,
            z contains the eigenvectors associated with the stored 
          eigenvalues.					
   ierr   set to zero for normal return, j if the j-th eigenvalue has
            not been determined after 30 iterations.		    


   Functions used
   --------------
   UTILITY	svd_fsign
   MISC		svd_pythag

 ***********************************************************************/
//typically, when called, nm==n, d and e are work arrays and z is the large array (s).

/*//This has now been placed in ritvec.c.
void imtql2(long nm, long n, double d[], double e[], double z[],dSpMem *sp, float *zf, fSpMem *spf)

{
  //d[] can be of size equal to number of ritz values stabilized.
  //modified by agb to include the option of using sparse memory - if z and zf is NULL. (can select between using floats and doubles depending on z or zf.
   long index, nnm, j, last, l, m, i, k, iteration, convergence, underflow;
   double b, test, g, r, s, c, p, f;
   clock_t starttime, endtime;
   //double *twocols=NULL;
   long k1,k2,k1off,k2off;
   int row;
   //long dmax=0;
   if (n == 1) return;
   ierr = 0;
   last = n - 1;
   for (i = 1; i < n; i++) e[i-1] = e[i];
   e[last] = 0.0;
   nnm = n * nm;
   printf("imtql2: n %ld, nm %ld, nnm %ld\n",n,nm,nnm);
   //if(z==NULL){
   //  twocols=(double*)malloc(sizeof(double)*2*nm);
   // }
   starttime = clock();
   for (l = 0; l < n; l++) {
      iteration = 0;
      printf("imtql2 computing evals/vecs %ld/%ld (gets a bit faster)       \r",l,n);
      fflush(NULL);
      // look for small sub-diagonal element 
      while (iteration <= 30) {
	 for (m = l; m < n; m++) {
	    convergence = FALSE;
	    if (m == last) break;
	    else {
	       test = fabs(d[m]) + fabs(d[m+1]);
	       if (test + fabs(e[m]) == test) convergence = TRUE;
	    }
	    if (convergence) break;
	 }
	 if (m != l) {

	   // set error -- no convergence to an eigenvalue after
	   // 30 iterations.
	    if (iteration == 30) {
	       ierr = l;
	       return;
	    }
	    p = d[l]; 
	    iteration += 1;

	    // form shift 
	    g = (d[l+1] - p) / (2.0 * e[l]);
	    r = svd_pythag(g, 1.0);
	    g = d[m] - p + e[l] / (g + svd_fsign(r, g));
	    s = 1.0;
	    c = 1.0;
	    p = 0.0;
	    underflow = FALSE;
	    i = m - 1;
	    while (underflow == FALSE && i >= l) {
	       f = s * e[i];
	       b = c * e[i];
	       r = svd_pythag(f, g);
	       e[i+1] = r;
	       if (r == 0.0)
		 underflow = TRUE;
	       else {
		  s = f / r;
		  c = g / r;
		  g = d[i+1] - p;
		  r = (d[i] - g) * s + 2.0 * c * b;
		  p = s * r;
		  d[i+1] = g + p;
		  g = c * r - b;

		  // form vector 
		  if(z!=NULL){//double array version
		    for (k = 0; k < nnm; k += n) {
		      index = k + i;
		      f = z[index+1];
		      z[index+1] = s * z[index] + c * f;
		      z[index] = c * z[index] - s * f;
		    } 
		  }else if(zf!=NULL){//use float version
		    for (k = 0; k < nnm; k += n) {
		      index = k + i;
		      f = (float)zf[index+1];
		      zf[index+1] = (float)(s * zf[index] + c * f);
		      zf[index] = (float)(c * zf[index] - s * f);
		    } 
		  }else if(sp!=NULL){
		    k1=sp->indptr[i];
		    k2=sp->indptr[i+1];
		    k1off=0;
		    k2off=0;
		    while(k1<sp->indptr[i+1] && k2<sp->indptr[i+2]){
		      //printf("col %ld, k1=%ld<%d, k2=%ld<%d row %d %d full %d\n",i,k1,sp->indptr[i+1],k2,sp->indptr[i+2],sp->rowind[k1],sp->rowind[k2],sp->cnt==sp->ndata); 
		      if(sp->rowind[k1]<sp->rowind[k2]){//only the first col contains a value for this row...
			test=sp->data[k1];
			//note, f==0 (=smGet(sp,k1,i+1)).
			row=sp->rowind[k1];//save incase if gets deleted when replacing the data (ie if data==0).
			k1off=smReplaceData(sp,k1,i,c*test);//return 1 if replaced, 0 if removed.
			//Need to take into account where the current min value is.  If this is <k1, or <k2, will need to adjust accordingly since this value will be removed.
			k2off=smInsert(sp,row,i+1,s*test);//increment k2 counter if a new value inserted (we don't want to use it).
			//We now find the new values for k1 and k2.  Note, it is not just as simple as incrementing k1.  There are many things to take account of - whether data has been removed or inserted, where the current minimum value was, whether this was removed (if the array was full) etc.
			k1=smGetIndxForRow(sp,row+1,i,k1);
			k2=smGetIndxForRow(sp,row+1,i+1,k2);
			//k1+=k1off;//only increment this if it hasn't been deleted.
			//k2+=k2off+k1off-1;//decrement k2 iff the element from col i was removed (ie if c*test==0).
		      }else if(sp->rowind[k1]>sp->rowind[k2]){//only 2nd col contains a value for this row...
			//note, test==0 (=smGet(sp,k2,i)).
			f=sp->data[k2];
			row=sp->rowind[k2];
			k2off=smReplaceData(sp,k2,i+1,c*f);
			k1off=smInsert(sp,row,i,-s*f);
			k1=smGetIndxForRow(sp,row+1,i,k1);
			k2=smGetIndxForRow(sp,row+1,i+1,k2);
			//k1+=k1off;
			//k2+=k2off+k1off;//increment if its own entry hasn't been deleted, and if a new value has been added to col i.
		      }else{//both cols contain a value for this row.
			test=sp->data[k1];
			f=sp->data[k2];
			row=sp->rowind[k1];
			k1off=smReplaceData(sp,k1,i,c*test-s*f);
			k2off=smReplaceData(sp,k2+k1off-1,i+1,s*test+c*f);
			k1=smGetIndxForRow(sp,row+1,i,k1);
			k2=smGetIndxForRow(sp,row+1,i+1,k2);
			//k1+=k1off;
			//k2+=k1off+k2off-1;
		      }
		    }
		    while(k1<sp->indptr[i+1]){//continue for k1...
		      //if we get here, means no entries in col i+1 are left
		      test=sp->data[k1];
		      row=sp->rowind[k1];
		      k1off=smReplaceData(sp,k1,i,c*test);
		      k2off=smInsert(sp,row,i+1,s*test);
		      k1=smGetIndxForRow(sp,row+1,i,k1);
		      k2=smGetIndxForRow(sp,row+1,i+1,k2);
		      //k1+=k1off;
		      //k2+=k2off+k1off-1;
		    }
		    while(k2<sp->indptr[i+2]){//continue for k2
		      //if get here, means no entries in col i were left.
		      f=sp->data[k2];
		      row=sp->rowind[k2];
		      k2off=smReplaceData(sp,k2,i+1,c*f);
		      k1off=smInsert(sp,row,i,-s*f);
		      k1=smGetIndxForRow(sp,row+1,i,k1);
		      k2=smGetIndxForRow(sp,row+1,i+1,k2);
		      //k1+=k1off;
		      //k2+=k2off+k1off;
		    }
		    
		    //for (k = 0; k < nm; k ++) {
		    //  f = smGet(sp,k,i+1);
		    //  test=smGet(sp,k,i);
		    //  smInsert(sp,k,i+1,s * test + c * f);
		    //  smInsert(sp,k,i  ,c * test - s * f);
		    //  } 
		    //for (k = 0; k < nnm; k += n) {
		    //  index = k + i;
		    //  f = smGet(sp,index+1,0);
		    //  test=smGet(sp,index,0);
		    //  smInsert(sp,index+1,0,s * test + c * f);
		    //  smInsert(sp,index,0, c * test - s * f);
		    //}
		  }else{
		    k1=spf->indptr[i];
		    k2=spf->indptr[i+1];
		    k1off=0;
		    k2off=0;
		    while(k1<spf->indptr[i+1] && k2<spf->indptr[i+2]){
		      //printf("col %ld, k1=%ld<%d, k2=%ld<%d row %d %d full %d\n",i,k1,sp->indptr[i+1],k2,sp->indptr[i+2],sp->rowind[k1],sp->rowind[k2],sp->cnt==sp->ndata); 
		      if(spf->rowind[k1]<spf->rowind[k2]){//only the first col contains a value for this row...
			test=(double)spf->data[k1];
			//note, f==0 (=smGet(sp,k1,i+1)).
			row=spf->rowind[k1];//save incase if gets deleted when replacing the data (ie if data==0).
			k1off=smReplaceDataFloat(spf,k1,i,(float)(c*test));//return 1 if replaced, 0 if removed.
			//Need to take into account where the current min value is.  If this is <k1, or <k2, will need to adjust accordingly since this value will be removed.
			k2off=smInsertFloat(spf,row,i+1,(float)(s*test));//increment k2 counter if a new value inserted (we don't want to use it).
			//We now find the new values for k1 and k2.  Note, it is not just as simple as incrementing k1.  There are many things to take account of - whether data has been removed or inserted, where the current minimum value was, whether this was removed (if the array was full) etc.
			k1=smGetIndxForRowFloat(spf,row+1,i,k1);
			k2=smGetIndxForRowFloat(spf,row+1,i+1,k2);
			//k1+=k1off;//only increment this if it hasn't been deleted.
			//k2+=k2off+k1off-1;//decrement k2 iff the element from col i was removed (ie if c*test==0).
		      }else if(spf->rowind[k1]>spf->rowind[k2]){//only 2nd col contains a value for this row...
			//note, test==0 (=smGet(sp,k2,i)).
			f=(double)spf->data[k2];
			row=spf->rowind[k2];
			k2off=smReplaceDataFloat(spf,k2,i+1,(float)(c*f));
			k1off=smInsertFloat(spf,row,i,(float)(-s*f));
			k1=smGetIndxForRowFloat(spf,row+1,i,k1);
			k2=smGetIndxForRowFloat(spf,row+1,i+1,k2);
			//k1+=k1off;
			//k2+=k2off+k1off;//increment if its own entry hasn't been deleted, and if a new value has been added to col i.
		      }else{//both cols contain a value for this row.
			test=(double)spf->data[k1];
			f=(double)spf->data[k2];
			row=spf->rowind[k1];
			k1off=smReplaceDataFloat(spf,k1,i,(float)(c*test-s*f));
			k2off=smReplaceDataFloat(spf,k2+k1off-1,i+1,(float)(s*test+c*f));
			k1=smGetIndxForRowFloat(spf,row+1,i,k1);
			k2=smGetIndxForRowFloat(spf,row+1,i+1,k2);
			//k1+=k1off;
			//k2+=k1off+k2off-1;
		      }
		    }
		    while(k1<spf->indptr[i+1]){//continue for k1...
		      //if we get here, means no entries in col i+1 are left
		      test=(double)spf->data[k1];
		      row=spf->rowind[k1];
		      k1off=smReplaceDataFloat(spf,k1,i,(float)(c*test));
		      k2off=smInsertFloat(spf,row,i+1,(float)(s*test));
		      k1=smGetIndxForRowFloat(spf,row+1,i,k1);
		      k2=smGetIndxForRowFloat(spf,row+1,i+1,k2);
		      //k1+=k1off;
		      //k2+=k2off+k1off-1;
		    }
		    while(k2<spf->indptr[i+2]){//continue for k2
		      //if get here, means no entries in col i were left.
		      f=(double)spf->data[k2];
		      row=spf->rowind[k2];
		      k2off=smReplaceDataFloat(spf,k2,i+1,(float)(c*f));
		      k1off=smInsertFloat(spf,row,i,(float)(-s*f));
		      k1=smGetIndxForRowFloat(spf,row+1,i,k1);
		      k2=smGetIndxForRowFloat(spf,row+1,i+1,k2);
		      //k1+=k1off;
		      //k2+=k2off+k1off;
		    }
		    


		  }
		  i--;
	       }
	    }   // end while (underflow != FALSE && i >= l) 
	    //........ recover from underflow .........
	    if (underflow) {
	       d[i+1] -= p;
	       e[m] = 0.0;
	    }
	    else {
	       d[l] -= p;
	       e[l] = g;
	       e[m] = 0.0;
	    }
	 }
	 else break;
      }		//...... end while (iteration <= 30) .........
   }		//...... end for (l=0; l<n; l++) .............
   endtime = clock();
   //SAFE_FREE(twocols);
   printf("imtql2 - eigenvec computation took %gs.  Now ordering\n",((double) (endtime - starttime)) / CLOCKS_PER_SEC);
   // order the eigenvalues 
   for (l = 1; l < n; l++) {
      i = l - 1;
      k = i;
      p = d[i];
      for (j = l; j < n; j++) {
	 if (d[j] < p) {
	    k = j;
	    p = d[j];
	 }
      }
      // ...and corresponding eigenvectors 
      if (k != i) {
	 d[k] = d[i];
	 d[i] = p;
	 if(z!=NULL){
	   for (j = 0; j < nnm; j += n) {
	     p = z[j+i];
	     z[j+i] = z[j+k];
	     z[j+k] = p;
	   }
	 }else if(zf!=NULL){
	   float pf;
	   for (j = 0; j < nnm; j += n) {
	     pf = zf[j+i];
	     zf[j+i] = zf[j+k];
	     zf[j+k] = pf;
	   }
	 }else if(sp!=NULL){
	   for (j = 0; j < nm; j ++) {
	     p = smGet(sp,j,i);
	     smInsert(sp,j,i, smGet(sp,j,k));
	     smInsert(sp,j,k, p);
	   }
	 }else{
	   float pf;
	   for (j = 0; j < nm; j ++) {
	     pf = smGetFloat(spf,j,i);
	     smInsertFloat(spf,j,i, smGetFloat(spf,j,k));
	     smInsertFloat(spf,j,k, pf);
	   }

	 }
      }   
   }
   return;
}		
*/


/***********************************************************************
 *                                                                     *
 *				machar()			       *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   This function is a partial translation of a Fortran-77 subroutine 
   written by W. J. Cody of Argonne National Laboratory.
   It dynamically determines the listed machine parameters of the
   floating-point arithmetic.  According to the documentation of
   the Fortran code, "the determination of the first three uses an
   extension of an algorithm due to M. Malcolm, ACM 15 (1972), 
   pp. 949-951, incorporating some, but not all, of the improvements
   suggested by M. Gentleman and S. Marovich, CACM 17 (1974), 
   pp. 276-277."  The complete Fortran version of this translation is
   documented in W. J. Cody, "Machar: a Subroutine to Dynamically 
   Determine Determine Machine Parameters," TOMS 14, December, 1988.


   Parameters reported 
   -------------------

   ibeta     the radix for the floating-point representation       
   it        the number of base ibeta digits in the floating-point
               significand					 
   irnd      0 if floating-point addition chops		      
             1 if floating-point addition rounds, but not in the 
                 ieee style					
             2 if floating-point addition rounds in the ieee style
             3 if floating-point addition chops, and there is    
                 partial underflow				
             4 if floating-point addition rounds, but not in the
                 ieee style, and there is partial underflow    
             5 if floating-point addition rounds in the ieee style,
                 and there is partial underflow                   
   machep    the largest negative integer such that              
                 1.0+float(ibeta)**machep .ne. 1.0, except that 
                 machep is bounded below by  -(it+3)          
   negeps    the largest negative integer such that          
                 1.0-float(ibeta)**negeps .ne. 1.0, except that 
                 negeps is bounded below by  -(it+3)	       

 ***********************************************************************/

void machar(long *ibeta, long *it, long *irnd, long *machep, long *negep) {

  volatile double beta, betain, betah, a, b, ZERO, ONE, TWO, temp, tempa,
    temp1;
  long i, itemp;
  
  ONE = (double) 1;
  TWO = ONE + ONE;
  ZERO = ONE - ONE;
  
  a = ONE;
  temp1 = ONE;
  while (temp1 - ONE == ZERO) {
    a = a + a;
    temp = a + ONE;
    temp1 = temp - a;
    b += a; /* to prevent icc compiler error */
  }
  b = ONE;
  itemp = 0;
  while (itemp == 0) {
    b = b + b;
    temp = a + b;
    itemp = (long)(temp - a);
  }
  *ibeta = itemp;
  beta = (double) *ibeta;
  
  *it = 0;
  b = ONE;
  temp1 = ONE;
  while (temp1 - ONE == ZERO) {
    *it = *it + 1;
    b = b * beta;
    temp = b + ONE;
    temp1 = temp - b;
  }
  *irnd = 0; 
  betah = beta / TWO; 
  temp = a + betah;
  if (temp - a != ZERO) *irnd = 1;
  tempa = a + beta;
  temp = tempa + betah;
  if ((*irnd == 0) && (temp - tempa != ZERO)) *irnd = 2;
  
  *negep = *it + 3;
  betain = ONE / beta;
  a = ONE;
  for (i = 0; i < *negep; i++) a = a * betain;
  b = a;
  temp = ONE - a;
  while (temp-ONE == ZERO) {
    a = a * beta;
    *negep = *negep - 1;
    temp = ONE - a;
  }
  *negep = -(*negep);
  
  *machep = -(*it) - 3;
  a = b;
  temp = ONE + a;
  while (temp - ONE == ZERO) {
    a = a * beta;
    *machep = *machep + 1;
    temp = ONE + a;
  }
  eps = a;
  return;
}

float macharfloat(long *ibeta, long *it, long *irnd, long *machep,long *negep){

  volatile float beta, betain, betah, a, b, ZERO, ONE, TWO, temp, tempa,
    temp1;
  long i, itemp;
  
  ONE = (float) 1;
  TWO = ONE + ONE;
  ZERO = ONE - ONE;
  
  a = ONE;
  temp1 = ONE;
  while (temp1 - ONE == ZERO) {
    a = a + a;
    temp = a + ONE;
    temp1 = temp - a;
    b += a; /* to prevent icc compiler error */
  }
  b = ONE;
  itemp = 0;
  while (itemp == 0) {
    b = b + b;
    temp = a + b;
    itemp = (long)(temp - a);
  }
  *ibeta = itemp;
  beta = (float) *ibeta;
  
  *it = 0;
  b = ONE;
  temp1 = ONE;
  while (temp1 - ONE == ZERO) {
    *it = *it + 1;
    b = b * beta;
    temp = b + ONE;
    temp1 = temp - b;
  }
  *irnd = 0; 
  betah = beta / TWO; 
  temp = a + betah;
  if (temp - a != ZERO) *irnd = 1;
  tempa = a + beta;
  temp = tempa + betah;
  if ((*irnd == 0) && (temp - tempa != ZERO)) *irnd = 2;
  
  *negep = *it + 3;
  betain = ONE / beta;
  a = ONE;
  for (i = 0; i < *negep; i++) a = a * betain;
  b = a;
  temp = ONE - a;
  while (temp-ONE == ZERO) {
    a = a * beta;
    *negep = *negep - 1;
    temp = ONE - a;
  }
  *negep = -(*negep);
  
  *machep = -(*it) - 3;
  a = b;
  temp = ONE + a;
  while (temp - ONE == ZERO) {
    a = a * beta;
    *machep = *machep + 1;
    temp = ONE + a;
  }
  //eps = a;
  return a;
}


/***********************************************************************
 *                                                                     *
 *                     store()                                         *
 *                                                                     *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   store() is a user-supplied function which, based on the input
   operation flag, stores to or retrieves from memory a vector.


   Arguments 
   ---------

   (input)
   n       length of vector to be stored or retrieved
   isw     operation flag:
	     isw = 1 request to store j-th Lanczos vector q(j)
	     isw = 2 request to retrieve j-th Lanczos vector q(j)
	     isw = 3 request to store q(j) for j = 0 or 1
	     isw = 4 request to retrieve q(j) for j = 0 or 1
   s	   contains the vector to be stored for a "store" request 

   (output)
   s	   contains the vector retrieved for a "retrieve" request 

   Functions used
   --------------

   BLAS		svd_dcopy

 ***********************************************************************/
//For EAGLE, might need to store these on disc...
void store(LanStoreStruct *lss,long n, long isw, long j, double *s) {
  //printf("called store %04ld %04ld %04ld\n", isw, j,n);
  if(lss->LanStore){//double...
    switch(isw) {
    case STORQ:
      if (!lss->LanStore[j + lss->maxll]) {
	//printf("Allocating LanStoreQ[%ld] size %ld doubles\n",j+MAXLL,n);
	if (!(lss->LanStore[j + lss->maxll] = svd_doubleArray(n, FALSE, "LanStore[j]")))
	  svd_fatalError("svdLAS2: failed to allocate LanStore[%d]", j + lss->maxll);
      }
      svd_dcopy(n, s, 1, lss->LanStore[j + lss->maxll], 1);
      break;
    case RETRQ:	
      if (!lss->LanStore[j + lss->maxll])
	svd_fatalError("svdLAS2: store (RETRQ) called on index %d (not allocated)", 
		       j + lss->maxll);
      svd_dcopy(n, lss->LanStore[j + lss->maxll], 1, s, 1);
      break;
    case STORP:	
      if (j >= lss->maxll) {
	svd_error("svdLAS2: store (STORP) called with j >= MAXLL");
	break;
      }
      if (!lss->LanStore[j]) {
	//printf("Allocating LanStoreP[%ld] size %ld doubles\n",j,n);
	if (!(lss->LanStore[j] = svd_doubleArray(n, FALSE, "LanStore[j]")))
	  svd_fatalError("svdLAS2: failed to allocate LanStore[%d]", j);
      }
      svd_dcopy(n, s, 1, lss->LanStore[j], 1);
      break;
    case RETRP:	
      if (j >= lss->maxll) {
	svd_error("svdLAS2: store (RETRP) called with j >= MAXLL");
	break;
      }
      if (!lss->LanStore[j])
	svd_fatalError("svdLAS2: store (RETRP) called on index %d (not allocated)", 
		       j);
      svd_dcopy(n, lss->LanStore[j], 1, s, 1);
      break;
    }
  }else{//storage in float format...
    switch(isw){
    case STORQ:
      if (!lss->LanStoreFloat[j + lss->maxll]) {
	//printf("Allocating LanStoreFloatQ[%ld] size %ld doubles\n",j+MAXLL,n);
	if (!(lss->LanStoreFloat[j + lss->maxll] = svd_floatArray(n, FALSE, "LanStoreFloat[j]")))
	  svd_fatalError("svdLAS2: failed to allocate LanStoreFloat[%d]", j + lss->maxll);
      }
      svd_dfcopy(n, s, 1, lss->LanStoreFloat[j + lss->maxll], 1);
      break;
    case RETRQ:	
      if (!lss->LanStoreFloat[j + lss->maxll])
	svd_fatalError("svdLAS2: store (RETRQ) called on index %d (not allocated)", 
		       j + lss->maxll);
      svd_fdcopy(n, lss->LanStoreFloat[j + lss->maxll], 1, s, 1);
      break;
    case STORP:	
      if (j >= lss->maxll) {
	svd_error("svdLAS2: store (STORP) called with j >= MAXLL");
	break;
      }
      if (!lss->LanStoreFloat[j]) {
	//printf("Allocating LanStoreFloatP[%ld] size %ld doubles\n",j,n);
	if (!(lss->LanStoreFloat[j] = svd_floatArray(n, FALSE, "LanStoreFloat[j]")))
	  svd_fatalError("svdLAS2: failed to allocate LanStoreFloat[%d]", j);
      }
      svd_dfcopy(n, s, 1, lss->LanStoreFloat[j], 1);
      break;
    case RETRP:	
      if (j >= lss->maxll) {
	svd_error("svdLAS2: store (RETRP) called with j >= MAXLL");
	break;
      }
      if (!lss->LanStoreFloat[j])
	svd_fatalError("svdLAS2: store (RETRP) called on index %d (not allocated)", 
		       j);
      svd_fdcopy(n, lss->LanStoreFloat[j], 1, s, 1);
      break;
    }
  }
  return;
}


void storeFloat(LanStoreStruct *lss,long n, long isw, long j, float *s) {
  //printf("called store %04ld %04ld %04ld\n", isw, j,n);
  if(lss->LanStore){//double...
    switch(isw) {
    case STORQ:
      if (!lss->LanStore[j + lss->maxll]) {
	//printf("Allocating LanStoreQ[%ld] size %ld doubles\n",j+MAXLL,n);
	if (!(lss->LanStore[j + lss->maxll] = svd_doubleArray(n, FALSE, "LanStore[j]")))
	  svd_fatalError("svdLAS2: failed to allocate LanStore[%d]", j + lss->maxll);
      }
      svd_fdcopy(n, s, 1, lss->LanStore[j + lss->maxll], 1);
      break;
    case RETRQ:	
      if (!lss->LanStore[j + lss->maxll])
	svd_fatalError("svdLAS2: store (RETRQ) called on index %d (not allocated)", 
		       j + lss->maxll);
      svd_dfcopy(n, lss->LanStore[j + lss->maxll], 1, s, 1);
      break;
    case STORP:	
      if (j >= lss->maxll) {
	svd_error("svdLAS2: store (STORP) called with j >= MAXLL");
	break;
      }
      if (!lss->LanStore[j]) {
	//printf("Allocating LanStoreP[%ld] size %ld doubles\n",j,n);
	if (!(lss->LanStore[j] = svd_doubleArray(n, FALSE, "LanStore[j]")))
	  svd_fatalError("svdLAS2: failed to allocate LanStore[%d]", j);
      }
      svd_fdcopy(n, s, 1, lss->LanStore[j], 1);
      break;
    case RETRP:	
      if (j >= lss->maxll) {
	svd_error("svdLAS2: store (RETRP) called with j >= MAXLL");
	break;
      }
      if (!lss->LanStore[j])
	svd_fatalError("svdLAS2: store (RETRP) called on index %d (not allocated)", 
		       j);
      svd_dfcopy(n, lss->LanStore[j], 1, s, 1);
      break;
    }
  }else{//storage in float format...
    switch(isw){
    case STORQ:
      if (!lss->LanStoreFloat[j + lss->maxll]) {
	//printf("Allocating LanStoreFloatQ[%ld] size %ld doubles\n",j+MAXLL,n);
	if (!(lss->LanStoreFloat[j + lss->maxll] = svd_floatArray(n, FALSE, "LanStoreFloat[j]")))
	  svd_fatalError("svdLAS2: failed to allocate LanStoreFloat[%d]", j + lss->maxll);
      }
      svd_fcopy(n, s, 1, lss->LanStoreFloat[j + lss->maxll], 1);
      break;
    case RETRQ:	
      if (!lss->LanStoreFloat[j + lss->maxll])
	svd_fatalError("svdLAS2: store (RETRQ) called on index %d (not allocated)", 
		       j + lss->maxll);
      svd_fcopy(n, lss->LanStoreFloat[j + lss->maxll], 1, s, 1);
      break;
    case STORP:	
      if (j >= lss->maxll) {
	svd_error("svdLAS2: store (STORP) called with j >= MAXLL");
	break;
      }
      if (!lss->LanStoreFloat[j]) {
	//printf("Allocating LanStoreFloatP[%ld] size %ld doubles\n",j,n);
	if (!(lss->LanStoreFloat[j] = svd_floatArray(n, FALSE, "LanStoreFloat[j]")))
	  svd_fatalError("svdLAS2: failed to allocate LanStoreFloat[%d]", j);
      }
      svd_fcopy(n, s, 1, lss->LanStoreFloat[j], 1);
      break;
    case RETRP:	
      if (j >= lss->maxll) {
	svd_error("svdLAS2: store (RETRP) called with j >= MAXLL");
	break;
      }
      if (!lss->LanStoreFloat[j])
	svd_fatalError("svdLAS2: store (RETRP) called on index %d (not allocated)", 
		       j);
      svd_fcopy(n, lss->LanStoreFloat[j], 1, s, 1);
      break;
    }
  }
  return;
}


/*
void storeFloat(long n, long isw, long j, float *s) {
  //printf("called store %04ld %04ld %04ld\n", isw, j,n);
  switch(isw) {
  case STORQ:
    if (!LanStoreFloat[j + MAXLL]) {
      //printf("Allocating LanStoreQ[%ld] size %ld doubles\n",j+MAXLL,n);
      if (!(LanStoreFloat[j + MAXLL] = svd_floatArray(n, FALSE, "LanStoreFloat[j]")))
        svd_fatalError("svdLAS2: failed to allocate LanStoreFloat[%d]", j + MAXLL);
    }
    svd_fcopy(n, s, 1, LanStoreFloat[j + MAXLL], 1);
    break;
  case RETRQ:	
    if (!LanStoreFloat[j + MAXLL])
      svd_fatalError("svdLAS2: storeFloat (RETRQ) called on index %d (not allocated)", 
                     j + MAXLL);
    svd_fcopy(n, LanStoreFloat[j + MAXLL], 1, s, 1);
    break;
  case STORP:	
    if (j >= MAXLL) {
      svd_error("svdLAS2: storeFloat (STORP) called with j >= MAXLL");
      break;
    }
    if (!LanStoreFloat[j]) {
      //printf("Allocating LanStoreP[%ld] size %ld doubles\n",j,n);
      if (!(LanStoreFloat[j] = svd_floatArray(n, FALSE, "LanStoreFloat[j]")))
        svd_fatalError("svdLAS2: failed to allocate LanStoreFloat[%d]", j);
    }
    svd_fcopy(n, s, 1, LanStoreFloat[j], 1);
    break;
  case RETRP:	
    if (j >= MAXLL) {
      svd_error("svdLAS2: storeFloat (RETRP) called with j >= MAXLL");
      break;
    }
    if (!LanStoreFloat[j])
      svd_fatalError("svdLAS2: storeFloat (RETRP) called on index %d (not allocated)", 
                     j);
    svd_fcopy(n, LanStoreFloat[j], 1, s, 1);
    break;
  }
  return;
}

void storeFloatUsingDouble(long n, long isw, long j, double *s) {
  //printf("called store %04ld %04ld %04ld\n", isw, j,n);
  switch(isw) {
  case STORQ:
    if (!LanStoreFloat[j + MAXLL]) {
      //printf("Allocating LanStoreQ[%ld] size %ld doubles\n",j+MAXLL,n);
      if (!(LanStoreFloat[j + MAXLL] = svd_floatArray(n, FALSE, "LanStoreFloat[j]")))
        svd_fatalError("svdLAS2: failed to allocate LanStoreFloat[%d]", j + MAXLL);
    }
    svd_dfcopy(n, s, 1, LanStoreFloat[j + MAXLL], 1);
    break;
  case RETRQ:	
    if (!LanStoreFloat[j + MAXLL])
      svd_fatalError("svdLAS2: storeFloat (RETRQ) called on index %d (not allocated)", 
                     j + MAXLL);
    svd_fdcopy(n, LanStoreFloat[j + MAXLL], 1, s, 1);
    break;
  case STORP:	
    if (j >= MAXLL) {
      svd_error("svdLAS2: storeFloat (STORP) called with j >= MAXLL");
      break;
    }
    if (!LanStoreFloat[j]) {
      //printf("Allocating LanStoreP[%ld] size %ld doubles\n",j,n);
      if (!(LanStoreFloat[j] = svd_floatArray(n, FALSE, "LanStoreFloat[j]")))
        svd_fatalError("svdLAS2: failed to allocate LanStoreFloat[%d]", j);
    }
    svd_dfcopy(n, s, 1, LanStoreFloat[j], 1);
    break;
  case RETRP:	
    if (j >= MAXLL) {
      svd_error("svdLAS2: storeFloat (RETRP) called with j >= MAXLL");
      break;
    }
    if (!LanStoreFloat[j])
      svd_fatalError("svdLAS2: storeFloat (RETRP) called on index %d (not allocated)", 
                     j);
    svd_fdcopy(n, LanStoreFloat[j], 1, s, 1);
    break;
  }
  return;
}

*/
