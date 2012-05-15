#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <fcntl.h>
#include <time.h>
#include <pthread.h>
#include "sparsemem.h"
#include "svdlib.h"
#include "svdutil.h"
#include "ritvec.h"

#ifdef SDOUBLE
  #define STYPE double
  #define SPTYPE dSpMem
#else
  #define STYPE float
  #define SPTYPE fSpMem

#endif

#ifdef STOREDOUBLE
  #define STORETYPE double
#else
  #define STORETYPE float
#endif

#ifdef VTDOUBLE
  #define WTYPE double
  #define COMPUTEGENINV computeGenInv_d
  #define COMPUTEGENINVTHREADED computeGenInv_d_t
  #define RUNGENINV runGenInv_d
#else
  #define WTYPE float
  #define COMPUTEGENINV computeGenInv_f
  #define COMPUTEGENINVTHREADED computeGenInv_f_t
  #define RUNGENINV runGenInv_f
#endif


int COMPUTEGENINV(DMat Ut, DMat Vt,float *S, long neig, long neigForGenInv, float fracEig, float minEig, ArrUnion *genInv,double minGIVal, int transposeGI, long evalStart);

extern void   store(LanStoreStruct *lss,long, long, long, double *);
extern void   storeFloat(LanStoreStruct *lss,long, long, long, float *);
extern void freeLanStore(LanStoreStruct *lss);
extern double *OPBTemp;

enum storeVals {STORQ = 1, RETRQ, STORP, RETRP};
#ifdef VTDOUBLE
 #ifdef STOREDOUBLE
  #ifdef SDOUBLE
   #define RITVEC ritvecddd
   #define IMTQL imtql2ddd
   #define CSSWAP csswapddd
   #define CSSWAPSTRUCT csswapstructddd
  #else
   #define RITVEC ritvecfdd
   #define IMTQL imtql2fdd
   #define CSSWAP csswapfdd
   #define CSSWAPSTRUCT csswapstructfdd
  #endif
 #else
  #ifdef SDOUBLE
   #define RITVEC ritvecdfd
   #define IMTQL imtql2dfd
   #define CSSWAP csswapdfd
   #define CSSWAPSTRUCT csswapstructdfd
  #else
   #define RITVEC ritvecffd
   #define IMTQL imtql2ffd
   #define CSSWAP csswapffd
   #define CSSWAPSTRUCT csswapstructffd
  #endif
 #endif
#else
 #ifdef STOREDOUBLE
  #ifdef SDOUBLE
   #define RITVEC ritvecddf
   #define IMTQL imtql2ddf
   #define CSSWAP csswapddf
   #define CSSWAPSTRUCT csswapstructddf
  #else
   #define RITVEC ritvecfdf
   #define IMTQL imtql2fdf
   #define CSSWAP csswapfdf
   #define CSSWAPSTRUCT csswapstructfdf
  #endif
 #else
  #ifdef SDOUBLE
   #define RITVEC ritvecdff
   #define IMTQL imtql2dff
   #define CSSWAP csswapdff
   #define CSSWAPSTRUCT csswapstructdff
  #else
   #define RITVEC ritvecfff
   #define IMTQL imtql2fff
   #define CSSWAP csswapfff
   #define CSSWAPSTRUCT csswapstructfff
  #endif
 #endif
#endif

typedef struct{
  STYPE *z;
  double c;
  double s;
  long start;
  long end;
  long nm;
}CSSWAPSTRUCT;

void CSSWAP(CSSWAPSTRUCT *css){//STYPE *z,double c,double s, long nm, long start,long end){
  long k;
  double f;
  for(k=css->start; k<css->end; k++){//transposed in memory - should give better performance...
    f = (double)css->z[k+css->nm];
    css->z[k+css->nm] =(STYPE)(css->s * css->z[k] + css->c * f);
    css->z[k] = (STYPE)(css->c * css->z[k] - css->s * f);
  }
}
void IMTQL(long nm, long n, double d[], STORETYPE *e, STYPE *z,SPTYPE *sp,int nthreads,long *ierr)

{
  //d[] can be of size equal to number of ritz values stabilized.
  //modified by agb to include the option of using sparse memory - if z and zf is NULL. (can select between using floats and doubles depending on z or zf.
   long nnm, j, last, l, m, i, k, iteration, convergence, underflow;
   long inm,knm;
   double b, test, g, r, s, c, p, f;
   clock_t starttime, endtime;
   //double *twocols=NULL;
   long k1,k2,k1off,k2off;
   int row;
   STYPE pp;//,sf,cf,ff;
   int ti;
   pthread_t *thread;
   CSSWAPSTRUCT *threadInfo;
   long nleft;
   thread=malloc(sizeof(pthread_t)*nthreads);
   threadInfo=malloc(sizeof(CSSWAPSTRUCT)*nthreads);
   for(ti=0; ti<nthreads; ti++){
     threadInfo[ti].z=z;
     threadInfo[ti].nm=nm;
   }
   //long dmax=0;
   if (n == 1) return;
   *ierr = 0;
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
      printf("imtql2 computing evals/vecs %ld/%ld        \r",l,n);
      fflush(NULL);
      /* look for small sub-diagonal element */
      while (iteration <= 30) {//number of iterations is typically 1-4, ie small.
	 for (m = l; m < n; m++) {
	    convergence = FALSE;
	    if (m == last) break;
	    else {
	       test = fabs(d[m]) + fabs(d[m+1]);
	       if (test + fabs((double)e[m]) == test) convergence = TRUE;
	    }
	    if (convergence) break;
	 }
	 if (m != l) {//m and l are constant in this block
	   //printf("iteration %ld\n",iteration);//The number of iterations required is typically 1-4.
	    /* set error -- no convergence to an eigenvalue after
	     * 30 iterations. */     
	    if (iteration == 30) {
	       *ierr = l;
	       return;
	    }
	    p = d[l]; 
	    iteration += 1;

	    /* form shift */
	    g = (d[l+1] - p) / (2.0 * (double)e[l]);
	    r = svd_pythag(g, 1.0);
	    g = d[m] - p + (double)e[l] / (g + svd_fsign(r, g));
	    s = 1.0;
	    c = 1.0;
	    p = 0.0;
	    underflow = FALSE;
	    i = m - 1;
	    while (underflow == FALSE && i >= l) {
	       //printf("i=%ld\n",i);
	       f = s * (double)e[i];
	       b = c * (double)e[i];
	       r = svd_pythag(f, g);
	       e[i+1] = (STORETYPE)r;
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

		  /* form vector */
		  if(z!=NULL){//double array version
/*		    for (k = 0; k < nnm; k += n) {
		      index = k + i;
		      f = (double)z[index+1];
		      z[index+1] =(STYPE)(s * z[index] + c * f);
		      z[index] = (STYPE)(c * z[index] - s * f);
		    } 
*/
		    //move it into an array so that can be parallelised.
		    threadInfo[0].c=c;
		    threadInfo[0].s=s;
		    threadInfo[0].start=i*nm;
		    threadInfo[0].end=threadInfo[0].start+n/nthreads;
		    nleft=n-(threadInfo[0].end-threadInfo[0].start);
		    for(ti=1; ti<nthreads; ti++){
		      threadInfo[ti].c=c;
		      threadInfo[ti].s=s;
		      threadInfo[ti].start=threadInfo[ti-1].end;
		      threadInfo[ti].end=threadInfo[ti].start+nleft/(nthreads-ti);
		      nleft-=threadInfo[ti].end-threadInfo[ti].start;
		      pthread_create(&thread[ti],NULL,(void*)CSSWAP,&threadInfo[ti]);
		    }
		    CSSWAP(&threadInfo[0]);
		    //CSSWAP(z,c,s,nm,i*nm,i*nm+n);
		    for(ti=1; ti<nthreads; ti++)
		      pthread_join(thread[ti],NULL);
		    /*
		    inm=i*nm+n;
		    for(k=i*nm; k<inm; k++){//transposed in memory - should give better performance...
		      //index=inm+k;
		      f = (double)z[k+nm];
		      z[k+nm] =(STYPE)(s * z[k] + c * f);
		      z[k] = (STYPE)(c * z[k] - s * f);
		      }*/
		  }else if(sp!=NULL){
		    /*
		    //expand the two columns into the vectors...
		    memset(twocols,0,sizeof(double)*2*nm);
		    for(k=sp->indptr[i];k<sp->indptr[i+1];k++)
		      twocols[sp->rowind[k]]=sp->data[k];
		    for(k=sp->indptr[i+1];k<sp->indptr[i+2];k++)
		      twocols[sp->rowind[k]+nm]=sp->data[k];
		    //perform the calc
		    for(k=0; k<nm; k++){
		      f=twocols[k+nm];
		      f = smGet(sp,k,i+1);
		      test=twocols[k];
		      twocols[k+nm]=s*test+c*f;
		      twocols[k]=c*test-s*f;
		    }
		    //now put back into the sparse memory.
		    for(k=0; k<nm; k++){
		      smInsert(sp,k,i+1,twocols[k+nm]);
		      smInsert(sp,k,i,twocols[k]);
		    }
		    */
		    k1=sp->indptr[i];
		    k2=sp->indptr[i+1];
		    k1off=0;
		    k2off=0;
		    while(k1<sp->indptr[i+1] && k2<sp->indptr[i+2]){
		      //printf("col %ld, k1=%ld<%d, k2=%ld<%d row %d %d full %d\n",i,k1,sp->indptr[i+1],k2,sp->indptr[i+2],sp->rowind[k1],sp->rowind[k2],sp->cnt==sp->ndata); 
		      if(sp->rowind[k1]<sp->rowind[k2]){//only the first col contains a value for this row...
			//note, f==0 (=smGet(sp,k1,i+1)).
			row=sp->rowind[k1];//save incase if gets deleted when replacing the data (ie if data==0).
#ifdef SDOUBLE
			test=sp->data[k1];
			k1off=smReplaceData(sp,k1,i,c*test);//return 1 if replaced, 0 if removed.
			k2off=smInsert(sp,row,i+1,s*test);//increment k2 counter if a new value inserted (we don't want to use it).
			//We now find the new values for k1 and k2.  Note, it is not just as simple as incrementing k1.  There are many things to take account of - whether data has been removed or inserted, where the current minimum value was, whether this was removed (if the array was full) etc.
			k1=smGetIndxForRow(sp,row+1,i,k1);
			k2=smGetIndxForRow(sp,row+1,i+1,k2);
#else
			test=(double)sp->data[k1];
			k1off=smReplaceDataFloat(sp,k1,i,(float)(c*test));//return 1 if replaced, 0 if removed.
			k2off=smInsertFloat(sp,row,i+1,(float)(s*test));//increment k2 counter if a new value inserted (we don't want to use it).
			k1=smGetIndxForRowFloat(sp,row+1,i,k1);
			k2=smGetIndxForRowFloat(sp,row+1,i+1,k2);
#endif
			//k1+=k1off;//only increment this if it hasn't been deleted.
			//k2+=k2off+k1off-1;//decrement k2 iff the element from col i was removed (ie if c*test==0).
		      }else if(sp->rowind[k1]>sp->rowind[k2]){//only 2nd col contains a value for this row...
			//note, test==0 (=smGet(sp,k2,i)).
			row=sp->rowind[k2];
#ifdef SDOUBLE
			f=sp->data[k2];
			k2off=smReplaceData(sp,k2,i+1,c*f);
			k1off=smInsert(sp,row,i,-s*f);
			k1=smGetIndxForRow(sp,row+1,i,k1);
			k2=smGetIndxForRow(sp,row+1,i+1,k2);
#else
			f=(double)sp->data[k2];
			k2off=smReplaceDataFloat(sp,k2,i+1,(float)(c*f));
			k1off=smInsertFloat(sp,row,i,(float)(-s*f));
			k1=smGetIndxForRowFloat(sp,row+1,i,k1);
			k2=smGetIndxForRowFloat(sp,row+1,i+1,k2);
#endif
			//k1+=k1off;
			//k2+=k2off+k1off;//increment if its own entry hasn't been deleted, and if a new value has been added to col i.
		      }else{//both cols contain a value for this row.
			row=sp->rowind[k1];
#ifdef SDOUBLE
			test=sp->data[k1];
			f=sp->data[k2];
			k1off=smReplaceData(sp,k1,i,c*test-s*f);
			k2off=smReplaceData(sp,k2+k1off-1,i+1,s*test+c*f);
			k1=smGetIndxForRow(sp,row+1,i,k1);
			k2=smGetIndxForRow(sp,row+1,i+1,k2);
#else
			test=(double)sp->data[k1];
			f=(double)sp->data[k2];
			k1off=smReplaceDataFloat(sp,k1,i,(float)(c*test-s*f));
			k2off=smReplaceDataFloat(sp,k2+k1off-1,i+1,(float)(s*test+c*f));
			k1=smGetIndxForRowFloat(sp,row+1,i,k1);
			k2=smGetIndxForRowFloat(sp,row+1,i+1,k2);
#endif
			//k1+=k1off;
			//k2+=k1off+k2off-1;
		      }
		    }
		    while(k1<sp->indptr[i+1]){//continue for k1...
		      //if we get here, means no entries in col i+1 are left
		      row=sp->rowind[k1];
#ifdef SDOUBLE
		      test=sp->data[k1];
		      k1off=smReplaceData(sp,k1,i,c*test);
		      k2off=smInsert(sp,row,i+1,s*test);
		      k1=smGetIndxForRow(sp,row+1,i,k1);
		      k2=smGetIndxForRow(sp,row+1,i+1,k2);
#else
		      test=(double)sp->data[k1];
		      k1off=smReplaceDataFloat(sp,k1,i,(float)(c*test));
		      k2off=smInsertFloat(sp,row,i+1,(float)(s*test));
		      k1=smGetIndxForRowFloat(sp,row+1,i,k1);
		      k2=smGetIndxForRowFloat(sp,row+1,i+1,k2);
#endif
		      //k1+=k1off;
		      //k2+=k2off+k1off-1;
		    }
		    while(k2<sp->indptr[i+2]){//continue for k2
		      //if get here, means no entries in col i were left.
		      row=sp->rowind[k2];
#ifdef SDOUBLE
		      f=sp->data[k2];
		      k2off=smReplaceData(sp,k2,i+1,c*f);
		      k1off=smInsert(sp,row,i,-s*f);
		      k1=smGetIndxForRow(sp,row+1,i,k1);
		      k2=smGetIndxForRow(sp,row+1,i+1,k2);
#else
		      f=(double)sp->data[k2];
		      k2off=smReplaceDataFloat(sp,k2,i+1,(float)(c*f));
		      k1off=smInsertFloat(sp,row,i,(float)(-s*f));
		      k1=smGetIndxForRowFloat(sp,row+1,i,k1);
		      k2=smGetIndxForRowFloat(sp,row+1,i+1,k2);
#endif
		      //k1+=k1off;
		      //k2+=k2off+k1off;
		    }
		    
		    /*
		    for (k = 0; k < nm; k ++) {
		      f = smGet(sp,k,i+1);
		      test=smGet(sp,k,i);
		      smInsert(sp,k,i+1,s * test + c * f);
		      smInsert(sp,k,i  ,c * test - s * f);
		      }*/ 
		    /*for (k = 0; k < nnm; k += n) {
		      index = k + i;
		      f = smGet(sp,index+1,0);
		      test=smGet(sp,index,0);
		      smInsert(sp,index+1,0,s * test + c * f);
		      smInsert(sp,index,0, c * test - s * f);
		    } */
		  }
		  i--;
	       }
	    }   /* end while (underflow != FALSE && i >= l) */
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
	 }
	 else break;
      }		/*...... end while (iteration <= 30) .........*/
   }		/*...... end for (l=0; l<n; l++) .............*/
   endtime = clock();
   //SAFE_FREE(twocols);
   printf("imtql2 - eigenvec computation took %gs.  Now ordering\n",((double) (endtime - starttime)) / CLOCKS_PER_SEC);
   SAFE_FREE(threadInfo);
   SAFE_FREE(thread);
   /* order the eigenvalues */
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
      /* ...and corresponding eigenvectors */
      if (k != i) {
	 d[k] = d[i];
	 d[i] = p;
	 if(z!=NULL){
/*	   for (j = 0; j < nnm; j += n) {
	     pp = z[j+i];
	     z[j+i] = z[j+k];
	     z[j+k] = pp;
	   }
*/
	   inm=i*nm;
	   knm=k*nm;
	   for(j=0; j<n; j++){//transposed version (should give better perofrmance).
	     pp=z[inm+j];
	     z[inm+j]=z[knm+j];
	     z[knm+j]=pp;
	   }
	 }else{// if(sp!=NULL){
	   for (j = 0; j < nm; j ++) {
#ifdef SDOUBLE
	     pp = smGet(sp,j,i);
	     smInsert(sp,j,i, smGet(sp,j,k));
	     smInsert(sp,j,k, pp);
#else
	     pp = smGetFloat(sp,j,i);
	     smInsertFloat(sp,j,i, smGetFloat(sp,j,k));
	     smInsertFloat(sp,j,k, pp);
#endif
	   }
	 }
      }   
   }
   return;
}		/*...... end main ............................*/


long RITVEC(long n, SMat A, SVDRec R, double kappa, double *ritz, double *bnd, 
            double *alf, double *bet, STORETYPE *w2, LanStoreStruct *lss,long steps, long neig, long size, ArrUnion *genInv,long neigForGenInv,float minEig,float fracEig,double minGIVal,int transposeGI,int prepareForGenInv,int considerSwap,int nthreads,long* ierr) {
  //size, added by agb is the maximum space that this array should take.  If the array needs to be larger than this, it will be implemented as a sparse memory matrix (sparsemem.h etc).
  //genInv, the dataspace for the generalised inverse added by agb.  Can be NULL.  If not NULL, V and U may be destroyed.  
  //neigForGenInv is the number of eigenvals to use in generalised inverse.
  //minEig is the value of minimum eigenval that will be accepted, or zero.
  //fracEig is the value of min eigenval as a fraction of max eigenval that will be accepted or zero.
  //If fracEig is set, this will overwrite minEig.
  //If minEig is set, this will overwrite neigForGenInv.
  //floatS, if set, will mean that float arrays are used for s or sp.  Otherwise, they will be double.
  //if considerSwap is set, will take care about writing to noncontiguous memory locations, possibly allocating extra memory if necessary.
  long js, jsq, i, k, /*size,*/ id2, tmp, nsig, x,xoffset;
  double *xv2, tmp0; //*w1 = R->Vt->value[0];
  int pagesize=0,ninpage=0;
  long xmod,xpos=0;
  float tmp1;
  WTYPE *w1;//double or float
  WTYPE *uttmp=NULL;//double or float
  WTYPE *tmpUtMem=NULL;
  STYPE *s=NULL;
  SPTYPE *sp=NULL;
  clock_t starttime;
  //dSpMem *sp=NULL;
  //fSpMem *spf=NULL;
  //float *sf=NULL;
  int xx,yy,allocatedUt=0;
  js = steps + 1;
  jsq = js * js;//agb - this could be large for EAGLE - approx nact**2. - ie about 8GB since double.
  /*size = sizeof(double) * n;*/
  printf("ritvec js=%ld, n=%ld\n",js,n);
  xoffset=n-neig;
  if(size==0 || jsq<size){
#ifdef SDOUBLE
    printf("Allocating large array s, size %ld elements (double)\n",jsq);
    s = svd_doubleArray(jsq, TRUE, "ritvec: s");//agb - large memory.
#else
    printf("Allocating large array s, size %ld elements (float)\n",jsq);
    s=svd_floatArray(jsq,TRUE,"ritvec: s float");
#endif
  }else{
    printf("Allocating sparse memory structure with %ld elements for a %ld x %ld array\n",size,js,js);
#ifdef SDOUBLE
    sp=smNewSparseMem(size,js,js);
#else
    sp=smNewSparseMemFloat(size,js,js);
#endif
  }
  xv2 = svd_doubleArray(n, FALSE, "ritvec: xv2");
  
  /* initialize s to an identity matrix */
  if(s!=NULL){
    for (i = 0; i < jsq; i+= (js+1))
      s[i] = 1.0;
  }else{
    for(i=0; i<js; i++){
#ifdef SDOUBLE
      smInsert(sp,i,i,1.0);
#else
      smInsertFloat(sp,i,i,1.);
#endif
    }
  }
  //svd_dcopy(js, alf, 1, w1, -1);//agb - w1=Vt->value[0].  copy from alf to w1 in reverse order (I think). alf is array of diagonal elements of the tridiagonal matrix T, w1 is the first column of Vt.  Agb replaced this with vecReverse.
  vecReverse(js,alf);//agb added this to use alf instead of w1 in imtql2.
#ifdef STOREDOUBLE
  svd_dcopy(steps, &bet[1], 1, &w2[1], -1);
#else
  svd_dfcopy(steps, &bet[1], 1, &w2[1], -1);
#endif
  /* on return from imtql2(), w1 contains eigenvalues in ascending 
   * order and s contains the corresponding eigenvectors */
  printf("Starting imtql2 (ierr=%ld)\n",*ierr);
  IMTQL(js, js, alf, w2, s,sp,nthreads,ierr);//w1=Vt->value[0] - agb - alf used to be w1.
  printf("Done imtql2 (ierr=%ld)\n",*ierr);
  if (*ierr) return 0;//agb - ierr is a global variable- so not thread safe.
  starttime=clock();
  /*fwrite((char *)&n, sizeof(n), 1, fp_out2);
    fwrite((char *)&js, sizeof(js), 1, fp_out2);
    fwrite((char *)&kappa, sizeof(kappa), 1, fp_out2);*/
  /*id = 0;*/
  nsig = 0;
  x = 0;
  //id2 = jsq - js;
  id2=js-1;//change for transposed version.
  xx=0;//added by agb for 2D s array...
  for (k = 0; k < js; k++) {
    printf("Iteration %ld/%ld\r",k,js);
    fflush(NULL);
    tmp = id2;
    yy=js-1;//added by agb for 2D s array...
    if (bnd[k] <= kappa * fabs(ritz[k]) && k > js-neig-1) {
      if (--x < 0)
	x = R->d - 1;
      //Note, there is no reason to have w1 as Vt.  It could be an allocated array, and the result at the end of each iteration saved to disk, in eg Vt_x.fits where x is the row number. (note in python Vt is Ut since transposed.
      if(x>=xoffset){
#ifdef VTDOUBLE
	w1 = R->Vt->value[x-xoffset];
#else
	w1 = R->Vt->valuef[x-xoffset];
#endif
      }else{
	printf("Error: las2.c - x<xoffset (%ld<%ld)\n",x,xoffset);
#ifdef VTDOUBLE
	w1=R->Vt->value[0];
#else
	w1=R->Vt->valuef[0];
#endif
      }
      memset(w1,0,sizeof(WTYPE)*n);//w1=Vt->value[x]

      //for (i = 0; i < n; i++) 
      //w1[i] = 0.0;//w1=Vt->value[x]
      if(s!=NULL){
	for (i = 0; i < js; i++) {
#ifdef STOREDOUBLE
	  store(lss,n, RETRQ, i, w2);//retrieve w2 from the store.  This could be slow for a disk store.
  #ifdef VTDOUBLE
	  svd_daxpy(n,(double)s[tmp], w2, 1, w1, 1);//w1=s[tmp]*w2+w1 (w1=Vt->value[x],w2 are vectors) All of s gets accessed, but some values are much smaller than others...
  #else
	  svd_daxpydf(n,(double)s[tmp], w2, 1, w1, 1);
  #endif
#else
	  storeFloat(lss,n, RETRQ, i, w2);//retrieve w2 from the store.  This could be slow for a disk store.
  #ifdef VTDOUBLE
	  svd_daxpyfd(n,(double)s[tmp], w2, 1, w1, 1);//w1=s[tmp]*w2+w1 (w1=Vt->value[x],w2 are vectors) All of s gets accessed, but some values are much smaller than others...
  #else
	  svd_daxpyff(n,(double)s[tmp], w2, 1, w1, 1);
  #endif
#endif
	  //tmp -= js;
	  tmp--;//changed for transposed version.
	}
      }else{
	long e=sp->indptr[xx+1];
	for(i=sp->indptr[xx]; i<e; i++){
	  //row=sp->rowind[i];
	  //val=sp->data[i];
#ifdef STOREDOUBLE
	  store(lss,n,RETRQ,js-1-sp->rowind[i],w2);
  #ifdef VTDOUBLE
	  svd_daxpy(n,(double)sp->data[i],w2,1,w1,1);//W1=Vt->value[x]
  #else
	  svd_daxpydf(n,(double)sp->data[i],w2,1,w1,1);//W1=Vt->value[x]
  #endif
#else
	  storeFloat(lss,n,RETRQ,js-1-sp->rowind[i],w2);
  #ifdef VTDOUBLE
	  svd_daxpyfd(n,(double)sp->data[i],w2,1,w1,1);//W1=Vt->value[x]
  #else
	  svd_daxpyff(n,(double)sp->data[i],w2,1,w1,1);//W1=Vt->value[x]
  #endif
#endif
	}
      }

      //agb - could we compute R->S (evals) here?
#ifdef VTDOUBLE
      svd_opb_t(A, w1, xv2, OPBTemp,nthreads);//compute Vt->value[x].AT.A, output in xv2.
      tmp0 = svd_ddot(n, w1, 1, xv2, 1);//dot Vt with xv2 (result is a scalar).
#else
      svd_opbf_t(A, w1, xv2, OPBTemp,nthreads);//compute Vt->value[x].AT.A, output in xv2.
      tmp0 = svd_fddot(n, w1, 1, xv2, 1);//dot Vt with xv2 (result is a scalar).
#endif
      tmp0=sqrt(tmp0);
      if(x>=xoffset)
	R->S[x-xoffset] = (float)tmp0;
      else
	printf("Error: las2.c - x<xoffset (%ld<%ld)\n",x,xoffset);
      //end R->S computation.
      nsig++;
    }
    //id2++;
    id2+=js;//changed for transposed version.
    xx++;//added by agb
  }
  printf("Finished intermediate section in %gs\n",((double)(clock()-starttime)/CLOCKS_PER_SEC));
  starttime=clock();
  SAFE_FREE(s);
#ifdef SDOUBLE
  smFreeSparseMem(sp);
#else
  smFreeSparseMemFloat(sp);
#endif
  freeLanStore(lss);//we now don't need to use the store any more...

  //Note, from this stage onwards, Vt is not changed - it could simply by loaded from disk as needed.

  //Also, note, this is the first point that UT is used.  So it could be allocated here - using eg the memory freed from LanStore... that is, if the user doesn't care about it at the end... ie they don't want to see it...
  if(R->Ut->value==NULL && R->Ut->valuef==NULL){
    //allocate it here...
    printf("Allocating space for Ut of size %ld x %ld\n",R->Ut->rows,R->Ut->cols);
    if(R->Vt->value==NULL){//as float
      R->Ut->valuef=(float **) malloc(R->Ut->rows * sizeof(float *));
      if (!R->Ut->valuef)
	printf("Fatal error: Unable to allocate Ut pointer array\n");
      R->Ut->valuef[0] = (float *) calloc(R->Ut->rows * R->Ut->cols, sizeof(float));
      if(!R->Ut->valuef[0])
	printf("Fatal error: Unable to allocate Ut array\n");
      for (i = 1; i < R->Ut->rows; i++)
	R->Ut->valuef[i]=R->Ut->valuef[i-1]+R->Ut->cols;
    }else{//as double
      R->Ut->value=(double **) malloc(R->Ut->rows * sizeof(double *));
      if (!R->Ut->value)
	printf("Fatal error: Unable to allocate Ut pointer array\n");
      R->Ut->value[0] = (double *) calloc(R->Ut->rows * R->Ut->cols, sizeof(double));
      if(!R->Ut->value[0])
	printf("Fatal error: Unable to allocate Ut array\n");
      for (i = 1; i < R->Ut->rows; i++)
	R->Ut->value[i]=R->Ut->value[i-1]+R->Ut->cols;
    }
    allocatedUt=1;
  }

  /* agb - removed Rotate the singular vectors and values. */
  /*rotateArray(R->Vt->value[0], R->Vt->rows * R->Vt->cols, x * R->Vt->cols);*/
  /* x is now the location of the highest singular value. */
  xoffset=0;
  R->d = svd_imin(R->d, nsig);
  //printf("rotating at x=%ld (this point moved to zero.  R->d now %ld\n",x,R->d);//577 downto 27.
  if((uttmp=malloc(sizeof(WTYPE)*R->Ut->cols))==NULL){//allocate a single row of Ut.
    printf("CRITICAL ERROR: unable to allocate uttmp\n");
  }
  if((genInv!=NULL || prepareForGenInv==1) && considerSwap==1){
    pagesize=getpagesize();
    if((tmpUtMem=malloc(pagesize*A->rows))==NULL){
      printf("unable to allocate tmpUtMem\n");
      considerSwap=0;
    }else{
      printf("Allocated tmpUtMem\n");
    }
    ninpage=pagesize/sizeof(WTYPE);//number of floats/doubles in a single page.
    xpos=0;
  }
  printf("Starting creation of U, pagesize=%d, ninpage=%d R->d=%d (considerSwap=%d)\n",pagesize,ninpage,R->d,considerSwap);
  for (x = 0; x < R->d; x++) {
    printf("Creating U %ld/%d\r",x,R->d);
    fflush(NULL);
    /* multiply by matrix B first */
/*    svd_opb(A, R->Vt->value[x+xoffset], xv2, OPBTemp);//compute Vt->value[x].AT.A, output in xv2.
    tmp0 = svd_ddot(n, R->Vt->value[x+xoffset], 1, xv2, 1);//dot Vt with xv2 (result is a scalar).
    //svd_daxpy(n, -tmp0, R->Vt->value[x+xoffset], 1, xv2, 1);//agb removed this.xv2=xv2+-tmp0*Vt->value[x]
    tmp0 = sqrt(tmp0);
    R->S[x] = tmp0;
*/    //xnorm = sqrt(svd_ddot(n, xv2, 1, xv2, 1));//agb removed this.
      
    /* multiply by matrix A to get (scaled) left s-vector */
    //This is the first place that U is used.
#ifdef VTDOUBLE
    w1=R->Vt->value[x+xoffset];
    svd_opa(A, w1,uttmp);// R->Ut->value[x]);//A.Vt->value[x] result stored in Ut->value[x].
    tmp1=(R->S[x]==0.?0.:1./R->S[x]);//R->S used to be tmp0.
    if(genInv!=NULL || prepareForGenInv==1){//note, this will destroy U - the scaling is for the first part of the generalised inverse, and stores rotated...
      tmp1*=tmp1;
      if(considerSwap==1){
	xmod=x%ninpage;
	svd_fddscalnew((long)A->rows,tmp1,uttmp,1,&(tmpUtMem[xmod]),ninpage);
	if(xmod==ninpage-1){
	  //we've filled the temporary memory, so now write it to the real array...
	  for(i=0; i<A->rows; i++){
	    memcpy(&(R->Ut->value[0][i*R->Ut->rows+xpos]),&(tmpUtMem[i*ninpage]),pagesize);
	  }
	  xpos+=ninpage;
	}
      }else{
	svd_fddscalnew((long)A->rows,tmp1,uttmp,1,&(R->Ut->value[0][x]),R->Ut->rows);
      }
    }else{
      svd_fddscalnew((long)A->rows,tmp1,uttmp,1,R->Ut->value[x],1);
    }
#else
    w1=R->Vt->valuef[x+xoffset];
    svd_opaff(A, w1,uttmp);// R->Ut->valuef[x]);//A.Vt->value[x] result stored in Ut->value[x].
    tmp1=(R->S[x]==0.?0.:1./R->S[x]);//R->S used to be tmp0.
    if(genInv!=NULL || prepareForGenInv==1){//note, this will destroy U - the scaling is for the first part of the generalised inverse.
      tmp1*=tmp1;
      if(considerSwap==1){
	xmod=x%ninpage;
	//printf("Doing scalnew xmod=%ld\n",xmod);
	svd_fffscalnew((long)A->rows,tmp1,uttmp,1,&(tmpUtMem[xmod]),ninpage);
	if(xmod==ninpage-1){
	  //we've filled the temporary memory, so now write it to the real array...
	  for(i=0; i<A->rows; i++){
	    //printf("memcpy %ld %ld %ld %ld\n",i,R->Ut->cols,xpos,xmod);
	    memcpy(&(R->Ut->valuef[0][i*R->Ut->rows+xpos]),&(tmpUtMem[i*ninpage]),pagesize);
	  }
	  xpos+=ninpage;
	}
      }else{
	svd_fffscalnew((long)A->rows,tmp1,uttmp,1,&(R->Ut->valuef[0][x]),R->Ut->rows);
      }
    }else{
      svd_fffscalnew((long)A->rows, tmp1,uttmp,1,R->Ut->valuef[x],1);
    }
#endif

    //if(genInv!=NULL){
      //perform the first part of the generalised inverse.  This is multiplying V by 1/S.
      //Note, this destroys V.
    //#ifdef VTDOUBLE
    //svd_fdscal(R->Vt->cols,tmp1,w1,1);
    //#else
    //svd_ffscal(R->Vt->cols,tmp1,w1,1);
    //#endif
    //}
    //xnorm *= tmp1;//agb removed this.
    //bnd[i] = xnorm;//agb removed this.
    //If wanted to compute the reconstruction matrix in place, here would be the place to do it.
      /*In python, the reconstructor is V.1/w.UT.
      However, here V and U are swapped around (since A is transposed).  So, I would want to do:
      U.1/w.VT
      The 1/w.VT simply scales each row [i] of VT by 1/w[i].  And R->S[x] is the w.
      This will have been done in the previous stage.
      After this, can then dot the cols of VT with the rows of U, ie the cols of UT.
      So, I could do something like (assuming Vt can be destroyed):
      svd_dscal(neig,1./tmp0,R->Vt->value[x+xoffset],1);
      save(Vt->value[x+xoffset]);//store it back in memory.
      //We need to dot every col of Ut by every col of Vt to get the reconmx (or rather, 
      //the neig first values of each col).
      */
  }
  //now, if considering swap, copy the last data into Ut...

  if((genInv!=NULL || prepareForGenInv==1) && considerSwap==1){
    if((R->d%ninpage)!=0){//otherwise it will have been done already...
      pagesize=(R->d%ninpage)*sizeof(WTYPE);//the number of bytes to copy for each row.
      for(i=0; i<A->rows; i++){
#ifdef VTDOUBLE
	memcpy(&(R->Ut->value[0][i*R->Ut->rows+xpos]),&(tmpUtMem[i*ninpage]),pagesize);
#else
	memcpy(&(R->Ut->valuef[0][i*R->Ut->rows+xpos]),&(tmpUtMem[i*ninpage]),pagesize);
#endif
      }
    }
    pagesize=getpagesize();
  }
  SAFE_FREE(uttmp);
  SAFE_FREE(tmpUtMem);
  if(genInv!=NULL || prepareForGenInv==1){
    printf("\nTransposing Vt for speed\n");
    for(i=0; i<R->Vt->rows; i++){
      printf("Done %ld/%ld\r",i,R->Vt->rows);
      fflush(NULL);
      for(k=i+1; k<R->Vt->cols; k++){
#ifdef VTDOUBLE
	tmp0=R->Vt->value[i][k];
	R->Vt->value[i][k]=R->Vt->value[k][i];
	R->Vt->value[k][i]=tmp0;
#else
	tmp1=R->Vt->valuef[i][k];
	R->Vt->valuef[i][k]=R->Vt->valuef[k][i];
	R->Vt->valuef[k][i]=tmp1;
#endif
      }
    }
  }

  //Initially, assuming all U and V stored in memory, for testing.  ie dot U with Vt.
  printf("Finished this section (%gs)\n",((double)(clock()-starttime)/CLOCKS_PER_SEC));
  if(genInv!=NULL){
    COMPUTEGENINV(R->Ut, R->Vt, R->S,neig,neigForGenInv,fracEig,minEig,genInv,minGIVal,transposeGI,0);
    /*
    starttime = clock();
    if(fracEig>0.)
      minEig=R->S[0]*fracEig;
    if(minEig>0.){
      i=0;
      while(i<neig && R->S[i]>=minEig)
	i++;
      if(i>1)
	neigForGenInv=i-1;
    }

    printf("Computing generalised inverse: ncent (ncols)=%ld, nact (nrows)=%ld, neig=%ld, neigForGenInv=%ld\n",R->Ut->cols,R->Vt->rows,neig,neigForGenInv);
    printf("Ut->cols %ld Ut->rows %ld Vt->cols %ld Vt->rows %ld R->d %d A->rows %u A->cols %u\n",R->Ut->cols,R->Ut->rows,R->Vt->cols,R->Vt->rows,R->d,A->rows,A->cols);
    //Note, Ut has been transposed previously, and so is now actually U.  This just affects how we access the array... (continuous instead of step).
    for(k=0; k<R->Ut->cols; k++){//ncent
      printf("Matrix multiply: done %ld/%ld      \r",k,R->Ut->cols);
      fflush(NULL);
#ifdef VTDOUBLE
      w1=&(R->Ut->value[0][R->Ut->rows*k]);//&(R->Ut->value[0][k]);
#else
      w1=&(R->Ut->valuef[0][R->Ut->rows*k]);//&(R->Ut->valuef[0][k]);
#endif
      for(i=0; i<R->Vt->cols; i++){//nact
	switch(genInv->typ){
	case 'd'://double array.
#ifdef VTDOUBLE
	  tmp0=(double)svd_ddot(neigForGenInv,w1,1,&(R->Vt->value[i][0]),1);
#else
	  tmp0=(double)svd_fdot(neigForGenInv,w1,1,&(R->Vt->valuef[i][0]),1);
#endif
	  if(fabs(tmp0)>minGIVal){
	    if(transposeGI)
	      genInv->data.ddata[k*R->Vt->cols+i]=tmp0;
	    else
	      genInv->data.ddata[i*R->Ut->cols+k]=tmp0;
	  }
	  break;
	case 'f'://float array
#ifdef VTDOUBLE
	  tmp1=(float)svd_ddot(neigForGenInv,w1,1,&(R->Vt->value[i][0]),1);
#else
	  tmp1=(float)svd_fdot(neigForGenInv,w1,1,&(R->Vt->valuef[i][0]),1);
#endif
	  if(fabsf(tmp1)>minGIVal){
	    if(transposeGI)
	      genInv->data.fdata[k*R->Vt->cols+i]=tmp1;
	    else
	      genInv->data.fdata[i*R->Ut->cols+k]=tmp1;
	  }
	  break;
	case 'D'://dSpMem object...
#ifdef VTDOUBLE
	  tmp0=(double)svd_ddot(neigForGenInv,w1,1,&(R->Vt->value[i][0]),1);
#else
	  tmp0=(double)svd_fdot(neigForGenInv,w1,1,&(R->Vt->valuef[i][0]),1);
#endif
	  if(fabs(tmp0)>minGIVal){
	    if(transposeGI)
	      smInsert(genInv->data.dSp,k,i,tmp0);
	    else
	      smInsert(genInv->data.dSp,i,k,tmp0);
	  }
	  break;
	case 'F':
#ifdef VTDOUBLE
	  tmp1=(float)svd_ddot(neigForGenInv,w1,1,&(R->Vt->value[i][0]),1);
#else
	  tmp1=(float)svd_fdot(neigForGenInv,w1,1,&(R->Vt->valuef[i][0]),1);
#endif
	  if(fabsf(tmp1)>minGIVal){
	    if(transposeGI)
	      smInsertFloat(genInv->data.fSp,k,i,tmp1);
	    else
	      smInsertFloat(genInv->data.fSp,i,k,tmp1);
	  }
	  break;
	}
      }
    }
    */
  }

  //for(i=0; i<neig; i++){//display Vt[:neig,9],Ut[:neig,9] - agrees with expected.
  // printf("%ld    %g    %g\n",i,R->Ut->value[0][9+i*R->Ut->cols],R->Vt->value[0][9+i*R->Vt->rows]);
  //}
  if(allocatedUt){
    if(R->Ut->value)
      SAFE_FREE(R->Ut->value[0]);
    if(R->Ut->valuef)
      SAFE_FREE(R->Ut->valuef[0]);
    SAFE_FREE(R->Ut->value);
    SAFE_FREE(R->Ut->valuef);
  }
  SAFE_FREE(xv2);
  return nsig;
}

//Only need 2 versions of computeGenInv (float or double for Vt).  So, only generate this in one setting
//of the other defines.
#ifdef SDOUBLE
#ifdef STOREDOUBLE
int COMPUTEGENINV(DMat Ut,DMat Vt,float *S, long neig, long neigForGenInv, float fracEig, float minEig, ArrUnion *genInv,double minGIVal, int transposeGI,long evalStart){
  //Previous to calling this, Ut and Vt should have been "prepared".  This means that Vt must have been transposed
  //and Ut must have had its rows scaled by the eigenvalues.  If the previous call to LAS2 had performed 
  //a generalised inverse, or if called with prepareGenInv set, then this will be prepared.  
  //transposeGI should probably be set if a dense array, but not set if a sparse array.
  //evalStart is the eigen value number to start at - all previous ones are ignored. - for testing, 0 by default.
  WTYPE *w1;//double or float
  double tmp0;
  float tmp1;
  long i,k;
  clock_t starttime;
  if(genInv!=NULL){
    starttime=clock();
    if(neigForGenInv<=0)
      neigForGenInv=neig;
    if(S!=NULL){
      if(fracEig>0.)
	minEig=S[0]*fracEig;
      if(minEig>0.){
	i=0;
	while(i<neig && S[i]>=minEig)
	  i++;
	if(i>1)
	  neigForGenInv=i-1;
      }
    }
    neigForGenInv-=evalStart;
    printf("Computing generalised inverse: ncent (ncols)=%ld, nact (nrows)=%ld, neig=%ld, neigForGenInv=%ld, evalStart=%ld\n",Ut->cols,Vt->rows,neig,neigForGenInv,evalStart);
    printf("Ut->cols %ld Ut->rows %ld Vt->cols %ld Vt->rows %ld\n",Ut->cols,Ut->rows,Vt->cols,Vt->rows);
    //Note, Ut has been transposed previously, and so is now actually U.  This just affects how we access the array... (continuous instead of step).
    for(k=0; k<Ut->cols; k++){//ncent
      printf("Matrix multiply: done %ld/%ld      \r",k,Ut->cols);
      fflush(NULL);
#ifdef VTDOUBLE
      w1=&(Ut->value[0][Ut->rows*k+evalStart]);//&(R->Ut->value[0][k]);
#else
      w1=&(Ut->valuef[0][Ut->rows*k+evalStart]);//&(R->Ut->valuef[0][k]);
#endif
      for(i=0; i<Vt->cols; i++){//nact
	switch(genInv->typ){
	case 'd'://double array.
#ifdef VTDOUBLE
	  tmp0=(double)svd_ddot(neigForGenInv,w1,1/*Ut->cols*/,&(Vt->value[i][evalStart]),1/*Vt->cols*/);
#else
	  tmp0=(double)svd_fdot(neigForGenInv,w1,1/*Ut->cols*/,&(Vt->valuef[i][evalStart]),1/*Vt->cols*/);
#endif
	  if(fabs(tmp0)>minGIVal){
	    if(transposeGI)
	      genInv->data.ddata[k*Vt->cols+i]=tmp0;
	    else
	      genInv->data.ddata[i*Ut->cols+k]=tmp0;
	  }
	  break;
	case 'f'://float array
#ifdef VTDOUBLE
	  tmp1=(float)svd_ddot(neigForGenInv,w1,1/*Ut->cols*/,&(Vt->value[i][evalStart]),1/*Vt->cols*/);
#else
	  tmp1=(float)svd_fdot(neigForGenInv,w1,1/*Ut->cols*/,&(Vt->valuef[i][evalStart]),1/*Vt->cols*/);
#endif
	  if(fabsf(tmp1)>minGIVal){
	    if(transposeGI)
	      genInv->data.fdata[k*Vt->cols+i]=tmp1;
	    else
	      genInv->data.fdata[i*Ut->cols+k]=tmp1;
	  }
	  break;
	case 'D'://dSpMem object...
#ifdef VTDOUBLE
	  tmp0=(double)svd_ddot(neigForGenInv,w1,1/*Ut->cols*/,&(Vt->value[i][evalStart]),1/*Vt->cols*/);
#else
	  tmp0=(double)svd_fdot(neigForGenInv,w1,1/*Ut->cols*/,&(Vt->valuef[i][evalStart]),1/*Vt->cols*/);
#endif
	  if(fabs(tmp0)>minGIVal){
	    if(transposeGI)
	      smInsert(genInv->data.dSp,k,i,tmp0);
	    else
	      smInsert(genInv->data.dSp,i,k,tmp0);
	  }
	  break;
	case 'F':
#ifdef VTDOUBLE
	  tmp1=(float)svd_ddot(neigForGenInv,w1,1/*Ut->cols*/,&(Vt->value[i][evalStart]),1/*Vt->cols*/);
#else
	  tmp1=(float)svd_fdot(neigForGenInv,w1,1/*Ut->cols*/,&(Vt->valuef[i][evalStart]),1/*Vt->cols*/);
#endif
	  if(fabsf(tmp1)>minGIVal){
	    if(transposeGI)
	      smInsertFloat(genInv->data.fSp,k,i,tmp1);
	    else
	      smInsertFloat(genInv->data.fSp,i,k,tmp1);
	  }
	  break;
	}
      }
    }
    printf("Finished computation of generalised inverse (took %gs)\n",((double)(clock()-starttime))/CLOCKS_PER_SEC);
  }
  return 0;
}
#endif //storedouble
#endif //sdouble

#ifdef SDOUBLE
#ifdef STOREDOUBLE//only define it twice...
//threaded versions
typedef struct{
  DMat Ut;
  DMat Vt;
  ArrUnion *genInv;
  double *minGIVal;
  long neigForGenInv;
  int transposeGI;
  long *cnt;
  int ncnt;
  long start;
  long end;
  int print;
  long evalStart;

}compGenInvStruct;

void RUNGENINV(compGenInvStruct *info){
  //Each thread calls this...
  WTYPE *w1;
  DMat Ut=info->Ut;
  DMat Vt=info->Vt;
  ArrUnion *genInv=info->genInv;
  long i,k;
  int j;
  double *minGIVal=info->minGIVal;
  long neigForGenInv=info->neigForGenInv;
  double tmp0;
  float tmp1;
  int transposeGI=info->transposeGI;
  long index;
  float *minGIValf;
  int ncnt=info->ncnt;
  long *cnt=info->cnt;
  long evalStart=info->evalStart;
  minGIValf=malloc(sizeof(float)*ncnt);
  for(j=0; j<ncnt; j++){
    cnt[j]=0;
    minGIValf[j]=(float)minGIVal[j];
  }
  //info->cnt=0;
  printf("ut->cols %ld ut->rows %ld vt->cols %ld vt->rows %ld start %ld end %ld\n",Ut->cols,Ut->rows,Vt->cols,Vt->rows,info->start,info->end);
  for(k=info->start; k<info->end; k++){//ncent
    if(info->print){
      printf("Matrix multiply: done %ld/%ld      \r",k,info->end);
      fflush(NULL);
    }
#ifdef VTDOUBLE
    w1=&(Ut->value[0][Ut->rows*k+evalStart]);//&(R->Ut->value[0][k]);
#else
    w1=&(Ut->valuef[0][Ut->rows*k+evalStart]);//&(R->Ut->valuef[0][k]);
#endif
    for(i=0; i<Vt->cols; i++){//nact
      if(transposeGI)
	index=k*Vt->cols+i;
      else
	index=i*Ut->cols+k;
      switch(genInv->typ){
	case 'd'://double array.
#ifdef VTDOUBLE
	  tmp0=(double)svd_ddot(neigForGenInv,w1,1/*Ut->cols*/,&(Vt->value[i][evalStart]),1/*Vt->cols*/);
#else
	  tmp0=(double)svd_fdot(neigForGenInv,w1,1/*Ut->cols*/,&(Vt->valuef[i][evalStart]),1/*Vt->cols*/);
#endif
	  genInv->data.ddata[index]=tmp0;
	  for(j=0; j<ncnt; j++){
	    if(fabs(tmp0)>minGIVal[j])
	      cnt[j]++;
	    else//skip the rest - should be in ascending order anyway.
	      j=ncnt;
	  }

	  //if(fabs(tmp0)>minGIVal){
	  //  info->cnt++;
	  //  genInv->data.ddata[index]=tmp0;
	  //}else
	  //  genInv->data.ddata[index]=0.;
	  
	  break;
	case 'f'://float array
#ifdef VTDOUBLE
	  tmp1=(float)svd_ddot(neigForGenInv,w1,1/*Ut->cols*/,&(Vt->value[i][evalStart]),1/*Vt->cols*/);
#else
	  tmp1=(float)svd_fdot(neigForGenInv,w1,1/*Ut->cols*/,&(Vt->valuef[i][evalStart]),1/*Vt->cols*/);
#endif
	  genInv->data.fdata[index]=tmp1;
	  for(j=0; j<ncnt; j++){
	    if(fabsf(tmp1)>minGIValf[j])
	      cnt[j]++;
	    else//skip the rest - should be in ascending order anyway.
	      j=ncnt;
	  }
	  
	  //if(fabsf(tmp1)>minGIVal){
	  //  genInv->data.fdata[index]=tmp1;
	  //  info->cnt++;
	  //}else
	  //  genInv->data.fdata[index]=0.;
	  //break;
      }
    }
  }
  free(minGIValf);
}

long COMPUTEGENINVTHREADED(DMat Ut,DMat Vt,float *S, long neig, long neigForGenInv, float fracEig, float minEig, ArrUnion *genInv,double *minGIVal,long *cnt, int ncnt,int transposeGI,int nthreads,long evalStart){
  //Note, that before calling this, U and V should have been prepared (ie U scaled by S, V transposed).

  //minGIVal is an array of values (ncnt values) for which we should
  //count the number of entries in the generalised inverse (to help
  //with sparsifying).  cnt is an empty array of same length,
  //which the results of these countings are placed into.  ncnt
  //is the number of different values to count for.
  //evalStart (for testing) is the number of initial evals/vecs to ignore...
  long i;
  clock_t starttime;
  pthread_t *thread;
  long nleft;
  //long cnt=0;
  int j;
  compGenInvStruct *threadInfo;
  if(genInv!=NULL && genInv->typ!='D' && genInv->typ!='F'){
    starttime=clock();
    if(neigForGenInv<=0)
      neigForGenInv=neig;
    if(S!=NULL){
      if(fracEig>0.)
	minEig=S[0]*fracEig;
      if(minEig>0.){
	i=0;
	while(i<neig && S[i]>=minEig)
	  i++;
	if(i>1)
	  neigForGenInv=i-1;
      }
    }
    neigForGenInv-=evalStart;
    printf("Computing generalised inverse: ncent (ncols)=%ld, nact (nrows)=%ld, neig=%ld, neigForGenInv=%ld, nthreads=%d, evalStart=%ld\n",Ut->cols,Vt->rows,neig,neigForGenInv,nthreads,evalStart);
    printf("Ut->cols %ld Ut->rows %ld Vt->cols %ld Vt->rows %ld\n",Ut->cols,Ut->rows,Vt->cols,Vt->rows);
    //Note, Ut has been transposed previously, and so is now actually U.  This just affects how we access the array... (continuous instead of step).

    thread=malloc(sizeof(pthread_t)*nthreads);
    threadInfo=malloc(sizeof(compGenInvStruct)*nthreads);
    threadInfo[0].print=1;//this thread prints the status...
    threadInfo[0].start=0;
    threadInfo[0].end=Ut->cols/nthreads;
    threadInfo[0].Ut=Ut;
    threadInfo[0].Vt=Vt;
    threadInfo[0].genInv=genInv;
    threadInfo[0].minGIVal=minGIVal;
    threadInfo[0].ncnt=ncnt;
    threadInfo[0].cnt=malloc(sizeof(long)*ncnt);
    threadInfo[0].neigForGenInv=neigForGenInv;
    threadInfo[0].transposeGI=transposeGI;
    threadInfo[0].evalStart=evalStart;
    nleft=Ut->cols-threadInfo[0].end;
    for(i=1; i<nthreads; i++){
      threadInfo[i].print=0;
      threadInfo[i].start=threadInfo[i-1].end;
      threadInfo[i].end=threadInfo[i].start+nleft/(nthreads-i);
      threadInfo[i].Ut=Ut;
      threadInfo[i].Vt=Vt;
      threadInfo[i].genInv=genInv;
      threadInfo[i].minGIVal=minGIVal;
      threadInfo[i].ncnt=ncnt;
      threadInfo[i].cnt=malloc(sizeof(long)*ncnt);
      threadInfo[i].neigForGenInv=neigForGenInv;
      threadInfo[i].transposeGI=transposeGI;
      threadInfo[i].evalStart=evalStart;
      nleft-=(threadInfo[i].end-threadInfo[i].start);
      printf("start %ld,end %ld\n",threadInfo[i].start,threadInfo[i].end);
      pthread_create(&thread[i],NULL,(void*)RUNGENINV,&threadInfo[i]);
    }
    RUNGENINV(&threadInfo[0]);
    printf("\nWaiting for threads\n");
    //cnt=threadInfo[0].cnt;
    for(j=0; j<ncnt; j++){
      cnt[j]=threadInfo[0].cnt[j];
      printf("cnt 0[%d]: %ld\n",j,threadInfo[0].cnt[j]);
    }
    for(i=1; i<nthreads; i++){
      pthread_join(thread[i],NULL);
      for(j=0; j<ncnt; j++){
	cnt[j]+=threadInfo[i].cnt[j];
	printf("cnt %ld[%d]: %ld\n",i,j,threadInfo[i].cnt[j]);
      }
    }
    printf("Done threads\n");
    for(j=0; j<ncnt; j++){
      printf("cnt[%d]=%ld\n",j,cnt[j]);
    }
    for(i=0; i<nthreads; i++){
      SAFE_FREE(threadInfo[i].cnt);
    }
    SAFE_FREE(thread);
    SAFE_FREE(threadInfo);
    printf("Cleaning up done\n");
  }
  return 0;
}

#endif //storedouble
#endif //sdouble

