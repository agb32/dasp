#ifndef RITVEC_H
#define RITVEC_H
/*This is the order...
#ifdef VTDOUBLE
 #ifdef STOREDOUBLE
  #ifdef SDOUBLE
   #define RITVEC ritvecddd
  #else
   #define RITVEC ritvecfdd
  #endif
 #else
  #ifdef SDOUBLE
   #define RITVEC ritvecdfd
  #else
   #define RITVEC ritvecffd
  #endif
 #endif
#else
 #ifdef STOREDOUBLE
  #ifdef SDOUBLE
   #define RITVEC ritvecddf
  #else
   #define RITVEC ritvecfdf
  #endif
 #else
  #ifdef SDOUBLE
   #define RITVEC ritvecdff
  #else
   #define RITVEC ritvecfff
  #endif
 #endif
#endif
*/


long   ritvecfff(long n, SMat A, SVDRec R, double kappa, double *ritz, 
		 double *bnd, double *alf, double *bet, float *w2, LanStoreStruct *lss,
              long steps, long neig,long size,ArrUnion *genInv,long neigForGenInv,float minEig,float fracEig,double minGIVal, int transposeGI,int prepareForGenInv,int considerSwap,int nthreads,long *ierr);
long   ritvecffd(long n, SMat A, SVDRec R, double kappa, double *ritz, 
              double *bnd, double *alf, double *bet, float *w2, LanStoreStruct *lss,
              long steps, long neig,long size,ArrUnion *genInv,long neigForGenInv,float minEig,float fracEig,double minGIVal, int transposeGI,int prepareForGenInv,int considerSwap,int nthreads,long *ierr);
long   ritvecdff(long n, SMat A, SVDRec R, double kappa, double *ritz, 
              double *bnd, double *alf, double *bet, float *w2, LanStoreStruct *lss,
              long steps, long neig,long size,ArrUnion *genInv,long neigForGenInv,float minEig,float fracEig,double minGIVal, int transposeGI,int prepareForGenInv,int considerSwap,int nthreads,long *ierr);
long   ritvecdfd(long n, SMat A, SVDRec R, double kappa, double *ritz, 
              double *bnd, double *alf, double *bet, float *w2, LanStoreStruct *lss,
              long steps, long neig,long size,ArrUnion *genInv,long neigForGenInv,float minEig,float fracEig,double minGIVal, int transposeGI,int prepareForGenInv,int considerSwap,int nthreads,long *ierr);
long   ritvecfdf(long n, SMat A, SVDRec R, double kappa, double *ritz, 
              double *bnd, double *alf, double *bet, double *w2, LanStoreStruct *lss,
              long steps, long neig,long size,ArrUnion *genInv,long neigForGenInv,float minEig,float fracEig,double minGIVal, int transposeGI,int prepareForGenInv,int considerSwap,int nthreads,long *ierr);
long   ritvecfdd(long n, SMat A, SVDRec R, double kappa, double *ritz, 
              double *bnd, double *alf, double *bet, double *w2, LanStoreStruct *lss,
              long steps, long neig,long size,ArrUnion *genInv,long neigForGenInv,float minEig,float fracEig,double minGIVal, int transposeGI,int prepareForGenInv,int considerSwap,int nthreads,long *ierr);
long   ritvecddf(long n, SMat A, SVDRec R, double kappa, double *ritz, 
              double *bnd, double *alf, double *bet, double *w2, LanStoreStruct *lss,
              long steps, long neig,long size,ArrUnion *genInv,long neigForGenInv,float minEig,float fracEig,double minGIVal, int transposeGI,int prepareForGenInv,int considerSwap,int nthreads,long *ierr);
long   ritvecddd(long n, SMat A, SVDRec R, double kappa, double *ritz, 
              double *bnd, double *alf, double *bet, double *w2, LanStoreStruct *lss,
              long steps, long neig,long size,ArrUnion *genInv,long neigForGenInv,float minEig,float fracEig,double minGIVal, int transposeGI,int prepareForGenInv,int considerSwap,int nthreads,long *ierr);


long computeGenInv_d_t(DMat Ut,DMat Vt,float *S, long neig, long neigForGenInv, float fracEig, float minEig, ArrUnion *genInv,double *minGIVal,long *cnt, int ncnt, int transposeGI,int nthreads,long evalStart);
long computeGenInv_f_t(DMat Ut,DMat Vt,float *S, long neig, long neigForGenInv, float fracEig, float minEig, ArrUnion *genInv,double *minGIVal,long *cnt, int ncnt, int transposeGI,int nthreads,long evalStart);
//void runGenInv_d(compGenInvStruct *info);
//void runGenInv_f(compGenInvStruct *info);
#endif //RITVEC_H
