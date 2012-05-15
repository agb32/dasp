#include <Python.h>
#include "numpy/arrayobject.h"//lib/python2.5/site-packages/numpy/core/include/numpy/arrayobject.h
//#include <Numeric/arrayobject.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "sparsemem.h"
#include "svdlib.h"
#include "ritvec.h"
#include "genInv.h"
#include <pthread.h>

//#define ISCONTIGUOUS(m) ((m)->flags & CONTIGUOUS)
static PyObject *SvdError;
#define SAFE_FREE(a) {if (a) {free(a); a = NULL;}}

int testCF(PyArrayObject *a){
  //test sparse float
  if(a->descr->type_num!=NPY_FLOAT || !PyArray_ISCONTIGUOUS(a)){
    return 0;
  }
  return 1;
}
int testCF1(PyArrayObject *a){
  //test sparse float 1D
  if(a->nd!=1 || a->descr->type_num!=NPY_FLOAT || !PyArray_ISCONTIGUOUS(a)){
    return 0;
  }
  return 1;
}
int testCF2(PyArrayObject *a){
  //test sparse float 1D
  if(a->nd!=2 || a->descr->type_num!=NPY_FLOAT || !PyArray_ISCONTIGUOUS(a)){
    return 0;
  }
  return 1;
}


int testCI1(PyArrayObject *a){
  //test sparse float 1D
  if(a->nd!=1 || a->descr->kind!='i' || a->descr->elsize!=sizeof(int) || !PyArray_ISCONTIGUOUS(a)){
    return 0;
  }
  return 1;
}
int testCL1(PyArrayObject *a){
  //test sparse float 1D
  if(a->nd!=1 || a->descr->kind!='i' || a->descr->elsize!=sizeof(long) || !PyArray_ISCONTIGUOUS(a)){
    return 0;
  }
  return 1;
}
int testCUI1(PyArrayObject *a){
  //test sparse float 1D
  if(a->nd!=1 || a->descr->type_num!=NPY_UINT || !PyArray_ISCONTIGUOUS(a)){
    return 0;
  }
  return 1;
}
int testCUL1(PyArrayObject *a){
  //test sparse float 1D
  if(a->nd!=1 || a->descr->type_num!=NPY_ULONG || !PyArray_ISCONTIGUOUS(a)){
    return 0;
  }
  return 1;
}

/*
static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
        (dmaxarg1) : (dmaxarg2))
static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
       (iminarg1) : (iminarg2))

#define NR_END 1
#define FREE_ARG char*
*/
/*double *dvector(long nl, long nh)
// allocate a double vector with subscript range v[nl..nh] 
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v){
	  printf("allocation failure in dvector()");
	  exit(0);
	}
	return v-nl+NR_END;
}

void free_dvector(double *v, long nl, long nh)
// free a double vector allocated with dvector() 
{
	free((FREE_ARG) (v+nl-NR_END));
}
*/
/*
double dpythag(double a, double b){
  //computes (a^2+b^2)^0.5 without destructive over or underflow.
  double absa,absb;
  absa=fabs(a);
  absb=fabs(b);
  if(absa>absb)
    return absa*sqrt(1.+DSQR(absb/absa));
  else
    return (absb==0.?0.:absb*sqrt(1.+DSQR(absa/absb)));
}

void dsvdcmp(double **a, int m, int n, double w[], double **v)
{
  //given a matrix a[1..m][1..n], this routine computes its SVD,
  //A=UWV^T.  The matrix U replaces a on output. The diagonal matrix
  //of singular values W is output as a vector w[1..n].  The matrix V
  //(not V^T) is output as v[1..n][1..n].
  //Here, n=ncents, m=nacts.
  int flag,i,its,j,jj,k,l,nm;
  double anorm,c,f,g,h,s,scale,x,y,z,*rv1;
  
  rv1=dvector(1,n);
  g=scale=anorm=0.0;
  //householder reduction to bidiagonal form.
  for (i=1;i<=n;i++) {
    l=i+1;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i <= m) {
      for (k=i;k<=m;k++) scale += fabs(a[k][i]);
      if (scale) {
	for (k=i;k<=m;k++) {
	  a[k][i] /= scale;
	  s += a[k][i]*a[k][i];
	}
	f=a[i][i];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	a[i][i]=f-g;
	for (j=l;j<=n;j++) {
	  for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
	  f=s/h;
	  for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
	}
	for (k=i;k<=m;k++) a[k][i] *= scale;
      }
    }
    w[i]=scale *g;
    g=s=scale=0.0;
    if (i <= m && i != n) {
      for (k=l;k<=n;k++) scale += fabs(a[i][k]);
      if (scale) {
	for (k=l;k<=n;k++) {
	  a[i][k] /= scale;
	  s += a[i][k]*a[i][k];
	}
	f=a[i][l];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	a[i][l]=f-g;
	for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
	for (j=l;j<=m;j++) {
	  for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
	  for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
	}
	for (k=l;k<=n;k++) a[i][k] *= scale;
      }
    }
    anorm=DMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
  }
  //accumulation of right-hand transformations.
  //here, l=n+1 I think
  for (i=n;i>=1;i--) {
    if (i < n) {
      if (g) {
	for (j=l;j<=n;j++)
	  v[j][i]=(a[i][j]/a[i][l])/g;//double division to avoid possible underflow.
	for (j=l;j<=n;j++) {
	  for (s=0.0,k=l;k<=n;k++)
	    s += a[i][k]*v[k][j];
	  for (k=l;k<=n;k++)
	    v[k][j] += s*v[k][i];
	}
      }
      for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
    }
    v[i][i]=1.0;
    g=rv1[i];
    l=i;
  }
  for (i=IMIN(m,n);i>=1;i--) {//accumulation of left-hand transformations.
    l=i+1;
    g=w[i];
    for (j=l;j<=n;j++) a[i][j]=0.0;
    if (g) {
      g=1.0/g;
      for (j=l;j<=n;j++) {
	for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
	f=(s/a[i][i])*g;
	for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
      }
      for (j=i;j<=m;j++) a[j][i] *= g;
    } else for (j=i;j<=m;j++) a[j][i]=0.0;
    ++a[i][i];
  }
  for (k=n;k>=1;k--) {//diagonalisation of the bidiagonal form:  Loop over singular values, and over allowed iterations.
    for (its=1;its<=30;its++) {
      flag=1;
      for (l=k;l>=1;l--) {//test for splitting. Note that rv1[1] is always zero.
	nm=l-1;
	if ((double)(fabs(rv1[l])+anorm) == anorm) {
	  flag=0;
	  break;
	}
	if ((double)(fabs(w[nm])+anorm) == anorm) break;
      }
      if (flag) {
	c=0.0;//cancellation of rv1[l] if l>1.
	s=1.0;
	for (i=l;i<=k;i++) {
	  f=s*rv1[i];
	  rv1[i]=c*rv1[i];
	  if ((double)(fabs(f)+anorm) == anorm) break;
	  g=w[i];
	  h=dpythag(f,g);
	  w[i]=h;
	  h=1.0/h;
	  c=g*h;
	  s = -f*h;
	  for (j=1;j<=m;j++) {
	    y=a[j][nm];
	    z=a[j][i];
	    a[j][nm]=y*c+z*s;
	    a[j][i]=z*c-y*s;
	  }
	}
      }
      z=w[k];
      if (l == k) {//convergence.
	if (z < 0.0) {//singular value is made nonnegative.
	  w[k] = -z;
	  for (j=1;j<=n;j++) v[j][k] = -v[j][k];
	}
	break;
      }
      if (its == 30){
	printf("no convergence in 30 dsvdcmp iterations");
	exit(0);
      }
      x=w[l];//shift from bottom 2 by 2 minor.
      nm=k-1;
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=dpythag(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
      c=s=1.0;//next QR transformation.
      for (j=l;j<=nm;j++) {
	i=j+1;
	g=rv1[i];
	y=w[i];
	h=s*g;
	g=c*g;
	z=dpythag(f,h);
	rv1[j]=z;
	c=f/z;
	s=h/z;
	f=x*c+g*s;
	g = g*c-x*s;
	h=y*s;
	y *= c;
	for (jj=1;jj<=n;jj++) {
	  x=v[jj][j];
	  z=v[jj][i];
	  v[jj][j]=x*c+z*s;
	  v[jj][i]=z*c-x*s;
	}
	z=dpythag(f,h);//rotation can be arbitrary if z=0.
	w[j]=z;
	if (z) {
	  z=1.0/z;
	  c=f*z;
	  s=h*z;
	}
	f=c*g+s*y;
	x=c*y-s*g;
	for (jj=1;jj<=m;jj++) {
	  y=a[jj][j];
	  z=a[jj][i];
	  a[jj][j]=y*c+z*s;
	  a[jj][i]=z*c-y*s;
	}
      }
      rv1[l]=0.0;
      rv1[k]=f;
      w[k]=x;
    }
  }
  free_dvector(rv1,1,n);
}

void dsvbksb(double **u, double w[], double **v, int m, int n, double b[], double x[]){
  //solves Ax=b where A is specified by u[1..m][1..n], w[1..n], v[1..n][1..n] as returned by
  //dsvdcmp.  m,n are the dimensions of a[1..m][1..n].  b[1..m] is the input right hand
  //side.  x[1..n] is the output solution vector.  No input quantities are
  //destroyed so the routine may be called many times.
  int jj,j,i;
  double s,*tmp;
  
  tmp=dvector(1,n);
  for (j=1;j<=n;j++) {
    s=0.0;
    if (w[j]) {
      for (i=1;i<=m;i++)
	s += u[i][j]*b[i];
      s /= w[j];
    }
    tmp[j]=s;
  }
  for (j=1;j<=n;j++) {
    s=0.0;
    for (jj=1;jj<=n;jj++) 
      s += v[j][jj]*tmp[jj];
    x[j]=s;
  }
  free_dvector(tmp,1,n);
}

*/
/*
int agbdsvdcmp(double *a, int m, int n, double *w, double *v){
  //given an array a (size m*n) -> [0..m-1][0..n-1], this routine computes its SVD,
  //A=UWV^T.  The matrix U replaces a on output. The diagonal matrix
  //of singular values W is output as an array w[0..n-1].  The matrix V
  //(not V^T) is output as an array v[0..n-1][0..n-1].
  //Here, n=ncents, m=nacts.
  int flag,i,its,j,jj,k,l,nm;
  double anorm,c,f,g,h,s,scale,x,y,z,*rv1;
  
  rv1=(double*)malloc(sizeof(double)*n);//dvector(1,n);
  if(rv1==NULL){
    return -1;
  }
  g=scale=anorm=0.0;
  //householder reduction to bidiagonal form.
  for (i=0;i<n;i++) {
    l=i+1;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i < m) {
      for (k=i;k<m;k++) scale += fabs(a[k*n+i]);
      if (scale) {
	for (k=i;k<m;k++) {
	  a[k*n+i] /= scale;
	  s += a[k*n+i]*a[k*n+i];
	}
	f=a[i*n+i];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	a[i*n+i]=f-g;
	for (j=l;j<n;j++) {
	  s=0.;
	  for (k=i;k<m;k++)
	    s += a[k*n+i]*a[k*n+j];
	  f=s/h;
	  for (k=i;k<m;k++)
	    a[k*n+j] += f*a[k*n+i];
	}
	for (k=i;k<m;k++) a[k*n+i] *= scale;
      }
    }
    w[i]=scale *g;
    g=s=scale=0.0;
    if (i < m && i != n-1) {
      for (k=l;k<n;k++)
	scale += fabs(a[i*n+k]);
      if (scale) {
	for (k=l;k<n;k++) {
	  a[i*n+k] /= scale;
	  s += a[i*n+k]*a[i*n+k];
	}
	f=a[i*n+l];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	a[i*n+l]=f-g;
	for (k=l;k<n;k++)
	  rv1[k]=a[i*n+k]/h;
	for (j=l;j<m;j++) {
	  s=0.;
	  for (k=l;k<n;k++)
	    s += a[j*n+k]*a[i*n+k];
	  for (k=l;k<n;k++)
	    a[j*n+k] += s*rv1[k];
	}
	for (k=l;k<n;k++)
	  a[i*n+k] *= scale;
      }
    }
    anorm=DMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
  }
  //accumulation of right-hand transformations.
  //here, l=n+1 I think
  for (i=n-1;i>=0;i--) {
    if (i < n-1) {
      if (g) {
	for (j=l;j<n;j++)
	  v[j*n+i]=(a[i*n+j]/a[i*n+l])/g;//double division to avoid possible underflow.
	for (j=l;j<n;j++) {
	  s=0.;
	  for (k=l;k<n;k++)
	    s += a[i*n+k]*v[k*n+j];
	  for (k=l;k<n;k++)
	    v[k*n+j] += s*v[k*n+i];
	}
      }
      for (j=l;j<n;j++)
	v[i*n+j]=v[j*n+i]=0.0;
    }
    v[i*n+i]=1.0;
    g=rv1[i];
    l=i;
  }
  for (i=IMIN(m,n)-1;i>=0;i--) {//accumulation of left-hand transformations.
    l=i+1;
    g=w[i];
    for (j=l;j<n;j++)
      a[i*n+j]=0.0;
    if (g) {
      g=1.0/g;
      for (j=l;j<n;j++) {
	s=0.;
	for (k=l;k<m;k++)
	  s += a[k*n+i]*a[k*n+j];
	f=(s/a[i*n+i])*g;
	for (k=i;k<m;k++)
	  a[k*n+j] += f*a[k*n+i];
      }
      for (j=i;j<m;j++)
	a[j*n+i] *= g;
    } else{
      for (j=i;j<m;j++)
	a[j*n+i]=0.0;
    }
    ++a[i*n+i];
  }
  for (k=n-1;k>=0;k--) {//diagonalisation of the bidiagonal form:  Loop over singular values, and over allowed iterations.
    for (its=1;its<=30;its++) {
      flag=1;
      for (l=k;l>=0;l--) {//test for splitting. Note that rv1[0] is always zero.
	nm=l-1;
	if ((double)(fabs(rv1[l])+anorm) == anorm) {
	  flag=0;
	  break;
	}
	if ((double)(fabs(w[nm])+anorm) == anorm) break;//should never get here when l=0...
      }
      if (flag) {
	c=0.0;//cancellation of rv1[l] if l>0.
	s=1.0;
	for (i=l;i<k;i++) {
	  f=s*rv1[i];
	  rv1[i]=c*rv1[i];
	  if ((double)(fabs(f)+anorm) == anorm) break;
	  g=w[i];
	  h=dpythag(f,g);
	  w[i]=h;
	  h=1.0/h;
	  c=g*h;
	  s = -f*h;
	  for (j=0;j<m;j++) {
	    y=a[j*n+nm];
	    z=a[j*n+i];
	    a[j*n+nm]=y*c+z*s;
	    a[j*n+i]=z*c-y*s;
	  }
	}
      }
      z=w[k];
      if (l == k) {//convergence.
	if (z < 0.0) {//singular value is made nonnegative.
	  w[k] = -z;
	  for (j=0;j<n;j++)
	    v[j*n+k] = -v[j*n+k];
	}
	break;
      }
      if (its == 30){
	printf("no convergence in 30 agbdsvdcmp iterations");
	free(rv1);//free_dvector(rv1,1,n);
	return -1;
      }
      x=w[l];//shift from bottom 2 by 2 minor.
      nm=k-1;//not sure what happens here when nm==-1...?!?  Actually, I think it works out okay!
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=dpythag(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
      c=s=1.0;//next QR transformation.
      for (j=l;j<nm;j++) {
	i=j+1;
	g=rv1[i];
	y=w[i];
	h=s*g;
	g=c*g;
	z=dpythag(f,h);
	rv1[j]=z;
	c=f/z;
	s=h/z;
	f=x*c+g*s;
	g = g*c-x*s;
	h=y*s;
	y *= c;
	for (jj=0;jj<n;jj++) {
	  x=v[jj*n+j];
	  z=v[jj*n+i];
	  v[jj*n+j]=x*c+z*s;
	  v[jj*n+i]=z*c-x*s;
	}
	z=dpythag(f,h);//rotation can be arbitrary if z=0.
	w[j]=z;
	if (z) {
	  z=1.0/z;
	  c=f*z;
	  s=h*z;
	}
	f=c*g+s*y;
	x=c*y-s*g;
	for (jj=0;jj<m;jj++) {
	  y=a[jj*n+j];
	  z=a[jj*n+i];
	  a[jj*n+j]=y*c+z*s;
	  a[jj*n+i]=z*c-y*s;
	}
      }
      rv1[l]=0.0;
      rv1[k]=f;
      w[k]=x;
    }
  }
  free(rv1);//free_dvector(rv1,1,n);
  return 0;
}
static PyObject* agbsvd(PyObject *self,PyObject *args){
  PyArrayObject *A,*W,*V;
  int m,n,err;
  if(!PyArg_ParseTuple(args,"O!O!O!",&PyArray_Type,&A,&PyArray_Type,&W,&PyArray_Type,&V)){
    printf("Usage: Input and U output, W, V arrays.\n");
    return NULL;
  }
  if(A->descr->type_num!=NPY_DOUBLE || W->descr->type_num!=NPY_DOUBLE || V->descr->type_num!=NPY_DOUBLE){
    printf("Arrays must be double\n");
    return NULL;
  }
  if(PyArray_ISCONTIGUOUS(A)==0 || PyArray_ISCONTIGUOUS(W)==0 || PyArray_ISCONTIGUOUS(V)==0){
    printf("Arrays must be contiguous\n");
    return NULL;
  }
  m=A->dimensions[0];
  n=A->dimensions[1];
  if(W->dimensions[0]!=n || V->dimensions[0]!=n || V->dimensions[0]!=n){
    printf("Arrays are the wrong shape.  A->m,n, W->n, V->n,n\n");
    return NULL;
  }
  if((err=agbdsvdcmp((double*)A->data, m, n, (double*)W->data, (double*)V->data))==-1){
    printf("Error during agbdsvdcmp %d\n",err);
    return NULL;
  }
  return Py_BuildValue("i",0);
}
*/
static PyObject *dolibsvd(PyObject *self,PyObject *args){
  //perform svd using libsvd... (which is specifically for sparse matricees).
  PyArrayObject *Adata,*Arowind,*Aindptr,*U,*W,*V;
  PyObject *Vobj;
  long m,n,size=0;
  //int i;
  SMat spA;
  SVDRec res,tmpres;
  dSpMem *genInvPtrD=NULL;
  double end[2] = {-1.0e-30, 1.0e-30};
  double kappa = 1e-6;
  int neig=0;
  PyObject *genInv=NULL;
  PyArrayObject *genInvArr;
  float minEig=0.,fracEig=0.;
  int useStoreFloat=0,useSFloat=0;
  double minGIVal=0.;
  int transposeGI=0;
  int nthreads=2;
  int prepareForGenInv=1;
  int considerSwap=0;
  ArrUnion *genInvUnion=NULL;
  if(!PyArg_ParseTuple(args,"llO!O!O!O!O!O|liOffiidiii",&m,&n,&PyArray_Type,&Adata,&PyArray_Type,&Arowind,&PyArray_Type,&Aindptr,&PyArray_Type,&U,&PyArray_Type,&W,&Vobj,&size,&neig,&genInv,&minEig,&fracEig,&useStoreFloat,&useSFloat,&minGIVal,&prepareForGenInv,&nthreads,&considerSwap)){
    printf("Usage: m(nacts),n(ncents),Input data, input rowind, input indptr (all from the sparse matrix csc format), U, W, V arrays (dense, V can be None), size (optional, max size for temporary array in algorithm. If 0, any size allowed.  If not, a sparse implementation may be used.  Size is the number of doubles to use, not memory.), neig (optional, the number of eigen values to use, defaults to zero, meaning all), genInv (optional, None or array (dense or sparse) to hold the generalised inverse), minEig (optional - value of eigenvalues to ignore), fracEig (optional - fractional value of eigenvalues to ignore), useStoreFloat (optional, default 0 - whether to use float32 for internal Store array), useSFloat (optional, default 0 - whether to use float32 for internal s array), minGIVal (optional - minimum value to put in generalised inverse),prepareForGenInv (optional, default 1), nthreads (optional, default 2), considerSwap (optional, can be used if Vt will be mmaped to speed things up, default=0).");
    return NULL;
  }
  if(Arowind->descr->kind!='i' || Arowind->descr->elsize!=sizeof(int) || Aindptr->descr->kind!='i' || Aindptr->descr->elsize!=sizeof(int)){
    printf("Arowind and Aindptr arrays must be type int\n");
    return NULL;
  }
  if(W->descr->type_num!=NPY_FLOAT){
    printf("W Array must be float\n");
    return NULL;
  }
  if(Adata->descr->type_num!=NPY_DOUBLE && Adata->descr->type_num!=NPY_FLOAT){
    printf("Adata array must be float or double\n");
    return NULL;
  }
  if(PyArray_Check(Vobj)){//its a Numeric array...
    V=(PyArrayObject*)Vobj;
  }else if(Vobj!=Py_None){//its not None
    printf("V must be None or a Numeric array\n");
    return NULL;
  }else{//its None...
    V=NULL;
  }

  if((U->descr->type_num!=NPY_DOUBLE && U->descr->type_num!=NPY_FLOAT)){
    printf("U must be double or float\n");
    return NULL;
  }
  if(V!=NULL && U->descr->type_num!=V->descr->type_num){
    printf("V must be of same type as U (or None)\n");
    return NULL;
  }
  
  if(PyArray_ISCONTIGUOUS(Adata)==0 || PyArray_ISCONTIGUOUS(Arowind)==0 || PyArray_ISCONTIGUOUS(Aindptr)==0 || PyArray_ISCONTIGUOUS(U)==0 || PyArray_ISCONTIGUOUS(W)==0 || (V!=NULL && PyArray_ISCONTIGUOUS(V)==0)){
    printf("Arrays must be contiguous\n");
    return NULL;
  }
  if(Aindptr->dimensions[0]!=n+1){
    printf("Aindptr must have dimensions n+1\n");
    return NULL;
  }
  if(Arowind->dimensions[0]!=Adata->dimensions[0]){
    printf("rowind and data must have same dimensions\n");
    return NULL;
  }
  if(W->nd!=1 || U->nd!=2 || W->dimensions[0]!=n || U->dimensions[0]!=m || U->dimensions[1]!=m || (V!=NULL && (V->nd!=2 || V->dimensions[0]!=m || V->dimensions[1]!=n))){
    printf("Output arrays are the wrong shape.  U->nacts,nacts, W->ncents, V->nacts,ncents\n");
    return NULL;
  }
  if(genInv!=NULL){
    genInvUnion=(ArrUnion*)malloc(sizeof(ArrUnion));
    genInvUnion->typ='n';
    if(PyArray_Check(genInv)){//its a Numeric array
      genInvArr=(PyArrayObject*)genInv;
      if((genInvArr->descr->type_num!=NPY_FLOAT && genInvArr->descr->type_num!=NPY_DOUBLE) || PyArray_ISCONTIGUOUS(genInvArr)==0 || genInvArr->nd!=2){
	printf("genInv must be a Float32/64 python array, contiguous, 2D with shape (%ld,%ld) or transposed\n",n,m);
	return NULL;
      }
      if(genInvArr->dimensions[0]==m && genInvArr->dimensions[1]==n){
	transposeGI=0;
      }else if(genInvArr->dimensions[0]==n && genInvArr->dimensions[1]==m){
	transposeGI=1;
      }else{
	printf("genInv shape must be %ld,%ld or %ld,%ld\n",n,m,m,n);
	return NULL;
      }
      if(genInvArr->descr->type_num==NPY_FLOAT){
	genInvUnion->typ='f';
	genInvUnion->data.fdata=(float*)genInvArr->data;
      }else{
	genInvUnion->typ='d';
	genInvUnion->data.ddata=(double*)genInvArr->data;
      }
    }else if(PyLong_Check(genInv) || PyInt_Check(genInv)){
      printf("Using sparse array for generalised inverse\n");
      if(PyLong_Check(genInv))
	genInvPtrD=(dSpMem*)PyLong_AsLong(genInv);
      else
	genInvPtrD=(dSpMem*)PyInt_AsLong(genInv);
      printf("Got address %p\n",genInvPtrD);
      if(genInvPtrD->typ=='f'){//its actually a float type...
	genInvUnion->typ='F';
	genInvUnion->data.fSp=(fSpMem*)genInvPtrD;
	if(genInvUnion->data.fSp->rows==m && genInvUnion->data.fSp->cols==n){
	  transposeGI=0;
	}else if(genInvUnion->data.fSp->rows==n && genInvUnion->data.fSp->cols==m){
	  transposeGI=1;
	}else{
	  printf("genInv must have shape %ld,%ld or %ld,%ld\n",n,m,m,n);
	  return NULL;
	}
      }else if(genInvPtrD->typ=='d'){//its a double type...
	genInvUnion->typ='D';
	genInvUnion->data.dSp=genInvPtrD;
	if(genInvUnion->data.dSp->rows==m && genInvUnion->data.dSp->cols==n){
	  transposeGI=0;
	}else if(genInvUnion->data.dSp->rows==n && genInvUnion->data.dSp->cols==m){
	  transposeGI=1;
	}else{
	  printf("genInv must have shape %ld,%ld or %ld,%ld\n",n,m,m,n);
	  return NULL;
	}

      }else{
	printf("svdmodule: Type of sparse array not known for generalised inverse.\n");
	return NULL;
      }
    }else if(genInv!=Py_None){//its not None
      printf("genInv must be None, a svd.sparseMatrix object or a Numeric array.\n");
      return NULL;
    }
    printf("Using generalised inverse as %c (from f,d,F,D)\n",genInvUnion->typ);
  }

  spA=(SMat)calloc(1,sizeof(struct smat));
  spA->rows=m;
  spA->cols=n;
  spA->vals=Adata->dimensions[0];
  spA->pointr=(unsigned int*)Aindptr->data;
  spA->rowind=(unsigned int*)Arowind->data;
  if(Adata->descr->type_num==NPY_DOUBLE){
    spA->value=(double*)Adata->data;
    spA->valuef=NULL;
  }else{
    spA->value=NULL;
    spA->valuef=(float*)Adata->data;
  }

  res=svdNewSVDRec();//(SVDRec)calloc(1,sizeof(struct svdrec));
  if(!res){
    printf("Allocation of svdRecord (for results) failed\n");
    free(spA);
    return NULL;
  }
  res->d=m;
  res->S=(float*)W->data;
  if (spA->cols >= spA->rows * 1.2){//A will be transposed for speed...
    if(U->descr->type_num==NPY_DOUBLE){
      if(V==NULL)
	res->Ut=svdNewDMatWithNoData(res->d,spA->cols);//allocated later...
      else
	res->Ut=svdNewDMatWithData(res->d,spA->cols,(double*)V->data);
      //this could be mmap'd if v large...
      res->Vt=svdNewDMatWithData(res->d,spA->rows,(double*)U->data);
    }else{//float
      if(V==NULL)
	res->Ut=svdNewDMatWithNoData(res->d,spA->cols);
      else
	res->Ut=svdNewDMatWithDataFloat(res->d,spA->cols,(float*)V->data);
      //this could be mmap'd if v large...
      res->Vt=svdNewDMatWithDataFloat(res->d,spA->rows,(float*)U->data);
    }
  }else{
    printf("Warning - svd has not been tested without the transpose\n");
    if(V==NULL){
      printf("ERROR: Not tested with V==None\n");
      return NULL;
    }
    if(U->descr->type_num==NPY_DOUBLE){
      res->Ut=svdNewDMatWithData(res->d,spA->rows,(double*)U->data);
      res->Vt=svdNewDMatWithData(res->d,spA->cols,(double*)V->data);
    }else{//float
      res->Ut=svdNewDMatWithDataFloat(res->d,spA->rows,(float*)U->data);
      res->Vt=svdNewDMatWithDataFloat(res->d,spA->cols,(float*)V->data);
    }
    
  }
  if (!res->Ut || !res->S || !res->Vt) {
    printf("Allocation of U W or V failed\n");
    free(spA);
    if(res->Ut){
      SAFE_FREE(res->Ut->value);
      SAFE_FREE(res->Ut->valuef);
      SAFE_FREE(res->Ut);
    }
    if(res->Vt){
      SAFE_FREE(res->Vt->value);
      SAFE_FREE(res->Vt->valuef);
      SAFE_FREE(res->Vt);
    }
    SAFE_FREE(genInvUnion);
    free(res);
    return NULL;
  }


  printf("Starting the svd\n");
  //compute the svd
  tmpres=svdLAS2(spA,0,0,end,kappa,res,size,neig,genInvUnion,minEig,fracEig,useStoreFloat,useSFloat,minGIVal,transposeGI,prepareForGenInv,considerSwap,nthreads);

  SAFE_FREE(genInvUnion);
  if(tmpres==NULL){
    printf("svdLAS2A failed\n");
    SAFE_FREE(res->Ut->value);
    SAFE_FREE(res->Ut->valuef);
    SAFE_FREE(res->Ut);
    SAFE_FREE(res->Vt->value);
    SAFE_FREE(res->Vt->valuef);
    SAFE_FREE(res->Vt);
    SAFE_FREE(res);
    free(spA);
    return NULL;
  }
  //now put results into the python arrays...
  if(res->err)
    printf("svdLAS2A Got error %d\n",res->err);
  printf("Results are of size: U %ld %ld, W %d, V %ld %ld\n",res->Ut->rows,res->Ut->cols,res->d,res->Vt->rows,res->Vt->cols);
  /*if(res->d<=n){
    printf("Copying evals\n");
    for(i=0; i<res->d; i++){
      ((double*)W->data)[i]=res->S[i];
    }
    }*/
  m=res->d;
  //now free stuff...
  SAFE_FREE(res->Ut->value);
  SAFE_FREE(res->Ut->valuef);
  SAFE_FREE(res->Ut);
  SAFE_FREE(res->Vt->value);
  SAFE_FREE(res->Vt->valuef);
  SAFE_FREE(res->Vt);
  SAFE_FREE(res);
  free(spA);
  return Py_BuildValue("i",m);//return the number of valid points of W.


}
static PyObject *computeDenseGenInv(PyObject *self,PyObject *args){
  //Compute the generalised inverse using Ut and Vt from previous.  Note, these must have been
  //altered for genInv, ie Ut scaled and Vt transposed.
  PyArrayObject *Ut,*Vt,*W,*genInv,*minGIValArr,*cntArr;
  long neig,neigForGenInv;
  float fracEig,minEig;
  double *minGIVal;
  long *cnt;
  int transposeGI,nthreads,ncnt;
  long evalStart=0;
  unsigned int nacts,ncents;
  //long n;
  DMat UtMat=NULL, VtMat=NULL;
  ArrUnion *genInvUnion;
  if(!PyArg_ParseTuple(args,"O!O!O!llffO!O!O!iil",&PyArray_Type,&Ut,&PyArray_Type,&Vt,&PyArray_Type,&W,&neig,&neigForGenInv,&fracEig,&minEig,&PyArray_Type,&genInv,&PyArray_Type,&minGIValArr,&PyArray_Type,&cntArr,&transposeGI,&nthreads,&evalStart)){
    printf("Usage: Ut, Vt, W, neig, neigForGenInv, fracEig, minEig, genInv, minGIVal, transposeGI, nthreads, evalStart \n");
    printf("Note, Ut/Vt should have been prepared (scaled and transposed) by the call to SVD decomposition.\n");
    return NULL;
  }
  if(W->descr->type_num!=NPY_FLOAT){
    printf("Error - W must be float32\n");
    return NULL;
  }
  if(W->nd!=1 || W->dimensions[0]<neig){
    printf("W is wrong shape\n");
    return NULL;
  }
  if(minGIValArr->descr->type_num!=NPY_DOUBLE || cntArr->descr->type_num!=NPY_LONG || minGIValArr->nd!=1 || cntArr->nd!=1 || minGIValArr->dimensions[0]!=cntArr->dimensions[0]){
    printf("minGIVal and cnt must be double and long respectively, 1D, and of same length.\n");
    return NULL;
  }
  ncnt=cntArr->dimensions[0];
  minGIVal=(double*)minGIValArr->data;
  cnt=(long*)cntArr->data;

  genInvUnion=malloc(sizeof(ArrUnion));
  if(genInv->descr->type_num==NPY_DOUBLE){
    genInvUnion->typ='d';
    genInvUnion->data.ddata=(double*)genInv->data;
  }else if(genInv->descr->type_num==NPY_FLOAT){
    genInvUnion->typ='f';
    genInvUnion->data.fdata=(float*)genInv->data;
  }else{
    printf("genInv must be float32 or 64\n");
    free(genInvUnion);
    return NULL;
  }
  //here, we assume that A was transposed - so we swap Vt and Ut around.
  if(Ut->descr->type_num==NPY_FLOAT && Vt->descr->type_num==NPY_FLOAT){
    VtMat=svdNewDMatWithDataFloat(Ut->dimensions[0],Ut->dimensions[1],(float*)Ut->data);
    UtMat=svdNewDMatWithDataFloat(Vt->dimensions[0],Vt->dimensions[1],(float*)Vt->data);
  }else if(Ut->descr->type_num==NPY_DOUBLE && Vt->descr->type_num==NPY_DOUBLE){
    VtMat=svdNewDMatWithData(Ut->dimensions[0],Ut->dimensions[1],(double*)Ut->data);
    UtMat=svdNewDMatWithData(Vt->dimensions[0],Vt->dimensions[1],(double*)Vt->data);
  }else{
    printf("Ut and Vt must be same type, and both either float32 or 64\n");
    return NULL;
  }
  //do some checks on shape
  nacts=VtMat->cols;
  ncents=UtMat->cols;
  if(transposeGI){
    if(genInv->nd!=2 || genInv->dimensions[0]!=ncents || genInv->dimensions[1]!=nacts){
      printf("genInv wrong shape\n");
      return NULL;
    }
  }else{
    if(genInv->nd!=2 || genInv->dimensions[1]!=ncents || genInv->dimensions[0]!=nacts){
      printf("genInv wrong shape\n");
      return NULL;
    }
  }
  if(Vt->descr->type_num==NPY_FLOAT)
    computeGenInv_f_t(UtMat,VtMat,(float*)W->data,neig,neigForGenInv,fracEig,minEig,genInvUnion,minGIVal,cnt,ncnt,transposeGI,nthreads,evalStart);
  else
    computeGenInv_d_t(UtMat,VtMat,(float*)W->data,neig,neigForGenInv,fracEig,minEig,genInvUnion,minGIVal,cnt,ncnt,transposeGI,nthreads,evalStart);
  printf("computeGenInv_f/d_t done\n");

  free(genInvUnion);
  SAFE_FREE(VtMat->value);
  SAFE_FREE(VtMat->valuef);
  SAFE_FREE(UtMat->value);
  SAFE_FREE(UtMat->valuef);
  SAFE_FREE(VtMat);
  SAFE_FREE(UtMat);
  return Py_BuildValue("l",0);
}

static PyObject *sparsifyGenInv(PyObject *self,PyObject *args){
  ArrUnion *genInvUnion,*spGenInvUnion;
  long ncols,nrows;
  PyArrayObject *genInv,*data,*rowind,*indptr;
  int transposeGI;
  dSpMem *dSp;
  fSpMem *fSp;
  long ndata;
  double minGIVal;
  int rtval=0;
  if(!PyArg_ParseTuple(args,"O!dilO!O!O!",&PyArray_Type,&genInv,&minGIVal,&transposeGI,&nrows,&PyArray_Type,&data,&PyArray_Type,&rowind,&PyArray_Type,&indptr)){
    printf("Usage: genInv, minGIVal, transposeGI, nrows,data,rowind,indptr\n");
    return NULL;
  }
  genInvUnion=malloc(sizeof(ArrUnion));
  spGenInvUnion=malloc(sizeof(ArrUnion));
  if(genInv->nd!=2 || data->nd!=1 || rowind->nd!=1 || indptr->nd!=1){
    printf("arrays are wrong shape\n");
    return NULL;
  }
  ncols=indptr->dimensions[0]-1;
  if((data->descr->type_num!=NPY_FLOAT && data->descr->type_num!=NPY_DOUBLE) || rowind->descr->kind!='i' || rowind->descr->elsize!=sizeof(int) || indptr->descr->kind!='i' || indptr->descr->elsize!=sizeof(int)){
    printf("csc arrays are wrong type\n");
    return NULL;
  }
  ndata=data->dimensions[0];
  if(genInv->descr->type_num==NPY_FLOAT){
    genInvUnion->typ='f';
    genInvUnion->data.fdata=(float*)genInv->data;
  }else if(genInv->descr->type_num==NPY_DOUBLE){
    genInvUnion->typ='d';
    genInvUnion->data.ddata=(double*)genInv->data;
  }else{
    printf("genInv must be double or float\n");
    return NULL;
  }
  if(data->descr->type_num==NPY_FLOAT){
    fSp=malloc(sizeof(fSpMem));
    fSp->data=(float*)data->data;
    fSp->rowind=(unsigned int*)rowind->data;
    fSp->indptr=(unsigned int*)indptr->data;
    fSp->typ='f';
    fSp->ndata=ndata;
    fSp->rows=nrows;
    fSp->cols=ncols;
    fSp->min=0.;
    fSp->rowmin=0;
    fSp->indmin=0;
    fSp->cnt=0;
    fSp->alloced=0;
    spGenInvUnion->typ='F';
    spGenInvUnion->data.fSp=fSp;
  }else{
    dSp=malloc(sizeof(dSpMem));
    dSp->data=(double*)data->data;
    dSp->rowind=(unsigned int*)rowind->data;
    dSp->indptr=(unsigned int*)indptr->data;
    dSp->typ='d';
    dSp->ndata=ndata;
    dSp->rows=nrows;
    dSp->cols=ncols;
    dSp->min=0.;
    dSp->rowmin=0;
    dSp->indmin=0;
    dSp->cnt=0;
    dSp->alloced=0;
    spGenInvUnion->typ='D';
    spGenInvUnion->data.dSp=dSp;
  }
  Py_BEGIN_ALLOW_THREADS;
  rtval=convertToSparse(genInvUnion,minGIVal,transposeGI,spGenInvUnion);
  Py_END_ALLOW_THREADS;
  if(spGenInvUnion->typ=='D'){
    SAFE_FREE(spGenInvUnion->data.dSp);
  }else{
    SAFE_FREE(spGenInvUnion->data.fSp);
  }
  SAFE_FREE(spGenInvUnion);
  SAFE_FREE(genInvUnion);
  return Py_BuildValue("i",rtval);
}

//now have some methods to do with the sparse matrix class that I've implemented...
static PyObject *sparseMatrixCreate(PyObject *self,PyObject *args){
  PyArrayObject *Adata,*Arowind,*Aindptr;
  unsigned int rows,cols,ndata;
  dSpMem *sp=NULL;
  fSpMem *fsp=NULL;
  //int i;
  if(!PyArg_ParseTuple(args,"IIIO!O!O!",&ndata,&rows,&cols,&PyArray_Type,&Adata,&PyArray_Type,&Arowind,&PyArray_Type,&Aindptr)){
    printf("Usage: ndata, rows, cols, data, rowind, indptr\n");
    return NULL;
  }
  if((Adata->descr->type_num!=NPY_DOUBLE && Adata->descr->type_num!=NPY_FLOAT) || PyArray_ISCONTIGUOUS(Adata)==0 || Adata->dimensions[0]!=ndata){
    printf("data must be float/double array, contiguous, of length ndata\n");
    return NULL;
  }
  if(Arowind->descr->kind!='i' || Arowind->descr->elsize!=sizeof(int) || PyArray_ISCONTIGUOUS(Arowind)==0 || Arowind->dimensions[0]!=ndata){
    printf("rowind must be int32 array, contiguous, of length ndata\n");
    return NULL;
  }
  if(Aindptr->descr->kind!='i' || Aindptr->descr->elsize!=sizeof(int) || PyArray_ISCONTIGUOUS(Aindptr)==0 || Aindptr->dimensions[0]!=cols+1){
    printf("indptr must be int32 array, contiguous, of length cols+1\n");
    return NULL;
  }
  if(Adata->descr->type_num==NPY_DOUBLE){
    sp=smNewFromExisting(ndata,rows,cols,(double*)Adata->data,(unsigned int*)Aindptr->data,(unsigned int*)Arowind->data);
    if(sp==NULL){
      printf("Error creating sparse memory object\n");
    }
  }else{
    fsp=smNewFromExistingFloat(ndata,rows,cols,(float*)Adata->data,(unsigned int*)Aindptr->data,(unsigned int*)Arowind->data);
    if(fsp==NULL){
      printf("Error creating sparse memory object\n");
    }
  }
  return Py_BuildValue("l",(long)(sp==NULL?(dSpMem*)fsp:sp));//return a pointer to it!
}
static PyObject *sparseMatrixFree(PyObject *self,PyObject *args){
  dSpMem *sp;
  //int i;
  if(!PyArg_ParseTuple(args,"l",&sp)){
    printf("Usage: handle to sp object\n");
    return NULL;
  }
  if(sp->typ=='d')
    smFreeSparseMem(sp);
  else if(sp->typ=='f')
    smFreeSparseMemFloat((fSpMem*)sp);
  else{
    printf("Sparse object invalid\n");
    return NULL;
  }
  return Py_BuildValue("i",0);
}
static PyObject *sparseMatrixInsert(PyObject *self,PyObject *args){
  //insert a value into the sparse matrix.
  dSpMem *sp;
  unsigned int row,col;
  double data;
  int val;
  //int i;
  if(!PyArg_ParseTuple(args,"lIId",&sp,&row,&col,&data)){
    printf("Usage: sparse object, row, col, data\n");
    return NULL;
  }
  if(sp->typ=='d')
    val=smInsert(sp,row,col,data);
  else if(sp->typ=='f')
    val=smInsertFloat((fSpMem*)sp,row,col,(float)data);
  else{
    printf("Sparse object invalid\n");
    return NULL;
  }

  if(val==-1){
    printf("Illegal position %d, %d for this sparse object\n",row,col);
    return NULL;
  }
  return Py_BuildValue("i",val);
}

static PyObject *sparseMatrixGet(PyObject *self,PyObject *args){
  //get a value from the sparse matrix.
  dSpMem *sp;
  unsigned int row,col;
  double data;
  if(!PyArg_ParseTuple(args,"lII",&sp,&row,&col)){
    printf("Usage: sparse object, row, col\n");
    return NULL;
  }
  if(sp->typ=='d')
    data=smGet(sp,row,col);
  else if(sp->typ=='f')
    data=smGetFloat((fSpMem*)sp,row,col);
  else{
    printf("Sparse object invalid\n");
    return NULL;
  }
  return Py_BuildValue("d",data);
}
static PyObject *sparseMatrixInfo(PyObject *self,PyObject *args){
  //get a value from the sparse matrix.
  dSpMem *sp;
  fSpMem *fsp;
  if(!PyArg_ParseTuple(args,"l",&sp)){
    printf("Usage: sparse object\n");
    return NULL;
  }
  if(sp->typ=='d')
    return Py_BuildValue("kkkdkkki",(long)sp->ndata,(long)sp->rows,(long)sp->cols,sp->min,(long)sp->rowmin,(long)sp->indmin,(long)sp->cnt,sp->alloced);
  else if(sp->typ=='f'){
    fsp=(fSpMem*)sp;
    return Py_BuildValue("kkkfkkki",(long)fsp->ndata,(long)fsp->rows,(long)fsp->cols,fsp->min,(long)fsp->rowmin,(long)fsp->indmin,(long)fsp->cnt,fsp->alloced);
  }else{
    printf("Sparse object invalid\n");
    return NULL;
  }
}

typedef struct{
  float *Adata;
  unsigned int *Acolind;
  unsigned int *Aindptr;
  float *Bdata;
  unsigned int *Browind;
  unsigned int *Bindptr;
  int BindptrDims;
  
}theDataStruct;

typedef struct{
  theDataStruct *AB;
  float *Cdata;
  unsigned int *Ccolind;
  unsigned int *Cindptr;
  int n;//thread number
  int astart;//starting row for this thread.
  int aend;//ending row for this thread.
  unsigned int Cdatasize;
  int err;
 
}threadInfoStruct;


int runCsrDotCsc(void *info){//this is run in a thread by csrDotCsc.
  threadInfoStruct *threadInfo;
  threadInfo=(threadInfoStruct*)info;
  threadInfo->err=sparseDotSparse(threadInfo->AB->Adata,
		  threadInfo->AB->Acolind,
		  &(threadInfo->AB->Aindptr[threadInfo->astart]),
		  threadInfo->aend-threadInfo->astart+1,
		  threadInfo->AB->Bdata,
		  threadInfo->AB->Browind,
		  threadInfo->AB->Bindptr,
		  threadInfo->AB->BindptrDims,
		  &threadInfo->Cdata,
		  &threadInfo->Ccolind,
		  threadInfo->Cindptr,
		  &threadInfo->Cdatasize,
		  (int)(threadInfo->n));
  return 0;
}


static PyObject *csrDotCsc(PyObject *self,PyObject *args){
  //do dot product of a sparse csr and csc matrix.
  //Converting to allow greater than 2G elements - ie using longs...
  PyArrayObject *Adata,*Acolind,*Aindptr,*Bdata,*Browind,*Bindptr;
  PyObject *CdataArr,*CcolindArr,*CindptrArr;
  float *Cdata;
  unsigned int *Ccolind;
  //unsigned long *Cindptr=NULL;
  //unsigned int *CindptrInt=NULL;
  unsigned int *Cindptr=NULL;
  threadInfoStruct *threadInfo;
  theDataStruct thedata;
  pthread_t *thread;
  int nthreads=1;
  npy_intp dims;
  int err=0;
  int arows,astart,aend,i,j,ndata;
  unsigned int indadd,totdata;
  unsigned int Cnzmax=0;
  long ltmp;
  //unsigned int CnzmaxInt;
  if(!PyArg_ParseTuple(args,"O!O!O!O!O!O!|i",&PyArray_Type,&Adata,&PyArray_Type,&Acolind,&PyArray_Type,&Aindptr,&PyArray_Type,&Bdata,&PyArray_Type,&Browind,&PyArray_Type,&Bindptr,/*&PyArray_Type,&Cdata,&PyArray_Type,&Ccolind,&PyArray_Type,&Cindptr,*/&nthreads)){
    printf("Usage: A.data,A.colind,A.indptr, B.data,B.rowind,B.indptr, nthreads (optional) where A is csr sparse format, B is csc sparse format and C is csr, with space for elements A.ndata + B.ndata.  This computes A dot B and puts result in C.\n");
    return NULL;
  }
  if(!testCF1(Adata) || !testCI1(Acolind) || !testCI1(Aindptr) ||
     !testCF1(Bdata) || !testCI1(Browind) || !testCI1(Bindptr)
     /*||  !testCF1(Cdata) || !testCI1(Ccolind) || !testCI1(Cindptr)*/
     ){
    printf("Arrays must be float (data) or int (ind and ptr), contiguous and 1D\n");
    return NULL;
  }
  //just check sizes are compatible.
  /*
  if(Aindptr->dimensions[0]!=Cindptr->dimensions[0] || Ccolind->dimensions[0]<Acolind->dimensions[0]+Browind->dimensions[0]){
    printf("C is wrong shape or doesn't have enough storage\n");
    return NULL;//note, we don't care about the x size of C (C.shape[1]), because this just depends on storage space... likewise, the x size of A and y size of B don't have to match here, and result will be okay, provided data is not present for the illegal areas.
    }*/
  //Cnzmax=Cdata->dimensions[0];
  if(nthreads==1){
    Cnzmax=Acolind->dimensions[0]+Browind->dimensions[0];//Cdata->dimensions[0];
    Ccolind=(unsigned int*)malloc(sizeof(int)*Cnzmax);
    Cdata=(float*)malloc(sizeof(float)*Cnzmax);
    //CindptrInt=(unsigned int*)malloc(sizeof(int)*Aindptr->dimensions[0]);
    //if(sparseDotSparse((float*)Adata->data,(unsigned int*)Acolind->data,(unsigned int*)Aindptr->data,(int)Aindptr->dimensions[0],(float*)Bdata->data,(unsigned int*)Browind->data,(unsigned int*)Bindptr->data,(int)Bindptr->dimensions[0],&Cdata,&Ccolind,CindptrInt,&CnzmaxInt,0)!=0)
    Cindptr=(unsigned int*)malloc(sizeof(int)*Aindptr->dimensions[0]);
    if(sparseDotSparse((float*)Adata->data,(unsigned int*)Acolind->data,(unsigned int*)Aindptr->data,(int)Aindptr->dimensions[0],(float*)Bdata->data,(unsigned int*)Browind->data,(unsigned int*)Bindptr->data,(int)Bindptr->dimensions[0],&Cdata,&Ccolind,Cindptr,&Cnzmax,0)!=0)
      return NULL;
  }else{//split the problem row wise between threads.
    thread=malloc(sizeof(pthread_t)*nthreads);

    threadInfo=(threadInfoStruct*)malloc(sizeof(threadInfoStruct)*nthreads);
    arows=Aindptr->dimensions[0]-1;
    thedata.Adata=(float*)Adata->data;
    thedata.Acolind=(unsigned int*)Acolind->data;
    thedata.Aindptr=(unsigned int*)Aindptr->data;
    thedata.Bdata=(float*)Bdata->data;
    thedata.Browind=(unsigned int*)Browind->data;
    thedata.Bindptr=(unsigned int*)Bindptr->data;
    thedata.BindptrDims=(int)Bindptr->dimensions[0];
    astart=0;
    for(i=0; i<nthreads; i++){
      aend=astart+(arows-astart)/(nthreads-i);
      if(i==nthreads-1)
	aend=arows;
      threadInfo[i].AB=&thedata;
      threadInfo[i].astart=astart;//starting row for this thread
      threadInfo[i].aend=aend;//ending row.
      threadInfo[i].Cdatasize=( ((unsigned int*)Bindptr->data)[Bindptr->dimensions[0]-1] + ((unsigned int*)Aindptr->data)[aend]-((unsigned int*)Aindptr->data)[astart]);
      printf("%d %d %d\n",astart,aend,threadInfo[i].Cdatasize);
      threadInfo[i].Cdata=(float*)malloc(sizeof(float)*threadInfo[i].Cdatasize);
      threadInfo[i].Ccolind=(unsigned int*)malloc(sizeof(int)*threadInfo[i].Cdatasize);
      threadInfo[i].Cindptr=(unsigned int*)malloc(sizeof(int)*(aend-astart+1));
      threadInfo[i].Cindptr[0]=0;
      threadInfo[i].n=i;
      pthread_create(&thread[i],NULL,(void*)runCsrDotCsc,&threadInfo[i]);
      astart=aend;
    }
    indadd=0;
    totdata=0;
    Cnzmax=0;
    for(i=0; i<nthreads; i++){
      pthread_join(thread[i],NULL);//wait for thread...
      printf("Joined %d\n",i);
      Cnzmax+=threadInfo[i].Cindptr[threadInfo[i].aend-threadInfo[i].astart];//datasize;
      err+=threadInfo[i].err;
    }
    if(err!=0){
      printf("Error non-zero\n");
      return NULL;
    }
    printf("All threads finished - now collating memory Cnzmax=%u\n",Cnzmax);
    //now allocate the memory...
    Ccolind=(unsigned int*)malloc(sizeof(int)*Cnzmax);
    Cdata=(float*)malloc(sizeof(float)*Cnzmax);
    Cindptr=(unsigned int*)malloc(sizeof(unsigned int)*Aindptr->dimensions[0]);
    if(Ccolind==NULL || Cdata==NULL || Cindptr==NULL){
      printf("Error allocating final memory: %p %p %p\n",Ccolind,Cdata,Cindptr);
      return NULL;
    }
    printf("Allocated okay\n");
    for(i=0; i<nthreads; i++){
      //collate C, and freem memory
      aend=threadInfo[i].aend;
      astart=threadInfo[i].astart;
      printf("astart %d aend %d\n",astart,aend);
      for(j=astart; j<=aend; j++){
	ltmp=(long)threadInfo[i].Cindptr[j-astart]+(long)indadd;
	if(ltmp>0xffffffff){
	  printf("Error - unsigned int overflow - too many elements for csrDotCsc\n");
	  return NULL;
	}
	Cindptr[j]=(unsigned int)ltmp;
      }
      printf("done for loop\n");
      ndata=threadInfo[i].Cindptr[aend-astart];
      printf("ndata %d\n",ndata);
      memcpy(&(Cdata[totdata]),threadInfo[i].Cdata,sizeof(float)*ndata);
      printf("done memcpy1\n");
      memcpy(&(Ccolind[totdata]),threadInfo[i].Ccolind,sizeof(float)*ndata);
      printf("done memcpy2\n");
      indadd=Cindptr[aend];
      totdata+=ndata;
      printf("indadd %u totdata %u, now freeing\n",indadd,totdata);
      free(threadInfo[i].Cdata);
      printf("free1\n");
      free(threadInfo[i].Ccolind);
      printf("free2\n");
      free(threadInfo[i].Cindptr);
      printf("Collated thread %d\n",i);
    }
    free(threadInfo);
    free(thread);
    //Cdata->dimensions[0]=((int*)(Cindptr->data))[];
  }
  dims=Cnzmax;
  CdataArr=PyArray_SimpleNewFromData(1,&dims,NPY_FLOAT,Cdata);
  CcolindArr=PyArray_SimpleNewFromData(1,&dims,NPY_INT,Ccolind);
  dims=Aindptr->dimensions[0];
  CindptrArr=PyArray_SimpleNewFromData(1,&dims,NPY_UINT,Cindptr);
  return Py_BuildValue("OOO",CdataArr,CcolindArr,CindptrArr);
}

typedef struct{
  float *Adata;
  int dim0;
  int dim1;
  float *bdata;
  unsigned int *browind;
  unsigned int *bindptr;
  int bdim;
  float *resdata;
  unsigned int *resrcind;
  unsigned int *resindptr;
  unsigned int resn;
  unsigned int resnmax;
  float valmin;
  int ctype;
  int end;
  int start;
  int n;
  int err;
}ddcStruct;//denseDotCsc struct

int runDenseDotCSC(void *info){//this is run in a thread by pyDenseDotCsc.
  ddcStruct *threadInfo;
  threadInfo=(ddcStruct*)info;
  //printf("thread started\n");

  threadInfo->err=denseDotCSC(threadInfo->Adata,threadInfo->dim0,threadInfo->dim1,threadInfo->bdata,threadInfo->browind,threadInfo->bindptr,threadInfo->bdim,&(threadInfo->resdata),&(threadInfo->resrcind),threadInfo->resindptr,&(threadInfo->resn),threadInfo->valmin,threadInfo->ctype,threadInfo->resnmax,threadInfo->n==0);


  
  return 0;
}


static PyObject *pyDenseDotCSC(PyObject *self,PyObject *args){
  //do dot product of a dense and csc matrix.
  PyArrayObject *A,*Bdata,*Bindptr,*Browind;
  PyObject *B,*shape;
  float *resdata;
  unsigned int *resrcind,*resindptr=NULL;
  ddcStruct *threadInfo;
  pthread_t *thread;
  int nthreads=1;
  npy_intp dims[2];
  float valmin=0;
  int ctype=0;
  unsigned int resn;
  int err=0;
  PyObject *CdataArr,*CcolindArr,*CindptrArr;
  int start,end,i,j,indadd,totdata,ndata,brows,bcols,rows;
  unsigned int resnmax=0;
  if(!PyArg_ParseTuple(args,"O!O|fiiI",&PyArray_Type,&A,&B,&valmin,&ctype,&nthreads,&resnmax)){
    printf("Usage: A (dense), B (csc), min value (optional), result type (optional, 0=csc, 1=csr, 2=dense), nthreads (optional), max number of elements (optional, only used for sparse types). This computes A dot B and returns result.\n");
    return NULL;
  }
  Bdata=(PyArrayObject*)PyObject_GetAttrString(B,"data");
  Bindptr=(PyArrayObject*)PyObject_GetAttrString(B,"indptr");
  Browind=(PyArrayObject*)PyObject_GetAttrString(B,"rowind");
  shape=PyObject_GetAttrString(B,"shape");
  if(Bdata==NULL || Bindptr==NULL || Browind==NULL || shape==NULL){
    printf("B mst be CSC\n");
    return NULL;
  }
  Py_DECREF(Bdata);//free it here, since we return at myriad points...
  Py_DECREF(Bindptr);
  Py_DECREF(Browind);
  Py_DECREF(shape);
  if(!PyArg_ParseTuple(shape,"ii",&brows,&bcols)){
    printf("Shape attribute of B not understood\n");
    return NULL;
  }

  if(!testCF2(A) || !testCF1(Bdata) || !(testCUI1(Bindptr)||testCI1(Bindptr)) || !(testCUI1(Browind)||testCI1(Browind))){
    printf("A must be 2D float contiguous, and B must be float32 with data, rowind and indptr contiguous\n");
    return NULL;
  }
  if(resnmax==0)
    resnmax=0xffffffff;
  if(nthreads==1){
    if(ctype==0 || ctype==1){//csc or csr
      //resn=Bdata->dimensions[0]/((float)(brows*bcols))*A->dimensions[0]*bcols;  //chosen such that sparsity of new matrix is same as B matrix.
      resn=(unsigned int)((Bdata->dimensions[0]/(float)brows)*A->dimensions[0]);
      if((resdata=malloc(sizeof(float)*resn))==NULL){
	printf("Unable to malloc resdata\n");
	return NULL;
      }
      if((resrcind=malloc(sizeof(int)*resn))==NULL){
	printf("Unable to malloc resrcind\n");
	return NULL;
      }
      if((resindptr=malloc(sizeof(int)*(ctype==0?Bindptr->dimensions[0]:A->dimensions[0]+1)))==NULL){
	printf("Unable to malloc resindptr\n");
	return NULL;
      }
    }else{
      if((resdata=malloc(sizeof(float)*A->dimensions[0]*bcols/*(Bindptr->dimensions[0]-1)*/))==NULL){
	printf("Unable to allocate memory for result\n");
	return NULL;
      }
      //memset(resdata,0,sizeof(float)*A->dimensions[0]*bcols);
      resn=0;
      resrcind=NULL;
      resindptr=NULL;
    }
    err=denseDotCSC((float*)A->data,(int)A->dimensions[0],(int)A->dimensions[1],(float*)Bdata->data,(unsigned int*)Browind->data,(unsigned int*)Bindptr->data,(int)Bindptr->dimensions[0],&resdata,&resrcind,resindptr,&resn,valmin,ctype,resnmax,A->dimensions[0]>1024);
  
  }else{//multi-threaded version

    thread=malloc(sizeof(pthread_t)*nthreads);

    threadInfo=(ddcStruct*)malloc(sizeof(ddcStruct)*nthreads);
    start=0;
    rows=(ctype==0?bcols:A->dimensions[0]);
    if(ctype!=0 && ctype!=1){
      if((resdata=malloc(sizeof(float)*A->dimensions[0]*bcols/*(Bindptr->dimensions[0]-1)*/))==NULL){
	printf("Unable to allocate memory for result\n");
	return NULL;
      }
      printf("Malloced resdata size %ld bytes 0x%p\n",sizeof(float)*A->dimensions[0]*bcols,resdata);
      //memset(resdata,0,sizeof(float)*A->dimensions[0]*bcols);
      resn=0;
      resrcind=NULL;
      resindptr=NULL;
	
    }
    for(i=0; i<nthreads; i++){
      end=start+(rows-start)/(nthreads-i);
      if(i==nthreads-1)
	end=rows;
      threadInfo[i].end=end;
      threadInfo[i].start=start;
      printf("%d %d\n",start,end);
      if(ctype==0){//result is csc
	threadInfo[i].Adata=(float*)A->data;
	threadInfo[i].dim0=(int)A->dimensions[0];
	threadInfo[i].dim1=(int)A->dimensions[1];
	threadInfo[i].bdata=(float*)Bdata->data;
	threadInfo[i].browind=(unsigned int*)Browind->data;
	threadInfo[i].bindptr=&(((unsigned int*)Bindptr->data)[start]);
	threadInfo[i].bdim=end-start+1;//(int)Bindptr->dimensions[0];
	threadInfo[i].resnmax=resnmax/nthreads;
	//resn=Bdata->dimensions[0]/((float)(brows*bcols))*A->dimensions[0]*(end-start);//chosen such that sparsity of new is same as B.
	resn=(unsigned int)((Bdata->dimensions[0]/(float)(brows*bcols))*A->dimensions[0]*(end-start));
	if((resdata=malloc(sizeof(float)*resn))==NULL){
	  printf("Unable to malloc resdata %ld\n",(long)resn);
	  return NULL;
	}
	if((resrcind=malloc(sizeof(int)*resn))==NULL){
	  printf("Unable to malloc resrcind\n");
	  return NULL;
	}
	if((resindptr=malloc(sizeof(int)*(end-start+1)))==NULL){
	  printf("Unable to malloc resindptr\n");
	  return NULL;
	}
	threadInfo[i].resdata=resdata;
	threadInfo[i].resrcind=resrcind;
	threadInfo[i].resindptr=resindptr;
	threadInfo[i].resn=resn;

	
      }else{//result is csr or dense...
	threadInfo[i].Adata=&(((float*)A->data)[start*A->dimensions[1]]);
	threadInfo[i].dim0=end-start;
	threadInfo[i].dim1=(int)A->dimensions[1];
	threadInfo[i].bdata=(float*)Bdata->data;
	threadInfo[i].browind=(unsigned int*)Browind->data;
	threadInfo[i].bindptr=(unsigned int*)Bindptr->data;
	threadInfo[i].bdim=(int)Bindptr->dimensions[0];
	threadInfo[i].resnmax=resnmax/nthreads;
	if(ctype==1){//csr result
	  //resn=Bdata->dimensions[0]/((float)(brows*bcols))*(end-start)*bcols;//chosen such that sparsity of new is same as B.
	  resn=(unsigned int)((Bdata->dimensions[0]/(float)brows)*(end-start));
	  //printf("resn: %u %ld %d %d %d %d\n",resn,(long)Bdata->dimensions[0],brows,bcols,end-start,bcols);
	  if((resdata=malloc(sizeof(float)*resn))==NULL){
	    printf("Unable to malloc resdata %ld\n",(long)resn);
	    return NULL;
	  }
	  if((resrcind=malloc(sizeof(int)*resn))==NULL){
	    printf("Unable to malloc resrcind\n");
	    return NULL;
	  }
	  if((resindptr=malloc(sizeof(int)*(end-start+1)))==NULL){
	    printf("Unable to malloc resindptr\n");
	    return NULL;
	  }
	  threadInfo[i].resdata=resdata;
	  threadInfo[i].resrcind=resrcind;
	  threadInfo[i].resindptr=resindptr;
	  threadInfo[i].resn=resn;
	  
	}else{//dense result
	  threadInfo[i].resdata=&(resdata[(long)start*(long)bcols]);
	  threadInfo[i].resrcind=NULL;
	  threadInfo[i].resindptr=NULL;
	  threadInfo[i].resn=0;
	}
      }
      threadInfo[i].valmin=valmin;
      threadInfo[i].ctype=ctype;
      threadInfo[i].n=i;
      start=end;
      pthread_create(&thread[i],NULL,(void*)runDenseDotCSC,&threadInfo[i]);
    }
    resn=0;
    for(i=0; i<nthreads; i++){
      pthread_join(thread[i],NULL);//wait for thread...
      if(ctype==0 || ctype==1)
	resn+=threadInfo[i].resindptr[threadInfo[i].end-threadInfo[i].start];//datasize;
      if(threadInfo[i].err<0)
	err=-1;
      else if(err>=0 && threadInfo[i].err>0)
	err=(err==0?threadInfo[i].err:(err<threadInfo[i].err?err:threadInfo[i].err));
    }
    resrcind=NULL;
    resdata=NULL;
    resindptr=NULL;
    if(err==0 && (ctype==0 || ctype==1)){//now allocate the csc/csr memory.
      indadd=0;
      totdata=0;
      resrcind=(unsigned int*)malloc(sizeof(unsigned int)*resn);
      resdata=(float*)malloc(sizeof(float)*resn);
      resindptr=(unsigned int*)malloc(sizeof(unsigned int)*(ctype==0?Bindptr->dimensions[0]:A->dimensions[0]+1));
      for(i=0; i<nthreads; i++){
	//collate C, and freem memory
	end=threadInfo[i].end;
	start=threadInfo[i].start;
	for(j=start; j<=end; j++){
	  resindptr[j]=threadInfo[i].resindptr[j-start]+indadd;
	}
	ndata=threadInfo[i].resindptr[end-start];
	memcpy(&(resdata[totdata]),threadInfo[i].resdata,sizeof(float)*ndata);
	memcpy(&(resrcind[totdata]),threadInfo[i].resrcind,sizeof(float)*ndata);
	indadd=resindptr[end];
	totdata+=ndata;
	free(threadInfo[i].resdata);
	free(threadInfo[i].resrcind);
	free(threadInfo[i].resindptr);
      }
    }
    free(thread);
    free(threadInfo);


  
    

  }
  if(err<0)//don't bother freeing, just raise an exception.
    return NULL;
  else if(err>0){
    //failed - reached max elements.  err is the number of rows that we'd got to.  Return this, which may be useful for the calling function (ie to change valmin)
    if(resdata!=NULL)free(resdata);
    if(resrcind!=NULL)free(resrcind);
    if(resindptr!=NULL)free(resindptr);
    return Py_BuildValue("i",err);
  }



  if(ctype==0 || ctype==1){
    dims[0]=resn;
    CdataArr=PyArray_SimpleNewFromData(1,dims,NPY_FLOAT,resdata);
    CcolindArr=PyArray_SimpleNewFromData(1,dims,NPY_UINT,resrcind);
    memset(&resrcind[resindptr[(ctype==0?bcols/*Bindptr->dimensions[0]-1*/:A->dimensions[0])]],0,sizeof(float)*(resn-resindptr[(ctype==0?bcols/*Bindptr->dimensions[0]-1*/:A->dimensions[0])]));
    dims[0]=(ctype==0?Bindptr->dimensions[0]:A->dimensions[0]+1);
    CindptrArr=PyArray_SimpleNewFromData(1,dims,NPY_UINT,resindptr);
    return Py_BuildValue("OOO",CdataArr,CcolindArr,CindptrArr);
  }else{
    dims[0]=A->dimensions[0];
    dims[1]=bcols/*Bindptr->dimensions[0]-1*/;
    CdataArr=PyArray_SimpleNewFromData(2,dims,NPY_FLOAT,resdata);
    return Py_BuildValue("O",CdataArr);
  }
}

static PyObject *pyCountInstances(PyObject *self,PyObject *args){
  //count number of values in array bigger than the selection.
  PyArrayObject *A,*vals,*cnt;
  int i;
  long size;
  if(!PyArg_ParseTuple(args,"O!O!O!",&PyArray_Type,&A,&PyArray_Type,&vals,&PyArray_Type,&cnt)){
    printf("Usage: A (dense), vals (sorted from small to large), cnt\n");
    return NULL;
  }
  if(!testCF(A) || !testCF1(vals) || !testCL1(cnt) || cnt->dimensions[0]!=vals->dimensions[0]){
    printf("Arrays must be contiguous and 2D (data) or 1D (vals and cnt)\n");
    return NULL;
  }
  size=1;
  for(i=0; i<A->nd; i++)
    size*=A->dimensions[i];
  countInstances((float*)A->data,size,(float*)vals->data,(long*)cnt->data, vals->dimensions[0]);
  Py_RETURN_NONE;
}
static PyObject *pySparsifyCsr(PyObject *self,PyObject *args){
  //sparsify dense to sparse.
  PyArrayObject *A,*data,*colind,*indptr;
  PyObject *csr,*shape;//,*tup;
  //PyObject *nnz;
  int brows,bcols;
  unsigned int nused;
  float val;
  if(!PyArg_ParseTuple(args,"O!Of",&PyArray_Type,&A,&csr,&val)){
    printf("Usage: Array, CSR object, min value allowed\nThe CSR object should be large enough to hold the sparsified data.\n");
    return NULL;
  }
  if(A->descr->type_num!=NPY_FLOAT || A->nd!=2){
    printf("A must be float32 and 2D\n");
    return NULL;
  }
  data=(PyArrayObject*)PyObject_GetAttrString(csr,"data");
  indptr=(PyArrayObject*)PyObject_GetAttrString(csr,"indptr");
  colind=(PyArrayObject*)PyObject_GetAttrString(csr,"colind");
  shape=PyObject_GetAttrString(csr,"shape");
  if(data==NULL || indptr==NULL || colind==NULL || shape==NULL){
    printf("Sparse array mst be CSR format\n");
    return NULL;
  }
  Py_DECREF(data);//free it here, since we return at myriad points...
  Py_DECREF(indptr);
  Py_DECREF(colind);
  Py_DECREF(shape);
  if(!PyArg_ParseTuple(shape,"ii",&brows,&bcols)){
    printf("Shape attribute of CSR not understood\n");
    return NULL;
  }
  if(!testCF1(data) ||  !(testCUI1(indptr)||testCI1(indptr)) || !(testCUI1(colind)||testCI1(colind)) || indptr->dimensions[0]!=A->dimensions[0]+1){
    printf("Wrong data types or not contiguous and 1D - csr.data,indptr,colind.\n");
    return NULL;
  }

  if(sparsifyCsr((float*)A->data,(int)A->dimensions[0],(int)A->dimensions[1],(int)A->strides[0]/sizeof(float),(int)A->strides[1]/sizeof(float),(float*)data->data,(unsigned int*)colind->data,(unsigned int*)indptr->data,(unsigned int)data->dimensions[0],val)!=0){
    printf("Sparsify failed - out of space\n");
    return NULL;
  }
  printf("Done sparsification... cleaning up\n");
  nused=((unsigned int*)indptr->data)[indptr->dimensions[0]-1];
  //zero the rest of colind.
  if(nused<colind->dimensions[0]){
    memset(&(((unsigned int*)colind->data)[nused]),0,sizeof(unsigned int*)*(colind->dimensions[0]-nused));
  }
  //printf("nused %d\n",nused);
  //tup=Py_BuildValue("()");
  //printf("pyobject _ call _check()\n");
  //PyObject_Call(shape=PyObject_GetAttrString(csr,"_check"),tup,NULL);
  //printf("decref tup\n");
  //Py_DECREF(tup);
  //printf("done\n");
  //Py_DECREF(shape);
  /*PyObject_SetAttrString(csr,"nnz");
  if(PyLong_Check(nnz)){
    printf("islong\n");
    ((PyLongObject*)nnz)->ob_ival=(long)nused;
    }*/
  Py_RETURN_NONE;
}

static PyObject *pyDensifyCsr(PyObject *self,PyObject *args){
  //sparsify dense to sparse.
  PyArrayObject *A,*data,*colind,*indptr;
  PyObject *csr,*shape;
  //PyObject *arrdata,*arrcolind,*arrindptr;
  //PyObject *nnz;
  int brows,bcols;
  if(!PyArg_ParseTuple(args,"O!O",&PyArray_Type,&A,&csr)){
    printf("Usage: Array (zero'd), CSR object\nThe array should have same shape as the CSR object\n");
    return NULL;
  }
  if(A->descr->type_num!=NPY_FLOAT || A->nd!=2){
    printf("A must be float32 and 2D\n");
    return NULL;
  }
  data=(PyArrayObject*)PyObject_GetAttrString(csr,"data");
  indptr=(PyArrayObject*)PyObject_GetAttrString(csr,"indptr");
  colind=(PyArrayObject*)PyObject_GetAttrString(csr,"colind");
  shape=PyObject_GetAttrString(csr,"shape");
  if(data==NULL || indptr==NULL || colind==NULL || shape==NULL){
    printf("Sparse array mst be CSR format\n");
    return NULL;
  }
  Py_DECREF(data);//free it here, since we return at myriad points...
  Py_DECREF(indptr);
  Py_DECREF(colind);
  Py_DECREF(shape);
  if(!PyArg_ParseTuple(shape,"ii",&brows,&bcols)){
    printf("Shape attribute of CSR not understood\n");
    return NULL;
  }
  if(!testCF1(data) ||  !(testCUI1(indptr)||testCI1(indptr)||testCL1(indptr)||testCUL1(indptr)) || !(testCUI1(colind)||testCI1(colind)) || indptr->dimensions[0]!=A->dimensions[0]+1){
    printf("Wrong data types or not contiguous and 1D - csr.data,indptr,colind.\n");
    return NULL;
  }
  if(indptr->descr->elsize==sizeof(int))
    densifyCsr((float*)A->data,(int)A->dimensions[0],(int)A->dimensions[1],(int)A->strides[0]/sizeof(float),A->strides[1]/sizeof(float),(float*)data->data,(unsigned int*)colind->data,(unsigned int*)indptr->data,(unsigned int)data->dimensions[0]);
  else
    densifyCsrL((float*)A->data,(int)A->dimensions[0],(int)A->dimensions[1],(int)A->strides[0]/sizeof(float),A->strides[1]/sizeof(float),(float*)data->data,(unsigned int*)colind->data,(unsigned long*)indptr->data,(unsigned int)data->dimensions[0]);
  Py_RETURN_NONE;
}


static PyMethodDef svdMethods[] = {
  //{"agbsvd",  agbsvd, METH_VARARGS,"Run the SVD algorithm."},
  {"dolibsvd",  dolibsvd, METH_VARARGS,"Run the SVD algorithm from libsvd."},
  {"computeDenseGenInv",computeDenseGenInv,METH_VARARGS,"Compute the generalised inverse"},
  {"sparsifyGenInv",sparsifyGenInv,METH_VARARGS,"Sparsify the dense generalised inverse"},
  {"sparseMatrixCreate",sparseMatrixCreate,METH_VARARGS,"Create a sparse matrix object"},
  {"sparseMatrixFree",sparseMatrixFree,METH_VARARGS,"Free the sparse matrix object"},
  {"sparseMatrixInsert",sparseMatrixInsert,METH_VARARGS,"Insert a value into the sparse matrix"},
  {"sparseMatrixGet",sparseMatrixGet,METH_VARARGS,"Get a value from the sparse matrix"},
  {"sparseMatrixInfo",sparseMatrixInfo,METH_VARARGS,"Get info from the sparse matrix struct"},
  {"csrDotCsc",csrDotCsc,METH_VARARGS,"Dot product of csr and csc matrix"},
  {"denseDotCsc",pyDenseDotCSC,METH_VARARGS,"Dot product of dense and csc matrix"},
  {"countInstances",pyCountInstances,METH_VARARGS,"Count number of occurrances greater than set values in a dense array"},
  {"sparsifyCsr",pySparsifyCsr,METH_VARARGS,"Sparsify a dense matrix to CSR"},
  {"densifyCsr",pyDensifyCsr,METH_VARARGS,"Densify a CSR matrix"},
  {NULL, NULL, 0, NULL}        /* Sentinel */
};
//PyMODINIT_FUNC 
void initsvd(void)
{
  PyObject *m;
  PyImport_AddModule("svd");
  m=Py_InitModule("svd", svdMethods);
  import_array();
  SvdError = PyErr_NewException("svd.error", NULL, NULL);
  Py_INCREF(SvdError);
  PyModule_AddObject(m, "error", SvdError);
}


int
main(int argc, char *argv[])
{
    // Pass argv[0] to the Python interpreter 
  Py_SetProgramName(argv[0]);
  
  // Initialize the Python interpreter.  Required. 
  Py_Initialize();
  
  // Add a static module 
  initsvd();
  return 0;
}


