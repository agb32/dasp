#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <netinet/in.h>
#include "sparsemem.h"
#include "svdlib.h"
#include "svdutil.h"
#include <pthread.h>
#define BUNZIP2  "bzip2 -d"
#define BZIP2    "bzip2 -1"
#define UNZIP    "gzip -d"
#define ZIP      "gzip -1"
#define COMPRESS "compress"

#define MAX_FILENAME 512
#define MAX_PIPES    64
static FILE *Pipe[MAX_PIPES];
static int numPipes = 0;

long *svd_longArray(long size, char empty, char *name) {
  long *a;
  if (empty) a = (long *) calloc(size, sizeof(long));
  else a = (long *) malloc(size * sizeof(long));
  if (!a) {
    perror(name);
    /* exit(errno); */
  }
  return a;
}
unsigned int *svd_uintArray(long size, char empty, char *name) {
  unsigned int *a;
  if (empty) a = (unsigned int *) calloc(size, sizeof(unsigned int));
  else a = (unsigned int *) malloc(size * sizeof(unsigned int));
  if (!a) {
    perror(name);
    /* exit(errno); */
  }
  return a;
}

double *svd_doubleArray(long size, char empty, char *name) {
  double *a;
  if (empty) a = (double *) calloc(size, sizeof(double));
  else a = (double *) malloc(size * sizeof(double));
  if (!a) {
    perror(name);
    /* exit(errno); */
  }
  return a;
}
float *svd_floatArray(long size, char empty, char *name) {
  float *a;
  if (empty) a = (float *) calloc(size, sizeof(float));
  else a = (float *) malloc(size * sizeof(float));
  if (!a) {
    perror(name);
    /* exit(errno); */
  }
  return a;
}


void svd_beep(void) {
  fputc('\a', stderr);
  fflush(stderr);
}

void svd_debug(char *fmt, ...) {
  va_list args;
  va_start(args, fmt);
  vfprintf(stderr, fmt, args);
  va_end(args);
}

void svd_error(char *fmt, ...) {
  va_list args;
  va_start(args, fmt);
  svd_beep();
  fprintf(stderr, "ERROR: ");
  vfprintf(stderr, fmt, args);
  fprintf(stderr, "\n");
  va_end(args);
}

void svd_fatalError(char *fmt, ...) {
  va_list args;
  va_start(args, fmt);
  svd_beep();
  fprintf(stderr, "ERROR: ");
  vfprintf(stderr, fmt, args);
  fprintf(stderr, "\a\n");
  va_end(args);
  exit(1);
}

static void registerPipe(FILE *p) {
  if (numPipes >= MAX_PIPES) svd_error("Too many pipes open");
  Pipe[numPipes++] = p;
}

static char isPipe(FILE *p) {
  int i;
  for (i = 0; i < numPipes && Pipe[i] != p; i++);
  if (i == numPipes) return FALSE;
  Pipe[i] = Pipe[--numPipes];
  return TRUE;
}

static FILE *openPipe(char *pipeName, char *mode) {
  FILE *pipe;
  fflush(stdout);
  if ((pipe = popen(pipeName, mode))) registerPipe(pipe);
  return pipe;
}

static FILE *readZippedFile(char *command, char *fileName) {
  char buf[MAX_FILENAME];
  sprintf(buf, "%s < %s 2>/dev/null", command, fileName);
  return openPipe(buf, "r");
}

FILE *svd_fatalReadFile(char *filename) {
  FILE *file;
  if (!(file = svd_readFile(filename)))
    svd_fatalError("couldn't read the file %s", filename);
  return file;
}

static int stringEndsIn(char *s, char *t) {
  int ls = strlen(s);
  int lt = strlen(t);
  if (ls < lt) return FALSE;
  return (strcmp(s + ls - lt, t)) ? FALSE : TRUE;
}

/* Will silently return NULL if file couldn't be opened */
FILE *svd_readFile(char *fileName) {
  char fileBuf[MAX_FILENAME];
  struct stat statbuf;

  /* Special file name */
  if (!strcmp(fileName, "-"))
    return stdin;
  
  /* If it is a pipe */
  if (fileName[0] == '|')
    return openPipe(fileName + 1, "r");

  /* Check if already ends in .gz or .Z and assume compressed */
  if (stringEndsIn(fileName, ".gz") || stringEndsIn(fileName, ".Z")) {
    if (!stat(fileName, &statbuf))
      return readZippedFile(UNZIP, fileName);
    return NULL;
  }
  /* Check if already ends in .bz or .bz2 and assume compressed */
  if (stringEndsIn(fileName, ".bz") || stringEndsIn(fileName, ".bz2")) {
    if (!stat(fileName, &statbuf))
      return readZippedFile(BUNZIP2, fileName);
    return NULL;
  }
  /* Try just opening normally */
  if (!stat(fileName, &statbuf))
    return fopen(fileName, "r");
  /* Try adding .gz */
  sprintf(fileBuf, "%s.gz", fileName);
  if (!stat(fileBuf, &statbuf))
    return readZippedFile(UNZIP, fileBuf);
  /* Try adding .Z */
  sprintf(fileBuf, "%s.Z", fileName);
  if (!stat(fileBuf, &statbuf))
    return readZippedFile(UNZIP, fileBuf);
  /* Try adding .bz2 */
  sprintf(fileBuf, "%s.bz2", fileName);
  if (!stat(fileBuf, &statbuf))
    return readZippedFile(BUNZIP2, fileBuf);
  /* Try adding .bz */
  sprintf(fileBuf, "%s.bz", fileName);
  if (!stat(fileBuf, &statbuf))
    return readZippedFile(BUNZIP2, fileBuf);

  return NULL;
}

static FILE *writeZippedFile(char *fileName, char append) {
  char buf[MAX_FILENAME];
  const char *op = (append) ? ">>" : ">";
  if (stringEndsIn(fileName, ".bz2") || stringEndsIn(fileName, ".bz"))
    sprintf(buf, "%s %s \"%s\"", BZIP2, op, fileName);
  else if (stringEndsIn(fileName, ".Z"))
    sprintf(buf, "%s %s \"%s\"", COMPRESS, op, fileName);
  else
    sprintf(buf, "%s %s \"%s\"", ZIP, op, fileName);
  return openPipe(buf, "w");
}

FILE *svd_writeFile(char *fileName, char append) {
  /* Special file name */
  if (!strcmp(fileName, "-"))
    return stdout;
  
  /* If it is a pipe */
  if (fileName[0] == '|')
    return openPipe(fileName + 1, "w");

  /* Check if ends in .gz, .Z, .bz, .bz2 */
  if (stringEndsIn(fileName, ".gz") || stringEndsIn(fileName, ".Z") ||
      stringEndsIn(fileName, ".bz") || stringEndsIn(fileName, ".bz2"))
    return writeZippedFile(fileName, append);
  return (append) ? fopen(fileName, "a") : fopen(fileName, "w");
}

/* Could be a file or a stream. */
void svd_closeFile(FILE *file) {
  if (file == stdin || file == stdout) return;
  if (isPipe(file)) pclose(file);
  else fclose(file);
}


char svd_readBinInt(FILE *file, int *val) {
  int x;
  if (fread(&x, sizeof(int), 1, file) == 1) {
    *val = ntohl(x);
    return FALSE;
  }
  return TRUE;
}

/* This reads a float in network order and converts to a real in host order. */
char svd_readBinFloat(FILE *file, float *val) {
  int x;
  float y;
  if (fread(&x, sizeof(int), 1, file) == 1) {
    x = ntohl(x);
    y = *((float *) &x);
    *val = y;
    return FALSE;
  }
  return TRUE;
}

char svd_writeBinInt(FILE *file, int x) {
  int y = htonl(x);
  if (fwrite(&y, sizeof(int), 1, file) != 1) return TRUE;
  return FALSE;
}

/* This takes a real in host order and writes a float in network order. */
char svd_writeBinFloat(FILE *file, float r) {
  int y = htonl(*((int *) &r));
  if (fwrite(&y, sizeof(int), 1, file) != 1) return TRUE;
  return FALSE;
}


/************************************************************** 
 * returns |a| if b is positive; else fsign returns -|a|      *
 **************************************************************/ 
double svd_fsign(double a, double b) {
  if ((a>=0.0 && b>=0.0) || (a<0.0 && b<0.0))return(a);
  else return -a;
}

/************************************************************** 
 * returns the larger of two double precision numbers         *
 **************************************************************/ 
double svd_dmax(double a, double b) {
   return (a > b) ? a : b;
}

/************************************************************** 
 * returns the smaller of two double precision numbers        *
 **************************************************************/ 
double svd_dmin(double a, double b) {
  return (a < b) ? a : b;
}

/************************************************************** 
 * returns the larger of two integers                         *
 **************************************************************/ 
long svd_imax(long a, long b) {
  return (a > b) ? a : b;
}

/************************************************************** 
 * returns the smaller of two integers                        *
 **************************************************************/ 
long svd_imin(long a, long b) {
  return (a < b) ? a : b;
}

/************************************************************** 
 * Function scales a vector by a constant.     		      *
 * Based on Fortran-77 routine from Linpack by J. Dongarra    *
 **************************************************************/ 
void svd_dscal(long n, double da, double *dx, long incx) {
  long i;
  
  if (n <= 0 || incx == 0) return;
  if (incx < 0) dx += (-n+1) * incx;
  for (i=0; i < n; i++) {
    *dx *= da;
    dx += incx;
  }
  return;
}

void svd_ffscal(long n, float da, float *dx, long incx) {
  long i;
  
  if (n <= 0 || incx == 0) return;
  if (incx < 0) dx += (-n+1) * incx;
  for (i=0; i < n; i++) {
    *dx *= da;
    dx += incx;
  }
  return;
}

void svd_fdscal(long n, float da, double *dx, long incx) {
  long i;
  
  if (n <= 0 || incx == 0) return;
  if (incx < 0) dx += (-n+1) * incx;
  for (i=0; i < n; i++) {
    *dx *= (double)da;
    dx += incx;
  }
  return;
}
void svd_dfscal(long n, double da, float *dx, long incx) {
  long i;
  
  if (n <= 0 || incx == 0) return;
  if (incx < 0) dx += (-n+1) * incx;
  for (i=0; i < n; i++) {
    *dx *= (float)da;
    dx += incx;
  }
  return;
}
void svd_fffscalnew(long n, float da, float *dx, long incx,float *dy,long incy) {
  long i;
  
  if (n <= 0 || incx == 0 || incy==0) return;
  if (incx < 0) dx += (-n+1) * incx;
  if(incy<0)dy+=(-n+1)*incy;
  for (i=0; i < n; i++) {
    *dy =*dx* da;
    dx+=incx;
    dy+=incy;
  }
  return;
}
void svd_fddscalnew(long n, float da, double *dx,long incx,double *dy,long incy) {
  long i;
  
  if (n <= 0 || incx == 0 || incy==0) return;
  if (incx < 0) dx += (-n+1) * incx;
  if(incy<0)dy+=(-n+1)*incy;
  for (i=0; i < n; i++) {
    *dy =*dx* da;
    dx+=incx;
    dy+=incy;
  }
  return;
}


/************************************************************** 
 * function scales a vector by a constant.	     	      *
 * Based on Fortran-77 routine from Linpack by J. Dongarra    *
 **************************************************************/ 
void svd_datx(long n, double da, double *dx, long incx, double *dy, long incy) {
  long i;
  
  if (n <= 0 || incx == 0 || incy == 0 || da == 0.0) return;
  if (incx == 1 && incy == 1) 
    for (i=0; i < n; i++) *dy++ = da * (*dx++);
  
  else {
    if (incx < 0) dx += (-n+1) * incx;
    if (incy < 0) dy += (-n+1) * incy;
    for (i=0; i < n; i++) {
      *dy = da * (*dx);
      dx += incx;
      dy += incy;
    }
  }
  return;
}

/************************************************************** 
 * Function copies a vector x to a vector y	     	      *
 * Based on Fortran-77 routine from Linpack by J. Dongarra    *
 **************************************************************/ 
void svd_dcopy(long n, double *dx, long incx, double *dy, long incy) {
  long i;
  
  if (n <= 0 || incx == 0 || incy == 0) return;
  if (incx == 1 && incy == 1) 
    for (i=0; i < n; i++) *dy++ = *dx++;
  
  else {
    if (incx < 0) dx += (-n+1) * incx;
    if (incy < 0) dy += (-n+1) * incy;
    for (i=0; i < n; i++) {
      *dy = *dx;
      dx += incx;
      dy += incy;
    }
  }
  return;
}

void svd_dfcopy(long n, double *dx, long incx, float *dy, long incy) {
  long i;
  
  if (n <= 0 || incx == 0 || incy == 0) return;
  if (incx == 1 && incy == 1) 
    for (i=0; i < n; i++) *dy++ = (float)*dx++;
  
  else {
    if (incx < 0) dx += (-n+1) * incx;
    if (incy < 0) dy += (-n+1) * incy;
    for (i=0; i < n; i++) {
      *dy = (float)*dx;
      dx += incx;
      dy += incy;
    }
  }
  return;
}

void svd_fcopy(long n, float *dx, long incx, float *dy, long incy) {
  long i;
  
  if (n <= 0 || incx == 0 || incy == 0) return;
  if (incx == 1 && incy == 1) 
    for (i=0; i < n; i++) *dy++ = *dx++;
  
  else {
    if (incx < 0) dx += (-n+1) * incx;
    if (incy < 0) dy += (-n+1) * incy;
    for (i=0; i < n; i++) {
      *dy = *dx;
      dx += incx;
      dy += incy;
    }
  }
  return;
}
void svd_fdcopy(long n, float *dx, long incx, double *dy, long incy) {
  long i;
  
  if (n <= 0 || incx == 0 || incy == 0) return;
  if (incx == 1 && incy == 1) 
    for (i=0; i < n; i++) *dy++ = (double)(*dx++);
  
  else {
    if (incx < 0) dx += (-n+1) * incx;
    if (incy < 0) dy += (-n+1) * incy;
    for (i=0; i < n; i++) {
      *dy = (double)(*dx);
      dx += incx;
      dy += incy;
    }
  }
  return;
}

//Reverse a vector...
void vecReverse(long n,double *v){
  //agb func to reverse a vector...
  double tmp;
  long i,e=n/2;
  for(i=0; i<e; i++){
    tmp=v[i];
    v[i]=v[n-1-i];
    v[n-1-i]=tmp;
  }
}

/************************************************************** 
 * Function forms the dot product of two vectors.      	      *
 * Based on Fortran-77 routine from Linpack by J. Dongarra    *
 **************************************************************/ 
double svd_ddot(long n, double *dx, long incx, double *dy, long incy) {
  long i;
  double dot_product;
  
  if (n <= 0 || incx == 0 || incy == 0) return(0.0);
  dot_product = 0.0;
  if (incx == 1 && incy == 1) 
    for (i=0; i < n; i++) dot_product += (*dx++) * (*dy++);
  else {
    if (incx < 0) dx += (-n+1) * incx;
    if (incy < 0) dy += (-n+1) * incy;
    for (i=0; i < n; i++) {
      dot_product += (*dx) * (*dy);
      dx += incx;
      dy += incy;
      }
  }
  return(dot_product);
}
double svd_fddot(long n, float *dx, long incx, double *dy, long incy) {
  long i;
  double dot_product;
  
  if (n <= 0 || incx == 0 || incy == 0) return(0.0);
  dot_product = 0.0;
  if (incx == 1 && incy == 1) 
    for (i=0; i < n; i++) dot_product += (*dx++) * (*dy++);
  else {
    if (incx < 0) dx += (-n+1) * incx;
    if (incy < 0) dy += (-n+1) * incy;
    for (i=0; i < n; i++) {
      dot_product += (*dx) * (*dy);
      dx += incx;
      dy += incy;
      }
  }
  return(dot_product);
}
double svd_fdot(long n, float *dx, long incx, float *dy, long incy) {
  long i;
  double dot_product;
  
  if (n <= 0 || incx == 0 || incy == 0) return(0.0);
  dot_product = 0.0;
  if (incx == 1 && incy == 1) 
    for (i=0; i < n; i++) dot_product += (double)((*dx++) * (*dy++));
  else {
    if (incx < 0) dx += (-n+1) * incx;
    if (incy < 0) dy += (-n+1) * incy;
    for (i=0; i < n; i++) {
      dot_product += (double)( (*dx) * (*dy));
      dx += incx;
      dy += incy;
      }
  }
  return(dot_product);
}
typedef struct{
  long n;
  float *dx;
  long incx;
  float *dy;
  long incy;
  double result;
}svd_fdot_struct;
void svd_fdot_run(svd_fdot_struct *st){
  long i;
  if (st->incx == 1 && st->incy == 1) 
    for (i=0; i < st->n; i++) st->result += (double)((*st->dx++) * (*st->dy++));
  else {
    for (i=0; i < st->n; i++) {
      st->result += (double)( (*st->dx) * (*st->dy));
      st->dx += st->incx;
      st->dy += st->incy;
    }
  }
}

double svd_fdot_t(long n, float *dx, long incx, float *dy, long incy,int nthreads) {
  //threaded version of fdot.  Actually slower, so don't use...
  //Still slower, even if n is 1000,000,000, even using 8 threads.  Actually 2 threads can be faster, if use n>=1e8
  int i;
  double dot_product;
  pthread_t *threadid;
  int nth=nthreads;
  svd_fdot_struct *st;
  if (n <= 0 || incx == 0 || incy == 0) return(0.0);
  if (incx < 0) dx += (-n+1) * incx;//move pointer to last element.
  if (incy < 0) dy += (-n+1) * incy;
  
  threadid=malloc(sizeof(pthread_t)*nthreads);
  st=malloc(sizeof(svd_fdot_struct)*nthreads);
  st[0].n=n/nthreads;
  st[0].dx=dx;
  st[0].dy=dy;
  st[0].incx=incx;
  st[0].incy=incy;
  st[0].result=0.;
  dx+=incx*st[0].n;
  dy+=incy*st[0].n;
  nth--;
  n-=st[0].n;
  for(i=1; i<nthreads; i++){
    st[i].n=n/nth;
    nth--;
    n-=st[i].n;
    st[i].dx=dx;
    dx+=incx*st[i].n;
    st[i].dy=dy;
    dy+=incy*st[i].n;
    st[i].incx=incx;
    st[i].incy=incy;
    st[i].result=0.;
    pthread_create(&threadid[i],NULL,(void*)svd_fdot_run,&st[i]);
  }
  svd_fdot_run(&st[0]);
  for(i=1; i<nthreads; i++){
    pthread_join(threadid[i],NULL);
    st[0].result+=st[i].result;
  }
  free(threadid);
  dot_product=st[0].result;
  free(st);
  return dot_product;
}

/************************************************************** 
 * Constant times a vector plus a vector     		      *
 * Based on Fortran-77 routine from Linpack by J. Dongarra    *
 **************************************************************/ 
void svd_daxpy (long n, double da, double *dx, long incx, double *dy, long incy) {
  long i;
  
  if (n <= 0 || incx == 0 || incy == 0 || da == 0.0) return;
  if (incx == 1 && incy == 1) 
    for (i=0; i < n; i++) {
      *dy += da * (*dx++);
      dy++;
    }
  else {
    if (incx < 0) dx += (-n+1) * incx;
    if (incy < 0) dy += (-n+1) * incy;
    for (i=0; i < n; i++) {
      *dy += da * (*dx);
      dx += incx;
      dy += incy;
    }
  }
  return;
}
void svd_daxpydf (long n, double da, double *dx, long incx, float *dy, long incy) {
  long i;
  
  if (n <= 0 || incx == 0 || incy == 0 || da == 0.0) return;
  if (incx == 1 && incy == 1) 
    for (i=0; i < n; i++) {
      *dy += (float)(da * (*dx++));
      dy++;
    }
  else {
    if (incx < 0) dx += (-n+1) * incx;
    if (incy < 0) dy += (-n+1) * incy;
    for (i=0; i < n; i++) {
      *dy +=(float)( da * (*dx));
      dx += incx;
      dy += incy;
    }
  }
  return;
}
void svd_daxpyff (long n, double da, float *dx, long incx, float *dy, long incy) {
  long i;
  
  if (n <= 0 || incx == 0 || incy == 0 || da == 0.0) return;
  if (incx == 1 && incy == 1) 
    for (i=0; i < n; i++) {
      *dy += (float)(da * (*dx++));
      dy++;
    }
  else {
    if (incx < 0) dx += (-n+1) * incx;
    if (incy < 0) dy += (-n+1) * incy;
    for (i=0; i < n; i++) {
      *dy +=(float)( da * (*dx));
      dx += incx;
      dy += incy;
    }
  }
  return;
}
void svd_daxpyfd (long n, double da, float *dx, long incx, double *dy, long incy) {
  long i;
  
  if (n <= 0 || incx == 0 || incy == 0 || da == 0.0) return;
  if (incx == 1 && incy == 1) 
    for (i=0; i < n; i++) {
      *dy += da * (*dx++);
      dy++;
    }
  else {
    if (incx < 0) dx += (-n+1) * incx;
    if (incy < 0) dy += (-n+1) * incy;
    for (i=0; i < n; i++) {
      *dy += da * (*dx);
      dx += incx;
      dy += incy;
    }
  }
  return;
}


/********************************************************************* 
 * Function sorts array1 and array2 into increasing order for array1 *
 *********************************************************************/
void svd_dsort2(long igap, long n, double *array1, double *array2) {
  double temp;
  long i, j, index;
  
  if (!igap) return;
  else {
    for (i = igap; i < n; i++) {
      j = i - igap;
      index = i;
      while (j >= 0 && array1[j] > array1[index]) {
        temp = array1[j];
        array1[j] = array1[index];
        array1[index] = temp;
        temp = array2[j];
        array2[j] = array2[index];
        array2[index] = temp;
        j -= igap;
        index = j + igap;
      }
    } 
  }
  svd_dsort2(igap/2,n,array1,array2);
}

/************************************************************** 
 * Function interchanges two vectors		     	      *
 * Based on Fortran-77 routine from Linpack by J. Dongarra    *
 **************************************************************/ 
void svd_dswap(long n, double *dx, long incx, double *dy, long incy) {
  long i;
  double dtemp;
  
  if (n <= 0 || incx == 0 || incy == 0) return;
  if (incx == 1 && incy == 1) {
    for (i=0; i < n; i++) {
      dtemp = *dy;
      *dy++ = *dx;
      *dx++ = dtemp;
    }	
  }
  else {
    if (incx < 0) dx += (-n+1) * incx;
    if (incy < 0) dy += (-n+1) * incy;
    for (i=0; i < n; i++) {
      dtemp = *dy;
      *dy = *dx;
      *dx = dtemp;
      dx += incx;
      dy += incy;
    }
  }
}

/***************************************************************** 
 * Function finds the index of element having max. absolute value*
 * based on FORTRAN 77 routine from Linpack by J. Dongarra       *
 *****************************************************************/ 
long svd_idamax(long n, double *dx, long incx) {
  long ix,i,imax;
  double dtemp, dmax;
  
  if (n < 1) return(-1);
  if (n == 1) return(0);
  if (incx == 0) return(-1);
  
  if (incx < 0) ix = (-n+1) * incx;
  else ix = 0;
  imax = ix;
  dx += ix;
  dmax = fabs(*dx);
  for (i=1; i < n; i++) {
    ix += incx;
    dx += incx;
    dtemp = fabs(*dx);
    if (dtemp > dmax) {
      dmax = dtemp;
      imax = ix;
    }
  }
  return(imax);
}

/**************************************************************
 * multiplication of matrix B by vector x, where B = A'A,     *
 * and A is nrow by ncol (nrow >> ncol). Hence, B is of order *
 * n = ncol (y stores product vector).		              *
 **************************************************************/
void svd_opb(SMat A, double *x, double *y, double *temp) {
  unsigned int i, j, end;
  unsigned int *pointr = A->pointr, *rowind = A->rowind;
  double *value = A->value;
  float *valuef=A->valuef;
  unsigned int n = A->cols;

  SVDCount[SVD_MXV] += 2;
  memset(y, 0, n * sizeof(double));
  memset(temp,0,A->rows*sizeof(double));
  //for (i = 0; i < A->rows; i++) temp[i] = 0.0;
  
  if(value!=NULL){//double version

    for (i = 0; i < A->cols; i++) {
      end = pointr[i+1];
      for (j = pointr[i]; j < end; j++) 
	temp[rowind[j]] += value[j] * (*x); 
      x++;
    }
    for (i = 0; i < A->cols; i++) {
      end = pointr[i+1];
      for (j = pointr[i]; j < end; j++) 
	*y += value[j] * temp[rowind[j]];
      y++;
    }
  }else{//float version
    for (i = 0; i < A->cols; i++) {
      end = pointr[i+1];
      for (j = pointr[i]; j < end; j++) 
	temp[rowind[j]] += (double)valuef[j] * (*x); 
      x++;
    }
    for (i = 0; i < A->cols; i++) {
      end = pointr[i+1];
      for (j = pointr[i]; j < end; j++) 
	*y += (double)valuef[j] * temp[rowind[j]];
      y++;
    }
  }
  return;
}
void svd_opbf(SMat A, float *x, double *y, double *temp) {
  unsigned int i, j, end;
  unsigned int *pointr = A->pointr, *rowind = A->rowind;
  double *value = A->value;
  float *valuef=A->valuef;
  unsigned int n = A->cols;

  SVDCount[SVD_MXV] += 2;
  memset(y, 0, n * sizeof(double));
  memset(temp,0,A->rows*sizeof(double));
  //for (i = 0; i < A->rows; i++) temp[i] = 0.0;

  if(value!=NULL){//double version
    for (i = 0; i < A->cols; i++) {
      end = pointr[i+1];
      for (j = pointr[i]; j < end; j++) 
	temp[rowind[j]] += value[j] * (*x); 
      x++;
    }
    for (i = 0; i < A->cols; i++) {
      end = pointr[i+1];
      for (j = pointr[i]; j < end; j++) 
	*y += value[j] * temp[rowind[j]];
      y++;
    }
  }else{//float version
    for (i = 0; i < A->cols; i++) {
      end = pointr[i+1];
      for (j = pointr[i]; j < end; j++) 
	temp[rowind[j]] += (double)(valuef[j] * (*x)); 
      x++;
    }
    for (i = 0; i < A->cols; i++) {
      end = pointr[i+1];
      for (j = pointr[i]; j < end; j++) 
	*y += (double)valuef[j] * temp[rowind[j]];
      y++;
    }
  }
  return;
}



typedef struct{
  double *value;
  float *valuef;
  double *temp;
  unsigned int *rowind;
  unsigned int *pointr;
  unsigned int colStart;
  unsigned int colEnd;
  double *y;
}svd_opb_struct;
void svd_opb_worker_d(svd_opb_struct *sos){
  unsigned int i,j,end;
  double *y;
  double *value=sos->value;
  double *temp=sos->temp;
  unsigned int *rowind=sos->rowind, *pointr=sos->pointr;

  y=&(sos->y[sos->colStart]);
  for(i=sos->colStart; i<sos->colEnd; i++){
    end=pointr[i+1];
    for(j=pointr[i]; j<end; j++){
      *y+=value[j]*temp[rowind[j]];
    }
    y++;
  } 
}
void svd_opb_worker_f(svd_opb_struct *sos){
  unsigned int i,j,end;
  double *y;
  float *valuef=sos->valuef;
  double *temp=sos->temp;
  unsigned int *rowind=sos->rowind, *pointr=sos->pointr;

  y=&(sos->y[sos->colStart]);
  for(i=sos->colStart; i<sos->colEnd; i++){
    end=pointr[i+1];
    for(j=pointr[i]; j<end; j++){
      *y+=(double)valuef[j]*temp[rowind[j]];
    }
    y++;
  } 
}

void svd_opb_t(SMat A, double *x, double *y, double *temp,int nthreads) {
  unsigned int i, j, end;
  unsigned int *pointr = A->pointr, *rowind = A->rowind;
  double *value = A->value;
  float *valuef=A->valuef;
  unsigned int n = A->cols;
  unsigned int nleft;
  pthread_t *thread;
  svd_opb_struct *threadInfo;
  SVDCount[SVD_MXV] += 2;
  memset(y, 0, n * sizeof(double));
  memset(temp,0,A->rows*sizeof(double));
  //for (i = 0; i < A->rows; i++) temp[i] = 0.0;
  
  threadInfo=malloc(sizeof(svd_opb_struct)*nthreads);
  thread=malloc(sizeof(pthread_t)*nthreads);
  for(i=0; i<nthreads; i++){
    threadInfo[i].y=y;
    threadInfo[i].temp=temp;
    threadInfo[i].pointr=pointr;
    threadInfo[i].value=value;
    threadInfo[i].valuef=valuef;
    threadInfo[i].rowind=rowind;
  }

  if(value!=NULL){//double version

    for (i = 0; i < A->cols; i++) {
      end = pointr[i+1];
      for (j = pointr[i]; j < end; j++) 
	temp[rowind[j]] += value[j] * (*x); 
      x++;
    }
    /*for (i = 0; i < A->cols; i++) {
      end = pointr[i+1];
      for (j = pointr[i]; j < end; j++) 
	*y += value[j] * temp[rowind[j]];
      y++;
      }*/
    threadInfo[0].colStart=0;
    threadInfo[0].colEnd=A->cols/nthreads;
    nleft=A->cols-threadInfo[0].colEnd;
    for(i=1; i<nthreads; i++){
      threadInfo[i].colStart=threadInfo[i-1].colEnd;
      threadInfo[i].colEnd=threadInfo[i].colStart+nleft/(nthreads-i);
      nleft-=threadInfo[i].colEnd-threadInfo[i].colStart;
      pthread_create(&thread[i],NULL,(void*)svd_opb_worker_d,&threadInfo[i]);
    }
    svd_opb_worker_d(&threadInfo[0]);
    for(i=1; i<nthreads; i++)
      pthread_join(thread[i],NULL);
  }else{//float version
    for (i = 0; i < A->cols; i++) {
      end = pointr[i+1];
      for (j = pointr[i]; j < end; j++) 
	temp[rowind[j]] += (double)valuef[j] * (*x); 
      x++;
    }
    /*for (i = 0; i < A->cols; i++) {
      end = pointr[i+1];
      for (j = pointr[i]; j < end; j++) 
	*y += (double)valuef[j] * temp[rowind[j]];
      y++;
    }
    */
    threadInfo[0].colStart=0;
    threadInfo[0].colEnd=A->cols/nthreads;
    nleft=A->cols-threadInfo[0].colEnd;
    for(i=1; i<nthreads; i++){
      threadInfo[i].colStart=threadInfo[i-1].colEnd;
      threadInfo[i].colEnd=threadInfo[i].colStart+nleft/(nthreads-i);
      nleft-=threadInfo[i].colEnd-threadInfo[i].colStart;
      pthread_create(&thread[i],NULL,(void*)svd_opb_worker_f,&threadInfo[i]);
    }
    svd_opb_worker_f(&threadInfo[0]);
    for(i=1; i<nthreads; i++)
      pthread_join(thread[i],NULL);
  }
  return;
}
void svd_opbf_t(SMat A, float *x, double *y, double *temp,int nthreads) {
  unsigned int i, j, end;
  unsigned int *pointr = A->pointr, *rowind = A->rowind;
  double *value = A->value;
  float *valuef=A->valuef;
  unsigned int n = A->cols;
  unsigned int nleft;
  pthread_t *thread;
  svd_opb_struct *threadInfo;

  SVDCount[SVD_MXV] += 2;
  memset(y, 0, n * sizeof(double));
  memset(temp,0,A->rows*sizeof(double));
  //for (i = 0; i < A->rows; i++) temp[i] = 0.0;

  threadInfo=malloc(sizeof(svd_opb_struct)*nthreads);
  thread=malloc(sizeof(pthread_t)*nthreads);
  for(i=0; i<nthreads; i++){
    threadInfo[i].y=y;
    threadInfo[i].temp=temp;
    threadInfo[i].pointr=pointr;
    threadInfo[i].value=value;
    threadInfo[i].valuef=valuef;
    threadInfo[i].rowind=rowind;
  }


  if(value!=NULL){//double version
    for (i = 0; i < A->cols; i++) {
      end = pointr[i+1];
      for (j = pointr[i]; j < end; j++) 
	temp[rowind[j]] += value[j] * (*x); 
      x++;
    }
    /*for (i = 0; i < A->cols; i++) {
      end = pointr[i+1];
      for (j = pointr[i]; j < end; j++) 
	*y += value[j] * temp[rowind[j]];
      y++;
    }
    */
    threadInfo[0].colStart=0;
    threadInfo[0].colEnd=A->cols/nthreads;
    nleft=A->cols-threadInfo[0].colEnd;
    for(i=1; i<nthreads; i++){
      threadInfo[i].colStart=threadInfo[i-1].colEnd;
      threadInfo[i].colEnd=threadInfo[i].colStart+nleft/(nthreads-i);
      nleft-=threadInfo[i].colEnd-threadInfo[i].colStart;
      pthread_create(&thread[i],NULL,(void*)svd_opb_worker_d,&threadInfo[i]);
    }
    svd_opb_worker_d(&threadInfo[0]);
    for(i=1; i<nthreads; i++)
      pthread_join(thread[i],NULL);

  }else{//float version
    for (i = 0; i < A->cols; i++) {
      end = pointr[i+1];
      for (j = pointr[i]; j < end; j++) 
	temp[rowind[j]] += (double)(valuef[j] * (*x)); 
      x++;
    }
    /*for (i = 0; i < A->cols; i++) {
      end = pointr[i+1];
      for (j = pointr[i]; j < end; j++) 
	*y += (double)valuef[j] * temp[rowind[j]];
      y++;
      }*/
    threadInfo[0].colStart=0;
    threadInfo[0].colEnd=A->cols/nthreads;
    nleft=A->cols-threadInfo[0].colEnd;
    for(i=1; i<nthreads; i++){
      threadInfo[i].colStart=threadInfo[i-1].colEnd;
      threadInfo[i].colEnd=threadInfo[i].colStart+nleft/(nthreads-i);
      nleft-=threadInfo[i].colEnd-threadInfo[i].colStart;
      pthread_create(&thread[i],NULL,(void*)svd_opb_worker_f,&threadInfo[i]);
    }
    svd_opb_worker_f(&threadInfo[0]);
    for(i=1; i<nthreads; i++)
      pthread_join(thread[i],NULL);
  }
  return;
}

/***********************************************************
 * multiplication of matrix A by vector x, where A is 	   *
 * nrow by ncol (nrow >> ncol).  y stores product vector.  *
 * Ax=y                                                    *
 ***********************************************************/
void svd_opa(SMat A, double *x, double *y) {
  unsigned int end, i, j;
  unsigned int *pointr = A->pointr, *rowind = A->rowind;
  double *value = A->value;
  float *valuef=A->valuef;
   
  SVDCount[SVD_MXV]++;
  memset(y, 0, A->rows * sizeof(double));
  if(value!=NULL){//double version
    for (i = 0; i < A->cols; i++) {
      end = pointr[i+1];
      for (j = pointr[i]; j < end; j++)
	y[rowind[j]] += value[j] * x[i]; 
    }
  }else{//float version
    for (i = 0; i < A->cols; i++) {
      end = pointr[i+1];
      for (j = pointr[i]; j < end; j++)
	y[rowind[j]] += (double)valuef[j] * x[i]; 
    }
  }

  return;
}
void svd_opaff(SMat A, float *x, float *y) {
  unsigned int end, i, j;
  unsigned int *pointr = A->pointr, *rowind = A->rowind;
  double *value = A->value;
  float *valuef=A->valuef;
  SVDCount[SVD_MXV]++;
  memset(y, 0, A->rows * sizeof(float));
  if(value!=NULL){//double version
    for (i = 0; i < A->cols; i++) {
      end = pointr[i+1];
      for (j = pointr[i]; j < end; j++)
	y[rowind[j]] += (float)(value[j] * x[i]); 
    }
  }else{//float version
    for (i = 0; i < A->cols; i++) {
      end = pointr[i+1];
      for (j = pointr[i]; j < end; j++)
	y[rowind[j]] += (float)(valuef[j] * x[i]); 
    }
  }
  return;
}


/***********************************************************************
 *                                                                     *
 *				random()                               *
 *                        (double precision)                           *
 ***********************************************************************/
/***********************************************************************

   Description
   -----------

   This is a translation of a Fortran-77 uniform random number
   generator.  The code is based  on  theory and suggestions  given in
   D. E. Knuth (1969),  vol  2.  The argument to the function should 
   be initialized to an arbitrary integer prior to the first call to 
   random.  The calling program should  not  alter  the value of the
   argument between subsequent calls to random.  Random returns values
   within the interval (0,1).


   Arguments 
   ---------

   (input)
   iy	   an integer seed whose value must not be altered by the caller
	   between subsequent calls

   (output)
   random  a double precision random number between (0,1)

 ***********************************************************************/
double svd_random2(long *iy) {
  //all these are only declared once, and never change.
   static long m2 = 0;
   static long ia, ic, mic;
   static double halfm, s;

   /* If first entry, compute (max int) / 2 */
   if (!m2) {
      m2 = 1 << (8 * (int)sizeof(int) - 2); 
      halfm = m2;

      /* compute multiplier and increment for linear congruential 
       * method */
      ia = 8 * (long)(halfm * atan(1.0) / 8.0) + 5;
      ic = 2 * (long)(halfm * (0.5 - sqrt(3.0)/6.0)) + 1;
      mic = (m2-ic) + m2;

      /* s is the scale factor for converting to floating point */
      s = 0.5 / halfm;
   }

   /* compute next random number */
   *iy = *iy * ia;

   /* for computers which do not allow integer overflow on addition */
   if (*iy > mic) *iy = (*iy - m2) - m2;

   *iy = *iy + ic;

   /* for computers whose word length for addition is greater than
    * for multiplication */
   if (*iy / 2 > m2) *iy = (*iy - m2) - m2;
  
   /* for computers whose integer overflow affects the sign bit */
   if (*iy < 0) *iy = (*iy + m2) + m2;

   return((double)(*iy) * s);
}

/************************************************************** 
 *							      *
 * Function finds sqrt(a^2 + b^2) without overflow or         *
 * destructive underflow.				      *
 *							      *
 **************************************************************/ 
/************************************************************** 

   Funtions used
   -------------

   UTILITY	dmax, dmin

 **************************************************************/ 
double svd_pythag(double a, double b) {
   double p, r, s, t, u, temp;

   p = svd_dmax(fabs(a), fabs(b));
   if (p != 0.0) {
      temp = svd_dmin(fabs(a), fabs(b)) / p;
      r = temp * temp; 
      t = 4.0 + r;
      while (t != 4.0) {
	 s = r / t;
	 u = 1.0 + 2.0 * s;
	 p *= u;
	 temp = s / u;
	 r *= temp * temp;
	 t = 4.0 + r;
      }
   }
   return(p);
}

