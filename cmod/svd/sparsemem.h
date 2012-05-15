#ifndef SPARSEMEM_H
#define SPARSEMEM_H


typedef struct{
  char typ;//should be 'd'
  double *data;
  unsigned int *rowind;
  unsigned int *indptr;
  unsigned int ndata;
  int rows;
  int cols;
  double min;//the current minimum value
  unsigned int rowmin;//position in rowind of this min
  unsigned int indmin;//position in indptr of this min.
  unsigned int cnt;//the current number of entries stored.
  int alloced;//whether we have allocated the memory, or was it passed in?
}dSpMem;
typedef struct{
  char typ;//should be 'f'
  float *data;
  unsigned int *rowind;
  unsigned int *indptr;
  unsigned int ndata;
  int rows;
  int cols;
  float min;//the current minimum value
  unsigned int rowmin;//position in rowind of this min
  unsigned int indmin;//position in indptr of this min.
  unsigned int cnt;//the current number of entries stored.
  int alloced;//whether we have allocated the memory, or was it passed in?
}fSpMem;

typedef struct{
  char typ;//can be f,d,F,D for float, double, sparse float, sparse double
  union {
    float *fdata;
    double *ddata;
    fSpMem *fSp;
    dSpMem *dSp;
  } data;
}ArrUnion;

unsigned int smGetIndxForRowFloat(fSpMem *sp,unsigned int row,unsigned int col, unsigned int suggestion);
int smReplaceDataFloat(fSpMem *sp,unsigned int rowind,unsigned int col, float data);
float smGetFloat(fSpMem *sp,unsigned int row, unsigned int col);
int smInsertFloat(fSpMem *sp,unsigned int row, unsigned int col, float data);
void smFreeSparseMemFloat(fSpMem *sp);
fSpMem *smNewSparseMemFloat(unsigned int ndata,unsigned int rows,unsigned int cols);
fSpMem *smNewFromExistingFloat(unsigned int ndata,unsigned int rows,unsigned int cols, float* data,unsigned int *indptr,unsigned int *rowind);

unsigned int smGetIndxForRow(dSpMem *sp,unsigned int row,unsigned int col,unsigned int suggestion);
int smReplaceData(dSpMem *sp,unsigned int rowind,unsigned int col, double data);
double smGet(dSpMem *sp,unsigned int row,unsigned int col);
int smInsert(dSpMem *sp,unsigned int row,unsigned int col, double data);
void smFreeSparseMem(dSpMem *sp);
dSpMem *smNewSparseMem(unsigned int ndata,unsigned int rows,unsigned int cols);
dSpMem *smNewFromExisting(unsigned int ndata,unsigned int rows,unsigned int cols, double* data,unsigned int *indptr,unsigned int *rowind);

int sparseDotSparse(float *dataa,unsigned int *colind,unsigned int *aindptr,int ninda,float *datab,unsigned int *rowind,unsigned int *bindptr,int nindb,float **resdataptr,unsigned int **rescolindptr,unsigned int* resindptr,unsigned int *resn,int pr);
int denseDotCSC(float *a,int arows, int acols, float *datab,unsigned int *browind, unsigned int *bindptr,int nindb,float **resdataptr,unsigned int **resrcindptr, unsigned int *resindptr,unsigned int *resn,float valmin,int resIsCsc,unsigned int resnmax,int pr);
int countInstances(float *a,long size,float *vals,long *cnt, int n);
int sparsifyCsr(float *a,int nrows,int ncols,int lda,int ldb,float *data,unsigned int *colind, unsigned int *indptr,unsigned int n,float val);
int densifyCsr(float *a, int nrows, int ncols, int lda, int ldb, float *data, unsigned int *colind,unsigned int *indptr,unsigned int n);
int densifyCsrL(float *a, int nrows, int ncols, int lda, int ldb, float *data, unsigned int *colind,unsigned long *indptr,unsigned int n);






#endif
