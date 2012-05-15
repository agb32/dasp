#ifndef SVDLIB_H
#define SVDLIB_H

#ifndef FALSE
#  define FALSE 0
#endif
#ifndef TRUE
#  define TRUE  1
#endif

/******************************** Structures *********************************/
typedef struct smat *SMat;
typedef struct dmat *DMat;
typedef struct svdrec *SVDRec;

/* Harwell-Boeing sparse matrix. */
struct smat {
  unsigned int rows;
  unsigned int cols;
  unsigned int vals;     /* Total non-zero entries. */
  unsigned int *pointr;  /* For each col (plus 1), index of first non-zero entry. */
  unsigned int *rowind;  /* For each nz entry, the row index. */
  double *value; /* For each nz entry, the value. */
  float *valuef;// For each nz entry the value is value is NULL.
};

/* Row-major dense matrix.  Rows are consecutive vectors. */
struct dmat {
  long rows;
  long cols;
  double **value; /* Accessed by [row][col]. Free value[0] and value to free.*/
  float **valuef;
};

struct svdrec {
  int d;      /* Dimensionality (rank) */
  DMat Ut;    /* Transpose of left singular vectors. (d by m)
                 The vectors are the rows of Ut. */
  float *S;  /* Array of singular values. (length d) */
  DMat Vt;    /* Transpose of right singular vectors. (d by n)
                 The vectors are the rows of Vt. */
  int err;
};

typedef struct{
  int iterations;
  int maxll;
  double **LanStore;
  float **LanStoreFloat;
}LanStoreStruct;


/******************************** Variables **********************************/

/* Version info */
extern char *SVDVersion;

/* How verbose is the package: 0, 1 (default), 2 */
extern long SVDVerbosity;

/* Counter(s) used to track how much work is done in computing the SVD. */
enum svdCounters {SVD_MXV, SVD_COUNTERS};
extern long SVDCount[SVD_COUNTERS];
extern void svdResetCounters(void);

enum svdFileFormats {SVD_F_STH, SVD_F_ST, SVD_F_SB, SVD_F_DT, SVD_F_DB};
/*
File formats:
SVD_F_STH: sparse text, SVDPACK-style
SVD_F_ST:  sparse text, SVDLIB-style
SVD_F_DT:  dense text
SVD_F_SB:  sparse binary
SVD_F_DB:  dense binary
*/

/* True if a file format is sparse: */
#define SVD_IS_SPARSE(format) ((format >= SVD_F_STH) && (format <= SVD_F_SB))


/******************************** Functions **********************************/

/* Creates an empty dense matrix. */
extern DMat svdNewDMat(int rows, int cols);
extern DMat svdNewDMatWithData(int rows, int cols,double *data);
extern DMat svdNewDMatWithNoData(int rows, int cols);
extern DMat svdNewDMatWithDataFloat(int rows, int cols,float *data);
/* Frees a dense matrix. */
extern void svdFreeDMat(DMat D);

/* Creates an empty sparse matrix. */
SMat svdNewSMat(int rows, int cols, unsigned int vals);
SMat svdNewSMatFloat(int rows, int cols, unsigned int vals);
/* Frees a sparse matrix. */
void svdFreeSMat(SMat S);

/* Creates an empty SVD record. */
SVDRec svdNewSVDRec(void);
/* Frees an svd rec and all its contents. */
void svdFreeSVDRec(SVDRec R);

/* Converts a sparse matrix to a dense one (without affecting former) */
DMat svdConvertStoD(SMat S);
/* Converts a dense matrix to a sparse one (without affecting former) */
SMat svdConvertDtoS(DMat D);

/* Transposes a dense matrix (returning a new one) */
DMat svdTransposeD(DMat D);
/* Transposes a sparse matrix (returning a new one) */
SMat svdTransposeS(SMat S);

/* Writes an array to a file. */
extern void svdWriteDenseArray(double *a, int n, char *filename, char binary);
extern void svdWriteDenseArrayFloat(float *a, int n, char *filename, char binary);
/* Reads an array from a file, storing its size in *np. */
extern double *svdLoadDenseArray(char *filename, int *np, char binary);

/* Loads a matrix file (in various formats) into a sparse matrix. */
extern SMat svdLoadSparseMatrix(char *filename, int format);
/* Loads a matrix file (in various formats) into a dense matrix. */
extern DMat svdLoadDenseMatrix(char *filename, int format);

/* Writes a dense matrix to a file in a given format. */
extern void svdWriteDenseMatrix(DMat A, char *filename, int format);
/* Writes a sparse matrix to a file in a given format. */
extern void svdWriteSparseMatrix(SMat A, char *filename, int format);


/* Performs the las2 SVD algorithm and returns the resulting Ut, S, and Vt. */
extern SVDRec svdLAS2(SMat A, long dimensions, long iterations, double end[2], 
                      double kappa,SVDRec R,long ritvecSize,int userNeig, ArrUnion *genInv, float minEig, float fracEig,int useStoreFloat, int useSFloat,double minGIVal,int transposeGI,int prepareForGenInv,int considerSwap,int nthreads);
/* Chooses default parameter values.  Set dimensions to 0 for all dimensions: */
extern SVDRec svdLAS2A(SMat A, long dimensions);

#endif /* SVDLIB_H */
