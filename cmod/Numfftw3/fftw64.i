typedef double fftw_complex[2]; 
typedef struct fftw_plan_s *fftw_plan;
typedef struct fftw_iodim_do_not_use_me fftw_iodim;
typedef enum fftw_r2r_kind_do_not_use_me fftw_r2r_kind;

/*****************************************************************************************************************/

/* Memory allocation functions */
void *fftw_malloc(size_t n);
void fftw_free(void *p);


/* Basic Initialisation functions */
void fftw_execute(const fftw_plan p);
void fftw_destroy_plan(fftw_plan p);
void fftw_cleanup(void);
void fftw_flops(const fftw_plan p, double *add, double *mul, double *fma);
void fftw_print_plan(const fftw_plan p);
void fftw_fprint_plan(const fftw_plan p, FILE *output_file);

/*****************************************************************************************************************/

/* Basic Complex DFT functions */
fftw_plan fftw_plan_dft_1d(int n, fftw_complex *in, fftw_complex *out, int sign,unsigned flags);
fftw_plan fftw_plan_dft_2d(int nx, int ny, fftw_complex *in, fftw_complex *out, int sign, unsigned flags);
fftw_plan fftw_plan_dft_3d(int nx, int ny, int nz,fftw_complex *in, fftw_complex *out, int sign, unsigned flags);
fftw_plan fftw_plan_dft(int rank, const int *n, fftw_complex *in, fftw_complex *out, int sign, unsigned flags);

/* Advanced Complex DFT */
fftw_plan fftw_plan_many_dft(int rank, const int *n, int howmany,fftw_complex *in, const int *inembed, int istride, int idist, fftw_complex *out, const int *onembed, int ostride, int odist, int sign, unsigned flags);

/* Guru Complex DFT functions */
fftw_plan fftw_plan_guru_dft(int rank, const fftw_iodim *dims, int
howmany_rank, const fftw_iodim *howmany_dims, fftw_complex *in,
fftw_complex *out, int sign, unsigned flags);
void fftw_execute_dft(const fftw_plan p, fftw_complex *in, fftw_complex *out);


/* Guru Split Complex DFT functions */
fftw_plan fftw_plan_guru_split_dft(int rank, const fftw_iodim *dims,
int howmany_rank, const fftw_iodim *howmany_dims, double *ri, double
*ii, double *ro, double *io, unsigned flags);
void fftw_execute_split_dft(const fftw_plan p, double *ri, double *ii, double *ro, double *io);


/*****************************************************************************************************************/

/* Basic Real DFT functions */
fftw_plan fftw_plan_dft_r2c_1d(int n,double *in,fftw_complex *out,unsigned flags);
fftw_plan fftw_plan_dft_r2c_2d(int nx, int ny, double *in, fftw_complex *out, unsigned flags);
fftw_plan fftw_plan_dft_r2c_3d(int nx, int ny, int nz, double *in, fftw_complex *out, unsigned flags);
fftw_plan fftw_plan_dft_r2c(int rank, const int *n, double *in, fftw_complex *out, unsigned flags);

/* Advanced Real DFT functions */
fftw_plan fftw_plan_many_dft_r2c(int rank, const int *n, int howmany,
double *in, const int *inembed, int istride, int idist, fftw_complex
*out, const int *onembed, int ostride, int odist, unsigned flags);

/* Guru Real DFT functions */
fftw_plan fftw_plan_guru_dft_r2c(int rank, const fftw_iodim *dims, int
howmany_rank, const fftw_iodim *howmany_dims, double *in, fftw_complex
*out, unsigned flags);
void fftw_execute_dft_r2c(const fftw_plan p, double *in, fftw_complex *out);

/* Guru Split Real DFT functions */
fftw_plan fftw_plan_guru_split_dft_r2c(int rank, const fftw_iodim
*dims, int howmany_rank, const fftw_iodim *howmany_dims, double *in,
double *ro, double *io,  unsigned flags);
void fftw_execute_split_dft_r2c(const fftw_plan p, double *in, double *ro, double *io);

/*****************************************************************************************************************/

/* Basic Inverse Real DFT functions */
fftw_plan fftw_plan_dft_c2r_1d(int n,fftw_complex *in,double *out,unsigned flags);
fftw_plan fftw_plan_dft_c2r_2d(int nx, int ny, fftw_complex *in, double *out, unsigned flags);
fftw_plan fftw_plan_dft_c2r_3d(int nx, int ny, int nz, fftw_complex *in, double *out, unsigned flags);
fftw_plan fftw_plan_dft_c2r(int rank, const int *n, fftw_complex *in, double *out, unsigned flags);

/* Advanced Inverse Real DFT functions */
fftw_plan fftw_plan_many_dft_c2r(int rank, const int *n, int howmany,
fftw_complex *in, const int *inembed, int istride, int idist, double
*out, const int *onembed, int ostride, int odist, unsigned flags);

/* Guru Inverse Real DFT functions */
fftw_plan fftw_plan_guru_dft_c2r(int rank, const fftw_iodim *dims, int
howmany_rank, const fftw_iodim *howmany_dims, fftw_complex *in, double
*out, unsigned flags);
void fftw_execute_dft_c2r(const fftw_plan p, fftw_complex *in, double *out);


/* Guru Split Inverse Real DFT functions */
fftw_plan fftw_plan_guru_split_dft_c2r(int rank, const fftw_iodim
*dims,  int howmany_rank,  const fftw_iodim *howmany_dims,  double
*ri, double *ii, double *out,  unsigned flags);
void fftw_execute_split_dft_c2r(const fftw_plan p, double *ri, double
*ii, double *out);

/*****************************************************************************************************************/

/* Basic Real to Real DFT functions */
fftw_plan fftw_plan_r2r_1d(int n, double *in, double *out,fftw_r2r_kind kind, unsigned flags);
fftw_plan fftw_plan_r2r_2d(int nx, int ny, double *in, double *out,fftw_r2r_kind kindx, fftw_r2r_kind kindy,unsigned flags);
fftw_plan fftw_plan_r2r_3d(int nx, int ny, int nz,double *in, double *out, fftw_r2r_kind kindx,fftw_r2r_kind kindy, fftw_r2r_kind kindz,unsigned flags);
fftw_plan fftw_plan_r2r(int rank, const int *n, double *in, double *out, const fftw_r2r_kind *kind, unsigned flags);

/* Advanced Real to Real DFT functions */
fftw_plan fftw_plan_many_r2r(int rank, const int *n,  int howmany,
double *in, const int *inembed,  int istride, int idist,  double *out,
const int *onembed,  int ostride, int odist,  const fftw_r2r_kind
*kind, unsigned flags);

/* Guru Real to Real DFT functions */
fftw_plan fftw_plan_guru_r2r(int rank, const fftw_iodim *dims,  int
howmany_rank,  const fftw_iodim *howmany_dims,  double *in, double
*out,  const fftw_r2r_kind *kind, unsigned flags);
void fftw_execute_r2r(const fftw_plan p, double *in, double *out);

/*****************************************************************************************************************/

/* Wisdom Export functions */
void fftw_export_wisdom_to_file(FILE *output_file);
char *fftw_export_wisdom_to_string(void);
void fftw_export_wisdom(void (*write_char)(char c, void *), void *data);

/* Wisdom Import functions */
int fftw_import_system_wisdom(void);
int fftw_import_wisdom_from_file(FILE *input_file);
int fftw_import_wisdom_from_string(const char *input_string);
int fftw_import_wisdom(int (*read_char)(void *), void *data);
void fftw_forget_wisdom(void);

/*****************************************************************************************************************/

/* Parallel FFTW functions */
/*int fftw_init_threads(void);
void fftw_plan_with_nthreads(int nthreads);
void fftw_cleanup_threads(void);*/

/*****************************************************************************************************************/

