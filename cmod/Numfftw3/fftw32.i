typedef float fftwf_complex[2]; 
typedef struct fftwf_plan_s *fftwf_plan;
typedef struct fftw_iodim_do_not_use_me fftwf_iodim;
typedef enum fftw_r2r_kind_do_not_use_me fftwf_r2r_kind;

/*****************************************************************************************************************/

/* Memory allocation functions */
void *fftwf_malloc(size_t n);
void fftwf_free(void *p);


/* Basic Initialisation functions */
void fftwf_execute(const fftwf_plan p);
void fftwf_destroy_plan(fftwf_plan p);
void fftwf_cleanup(void);
void fftwf_flops(const fftwf_plan p, double *add, double *mul, double *fma);
void fftwf_print_plan(const fftwf_plan p);
void fftwf_fprint_plan(const fftwf_plan p, FILE *output_file);

/*****************************************************************************************************************/

/* Basic Complex DFT functions */
fftwf_plan fftwf_plan_dft_1d(int n, fftwf_complex *in, fftwf_complex *out, int sign,unsigned flags);
fftwf_plan fftwf_plan_dft_2d(int nx, int ny, fftwf_complex *in, fftwf_complex *out, int sign, unsigned flags);
fftwf_plan fftwf_plan_dft_3d(int nx, int ny, int nz,fftwf_complex *in, fftwf_complex *out, int sign, unsigned flags);
fftwf_plan fftwf_plan_dft(int rank, const int *n, fftwf_complex *in, fftwf_complex *out, int sign, unsigned flags);

/* Advanced Complex DFT
fftwf_plan fftwf_plan_many_dft(int rank, const int *n, int
howmany,fftwf_complex *in, const int *inembed, int istride, int idist,
fftwf_complex *out, const int *onembed, int ostride, int odist, int sign, unsigned flags);

/* Guru Complex DFT functions */
fftwf_plan fftwf_plan_guru_dft(int rank, const fftwf_iodim *dims, int
howmany_rank, const fftwf_iodim *howmany_dims, fftwf_complex *in,
fftwf_complex *out, int sign, unsigned flags);
void fftwf_execute_dft(const fftwf_plan p, fftwf_complex *in, fftwf_complex *out);

/* Guru Split Complex DFT functions */
fftwf_plan fftwf_plan_guru_split_dft(int rank, const fftwf_iodim *dims,
int howmany_rank, const fftwf_iodim *howmany_dims, float *ri, float
*ii, float *ro, float *io, unsigned flags);
void fftwf_execute_split_dft(const fftwf_plan p, float *ri, float *ii, float *ro, float *io);

/*****************************************************************************************************************/

/* Basic Real DFT functions */
fftwf_plan fftwf_plan_dft_r2c_1d(int n,float *in,fftwf_complex *out,unsigned flags);
fftwf_plan fftwf_plan_dft_r2c_2d(int nx, int ny, float *in, fftwf_complex *out, unsigned flags);
fftwf_plan fftwf_plan_dft_r2c_3d(int nx, int ny, int nz, float *in, fftwf_complex *out, unsigned flags);
fftwf_plan fftwf_plan_dft_r2c(int rank, const int *n, float *in, fftwf_complex *out, unsigned flags);

/* Advanced Real DFT functions */
fftwf_plan fftwf_plan_many_dft_r2c(int rank, const int *n, int howmany,
float *in, const int *inembed, int istride, int idist, fftwf_complex
*out, const int *onembed, int ostride, int odist, unsigned flags);

/* Guru Real DFT functions */
fftwf_plan fftwf_plan_guru_dft_r2c(int rank, const fftwf_iodim *dims, int
howmany_rank, const fftwf_iodim *howmany_dims, float *in, fftwf_complex
*out, unsigned flags);
void fftwf_execute_dft_r2c(const fftwf_plan p, float *in, fftwf_complex *out);

/* Guru Split Real DFT functions */
fftwf_plan fftwf_plan_guru_split_dft_r2c(int rank, const fftwf_iodim
*dims, int howmany_rank, const fftwf_iodim *howmany_dims, float *in,
float *ro, float *io,  unsigned flags);
void fftwf_execute_split_dft_r2c(const fftwf_plan p, float *in, float *ro, float *io);

/*****************************************************************************************************************/

/* Basic Inverse Real DFT functions */
fftwf_plan fftwf_plan_dft_c2r_1d(int n,fftwf_complex *in,float *out,unsigned flags);
fftwf_plan fftwf_plan_dft_c2r_2d(int nx, int ny, fftwf_complex *in, float *out, unsigned flags);
fftwf_plan fftwf_plan_dft_c2r_3d(int nx, int ny, int nz, fftwf_complex *in, float *out, unsigned flags);
fftwf_plan fftwf_plan_dft_c2r(int rank, const int *n, fftwf_complex *in, float *out, unsigned flags);

/* Advanced Inverse Real DFT functions */
fftwf_plan fftwf_plan_many_dft_c2r(int rank, const int *n, int howmany,
fftwf_complex *in, const int *inembed, int istride, int idist, float
*out, const int *onembed, int ostride, int odist, unsigned flags);

/* Guru Inverse Real DFT functions */
fftwf_plan fftwf_plan_guru_dft_c2r(int rank, const fftwf_iodim *dims, int
howmany_rank, const fftwf_iodim *howmany_dims, fftwf_complex *in, float
*out, unsigned flags);
void fftwf_execute_dft_c2r(const fftwf_plan p, fftwf_complex *in, float *out);

/* Guru Split Inverse Real DFT functions */
fftwf_plan fftwf_plan_guru_split_dft_c2r(int rank, const fftwf_iodim
*dims,  int howmany_rank,  const fftwf_iodim *howmany_dims,  float
*ri, float *ii, float *out,  unsigned flags);
void fftwf_execute_split_dft_c2r(const fftwf_plan p, float *ri, float
*ii, float *out);

/*****************************************************************************************************************/

/* Basic Real to Real DFT functions */
fftwf_plan fftwf_plan_r2r_1d(int n, float *in, float
*out,fftwf_r2r_kind kind, unsigned flags);
fftwf_plan fftwf_plan_r2r_2d(int nx, int ny, float *in, float
*out,fftwf_r2r_kind kindx, fftwf_r2r_kind kindy,unsigned flags);
fftwf_plan fftwf_plan_r2r_3d(int nx, int ny, int nz,float *in, float
*out, fftwf_r2r_kind kindx,fftwf_r2r_kind kindy, fftwf_r2r_kind
kindz,unsigned flags);
fftwf_plan fftwf_plan_r2r(int rank, const int *n, float *in, float
*out, const fftwf_r2r_kind *kind, unsigned flags);

/* Advanced Real to Real DFT functions */
fftwf_plan fftwf_plan_many_r2r(int rank, const int *n,  int howmany,
float *in, const int *inembed,  int istride, int idist,  float *out,
const int *onembed,  int ostride, int odist,  const fftwf_r2r_kind
*kind, unsigned flags);


/* Guru Real to Real DFT functions */
fftwf_plan fftwf_plan_guru_r2r(int rank, const fftwf_iodim *dims,  int
howmany_rank,  const fftwf_iodim *howmany_dims,  float *in, float
*out,  const fftwf_r2r_kind *kind, unsigned flags);
void fftwf_execute_r2r(const fftwf_plan p, float *in, float *out);

/*****************************************************************************************************************/

/* Wisdom Export functions */
void fftwf_export_wisdom_to_file(FILE *output_file);
char *fftwf_export_wisdom_to_string(void);
void fftwf_export_wisdom(void (*write_char)(char c, void *), void *data);

/* Wisdom Import functions */
int fftwf_import_system_wisdom(void);
int fftwf_import_wisdom_from_file(FILE *input_file);
int fftwf_import_wisdom_from_string(const char *input_string);
int fftwf_import_wisdom(int (*read_char)(void *), void *data);
void fftwf_forget_wisdom(void);


/*****************************************************************************************************************/

/* Parallel FFTW functions */
/*int fftwf_init_threads(void);
void fftwf_plan_with_nthreads(int nthreads);
void fftwf_cleanup_threads(void);*/

/*****************************************************************************************************************/

