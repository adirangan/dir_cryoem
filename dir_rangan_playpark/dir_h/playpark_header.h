#include <fcntl.h>
#include <math.h>
#include <pthread.h>
#include <omp.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <time.h>
#include <unistd.h>
#include <errno.h>
#include <assert.h>
#include <emmintrin.h>
#include <immintrin.h>
#include <cblas.h>
#include <complex.h>
#include <fftw3.h>
#include "finufft.h"

#define PI 3.141592653589793
#define CACHELINE 64 // assumed number of bytes per cache line, for alignment
#define CACHELINE_DBL (CACHELINE / sizeof(double))
#define MAX_THREADS 144
#define MAX_THREADS_P1 145
#define RSNSIZE 64
#define FNAMESIZE (4096+256)
#define PNAMESIZE (8192+512)
#define GLOBAL_MEMORY_GB_DEFAULT 4
#define MULTIPLEX_DIST 960
#define MULTIPLEX_2DIST (MULTIPLEX_DIST * 2)
#define BIT8 8
#define BITJ 16
#define bget_off(bmr,nr) (((bmr[(nr)/BIT8] >> (7-((nr)%BIT8))) & 1) ^ 1)
#define bget_0on(bmr,nr) (((bmr[(nr)/BIT8] >> (7-((nr)%BIT8))) & 1) ^ 0)
#define bget____(bmr,nr) ((int)(((bmr[(nr)/BIT8] >> (7-((nr)%BIT8))) & 1) << 1)-(int)1)
#define bset_off(bmr,nr) (bmr[(nr)/BIT8] &= ~(1 << (7 - ((nr)%BIT8))))
#define bset_0on(bmr,nr) (bmr[(nr)/BIT8] |=  (1 << (7 - ((nr)%BIT8))))
#define uchar_mask_size(A) (((A) + ((BIT8 - ((A) % BIT8)) % BIT8))/BIT8)

/* 5 --> 0101 */
#define FIVEMASK 0x5555555555555555LU
/* cachealign(val) returns the multiple of 64 geq to val */
#define CACHEALIGN(val) ((val + (CACHELINE - 1)) & (~(CACHELINE - 1)))
/* cachealign_dbl(val) returns the multiple of 8 geq to val */
#define CACHEALIGN_DBL(val) ((val + (CACHELINE_DBL - 1)) & (~(CACHELINE_DBL - 1)))

#define RET_SUCCESS 0
#define RET_NOMEM 1
#define RET_READ_FAIL 2
#define RET_INVALID_FORMAT 3
#define RET_CALC_NOT_YET_SUPPORTED 4
#define RET_INVALID_CMDLINE 5

/* ---------------------------------------------------------------- */

#define GLOBAL_omp_unused -1 /* unused */
#define GLOBAL_omp_off 0 /* parallelize manually */
#define GLOBAL_omp_0on 1 /* use omp in addition */

/* global variables used for timing */
#define GLOBAL_NTICKS 8
extern clock_t GLOBAL_t_start[GLOBAL_NTICKS];
extern clock_t GLOBAL_t_final[GLOBAL_NTICKS];
extern struct timeval GLOBAL_d_start[GLOBAL_NTICKS];
extern struct timeval GLOBAL_d_final[GLOBAL_NTICKS];
extern long GLOBAL_l_msec[GLOBAL_NTICKS],GLOBAL_l_ssec[GLOBAL_NTICKS],GLOBAL_l_usec[GLOBAL_NTICKS]; 
extern double GLOBAL_elct[GLOBAL_NTICKS],GLOBAL_elrt[GLOBAL_NTICKS];

/* global variables used for most routines */
extern char GLOBAL_CWD[FNAMESIZE];
extern char GLOBAL_dir_base[FNAMESIZE];
extern char GLOBAL_dir_name[FNAMESIZE];
extern char GLOBAL_dir_xpre[FNAMESIZE];
extern char GLOBAL_out_name[FNAMESIZE];
extern int GLOBAL_verbose;
extern int GLOBAL_thread_count;
extern int GLOBAL_omp_type;
extern int GLOBAL_1_icache_linesize;
extern int GLOBAL_1_icache_size;
extern int GLOBAL_1_icache_assoc;
extern int GLOBAL_1_dcache_linesize;
extern int GLOBAL_1_dcache_size;
extern int GLOBAL_1_dcache_assoc;
extern int GLOBAL_2_cache_linesize;
extern int GLOBAL_2_cache_size;
extern int GLOBAL_2_cache_assoc;
extern int GLOBAL_3_cache_linesize;
extern int GLOBAL_3_cache_size;
extern int GLOBAL_3_cache_assoc;
extern int GLOBAL_4_cache_linesize;
extern int GLOBAL_4_cache_size;
extern int GLOBAL_4_cache_assoc;
extern double GLOBAL_tolerance;
extern unsigned int GLOBAL_recursion_limit;
extern int addressable_1;
extern int addressable_0;
extern int addressable_int_length;
extern int addressable_int[128];

#define rup(A,B) ((A) + !!((A)%(B))*((B) - ((A)%(B))))
#define maximum(A,B) ((A) > (B) ? (A) : (B))
#define minimum(A,B) ((A) < (B) ? (A) : (B))
#define periodize(A,B,C) ((A) < (B) ? (A) + (C) - (B) : ((A) >= (C) ? (A) - (C) + (B) : (A)))
#define crop(A,B,C) ((A) < (B) ? (B) : ((A) > (C) ? (C) : (A)))
#define rand01 ((double)rand()/(double)RAND_MAX)
#define bsize(A) ((rup((A) + ((BITJ - ((A) % BITJ)) % BITJ),POPLENGTH))/BIT8)
#define psize(A) ((rup((A) + ((BITJ - ((A) % BITJ)) % BITJ),POPLENGTH))/POPLENGTH)

#define RGB3 3 
#define RGBA 4 
#define UCHAR_MAX 255 
#define POPLENGTH (MULTIPLEX_2DIST / 128 * sizeof(__m128i) * BIT8)

/* thread management */
extern int GLOBAL_nf; extern int GLOBAL_nf_cur; extern int GLOBAL_nf_ind; extern int GLOBAL_nf_opn;
extern int GLOBAL_tint[MAX_THREADS];
extern void *GLOBAL_tvp[MAX_THREADS][128];
extern pthread_t GLOBAL_threads[MAX_THREADS];
extern unsigned long long int GLOBAL_ops_f_[MAX_THREADS];
extern unsigned long long int GLOBAL_ops_f_sum;
extern unsigned long long int GLOBAL_ops_b_[MAX_THREADS];
extern unsigned long long int GLOBAL_ops_b_sum;

/* ---------------------------------------------------------------- */
unsigned long long int GLOBAL_n_malloc1;
int GLOBAL_n_malloc1_[GLOBAL_NTICKS];
/* manually managed memory stack */
extern unsigned char* wkspace;
extern unsigned char* wkspace_base;
extern int GLOBAL_wkspace_point;
struct wkspace_point
{
  struct wkspace_point * parent;
  struct wkspace_point * child;
  long long int check;
  long long int *point;
};
extern struct wkspace_point *wkspace_point_0,*wkspace_point_t;
extern unsigned long long int GLOBAL_memory_kb;
extern unsigned long long int GLOBAL_memory_mb;
extern unsigned long long int GLOBAL_memory_gb;
extern long long int wkspace_left;
extern long long int wkspace_used;

/* ---------------------------------------------------------------- */

/* RAND functions */
extern unsigned long int POW2RPOWPLUSRADD;
extern unsigned long int POW22RPOWMINUSONE;
extern int RCYCLENUM;
unsigned long int lrand();
double randn();
unsigned long int RGET(unsigned long int *);
double R01GET(unsigned long int *);
void RSEED_avd8(unsigned long int *);
double RNGET(unsigned long int *);
double RISIGET(unsigned long int *,double);

/* ---------------------------------------------------------------- */

#ifndef _MONOLITH
void fill_uchar_zero(unsigned char* iarr, size_t size);
void fill_uchar_ones(unsigned char* iarr, size_t size);
void fill_long_zero(long* larr, size_t size);
#endif /* _MONOLITH */
void ping();
void pong();
void MDA_io_test();
void MDA_write_i4(int n_d,int *d_,int *i4_,char *fname);
void MDA_read_i4(int *n_d_p,int **d_p,int **i4_p,char *fname);
void MDA_write_r8(int n_d,int *d_,double *r8_,char *fname);
void MDA_read_r8(int *n_d_p,int **d_p,double **r8_p,char *fname);
void MDA_write_ulli(int n_d,int *d_,unsigned long long int *ulli_,char *fname);
void MDA_read_ulli(int *n_d_p,int **d_p,unsigned long long int **ulli_p,char *fname);
void GLOBAL_ops_reset_one(int nt);
void GLOBAL_ops_count_one(int nt,unsigned long long int f,unsigned long long int b);
void GLOBAL_ops_addup_all();
void GLOBAL_ops_addup_one(int nt);
void GLOBAL_ops_printf_all(int verbose,char *prefix);
void GLOBAL_tic(int nx);
void GLOBAL_toc(int nx,int verbose,char *prefix);
void GLOBAL_ops_toc(int nt,int nx,int verbose,char *prefix);
void GLOBAL_pthread_tic();
void GLOBAL_pthread_toc();
void GLOBAL_pthread_tuc();
void hsv2rgb(double h,double s,double v,double *r,double *g,double *b);
void colorscale(double val,double valmin,double valmax,double *rcolor,double *gcolor,double *bcolor);
int WritePNMfile_color(double *ra,int rows,int cols,double min,double max,char *filename);
void updateglobals(char *vname);
void read_input();
void* malloc1(size_t size);
void free1(void **vp);
unsigned char *wkspace_alloc_nocheck(unsigned long long int size);
struct wkspace_point * wkspace_make_point();
struct wkspace_point * wkspace_set_point(struct wkspace_point *w);
void wkspace_printf_point(struct wkspace_point *w);
long long int wkspace_check_point(struct wkspace_point *w);
unsigned char* wkspace_alloc(unsigned long long int size);
unsigned char* wkspace_all0c(unsigned long long int size);
void wkspace_reset(void* new_base);
void wkspace_printf();
/* %%%%%%%%%%%%%%%% */
float erfcln_single_f(float f_0in);
void erfcln_f(int n_r,float *f_0in_,float **f_out_p_);
void erfcln_f_test();
double erfcln_single_d(double d_0in);
void erfcln_d(int n_r,double *d_0in_,double **d_out_p_);
void erfcln_d_test();
double z_to_lp_single_d(double d_0in);
void z_to_lp_d(int n_r,double *d_0in_,double **d_out_p_);
void z_to_lp_d_test();
/* %%%%%%%%%%%%%%%% */
void array_printf(void *v_,char *type,int n_row,int n_col,char *prefix);
void array_fprintf(char *fname,void *v_,char *type,int n_row,int n_col,char *prefix);
void bitstring_from_uchar(unsigned char *w_, char *str_, int k);
void bitstring_from_uchar_printf(unsigned char *w_,int nrows,int ncols,char *prefix);
/* %%%%%%%%%%%%%%%% */
void icumsum(unsigned long long int ulli_length,int *i_0in_,int **i_out_p_);
double ifnormn(unsigned long long int ulli_length,int *i_0_,int *i_1_);
float ffnormn(unsigned long long int ulli_length,float *f_0_,float *f_1_);
float ffnorm(unsigned long long int ulli_length,float *f_0_,float *f_1_);
double dfnormn(unsigned long long int ulli_length,double *d_0_,double *d_1_);
double dfnorm(unsigned long long int ulli_length,double *d_0_,double *d_1_);
double cfnorm(unsigned long long int ulli_length,float complex *c_0_,float complex *c_1_);
double zfnorm(unsigned long long int ulli_length,double complex *z_0_,double complex *z_1_);
void array_stats(void *v_,char *type,unsigned long long int ulli_length,void *max_p,void *min_p,double *mean_p,double *stdev_p);
/* %%%%%%%%%%%%%%%% */
int is_internal_maximum(int n_Z,double *Z_,int index);
void find_local_maxima(int n_Z,double *Z_,int *n_index_p,int **index_p_);
void find_internal_maximum(int verbose,int n_Z,double *Z_,double Z_min,double *zone_max_p,int *zone_max_index_p);
void find_internal_maximum_test();
/* %%%%%%%%%%%%%%%% */
void dtranspose_bruteforce(int n_row_A,int n_col_A,double *d_0in__,double *d_out__);
void dtranspose_block_AtoB(const int n_row_A,const int n_col_A,const double* A_,double* B_,const int block_size);
void dtranspose_block_BtoA(const int n_row_A,const int n_col_A,const double* A_,double* B_,const int block_size);
void dtranspose(const int n_row_A,const int n_col_A,const double* A_,const double* B_);
void dtranspose_test();
void cntranspose_bruteforce(int n_row_A,int n_col_A,float complex *c_0in__,float complex *c_out__);
void cntranspose_block_AtoB(const int n_row_A,const int n_col_A,const float complex* A_,float complex* B_,const int block_size);
void cntranspose_block_BtoA(const int n_row_A,const int n_col_A,const float complex* A_,float complex* B_,const int block_size);
void cntranspose(const int n_row_A,const int n_col_A,const float complex* A_,float complex* B_);
void cntranspose_test();
void cctranspose_bruteforce(int n_row_A,int n_col_A,float complex *c_0in__,float complex *c_out__);
void cctranspose_block_AtoB(const int n_row_A,const int n_col_A,const float complex* A_,float complex* B_,const int block_size);
void cctranspose_block_BtoA(const int n_row_A,const int n_col_A,const float complex* A_,float complex* B_,const int block_size);
void cctranspose(const int n_row_A,const int n_col_A,const float complex* A_,float complex* B_);
void cctranspose_test();
/* %%%%%%%%%%%%%%%% */
void fftf_1d_bruteforce(int n_row,float complex *c_0in_,float complex *c_out_,int sgn);
void fftf__1d_bruteforce(int n_row,int n_col,float complex *c_0in_,float complex *c_out_,int sgn);
void fftw3f__1d_andplan(int n_row,int n_col,float complex *c_0in_,float complex *c_out_,int sgn);
void fftw3f__1d_test();
/* %%%%%%%%%%%%%%%% */
int finufft_1d_test();
/* %%%%%%%%%%%%%%%% */
void dp_ps_bruteforce(int n_col_X,float *f_A_,float *f_B_,float *f_C_);
void dp_pd_immintrin_loadu(int n_col_X,double *d_A_,double *d_B_,double *d_C_);
void dp_ps_immintrin_loadu(int n_col_X,float *f_A_,float *f_B_,float *f_C_);
void dp_ps_mult_immintrin_loadu(int n_row_A,int n_col_X,float *f_A_trn__,int n_row_B,float *f_B_trn__,float **f_C_p__);
void dp_ps_mult_bruteforce(int n_row_A,int n_col_X,float *f_A_trn__,int n_row_B,float *f_B_trn__,float **f_C_p__);
void dp_ps_mult_immintrin_load1(int n_row_A,int n_col_X,float *f_A_trn__,int n_row_B,float *f_B_trn__,float **f_C_p__);
void dp_ps_mult_immintrin_fma2(int n_row_A,int n_col_X,__m256 *ps_A_trn__,int n_row_B,__m256 *ps_B_trn__,float **f_C_p__);
void dp_ps_mult_immintrin_fma(int n_row_A,int n_col_X,__m256 *ps_A_trn__,int n_row_B,__m256 *ps_B_trn__,float **f_C_p__);
void dp_ps_mult_immintrin_avx(int n_row_A,int n_col_X,__m256 *ps_A_trn__,int n_row_B,__m256 *ps_B_trn__,float **f_C_p__);
void dp_ps_mult_immintrin_test();
void dp_ps_single_test();
/* %%%%%%%%%%%%%%%% */
int iquicksort_partition_index(int *i_,int stride,int *index_,int l,int r);
unsigned int iquicksort_index(unsigned int recursion_level,int *i_,int stride,int *index_,int l,int r);
void iquicksort_index_driver(int n_i,int *i_,int stride,int *i_workspace_,int *index_);
void iquicksort_index_driver_test();
void iquicksort_index_index_driver(int n_i,int *i_,int stride,int *i_workspace_,int *index_orig_from_sort_,int *index_sort_from_orig_,int *index_workspace_);
void iquicksort_index_index_driver_test();
void irandperm(int n_i,int **index_p_,unsigned long long int *rseed);
void irandperm_test();
int fquicksort_partition_index(float *f_,int stride,int *index_,int l,int r);
unsigned int fquicksort_index(unsigned int recursion_level,float *f_,int stride,int *index_,int l,int r);
void fquicksort_index_driver(int n_f,float *f_,int stride,float *f_workspace_,int *index_);
void fquicksort_index_driver_test();
void fquicksort_index_index_driver(int n_f,float *f_,int stride,float *f_workspace_,int *index_orig_from_sort_,int *index_sort_from_orig_,int *index_workspace_);
void fquicksort_index_index_driver_test();
int dquicksort_partition_index(double *d_,int stride,int *index_,int l,int r);
unsigned int dquicksort_index(unsigned int recursion_level,double *d_,int stride,int *index_,int l,int r);
void dquicksort_index_driver(int n_d,double *d_,int stride,double *d_workspace_,int *index_);
void dquicksort_index_driver_test();
void dquicksort_index_index_driver(int n_d,double *d_,int stride,double *d_workspace_,int *index_orig_from_sort_,int *index_sort_from_orig_,int *index_workspace_);
void dquicksort_index_index_driver_test();
/* %%%%%%%%%%%%%%%% */
void get_xdrop_logscale(double lrij,double lcij,double gamma,int *rdrop_p,int *cdrop_p);
int get_xdrop_logscale_length(double lrij,double lcij,double gamma);
void get_xdrop_logscale_array(double lrij,double lcij,double gamma,int *length_p,int **rdrop_p_,int **cdrop_p_,int **rkeep_p_,int **ckeep_p_);
void get_xdrop_logscale_array_test();
/* %%%%%%%%%%%%%%%% */
void gumbel_nll(int n_x,double *x_,double *g_,double tol,double **nll_p_,double *nll_sum_p);
void gumbel_nll_wrap(int n_g,double *g_,double *nll_sum_p,void *args);
void gumbel_nll_test();
/* %%%%%%%%%%%%%%%% */
void dexcluster_nonbinary_f_recursive_helper_QR__
(
  int verbose
 ,int flag_rdrop_vs_rcdrop
 ,int flag_force_create
 ,int n_r
 ,int n_c
 ,float * E_base_rc__
 ,int n_r_index
 ,int *r_index_
 ,int n_c_index
 ,int *c_index_
 ,double gamma
 ,int n_shuffle
 ,char *fname_trace__
 ,char *fname_xdrop__
 ,char *fname_QR__
);
void dexcluster_nonbinary_recursive_helper_ZR__
(
 int verbose
 ,int flag_rdrop_vs_rcdrop
 ,char *fname_trace__
 ,char *fname_xdrop__
 ,char *fname_QR__
 ,double p_use
 ,int n_member_lob
 ,double *nlp_ZR_max_p
 ,int *nlp_ZR_index_p
 ,int *n_r_rtn_index_p
 ,int **r_rtn_index_p_
 ,int *n_r_rmv_index_p
 ,int **r_rmv_index_p_
 ,int *n_c_rtn_index_p
 ,int **c_rtn_index_p_
 ,int *n_c_rmv_index_p
 ,int **c_rmv_index_p_
 ,double *nlp_gumb_opt_p
 ,double *nlp_gumb_emp_p
);
void dexcluster_nonbinary_f_recursive
(
 char *dir_trunk_0in
 ,char *dir_out_0in
 ,char *prefix_base_0in
 ,int n_r
 ,int n_c
 ,float *E_base_rc__
 ,int n_r_index_0in
 ,int *r_index_0in_
 ,int n_c_index_0in
 ,int *c_index_0in_
 ,double gamma_0in
 ,int n_shuffle_0in
 ,double p_set_0in
 ,int n_member_lob_0in
 ,double p_prev_0in
 ,int flag_force_create_0in
 ,char ***output_label_p_
 ,char ***lpFmax_label_p_
 ,char ***lpnext_label_p_
);
void dexcluster_nonbinary_f_recursive_test();




