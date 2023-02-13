#ifndef _MONOLITH
#include "ampm_header.h"
#endif /* _MONOLITH */

/* %%%%%%%% */
/* testing transposes found at: */
/* https://stackoverflow.com/questions/16737298/what-is-the-fastest-way-to-transpose-a-matrix-in-c */
/* %%%%%%%% */
void transpose_float_rc_bruteforce(float *A__, float *B__, const int n_row, const int n_col)
{
  int nr=0,nc=0,na=0,nb=0;
  for (nr=0;nr<n_row;nr++){
    nb=nr;
    for (nc=0;nc<n_col;nc++){
      B__[na++] = A__[nb]; nb+=n_row;
      /* for (nc=0;nc<n_col;nc++){ } */}
    /* for (nr=0;nr<n_row;nr++){ } */}
}

void transpose_float_cr_bruteforce(float *A__, float *B__, const int n_row, const int n_col)
{
  int nr=0,nc=0,na=0,nb=0;
  for (nc=0;nc<n_col;nc++){
    nb=nc;
    for (nr=0;nr<n_row;nr++){
      B__[nb] = A__[na++]; nb+=n_col;
      /* for (nr=0;nr<n_row;nr++){ } */}
    /* for (nc=0;nc<n_col;nc++){ } */}
}

void transpose_float_block0(float *A__, float *B__, const int lda, const int ldb, const int block_size)
{
  int i=0,j=0;
  for(i=0; i<block_size; i++) {
    for(j=0; j<block_size; j++) {
      B__[j*ldb + i] = A__[i*lda +j];
      /* for(j=0; j<block_size; j++) { } */}
    /* for(i=0; i<block_size; i++) { } */}
}

void transpose_float_block1(float *A__, float *B__, const int m, const int n, const int lda, const int ldb, const int block_size)
{
  int i=0,j=0;
  for(i=0; i<n; i+=block_size) {
    for(j=0; j<m; j+=block_size) {
      transpose_float_block0(&A__[i*lda +j], &B__[j*ldb + i], lda, ldb, block_size);
      /* for(j=0; j<m; j+=block_size) { } */}
    /* for(i=0; i<n; i+=block_size) { } */}
}

void transpose_ps_block0(float *A, float *B, const int lda, const int ldb)
{
  __m128 row1 = _mm_load_ps(&A[0*lda]);
  __m128 row2 = _mm_load_ps(&A[1*lda]);
  __m128 row3 = _mm_load_ps(&A[2*lda]);
  __m128 row4 = _mm_load_ps(&A[3*lda]);
  _MM_TRANSPOSE4_PS(row1, row2, row3, row4);
  _mm_store_ps(&B[0*ldb], row1);
  _mm_store_ps(&B[1*ldb], row2);
  _mm_store_ps(&B[2*ldb], row3);
  _mm_store_ps(&B[3*ldb], row4);
}

void transpose_ps_block1(float *A, float *B, const int m, const int n, const int lda, const int ldb ,const int block_size)
{
  int i=0,j=0,max_i2=0,max_j2=0,i2=0,j2=0;
  for(i=0; i<n; i+=block_size) {
    for(j=0; j<m; j+=block_size) {
      max_i2 = i+block_size < n ? i + block_size : n;
      max_j2 = j+block_size < m ? j + block_size : m;
      for(i2=i; i2<max_i2; i2+=4) {
	for(j2=j; j2<max_j2; j2+=4) {
	  transpose_ps_block0(&A[i2*lda +j2], &B[j2*ldb + i2], lda, ldb);
	  /* for(j2=j; j2<max_j2; j2+=4) { } */}
	/* for(i2=i; i2<max_i2; i2+=4) { } */}
      /* for(j=0; j<m; j+=block_size) { } */}
    /* for(i=0; i<n; i+=block_size) { } */}
}

void transpose_float_test()
{
  int verbose=1;
  const int n_row = 17001;
  const int n_col = 15001;
  int n_row_rup = rup(n_row, 16);
  int n_col_rup = rup(n_col, 16);
  unsigned long long int ulli_total = (unsigned long long int)n_row_rup*(unsigned long long int)n_col_rup;
  int n_row_sub = 5;
  int n_col_sub = 6;
  unsigned long long int tab=0;
  float *A__ = NULL;
  float *B_bf__ = NULL;
  float *B_al__ = NULL;
  /* %%%% */
  int nr=0,nc=0;
  /* %%%% */
  int ntick=0,n_tick = 8;
  clock_t t_start_[n_tick];
  clock_t t_final_[n_tick];
  struct timeval d_start_[n_tick];
  struct timeval d_final_[n_tick];
  long l_msec_[n_tick],l_ssec_[n_tick],l_usec_[n_tick];
  double elct_[n_tick],elrt_[n_tick];
  /* %%%% */
  if (verbose){ printf(" %% [entering test_transpose_float]\n");}
  /* %%%% */
  A__ = (float*)_mm_malloc(sizeof(float)*ulli_total, 64);
  for (nr=0;nr<n_row;nr++){
    for (nc=0;nc<n_col;nc++){
      A__[nr + nc*n_row_rup] = nr*100 + nc;
      /* for (nc=0;nc<n_col;nc++){ } */}
    /* for (nr=0;nr<n_row;nr++){ } */}
  if (verbose>1){ array_sub_printf(A__,"float",n_row_rup,n_row_sub,n_col_rup,n_col_sub," %% A_sub__: ");}
  B_bf__ = (float*)_mm_malloc(sizeof(float)*ulli_total, 64);
  B_al__ = (float*)_mm_malloc(sizeof(float)*ulli_total, 64);
  /* %%%% */
  memset(B_bf__,0,ulli_total*sizeof(float));
  local_tic(0,t_start_,d_start_);
  tab = ulli_total;
  transpose_float_rc_bruteforce(A__,B_bf__,n_row_rup,n_col_rup);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose,"transpose_float_rc_bruteforce: ");
  if (verbose>1){ array_sub_printf(B_bf__,"float",n_col_rup,n_col_sub,n_row_rup,n_row_sub," %% B_sub_bf__: ");}
  /* %%%% */
  memset(B_al__,0,ulli_total*sizeof(float));
  local_tic(0,t_start_,d_start_);
  tab = ulli_total;
  transpose_float_cr_bruteforce(A__,B_al__,n_row_rup,n_col_rup);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose,"transpose_float_cr_bruteforce: ");
  if (verbose){ printf(" %% error: %0.16f\n",ffnormn(ulli_total,B_bf__,B_al__));}
  /* %%%% */
  memset(B_al__,0,ulli_total*sizeof(float));
  local_tic(0,t_start_,d_start_);
  tab = ulli_total;
  transpose_float_block1(A__,B_al__,n_row_rup,n_col_rup,n_row_rup,n_col_rup,4);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose,"transpose_float_block1(4): ");
  if (verbose){ printf(" %% error: %0.16f\n",ffnormn(ulli_total,B_bf__,B_al__));}
  /* %%%% */
  memset(B_al__,0,ulli_total*sizeof(float));
  local_tic(0,t_start_,d_start_);
  tab = ulli_total;
  transpose_float_block1(A__,B_al__,n_row_rup,n_col_rup,n_row_rup,n_col_rup,8);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose,"transpose_float_block1(8): ");
  if (verbose){ printf(" %% error: %0.16f\n",ffnormn(ulli_total,B_bf__,B_al__));}
  /* %%%% */
  memset(B_al__,0,ulli_total*sizeof(float));
  local_tic(0,t_start_,d_start_);
  tab = ulli_total;
  transpose_float_block1(A__,B_al__,n_row_rup,n_col_rup,n_row_rup,n_col_rup,16);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose,"transpose_float_block1(16): ");
  if (verbose){ printf(" %% error: %0.16f\n",ffnormn(ulli_total,B_bf__,B_al__));}
  /* %%%% */
  memset(B_al__,0,ulli_total*sizeof(float));
  local_tic(0,t_start_,d_start_);
  tab = ulli_total;
  transpose_ps_block1(A__,B_al__,n_row_rup,n_col_rup,n_row_rup,n_col_rup,16);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose,"transpose_ps_block1(16): ");
  if (verbose){ printf(" %% error: %0.16f\n",ffnormn(ulli_total,B_bf__,B_al__));}
  /* %%%% */
  _mm_free(A__); A__=NULL;
  _mm_free(B_bf__); B_bf__=NULL;
  _mm_free(B_al__); B_al__=NULL;
  /* %%%% */
  if (verbose){ printf(" %% [finished test_transpose_float]\n");}
}

