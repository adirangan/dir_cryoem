#ifndef _MONOLITH
#include "ampm_header.h"
#endif /* _MONOLITH */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* testing matrix multiplication */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
void dp_ps_bruteforce(int n_col_X,float *f_A_,float *f_B_,float *f_C_)
{
  int ncol_X=0;
  float *f_A_point0_=NULL,*f_B_point0_=NULL;
  float f_ac0=0;
  f_A_point0_ = f_A_;
  f_B_point0_ = f_B_;
  for (ncol_X=0;ncol_X<n_col_X;ncol_X++){
    f_ac0 += *(f_A_point0_++) * *(f_B_point0_++);
    /* for (ncol_X=0;ncol_X<n_col_X;ncol_X++){ } */}
  *f_C_ = f_ac0;
}

void dp_ps_mult_bruteforce(int n_row_A,int n_col_X,float *f_A_trn__,int n_row_B,float *f_B_trn__,float **f_C_p_)
{
  int n_col_X_rup = rup(n_col_X,8);
  int nrow_A=0,nrow_B=0,ncol_X=0;
  float *f_C__=NULL;
  f_C__=NULL;
  if (f_C_p_!=NULL){
    if ( (*f_C_p_)==NULL ){ (*f_C_p_) = (float *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
    f_C__ = *f_C_p_;
    /* if (f_C_p_!=NULL){ } */}
  if (f_C__!=NULL){
    for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
	for (ncol_X=0;ncol_X<n_col_X;ncol_X++){
	  f_C__[nrow_A + nrow_B*n_row_A] += f_A_trn__[ncol_X + nrow_A*n_col_X_rup] * f_B_trn__[ncol_X + nrow_B*n_col_X_rup];
	  /* for (ncol_X=0;ncol_X<n_col_X;ncol_X++){ } */}
	/* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ } */}
    /* if (f_C__!=NULL){ } */}
}

void dp_ps_mult_cblas_sgemm(int n_row_A,int n_col_X,float *f_A_trn__,int n_row_B,float *f_B_trn__,float **f_C_p_)
{
  int n_col_X_rup = rup(n_col_X,8);
  float *f_C__=NULL;
  f_C__=NULL;
  if (f_C_p_!=NULL){
    if ( (*f_C_p_)==NULL ){ (*f_C_p_) = (float *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
    f_C__ = *f_C_p_;
    /* if (f_C_p_!=NULL){ } */}
  if (f_C__!=NULL){
#ifdef _CBLAS
    cblas_sgemm(CblasColMajor,CblasTrans,CblasNoTrans,n_row_A,n_row_B,n_col_X,1.0,f_A_trn__,n_col_X_rup,f_B_trn__,n_col_X_rup,0.0,f_C__,n_row_A);
#endif /* _CBLAS */
#ifndef _CBLAS
    printf(" %% Warning, _CBLAS not defined in dp_ps_mult_cblas_sgemm, using bruteforce instead.\n");
    dp_ps_mult_bruteforce(n_row_A,n_col_X,f_A_trn__,n_row_B,f_B_trn__,f_C_p_);
#endif /* _CBLAS */
  /* if (f_C__!=NULL){ } */}
}

void dp_ps_mult_cblas_sdot(int n_row_A,int n_col_X,float *f_A_trn__,int n_row_B,float *f_B_trn__,float **f_C_p_)
{
  int n_col_X_rup = rup(n_col_X,8);
  int nrow_A=0,nrow_B=0,ncol_X=0;
  float *f_C__=NULL;
  f_C__=NULL;
  if (f_C_p_!=NULL){
    if ( (*f_C_p_)==NULL ){ (*f_C_p_) = (float *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
    f_C__ = *f_C_p_;
    /* if (f_C_p_!=NULL){ } */}
  if (f_C__!=NULL){
#ifdef _CBLAS
    for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
	f_C__[nrow_A + nrow_B*n_row_A] = cblas_sdot(n_col_X,f_A_trn__+nrow_A*n_col_X_rup,1,f_B_trn__+nrow_B*n_col_X_rup,1);
	/* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ } */}
#endif /* _CBLAS */
#ifndef _CBLAS
    printf(" %% Warning, _CBLAS not defined in dp_ps_mult_cblas_sdot, using bruteforce instead.\n");
    dp_ps_mult_bruteforce(n_row_A,n_col_X,f_A_trn__,n_row_B,f_B_trn__,f_C_p_);
#endif /* _CBLAS */
    /* if (f_C__!=NULL){ } */}
}

void dp_ps_immintrin_loadu_avx(int n_col_X,float *f_A_,float *f_B_,float *f_C_)
{
  /* does not assume alignment */
  int n_col_X_rup = rup(n_col_X,8);
  int n_col_X_256 = n_col_X_rup/8; //%<-- 8 floats per __m256. ;
  int ncol_X=0,ncol_X_256=0;
  float *f_A_point0_=NULL,*f_B_point0_=NULL;
  __m256 ps_A0,ps_B0;
  __m256 ps_ac0;
  float *f_ac0_ = (float *) &ps_ac0;
  float f_ac0=0;
  f_A_point0_ = f_A_;
  f_B_point0_ = f_B_;
  ps_ac0  = _mm256_set1_ps((float)0.0);
  for (ncol_X_256=0;ncol_X_256<n_col_X_256-1;ncol_X_256++){
    ps_A0 = _mm256_loadu_ps(f_A_point0_);
    ps_B0 = _mm256_loadu_ps(f_B_point0_);
    ps_ac0 = _mm256_dp_ps(ps_A0,ps_B0,0xF1);
    f_ac0 += f_ac0_[0] + f_ac0_[4];
    f_A_point0_+=8;f_B_point0_+=8;
    /* for (ncol_X_256=0;ncol_X_256<n_col_X_256-1;ncol_X_256++){ } */}
  f_A_point0_ = &(f_A_[8*ncol_X_256]);
  f_B_point0_ = &(f_B_[8*ncol_X_256]);
  for (ncol_X=8*ncol_X_256;ncol_X<n_col_X;ncol_X++){
    f_ac0 += *(f_A_point0_++) * *(f_B_point0_++);
    /* for (ncol_X=8*ncol_X_256;ncol_X<n_col_X;ncol_X++){ } */}
  *f_C_ = f_ac0;
}

void dp_ps_mult_immintrin_loadu_avx(int n_row_A,int n_col_X,float *f_A_trn__,int n_row_B,float *f_B_trn__,float **f_C_p_)
{
  /* does not assume alignment */
  int nrow_A=0,nrow_B=0;
  float *tmp_f_A_point0_=NULL;
  float *tmp_f_B_point0_=NULL;
  float *f_C__=NULL;
  f_C__=NULL;
  if (f_C_p_!=NULL){
    if ( (*f_C_p_)==NULL ){ (*f_C_p_) = (float *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
    f_C__ = *f_C_p_;
    /* if (f_C_p_!=NULL){ } */}
  if (f_C__!=NULL){
    for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
	tmp_f_A_point0_ = &(f_A_trn__[nrow_A*n_col_X]);
	tmp_f_B_point0_ = &(f_B_trn__[nrow_B*n_col_X]);
	dp_ps_immintrin_loadu_avx(n_col_X,tmp_f_A_point0_,tmp_f_B_point0_,&(f_C__[nrow_A + nrow_B*n_row_A]));
	/* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ } */}
    /* if (f_C__!=NULL){ } */}
}

void dp_ps_immintrin_loadu_fma(int n_col_X,float *f_A_,float *f_B_,float *f_C_)
{
  /* does not assume alignment */
  int n_col_X_rup = rup(n_col_X,8);
  int n_col_X_256 = n_col_X_rup/8; //%<-- 8 floats per __m256. ;
  int ncol_X=0,ncol_X_256=0;
  float *f_A_point0_=NULL,*f_B_point0_=NULL;
  __m256 ps_A0,ps_B0;
  __m256 ps_ac0;
  float *f_ac0_ = (float *) &ps_ac0;
  float f_ac0=0;
  f_A_point0_ = f_A_;
  f_B_point0_ = f_B_;
  ps_ac0  = _mm256_set1_ps((float)0.0);
  for (ncol_X_256=0;ncol_X_256<n_col_X_256-1;ncol_X_256++){
    ps_A0 = _mm256_loadu_ps(f_A_point0_);
    ps_B0 = _mm256_loadu_ps(f_B_point0_);
#ifdef _FMA
    ps_ac0 = _mm256_fmadd_ps(ps_A0,ps_B0,ps_ac0);
#endif /* _FMA */
    f_A_point0_+=8;f_B_point0_+=8;
    /* for (ncol_X_256=0;ncol_X_256<n_col_X_256-1;ncol_X_256++){ } */}
  f_ac0 = f_ac0_[0] + f_ac0_[1] + f_ac0_[2] + f_ac0_[3] + f_ac0_[4] + f_ac0_[5] + f_ac0_[6] + f_ac0_[7];
  f_A_point0_ = &(f_A_[8*ncol_X_256]);
  f_B_point0_ = &(f_B_[8*ncol_X_256]);
  for (ncol_X=8*ncol_X_256;ncol_X<n_col_X;ncol_X++){
    f_ac0 += *(f_A_point0_++) * *(f_B_point0_++);
    /* for (ncol_X=8*ncol_X_256;ncol_X<n_col_X;ncol_X++){ } */}
  *f_C_ = f_ac0;
}

void dp_ps_mult_immintrin_loadu_fma(int n_row_A,int n_col_X,float *f_A_trn__,int n_row_B,float *f_B_trn__,float **f_C_p_)
{
  /* does not assume alignment */
  int nrow_A=0,nrow_B=0;
  float *tmp_f_A_point0_=NULL;
  float *tmp_f_B_point0_=NULL;
  float *f_C__=NULL;
  f_C__=NULL;
  if (f_C_p_!=NULL){
    if ( (*f_C_p_)==NULL ){ (*f_C_p_) = (float *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
    f_C__ = *f_C_p_;
    /* if (f_C_p_!=NULL){ } */}
  if (f_C__!=NULL){
    for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
	tmp_f_A_point0_ = &(f_A_trn__[nrow_A*n_col_X]);
	tmp_f_B_point0_ = &(f_B_trn__[nrow_B*n_col_X]);
	dp_ps_immintrin_loadu_fma(n_col_X,tmp_f_A_point0_,tmp_f_B_point0_,&(f_C__[nrow_A + nrow_B*n_row_A]));
	/* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ } */}
    /* if (f_C__!=NULL){ } */}
}

void dp_ps_mult_immintrin_load1_fma(int n_row_A,int n_col_X,float *f_A_trn__,int n_row_B,float *f_B_trn__,float **f_C_p_)
{
  int n_col_X_rup = rup(n_col_X,8);
  int n_col_X_256 = n_col_X_rup/8; //%<-- 8 floats per __m256. ;
  int nrow_A=0,nrow_B=0,ncol_X=0;
  float *tmp_f_A_point0_=NULL;
  float *tmp_f_B_point0_=NULL;
  __m256 ps_A0,ps_B0;
  __m256 ps_ac0;
  float *tmp_f_ac0_ = (float *) &ps_ac0;
  float tmp_f_ac0=0;
  float *f_C__=NULL;
  f_C__=NULL;
  if (f_C_p_!=NULL){
    if ( (*f_C_p_)==NULL ){ (*f_C_p_) = (float *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
    f_C__ = *f_C_p_;
    /* if (f_C_p_!=NULL){ } */}
  if (f_C__!=NULL){
    for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
	tmp_f_A_point0_ = &(f_A_trn__[nrow_A*n_col_X_rup]);
	tmp_f_B_point0_ = &(f_B_trn__[nrow_B*n_col_X_rup]);
	ps_ac0 = _mm256_set1_ps((float)0.0);
	for (ncol_X=0;ncol_X<n_col_X_256;ncol_X++){
	  ps_A0 = _mm256_load_ps(tmp_f_A_point0_);
	  ps_B0 = _mm256_load_ps(tmp_f_B_point0_);
#ifdef _FMA
	  ps_ac0 = _mm256_fmadd_ps(ps_A0,ps_B0,ps_ac0);
#endif /* _FMA */
#ifdef _FMA
	  tmp_f_A_point0_+=8;tmp_f_B_point0_+=8;
#endif /* _FMA */
	  /* for (ncol_X=0;ncol_X<n_col_X_256;ncol_X++){ } */}
	tmp_f_ac0 = tmp_f_ac0_[0] + tmp_f_ac0_[1] + tmp_f_ac0_[2] + tmp_f_ac0_[3] + tmp_f_ac0_[4] + tmp_f_ac0_[5] + tmp_f_ac0_[6] + tmp_f_ac0_[7];
	f_C__[nrow_A + nrow_B*n_row_A] = tmp_f_ac0;
	/* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ } */}
    /* if (f_C__!=NULL){ } */}
}

void dp_ps_mult_immintrin_fma2(int n_row_A,int n_col_X,__m256 *ps_A_trn__,int n_row_B,__m256 *ps_B_trn__,float **f_C_p_)
{
  int n_col_X_rup = rup(n_col_X,8);
  int n_col_X_256 = n_col_X_rup/8; //%<-- 8 floats per __m256. ;
  int nrow_A=0,nrow_B=0,ncol_X=0;
  __m256 *tmp_ps_A_trn__=NULL;
  __m256 *tmp_ps_B_trn__=NULL;
  __m256 *tmp_ps_A_trn_=NULL;
  __m256 *tmp_ps_B_trn_=NULL;
  __m256 tmp_ps_A_trn,tmp_ps_B_trn;
  __m256 ps_ac0,ps_ac1;
  float *tmp_f_ac0_ = (float *) &ps_ac0;
  float *tmp_f_ac1_ = (float *) &ps_ac1;
  float tmp_f_ac0=0,tmp_f_ac1=0; 
  float *f_C__=NULL;
  f_C__=NULL;
  if (f_C_p_!=NULL){
    if ( (*f_C_p_)==NULL ){ (*f_C_p_) = (float *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
    f_C__ = *f_C_p_;
    /* if (f_C_p_!=NULL){ } */}
  if (f_C__!=NULL){
    printf(" n_col_X %d n_col_X_rup %d n_col_X_256 %d \n",n_col_X,n_col_X_rup,n_col_X_256);
    tmp_ps_B_trn__ = ps_B_trn__;
    for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      tmp_ps_A_trn__ = ps_A_trn__;
      for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
	tmp_ps_A_trn_ = tmp_ps_A_trn__;
	tmp_ps_B_trn_ = tmp_ps_B_trn__;
	ps_ac0 = _mm256_set1_ps((float)0.0);
	ps_ac1 = _mm256_set1_ps((float)0.0);
	for (ncol_X=0;ncol_X<n_col_X_256-1;ncol_X+=2){
#ifdef _FMA
	  ps_ac0 = _mm256_fmadd_ps(*(tmp_ps_A_trn_),*(tmp_ps_B_trn_),ps_ac0);
#endif /* _FMA */
#ifdef _FMA
	  ps_ac1 = _mm256_fmadd_ps(*(tmp_ps_A_trn_+1),*(tmp_ps_B_trn_+1),ps_ac1);
#endif /* _FMA */
	  tmp_ps_A_trn_+=2;
	  tmp_ps_B_trn_+=2;
	  /* for (ncol_X=0;ncol_X<n_col_X_256-1;ncol_X+=2){ } */}
	if ( (nrow_A==0) && (nrow_B==0) ){ printf(" ncol_X %d\n",ncol_X);}
#ifdef _FMA
	if (n_col_X_256%2==1){ ps_ac0 = _mm256_fmadd_ps(*(tmp_ps_A_trn_++),*(tmp_ps_B_trn_++),ps_ac0);}
#endif /* _FMA */
	tmp_f_ac0 = tmp_f_ac0_[0] + tmp_f_ac0_[1] + tmp_f_ac0_[2] + tmp_f_ac0_[3] + tmp_f_ac0_[4] + tmp_f_ac0_[5] + tmp_f_ac0_[6] + tmp_f_ac0_[7];
	tmp_f_ac1 = tmp_f_ac1_[0] + tmp_f_ac1_[1] + tmp_f_ac1_[2] + tmp_f_ac1_[3] + tmp_f_ac1_[4] + tmp_f_ac1_[5] + tmp_f_ac1_[6] + tmp_f_ac1_[7];
	f_C__[nrow_A + nrow_B*n_row_A] = tmp_f_ac0 + tmp_f_ac1;
	tmp_ps_A_trn__ += n_col_X_256;
	/* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
      tmp_ps_B_trn__ += n_col_X_256;
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ } */}
    /* if (f_C__!=NULL){ } */}
}

void dp_ps_mult_immintrin_fma(int n_row_A,int n_col_X,__m256 *ps_A_trn__,int n_row_B,__m256 *ps_B_trn__,float **f_C_p_)
{
  int n_col_X_rup = rup(n_col_X,8);
  int n_col_X_256 = n_col_X_rup/8; //%<-- 8 floats per __m256. ;
  int nrow_A=0,nrow_B=0,ncol_X=0;
  __m256 *tmp_ps_A_trn__=NULL;
  __m256 *tmp_ps_B_trn__=NULL;
  __m256 *tmp_ps_A_trn_=NULL;
  __m256 *tmp_ps_B_trn_=NULL;
  __m256 tmp_ps_A_trn,tmp_ps_B_trn;
  __m256 ps_acc;
  float *tmp_f_ = (float *) &ps_acc;
  float tmp_f=0; 
  float *f_C__=NULL;
  f_C__=NULL;
  if (f_C_p_!=NULL){
    if ( (*f_C_p_)==NULL ){ (*f_C_p_) = (float *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
    f_C__ = *f_C_p_;
    /* if (f_C_p_!=NULL){ } */}
  if (f_C__!=NULL){
    tmp_ps_B_trn__ = ps_B_trn__;
    for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      tmp_ps_A_trn__ = ps_A_trn__;
      for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
	ps_acc = _mm256_set1_ps((float)0.0);
	tmp_ps_A_trn_ = tmp_ps_A_trn__;
	tmp_ps_B_trn_ = tmp_ps_B_trn__;
	for (ncol_X=0;ncol_X<n_col_X_256;ncol_X++){
#ifdef _FMA
	  ps_acc = _mm256_fmadd_ps(*(tmp_ps_A_trn_++),*(tmp_ps_B_trn_++),ps_acc);
#endif /* _FMA */
	  /* for (ncol_X=0;ncol_X<n_col_X_256;ncol_X++){ } */}
	tmp_f = tmp_f_[0] + tmp_f_[1] + tmp_f_[2] + tmp_f_[3] + tmp_f_[4] + tmp_f_[5] + tmp_f_[6] + tmp_f_[7];
	f_C__[nrow_A + nrow_B*n_row_A] = tmp_f;
	tmp_ps_A_trn__ += n_col_X_256;
	/* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
      tmp_ps_B_trn__ += n_col_X_256;
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ } */}
    /* if (f_C__!=NULL){ } */}
}

void dp_ps_mult_immintrin_avx(int n_row_A,int n_col_X,__m256 *ps_A_trn__,int n_row_B,__m256 *ps_B_trn__,float **f_C_p_)
{
  int n_col_X_rup = rup(n_col_X,8);
  int n_col_X_256 = n_col_X_rup/8; //%<-- 8 floats per __m256. ;
  int nrow_A=0,nrow_B=0,ncol_X=0;
  __m256 *tmp_ps_A_trn__=NULL;
  __m256 *tmp_ps_B_trn__=NULL;
  __m256 *tmp_ps_A_trn_=NULL;
  __m256 *tmp_ps_B_trn_=NULL;
  __m256 tmp_ps,tmp_ps_A_trn,tmp_ps_B_trn;
  float *tmp_f_ = (float *) &tmp_ps;
  float tmp_f=0; 
  float *f_C__=NULL;
  f_C__=NULL;
  if (f_C_p_!=NULL){
    if ( (*f_C_p_)==NULL ){ (*f_C_p_) = (float *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
    f_C__ = *f_C_p_;
    /* if (f_C_p_!=NULL){ } */}
  if (f_C__!=NULL){
    tmp_ps_B_trn__ = ps_B_trn__;
    for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      tmp_ps_A_trn__ = ps_A_trn__;
      for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
	tmp_f=0;
	tmp_ps_A_trn_ = tmp_ps_A_trn__;
	tmp_ps_B_trn_ = tmp_ps_B_trn__;
	for (ncol_X=0;ncol_X<n_col_X_256;ncol_X++){
	  tmp_ps = _mm256_dp_ps(*(tmp_ps_A_trn_++),*(tmp_ps_B_trn_++),0xF1);
	  tmp_f += tmp_f_[0] + tmp_f_[4];
	  /* for (ncol_X=0;ncol_X<n_col_X_256;ncol_X++){ } */}
	f_C__[nrow_A + nrow_B*n_row_A] = tmp_f;
	tmp_ps_A_trn__ += n_col_X_256;
	/* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
      tmp_ps_B_trn__ += n_col_X_256;
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ } */}
    /* if (f_C__!=NULL){ } */}
}

void dp_ps_immintrin_loadu_wrap(int n_col_X,float *f_A_,float *f_B_,float *f_C_)
{
  /* does not assume alignment */
#ifdef _FMA
  dp_ps_immintrin_loadu_fma(n_col_X,f_A_,f_B_,f_C_);
#endif /* _FMA */
#ifndef _FMA
  dp_ps_immintrin_loadu_avx(n_col_X,f_A_,f_B_,f_C_);
#endif /* _FMA */
}

void dp_ps_mult_immintrin_loadu_wrap(int n_row_A,int n_col_X,float *f_A_trn__,int n_row_B,float *f_B_trn__,float **f_C_p_)
{
  /* does not assume alignment */
#ifdef _FMA
  dp_ps_mult_immintrin_loadu_fma(n_row_A,n_col_X,f_A_trn__,n_row_B,f_B_trn__,f_C_p_);
#endif /* _FMA */
#ifndef _FMA
  dp_ps_mult_immintrin_loadu_avx(n_row_A,n_col_X,f_A_trn__,n_row_B,f_B_trn__,f_C_p_);
#endif /* _FMA */
}

void dp_ps_mult_immintrin_test()
{
  int verbose=1;
  /* int n_row_A = 1024*1 + 723; */
  /* int n_col_X = 1024*6 + 817; */
  /* int n_row_B = 1024*1 + 511; */
  int n_row_A = 531;
  int n_col_X = 15039;
  int n_col_X_rup = rup(n_col_X,8);
  int n_col_X_256 = n_col_X_rup/8; //%<-- 8 floats per __m256. ;
  int n_row_B = 430;
  int n_row_A_sub = 5;
  int n_col_X_sub = 15;
  int n_row_B_sub = 4;
  int nrow_A=0,ncol_X=0,nrow_B=0;
  unsigned long long int ulli_A_total = (unsigned long long int) n_row_A * (unsigned long long int) n_col_X_rup ;
  unsigned long long int ulli_B_total = (unsigned long long int) n_row_B * (unsigned long long int) n_col_X_rup ;
  unsigned long long int ulli_C_total = (unsigned long long int) n_row_A * (unsigned long long int) n_row_B ;
  unsigned long long int ulli=0;
  unsigned long long int tab=0;
  __m256 *ps_A_trn__ = NULL;
  __m256 *ps_B_trn__ = NULL;
  float *f_A_trn__ = NULL;
  float *f_A_u_trn__ = NULL;
  float *f_B_trn__ = NULL;
  float *f_B_u_trn__ = NULL;
  float *f_C_bf__ = NULL;
  float *f_C_ps__ = NULL;
  float *f_A_sub__ = NULL;
  float *f_B_sub__ = NULL;
  float *f_C_sub__ = NULL;
  float ferror=0;
  /* %%%% */
  int ntick=0,n_tick = 8;
  clock_t t_start_[n_tick];
  clock_t t_final_[n_tick];
  struct timeval d_start_[n_tick];
  struct timeval d_final_[n_tick];
  long l_msec_[n_tick],l_ssec_[n_tick],l_usec_[n_tick];
  double elct_[n_tick],elrt_[n_tick];
  /* %%%% */
  ps_A_trn__ = (__m256 *) _mm_malloc(n_col_X_256*n_row_A*sizeof(__m256),32);
  ps_B_trn__ = (__m256 *) _mm_malloc(n_col_X_256*n_row_B*sizeof(__m256),32);
  f_A_trn__ = (float *) ps_A_trn__;
  f_B_trn__ = (float *) ps_B_trn__;
  f_A_u_trn__ = (float *) malloc(n_col_X*n_row_A*sizeof(float));
  f_B_u_trn__ = (float *) malloc(n_col_X*n_row_B*sizeof(float));
  f_C_bf__ = (float *) malloc(n_row_A*n_row_B*sizeof(float));
  f_C_ps__ = (float *) malloc(n_row_A*n_row_B*sizeof(float));
  f_A_sub__ = (float *) malloc(n_row_A_sub*n_col_X_sub*sizeof(float));
  f_B_sub__ = (float *) malloc(n_row_B_sub*n_col_X_sub*sizeof(float));
  f_C_sub__ = (float *) malloc(n_row_A_sub*n_row_B_sub*sizeof(float));
  tab = (unsigned long long int) n_row_A * (unsigned long long int) n_row_B * (unsigned long long int) n_col_X_rup ;
  local_tic(0,t_start_,d_start_);
  ulli=0;
  for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ for (ncol_X=0;ncol_X<n_col_X;ncol_X++){
      f_A_trn__[ncol_X + nrow_A*n_col_X_rup] = (float)(((int)ulli%7)-3); 
      f_A_u_trn__[ncol_X + nrow_A*n_col_X] = (float)(((int)ulli%7)-3); 
      ulli++;
      /* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ for (ncol_X=0;ncol_X<n_col_X;ncol_X++){ }} */}}
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (ncol_X=0;ncol_X<n_col_X_sub;ncol_X++){ f_A_sub__[nrow_A+ncol_X*n_row_A_sub] = f_A_trn__[ncol_X + nrow_A*n_col_X_rup];}}
  if (verbose>1){
    printf(" %% upper corner of f_A_trn__: \n");
    array_printf(f_A_sub__,"float",n_row_A_sub,n_col_X_sub," % f_A_sub__: ");
    /* if (verbose>1){ } */}
  ulli=0;
  for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ for (ncol_X=0;ncol_X<n_col_X;ncol_X++){
      f_B_trn__[ncol_X + nrow_B*n_col_X_rup] = (float)(((int)ulli%5)-2);
      f_B_u_trn__[ncol_X + nrow_B*n_col_X] = (float)(((int)ulli%5)-2);
      ulli++;
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ for (ncol_X=0;ncol_X<n_col_X;ncol_X++){ }} */}}
  for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ for (ncol_X=0;ncol_X<n_col_X_sub;ncol_X++){ f_B_sub__[nrow_B+ncol_X*n_row_B_sub] = f_B_trn__[ncol_X + nrow_B*n_col_X_rup];}}
  if (verbose>1){
    printf(" %% upper corner of f_B_trn__: \n");
    array_printf(f_B_sub__,"float",n_row_B_sub,n_col_X_sub," % f_B_sub__: ");
    /* if (verbose>1){ } */}
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," initialize: ");
  /* %%%%%%%% */
  local_tic(0,t_start_,d_start_);
  memset(f_C_bf__,0,ulli_C_total*sizeof(float));
  dp_ps_mult_bruteforce(n_row_A,n_col_X,f_A_trn__,n_row_B,f_B_trn__,&f_C_bf__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," dp_ps_mult_bruteforce: ");
  printf("Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/elrt_[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ f_C_sub__[nrow_A+nrow_B*n_row_A_sub] = f_C_bf__[nrow_A+nrow_B*n_row_A];}}
  if (verbose>1){
    printf(" %% upper corner of f_C_bf__: \n");
    array_printf(f_C_sub__,"float",n_row_A_sub,n_row_B_sub," % f_C_bf__: ");
    /* if (verbose>1){ } */}
  ferror = ffnorm(ulli_C_total,f_C_bf__,f_C_bf__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
#ifdef _CBLAS
  local_tic(0,t_start_,d_start_);
  memset(f_C_ps__,0,ulli_C_total*sizeof(float));
  dp_ps_mult_cblas_sgemm(n_row_A,n_col_X,f_A_trn__,n_row_B,f_B_trn__,&f_C_ps__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," dp_ps_mult_cblas_sgemm: ");
  printf("Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/elrt_[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ f_C_sub__[nrow_A+nrow_B*n_row_A_sub] = f_C_ps__[nrow_A+nrow_B*n_row_A];}}
  if (verbose>1){
    printf(" %% upper corner of f_C_ps__: \n");
    array_printf(f_C_sub__,"float",n_row_A_sub,n_row_B_sub," % f_C_ps__: ");
    /* if (verbose>1){ } */}
  ferror = ffnorm(ulli_C_total,f_C_bf__,f_C_ps__);
  printf(" %% ferror %0.16f\n",ferror);
#endif /* _CBLAS */
  /* %%%%%%%% */
#ifdef _CBLAS
  local_tic(0,t_start_,d_start_);
  memset(f_C_ps__,0,ulli_C_total*sizeof(float));
  dp_ps_mult_cblas_sdot(n_row_A,n_col_X,f_A_trn__,n_row_B,f_B_trn__,&f_C_ps__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," dp_ps_mult_cblas_sdot: ");
  printf("Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/elrt_[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ f_C_sub__[nrow_A+nrow_B*n_row_A_sub] = f_C_ps__[nrow_A+nrow_B*n_row_A];}}
  if (verbose>1){
    printf(" %% upper corner of f_C_ps__: \n");
    array_printf(f_C_sub__,"float",n_row_A_sub,n_row_B_sub," % f_C_ps__: ");
    /* if (verbose>1){ } */}
  ferror = ffnorm(ulli_C_total,f_C_bf__,f_C_ps__);
  printf(" %% ferror %0.16f\n",ferror);
#endif /* _CBLAS */
  /* %%%%%%%% */
  local_tic(0,t_start_,d_start_);
  memset(f_C_ps__,0,ulli_C_total*sizeof(float));
  dp_ps_mult_immintrin_avx(n_row_A,n_col_X,ps_A_trn__,n_row_B,ps_B_trn__,&f_C_ps__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," dp_ps_mult_immintrin_avx: ");
  printf("Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/elrt_[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ f_C_sub__[nrow_A+nrow_B*n_row_A_sub] = f_C_ps__[nrow_A+nrow_B*n_row_A];}}
  if (verbose>1){
    printf(" %% upper corner of f_C_ps__: \n");
    array_printf(f_C_sub__,"float",n_row_A_sub,n_row_B_sub," % f_C_ps__: ");
    /* if (verbose>1){ } */}
  ferror = ffnorm(ulli_C_total,f_C_bf__,f_C_ps__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
  local_tic(0,t_start_,d_start_);
  memset(f_C_ps__,0,ulli_C_total*sizeof(float));
  dp_ps_mult_immintrin_fma(n_row_A,n_col_X,ps_A_trn__,n_row_B,ps_B_trn__,&f_C_ps__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," dp_ps_mult_immintrin_fma: ");
  printf("Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/elrt_[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ f_C_sub__[nrow_A+nrow_B*n_row_A_sub] = f_C_ps__[nrow_A+nrow_B*n_row_A];}}
  if (verbose>1){
    printf(" %% upper corner of f_C_ps__: \n");
    array_printf(f_C_sub__,"float",n_row_A_sub,n_row_B_sub," % f_C_ps__: ");
    /* if (verbose>1){ } */}
  ferror = ffnorm(ulli_C_total,f_C_bf__,f_C_ps__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
  local_tic(0,t_start_,d_start_);
  memset(f_C_ps__,0,ulli_C_total*sizeof(float));
  dp_ps_mult_immintrin_fma2(n_row_A,n_col_X,ps_A_trn__,n_row_B,ps_B_trn__,&f_C_ps__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," dp_ps_mult_immintrin_fma2: ");
  printf("Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/elrt_[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ f_C_sub__[nrow_A+nrow_B*n_row_A_sub] = f_C_ps__[nrow_A+nrow_B*n_row_A];}}
  if (verbose>1){
    printf(" %% upper corner of f_C_ps__: \n");
    array_printf(f_C_sub__,"float",n_row_A_sub,n_row_B_sub," % f_C_ps__: ");
    /* if (verbose>1){ } */}
  ferror = ffnorm(ulli_C_total,f_C_bf__,f_C_ps__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
  local_tic(0,t_start_,d_start_);
  memset(f_C_ps__,0,ulli_C_total*sizeof(float));
  dp_ps_mult_immintrin_load1_fma(n_row_A,n_col_X,f_A_trn__,n_row_B,f_B_trn__,&f_C_ps__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," dp_ps_mult_immintrin_load1_fma: ");
  printf("Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/elrt_[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ f_C_sub__[nrow_A+nrow_B*n_row_A_sub] = f_C_ps__[nrow_A+nrow_B*n_row_A];}}
  if (verbose>1){
    printf(" %% upper corner of f_C_ps__: \n");
    array_printf(f_C_sub__,"float",n_row_A_sub,n_row_B_sub," % f_C_ps__: ");
    /* if (verbose>1){ } */}
  ferror = ffnorm(ulli_C_total,f_C_bf__,f_C_ps__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
  local_tic(0,t_start_,d_start_);
  memset(f_C_ps__,0,ulli_C_total*sizeof(float));
  dp_ps_mult_immintrin_loadu_fma(n_row_A,n_col_X,f_A_u_trn__,n_row_B,f_B_u_trn__,&f_C_ps__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," dp_ps_mult_immintrin_loadu_fma: ");
  printf("Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X/elrt_[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ f_C_sub__[nrow_A+nrow_B*n_row_A_sub] = f_C_ps__[nrow_A+nrow_B*n_row_A];}}
  if (verbose>1){
    printf(" %% upper corner of f_C_ps__: \n");
    array_printf(f_C_sub__,"float",n_row_A_sub,n_row_B_sub," % f_C_ps__: ");
    /* if (verbose>1){ } */}
  ferror = ffnorm(ulli_C_total,f_C_bf__,f_C_ps__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
  local_tic(0,t_start_,d_start_);
  memset(f_C_ps__,0,ulli_C_total*sizeof(float));
  dp_ps_mult_immintrin_loadu_avx(n_row_A,n_col_X,f_A_u_trn__,n_row_B,f_B_u_trn__,&f_C_ps__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," dp_ps_mult_immintrin_loadu_avx: ");
  printf("Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X/elrt_[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ f_C_sub__[nrow_A+nrow_B*n_row_A_sub] = f_C_ps__[nrow_A+nrow_B*n_row_A];}}
  if (verbose>1){
    printf(" %% upper corner of f_C_ps__: \n");
    array_printf(f_C_sub__,"float",n_row_A_sub,n_row_B_sub," % f_C_ps__: ");
    /* if (verbose>1){ } */}
  ferror = ffnorm(ulli_C_total,f_C_bf__,f_C_ps__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
  local_tic(0,t_start_,d_start_);
  memset(f_C_ps__,0,ulli_C_total*sizeof(float));
  dp_ps_mult_immintrin_loadu_wrap(n_row_A,n_col_X,f_A_u_trn__,n_row_B,f_B_u_trn__,&f_C_ps__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," dp_ps_mult_immintrin_loadu_wrap: ");
  printf("Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X/elrt_[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ f_C_sub__[nrow_A+nrow_B*n_row_A_sub] = f_C_ps__[nrow_A+nrow_B*n_row_A];}}
  if (verbose>1){
    printf(" %% upper corner of f_C_ps__: \n");
    array_printf(f_C_sub__,"float",n_row_A_sub,n_row_B_sub," % f_C_ps__: ");
    /* if (verbose>1){ } */}
  ferror = ffnorm(ulli_C_total,f_C_bf__,f_C_ps__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
  _mm_free(ps_A_trn__); ps_A_trn__ = NULL;
  _mm_free(ps_B_trn__); ps_B_trn__ = NULL;
  free(f_A_u_trn__); f_A_u_trn__=NULL;
  free(f_B_u_trn__); f_B_u_trn__=NULL;
  free(f_C_bf__); f_C_bf__=NULL;
  free(f_C_ps__); f_C_ps__=NULL;
  free(f_A_sub__); f_A_sub__=NULL;
  free(f_B_sub__); f_B_sub__=NULL;
  free(f_C_sub__); f_C_sub__=NULL;
  //wkspace_printf();
}
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void hp_interleave_mult_bruteforce(int n_row_A,int n_col_X,float complex *c_A_trn__,int n_row_B,float complex *c_B_trn__,float complex **c_C_p_)
{
  /* unaligned */
  int nrow_A=0,nrow_B=0,ncol_X=0;
  float complex *c_C__=NULL;
  c_C__=NULL;
  if (c_C_p_!=NULL){
    if ( (*c_C_p_)==NULL ){ (*c_C_p_) = (float complex *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float complex));}
    c_C__ = *c_C_p_;
    /* if (c_C_p_!=NULL){ } */}
  if (c_C__!=NULL){
    for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
	for (ncol_X=0;ncol_X<n_col_X;ncol_X++){
	  c_C__[nrow_A + nrow_B*n_row_A] += conjf(c_A_trn__[ncol_X + nrow_A*n_col_X]) * c_B_trn__[ncol_X + nrow_B*n_col_X];
	  /* for (ncol_X=0;ncol_X<n_col_X;ncol_X++){ } */}
	/* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ } */}
    /* if (c_C__!=NULL){ } */}
}

void hp_interleave_mult_cblas_cgemm(int n_row_A,int n_col_X,float complex *c_A_trn__,int n_row_B,float complex *c_B_trn__,float complex **c_C_p_)
{
  /* unaligned */
  float complex *c_C__=NULL;
  float complex calpha = (float complex)1.0;
  float complex cbeta = (float complex)0.0;
  c_C__=NULL;
  if (c_C_p_!=NULL){
    if ( (*c_C_p_)==NULL ){ (*c_C_p_) = (float complex *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float complex));}
    c_C__ = *c_C_p_;
    /* if (c_C_p_!=NULL){ } */}
  if (c_C__!=NULL){
#ifdef _CBLAS
    cblas_cgemm(CblasColMajor,CblasConjTrans,CblasNoTrans,n_row_A,n_row_B,n_col_X,&calpha,c_A_trn__,n_col_X,c_B_trn__,n_col_X,&cbeta,c_C__,n_row_A);
#endif /* _CBLAS */
#ifndef _CBLAS
    printf(" %% Warning, _CBLAS not defined in hp_interleave_mult_cblas_cgemm, using bruteforce instead.\n");
    hp_interleave_mult_bruteforce(n_row_A,n_col_X,c_A_trn__,n_row_B,c_B_trn__,c_C_p_);
#endif /* _CBLAS */
  /* if (c_C__!=NULL){ } */}
}

void hp_segregated_mult_bruteforce(int n_row_A,int n_col_X,float *f_AR_trn__,float *f_AI_trn__,int n_row_B,float *f_BR_trn__,float *f_BI_trn__,float complex **c_C_p_)
{
  /* unaligned */
  int nrow_A=0,nrow_B=0,ncol_X=0;
  float complex *c_C__=NULL;
  c_C__=NULL;
  if (c_C_p_!=NULL){
    if ( (*c_C_p_)==NULL ){ (*c_C_p_) = (float complex *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float complex));}
    c_C__ = *c_C_p_;
    /* if (c_C_p_!=NULL){ } */}
  if (c_C__!=NULL){
    for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
	c_C__[nrow_A + nrow_B*n_row_A] = (float complex) 0.0;
	for (ncol_X=0;ncol_X<n_col_X;ncol_X++){
	  c_C__[nrow_A + nrow_B*n_row_A] +=
	    ((float complex)f_AR_trn__[ncol_X + nrow_A*n_col_X] - _Complex_I * (float complex)f_AI_trn__[ncol_X + nrow_A*n_col_X]) *
	    ((float complex)f_BR_trn__[ncol_X + nrow_B*n_col_X] + _Complex_I * (float complex)f_BI_trn__[ncol_X + nrow_B*n_col_X])
	    ;
	  /* for (ncol_X=0;ncol_X<n_col_X;ncol_X++){ } */}
	/* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ } */}
    /* if (c_C__!=NULL){ } */}
}

void hp_segregated_to_interleaved_mult_immintrin_loadu_fma(int n_row_A,int n_col_X,float *f_AR_trn__,float *f_AI_trn__,int n_row_B,float *f_BR_trn__,float *f_BI_trn__,float complex **c_C_p_)
{
  /* does not assume alignment */
  int nrow_A=0,nrow_B=0;
  float *tmp_f_AR_point0_=NULL;
  float *tmp_f_AI_point0_=NULL;
  float *tmp_f_BR_point0_=NULL;
  float *tmp_f_BI_point0_=NULL;
  float complex *c_C__=NULL;
  float f_ARBR=0.0,f_ARBI=0.0,f_AIBR=0.0,f_AIBI=0.0;
  int na=0;
  c_C__=NULL;
  if (c_C_p_!=NULL){
    if ( (*c_C_p_)==NULL ){ (*c_C_p_) = (float complex *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float complex));}
    c_C__ = *c_C_p_;
    /* if (c_C_p_!=NULL){ } */}
  if (c_C__!=NULL){
    na=0;
    for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
	tmp_f_AR_point0_ = &(f_AR_trn__[nrow_A*n_col_X]);
	tmp_f_AI_point0_ = &(f_AI_trn__[nrow_A*n_col_X]);
	tmp_f_BR_point0_ = &(f_BR_trn__[nrow_B*n_col_X]);
	tmp_f_BI_point0_ = &(f_BI_trn__[nrow_B*n_col_X]);
	dp_ps_immintrin_loadu_fma(n_col_X,tmp_f_AR_point0_,tmp_f_BR_point0_,&f_ARBR);
	dp_ps_immintrin_loadu_fma(n_col_X,tmp_f_AR_point0_,tmp_f_BI_point0_,&f_ARBI);
	dp_ps_immintrin_loadu_fma(n_col_X,tmp_f_AI_point0_,tmp_f_BI_point0_,&f_AIBI);
	dp_ps_immintrin_loadu_fma(n_col_X,tmp_f_AI_point0_,tmp_f_BR_point0_,&f_AIBR);
	c_C__[na] = (float complex)(f_ARBR + f_AIBI) + _Complex_I * (float complex)(f_ARBI - f_AIBR);
	na += 1;
	/* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ } */}
    /* if (c_C__!=NULL){ } */}
}

void hp_segregated_to_segregated_mult_immintrin_loadu_fma(int n_row_A,int n_col_X,float *f_AR_trn__,float *f_AI_trn__,int n_row_B,float *f_BR_trn__,float *f_BI_trn__,float **f_CR_p_,float **f_CI_p_)
{
  /* does not assume alignment */
  int nrow_A=0,nrow_B=0;
  float *tmp_f_AR_point0_=NULL;
  float *tmp_f_AI_point0_=NULL;
  float *tmp_f_BR_point0_=NULL;
  float *tmp_f_BI_point0_=NULL;
  float *f_CR__=NULL;
  float *f_CI__=NULL;
  float f_ARBR=0.0,f_ARBI=0.0,f_AIBR=0.0,f_AIBI=0.0;
  int na=0;
  f_CR__=NULL;
  if (f_CR_p_!=NULL){
    if ( (*f_CR_p_)==NULL ){ (*f_CR_p_) = (float *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
    f_CR__ = *f_CR_p_;
    /* if (f_CR_p_!=NULL){ } */}
  f_CI__=NULL;
  if (f_CI_p_!=NULL){
    if ( (*f_CI_p_)==NULL ){ (*f_CI_p_) = (float *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
    f_CI__ = *f_CI_p_;
    /* if (f_CI_p_!=NULL){ } */}
  if ((f_CR__!=NULL) && (f_CI__!=NULL)){
    na=0;
    for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
	tmp_f_AR_point0_ = &(f_AR_trn__[nrow_A*n_col_X]);
	tmp_f_AI_point0_ = &(f_AI_trn__[nrow_A*n_col_X]);
	tmp_f_BR_point0_ = &(f_BR_trn__[nrow_B*n_col_X]);
	tmp_f_BI_point0_ = &(f_BI_trn__[nrow_B*n_col_X]);
	dp_ps_immintrin_loadu_fma(n_col_X,tmp_f_AR_point0_,tmp_f_BR_point0_,&f_ARBR);
	dp_ps_immintrin_loadu_fma(n_col_X,tmp_f_AR_point0_,tmp_f_BI_point0_,&f_ARBI);
	dp_ps_immintrin_loadu_fma(n_col_X,tmp_f_AI_point0_,tmp_f_BI_point0_,&f_AIBI);
	dp_ps_immintrin_loadu_fma(n_col_X,tmp_f_AI_point0_,tmp_f_BR_point0_,&f_AIBR);
	f_CR__[na] = (f_ARBR + f_AIBI);
	f_CI__[na] = (f_ARBI - f_AIBR);
	na += 1;
	/* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ } */}
    /* if ((f_CR__!=NULL) && (f_CI__!=NULL)){ } */}
}

void hp_segregated_to_segregated_mult_immintrin_load1_fma(int n_row_A,int n_col_X,float *f_AR_trn__,float *f_AI_trn__,int n_row_B,float *f_BR_trn__,float *f_BI_trn__,float **f_CR_p_,float **f_CI_p_)
{
  /* assumes alignment */
  int n_col_X_rup = rup(n_col_X,8);
  int n_col_X_256 = n_col_X_rup/8; //%<-- 8 floats per __m256. ;
  int nrow_A=0,nrow_B=0,ncol_X=0;
  float *tmp_f_AR_point0_=NULL;
  float *tmp_f_AI_point0_=NULL;
  float *tmp_f_BR_point0_=NULL;
  float *tmp_f_BI_point0_=NULL;
  __m256 ps_AR0,ps_BR0;
  __m256 ps_AI0,ps_BI0;
  __m256 ps_acARBR0;
  __m256 ps_acARBI0;
  __m256 ps_acAIBR0;
  __m256 ps_acAIBI0;
  float *tmp_f_acARBR0_ = (float *) &ps_acARBR0;
  float *tmp_f_acARBI0_ = (float *) &ps_acARBI0;
  float *tmp_f_acAIBR0_ = (float *) &ps_acAIBR0;
  float *tmp_f_acAIBI0_ = (float *) &ps_acAIBI0;
  float tmp_f_acARBR0=0;
  float tmp_f_acARBI0=0;
  float tmp_f_acAIBR0=0;
  float tmp_f_acAIBI0=0;
  float *f_CR__=NULL;
  float *f_CI__=NULL;
  int na=0;
  f_CR__=NULL;
  if (f_CR_p_!=NULL){
    if ( (*f_CR_p_)==NULL ){ (*f_CR_p_) = (float *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
    f_CR__ = *f_CR_p_;
    /* if (f_CR_p_!=NULL){ } */}
  f_CI__=NULL;
  if (f_CI_p_!=NULL){
    if ( (*f_CI_p_)==NULL ){ (*f_CI_p_) = (float *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
    f_CI__ = *f_CI_p_;
    /* if (f_CI_p_!=NULL){ } */}
  if ((f_CR__!=NULL) && (f_CI__!=NULL)){
    na=0;
    for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
	tmp_f_AR_point0_ = &(f_AR_trn__[nrow_A*n_col_X_rup]);
	tmp_f_AI_point0_ = &(f_AI_trn__[nrow_A*n_col_X_rup]);
	tmp_f_BR_point0_ = &(f_BR_trn__[nrow_B*n_col_X_rup]);
	tmp_f_BI_point0_ = &(f_BI_trn__[nrow_B*n_col_X_rup]);
	ps_acARBR0 = _mm256_set1_ps((float)0.0);
	ps_acARBI0 = _mm256_set1_ps((float)0.0);
	ps_acAIBR0 = _mm256_set1_ps((float)0.0);
	ps_acAIBI0 = _mm256_set1_ps((float)0.0);
	for (ncol_X=0;ncol_X<n_col_X_256;ncol_X++){
	  ps_AR0 = _mm256_load_ps(tmp_f_AR_point0_);
	  ps_BR0 = _mm256_load_ps(tmp_f_BR_point0_);
	  ps_AI0 = _mm256_load_ps(tmp_f_AI_point0_);
	  ps_BI0 = _mm256_load_ps(tmp_f_BI_point0_);
#ifdef _FMA
	  ps_acARBR0 = _mm256_fmadd_ps(ps_AR0,ps_BR0,ps_acARBR0);
	  ps_acARBI0 = _mm256_fmadd_ps(ps_AR0,ps_BI0,ps_acARBI0);
	  ps_acAIBR0 = _mm256_fmadd_ps(ps_AI0,ps_BR0,ps_acAIBR0);
	  ps_acAIBI0 = _mm256_fmadd_ps(ps_AI0,ps_BI0,ps_acAIBI0);
#endif /* _FMA */
#ifdef _FMA
	  tmp_f_AR_point0_+=8;tmp_f_BR_point0_+=8;
	  tmp_f_AI_point0_+=8;tmp_f_BI_point0_+=8;
#endif /* _FMA */
	  /* for (ncol_X=0;ncol_X<n_col_X_256;ncol_X++){ } */}
	tmp_f_acARBR0 = tmp_f_acARBR0_[0] + tmp_f_acARBR0_[1] + tmp_f_acARBR0_[2] + tmp_f_acARBR0_[3] + tmp_f_acARBR0_[4] + tmp_f_acARBR0_[5] + tmp_f_acARBR0_[6] + tmp_f_acARBR0_[7];
	tmp_f_acARBI0 = tmp_f_acARBI0_[0] + tmp_f_acARBI0_[1] + tmp_f_acARBI0_[2] + tmp_f_acARBI0_[3] + tmp_f_acARBI0_[4] + tmp_f_acARBI0_[5] + tmp_f_acARBI0_[6] + tmp_f_acARBI0_[7];
	tmp_f_acAIBR0 = tmp_f_acAIBR0_[0] + tmp_f_acAIBR0_[1] + tmp_f_acAIBR0_[2] + tmp_f_acAIBR0_[3] + tmp_f_acAIBR0_[4] + tmp_f_acAIBR0_[5] + tmp_f_acAIBR0_[6] + tmp_f_acAIBR0_[7];
	tmp_f_acAIBI0 = tmp_f_acAIBI0_[0] + tmp_f_acAIBI0_[1] + tmp_f_acAIBI0_[2] + tmp_f_acAIBI0_[3] + tmp_f_acAIBI0_[4] + tmp_f_acAIBI0_[5] + tmp_f_acAIBI0_[6] + tmp_f_acAIBI0_[7];
	f_CR__[na] = tmp_f_acARBR0 + tmp_f_acAIBI0;
	f_CI__[na] = tmp_f_acARBI0 - tmp_f_acAIBR0;
	na += 1;
	/* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ } */}
    /* if ((f_CR__!=NULL) && (f_CI__!=NULL)){ } */}
}

void hp_ps_mult_immintrin_test()
{
  int verbose = 1;
  int n_row_A = 531;
  int n_col_X = 15039;
  int n_col_X_rup = rup(n_col_X,8);
  int n_col_X_256 = n_col_X_rup/8; //%<-- 8 floats per __m256. ;
  int n_row_B = 430;
  int n_row_A_sub = 5;
  int n_col_X_sub = 15;
  int n_row_B_sub = 4;
  int nrow_A=0,ncol_X=0,nrow_B=0;
  unsigned long long int ulli_A_total = (unsigned long long int) n_row_A * (unsigned long long int) n_col_X_rup ;
  unsigned long long int ulli_B_total = (unsigned long long int) n_row_B * (unsigned long long int) n_col_X_rup ;
  unsigned long long int ulli_C_total = (unsigned long long int) n_row_A * (unsigned long long int) n_row_B ;
  unsigned long long int ulli=0;
  unsigned long long int tab=0;
  __m256 *ps_AR_trn__ = NULL;
  __m256 *ps_AI_trn__ = NULL;
  __m256 *ps_BR_trn__ = NULL;
  __m256 *ps_BI_trn__ = NULL;
  float *f_AR_trn__ = NULL;
  float *f_AI_trn__ = NULL;
  float *f_AR_u_trn__ = NULL;
  float *f_AI_u_trn__ = NULL;
  float *f_BR_trn__ = NULL;
  float *f_BI_trn__ = NULL;
  float *f_BR_u_trn__ = NULL;
  float *f_BI_u_trn__ = NULL;
  float *f_CR_bf__ = NULL;
  float *f_CI_bf__ = NULL;
  float *f_CR_al__ = NULL;
  float *f_CI_al__ = NULL;
  float *f_AR_sub__ = NULL;
  float *f_AI_sub__ = NULL;
  float *f_BR_sub__ = NULL;
  float *f_BI_sub__ = NULL;
  float *f_CR_sub__ = NULL;
  float *f_CI_sub__ = NULL;
  float complex *c_A_u_trn__ = NULL;
  float complex *c_B_u_trn__ = NULL;
  float complex *c_C_bf__ = NULL;
  float complex *c_C_al__ = NULL;
  float complex *c_A_sub__ = NULL;
  float complex *c_B_sub__ = NULL;
  float complex *c_C_sub__ = NULL;
  float ferror=0;
  /* %%%% */
  int ntick=0,n_tick = 8;
  clock_t t_start_[n_tick];
  clock_t t_final_[n_tick];
  struct timeval d_start_[n_tick];
  struct timeval d_final_[n_tick];
  long l_msec_[n_tick],l_ssec_[n_tick],l_usec_[n_tick];
  double elct_[n_tick],elrt_[n_tick];
  /* %%%% */
  ps_AR_trn__ = (__m256 *) _mm_malloc(n_col_X_256*n_row_A*sizeof(__m256),32);
  ps_AI_trn__ = (__m256 *) _mm_malloc(n_col_X_256*n_row_A*sizeof(__m256),32);
  ps_BR_trn__ = (__m256 *) _mm_malloc(n_col_X_256*n_row_B*sizeof(__m256),32);
  ps_BI_trn__ = (__m256 *) _mm_malloc(n_col_X_256*n_row_B*sizeof(__m256),32);
  f_AR_trn__ = (float *) ps_AR_trn__;
  f_AI_trn__ = (float *) ps_AI_trn__;
  f_BR_trn__ = (float *) ps_BR_trn__;
  f_BI_trn__ = (float *) ps_BI_trn__;
  f_AR_u_trn__ = (float *) malloc(n_col_X*n_row_A*sizeof(float));
  f_AI_u_trn__ = (float *) malloc(n_col_X*n_row_A*sizeof(float));
  f_BR_u_trn__ = (float *) malloc(n_col_X*n_row_B*sizeof(float));
  f_BI_u_trn__ = (float *) malloc(n_col_X*n_row_B*sizeof(float));
  f_CR_bf__ = (float *) malloc(n_row_A*n_row_B*sizeof(float));
  f_CI_bf__ = (float *) malloc(n_row_A*n_row_B*sizeof(float));
  f_CR_al__ = (float *) malloc(n_row_A*n_row_B*sizeof(float));
  f_CI_al__ = (float *) malloc(n_row_A*n_row_B*sizeof(float));
  f_AR_sub__ = (float *) malloc(n_row_A_sub*n_col_X_sub*sizeof(float));
  f_AI_sub__ = (float *) malloc(n_row_A_sub*n_col_X_sub*sizeof(float));
  f_BR_sub__ = (float *) malloc(n_row_B_sub*n_col_X_sub*sizeof(float));
  f_BI_sub__ = (float *) malloc(n_row_B_sub*n_col_X_sub*sizeof(float));
  f_CR_sub__ = (float *) malloc(n_row_A_sub*n_row_B_sub*sizeof(float));
  f_CI_sub__ = (float *) malloc(n_row_A_sub*n_row_B_sub*sizeof(float));
  c_A_u_trn__ = (float complex *) malloc(n_col_X*n_row_A*sizeof(float complex));
  c_B_u_trn__ = (float complex *) malloc(n_col_X*n_row_B*sizeof(float complex));
  c_C_bf__ = (float complex *) malloc(n_row_A*n_row_B*sizeof(float complex));
  c_C_al__ = (float complex *) malloc(n_row_A*n_row_B*sizeof(float complex));
  c_A_sub__ = (float complex *) malloc(n_row_A_sub*n_col_X_sub*sizeof(float complex));
  c_B_sub__ = (float complex *) malloc(n_row_B_sub*n_col_X_sub*sizeof(float complex));
  c_C_sub__ = (float complex *) malloc(n_row_A_sub*n_row_B_sub*sizeof(float complex));
  tab = (unsigned long long int) n_row_A * (unsigned long long int) n_row_B * (unsigned long long int) n_col_X_rup;
  local_tic(0,t_start_,d_start_);
  ulli=0;
  for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ for (ncol_X=0;ncol_X<n_col_X;ncol_X++){
      f_AR_trn__[ncol_X + nrow_A*n_col_X_rup] = (float)(((int)ulli%7)-3);
      f_AI_trn__[ncol_X + nrow_A*n_col_X_rup] = (float)(((int)ulli%7)-2);
      f_AR_u_trn__[ncol_X + nrow_A*n_col_X] = (float)(((int)ulli%7)-3);
      f_AI_u_trn__[ncol_X + nrow_A*n_col_X] = (float)(((int)ulli%7)-2);
      c_A_u_trn__[ncol_X + nrow_A*n_col_X] = (float complex) f_AR_u_trn__[ncol_X + nrow_A*n_col_X] + _Complex_I * (float complex) f_AI_u_trn__[ncol_X + nrow_A*n_col_X];
      ulli++;
      /* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ for (ncol_X=0;ncol_X<n_col_X;ncol_X++){ }} */}}  
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (ncol_X=0;ncol_X<n_col_X_sub;ncol_X++){ 
      f_AR_sub__[nrow_A+ncol_X*n_row_A_sub] = f_AR_trn__[ncol_X + nrow_A*n_col_X_rup];
      f_AI_sub__[nrow_A+ncol_X*n_row_A_sub] = f_AI_trn__[ncol_X + nrow_A*n_col_X_rup];
      c_A_sub__[nrow_A+ncol_X*n_row_A_sub] = (float complex) f_AR_sub__[nrow_A+ncol_X*n_row_A_sub] + _Complex_I * (float complex) f_AI_sub__[nrow_A+ncol_X*n_row_A_sub];
      /* for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (ncol_X=0;ncol_X<n_col_X_sub;ncol_X++){ }} */}}
  if (verbose>1){
    printf(" %% upper corner of f_AR_trn__: \n");
    array_printf(f_AR_sub__,"float",n_row_A_sub,n_col_X_sub," % f_AR_sub__: ");
    printf(" %% upper corner of f_AI_trn__: \n");
    array_printf(f_AI_sub__,"float",n_row_A_sub,n_col_X_sub," % f_AI_sub__: ");
    printf(" %% upper corner of c_A_trn__: \n");
    array_printf(c_A_sub__,"float complex",n_row_A_sub,n_col_X_sub," % c_A_sub__: ");
    /* if (verbose>1){ } */}
  ulli=0;
  for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ for (ncol_X=0;ncol_X<n_col_X;ncol_X++){
      f_BR_trn__[ncol_X + nrow_B*n_col_X_rup] = (float)(((int)ulli%11)-2);
      f_BI_trn__[ncol_X + nrow_B*n_col_X_rup] = (float)(((int)ulli%11)-1);
      f_BR_u_trn__[ncol_X + nrow_B*n_col_X] = (float)(((int)ulli%11)-2);
      f_BI_u_trn__[ncol_X + nrow_B*n_col_X] = (float)(((int)ulli%11)-1);
      c_B_u_trn__[ncol_X + nrow_B*n_col_X] = (float complex) f_BR_u_trn__[ncol_X + nrow_B*n_col_X] + _Complex_I * (float complex) f_BI_u_trn__[ncol_X + nrow_B*n_col_X];
      ulli++;
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ for (ncol_X=0;ncol_X<n_col_X;ncol_X++){ }} */}}
  if (verbose>1){
    for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ for (ncol_X=0;ncol_X<n_col_X_sub;ncol_X++){
	f_BR_sub__[nrow_B+ncol_X*n_row_B_sub] = f_BR_trn__[ncol_X + nrow_B*n_col_X_rup];
	f_BI_sub__[nrow_B+ncol_X*n_row_B_sub] = f_BI_trn__[ncol_X + nrow_B*n_col_X_rup];
	c_B_sub__[nrow_B+ncol_X*n_row_B_sub] = (float complex) f_BR_sub__[nrow_B+ncol_X*n_row_B_sub] + _Complex_I * (float complex) f_BI_sub__[nrow_B+ncol_X*n_row_B_sub];
	/* for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ for (ncol_X=0;ncol_X<n_col_X_sub;ncol_X++){ }} */}}
    printf(" %% upper corner of f_BR_trn__: \n");
    array_printf(f_BR_sub__,"float",n_row_B_sub,n_col_X_sub," % f_BR_sub__: ");
    printf(" %% upper corner of f_BI_trn__: \n");
    array_printf(f_BI_sub__,"float",n_row_B_sub,n_col_X_sub," % f_BI_sub__: ");
    printf(" %% upper corner of c_B_trn__: \n");
    array_printf(c_B_sub__,"float complex",n_row_B_sub,n_col_X_sub," % c_B_sub__: ");
    /* if (verbose>1){ } */}
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," initialize: ");
  /* %%%%%%%% */
  local_tic(0,t_start_,d_start_);
  memset(c_C_bf__,0,ulli_C_total*sizeof(float complex));
  hp_interleave_mult_bruteforce(n_row_A,n_col_X,c_A_u_trn__,n_row_B,c_B_u_trn__,&c_C_bf__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," hp_ps_mult_bruteforce: ");
  printf(" %% Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/elrt_[0]/1e9);
  if (verbose>1){
    for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){
	c_C_sub__[nrow_A+nrow_B*n_row_A_sub] = c_C_bf__[nrow_A+nrow_B*n_row_A];
	/* for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ }} */}}
    printf(" %% upper corner of c_C_bf__: \n");
    array_printf(c_C_sub__,"float complex",n_row_A_sub,n_row_B_sub," % c_C_bf__: ");
    /* if (verbose>1){ } */}
  ferror = cfnorm(ulli_C_total,c_C_bf__,c_C_bf__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
  local_tic(0,t_start_,d_start_);
  memset(c_C_al__,0,ulli_C_total*sizeof(float complex));
  hp_segregated_mult_bruteforce(n_row_A,n_col_X,f_AR_u_trn__,f_AI_u_trn__,n_row_B,f_BR_u_trn__,f_BI_u_trn__,&c_C_al__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," hp_segregated_mult_bruteforce: ");
  printf(" %% Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/elrt_[0]/1e9);
  if (verbose>1){
    for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){
	c_C_sub__[nrow_A+nrow_B*n_row_A_sub] = c_C_al__[nrow_A+nrow_B*n_row_A];
	/* for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ }} */}}
    printf(" %% upper corner of c_C_al__: \n");
    array_printf(c_C_sub__,"float complex",n_row_A_sub,n_row_B_sub," % c_C_al__: ");
    /* if (verbose>1){ } */}
  ferror = cfnorm(ulli_C_total,c_C_bf__,c_C_al__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
  local_tic(0,t_start_,d_start_);
  memset(c_C_al__,0,ulli_C_total*sizeof(float complex));
  hp_segregated_to_interleaved_mult_immintrin_loadu_fma(n_row_A,n_col_X,f_AR_u_trn__,f_AI_u_trn__,n_row_B,f_BR_u_trn__,f_BI_u_trn__,&c_C_al__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," hp_segregated_to_interleaved_mult_immintrin_loadu_fma: ");
  printf(" %% Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/elrt_[0]/1e9);
  if (verbose>1){
    for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){
	c_C_sub__[nrow_A+nrow_B*n_row_A_sub] = c_C_al__[nrow_A+nrow_B*n_row_A];
	/* for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ }} */}}
    printf(" %% upper corner of c_C_al__: \n");
    array_printf(c_C_sub__,"float complex",n_row_A_sub,n_row_B_sub," % c_C_al__: ");
    /* if (verbose>1){ } */}
  ferror = cfnorm(ulli_C_total,c_C_bf__,c_C_al__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
  local_tic(0,t_start_,d_start_);
  memset(f_CR_al__,0,ulli_C_total*sizeof(float));
  memset(f_CI_al__,0,ulli_C_total*sizeof(float));
  hp_segregated_to_segregated_mult_immintrin_loadu_fma(n_row_A,n_col_X,f_AR_u_trn__,f_AI_u_trn__,n_row_B,f_BR_u_trn__,f_BI_u_trn__,&f_CR_al__,&f_CI_al__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," hp_segregated_to_segregated_mult_immintrin_loadu_fma: ");
  printf(" %% Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/elrt_[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      c_C_al__[nrow_A+nrow_B*n_row_A] = f_CR_al__[nrow_A+nrow_B*n_row_A] + _Complex_I * (float complex) f_CI_al__[nrow_A+nrow_B*n_row_A];
      /* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ }} */}}
  if (verbose>1){
    for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){
	c_C_sub__[nrow_A+nrow_B*n_row_A_sub] = c_C_al__[nrow_A+nrow_B*n_row_A];
	/* for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ }} */}}
    printf(" %% upper corner of c_C_al__: \n");
    array_printf(c_C_sub__,"float complex",n_row_A_sub,n_row_B_sub," % c_C_al__: ");
    /* if (verbose>1){ } */}
  ferror = cfnorm(ulli_C_total,c_C_bf__,c_C_al__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
  local_tic(0,t_start_,d_start_);
  memset(f_CR_al__,0,ulli_C_total*sizeof(float));
  memset(f_CI_al__,0,ulli_C_total*sizeof(float));
  hp_segregated_to_segregated_mult_immintrin_load1_fma(n_row_A,n_col_X,f_AR_trn__,f_AI_trn__,n_row_B,f_BR_trn__,f_BI_trn__,&f_CR_al__,&f_CI_al__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," hp_segregated_to_segregated_mult_immintrin_load1_fma: ");
  printf(" %% Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/elrt_[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      c_C_al__[nrow_A+nrow_B*n_row_A] = f_CR_al__[nrow_A+nrow_B*n_row_A] + _Complex_I * (float complex) f_CI_al__[nrow_A+nrow_B*n_row_A];
      /* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ }} */}}
  if (verbose>1){
    for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){
	c_C_sub__[nrow_A+nrow_B*n_row_A_sub] = c_C_al__[nrow_A+nrow_B*n_row_A];
	/* for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ }} */}}
    printf(" %% upper corner of c_C_al__: \n");
    array_printf(c_C_sub__,"float complex",n_row_A_sub,n_row_B_sub," % c_C_al__: ");
    /* if (verbose>1){ } */}
  ferror = cfnorm(ulli_C_total,c_C_bf__,c_C_al__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
#ifdef _CBLAS
  local_tic(0,t_start_,d_start_);
  memset(c_C_al__,0,ulli_C_total*sizeof(float complex));
  hp_interleave_mult_cblas_cgemm(n_row_A,n_col_X,c_A_u_trn__,n_row_B,c_B_u_trn__,&c_C_al__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," hp_interleave_mult_cblas_cgemm: ");
  printf(" %% Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/elrt_[0]/1e9);
  if (verbose>1){
    for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){
	c_C_sub__[nrow_A+nrow_B*n_row_A_sub] = c_C_al__[nrow_A+nrow_B*n_row_A];
	/* for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ }} */}}
    printf(" %% upper corner of c_C_al__: \n");
    array_printf(c_C_sub__,"float complex",n_row_A_sub,n_row_B_sub," % c_C_al__: ");
    /* if (verbose>1){ } */}
  ferror = cfnorm(ulli_C_total,c_C_bf__,c_C_al__);
  printf(" %% ferror %0.16f\n",ferror);
#endif /* _CBLAS */
  /* %%%%%%%% */
  _mm_free(ps_AR_trn__); ps_AR_trn__ = NULL;
  _mm_free(ps_AI_trn__); ps_AI_trn__ = NULL;
  _mm_free(ps_BR_trn__); ps_BR_trn__ = NULL;
  _mm_free(ps_BI_trn__); ps_BI_trn__ = NULL;
  free(f_AR_u_trn__); f_AR_u_trn__=NULL;
  free(f_AI_u_trn__); f_AI_u_trn__=NULL;
  free(f_BR_u_trn__); f_BR_u_trn__=NULL;
  free(f_BI_u_trn__); f_BI_u_trn__=NULL;
  free(f_CR_bf__); f_CR_bf__=NULL;
  free(f_CI_bf__); f_CI_bf__=NULL;
  free(f_CR_al__); f_CR_al__=NULL;
  free(f_CI_al__); f_CI_al__=NULL;
  free(f_AR_sub__); f_AR_sub__=NULL;
  free(f_AI_sub__); f_AI_sub__=NULL;
  free(f_BR_sub__); f_BR_sub__=NULL;
  free(f_BI_sub__); f_BI_sub__=NULL;
  free(f_CR_sub__); f_CR_sub__=NULL;
  free(f_CI_sub__); f_CI_sub__=NULL;
  free(c_A_u_trn__); c_A_u_trn__=NULL;
  free(c_B_u_trn__); c_B_u_trn__=NULL;
  free(c_C_bf__); c_C_bf__=NULL;
  free(c_C_al__); c_C_al__=NULL;
  free(c_A_sub__); c_A_sub__=NULL;
  free(c_B_sub__); c_B_sub__=NULL;
  free(c_C_sub__); c_C_sub__=NULL;
  //wkspace_printf();
}
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void cblas_test()
{
  int verbose=1;
  int M=0,N=0,K=0,T=0,nT=0;
  float *sA__=NULL;
  float *sB__=NULL;
  float *sC__=NULL;
  double *dA__=NULL;
  double *dB__=NULL;
  double *dC__=NULL;
  complex *cA__=NULL;
  complex *cB__=NULL;
  complex *cC__=NULL;
  double complex *zA__=NULL;
  double complex *zB__=NULL;
  double complex *zC__=NULL;
  unsigned long long int tab=0;
  /* %%%% */
  float cblas_salpha = (float) 1.0;
  float cblas_sbeta  = (float) 0.0;
  double cblas_dalpha = (double) 1.0;
  double cblas_dbeta  = (double) 0.0;
  complex cblas_calpha = (complex) 1.0;
  complex cblas_cbeta  = (complex) 0.0;
  double complex cblas_zalpha = (double complex) 1.0;
  double complex cblas_zbeta  = (double complex) 0.0;
  /* %%%% */
  int ntick=0,n_tick = 8;
  clock_t t_start_[n_tick];
  clock_t t_final_[n_tick];
  struct timeval d_start_[n_tick];
  struct timeval d_final_[n_tick];
  long l_msec_[n_tick],l_ssec_[n_tick],l_usec_[n_tick];
  double elct_[n_tick],elrt_[n_tick];
  if (verbose){ printf(" %% [entering cblas_test]\n");}
  /* %%%% */
  M = 130; N = 310; K = 160; T = 200;
  sA__ = (float *) calloc(M*K,sizeof(float)); sB__ = (float *) calloc(N*K,sizeof(float)); sC__ = (float *) calloc(M*N,sizeof(float));
  dA__ = (double *) calloc(M*K,sizeof(double)); dB__ = (double *) calloc(N*K,sizeof(double)); dC__ = (double *) calloc(M*N,sizeof(double));
  cA__ = (complex *) calloc(M*K,sizeof(complex)); cB__ = (complex *) calloc(N*K,sizeof(complex)); cC__ = (complex *) calloc(M*N,sizeof(complex));
  zA__ = (double complex *) calloc(M*K,sizeof(double complex)); zB__ = (double complex *) calloc(N*K,sizeof(double complex)); zC__ = (double complex *) calloc(M*N,sizeof(double complex));
  /* %%%% */
  local_tic(0,t_start_,d_start_);
  tab = (unsigned long long int)(M*N*K)*(unsigned long long int)T;
  for (nT=0;nT<T;nT++){
    cblas_sgemm(
		CblasColMajor
		,CblasTrans
		,CblasNoTrans
		,M
		,N
		,K
		,cblas_salpha
		,sA__
		,K
		,sB__
		,K
		,cblas_sbeta
		,sC__
		,M
		);
    /* for (nT=0;nT<T;nT++){ } */}
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose,"sgemm: ");  
  /* %%%% */
  local_tic(0,t_start_,d_start_);
  tab = (unsigned long long int)(M*N*K)*(unsigned long long int)T;
  for (nT=0;nT<T;nT++){
    cblas_dgemm(
		CblasColMajor
		,CblasTrans
		,CblasNoTrans
		,M
		,N
		,K
		,cblas_dalpha
		,dA__
		,K
		,dB__
		,K
		,cblas_dbeta
		,dC__
		,M
		);
    /* for (nT=0;nT<T;nT++){ } */}
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose,"dgemm: ");  
  /* %%%% */
  local_tic(0,t_start_,d_start_);
  tab = (unsigned long long int)(M*N*K)*(unsigned long long int)T;
  for (nT=0;nT<T;nT++){
    cblas_cgemm(
		CblasColMajor
		,CblasConjTrans
		,CblasNoTrans
		,M
		,N
		,K
		,&cblas_calpha
		,cA__
		,K
		,cB__
		,K
		,&cblas_cbeta
		,cC__
		,M
		);
    /* for (nT=0;nT<T;nT++){ } */}
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose,"cgemm: ");  
  /* %%%% */
  local_tic(0,t_start_,d_start_);
  tab = (unsigned long long int)(M*N*K)*(unsigned long long int)T;
  for (nT=0;nT<T;nT++){
    cblas_zgemm(
		CblasColMajor
		,CblasConjTrans
		,CblasNoTrans
		,M
		,N
		,K
		,&cblas_zalpha
		,zA__
		,K
		,zB__
		,K
		,&cblas_zbeta
		,zC__
		,M
		);
    /* for (nT=0;nT<T;nT++){ } */}
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose,"zgemm: ");  
  /* %%%% */
  free(sA__);free(sB__);free(sC__);
  free(dA__);free(dB__);free(dC__);
  free(cA__);free(cB__);free(cC__);
  free(zA__);free(zB__);free(zC__);
  /* %%%% */
  if (verbose){ printf(" %% [finished cblas_test]\n");}  
}

void nhp_segregated_to_segregated_mult_immintrin_load1_fma(int n_row_A,int n_col_X,float *f_AR_trn__,float *f_AI_trn__,int n_row_B,float *f_BR_trn__,float *f_BI_trn__,float **f_CR_p_,float **f_CI_p_)
{
  /* non hermitian product */
  /* assumes alignment */
  int n_col_X_rup = rup(n_col_X,8);
  int n_col_X_256 = n_col_X_rup/8; //%<-- 8 floats per __m256. ;
  int nrow_A=0,nrow_B=0,ncol_X=0;
  float *tmp_f_AR_point0_=NULL;
  float *tmp_f_AI_point0_=NULL;
  float *tmp_f_BR_point0_=NULL;
  float *tmp_f_BI_point0_=NULL;
  __m256 ps_AR0,ps_BR0;
  __m256 ps_AI0,ps_BI0;
  __m256 ps_acARBR0;
  __m256 ps_acARBI0;
  __m256 ps_acAIBR0;
  __m256 ps_acAIBI0;
  float *tmp_f_acARBR0_ = (float *) &ps_acARBR0;
  float *tmp_f_acARBI0_ = (float *) &ps_acARBI0;
  float *tmp_f_acAIBR0_ = (float *) &ps_acAIBR0;
  float *tmp_f_acAIBI0_ = (float *) &ps_acAIBI0;
  float tmp_f_acARBR0=0;
  float tmp_f_acARBI0=0;
  float tmp_f_acAIBR0=0;
  float tmp_f_acAIBI0=0;
  float *f_CR__=NULL;
  float *f_CI__=NULL;
  int na=0;
  f_CR__=NULL;
  if (f_CR_p_!=NULL){
    if ( (*f_CR_p_)==NULL ){ (*f_CR_p_) = (float *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
    f_CR__ = *f_CR_p_;
    /* if (f_CR_p_!=NULL){ } */}
  f_CI__=NULL;
  if (f_CI_p_!=NULL){
    if ( (*f_CI_p_)==NULL ){ (*f_CI_p_) = (float *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
    f_CI__ = *f_CI_p_;
    /* if (f_CI_p_!=NULL){ } */}
  if ((f_CR__!=NULL) && (f_CI__!=NULL)){
    na=0;
    for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
	tmp_f_AR_point0_ = &(f_AR_trn__[nrow_A*n_col_X_rup]);
	tmp_f_AI_point0_ = &(f_AI_trn__[nrow_A*n_col_X_rup]);
	tmp_f_BR_point0_ = &(f_BR_trn__[nrow_B*n_col_X_rup]);
	tmp_f_BI_point0_ = &(f_BI_trn__[nrow_B*n_col_X_rup]);
	ps_acARBR0 = _mm256_set1_ps((float)0.0);
	ps_acARBI0 = _mm256_set1_ps((float)0.0);
	ps_acAIBR0 = _mm256_set1_ps((float)0.0);
	ps_acAIBI0 = _mm256_set1_ps((float)0.0);
	for (ncol_X=0;ncol_X<n_col_X_256;ncol_X++){
	  ps_AR0 = _mm256_load_ps(tmp_f_AR_point0_);
	  ps_BR0 = _mm256_load_ps(tmp_f_BR_point0_);
	  ps_AI0 = _mm256_load_ps(tmp_f_AI_point0_);
	  ps_BI0 = _mm256_load_ps(tmp_f_BI_point0_);
#ifdef _FMA
	  ps_acARBR0 = _mm256_fmadd_ps(ps_AR0,ps_BR0,ps_acARBR0);
	  ps_acARBI0 = _mm256_fmadd_ps(ps_AR0,ps_BI0,ps_acARBI0);
	  ps_acAIBR0 = _mm256_fmadd_ps(ps_AI0,ps_BR0,ps_acAIBR0);
	  ps_acAIBI0 = _mm256_fmadd_ps(ps_AI0,ps_BI0,ps_acAIBI0);
#endif /* _FMA */
#ifdef _FMA
	  tmp_f_AR_point0_+=8;tmp_f_BR_point0_+=8;
	  tmp_f_AI_point0_+=8;tmp_f_BI_point0_+=8;
#endif /* _FMA */
	  /* for (ncol_X=0;ncol_X<n_col_X_256;ncol_X++){ } */}
	tmp_f_acARBR0 = tmp_f_acARBR0_[0] + tmp_f_acARBR0_[1] + tmp_f_acARBR0_[2] + tmp_f_acARBR0_[3] + tmp_f_acARBR0_[4] + tmp_f_acARBR0_[5] + tmp_f_acARBR0_[6] + tmp_f_acARBR0_[7];
	tmp_f_acARBI0 = tmp_f_acARBI0_[0] + tmp_f_acARBI0_[1] + tmp_f_acARBI0_[2] + tmp_f_acARBI0_[3] + tmp_f_acARBI0_[4] + tmp_f_acARBI0_[5] + tmp_f_acARBI0_[6] + tmp_f_acARBI0_[7];
	tmp_f_acAIBR0 = tmp_f_acAIBR0_[0] + tmp_f_acAIBR0_[1] + tmp_f_acAIBR0_[2] + tmp_f_acAIBR0_[3] + tmp_f_acAIBR0_[4] + tmp_f_acAIBR0_[5] + tmp_f_acAIBR0_[6] + tmp_f_acAIBR0_[7];
	tmp_f_acAIBI0 = tmp_f_acAIBI0_[0] + tmp_f_acAIBI0_[1] + tmp_f_acAIBI0_[2] + tmp_f_acAIBI0_[3] + tmp_f_acAIBI0_[4] + tmp_f_acAIBI0_[5] + tmp_f_acAIBI0_[6] + tmp_f_acAIBI0_[7];
	f_CR__[na] = tmp_f_acARBR0 - tmp_f_acAIBI0;
	f_CI__[na] = tmp_f_acARBI0 + tmp_f_acAIBR0;
	na += 1;
	/* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ } */}
    /* if ((f_CR__!=NULL) && (f_CI__!=NULL)){ } */}
}

void nhpr_segregated_to_segregated_mult_immintrin_load1_fma(int n_row_A,int n_col_X,float *f_AR_trn__,float *f_AI_trn__,int n_row_B,float *f_BR_trn__,float *f_BI_trn__,float **f_CR_p_,float **f_CI_p_)
{
  /* non hermitian product, real output only */
  /* assumes alignment */
  int n_col_X_rup = rup(n_col_X,8);
  int n_col_X_256 = n_col_X_rup/8; //%<-- 8 floats per __m256. ;
  int nrow_A=0,nrow_B=0,ncol_X=0;
  float *tmp_f_AR_point0_=NULL;
  float *tmp_f_AI_point0_=NULL;
  float *tmp_f_BR_point0_=NULL;
  float *tmp_f_BI_point0_=NULL;
  __m256 ps_AR0,ps_BR0;
  __m256 ps_AI0,ps_BI0;
  __m256 ps_acARBR0;
  __m256 ps_acAIBI0;
  float *tmp_f_acARBR0_ = (float *) &ps_acARBR0;
  float *tmp_f_acAIBI0_ = (float *) &ps_acAIBI0;
  float tmp_f_acARBR0=0;
  float tmp_f_acAIBI0=0;
  float *f_CR__=NULL;
  float *f_CI__=NULL;
  int na=0;
  f_CR__=NULL;
  if (f_CR_p_!=NULL){
    if ( (*f_CR_p_)==NULL ){ (*f_CR_p_) = (float *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
    f_CR__ = *f_CR_p_;
    /* if (f_CR_p_!=NULL){ } */}
  f_CI__=NULL;
  if (f_CI_p_!=NULL){
    if ( (*f_CI_p_)==NULL ){ (*f_CI_p_) = (float *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
    f_CI__ = *f_CI_p_;
    /* if (f_CI_p_!=NULL){ } */}
  if ((f_CR__!=NULL) && (f_CI__!=NULL)){
    na=0;
    for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
	tmp_f_AR_point0_ = &(f_AR_trn__[nrow_A*n_col_X_rup]);
	tmp_f_AI_point0_ = &(f_AI_trn__[nrow_A*n_col_X_rup]);
	tmp_f_BR_point0_ = &(f_BR_trn__[nrow_B*n_col_X_rup]);
	tmp_f_BI_point0_ = &(f_BI_trn__[nrow_B*n_col_X_rup]);
	ps_acARBR0 = _mm256_set1_ps((float)0.0);
	ps_acAIBI0 = _mm256_set1_ps((float)0.0);
	for (ncol_X=0;ncol_X<n_col_X_256;ncol_X++){
	  ps_AR0 = _mm256_load_ps(tmp_f_AR_point0_);
	  ps_BR0 = _mm256_load_ps(tmp_f_BR_point0_);
	  ps_AI0 = _mm256_load_ps(tmp_f_AI_point0_);
	  ps_BI0 = _mm256_load_ps(tmp_f_BI_point0_);
#ifdef _FMA
	  ps_acARBR0 = _mm256_fmadd_ps(ps_AR0,ps_BR0,ps_acARBR0);
	  ps_acAIBI0 = _mm256_fmadd_ps(ps_AI0,ps_BI0,ps_acAIBI0);
#endif /* _FMA */
#ifdef _FMA
	  tmp_f_AR_point0_+=8;tmp_f_BR_point0_+=8;
	  tmp_f_AI_point0_+=8;tmp_f_BI_point0_+=8;
#endif /* _FMA */
	  /* for (ncol_X=0;ncol_X<n_col_X_256;ncol_X++){ } */}
	tmp_f_acARBR0 = tmp_f_acARBR0_[0] + tmp_f_acARBR0_[1] + tmp_f_acARBR0_[2] + tmp_f_acARBR0_[3] + tmp_f_acARBR0_[4] + tmp_f_acARBR0_[5] + tmp_f_acARBR0_[6] + tmp_f_acARBR0_[7];
	tmp_f_acAIBI0 = tmp_f_acAIBI0_[0] + tmp_f_acAIBI0_[1] + tmp_f_acAIBI0_[2] + tmp_f_acAIBI0_[3] + tmp_f_acAIBI0_[4] + tmp_f_acAIBI0_[5] + tmp_f_acAIBI0_[6] + tmp_f_acAIBI0_[7];
	f_CR__[na] = tmp_f_acARBR0 - tmp_f_acAIBI0;
	na += 1;
	/* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ } */}
    /* if ((f_CR__!=NULL) && (f_CI__!=NULL)){ } */}
}

void rcp_segregated_to_segregated_mult_immintrin_load1_fma(int n_row_A,int n_col_X,float *f_AR_trn__,float *f_AI_trn__,int n_row_B,float *f_BR_trn__,float *f_BI_trn__,float **f_CR_p_,float **f_CI_p_)
{
  /* non hermitian product */
  /* assumes A is real */
  /* assumes alignment */
  int n_col_X_rup = rup(n_col_X,8);
  int n_col_X_256 = n_col_X_rup/8; //%<-- 8 floats per __m256. ;
  int nrow_A=0,nrow_B=0,ncol_X=0;
  float *tmp_f_AR_point0_=NULL;
  float *tmp_f_BR_point0_=NULL;
  float *tmp_f_BI_point0_=NULL;
  __m256 ps_AR0,ps_BR0;
  __m256 ps_BI0;
  __m256 ps_acARBR0;
  __m256 ps_acARBI0;
  float *tmp_f_acARBR0_ = (float *) &ps_acARBR0;
  float *tmp_f_acARBI0_ = (float *) &ps_acARBI0;
  float tmp_f_acARBR0=0;
  float tmp_f_acARBI0=0;
  float *f_CR__=NULL;
  float *f_CI__=NULL;
  int na=0;
  f_CR__=NULL;
  if (f_CR_p_!=NULL){
    if ( (*f_CR_p_)==NULL ){ (*f_CR_p_) = (float *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
    f_CR__ = *f_CR_p_;
    /* if (f_CR_p_!=NULL){ } */}
  f_CI__=NULL;
  if (f_CI_p_!=NULL){
    if ( (*f_CI_p_)==NULL ){ (*f_CI_p_) = (float *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
    f_CI__ = *f_CI_p_;
    /* if (f_CI_p_!=NULL){ } */}
  if ((f_CR__!=NULL) && (f_CI__!=NULL)){
    na=0;
    for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
	tmp_f_AR_point0_ = &(f_AR_trn__[nrow_A*n_col_X_rup]);
	tmp_f_BR_point0_ = &(f_BR_trn__[nrow_B*n_col_X_rup]);
	tmp_f_BI_point0_ = &(f_BI_trn__[nrow_B*n_col_X_rup]);
	ps_acARBR0 = _mm256_set1_ps((float)0.0);
	ps_acARBI0 = _mm256_set1_ps((float)0.0);
	for (ncol_X=0;ncol_X<n_col_X_256;ncol_X++){
	  ps_AR0 = _mm256_load_ps(tmp_f_AR_point0_);
	  ps_BR0 = _mm256_load_ps(tmp_f_BR_point0_);
	  ps_BI0 = _mm256_load_ps(tmp_f_BI_point0_);
#ifdef _FMA
	  ps_acARBR0 = _mm256_fmadd_ps(ps_AR0,ps_BR0,ps_acARBR0);
	  ps_acARBI0 = _mm256_fmadd_ps(ps_AR0,ps_BI0,ps_acARBI0);
#endif /* _FMA */
#ifdef _FMA
	  tmp_f_AR_point0_+=8;tmp_f_BR_point0_+=8;
	  tmp_f_BI_point0_+=8;
#endif /* _FMA */
	  /* for (ncol_X=0;ncol_X<n_col_X_256;ncol_X++){ } */}
	tmp_f_acARBR0 = tmp_f_acARBR0_[0] + tmp_f_acARBR0_[1] + tmp_f_acARBR0_[2] + tmp_f_acARBR0_[3] + tmp_f_acARBR0_[4] + tmp_f_acARBR0_[5] + tmp_f_acARBR0_[6] + tmp_f_acARBR0_[7];
	tmp_f_acARBI0 = tmp_f_acARBI0_[0] + tmp_f_acARBI0_[1] + tmp_f_acARBI0_[2] + tmp_f_acARBI0_[3] + tmp_f_acARBI0_[4] + tmp_f_acARBI0_[5] + tmp_f_acARBI0_[6] + tmp_f_acARBI0_[7];
	f_CR__[na] = tmp_f_acARBR0 ;
	f_CI__[na] = tmp_f_acARBI0 ;
	na += 1;
	/* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ } */}
    /* if ((f_CR__!=NULL) && (f_CI__!=NULL)){ } */}
}

void rcp_ps_mult_immintrin_test()
{
  int verbose = 1;
  int n_row_A = 531;
  int n_col_X = 15039;
  int n_col_X_rup = rup(n_col_X,8);
  int n_col_X_256 = n_col_X_rup/8; //%<-- 8 floats per __m256. ;
  int n_row_B = 430;
  int n_row_A_sub = 5;
  int n_col_X_sub = 15;
  int n_row_B_sub = 4;
  int nrow_A=0,ncol_X=0,nrow_B=0;
  unsigned long long int ulli_A_total = (unsigned long long int) n_row_A * (unsigned long long int) n_col_X_rup ;
  unsigned long long int ulli_B_total = (unsigned long long int) n_row_B * (unsigned long long int) n_col_X_rup ;
  unsigned long long int ulli_C_total = (unsigned long long int) n_row_A * (unsigned long long int) n_row_B ;
  unsigned long long int ulli=0;
  unsigned long long int tab=0;
  __m256 *ps_AR_trn__ = NULL;
  __m256 *ps_AI_trn__ = NULL;
  __m256 *ps_AO_trn__ = NULL;
  __m256 *ps_BR_trn__ = NULL;
  __m256 *ps_BI_trn__ = NULL;
  float *f_AR_trn__ = NULL;
  float *f_AI_trn__ = NULL;
  float *f_AO_trn__ = NULL;
  float *f_AR_u_trn__ = NULL;
  float *f_AI_u_trn__ = NULL;
  float *f_AO_u_trn__ = NULL;
  float *f_BR_trn__ = NULL;
  float *f_BI_trn__ = NULL;
  float *f_BR_u_trn__ = NULL;
  float *f_BI_u_trn__ = NULL;
  float *f_CR_bf__ = NULL;
  float *f_CI_bf__ = NULL;
  float *f_CR_al__ = NULL;
  float *f_CI_al__ = NULL;
  float *f_AR_sub__ = NULL;
  float *f_AI_sub__ = NULL;
  float *f_AO_sub__ = NULL;
  float *f_BR_sub__ = NULL;
  float *f_BI_sub__ = NULL;
  float *f_CR_sub__ = NULL;
  float *f_CI_sub__ = NULL;
  float complex *c_A_u_trn__ = NULL;
  float complex *c_B_u_trn__ = NULL;
  float complex *c_C_bf__ = NULL;
  float complex *c_C_al__ = NULL;
  float complex *c_A_sub__ = NULL;
  float complex *c_B_sub__ = NULL;
  float complex *c_C_sub__ = NULL;
  float ferror=0;
  /* %%%% */
  int ntick=0,n_tick = 8;
  clock_t t_start_[n_tick];
  clock_t t_final_[n_tick];
  struct timeval d_start_[n_tick];
  struct timeval d_final_[n_tick];
  long l_msec_[n_tick],l_ssec_[n_tick],l_usec_[n_tick];
  double elct_[n_tick],elrt_[n_tick];
  /* %%%% */
  ps_AR_trn__ = (__m256 *) _mm_malloc(n_col_X_256*n_row_A*sizeof(__m256),32);
  ps_AI_trn__ = (__m256 *) _mm_malloc(n_col_X_256*n_row_A*sizeof(__m256),32);
  ps_AO_trn__ = (__m256 *) _mm_malloc(n_col_X_256*n_row_A*sizeof(__m256),32);
  ps_BR_trn__ = (__m256 *) _mm_malloc(n_col_X_256*n_row_B*sizeof(__m256),32);
  ps_BI_trn__ = (__m256 *) _mm_malloc(n_col_X_256*n_row_B*sizeof(__m256),32);
  f_AR_trn__ = (float *) ps_AR_trn__;
  f_AI_trn__ = (float *) ps_AI_trn__;
  f_AO_trn__ = (float *) ps_AO_trn__;
  f_BR_trn__ = (float *) ps_BR_trn__;
  f_BI_trn__ = (float *) ps_BI_trn__;
  f_AR_u_trn__ = (float *) malloc(n_col_X*n_row_A*sizeof(float));
  f_AI_u_trn__ = (float *) malloc(n_col_X*n_row_A*sizeof(float));
  f_AO_u_trn__ = (float *) malloc(n_col_X*n_row_A*sizeof(float));
  f_BR_u_trn__ = (float *) malloc(n_col_X*n_row_B*sizeof(float));
  f_BI_u_trn__ = (float *) malloc(n_col_X*n_row_B*sizeof(float));
  f_CR_bf__ = (float *) malloc(n_row_A*n_row_B*sizeof(float));
  f_CI_bf__ = (float *) malloc(n_row_A*n_row_B*sizeof(float));
  f_CR_al__ = (float *) malloc(n_row_A*n_row_B*sizeof(float));
  f_CI_al__ = (float *) malloc(n_row_A*n_row_B*sizeof(float));
  f_AR_sub__ = (float *) malloc(n_row_A_sub*n_col_X_sub*sizeof(float));
  f_AI_sub__ = (float *) malloc(n_row_A_sub*n_col_X_sub*sizeof(float));
  f_AO_sub__ = (float *) malloc(n_row_A_sub*n_col_X_sub*sizeof(float));
  f_BR_sub__ = (float *) malloc(n_row_B_sub*n_col_X_sub*sizeof(float));
  f_BI_sub__ = (float *) malloc(n_row_B_sub*n_col_X_sub*sizeof(float));
  f_CR_sub__ = (float *) malloc(n_row_A_sub*n_row_B_sub*sizeof(float));
  f_CI_sub__ = (float *) malloc(n_row_A_sub*n_row_B_sub*sizeof(float));
  c_A_u_trn__ = (float complex *) malloc(n_col_X*n_row_A*sizeof(float complex));
  c_B_u_trn__ = (float complex *) malloc(n_col_X*n_row_B*sizeof(float complex));
  c_C_bf__ = (float complex *) malloc(n_row_A*n_row_B*sizeof(float complex));
  c_C_al__ = (float complex *) malloc(n_row_A*n_row_B*sizeof(float complex));
  c_A_sub__ = (float complex *) malloc(n_row_A_sub*n_col_X_sub*sizeof(float complex));
  c_B_sub__ = (float complex *) malloc(n_row_B_sub*n_col_X_sub*sizeof(float complex));
  c_C_sub__ = (float complex *) malloc(n_row_A_sub*n_row_B_sub*sizeof(float complex));
  tab = (unsigned long long int) n_row_A * (unsigned long long int) n_row_B * (unsigned long long int) n_col_X_rup;
  local_tic(0,t_start_,d_start_);
  ulli=0;
  for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ for (ncol_X=0;ncol_X<n_col_X;ncol_X++){
      f_AR_trn__[ncol_X + nrow_A*n_col_X_rup] = (float)(((int)ulli%7)-3);
      f_AI_trn__[ncol_X + nrow_A*n_col_X_rup] = (float)(((int)ulli%7)-2);
      f_AO_trn__[ncol_X + nrow_A*n_col_X_rup] = (float)(0);
      f_AR_u_trn__[ncol_X + nrow_A*n_col_X] = (float)(((int)ulli%7)-3);
      f_AI_u_trn__[ncol_X + nrow_A*n_col_X] = (float)(((int)ulli%7)-2);
      f_AO_u_trn__[ncol_X + nrow_A*n_col_X] = (float)(0);
      c_A_u_trn__[ncol_X + nrow_A*n_col_X] = (float complex) f_AR_u_trn__[ncol_X + nrow_A*n_col_X] + _Complex_I * (float complex) f_AO_u_trn__[ncol_X + nrow_A*n_col_X];
      ulli++;
      /* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ for (ncol_X=0;ncol_X<n_col_X;ncol_X++){ }} */}}  
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (ncol_X=0;ncol_X<n_col_X_sub;ncol_X++){ 
      f_AR_sub__[nrow_A+ncol_X*n_row_A_sub] = f_AR_trn__[ncol_X + nrow_A*n_col_X_rup];
      f_AI_sub__[nrow_A+ncol_X*n_row_A_sub] = f_AI_trn__[ncol_X + nrow_A*n_col_X_rup];
      f_AO_sub__[nrow_A+ncol_X*n_row_A_sub] = f_AO_trn__[ncol_X + nrow_A*n_col_X_rup];
      c_A_sub__[nrow_A+ncol_X*n_row_A_sub] = (float complex) f_AR_sub__[nrow_A+ncol_X*n_row_A_sub] + _Complex_I * (float complex) f_AO_sub__[nrow_A+ncol_X*n_row_A_sub];
      /* for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (ncol_X=0;ncol_X<n_col_X_sub;ncol_X++){ }} */}}
  if (verbose>1){
    printf(" %% upper corner of f_AR_trn__: \n");
    array_printf(f_AR_sub__,"float",n_row_A_sub,n_col_X_sub," % f_AR_sub__: ");
    printf(" %% upper corner of f_AI_trn__: \n");
    array_printf(f_AI_sub__,"float",n_row_A_sub,n_col_X_sub," % f_AI_sub__: ");
    printf(" %% upper corner of f_AO_trn__: \n");
    array_printf(f_AO_sub__,"float",n_row_A_sub,n_col_X_sub," % f_AO_sub__: ");
    printf(" %% upper corner of c_A_trn__: \n");
    array_printf(c_A_sub__,"float complex",n_row_A_sub,n_col_X_sub," % c_A_sub__: ");
    /* if (verbose>1){ } */}
  ulli=0;
  for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ for (ncol_X=0;ncol_X<n_col_X;ncol_X++){
      f_BR_trn__[ncol_X + nrow_B*n_col_X_rup] = (float)(((int)ulli%11)-2);
      f_BI_trn__[ncol_X + nrow_B*n_col_X_rup] = (float)(((int)ulli%11)-1);
      f_BR_u_trn__[ncol_X + nrow_B*n_col_X] = (float)(((int)ulli%11)-2);
      f_BI_u_trn__[ncol_X + nrow_B*n_col_X] = (float)(((int)ulli%11)-1);
      c_B_u_trn__[ncol_X + nrow_B*n_col_X] = (float complex) f_BR_u_trn__[ncol_X + nrow_B*n_col_X] + _Complex_I * (float complex) f_BI_u_trn__[ncol_X + nrow_B*n_col_X];
      ulli++;
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ for (ncol_X=0;ncol_X<n_col_X;ncol_X++){ }} */}}
  if (verbose>1){
    for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ for (ncol_X=0;ncol_X<n_col_X_sub;ncol_X++){
	f_BR_sub__[nrow_B+ncol_X*n_row_B_sub] = f_BR_trn__[ncol_X + nrow_B*n_col_X_rup];
	f_BI_sub__[nrow_B+ncol_X*n_row_B_sub] = f_BI_trn__[ncol_X + nrow_B*n_col_X_rup];
	c_B_sub__[nrow_B+ncol_X*n_row_B_sub] = (float complex) f_BR_sub__[nrow_B+ncol_X*n_row_B_sub] + _Complex_I * (float complex) f_BI_sub__[nrow_B+ncol_X*n_row_B_sub];
	/* for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ for (ncol_X=0;ncol_X<n_col_X_sub;ncol_X++){ }} */}}
    printf(" %% upper corner of f_BR_trn__: \n");
    array_printf(f_BR_sub__,"float",n_row_B_sub,n_col_X_sub," % f_BR_sub__: ");
    printf(" %% upper corner of f_BI_trn__: \n");
    array_printf(f_BI_sub__,"float",n_row_B_sub,n_col_X_sub," % f_BI_sub__: ");
    printf(" %% upper corner of c_B_trn__: \n");
    array_printf(c_B_sub__,"float complex",n_row_B_sub,n_col_X_sub," % c_B_sub__: ");
    /* if (verbose>1){ } */}
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," initialize: ");
  /* %%%%%%%% */
  local_tic(0,t_start_,d_start_);
  memset(c_C_bf__,0,ulli_C_total*sizeof(float complex));
  hp_interleave_mult_bruteforce(n_row_A,n_col_X,c_A_u_trn__,n_row_B,c_B_u_trn__,&c_C_bf__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," hp_ps_mult_bruteforce: ");
  printf(" %% Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/elrt_[0]/1e9);
  if (verbose>1){
    for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){
	c_C_sub__[nrow_A+nrow_B*n_row_A_sub] = c_C_bf__[nrow_A+nrow_B*n_row_A];
	/* for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ }} */}}
    printf(" %% upper corner of c_C_bf__: \n");
    array_printf(c_C_sub__,"float complex",n_row_A_sub,n_row_B_sub," % c_C_bf__: ");
    /* if (verbose>1){ } */}
  ferror = cfnorm(ulli_C_total,c_C_bf__,c_C_bf__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
  local_tic(0,t_start_,d_start_);
  memset(c_C_al__,0,ulli_C_total*sizeof(float complex));
  hp_segregated_mult_bruteforce(n_row_A,n_col_X,f_AR_u_trn__,f_AO_u_trn__,n_row_B,f_BR_u_trn__,f_BI_u_trn__,&c_C_al__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," hp_segregated_mult_bruteforce: ");
  printf(" %% Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/elrt_[0]/1e9);
  if (verbose>1){
    for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){
	c_C_sub__[nrow_A+nrow_B*n_row_A_sub] = c_C_al__[nrow_A+nrow_B*n_row_A];
	/* for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ }} */}}
    printf(" %% upper corner of c_C_al__: \n");
    array_printf(c_C_sub__,"float complex",n_row_A_sub,n_row_B_sub," % c_C_al__: ");
    /* if (verbose>1){ } */}
  ferror = cfnorm(ulli_C_total,c_C_bf__,c_C_al__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
  local_tic(0,t_start_,d_start_);
  memset(c_C_al__,0,ulli_C_total*sizeof(float complex));
  hp_segregated_to_interleaved_mult_immintrin_loadu_fma(n_row_A,n_col_X,f_AR_u_trn__,f_AO_u_trn__,n_row_B,f_BR_u_trn__,f_BI_u_trn__,&c_C_al__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," hp_segregated_to_interleaved_mult_immintrin_loadu_fma: ");
  printf(" %% Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/elrt_[0]/1e9);
  if (verbose>1){
    for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){
	c_C_sub__[nrow_A+nrow_B*n_row_A_sub] = c_C_al__[nrow_A+nrow_B*n_row_A];
	/* for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ }} */}}
    printf(" %% upper corner of c_C_al__: \n");
    array_printf(c_C_sub__,"float complex",n_row_A_sub,n_row_B_sub," % c_C_al__: ");
    /* if (verbose>1){ } */}
  ferror = cfnorm(ulli_C_total,c_C_bf__,c_C_al__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
  local_tic(0,t_start_,d_start_);
  memset(f_CR_al__,0,ulli_C_total*sizeof(float));
  memset(f_CI_al__,0,ulli_C_total*sizeof(float));
  hp_segregated_to_segregated_mult_immintrin_loadu_fma(n_row_A,n_col_X,f_AR_u_trn__,f_AO_u_trn__,n_row_B,f_BR_u_trn__,f_BI_u_trn__,&f_CR_al__,&f_CI_al__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," hp_segregated_to_segregated_mult_immintrin_loadu_fma: ");
  printf(" %% Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/elrt_[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      c_C_al__[nrow_A+nrow_B*n_row_A] = f_CR_al__[nrow_A+nrow_B*n_row_A] + _Complex_I * (float complex) f_CI_al__[nrow_A+nrow_B*n_row_A];
      /* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ }} */}}
  if (verbose>1){
    for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){
	c_C_sub__[nrow_A+nrow_B*n_row_A_sub] = c_C_al__[nrow_A+nrow_B*n_row_A];
	/* for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ }} */}}
    printf(" %% upper corner of c_C_al__: \n");
    array_printf(c_C_sub__,"float complex",n_row_A_sub,n_row_B_sub," % c_C_al__: ");
    /* if (verbose>1){ } */}
  ferror = cfnorm(ulli_C_total,c_C_bf__,c_C_al__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
  local_tic(0,t_start_,d_start_);
  memset(f_CR_al__,0,ulli_C_total*sizeof(float));
  memset(f_CI_al__,0,ulli_C_total*sizeof(float));
  hp_segregated_to_segregated_mult_immintrin_load1_fma(n_row_A,n_col_X,f_AR_trn__,f_AO_trn__,n_row_B,f_BR_trn__,f_BI_trn__,&f_CR_al__,&f_CI_al__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," hp_segregated_to_segregated_mult_immintrin_load1_fma: ");
  printf(" %% Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/elrt_[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      c_C_al__[nrow_A+nrow_B*n_row_A] = f_CR_al__[nrow_A+nrow_B*n_row_A] + _Complex_I * (float complex) f_CI_al__[nrow_A+nrow_B*n_row_A];
      /* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ }} */}}
  if (verbose>1){
    for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){
	c_C_sub__[nrow_A+nrow_B*n_row_A_sub] = c_C_al__[nrow_A+nrow_B*n_row_A];
	/* for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ }} */}}
    printf(" %% upper corner of c_C_al__: \n");
    array_printf(c_C_sub__,"float complex",n_row_A_sub,n_row_B_sub," % c_C_al__: ");
    /* if (verbose>1){ } */}
  ferror = cfnorm(ulli_C_total,c_C_bf__,c_C_al__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
  local_tic(0,t_start_,d_start_);
  memset(f_CR_al__,0,ulli_C_total*sizeof(float));
  memset(f_CI_al__,0,ulli_C_total*sizeof(float));
  rcp_segregated_to_segregated_mult_immintrin_load1_fma(n_row_A,n_col_X,f_AR_trn__,NULL,n_row_B,f_BR_trn__,f_BI_trn__,&f_CR_al__,&f_CI_al__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," rcp_segregated_to_segregated_mult_immintrin_load1_fma: ");
  printf(" %% Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/elrt_[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      c_C_al__[nrow_A+nrow_B*n_row_A] = f_CR_al__[nrow_A+nrow_B*n_row_A] + _Complex_I * (float complex) f_CI_al__[nrow_A+nrow_B*n_row_A];
      /* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ }} */}}
  if (verbose>1){
    for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){
	c_C_sub__[nrow_A+nrow_B*n_row_A_sub] = c_C_al__[nrow_A+nrow_B*n_row_A];
	/* for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ }} */}}
    printf(" %% upper corner of c_C_al__: \n");
    array_printf(c_C_sub__,"float complex",n_row_A_sub,n_row_B_sub," % c_C_al__: ");
    /* if (verbose>1){ } */}
  ferror = cfnorm(ulli_C_total,c_C_bf__,c_C_al__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
#ifdef _CBLAS
  local_tic(0,t_start_,d_start_);
  memset(c_C_al__,0,ulli_C_total*sizeof(float complex));
  hp_interleave_mult_cblas_cgemm(n_row_A,n_col_X,c_A_u_trn__,n_row_B,c_B_u_trn__,&c_C_al__);
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," hp_interleave_mult_cblas_cgemm: ");
  printf(" %% Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/elrt_[0]/1e9);
  if (verbose>1){
    for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){
	c_C_sub__[nrow_A+nrow_B*n_row_A_sub] = c_C_al__[nrow_A+nrow_B*n_row_A];
	/* for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ }} */}}
    printf(" %% upper corner of c_C_al__: \n");
    array_printf(c_C_sub__,"float complex",n_row_A_sub,n_row_B_sub," % c_C_al__: ");
    /* if (verbose>1){ } */}
  ferror = cfnorm(ulli_C_total,c_C_bf__,c_C_al__);
  printf(" %% ferror %0.16f\n",ferror);
#endif /* _CBLAS */
  /* %%%%%%%% */
  _mm_free(ps_AR_trn__); ps_AR_trn__ = NULL;
  _mm_free(ps_AI_trn__); ps_AI_trn__ = NULL;
  _mm_free(ps_AO_trn__); ps_AO_trn__ = NULL;
  _mm_free(ps_BR_trn__); ps_BR_trn__ = NULL;
  _mm_free(ps_BI_trn__); ps_BI_trn__ = NULL;
  free(f_AR_u_trn__); f_AR_u_trn__=NULL;
  free(f_AI_u_trn__); f_AI_u_trn__=NULL;
  free(f_AO_u_trn__); f_AO_u_trn__=NULL;
  free(f_BR_u_trn__); f_BR_u_trn__=NULL;
  free(f_BI_u_trn__); f_BI_u_trn__=NULL;
  free(f_CR_bf__); f_CR_bf__=NULL;
  free(f_CI_bf__); f_CI_bf__=NULL;
  free(f_CR_al__); f_CR_al__=NULL;
  free(f_CI_al__); f_CI_al__=NULL;
  free(f_AR_sub__); f_AR_sub__=NULL;
  free(f_AI_sub__); f_AI_sub__=NULL;
  free(f_AO_sub__); f_AO_sub__=NULL;
  free(f_BR_sub__); f_BR_sub__=NULL;
  free(f_BI_sub__); f_BI_sub__=NULL;
  free(f_CR_sub__); f_CR_sub__=NULL;
  free(f_CI_sub__); f_CI_sub__=NULL;
  free(c_A_u_trn__); c_A_u_trn__=NULL;
  free(c_B_u_trn__); c_B_u_trn__=NULL;
  free(c_C_bf__); c_C_bf__=NULL;
  free(c_C_al__); c_C_al__=NULL;
  free(c_A_sub__); c_A_sub__=NULL;
  free(c_B_sub__); c_B_sub__=NULL;
  free(c_C_sub__); c_C_sub__=NULL;
  //wkspace_printf();
}
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
