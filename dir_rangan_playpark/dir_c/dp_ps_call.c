#ifndef _MONOLITH
#include "playpark_header.h"
#endif /* _MONOLITH */

void dp_pd_bruteforce(int n_col_X,double *d_A_,double *d_B_,double *d_C_)
{
  int ncol_X=0;
  double *d_A_point0_=NULL,*d_B_point0_=NULL;
  double d_ac0=0;
  d_A_point0_ = d_A_;
  d_B_point0_ = d_B_;
  for (ncol_X=0;ncol_X<n_col_X;ncol_X++){
    d_ac0 += *(d_A_point0_++) * *(d_B_point0_++);
    /* for (ncol_X=0;ncol_X<n_col_X;ncol_X++){ } */}
  *d_C_ = d_ac0;
}

void dp_pd_mult_bruteforce(int n_row_A,int n_col_X,double *d_A_trn__,int n_row_B,double *d_B_trn__,double **d_C_p_)
{
  int n_col_X_rup = rup(n_col_X,8);
  int nrow_A=0,nrow_B=0,ncol_X=0;
  double *d_C__=NULL;
  d_C__=NULL;
  if (d_C_p_!=NULL){
    if ( (*d_C_p_)==NULL ){ (*d_C_p_) = (double *) malloc1((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(double));}
    d_C__ = *d_C_p_;
    /* if (d_C_p_!=NULL){ } */}
  if (d_C__!=NULL){
    for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
	for (ncol_X=0;ncol_X<n_col_X;ncol_X++){
	  d_C__[nrow_A + nrow_B*n_row_A] += d_A_trn__[ncol_X + nrow_A*n_col_X_rup] * d_B_trn__[ncol_X + nrow_B*n_col_X_rup];
	  /* for (ncol_X=0;ncol_X<n_col_X;ncol_X++){ } */}
	/* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ } */}
    /* if (d_C__!=NULL){ } */}
}

void dp_pd_immintrin_loadu(int n_col_X,double *d_A_,double *d_B_,double *d_C_)
{
  /* does not assume alignment */
  int n_col_X_rup = rup(n_col_X,4);
  int n_col_X_256d = n_col_X_rup/4; //%<-- 4 doubles per __m256d. ;
  int ncol_X=0,ncol_X_256d=0;
  double *d_A_point0_=NULL,*d_B_point0_=NULL;
  __m256d pd_A0,pd_B0;
  __m256d pd_ac0;
  double *d_ac0_ = (double *) &pd_ac0;
  double d_ac0=0;
  d_A_point0_ = d_A_;
  d_B_point0_ = d_B_;
  pd_ac0  = _mm256_set1_pd((double)0.0);
  for (ncol_X_256d=0;ncol_X_256d<n_col_X_256d-1;ncol_X_256d++){
    pd_A0 = _mm256_loadu_pd(d_A_point0_);
    pd_B0 = _mm256_loadu_pd(d_B_point0_);
    pd_ac0 = _mm256_fmadd_pd(pd_A0,pd_B0,pd_ac0);
    d_A_point0_+=4;d_B_point0_+=4;
    /* for (ncol_X_256d=0;ncol_X_256d<n_col_X_256d-1;ncol_X_256d++){ } */}
  d_ac0 = d_ac0_[0] + d_ac0_[1] + d_ac0_[2] + d_ac0_[3];
  d_A_point0_ = &(d_A_[4*ncol_X_256d]);
  d_B_point0_ = &(d_B_[4*ncol_X_256d]);
  for (ncol_X=4*ncol_X_256d;ncol_X<n_col_X;ncol_X++){
    d_ac0 += *(d_A_point0_++) * *(d_B_point0_++);
    /* for (ncol_X=4*ncol_X_256d;ncol_X<n_col_X;ncol_X++){ } */}
  *d_C_ = d_ac0;
}

void dp_pd_mult_immintrin_loadu(int n_row_A,int n_col_X,double *d_A_trn__,int n_row_B,double *d_B_trn__,double **d_C_p_)
{
  /* does not assume alignment */
  int nrow_A=0,nrow_B=0;
  double *tmp_d_A_point0_=NULL;
  double *tmp_d_B_point0_=NULL;
  double *d_C__=NULL;
  d_C__=NULL;
  if (d_C_p_!=NULL){
    if ( (*d_C_p_)==NULL ){ (*d_C_p_) = (double *) malloc1((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(double));}
    d_C__ = *d_C_p_;
    /* if (d_C_p_!=NULL){ } */}
  if (d_C__!=NULL){
    for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
	tmp_d_A_point0_ = &(d_A_trn__[nrow_A*n_col_X]);
	tmp_d_B_point0_ = &(d_B_trn__[nrow_B*n_col_X]);
	dp_pd_immintrin_loadu(n_col_X,tmp_d_A_point0_,tmp_d_B_point0_,&(d_C__[nrow_A + nrow_B*n_row_A]));
	/* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ } */}
    /* if (d_C__!=NULL){ } */}
}

void dp_pd_mult_immintrin_test()
{
  /* int n_row_A = 1024*1 + 723; */
  /* int n_col_X = 1024*6 + 817; */
  /* int n_row_B = 1024*1 + 511; */
  int n_row_A = 531;
  int n_col_X = 15039;
  int n_col_X_rup = rup(n_col_X,8);
  int n_col_X_256d = n_col_X_rup/4; //%<-- 4 doubles per __m256d. ;
  int n_row_B = 430;
  int n_row_A_sub = 5;
  int n_col_X_sub = 15;
  int n_row_B_sub = 4;
  int nrow_A=0,ncol_X=0,nrow_B=0;
  unsigned long long int ulli_A_total = (unsigned long long int) n_row_A * (unsigned long long int) n_col_X_rup ;
  unsigned long long int ulli_B_total = (unsigned long long int) n_row_B * (unsigned long long int) n_col_X_rup ;
  unsigned long long int ulli_C_total = (unsigned long long int) n_row_A * (unsigned long long int) n_row_B ;
  unsigned long long int ulli=0;
  __m256d *pd_A_trn__ = NULL;
  __m256d *pd_B_trn__ = NULL;
  double *d_A_trn__ = NULL;
  double *d_A_u_trn__ = NULL;
  double *d_B_trn__ = NULL;
  double *d_B_u_trn__ = NULL;
  double *d_C_bf__ = NULL;
  double *d_C_pd__ = NULL;
  double *d_A_sub__ = NULL;
  double *d_B_sub__ = NULL;
  double *d_C_sub__ = NULL;
  double derror=0;
  pd_A_trn__ = (__m256d *) _mm_malloc(n_col_X_256d*n_row_A*sizeof(__m256d),64);
  pd_B_trn__ = (__m256d *) _mm_malloc(n_col_X_256d*n_row_A*sizeof(__m256d),64);
  d_A_trn__ = (double *) pd_A_trn__;
  d_B_trn__ = (double *) pd_B_trn__;
  d_A_u_trn__ = (double *) malloc1(n_col_X*n_row_A*sizeof(double));
  d_B_u_trn__ = (double *) malloc1(n_col_X*n_row_A*sizeof(double));
  d_C_bf__ = (double *) malloc1(n_row_A*n_row_B*sizeof(double));
  d_C_pd__ = (double *) malloc1(n_row_A*n_row_B*sizeof(double));
  d_A_sub__ = (double *) malloc1(n_row_A_sub*n_col_X_sub*sizeof(double));
  d_B_sub__ = (double *) malloc1(n_row_B_sub*n_col_X_sub*sizeof(double));
  d_C_sub__ = (double *) malloc1(n_row_A_sub*n_row_B_sub*sizeof(double));
  GLOBAL_tic(0);
  ulli=0;
  for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ for (ncol_X=0;ncol_X<n_col_X;ncol_X++){
      d_A_trn__[ncol_X + nrow_A*n_col_X_rup] = (double)(((int)ulli%7)-3); 
      d_A_u_trn__[ncol_X + nrow_A*n_col_X] = (double)(((int)ulli%7)-3); 
      ulli++;
      /* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ for (ncol_X=0;ncol_X<n_col_X;ncol_X++){ }} */}}
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (ncol_X=0;ncol_X<n_col_X_sub;ncol_X++){ d_A_sub__[nrow_A+ncol_X*n_row_A_sub] = d_A_trn__[ncol_X + nrow_A*n_col_X_rup];}}
  printf(" %% upper corner of d_A_trn__: \n");
  array_printf(d_A_sub__,"double",n_row_A_sub,n_col_X_sub," % d_A_sub__: ");
  ulli=0;
  for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ for (ncol_X=0;ncol_X<n_col_X;ncol_X++){
      d_B_trn__[ncol_X + nrow_B*n_col_X_rup] = (double)(((int)ulli%5)-2);
      d_B_u_trn__[ncol_X + nrow_B*n_col_X] = (double)(((int)ulli%5)-2);
      ulli++;
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ for (ncol_X=0;ncol_X<n_col_X;ncol_X++){ }} */}}
  for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ for (ncol_X=0;ncol_X<n_col_X_sub;ncol_X++){ d_B_sub__[nrow_B+ncol_X*n_row_B_sub] = d_B_trn__[ncol_X + nrow_B*n_col_X_rup];}}
  printf(" %% upper corner of d_B_trn__: \n");
  array_printf(d_B_sub__,"double",n_row_B_sub,n_col_X_sub," % d_B_sub__: ");
  GLOBAL_toc(0,1," initialize: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  memset(d_C_bf__,0,ulli_C_total*sizeof(double));
  dp_pd_mult_bruteforce(n_row_A,n_col_X,d_A_trn__,n_row_B,d_B_trn__,&d_C_bf__);
  GLOBAL_toc(0,1," dp_pd_mult_bruteforce: ");
  printf("Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/GLOBAL_elrt[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ d_C_sub__[nrow_A+nrow_B*n_row_A_sub] = d_C_bf__[nrow_A+nrow_B*n_row_A];}}
  printf(" %% upper corner of d_C_bf__: \n");
  array_printf(d_C_sub__,"double",n_row_A_sub,n_row_B_sub," % d_C_bf__: ");
  derror = dfnorm(ulli_C_total,d_C_bf__,d_C_bf__);
  printf(" %% derror %0.16f\n",derror);
  /* %%%%%%%% */
  GLOBAL_tic(0);
  memset(d_C_pd__,0,ulli_C_total*sizeof(double));
  dp_pd_mult_immintrin_loadu(n_row_A,n_col_X,d_A_u_trn__,n_row_B,d_B_u_trn__,&d_C_pd__);
  GLOBAL_toc(0,1," dp_pd_mult_immintrin_loadu: ");
  printf("Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X/GLOBAL_elrt[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ d_C_sub__[nrow_A+nrow_B*n_row_A_sub] = d_C_pd__[nrow_A+nrow_B*n_row_A];}}
  printf(" %% upper corner of d_C_pd__: \n");
  array_printf(d_C_sub__,"double",n_row_A_sub,n_row_B_sub," % d_C_pd__: ");
  derror = dfnorm(ulli_C_total,d_C_bf__,d_C_pd__);
  printf(" %% derror %0.16f\n",derror);
  /* %%%%%%%% */
  _mm_free(pd_A_trn__); pd_A_trn__ = NULL;
  _mm_free(pd_B_trn__); pd_B_trn__ = NULL;
  free1(&d_A_u_trn__);
  free1(&d_B_u_trn__);
  free1(&d_C_bf__);
  free1(&d_C_pd__);
  free1(&d_A_sub__);
  free1(&d_B_sub__);
  free1(&d_C_sub__);
  //wkspace_printf();
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void dp_ps_mult_cblas_sgemm(int n_row_A,int n_col_X,float *f_A_trn__,int n_row_B,float *f_B_trn__,float **f_C_p_)
{
  int n_col_X_rup = rup(n_col_X,8);
  float *f_C__=NULL;
  f_C__=NULL;
  if (f_C_p_!=NULL){
    if ( (*f_C_p_)==NULL ){ (*f_C_p_) = (float *) malloc1((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
    f_C__ = *f_C_p_;
    /* if (f_C_p_!=NULL){ } */}
  if (f_C__!=NULL){
    cblas_sgemm(CblasColMajor,CblasTrans,CblasNoTrans,n_row_A,n_row_B,n_col_X,1.0,f_A_trn__,n_col_X_rup,f_B_trn__,n_col_X_rup,0.0,f_C__,n_row_A);
  /* if (f_C__!=NULL){ } */}
}

void dp_ps_mult_cblas_sdot(int n_row_A,int n_col_X,float *f_A_trn__,int n_row_B,float *f_B_trn__,float **f_C_p_)
{
  int n_col_X_rup = rup(n_col_X,8);
  int nrow_A=0,nrow_B=0,ncol_X=0;
  float *f_C__=NULL;
  f_C__=NULL;
  if (f_C_p_!=NULL){
    if ( (*f_C_p_)==NULL ){ (*f_C_p_) = (float *) malloc1((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
    f_C__ = *f_C_p_;
    /* if (f_C_p_!=NULL){ } */}
  if (f_C__!=NULL){
    for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
	f_C__[nrow_A + nrow_B*n_row_A] = cblas_sdot(n_col_X,f_A_trn__+nrow_A*n_col_X_rup,1,f_B_trn__+nrow_B*n_col_X_rup,1);
	/* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ } */}
    /* if (f_C__!=NULL){ } */}
}

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
    if ( (*f_C_p_)==NULL ){ (*f_C_p_) = (float *) malloc1((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
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

void dp_ps_immintrin_loadu(int n_col_X,float *f_A_,float *f_B_,float *f_C_)
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
    ps_ac0 = _mm256_fmadd_ps(ps_A0,ps_B0,ps_ac0);
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

void dp_ps_mult_immintrin_loadu(int n_row_A,int n_col_X,float *f_A_trn__,int n_row_B,float *f_B_trn__,float **f_C_p_)
{
  /* does not assume alignment */
  int nrow_A=0,nrow_B=0;
  float *tmp_f_A_point0_=NULL;
  float *tmp_f_B_point0_=NULL;
  float *f_C__=NULL;
  f_C__=NULL;
  if (f_C_p_!=NULL){
    if ( (*f_C_p_)==NULL ){ (*f_C_p_) = (float *) malloc1((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
    f_C__ = *f_C_p_;
    /* if (f_C_p_!=NULL){ } */}
  if (f_C__!=NULL){
    for (nrow_B=0;nrow_B<n_row_B;nrow_B++){
      for (nrow_A=0;nrow_A<n_row_A;nrow_A++){
	tmp_f_A_point0_ = &(f_A_trn__[nrow_A*n_col_X]);
	tmp_f_B_point0_ = &(f_B_trn__[nrow_B*n_col_X]);
	dp_ps_immintrin_loadu(n_col_X,tmp_f_A_point0_,tmp_f_B_point0_,&(f_C__[nrow_A + nrow_B*n_row_A]));
	/* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ } */}
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ } */}
    /* if (f_C__!=NULL){ } */}
}

void dp_ps_mult_immintrin_load1(int n_row_A,int n_col_X,float *f_A_trn__,int n_row_B,float *f_B_trn__,float **f_C_p_)
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
    if ( (*f_C_p_)==NULL ){ (*f_C_p_) = (float *) malloc1((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
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
	  ps_ac0 = _mm256_fmadd_ps(ps_A0,ps_B0,ps_ac0);
	  tmp_f_A_point0_+=8;tmp_f_B_point0_+=8;
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
    if ( (*f_C_p_)==NULL ){ (*f_C_p_) = (float *) malloc1((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
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
	  ps_ac0 = _mm256_fmadd_ps(*(tmp_ps_A_trn_),*(tmp_ps_B_trn_),ps_ac0);
	  ps_ac1 = _mm256_fmadd_ps(*(tmp_ps_A_trn_+1),*(tmp_ps_B_trn_+1),ps_ac1);
	  tmp_ps_A_trn_+=2;
	  tmp_ps_B_trn_+=2;
	  /* for (ncol_X=0;ncol_X<n_col_X_256-1;ncol_X+=2){ } */}
	if ( (nrow_A==0) && (nrow_B==0) ){ printf(" ncol_X %d\n",ncol_X);}
	if (n_col_X_256%2==1){ ps_ac0 = _mm256_fmadd_ps(*(tmp_ps_A_trn_++),*(tmp_ps_B_trn_++),ps_ac0);}
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
    if ( (*f_C_p_)==NULL ){ (*f_C_p_) = (float *) malloc1((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
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
	  ps_acc = _mm256_fmadd_ps(*(tmp_ps_A_trn_++),*(tmp_ps_B_trn_++),ps_acc);
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
    if ( (*f_C_p_)==NULL ){ (*f_C_p_) = (float *) malloc1((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
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

void dp_ps_mult_immintrin_test()
{
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
  ps_A_trn__ = (__m256 *) _mm_malloc(n_col_X_256*n_row_A*sizeof(__m256),32);
  ps_B_trn__ = (__m256 *) _mm_malloc(n_col_X_256*n_row_A*sizeof(__m256),32);
  f_A_trn__ = (float *) ps_A_trn__;
  f_B_trn__ = (float *) ps_B_trn__;
  f_A_u_trn__ = (float *) malloc1(n_col_X*n_row_A*sizeof(float));
  f_B_u_trn__ = (float *) malloc1(n_col_X*n_row_A*sizeof(float));
  f_C_bf__ = (float *) malloc1(n_row_A*n_row_B*sizeof(float));
  f_C_ps__ = (float *) malloc1(n_row_A*n_row_B*sizeof(float));
  f_A_sub__ = (float *) malloc1(n_row_A_sub*n_col_X_sub*sizeof(float));
  f_B_sub__ = (float *) malloc1(n_row_B_sub*n_col_X_sub*sizeof(float));
  f_C_sub__ = (float *) malloc1(n_row_A_sub*n_row_B_sub*sizeof(float));
  GLOBAL_tic(0);
  ulli=0;
  for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ for (ncol_X=0;ncol_X<n_col_X;ncol_X++){
      f_A_trn__[ncol_X + nrow_A*n_col_X_rup] = (float)(((int)ulli%7)-3); 
      f_A_u_trn__[ncol_X + nrow_A*n_col_X] = (float)(((int)ulli%7)-3); 
      ulli++;
      /* for (nrow_A=0;nrow_A<n_row_A;nrow_A++){ for (ncol_X=0;ncol_X<n_col_X;ncol_X++){ }} */}}
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (ncol_X=0;ncol_X<n_col_X_sub;ncol_X++){ f_A_sub__[nrow_A+ncol_X*n_row_A_sub] = f_A_trn__[ncol_X + nrow_A*n_col_X_rup];}}
  printf(" %% upper corner of f_A_trn__: \n");
  array_printf(f_A_sub__,"float",n_row_A_sub,n_col_X_sub," % f_A_sub__: ");
  ulli=0;
  for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ for (ncol_X=0;ncol_X<n_col_X;ncol_X++){
      f_B_trn__[ncol_X + nrow_B*n_col_X_rup] = (float)(((int)ulli%5)-2);
      f_B_u_trn__[ncol_X + nrow_B*n_col_X] = (float)(((int)ulli%5)-2);
      ulli++;
      /* for (nrow_B=0;nrow_B<n_row_B;nrow_B++){ for (ncol_X=0;ncol_X<n_col_X;ncol_X++){ }} */}}
  for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ for (ncol_X=0;ncol_X<n_col_X_sub;ncol_X++){ f_B_sub__[nrow_B+ncol_X*n_row_B_sub] = f_B_trn__[ncol_X + nrow_B*n_col_X_rup];}}
  printf(" %% upper corner of f_B_trn__: \n");
  array_printf(f_B_sub__,"float",n_row_B_sub,n_col_X_sub," % f_B_sub__: ");
  GLOBAL_toc(0,1," initialize: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);
  memset(f_C_bf__,0,ulli_C_total*sizeof(float));
  dp_ps_mult_bruteforce(n_row_A,n_col_X,f_A_trn__,n_row_B,f_B_trn__,&f_C_bf__);
  GLOBAL_toc(0,1," dp_ps_mult_bruteforce: ");
  printf("Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/GLOBAL_elrt[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ f_C_sub__[nrow_A+nrow_B*n_row_A_sub] = f_C_bf__[nrow_A+nrow_B*n_row_A];}}
  printf(" %% upper corner of f_C_bf__: \n");
  array_printf(f_C_sub__,"float",n_row_A_sub,n_row_B_sub," % f_C_bf__: ");
  ferror = ffnorm(ulli_C_total,f_C_bf__,f_C_bf__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
  GLOBAL_tic(0);
  memset(f_C_ps__,0,ulli_C_total*sizeof(float));
  dp_ps_mult_cblas_sgemm(n_row_A,n_col_X,ps_A_trn__,n_row_B,ps_B_trn__,&f_C_ps__);
  GLOBAL_toc(0,1," dp_ps_mult_cblas_sgemm: ");
  printf("Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/GLOBAL_elrt[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ f_C_sub__[nrow_A+nrow_B*n_row_A_sub] = f_C_ps__[nrow_A+nrow_B*n_row_A];}}
  printf(" %% upper corner of f_C_ps__: \n");
  array_printf(f_C_sub__,"float",n_row_A_sub,n_row_B_sub," % f_C_ps__: ");
  ferror = ffnorm(ulli_C_total,f_C_bf__,f_C_ps__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
  GLOBAL_tic(0);
  memset(f_C_ps__,0,ulli_C_total*sizeof(float));
  dp_ps_mult_cblas_sdot(n_row_A,n_col_X,ps_A_trn__,n_row_B,ps_B_trn__,&f_C_ps__);
  GLOBAL_toc(0,1," dp_ps_mult_cblas_sdot: ");
  printf("Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/GLOBAL_elrt[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ f_C_sub__[nrow_A+nrow_B*n_row_A_sub] = f_C_ps__[nrow_A+nrow_B*n_row_A];}}
  printf(" %% upper corner of f_C_ps__: \n");
  array_printf(f_C_sub__,"float",n_row_A_sub,n_row_B_sub," % f_C_ps__: ");
  ferror = ffnorm(ulli_C_total,f_C_bf__,f_C_ps__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
  GLOBAL_tic(0);
  memset(f_C_ps__,0,ulli_C_total*sizeof(float));
  dp_ps_mult_immintrin_avx(n_row_A,n_col_X,ps_A_trn__,n_row_B,ps_B_trn__,&f_C_ps__);
  GLOBAL_toc(0,1," dp_ps_mult_immintrin_avx: ");
  printf("Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/GLOBAL_elrt[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ f_C_sub__[nrow_A+nrow_B*n_row_A_sub] = f_C_ps__[nrow_A+nrow_B*n_row_A];}}
  printf(" %% upper corner of f_C_ps__: \n");
  array_printf(f_C_sub__,"float",n_row_A_sub,n_row_B_sub," % f_C_ps__: ");
  ferror = ffnorm(ulli_C_total,f_C_bf__,f_C_ps__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
  GLOBAL_tic(0);
  memset(f_C_ps__,0,ulli_C_total*sizeof(float));
  dp_ps_mult_immintrin_fma(n_row_A,n_col_X,ps_A_trn__,n_row_B,ps_B_trn__,&f_C_ps__);
  GLOBAL_toc(0,1," dp_ps_mult_immintrin_fma: ");
  printf("Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/GLOBAL_elrt[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ f_C_sub__[nrow_A+nrow_B*n_row_A_sub] = f_C_ps__[nrow_A+nrow_B*n_row_A];}}
  printf(" %% upper corner of f_C_ps__: \n");
  array_printf(f_C_sub__,"float",n_row_A_sub,n_row_B_sub," % f_C_ps__: ");
  ferror = ffnorm(ulli_C_total,f_C_bf__,f_C_ps__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
  GLOBAL_tic(0);
  memset(f_C_ps__,0,ulli_C_total*sizeof(float));
  dp_ps_mult_immintrin_fma2(n_row_A,n_col_X,ps_A_trn__,n_row_B,ps_B_trn__,&f_C_ps__);
  GLOBAL_toc(0,1," dp_ps_mult_immintrin_fma2: ");
  printf("Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/GLOBAL_elrt[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ f_C_sub__[nrow_A+nrow_B*n_row_A_sub] = f_C_ps__[nrow_A+nrow_B*n_row_A];}}
  printf(" %% upper corner of f_C_ps__: \n");
  array_printf(f_C_sub__,"float",n_row_A_sub,n_row_B_sub," % f_C_ps__: ");
  ferror = ffnorm(ulli_C_total,f_C_bf__,f_C_ps__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
  GLOBAL_tic(0);
  memset(f_C_ps__,0,ulli_C_total*sizeof(float));
  dp_ps_mult_immintrin_load1(n_row_A,n_col_X,f_A_trn__,n_row_B,f_B_trn__,&f_C_ps__);
  GLOBAL_toc(0,1," dp_ps_mult_immintrin_load1: ");
  printf("Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X_rup/GLOBAL_elrt[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ f_C_sub__[nrow_A+nrow_B*n_row_A_sub] = f_C_ps__[nrow_A+nrow_B*n_row_A];}}
  printf(" %% upper corner of f_C_ps__: \n");
  array_printf(f_C_sub__,"float",n_row_A_sub,n_row_B_sub," % f_C_ps__: ");
  ferror = ffnorm(ulli_C_total,f_C_bf__,f_C_ps__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
  GLOBAL_tic(0);
  memset(f_C_ps__,0,ulli_C_total*sizeof(float));
  dp_ps_mult_immintrin_loadu(n_row_A,n_col_X,f_A_u_trn__,n_row_B,f_B_u_trn__,&f_C_ps__);
  GLOBAL_toc(0,1," dp_ps_mult_immintrin_loadu: ");
  printf("Gops %0.6f\n",(double)n_row_A*(double)n_row_B*(double)n_col_X/GLOBAL_elrt[0]/1e9);
  for (nrow_A=0;nrow_A<n_row_A_sub;nrow_A++){ for (nrow_B=0;nrow_B<n_row_B_sub;nrow_B++){ f_C_sub__[nrow_A+nrow_B*n_row_A_sub] = f_C_ps__[nrow_A+nrow_B*n_row_A];}}
  printf(" %% upper corner of f_C_ps__: \n");
  array_printf(f_C_sub__,"float",n_row_A_sub,n_row_B_sub," % f_C_ps__: ");
  ferror = ffnorm(ulli_C_total,f_C_bf__,f_C_ps__);
  printf(" %% ferror %0.16f\n",ferror);
  /* %%%%%%%% */
  _mm_free(ps_A_trn__); ps_A_trn__ = NULL;
  _mm_free(ps_B_trn__); ps_B_trn__ = NULL;
  free1(&f_A_u_trn__);
  free1(&f_B_u_trn__);
  free1(&f_C_bf__);
  free1(&f_C_ps__);
  free1(&f_A_sub__);
  free1(&f_B_sub__);
  free1(&f_C_sub__);
  //wkspace_printf();
}

void dp_ps_single_test()
{
  __m256 ps0 = _mm256_setr_ps((float)10.0,(float)20.0,(float)30.0,(float)40.0,(float)50.0,(float)60.0,(float)70.0,(float)80.0);
  __m256 ps1 = _mm256_setr_ps((float)11.0,(float)21.0,(float)31.0,(float)41.0,(float)51.0,(float)61.0,(float)71.0,(float)81.0);
  unsigned char uchar_mask = 0x01; //%<-- mask 01, 11, 21, 41, 81 correspond to none, first, second, third and fourth entries, respectively. ;
  /* Note also that uchar_mask must be 8 bit immediate, available as compile-time constant. */
  __m256 pso;
  float *f_ = NULL;
  printf(" %% [entering dp_ps_test]\n");
  f_ = (float *) &ps0; array_printf(f_,"float",1,8,"ps0: ");
  f_ = (float *) &ps1; array_printf(f_,"float",1,8,"ps1: ");
  uchar_mask = 0x01; bitstring_from_uchar_printf(&uchar_mask,1,BIT8,"uchar_mask: ");
  pso = _mm256_dp_ps(ps0,ps1,0x01); f_ = (float *) &pso; array_printf(f_,"float",1,8,"pso: ");
  uchar_mask = 0x11; bitstring_from_uchar_printf(&uchar_mask,1,BIT8,"uchar_mask: ");
  pso = _mm256_dp_ps(ps0,ps1,0x11); f_ = (float *) &pso; array_printf(f_,"float",1,8,"pso: ");
  uchar_mask = 0x21; bitstring_from_uchar_printf(&uchar_mask,1,BIT8,"uchar_mask: ");
  pso = _mm256_dp_ps(ps0,ps1,0x21); f_ = (float *) &pso; array_printf(f_,"float",1,8,"pso: ");
  uchar_mask = 0x41; bitstring_from_uchar_printf(&uchar_mask,1,BIT8,"uchar_mask: ");
  pso = _mm256_dp_ps(ps0,ps1,0x41); f_ = (float *) &pso; array_printf(f_,"float",1,8,"pso: ");
  uchar_mask = 0x81; bitstring_from_uchar_printf(&uchar_mask,1,BIT8,"uchar_mask: ");
  pso = _mm256_dp_ps(ps0,ps1,0x81); f_ = (float *) &pso; array_printf(f_,"float",1,8,"pso: ");
  uchar_mask = 0x31; bitstring_from_uchar_printf(&uchar_mask,1,BIT8,"uchar_mask: ");
  pso = _mm256_dp_ps(ps0,ps1,0x31); f_ = (float *) &pso; array_printf(f_,"float",1,8,"pso: ");
  printf(" %% [finished dp_ps_test]\n");
}


