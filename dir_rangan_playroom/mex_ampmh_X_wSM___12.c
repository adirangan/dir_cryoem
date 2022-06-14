/* The computational routine (some optimization) */

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

float* float_m256_malloc_from_double_(int n_a,double *d_)
{
  int na=0;
  int n_a_rup = rup(n_a,8);
  int n_a_256 = n_a_rup/8; //%<-- 8 floats per __m256. ;
  float *f_=NULL;
  f_ = (float *) _mm_malloc(n_a_256*sizeof(__m256),32);
  memset(f_,0,n_a_256*sizeof(float));
  if (d_!=NULL){
    for (na=0;na<n_a;na++){ f_[na] = (float)(d_[na]);}
  /* if (d_!=NULL){ } */}
  return f_;
}

float* float__m256_malloc_from_double__(int n_a,int n_b,double *d__)
{
  int na=0;
  int n_a_rup = rup(n_a,8);
  int n_a_256 = n_a_rup/8; //%<-- 8 floats per __m256. ;
  int nb=0;
  int n_b_rup = rup(n_b,8);
  int n_b_256 = n_b_rup/8; //%<-- 8 floats per __m256. ;
  unsigned long long int ulli=0;
  float *f__=NULL;
  ulli = (unsigned long long int)n_a_256*(unsigned long long int)n_b_rup;
  f__ = (float *) _mm_malloc(ulli*sizeof(__m256),32);
  memset(f__,0,ulli*sizeof(float));
  if (d__!=NULL){
    for (nb=0;nb<n_b;nb++){
      for (na=0;na<n_a;na++){
	f__[na + nb*n_a_rup] = (float)(d__[na + nb*n_a]);
	/* for (na=0;na<n_a;na++){  } */}
      /* for (nb=0;nb<n_b;nb++){  } */}
    /* if (d__!=NULL){ } */}
  return f__;
}

float* float___m256_malloc_from_double___(int n_a,int n_b,int n_c,double *d___)
{
  int na=0;
  int n_a_rup = rup(n_a,8);
  int n_a_256 = n_a_rup/8; //%<-- 8 floats per __m256. ;
  int nb=0;
  int n_b_rup = rup(n_b,8);
  int n_b_256 = n_b_rup/8; //%<-- 8 floats per __m256. ;
  int nc=0;
  int n_c_rup = rup(n_c,8);
  int n_c_256 = n_c_rup/8; //%<-- 8 floats per __m256. ;
  unsigned long long int ulli=0;
  float *f___=NULL;
  ulli = (unsigned long long int)n_a_256*(unsigned long long int)n_b_rup*(unsigned long long int)n_c_rup;
  f___ = (float *) _mm_malloc(ulli*sizeof(__m256),32);
  memset(f___,0,ulli*sizeof(float));
  if (d___!=NULL){
    for (nc=0;nc<n_c;nc++){
      for (nb=0;nb<n_b;nb++){
	for (na=0;na<n_a;na++){
	  f___[na + (nb + nc*n_b_rup)*n_a_rup] = (float)(d___[na + (nb + nc*n_b)*n_a]);
	  /* for (na=0;na<n_a;na++){  } */}
	/* for (nb=0;nb<n_b;nb++){  } */}
      /* for (nc=0;nc<n_c;nc++){  } */}
    /* if (d___!=NULL){ } */}
  return f___;
}

float* float____m256_malloc_from_double____(int n_a,int n_b,int n_c,int n_d,double *d____)
{
  int na=0;
  int n_a_rup = rup(n_a,8);
  int n_a_256 = n_a_rup/8; //%<-- 8 floats per __m256. ;
  int nb=0;
  int n_b_rup = rup(n_b,8);
  int n_b_256 = n_b_rup/8; //%<-- 8 floats per __m256. ;
  int nc=0;
  int n_c_rup = rup(n_c,8);
  int n_c_256 = n_c_rup/8; //%<-- 8 floats per __m256. ;
  int nd=0;
  int n_d_rup = rup(n_d,8);
  int n_d_256 = n_d_rup/8; //%<-- 8 floats per __m256. ;
  unsigned long long int ulli=0;
  float *f____=NULL;
  ulli = (unsigned long long int)n_a_256*(unsigned long long int)n_b_rup*(unsigned long long int)n_c_rup*(unsigned long long int)n_d_rup;
  f____ = (float *) _mm_malloc(ulli*sizeof(__m256),32);
  memset(f____,0,ulli*sizeof(float));
  if (d____!=NULL){
    for (nd=0;nd<n_d;nd++){
      for (nc=0;nc<n_c;nc++){
	for (nb=0;nb<n_b;nb++){
	  for (na=0;na<n_a;na++){
	    f____[na + (nb + (nc + nd*n_c_rup)*n_b_rup)*n_a_rup] = (float)(d____[na + (nb + (nc + nd*n_c)*n_b)*n_a]);
	    /* for (na=0;na<n_a;na++){  } */}
	  /* for (nb=0;nb<n_b;nb++){  } */}
	/* for (nc=0;nc<n_c;nc++){  } */}
      /* for (nd=0;nd<n_d;nd++){  } */}
    /* if (d____!=NULL){ } */}
  return f____;
}

void mex_ampmh_X_wSM___12
(
  int n_M_per_Mbatch
 ,int n_S_per_Sbatch
 ,int flag_optimize_over_gamma_z
 ,int flag_compute_I_value
 ,double tolerance_master
 ,int FTK_n_svd_l
 ,int FTK_n_delta_v
 ,double *FTK_svd_U_d_expiw_s_real__
 ,double *FTK_svd_U_d_expiw_s_imag__
 ,double *FTK_delta_x_
 ,double *FTK_delta_y_
 ,int n_w_max
 ,int pm_n_UX_rank
 ,int n_S
 ,double *CTF_UX_S_k_q_wnS_real__
 ,double *CTF_UX_S_k_q_wnS_imag__
 ,double *CTF_UX_S_l2_
 ,int n_M
 ,double *svd_VUXM_lwnM_real____
 ,double *svd_VUXM_lwnM_imag____
 ,double *UX_M_l2_dM__
 ,double *X_wSM___
 ,double *delta_x_wSM___
 ,double *delta_y_wSM___
 ,double *gamma_z_wSM___
 ,double *I_value_wSM___
 )
{
  int verbose=1;
  int flag_disp=0;
  int flag_dump=0;
  /* %%%%; */
  int n_M_per_Mbatch_rup=0,n_M_per_Mbatch_256=0;
  int n_S_per_Sbatch_rup=0,n_S_per_Sbatch_256=0;
  int FTK_n_svd_l_rup=0,FTK_n_svd_l_256=0;
  int FTK_n_delta_v_rup=0,FTK_n_delta_v_256=0;
  int n_w_max_rup=0,n_w_max_256=0;
  int pm_n_UX_rank_rup=0,pm_n_UX_rank_256=0;
  int n_S_rup=0,n_S_256=0;
  int n_S_sub_rup=0,n_S_sub_256=0;
  int n_M_rup=0,n_M_256=0;
  int n_M_sub_rup=0,n_M_sub_256=0;
  /* %%%%; */
  double complex *FTK_svd_U_d_expiw_s__=NULL;
  float *f_FTK_svd_U_d_expiw_s_real__=NULL;
  float *f_FTK_svd_U_d_expiw_s_imag__=NULL;
  float *f_FTK_svd_U_d_expiw_s_tran_real__=NULL;
  float *f_FTK_svd_U_d_expiw_s_tran_imag__=NULL;
  double complex *CTF_UX_S_k_q_wnS__=NULL;
  float *f_CTF_UX_S_k_q_wnS_real___=NULL;
  float *f_CTF_UX_S_k_q_wnS_imag___=NULL;
  float *f_CTF_UX_S_k_q_nSw_real___=NULL;
  float *f_CTF_UX_S_k_q_nSw_imag___=NULL;
  double complex *svd_VUXM_lwnM____=NULL;
  float *f_svd_VUXM_lwnM_real____=NULL;
  float *f_svd_VUXM_lwnM_imag____=NULL;
  float *f_svd_VUXM_nMlw_real____=NULL;
  float *f_svd_VUXM_nMlw_imag____=NULL;
  unsigned long long int na=0,n_a=0;
  /* %%%% */
  int nMbatch=0,n_Mbatch=0,nSbatch=0,n_Sbatch=0;
  int *index_M_in_Mbatch_=NULL;
  int *index_S_in_Sbatch_=NULL;
  int n_M_sub=0,nM_sub=0;
  int n_S_sub=0,nS_sub=0;
  int nl=0,ndelta_v=0,nw=0,ndw=0,nS=0,nM=0;
  int pm_nUX_rank=0;
  unsigned long long int tab=0,tabA=0,tabB=0,tabC=0;
  double complex *CTF_UX_S_k_q_nSw___=NULL;
  double complex *svd_VUXM_nMwl____=NULL;
  double complex *svd_VUXM_nMlw____=NULL;
  double complex *svd_SVUXM_SMwl____=NULL;
  double complex *svd_SVUXM_SMlw____=NULL;
  float *f_svd_SVUXM_SMlw_real____=NULL;
  float *f_svd_SVUXM_SMlw_imag____=NULL;
  float *f_svd_SVUXM_wSMl_real____=NULL;
  float *f_svd_SVUXM_wSMl_imag____=NULL;
  fftw_complex *svd_SVUXM_0in_wlSM____;
  fftw_complex *svd_SVUXM_out_wlSM____;
  fftwf_complex *f_svd_SVUXM_0in_wSMl____;
  fftwf_complex *f_svd_SVUXM_out_wSMl____;
  double complex *svd_SVUXM_lwSM____=NULL;
  double complex *svd_SVUXM_lwsM____=NULL;
  float *f_svd_SVUXM_lwSM_real____=NULL;
  float *f_svd_SVUXM_lwSM_imag____=NULL;
  float *f_svd_SVUXM_lwsM_real____=NULL;
  float *f_svd_SVUXM_lwsM_imag____=NULL;
  double complex *svd_USESVUXM_dwSM____=NULL;
  float *f_svd_USESVUXM_dwSM_real____=NULL;
  float *f_svd_USESVUXM_dwSM_imag____=NULL;
  double complex cblas_alpha = (double complex) 1.0;
  double complex cblas_beta  = (double complex) 0.0;
  double *l2_dSM___=NULL;
  double *n2_dSM___=NULL;
  double *f2_dSM___=NULL;
  double *X_sub_dwSM____=NULL;
  double *d_X_sub_dwSM____=NULL;
  double *I_value_use_dwSM____=NULL;
  double *gamma_z_=NULL;
  float *f_CR_=NULL,*f_CI_=NULL,*f_C_=NULL;
  /* %%%% */
  fftw_plan fftw_plan_many_plan;
  int fftw_plan_many_rank=0;
  int fftw_plan_many_n[1]={0};
  int fftw_plan_many_howmany=0;
  int fftw_plan_many_istride=0;
  int fftw_plan_many_idist=0;
  int fftw_plan_many_ostride=0;
  int fftw_plan_many_odist=0;
  int fftw_plan_many_sign=0;
  unsigned int fftw_plan_many_flags=0;
  /* %%%% */
  fftwf_plan fftwf_plan_many_plan;
  int fftwf_plan_many_rank=0;
  int fftwf_plan_many_n[1]={0};
  int fftwf_plan_many_howmany=0;
  int fftwf_plan_many_istride=0;
  int fftwf_plan_many_idist=0;
  int fftwf_plan_many_ostride=0;
  int fftwf_plan_many_odist=0;
  int fftwf_plan_many_sign=0;
  unsigned int fftwf_plan_many_flags=0;
  /* %%%% */
  int MDA_n_dim=0,MDA_ndim = 0;
  int MDA_dim_[8];
  char MDA_fname[32];
  /* %%%% */
  int ntick=0,n_tick = 8;
  clock_t t_start_[n_tick];
  clock_t t_final_[n_tick];
  struct timeval d_start_[n_tick];
  struct timeval d_final_[n_tick];
  long l_msec_[n_tick],l_ssec_[n_tick],l_usec_[n_tick];
  double elct_[n_tick],elrt_[n_tick];
  double elrt_sum_CTF_UX_S_k_q_nSw___=0;
  double n_op_sum_CTF_UX_S_k_q_nSw___=0;
  double elrt_sum_f_CTF_UX_S_k_q_nSw___=0;
  double n_op_sum_f_CTF_UX_S_k_q_nSw___=0;
  double elrt_sum_svd_VUXM_nMwl____=0;
  double n_op_sum_svd_VUXM_nMwl____=0;
  double elrt_sum_f_svd_VUXM_nMlw____=0;
  double n_op_sum_f_svd_VUXM_nMlw____=0;
  double elrt_sum_svd_SVUXM_SMwl____=0;
  double n_op_sum_svd_SVUXM_SMwl____=0;
  double elrt_sum_svd_SVUXM_SMlw____=0;
  double n_op_sum_svd_SVUXM_SMlw____=0;
  double elrt_sum_f_svd_SVUXM_SMlw____=0;
  double n_op_sum_f_svd_SVUXM_SMlw____=0;
  double elrt_sum_svd_SVUXM_0in_wlSM____=0;
  double n_op_sum_svd_SVUXM_0in_wlSM____=0;
  double elrt_sum_f_svd_SVUXM_0in_wSMl____=0;
  double n_op_sum_f_svd_SVUXM_0in_wSMl____=0;
  double elrt_sum_svd_SVUXM_out_wlSM____=0;
  double n_op_sum_svd_SVUXM_out_wlSM____=0;
  double elrt_sum_f_svd_SVUXM_out_wSMl____=0;
  double n_op_sum_f_svd_SVUXM_out_wSMl____=0;
  double elrt_sum_svd_SVUXM_lwSM____=0;
  double n_op_sum_svd_SVUXM_lwSM____=0;
  double elrt_sum_f_svd_SVUXM_lwSM____=0;
  double n_op_sum_f_svd_SVUXM_lwSM____=0;
  double elrt_sum_svd_SVUXM_lwsM____=0;
  double n_op_sum_svd_SVUXM_lwsM____=0;
  double elrt_sum_f_svd_SVUXM_lwsM____=0;
  double n_op_sum_f_svd_SVUXM_lwsM____=0;
  double elrt_sum_svd_USESVUXM_dwSM____=0;
  double n_op_sum_svd_USESVUXM_dwSM____=0;
  double elrt_sum_f_svd_USESVUXM_dwSM____=0;
  double n_op_sum_f_svd_USESVUXM_dwSM____=0;
  double elrt_sum_l2_dSM___=0;
  double n_op_sum_l2_dSM___=0;
  double elrt_sum_X_sub_dwSM____=0;
  double n_op_sum_X_sub_dwSM____=0;
  double elrt_sum_d_X_sub_dwSM____=0;
  double n_op_sum_d_X_sub_dwSM____=0;
  double elrt_sum_I_value_use_dwSM____=0;
  double n_op_sum_I_value_use_dwSM____=0;
  double elrt_sum_X_wSM___=0;
  double n_op_sum_X_wSM___=0;
  /* %%%% */
  if (verbose){ printf(" %% [entering mex_ampmh_X_wSM___12]\n");}
  /* %%%%%%%%%%%%%%%% */
  n_a = FTK_n_delta_v*FTK_n_svd_l;
  FTK_svd_U_d_expiw_s__ = double_complex_malloc_and_interleave(n_a,FTK_svd_U_d_expiw_s_real__,FTK_svd_U_d_expiw_s_imag__);
  n_a = n_w_max*pm_n_UX_rank*n_S;
  CTF_UX_S_k_q_wnS__ = double_complex_malloc_and_interleave(n_a,CTF_UX_S_k_q_wnS_real__,CTF_UX_S_k_q_wnS_imag__);
  n_a = FTK_n_svd_l*n_w_max*pm_n_UX_rank*n_M;
  svd_VUXM_lwnM____ = double_complex_malloc_and_interleave(n_a,svd_VUXM_lwnM_real____,svd_VUXM_lwnM_imag____);
  /* %%%%%%%%%%%%%%%% */
  FTK_n_delta_v_rup = rup(FTK_n_delta_v,8); FTK_n_delta_v_256 = FTK_n_delta_v_rup/8;
  FTK_n_svd_l_rup = rup(FTK_n_svd_l,8); FTK_n_svd_l_256 = FTK_n_svd_l_rup/8;
  n_w_max_rup = rup(n_w_max,8); n_w_max_256 = n_w_max_rup/8;
  pm_n_UX_rank_rup = rup(pm_n_UX_rank,8); pm_n_UX_rank_256 = pm_n_UX_rank_rup/8;
  n_S_rup = rup(n_S,8); n_S_256 = n_S_rup/8;
  n_S_per_Sbatch_rup = rup(n_S_per_Sbatch,8); n_S_per_Sbatch_256 = n_S_per_Sbatch_rup/8;
  n_M_rup = rup(n_M,8); n_M_256 = n_M_rup/8;
  n_M_per_Mbatch_rup = rup(n_M_per_Mbatch,8); n_M_per_Mbatch_256 = n_M_per_Mbatch_rup/8;
  f_FTK_svd_U_d_expiw_s_real__ = float__m256_malloc_from_double__(FTK_n_delta_v,FTK_n_svd_l,FTK_svd_U_d_expiw_s_real__);
  f_FTK_svd_U_d_expiw_s_imag__ = float__m256_malloc_from_double__(FTK_n_delta_v,FTK_n_svd_l,FTK_svd_U_d_expiw_s_imag__);
  f_FTK_svd_U_d_expiw_s_tran_real__ = float__m256_malloc_from_double__(FTK_n_svd_l,FTK_n_delta_v,NULL);
  f_FTK_svd_U_d_expiw_s_tran_imag__ = float__m256_malloc_from_double__(FTK_n_svd_l,FTK_n_delta_v,NULL);
  transpose_ps_block1(
		      f_FTK_svd_U_d_expiw_s_real__
		      ,f_FTK_svd_U_d_expiw_s_tran_real__
		      ,FTK_n_delta_v_rup
		      ,FTK_n_svd_l_rup
		      ,FTK_n_delta_v_rup
		      ,FTK_n_svd_l_rup
		      ,16
		      );
  transpose_ps_block1(
		      f_FTK_svd_U_d_expiw_s_imag__
		      ,f_FTK_svd_U_d_expiw_s_tran_imag__
		      ,FTK_n_delta_v_rup
		      ,FTK_n_svd_l_rup
		      ,FTK_n_delta_v_rup
		      ,FTK_n_svd_l_rup
		      ,16
		      );
  if (verbose>2){
    array_sub_printf(f_FTK_svd_U_d_expiw_s_real__     ,"float" ,FTK_n_delta_v_rup,3,FTK_n_svd_l_rup,4," %% f_FTK_svd_U_d_expiw_s_real__     : ");
    array_sub_printf(f_FTK_svd_U_d_expiw_s_tran_real__,"float" ,FTK_n_svd_l_rup,4,FTK_n_delta_v_rup,3," %% f_FTK_svd_U_d_expiw_s_tran_real__: ");
    array_sub_printf(f_FTK_svd_U_d_expiw_s_imag__     ,"float" ,FTK_n_delta_v_rup,3,FTK_n_svd_l_rup,4," %% f_FTK_svd_U_d_expiw_s_imag__     : ");
    array_sub_printf(f_FTK_svd_U_d_expiw_s_tran_imag__,"float" ,FTK_n_svd_l_rup,4,FTK_n_delta_v_rup,3," %% f_FTK_svd_U_d_expiw_s_tran_imag__: ");
    exit(0);
    /* if (verbose>2){ } */}
  f_CTF_UX_S_k_q_wnS_real___ = float___m256_malloc_from_double___(n_w_max,pm_n_UX_rank,n_S,CTF_UX_S_k_q_wnS_real__);
  f_CTF_UX_S_k_q_wnS_imag___ = float___m256_malloc_from_double___(n_w_max,pm_n_UX_rank,n_S,CTF_UX_S_k_q_wnS_imag__);
  f_CTF_UX_S_k_q_nSw_real___ = float___m256_malloc_from_double___(pm_n_UX_rank,n_S,n_w_max,NULL);
  f_CTF_UX_S_k_q_nSw_imag___ = float___m256_malloc_from_double___(pm_n_UX_rank,n_S,n_w_max,NULL);
  f_svd_VUXM_lwnM_real____ = float____m256_malloc_from_double____(FTK_n_svd_l,n_w_max,pm_n_UX_rank,n_M,svd_VUXM_lwnM_real____);
  f_svd_VUXM_lwnM_imag____ = float____m256_malloc_from_double____(FTK_n_svd_l,n_w_max,pm_n_UX_rank,n_M,svd_VUXM_lwnM_imag____);
  f_svd_VUXM_nMlw_real____ = float____m256_malloc_from_double____(FTK_n_svd_l,n_w_max,pm_n_UX_rank,n_M_per_Mbatch,NULL);
  f_svd_VUXM_nMlw_imag____ = float____m256_malloc_from_double____(FTK_n_svd_l,n_w_max,pm_n_UX_rank,n_M_per_Mbatch,NULL);
  if (verbose>2){
    /* %%%% */
    array_sub_printf(  FTK_svd_U_d_expiw_s_real__,"double",FTK_n_delta_v    ,3,FTK_n_svd_l    ,4," %%   FTK_svd_U_d_expiw_s_real__: ");
    array_sub_printf(f_FTK_svd_U_d_expiw_s_real__,"float" ,FTK_n_delta_v_rup,3,FTK_n_svd_l_rup,4," %% f_FTK_svd_U_d_expiw_s_real__: ");
    array_sub_printf(  FTK_svd_U_d_expiw_s_imag__,"double",FTK_n_delta_v    ,3,FTK_n_svd_l    ,4," %%   FTK_svd_U_d_expiw_s_imag__: ");
    array_sub_printf(f_FTK_svd_U_d_expiw_s_imag__,"float" ,FTK_n_delta_v_rup,3,FTK_n_svd_l_rup,4," %% f_FTK_svd_U_d_expiw_s_imag__: ");
    /* %%%% */
    array_sub_printf(  CTF_UX_S_k_q_wnS_real__ ,"double",n_w_max    *pm_n_UX_rank    ,3,n_S    ,4," %%   CTF_UX_S_k_q_wnS_real__ : ");
    array_sub_printf(f_CTF_UX_S_k_q_wnS_real___,"float" ,n_w_max_rup*pm_n_UX_rank_rup,3,n_S_rup,4," %% f_CTF_UX_S_k_q_wnS_real___: ");
    array_sub_printf(  CTF_UX_S_k_q_wnS_imag__ ,"double",n_w_max    *pm_n_UX_rank    ,3,n_S    ,4," %%   CTF_UX_S_k_q_wnS_imag__ : ");
    array_sub_printf(f_CTF_UX_S_k_q_wnS_imag___,"float" ,n_w_max_rup*pm_n_UX_rank_rup,3,n_S_rup,4," %% f_CTF_UX_S_k_q_wnS_imag___: ");
    /* %%%% */
    array_sub_printf(  svd_VUXM_lwnM_real____,"double",FTK_n_svd_l    *n_w_max    *pm_n_UX_rank    ,3,n_M    ,4," %%   svd_VUXM_lwnM_real____: ");
    array_sub_printf(f_svd_VUXM_lwnM_real____,"float" ,FTK_n_svd_l_rup*n_w_max_rup*pm_n_UX_rank_rup,3,n_M_rup,4," %% f_svd_VUXM_lwnM_real____: ");
    array_sub_printf(  svd_VUXM_lwnM_imag____,"double",FTK_n_svd_l    *n_w_max    *pm_n_UX_rank    ,3,n_M    ,4," %%   svd_VUXM_lwnM_imag____: ");
    array_sub_printf(f_svd_VUXM_lwnM_imag____,"float" ,FTK_n_svd_l_rup*n_w_max_rup*pm_n_UX_rank_rup,3,n_M_rup,4," %% f_svd_VUXM_lwnM_imag____: ");
    /* %%%% */
    array_sub_printf(  svd_VUXM_lwnM_real____,"double",FTK_n_svd_l    *n_w_max    ,3,pm_n_UX_rank    *n_M    ,4," %%   svd_VUXM_lwnM_real____: ");
    array_sub_printf(f_svd_VUXM_lwnM_real____,"float" ,FTK_n_svd_l_rup*n_w_max_rup,3,pm_n_UX_rank_rup*n_M_rup,4," %% f_svd_VUXM_lwnM_real____: ");
    array_sub_printf(  svd_VUXM_lwnM_imag____,"double",FTK_n_svd_l    *n_w_max    ,3,pm_n_UX_rank    *n_M    ,4," %%   svd_VUXM_lwnM_imag____: ");
    array_sub_printf(f_svd_VUXM_lwnM_imag____,"float" ,FTK_n_svd_l_rup*n_w_max_rup,3,pm_n_UX_rank_rup*n_M_rup,4," %% f_svd_VUXM_lwnM_imag____: ");
    /* %%%% */
    exit(0);
    /*if (verbose>2){ } */}
  /* %%%%%%%%%%%%%%%% */
  if (verbose>1){ printf(" %% n_M_per_Mbatch %d\n",n_M_per_Mbatch);}
  if (verbose>1){ printf(" %% n_S_per_Sbatch %d\n",n_S_per_Sbatch);}
  if (verbose>1){ printf(" %% flag_optimize_over_gamma_z %d\n",flag_optimize_over_gamma_z);}
  if (verbose>1){ printf(" %% flag_compute_I_value %d\n",flag_compute_I_value);}
  if (verbose>1){ printf(" %% tolerance_master %0.16f\n",tolerance_master);}
  if (verbose>1){ printf(" %% FTK_n_svd_l %d\n",FTK_n_svd_l);}
  if (verbose>1){ printf(" %% FTK_n_delta_v %d\n",FTK_n_delta_v);}
  n_a = FTK_n_delta_v*FTK_n_svd_l;
  if (verbose>1){ printf(" %% FTK_svd_U_d_expiw_s__ %+0.16f %+0.16f*i --> %+0.16f %+0.16f*i\n",creal(FTK_svd_U_d_expiw_s__[0]),cimag(FTK_svd_U_d_expiw_s__[0]),creal(FTK_svd_U_d_expiw_s__[n_a-1]),cimag(FTK_svd_U_d_expiw_s__[n_a-1]));}
  if (verbose>1){ printf(" %% FTK_delta_x_ %0.16f --> %0.16f\n",FTK_delta_x_[0],FTK_delta_x_[FTK_n_delta_v-1]);}
  if (verbose>1){ printf(" %% FTK_delta_y_ %0.16f --> %0.16f\n",FTK_delta_y_[0],FTK_delta_y_[FTK_n_delta_v-1]);}
  if (verbose>1){ printf(" %% n_w_max %d\n",n_w_max);}
  if (verbose>1){ printf(" %% pm_n_UX_rank %d\n",pm_n_UX_rank);}
  if (verbose>1){ printf(" %% n_S %d\n",n_S);}
  n_a = n_w_max*pm_n_UX_rank*n_S;
  if (verbose>1){ printf(" %% CTF_UX_S_k_q_wnS__ %+0.16f %+0.16f*i --> %+0.16f %+0.16f*i\n",creal(CTF_UX_S_k_q_wnS__[0]),cimag(CTF_UX_S_k_q_wnS__[0]),creal(CTF_UX_S_k_q_wnS__[n_a-1]),cimag(CTF_UX_S_k_q_wnS__[n_a-1]));}
  if (verbose>1){ printf(" %% CTF_UX_S_l2_ %0.16f --> %0.16f\n",CTF_UX_S_l2_[0],CTF_UX_S_l2_[n_S-1]);}
  if (verbose>1){ printf(" %% n_M %d\n",n_M);}
  n_a = FTK_n_svd_l*n_w_max*pm_n_UX_rank*n_M;
  if (verbose>1){ printf(" %% svd_VUXM_lwnM____ %+0.16f %+0.16f*i --> %+0.16f %+0.16f*i\n",creal(svd_VUXM_lwnM____[0]),cimag(svd_VUXM_lwnM____[0]),creal(svd_VUXM_lwnM____[n_a-1]),cimag(svd_VUXM_lwnM____[n_a-1]));}
  if (verbose>1){ printf(" %% UX_M_l2_dM__ %0.16f --> %0.16f\n",UX_M_l2_dM__[0],UX_M_l2_dM__[FTK_n_delta_v*n_M-1]);}
  /* %%%% */
  /* if ( (nargout>5) & (isempty(I_value_wSM___)) ); I_value_wSM___ = ones(size(X_wSM___)); end; */
  /* %%%% */
  if (flag_compute_I_value){
    if (flag_optimize_over_gamma_z==0){ n_a = (unsigned long long int)n_M*(unsigned long long int)n_S*(unsigned long long int)n_w_max;}
    if (flag_optimize_over_gamma_z==1){ n_a = (unsigned long long int)n_M*(unsigned long long int)n_S;}
    for (na=0;na<n_a;na++){ I_value_wSM___[na] = 1.0;}
    /* if (flag_compute_I_value){ } */}
  /* %%%% */
  /* CTF_UX_S_k_q_nSw___ = permute(CTF_UX_S_k_q_wnS___,[2,3,1]); */
  /* %%%% */
  local_tic(0,t_start_,d_start_);
  tab = (unsigned long long int)pm_n_UX_rank*(unsigned long long int)n_S*(unsigned long long int)n_w_max;
  CTF_UX_S_k_q_nSw___ = (double complex *) malloc(tab*sizeof(double complex));
  na=0;
  for (nS=0;nS<n_S;nS++){
    for (pm_nUX_rank=0;pm_nUX_rank<pm_n_UX_rank;pm_nUX_rank++){
      for (nw=0;nw<n_w_max;nw++){
	/* tabA = */
	/*   (unsigned long long int)pm_nUX_rank */
	/*   + (unsigned long long int)nS*(unsigned long long int)pm_n_UX_rank */
	/*   + (unsigned long long int)nw*(unsigned long long int)n_S*(unsigned long long int)pm_n_UX_rank; */
	tabA = (unsigned long long int)pm_nUX_rank + ((unsigned long long int)nS + (unsigned long long int)nw*(unsigned long long int)n_S)*(unsigned long long int)pm_n_UX_rank;
	CTF_UX_S_k_q_nSw___[tabA] = CTF_UX_S_k_q_wnS__[na]; na++;
	/* for (nw=0;nw<n_w_max;nw++){ } */}
      /* for (pm_nUX_rank=0;pm_nUX_rank<pm_n_UX_rank;pm_nUX_rank++){ } */}
    /* for (nS=0;nS<n_S;nS++){ } */}
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"CTF_UX_S_k_q_nSw___: ");
  elrt_sum_CTF_UX_S_k_q_nSw___ += elrt_[0];
  n_op_sum_CTF_UX_S_k_q_nSw___ += tab;
  /* %%%%%%%%%%%%%%%% */
  local_tic(0,t_start_,d_start_);
  tab = (unsigned long long int)pm_n_UX_rank*(unsigned long long int)n_S*(unsigned long long int)n_w_max;
  transpose_ps_block1(
		      f_CTF_UX_S_k_q_wnS_real___
		      ,f_CTF_UX_S_k_q_nSw_real___
		      ,n_w_max_rup
		      ,pm_n_UX_rank_rup*n_S_rup
		      ,n_w_max_rup
		      ,pm_n_UX_rank_rup*n_S_rup
		      ,16
		      );
  transpose_ps_block1(
		      f_CTF_UX_S_k_q_wnS_imag___
		      ,f_CTF_UX_S_k_q_nSw_imag___
		      ,n_w_max_rup
		      ,pm_n_UX_rank_rup*n_S_rup
		      ,n_w_max_rup
		      ,pm_n_UX_rank_rup*n_S_rup
		      ,16
		      );
  if (verbose>2){
    array_sub_printf(  CTF_UX_S_k_q_nSw___ ,"double complex",pm_n_UX_rank    *n_S    ,4,n_w_max    ,3," %%   CTF_UX_S_k_q_nSw___ : ");
    array_sub_printf(f_CTF_UX_S_k_q_nSw_real___,"float" ,pm_n_UX_rank_rup*n_S_rup,4,n_w_max_rup,3," %% f_CTF_UX_S_k_q_nSw_real___: ");
    array_sub_printf(f_CTF_UX_S_k_q_nSw_imag___,"float" ,pm_n_UX_rank_rup*n_S_rup,4,n_w_max_rup,3," %% f_CTF_UX_S_k_q_nSw_imag___: ");
    exit(0);
    /* if (verbose>2){ } */}  
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"f_CTF_UX_S_k_q_nSw___: ");
  elrt_sum_f_CTF_UX_S_k_q_nSw___ += elrt_[0];
  n_op_sum_f_CTF_UX_S_k_q_nSw___ += tab;
  /* %%%%%%%%%%%%%%%% */
  tab = (unsigned long long int)pm_n_UX_rank*(unsigned long long int)n_M_per_Mbatch*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_svd_l;
  svd_VUXM_nMwl____ = (double complex *) malloc(tab*sizeof(double complex));
  svd_VUXM_nMlw____ = (double complex *) malloc(tab*sizeof(double complex));
  tab = (unsigned long long int)n_S*(unsigned long long int)n_M_per_Mbatch*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_svd_l;
  svd_SVUXM_SMwl____ = (double complex *) malloc(tab*sizeof(double complex));
  svd_SVUXM_SMlw____ = (double complex *) malloc(tab*sizeof(double complex));
  f_svd_SVUXM_SMlw_real____ = float____m256_malloc_from_double____(n_S,n_M_per_Mbatch,n_w_max,FTK_n_svd_l,NULL);
  f_svd_SVUXM_SMlw_imag____ = float____m256_malloc_from_double____(n_S,n_M_per_Mbatch,n_w_max,FTK_n_svd_l,NULL);
  f_svd_SVUXM_wSMl_real____ = float____m256_malloc_from_double____(n_S,n_M_per_Mbatch,n_w_max,FTK_n_svd_l,NULL);
  f_svd_SVUXM_wSMl_imag____ = float____m256_malloc_from_double____(n_S,n_M_per_Mbatch,n_w_max,FTK_n_svd_l,NULL);
  /* %%%% */
  tab = (unsigned long long int)n_w_max*(unsigned long long int)FTK_n_svd_l*(unsigned long long int)n_S*(unsigned long long int)n_M_per_Mbatch;
  svd_SVUXM_0in_wlSM____ = (fftw_complex *) fftw_malloc(tab*sizeof(fftw_complex));
  svd_SVUXM_out_wlSM____ = (fftw_complex *) fftw_malloc(tab*sizeof(fftw_complex));
  fftw_plan_many_rank = 1;
  fftw_plan_many_n[0] = n_w_max;
  fftw_plan_many_howmany = (unsigned long long int)FTK_n_svd_l*(unsigned long long int)n_S*(unsigned long long int)n_M_per_Mbatch;
  fftw_plan_many_istride = 1;
  fftw_plan_many_idist = n_w_max;
  fftw_plan_many_ostride = 1;
  fftw_plan_many_odist = n_w_max;
  fftw_plan_many_sign = FFTW_BACKWARD;
  fftw_plan_many_flags = FFTW_ESTIMATE;
  fftw_plan_many_plan = fftw_plan_many_dft(
					   fftw_plan_many_rank
					   ,fftw_plan_many_n
					   ,fftw_plan_many_howmany
					   ,svd_SVUXM_0in_wlSM____
					   ,NULL
					   ,fftw_plan_many_istride
					   ,fftw_plan_many_idist
					   ,svd_SVUXM_out_wlSM____
					   ,NULL
					   ,fftw_plan_many_ostride
					   ,fftw_plan_many_odist
					   ,fftw_plan_many_sign
					   ,fftw_plan_many_flags
					   );
  /* %%%% */
  tab = (unsigned long long int)n_w_max_rup*(unsigned long long int)n_S_rup*(unsigned long long int)n_M_per_Mbatch_rup*(unsigned long long int)FTK_n_svd_l_rup;
  f_svd_SVUXM_0in_wSMl____ = (fftwf_complex *) fftwf_malloc(tab*sizeof(fftwf_complex));
  f_svd_SVUXM_out_wSMl____ = (fftwf_complex *) fftwf_malloc(tab*sizeof(fftwf_complex));
  fftwf_plan_many_rank = 1;
  fftwf_plan_many_n[0] = n_w_max;
  fftwf_plan_many_howmany = (unsigned long long int)n_S_rup*(unsigned long long int)n_M_per_Mbatch_rup*(unsigned long long int)FTK_n_svd_l_rup;
  fftwf_plan_many_istride = 1;
  fftwf_plan_many_idist = n_w_max_rup;
  fftwf_plan_many_ostride = 1;
  fftwf_plan_many_odist = n_w_max_rup;
  fftwf_plan_many_sign = FFTW_BACKWARD;
  fftwf_plan_many_flags = FFTW_ESTIMATE;
  fftwf_plan_many_plan = fftwf_plan_many_dft(
					   fftwf_plan_many_rank
					   ,fftwf_plan_many_n
					   ,fftwf_plan_many_howmany
					   ,f_svd_SVUXM_0in_wSMl____
					   ,NULL
					   ,fftwf_plan_many_istride
					   ,fftwf_plan_many_idist
					   ,f_svd_SVUXM_out_wSMl____
					   ,NULL
					   ,fftwf_plan_many_ostride
					   ,fftwf_plan_many_odist
					   ,fftwf_plan_many_sign
					   ,fftwf_plan_many_flags
					   );
  /* %%%% */
  tab = (unsigned long long int)FTK_n_svd_l*(unsigned long long int)n_w_max*(unsigned long long int)n_S*(unsigned long long int)n_M_per_Mbatch;
  svd_SVUXM_lwSM____ = (double complex *) malloc(tab*sizeof(double complex));
  f_svd_SVUXM_lwSM_real____ = float____m256_malloc_from_double____(FTK_n_svd_l,n_w_max,n_S,n_M_per_Mbatch,NULL);
  f_svd_SVUXM_lwSM_imag____ = float____m256_malloc_from_double____(FTK_n_svd_l,n_w_max,n_S,n_M_per_Mbatch,NULL);
  tab = (unsigned long long int)n_S_per_Sbatch*(unsigned long long int)n_M_per_Mbatch*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_svd_l;
  svd_SVUXM_lwsM____ = (double complex *) malloc(tab*sizeof(double complex));
  f_svd_SVUXM_lwsM_real____ = float____m256_malloc_from_double____(FTK_n_svd_l,n_w_max,n_S_per_Sbatch,n_M_per_Mbatch,NULL);
  f_svd_SVUXM_lwsM_imag____ = float____m256_malloc_from_double____(FTK_n_svd_l,n_w_max,n_S_per_Sbatch,n_M_per_Mbatch,NULL);
  tab = (unsigned long long int)n_S_per_Sbatch*(unsigned long long int)n_M_per_Mbatch*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_delta_v;
  svd_USESVUXM_dwSM____ = (double complex *) malloc(tab*sizeof(double complex));
  f_svd_USESVUXM_dwSM_real____ = float____m256_malloc_from_double____(FTK_n_delta_v,n_w_max,n_S_per_Sbatch,n_M_per_Mbatch,NULL);
  f_svd_USESVUXM_dwSM_imag____ = float____m256_malloc_from_double____(FTK_n_delta_v,n_w_max,n_S_per_Sbatch,n_M_per_Mbatch,NULL);
  tab = (unsigned long long int)n_S_per_Sbatch*(unsigned long long int)n_M_per_Mbatch*(unsigned long long int)FTK_n_delta_v;
  l2_dSM___ = (double *) malloc(tab*sizeof(double));
  n2_dSM___ = (double *) malloc(tab*sizeof(double));
  f2_dSM___ = (double *) malloc(tab*sizeof(double));
  tab = (unsigned long long int)n_S_per_Sbatch*(unsigned long long int)n_M_per_Mbatch*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_delta_v;
  X_sub_dwSM____ = (double *) malloc(tab*sizeof(double));
  d_X_sub_dwSM____ = (double *) malloc(tab*sizeof(double));
  if (flag_compute_I_value){ I_value_use_dwSM____ = (double *) malloc(tab*sizeof(double));}
  tab = (unsigned long long int)n_w_max;
  gamma_z_ = (double *) malloc(tab*sizeof(double));
  for (nw=0;nw<n_w_max;nw++){ gamma_z_[nw] = 2*PI*(double)nw/(double)n_w_max;}
  /* %%%%%%%%%%%%%%%% */
  n_Mbatch = ceil((double)n_M/(double)n_M_per_Mbatch);
  if (verbose>1){ printf(" %% n_Mbatch %d\n",n_Mbatch);}
  index_M_in_Mbatch_ = (int *) malloc(n_M_per_Mbatch*sizeof(int));
  n_Sbatch = ceil((double)n_S/(double)n_S_per_Sbatch);
  if (verbose>1){ printf(" %% n_Sbatch %d\n",n_Sbatch);}
  index_S_in_Sbatch_ = (int *) malloc(n_S_per_Sbatch*sizeof(int));
  for (nMbatch=0;nMbatch<n_Mbatch;nMbatch++){
    n_M_sub = n_M_per_Mbatch; if (nMbatch==n_Mbatch-1){ n_M_sub = n_M - n_M_per_Mbatch*nMbatch;}
    for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ index_M_in_Mbatch_[nM_sub] = n_M_per_Mbatch*nMbatch + nM_sub;}
    if (verbose>0){ printf(" %% nMbatch %d/%d: index_M_in_Mbatch_ %d --> %d\n",nMbatch,n_Mbatch,index_M_in_Mbatch_[0],index_M_in_Mbatch_[n_M_sub-1]);}
    if (n_M_sub>0){
      n_M_sub_rup = rup(n_M_sub,8);
      n_M_sub_256 = n_M_sub_rup/8;
      /* %%%% */
      /* svd_VUXM_nMwl____ = permute(svd_VUXM_lwnM____(:,:,:,1+index_M_in_Mbatch_),[3,4,2,1]); */
      /* %%%% */      
      local_tic(0,t_start_,d_start_);
      tab = (unsigned long long int)FTK_n_svd_l*(unsigned long long int)n_w_max*(unsigned long long int)n_M_sub*(unsigned long long int)pm_n_UX_rank;
      na = (unsigned long long int)FTK_n_svd_l*(unsigned long long int)n_w_max*(unsigned long long int)pm_n_UX_rank*(unsigned long long int)(n_M_per_Mbatch*nMbatch);
      for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
	for (pm_nUX_rank=0;pm_nUX_rank<pm_n_UX_rank;pm_nUX_rank++){
	  for (nw=0;nw<n_w_max;nw++){
	    for (nl=0;nl<FTK_n_svd_l;nl++){
	      /* tabA = */
	      /* 	(unsigned long long int)pm_nUX_rank */
	      /* 	+ (unsigned long long int)nM_sub*(unsigned long long int)pm_n_UX_rank */
	      /* 	+ (unsigned long long int)nw*(unsigned long long int)n_M_sub*(unsigned long long int)pm_n_UX_rank */
	      /* 	+ (unsigned long long int)nl*(unsigned long long int)n_w_max*(unsigned long long int)n_M_sub*(unsigned long long int)pm_n_UX_rank; */
	      tabA = (unsigned long long int)pm_nUX_rank + ((unsigned long long int)nM_sub + ((unsigned long long int)nw + (unsigned long long int)nl*(unsigned long long int)n_w_max)*(unsigned long long int)n_M_sub)*(unsigned long long int)pm_n_UX_rank;
	      tabB = (unsigned long long int)pm_nUX_rank + ((unsigned long long int)nM_sub + ((unsigned long long int)nl + (unsigned long long int)nw*(unsigned long long int)FTK_n_svd_l)*(unsigned long long int)n_M_sub)*(unsigned long long int)pm_n_UX_rank;
	      svd_VUXM_nMwl____[tabA] = svd_VUXM_lwnM____[na];
	      svd_VUXM_nMlw____[tabB] = svd_VUXM_lwnM____[na];
	      na++;
	      /* for (nl=0;nl<FTK_n_svd_l;nl++){ } */}
	    /* for (nw=0;nw<n_w_max;nw++){ } */}
	  /* for (pm_nUX_rank=0;pm_nUX_rank<pm_n_UX_rank;pm_nUX_rank++){ } */}
	/* for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ } */}
      local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"svd_VUXM_nMwl____: ");
      elrt_sum_svd_VUXM_nMwl____ += elrt_[0];
      n_op_sum_svd_VUXM_nMwl____ += tab;
      /* %%%% */
      local_tic(0,t_start_,d_start_);
      tab = (unsigned long long int)FTK_n_svd_l*(unsigned long long int)n_w_max*(unsigned long long int)pm_n_UX_rank*(unsigned long long int)n_M_sub;
      tabA = (unsigned long long int)FTK_n_svd_l_rup*(unsigned long long int)n_w_max_rup*(unsigned long long int)pm_n_UX_rank_rup*(unsigned long long int)n_M_sub;
      memset(f_svd_VUXM_nMlw_real____,0,tabA*sizeof(float));
      memset(f_svd_VUXM_nMlw_imag____,0,tabA*sizeof(float));
      tabB = (unsigned long long int)FTK_n_svd_l_rup*(unsigned long long int)n_w_max_rup*(unsigned long long int)pm_n_UX_rank_rup*(unsigned long long int)(n_M_per_Mbatch*nMbatch);
      transpose_ps_block1(
      			  f_svd_VUXM_lwnM_real____ + tabB
      			  ,f_svd_VUXM_nMlw_real____
      			  ,FTK_n_svd_l_rup*n_w_max_rup
      			  ,pm_n_UX_rank_rup*n_M_sub_rup
      			  ,FTK_n_svd_l_rup*n_w_max_rup
      			  ,pm_n_UX_rank_rup*n_M_sub_rup
      			  ,16
      			  );
      transpose_ps_block1(
      			  f_svd_VUXM_lwnM_imag____ + tabB
      			  ,f_svd_VUXM_nMlw_imag____
      			  ,FTK_n_svd_l_rup*n_w_max_rup
      			  ,pm_n_UX_rank_rup*n_M_sub_rup
      			  ,FTK_n_svd_l_rup*n_w_max_rup
      			  ,pm_n_UX_rank_rup*n_M_sub_rup
      			  ,16
      			  );
      if (verbose>2){
	array_sub_printf(  svd_VUXM_nMlw____ ,"double complex",pm_n_UX_rank    *n_M_sub    ,4,FTK_n_svd_l    *n_w_max    ,3," %%   svd_VUXM_nMlw____ : ");
	array_sub_printf(f_svd_VUXM_nMlw_real____,"float" ,pm_n_UX_rank_rup*n_M_sub_rup,4,FTK_n_svd_l_rup*n_w_max_rup,3," %% f_svd_VUXM_nMlw_real____: ");
	array_sub_printf(f_svd_VUXM_nMlw_imag____,"float" ,pm_n_UX_rank_rup*n_M_sub_rup,4,FTK_n_svd_l_rup*n_w_max_rup,3," %% f_svd_VUXM_nMlw_imag____: ");
	exit(0);
	/* if (verbose>2){ } */}  
      local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"f_svd_VUXM_nMlw____: ");
      elrt_sum_f_svd_VUXM_nMlw____ += elrt_[0];
      n_op_sum_f_svd_VUXM_nMlw____ += tab;
      /* %%%% */
      /* 
	 svd_SVUXM_SMwl____ = zeros(n_S,n_M_sub,n_w_max,FTK.n_svd_l);
	 for nl=0:FTK.n_svd_l-1;
	 for nw=0:n_w_max-1;
	 svd_SVUXM_SMwl____(:,:,1+nw,1+nl) = ctranspose(CTF_UX_S_k_q_nSw___(:,:,1+nw))*svd_VUXM_nMwl____(:,:,1+nw,1+nl);
	 end;%for nw=0:n_w_max-1;
	 end;%for nl=0:FTK.n_svd_l-1;
      */
      /* %%%% */
      local_tic(0,t_start_,d_start_);
      tab = (unsigned long long int)(FTK_n_svd_l*n_w_max)*(unsigned long long int)(n_S*n_M_sub)*(unsigned long long int)pm_n_UX_rank;
      for (nl=0;nl<FTK_n_svd_l;nl++){
	for (nw=0;nw<n_w_max;nw++){
	  tabA = (unsigned long long int)nw*(unsigned long long int)pm_n_UX_rank*(unsigned long long int)n_S;
	  tabB = (unsigned long long int)(nw+nl*n_w_max)*(unsigned long long int)pm_n_UX_rank*(unsigned long long int)n_M_sub;
	  tabC = (unsigned long long int)(nw+nl*n_w_max)*(unsigned long long int)n_S*(unsigned long long int)n_M_sub;
	  cblas_zgemm(
		      CblasColMajor
		      ,CblasConjTrans
		      ,CblasNoTrans
		      ,n_S
		      ,n_M_sub
		      ,pm_n_UX_rank
		      ,&cblas_alpha
		      ,CTF_UX_S_k_q_nSw___ + tabA
		      ,pm_n_UX_rank
		      ,svd_VUXM_nMwl____ + tabB
		      ,pm_n_UX_rank
		      ,&cblas_beta
		      ,svd_SVUXM_SMwl____ + tabC
		      ,n_S
		      );
	  /* for (nw=0;nw<n_w_max;nw++){ } */}
	/* for (nl=0;nl<FTK_n_svd_l;nl++){ } */}
      local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"svd_SVUXM_SMwl____: ");
      elrt_sum_svd_SVUXM_SMwl____ += elrt_[0];
      n_op_sum_svd_SVUXM_SMwl____ += tab;
      /* %%%% */
      local_tic(0,t_start_,d_start_);
      tab = (unsigned long long int)(FTK_n_svd_l*n_w_max)*(unsigned long long int)(n_S*n_M_sub)*(unsigned long long int)pm_n_UX_rank;
      for (nl=0;nl<FTK_n_svd_l;nl++){
	for (nw=0;nw<n_w_max;nw++){
	  tabA = (unsigned long long int)nw*(unsigned long long int)pm_n_UX_rank*(unsigned long long int)n_S;
	  tabB = (unsigned long long int)(nl+nw*FTK_n_svd_l)*(unsigned long long int)pm_n_UX_rank*(unsigned long long int)n_M_sub;
	  tabC = (unsigned long long int)(nl+nw*FTK_n_svd_l)*(unsigned long long int)n_S*(unsigned long long int)n_M_sub;
	  cblas_zgemm(
		      CblasColMajor
		      ,CblasConjTrans
		      ,CblasNoTrans
		      ,n_S
		      ,n_M_sub
		      ,pm_n_UX_rank
		      ,&cblas_alpha
		      ,CTF_UX_S_k_q_nSw___ + tabA
		      ,pm_n_UX_rank
		      ,svd_VUXM_nMlw____ + tabB
		      ,pm_n_UX_rank
		      ,&cblas_beta
		      ,svd_SVUXM_SMlw____ + tabC
		      ,n_S
		      );
	  /* for (nw=0;nw<n_w_max;nw++){ } */}
	/* for (nl=0;nl<FTK_n_svd_l;nl++){ } */}
      local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"svd_SVUXM_SMlw____: ");
      elrt_sum_svd_SVUXM_SMlw____ += elrt_[0];
      n_op_sum_svd_SVUXM_SMlw____ += tab;
      /* %%%% */
      local_tic(0,t_start_,d_start_);
      tab = (unsigned long long int)(FTK_n_svd_l*n_w_max)*(unsigned long long int)(n_S*n_M_sub)*(unsigned long long int)pm_n_UX_rank;
      tabA = (unsigned long long int)(FTK_n_svd_l_rup*n_w_max_rup)*(unsigned long long int)(n_S_rup*n_M_sub_rup);
      memset(f_svd_SVUXM_SMlw_real____,0,tabA*sizeof(float));
      memset(f_svd_SVUXM_SMlw_imag____,0,tabA*sizeof(float));
      for (nl=0;nl<FTK_n_svd_l;nl++){
	for (nw=0;nw<n_w_max;nw++){
	  tabA = (unsigned long long int)nw*(unsigned long long int)pm_n_UX_rank_rup*(unsigned long long int)n_S_rup;
	  tabB = (unsigned long long int)(nl+nw*FTK_n_svd_l_rup)*(unsigned long long int)pm_n_UX_rank_rup*(unsigned long long int)n_M_sub_rup;
	  tabC = (unsigned long long int)(nl+nw*FTK_n_svd_l_rup)*(unsigned long long int)n_S_rup*(unsigned long long int)n_M_sub_rup;
	  f_CR_ = f_svd_SVUXM_SMlw_real____ + tabC;
	  f_CI_ = f_svd_SVUXM_SMlw_imag____ + tabC;
	  hp_segregated_to_segregated_mult_immintrin_load1_fma(
							       n_S_rup
							       ,pm_n_UX_rank_rup
							       ,f_CTF_UX_S_k_q_nSw_real___ + tabA
							       ,f_CTF_UX_S_k_q_nSw_imag___ + tabA
							       ,n_M_sub_rup
							       ,f_svd_VUXM_nMlw_real____ + tabB
							       ,f_svd_VUXM_nMlw_imag____ + tabB
							       ,&f_CR_
							       ,&f_CI_
							       );
	  /* for (nw=0;nw<n_w_max;nw++){ } */}
	/* for (nl=0;nl<FTK_n_svd_l;nl++){ } */}
      local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"svd_f_SVUXM_SMlw____: ");
      elrt_sum_f_svd_SVUXM_SMlw____ += elrt_[0];
      n_op_sum_f_svd_SVUXM_SMlw____ += tab;
      if (verbose>2){
	array_sub_printf(  svd_SVUXM_SMlw____ ,"double complex",n_S    *n_M_sub    ,4,FTK_n_svd_l    *n_w_max    ,3," %%   svd_SVUXM_SMlw____ : ");
	array_sub_printf(f_svd_SVUXM_SMlw_real____,"float" ,n_S_rup*n_M_sub_rup,4,FTK_n_svd_l_rup*n_w_max_rup,3," %% f_svd_SVUXM_SMlw_real____: ");
	array_sub_printf(f_svd_SVUXM_SMlw_imag____,"float" ,n_S_rup*n_M_sub_rup,4,FTK_n_svd_l_rup*n_w_max_rup,3," %% f_svd_SVUXM_SMlw_imag____: ");
	exit(0);
	/* if (verbose>2){ } */}
      /* %%%% */
      /* svd_SVUXM_lwSM____ = permute(ifft(permute(svd_SVUXM_SMwl____,[3,4,1,2]),[],1)*n_w_max,[2,1,3,4]); */
      /* %%%% */      
      local_tic(0,t_start_,d_start_);
      tab = (unsigned long long int)n_S*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_svd_l;
      memset(svd_SVUXM_0in_wlSM____,0,tab*sizeof(fftw_complex));
      na=0; tabA=0;
      for (nl=0;nl<FTK_n_svd_l;nl++){
	for (nw=0;nw<n_w_max;nw++){
	  for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
	    for (nS=0;nS<n_S;nS++){
	      /* tabA = */
	      /* 	(unsigned long long int)nw */
	      /* 	+ (unsigned long long int)nl*(unsigned long long int)n_w_max */
	      /* 	+ (unsigned long long int)nS*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_svd_l */
	      /* 	+ (unsigned long long int)nM_sub*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_svd_l*(unsigned long long int)n_S; */
	      tabA = (unsigned long long int)nw + ((unsigned long long int)nl + ((unsigned long long int)nS + (unsigned long long int)nM_sub*(unsigned long long int)n_S)*(unsigned long long int)FTK_n_svd_l)*(unsigned long long int)n_w_max;
	      svd_SVUXM_0in_wlSM____[tabA] = (fftw_complex) (svd_SVUXM_SMwl____[na]); na++;
	      /* for (nS=0;nS<n_S;nS++){ } */}
	    /* for (nM_sub=0;nM_sub<n_M_sub;nS++){ } */}
	  /* for (nw=0;nw<n_w_max;nw++){ } */}
	/* for (nl=0;nl<FTK_n_svd_l;nl++){ } */}
      local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"svd_SVUXM_0in_wlSM____: ");
      elrt_sum_svd_SVUXM_0in_wlSM____ += elrt_[0];
      n_op_sum_svd_SVUXM_0in_wlSM____ += tab;
      local_tic(0,t_start_,d_start_);
      fftw_execute(fftw_plan_many_plan);
      local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"svd_SVUXM_out_wlSM____: ");
      elrt_sum_svd_SVUXM_out_wlSM____ += elrt_[0];
      n_op_sum_svd_SVUXM_out_wlSM____ += tab;
      local_tic(0,t_start_,d_start_);
      tab = (unsigned long long int)n_S*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_svd_l;
      memset(svd_SVUXM_lwSM____,0,tab*sizeof(double complex));
      na=0; tabA=0;
      for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
	for (nS=0;nS<n_S;nS++){
	  for (nl=0;nl<FTK_n_svd_l;nl++){
	    for (nw=0;nw<n_w_max;nw++){
	      /* tabA = */
	      /* 	(unsigned long long int)nl */
	      /* 	+ (unsigned long long int)nw*(unsigned long long int)FTK_n_svd_l */
	      /* 	+ (unsigned long long int)nS*(unsigned long long int)FTK_n_svd_l*(unsigned long long int)n_w_max */
	      /* 	+ (unsigned long long int)nM_sub*(unsigned long long int)FTK_n_svd_l*(unsigned long long int)n_w_max*(unsigned long long int)n_S; */
	      tabA = (unsigned long long int)nl + ((unsigned long long int)nw + ((unsigned long long int)nS + (unsigned long long int)nM_sub*(unsigned long long int)n_S)*(unsigned long long int)n_w_max)*(unsigned long long int)FTK_n_svd_l;
	      svd_SVUXM_lwSM____[tabA] = (double complex) (svd_SVUXM_out_wlSM____[na]); na++;
	      /* for (nw=0;nw<n_w_max;nw++){ } */}
	    /* for (nl=0;nl<FTK_n_svd_l;nl++){ } */}
	  /* for (nS=0;nS<n_S;nS++){ } */}
	/* for (nM_sub=0;nM_sub<n_M_sub;nS++){ } */}
      local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"svd_SVUXM_lwSM____: ");
      elrt_sum_svd_SVUXM_lwSM____ += elrt_[0];
      n_op_sum_svd_SVUXM_lwSM____ += tab;
      /* %%%% */
      local_tic(0,t_start_,d_start_);
      tab = (unsigned long long int)n_w_max*(unsigned long long int)n_S*(unsigned long long int)n_M_sub*(unsigned long long int)FTK_n_svd_l;
      tabA = (unsigned long long int)n_w_max_rup*(unsigned long long int)n_S_rup*(unsigned long long int)n_M_sub_rup*(unsigned long long int)FTK_n_svd_l_rup;
      memset(f_svd_SVUXM_0in_wSMl____,0,tabA*sizeof(fftwf_complex));
      f_C_ = (float *) f_svd_SVUXM_0in_wSMl____;
      for (nl=0;nl<FTK_n_svd_l;nl++){
	for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
	  for (nS=0;nS<n_S;nS++){
	    for (nw=0;nw<n_w_max;nw++){
	      tabA = (unsigned long long int)nS + ((unsigned long long int)nM_sub + ((unsigned long long int)nl + (unsigned long long int)nw*(unsigned long long int)FTK_n_svd_l_rup)*(unsigned long long int)n_M_sub_rup)*(unsigned long long int)n_S_rup;
	      tabB = (unsigned long long int)nw + ((unsigned long long int)nS + ((unsigned long long int)nM_sub + (unsigned long long int)nl*(unsigned long long int)n_M_sub_rup)*(unsigned long long int)n_S_rup)*(unsigned long long int)n_w_max_rup;
	      f_C_[2*tabB + 0] = f_svd_SVUXM_SMlw_real____[tabA]; 
	      f_C_[2*tabB + 1] = f_svd_SVUXM_SMlw_imag____[tabA]; 
	      na++;
	      /* for (nw=0;nw<n_w_max;nw++){ } */}
	    /* for (nS=0;nS<n_S;nS++){ } */}
	  /* for (nM_sub=0;nM_sub<n_M_sub;nS++){ } */}
	/* for (nl=0;nl<FTK_n_svd_l;nl++){ } */}
      local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"f_svd_SVUXM_0in_wSMl____: ");
      elrt_sum_f_svd_SVUXM_0in_wSMl____ += elrt_[0];
      n_op_sum_f_svd_SVUXM_0in_wSMl____ += tab;
      local_tic(0,t_start_,d_start_);
      fftwf_execute(fftwf_plan_many_plan);
      local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"f_svd_SVUXM_out_wSMl____: ");
      elrt_sum_f_svd_SVUXM_out_wSMl____ += elrt_[0];
      n_op_sum_f_svd_SVUXM_out_wSMl____ += tab;
      local_tic(0,t_start_,d_start_);
      tab = (unsigned long long int)FTK_n_svd_l*(unsigned long long int)n_w_max*(unsigned long long int)n_S*(unsigned long long int)n_M_sub;
      tabA = (unsigned long long int)FTK_n_svd_l_rup*(unsigned long long int)n_w_max_rup*(unsigned long long int)n_S_rup*(unsigned long long int)n_M_sub_rup;
      memset(f_svd_SVUXM_lwSM_real____,0,tabA*sizeof(float));
      memset(f_svd_SVUXM_lwSM_imag____,0,tabA*sizeof(float));
      f_C_ = (float *) f_svd_SVUXM_out_wSMl____;
      for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
	for (nS=0;nS<n_S;nS++){
	  for (nw=0;nw<n_w_max;nw++){
	    for (nl=0;nl<FTK_n_svd_l;nl++){
      	      tabA = (unsigned long long int)nw + ((unsigned long long int)nS + ((unsigned long long int)nM_sub + (unsigned long long int)nl*(unsigned long long int)n_M_sub_rup)*(unsigned long long int)n_S_rup)*(unsigned long long int)n_w_max_rup;
      	      tabB = (unsigned long long int)nl + ((unsigned long long int)nw + ((unsigned long long int)nS + (unsigned long long int)nM_sub*(unsigned long long int)n_S_rup)*(unsigned long long int)n_w_max_rup)*(unsigned long long int)FTK_n_svd_l_rup;
      	      f_svd_SVUXM_lwSM_real____[tabB] = f_C_[2*tabA+0];
      	      f_svd_SVUXM_lwSM_imag____[tabB] = f_C_[2*tabA+1];
	      /* for (nl=0;nl<FTK_n_svd_l;nl++){ } */}
	    /* for (nw=0;nw<n_w_max;nw++){ } */}
	  /* for (nS=0;nS<n_S;nS++){ } */}
	/* for (nM_sub=0;nM_sub<n_M_sub;nS++){ } */}
      local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"f_svd_SVUXM_lwSM____: ");
      elrt_sum_f_svd_SVUXM_lwSM____ += elrt_[0];
      n_op_sum_f_svd_SVUXM_lwSM____ += tab;
      if (verbose>2){
	array_sub_printf(  svd_SVUXM_lwSM____ ,"double complex",FTK_n_svd_l    *n_w_max    ,3,n_S    *n_M_sub    ,4," %%   svd_SVUXM_lwSM____ : ");
	array_sub_printf(f_svd_SVUXM_lwSM_real____,"float" ,FTK_n_svd_l_rup*n_w_max_rup,3,n_S_rup*n_M_sub_rup,4," %% f_svd_SVUXM_lwSM_real____: ");
	array_sub_printf(f_svd_SVUXM_lwSM_imag____,"float" ,FTK_n_svd_l_rup*n_w_max_rup,3,n_S_rup*n_M_sub_rup,4," %% f_svd_SVUXM_lwSM_imag____: ");
	exit(0);
	/* if (verbose>2){ } */}      
      /* %%%% */
      /* display output */
      /* %%%% */
      if (flag_disp){
	printf(" %% nMbatch %d/%d\n",nMbatch,n_Mbatch);
	tab = minimum(8,(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_svd_l*(unsigned long long int)n_S*(unsigned long long int)n_M_sub);
	array_printf(svd_SVUXM_SMwl____,"double complex",1,tab," %% svd_SVUXM_SMwl____: ");
	array_printf(svd_SVUXM_lwSM____,"double complex",1,tab," %% svd_SVUXM_lwSM____: ");
	/* if (flag_disp){ } */}
      if (flag_dump){
	MDA_ndim=0;
	MDA_dim_[MDA_ndim++] = n_S;
	MDA_dim_[MDA_ndim++] = n_M_sub;
	MDA_dim_[MDA_ndim++] = n_w_max;
	MDA_dim_[MDA_ndim++] = FTK_n_svd_l;
	MDA_n_dim = MDA_ndim;
	sprintf(MDA_fname,"mex_svd_SVUXM_SMwl____.mat");
	MDA_write_c16(MDA_n_dim,MDA_dim_,svd_SVUXM_SMwl____,MDA_fname);
	MDA_ndim=0;
	MDA_dim_[MDA_ndim++] = FTK_n_svd_l;
	MDA_dim_[MDA_ndim++] = n_w_max;
	MDA_dim_[MDA_ndim++] = n_S;
	MDA_dim_[MDA_ndim++] = n_M_sub;
	MDA_n_dim = MDA_ndim;
	sprintf(MDA_fname,"mex_svd_SVUXM_lwSM____.mat");
	MDA_write_c16(MDA_n_dim,MDA_dim_,svd_SVUXM_lwSM____,MDA_fname);
	/* if (flag_dump){ } */}
      for (nSbatch=0;nSbatch<n_Sbatch;nSbatch++){
	n_S_sub = n_S_per_Sbatch; if (nSbatch==n_Sbatch-1){ n_S_sub = n_S - n_S_per_Sbatch*nSbatch;}
	for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){ index_S_in_Sbatch_[nS_sub] = n_S_per_Sbatch*nSbatch + nS_sub;}
	if (verbose>1){ printf(" %% nSbatch %d/%d: index_S_in_Sbatch_ %d --> %d\n",nSbatch,n_Sbatch,index_S_in_Sbatch_[0],index_S_in_Sbatch_[n_S_sub-1]);}
	if (n_S_sub>0){
	  n_S_sub_rup  = rup(n_S_sub,8);
	  n_S_sub_256 = n_S_sub_rup/8;
	  /* %%%% */
	  /* svd_SVUXM_lwsM____ = svd_SVUXM_lwSM____(:,:,1+index_S_in_Sbatch_,:); */
	  /* %%%% */	  
	  local_tic(0,t_start_,d_start_);
	  tab = (unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_svd_l;
	  memset(svd_SVUXM_lwsM____,0,tab*sizeof(double complex));
	  for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
	    for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){
	      nS = index_S_in_Sbatch_[nS_sub];
	      tabA = (((unsigned long long int)nS_sub + (unsigned long long int)nM_sub*(unsigned long long int)n_S_sub)*(unsigned long long int)n_w_max)*(unsigned long long int)FTK_n_svd_l;
	      tabB = (((unsigned long long int)nS     + (unsigned long long int)nM_sub*(unsigned long long int)n_S    )*(unsigned long long int)n_w_max)*(unsigned long long int)FTK_n_svd_l;
	      tabC = (unsigned long long int)n_w_max*(unsigned long long int)FTK_n_svd_l;
	      memcpy( svd_SVUXM_lwsM____ + tabA , svd_SVUXM_lwSM____ + tabB , tabC*sizeof(double complex));
	      /* for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){ } */}
	    /* for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ } */}
	  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"svd_SVUXM_lwsM____: ");
	  elrt_sum_svd_SVUXM_lwsM____ += elrt_[0];
	  n_op_sum_svd_SVUXM_lwsM____ += tab;
	  /* %%%% */
	  local_tic(0,t_start_,d_start_);
	  tab = (unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_svd_l;
	  tabA = (unsigned long long int)n_S_sub_rup*(unsigned long long int)n_M_sub_rup*(unsigned long long int)n_w_max_rup*(unsigned long long int)FTK_n_svd_l_rup;
	  memset(f_svd_SVUXM_lwsM_real____,0,tabA*sizeof(float));
	  memset(f_svd_SVUXM_lwsM_imag____,0,tabA*sizeof(float));
	  for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
	    for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){
	      nS = index_S_in_Sbatch_[nS_sub];
	      tabA = (((unsigned long long int)nS_sub + (unsigned long long int)nM_sub*(unsigned long long int)n_S_sub_rup)*(unsigned long long int)n_w_max_rup)*(unsigned long long int)FTK_n_svd_l_rup;
	      tabB = (((unsigned long long int)nS     + (unsigned long long int)nM_sub*(unsigned long long int)n_S_rup    )*(unsigned long long int)n_w_max_rup)*(unsigned long long int)FTK_n_svd_l_rup;
	      tabC = (unsigned long long int)n_w_max_rup*(unsigned long long int)FTK_n_svd_l_rup;
	      memcpy( f_svd_SVUXM_lwsM_real____ + tabA , f_svd_SVUXM_lwSM_real____ + tabB , tabC*sizeof(float));
	      memcpy( f_svd_SVUXM_lwsM_imag____ + tabA , f_svd_SVUXM_lwSM_imag____ + tabB , tabC*sizeof(float));
	      /* for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){ } */}
	    /* for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ } */}
	  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"f_svd_SVUXM_lwsM____: ");
	  elrt_sum_f_svd_SVUXM_lwsM____ += elrt_[0];
	  n_op_sum_f_svd_SVUXM_lwsM____ += tab;
	  if (verbose>2){
	    array_sub_printf(  svd_SVUXM_lwsM____ ,"double complex",FTK_n_svd_l    *n_w_max    ,3,n_S_sub    *n_M_sub    ,4," %%   svd_SVUXM_lwsM____ : ");
	    array_sub_printf(f_svd_SVUXM_lwsM_real____,"float" ,FTK_n_svd_l_rup*n_w_max_rup,3,n_S_sub_rup*n_M_sub_rup,4," %% f_svd_SVUXM_lwsM_real____: ");
	    array_sub_printf(f_svd_SVUXM_lwsM_imag____,"float" ,FTK_n_svd_l_rup*n_w_max_rup,3,n_S_sub_rup*n_M_sub_rup,4," %% f_svd_SVUXM_lwsM_imag____: ");
	    exit(0);
	    /* if (verbose>2){ } */}      
	  /* %%%% */	  
	  if (flag_dump){
	    MDA_ndim=0;
	    MDA_dim_[MDA_ndim++] = FTK_n_svd_l;
	    MDA_dim_[MDA_ndim++] = n_w_max;
	    MDA_dim_[MDA_ndim++] = n_S_sub;
	    MDA_dim_[MDA_ndim++] = n_M_sub;
	    MDA_n_dim = MDA_ndim;
	    sprintf(MDA_fname,"mex_svd_SVUXM_lwsM____.mat");
	    MDA_write_c16(MDA_n_dim,MDA_dim_,svd_SVUXM_lwsM____,MDA_fname);
	    /* if (flag_dump){ } */}
	  /* %%%% */
	  /* svd_USESVUXM_dwSM____ = real(reshape(FTK.svd_U_d_expiw_s__*reshape(svd_SVUXM_lwsM____,[FTK.n_svd_l,n_w_max*n_S_sub*n_M_sub]),[FTK.n_delta_v,n_w_max,n_S_sub,n_M_sub])); */
	  /* %%%% */
	  local_tic(0,t_start_,d_start_);
	  tab = (unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_delta_v;
	  memset(svd_USESVUXM_dwSM____,0,tab*sizeof(double complex));
	  tabC = (unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max;
	  tab = (unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_delta_v*(unsigned long long int)FTK_n_svd_l;
	  cblas_zgemm(
		      CblasColMajor
		      ,CblasNoTrans
		      ,CblasNoTrans
		      ,FTK_n_delta_v
		      ,tabC
		      ,FTK_n_svd_l
		      ,&cblas_alpha
		      ,FTK_svd_U_d_expiw_s__
		      ,FTK_n_delta_v
		      ,svd_SVUXM_lwsM____
		      ,FTK_n_svd_l
		      ,&cblas_beta
		      ,svd_USESVUXM_dwSM____
		      ,FTK_n_delta_v
		      );
	  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"svd_USESVUXM_dwSM____: ");
	  elrt_sum_svd_USESVUXM_dwSM____ += elrt_[0];
	  n_op_sum_svd_USESVUXM_dwSM____ += tab;
	  /* %%%% */
	  local_tic(0,t_start_,d_start_);
	  tab = (unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_delta_v*(unsigned long long int)FTK_n_svd_l;
	  tabA = (unsigned long long int)n_S_sub_rup*(unsigned long long int)n_M_sub_rup*(unsigned long long int)n_w_max_rup*(unsigned long long int)FTK_n_delta_v_rup;
	  memset(f_svd_USESVUXM_dwSM_real____,0,tabA*sizeof(float));
	  memset(f_svd_USESVUXM_dwSM_imag____,0,tabA*sizeof(float));
	  tabC = (unsigned long long int)n_S_sub_rup*(unsigned long long int)n_M_sub_rup*(unsigned long long int)n_w_max_rup;
	  nhpr_segregated_to_segregated_mult_immintrin_load1_fma(
							       FTK_n_delta_v_rup
							       ,FTK_n_svd_l_rup
							       ,f_FTK_svd_U_d_expiw_s_tran_real__
							       ,f_FTK_svd_U_d_expiw_s_tran_imag__
							       ,tabC
							       ,f_svd_SVUXM_lwsM_real____
							       ,f_svd_SVUXM_lwsM_imag____
							       ,&f_svd_USESVUXM_dwSM_real____
							       ,&f_svd_USESVUXM_dwSM_imag____
							       );
	  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"f_svd_USESVUXM_dwSM____: ");
	  elrt_sum_f_svd_USESVUXM_dwSM____ += elrt_[0];
	  n_op_sum_f_svd_USESVUXM_dwSM____ += tab;
	  if (verbose>2){
	    array_sub_printf(  svd_USESVUXM_dwSM____ ,"double complex",FTK_n_delta_v    *n_w_max    ,3,n_S_sub    *n_M_sub    ,4," %%   svd_USESVUXM_dwSM____ : ");
	    array_sub_printf(f_svd_USESVUXM_dwSM_real____,"float" ,FTK_n_delta_v_rup*n_w_max_rup,3,n_S_sub_rup*n_M_sub_rup,4," %% f_svd_USESVUXM_dwSM_real____: ");
	    array_sub_printf(f_svd_USESVUXM_dwSM_imag____,"float" ,FTK_n_delta_v_rup*n_w_max_rup,3,n_S_sub_rup*n_M_sub_rup,4," %% f_svd_USESVUXM_dwSM_imag____: ");
	    exit(0);
	    /* if (verbose>2){ } */}      
	  /* %%%% */
	  if (flag_dump){
	    MDA_ndim=0;
	    MDA_dim_[MDA_ndim++] = FTK_n_delta_v;
	    MDA_dim_[MDA_ndim++] = n_w_max;
	    MDA_dim_[MDA_ndim++] = n_S_sub;
	    MDA_dim_[MDA_ndim++] = n_M_sub;
	    MDA_n_dim = MDA_ndim;
	    sprintf(MDA_fname,"mex_svd_USESVUXM_dwSM____.mat");
	    MDA_write_c16(MDA_n_dim,MDA_dim_,svd_USESVUXM_dwSM____,MDA_fname);
	    /* if (flag_dump){ } */}
	  /* %%%% */
	  /* 
	     l2_dSM___ = permute(reshape(reshape(sqrt(CTF_UX_S_l2_use_(1+index_S_in_Sbatch_)),[n_S_sub,1])*reshape(sqrt(UX_M_l2_dM_use__(:,1+index_M_in_Mbatch_)),[1,FTK.n_delta_v*n_M_sub]),[n_S_sub,FTK.n_delta_v,n_M_sub]),[2,1,3]);
	     n2_dSM___ = 1./max(1e-14,l2_dSM___);
	     f2_dSM___ = permute(reshape(reshape(sqrt(CTF_UX_S_l2_use_(1+index_S_in_Sbatch_)),[n_S_sub,1])*reshape(1./max(1e-14,sqrt(UX_M_l2_dM_use__(:,1+index_M_in_Mbatch_))),[1,FTK.n_delta_v*n_M_sub]),[n_S_sub,FTK.n_delta_v,n_M_sub]),[2,1,3]);
	     ss_S_ = reshape(CTF_UX_S_l2_use_(1+index_S_in_Sbatch_),[n_S_sub,1]);
	  */
	  /* %%%% */
	  local_tic(0,t_start_,d_start_);
	  tab = (unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub*(unsigned long long int)FTK_n_delta_v;
	  memset(l2_dSM___,0,tab*sizeof(double));
	  memset(n2_dSM___,0,tab*sizeof(double));
	  memset(f2_dSM___,0,tab*sizeof(double));
	  for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){
	    nS = index_S_in_Sbatch_[nS_sub];
	    for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
	      nM = index_M_in_Mbatch_[nM_sub];
	      for (ndelta_v=0;ndelta_v<FTK_n_delta_v;ndelta_v++){
		tabA = (unsigned long long int)ndelta_v + (unsigned long long int)nM*(unsigned long long int)FTK_n_delta_v;
		tabB = (unsigned long long int)ndelta_v + ((unsigned long long int)nS_sub + (unsigned long long int)nM_sub*(unsigned long long int)n_S_sub)*(unsigned long long int)FTK_n_delta_v;
		l2_dSM___[tabB] = sqrt(CTF_UX_S_l2_[nS]) * sqrt(UX_M_l2_dM__[tabA]);
		n2_dSM___[tabB] = 1.0 / maximum(1.0e-14,l2_dSM___[tabB]);
		f2_dSM___[tabB] = sqrt(CTF_UX_S_l2_[nS]) / maximum(1.0e-14,sqrt(UX_M_l2_dM__[tabA]));
		/* for (ndelta_v=0;ndelta_v<FTK_n_delta_v;ndelta_v++){ } */}
	      /* for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ } */}
	    /* for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){ } */}
	  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"l2_dSM___: ");
	  elrt_sum_l2_dSM___ += elrt_[0];
	  n_op_sum_l2_dSM___ += tab;
	  /* %%%% */
	  if (flag_dump){
	    MDA_ndim=0;
	    MDA_dim_[MDA_ndim++] = FTK_n_delta_v;
	    MDA_dim_[MDA_ndim++] = n_S_sub;
	    MDA_dim_[MDA_ndim++] = n_M_sub;
	    MDA_n_dim = MDA_ndim;
	    sprintf(MDA_fname,"mex_l2_dSM___.mat");
	    MDA_write_r8(MDA_n_dim,MDA_dim_,l2_dSM___,MDA_fname);
	    sprintf(MDA_fname,"mex_n2_dSM___.mat");
	    MDA_write_r8(MDA_n_dim,MDA_dim_,n2_dSM___,MDA_fname);
	    sprintf(MDA_fname,"mex_f2_dSM___.mat");
	    MDA_write_r8(MDA_n_dim,MDA_dim_,f2_dSM___,MDA_fname);
	    /* if (flag_dump){ } */}
	  /* %%%% */
	  /* X_sub_dwSM____ = repmat(reshape(n2_dSM___,[FTK.n_delta_v,1,n_S_sub,n_M_sub]),[1,n_w_max,1,1]) .* svd_USESVUXM_dwSM____; %<-- correlation. ; */
	  /* %%%% */
	  local_tic(0,t_start_,d_start_);
	  tab = (unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_delta_v;
	  memset(X_sub_dwSM____,0,tab*sizeof(double));
	  na=0;
	  for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
	    for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){
	      for (nw=0;nw<n_w_max;nw++){
		for (ndelta_v=0;ndelta_v<FTK_n_delta_v;ndelta_v++){
		  tabA = (unsigned long long int)ndelta_v + ((unsigned long long int)nS_sub + (unsigned long long int)nM_sub*(unsigned long long int)n_S_sub)*FTK_n_delta_v;
		  X_sub_dwSM____[na] = n2_dSM___[tabA] * creal(svd_USESVUXM_dwSM____[na]) ;
		  na += 1;
		  /* for (ndelta_v=0;ndelta_v<FTK_n_delta_v;ndelta_v++){ } */}
		/* for (nw=0;nw<n_w_max;nw++){ } */}
	      /* for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){ } */}
	    /* for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ } */}
	  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"X_sub_dwSM____: ");
	  elrt_sum_X_sub_dwSM____ += elrt_[0];
	  n_op_sum_X_sub_dwSM____ += tab;
	  /* %%%% */
	  local_tic(0,t_start_,d_start_);
	  tab = (unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_delta_v;
	  memset(d_X_sub_dwSM____,0,tab*sizeof(double));
	  na=0;
	  for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
	    for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){
	      for (nw=0;nw<n_w_max;nw++){
		for (ndelta_v=0;ndelta_v<FTK_n_delta_v;ndelta_v++){
		  tabA = (unsigned long long int)ndelta_v + ((unsigned long long int)nS_sub + (unsigned long long int)nM_sub*(unsigned long long int)n_S_sub)*FTK_n_delta_v;
		  tabB = (unsigned long long int)ndelta_v + ((unsigned long long int)nw + ((unsigned long long int)nS_sub + (unsigned long long int)nM_sub*(unsigned long long int)n_S_sub_rup)*(unsigned long long int)n_w_max_rup)*(unsigned long long int)FTK_n_delta_v_rup;
		  d_X_sub_dwSM____[na] = n2_dSM___[tabA] * f_svd_USESVUXM_dwSM_real____[tabB] ;
		  na += 1;
		  /* for (ndelta_v=0;ndelta_v<FTK_n_delta_v;ndelta_v++){ } */}
		/* for (nw=0;nw<n_w_max;nw++){ } */}
	      /* for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){ } */}
	    /* for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ } */}
	  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"d_X_sub_dwSM____: ");
	  elrt_sum_d_X_sub_dwSM____ += elrt_[0];
	  n_op_sum_d_X_sub_dwSM____ += tab;
	  if (verbose>2){
	    array_sub_printf(  X_sub_dwSM____ ,"double complex",FTK_n_delta_v*n_w_max,3,n_S_sub*n_M_sub,4," %%   X_sub_dwSM____ : ");
	    array_sub_printf(d_X_sub_dwSM____ ,"double complex",FTK_n_delta_v*n_w_max,3,n_S_sub*n_M_sub,4," %% d_X_sub_dwSM____ : ");
	    exit(0);
	    /* if (verbose>2){ } */}
	  /* %%%% */
	  tab = (unsigned long long int)FTK_n_delta_v*(unsigned long long int)n_w_max*(unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub;
	  memcpy( X_sub_dwSM____ , d_X_sub_dwSM____ , tab*sizeof(float));
	  /* %%%% */
	  if (flag_dump){
	    MDA_ndim=0;
	    MDA_dim_[MDA_ndim++] = FTK_n_delta_v;
	    MDA_dim_[MDA_ndim++] = n_w_max;
	    MDA_dim_[MDA_ndim++] = n_S_sub;
	    MDA_dim_[MDA_ndim++] = n_M_sub;
	    MDA_n_dim = MDA_ndim;
	    sprintf(MDA_fname,"mex_X_sub_dwSM____.mat");
	    MDA_write_r8(MDA_n_dim,MDA_dim_,X_sub_dwSM____,MDA_fname);
	    /* if (flag_dump){ } */}
	  /* %%%% */
	  /* 
	     I_value_sub_dwSM____ = repmat(reshape(f2_dSM___,[FTK.n_delta_v,1,n_S_sub,n_M_sub]),[1,n_w_max,1,1]) .* X_sub_dwSM____; %<-- I_value. ;
	     I_value_use_dwSM____ = max(0,real(I_value_sub_dwSM____));
	  */
	  /* %%%% */
	  if (flag_compute_I_value){
	    local_tic(0,t_start_,d_start_);
	    tab = (unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_delta_v;
	    memset(I_value_use_dwSM____,0,tab*sizeof(double));
	    na=0;
	    for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
	      for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){
		for (nw=0;nw<n_w_max;nw++){
		  for (ndelta_v=0;ndelta_v<FTK_n_delta_v;ndelta_v++){
		    tabA = (unsigned long long int)ndelta_v + ((unsigned long long int)nS_sub + (unsigned long long int)nM_sub*(unsigned long long int)n_S_sub)*FTK_n_delta_v;
		    I_value_use_dwSM____[na] = maximum( 0.0 , f2_dSM___[tabA] * X_sub_dwSM____[na] ) ;
		    na += 1;
		    /* for (ndelta_v=0;ndelta_v<FTK_n_delta_v;ndelta_v++){ } */}
		  /* for (nw=0;nw<n_w_max;nw++){ } */}
	      /* for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){ } */}
	    /* for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ } */}
	    local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"I_value_use_dwSM____: ");
	    elrt_sum_I_value_use_dwSM____ += elrt_[0];
	    n_op_sum_I_value_use_dwSM____ += tab;
	    /* %%%% */
	    if (flag_dump){
	      MDA_ndim=0;
	      MDA_dim_[MDA_ndim++] = FTK_n_delta_v;
	      MDA_dim_[MDA_ndim++] = n_w_max;
	      MDA_dim_[MDA_ndim++] = n_S_sub;
	      MDA_dim_[MDA_ndim++] = n_M_sub;
	      MDA_n_dim = MDA_ndim;
	      sprintf(MDA_fname,"mex_I_value_use_dwSM____.mat");
	      MDA_write_r8(MDA_n_dim,MDA_dim_,I_value_use_dwSM____,MDA_fname);
	      /* if (flag_dump){ } */}
	    /* if (flag_compute_I_value){ } */}
	  /* %%%% */
	  /*
	    if (flag_optimize_over_gamma_z == 0);
	    tmp_t = tic(); nop=0;
	    [tmp_X_wSM___,tmp_delta_ij___] = max(reshape(X_sub_dwSM____,[FTK.n_delta_v,n_w_max*n_S_sub*n_M_sub]),[],1); %<-- maximize correlation. ;
	    assert(min(tmp_delta_ij___)>=1); assert(max(tmp_delta_ij___)<=FTK.n_delta_v);
	    tmp_X_wSM___ = reshape(tmp_X_wSM___,[n_w_max,n_S_sub,n_M_sub]);
	    tmp_delta_ij___ = reshape(tmp_delta_ij___,[n_w_max,n_S_sub,n_M_sub]);
	    tmp_delta_x___ = FTK.delta_x_(tmp_delta_ij___);
	    tmp_delta_y___ = FTK.delta_y_(tmp_delta_ij___);
	    tmp_gamma_z___ = 2*pi*(0:n_w_max-1)/n_w_max;
	    X_wSM___(:,1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = tmp_X_wSM___;
	    delta_x_wSM___(:,1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = reshape(tmp_delta_x___,[n_w_max,n_S_sub,n_M_sub]);
	    delta_y_wSM___(:,1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = reshape(tmp_delta_y___,[n_w_max,n_S_sub,n_M_sub]);
	    gamma_z_wSM___(:,1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = repmat(tmp_gamma_z___(:),[1,n_S_sub,n_M_sub]);
	    if (flag_compute_I_value);
	    tmp_I_value_use_dwSM__ = reshape(I_value_use_dwSM____,[FTK.n_delta_v,n_w_max*n_S_sub*n_M_sub]);
	    tmp_I_value_use_wSM_ = zeros(n_w_max*n_S_sub*n_M_sub,1);
	    tmp_t2=tic();
	    for nl=0:n_w_max*n_S_sub*n_M_sub-1;
	    tmp_I_value_use_wSM_(1+nl) = tmp_I_value_use_dwSM__(tmp_delta_ij___(1+nl),1+nl);
	    end;%for nl=0:n_w_max*n_S_sub*n_M_sub-1;
	    I_value_wSM___(:,1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = reshape(tmp_I_value_use_wSM_,[n_w_max,n_S_sub,n_M_sub]);
	    tmp_t2 = toc(tmp_t2); if (verbose>1); disp(sprintf(' %% I_value_wSM___ %0.6fs',tmp_t2)); end;
	    end;%if (flag_compute_I_value);
	    tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_wSM___: %0.6f',tmp_t)); end;
	    parameter = parameter_timing_update(parameter,'ampmh_X_wSM_mex___8: X_wSM___',tmp_t,1,nop);
	    end;%if (flag_optimize_over_gamma_z == 0);
	  */
	  /* %%%% */
	  if (flag_optimize_over_gamma_z==0){
	    local_tic(0,t_start_,d_start_);
	    tab = (unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_delta_v;
	    na=0;
	    for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
	      nM = index_M_in_Mbatch_[nM_sub];
	      for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){
		nS = index_S_in_Sbatch_[nS_sub];
		for (nw=0;nw<n_w_max;nw++){
		  tabA = (unsigned long long int)nw + ((unsigned long long int)nS_sub + (unsigned long long int)nM_sub*(unsigned long long int)n_S_sub)*(unsigned long long int)n_w_max;
		  ndelta_v = dmax_index(FTK_n_delta_v,X_sub_dwSM____ + tabA*(unsigned long long int)FTK_n_delta_v);
		  tabB = (unsigned long long int)nw + ((unsigned long long int)nS     + (unsigned long long int)nM    *(unsigned long long int)n_S    )*(unsigned long long int)n_w_max;
		  tabC = ndelta_v + tabA*(unsigned long long int)FTK_n_delta_v;
		  X_wSM___[tabB] = X_sub_dwSM____[tabC];
		  delta_x_wSM___[tabB] = FTK_delta_x_[ndelta_v];
		  delta_y_wSM___[tabB] = FTK_delta_y_[ndelta_v];
		  gamma_z_wSM___[tabB] = gamma_z_[nw];
		  if (flag_compute_I_value){ I_value_wSM___[tabB] = I_value_use_dwSM____[tabC];}
		  na += FTK_n_delta_v;
		  /* for (nw=0;nw<n_w_max;nw++){ } */}
		/* for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){ } */}
	      /* for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ } */}
	    local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"X_wSM___: ");
	    elrt_sum_X_wSM___ += elrt_[0];
	    n_op_sum_X_wSM___ += tab;
	    /* if (flag_optimize_over_gamma_z==0){ } */}
	  /* %%%% */
	  /*
	    if (flag_optimize_over_gamma_z == 1);
	    tmp_t = tic(); nop=0;
	    [tmp_X_SM__,tmp_dw_ij__] = max(reshape(X_sub_dwSM____,[FTK.n_delta_v*n_w_max,n_S_sub*n_M_sub]),[],1); %<-- maximize correlation. ;
	    [tmp_delta_ij__,tmp_gamma_ij__] = ind2sub([FTK.n_delta_v,n_w_max],tmp_dw_ij__);
	    assert(min(tmp_delta_ij__)>=1); assert(max(tmp_delta_ij__)<=FTK.n_delta_v);
	    assert(min(tmp_gamma_ij__)>=1); assert(max(tmp_gamma_ij__)<=n_w_max);
	    tmp_X_SM__ = reshape(tmp_X_SM__,[n_S_sub,n_M_sub]);
	    tmp_delta_ij__ = reshape(tmp_delta_ij__,[n_S_sub,n_M_sub]);
	    tmp_gamma_ij__ = reshape(tmp_gamma_ij__,[n_S_sub,n_M_sub]);
	    tmp_delta_x__ = FTK.delta_x_(tmp_delta_ij__);
	    tmp_delta_y__ = FTK.delta_y_(tmp_delta_ij__);
	    tmp_gamma_z__ = 2*pi*(tmp_gamma_ij__-1)/n_w_max;
	    X_SM__(1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = tmp_X_SM__;
	    delta_x_SM__(1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = reshape(tmp_delta_x__,[n_S_sub,n_M_sub]);
	    delta_y_SM__(1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = reshape(tmp_delta_y__,[n_S_sub,n_M_sub]);
	    gamma_z_SM__(1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = reshape(tmp_gamma_z__,[n_S_sub,n_M_sub]);
	    if (flag_compute_I_value);
	    tmp_I_value_use_dwSM___ = reshape(I_value_use_dwSM____,[FTK.n_delta_v,n_w_max,n_S_sub*n_M_sub]);
	    tmp_I_value_use_SM_ = zeros(n_S_sub*n_M_sub,1);
	    tmp_t2=tic();
	    for nl=0:n_S_sub*n_M_sub-1;
	    tmp_I_value_use_SM_(1+nl) = tmp_I_value_use_dwSM___(tmp_delta_ij__(1+nl),tmp_gamma_ij__(1+nl),1+nl);
	    end;%for nl=0:n_S_sub*n_M_sub-1;
	    I_value_SM__(1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = reshape(tmp_I_value_use_SM_,[n_S_sub,n_M_sub]);
	    tmp_t2 = toc(tmp_t2); if (verbose>1); disp(sprintf(' %% I_value_SM__ %0.6fs',tmp_t2)); end;
	    end;%if (flag_compute_I_value);
	    tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_wSM___: %0.6f',tmp_t)); end;
	    parameter = parameter_timing_update(parameter,'ampmh_X_wSM_mex___8: X_wSM___',tmp_t,1,nop);
	    end;%if (flag_optimize_over_gamma_z == 1);	    
	  */
	  /* %%%% */
	  if (flag_optimize_over_gamma_z==1){
	    local_tic(0,t_start_,d_start_);
	    tab = (unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_delta_v;
	    na=0;
	    for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
	      nM = index_M_in_Mbatch_[nM_sub];
	      for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){
		nS = index_S_in_Sbatch_[nS_sub];
		tabA = (unsigned long long int)nS_sub + (unsigned long long int)nM_sub*(unsigned long long int)n_S_sub;
		ndw = dmax_index((unsigned long long int)FTK_n_delta_v*(unsigned long long int)n_w_max,X_sub_dwSM____ + tabA*(unsigned long long int)FTK_n_delta_v*(unsigned long long int)n_w_max);
		nw = ndw / maximum(1,FTK_n_delta_v);
		ndelta_v = ndw % FTK_n_delta_v;
		tabB = (unsigned long long int)nS     + (unsigned long long int)nM    *(unsigned long long int)n_S    ;
		tabC = ndw + tabA*(unsigned long long int)FTK_n_delta_v*(unsigned long long int)n_w_max;
		X_wSM___[tabB] = X_sub_dwSM____[tabC];
		delta_x_wSM___[tabB] = FTK_delta_x_[ndelta_v];
		delta_y_wSM___[tabB] = FTK_delta_y_[ndelta_v];
		gamma_z_wSM___[tabB] = gamma_z_[nw];
		if (flag_compute_I_value){ I_value_wSM___[tabB] = I_value_use_dwSM____[tabC];}
		na += (unsigned long long int)FTK_n_delta_v*(unsigned long long int)n_w_max;
		/* for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){ } */}
	      /* for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ } */}
	    local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"X_wSM___: ");
	    elrt_sum_X_wSM___ += elrt_[0];
	    n_op_sum_X_wSM___ += tab;
	    /* if (flag_optimize_over_gamma_z==1){ } */}
	  /* %%%% */
	  if (flag_disp){
	    printf(" %% nMbatch %d/%d nSbatch %d/%d\n",nMbatch,n_Mbatch,nSbatch,n_Sbatch);
	    if (flag_optimize_over_gamma_z==0){ tab = minimum(8,(unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max);}
	    if (flag_optimize_over_gamma_z==1){ tab = minimum(8,(unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub);}
	    array_printf(X_wSM___,"double",1,tab," %% X_wSM___: ");
	    array_printf(delta_x_wSM___,"double",1,tab," %% delta_x_wSM___: ");
	    array_printf(delta_y_wSM___,"double",1,tab," %% delta_y_wSM___: ");
	    array_printf(gamma_z_wSM___,"double",1,tab," %% gamma_z_wSM___: ");
	    if (flag_compute_I_value){ array_printf(I_value_wSM___,"double",1,tab," %% I_value_wSM___: ");}
	    /* if (flag_disp){ } */}
	  if (flag_dump){
	    MDA_ndim=0;
	    if (flag_optimize_over_gamma_z==0){ MDA_dim_[MDA_ndim++] = n_w_max;}
	    MDA_dim_[MDA_ndim++] = n_S;
	    MDA_dim_[MDA_ndim++] = n_M;
	    MDA_n_dim = MDA_ndim;
	    sprintf(MDA_fname,"mex_X_wSM___.mat");MDA_write_r8(MDA_n_dim,MDA_dim_,X_wSM___,MDA_fname);
	    sprintf(MDA_fname,"mex_delta_x_wSM___.mat");MDA_write_r8(MDA_n_dim,MDA_dim_,delta_x_wSM___,MDA_fname);
	    sprintf(MDA_fname,"mex_delta_y_wSM___.mat");MDA_write_r8(MDA_n_dim,MDA_dim_,delta_y_wSM___,MDA_fname);
	    sprintf(MDA_fname,"mex_gamma_z_wSM___.mat");MDA_write_r8(MDA_n_dim,MDA_dim_,gamma_z_wSM___,MDA_fname);
	    if (flag_compute_I_value){ sprintf(MDA_fname,"mex_I_value_wSM___.mat");MDA_write_r8(MDA_n_dim,MDA_dim_,I_value_wSM___,MDA_fname);}
	    /* if (flag_dump){ } */}
	  /* if (n_S_sub>0){ } */}
	/* for (nSbatch=0;nSbatch<n_Sbatch;nSbatch++){ } */}
      /* if (n_M_sub>0){ } */}
    /* for (nMbatch=0;nMbatch<n_Mbatch;nMbatch++){ } */}
  /* %%%%%%%%%%%%%%%% */
  if (verbose){
    printf(" %% elrt_sum_CTF_UX_S_k_q_nSw___     : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_CTF_UX_S_k_q_nSw___,n_op_sum_CTF_UX_S_k_q_nSw___/elrt_sum_CTF_UX_S_k_q_nSw___/1e6,n_op_sum_CTF_UX_S_k_q_nSw___/elrt_sum_CTF_UX_S_k_q_nSw___/1e9);
    printf(" %% elrt_sum_f_CTF_UX_S_k_q_nSw___   : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_f_CTF_UX_S_k_q_nSw___,n_op_sum_f_CTF_UX_S_k_q_nSw___/elrt_sum_f_CTF_UX_S_k_q_nSw___/1e6,n_op_sum_f_CTF_UX_S_k_q_nSw___/elrt_sum_f_CTF_UX_S_k_q_nSw___/1e9);
    printf(" %% elrt_sum_svd_VUXM_nMwl____       : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_svd_VUXM_nMwl____,n_op_sum_svd_VUXM_nMwl____/elrt_sum_svd_VUXM_nMwl____/1e6,n_op_sum_svd_VUXM_nMwl____/elrt_sum_svd_VUXM_nMwl____/1e9);
    printf(" %% elrt_sum_f_svd_VUXM_nMlw____     : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_f_svd_VUXM_nMlw____,n_op_sum_f_svd_VUXM_nMlw____/elrt_sum_f_svd_VUXM_nMlw____/1e6,n_op_sum_f_svd_VUXM_nMlw____/elrt_sum_f_svd_VUXM_nMlw____/1e9);
    printf(" %% elrt_sum_svd_SVUXM_SMwl____      : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_svd_SVUXM_SMwl____,n_op_sum_svd_SVUXM_SMwl____/elrt_sum_svd_SVUXM_SMwl____/1e6,n_op_sum_svd_SVUXM_SMwl____/elrt_sum_svd_SVUXM_SMwl____/1e9);
    printf(" %% elrt_sum_svd_SVUXM_SMlw____      : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_svd_SVUXM_SMlw____,n_op_sum_svd_SVUXM_SMlw____/elrt_sum_svd_SVUXM_SMlw____/1e6,n_op_sum_svd_SVUXM_SMlw____/elrt_sum_svd_SVUXM_SMlw____/1e9);
    printf(" %% elrt_sum_f_svd_SVUXM_SMlw____    : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_f_svd_SVUXM_SMlw____,n_op_sum_f_svd_SVUXM_SMlw____/elrt_sum_f_svd_SVUXM_SMlw____/1e6,n_op_sum_f_svd_SVUXM_SMlw____/elrt_sum_f_svd_SVUXM_SMlw____/1e9);
    printf(" %% elrt_sum_svd_SVUXM_0in_wlSM____  : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_svd_SVUXM_0in_wlSM____,n_op_sum_svd_SVUXM_0in_wlSM____/elrt_sum_svd_SVUXM_0in_wlSM____/1e6,n_op_sum_svd_SVUXM_0in_wlSM____/elrt_sum_svd_SVUXM_0in_wlSM____/1e9);
    printf(" %% elrt_sum_f_svd_SVUXM_0in_wSMl____: %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_f_svd_SVUXM_0in_wSMl____,n_op_sum_f_svd_SVUXM_0in_wSMl____/elrt_sum_f_svd_SVUXM_0in_wSMl____/1e6,n_op_sum_f_svd_SVUXM_0in_wSMl____/elrt_sum_f_svd_SVUXM_0in_wSMl____/1e9);
    printf(" %% elrt_sum_svd_SVUXM_out_wlSM____  : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_svd_SVUXM_out_wlSM____,n_op_sum_svd_SVUXM_out_wlSM____/elrt_sum_svd_SVUXM_out_wlSM____/1e6,n_op_sum_svd_SVUXM_out_wlSM____/elrt_sum_svd_SVUXM_out_wlSM____/1e9);
    printf(" %% elrt_sum_f_svd_SVUXM_out_wSMl____: %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_f_svd_SVUXM_out_wSMl____,n_op_sum_f_svd_SVUXM_out_wSMl____/elrt_sum_f_svd_SVUXM_out_wSMl____/1e6,n_op_sum_f_svd_SVUXM_out_wSMl____/elrt_sum_f_svd_SVUXM_out_wSMl____/1e9);
    printf(" %% elrt_sum_svd_SVUXM_lwSM____      : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_svd_SVUXM_lwSM____,n_op_sum_svd_SVUXM_lwSM____/elrt_sum_svd_SVUXM_lwSM____/1e6,n_op_sum_svd_SVUXM_lwSM____/elrt_sum_svd_SVUXM_lwSM____/1e9);
    printf(" %% elrt_sum_f_svd_SVUXM_lwSM____    : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_f_svd_SVUXM_lwSM____,n_op_sum_f_svd_SVUXM_lwSM____/elrt_sum_f_svd_SVUXM_lwSM____/1e6,n_op_sum_f_svd_SVUXM_lwSM____/elrt_sum_f_svd_SVUXM_lwSM____/1e9);
    printf(" %% elrt_sum_svd_SVUXM_lwsM____      : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_svd_SVUXM_lwsM____,n_op_sum_svd_SVUXM_lwsM____/elrt_sum_svd_SVUXM_lwsM____/1e6,n_op_sum_svd_SVUXM_lwsM____/elrt_sum_svd_SVUXM_lwsM____/1e9);
    printf(" %% elrt_sum_f_svd_SVUXM_lwsM____    : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_f_svd_SVUXM_lwsM____,n_op_sum_f_svd_SVUXM_lwsM____/elrt_sum_f_svd_SVUXM_lwsM____/1e6,n_op_sum_f_svd_SVUXM_lwsM____/elrt_sum_f_svd_SVUXM_lwsM____/1e9);
    printf(" %% elrt_sum_svd_USESVUXM_dwSM____   : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_svd_USESVUXM_dwSM____,n_op_sum_svd_USESVUXM_dwSM____/elrt_sum_svd_USESVUXM_dwSM____/1e6,n_op_sum_svd_USESVUXM_dwSM____/elrt_sum_svd_USESVUXM_dwSM____/1e9);
    printf(" %% elrt_sum_f_svd_USESVUXM_dwSM____ : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_f_svd_USESVUXM_dwSM____,n_op_sum_f_svd_USESVUXM_dwSM____/elrt_sum_f_svd_USESVUXM_dwSM____/1e6,n_op_sum_f_svd_USESVUXM_dwSM____/elrt_sum_f_svd_USESVUXM_dwSM____/1e9);
    printf(" %% elrt_sum_l2_dSM___               : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_l2_dSM___,n_op_sum_l2_dSM___/elrt_sum_l2_dSM___/1e6,n_op_sum_l2_dSM___/elrt_sum_l2_dSM___/1e9);
    printf(" %% elrt_sum_X_sub_dwSM____          : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_X_sub_dwSM____,n_op_sum_X_sub_dwSM____/elrt_sum_X_sub_dwSM____/1e6,n_op_sum_X_sub_dwSM____/elrt_sum_X_sub_dwSM____/1e9);
    printf(" %% elrt_sum_d_X_sub_dwSM____        : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_d_X_sub_dwSM____,n_op_sum_d_X_sub_dwSM____/elrt_sum_d_X_sub_dwSM____/1e6,n_op_sum_d_X_sub_dwSM____/elrt_sum_d_X_sub_dwSM____/1e9);
    printf(" %% elrt_sum_I_value_use_dwSM____    : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_I_value_use_dwSM____,n_op_sum_I_value_use_dwSM____/elrt_sum_I_value_use_dwSM____/1e6,n_op_sum_I_value_use_dwSM____/elrt_sum_I_value_use_dwSM____/1e9);
    printf(" %% elrt_sum_X_wSM___                : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_X_wSM___,n_op_sum_X_wSM___/elrt_sum_X_wSM___/1e6,n_op_sum_X_wSM___/elrt_sum_X_wSM___/1e9);
    /* if (verbose){ } */}
  /* %%%%%%%%%%%%%%%% */
  free(gamma_z_); gamma_z_=NULL;
  if (flag_compute_I_value){ free(I_value_use_dwSM____); I_value_use_dwSM____=NULL;}
  free(X_sub_dwSM____); X_sub_dwSM____=NULL;
  free(l2_dSM___); l2_dSM___=NULL;
  free(n2_dSM___); n2_dSM___=NULL;
  free(f2_dSM___); f2_dSM___=NULL;
  free(svd_USESVUXM_dwSM____); svd_USESVUXM_dwSM____=NULL;
  free(svd_SVUXM_lwsM____); svd_SVUXM_lwsM____=NULL;
  free(svd_SVUXM_lwSM____); svd_SVUXM_lwSM____=NULL;
  fftw_free(svd_SVUXM_0in_wlSM____); svd_SVUXM_0in_wlSM____=NULL;
  fftw_free(svd_SVUXM_out_wlSM____); svd_SVUXM_out_wlSM____=NULL;
  fftw_destroy_plan(fftw_plan_many_plan);
  fftwf_free(f_svd_SVUXM_0in_wSMl____); f_svd_SVUXM_0in_wSMl____=NULL;
  fftwf_free(f_svd_SVUXM_out_wSMl____); f_svd_SVUXM_out_wSMl____=NULL;
  fftwf_destroy_plan(fftwf_plan_many_plan);
  free(svd_SVUXM_SMwl____); svd_SVUXM_SMwl____=NULL;
  free(svd_SVUXM_SMlw____); svd_SVUXM_SMlw____=NULL;
  free(svd_VUXM_nMwl____); svd_VUXM_nMwl____=NULL;
  free(svd_VUXM_nMlw____); svd_VUXM_nMlw____=NULL;
  free(index_S_in_Sbatch_); index_S_in_Sbatch_=NULL;
  free(index_M_in_Mbatch_); index_M_in_Mbatch_=NULL;
  free(CTF_UX_S_k_q_nSw___); CTF_UX_S_k_q_nSw___=NULL;
  free(FTK_svd_U_d_expiw_s__); FTK_svd_U_d_expiw_s__=NULL;
  free(CTF_UX_S_k_q_wnS__); CTF_UX_S_k_q_wnS__=NULL;
  free(svd_VUXM_lwnM____); svd_VUXM_lwnM____=NULL;
  _mm_free(f_FTK_svd_U_d_expiw_s_real__); f_FTK_svd_U_d_expiw_s_real__=NULL;
  _mm_free(f_FTK_svd_U_d_expiw_s_imag__); f_FTK_svd_U_d_expiw_s_imag__=NULL;
  _mm_free(f_FTK_svd_U_d_expiw_s_tran_real__); f_FTK_svd_U_d_expiw_s_tran_real__=NULL;
  _mm_free(f_FTK_svd_U_d_expiw_s_tran_imag__); f_FTK_svd_U_d_expiw_s_tran_imag__=NULL;
  _mm_free(f_CTF_UX_S_k_q_nSw_real___); f_CTF_UX_S_k_q_nSw_real___=NULL;
  _mm_free(f_CTF_UX_S_k_q_nSw_imag___); f_CTF_UX_S_k_q_nSw_imag___=NULL;
  _mm_free(f_CTF_UX_S_k_q_wnS_real___); f_CTF_UX_S_k_q_wnS_real___=NULL;
  _mm_free(f_CTF_UX_S_k_q_wnS_imag___); f_CTF_UX_S_k_q_wnS_imag___=NULL;
  _mm_free(f_svd_VUXM_lwnM_real____); f_svd_VUXM_lwnM_real____=NULL;
  _mm_free(f_svd_VUXM_lwnM_imag____); f_svd_VUXM_lwnM_imag____=NULL;  
  _mm_free(f_svd_VUXM_nMlw_real____); f_svd_VUXM_nMlw_real____=NULL;
  _mm_free(f_svd_VUXM_nMlw_imag____); f_svd_VUXM_nMlw_imag____=NULL;  
  _mm_free(f_svd_SVUXM_SMlw_real____); f_svd_SVUXM_SMlw_real____=NULL;
  _mm_free(f_svd_SVUXM_SMlw_imag____); f_svd_SVUXM_SMlw_imag____=NULL;
  _mm_free(f_svd_SVUXM_wSMl_real____); f_svd_SVUXM_wSMl_real____=NULL;
  _mm_free(f_svd_SVUXM_wSMl_imag____); f_svd_SVUXM_wSMl_imag____=NULL;
  _mm_free(f_svd_SVUXM_lwsM_real____); f_svd_SVUXM_lwsM_real____=NULL;
  _mm_free(f_svd_SVUXM_lwsM_imag____); f_svd_SVUXM_lwsM_imag____=NULL;
  _mm_free(f_svd_SVUXM_lwSM_real____); f_svd_SVUXM_lwSM_real____=NULL;
  _mm_free(f_svd_SVUXM_lwSM_imag____); f_svd_SVUXM_lwSM_imag____=NULL;
  _mm_free(f_svd_USESVUXM_dwSM_real____); f_svd_USESVUXM_dwSM_real____=NULL;
  _mm_free(f_svd_USESVUXM_dwSM_imag____); f_svd_USESVUXM_dwSM_imag____=NULL;
  free(d_X_sub_dwSM____); d_X_sub_dwSM____=NULL;
  /* %%%%%%%%%%%%%%%% */
  if (verbose){ printf(" %% [finished mex_ampmh_X_wSM___12]\n");}
}

void mex_ampmh_X_wSM___13
(
  int n_M_per_Mbatch
 ,int n_S_per_Sbatch
 ,int flag_optimize_over_gamma_z
 ,int flag_compute_I_value
 ,double tolerance_master
 ,int FTK_n_svd_l
 ,int FTK_n_delta_v
 ,double *FTK_svd_U_d_expiw_s_real__
 ,double *FTK_svd_U_d_expiw_s_imag__
 ,double *FTK_delta_x_
 ,double *FTK_delta_y_
 ,int n_w_max
 ,int pm_n_UX_rank
 ,int n_S
 ,double *CTF_UX_S_k_q_wnS_real__
 ,double *CTF_UX_S_k_q_wnS_imag__
 ,double *CTF_UX_S_l2_
 ,int n_M
 ,double *svd_VUXM_lwnM_real____
 ,double *svd_VUXM_lwnM_imag____
 ,double *UX_M_l2_dM__
 ,double *X_wSM___
 ,double *delta_x_wSM___
 ,double *delta_y_wSM___
 ,double *gamma_z_wSM___
 ,double *I_value_wSM___
 )
{
  int verbose=1;
  int flag_disp=0;
  int flag_dump=0;
  int flag_slow=0;
  /* %%%%; */
  int n_M_per_Mbatch_rup=0,n_M_per_Mbatch_256=0;
  int n_S_per_Sbatch_rup=0,n_S_per_Sbatch_256=0;
  int FTK_n_svd_l_rup=0,FTK_n_svd_l_256=0;
  int FTK_n_delta_v_rup=0,FTK_n_delta_v_256=0;
  int n_w_max_rup=0,n_w_max_256=0;
  int pm_n_UX_rank_rup=0,pm_n_UX_rank_256=0;
  int n_S_rup=0,n_S_256=0;
  int n_S_sub_rup=0,n_S_sub_256=0;
  int n_M_rup=0,n_M_256=0;
  int n_M_sub_rup=0,n_M_sub_256=0;
  /* %%%%; */
  double complex *FTK_svd_U_d_expiw_s__=NULL;
  float *f_FTK_svd_U_d_expiw_s_real__=NULL;
  float *f_FTK_svd_U_d_expiw_s_imag__=NULL;
  float *f_FTK_svd_U_d_expiw_s_tran_real__=NULL;
  float *f_FTK_svd_U_d_expiw_s_tran_imag__=NULL;
  double complex *CTF_UX_S_k_q_wnS__=NULL;
  float *f_CTF_UX_S_k_q_wnS_real___=NULL;
  float *f_CTF_UX_S_k_q_wnS_imag___=NULL;
  float *f_CTF_UX_S_k_q_nSw_real___=NULL;
  float *f_CTF_UX_S_k_q_nSw_imag___=NULL;
  double complex *svd_VUXM_lwnM____=NULL;
  float *f_svd_VUXM_lwnM_real____=NULL;
  float *f_svd_VUXM_lwnM_imag____=NULL;
  float *f_svd_VUXM_nMlw_real____=NULL;
  float *f_svd_VUXM_nMlw_imag____=NULL;
  unsigned long long int na=0,n_a=0;
  /* %%%% */
  int nMbatch=0,n_Mbatch=0,nSbatch=0,n_Sbatch=0;
  int *index_M_in_Mbatch_=NULL;
  int *index_S_in_Sbatch_=NULL;
  int n_M_sub=0,nM_sub=0;
  int n_S_sub=0,nS_sub=0;
  int nl=0,ndelta_v=0,nw=0,ndw=0,nS=0,nM=0;
  int pm_nUX_rank=0;
  unsigned long long int tab=0,tabA=0,tabB=0,tabC=0;
  double complex *CTF_UX_S_k_q_nSw___=NULL;
  double complex *svd_VUXM_nMwl____=NULL;
  double complex *svd_VUXM_nMlw____=NULL;
  double complex *svd_SVUXM_SMwl____=NULL;
  double complex *svd_SVUXM_SMlw____=NULL;
  float *f_svd_SVUXM_SMlw_real____=NULL;
  float *f_svd_SVUXM_SMlw_imag____=NULL;
  float *f_svd_SVUXM_wSMl_real____=NULL;
  float *f_svd_SVUXM_wSMl_imag____=NULL;
  fftw_complex *svd_SVUXM_0in_wlSM____;
  fftw_complex *svd_SVUXM_out_wlSM____;
  fftwf_complex *f_svd_SVUXM_0in_wSMl____;
  fftwf_complex *f_svd_SVUXM_out_wSMl____;
  double complex *svd_SVUXM_lwSM____=NULL;
  double complex *svd_SVUXM_lwsM____=NULL;
  float *f_svd_SVUXM_lwSM_real____=NULL;
  float *f_svd_SVUXM_lwSM_imag____=NULL;
  float *f_svd_SVUXM_lwsM_real____=NULL;
  float *f_svd_SVUXM_lwsM_imag____=NULL;
  double complex *svd_USESVUXM_dwSM____=NULL;
  float *f_svd_USESVUXM_dwSM_real____=NULL;
  float *f_svd_USESVUXM_dwSM_imag____=NULL;
  double complex cblas_alpha = (double complex) 1.0;
  double complex cblas_beta  = (double complex) 0.0;
  double *l2_dSM___=NULL;
  double *n2_dSM___=NULL;
  double *f2_dSM___=NULL;
  double *X_sub_dwSM____=NULL;
  double *d_X_sub_dwSM____=NULL;
  double *I_value_use_dwSM____=NULL;
  double *gamma_z_=NULL;
  float *f_CR_=NULL,*f_CI_=NULL,*f_C_=NULL;
  /* %%%% */
  fftw_plan fftw_plan_many_plan;
  int fftw_plan_many_rank=0;
  int fftw_plan_many_n[1]={0};
  int fftw_plan_many_howmany=0;
  int fftw_plan_many_istride=0;
  int fftw_plan_many_idist=0;
  int fftw_plan_many_ostride=0;
  int fftw_plan_many_odist=0;
  int fftw_plan_many_sign=0;
  unsigned int fftw_plan_many_flags=0;
  /* %%%% */
  fftwf_plan fftwf_plan_many_plan;
  int fftwf_plan_many_rank=0;
  int fftwf_plan_many_n[1]={0};
  int fftwf_plan_many_howmany=0;
  int fftwf_plan_many_istride=0;
  int fftwf_plan_many_idist=0;
  int fftwf_plan_many_ostride=0;
  int fftwf_plan_many_odist=0;
  int fftwf_plan_many_sign=0;
  unsigned int fftwf_plan_many_flags=0;
  /* %%%% */
  int MDA_n_dim=0,MDA_ndim = 0;
  int MDA_dim_[8];
  char MDA_fname[32];
  /* %%%% */
  int ntick=0,n_tick = 8;
  clock_t t_start_[n_tick];
  clock_t t_final_[n_tick];
  struct timeval d_start_[n_tick];
  struct timeval d_final_[n_tick];
  long l_msec_[n_tick],l_ssec_[n_tick],l_usec_[n_tick];
  double elct_[n_tick],elrt_[n_tick];
  double elrt_sum_CTF_UX_S_k_q_nSw___=0;
  double n_op_sum_CTF_UX_S_k_q_nSw___=0;
  double elrt_sum_f_CTF_UX_S_k_q_nSw___=0;
  double n_op_sum_f_CTF_UX_S_k_q_nSw___=0;
  double elrt_sum_svd_VUXM_nMwl____=0;
  double n_op_sum_svd_VUXM_nMwl____=0;
  double elrt_sum_f_svd_VUXM_nMlw____=0;
  double n_op_sum_f_svd_VUXM_nMlw____=0;
  double elrt_sum_svd_SVUXM_SMwl____=0;
  double n_op_sum_svd_SVUXM_SMwl____=0;
  double elrt_sum_svd_SVUXM_SMlw____=0;
  double n_op_sum_svd_SVUXM_SMlw____=0;
  double elrt_sum_f_svd_SVUXM_SMlw____=0;
  double n_op_sum_f_svd_SVUXM_SMlw____=0;
  double elrt_sum_svd_SVUXM_0in_wlSM____=0;
  double n_op_sum_svd_SVUXM_0in_wlSM____=0;
  double elrt_sum_f_svd_SVUXM_0in_wSMl____=0;
  double n_op_sum_f_svd_SVUXM_0in_wSMl____=0;
  double elrt_sum_svd_SVUXM_out_wlSM____=0;
  double n_op_sum_svd_SVUXM_out_wlSM____=0;
  double elrt_sum_f_svd_SVUXM_out_wSMl____=0;
  double n_op_sum_f_svd_SVUXM_out_wSMl____=0;
  double elrt_sum_svd_SVUXM_lwSM____=0;
  double n_op_sum_svd_SVUXM_lwSM____=0;
  double elrt_sum_f_svd_SVUXM_lwSM____=0;
  double n_op_sum_f_svd_SVUXM_lwSM____=0;
  double elrt_sum_svd_SVUXM_lwsM____=0;
  double n_op_sum_svd_SVUXM_lwsM____=0;
  double elrt_sum_f_svd_SVUXM_lwsM____=0;
  double n_op_sum_f_svd_SVUXM_lwsM____=0;
  double elrt_sum_svd_USESVUXM_dwSM____=0;
  double n_op_sum_svd_USESVUXM_dwSM____=0;
  double elrt_sum_f_svd_USESVUXM_dwSM____=0;
  double n_op_sum_f_svd_USESVUXM_dwSM____=0;
  double elrt_sum_l2_dSM___=0;
  double n_op_sum_l2_dSM___=0;
  double elrt_sum_X_sub_dwSM____=0;
  double n_op_sum_X_sub_dwSM____=0;
  double elrt_sum_d_X_sub_dwSM____=0;
  double n_op_sum_d_X_sub_dwSM____=0;
  double elrt_sum_I_value_use_dwSM____=0;
  double n_op_sum_I_value_use_dwSM____=0;
  double elrt_sum_X_wSM___=0;
  double n_op_sum_X_wSM___=0;
  /* %%%% */
  if (verbose){ printf(" %% [entering mex_ampmh_X_wSM___13]\n");}
  /* %%%%%%%%%%%%%%%% */
  if (flag_slow){
    n_a = FTK_n_delta_v*FTK_n_svd_l;
    FTK_svd_U_d_expiw_s__ = double_complex_malloc_and_interleave(n_a,FTK_svd_U_d_expiw_s_real__,FTK_svd_U_d_expiw_s_imag__);
    n_a = n_w_max*pm_n_UX_rank*n_S;
    CTF_UX_S_k_q_wnS__ = double_complex_malloc_and_interleave(n_a,CTF_UX_S_k_q_wnS_real__,CTF_UX_S_k_q_wnS_imag__);
    n_a = FTK_n_svd_l*n_w_max*pm_n_UX_rank*n_M;
    svd_VUXM_lwnM____ = double_complex_malloc_and_interleave(n_a,svd_VUXM_lwnM_real____,svd_VUXM_lwnM_imag____);
    /* if (flag_slow){ } */}
  /* %%%%%%%%%%%%%%%% */
  FTK_n_delta_v_rup = rup(FTK_n_delta_v,8); FTK_n_delta_v_256 = FTK_n_delta_v_rup/8;
  FTK_n_svd_l_rup = rup(FTK_n_svd_l,8); FTK_n_svd_l_256 = FTK_n_svd_l_rup/8;
  n_w_max_rup = rup(n_w_max,8); n_w_max_256 = n_w_max_rup/8;
  pm_n_UX_rank_rup = rup(pm_n_UX_rank,8); pm_n_UX_rank_256 = pm_n_UX_rank_rup/8;
  n_S_rup = rup(n_S,8); n_S_256 = n_S_rup/8;
  n_S_per_Sbatch_rup = rup(n_S_per_Sbatch,8); n_S_per_Sbatch_256 = n_S_per_Sbatch_rup/8;
  n_M_rup = rup(n_M,8); n_M_256 = n_M_rup/8;
  n_M_per_Mbatch_rup = rup(n_M_per_Mbatch,8); n_M_per_Mbatch_256 = n_M_per_Mbatch_rup/8;
  f_FTK_svd_U_d_expiw_s_real__ = float__m256_malloc_from_double__(FTK_n_delta_v,FTK_n_svd_l,FTK_svd_U_d_expiw_s_real__);
  f_FTK_svd_U_d_expiw_s_imag__ = float__m256_malloc_from_double__(FTK_n_delta_v,FTK_n_svd_l,FTK_svd_U_d_expiw_s_imag__);
  f_FTK_svd_U_d_expiw_s_tran_real__ = float__m256_malloc_from_double__(FTK_n_svd_l,FTK_n_delta_v,NULL);
  f_FTK_svd_U_d_expiw_s_tran_imag__ = float__m256_malloc_from_double__(FTK_n_svd_l,FTK_n_delta_v,NULL);
  transpose_ps_block1(
		      f_FTK_svd_U_d_expiw_s_real__
		      ,f_FTK_svd_U_d_expiw_s_tran_real__
		      ,FTK_n_delta_v_rup
		      ,FTK_n_svd_l_rup
		      ,FTK_n_delta_v_rup
		      ,FTK_n_svd_l_rup
		      ,16
		      );
  transpose_ps_block1(
		      f_FTK_svd_U_d_expiw_s_imag__
		      ,f_FTK_svd_U_d_expiw_s_tran_imag__
		      ,FTK_n_delta_v_rup
		      ,FTK_n_svd_l_rup
		      ,FTK_n_delta_v_rup
		      ,FTK_n_svd_l_rup
		      ,16
		      );
  if (verbose>2){
    array_sub_printf(f_FTK_svd_U_d_expiw_s_real__     ,"float" ,FTK_n_delta_v_rup,3,FTK_n_svd_l_rup,4," %% f_FTK_svd_U_d_expiw_s_real__     : ");
    array_sub_printf(f_FTK_svd_U_d_expiw_s_tran_real__,"float" ,FTK_n_svd_l_rup,4,FTK_n_delta_v_rup,3," %% f_FTK_svd_U_d_expiw_s_tran_real__: ");
    array_sub_printf(f_FTK_svd_U_d_expiw_s_imag__     ,"float" ,FTK_n_delta_v_rup,3,FTK_n_svd_l_rup,4," %% f_FTK_svd_U_d_expiw_s_imag__     : ");
    array_sub_printf(f_FTK_svd_U_d_expiw_s_tran_imag__,"float" ,FTK_n_svd_l_rup,4,FTK_n_delta_v_rup,3," %% f_FTK_svd_U_d_expiw_s_tran_imag__: ");
    exit(0);
    /* if (verbose>2){ } */}
  f_CTF_UX_S_k_q_wnS_real___ = float___m256_malloc_from_double___(n_w_max,pm_n_UX_rank,n_S,CTF_UX_S_k_q_wnS_real__);
  f_CTF_UX_S_k_q_wnS_imag___ = float___m256_malloc_from_double___(n_w_max,pm_n_UX_rank,n_S,CTF_UX_S_k_q_wnS_imag__);
  f_CTF_UX_S_k_q_nSw_real___ = float___m256_malloc_from_double___(pm_n_UX_rank,n_S,n_w_max,NULL);
  f_CTF_UX_S_k_q_nSw_imag___ = float___m256_malloc_from_double___(pm_n_UX_rank,n_S,n_w_max,NULL);
  f_svd_VUXM_lwnM_real____ = float____m256_malloc_from_double____(FTK_n_svd_l,n_w_max,pm_n_UX_rank,n_M,svd_VUXM_lwnM_real____);
  f_svd_VUXM_lwnM_imag____ = float____m256_malloc_from_double____(FTK_n_svd_l,n_w_max,pm_n_UX_rank,n_M,svd_VUXM_lwnM_imag____);
  f_svd_VUXM_nMlw_real____ = float____m256_malloc_from_double____(FTK_n_svd_l,n_w_max,pm_n_UX_rank,n_M_per_Mbatch,NULL);
  f_svd_VUXM_nMlw_imag____ = float____m256_malloc_from_double____(FTK_n_svd_l,n_w_max,pm_n_UX_rank,n_M_per_Mbatch,NULL);
  if (verbose>2){
    /* %%%% */
    array_sub_printf(  FTK_svd_U_d_expiw_s_real__,"double",FTK_n_delta_v    ,3,FTK_n_svd_l    ,4," %%   FTK_svd_U_d_expiw_s_real__: ");
    array_sub_printf(f_FTK_svd_U_d_expiw_s_real__,"float" ,FTK_n_delta_v_rup,3,FTK_n_svd_l_rup,4," %% f_FTK_svd_U_d_expiw_s_real__: ");
    array_sub_printf(  FTK_svd_U_d_expiw_s_imag__,"double",FTK_n_delta_v    ,3,FTK_n_svd_l    ,4," %%   FTK_svd_U_d_expiw_s_imag__: ");
    array_sub_printf(f_FTK_svd_U_d_expiw_s_imag__,"float" ,FTK_n_delta_v_rup,3,FTK_n_svd_l_rup,4," %% f_FTK_svd_U_d_expiw_s_imag__: ");
    /* %%%% */
    array_sub_printf(  CTF_UX_S_k_q_wnS_real__ ,"double",n_w_max    *pm_n_UX_rank    ,3,n_S    ,4," %%   CTF_UX_S_k_q_wnS_real__ : ");
    array_sub_printf(f_CTF_UX_S_k_q_wnS_real___,"float" ,n_w_max_rup*pm_n_UX_rank_rup,3,n_S_rup,4," %% f_CTF_UX_S_k_q_wnS_real___: ");
    array_sub_printf(  CTF_UX_S_k_q_wnS_imag__ ,"double",n_w_max    *pm_n_UX_rank    ,3,n_S    ,4," %%   CTF_UX_S_k_q_wnS_imag__ : ");
    array_sub_printf(f_CTF_UX_S_k_q_wnS_imag___,"float" ,n_w_max_rup*pm_n_UX_rank_rup,3,n_S_rup,4," %% f_CTF_UX_S_k_q_wnS_imag___: ");
    /* %%%% */
    array_sub_printf(  svd_VUXM_lwnM_real____,"double",FTK_n_svd_l    *n_w_max    *pm_n_UX_rank    ,3,n_M    ,4," %%   svd_VUXM_lwnM_real____: ");
    array_sub_printf(f_svd_VUXM_lwnM_real____,"float" ,FTK_n_svd_l_rup*n_w_max_rup*pm_n_UX_rank_rup,3,n_M_rup,4," %% f_svd_VUXM_lwnM_real____: ");
    array_sub_printf(  svd_VUXM_lwnM_imag____,"double",FTK_n_svd_l    *n_w_max    *pm_n_UX_rank    ,3,n_M    ,4," %%   svd_VUXM_lwnM_imag____: ");
    array_sub_printf(f_svd_VUXM_lwnM_imag____,"float" ,FTK_n_svd_l_rup*n_w_max_rup*pm_n_UX_rank_rup,3,n_M_rup,4," %% f_svd_VUXM_lwnM_imag____: ");
    /* %%%% */
    array_sub_printf(  svd_VUXM_lwnM_real____,"double",FTK_n_svd_l    *n_w_max    ,3,pm_n_UX_rank    *n_M    ,4," %%   svd_VUXM_lwnM_real____: ");
    array_sub_printf(f_svd_VUXM_lwnM_real____,"float" ,FTK_n_svd_l_rup*n_w_max_rup,3,pm_n_UX_rank_rup*n_M_rup,4," %% f_svd_VUXM_lwnM_real____: ");
    array_sub_printf(  svd_VUXM_lwnM_imag____,"double",FTK_n_svd_l    *n_w_max    ,3,pm_n_UX_rank    *n_M    ,4," %%   svd_VUXM_lwnM_imag____: ");
    array_sub_printf(f_svd_VUXM_lwnM_imag____,"float" ,FTK_n_svd_l_rup*n_w_max_rup,3,pm_n_UX_rank_rup*n_M_rup,4," %% f_svd_VUXM_lwnM_imag____: ");
    /* %%%% */
    exit(0);
    /*if (verbose>2){ } */}
  /* %%%%%%%%%%%%%%%% */
  if (verbose>1){ printf(" %% n_M_per_Mbatch %d\n",n_M_per_Mbatch);}
  if (verbose>1){ printf(" %% n_S_per_Sbatch %d\n",n_S_per_Sbatch);}
  if (verbose>1){ printf(" %% flag_optimize_over_gamma_z %d\n",flag_optimize_over_gamma_z);}
  if (verbose>1){ printf(" %% flag_compute_I_value %d\n",flag_compute_I_value);}
  if (verbose>1){ printf(" %% tolerance_master %0.16f\n",tolerance_master);}
  if (verbose>1){ printf(" %% FTK_n_svd_l %d\n",FTK_n_svd_l);}
  if (verbose>1){ printf(" %% FTK_n_delta_v %d\n",FTK_n_delta_v);}
  n_a = FTK_n_delta_v*FTK_n_svd_l;
  if (flag_slow){
    if (verbose>1){ printf(" %% FTK_svd_U_d_expiw_s__ %+0.16f %+0.16f*i --> %+0.16f %+0.16f*i\n",creal(FTK_svd_U_d_expiw_s__[0]),cimag(FTK_svd_U_d_expiw_s__[0]),creal(FTK_svd_U_d_expiw_s__[n_a-1]),cimag(FTK_svd_U_d_expiw_s__[n_a-1]));}
    /* if (flag_slow){ } */}
  if (verbose>1){ printf(" %% FTK_delta_x_ %0.16f --> %0.16f\n",FTK_delta_x_[0],FTK_delta_x_[FTK_n_delta_v-1]);}
  if (verbose>1){ printf(" %% FTK_delta_y_ %0.16f --> %0.16f\n",FTK_delta_y_[0],FTK_delta_y_[FTK_n_delta_v-1]);}
  if (verbose>1){ printf(" %% n_w_max %d\n",n_w_max);}
  if (verbose>1){ printf(" %% pm_n_UX_rank %d\n",pm_n_UX_rank);}
  if (verbose>1){ printf(" %% n_S %d\n",n_S);}
  if (flag_slow){
    n_a = n_w_max*pm_n_UX_rank*n_S;
    if (verbose>1){ printf(" %% CTF_UX_S_k_q_wnS__ %+0.16f %+0.16f*i --> %+0.16f %+0.16f*i\n",creal(CTF_UX_S_k_q_wnS__[0]),cimag(CTF_UX_S_k_q_wnS__[0]),creal(CTF_UX_S_k_q_wnS__[n_a-1]),cimag(CTF_UX_S_k_q_wnS__[n_a-1]));}
    /* if (flag_slow){ } */}
  if (verbose>1){ printf(" %% CTF_UX_S_l2_ %0.16f --> %0.16f\n",CTF_UX_S_l2_[0],CTF_UX_S_l2_[n_S-1]);}
  if (verbose>1){ printf(" %% n_M %d\n",n_M);}
  if (flag_slow){
    n_a = FTK_n_svd_l*n_w_max*pm_n_UX_rank*n_M;
    if (verbose>1){ printf(" %% svd_VUXM_lwnM____ %+0.16f %+0.16f*i --> %+0.16f %+0.16f*i\n",creal(svd_VUXM_lwnM____[0]),cimag(svd_VUXM_lwnM____[0]),creal(svd_VUXM_lwnM____[n_a-1]),cimag(svd_VUXM_lwnM____[n_a-1]));}
    /* if (flag_slow){ } */}
  if (verbose>1){ printf(" %% UX_M_l2_dM__ %0.16f --> %0.16f\n",UX_M_l2_dM__[0],UX_M_l2_dM__[FTK_n_delta_v*n_M-1]);}
  /* %%%% */
  /* if ( (nargout>5) & (isempty(I_value_wSM___)) ); I_value_wSM___ = ones(size(X_wSM___)); end; */
  /* %%%% */
  if (flag_compute_I_value){
    if (flag_optimize_over_gamma_z==0){ n_a = (unsigned long long int)n_M*(unsigned long long int)n_S*(unsigned long long int)n_w_max;}
    if (flag_optimize_over_gamma_z==1){ n_a = (unsigned long long int)n_M*(unsigned long long int)n_S;}
    for (na=0;na<n_a;na++){ I_value_wSM___[na] = 1.0;}
    /* if (flag_compute_I_value){ } */}
  /* %%%% */
  /* CTF_UX_S_k_q_nSw___ = permute(CTF_UX_S_k_q_wnS___,[2,3,1]); */
  /* %%%% */
  if (flag_slow){
    local_tic(0,t_start_,d_start_);
    tab = (unsigned long long int)pm_n_UX_rank*(unsigned long long int)n_S*(unsigned long long int)n_w_max;
    CTF_UX_S_k_q_nSw___ = (double complex *) malloc(tab*sizeof(double complex));
    na=0;
    for (nS=0;nS<n_S;nS++){
      for (pm_nUX_rank=0;pm_nUX_rank<pm_n_UX_rank;pm_nUX_rank++){
	for (nw=0;nw<n_w_max;nw++){
	  /* tabA = */
	  /*   (unsigned long long int)pm_nUX_rank */
	  /*   + (unsigned long long int)nS*(unsigned long long int)pm_n_UX_rank */
	  /*   + (unsigned long long int)nw*(unsigned long long int)n_S*(unsigned long long int)pm_n_UX_rank; */
	  tabA = (unsigned long long int)pm_nUX_rank + ((unsigned long long int)nS + (unsigned long long int)nw*(unsigned long long int)n_S)*(unsigned long long int)pm_n_UX_rank;
	  CTF_UX_S_k_q_nSw___[tabA] = CTF_UX_S_k_q_wnS__[na]; na++;
	  /* for (nw=0;nw<n_w_max;nw++){ } */}
	/* for (pm_nUX_rank=0;pm_nUX_rank<pm_n_UX_rank;pm_nUX_rank++){ } */}
      /* for (nS=0;nS<n_S;nS++){ } */}
    local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"CTF_UX_S_k_q_nSw___: ");
    elrt_sum_CTF_UX_S_k_q_nSw___ += elrt_[0];
    n_op_sum_CTF_UX_S_k_q_nSw___ += tab;
    /* if (flag_slow){ } */}
  /* %%%%%%%%%%%%%%%% */
  local_tic(0,t_start_,d_start_);
  tab = (unsigned long long int)pm_n_UX_rank*(unsigned long long int)n_S*(unsigned long long int)n_w_max;
  transpose_ps_block1(
		      f_CTF_UX_S_k_q_wnS_real___
		      ,f_CTF_UX_S_k_q_nSw_real___
		      ,n_w_max_rup
		      ,pm_n_UX_rank_rup*n_S_rup
		      ,n_w_max_rup
		      ,pm_n_UX_rank_rup*n_S_rup
		      ,16
		      );
  transpose_ps_block1(
		      f_CTF_UX_S_k_q_wnS_imag___
		      ,f_CTF_UX_S_k_q_nSw_imag___
		      ,n_w_max_rup
		      ,pm_n_UX_rank_rup*n_S_rup
		      ,n_w_max_rup
		      ,pm_n_UX_rank_rup*n_S_rup
		      ,16
		      );
  if (verbose>2 && flag_slow){
    array_sub_printf(  CTF_UX_S_k_q_nSw___ ,"double complex",pm_n_UX_rank    *n_S    ,4,n_w_max    ,3," %%   CTF_UX_S_k_q_nSw___ : ");
    array_sub_printf(f_CTF_UX_S_k_q_nSw_real___,"float" ,pm_n_UX_rank_rup*n_S_rup,4,n_w_max_rup,3," %% f_CTF_UX_S_k_q_nSw_real___: ");
    array_sub_printf(f_CTF_UX_S_k_q_nSw_imag___,"float" ,pm_n_UX_rank_rup*n_S_rup,4,n_w_max_rup,3," %% f_CTF_UX_S_k_q_nSw_imag___: ");
    exit(0);
    /* if (verbose>2){ } */}  
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"f_CTF_UX_S_k_q_nSw___: ");
  elrt_sum_f_CTF_UX_S_k_q_nSw___ += elrt_[0];
  n_op_sum_f_CTF_UX_S_k_q_nSw___ += tab;
  /* %%%%%%%%%%%%%%%% */
  if (flag_slow){
    tab = (unsigned long long int)pm_n_UX_rank*(unsigned long long int)n_M_per_Mbatch*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_svd_l;
    svd_VUXM_nMwl____ = (double complex *) malloc(tab*sizeof(double complex));
    svd_VUXM_nMlw____ = (double complex *) malloc(tab*sizeof(double complex));
    tab = (unsigned long long int)n_S*(unsigned long long int)n_M_per_Mbatch*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_svd_l;
    svd_SVUXM_SMwl____ = (double complex *) malloc(tab*sizeof(double complex));
    svd_SVUXM_SMlw____ = (double complex *) malloc(tab*sizeof(double complex));
    /* if (flag_slow){ } */}
  f_svd_SVUXM_SMlw_real____ = float____m256_malloc_from_double____(n_S,n_M_per_Mbatch,n_w_max,FTK_n_svd_l,NULL);
  f_svd_SVUXM_SMlw_imag____ = float____m256_malloc_from_double____(n_S,n_M_per_Mbatch,n_w_max,FTK_n_svd_l,NULL);
  f_svd_SVUXM_wSMl_real____ = float____m256_malloc_from_double____(n_S,n_M_per_Mbatch,n_w_max,FTK_n_svd_l,NULL);
  f_svd_SVUXM_wSMl_imag____ = float____m256_malloc_from_double____(n_S,n_M_per_Mbatch,n_w_max,FTK_n_svd_l,NULL);
  /* %%%% */
  if (flag_slow){
    tab = (unsigned long long int)n_w_max*(unsigned long long int)FTK_n_svd_l*(unsigned long long int)n_S*(unsigned long long int)n_M_per_Mbatch;
    svd_SVUXM_0in_wlSM____ = (fftw_complex *) fftw_malloc(tab*sizeof(fftw_complex));
    svd_SVUXM_out_wlSM____ = (fftw_complex *) fftw_malloc(tab*sizeof(fftw_complex));
    fftw_plan_many_rank = 1;
    fftw_plan_many_n[0] = n_w_max;
    fftw_plan_many_howmany = (unsigned long long int)FTK_n_svd_l*(unsigned long long int)n_S*(unsigned long long int)n_M_per_Mbatch;
    fftw_plan_many_istride = 1;
    fftw_plan_many_idist = n_w_max;
    fftw_plan_many_ostride = 1;
    fftw_plan_many_odist = n_w_max;
    fftw_plan_many_sign = FFTW_BACKWARD;
    fftw_plan_many_flags = FFTW_ESTIMATE;
    fftw_plan_many_plan = fftw_plan_many_dft(
					     fftw_plan_many_rank
					     ,fftw_plan_many_n
					     ,fftw_plan_many_howmany
					     ,svd_SVUXM_0in_wlSM____
					     ,NULL
					     ,fftw_plan_many_istride
					     ,fftw_plan_many_idist
					     ,svd_SVUXM_out_wlSM____
					     ,NULL
					     ,fftw_plan_many_ostride
					     ,fftw_plan_many_odist
					     ,fftw_plan_many_sign
					     ,fftw_plan_many_flags
					     );
    /* if (flag_slow){ } */}
  /* %%%% */
  tab = (unsigned long long int)n_w_max_rup*(unsigned long long int)n_S_rup*(unsigned long long int)n_M_per_Mbatch_rup*(unsigned long long int)FTK_n_svd_l_rup;
  f_svd_SVUXM_0in_wSMl____ = (fftwf_complex *) fftwf_malloc(tab*sizeof(fftwf_complex));
  f_svd_SVUXM_out_wSMl____ = (fftwf_complex *) fftwf_malloc(tab*sizeof(fftwf_complex));
  fftwf_plan_many_rank = 1;
  fftwf_plan_many_n[0] = n_w_max;
  fftwf_plan_many_howmany = (unsigned long long int)n_S_rup*(unsigned long long int)n_M_per_Mbatch_rup*(unsigned long long int)FTK_n_svd_l_rup;
  fftwf_plan_many_istride = 1;
  fftwf_plan_many_idist = n_w_max_rup;
  fftwf_plan_many_ostride = 1;
  fftwf_plan_many_odist = n_w_max_rup;
  fftwf_plan_many_sign = FFTW_BACKWARD;
  fftwf_plan_many_flags = FFTW_ESTIMATE;
  fftwf_plan_many_plan = fftwf_plan_many_dft(
					   fftwf_plan_many_rank
					   ,fftwf_plan_many_n
					   ,fftwf_plan_many_howmany
					   ,f_svd_SVUXM_0in_wSMl____
					   ,NULL
					   ,fftwf_plan_many_istride
					   ,fftwf_plan_many_idist
					   ,f_svd_SVUXM_out_wSMl____
					   ,NULL
					   ,fftwf_plan_many_ostride
					   ,fftwf_plan_many_odist
					   ,fftwf_plan_many_sign
					   ,fftwf_plan_many_flags
					   );
  /* %%%% */
  if (flag_slow){
    tab = (unsigned long long int)FTK_n_svd_l*(unsigned long long int)n_w_max*(unsigned long long int)n_S*(unsigned long long int)n_M_per_Mbatch;
    svd_SVUXM_lwSM____ = (double complex *) malloc(tab*sizeof(double complex));
    /* if (flag_slow){ } */}
  f_svd_SVUXM_lwSM_real____ = float____m256_malloc_from_double____(FTK_n_svd_l,n_w_max,n_S,n_M_per_Mbatch,NULL);
  f_svd_SVUXM_lwSM_imag____ = float____m256_malloc_from_double____(FTK_n_svd_l,n_w_max,n_S,n_M_per_Mbatch,NULL);
  if (flag_slow){
    tab = (unsigned long long int)n_S_per_Sbatch*(unsigned long long int)n_M_per_Mbatch*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_svd_l;
    svd_SVUXM_lwsM____ = (double complex *) malloc(tab*sizeof(double complex));
    /* if (flag_slow){ } */}
  f_svd_SVUXM_lwsM_real____ = float____m256_malloc_from_double____(FTK_n_svd_l,n_w_max,n_S_per_Sbatch,n_M_per_Mbatch,NULL);
  f_svd_SVUXM_lwsM_imag____ = float____m256_malloc_from_double____(FTK_n_svd_l,n_w_max,n_S_per_Sbatch,n_M_per_Mbatch,NULL);
  if (flag_slow){
    tab = (unsigned long long int)n_S_per_Sbatch*(unsigned long long int)n_M_per_Mbatch*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_delta_v;
    svd_USESVUXM_dwSM____ = (double complex *) malloc(tab*sizeof(double complex));
    /* if (flag_slow){ } */}
  f_svd_USESVUXM_dwSM_real____ = float____m256_malloc_from_double____(FTK_n_delta_v,n_w_max,n_S_per_Sbatch,n_M_per_Mbatch,NULL);
  f_svd_USESVUXM_dwSM_imag____ = float____m256_malloc_from_double____(FTK_n_delta_v,n_w_max,n_S_per_Sbatch,n_M_per_Mbatch,NULL);
  tab = (unsigned long long int)n_S_per_Sbatch*(unsigned long long int)n_M_per_Mbatch*(unsigned long long int)FTK_n_delta_v;
  l2_dSM___ = (double *) malloc(tab*sizeof(double));
  n2_dSM___ = (double *) malloc(tab*sizeof(double));
  f2_dSM___ = (double *) malloc(tab*sizeof(double));
  tab = (unsigned long long int)n_S_per_Sbatch*(unsigned long long int)n_M_per_Mbatch*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_delta_v;
  X_sub_dwSM____ = (double *) malloc(tab*sizeof(double));
  d_X_sub_dwSM____ = (double *) malloc(tab*sizeof(double));
  if (flag_compute_I_value){ I_value_use_dwSM____ = (double *) malloc(tab*sizeof(double));}
  tab = (unsigned long long int)n_w_max;
  gamma_z_ = (double *) malloc(tab*sizeof(double));
  for (nw=0;nw<n_w_max;nw++){ gamma_z_[nw] = 2*PI*(double)nw/(double)n_w_max;}
  /* %%%%%%%%%%%%%%%% */
  n_Mbatch = ceil((double)n_M/(double)n_M_per_Mbatch);
  if (verbose>1){ printf(" %% n_Mbatch %d\n",n_Mbatch);}
  index_M_in_Mbatch_ = (int *) malloc(n_M_per_Mbatch*sizeof(int));
  n_Sbatch = ceil((double)n_S/(double)n_S_per_Sbatch);
  if (verbose>1){ printf(" %% n_Sbatch %d\n",n_Sbatch);}
  index_S_in_Sbatch_ = (int *) malloc(n_S_per_Sbatch*sizeof(int));
  for (nMbatch=0;nMbatch<n_Mbatch;nMbatch++){
    n_M_sub = n_M_per_Mbatch; if (nMbatch==n_Mbatch-1){ n_M_sub = n_M - n_M_per_Mbatch*nMbatch;}
    for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ index_M_in_Mbatch_[nM_sub] = n_M_per_Mbatch*nMbatch + nM_sub;}
    if (verbose>0){ printf(" %% nMbatch %d/%d: index_M_in_Mbatch_ %d --> %d\n",nMbatch,n_Mbatch,index_M_in_Mbatch_[0],index_M_in_Mbatch_[n_M_sub-1]);}
    if (n_M_sub>0){
      n_M_sub_rup = rup(n_M_sub,8);
      n_M_sub_256 = n_M_sub_rup/8;
      /* %%%% */
      /* svd_VUXM_nMwl____ = permute(svd_VUXM_lwnM____(:,:,:,1+index_M_in_Mbatch_),[3,4,2,1]); */
      /* %%%% */
      if (flag_slow){
	local_tic(0,t_start_,d_start_);
	tab = (unsigned long long int)FTK_n_svd_l*(unsigned long long int)n_w_max*(unsigned long long int)n_M_sub*(unsigned long long int)pm_n_UX_rank;
	na = (unsigned long long int)FTK_n_svd_l*(unsigned long long int)n_w_max*(unsigned long long int)pm_n_UX_rank*(unsigned long long int)(n_M_per_Mbatch*nMbatch);
	for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
	  for (pm_nUX_rank=0;pm_nUX_rank<pm_n_UX_rank;pm_nUX_rank++){
	    for (nw=0;nw<n_w_max;nw++){
	      for (nl=0;nl<FTK_n_svd_l;nl++){
		/* tabA = */
		/* 	(unsigned long long int)pm_nUX_rank */
		/* 	+ (unsigned long long int)nM_sub*(unsigned long long int)pm_n_UX_rank */
		/* 	+ (unsigned long long int)nw*(unsigned long long int)n_M_sub*(unsigned long long int)pm_n_UX_rank */
		/* 	+ (unsigned long long int)nl*(unsigned long long int)n_w_max*(unsigned long long int)n_M_sub*(unsigned long long int)pm_n_UX_rank; */
		tabA = (unsigned long long int)pm_nUX_rank + ((unsigned long long int)nM_sub + ((unsigned long long int)nw + (unsigned long long int)nl*(unsigned long long int)n_w_max)*(unsigned long long int)n_M_sub)*(unsigned long long int)pm_n_UX_rank;
		tabB = (unsigned long long int)pm_nUX_rank + ((unsigned long long int)nM_sub + ((unsigned long long int)nl + (unsigned long long int)nw*(unsigned long long int)FTK_n_svd_l)*(unsigned long long int)n_M_sub)*(unsigned long long int)pm_n_UX_rank;
		svd_VUXM_nMwl____[tabA] = svd_VUXM_lwnM____[na];
		svd_VUXM_nMlw____[tabB] = svd_VUXM_lwnM____[na];
		na++;
		/* for (nl=0;nl<FTK_n_svd_l;nl++){ } */}
	      /* for (nw=0;nw<n_w_max;nw++){ } */}
	    /* for (pm_nUX_rank=0;pm_nUX_rank<pm_n_UX_rank;pm_nUX_rank++){ } */}
	  /* for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ } */}
	local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"svd_VUXM_nMwl____: ");
	elrt_sum_svd_VUXM_nMwl____ += elrt_[0];
	n_op_sum_svd_VUXM_nMwl____ += tab;
	/* if (flag_slow){ } */}
      /* %%%% */
      local_tic(0,t_start_,d_start_);
      tab = (unsigned long long int)FTK_n_svd_l*(unsigned long long int)n_w_max*(unsigned long long int)pm_n_UX_rank*(unsigned long long int)n_M_sub;
      tabA = (unsigned long long int)FTK_n_svd_l_rup*(unsigned long long int)n_w_max_rup*(unsigned long long int)pm_n_UX_rank_rup*(unsigned long long int)n_M_sub;
      memset(f_svd_VUXM_nMlw_real____,0,tabA*sizeof(float));
      memset(f_svd_VUXM_nMlw_imag____,0,tabA*sizeof(float));
      tabB = (unsigned long long int)FTK_n_svd_l_rup*(unsigned long long int)n_w_max_rup*(unsigned long long int)pm_n_UX_rank_rup*(unsigned long long int)(n_M_per_Mbatch*nMbatch);
      transpose_ps_block1(
      			  f_svd_VUXM_lwnM_real____ + tabB
      			  ,f_svd_VUXM_nMlw_real____
      			  ,FTK_n_svd_l_rup*n_w_max_rup
      			  ,pm_n_UX_rank_rup*n_M_sub_rup
      			  ,FTK_n_svd_l_rup*n_w_max_rup
      			  ,pm_n_UX_rank_rup*n_M_sub_rup
      			  ,16
      			  );
      transpose_ps_block1(
      			  f_svd_VUXM_lwnM_imag____ + tabB
      			  ,f_svd_VUXM_nMlw_imag____
      			  ,FTK_n_svd_l_rup*n_w_max_rup
      			  ,pm_n_UX_rank_rup*n_M_sub_rup
      			  ,FTK_n_svd_l_rup*n_w_max_rup
      			  ,pm_n_UX_rank_rup*n_M_sub_rup
      			  ,16
      			  );
      if (verbose>2 && flag_slow){
	array_sub_printf(  svd_VUXM_nMlw____ ,"double complex",pm_n_UX_rank    *n_M_sub    ,4,FTK_n_svd_l    *n_w_max    ,3," %%   svd_VUXM_nMlw____ : ");
	array_sub_printf(f_svd_VUXM_nMlw_real____,"float" ,pm_n_UX_rank_rup*n_M_sub_rup,4,FTK_n_svd_l_rup*n_w_max_rup,3," %% f_svd_VUXM_nMlw_real____: ");
	array_sub_printf(f_svd_VUXM_nMlw_imag____,"float" ,pm_n_UX_rank_rup*n_M_sub_rup,4,FTK_n_svd_l_rup*n_w_max_rup,3," %% f_svd_VUXM_nMlw_imag____: ");
	exit(0);
	/* if (verbose>2){ } */}  
      local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"f_svd_VUXM_nMlw____: ");
      elrt_sum_f_svd_VUXM_nMlw____ += elrt_[0];
      n_op_sum_f_svd_VUXM_nMlw____ += tab;
      /* %%%% */
      /* 
	 svd_SVUXM_SMwl____ = zeros(n_S,n_M_sub,n_w_max,FTK.n_svd_l);
	 for nl=0:FTK.n_svd_l-1;
	 for nw=0:n_w_max-1;
	 svd_SVUXM_SMwl____(:,:,1+nw,1+nl) = ctranspose(CTF_UX_S_k_q_nSw___(:,:,1+nw))*svd_VUXM_nMwl____(:,:,1+nw,1+nl);
	 end;%for nw=0:n_w_max-1;
	 end;%for nl=0:FTK.n_svd_l-1;
      */
      /* %%%% */
      if (flag_slow){
	local_tic(0,t_start_,d_start_);
	tab = (unsigned long long int)(FTK_n_svd_l*n_w_max)*(unsigned long long int)(n_S*n_M_sub)*(unsigned long long int)pm_n_UX_rank;
	for (nl=0;nl<FTK_n_svd_l;nl++){
	  for (nw=0;nw<n_w_max;nw++){
	    tabA = (unsigned long long int)nw*(unsigned long long int)pm_n_UX_rank*(unsigned long long int)n_S;
	    tabB = (unsigned long long int)(nw+nl*n_w_max)*(unsigned long long int)pm_n_UX_rank*(unsigned long long int)n_M_sub;
	    tabC = (unsigned long long int)(nw+nl*n_w_max)*(unsigned long long int)n_S*(unsigned long long int)n_M_sub;
	    cblas_zgemm(
			CblasColMajor
			,CblasConjTrans
			,CblasNoTrans
			,n_S
			,n_M_sub
			,pm_n_UX_rank
			,&cblas_alpha
			,CTF_UX_S_k_q_nSw___ + tabA
			,pm_n_UX_rank
			,svd_VUXM_nMwl____ + tabB
			,pm_n_UX_rank
			,&cblas_beta
			,svd_SVUXM_SMwl____ + tabC
			,n_S
			);
	    /* for (nw=0;nw<n_w_max;nw++){ } */}
	  /* for (nl=0;nl<FTK_n_svd_l;nl++){ } */}
	local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"svd_SVUXM_SMwl____: ");
	elrt_sum_svd_SVUXM_SMwl____ += elrt_[0];
	n_op_sum_svd_SVUXM_SMwl____ += tab;
	/* if (flag_slow){ } */}
      /* %%%% */
      if (flag_slow){
	local_tic(0,t_start_,d_start_);
	tab = (unsigned long long int)(FTK_n_svd_l*n_w_max)*(unsigned long long int)(n_S*n_M_sub)*(unsigned long long int)pm_n_UX_rank;
	for (nl=0;nl<FTK_n_svd_l;nl++){
	  for (nw=0;nw<n_w_max;nw++){
	    tabA = (unsigned long long int)nw*(unsigned long long int)pm_n_UX_rank*(unsigned long long int)n_S;
	    tabB = (unsigned long long int)(nl+nw*FTK_n_svd_l)*(unsigned long long int)pm_n_UX_rank*(unsigned long long int)n_M_sub;
	    tabC = (unsigned long long int)(nl+nw*FTK_n_svd_l)*(unsigned long long int)n_S*(unsigned long long int)n_M_sub;
	    cblas_zgemm(
			CblasColMajor
			,CblasConjTrans
			,CblasNoTrans
			,n_S
			,n_M_sub
			,pm_n_UX_rank
			,&cblas_alpha
			,CTF_UX_S_k_q_nSw___ + tabA
			,pm_n_UX_rank
			,svd_VUXM_nMlw____ + tabB
			,pm_n_UX_rank
			,&cblas_beta
			,svd_SVUXM_SMlw____ + tabC
			,n_S
			);
	    /* for (nw=0;nw<n_w_max;nw++){ } */}
	  /* for (nl=0;nl<FTK_n_svd_l;nl++){ } */}
	local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"svd_SVUXM_SMlw____: ");
	elrt_sum_svd_SVUXM_SMlw____ += elrt_[0];
	n_op_sum_svd_SVUXM_SMlw____ += tab;
	/* if (flag_slow){ } */}
      /* %%%% */
      local_tic(0,t_start_,d_start_);
      tab = (unsigned long long int)(FTK_n_svd_l*n_w_max)*(unsigned long long int)(n_S*n_M_sub)*(unsigned long long int)pm_n_UX_rank;
      tabA = (unsigned long long int)(FTK_n_svd_l_rup*n_w_max_rup)*(unsigned long long int)(n_S_rup*n_M_sub_rup);
      memset(f_svd_SVUXM_SMlw_real____,0,tabA*sizeof(float));
      memset(f_svd_SVUXM_SMlw_imag____,0,tabA*sizeof(float));
      for (nl=0;nl<FTK_n_svd_l;nl++){
	for (nw=0;nw<n_w_max;nw++){
	  tabA = (unsigned long long int)nw*(unsigned long long int)pm_n_UX_rank_rup*(unsigned long long int)n_S_rup;
	  tabB = (unsigned long long int)(nl+nw*FTK_n_svd_l_rup)*(unsigned long long int)pm_n_UX_rank_rup*(unsigned long long int)n_M_sub_rup;
	  tabC = (unsigned long long int)(nl+nw*FTK_n_svd_l_rup)*(unsigned long long int)n_S_rup*(unsigned long long int)n_M_sub_rup;
	  f_CR_ = f_svd_SVUXM_SMlw_real____ + tabC;
	  f_CI_ = f_svd_SVUXM_SMlw_imag____ + tabC;
	  hp_segregated_to_segregated_mult_immintrin_load1_fma(
							       n_S_rup
							       ,pm_n_UX_rank_rup
							       ,f_CTF_UX_S_k_q_nSw_real___ + tabA
							       ,f_CTF_UX_S_k_q_nSw_imag___ + tabA
							       ,n_M_sub_rup
							       ,f_svd_VUXM_nMlw_real____ + tabB
							       ,f_svd_VUXM_nMlw_imag____ + tabB
							       ,&f_CR_
							       ,&f_CI_
							       );
	  /* for (nw=0;nw<n_w_max;nw++){ } */}
	/* for (nl=0;nl<FTK_n_svd_l;nl++){ } */}
      local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"svd_f_SVUXM_SMlw____: ");
      elrt_sum_f_svd_SVUXM_SMlw____ += elrt_[0];
      n_op_sum_f_svd_SVUXM_SMlw____ += tab;
      if (verbose>2 && flag_slow){
	array_sub_printf(  svd_SVUXM_SMlw____ ,"double complex",n_S    *n_M_sub    ,4,FTK_n_svd_l    *n_w_max    ,3," %%   svd_SVUXM_SMlw____ : ");
	array_sub_printf(f_svd_SVUXM_SMlw_real____,"float" ,n_S_rup*n_M_sub_rup,4,FTK_n_svd_l_rup*n_w_max_rup,3," %% f_svd_SVUXM_SMlw_real____: ");
	array_sub_printf(f_svd_SVUXM_SMlw_imag____,"float" ,n_S_rup*n_M_sub_rup,4,FTK_n_svd_l_rup*n_w_max_rup,3," %% f_svd_SVUXM_SMlw_imag____: ");
	exit(0);
	/* if (verbose>2){ } */}
      /* %%%% */
      /* svd_SVUXM_lwSM____ = permute(ifft(permute(svd_SVUXM_SMwl____,[3,4,1,2]),[],1)*n_w_max,[2,1,3,4]); */
      /* %%%% */
      if (flag_slow){
	local_tic(0,t_start_,d_start_);
	tab = (unsigned long long int)n_S*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_svd_l;
	memset(svd_SVUXM_0in_wlSM____,0,tab*sizeof(fftw_complex));
	na=0; tabA=0;
	for (nl=0;nl<FTK_n_svd_l;nl++){
	  for (nw=0;nw<n_w_max;nw++){
	    for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
	      for (nS=0;nS<n_S;nS++){
		/* tabA = */
		/* 	(unsigned long long int)nw */
		/* 	+ (unsigned long long int)nl*(unsigned long long int)n_w_max */
		/* 	+ (unsigned long long int)nS*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_svd_l */
		/* 	+ (unsigned long long int)nM_sub*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_svd_l*(unsigned long long int)n_S; */
		tabA = (unsigned long long int)nw + ((unsigned long long int)nl + ((unsigned long long int)nS + (unsigned long long int)nM_sub*(unsigned long long int)n_S)*(unsigned long long int)FTK_n_svd_l)*(unsigned long long int)n_w_max;
		svd_SVUXM_0in_wlSM____[tabA] = (fftw_complex) (svd_SVUXM_SMwl____[na]); na++;
		/* for (nS=0;nS<n_S;nS++){ } */}
	      /* for (nM_sub=0;nM_sub<n_M_sub;nS++){ } */}
	    /* for (nw=0;nw<n_w_max;nw++){ } */}
	  /* for (nl=0;nl<FTK_n_svd_l;nl++){ } */}
	local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"svd_SVUXM_0in_wlSM____: ");
	elrt_sum_svd_SVUXM_0in_wlSM____ += elrt_[0];
	n_op_sum_svd_SVUXM_0in_wlSM____ += tab;
	local_tic(0,t_start_,d_start_);
	fftw_execute(fftw_plan_many_plan);
	local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"svd_SVUXM_out_wlSM____: ");
	elrt_sum_svd_SVUXM_out_wlSM____ += elrt_[0];
	n_op_sum_svd_SVUXM_out_wlSM____ += tab;
	local_tic(0,t_start_,d_start_);
	tab = (unsigned long long int)n_S*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_svd_l;
	memset(svd_SVUXM_lwSM____,0,tab*sizeof(double complex));
	na=0; tabA=0;
	for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
	  for (nS=0;nS<n_S;nS++){
	    for (nl=0;nl<FTK_n_svd_l;nl++){
	      for (nw=0;nw<n_w_max;nw++){
		/* tabA = */
		/* 	(unsigned long long int)nl */
		/* 	+ (unsigned long long int)nw*(unsigned long long int)FTK_n_svd_l */
		/* 	+ (unsigned long long int)nS*(unsigned long long int)FTK_n_svd_l*(unsigned long long int)n_w_max */
		/* 	+ (unsigned long long int)nM_sub*(unsigned long long int)FTK_n_svd_l*(unsigned long long int)n_w_max*(unsigned long long int)n_S; */
		tabA = (unsigned long long int)nl + ((unsigned long long int)nw + ((unsigned long long int)nS + (unsigned long long int)nM_sub*(unsigned long long int)n_S)*(unsigned long long int)n_w_max)*(unsigned long long int)FTK_n_svd_l;
		svd_SVUXM_lwSM____[tabA] = (double complex) (svd_SVUXM_out_wlSM____[na]); na++;
		/* for (nw=0;nw<n_w_max;nw++){ } */}
	      /* for (nl=0;nl<FTK_n_svd_l;nl++){ } */}
	    /* for (nS=0;nS<n_S;nS++){ } */}
	  /* for (nM_sub=0;nM_sub<n_M_sub;nS++){ } */}
	local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"svd_SVUXM_lwSM____: ");
	elrt_sum_svd_SVUXM_lwSM____ += elrt_[0];
	n_op_sum_svd_SVUXM_lwSM____ += tab;
	/* if (flag_slow){ } */}
      /* %%%% */
      local_tic(0,t_start_,d_start_);
      tab = (unsigned long long int)n_w_max*(unsigned long long int)n_S*(unsigned long long int)n_M_sub*(unsigned long long int)FTK_n_svd_l;
      tabA = (unsigned long long int)n_w_max_rup*(unsigned long long int)n_S_rup*(unsigned long long int)n_M_sub_rup*(unsigned long long int)FTK_n_svd_l_rup;
      memset(f_svd_SVUXM_0in_wSMl____,0,tabA*sizeof(fftwf_complex));
      f_C_ = (float *) f_svd_SVUXM_0in_wSMl____;
      for (nl=0;nl<FTK_n_svd_l;nl++){
	for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
	  for (nS=0;nS<n_S;nS++){
	    for (nw=0;nw<n_w_max;nw++){
	      tabA = (unsigned long long int)nS + ((unsigned long long int)nM_sub + ((unsigned long long int)nl + (unsigned long long int)nw*(unsigned long long int)FTK_n_svd_l_rup)*(unsigned long long int)n_M_sub_rup)*(unsigned long long int)n_S_rup;
	      tabB = (unsigned long long int)nw + ((unsigned long long int)nS + ((unsigned long long int)nM_sub + (unsigned long long int)nl*(unsigned long long int)n_M_sub_rup)*(unsigned long long int)n_S_rup)*(unsigned long long int)n_w_max_rup;
	      f_C_[2*tabB + 0] = f_svd_SVUXM_SMlw_real____[tabA]; 
	      f_C_[2*tabB + 1] = f_svd_SVUXM_SMlw_imag____[tabA]; 
	      na++;
	      /* for (nw=0;nw<n_w_max;nw++){ } */}
	    /* for (nS=0;nS<n_S;nS++){ } */}
	  /* for (nM_sub=0;nM_sub<n_M_sub;nS++){ } */}
	/* for (nl=0;nl<FTK_n_svd_l;nl++){ } */}
      local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"f_svd_SVUXM_0in_wSMl____: ");
      elrt_sum_f_svd_SVUXM_0in_wSMl____ += elrt_[0];
      n_op_sum_f_svd_SVUXM_0in_wSMl____ += tab;
      local_tic(0,t_start_,d_start_);
      fftwf_execute(fftwf_plan_many_plan);
      local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"f_svd_SVUXM_out_wSMl____: ");
      elrt_sum_f_svd_SVUXM_out_wSMl____ += elrt_[0];
      n_op_sum_f_svd_SVUXM_out_wSMl____ += tab;
      local_tic(0,t_start_,d_start_);
      tab = (unsigned long long int)FTK_n_svd_l*(unsigned long long int)n_w_max*(unsigned long long int)n_S*(unsigned long long int)n_M_sub;
      tabA = (unsigned long long int)FTK_n_svd_l_rup*(unsigned long long int)n_w_max_rup*(unsigned long long int)n_S_rup*(unsigned long long int)n_M_sub_rup;
      memset(f_svd_SVUXM_lwSM_real____,0,tabA*sizeof(float));
      memset(f_svd_SVUXM_lwSM_imag____,0,tabA*sizeof(float));
      f_C_ = (float *) f_svd_SVUXM_out_wSMl____;
      for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
	for (nS=0;nS<n_S;nS++){
	  for (nw=0;nw<n_w_max;nw++){
	    for (nl=0;nl<FTK_n_svd_l;nl++){
      	      tabA = (unsigned long long int)nw + ((unsigned long long int)nS + ((unsigned long long int)nM_sub + (unsigned long long int)nl*(unsigned long long int)n_M_sub_rup)*(unsigned long long int)n_S_rup)*(unsigned long long int)n_w_max_rup;
      	      tabB = (unsigned long long int)nl + ((unsigned long long int)nw + ((unsigned long long int)nS + (unsigned long long int)nM_sub*(unsigned long long int)n_S_rup)*(unsigned long long int)n_w_max_rup)*(unsigned long long int)FTK_n_svd_l_rup;
      	      f_svd_SVUXM_lwSM_real____[tabB] = f_C_[2*tabA+0];
      	      f_svd_SVUXM_lwSM_imag____[tabB] = f_C_[2*tabA+1];
	      /* for (nl=0;nl<FTK_n_svd_l;nl++){ } */}
	    /* for (nw=0;nw<n_w_max;nw++){ } */}
	  /* for (nS=0;nS<n_S;nS++){ } */}
	/* for (nM_sub=0;nM_sub<n_M_sub;nS++){ } */}
      local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"f_svd_SVUXM_lwSM____: ");
      elrt_sum_f_svd_SVUXM_lwSM____ += elrt_[0];
      n_op_sum_f_svd_SVUXM_lwSM____ += tab;
      if (verbose>2 && flag_slow){
	array_sub_printf(  svd_SVUXM_lwSM____ ,"double complex",FTK_n_svd_l    *n_w_max    ,3,n_S    *n_M_sub    ,4," %%   svd_SVUXM_lwSM____ : ");
	array_sub_printf(f_svd_SVUXM_lwSM_real____,"float" ,FTK_n_svd_l_rup*n_w_max_rup,3,n_S_rup*n_M_sub_rup,4," %% f_svd_SVUXM_lwSM_real____: ");
	array_sub_printf(f_svd_SVUXM_lwSM_imag____,"float" ,FTK_n_svd_l_rup*n_w_max_rup,3,n_S_rup*n_M_sub_rup,4," %% f_svd_SVUXM_lwSM_imag____: ");
	exit(0);
	/* if (verbose>2){ } */}      
      /* %%%% */
      /* display output */
      /* %%%% */
      if (flag_disp && flag_slow){
	printf(" %% nMbatch %d/%d\n",nMbatch,n_Mbatch);
	tab = minimum(8,(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_svd_l*(unsigned long long int)n_S*(unsigned long long int)n_M_sub);
	array_printf(svd_SVUXM_SMwl____,"double complex",1,tab," %% svd_SVUXM_SMwl____: ");
	array_printf(svd_SVUXM_lwSM____,"double complex",1,tab," %% svd_SVUXM_lwSM____: ");
	/* if (flag_disp){ } */}
      if (flag_dump && flag_slow){
	MDA_ndim=0;
	MDA_dim_[MDA_ndim++] = n_S;
	MDA_dim_[MDA_ndim++] = n_M_sub;
	MDA_dim_[MDA_ndim++] = n_w_max;
	MDA_dim_[MDA_ndim++] = FTK_n_svd_l;
	MDA_n_dim = MDA_ndim;
	sprintf(MDA_fname,"mex_svd_SVUXM_SMwl____.mat");
	MDA_write_c16(MDA_n_dim,MDA_dim_,svd_SVUXM_SMwl____,MDA_fname);
	MDA_ndim=0;
	MDA_dim_[MDA_ndim++] = FTK_n_svd_l;
	MDA_dim_[MDA_ndim++] = n_w_max;
	MDA_dim_[MDA_ndim++] = n_S;
	MDA_dim_[MDA_ndim++] = n_M_sub;
	MDA_n_dim = MDA_ndim;
	sprintf(MDA_fname,"mex_svd_SVUXM_lwSM____.mat");
	MDA_write_c16(MDA_n_dim,MDA_dim_,svd_SVUXM_lwSM____,MDA_fname);
	/* if (flag_dump){ } */}
      for (nSbatch=0;nSbatch<n_Sbatch;nSbatch++){
	n_S_sub = n_S_per_Sbatch; if (nSbatch==n_Sbatch-1){ n_S_sub = n_S - n_S_per_Sbatch*nSbatch;}
	for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){ index_S_in_Sbatch_[nS_sub] = n_S_per_Sbatch*nSbatch + nS_sub;}
	if (verbose>1){ printf(" %% nSbatch %d/%d: index_S_in_Sbatch_ %d --> %d\n",nSbatch,n_Sbatch,index_S_in_Sbatch_[0],index_S_in_Sbatch_[n_S_sub-1]);}
	if (n_S_sub>0){
	  n_S_sub_rup  = rup(n_S_sub,8);
	  n_S_sub_256 = n_S_sub_rup/8;
	  /* %%%% */
	  /* svd_SVUXM_lwsM____ = svd_SVUXM_lwSM____(:,:,1+index_S_in_Sbatch_,:); */
	  /* %%%% */
	  if (flag_slow){
	    local_tic(0,t_start_,d_start_);
	    tab = (unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_svd_l;
	    memset(svd_SVUXM_lwsM____,0,tab*sizeof(double complex));
	    for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
	      for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){
		nS = index_S_in_Sbatch_[nS_sub];
		tabA = (((unsigned long long int)nS_sub + (unsigned long long int)nM_sub*(unsigned long long int)n_S_sub)*(unsigned long long int)n_w_max)*(unsigned long long int)FTK_n_svd_l;
		tabB = (((unsigned long long int)nS     + (unsigned long long int)nM_sub*(unsigned long long int)n_S    )*(unsigned long long int)n_w_max)*(unsigned long long int)FTK_n_svd_l;
		tabC = (unsigned long long int)n_w_max*(unsigned long long int)FTK_n_svd_l;
		memcpy( svd_SVUXM_lwsM____ + tabA , svd_SVUXM_lwSM____ + tabB , tabC*sizeof(double complex));
		/* for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){ } */}
	      /* for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ } */}
	    local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"svd_SVUXM_lwsM____: ");
	    elrt_sum_svd_SVUXM_lwsM____ += elrt_[0];
	    n_op_sum_svd_SVUXM_lwsM____ += tab;
	    /* if (flag_slow){ } */}
	  /* %%%% */
	  local_tic(0,t_start_,d_start_);
	  tab = (unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_svd_l;
	  tabA = (unsigned long long int)n_S_sub_rup*(unsigned long long int)n_M_sub_rup*(unsigned long long int)n_w_max_rup*(unsigned long long int)FTK_n_svd_l_rup;
	  memset(f_svd_SVUXM_lwsM_real____,0,tabA*sizeof(float));
	  memset(f_svd_SVUXM_lwsM_imag____,0,tabA*sizeof(float));
	  for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
	    for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){
	      nS = index_S_in_Sbatch_[nS_sub];
	      tabA = (((unsigned long long int)nS_sub + (unsigned long long int)nM_sub*(unsigned long long int)n_S_sub_rup)*(unsigned long long int)n_w_max_rup)*(unsigned long long int)FTK_n_svd_l_rup;
	      tabB = (((unsigned long long int)nS     + (unsigned long long int)nM_sub*(unsigned long long int)n_S_rup    )*(unsigned long long int)n_w_max_rup)*(unsigned long long int)FTK_n_svd_l_rup;
	      tabC = (unsigned long long int)n_w_max_rup*(unsigned long long int)FTK_n_svd_l_rup;
	      memcpy( f_svd_SVUXM_lwsM_real____ + tabA , f_svd_SVUXM_lwSM_real____ + tabB , tabC*sizeof(float));
	      memcpy( f_svd_SVUXM_lwsM_imag____ + tabA , f_svd_SVUXM_lwSM_imag____ + tabB , tabC*sizeof(float));
	      /* for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){ } */}
	    /* for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ } */}
	  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"f_svd_SVUXM_lwsM____: ");
	  elrt_sum_f_svd_SVUXM_lwsM____ += elrt_[0];
	  n_op_sum_f_svd_SVUXM_lwsM____ += tab;
	  if (verbose>2 && flag_slow){
	    array_sub_printf(  svd_SVUXM_lwsM____ ,"double complex",FTK_n_svd_l    *n_w_max    ,3,n_S_sub    *n_M_sub    ,4," %%   svd_SVUXM_lwsM____ : ");
	    array_sub_printf(f_svd_SVUXM_lwsM_real____,"float" ,FTK_n_svd_l_rup*n_w_max_rup,3,n_S_sub_rup*n_M_sub_rup,4," %% f_svd_SVUXM_lwsM_real____: ");
	    array_sub_printf(f_svd_SVUXM_lwsM_imag____,"float" ,FTK_n_svd_l_rup*n_w_max_rup,3,n_S_sub_rup*n_M_sub_rup,4," %% f_svd_SVUXM_lwsM_imag____: ");
	    exit(0);
	    /* if (verbose>2){ } */}      
	  /* %%%% */	  
	  if (flag_dump && flag_slow){
	    MDA_ndim=0;
	    MDA_dim_[MDA_ndim++] = FTK_n_svd_l;
	    MDA_dim_[MDA_ndim++] = n_w_max;
	    MDA_dim_[MDA_ndim++] = n_S_sub;
	    MDA_dim_[MDA_ndim++] = n_M_sub;
	    MDA_n_dim = MDA_ndim;
	    sprintf(MDA_fname,"mex_svd_SVUXM_lwsM____.mat");
	    MDA_write_c16(MDA_n_dim,MDA_dim_,svd_SVUXM_lwsM____,MDA_fname);
	    /* if (flag_dump){ } */}
	  /* %%%% */
	  /* svd_USESVUXM_dwSM____ = real(reshape(FTK.svd_U_d_expiw_s__*reshape(svd_SVUXM_lwsM____,[FTK.n_svd_l,n_w_max*n_S_sub*n_M_sub]),[FTK.n_delta_v,n_w_max,n_S_sub,n_M_sub])); */
	  /* %%%% */
	  if (flag_slow){
	    local_tic(0,t_start_,d_start_);
	    tab = (unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_delta_v;
	    memset(svd_USESVUXM_dwSM____,0,tab*sizeof(double complex));
	    tabC = (unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max;
	    tab = (unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_delta_v*(unsigned long long int)FTK_n_svd_l;
	    cblas_zgemm(
			CblasColMajor
			,CblasNoTrans
			,CblasNoTrans
			,FTK_n_delta_v
			,tabC
			,FTK_n_svd_l
			,&cblas_alpha
			,FTK_svd_U_d_expiw_s__
			,FTK_n_delta_v
			,svd_SVUXM_lwsM____
			,FTK_n_svd_l
			,&cblas_beta
			,svd_USESVUXM_dwSM____
			,FTK_n_delta_v
			);
	    local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"svd_USESVUXM_dwSM____: ");
	    elrt_sum_svd_USESVUXM_dwSM____ += elrt_[0];
	    n_op_sum_svd_USESVUXM_dwSM____ += tab;
	    /* if (flag_slow){ } */}
	  /* %%%% */
	  local_tic(0,t_start_,d_start_);
	  tab = (unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_delta_v*(unsigned long long int)FTK_n_svd_l;
	  tabA = (unsigned long long int)n_S_sub_rup*(unsigned long long int)n_M_sub_rup*(unsigned long long int)n_w_max_rup*(unsigned long long int)FTK_n_delta_v_rup;
	  memset(f_svd_USESVUXM_dwSM_real____,0,tabA*sizeof(float));
	  memset(f_svd_USESVUXM_dwSM_imag____,0,tabA*sizeof(float));
	  tabC = (unsigned long long int)n_S_sub_rup*(unsigned long long int)n_M_sub_rup*(unsigned long long int)n_w_max_rup;
	  nhpr_segregated_to_segregated_mult_immintrin_load1_fma(
							       FTK_n_delta_v_rup
							       ,FTK_n_svd_l_rup
							       ,f_FTK_svd_U_d_expiw_s_tran_real__
							       ,f_FTK_svd_U_d_expiw_s_tran_imag__
							       ,tabC
							       ,f_svd_SVUXM_lwsM_real____
							       ,f_svd_SVUXM_lwsM_imag____
							       ,&f_svd_USESVUXM_dwSM_real____
							       ,&f_svd_USESVUXM_dwSM_imag____
							       );
	  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"f_svd_USESVUXM_dwSM____: ");
	  elrt_sum_f_svd_USESVUXM_dwSM____ += elrt_[0];
	  n_op_sum_f_svd_USESVUXM_dwSM____ += tab;
	  if (verbose>2 && flag_slow){
	    array_sub_printf(  svd_USESVUXM_dwSM____ ,"double complex",FTK_n_delta_v    *n_w_max    ,3,n_S_sub    *n_M_sub    ,4," %%   svd_USESVUXM_dwSM____ : ");
	    array_sub_printf(f_svd_USESVUXM_dwSM_real____,"float" ,FTK_n_delta_v_rup*n_w_max_rup,3,n_S_sub_rup*n_M_sub_rup,4," %% f_svd_USESVUXM_dwSM_real____: ");
	    array_sub_printf(f_svd_USESVUXM_dwSM_imag____,"float" ,FTK_n_delta_v_rup*n_w_max_rup,3,n_S_sub_rup*n_M_sub_rup,4," %% f_svd_USESVUXM_dwSM_imag____: ");
	    exit(0);
	    /* if (verbose>2){ } */}      
	  /* %%%% */
	  if (flag_dump && flag_slow){
	    MDA_ndim=0;
	    MDA_dim_[MDA_ndim++] = FTK_n_delta_v;
	    MDA_dim_[MDA_ndim++] = n_w_max;
	    MDA_dim_[MDA_ndim++] = n_S_sub;
	    MDA_dim_[MDA_ndim++] = n_M_sub;
	    MDA_n_dim = MDA_ndim;
	    sprintf(MDA_fname,"mex_svd_USESVUXM_dwSM____.mat");
	    MDA_write_c16(MDA_n_dim,MDA_dim_,svd_USESVUXM_dwSM____,MDA_fname);
	    /* if (flag_dump){ } */}
	  /* %%%% */
	  /* 
	     l2_dSM___ = permute(reshape(reshape(sqrt(CTF_UX_S_l2_use_(1+index_S_in_Sbatch_)),[n_S_sub,1])*reshape(sqrt(UX_M_l2_dM_use__(:,1+index_M_in_Mbatch_)),[1,FTK.n_delta_v*n_M_sub]),[n_S_sub,FTK.n_delta_v,n_M_sub]),[2,1,3]);
	     n2_dSM___ = 1./max(1e-14,l2_dSM___);
	     f2_dSM___ = permute(reshape(reshape(sqrt(CTF_UX_S_l2_use_(1+index_S_in_Sbatch_)),[n_S_sub,1])*reshape(1./max(1e-14,sqrt(UX_M_l2_dM_use__(:,1+index_M_in_Mbatch_))),[1,FTK.n_delta_v*n_M_sub]),[n_S_sub,FTK.n_delta_v,n_M_sub]),[2,1,3]);
	     ss_S_ = reshape(CTF_UX_S_l2_use_(1+index_S_in_Sbatch_),[n_S_sub,1]);
	  */
	  /* %%%% */
	  local_tic(0,t_start_,d_start_);
	  tab = (unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub*(unsigned long long int)FTK_n_delta_v;
	  memset(l2_dSM___,0,tab*sizeof(double));
	  memset(n2_dSM___,0,tab*sizeof(double));
	  memset(f2_dSM___,0,tab*sizeof(double));
	  for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){
	    nS = index_S_in_Sbatch_[nS_sub];
	    for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
	      nM = index_M_in_Mbatch_[nM_sub];
	      for (ndelta_v=0;ndelta_v<FTK_n_delta_v;ndelta_v++){
		tabA = (unsigned long long int)ndelta_v + (unsigned long long int)nM*(unsigned long long int)FTK_n_delta_v;
		tabB = (unsigned long long int)ndelta_v + ((unsigned long long int)nS_sub + (unsigned long long int)nM_sub*(unsigned long long int)n_S_sub)*(unsigned long long int)FTK_n_delta_v;
		l2_dSM___[tabB] = sqrt(CTF_UX_S_l2_[nS]) * sqrt(UX_M_l2_dM__[tabA]);
		n2_dSM___[tabB] = 1.0 / maximum(1.0e-14,l2_dSM___[tabB]);
		f2_dSM___[tabB] = sqrt(CTF_UX_S_l2_[nS]) / maximum(1.0e-14,sqrt(UX_M_l2_dM__[tabA]));
		/* for (ndelta_v=0;ndelta_v<FTK_n_delta_v;ndelta_v++){ } */}
	      /* for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ } */}
	    /* for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){ } */}
	  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"l2_dSM___: ");
	  elrt_sum_l2_dSM___ += elrt_[0];
	  n_op_sum_l2_dSM___ += tab;
	  /* %%%% */
	  if (flag_dump){
	    MDA_ndim=0;
	    MDA_dim_[MDA_ndim++] = FTK_n_delta_v;
	    MDA_dim_[MDA_ndim++] = n_S_sub;
	    MDA_dim_[MDA_ndim++] = n_M_sub;
	    MDA_n_dim = MDA_ndim;
	    sprintf(MDA_fname,"mex_l2_dSM___.mat");
	    MDA_write_r8(MDA_n_dim,MDA_dim_,l2_dSM___,MDA_fname);
	    sprintf(MDA_fname,"mex_n2_dSM___.mat");
	    MDA_write_r8(MDA_n_dim,MDA_dim_,n2_dSM___,MDA_fname);
	    sprintf(MDA_fname,"mex_f2_dSM___.mat");
	    MDA_write_r8(MDA_n_dim,MDA_dim_,f2_dSM___,MDA_fname);
	    /* if (flag_dump){ } */}
	  /* %%%% */
	  /* X_sub_dwSM____ = repmat(reshape(n2_dSM___,[FTK.n_delta_v,1,n_S_sub,n_M_sub]),[1,n_w_max,1,1]) .* svd_USESVUXM_dwSM____; %<-- correlation. ; */
	  /* %%%% */
	  if (flag_slow){
	    local_tic(0,t_start_,d_start_);
	    tab = (unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_delta_v;
	    memset(X_sub_dwSM____,0,tab*sizeof(double));
	    na=0;
	    for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
	      for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){
		for (nw=0;nw<n_w_max;nw++){
		  for (ndelta_v=0;ndelta_v<FTK_n_delta_v;ndelta_v++){
		    tabA = (unsigned long long int)ndelta_v + ((unsigned long long int)nS_sub + (unsigned long long int)nM_sub*(unsigned long long int)n_S_sub)*FTK_n_delta_v;
		    X_sub_dwSM____[na] = n2_dSM___[tabA] * creal(svd_USESVUXM_dwSM____[na]) ;
		    na += 1;
		    /* for (ndelta_v=0;ndelta_v<FTK_n_delta_v;ndelta_v++){ } */}
		  /* for (nw=0;nw<n_w_max;nw++){ } */}
		/* for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){ } */}
	      /* for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ } */}
	    local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"X_sub_dwSM____: ");
	    elrt_sum_X_sub_dwSM____ += elrt_[0];
	    n_op_sum_X_sub_dwSM____ += tab;
	    /* if (flag_slow){ } */}
	  /* %%%% */
	  local_tic(0,t_start_,d_start_);
	  tab = (unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_delta_v;
	  memset(d_X_sub_dwSM____,0,tab*sizeof(double));
	  na=0;
	  for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
	    for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){
	      for (nw=0;nw<n_w_max;nw++){
		for (ndelta_v=0;ndelta_v<FTK_n_delta_v;ndelta_v++){
		  tabA = (unsigned long long int)ndelta_v + ((unsigned long long int)nS_sub + (unsigned long long int)nM_sub*(unsigned long long int)n_S_sub)*FTK_n_delta_v;
		  tabB = (unsigned long long int)ndelta_v + ((unsigned long long int)nw + ((unsigned long long int)nS_sub + (unsigned long long int)nM_sub*(unsigned long long int)n_S_sub_rup)*(unsigned long long int)n_w_max_rup)*(unsigned long long int)FTK_n_delta_v_rup;
		  d_X_sub_dwSM____[na] = n2_dSM___[tabA] * f_svd_USESVUXM_dwSM_real____[tabB] ;
		  na += 1;
		  /* for (ndelta_v=0;ndelta_v<FTK_n_delta_v;ndelta_v++){ } */}
		/* for (nw=0;nw<n_w_max;nw++){ } */}
	      /* for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){ } */}
	    /* for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ } */}
	  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"d_X_sub_dwSM____: ");
	  elrt_sum_d_X_sub_dwSM____ += elrt_[0];
	  n_op_sum_d_X_sub_dwSM____ += tab;
	  if (verbose>2 && flag_slow){
	    array_sub_printf(  X_sub_dwSM____ ,"double complex",FTK_n_delta_v*n_w_max,3,n_S_sub*n_M_sub,4," %%   X_sub_dwSM____ : ");
	    array_sub_printf(d_X_sub_dwSM____ ,"double complex",FTK_n_delta_v*n_w_max,3,n_S_sub*n_M_sub,4," %% d_X_sub_dwSM____ : ");
	    exit(0);
	    /* if (verbose>2){ } */}
	  /* %%%% */
	  tab = (unsigned long long int)FTK_n_delta_v*(unsigned long long int)n_w_max*(unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub;
	  memcpy( X_sub_dwSM____ , d_X_sub_dwSM____ , tab*sizeof(float));
	  /* %%%% */
	  if (flag_dump){
	    MDA_ndim=0;
	    MDA_dim_[MDA_ndim++] = FTK_n_delta_v;
	    MDA_dim_[MDA_ndim++] = n_w_max;
	    MDA_dim_[MDA_ndim++] = n_S_sub;
	    MDA_dim_[MDA_ndim++] = n_M_sub;
	    MDA_n_dim = MDA_ndim;
	    sprintf(MDA_fname,"mex_X_sub_dwSM____.mat");
	    MDA_write_r8(MDA_n_dim,MDA_dim_,X_sub_dwSM____,MDA_fname);
	    /* if (flag_dump){ } */}
	  /* %%%% */
	  /* 
	     I_value_sub_dwSM____ = repmat(reshape(f2_dSM___,[FTK.n_delta_v,1,n_S_sub,n_M_sub]),[1,n_w_max,1,1]) .* X_sub_dwSM____; %<-- I_value. ;
	     I_value_use_dwSM____ = max(0,real(I_value_sub_dwSM____));
	  */
	  /* %%%% */
	  if (flag_compute_I_value){
	    local_tic(0,t_start_,d_start_);
	    tab = (unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_delta_v;
	    memset(I_value_use_dwSM____,0,tab*sizeof(double));
	    na=0;
	    for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
	      for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){
		for (nw=0;nw<n_w_max;nw++){
		  for (ndelta_v=0;ndelta_v<FTK_n_delta_v;ndelta_v++){
		    tabA = (unsigned long long int)ndelta_v + ((unsigned long long int)nS_sub + (unsigned long long int)nM_sub*(unsigned long long int)n_S_sub)*FTK_n_delta_v;
		    I_value_use_dwSM____[na] = maximum( 0.0 , f2_dSM___[tabA] * X_sub_dwSM____[na] ) ;
		    na += 1;
		    /* for (ndelta_v=0;ndelta_v<FTK_n_delta_v;ndelta_v++){ } */}
		  /* for (nw=0;nw<n_w_max;nw++){ } */}
	      /* for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){ } */}
	    /* for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ } */}
	    local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"I_value_use_dwSM____: ");
	    elrt_sum_I_value_use_dwSM____ += elrt_[0];
	    n_op_sum_I_value_use_dwSM____ += tab;
	    /* %%%% */
	    if (flag_dump){
	      MDA_ndim=0;
	      MDA_dim_[MDA_ndim++] = FTK_n_delta_v;
	      MDA_dim_[MDA_ndim++] = n_w_max;
	      MDA_dim_[MDA_ndim++] = n_S_sub;
	      MDA_dim_[MDA_ndim++] = n_M_sub;
	      MDA_n_dim = MDA_ndim;
	      sprintf(MDA_fname,"mex_I_value_use_dwSM____.mat");
	      MDA_write_r8(MDA_n_dim,MDA_dim_,I_value_use_dwSM____,MDA_fname);
	      /* if (flag_dump){ } */}
	    /* if (flag_compute_I_value){ } */}
	  /* %%%% */
	  /*
	    if (flag_optimize_over_gamma_z == 0);
	    tmp_t = tic(); nop=0;
	    [tmp_X_wSM___,tmp_delta_ij___] = max(reshape(X_sub_dwSM____,[FTK.n_delta_v,n_w_max*n_S_sub*n_M_sub]),[],1); %<-- maximize correlation. ;
	    assert(min(tmp_delta_ij___)>=1); assert(max(tmp_delta_ij___)<=FTK.n_delta_v);
	    tmp_X_wSM___ = reshape(tmp_X_wSM___,[n_w_max,n_S_sub,n_M_sub]);
	    tmp_delta_ij___ = reshape(tmp_delta_ij___,[n_w_max,n_S_sub,n_M_sub]);
	    tmp_delta_x___ = FTK.delta_x_(tmp_delta_ij___);
	    tmp_delta_y___ = FTK.delta_y_(tmp_delta_ij___);
	    tmp_gamma_z___ = 2*pi*(0:n_w_max-1)/n_w_max;
	    X_wSM___(:,1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = tmp_X_wSM___;
	    delta_x_wSM___(:,1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = reshape(tmp_delta_x___,[n_w_max,n_S_sub,n_M_sub]);
	    delta_y_wSM___(:,1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = reshape(tmp_delta_y___,[n_w_max,n_S_sub,n_M_sub]);
	    gamma_z_wSM___(:,1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = repmat(tmp_gamma_z___(:),[1,n_S_sub,n_M_sub]);
	    if (flag_compute_I_value);
	    tmp_I_value_use_dwSM__ = reshape(I_value_use_dwSM____,[FTK.n_delta_v,n_w_max*n_S_sub*n_M_sub]);
	    tmp_I_value_use_wSM_ = zeros(n_w_max*n_S_sub*n_M_sub,1);
	    tmp_t2=tic();
	    for nl=0:n_w_max*n_S_sub*n_M_sub-1;
	    tmp_I_value_use_wSM_(1+nl) = tmp_I_value_use_dwSM__(tmp_delta_ij___(1+nl),1+nl);
	    end;%for nl=0:n_w_max*n_S_sub*n_M_sub-1;
	    I_value_wSM___(:,1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = reshape(tmp_I_value_use_wSM_,[n_w_max,n_S_sub,n_M_sub]);
	    tmp_t2 = toc(tmp_t2); if (verbose>1); disp(sprintf(' %% I_value_wSM___ %0.6fs',tmp_t2)); end;
	    end;%if (flag_compute_I_value);
	    tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_wSM___: %0.6f',tmp_t)); end;
	    parameter = parameter_timing_update(parameter,'ampmh_X_wSM_mex___8: X_wSM___',tmp_t,1,nop);
	    end;%if (flag_optimize_over_gamma_z == 0);
	  */
	  /* %%%% */
	  if (flag_optimize_over_gamma_z==0){
	    local_tic(0,t_start_,d_start_);
	    tab = (unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_delta_v;
	    na=0;
	    for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
	      nM = index_M_in_Mbatch_[nM_sub];
	      for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){
		nS = index_S_in_Sbatch_[nS_sub];
		for (nw=0;nw<n_w_max;nw++){
		  tabA = (unsigned long long int)nw + ((unsigned long long int)nS_sub + (unsigned long long int)nM_sub*(unsigned long long int)n_S_sub)*(unsigned long long int)n_w_max;
		  ndelta_v = dmax_index(FTK_n_delta_v,X_sub_dwSM____ + tabA*(unsigned long long int)FTK_n_delta_v);
		  tabB = (unsigned long long int)nw + ((unsigned long long int)nS     + (unsigned long long int)nM    *(unsigned long long int)n_S    )*(unsigned long long int)n_w_max;
		  tabC = ndelta_v + tabA*(unsigned long long int)FTK_n_delta_v;
		  X_wSM___[tabB] = X_sub_dwSM____[tabC];
		  delta_x_wSM___[tabB] = FTK_delta_x_[ndelta_v];
		  delta_y_wSM___[tabB] = FTK_delta_y_[ndelta_v];
		  gamma_z_wSM___[tabB] = gamma_z_[nw];
		  if (flag_compute_I_value){ I_value_wSM___[tabB] = I_value_use_dwSM____[tabC];}
		  na += FTK_n_delta_v;
		  /* for (nw=0;nw<n_w_max;nw++){ } */}
		/* for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){ } */}
	      /* for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ } */}
	    local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"X_wSM___: ");
	    elrt_sum_X_wSM___ += elrt_[0];
	    n_op_sum_X_wSM___ += tab;
	    /* if (flag_optimize_over_gamma_z==0){ } */}
	  /* %%%% */
	  /*
	    if (flag_optimize_over_gamma_z == 1);
	    tmp_t = tic(); nop=0;
	    [tmp_X_SM__,tmp_dw_ij__] = max(reshape(X_sub_dwSM____,[FTK.n_delta_v*n_w_max,n_S_sub*n_M_sub]),[],1); %<-- maximize correlation. ;
	    [tmp_delta_ij__,tmp_gamma_ij__] = ind2sub([FTK.n_delta_v,n_w_max],tmp_dw_ij__);
	    assert(min(tmp_delta_ij__)>=1); assert(max(tmp_delta_ij__)<=FTK.n_delta_v);
	    assert(min(tmp_gamma_ij__)>=1); assert(max(tmp_gamma_ij__)<=n_w_max);
	    tmp_X_SM__ = reshape(tmp_X_SM__,[n_S_sub,n_M_sub]);
	    tmp_delta_ij__ = reshape(tmp_delta_ij__,[n_S_sub,n_M_sub]);
	    tmp_gamma_ij__ = reshape(tmp_gamma_ij__,[n_S_sub,n_M_sub]);
	    tmp_delta_x__ = FTK.delta_x_(tmp_delta_ij__);
	    tmp_delta_y__ = FTK.delta_y_(tmp_delta_ij__);
	    tmp_gamma_z__ = 2*pi*(tmp_gamma_ij__-1)/n_w_max;
	    X_SM__(1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = tmp_X_SM__;
	    delta_x_SM__(1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = reshape(tmp_delta_x__,[n_S_sub,n_M_sub]);
	    delta_y_SM__(1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = reshape(tmp_delta_y__,[n_S_sub,n_M_sub]);
	    gamma_z_SM__(1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = reshape(tmp_gamma_z__,[n_S_sub,n_M_sub]);
	    if (flag_compute_I_value);
	    tmp_I_value_use_dwSM___ = reshape(I_value_use_dwSM____,[FTK.n_delta_v,n_w_max,n_S_sub*n_M_sub]);
	    tmp_I_value_use_SM_ = zeros(n_S_sub*n_M_sub,1);
	    tmp_t2=tic();
	    for nl=0:n_S_sub*n_M_sub-1;
	    tmp_I_value_use_SM_(1+nl) = tmp_I_value_use_dwSM___(tmp_delta_ij__(1+nl),tmp_gamma_ij__(1+nl),1+nl);
	    end;%for nl=0:n_S_sub*n_M_sub-1;
	    I_value_SM__(1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = reshape(tmp_I_value_use_SM_,[n_S_sub,n_M_sub]);
	    tmp_t2 = toc(tmp_t2); if (verbose>1); disp(sprintf(' %% I_value_SM__ %0.6fs',tmp_t2)); end;
	    end;%if (flag_compute_I_value);
	    tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_wSM___: %0.6f',tmp_t)); end;
	    parameter = parameter_timing_update(parameter,'ampmh_X_wSM_mex___8: X_wSM___',tmp_t,1,nop);
	    end;%if (flag_optimize_over_gamma_z == 1);	    
	  */
	  /* %%%% */
	  if (flag_optimize_over_gamma_z==1){
	    local_tic(0,t_start_,d_start_);
	    tab = (unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_delta_v;
	    na=0;
	    for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
	      nM = index_M_in_Mbatch_[nM_sub];
	      for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){
		nS = index_S_in_Sbatch_[nS_sub];
		tabA = (unsigned long long int)nS_sub + (unsigned long long int)nM_sub*(unsigned long long int)n_S_sub;
		ndw = dmax_index((unsigned long long int)FTK_n_delta_v*(unsigned long long int)n_w_max,X_sub_dwSM____ + tabA*(unsigned long long int)FTK_n_delta_v*(unsigned long long int)n_w_max);
		nw = ndw / maximum(1,FTK_n_delta_v);
		ndelta_v = ndw % FTK_n_delta_v;
		tabB = (unsigned long long int)nS     + (unsigned long long int)nM    *(unsigned long long int)n_S    ;
		tabC = ndw + tabA*(unsigned long long int)FTK_n_delta_v*(unsigned long long int)n_w_max;
		X_wSM___[tabB] = X_sub_dwSM____[tabC];
		delta_x_wSM___[tabB] = FTK_delta_x_[ndelta_v];
		delta_y_wSM___[tabB] = FTK_delta_y_[ndelta_v];
		gamma_z_wSM___[tabB] = gamma_z_[nw];
		if (flag_compute_I_value){ I_value_wSM___[tabB] = I_value_use_dwSM____[tabC];}
		na += (unsigned long long int)FTK_n_delta_v*(unsigned long long int)n_w_max;
		/* for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){ } */}
	      /* for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ } */}
	    local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"X_wSM___: ");
	    elrt_sum_X_wSM___ += elrt_[0];
	    n_op_sum_X_wSM___ += tab;
	    /* if (flag_optimize_over_gamma_z==1){ } */}
	  /* %%%% */
	  if (flag_disp){
	    printf(" %% nMbatch %d/%d nSbatch %d/%d\n",nMbatch,n_Mbatch,nSbatch,n_Sbatch);
	    if (flag_optimize_over_gamma_z==0){ tab = minimum(8,(unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max);}
	    if (flag_optimize_over_gamma_z==1){ tab = minimum(8,(unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub);}
	    array_printf(X_wSM___,"double",1,tab," %% X_wSM___: ");
	    array_printf(delta_x_wSM___,"double",1,tab," %% delta_x_wSM___: ");
	    array_printf(delta_y_wSM___,"double",1,tab," %% delta_y_wSM___: ");
	    array_printf(gamma_z_wSM___,"double",1,tab," %% gamma_z_wSM___: ");
	    if (flag_compute_I_value){ array_printf(I_value_wSM___,"double",1,tab," %% I_value_wSM___: ");}
	    /* if (flag_disp){ } */}
	  if (flag_dump){
	    MDA_ndim=0;
	    if (flag_optimize_over_gamma_z==0){ MDA_dim_[MDA_ndim++] = n_w_max;}
	    MDA_dim_[MDA_ndim++] = n_S;
	    MDA_dim_[MDA_ndim++] = n_M;
	    MDA_n_dim = MDA_ndim;
	    sprintf(MDA_fname,"mex_X_wSM___.mat");MDA_write_r8(MDA_n_dim,MDA_dim_,X_wSM___,MDA_fname);
	    sprintf(MDA_fname,"mex_delta_x_wSM___.mat");MDA_write_r8(MDA_n_dim,MDA_dim_,delta_x_wSM___,MDA_fname);
	    sprintf(MDA_fname,"mex_delta_y_wSM___.mat");MDA_write_r8(MDA_n_dim,MDA_dim_,delta_y_wSM___,MDA_fname);
	    sprintf(MDA_fname,"mex_gamma_z_wSM___.mat");MDA_write_r8(MDA_n_dim,MDA_dim_,gamma_z_wSM___,MDA_fname);
	    if (flag_compute_I_value){ sprintf(MDA_fname,"mex_I_value_wSM___.mat");MDA_write_r8(MDA_n_dim,MDA_dim_,I_value_wSM___,MDA_fname);}
	    /* if (flag_dump){ } */}
	  /* if (n_S_sub>0){ } */}
	/* for (nSbatch=0;nSbatch<n_Sbatch;nSbatch++){ } */}
      /* if (n_M_sub>0){ } */}
    /* for (nMbatch=0;nMbatch<n_Mbatch;nMbatch++){ } */}
  /* %%%%%%%%%%%%%%%% */
  if (verbose){
    printf(" %% elrt_sum_CTF_UX_S_k_q_nSw___     : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_CTF_UX_S_k_q_nSw___,n_op_sum_CTF_UX_S_k_q_nSw___/elrt_sum_CTF_UX_S_k_q_nSw___/1e6,n_op_sum_CTF_UX_S_k_q_nSw___/elrt_sum_CTF_UX_S_k_q_nSw___/1e9);
    printf(" %% elrt_sum_f_CTF_UX_S_k_q_nSw___   : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_f_CTF_UX_S_k_q_nSw___,n_op_sum_f_CTF_UX_S_k_q_nSw___/elrt_sum_f_CTF_UX_S_k_q_nSw___/1e6,n_op_sum_f_CTF_UX_S_k_q_nSw___/elrt_sum_f_CTF_UX_S_k_q_nSw___/1e9);
    printf(" %% elrt_sum_svd_VUXM_nMwl____       : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_svd_VUXM_nMwl____,n_op_sum_svd_VUXM_nMwl____/elrt_sum_svd_VUXM_nMwl____/1e6,n_op_sum_svd_VUXM_nMwl____/elrt_sum_svd_VUXM_nMwl____/1e9);
    printf(" %% elrt_sum_f_svd_VUXM_nMlw____     : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_f_svd_VUXM_nMlw____,n_op_sum_f_svd_VUXM_nMlw____/elrt_sum_f_svd_VUXM_nMlw____/1e6,n_op_sum_f_svd_VUXM_nMlw____/elrt_sum_f_svd_VUXM_nMlw____/1e9);
    printf(" %% elrt_sum_svd_SVUXM_SMwl____      : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_svd_SVUXM_SMwl____,n_op_sum_svd_SVUXM_SMwl____/elrt_sum_svd_SVUXM_SMwl____/1e6,n_op_sum_svd_SVUXM_SMwl____/elrt_sum_svd_SVUXM_SMwl____/1e9);
    printf(" %% elrt_sum_svd_SVUXM_SMlw____      : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_svd_SVUXM_SMlw____,n_op_sum_svd_SVUXM_SMlw____/elrt_sum_svd_SVUXM_SMlw____/1e6,n_op_sum_svd_SVUXM_SMlw____/elrt_sum_svd_SVUXM_SMlw____/1e9);
    printf(" %% elrt_sum_f_svd_SVUXM_SMlw____    : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_f_svd_SVUXM_SMlw____,n_op_sum_f_svd_SVUXM_SMlw____/elrt_sum_f_svd_SVUXM_SMlw____/1e6,n_op_sum_f_svd_SVUXM_SMlw____/elrt_sum_f_svd_SVUXM_SMlw____/1e9);
    printf(" %% elrt_sum_svd_SVUXM_0in_wlSM____  : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_svd_SVUXM_0in_wlSM____,n_op_sum_svd_SVUXM_0in_wlSM____/elrt_sum_svd_SVUXM_0in_wlSM____/1e6,n_op_sum_svd_SVUXM_0in_wlSM____/elrt_sum_svd_SVUXM_0in_wlSM____/1e9);
    printf(" %% elrt_sum_f_svd_SVUXM_0in_wSMl____: %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_f_svd_SVUXM_0in_wSMl____,n_op_sum_f_svd_SVUXM_0in_wSMl____/elrt_sum_f_svd_SVUXM_0in_wSMl____/1e6,n_op_sum_f_svd_SVUXM_0in_wSMl____/elrt_sum_f_svd_SVUXM_0in_wSMl____/1e9);
    printf(" %% elrt_sum_svd_SVUXM_out_wlSM____  : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_svd_SVUXM_out_wlSM____,n_op_sum_svd_SVUXM_out_wlSM____/elrt_sum_svd_SVUXM_out_wlSM____/1e6,n_op_sum_svd_SVUXM_out_wlSM____/elrt_sum_svd_SVUXM_out_wlSM____/1e9);
    printf(" %% elrt_sum_f_svd_SVUXM_out_wSMl____: %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_f_svd_SVUXM_out_wSMl____,n_op_sum_f_svd_SVUXM_out_wSMl____/elrt_sum_f_svd_SVUXM_out_wSMl____/1e6,n_op_sum_f_svd_SVUXM_out_wSMl____/elrt_sum_f_svd_SVUXM_out_wSMl____/1e9);
    printf(" %% elrt_sum_svd_SVUXM_lwSM____      : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_svd_SVUXM_lwSM____,n_op_sum_svd_SVUXM_lwSM____/elrt_sum_svd_SVUXM_lwSM____/1e6,n_op_sum_svd_SVUXM_lwSM____/elrt_sum_svd_SVUXM_lwSM____/1e9);
    printf(" %% elrt_sum_f_svd_SVUXM_lwSM____    : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_f_svd_SVUXM_lwSM____,n_op_sum_f_svd_SVUXM_lwSM____/elrt_sum_f_svd_SVUXM_lwSM____/1e6,n_op_sum_f_svd_SVUXM_lwSM____/elrt_sum_f_svd_SVUXM_lwSM____/1e9);
    printf(" %% elrt_sum_svd_SVUXM_lwsM____      : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_svd_SVUXM_lwsM____,n_op_sum_svd_SVUXM_lwsM____/elrt_sum_svd_SVUXM_lwsM____/1e6,n_op_sum_svd_SVUXM_lwsM____/elrt_sum_svd_SVUXM_lwsM____/1e9);
    printf(" %% elrt_sum_f_svd_SVUXM_lwsM____    : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_f_svd_SVUXM_lwsM____,n_op_sum_f_svd_SVUXM_lwsM____/elrt_sum_f_svd_SVUXM_lwsM____/1e6,n_op_sum_f_svd_SVUXM_lwsM____/elrt_sum_f_svd_SVUXM_lwsM____/1e9);
    printf(" %% elrt_sum_svd_USESVUXM_dwSM____   : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_svd_USESVUXM_dwSM____,n_op_sum_svd_USESVUXM_dwSM____/elrt_sum_svd_USESVUXM_dwSM____/1e6,n_op_sum_svd_USESVUXM_dwSM____/elrt_sum_svd_USESVUXM_dwSM____/1e9);
    printf(" %% elrt_sum_f_svd_USESVUXM_dwSM____ : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_f_svd_USESVUXM_dwSM____,n_op_sum_f_svd_USESVUXM_dwSM____/elrt_sum_f_svd_USESVUXM_dwSM____/1e6,n_op_sum_f_svd_USESVUXM_dwSM____/elrt_sum_f_svd_USESVUXM_dwSM____/1e9);
    printf(" %% elrt_sum_l2_dSM___               : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_l2_dSM___,n_op_sum_l2_dSM___/elrt_sum_l2_dSM___/1e6,n_op_sum_l2_dSM___/elrt_sum_l2_dSM___/1e9);
    printf(" %% elrt_sum_X_sub_dwSM____          : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_X_sub_dwSM____,n_op_sum_X_sub_dwSM____/elrt_sum_X_sub_dwSM____/1e6,n_op_sum_X_sub_dwSM____/elrt_sum_X_sub_dwSM____/1e9);
    printf(" %% elrt_sum_d_X_sub_dwSM____        : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_d_X_sub_dwSM____,n_op_sum_d_X_sub_dwSM____/elrt_sum_d_X_sub_dwSM____/1e6,n_op_sum_d_X_sub_dwSM____/elrt_sum_d_X_sub_dwSM____/1e9);
    printf(" %% elrt_sum_I_value_use_dwSM____    : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_I_value_use_dwSM____,n_op_sum_I_value_use_dwSM____/elrt_sum_I_value_use_dwSM____/1e6,n_op_sum_I_value_use_dwSM____/elrt_sum_I_value_use_dwSM____/1e9);
    printf(" %% elrt_sum_X_wSM___                : %0.6f <-- %0.2fMhz %0.2fGhz \n",elrt_sum_X_wSM___,n_op_sum_X_wSM___/elrt_sum_X_wSM___/1e6,n_op_sum_X_wSM___/elrt_sum_X_wSM___/1e9);
    /* if (verbose){ } */}
  /* %%%%%%%%%%%%%%%% */
  free(gamma_z_); gamma_z_=NULL;
  if (flag_compute_I_value){ free(I_value_use_dwSM____); I_value_use_dwSM____=NULL;}
  free(X_sub_dwSM____); X_sub_dwSM____=NULL;
  free(l2_dSM___); l2_dSM___=NULL;
  free(n2_dSM___); n2_dSM___=NULL;
  free(f2_dSM___); f2_dSM___=NULL;
  if (flag_slow){
    free(svd_USESVUXM_dwSM____); svd_USESVUXM_dwSM____=NULL;
    free(svd_SVUXM_lwsM____); svd_SVUXM_lwsM____=NULL;
    free(svd_SVUXM_lwSM____); svd_SVUXM_lwSM____=NULL;
    fftw_free(svd_SVUXM_0in_wlSM____); svd_SVUXM_0in_wlSM____=NULL;
    fftw_free(svd_SVUXM_out_wlSM____); svd_SVUXM_out_wlSM____=NULL;
    fftw_destroy_plan(fftw_plan_many_plan);
    /* if (flag_slow){ } */}
  fftwf_free(f_svd_SVUXM_0in_wSMl____); f_svd_SVUXM_0in_wSMl____=NULL;
  fftwf_free(f_svd_SVUXM_out_wSMl____); f_svd_SVUXM_out_wSMl____=NULL;
  fftwf_destroy_plan(fftwf_plan_many_plan);
  if (flag_slow){
    free(svd_SVUXM_SMwl____); svd_SVUXM_SMwl____=NULL;
    free(svd_SVUXM_SMlw____); svd_SVUXM_SMlw____=NULL;
    free(svd_VUXM_nMwl____); svd_VUXM_nMwl____=NULL;
    free(svd_VUXM_nMlw____); svd_VUXM_nMlw____=NULL;
    /* if (flag_slow){ } */}
  free(index_S_in_Sbatch_); index_S_in_Sbatch_=NULL;
  free(index_M_in_Mbatch_); index_M_in_Mbatch_=NULL;
  if (flag_slow){
    free(CTF_UX_S_k_q_nSw___); CTF_UX_S_k_q_nSw___=NULL;
    free(FTK_svd_U_d_expiw_s__); FTK_svd_U_d_expiw_s__=NULL;
    free(CTF_UX_S_k_q_wnS__); CTF_UX_S_k_q_wnS__=NULL;
    free(svd_VUXM_lwnM____); svd_VUXM_lwnM____=NULL;
    /* if (flag_slow){ } */}
  _mm_free(f_FTK_svd_U_d_expiw_s_real__); f_FTK_svd_U_d_expiw_s_real__=NULL;
  _mm_free(f_FTK_svd_U_d_expiw_s_imag__); f_FTK_svd_U_d_expiw_s_imag__=NULL;
  _mm_free(f_FTK_svd_U_d_expiw_s_tran_real__); f_FTK_svd_U_d_expiw_s_tran_real__=NULL;
  _mm_free(f_FTK_svd_U_d_expiw_s_tran_imag__); f_FTK_svd_U_d_expiw_s_tran_imag__=NULL;
  _mm_free(f_CTF_UX_S_k_q_nSw_real___); f_CTF_UX_S_k_q_nSw_real___=NULL;
  _mm_free(f_CTF_UX_S_k_q_nSw_imag___); f_CTF_UX_S_k_q_nSw_imag___=NULL;
  _mm_free(f_CTF_UX_S_k_q_wnS_real___); f_CTF_UX_S_k_q_wnS_real___=NULL;
  _mm_free(f_CTF_UX_S_k_q_wnS_imag___); f_CTF_UX_S_k_q_wnS_imag___=NULL;
  _mm_free(f_svd_VUXM_lwnM_real____); f_svd_VUXM_lwnM_real____=NULL;
  _mm_free(f_svd_VUXM_lwnM_imag____); f_svd_VUXM_lwnM_imag____=NULL;  
  _mm_free(f_svd_VUXM_nMlw_real____); f_svd_VUXM_nMlw_real____=NULL;
  _mm_free(f_svd_VUXM_nMlw_imag____); f_svd_VUXM_nMlw_imag____=NULL;  
  _mm_free(f_svd_SVUXM_SMlw_real____); f_svd_SVUXM_SMlw_real____=NULL;
  _mm_free(f_svd_SVUXM_SMlw_imag____); f_svd_SVUXM_SMlw_imag____=NULL;
  _mm_free(f_svd_SVUXM_wSMl_real____); f_svd_SVUXM_wSMl_real____=NULL;
  _mm_free(f_svd_SVUXM_wSMl_imag____); f_svd_SVUXM_wSMl_imag____=NULL;
  _mm_free(f_svd_SVUXM_lwsM_real____); f_svd_SVUXM_lwsM_real____=NULL;
  _mm_free(f_svd_SVUXM_lwsM_imag____); f_svd_SVUXM_lwsM_imag____=NULL;
  _mm_free(f_svd_SVUXM_lwSM_real____); f_svd_SVUXM_lwSM_real____=NULL;
  _mm_free(f_svd_SVUXM_lwSM_imag____); f_svd_SVUXM_lwSM_imag____=NULL;
  _mm_free(f_svd_USESVUXM_dwSM_real____); f_svd_USESVUXM_dwSM_real____=NULL;
  _mm_free(f_svd_USESVUXM_dwSM_imag____); f_svd_USESVUXM_dwSM_imag____=NULL;
  free(d_X_sub_dwSM____); d_X_sub_dwSM____=NULL;
  /* %%%%%%%%%%%%%%%% */
  if (verbose){ printf(" %% [finished mex_ampmh_X_wSM___13]\n");}
}
