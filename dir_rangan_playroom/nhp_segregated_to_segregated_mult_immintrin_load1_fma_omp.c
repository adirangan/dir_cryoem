/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* The helper function */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
void nhp_segregated_to_segregated_mult_immintrin_load1_fma_omp_helper
(
 int n_row_A
 ,int n_col_X
 ,float *f_AR_trn__
 ,float *f_AI_trn__
 ,int n_row_B
 ,float *f_BR_trn__
 ,float *f_BI_trn__
 ,float **f_CR_p_
 ,float **f_CI_p_
  /* %%%%; */
 ,int nrbatch
 ,int n_rbatch
 )
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
  int nrow_B_per=0,nrow_B_start=0,nrow_B_final=0;
  unsigned long long int na=0;
  unsigned long long int na_start=0,na_final=0;
  char str_thisfunction[] = "nhp_segregated_to_segregated_mult_immintrin_load1_fma_omp_helper";
  f_CR__=NULL;
  if (f_CR_p_!=NULL){
    if ( (*f_CR_p_)==NULL ){ printf(" %% Warning, (*f_CR_p_)==NULL in %s\n",str_thisfunction);}
    f_CR__ = *f_CR_p_;
    /* if (f_CR_p_!=NULL){ } */}
  f_CI__=NULL;
  if (f_CI_p_!=NULL){
    if ( (*f_CI_p_)==NULL ){ printf(" %% Warning, (*f_CI_p_)==NULL in %s\n",str_thisfunction);}
    f_CI__ = *f_CI_p_;
    /* if (f_CI_p_!=NULL){ } */}
  if ((f_CR__!=NULL) && (f_CI__!=NULL)){
    nrow_B_per = (int)ceil((double)n_row_B/maximum(1,(double)n_rbatch));
    nrow_B_start = maximum(0,minimum(n_row_B-1,(nrbatch+0)*nrow_B_per));
    nrow_B_final = maximum(0,minimum(n_row_B-1,(nrbatch+1)*nrow_B_per));
    na=(unsigned long long int)nrow_B_start*(unsigned long long int)n_row_A;
    for (nrow_B=nrow_B_start;nrow_B<=nrow_B_final;nrow_B++){
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
      /* for (nrow_B=nrow_B_start;nrow_B<=nrow_B_final;nrow_B++){ } */}
    /* if ((f_CR__!=NULL) && (f_CI__!=NULL)){ } */}
}

void nhp_segregated_to_segregated_mult_immintrin_load1_fma_omp
(
 int n_row_A
 ,int n_col_X
 ,float *f_AR_trn__
 ,float *f_AI_trn__
 ,int n_row_B
 ,float *f_BR_trn__
 ,float *f_BI_trn__
 ,float **f_CR_p_
 ,float **f_CI_p_
 )
{
  /* non hermitian product */
  /* assumes alignment */
  int flag_omp = 1;
  float *f_CR__=NULL;
  float *f_CI__=NULL;
  int nrbatch=0,n_r_per_rbatch=0,n_rbatch=0;
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
    n_r_per_rbatch = 64;
    n_rbatch = ceil((double)n_row_B/(double)n_r_per_rbatch);
    if (flag_omp==0){
      for (nrbatch=0;nrbatch<n_rbatch;nrbatch++){
	nhp_segregated_to_segregated_mult_immintrin_load1_fma_omp_helper
	  (
	   n_row_A
	   ,n_col_X
	   ,f_AR_trn__
	   ,f_AI_trn__
	   ,n_row_B
	   ,f_BR_trn__
	   ,f_BI_trn__
	   ,f_CR_p_
	   ,f_CI_p_
	   ,nrbatch
	   ,n_rbatch
	   );
	/* for (nrbatch=0;nrbatch<n_rbatch;nrbatch++){ } */}
      /* if (flag_omp==0){ } */}
    if (flag_omp==1){
#pragma omp parallel private(nrbatch)
      { /* begin omp parallel */
	nrbatch = 0;
#pragma omp for schedule(dynamic)
	for (nrbatch=0;nrbatch<n_rbatch;nrbatch++){
	  nhp_segregated_to_segregated_mult_immintrin_load1_fma_omp_helper
	    (
	     n_row_A
	     ,n_col_X
	     ,f_AR_trn__
	     ,f_AI_trn__
	     ,n_row_B
	     ,f_BR_trn__
	     ,f_BI_trn__
	     ,f_CR_p_
	     ,f_CI_p_
	     ,nrbatch
	     ,n_rbatch
	     );
	  /* for (nrbatch=0;nrbatch<n_rbatch;nrbatch++){ } */}
	/* end omp parallel */}
      /* if (flag_omp==1){ } */}
    /* if ((f_CR__!=NULL) && (f_CI__!=NULL)){ } */}
}
