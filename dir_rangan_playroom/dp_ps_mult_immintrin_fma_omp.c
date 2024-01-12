/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* The helper function */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
void dp_ps_mult_immintrin_fma_omp_helper
(
 int n_row_A
 ,int n_col_X
 ,__m256 *ps_A_trn__
 ,int n_row_B
 ,__m256 *ps_B_trn__
 ,float **f_C_p_
 /* %%%% */
 ,int nrbatch
 ,int n_rbatch
 )
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
  int nrow_B_per=0,nrow_B_start=0,nrow_B_final=0;
  char str_thisfunction[] = "hp_segregated_to_segregated_mult_immintrin_load1_fma_omp_helper";
  f_C__=NULL;
  if (f_C_p_!=NULL){
    if ( (*f_C_p_)==NULL ){ printf(" %% Warning, (*f_C_p_)==NULL in %s\n",str_thisfunction);}
    f_C__ = *f_C_p_;
    /* if (f_C_p_!=NULL){ } */}
  if (f_C__!=NULL){
    tmp_ps_B_trn__ = ps_B_trn__;
    nrow_B_per = (int)ceil((double)n_row_B/maximum(1,(double)n_rbatch));
    nrow_B_start = maximum(0,minimum(n_row_B-1,(nrbatch+0)*nrow_B_per));
    nrow_B_final = maximum(0,minimum(n_row_B-1,(nrbatch+1)*nrow_B_per));
    tmp_ps_B_trn__ += (unsigned long long int)n_col_X_256*(unsigned long long int)nrow_B_start;
    for (nrow_B=nrow_B_start;nrow_B<=nrow_B_final;nrow_B++){
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
      /* for (nrow_B=nrow_B_start;nrow_B<=nrow_B_final;nrow_B++){ } */}
    /* if (f_C__!=NULL){ } */}
}

void dp_ps_mult_immintrin_fma_omp
(
 int n_row_A
 ,int n_col_X
 ,__m256 *ps_A_trn__
 ,int n_row_B
 ,__m256 *ps_B_trn__
 ,float **f_C_p_
 )
{
  int flag_omp = 1;
  float *f_C__=NULL;
  int nrbatch=0,n_r_per_rbatch=0,n_rbatch=0;
  f_C__=NULL;
  if (f_C_p_!=NULL){
    if ( (*f_C_p_)==NULL ){ (*f_C_p_) = (float *) malloc((unsigned long long int)n_row_A*(unsigned long long int)n_row_B*sizeof(float));}
    f_C__ = *f_C_p_;
    /* if (f_C_p_!=NULL){ } */}
  if (f_C__!=NULL){
    n_r_per_rbatch = 64;
    n_rbatch = ceil((double)n_row_B/(double)n_r_per_rbatch);
    if (flag_omp==0){
      for (nrbatch=0;nrbatch<n_rbatch;nrbatch++){
	dp_ps_mult_immintrin_fma_omp_helper
	  (
	   n_row_A
	   ,n_col_X
	   ,ps_A_trn__
	   ,n_row_B
	   ,ps_B_trn__
	   ,f_C_p_
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
	  dp_ps_mult_immintrin_fma_omp_helper
	    (
	     n_row_A
	     ,n_col_X
	     ,ps_A_trn__
	     ,n_row_B
	     ,ps_B_trn__
	     ,f_C_p_
	     ,nrbatch
	     ,n_rbatch
	     );
	  /* for (nrbatch=0;nrbatch<n_rbatch;nrbatch++){ } */}
	/* end omp parallel */}
      /* if (flag_omp==1){ } */}
    /* if (f_C__!=NULL){ } */}
}
