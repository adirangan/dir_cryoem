/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* The helper function */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
void mex_ampmh_X_wSM___fftwf_fast_helper
(
 int verbose
 ,int n_M_per_Mbatch
 ,int n_M_sub
 ,int n_S
 ,int FTK_n_svd_l
 ,int n_w_max
 ,float *f_svd_SVUXM_0in_wSMl_real____
 ,float *f_svd_SVUXM_0in_wSMl_imag____
 ,float *f_svd_SVUXM_out_wSMl_real____
 ,float *f_svd_SVUXM_out_wSMl_imag____
 ,fftwf_plan fftwf_plan_guru_split_dft_plan
 )
{
  /* %%%%; */
  unsigned long long int tabA=0,tabB=0,tabC=0;
  int n_M_per_Mbatch_rup=0,n_M_per_Mbatch_256=0;
  int n_M_sub_rup=0,n_M_sub_256=0;
  int FTK_n_svd_l_rup=0,FTK_n_svd_l_256=0;
  int n_w_max_rup=0,n_w_max_256=0;
  int n_S_rup=0,n_S_256=0;
  /* %%%% */
  if (verbose){ printf(" %% [entering mex_ampmh_X_wSM___fftwf_fast_helper]\n");}
  FTK_n_svd_l_rup = rup(FTK_n_svd_l,8); FTK_n_svd_l_256 = FTK_n_svd_l_rup/8;
  n_w_max_rup = rup(n_w_max,8); n_w_max_256 = n_w_max_rup/8;
  n_S_rup = rup(n_S,8); n_S_256 = n_S_rup/8;
  n_M_per_Mbatch_rup = rup(n_M_per_Mbatch,8); n_M_per_Mbatch_256 = n_M_per_Mbatch_rup/8;
  n_M_sub_rup = rup(n_M_sub,8); n_M_sub_256 = n_M_sub_rup/8;
  /* %%%%%%%%%%%%%%%% */
  if (n_M_sub>0){
    /* %%%% */
    /* svd_SVUXM_wSMl____ = ifft(svd_SVUXM_wSMl____,[],1); */
    /* %%%% */
    tabA =
      (unsigned long long int)n_w_max_rup
      *(unsigned long long int)n_S_rup
      *(unsigned long long int)n_M_per_Mbatch_rup
      *(unsigned long long int)FTK_n_svd_l_rup;
    memset(f_svd_SVUXM_out_wSMl_real____,0,tabA*sizeof(float));
    memset(f_svd_SVUXM_out_wSMl_imag____,0,tabA*sizeof(float));
    tabB =
      (unsigned long long int)n_S_rup
      *(unsigned long long int)n_M_sub_rup
      *(unsigned long long int)FTK_n_svd_l_rup;
    for (tabA=0;tabA<tabB;tabA++){
      tabC = tabA*(unsigned long long int)n_w_max_rup;
      fftwf_execute_split_dft(
			      fftwf_plan_guru_split_dft_plan
			      ,f_svd_SVUXM_0in_wSMl_imag____ + tabC
			      ,f_svd_SVUXM_0in_wSMl_real____ + tabC
			      ,f_svd_SVUXM_out_wSMl_imag____ + tabC
			      ,f_svd_SVUXM_out_wSMl_real____ + tabC
			      ); /* %<-- note, for backwards dft we swap real and imag. ; */
      /* for (tabA=0;tabA<tabB;tabA++){ } */}
    /* %%%% */
    /* if (n_M_sub>0){ } */}
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  if (verbose){ printf(" %% [finished mex_ampmh_X_wSM___fftwf_fast_helper]\n");}
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* fast fftwf */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
void mex_ampmh_X_wSM___fftwf_fast_test()
{
  int verbose = 1;
  int n_M_per_Mbatch = 24;
  int n_M_sub = n_M_per_Mbatch-1; //%<-- just to check dimensions. ;
  int n_S = 993;
  int FTK_n_svd_l = 24;
  int n_w_max = 98;
  int n_M_per_Mbatch_rup=0,n_M_per_Mbatch_256=0;
  int n_M_sub_rup=0,n_M_sub_256=0;
  int FTK_n_svd_l_rup=0,FTK_n_svd_l_256=0;
  int n_w_max_rup=0,n_w_max_256=0;
  int n_S_rup=0,n_S_256=0;
  /* %%%% */
  fftwf_plan fftwf_plan_guru_split_dft_plan;
  fftwf_iodim fftwf_iodim_use;
  float * f_svd_SVUXM_0inout_wSMl_realimag____=NULL;
  float * f_svd_SVUXM_0in_wSMl_real____=NULL;
  float * f_svd_SVUXM_0in_wSMl_imag____=NULL;
  float * f_svd_SVUXM_out_wSMl_real____=NULL;
  float * f_svd_SVUXM_out_wSMl_imag____=NULL;
  /* %%%% */
  unsigned long long int tab=0,tabA=0,tabB=0,ulli=0;
  /* %%%% */
  int ntick=0,n_tick = 8;
  clock_t t_start_[n_tick];
  clock_t t_final_[n_tick];
  struct timeval d_start_[n_tick];
  struct timeval d_final_[n_tick];
  long l_msec_[n_tick],l_ssec_[n_tick],l_usec_[n_tick];
  double elct_[n_tick],elrt_[n_tick];
  /* %%%% */
  if (verbose){ printf(" %% [entering mex_ampmh_X_wSM___fftwf_fast_helper]\n");}
  /* %%%%%%%%%%%%%%%% */
  FTK_n_svd_l_rup = rup(FTK_n_svd_l,8); FTK_n_svd_l_256 = FTK_n_svd_l_rup/8;
  n_w_max_rup = rup(n_w_max,8); n_w_max_256 = n_w_max_rup/8;
  n_S_rup = rup(n_S,8); n_S_256 = n_S_rup/8;
  n_M_per_Mbatch_rup = rup(n_M_per_Mbatch,8); n_M_per_Mbatch_256 = n_M_per_Mbatch_rup/8;
  n_M_sub_rup = rup(n_M_sub,8); n_M_sub_256 = n_M_sub_rup/8;
  tab = 
    (unsigned long long int)n_w_max_rup
    *(unsigned long long int)n_S_rup
    *(unsigned long long int)n_M_per_Mbatch_rup
    *(unsigned long long int)FTK_n_svd_l_rup;
  f_svd_SVUXM_0inout_wSMl_realimag____ = (float *) _mm_malloc(tab/2*sizeof(__m256),32); //%<-- 8 floats per __m256 --> 4 float complex per __m256 ;
  f_svd_SVUXM_0in_wSMl_real____ = 0*tab + (float *) f_svd_SVUXM_0inout_wSMl_realimag____;
  f_svd_SVUXM_0in_wSMl_imag____ = 1*tab + (float *) f_svd_SVUXM_0inout_wSMl_realimag____;
  f_svd_SVUXM_out_wSMl_real____ = 2*tab + (float *) f_svd_SVUXM_0inout_wSMl_realimag____;
  f_svd_SVUXM_out_wSMl_imag____ = 3*tab + (float *) f_svd_SVUXM_0inout_wSMl_realimag____;
  fftwf_iodim_use.n = n_w_max;
  fftwf_iodim_use.is = 1;
  fftwf_iodim_use.os = 1;
  fftwf_plan_guru_split_dft_plan = fftwf_plan_guru_split_dft(
							     1
							     ,&fftwf_iodim_use
							     ,0
							     ,NULL
							     ,f_svd_SVUXM_0in_wSMl_imag____
							     ,f_svd_SVUXM_0in_wSMl_real____
							     ,f_svd_SVUXM_out_wSMl_imag____
							     ,f_svd_SVUXM_out_wSMl_real____
							     ,FFTW_MEASURE
							     ); /* %<-- note, for backwards dft we swap real and imag. ; */
  memset(f_svd_SVUXM_0in_wSMl_real____,0,tab*sizeof(float));
  memset(f_svd_SVUXM_0in_wSMl_imag____,0,tab*sizeof(float));
  local_tic(0,t_start_,d_start_);
  for (tabB=0;tabB<FTK_n_svd_l*n_S*n_M_per_Mbatch;tabB++){
    for (tabA=0;tabA<n_w_max;tabA++){
      ulli = tabA + tabB*n_w_max_rup;
      f_svd_SVUXM_0in_wSMl_real____[ulli] = (float) tabA;
      f_svd_SVUXM_0in_wSMl_imag____[ulli] = (float) tabB;
      /* for (tabA=0;tabA<n_w_max;tabA++){ } */}
    /* for (tabB=0;tabB<FTK_n_svd_l*n_S*n_M_per_Mbatch;tabB++){ } */}
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," initialization: ");
  local_tic(0,t_start_,d_start_);
  mex_ampmh_X_wSM___fftwf_fast_helper
    (
     verbose
     ,n_M_per_Mbatch
     ,n_M_sub
     ,n_S
     ,FTK_n_svd_l
     ,n_w_max
     ,f_svd_SVUXM_0in_wSMl_real____
     ,f_svd_SVUXM_0in_wSMl_imag____
     ,f_svd_SVUXM_out_wSMl_real____
     ,f_svd_SVUXM_out_wSMl_imag____
     ,fftwf_plan_guru_split_dft_plan
     );
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," mex_ampmh_X_wSM___fftwf_fast_helper: ");
  if (verbose){
    array_sub_printf
      (
       f_svd_SVUXM_0in_wSMl_real____
       ,"float"
       ,n_w_max_rup
       ,4
       ,n_S_rup*n_M_sub_rup*FTK_n_svd_l_rup
       ,5
       ," %% f_svd_SVUXM_0in_wSMl_real____: "
       );
      /* if (verbose){ } */
    array_sub_printf
      (
       f_svd_SVUXM_0in_wSMl_imag____
       ,"float"
       ,n_w_max_rup
       ,4
       ,n_S_rup*n_M_sub_rup*FTK_n_svd_l_rup
       ,5
       ," %% f_svd_SVUXM_0in_wSMl_imag____: "
       );
    array_sub_printf
      (
       f_svd_SVUXM_out_wSMl_real____
       ,"float"
       ,n_w_max_rup
       ,4
       ,n_S_rup*n_M_sub_rup*FTK_n_svd_l_rup
       ,5
       ," %% f_svd_SVUXM_out_wSMl_real____: "
       );
      /* if (verbose){ } */
    array_sub_printf
      (
       f_svd_SVUXM_out_wSMl_imag____
       ,"float"
       ,n_w_max_rup
       ,4
       ,n_S_rup*n_M_sub_rup*FTK_n_svd_l_rup
       ,5
       ," %% f_svd_SVUXM_out_wSMl_imag____: "
       );
      /* if (verbose){ } */}
  _mm_free(f_svd_SVUXM_0inout_wSMl_realimag____); f_svd_SVUXM_0inout_wSMl_realimag____=NULL;
  fftwf_destroy_plan(fftwf_plan_guru_split_dft_plan);
  /* %%%%%%%%%%%%%%%% */
  if (verbose){ printf(" %% [finished mex_ampmh_X_wSM___fftwf_fast_helper]\n");}
}


