/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* The helper function */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
void mex_ampmh_X_wSM___fftw_slow_helper
(
 int verbose
 ,int n_M_per_Mbatch
 ,int n_M_sub
 ,int n_S
 ,int FTK_n_svd_l
 ,int n_w_max
 ,fftw_complex *svd_SVUXM_0in_wlSM____
 ,fftw_complex *svd_SVUXM_out_wlSM____
 ,fftw_plan fftw_plan_guru_split_dft_plan
 )
{
  /* %%%%; */
  unsigned long long int tab=0,tabA=0,tabB=0,tabC=0;
  /* %%%% */
  if (verbose){ printf(" %% [entering mex_ampmh_X_wSM___fftw_slow_helper]\n");}
  /* %%%%%%%%%%%%%%%% */
  if (n_M_sub>0){
    /* %%%% */
    /* svd_SVUXM_wlSM____ = ifft(svd_SVUXM_wlSM____,[],1); */
    /* %%%% */
    tab =
      (unsigned long long int)n_w_max
      *(unsigned long long int)n_S
      *(unsigned long long int)n_M_per_Mbatch
      *(unsigned long long int)FTK_n_svd_l;
    memset(svd_SVUXM_out_wlSM____,0,tab*sizeof(fftw_complex));
    tabB =
      (unsigned long long int)n_S
      *(unsigned long long int)n_M_sub
      *(unsigned long long int)FTK_n_svd_l;
    for (tabA=0;tabA<tabB;tabA++){
      tabC = tabA*(unsigned long long int)n_w_max;
      fftw_execute_split_dft(
			     fftw_plan_guru_split_dft_plan
			     ,1 + (double *)(svd_SVUXM_0in_wlSM____ + tabC)
			     ,0 + (double *)(svd_SVUXM_0in_wlSM____ + tabC)
			     ,1 + (double *)(svd_SVUXM_out_wlSM____ + tabC)
			     ,0 + (double *)(svd_SVUXM_out_wlSM____ + tabC)
			     ); /* %<-- note, for backwards dft we swap real and imag. ; */
      /* for (tabA=0;tabA<tabB;tabA++){ } */}
    /* %%%% */
    /* if (n_M_sub>0){ } */}
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  if (verbose){ printf(" %% [finished mex_ampmh_X_wSM___fftw_slow_helper]\n");}
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* slow fftw */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
void mex_ampmh_X_wSM___fftw_slow_test()
{
  int verbose = 1;
  int n_M_per_Mbatch = 24;
  int n_M_sub = n_M_per_Mbatch-1; //%<-- just to check dimensions. ;
  int n_S = 993;
  int FTK_n_svd_l = 24;
  int n_w_max = 98;
  int n_M_per_Mbatch_rup=0,n_M_per_Mbatch_256=0;
  int FTK_n_svd_l_rup=0,FTK_n_svd_l_256=0;
  int n_w_max_rup=0,n_w_max_256=0;
  int n_S_rup=0,n_S_256=0;
  fftw_plan fftw_plan_guru_split_dft_plan;
  fftw_iodim fftw_iodim_use;
  fftw_complex * svd_SVUXM_0inout_wlSM____=NULL;
  fftw_complex * svd_SVUXM_0in_wlSM____=NULL;
  fftw_complex * svd_SVUXM_out_wlSM____=NULL;
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
  if (verbose){ printf(" %% [entering mex_ampmh_X_wSM___fftw_slow_helper]\n");}
  /* %%%%%%%%%%%%%%%% */
  FTK_n_svd_l_rup = rup(FTK_n_svd_l,8); FTK_n_svd_l_256 = FTK_n_svd_l_rup/8;
  n_w_max_rup = rup(n_w_max,8); n_w_max_256 = n_w_max_rup/8;
  n_S_rup = rup(n_S,8); n_S_256 = n_S_rup/8;
  n_M_per_Mbatch_rup = rup(n_M_per_Mbatch,8); n_M_per_Mbatch_256 = n_M_per_Mbatch_rup/8;
  tab = 
    (unsigned long long int)n_w_max_rup
    *(unsigned long long int)FTK_n_svd_l_rup
    *(unsigned long long int)n_S_rup
    *(unsigned long long int)n_M_per_Mbatch_rup;
  svd_SVUXM_0inout_wlSM____ = (fftw_complex *) _mm_malloc(tab/1*sizeof(__m256),32); //%<-- 8 floats per __m256 --> 2 double complex per __m256 ;
  svd_SVUXM_0in_wlSM____ = 0*tab + (fftw_complex *) svd_SVUXM_0inout_wlSM____;
  svd_SVUXM_out_wlSM____ = 1*tab + (fftw_complex *) svd_SVUXM_0inout_wlSM____;
  fftw_iodim_use.n = n_w_max;
  fftw_iodim_use.is = 2;
  fftw_iodim_use.os = 2;
  fftw_plan_guru_split_dft_plan = fftw_plan_guru_split_dft(
							   1
							   ,&fftw_iodim_use
							   ,0
							   ,NULL
							   ,1 + (double *)svd_SVUXM_0in_wlSM____
							   ,0 + (double *)svd_SVUXM_0in_wlSM____
							   ,1 + (double *)svd_SVUXM_out_wlSM____
							   ,0 + (double *)svd_SVUXM_out_wlSM____
							   ,FFTW_MEASURE
							   ); /* %<-- note, for backwards dft we swap real and imag. ; */
  memset(svd_SVUXM_0in_wlSM____,0,tab*sizeof(fftw_complex));
  local_tic(0,t_start_,d_start_);
  ulli=0;
  for (tabA=0;tabA<tab;tabA++){
    svd_SVUXM_0in_wlSM____[tabA] = (fftw_complex) ((float complex) (((int)ulli%7)-3) + _Complex_I * (float complex) (((int)ulli%7)-2));
    ulli++;
    /* for (tabA=0;tabA<tab;tabA++){ } */}
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," initialization: ");
  local_tic(0,t_start_,d_start_);
  mex_ampmh_X_wSM___fftw_slow_helper
    (
     verbose
     ,n_M_per_Mbatch
     ,n_M_sub
     ,n_S
     ,FTK_n_svd_l
     ,n_w_max
     ,svd_SVUXM_0in_wlSM____
     ,svd_SVUXM_out_wlSM____
     ,fftw_plan_guru_split_dft_plan
     );
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose," mex_ampmh_X_wSM___fftw_slow_helper: ");
  _mm_free(svd_SVUXM_0inout_wlSM____); svd_SVUXM_0inout_wlSM____=NULL;
  fftw_destroy_plan(fftw_plan_guru_split_dft_plan);
  /* %%%%%%%%%%%%%%%% */
  if (verbose){ printf(" %% [finished mex_ampmh_X_wSM___fftw_slow_helper]\n");}
}


