/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* The helper function */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
void mex_ampmh_X_wSM___fft_mkl_n1po_helper
(
 int verbose
 ,int n_M
 ,int n_M_sub
 ,int n_S
 ,int FTK_n_svd_l
 ,int n_w_max
 ,float complex *f0_svd_SVUXM_0in_wlSM____
 ,float complex *f0_svd_SVUXM_out_wlSM____
 )
{
  MKL_LONG mkll_n_dfti=0; //%<-- DFTI_NUMBER_OF_TRANSFORMS ;
  MKL_LONG mkll_i_dfti=0; //%<-- DFTI_INPUT_DISTANCE ;
  MKL_LONG mkll_err=0;
  DFTI_DESCRIPTOR *mkl_dfti_fft1d_p_ = NULL;  
  /* %%%%; */
  unsigned long long int tab=0,tabA=0,tabB=0,tabC=0;
  /* %%%% */
  if (verbose){ printf(" %% [entering mex_ampmh_X_wSM___fft_mkl_n1po_helper]\n");}
  /* %%%%%%%%%%%%%%%% */
  if (n_M_sub>0){
    /* %%%% */
    /* svd_SVUXM_wlSM____ = ifft(svd_SVUXM_wlSM____,[],1); */
    /* %%%% */
    tab =
      (unsigned long long int)n_w_max
      *(unsigned long long int)FTK_n_svd_l
      *(unsigned long long int)n_S
      *(unsigned long long int)n_M;
    memset(f0_svd_SVUXM_out_wlSM____,0,tab*sizeof(float complex));
    mkll_i_dfti = n_w_max; mkll_n_dfti = 1; //%<-- only do one at a time. ;
    mkll_err = DftiCreateDescriptor(&mkl_dfti_fft1d_p_,DFTI_SINGLE,DFTI_COMPLEX,1,mkll_i_dfti);
    mkll_err = DftiSetValue(mkl_dfti_fft1d_p_,DFTI_NUMBER_OF_TRANSFORMS,mkll_n_dfti);
    mkll_err = DftiSetValue(mkl_dfti_fft1d_p_,DFTI_PLACEMENT,DFTI_NOT_INPLACE);
    mkll_err = DftiCommitDescriptor(mkl_dfti_fft1d_p_);
    tabB =
      (unsigned long long int)FTK_n_svd_l
      *(unsigned long long int)n_S
      *(unsigned long long int)n_M_sub;
    for (tabA=0;tabA<tabB;tabA++){
      tabC = tabA*(unsigned long long int)n_w_max;
      mkll_err = DftiComputeBackward
      	(
      	 mkl_dfti_fft1d_p_
      	 ,f0_svd_SVUXM_0in_wlSM____ + tabC
      	 ,f0_svd_SVUXM_out_wlSM____ + tabC
      	 );
      /* for (tabA=0;tabA<tabB;tabA++){ } */}
    mkll_err = DftiFreeDescriptor(&mkl_dfti_fft1d_p_);
    /* %%%% */
    /* if (n_M_sub>0){ } */}
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  if (verbose){ printf(" %% [finished mex_ampmh_X_wSM___fft_mkl_n1po_helper]\n");}
}

void mex_ampmh_X_wSM___fft_mkl_nxpo_helper
(
 int verbose
 ,int n_M
 ,int n_M_sub
 ,int n_S
 ,int FTK_n_svd_l
 ,int n_w_max
 ,float complex *f0_svd_SVUXM_0in_wlSM____
 ,float complex *f0_svd_SVUXM_out_wlSM____
 )
{
  MKL_LONG mkll_n_dfti=0; //%<-- DFTI_NUMBER_OF_TRANSFORMS ;
  MKL_LONG mkll_i_dfti=0; //%<-- DFTI_INPUT_DISTANCE ;
  MKL_LONG mkll_err=0;
  DFTI_DESCRIPTOR *mkl_dfti_fft1d_p_ = NULL;  
  /* %%%%; */
  unsigned long long int tabC=0;
  /* %%%% */
  if (verbose){ printf(" %% [entering mex_ampmh_X_wSM___fft_mkl_nxpo_helper]\n");}
  /* %%%%%%%%%%%%%%%% */
  if (n_M_sub>0){
    /* %%%% */
    /* svd_SVUXM_wlSM____ = ifft(svd_SVUXM_wlSM____,[],1); */
    /* %%%% */
    tabC =
      (unsigned long long int)n_w_max
      *(unsigned long long int)FTK_n_svd_l
      *(unsigned long long int)n_S
      *(unsigned long long int)n_M;
    memset(f0_svd_SVUXM_out_wlSM____,0,tabC*sizeof(float complex));
    mkll_i_dfti = n_w_max; mkll_n_dfti = FTK_n_svd_l*n_S*n_M_sub; //%<-- do all at once. ;
    mkll_err = DftiCreateDescriptor(&mkl_dfti_fft1d_p_,DFTI_SINGLE,DFTI_COMPLEX,1,mkll_i_dfti);
    mkll_err = DftiSetValue(mkl_dfti_fft1d_p_,DFTI_NUMBER_OF_TRANSFORMS,mkll_n_dfti);
    mkll_err = DftiSetValue(mkl_dfti_fft1d_p_,DFTI_INPUT_DISTANCE,mkll_i_dfti);
    mkll_err = DftiSetValue(mkl_dfti_fft1d_p_,DFTI_OUTPUT_DISTANCE,mkll_i_dfti);
    mkll_err = DftiSetValue(mkl_dfti_fft1d_p_,DFTI_PLACEMENT,DFTI_NOT_INPLACE);
    mkll_err = DftiCommitDescriptor(mkl_dfti_fft1d_p_);
    mkll_err = DftiComputeBackward
      (
       mkl_dfti_fft1d_p_
       ,f0_svd_SVUXM_0in_wlSM____
       ,f0_svd_SVUXM_out_wlSM____
       );
    mkll_err = DftiFreeDescriptor(&mkl_dfti_fft1d_p_);
    /* %%%% */
    /* if (n_M_sub>0){ } */}
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  if (verbose){ printf(" %% [finished mex_ampmh_X_wSM___fft_mkl_nxpo_helper]\n");}
}

void mex_ampmh_X_wSM___fft_mkl_nxpo_omp_helper
(
 int verbose
 ,int n_M_per_Mbatch
 ,int n_M
 ,int n_S
 ,int FTK_n_svd_l
 ,int n_w_max
 ,float complex *f0_svd_SVUXM_0in_wlSM____
 ,float complex *f0_svd_SVUXM_out_wlSM____
/* %%%% */
 ,int nMbatch
 ,int n_Mbatch
 )
{
  MKL_LONG mkll_n_dfti=0; //%<-- DFTI_NUMBER_OF_TRANSFORMS ;
  MKL_LONG mkll_i_dfti=0; //%<-- DFTI_INPUT_DISTANCE ;
  MKL_LONG mkll_err=0;
  DFTI_DESCRIPTOR *mkl_dfti_fft1d_p_ = NULL;  
  /* %%%%; */
  unsigned long long int tab=0,tabC=0;
  int n_M_sub=0,nM_sub=0;
  /* %%%% */
  if (verbose){ printf(" %% [entering mex_ampmh_X_wSM___fft_mkl_nxpo_omp_helper]\n");}
  n_M_sub = n_M_per_Mbatch; if (nMbatch==n_Mbatch-1){ n_M_sub = n_M - n_M_per_Mbatch*nMbatch;}
  /* %%%%%%%%%%%%%%%% */
  if (n_M_sub>0){
    /* %%%% */
    /* svd_SVUXM_wlSM____ = ifft(svd_SVUXM_wlSM____,[],1); */
    /* %%%% */
    tabC =
      (unsigned long long int)n_w_max
      *(unsigned long long int)FTK_n_svd_l
      *(unsigned long long int)n_S
      *(unsigned long long int)n_M_sub;
    tab =
      (unsigned long long int)n_w_max
      *(unsigned long long int)FTK_n_svd_l
      *(unsigned long long int)n_S
      *(unsigned long long int)n_M_per_Mbatch
      *(unsigned long long int)nMbatch;
    memset(f0_svd_SVUXM_out_wlSM____+tab,0,tabC*sizeof(float complex));
    mkll_i_dfti = n_w_max; mkll_n_dfti = FTK_n_svd_l*n_S*n_M_sub; //%<-- do all at once. ;
    mkll_err = DftiCreateDescriptor(&mkl_dfti_fft1d_p_,DFTI_SINGLE,DFTI_COMPLEX,1,mkll_i_dfti);
    mkll_err = DftiSetValue(mkl_dfti_fft1d_p_,DFTI_NUMBER_OF_TRANSFORMS,mkll_n_dfti);
    mkll_err = DftiSetValue(mkl_dfti_fft1d_p_,DFTI_INPUT_DISTANCE,mkll_i_dfti);
    mkll_err = DftiSetValue(mkl_dfti_fft1d_p_,DFTI_OUTPUT_DISTANCE,mkll_i_dfti);
    mkll_err = DftiSetValue(mkl_dfti_fft1d_p_,DFTI_PLACEMENT,DFTI_NOT_INPLACE);
    mkll_err = DftiCommitDescriptor(mkl_dfti_fft1d_p_);
    mkll_err = DftiComputeBackward(mkl_dfti_fft1d_p_,f0_svd_SVUXM_0in_wlSM____ + tab,f0_svd_SVUXM_out_wlSM____ + tab);
    mkll_err = DftiFreeDescriptor(&mkl_dfti_fft1d_p_);
    /* %%%% */
    /* if (n_M_sub>0){ } */}
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  if (verbose){ printf(" %% [finished mex_ampmh_X_wSM___fft_mkl_nxpo_omp_helper]\n");}
}

void mex_ampmh_X_wSM___fft_mkl_n1pi_helper
(
 int verbose
 ,int n_M
 ,int n_M_sub
 ,int n_S
 ,int FTK_n_svd_l
 ,int n_w_max
 ,float complex *f0_svd_SVUXM_0in_wlSM____
 )
{
  MKL_LONG mkll_n_dfti=0; //%<-- DFTI_NUMBER_OF_TRANSFORMS ;
  MKL_LONG mkll_i_dfti=0; //%<-- DFTI_INPUT_DISTANCE ;
  MKL_LONG mkll_err=0;
  DFTI_DESCRIPTOR *mkl_dfti_fft1d_p_ = NULL;  
  /* %%%%; */
  unsigned long long int tab=0,tabA=0,tabB=0,tabC=0;
  /* %%%% */
  if (verbose){ printf(" %% [entering mex_ampmh_X_wSM___fft_mkl_n1pi_helper]\n");}
  /* %%%%%%%%%%%%%%%% */
  if (n_M_sub>0){
    /* %%%% */
    /* svd_SVUXM_wlSM____ = ifft(svd_SVUXM_wlSM____,[],1); */
    /* %%%% */
    tab =
      (unsigned long long int)n_w_max
      *(unsigned long long int)FTK_n_svd_l
      *(unsigned long long int)n_S
      *(unsigned long long int)n_M;
    mkll_i_dfti = n_w_max; mkll_n_dfti = 1; //%<-- only do one at a time. ;
    mkll_err = DftiCreateDescriptor(&mkl_dfti_fft1d_p_,DFTI_SINGLE,DFTI_COMPLEX,1,mkll_i_dfti);
    mkll_err = DftiSetValue(mkl_dfti_fft1d_p_,DFTI_NUMBER_OF_TRANSFORMS,mkll_n_dfti);
    mkll_err = DftiSetValue(mkl_dfti_fft1d_p_,DFTI_PLACEMENT,DFTI_INPLACE);
    mkll_err = DftiCommitDescriptor(mkl_dfti_fft1d_p_);
    tabB =
      (unsigned long long int)FTK_n_svd_l
      *(unsigned long long int)n_S
      *(unsigned long long int)n_M_sub;
    for (tabA=0;tabA<tabB;tabA++){
      tabC = tabA*(unsigned long long int)n_w_max;
      mkll_err = DftiComputeBackward
      	(
      	 mkl_dfti_fft1d_p_
      	 ,f0_svd_SVUXM_0in_wlSM____ + tabC
      	 );
      /* for (tabA=0;tabA<tabB;tabA++){ } */}
    mkll_err = DftiFreeDescriptor(&mkl_dfti_fft1d_p_);
    /* %%%% */
    /* if (n_M_sub>0){ } */}
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  if (verbose){ printf(" %% [finished mex_ampmh_X_wSM___fft_mkl_n1pi_helper]\n");}
}

void mex_ampmh_X_wSM___fft_mkl_nxpi_helper
(
 int verbose
 ,int n_M
 ,int n_M_sub
 ,int n_S
 ,int FTK_n_svd_l
 ,int n_w_max
 ,float complex *f0_svd_SVUXM_0in_wlSM____
 )
{
  MKL_LONG mkll_n_dfti=0; //%<-- DFTI_NUMBER_OF_TRANSFORMS ;
  MKL_LONG mkll_i_dfti=0; //%<-- DFTI_INPUT_DISTANCE ;
  MKL_LONG mkll_err=0;
  DFTI_DESCRIPTOR *mkl_dfti_fft1d_p_ = NULL;  
  /* %%%%; */
  if (verbose){ printf(" %% [entering mex_ampmh_X_wSM___fft_mkl_nxpi_helper]\n");}
  /* %%%%%%%%%%%%%%%% */
  if (n_M_sub>0){
    /* %%%% */
    /* svd_SVUXM_wlSM____ = ifft(svd_SVUXM_wlSM____,[],1); */
    /* %%%% */
    mkll_i_dfti = n_w_max; mkll_n_dfti = FTK_n_svd_l*n_S*n_M_sub; //%<-- do all at once. ;
    mkll_err = DftiCreateDescriptor(&mkl_dfti_fft1d_p_,DFTI_SINGLE,DFTI_COMPLEX,1,mkll_i_dfti);
    mkll_err = DftiSetValue(mkl_dfti_fft1d_p_,DFTI_NUMBER_OF_TRANSFORMS,mkll_n_dfti);
    mkll_err = DftiSetValue(mkl_dfti_fft1d_p_,DFTI_INPUT_DISTANCE,mkll_i_dfti);
    mkll_err = DftiSetValue(mkl_dfti_fft1d_p_,DFTI_OUTPUT_DISTANCE,mkll_i_dfti);
    mkll_err = DftiSetValue(mkl_dfti_fft1d_p_,DFTI_PLACEMENT,DFTI_INPLACE);
    mkll_err = DftiCommitDescriptor(mkl_dfti_fft1d_p_);
    mkll_err = DftiComputeBackward(mkl_dfti_fft1d_p_,f0_svd_SVUXM_0in_wlSM____);
    mkll_err = DftiFreeDescriptor(&mkl_dfti_fft1d_p_);
    /* %%%% */
    /* if (n_M_sub>0){ } */}
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  if (verbose){ printf(" %% [finished mex_ampmh_X_wSM___fft_mkl_nxpi_helper]\n");}
}

void mex_ampmh_X_wSM___fft_mkl_nxpi_omp_helper
(
 int verbose
 ,int n_M_per_Mbatch
 ,int n_M
 ,int n_S
 ,int FTK_n_svd_l
 ,int n_w_max
 ,float complex *f0_svd_SVUXM_0in_wlSM____
/* %%%% */
 ,int nMbatch
 ,int n_Mbatch
 )
{
  MKL_LONG mkll_n_dfti=0; //%<-- DFTI_NUMBER_OF_TRANSFORMS ;
  MKL_LONG mkll_i_dfti=0; //%<-- DFTI_INPUT_DISTANCE ;
  MKL_LONG mkll_err=0;
  DFTI_DESCRIPTOR *mkl_dfti_fft1d_p_ = NULL;  
  /* %%%%; */
  unsigned long long int tab=0;
  int n_M_sub=0,nM_sub=0;
  /* %%%% */
  if (verbose){ printf(" %% [entering mex_ampmh_X_wSM___fft_mkl_nxpi_omp_helper]\n");}
  n_M_sub = n_M_per_Mbatch; if (nMbatch==n_Mbatch-1){ n_M_sub = n_M - n_M_per_Mbatch*nMbatch;}
  /* %%%%%%%%%%%%%%%% */
  if (n_M_sub>0){
    /* %%%% */
    /* svd_SVUXM_wlSM____ = ifft(svd_SVUXM_wlSM____,[],1); */
    /* %%%% */
    tab =
      (unsigned long long int)n_w_max
      *(unsigned long long int)FTK_n_svd_l
      *(unsigned long long int)n_S
      *(unsigned long long int)n_M_per_Mbatch
      *(unsigned long long int)nMbatch;
    mkll_i_dfti = n_w_max; mkll_n_dfti = FTK_n_svd_l*n_S*n_M_sub; //%<-- do all at once. ;
    mkll_err = DftiCreateDescriptor(&mkl_dfti_fft1d_p_,DFTI_SINGLE,DFTI_COMPLEX,1,mkll_i_dfti);
    mkll_err = DftiSetValue(mkl_dfti_fft1d_p_,DFTI_NUMBER_OF_TRANSFORMS,mkll_n_dfti);
    mkll_err = DftiSetValue(mkl_dfti_fft1d_p_,DFTI_INPUT_DISTANCE,mkll_i_dfti);
    mkll_err = DftiSetValue(mkl_dfti_fft1d_p_,DFTI_OUTPUT_DISTANCE,mkll_i_dfti);
    mkll_err = DftiSetValue(mkl_dfti_fft1d_p_,DFTI_PLACEMENT,DFTI_INPLACE);
    mkll_err = DftiCommitDescriptor(mkl_dfti_fft1d_p_);
    mkll_err = DftiComputeBackward(mkl_dfti_fft1d_p_,f0_svd_SVUXM_0in_wlSM____ + tab);
    mkll_err = DftiFreeDescriptor(&mkl_dfti_fft1d_p_);
    /* %%%% */
    /* if (n_M_sub>0){ } */}
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  if (verbose){ printf(" %% [finished mex_ampmh_X_wSM___fft_mkl_nxpi_omp_helper]\n");}
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* mkl fft */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
void mex_ampmh_X_wSM___fft_mkl_test()
{
  /* modified from example C18 on page 3001 of the Intel MKL Reference Manual: */
  /* https://sep.stanford.edu/sep/claudio/Research/Prst_ExpRefl/ShtPSPI/intel/mkl/10.0.3.020/doc/mklman.pdf */
  int verbose = 1;
  int n_M = 24*3+3;
  int n_M_use = n_M-0; //%<-- set to less than n_M to check partial fft. ;
  int n_M_sub = 0;
  int n_M_per_Mbatch_use = 12;
  int n_S = 993;
  int FTK_n_svd_l = 24;
  int n_w_max = 98;
  int n_M_rup=0,n_M_256=0;
  int FTK_n_svd_l_rup=0,FTK_n_svd_l_256=0;
  int n_w_max_rup=0,n_w_max_256=0;
  int n_S_rup=0,n_S_256=0;
  fftw_plan fftw_plan_guru_split_dft_plan;
  fftw_iodim fftw_iodim_use;
  fftw_complex * f1_svd_SVUXM_0inout_wlSM____=NULL;
  fftw_complex * f1_svd_SVUXM_0in_wlSM____=NULL;
  fftw_complex * f1_svd_SVUXM_out_wlSM____=NULL;
  float complex * f0_svd_SVUXM_bkp_wlSM____=NULL;
  float complex * f0_svd_SVUXM_0in_wlSM____=NULL;
  float complex * f0_svd_SVUXM_out_wlSM____=NULL;
  unsigned long long int ulli_numel=0,ulli_total=0,tabA=0,tabB=0,ulli=0;
  float rerror=0;
  /* %%%% */
  int ntick=0,n_tick = 8;
  clock_t t_start_[n_tick];
  clock_t t_final_[n_tick];
  struct timeval d_start_[n_tick];
  struct timeval d_final_[n_tick];
  long l_msec_[n_tick],l_ssec_[n_tick],l_usec_[n_tick];
  double elct_[n_tick],elrt_[n_tick];
  /* %%%% */
  int n_M_per_Mbatch=0,nMbatch=0,n_Mbatch=0;
  /* %%%% */
  if (verbose){ printf(" %% [entering mex_ampmh_X_wSM___fft_mkl_test]\n");}
  /* %%%%%%%%%%%%%%%% */
  FTK_n_svd_l_rup = rup(FTK_n_svd_l,8); FTK_n_svd_l_256 = FTK_n_svd_l_rup/8;
  n_w_max_rup = rup(n_w_max,8); n_w_max_256 = n_w_max_rup/8;
  n_S_rup = rup(n_S,8); n_S_256 = n_S_rup/8;
  n_M_rup = rup(n_M,8); n_M_256 = n_M_rup/8;
  if (verbose>1){ printf(" %% FTK_n_svd_l %d(%d,%d*8)\n",FTK_n_svd_l,FTK_n_svd_l_rup,FTK_n_svd_l_256);}
  if (verbose>1){ printf(" %% n_w_max %d(%d,%d*8)\n",n_w_max,n_w_max_rup,n_w_max_256);}
  if (verbose>1){ printf(" %% n_S %d(%d,%d*8)\n",n_S,n_S_rup,n_S_256);}
  if (verbose>1){ printf(" %% n_M %d(%d,%d*8)\n",n_M,n_M_rup,n_M_256);}
  ulli_numel = 
    (unsigned long long int)n_w_max
    *(unsigned long long int)FTK_n_svd_l
    *(unsigned long long int)n_S
    *(unsigned long long int)n_M_use;
  ulli_total = 
    (unsigned long long int)n_w_max_rup
    *(unsigned long long int)FTK_n_svd_l_rup
    *(unsigned long long int)n_S_rup
    *(unsigned long long int)n_M_rup;
  if (verbose>1){ printf(" %% ulli_numel %lld --> %0.2fGB\n",ulli_numel,(double)ulli_numel*sizeof(__m256)/1.0e9);}
  if (verbose>1){ printf(" %% ulli_total %lld --> %0.2fGB\n",ulli_total,(double)ulli_total*sizeof(__m256)/1.0e9);}
  /* %%%% */
  f1_svd_SVUXM_0inout_wlSM____ = (fftw_complex *) _mm_malloc(ulli_total/1*sizeof(__m256),32); //%<-- 8 floats per __m256 --> 2 double complex per __m256 ;
  f1_svd_SVUXM_0in_wlSM____ = 0*ulli_total + (fftw_complex *) f1_svd_SVUXM_0inout_wlSM____;
  f1_svd_SVUXM_out_wlSM____ = 1*ulli_total + (fftw_complex *) f1_svd_SVUXM_0inout_wlSM____;
  fftw_iodim_use.n = n_w_max;
  fftw_iodim_use.is = 2;
  fftw_iodim_use.os = 2;
  fftw_plan_guru_split_dft_plan = fftw_plan_guru_split_dft(
							   1
							   ,&fftw_iodim_use
							   ,0
							   ,NULL
							   ,1 + (double *)f1_svd_SVUXM_0in_wlSM____
							   ,0 + (double *)f1_svd_SVUXM_0in_wlSM____
							   ,1 + (double *)f1_svd_SVUXM_out_wlSM____
							   ,0 + (double *)f1_svd_SVUXM_out_wlSM____
							   ,FFTW_MEASURE
							   ); /* %<-- note, for backwards dft we swap real and imag. ; */
  memset(f1_svd_SVUXM_0in_wlSM____,0,ulli_total*sizeof(fftw_complex));
  memset(f1_svd_SVUXM_out_wlSM____,0,ulli_total*sizeof(fftw_complex));
  local_tic(0,t_start_,d_start_);
  for (tabB=0;tabB<FTK_n_svd_l*n_S*n_M;tabB++){
    for (tabA=0;tabA<n_w_max;tabA++){
      ulli = tabA + tabB*n_w_max;
      f1_svd_SVUXM_0in_wlSM____[ulli] = (fftw_complex) (tabA + _Complex_I * (tabB%13));
      /* for (tabA=0;tabA<n_w_max;tabA++){ } */}
    /* for (tabB=0;tabB<FTK_n_svd_l*n_S*n_M;tabB++){ } */}
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,ulli_total,verbose," initialization: ");
  /* %%%% */
  local_tic(0,t_start_,d_start_);
  mex_ampmh_X_wSM___fftw_slow_helper
    (
     verbose
     ,n_M
     ,n_M_use
     ,n_S
     ,FTK_n_svd_l
     ,n_w_max
     ,f1_svd_SVUXM_0in_wlSM____
     ,f1_svd_SVUXM_out_wlSM____
     ,fftw_plan_guru_split_dft_plan
     );
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,ulli_total,verbose," mex_ampmh_X_wSM___fftw_slow_helper: ");
  /* %%%% */
  f0_svd_SVUXM_bkp_wlSM____ = (float complex *) _mm_malloc(ulli_total*sizeof(float complex),64);
  memset(f0_svd_SVUXM_bkp_wlSM____,0,ulli_total*sizeof(float complex));
  for (tabB=0;tabB<FTK_n_svd_l*n_S*n_M;tabB++){
    for (tabA=0;tabA<n_w_max;tabA++){
      ulli = tabA + tabB*n_w_max;
      f0_svd_SVUXM_bkp_wlSM____[ulli] = (float complex) f1_svd_SVUXM_out_wlSM____[ulli];
      /* for (tabA=0;tabA<n_w_max;tabA++){ } */}
    /* for (tabB=0;tabB<FTK_n_svd_l*n_S*n_M;tabB++){ } */}
  /* %%%% */
  if (verbose>1){
    array_sub_printf
      (
       f1_svd_SVUXM_0in_wlSM____
       ,"double complex"
       ,n_w_max
       ,4
       ,FTK_n_svd_l*n_S*n_M_use
       ,5
       ," %% f1_svd_SVUXM_0in_wlSM____: "
       );
    array_sub_printf
      (
       f1_svd_SVUXM_out_wlSM____
       ,"double complex"
       ,n_w_max
       ,4
       ,FTK_n_svd_l*n_S*n_M_use
       ,5
       ," %% f1_svd_SVUXM_out_wlSM____: "
       );
    array_sub_printf
      (
       f0_svd_SVUXM_bkp_wlSM____
       ,"float complex"
       ,n_w_max
       ,4
       ,FTK_n_svd_l*n_S*n_M_use
       ,5
       ," %% f1_svd_SVUXM_bkp_wlSM____: "
       );
    array_printf_margin(f0_svd_SVUXM_bkp_wlSM____,"float complex",n_w_max,FTK_n_svd_l*n_S*n_M_use," %% f0_svd_SVUXM_bkp_wlSM____: ");
    /* if (verbose){ } */}
  /* %%%% */
  f0_svd_SVUXM_0in_wlSM____ = (float complex *) _mm_malloc(ulli_total*sizeof(float complex),64);
  f0_svd_SVUXM_out_wlSM____ = (float complex *) _mm_malloc(ulli_total*sizeof(float complex),64);
  memset(f0_svd_SVUXM_0in_wlSM____,0,ulli_total*sizeof(float complex));
  local_tic(0,t_start_,d_start_);
  for (tabB=0;tabB<FTK_n_svd_l*n_S*n_M;tabB++){
    for (tabA=0;tabA<n_w_max;tabA++){
      ulli = tabA + tabB*n_w_max;
      f0_svd_SVUXM_0in_wlSM____[ulli] = (float complex) (tabA + _Complex_I * (tabB%13));
      /* for (tabA=0;tabA<n_w_max;tabA++){ } */}
    /* for (tabB=0;tabB<FTK_n_svd_l*n_S*n_M;tabB++){ } */}
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,ulli_total,verbose," initialization: ");
  local_tic(0,t_start_,d_start_);
  mex_ampmh_X_wSM___fft_mkl_n1po_helper
    (
     verbose
     ,n_M
     ,n_M_use
     ,n_S
     ,FTK_n_svd_l
     ,n_w_max
     ,f0_svd_SVUXM_0in_wlSM____
     ,f0_svd_SVUXM_out_wlSM____
     );
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,ulli_total,verbose," mex_ampmh_X_wSM___fft_mkl_n1po_helper: ");
  rerror = cfnormn(ulli_numel,f0_svd_SVUXM_bkp_wlSM____,f0_svd_SVUXM_out_wlSM____);
  printf(" %% n1po relative error: %0.16f\n",rerror);
  if (verbose>1){
    array_sub_printf
      (
       f0_svd_SVUXM_0in_wlSM____
       ,"float complex"
       ,n_w_max
       ,4
       ,FTK_n_svd_l*n_S*n_M_use
       ,5
       ," %% f0_svd_SVUXM_0in_wlSM____: "
       );
    array_sub_printf
      (
       f0_svd_SVUXM_out_wlSM____
       ,"float complex"
       ,n_w_max
       ,4
       ,FTK_n_svd_l*n_S*n_M_use
       ,5
       ," %% f0_svd_SVUXM_out_wlSM____: "
       );
    array_printf_margin(f0_svd_SVUXM_out_wlSM____,"float complex",n_w_max,FTK_n_svd_l*n_S*n_M_use," %% f0_svd_SVUXM_out_wlSM____: ");
    /* if (verbose){ } */}
  _mm_free(f0_svd_SVUXM_0in_wlSM____); f0_svd_SVUXM_0in_wlSM____=NULL;
  _mm_free(f0_svd_SVUXM_out_wlSM____); f0_svd_SVUXM_out_wlSM____=NULL;
  /* %%%% */
  f0_svd_SVUXM_0in_wlSM____ = (float complex *) _mm_malloc(ulli_total*sizeof(float complex),64);
  memset(f0_svd_SVUXM_0in_wlSM____,0,ulli_total*sizeof(float complex));
  local_tic(0,t_start_,d_start_);
  for (tabB=0;tabB<FTK_n_svd_l*n_S*n_M;tabB++){
    for (tabA=0;tabA<n_w_max;tabA++){
      ulli = tabA + tabB*n_w_max;
      f0_svd_SVUXM_0in_wlSM____[ulli] = (float complex) (tabA + _Complex_I * (tabB%13));
      /* for (tabA=0;tabA<n_w_max;tabA++){ } */}
    /* for (tabB=0;tabB<FTK_n_svd_l*n_S*n_M;tabB++){ } */}
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,ulli_total,verbose," initialization: ");
  local_tic(0,t_start_,d_start_);
  mex_ampmh_X_wSM___fft_mkl_n1pi_helper
    (
     verbose
     ,n_M
     ,n_M_use
     ,n_S
     ,FTK_n_svd_l
     ,n_w_max
     ,f0_svd_SVUXM_0in_wlSM____
     );
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,ulli_total,verbose," mex_ampmh_X_wSM___fft_mkl_n1pi_helper: ");
  rerror = cfnormn(ulli_numel,f0_svd_SVUXM_bkp_wlSM____,f0_svd_SVUXM_0in_wlSM____);
  printf(" %% n1pi relative error: %0.16f\n",rerror);
  if (verbose>1){ array_printf_margin(f0_svd_SVUXM_0in_wlSM____,"float complex",n_w_max,FTK_n_svd_l*n_S*n_M_use," %% f0_svd_SVUXM_0in_wlSM____: "); /* if (verbose){ } */}
  _mm_free(f0_svd_SVUXM_0in_wlSM____); f0_svd_SVUXM_0in_wlSM____=NULL;
  /* %%%% */
  f0_svd_SVUXM_0in_wlSM____ = (float complex *) _mm_malloc(ulli_total*sizeof(float complex),64);
  f0_svd_SVUXM_out_wlSM____ = (float complex *) _mm_malloc(ulli_total*sizeof(float complex),64);
  memset(f0_svd_SVUXM_0in_wlSM____,0,ulli_total*sizeof(float complex));
  local_tic(0,t_start_,d_start_);
  for (tabB=0;tabB<FTK_n_svd_l*n_S*n_M;tabB++){
    for (tabA=0;tabA<n_w_max;tabA++){
      ulli = tabA + tabB*n_w_max;
      f0_svd_SVUXM_0in_wlSM____[ulli] = (float complex) (tabA + _Complex_I * (tabB%13));
      /* for (tabA=0;tabA<n_w_max;tabA++){ } */}
    /* for (tabB=0;tabB<FTK_n_svd_l*n_S*n_M;tabB++){ } */}
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,ulli_total,verbose," initialization: ");
  local_tic(0,t_start_,d_start_);
  mex_ampmh_X_wSM___fft_mkl_nxpo_helper
    (
     verbose
     ,n_M
     ,n_M_use
     ,n_S
     ,FTK_n_svd_l
     ,n_w_max
     ,f0_svd_SVUXM_0in_wlSM____
     ,f0_svd_SVUXM_out_wlSM____
     );
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,ulli_total,verbose," mex_ampmh_X_wSM___fft_mkl_nxpo_helper: ");
  rerror = cfnormn(ulli_numel,f0_svd_SVUXM_bkp_wlSM____,f0_svd_SVUXM_out_wlSM____);
  printf(" %% nxpo relative error: %0.16f\n",rerror);
  if (verbose>1){ array_printf_margin(f0_svd_SVUXM_out_wlSM____,"float complex",n_w_max,FTK_n_svd_l*n_S*n_M_use," %% f0_svd_SVUXM_out_wlSM____: "); /* if (verbose){ } */}
  _mm_free(f0_svd_SVUXM_0in_wlSM____); f0_svd_SVUXM_0in_wlSM____=NULL;
  _mm_free(f0_svd_SVUXM_out_wlSM____); f0_svd_SVUXM_out_wlSM____=NULL;
  /* %%%% */
  f0_svd_SVUXM_0in_wlSM____ = (float complex *) _mm_malloc(ulli_total*sizeof(float complex),64);
  f0_svd_SVUXM_out_wlSM____ = (float complex *) _mm_malloc(ulli_total*sizeof(float complex),64);
  memset(f0_svd_SVUXM_0in_wlSM____,0,ulli_total*sizeof(float complex));
  memset(f0_svd_SVUXM_out_wlSM____,0,ulli_total*sizeof(float complex));
  local_tic(0,t_start_,d_start_);
  for (tabB=0;tabB<FTK_n_svd_l*n_S*n_M;tabB++){
    for (tabA=0;tabA<n_w_max;tabA++){
      ulli = tabA + tabB*n_w_max;
      f0_svd_SVUXM_0in_wlSM____[ulli] = (float complex) (tabA + _Complex_I * (tabB%13));
      /* for (tabA=0;tabA<n_w_max;tabA++){ } */}
    /* for (tabB=0;tabB<FTK_n_svd_l*n_S*n_M;tabB++){ } */}
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,ulli_total,verbose," initialization: ");
  local_tic(0,t_start_,d_start_);
  n_M_per_Mbatch = n_M_per_Mbatch_use;
  n_Mbatch = ceil((double)n_M/(double)n_M_per_Mbatch);
  if (verbose>1){ printf(" %% n_Mbatch %d\n",n_Mbatch);}
#pragma omp parallel private(nMbatch)
  { /* begin omp parallel */
    nMbatch = 0;
#pragma omp for schedule(dynamic)
    for (nMbatch=0;nMbatch<n_Mbatch;nMbatch++){
      mex_ampmh_X_wSM___fft_mkl_nxpo_omp_helper
	(
	 0*verbose
	 ,n_M_per_Mbatch
	 ,n_M
	 ,n_S
	 ,FTK_n_svd_l
	 ,n_w_max
	 ,f0_svd_SVUXM_0in_wlSM____
	 ,f0_svd_SVUXM_out_wlSM____
	 ,nMbatch
	 ,n_Mbatch
	 );
      /* for (nMbatch=0;nMbatch<n_Mbatch;nMbatch++){ } */}
    /* end omp parallel */}
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,ulli_total,verbose," mex_ampmh_X_wSM___fft_mkl_nxpo_omp_helper: ");
  rerror = cfnormn(ulli_numel,f0_svd_SVUXM_bkp_wlSM____,f0_svd_SVUXM_out_wlSM____);
  printf(" %% nxpo_omp relative error: %0.16f\n",rerror);
  if (verbose>1){ array_printf_margin(f0_svd_SVUXM_out_wlSM____,"float complex",n_w_max,FTK_n_svd_l*n_S*n_M_use," %% f0_svd_SVUXM_out_wlSM____: "); /* if (verbose){ } */}
  _mm_free(f0_svd_SVUXM_0in_wlSM____); f0_svd_SVUXM_0in_wlSM____=NULL;
  _mm_free(f0_svd_SVUXM_out_wlSM____); f0_svd_SVUXM_out_wlSM____=NULL;
  /* %%%% */
  f0_svd_SVUXM_0in_wlSM____ = (float complex *) _mm_malloc(ulli_total*sizeof(float complex),64);
  memset(f0_svd_SVUXM_0in_wlSM____,0,ulli_total*sizeof(float complex));
  local_tic(0,t_start_,d_start_);
  for (tabB=0;tabB<FTK_n_svd_l*n_S*n_M;tabB++){
    for (tabA=0;tabA<n_w_max;tabA++){
      ulli = tabA + tabB*n_w_max;
      f0_svd_SVUXM_0in_wlSM____[ulli] = (float complex) (tabA + _Complex_I * (tabB%13));
      /* for (tabA=0;tabA<n_w_max;tabA++){ } */}
    /* for (tabB=0;tabB<FTK_n_svd_l*n_S*n_M;tabB++){ } */}
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,ulli_total,verbose," initialization: ");
  local_tic(0,t_start_,d_start_);
  mex_ampmh_X_wSM___fft_mkl_nxpi_helper
    (
     verbose
     ,n_M
     ,n_M_use
     ,n_S
     ,FTK_n_svd_l
     ,n_w_max
     ,f0_svd_SVUXM_0in_wlSM____
     );
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,ulli_total,verbose," mex_ampmh_X_wSM___fft_mkl_nxpi_helper: ");
  rerror = cfnormn(ulli_numel,f0_svd_SVUXM_bkp_wlSM____,f0_svd_SVUXM_0in_wlSM____);
  printf(" %% nxpi relative error: %0.16f\n",rerror);
  if (verbose>1){ array_printf_margin(f0_svd_SVUXM_0in_wlSM____,"float complex",n_w_max,FTK_n_svd_l*n_S*n_M_use," %% f0_svd_SVUXM_0in_wlSM____: "); /* if (verbose){ } */}
  _mm_free(f0_svd_SVUXM_0in_wlSM____); f0_svd_SVUXM_0in_wlSM____=NULL;
  /* %%%% */
  f0_svd_SVUXM_0in_wlSM____ = (float complex *) _mm_malloc(ulli_total*sizeof(float complex),64);
  memset(f0_svd_SVUXM_0in_wlSM____,0,ulli_total*sizeof(float complex));
  local_tic(0,t_start_,d_start_);
  for (tabB=0;tabB<FTK_n_svd_l*n_S*n_M;tabB++){
    for (tabA=0;tabA<n_w_max;tabA++){
      ulli = tabA + tabB*n_w_max;
      f0_svd_SVUXM_0in_wlSM____[ulli] = (float complex) (tabA + _Complex_I * (tabB%13));
      /* for (tabA=0;tabA<n_w_max;tabA++){ } */}
    /* for (tabB=0;tabB<FTK_n_svd_l*n_S*n_M;tabB++){ } */}
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,ulli_total,verbose," initialization: ");
  local_tic(0,t_start_,d_start_);
  n_M_per_Mbatch = n_M_per_Mbatch_use;
  n_Mbatch = ceil((double)n_M/(double)n_M_per_Mbatch);
  if (verbose>1){ printf(" %% n_Mbatch %d\n",n_Mbatch);}
#pragma omp parallel private(nMbatch)
  { /* begin omp parallel */
    nMbatch = 0;
#pragma omp for schedule(dynamic)
    for (nMbatch=0;nMbatch<n_Mbatch;nMbatch++){
      mex_ampmh_X_wSM___fft_mkl_nxpi_omp_helper
	(
	 0*verbose
	 ,n_M_per_Mbatch
	 ,n_M
	 ,n_S
	 ,FTK_n_svd_l
	 ,n_w_max
	 ,f0_svd_SVUXM_0in_wlSM____
	 ,nMbatch
	 ,n_Mbatch
	 );
      /* for (nMbatch=0;nMbatch<n_Mbatch;nMbatch++){ } */}
    /* end omp parallel */}
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,ulli_total,verbose," mex_ampmh_X_wSM___fft_mkl_nxpi_omp_helper: ");
  rerror = cfnormn(ulli_numel,f0_svd_SVUXM_bkp_wlSM____,f0_svd_SVUXM_0in_wlSM____);
  printf(" %% nxpi_omp relative error: %0.16f\n",rerror);
  if (verbose>1){ array_printf_margin(f0_svd_SVUXM_0in_wlSM____,"float complex",n_w_max,FTK_n_svd_l*n_S*n_M_use," %% f0_svd_SVUXM_0in_wlSM____: "); /* if (verbose){ } */}
  _mm_free(f0_svd_SVUXM_0in_wlSM____); f0_svd_SVUXM_0in_wlSM____=NULL;
  /* %%%%%%%%%%%%%%%% */
  _mm_free(f1_svd_SVUXM_0inout_wlSM____); f1_svd_SVUXM_0inout_wlSM____=NULL;
  fftw_destroy_plan(fftw_plan_guru_split_dft_plan);
  _mm_free(f0_svd_SVUXM_bkp_wlSM____); f0_svd_SVUXM_bkp_wlSM____=NULL;
  if (verbose){ printf(" %% [finished mex_ampmh_X_wSM___fft_mkl_test]\n");}
}


