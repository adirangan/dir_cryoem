/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* The helper function */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
void mex_ampmh_X_wSM___15_omp_helper
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
  /* %%%%; */
 ,int verbose
 ,int flag_slow
 ,int flag_fast
  /* %%%%; */
 ,double complex *FTK_svd_U_d_expiw_s__
 ,float *f_FTK_svd_U_d_expiw_s_dl_real__
 ,float *f_FTK_svd_U_d_expiw_s_dl_imag__
 ,float *f_FTK_svd_U_d_expiw_s_ld_real__
 ,float *f_FTK_svd_U_d_expiw_s_ld_imag__
 ,double complex *CTF_UX_S_k_q_wnS__
 ,double complex *CTF_UX_S_k_q_nSw___
 ,float *f_CTF_UX_S_k_q_wnS_real___
 ,float *f_CTF_UX_S_k_q_wnS_imag___
 ,float *f_CTF_UX_S_k_q_nSw_real___
 ,float *f_CTF_UX_S_k_q_nSw_imag___
 ,double complex *svd_VUXM_lwnM____
 ,float *f_svd_VUXM_lwnM_real____
 ,float *f_svd_VUXM_lwnM_imag____
 ,double *gamma_z_
  /* %%%%; */
 ,int nMbatch
 ,int n_Mbatch
 ,fftw_plan fftw_plan_guru_split_dft_plan
 ,fftwf_plan fftwf_plan_guru_split_dft_plan
 )
{
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
  float *f_svd_VUXM_nMlw_real____=NULL;
  float *f_svd_VUXM_nMlw_imag____=NULL;
  unsigned long long int na=0,n_a=0;
  /* %%%% */
  int nSbatch=0,n_Sbatch=0;
  int *index_M_in_Mbatch_=NULL;
  int *index_S_in_Sbatch_=NULL;
  int n_M_sub=0,nM_sub=0;
  int n_S_sub=0,nS_sub=0;
  int nl=0,ndelta_v=0,nw=0,ndw=0,nS=0,nM=0;
  int pm_nUX_rank=0;
  unsigned long long int tab=0,tabA=0,tabB=0,tabC=0;
  double complex *svd_VUXM_nMwl____=NULL;
  double complex *svd_VUXM_nMlw____=NULL;
  double complex *svd_SVUXM_SMwl____=NULL;
  double complex *svd_SVUXM_SMlw____=NULL;
  float *f_svd_SVUXM_SMlw_real____=NULL;
  float *f_svd_SVUXM_SMlw_imag____=NULL;
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
  float *f_CR_=NULL,*f_CI_=NULL,*f_C_=NULL;
  double d_l2_R=0,d_l2_I=0,d_l2=0,d_l2_iso=0;
  /* %%%% */
  fftw_complex * svd_SVUXM_0inout_wlSM____=NULL;
  fftw_complex * svd_SVUXM_0in_wlSM____=NULL;
  fftw_complex * svd_SVUXM_out_wlSM____=NULL;
  /* %%%% */
  float * f_svd_SVUXM_0inout_wSMl_realimag____=NULL;
  float * f_svd_SVUXM_0in_wSMl_real____=NULL;
  float * f_svd_SVUXM_0in_wSMl_imag____=NULL;
  float * f_svd_SVUXM_out_wSMl_real____=NULL;
  float * f_svd_SVUXM_out_wSMl_imag____=NULL;
  /* %%%% */
  if (verbose){ printf(" %% [entering mex_ampmh_X_wSM___15_omp_helper]\n");}
  /* %%%%%%%%%%%%%%%% */
  FTK_n_delta_v_rup = rup(FTK_n_delta_v,8); FTK_n_delta_v_256 = FTK_n_delta_v_rup/8;
  FTK_n_svd_l_rup = rup(FTK_n_svd_l,8); FTK_n_svd_l_256 = FTK_n_svd_l_rup/8;
  n_w_max_rup = rup(n_w_max,8); n_w_max_256 = n_w_max_rup/8;
  pm_n_UX_rank_rup = rup(pm_n_UX_rank,8); pm_n_UX_rank_256 = pm_n_UX_rank_rup/8;
  n_S_rup = rup(n_S,8); n_S_256 = n_S_rup/8;
  n_S_per_Sbatch_rup = rup(n_S_per_Sbatch,8); n_S_per_Sbatch_256 = n_S_per_Sbatch_rup/8;
  n_M_rup = rup(n_M,8); n_M_256 = n_M_rup/8;
  n_M_per_Mbatch_rup = rup(n_M_per_Mbatch,8); n_M_per_Mbatch_256 = n_M_per_Mbatch_rup/8;
  /* %%%%; */
  if (flag_slow){
    tab = 
      (unsigned long long int)pm_n_UX_rank
      *(unsigned long long int)n_M_per_Mbatch
      *(unsigned long long int)n_w_max
      *(unsigned long long int)FTK_n_svd_l;
    svd_VUXM_nMwl____ = (double complex *) malloc(tab*sizeof(double complex));
    svd_VUXM_nMlw____ = (double complex *) malloc(tab*sizeof(double complex));
    tab = 
      (unsigned long long int)n_S
      *(unsigned long long int)n_M_per_Mbatch
      *(unsigned long long int)n_w_max
      *(unsigned long long int)FTK_n_svd_l;
    svd_SVUXM_SMwl____ = (double complex *) malloc(tab*sizeof(double complex));
    svd_SVUXM_SMlw____ = (double complex *) malloc(tab*sizeof(double complex));
    /* if (flag_slow){ } */}
  if (flag_fast){
    f_svd_VUXM_nMlw_real____ = float____m256_malloc_from_double____(pm_n_UX_rank,n_M_per_Mbatch,FTK_n_svd_l,n_w_max,NULL);
    f_svd_VUXM_nMlw_imag____ = float____m256_malloc_from_double____(pm_n_UX_rank,n_M_per_Mbatch,FTK_n_svd_l,n_w_max,NULL);
    f_svd_SVUXM_SMlw_real____ = float____m256_malloc_from_double____(FTK_n_svd_l,n_w_max,n_S,n_M_per_Mbatch,NULL);
    f_svd_SVUXM_SMlw_imag____ = float____m256_malloc_from_double____(FTK_n_svd_l,n_w_max,n_S,n_M_per_Mbatch,NULL);
    /* if (flag_fast){ } */}
  /* %%%% */
  if (flag_slow){
    tab = 
      (unsigned long long int)n_w_max_rup
      *(unsigned long long int)FTK_n_svd_l_rup
      *(unsigned long long int)n_S_rup
      *(unsigned long long int)n_M_per_Mbatch_rup;
    svd_SVUXM_0inout_wlSM____ = (fftw_complex *) _mm_malloc(tab/1*sizeof(__m256),32); //%<-- 8 floats per __m256 --> 2 double complex per __m256 ;
    svd_SVUXM_0in_wlSM____ = 0*tab + (fftw_complex *) svd_SVUXM_0inout_wlSM____;
    svd_SVUXM_out_wlSM____ = 1*tab + (fftw_complex *) svd_SVUXM_0inout_wlSM____;
    /* if (flag_slow){ } */}
  /* %%%% */
  if (flag_fast){
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
    /* if (flag_fast){ } */}
  /* %%%% */
  if (flag_slow){
    tab = 
      (unsigned long long int)FTK_n_svd_l
      *(unsigned long long int)n_w_max
      *(unsigned long long int)n_S
      *(unsigned long long int)n_M_per_Mbatch;
    svd_SVUXM_lwSM____ = (double complex *) malloc(tab*sizeof(double complex));
    /* if (flag_slow){ } */}
  if (flag_fast){
    f_svd_SVUXM_lwSM_real____ = float____m256_malloc_from_double____(FTK_n_svd_l,n_w_max,n_S,n_M_per_Mbatch,NULL);
    f_svd_SVUXM_lwSM_imag____ = float____m256_malloc_from_double____(FTK_n_svd_l,n_w_max,n_S,n_M_per_Mbatch,NULL);
    /* if (flag_fast){ } */}
  if (flag_slow){
    tab = 
      (unsigned long long int)n_S_per_Sbatch
      *(unsigned long long int)n_M_per_Mbatch
      *(unsigned long long int)n_w_max
      *(unsigned long long int)FTK_n_svd_l;
    svd_SVUXM_lwsM____ = (double complex *) malloc(tab*sizeof(double complex));
    /* if (flag_slow){ } */}
  if (flag_fast){
    f_svd_SVUXM_lwsM_real____ = float____m256_malloc_from_double____(FTK_n_svd_l,n_w_max,n_S_per_Sbatch,n_M_per_Mbatch,NULL);
    f_svd_SVUXM_lwsM_imag____ = float____m256_malloc_from_double____(FTK_n_svd_l,n_w_max,n_S_per_Sbatch,n_M_per_Mbatch,NULL);
    /* if (flag_fast){ } */}
  if (flag_slow){
    tab = 
      (unsigned long long int)n_S_per_Sbatch
      *(unsigned long long int)n_M_per_Mbatch
      *(unsigned long long int)n_w_max
      *(unsigned long long int)FTK_n_delta_v;
    svd_USESVUXM_dwSM____ = (double complex *) malloc(tab*sizeof(double complex));
    /* if (flag_slow){ } */}
  if (flag_fast){
    f_svd_USESVUXM_dwSM_real____ = float____m256_malloc_from_double____(FTK_n_delta_v,n_w_max,n_S_per_Sbatch,n_M_per_Mbatch,NULL);
    f_svd_USESVUXM_dwSM_imag____ = float____m256_malloc_from_double____(FTK_n_delta_v,n_w_max,n_S_per_Sbatch,n_M_per_Mbatch,NULL);
    /* if (flag_fast){ } */}
  tab = (unsigned long long int)n_S_per_Sbatch*(unsigned long long int)n_M_per_Mbatch*(unsigned long long int)FTK_n_delta_v;
  l2_dSM___ = (double *) malloc(tab*sizeof(double));
  n2_dSM___ = (double *) malloc(tab*sizeof(double));
  f2_dSM___ = (double *) malloc(tab*sizeof(double));
  tab = 
    (unsigned long long int)n_S_per_Sbatch
    *(unsigned long long int)n_M_per_Mbatch
    *(unsigned long long int)n_w_max
    *(unsigned long long int)FTK_n_delta_v;
  X_sub_dwSM____ = (double *) malloc(tab*sizeof(double));
  d_X_sub_dwSM____ = (double *) malloc(tab*sizeof(double));
  if (flag_compute_I_value){ I_value_use_dwSM____ = (double *) malloc(tab*sizeof(double));}
  /* %%%%%%%%%%%%%%%% */
  n_Mbatch = ceil((double)n_M/(double)n_M_per_Mbatch);
  if (verbose>1){ printf(" %% n_Mbatch %d\n",n_Mbatch);}
  index_M_in_Mbatch_ = (int *) malloc(n_M_per_Mbatch*sizeof(int));
  n_Sbatch = ceil((double)n_S/(double)n_S_per_Sbatch);
  if (verbose>1){ printf(" %% n_Sbatch %d\n",n_Sbatch);}
  index_S_in_Sbatch_ = (int *) malloc(n_S_per_Sbatch*sizeof(int));
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
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
      tab =
  	(unsigned long long int)FTK_n_svd_l
  	*(unsigned long long int)n_w_max
  	*(unsigned long long int)n_M_sub
  	*(unsigned long long int)pm_n_UX_rank;
      na =
  	(unsigned long long int)FTK_n_svd_l
  	*(unsigned long long int)n_w_max
  	*(unsigned long long int)pm_n_UX_rank
  	*(unsigned long long int)(n_M_per_Mbatch*nMbatch);
      for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
  	for (pm_nUX_rank=0;pm_nUX_rank<pm_n_UX_rank;pm_nUX_rank++){
  	  for (nw=0;nw<n_w_max;nw++){
  	    for (nl=0;nl<FTK_n_svd_l;nl++){
  	      /* tabA = */
  	      /* 	(unsigned long long int)pm_nUX_rank */
  	      /* 	+ (unsigned long long int)nM_sub*(unsigned long long int)pm_n_UX_rank */
  	      /* 	+ (unsigned long long int)nw*(unsigned long long int)n_M_sub*(unsigned long long int)pm_n_UX_rank */
  	      /* 	+ (unsigned long long int)nl*(unsigned long long int)n_w_max*(unsigned long long int)n_M_sub*(unsigned long long int)pm_n_UX_rank; */
  	      tabA =
  		(unsigned long long int)pm_nUX_rank +
  		((unsigned long long int)nM_sub +
  		 ((unsigned long long int)nw +
  		  (unsigned long long int)nl
  		  *(unsigned long long int)n_w_max)
  		 *(unsigned long long int)n_M_sub)
  		*(unsigned long long int)pm_n_UX_rank;
  	      tabB =
  		(unsigned long long int)pm_nUX_rank +
  		((unsigned long long int)nM_sub +
  		 ((unsigned long long int)nl +
  		  (unsigned long long int)nw
  		  *(unsigned long long int)FTK_n_svd_l)
  		 *(unsigned long long int)n_M_sub)
  		*(unsigned long long int)pm_n_UX_rank;
  	      svd_VUXM_nMwl____[tabA] = svd_VUXM_lwnM____[na];
  	      svd_VUXM_nMlw____[tabB] = svd_VUXM_lwnM____[na];
  	      na++;
  	      /* for (nl=0;nl<FTK_n_svd_l;nl++){ } */}
  	    /* for (nw=0;nw<n_w_max;nw++){ } */}
  	  /* for (pm_nUX_rank=0;pm_nUX_rank<pm_n_UX_rank;pm_nUX_rank++){ } */}
  	/* for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ } */}
      /* if (flag_slow){ } */}
    /* %%%% */
    if (flag_fast){
      tab =
  	(unsigned long long int)FTK_n_svd_l
  	*(unsigned long long int)n_w_max
  	*(unsigned long long int)pm_n_UX_rank
  	*(unsigned long long int)n_M_sub;
      tabA =
  	(unsigned long long int)FTK_n_svd_l_rup
  	*(unsigned long long int)n_w_max_rup
  	*(unsigned long long int)pm_n_UX_rank_rup
  	*(unsigned long long int)n_M_sub_rup;
      memset(f_svd_VUXM_nMlw_real____,0,tabA*sizeof(float));
      memset(f_svd_VUXM_nMlw_imag____,0,tabA*sizeof(float));
      tabB =
  	(unsigned long long int)FTK_n_svd_l_rup
  	*(unsigned long long int)n_w_max_rup
  	*(unsigned long long int)pm_n_UX_rank_rup
  	*(unsigned long long int)(n_M_per_Mbatch*nMbatch);
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
      /* if (flag_fast){ } */}
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
      /* if (flag_slow){ } */}
    /* %%%% */
    if (flag_slow){
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
      /* if (flag_slow){ } */}
    /* %%%% */
    if (flag_fast){
      tab = (unsigned long long int)(FTK_n_svd_l*n_w_max)*(unsigned long long int)(n_S*n_M_sub)*(unsigned long long int)pm_n_UX_rank;
      tabA = (unsigned long long int)(FTK_n_svd_l_rup*n_w_max_rup)*(unsigned long long int)(n_S_rup*n_M_sub_rup);
      memset(f_svd_SVUXM_SMlw_real____,0,tabA*sizeof(float));
      memset(f_svd_SVUXM_SMlw_imag____,0,tabA*sizeof(float));
      for (nl=0;nl<FTK_n_svd_l;nl++){
  	for (nw=0;nw<n_w_max;nw++){
  	  tabA =
  	    (unsigned long long int)nw
  	    *(unsigned long long int)pm_n_UX_rank_rup
  	    *(unsigned long long int)n_S_rup;
  	  tabB =
  	    (unsigned long long int)(nl+nw*FTK_n_svd_l_rup)
  	    *(unsigned long long int)pm_n_UX_rank_rup
  	    *(unsigned long long int)n_M_sub_rup;
  	  tabC =
  	    (unsigned long long int)(nl+nw*FTK_n_svd_l_rup)
  	    *(unsigned long long int)n_S_rup
  	    *(unsigned long long int)n_M_sub_rup;
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
      /* if (flag_fast){ } */}
    /* %%%% */
    /* svd_SVUXM_lwSM____ = permute(ifft(permute(svd_SVUXM_SMwl____,[3,4,1,2]),[],1)*n_w_max,[2,1,3,4]); */
    /* %%%% */
    if (flag_slow){
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
  	      tabA =
  		(unsigned long long int)nw +
  		((unsigned long long int)nl +
  		 ((unsigned long long int)nS +
  		  (unsigned long long int)nM_sub
  		  *(unsigned long long int)n_S)
  		 *(unsigned long long int)FTK_n_svd_l)
  		*(unsigned long long int)n_w_max;
  	      svd_SVUXM_0in_wlSM____[tabA] = (fftw_complex) (svd_SVUXM_SMwl____[na]); na++;
  	      /* for (nS=0;nS<n_S;nS++){ } */}
  	    /* for (nM_sub=0;nM_sub<n_M_sub;nS++){ } */}
  	  /* for (nw=0;nw<n_w_max;nw++){ } */}
  	/* for (nl=0;nl<FTK_n_svd_l;nl++){ } */}
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
  	      tabA =
  		(unsigned long long int)nl +
  		((unsigned long long int)nw +
  		 ((unsigned long long int)nS +
  		  (unsigned long long int)nM_sub
  		  *(unsigned long long int)n_S)
  		 *(unsigned long long int)n_w_max)
  		*(unsigned long long int)FTK_n_svd_l;
  	      svd_SVUXM_lwSM____[tabA] = (double complex) (svd_SVUXM_out_wlSM____[na]); na++;
  	      /* for (nw=0;nw<n_w_max;nw++){ } */}
  	    /* for (nl=0;nl<FTK_n_svd_l;nl++){ } */}
  	  /* for (nS=0;nS<n_S;nS++){ } */}
  	/* for (nM_sub=0;nM_sub<n_M_sub;nS++){ } */}
      /* if (flag_slow){ } */}
    /* %%%% */
    if (flag_fast){
      tab =
  	(unsigned long long int)n_w_max
  	*(unsigned long long int)n_S
  	*(unsigned long long int)n_M_sub
  	*(unsigned long long int)FTK_n_svd_l;
      tabA =
  	(unsigned long long int)n_w_max_rup
  	*(unsigned long long int)n_S_rup
  	*(unsigned long long int)n_M_sub_rup
  	*(unsigned long long int)FTK_n_svd_l_rup;
      memset(f_svd_SVUXM_0in_wSMl_real____,0,tabA*sizeof(float));
      memset(f_svd_SVUXM_0in_wSMl_imag____,0,tabA*sizeof(float));
      tabB =
  	(unsigned long long int)n_S_rup
  	*(unsigned long long int)n_M_sub_rup
  	*(unsigned long long int)FTK_n_svd_l_rup;
      tabC =
  	(unsigned long long int)n_w_max_rup;
      transpose_ps_block1(
  			  f_svd_SVUXM_SMlw_real____
  			  ,f_svd_SVUXM_0in_wSMl_real____
  			  ,tabB
  			  ,tabC
  			  ,tabB
  			  ,tabC
  			  ,16
  			  );
      transpose_ps_block1(
  			  f_svd_SVUXM_SMlw_imag____
  			  ,f_svd_SVUXM_0in_wSMl_imag____
  			  ,tabB
  			  ,tabC
  			  ,tabB
  			  ,tabC
  			  ,16
  			  );
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
      tab =
  	(unsigned long long int)FTK_n_svd_l
  	*(unsigned long long int)n_w_max
  	*(unsigned long long int)n_S
  	*(unsigned long long int)n_M_sub;
      tabA =
  	(unsigned long long int)FTK_n_svd_l_rup
  	*(unsigned long long int)n_w_max_rup
  	*(unsigned long long int)n_S_rup
  	*(unsigned long long int)n_M_sub_rup;
      memset(f_svd_SVUXM_lwSM_real____,0,tabA*sizeof(float));
      memset(f_svd_SVUXM_lwSM_imag____,0,tabA*sizeof(float));
      tabB =
  	(unsigned long long int)n_w_max_rup
  	*(unsigned long long int)n_S_rup
  	*(unsigned long long int)n_M_sub_rup;
      tabC =
  	(unsigned long long int)FTK_n_svd_l_rup;
      transpose_ps_block1(
  			  f_svd_SVUXM_out_wSMl_real____
  			  ,f_svd_SVUXM_lwSM_real____
  			  ,tabB
  			  ,tabC
  			  ,tabB
  			  ,tabC
  			  ,16
  			  );
      transpose_ps_block1(
  			  f_svd_SVUXM_out_wSMl_imag____
  			  ,f_svd_SVUXM_lwSM_imag____
  			  ,tabB
  			  ,tabC
  			  ,tabB
  			  ,tabC
  			  ,16
  			  );
      /* if (flag_fast){ } */}
    /* %%%% */
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
  	  tab = (unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub*(unsigned long long int)n_w_max*(unsigned long long int)FTK_n_svd_l;
  	  memset(svd_SVUXM_lwsM____,0,tab*sizeof(double complex));
  	  for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
  	    for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){
  	      nS = index_S_in_Sbatch_[nS_sub];
  	      tabA =
  		(((unsigned long long int)nS_sub +
  		  (unsigned long long int)nM_sub
  		  *(unsigned long long int)n_S_sub)
  		 *(unsigned long long int)n_w_max)
  		*(unsigned long long int)FTK_n_svd_l;
  	      tabB =
  		(((unsigned long long int)nS     +
  		  (unsigned long long int)nM_sub
  		  *(unsigned long long int)n_S    )
  		 *(unsigned long long int)n_w_max)
  		*(unsigned long long int)FTK_n_svd_l;
  	      tabC = (unsigned long long int)n_w_max*(unsigned long long int)FTK_n_svd_l;
  	      memcpy( svd_SVUXM_lwsM____ + tabA , svd_SVUXM_lwSM____ + tabB , tabC*sizeof(double complex));
  	      /* for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){ } */}
  	    /* for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ } */}
  	  /* if (flag_slow){ } */}
  	/* %%%% */
  	if (flag_fast){
  	  tab =
  	    (unsigned long long int)n_S_sub
  	    *(unsigned long long int)n_M_sub
  	    *(unsigned long long int)n_w_max
  	    *(unsigned long long int)FTK_n_svd_l;
  	  tabA =
  	    (unsigned long long int)n_S_sub_rup
  	    *(unsigned long long int)n_M_sub_rup
  	    *(unsigned long long int)n_w_max_rup
  	    *(unsigned long long int)FTK_n_svd_l_rup;
  	  memset(f_svd_SVUXM_lwsM_real____,12,tabA*sizeof(float));
  	  memset(f_svd_SVUXM_lwsM_imag____,13,tabA*sizeof(float));
  	  for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
  	    for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){
  	      nS = index_S_in_Sbatch_[nS_sub];
  	      tabA =
  		(((unsigned long long int)nS_sub +
  		  (unsigned long long int)nM_sub
  		  *(unsigned long long int)n_S_sub_rup)
  		 *(unsigned long long int)n_w_max_rup)
  		*(unsigned long long int)FTK_n_svd_l_rup;
  	      tabB =
  		(((unsigned long long int)nS     +
  		  (unsigned long long int)nM_sub
  		  *(unsigned long long int)n_S_rup    )
  		 *(unsigned long long int)n_w_max_rup)
  		*(unsigned long long int)FTK_n_svd_l_rup;
  	      tabC =
  		(unsigned long long int)n_w_max_rup
  		*(unsigned long long int)FTK_n_svd_l_rup;
  	      memcpy( f_svd_SVUXM_lwsM_real____ + tabA , f_svd_SVUXM_lwSM_real____ + tabB , tabC*sizeof(float));
  	      memcpy( f_svd_SVUXM_lwsM_imag____ + tabA , f_svd_SVUXM_lwSM_imag____ + tabB , tabC*sizeof(float));
  	      /* for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){ } */}
  	    /* for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ } */}
  	  /* if (flag_fast){ } */}
  	/* %%%% */
  	/* svd_USESVUXM_dwSM____ = real(reshape(FTK.svd_U_d_expiw_s__*reshape(svd_SVUXM_lwsM____,[FTK.n_svd_l,n_w_max*n_S_sub*n_M_sub]),[FTK.n_delta_v,n_w_max,n_S_sub,n_M_sub])); */
  	/* %%%% */
  	if (flag_slow){
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
  	  /* if (flag_slow){ } */}
  	/* %%%% */
  	if (flag_fast){
  	  tab =
  	    (unsigned long long int)n_S_sub
  	    *(unsigned long long int)n_M_sub
  	    *(unsigned long long int)n_w_max
  	    *(unsigned long long int)FTK_n_delta_v
  	    *(unsigned long long int)FTK_n_svd_l;
  	  tabA =
  	    (unsigned long long int)n_S_sub_rup
  	    *(unsigned long long int)n_M_sub_rup
  	    *(unsigned long long int)n_w_max_rup
  	    *(unsigned long long int)FTK_n_delta_v_rup;
  	  memset(f_svd_USESVUXM_dwSM_real____,0,tabA*sizeof(float));
  	  memset(f_svd_USESVUXM_dwSM_imag____,0,tabA*sizeof(float));
  	  tabC = (unsigned long long int)n_S_sub_rup*(unsigned long long int)n_M_sub_rup*(unsigned long long int)n_w_max_rup;
  	  nhpr_segregated_to_segregated_mult_immintrin_load1_fma(
  								 FTK_n_delta_v_rup
  								 ,FTK_n_svd_l_rup
  								 ,f_FTK_svd_U_d_expiw_s_ld_real__
  								 ,f_FTK_svd_U_d_expiw_s_ld_imag__
  								 ,tabC
  								 ,f_svd_SVUXM_lwsM_real____
  								 ,f_svd_SVUXM_lwsM_imag____
  								 ,&f_svd_USESVUXM_dwSM_real____
  								 ,&f_svd_USESVUXM_dwSM_imag____
  								 );
  	  /* if (flag_fast){ } */}
  	/* %%%% */
  	/* %%%% */
  	/*
  	   l2_dSM___ = permute(reshape(reshape(sqrt(CTF_UX_S_l2_use_(1+index_S_in_Sbatch_)),[n_S_sub,1])*reshape(sqrt(UX_M_l2_dM_use__(:,1+index_M_in_Mbatch_)),[1,FTK.n_delta_v*n_M_sub]),[n_S_sub,FTK.n_delta_v,n_M_sub]),[2,1,3]);
  	   n2_dSM___ = 1./max(1e-14,l2_dSM___);
  	   f2_dSM___ = permute(reshape(reshape(sqrt(CTF_UX_S_l2_use_(1+index_S_in_Sbatch_)),[n_S_sub,1])*reshape(1./max(1e-14,sqrt(UX_M_l2_dM_use__(:,1+index_M_in_Mbatch_))),[1,FTK.n_delta_v*n_M_sub]),[n_S_sub,FTK.n_delta_v,n_M_sub]),[2,1,3]);
  	   ss_S_ = reshape(CTF_UX_S_l2_use_(1+index_S_in_Sbatch_),[n_S_sub,1]);
  	*/
  	/* %%%% */
  	tab =
  	  (unsigned long long int)n_S_sub
  	  *(unsigned long long int)n_M_sub
  	  *(unsigned long long int)FTK_n_delta_v;
  	memset(l2_dSM___,0,tab*sizeof(double));
  	memset(n2_dSM___,0,tab*sizeof(double));
  	memset(f2_dSM___,0,tab*sizeof(double));
  	for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){
  	  nS = index_S_in_Sbatch_[nS_sub];
  	  for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
  	    nM = index_M_in_Mbatch_[nM_sub];
  	    for (ndelta_v=0;ndelta_v<FTK_n_delta_v;ndelta_v++){
  	      tabA =
  		(unsigned long long int)ndelta_v +
  		(unsigned long long int)nM
  		*(unsigned long long int)FTK_n_delta_v;
  	      tabB =
  		(unsigned long long int)ndelta_v +
  		((unsigned long long int)nS_sub +
  		 (unsigned long long int)nM_sub
  		 *(unsigned long long int)n_S_sub)
  		*(unsigned long long int)FTK_n_delta_v;
  	      l2_dSM___[tabB] = sqrt(CTF_UX_S_l2_[nS]) * sqrt(UX_M_l2_dM__[tabA]);
  	      n2_dSM___[tabB] = 1.0 / maximum(1.0e-14,l2_dSM___[tabB]);
  	      f2_dSM___[tabB] = sqrt(CTF_UX_S_l2_[nS]) / maximum(1.0e-14,sqrt(UX_M_l2_dM__[tabA]));
  	      /* for (ndelta_v=0;ndelta_v<FTK_n_delta_v;ndelta_v++){ } */}
  	    /* for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ } */}
  	  /* for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){ } */}
  	/* %%%% */
  	/* X_sub_dwSM____ = repmat(reshape(n2_dSM___,[FTK.n_delta_v,1,n_S_sub,n_M_sub]),[1,n_w_max,1,1]) .* svd_USESVUXM_dwSM____; %<-- correlation. ; */
  	/* %%%% */
  	if (flag_slow){
  	  tab =
  	    (unsigned long long int)n_S_sub
  	    *(unsigned long long int)n_M_sub
  	    *(unsigned long long int)n_w_max
  	    *(unsigned long long int)FTK_n_delta_v;
  	  memset(X_sub_dwSM____,0,tab*sizeof(double));
  	  na=0;
  	  for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
  	    for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){
  	      for (nw=0;nw<n_w_max;nw++){
  		for (ndelta_v=0;ndelta_v<FTK_n_delta_v;ndelta_v++){
  		  tabA =
  		    (unsigned long long int)ndelta_v +
  		    ((unsigned long long int)nS_sub +
  		     (unsigned long long int)nM_sub
  		     *(unsigned long long int)n_S_sub)
  		    *(unsigned long long int)FTK_n_delta_v;
  		  X_sub_dwSM____[na] = n2_dSM___[tabA] * creal(svd_USESVUXM_dwSM____[na]) ;
  		  na += 1;
  		  /* for (ndelta_v=0;ndelta_v<FTK_n_delta_v;ndelta_v++){ } */}
  		/* for (nw=0;nw<n_w_max;nw++){ } */}
  	      /* for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){ } */}
  	    /* for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ } */}
  	  /* if (flag_slow){ } */}
  	/* %%%% */
  	if (flag_fast){
  	  tab =
  	    (unsigned long long int)n_S_sub
  	    *(unsigned long long int)n_M_sub
  	    *(unsigned long long int)n_w_max
  	    *(unsigned long long int)FTK_n_delta_v;
  	  memset(d_X_sub_dwSM____,0,tab*sizeof(double));
  	  na=0;
  	  for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){
  	    for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){
  	      for (nw=0;nw<n_w_max;nw++){
  		for (ndelta_v=0;ndelta_v<FTK_n_delta_v;ndelta_v++){
  		  tabA =
  		    (unsigned long long int)ndelta_v +
  		    ((unsigned long long int)nS_sub +
  		     (unsigned long long int)nM_sub
  		     *(unsigned long long int)n_S_sub)
  		    *(unsigned long long int)FTK_n_delta_v;
  		  tabB =
  		    (unsigned long long int)ndelta_v +
  		    ((unsigned long long int)nw +
  		     ((unsigned long long int)nS_sub +
  		      (unsigned long long int)nM_sub
  		      *(unsigned long long int)n_S_sub_rup)
  		     *(unsigned long long int)n_w_max_rup)
  		    *(unsigned long long int)FTK_n_delta_v_rup;
  		  d_X_sub_dwSM____[na] = n2_dSM___[tabA] * f_svd_USESVUXM_dwSM_real____[tabB] ;
  		  na += 1;
  		  /* for (ndelta_v=0;ndelta_v<FTK_n_delta_v;ndelta_v++){ } */}
  		/* for (nw=0;nw<n_w_max;nw++){ } */}
  	      /* for (nS_sub=0;nS_sub<n_S_sub;nS_sub++){ } */}
  	    /* for (nM_sub=0;nM_sub<n_M_sub;nM_sub++){ } */}
  	  /* if (flag_fast){ } */}
  	/* %%%% */
  	tab = (unsigned long long int)FTK_n_delta_v*(unsigned long long int)n_w_max*(unsigned long long int)n_S_sub*(unsigned long long int)n_M_sub;
  	if (!flag_slow && flag_fast){ memcpy( X_sub_dwSM____ , d_X_sub_dwSM____ , tab*sizeof(double));}
  	/* %%%% */
  	/* %%%% */
  	/*
  	   I_value_sub_dwSM____ = repmat(reshape(f2_dSM___,[FTK.n_delta_v,1,n_S_sub,n_M_sub]),[1,n_w_max,1,1]) .* X_sub_dwSM____; %<-- I_value. ;
  	   I_value_use_dwSM____ = max(0,real(I_value_sub_dwSM____));
  	*/
  	/* %%%% */
  	if (flag_compute_I_value){
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
  	  /* %%%% */
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
  	  /* if (flag_optimize_over_gamma_z==1){ } */}
  	/* if (n_S_sub>0){ } */}
      /* for (nSbatch=0;nSbatch<n_Sbatch;nSbatch++){ } */}
    /* if (n_M_sub>0){ } */}
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  if (flag_compute_I_value){ free(I_value_use_dwSM____); I_value_use_dwSM____=NULL;}
  free(X_sub_dwSM____); X_sub_dwSM____=NULL;
  free(d_X_sub_dwSM____); d_X_sub_dwSM____=NULL;
  free(l2_dSM___); l2_dSM___=NULL;
  free(n2_dSM___); n2_dSM___=NULL;
  free(f2_dSM___); f2_dSM___=NULL;
  if (flag_slow){
    free(svd_USESVUXM_dwSM____); svd_USESVUXM_dwSM____=NULL;
    free(svd_SVUXM_lwsM____); svd_SVUXM_lwsM____=NULL;
    free(svd_SVUXM_lwSM____); svd_SVUXM_lwSM____=NULL;
    _mm_free(svd_SVUXM_0inout_wlSM____); svd_SVUXM_0inout_wlSM____=NULL;
    /* if (flag_slow){ } */}
  if (flag_fast){
    _mm_free(f_svd_SVUXM_0inout_wSMl_realimag____); f_svd_SVUXM_0inout_wSMl_realimag____=NULL;
    /* if (flag_fast){ } */}
  if (flag_slow){
    free(svd_SVUXM_SMwl____); svd_SVUXM_SMwl____=NULL;
    free(svd_SVUXM_SMlw____); svd_SVUXM_SMlw____=NULL;
    free(svd_VUXM_nMwl____); svd_VUXM_nMwl____=NULL;
    free(svd_VUXM_nMlw____); svd_VUXM_nMlw____=NULL;
    /* if (flag_slow){ } */}
  free(index_S_in_Sbatch_); index_S_in_Sbatch_=NULL;
  free(index_M_in_Mbatch_); index_M_in_Mbatch_=NULL;
  if (flag_fast){
    _mm_free(f_svd_VUXM_nMlw_real____); f_svd_VUXM_nMlw_real____=NULL;
    _mm_free(f_svd_VUXM_nMlw_imag____); f_svd_VUXM_nMlw_imag____=NULL;  
    _mm_free(f_svd_SVUXM_SMlw_real____); f_svd_SVUXM_SMlw_real____=NULL;
    _mm_free(f_svd_SVUXM_SMlw_imag____); f_svd_SVUXM_SMlw_imag____=NULL;
    _mm_free(f_svd_SVUXM_lwsM_real____); f_svd_SVUXM_lwsM_real____=NULL;
    _mm_free(f_svd_SVUXM_lwsM_imag____); f_svd_SVUXM_lwsM_imag____=NULL;
    _mm_free(f_svd_SVUXM_lwSM_real____); f_svd_SVUXM_lwSM_real____=NULL;
    _mm_free(f_svd_SVUXM_lwSM_imag____); f_svd_SVUXM_lwSM_imag____=NULL;
    _mm_free(f_svd_USESVUXM_dwSM_real____); f_svd_USESVUXM_dwSM_real____=NULL;
    _mm_free(f_svd_USESVUXM_dwSM_imag____); f_svd_USESVUXM_dwSM_imag____=NULL;
    /* if (flag_fast){ } */}
  /* %%%%%%%%%%%%%%%% */
  if (verbose){ printf(" %% [finished mex_ampmh_X_wSM___15_omp_helper]\n");}
}  

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* The computational routine (slow vs fast) */
/* attempting omp */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
void mex_ampmh_X_wSM___15_omp
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
  int flag_omp=1;
  int verbose=2; // %<-- can use printf within parallel block more generally, since printf is thread safe. ;
#ifndef WITHOUT_MEX
  verbose = verbose && !flag_omp; // %<-- cannot use printf within parallel block in mex, since printf is sent to mexPrintf, which is not thread safe. ;
#endif
  int flag_slow=0;
  int flag_fast=1;
  /* %%%%; */
  int n_M_per_Mbatch_rup=0,n_M_per_Mbatch_256=0;
  int n_S_per_Sbatch_rup=0,n_S_per_Sbatch_256=0;
  int FTK_n_svd_l_rup=0,FTK_n_svd_l_256=0;
  int FTK_n_delta_v_rup=0,FTK_n_delta_v_256=0;
  int n_w_max_rup=0,n_w_max_256=0;
  int pm_n_UX_rank_rup=0,pm_n_UX_rank_256=0;
  int n_S_rup=0,n_S_256=0;
  int n_M_rup=0,n_M_256=0;
  /* %%%%; */
  double complex *FTK_svd_U_d_expiw_s__=NULL;
  float *f_FTK_svd_U_d_expiw_s_dl_real__=NULL;
  float *f_FTK_svd_U_d_expiw_s_dl_imag__=NULL;
  float *f_FTK_svd_U_d_expiw_s_ld_real__=NULL;
  float *f_FTK_svd_U_d_expiw_s_ld_imag__=NULL;
  double complex *CTF_UX_S_k_q_wnS__=NULL;
  double complex *CTF_UX_S_k_q_nSw___=NULL;
  float *f_CTF_UX_S_k_q_wnS_real___=NULL;
  float *f_CTF_UX_S_k_q_wnS_imag___=NULL;
  float *f_CTF_UX_S_k_q_nSw_real___=NULL;
  float *f_CTF_UX_S_k_q_nSw_imag___=NULL;
  double complex *svd_VUXM_lwnM____=NULL;
  float *f_svd_VUXM_lwnM_real____=NULL;
  float *f_svd_VUXM_lwnM_imag____=NULL;
  unsigned long long int na=0,n_a=0;
  /* %%%% */
  int nMbatch=0,n_Mbatch=0;
  int nw=0,nS=0;
  int pm_nUX_rank=0;
  unsigned long long int tab=0,tabA=0,tabB=0,tabC=0;
  double *gamma_z_=NULL;
  /* %%%% */
  fftw_plan fftw_plan_guru_split_dft_plan;
  fftw_iodim fftw_iodim_use;
  fftw_complex * svd_SVUXM_0inout_wlSM____=NULL;
  fftw_complex * svd_SVUXM_0in_wlSM____=NULL;
  fftw_complex * svd_SVUXM_out_wlSM____=NULL;
  /* %%%% */
  fftwf_plan fftwf_plan_guru_split_dft_plan;
  fftwf_iodim fftwf_iodim_use;
  float * f_svd_SVUXM_0inout_wSMl_realimag____=NULL;
  float * f_svd_SVUXM_0in_wSMl_real____=NULL;
  float * f_svd_SVUXM_0in_wSMl_imag____=NULL;
  float * f_svd_SVUXM_out_wSMl_real____=NULL;
  float * f_svd_SVUXM_out_wSMl_imag____=NULL;
  /* %%%% */
  int ntick=0,n_tick = 8;
  clock_t t_start_[n_tick];
  clock_t t_final_[n_tick];
  struct timeval d_start_[n_tick];
  struct timeval d_final_[n_tick];
  long l_msec_[n_tick],l_ssec_[n_tick],l_usec_[n_tick];
  double elct_[n_tick],elrt_[n_tick];
  /* %%%% */
  if (verbose){ printf(" %% [entering mex_ampmh_X_wSM___15_omp]\n");}
  /* %%%%%%%%%%%%%%%% */
  tab = (unsigned long long int)n_w_max;
  gamma_z_ = (double *) malloc(tab*sizeof(double));
  for (nw=0;nw<n_w_max;nw++){ gamma_z_[nw] = 2*PI*(double)nw/(double)n_w_max;}
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
  if (flag_fast){
    f_FTK_svd_U_d_expiw_s_dl_real__ = float__m256_malloc_from_double__(FTK_n_delta_v,FTK_n_svd_l,FTK_svd_U_d_expiw_s_real__);
    f_FTK_svd_U_d_expiw_s_dl_imag__ = float__m256_malloc_from_double__(FTK_n_delta_v,FTK_n_svd_l,FTK_svd_U_d_expiw_s_imag__);
    f_FTK_svd_U_d_expiw_s_ld_real__ = float__m256_malloc_from_double__(FTK_n_svd_l,FTK_n_delta_v,NULL);
    f_FTK_svd_U_d_expiw_s_ld_imag__ = float__m256_malloc_from_double__(FTK_n_svd_l,FTK_n_delta_v,NULL);
    transpose_ps_block1(
			f_FTK_svd_U_d_expiw_s_dl_real__
			,f_FTK_svd_U_d_expiw_s_ld_real__
			,FTK_n_delta_v_rup
			,FTK_n_svd_l_rup
			,FTK_n_delta_v_rup
			,FTK_n_svd_l_rup
			,16
			);
    transpose_ps_block1(
			f_FTK_svd_U_d_expiw_s_dl_imag__
			,f_FTK_svd_U_d_expiw_s_ld_imag__
			,FTK_n_delta_v_rup
			,FTK_n_svd_l_rup
			,FTK_n_delta_v_rup
			,FTK_n_svd_l_rup
			,16
			);
    /* if (flag_fast){ } */}
  if (flag_fast){
    f_CTF_UX_S_k_q_wnS_real___ = float___m256_malloc_from_double___(n_w_max,pm_n_UX_rank,n_S,CTF_UX_S_k_q_wnS_real__);
    f_CTF_UX_S_k_q_wnS_imag___ = float___m256_malloc_from_double___(n_w_max,pm_n_UX_rank,n_S,CTF_UX_S_k_q_wnS_imag__);
    f_CTF_UX_S_k_q_nSw_real___ = float___m256_malloc_from_double___(pm_n_UX_rank,n_S,n_w_max,NULL);
    f_CTF_UX_S_k_q_nSw_imag___ = float___m256_malloc_from_double___(pm_n_UX_rank,n_S,n_w_max,NULL);
    f_svd_VUXM_lwnM_real____ = float____m256_malloc_from_double____(FTK_n_svd_l,n_w_max,pm_n_UX_rank,n_M,svd_VUXM_lwnM_real____);
    f_svd_VUXM_lwnM_imag____ = float____m256_malloc_from_double____(FTK_n_svd_l,n_w_max,pm_n_UX_rank,n_M,svd_VUXM_lwnM_imag____);
    /* if (flag_fast){ } */}
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
	  tabA = 
	    (unsigned long long int)pm_nUX_rank + 
	    ((unsigned long long int)nS + 
	     (unsigned long long int)nw
	     *(unsigned long long int)n_S)
	    *(unsigned long long int)pm_n_UX_rank;
	  CTF_UX_S_k_q_nSw___[tabA] = CTF_UX_S_k_q_wnS__[na]; na++;
	  /* for (nw=0;nw<n_w_max;nw++){ } */}
	/* for (pm_nUX_rank=0;pm_nUX_rank<pm_n_UX_rank;pm_nUX_rank++){ } */}
      /* for (nS=0;nS<n_S;nS++){ } */}
    local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"CTF_UX_S_k_q_nSw___: ");
    /* if (flag_slow){ } */}
  /* %%%%%%%%%%%%%%%% */
  if (flag_fast){
    local_tic(0,t_start_,d_start_);
    tab = 
      (unsigned long long int)pm_n_UX_rank
      *(unsigned long long int)n_S
      *(unsigned long long int)n_w_max;
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
    local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"f_CTF_UX_S_k_q_nSw___: ");
    /* if (flag_fast){ } */}
  /* %%%%%%%%%%%%%%%% */
  if (flag_slow){
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
    _mm_free(svd_SVUXM_0inout_wlSM____); svd_SVUXM_0inout_wlSM____=NULL;
    /* if (flag_slow){ } */}
  /* %%%% */
  if (flag_fast){
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
    _mm_free(f_svd_SVUXM_0inout_wSMl_realimag____); f_svd_SVUXM_0inout_wSMl_realimag____=NULL;
    /* if (flag_fast){ } */}
  /* %%%%%%%%%%%%%%%% */
  local_tic(0,t_start_,d_start_);
  n_Mbatch = ceil((double)n_M/(double)n_M_per_Mbatch);
  if (verbose>1){ printf(" %% n_Mbatch %d\n",n_Mbatch);}
  if (flag_omp){
#pragma omp parallel private(nMbatch)
  { /* begin omp parallel */
    nMbatch = 0;
#pragma omp for schedule(dynamic)
    for (nMbatch=0;nMbatch<n_Mbatch;nMbatch++){
      mex_ampmh_X_wSM___15_omp_helper
	(
	 n_M_per_Mbatch
	 ,n_S_per_Sbatch
	 ,flag_optimize_over_gamma_z
	 ,flag_compute_I_value
	 ,tolerance_master
	 ,FTK_n_svd_l
	 ,FTK_n_delta_v
	 ,FTK_svd_U_d_expiw_s_real__
	 ,FTK_svd_U_d_expiw_s_imag__
	 ,FTK_delta_x_
	 ,FTK_delta_y_
	 ,n_w_max
	 ,pm_n_UX_rank
	 ,n_S
	 ,CTF_UX_S_k_q_wnS_real__
	 ,CTF_UX_S_k_q_wnS_imag__
	 ,CTF_UX_S_l2_
	 ,n_M
	 ,svd_VUXM_lwnM_real____
	 ,svd_VUXM_lwnM_imag____
	 ,UX_M_l2_dM__
	 ,X_wSM___
	 ,delta_x_wSM___
	 ,delta_y_wSM___
	 ,gamma_z_wSM___
	 ,I_value_wSM___
	 /* %%%%; */
	 ,verbose
	 ,flag_slow
	 ,flag_fast
	 /* %%%%; */
	 ,FTK_svd_U_d_expiw_s__
	 ,f_FTK_svd_U_d_expiw_s_dl_real__
	 ,f_FTK_svd_U_d_expiw_s_dl_imag__
	 ,f_FTK_svd_U_d_expiw_s_ld_real__
	 ,f_FTK_svd_U_d_expiw_s_ld_imag__
	 ,CTF_UX_S_k_q_wnS__
	 ,CTF_UX_S_k_q_nSw___
	 ,f_CTF_UX_S_k_q_wnS_real___
	 ,f_CTF_UX_S_k_q_wnS_imag___
	 ,f_CTF_UX_S_k_q_nSw_real___
	 ,f_CTF_UX_S_k_q_nSw_imag___
	 ,svd_VUXM_lwnM____
	 ,f_svd_VUXM_lwnM_real____
	 ,f_svd_VUXM_lwnM_imag____
	 ,gamma_z_
	 /* %%%%; */
	 ,nMbatch
	 ,n_Mbatch
	 ,fftw_plan_guru_split_dft_plan
	 ,fftwf_plan_guru_split_dft_plan
	 );
      /* for (nMbatch=0;nMbatch<n_Mbatch;nMbatch++){ } */}
    /* end omp parallel */}
  /* if (flag_omp){ } */}
  if (!flag_omp){
    for (nMbatch=0;nMbatch<n_Mbatch;nMbatch++){
      mex_ampmh_X_wSM___15_omp_helper
	(
	 n_M_per_Mbatch
	 ,n_S_per_Sbatch
	 ,flag_optimize_over_gamma_z
	 ,flag_compute_I_value
	 ,tolerance_master
	 ,FTK_n_svd_l
	 ,FTK_n_delta_v
	 ,FTK_svd_U_d_expiw_s_real__
	 ,FTK_svd_U_d_expiw_s_imag__
	 ,FTK_delta_x_
	 ,FTK_delta_y_
	 ,n_w_max
	 ,pm_n_UX_rank
	 ,n_S
	 ,CTF_UX_S_k_q_wnS_real__
	 ,CTF_UX_S_k_q_wnS_imag__
	 ,CTF_UX_S_l2_
	 ,n_M
	 ,svd_VUXM_lwnM_real____
	 ,svd_VUXM_lwnM_imag____
	 ,UX_M_l2_dM__
	 ,X_wSM___
	 ,delta_x_wSM___
	 ,delta_y_wSM___
	 ,gamma_z_wSM___
	 ,I_value_wSM___
	 /* %%%%; */
	 ,verbose
	 ,flag_slow
	 ,flag_fast
	 /* %%%%; */
	 ,FTK_svd_U_d_expiw_s__
	 ,f_FTK_svd_U_d_expiw_s_dl_real__
	 ,f_FTK_svd_U_d_expiw_s_dl_imag__
	 ,f_FTK_svd_U_d_expiw_s_ld_real__
	 ,f_FTK_svd_U_d_expiw_s_ld_imag__
	 ,CTF_UX_S_k_q_wnS__
	 ,CTF_UX_S_k_q_nSw___
	 ,f_CTF_UX_S_k_q_wnS_real___
	 ,f_CTF_UX_S_k_q_wnS_imag___
	 ,f_CTF_UX_S_k_q_nSw_real___
	 ,f_CTF_UX_S_k_q_nSw_imag___
	 ,svd_VUXM_lwnM____
	 ,f_svd_VUXM_lwnM_real____
	 ,f_svd_VUXM_lwnM_imag____
	 ,gamma_z_
	 /* %%%%; */
	 ,nMbatch
	 ,n_Mbatch
	 ,fftw_plan_guru_split_dft_plan
	 ,fftwf_plan_guru_split_dft_plan
	 );    
      /* for (nMbatch=0;nMbatch<n_Mbatch;nMbatch++){ } */}
    /* if (!flag_omp){ } */}
  local_toc(0,t_start_,t_final_,d_start_,d_final_,l_msec_,l_ssec_,l_usec_,elct_,elrt_,tab,verbose-1,"f_CTF_UX_S_k_q_nSw___: ");
  /* %%%%%%%%%%%%%%%% */
  free(gamma_z_); gamma_z_=NULL;
  if (flag_slow){
    fftw_destroy_plan(fftw_plan_guru_split_dft_plan);
    /* if (flag_slow){ } */}
  if (flag_fast){
    fftwf_destroy_plan(fftwf_plan_guru_split_dft_plan);
    /* if (flag_fast){ } */}
  if (flag_slow){
    free(CTF_UX_S_k_q_nSw___); CTF_UX_S_k_q_nSw___=NULL;
    free(FTK_svd_U_d_expiw_s__); FTK_svd_U_d_expiw_s__=NULL;
    free(CTF_UX_S_k_q_wnS__); CTF_UX_S_k_q_wnS__=NULL;
    free(svd_VUXM_lwnM____); svd_VUXM_lwnM____=NULL;
    /* if (flag_slow){ } */}
  if (flag_fast){
    _mm_free(f_FTK_svd_U_d_expiw_s_dl_real__); f_FTK_svd_U_d_expiw_s_dl_real__=NULL;
    _mm_free(f_FTK_svd_U_d_expiw_s_dl_imag__); f_FTK_svd_U_d_expiw_s_dl_imag__=NULL;
    _mm_free(f_FTK_svd_U_d_expiw_s_ld_real__); f_FTK_svd_U_d_expiw_s_ld_real__=NULL;
    _mm_free(f_FTK_svd_U_d_expiw_s_ld_imag__); f_FTK_svd_U_d_expiw_s_ld_imag__=NULL;
    _mm_free(f_CTF_UX_S_k_q_nSw_real___); f_CTF_UX_S_k_q_nSw_real___=NULL;
    _mm_free(f_CTF_UX_S_k_q_nSw_imag___); f_CTF_UX_S_k_q_nSw_imag___=NULL;
    _mm_free(f_CTF_UX_S_k_q_wnS_real___); f_CTF_UX_S_k_q_wnS_real___=NULL;
    _mm_free(f_CTF_UX_S_k_q_wnS_imag___); f_CTF_UX_S_k_q_wnS_imag___=NULL;
    _mm_free(f_svd_VUXM_lwnM_real____); f_svd_VUXM_lwnM_real____=NULL;
    _mm_free(f_svd_VUXM_lwnM_imag____); f_svd_VUXM_lwnM_imag____=NULL;  
    /* if (flag_fast){ } */}
  /* %%%%%%%%%%%%%%%% */
  if (verbose){ printf(" %% [finished mex_ampmh_X_wSM___15_omp]\n");}
}
