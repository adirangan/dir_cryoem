void dexcluster_nonbinary_rdrop_d_bkp(int n_r,int n_c,double *A1_rc__,double *A1_cr__,double gamma,int **out_xdrop_p_,int *n_iteration_p,double **trace_p_)
{
  int verbose=1;
  unsigned long long int ulli_rc=(unsigned long long int)n_r*(unsigned long long int)n_c;
  unsigned long long int ulli_r1=0;
  unsigned long long int ulli_1c=0;
  int n_iteration=0;
  int *rdrop_=NULL,*cdrop_=NULL;
  int *rkeep_=NULL,*ckeep_=NULL;
  int n_index_r_rem=0,n_index_c_rem=0;
  int *index_r_rem_=NULL,*index_c_rem_=NULL;
  int nr=0,nc=0;
  int *out_xdrop_trn__=NULL;
  double *trace_trn__=NULL;
  double *A1_rc_=NULL,*A1_cr_=NULL;
  double *er_=NULL,*ec_=NULL,*etA1n_=NULL,*etA2n_=NULL;
  double *QR_pre_=NULL,*QC_pre_=NULL,*QR_pos_=NULL,*QC_pos_=NULL;
  double *QR_pre_r_rem_=NULL,*QC_pre_c_rem_=NULL;
  double *QR_workspace_=NULL,*QC_workspace_=NULL;
  double *etx1n_=NULL,*etx2n_=NULL,*etx1n_stretch_=NULL;
  double *etw1n_=NULL,*etw2n_=NULL;
  double *ety1n_=NULL;
  double *etz1n_=NULL;
  double *ynyten_=NULL,*wnxten_=NULL;
  double *ynzten_=NULL,*y2ny2n_=NULL;
  int *tmp_index_r_=NULL,*tmp_index_c_=NULL;
  int n_r_rmv=0,n_r_rtn=0;
  int *r_rmv_=NULL,*r_rtn_=NULL;
  int n_c_rmv=0,n_c_rtn=0;
  int *c_rmv_=NULL,*c_rtn_=NULL;
  int nx=0,niteration=0;
  double tmp_d_0=0,tmp_d_1=0;
  if (verbose){ printf(" %% [entering dexcluster_nonbinary_rdrop_d] n_r %d n_c %d gamma %0.2f\n",n_r,n_c,gamma);}
  if (verbose>1){ darray_printf_margin(A1_rc__,n_r,n_c," % A1_rc__: ");printf("\n");}
  if (verbose>1){ darray_printf_margin(A1_cr__,n_c,n_r," % A1_cr__: ");printf("\n");}
  /* %%%%%%%%%%%%%%%% */  
  get_xdrop_logscale_array(n_r,n_c,gamma,&n_iteration,&rdrop_,&cdrop_,&rkeep_,&ckeep_);
  memset(cdrop_,0,niteration*sizeof(int));
  if (n_iteration_p!=NULL){ *n_iteration_p = n_iteration;}
  if (verbose>1){ printf(" %% n_iteration %d\n",n_iteration);}
  if (verbose>1){ array_printf(rdrop_,"int",1,n_iteration," % rdrop_: ");}
  if (verbose>1){ array_printf(cdrop_,"int",1,n_iteration," % cdrop_: ");}
  ulli_r1 = (unsigned long long int)n_r*(unsigned long long int)(1+cdrop_[0]);
  ulli_1c = (unsigned long long int)n_c*(unsigned long long int)(1+rdrop_[0]);
  n_index_r_rem = n_r; n_index_c_rem = n_c;
  index_r_rem_ = (int *) malloc1(n_index_r_rem*sizeof(int)); for (nr=0;nr<n_r;nr++){ index_r_rem_[nr]=nr;}
  index_c_rem_ = (int *) malloc1(n_index_c_rem*sizeof(int)); for (nc=0;nc<n_c;nc++){ index_c_rem_[nc]=nc;}
  /* %%%%%%%%%%%%%%%% */  
  if (*out_xdrop_p_==NULL){ (*out_xdrop_p_) = (int *) malloc1(2*n_r*sizeof(int));}
  if (*trace_p_==NULL){ (*trace_p_) = (double *) malloc1(6*n_iteration*sizeof(double));}
  out_xdrop_trn__ = *out_xdrop_p_; trace_trn__ = *trace_p_;
  /* %%%%%%%%%%%%%%%% */  
  index_r_rem_ = (int *) malloc1(n_index_r_rem*sizeof(int)); for (nr=0;nr<n_r;nr++){ index_r_rem_[nr]=nr;}
  index_c_rem_ = (int *) malloc1(n_index_c_rem*sizeof(int)); for (nc=0;nc<n_c;nc++){ index_c_rem_[nc]=nc;}
  etx1n_ = (double *) malloc1(n_c*sizeof(double));
  etx2n_ = (double *) malloc1(n_c*sizeof(double));
  etx1n_stretch_ = (double *) malloc1(n_c*sizeof(double));
  etw1n_ = (double *) malloc1(n_c*sizeof(double));
  etw2n_ = (double *) malloc1(n_c*sizeof(double));
  ety1n_ = (double *) malloc1(n_c*sizeof(double));
  etz1n_ = (double *) malloc1(n_c*sizeof(double));
  ynyten_ = (double *) malloc1(n_r*sizeof(double));
  wnxten_ = (double *) malloc1(n_r*sizeof(double));
  ynzten_ = (double *) malloc1(n_r*sizeof(double));
  y2ny2n_ = (double *) malloc1(n_r*sizeof(double));
  er_ = (double *) malloc1(n_r*sizeof(double));
  ec_ = (double *) malloc1(n_c*sizeof(double));
  etA1n_ = (double *) malloc1(n_c*sizeof(double));
  etA2n_ = (double *) malloc1(n_c*sizeof(double));
  QR_pre_ = (double *) malloc1(n_r*sizeof(double));
  QR_pre_r_rem_ = (double *) malloc1(n_r*sizeof(double));
  QC_pre_ = (double *) malloc1(n_c*sizeof(double));
  QC_pre_c_rem_ = (double *) malloc1(n_c*sizeof(double));
  QR_pos_ = (double *) malloc1(n_r*sizeof(double));
  QC_pos_ = (double *) malloc1(n_c*sizeof(double));
  QR_workspace_ = (double *) malloc1(n_r*sizeof(double));
  QC_workspace_ = (double *) malloc1(n_c*sizeof(double));
  tmp_index_r_ = (int *) malloc1(n_r*sizeof(int));
  tmp_index_c_ = (int *) malloc1(n_c*sizeof(int));
  r_rmv_ = (int *) malloc1(n_r*sizeof(int));
  r_rtn_ = (int *) malloc1(n_r*sizeof(int));
  c_rmv_ = (int *) malloc1(n_c*sizeof(int));
  c_rtn_ = (int *) malloc1(n_c*sizeof(int));
  /* %%%%%%%%%%%%%%%% */  
  for (nr=0;nr<n_r;nr++){ er_[nr] = (double)1.0;}
  for (nc=0;nc<n_c;nc++){ ec_[nc] = (double)1.0;}
  for (nc=0;nc<n_c;nc++){
    A1_rc_ = A1_rc__ + (unsigned long long int)nc*(unsigned long long int)n_r;
    dp_pd_immintrin_loadu(n_r,er_,A1_rc_,&(etA1n_[nc]));
    dp_pd_immintrin_loadu(n_r,A1_rc_,A1_rc_,&(etA2n_[nc]));
    QC_pre_[nc] = etA1n_[nc]*etA1n_[nc] - etA2n_[nc];
    /* for (nc=0;nc<n_c;nc++){ } */}
  for (nr=0;nr<n_r;nr++){
    A1_cr_ = A1_cr__ + (unsigned long long int)nr*(unsigned long long int)n_c;
    dp_pd_immintrin_loadu(n_c,etA1n_,A1_cr_,&(tmp_d_0));
    dp_pd_immintrin_loadu(n_c,A1_cr_,A1_cr_,&(tmp_d_1));
    QR_pre_[nr] = tmp_d_0 - tmp_d_1;
    /* for (nr=0;nr<n_r;nr++){ } */}
  if (verbose>2){ array_printf(etA1n_,"double",1,n_c," % etA1n_: ");}
  if (verbose>2){ array_printf(etA2n_,"double",1,n_c," % etA2n_: ");}
  if (verbose>2){ array_printf(QR_pre_,"double",1,n_r," % QR_pre_: ");}
  if (verbose>2){ array_printf(QC_pre_,"double",1,n_c," % QC_pre_: ");}
  if (verbose>2){ printf("\n");}
  /* %%%%%%%%%%%%%%%% */  
  nx=0;
  for (niteration=0;niteration<n_iteration;niteration++){
    /* %%%%%%%% */
    trace_trn__[0+6*niteration] = 1+niteration;
    trace_trn__[1+6*niteration] = n_index_r_rem;
    trace_trn__[2+6*niteration] = n_index_c_rem;
    dp_pd_immintrin_loadu(n_r,er_,QR_pre_,&(tmp_d_0));
    trace_trn__[3+6*niteration] = tmp_d_0/n_r;
    dp_pd_immintrin_loadu(n_c,ec_,QC_pre_,&(tmp_d_1));
    trace_trn__[4+6*niteration] = tmp_d_1/n_c;
    trace_trn__[5+6*niteration] = 1.0;
    if (verbose>2){ array_printf(&(trace_trn__[6*niteration]),"double",1,6," % trace_trn__: ");}
    /* %%%%%%%% */
    for (nr=0;nr<n_index_r_rem;nr++){ QR_pre_r_rem_[nr] = QR_pre_[index_r_rem_[nr]];}
    dquicksort_index_driver(n_index_r_rem,QR_pre_r_rem_,1,QR_workspace_,tmp_index_r_);
    if (verbose>2){ array_printf(tmp_index_r_,"int",1,n_index_r_rem," % tmp_index_r_: ");}
    n_r_rmv = rdrop_[niteration]; n_r_rtn = n_index_r_rem - n_r_rmv;
    for (nr=0;nr<n_r_rmv;nr++){ r_rmv_[nr] = index_r_rem_[tmp_index_r_[nr]];}
    for (nr=0;nr<n_r_rtn;nr++){ r_rtn_[nr] = index_r_rem_[tmp_index_r_[nr+n_r_rmv]];}
    if (verbose>2){ array_printf(r_rmv_,"int",1,n_r_rmv," % r_rmv_: ");}
    if (verbose>2){ array_printf(r_rtn_,"int",1,n_r_rtn," % r_rtn_: ");}
    /* %%%%%%%% */
    memset(etx1n_,0,n_c*sizeof(double));
    memset(etx2n_,0,n_c*sizeof(double));
    memset(etw1n_,0,n_c*sizeof(double));
    memset(etw2n_,0,n_c*sizeof(double));
    memset(QC_pos_,0,n_c*sizeof(double));
    for (nc=0;nc<n_index_c_rem;nc++){
      etx1n_[nc] = 0; etx2n_[nc] = 0; etw1n_[nc] = 0; etw2n_[nc] = 0;
      A1_rc_ = A1_rc__ + (unsigned long long int)(index_c_rem_[nc]) * (unsigned long long int)n_r;
      for (nr=0;nr<n_r_rmv;nr++){
	tmp_d_0 = A1_rc_[r_rmv_[nr]];
	etx1n_[nc] += tmp_d_0;
	etx2n_[nc] += tmp_d_0*tmp_d_0;
	/* for (nr=0;nr<n_r_rmv;nr++){ } */}
      etw1n_[nc] = etA1n_[index_c_rem_[nc]] - etx1n_[nc];
      etw2n_[nc] = etA2n_[index_c_rem_[nc]] - etx2n_[nc];
      etA1n_[index_c_rem_[nc]] -= etx1n_[nc];
      etA2n_[index_c_rem_[nc]] -= etx2n_[nc];
      QC_pos_[index_c_rem_[nc]] = QC_pre_[index_c_rem_[nc]] - (2*etw1n_[nc]*etx1n_[nc] + etx1n_[nc]*etx1n_[nc] - etx2n_[nc]);
      /* for (nc=0;nc<n_index_c_rem;nc++){ } */}
    memcpy(QC_pre_,QC_pos_,n_c*sizeof(double));
    if (verbose>2){ array_printf(etx1n_,"double",1,n_index_c_rem," % etx1n_: ");}
    if (verbose>2){ array_printf(etx2n_,"double",1,n_index_c_rem," % etx2n_: ");}
    if (verbose>2){ array_printf(etw1n_,"double",1,n_index_c_rem," % etw1n_: ");}
    if (verbose>2){ array_printf(etw2n_,"double",1,n_index_c_rem," % etw2n_: ");}
    if (verbose>2){ array_printf(QC_pos_,"double",1,n_c," % QC_pos_: ");}
    if (verbose>2){ array_printf(QC_pre_,"double",1,n_c," % QC_pre_: ");}
    /* %%%%%%%% */
    for (nc=0;nc<n_index_c_rem;nc++){ QC_pre_c_rem_[nc] = QC_pre_[index_c_rem_[nc]];}
    dquicksort_index_driver(n_index_c_rem,QC_pre_c_rem_,1,QC_workspace_,tmp_index_c_);
    if (verbose>2){ array_printf(tmp_index_c_,"int",1,n_index_c_rem," % tmp_index_c_: ");}
    n_c_rmv = 0; n_c_rtn = n_index_c_rem;
    memcpy(c_rtn_,index_c_rem_,n_c*sizeof(int));
    /* %%%%%%%% */
    memset(etx1n_,0,n_c*sizeof(double));
    memset(etx1n_stretch_,0,n_c*sizeof(double));
    for (nc=0;nc<n_c_rtn;nc++){
      etx1n_[nc] = 0;
      A1_rc_ = A1_rc__ + (unsigned long long int)(c_rtn_[nc]) * (unsigned long long int)n_r;
      for (nr=0;nr<n_r_rmv;nr++){
	tmp_d_0 = A1_rc_[r_rmv_[nr]];
	etx1n_[nc] += tmp_d_0;
	/* for (nr=0;nr<n_r_rmv;nr++){ } */}
      etx1n_stretch_[c_rtn_[nc]] = etx1n_[nc];
      /* for (nc=0;nc<n_c_rtn;nc++){ } */}
    memset(etz1n_,0,n_c*sizeof(double));
    for (nc=0;nc<n_c_rmv;nc++){
      etz1n_[nc] = 0;
      A1_rc_ = A1_rc__ + (unsigned long long int)(c_rmv_[nc]) * (unsigned long long int)n_r;
      for (nr=0;nr<n_r_rmv;nr++){
	tmp_d_0 = A1_rc_[r_rmv_[nr]];
	etz1n_[nc] += tmp_d_0;
	/* for (nr=0;nr<n_r_rmv;nr++){ } */}
      /* for (nc=0;nc<n_c_rmv;nc++){ } */}
    memset(ety1n_,0,n_c*sizeof(double));
    for (nc=0;nc<n_c_rmv;nc++){
      ety1n_[nc] = 0;
      A1_rc_ = A1_rc__ + (unsigned long long int)(c_rmv_[nc]) * (unsigned long long int)n_r;
      for (nr=0;nr<n_r_rtn;nr++){
	tmp_d_0 = A1_rc_[r_rtn_[nr]];
	ety1n_[nc] += tmp_d_0;
	/* for (nr=0;nr<n_r_rtn;nr++){ } */}
      /* for (nc=0;nc<n_c_rmv;nc++){ } */}
    memset(ynyten_,0,n_r*sizeof(double));
    memset(wnxten_,0,n_r*sizeof(double));
    memset(ynzten_,0,n_r*sizeof(double));
    memset(y2ny2n_,0,n_r*sizeof(double));
    for (nr=0;nr<n_r_rtn;nr++){
      ynyten_[nr] = 0; ynzten_[nr] = 0; y2ny2n_[nr] = 0;
      A1_cr_ = A1_cr__ + (unsigned long long int)(r_rtn_[nr]) * (unsigned long long int)n_c;
      for (nc=0;nc<n_c_rmv;nc++){
	tmp_d_0 = A1_cr_[c_rmv_[nc]];
	ynyten_[nr] += tmp_d_0*ety1n_[nc];
	ynzten_[nr] += tmp_d_0*etz1n_[nc];
	y2ny2n_[nr] += tmp_d_0*tmp_d_0;
	/* for (nc=0;nc<n_c_rmv;nc++){ } */}
      if (n_c_rtn>n_c/4){ dp_pd_immintrin_loadu(n_c,A1_cr_,etx1n_stretch_,&(wnxten_[nr]));}
      else{
	wnxten_[nr] = 0;
	for (nc=0;nc<n_c_rtn;nc++){
	  tmp_d_0 = A1_cr_[c_rtn_[nc]];
	  wnxten_[nr] += tmp_d_0*etx1n_[nc];
	  /* for (nc=0;nc<n_c_rtn;nc++){ } */}
	/* use dp_pd */}
      /* for (nr=0;nr<n_r_rtn;nr++){ } */}
    memset(QR_pos_,0,n_r*sizeof(double));
    for (nr=0;nr<n_r_rtn;nr++){
      QR_pos_[r_rtn_[nr]] = QR_pre_[r_rtn_[nr]] - (ynyten_[nr] + wnxten_[nr] + ynzten_[nr] - y2ny2n_[nr]);
      /* for (nr=0;nr<n_r_rtn;nr++){ } */}
    memcpy(QR_pre_,QR_pos_,n_r*sizeof(double));
    if (verbose>2){ array_printf(etx1n_,"double",1,n_c_rtn," % etx1n_: ");}
    if (verbose>2){ array_printf(etz1n_,"double",1,n_c_rmv," % etz1n_: ");}
    if (verbose>2){ array_printf(ety1n_,"double",1,n_c_rmv," % ety1n_: ");}
    if (verbose>2){ array_printf(ynyten_,"double",1,n_r_rtn," % ynyten_: ");}
    if (verbose>2){ array_printf(ynzten_,"double",1,n_r_rtn," % ynzten_: ");}
    if (verbose>2){ array_printf(y2ny2n_,"double",1,n_r_rtn," % y2ny2n_: ");}
    if (verbose>2){ array_printf(wnxten_,"double",1,n_r_rtn," % wnxten_: ");}
    if (verbose>2){ array_printf(QR_pos_,"double",1,n_r," % QR_pos_: ");}
    if (verbose>2){ array_printf(QR_pre_,"double",1,n_r," % QR_pre_: ");}
    /* %%%%%%%% */
    for (nc=0;nc<n_c_rmv;nc++){ etA1n_[c_rmv_[nc]] = 0; etA2n_[c_rmv_[nc]] = 0; QC_pre_[c_rmv_[nc]] = 0;}
    if (verbose>2){ array_printf(etA1n_,"double",1,n_c," % etA1n_: ");}
    if (verbose>2){ array_printf(etA2n_,"double",1,n_c," % etA2n_: ");}
    if (verbose>2){ array_printf(QC_pre_,"double",1,n_c," % QC_pre_: ");}
    iquicksort_index_driver(n_r_rtn,r_rtn_,1,index_r_rem_,tmp_index_r_); n_index_r_rem = n_r_rtn;
    if (verbose>2){ array_printf(index_r_rem_,"int",1,n_index_r_rem," % index_r_rem_: ");}
    for (nr=0;nr<n_r_rmv;nr++){
      out_xdrop_trn__[ 0 + nx*2 ] = r_rmv_[nr]; out_xdrop_trn__[ 1 + nx*2 ] = -1;
      if (verbose>2){ array_printf(&(out_xdrop_trn__[2*nx]),"int",1,2," % out_xdrop_trn__: ");}
      nx += 1;
      /* for (nr=0;nr<n_r_rmv;nr++){ } */}
    if (verbose>2){ printf(" %% nx %d\n",nx);}
    if (verbose>2){ printf(" %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%\n");}
    /* %%%%%%%% */
    /* for (niteration=0;niteration<n_iteration;niteration++){ } */}
  /* %%%%%%%%%%%%%%%% */  
  free1(index_r_rem_);index_r_rem_=NULL;
  free1(index_c_rem_);index_c_rem_=NULL;
  free1(etx1n_);etx1n_=NULL;
  free1(etx2n_);etx2n_=NULL;
  free1(etx1n_stretch_);etx1n_stretch_=NULL;
  free1(etw1n_);etw1n_=NULL;
  free1(etw2n_);etw2n_=NULL;
  free1(ety1n_);ety1n_=NULL;
  free1(etz1n_);etz1n_=NULL;
  free1(ynyten_);ynyten_=NULL;
  free1(wnxten_);wnxten_=NULL;
  free1(ynzten_);ynzten_=NULL;
  free1(y2ny2n_);y2ny2n_=NULL;
  free1(er_);er_=NULL;
  free1(ec_);ec_=NULL;
  free1(etA1n_);etA1n_=NULL;
  free1(etA2n_);etA2n_=NULL;
  free1(QR_pre_);QR_pre_=NULL;
  free1(QR_pre_r_rem_);QR_pre_r_rem_=NULL;
  free1(QC_pre_);QC_pre_=NULL;
  free1(QC_pre_c_rem_);QC_pre_c_rem_=NULL;
  free1(QR_pos_);QR_pos_=NULL;
  free1(QC_pos_);QC_pos_=NULL;
  free1(QR_workspace_);QR_workspace_=NULL;
  free1(QC_workspace_);QC_workspace_=NULL;
  free1(tmp_index_r_);tmp_index_r_=NULL;
  free1(tmp_index_c_);tmp_index_c_=NULL;
  free1(r_rmv_);r_rmv_=NULL;
  free1(r_rtn_);r_rtn_=NULL;
  free1(c_rmv_);c_rmv_=NULL;
  free1(c_rtn_);c_rtn_=NULL;
  /* %%%% */
  free1(rdrop_);rdrop_=NULL;
  free1(cdrop_);cdrop_=NULL;
  free1(rkeep_);rkeep_=NULL;
  free1(ckeep_);ckeep_=NULL;
  /* %%%%%%%%%%%%%%%% */  
  if (verbose){ printf(" %% [finished dexcluster_nonbinary_rdrop_d] n_r %d n_c %d gamma %0.2f\n",n_r,n_c,gamma);}
}

void dexcluster_nonbinary_rdrop_d(int n_r,int n_c,double *A1_rc__,double *A1_cr__,double gamma,int **out_xdrop_p_,int *n_iteration_p,double **trace_p_)
{
  int verbose=1;
  unsigned long long int ulli_rc=(unsigned long long int)n_r*(unsigned long long int)n_c;
  unsigned long long int ulli_r1=0;
  unsigned long long int ulli_1c=0;
  int n_iteration=0;
  int *rdrop_=NULL,*cdrop_=NULL;
  int *rkeep_=NULL,*ckeep_=NULL;
  int n_index_r_rem=0;
  int *index_r_rem_=NULL;
  int nr=0,nc=0;
  int *out_xdrop_trn__=NULL;
  double *trace_trn__=NULL;
  double *A1_rc_=NULL,*A1_cr_=NULL;
  double *er_=NULL,*etA1n_=NULL;
  double *QR_pre_=NULL,*QR_pos_=NULL;
  double *QR_pre_r_rem_=NULL;
  double *QR_workspace_=NULL;
  double *etx1n_=NULL;
  double *wnxten_=NULL;
  int *tmp_index_r_=NULL,*tmp_index_c_=NULL;
  int n_r_rmv=0,n_r_rtn=0;
  int *r_rmv_=NULL,*r_rtn_=NULL;
  int nx=0,niteration=0;
  double tmp_d_0=0,tmp_d_1=0;
  if (verbose){ printf(" %% [entering dexcluster_nonbinary_rdrop_d] n_r %d n_c %d gamma %0.2f\n",n_r,n_c,gamma);}
  if (verbose>1){ darray_printf_margin(A1_rc__,n_r,n_c," % A1_rc__: ");printf("\n");}
  if (verbose>1){ darray_printf_margin(A1_cr__,n_c,n_r," % A1_cr__: ");printf("\n");}
  /* %%%%%%%%%%%%%%%% */  
  get_xdrop_logscale_array(n_r,n_c,gamma,&n_iteration,&rdrop_,&cdrop_,&rkeep_,&ckeep_);
  memset(cdrop_,0,niteration*sizeof(int));
  if (n_iteration_p!=NULL){ *n_iteration_p = n_iteration;}
  if (verbose>1){ printf(" %% n_iteration %d\n",n_iteration);}
  if (verbose>1){ array_printf(rdrop_,"int",1,n_iteration," % rdrop_: ");}
  if (verbose>1){ array_printf(cdrop_,"int",1,n_iteration," % cdrop_: ");}
  ulli_r1 = (unsigned long long int)n_r*(unsigned long long int)(1+cdrop_[0]);
  ulli_1c = (unsigned long long int)n_c*(unsigned long long int)(1+rdrop_[0]);
  n_index_r_rem = n_r;
  index_r_rem_ = (int *) malloc1(n_index_r_rem*sizeof(int)); for (nr=0;nr<n_r;nr++){ index_r_rem_[nr]=nr;}
  /* %%%%%%%%%%%%%%%% */  
  if (*out_xdrop_p_==NULL){ (*out_xdrop_p_) = (int *) malloc1(2*n_r*sizeof(int));}
  if (*trace_p_==NULL){ (*trace_p_) = (double *) malloc1(6*n_iteration*sizeof(double));}
  out_xdrop_trn__ = *out_xdrop_p_; trace_trn__ = *trace_p_;
  /* %%%%%%%%%%%%%%%% */  
  index_r_rem_ = (int *) malloc1(n_index_r_rem*sizeof(int)); for (nr=0;nr<n_r;nr++){ index_r_rem_[nr]=nr;}
  etx1n_ = (double *) malloc1(n_c*sizeof(double));
  wnxten_ = (double *) malloc1(n_r*sizeof(double));
  er_ = (double *) malloc1(n_r*sizeof(double));
  etA1n_ = (double *) malloc1(n_c*sizeof(double));
  QR_pre_ = (double *) malloc1(n_r*sizeof(double));
  QR_pre_r_rem_ = (double *) malloc1(n_r*sizeof(double));
  QR_pos_ = (double *) malloc1(n_r*sizeof(double));
  QR_workspace_ = (double *) malloc1(n_r*sizeof(double));
  tmp_index_r_ = (int *) malloc1(n_r*sizeof(int));
  r_rmv_ = (int *) malloc1(n_r*sizeof(int));
  r_rtn_ = (int *) malloc1(n_r*sizeof(int));
  /* %%%%%%%%%%%%%%%% */  
  for (nr=0;nr<n_r;nr++){ er_[nr] = (double)1.0;}
  for (nc=0;nc<n_c;nc++){
    A1_rc_ = A1_rc__ + (unsigned long long int)nc*(unsigned long long int)n_r;
    dp_pd_immintrin_loadu(n_r,er_,A1_rc_,&(etA1n_[nc]));
    /* for (nc=0;nc<n_c;nc++){ } */}
  for (nr=0;nr<n_r;nr++){
    A1_cr_ = A1_cr__ + (unsigned long long int)nr*(unsigned long long int)n_c;
    dp_pd_immintrin_loadu(n_c,etA1n_,A1_cr_,&(tmp_d_0));
    dp_pd_immintrin_loadu(n_c,A1_cr_,A1_cr_,&(tmp_d_1));
    QR_pre_[nr] = tmp_d_0 - tmp_d_1;
    /* for (nr=0;nr<n_r;nr++){ } */}
  if (verbose>2){ array_printf(etA1n_,"double",1,n_c," % etA1n_: ");}
  if (verbose>2){ array_printf(QR_pre_,"double",1,n_r," % QR_pre_: ");}
  if (verbose>2){ printf("\n");}
  /* %%%%%%%%%%%%%%%% */  
  nx=0;
  for (niteration=0;niteration<n_iteration;niteration++){
    /* %%%%%%%% */
    trace_trn__[0+6*niteration] = 1+niteration;
    trace_trn__[1+6*niteration] = n_index_r_rem;
    trace_trn__[2+6*niteration] = n_c;
    dp_pd_immintrin_loadu(n_r,er_,QR_pre_,&(tmp_d_0));
    trace_trn__[3+6*niteration] = tmp_d_0/n_r;
    trace_trn__[4+6*niteration] = 0.0;
    trace_trn__[5+6*niteration] = 1.0;
    if (verbose>2){ array_printf(&(trace_trn__[6*niteration]),"double",1,6," % trace_trn__: ");}
    /* %%%%%%%% */
    for (nr=0;nr<n_index_r_rem;nr++){ QR_pre_r_rem_[nr] = QR_pre_[index_r_rem_[nr]];}
    dquicksort_index_driver(n_index_r_rem,QR_pre_r_rem_,1,QR_workspace_,tmp_index_r_);
    if (verbose>2){ array_printf(tmp_index_r_,"int",1,n_index_r_rem," % tmp_index_r_: ");}
    n_r_rmv = rdrop_[niteration]; n_r_rtn = n_index_r_rem - n_r_rmv;
    for (nr=0;nr<n_r_rmv;nr++){ r_rmv_[nr] = index_r_rem_[tmp_index_r_[nr]];}
    for (nr=0;nr<n_r_rtn;nr++){ r_rtn_[nr] = index_r_rem_[tmp_index_r_[nr+n_r_rmv]];}
    if (verbose>2){ array_printf(r_rmv_,"int",1,n_r_rmv," % r_rmv_: ");}
    if (verbose>2){ array_printf(r_rtn_,"int",1,n_r_rtn," % r_rtn_: ");}
    /* %%%%%%%% */
    if (n_r_rmv==1){
      A1_cr_ = A1_cr__ + (unsigned long long int)(r_rmv_[0]) * (unsigned long long int)n_c;
      memcpy(etx1n_,A1_cr_,n_c*sizeof(double));
      /* if (n_r_rmv==1){ } */}
    else /* if (n_r_rmv!=1) */{
      memset(etx1n_,0,n_c*sizeof(double));
      for (nr=0;nr<n_r_rmv;nr++){
	A1_cr_ = A1_cr__ + (unsigned long long int)(r_rmv_[nr]) * (unsigned long long int)n_c;
	for (nc=0;nc<n_c;nc++){
	  etx1n_[nc] += A1_cr__[nc];
	  /* for (nc=0;nc<n_c;nc++){ } */}
	/* for (nr=0;nr<n_r_rmv;nr++){ } */}
      /* else if (n_r_rmv!=1){ } */}
    if (0){ /* old version */
      memset(etx1n_,0,n_c*sizeof(double));
      for (nc=0;nc<n_c;nc++){
	etx1n_[nc] = 0;
	A1_rc_ = A1_rc__ + (unsigned long long int)(nc) * (unsigned long long int)n_r;
	for (nr=0;nr<n_r_rmv;nr++){
	  tmp_d_0 = A1_rc_[r_rmv_[nr]];
	  etx1n_[nc] += tmp_d_0;
	  /* for (nr=0;nr<n_r_rmv;nr++){ } */}
	/* for (nc=0;nc<n_c;nc++){ } */}
      /* if (0){ } */}
    memset(wnxten_,0,n_r*sizeof(double));
    for (nr=0;nr<n_r_rtn;nr++){
      A1_cr_ = A1_cr__ + (unsigned long long int)(r_rtn_[nr]) * (unsigned long long int)n_c;
      dp_pd_immintrin_loadu(n_c,A1_cr_,etx1n_,&(wnxten_[nr]));
      /* for (nr=0;nr<n_r_rtn;nr++){ } */}
    memset(QR_pos_,0,n_r*sizeof(double));
    for (nr=0;nr<n_r_rtn;nr++){
      QR_pos_[r_rtn_[nr]] = QR_pre_[r_rtn_[nr]] - wnxten_[nr];
      /* for (nr=0;nr<n_r_rtn;nr++){ } */}
    memcpy(QR_pre_,QR_pos_,n_r*sizeof(double));
    if (verbose>2){ array_printf(etx1n_,"double",1,n_c," % etx1n_: ");}
    if (verbose>2){ array_printf(wnxten_,"double",1,n_r_rtn," % wnxten_: ");}
    if (verbose>2){ array_printf(QR_pos_,"double",1,n_r," % QR_pos_: ");}
    if (verbose>2){ array_printf(QR_pre_,"double",1,n_r," % QR_pre_: ");}
    /* %%%%%%%% */
    iquicksort_index_driver(n_r_rtn,r_rtn_,1,index_r_rem_,tmp_index_r_); n_index_r_rem = n_r_rtn;
    if (verbose>2){ array_printf(index_r_rem_,"int",1,n_index_r_rem," % index_r_rem_: ");}
    for (nr=0;nr<n_r_rmv;nr++){
      out_xdrop_trn__[ 0 + nx*2 ] = r_rmv_[nr]; out_xdrop_trn__[ 1 + nx*2 ] = -1;
      if (verbose>2){ array_printf(&(out_xdrop_trn__[2*nx]),"int",1,2," % out_xdrop_trn__: ");}
      nx += 1;
      /* for (nr=0;nr<n_r_rmv;nr++){ } */}
    if (verbose>2){ printf(" %% nx %d\n",nx);}
    if (verbose>2){ printf(" %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%\n");}
    /* %%%%%%%% */
    /* for (niteration=0;niteration<n_iteration;niteration++){ } */}
  /* %%%%%%%%%%%%%%%% */  
  free1(index_r_rem_);index_r_rem_=NULL;
  free1(etx1n_);etx1n_=NULL;
  free1(wnxten_);wnxten_=NULL;
  free1(er_);er_=NULL;
  free1(etA1n_);etA1n_=NULL;
  free1(QR_pre_);QR_pre_=NULL;
  free1(QR_pre_r_rem_);QR_pre_r_rem_=NULL;
  free1(QR_pos_);QR_pos_=NULL;
  free1(QR_workspace_);QR_workspace_=NULL;
  free1(tmp_index_r_);tmp_index_r_=NULL;
  free1(r_rmv_);r_rmv_=NULL;
  free1(r_rtn_);r_rtn_=NULL;
  /* %%%% */
  free1(rdrop_);rdrop_=NULL;
  free1(cdrop_);cdrop_=NULL;
  free1(rkeep_);rkeep_=NULL;
  free1(ckeep_);ckeep_=NULL;
  /* %%%%%%%%%%%%%%%% */  
  if (verbose){ printf(" %% [finished dexcluster_nonbinary_rdrop_d] n_r %d n_c %d gamma %0.2f\n",n_r,n_c,gamma);}
}

