#ifndef _MONOLITH
#include "playpark_header.h"
#endif /* _MONOLITH */

void dexcluster_nonbinary_f(int n_r,int n_c,float *A1_rc__,float *A1_cr__,double gamma,int **out_xdrop_p_,int *n_iteration_p,double **trace_p_)
{
  int verbose=0;
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
  float *A1_rc_=NULL,*A1_cr_=NULL;
  float *er_=NULL,*ec_=NULL,*etA1n_=NULL,*etA2n_=NULL;
  float *QR_pre_=NULL,*QC_pre_=NULL,*QR_pos_=NULL,*QC_pos_=NULL;
  float *QR_pre_r_rem_=NULL,*QC_pre_c_rem_=NULL;
  float *QR_workspace_=NULL,*QC_workspace_=NULL;
  float *etx1n_=NULL,*etx2n_=NULL,*etx1n_stretch_=NULL;
  float *etw1n_=NULL,*etw2n_=NULL;
  float *ety1n_=NULL;
  float *etz1n_=NULL;
  float *ynyten_=NULL,*wnxten_=NULL;
  float *ynzten_=NULL,*y2ny2n_=NULL;
  int *tmp_index_r_=NULL,*tmp_index_c_=NULL;
  int n_r_rmv=0,n_r_rtn=0;
  int *r_rmv_=NULL,*r_rtn_=NULL;
  int n_c_rmv=0,n_c_rtn=0;
  int *c_rmv_=NULL,*c_rtn_=NULL;
  int nx=0,niteration=0;
  float tmp_f_0=0,tmp_f_1=0;
  if (verbose){ printf(" %% [entering dexcluster_nonbinary_f] n_r %d n_c %d gamma %0.2f\n",n_r,n_c,gamma);}
  if (verbose>1){ farray_printf_margin(A1_rc__,n_r,n_c," % A1_rc__: ");printf("\n");}
  if (verbose>1){ farray_printf_margin(A1_cr__,n_c,n_r," % A1_cr__: ");printf("\n");}
  /* %%%%%%%%%%%%%%%% */  
  get_xdrop_logscale_array(n_r,n_c,gamma,&n_iteration,&rdrop_,&cdrop_,&rkeep_,&ckeep_);
  if (n_iteration_p!=NULL){ *n_iteration_p = n_iteration;}
  if (verbose>1){ printf(" %% n_iteration %d\n",n_iteration);}
  if (verbose>1){ array_printf(rdrop_,"int",1,n_iteration," % rdrop_: ");}
  if (verbose>1){ array_printf(cdrop_,"int",1,n_iteration," % cdrop_: ");}
  ulli_r1 = (unsigned long long int)n_r*(unsigned long long int)(1+cdrop_[0]);
  ulli_1c = (unsigned long long int)n_c*(unsigned long long int)(1+rdrop_[0]);
  /* %%%%%%%%%%%%%%%% */  
  if (*out_xdrop_p_==NULL){ (*out_xdrop_p_) = (int *) malloc1(2*(n_r+n_c)*sizeof(int));}
  if (*trace_p_==NULL){ (*trace_p_) = (double *) malloc1(6*n_iteration*sizeof(double));}
  out_xdrop_trn__ = *out_xdrop_p_; trace_trn__ = *trace_p_;
  /* %%%%%%%%%%%%%%%% */  
  n_index_r_rem = n_r; n_index_c_rem = n_c;
  index_r_rem_ = (int *) malloc1(n_index_r_rem*sizeof(int)); for (nr=0;nr<n_r;nr++){ index_r_rem_[nr]=nr;}
  index_c_rem_ = (int *) malloc1(n_index_c_rem*sizeof(int)); for (nc=0;nc<n_c;nc++){ index_c_rem_[nc]=nc;}
  etx1n_ = (float *) malloc1(n_c*sizeof(float));
  etx2n_ = (float *) malloc1(n_c*sizeof(float));
  etx1n_stretch_ = (float *) malloc1(n_c*sizeof(float));
  etw1n_ = (float *) malloc1(n_c*sizeof(float));
  etw2n_ = (float *) malloc1(n_c*sizeof(float));
  ety1n_ = (float *) malloc1(n_c*sizeof(float));
  etz1n_ = (float *) malloc1(n_c*sizeof(float));
  ynyten_ = (float *) malloc1(n_r*sizeof(float));
  wnxten_ = (float *) malloc1(n_r*sizeof(float));
  ynzten_ = (float *) malloc1(n_r*sizeof(float));
  y2ny2n_ = (float *) malloc1(n_r*sizeof(float));
  er_ = (float *) malloc1(n_r*sizeof(float));
  ec_ = (float *) malloc1(n_c*sizeof(float));
  etA1n_ = (float *) malloc1(n_c*sizeof(float));
  etA2n_ = (float *) malloc1(n_c*sizeof(float));
  QR_pre_ = (float *) malloc1(n_r*sizeof(float));
  QR_pre_r_rem_ = (float *) malloc1(n_r*sizeof(float));
  QC_pre_ = (float *) malloc1(n_c*sizeof(float));
  QC_pre_c_rem_ = (float *) malloc1(n_c*sizeof(float));
  QR_pos_ = (float *) malloc1(n_r*sizeof(float));
  QC_pos_ = (float *) malloc1(n_c*sizeof(float));
  QR_workspace_ = (float *) malloc1(n_r*sizeof(float));
  QC_workspace_ = (float *) malloc1(n_c*sizeof(float));
  tmp_index_r_ = (int *) malloc1(n_r*sizeof(int));
  tmp_index_c_ = (int *) malloc1(n_c*sizeof(int));
  r_rmv_ = (int *) malloc1(n_r*sizeof(int));
  r_rtn_ = (int *) malloc1(n_r*sizeof(int));
  c_rmv_ = (int *) malloc1(n_c*sizeof(int));
  c_rtn_ = (int *) malloc1(n_c*sizeof(int));
  /* %%%%%%%%%%%%%%%% */  
  for (nr=0;nr<n_r;nr++){ er_[nr] = (float)1.0;}
  for (nc=0;nc<n_c;nc++){ ec_[nc] = (float)1.0;}
  for (nc=0;nc<n_c;nc++){
    A1_rc_ = A1_rc__ + (unsigned long long int)nc*(unsigned long long int)n_r;
    dp_ps_immintrin_loadu(n_r,er_,A1_rc_,&(etA1n_[nc]));
    dp_ps_immintrin_loadu(n_r,A1_rc_,A1_rc_,&(etA2n_[nc]));
    QC_pre_[nc] = etA1n_[nc]*etA1n_[nc] - etA2n_[nc];
    /* for (nc=0;nc<n_c;nc++){ } */}
  for (nr=0;nr<n_r;nr++){
    A1_cr_ = A1_cr__ + (unsigned long long int)nr*(unsigned long long int)n_c;
    dp_ps_immintrin_loadu(n_c,etA1n_,A1_cr_,&(tmp_f_0));
    dp_ps_immintrin_loadu(n_c,A1_cr_,A1_cr_,&(tmp_f_1));
    QR_pre_[nr] = tmp_f_0 - tmp_f_1;
    /* for (nr=0;nr<n_r;nr++){ } */}
  if (verbose>2){ array_printf(etA1n_,"float",1,n_c," % etA1n_: ");}
  if (verbose>2){ array_printf(etA2n_,"float",1,n_c," % etA2n_: ");}
  if (verbose>2){ array_printf(QR_pre_,"float",1,n_r," % QR_pre_: ");}
  if (verbose>2){ array_printf(QC_pre_,"float",1,n_c," % QC_pre_: ");}
  if (verbose>2){ printf("\n");}
  /* %%%%%%%%%%%%%%%% */  
  nx=0;
  for (niteration=0;niteration<n_iteration;niteration++){
    /* %%%%%%%% */
    trace_trn__[0+6*niteration] = 1+niteration;
    trace_trn__[1+6*niteration] = n_index_r_rem;
    trace_trn__[2+6*niteration] = n_index_c_rem;
    dp_ps_immintrin_loadu(n_r,er_,QR_pre_,&(tmp_f_0));
    trace_trn__[3+6*niteration] = tmp_f_0/n_r;
    dp_ps_immintrin_loadu(n_c,ec_,QC_pre_,&(tmp_f_1));
    trace_trn__[4+6*niteration] = tmp_f_1/n_c;
    trace_trn__[5+6*niteration] = 1.0;
    if (verbose>2){ array_printf(&(trace_trn__[6*niteration]),"double",1,6," % trace_trn__: ");}
    /* %%%%%%%% */
    for (nr=0;nr<n_index_r_rem;nr++){ QR_pre_r_rem_[nr] = QR_pre_[index_r_rem_[nr]];}
    fquicksort_index_driver(n_index_r_rem,QR_pre_r_rem_,1,QR_workspace_,tmp_index_r_);
    if (verbose>2){ array_printf(tmp_index_r_,"int",1,n_index_r_rem," % tmp_index_r_: ");}
    n_r_rmv = rdrop_[niteration]; n_r_rtn = n_index_r_rem - n_r_rmv;
    for (nr=0;nr<n_r_rmv;nr++){ r_rmv_[nr] = index_r_rem_[tmp_index_r_[nr]];}
    for (nr=0;nr<n_r_rtn;nr++){ r_rtn_[nr] = index_r_rem_[tmp_index_r_[nr+n_r_rmv]];}
    if (verbose>2){ array_printf(r_rmv_,"int",1,n_r_rmv," % r_rmv_: ");}
    if (verbose>2){ array_printf(r_rtn_,"int",1,n_r_rtn," % r_rtn_: ");}
    /* %%%%%%%% */
    memset(etx1n_,0,n_c*sizeof(float));
    memset(etx2n_,0,n_c*sizeof(float));
    memset(etw1n_,0,n_c*sizeof(float));
    memset(etw2n_,0,n_c*sizeof(float));
    memset(QC_pos_,0,n_c*sizeof(float));
    for (nc=0;nc<n_index_c_rem;nc++){
      etx1n_[nc] = 0; etx2n_[nc] = 0; etw1n_[nc] = 0; etw2n_[nc] = 0;
      A1_rc_ = A1_rc__ + (unsigned long long int)(index_c_rem_[nc]) * (unsigned long long int)n_r;
      for (nr=0;nr<n_r_rmv;nr++){
	tmp_f_0 = A1_rc_[r_rmv_[nr]];
	etx1n_[nc] += tmp_f_0;
	etx2n_[nc] += tmp_f_0*tmp_f_0;
	/* for (nr=0;nr<n_r_rmv;nr++){ } */}
      etw1n_[nc] = etA1n_[index_c_rem_[nc]] - etx1n_[nc];
      etw2n_[nc] = etA2n_[index_c_rem_[nc]] - etx2n_[nc];
      etA1n_[index_c_rem_[nc]] -= etx1n_[nc];
      etA2n_[index_c_rem_[nc]] -= etx2n_[nc];
      QC_pos_[index_c_rem_[nc]] = QC_pre_[index_c_rem_[nc]] - (2*etw1n_[nc]*etx1n_[nc] + etx1n_[nc]*etx1n_[nc] - etx2n_[nc]);
      /* for (nc=0;nc<n_index_c_rem;nc++){ } */}
    memcpy(QC_pre_,QC_pos_,n_c*sizeof(float));
    if (verbose>2){ array_printf(etx1n_,"float",1,n_index_c_rem," % etx1n_: ");}
    if (verbose>2){ array_printf(etx2n_,"float",1,n_index_c_rem," % etx2n_: ");}
    if (verbose>2){ array_printf(etw1n_,"float",1,n_index_c_rem," % etw1n_: ");}
    if (verbose>2){ array_printf(etw2n_,"float",1,n_index_c_rem," % etw2n_: ");}
    if (verbose>2){ array_printf(QC_pos_,"float",1,n_c," % QC_pos_: ");}
    if (verbose>2){ array_printf(QC_pre_,"float",1,n_c," % QC_pre_: ");}
    /* %%%%%%%% */
    for (nc=0;nc<n_index_c_rem;nc++){ QC_pre_c_rem_[nc] = QC_pre_[index_c_rem_[nc]];}
    fquicksort_index_driver(n_index_c_rem,QC_pre_c_rem_,1,QC_workspace_,tmp_index_c_);
    if (verbose>2){ array_printf(tmp_index_c_,"int",1,n_index_c_rem," % tmp_index_c_: ");}
    n_c_rmv = cdrop_[niteration]; n_c_rtn = n_index_c_rem - n_c_rmv;
    for (nc=0;nc<n_c_rmv;nc++){ c_rmv_[nc] = index_c_rem_[tmp_index_c_[nc]];}
    for (nc=0;nc<n_c_rtn;nc++){ c_rtn_[nc] = index_c_rem_[tmp_index_c_[nc+n_c_rmv]];}
    if (verbose>2){ array_printf(c_rmv_,"int",1,n_c_rmv," % c_rmv_: ");}
    if (verbose>2){ array_printf(c_rtn_,"int",1,n_c_rtn," % c_rtn_: ");}
    /* %%%%%%%% */
    memset(etx1n_,0,n_c*sizeof(float));
    memset(etx1n_stretch_,0,n_c*sizeof(float));
    for (nc=0;nc<n_c_rtn;nc++){
      etx1n_[nc] = 0;
      A1_rc_ = A1_rc__ + (unsigned long long int)(c_rtn_[nc]) * (unsigned long long int)n_r;
      for (nr=0;nr<n_r_rmv;nr++){
	tmp_f_0 = A1_rc_[r_rmv_[nr]];
	etx1n_[nc] += tmp_f_0;
	/* for (nr=0;nr<n_r_rmv;nr++){ } */}
      etx1n_stretch_[c_rtn_[nc]] = etx1n_[nc];
      /* for (nc=0;nc<n_c_rtn;nc++){ } */}
    memset(etz1n_,0,n_c*sizeof(float));
    for (nc=0;nc<n_c_rmv;nc++){
      etz1n_[nc] = 0;
      A1_rc_ = A1_rc__ + (unsigned long long int)(c_rmv_[nc]) * (unsigned long long int)n_r;
      for (nr=0;nr<n_r_rmv;nr++){
	tmp_f_0 = A1_rc_[r_rmv_[nr]];
	etz1n_[nc] += tmp_f_0;
	/* for (nr=0;nr<n_r_rmv;nr++){ } */}
      /* for (nc=0;nc<n_c_rmv;nc++){ } */}
    memset(ety1n_,0,n_c*sizeof(float));
    for (nc=0;nc<n_c_rmv;nc++){
      ety1n_[nc] = 0;
      A1_rc_ = A1_rc__ + (unsigned long long int)(c_rmv_[nc]) * (unsigned long long int)n_r;
      for (nr=0;nr<n_r_rtn;nr++){
	tmp_f_0 = A1_rc_[r_rtn_[nr]];
	ety1n_[nc] += tmp_f_0;
	/* for (nr=0;nr<n_r_rtn;nr++){ } */}
      /* for (nc=0;nc<n_c_rmv;nc++){ } */}
    memset(ynyten_,0,n_r*sizeof(float));
    memset(wnxten_,0,n_r*sizeof(float));
    memset(ynzten_,0,n_r*sizeof(float));
    memset(y2ny2n_,0,n_r*sizeof(float));
    for (nr=0;nr<n_r_rtn;nr++){
      ynyten_[nr] = 0; ynzten_[nr] = 0; y2ny2n_[nr] = 0;
      A1_cr_ = A1_cr__ + (unsigned long long int)(r_rtn_[nr]) * (unsigned long long int)n_c;
      for (nc=0;nc<n_c_rmv;nc++){
	tmp_f_0 = A1_cr_[c_rmv_[nc]];
	ynyten_[nr] += tmp_f_0*ety1n_[nc];
	ynzten_[nr] += tmp_f_0*etz1n_[nc];
	y2ny2n_[nr] += tmp_f_0*tmp_f_0;
	/* for (nc=0;nc<n_c_rmv;nc++){ } */}
      if (n_c_rtn>n_c/4){ dp_ps_immintrin_loadu(n_c,A1_cr_,etx1n_stretch_,&(wnxten_[nr]));}
      else{
	wnxten_[nr] = 0;
	for (nc=0;nc<n_c_rtn;nc++){
	  tmp_f_0 = A1_cr_[c_rtn_[nc]];
	  wnxten_[nr] += tmp_f_0*etx1n_[nc];
	  /* for (nc=0;nc<n_c_rtn;nc++){ } */}
	/* use dp_ps */}
      /* for (nr=0;nr<n_r_rtn;nr++){ } */}
    memset(QR_pos_,0,n_r*sizeof(float));
    for (nr=0;nr<n_r_rtn;nr++){
      QR_pos_[r_rtn_[nr]] = QR_pre_[r_rtn_[nr]] - (ynyten_[nr] + wnxten_[nr] + ynzten_[nr] - y2ny2n_[nr]);
      /* for (nr=0;nr<n_r_rtn;nr++){ } */}
    memcpy(QR_pre_,QR_pos_,n_r*sizeof(float));
    if (verbose>2){ array_printf(etx1n_,"float",1,n_c_rtn," % etx1n_: ");}
    if (verbose>2){ array_printf(etz1n_,"float",1,n_c_rmv," % etz1n_: ");}
    if (verbose>2){ array_printf(ety1n_,"float",1,n_c_rmv," % ety1n_: ");}
    if (verbose>2){ array_printf(ynyten_,"float",1,n_r_rtn," % ynyten_: ");}
    if (verbose>2){ array_printf(ynzten_,"float",1,n_r_rtn," % ynzten_: ");}
    if (verbose>2){ array_printf(y2ny2n_,"float",1,n_r_rtn," % y2ny2n_: ");}
    if (verbose>2){ array_printf(wnxten_,"float",1,n_r_rtn," % wnxten_: ");}
    if (verbose>2){ array_printf(QR_pos_,"float",1,n_r," % QR_pos_: ");}
    if (verbose>2){ array_printf(QR_pre_,"float",1,n_r," % QR_pre_: ");}
    /* %%%%%%%% */
    for (nc=0;nc<n_c_rmv;nc++){ etA1n_[c_rmv_[nc]] = 0; etA2n_[c_rmv_[nc]] = 0; QC_pre_[c_rmv_[nc]] = 0;}
    if (verbose>2){ array_printf(etA1n_,"float",1,n_c," % etA1n_: ");}
    if (verbose>2){ array_printf(etA2n_,"float",1,n_c," % etA2n_: ");}
    if (verbose>2){ array_printf(QC_pre_,"float",1,n_c," % QC_pre_: ");}
    iquicksort_index_driver(n_r_rtn,r_rtn_,1,index_r_rem_,tmp_index_r_); n_index_r_rem = n_r_rtn;
    iquicksort_index_driver(n_c_rtn,c_rtn_,1,index_c_rem_,tmp_index_c_); n_index_c_rem = n_c_rtn;
    if (verbose>2){ array_printf(index_r_rem_,"int",1,n_index_r_rem," % index_r_rem_: ");}
    if (verbose>2){ array_printf(index_c_rem_,"int",1,n_index_c_rem," % index_c_rem_: ");}
    for (nr=0;nr<n_r_rmv;nr++){
      out_xdrop_trn__[ 0 + nx*2 ] = r_rmv_[nr]; out_xdrop_trn__[ 1 + nx*2 ] = -1;
      if (verbose>2){ array_printf(&(out_xdrop_trn__[2*nx]),"int",1,2," % out_xdrop_trn__: ");}
      nx += 1;
      /* for (nr=0;nr<n_r_rmv;nr++){ } */}
    for (nc=0;nc<n_c_rmv;nc++){
      out_xdrop_trn__[ 0 + nx*2 ] = -1; out_xdrop_trn__[ 1 + nx*2 ] = c_rmv_[nc];
      if (verbose>2){ array_printf(&(out_xdrop_trn__[2*nx]),"int",1,2," % out_xdrop_trn__: ");}
      nx += 1;
      /* for (nc=0;nc<n_c_rmv;nc++){ } */}
    if (verbose>2){ printf(" %% nx %d\n",nx);}
    if (verbose>2){ printf(" %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%\n");}
    /* %%%%%%%% */
    /* for (niteration=0;niteration<n_iteration;niteration++){ } */}
  /* %%%%%%%%%%%%%%%% */  
  free1(&index_r_rem_);
  free1(&index_c_rem_);
  free1(&etx1n_);
  free1(&etx2n_);
  free1(&etx1n_stretch_);
  free1(&etw1n_);
  free1(&etw2n_);
  free1(&ety1n_);
  free1(&etz1n_);
  free1(&ynyten_);
  free1(&wnxten_);
  free1(&ynzten_);
  free1(&y2ny2n_);
  free1(&er_);
  free1(&ec_);
  free1(&etA1n_);
  free1(&etA2n_);
  free1(&QR_pre_);
  free1(&QR_pre_r_rem_);
  free1(&QC_pre_);
  free1(&QC_pre_c_rem_);
  free1(&QR_pos_);
  free1(&QC_pos_);
  free1(&QR_workspace_);
  free1(&QC_workspace_);
  free1(&tmp_index_r_);
  free1(&tmp_index_c_);
  free1(&r_rmv_);
  free1(&r_rtn_);
  free1(&c_rmv_);
  free1(&c_rtn_);
  /* %%%% */
  free1(&rdrop_);
  free1(&cdrop_);
  free1(&rkeep_);
  free1(&ckeep_);
  /* %%%%%%%%%%%%%%%% */  
  if (verbose){ printf(" %% [finished dexcluster_nonbinary_f] n_r %d n_c %d gamma %0.2f\n",n_r,n_c,gamma);}
}

void dexcluster_nonbinary_d(int n_r,int n_c,double *B1_rc__,double *B1_cr__,double gamma,int **out_xdrop_p_,int *n_iteration_p,double **trace_p_)
{
  int verbose=0;
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
  double *B1_rc_=NULL,*B1_cr_=NULL;
  double *er_=NULL,*ec_=NULL,*etB1n_=NULL,*etB2n_=NULL;
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
  if (verbose){ printf(" %% [entering dexcluster_nonbinary_d] n_r %d n_c %d gamma %0.2f\n",n_r,n_c,gamma);}
  if (verbose>1){ darray_printf_margin(B1_rc__,n_r,n_c," % B1_rc__: ");printf("\n");}
  if (verbose>1){ darray_printf_margin(B1_cr__,n_c,n_r," % B1_cr__: ");printf("\n");}
  /* %%%%%%%%%%%%%%%% */  
  get_xdrop_logscale_array(n_r,n_c,gamma,&n_iteration,&rdrop_,&cdrop_,&rkeep_,&ckeep_);
  if (n_iteration_p!=NULL){ *n_iteration_p = n_iteration;}
  if (verbose>1){ printf(" %% n_iteration %d\n",n_iteration);}
  if (verbose>1){ array_printf(rdrop_,"int",1,n_iteration," % rdrop_: ");}
  if (verbose>1){ array_printf(cdrop_,"int",1,n_iteration," % cdrop_: ");}
  ulli_r1 = (unsigned long long int)n_r*(unsigned long long int)(1+cdrop_[0]);
  ulli_1c = (unsigned long long int)n_c*(unsigned long long int)(1+rdrop_[0]);
  /* %%%%%%%%%%%%%%%% */  
  if (*out_xdrop_p_==NULL){ (*out_xdrop_p_) = (int *) malloc1(2*(n_r+n_c)*sizeof(int));}
  if (*trace_p_==NULL){ (*trace_p_) = (double *) malloc1(6*n_iteration*sizeof(double));}
  out_xdrop_trn__ = *out_xdrop_p_; trace_trn__ = *trace_p_;
  /* %%%%%%%%%%%%%%%% */  
  n_index_r_rem = n_r; n_index_c_rem = n_c;
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
  etB1n_ = (double *) malloc1(n_c*sizeof(double));
  etB2n_ = (double *) malloc1(n_c*sizeof(double));
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
    B1_rc_ = B1_rc__ + (unsigned long long int)nc*(unsigned long long int)n_r;
    dp_pd_immintrin_loadu(n_r,er_,B1_rc_,&(etB1n_[nc]));
    dp_pd_immintrin_loadu(n_r,B1_rc_,B1_rc_,&(etB2n_[nc]));
    QC_pre_[nc] = etB1n_[nc]*etB1n_[nc] - etB2n_[nc];
    /* for (nc=0;nc<n_c;nc++){ } */}
  for (nr=0;nr<n_r;nr++){
    B1_cr_ = B1_cr__ + (unsigned long long int)nr*(unsigned long long int)n_c;
    dp_pd_immintrin_loadu(n_c,etB1n_,B1_cr_,&(tmp_d_0));
    dp_pd_immintrin_loadu(n_c,B1_cr_,B1_cr_,&(tmp_d_1));
    QR_pre_[nr] = tmp_d_0 - tmp_d_1;
    /* for (nr=0;nr<n_r;nr++){ } */}
  if (verbose>2){ array_printf(etB1n_,"double",1,n_c," % etB1n_: ");}
  if (verbose>2){ array_printf(etB2n_,"double",1,n_c," % etB2n_: ");}
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
      B1_rc_ = B1_rc__ + (unsigned long long int)(index_c_rem_[nc]) * (unsigned long long int)n_r;
      for (nr=0;nr<n_r_rmv;nr++){
	tmp_d_0 = B1_rc_[r_rmv_[nr]];
	etx1n_[nc] += tmp_d_0;
	etx2n_[nc] += tmp_d_0*tmp_d_0;
	/* for (nr=0;nr<n_r_rmv;nr++){ } */}
      etw1n_[nc] = etB1n_[index_c_rem_[nc]] - etx1n_[nc];
      etw2n_[nc] = etB2n_[index_c_rem_[nc]] - etx2n_[nc];
      etB1n_[index_c_rem_[nc]] -= etx1n_[nc];
      etB2n_[index_c_rem_[nc]] -= etx2n_[nc];
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
    n_c_rmv = cdrop_[niteration]; n_c_rtn = n_index_c_rem - n_c_rmv;
    for (nc=0;nc<n_c_rmv;nc++){ c_rmv_[nc] = index_c_rem_[tmp_index_c_[nc]];}
    for (nc=0;nc<n_c_rtn;nc++){ c_rtn_[nc] = index_c_rem_[tmp_index_c_[nc+n_c_rmv]];}
    if (verbose>2){ array_printf(c_rmv_,"int",1,n_c_rmv," % c_rmv_: ");}
    if (verbose>2){ array_printf(c_rtn_,"int",1,n_c_rtn," % c_rtn_: ");}
    /* %%%%%%%% */
    memset(etx1n_,0,n_c*sizeof(double));
    memset(etx1n_stretch_,0,n_c*sizeof(double));
    for (nc=0;nc<n_c_rtn;nc++){
      etx1n_[nc] = 0;
      B1_rc_ = B1_rc__ + (unsigned long long int)(c_rtn_[nc]) * (unsigned long long int)n_r;
      for (nr=0;nr<n_r_rmv;nr++){
	tmp_d_0 = B1_rc_[r_rmv_[nr]];
	etx1n_[nc] += tmp_d_0;
	/* for (nr=0;nr<n_r_rmv;nr++){ } */}
      etx1n_stretch_[c_rtn_[nc]] = etx1n_[nc];
      /* for (nc=0;nc<n_c_rtn;nc++){ } */}
    memset(etz1n_,0,n_c*sizeof(double));
    for (nc=0;nc<n_c_rmv;nc++){
      etz1n_[nc] = 0;
      B1_rc_ = B1_rc__ + (unsigned long long int)(c_rmv_[nc]) * (unsigned long long int)n_r;
      for (nr=0;nr<n_r_rmv;nr++){
	tmp_d_0 = B1_rc_[r_rmv_[nr]];
	etz1n_[nc] += tmp_d_0;
	/* for (nr=0;nr<n_r_rmv;nr++){ } */}
      /* for (nc=0;nc<n_c_rmv;nc++){ } */}
    memset(ety1n_,0,n_c*sizeof(double));
    for (nc=0;nc<n_c_rmv;nc++){
      ety1n_[nc] = 0;
      B1_rc_ = B1_rc__ + (unsigned long long int)(c_rmv_[nc]) * (unsigned long long int)n_r;
      for (nr=0;nr<n_r_rtn;nr++){
	tmp_d_0 = B1_rc_[r_rtn_[nr]];
	ety1n_[nc] += tmp_d_0;
	/* for (nr=0;nr<n_r_rtn;nr++){ } */}
      /* for (nc=0;nc<n_c_rmv;nc++){ } */}
    memset(ynyten_,0,n_r*sizeof(double));
    memset(wnxten_,0,n_r*sizeof(double));
    memset(ynzten_,0,n_r*sizeof(double));
    memset(y2ny2n_,0,n_r*sizeof(double));
    for (nr=0;nr<n_r_rtn;nr++){
      ynyten_[nr] = 0; ynzten_[nr] = 0; y2ny2n_[nr] = 0;
      B1_cr_ = B1_cr__ + (unsigned long long int)(r_rtn_[nr]) * (unsigned long long int)n_c;
      for (nc=0;nc<n_c_rmv;nc++){
	tmp_d_0 = B1_cr_[c_rmv_[nc]];
	ynyten_[nr] += tmp_d_0*ety1n_[nc];
	ynzten_[nr] += tmp_d_0*etz1n_[nc];
	y2ny2n_[nr] += tmp_d_0*tmp_d_0;
	/* for (nc=0;nc<n_c_rmv;nc++){ } */}
      if (n_c_rtn>n_c/4){ dp_pd_immintrin_loadu(n_c,B1_cr_,etx1n_stretch_,&(wnxten_[nr]));}
      else{
	wnxten_[nr] = 0;
	for (nc=0;nc<n_c_rtn;nc++){
	  tmp_d_0 = B1_cr_[c_rtn_[nc]];
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
    for (nc=0;nc<n_c_rmv;nc++){ etB1n_[c_rmv_[nc]] = 0; etB2n_[c_rmv_[nc]] = 0; QC_pre_[c_rmv_[nc]] = 0;}
    if (verbose>2){ array_printf(etB1n_,"double",1,n_c," % etB1n_: ");}
    if (verbose>2){ array_printf(etB2n_,"double",1,n_c," % etB2n_: ");}
    if (verbose>2){ array_printf(QC_pre_,"double",1,n_c," % QC_pre_: ");}
    iquicksort_index_driver(n_r_rtn,r_rtn_,1,index_r_rem_,tmp_index_r_); n_index_r_rem = n_r_rtn;
    iquicksort_index_driver(n_c_rtn,c_rtn_,1,index_c_rem_,tmp_index_c_); n_index_c_rem = n_c_rtn;
    if (verbose>2){ array_printf(index_r_rem_,"int",1,n_index_r_rem," % index_r_rem_: ");}
    if (verbose>2){ array_printf(index_c_rem_,"int",1,n_index_c_rem," % index_c_rem_: ");}
    for (nr=0;nr<n_r_rmv;nr++){
      out_xdrop_trn__[ 0 + nx*2 ] = r_rmv_[nr]; out_xdrop_trn__[ 1 + nx*2 ] = -1;
      if (verbose>2){ array_printf(&(out_xdrop_trn__[2*nx]),"int",1,2," % out_xdrop_trn__: ");}
      nx += 1;
      /* for (nr=0;nr<n_r_rmv;nr++){ } */}
    for (nc=0;nc<n_c_rmv;nc++){
      out_xdrop_trn__[ 0 + nx*2 ] = -1; out_xdrop_trn__[ 1 + nx*2 ] = c_rmv_[nc];
      if (verbose>2){ array_printf(&(out_xdrop_trn__[2*nx]),"int",1,2," % out_xdrop_trn__: ");}
      nx += 1;
      /* for (nc=0;nc<n_c_rmv;nc++){ } */}
    if (verbose>2){ printf(" %% nx %d\n",nx);}
    if (verbose>2){ printf(" %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%\n");}
    /* %%%%%%%% */
    /* for (niteration=0;niteration<n_iteration;niteration++){ } */}
  /* %%%%%%%%%%%%%%%% */  
  free1(&index_r_rem_);
  free1(&index_c_rem_);
  free1(&etx1n_);
  free1(&etx2n_);
  free1(&etx1n_stretch_);
  free1(&etw1n_);
  free1(&etw2n_);
  free1(&ety1n_);
  free1(&etz1n_);
  free1(&ynyten_);
  free1(&wnxten_);
  free1(&ynzten_);
  free1(&y2ny2n_);
  free1(&er_);
  free1(&ec_);
  free1(&etB1n_);
  free1(&etB2n_);
  free1(&QR_pre_);
  free1(&QR_pre_r_rem_);
  free1(&QC_pre_);
  free1(&QC_pre_c_rem_);
  free1(&QR_pos_);
  free1(&QC_pos_);
  free1(&QR_workspace_);
  free1(&QC_workspace_);
  free1(&tmp_index_r_);
  free1(&tmp_index_c_);
  free1(&r_rmv_);
  free1(&r_rtn_);
  free1(&c_rmv_);
  free1(&c_rtn_);
  /* %%%% */
  free1(&rdrop_);
  free1(&cdrop_);
  free1(&rkeep_);
  free1(&ckeep_);
  /* %%%%%%%%%%%%%%%% */  
  if (verbose){ printf(" %% [finished dexcluster_nonbinary_d] n_r %d n_c %d gamma %0.2f\n",n_r,n_c,gamma);}
}

void dexcluster_nonbinary_test_error()
{
  int verbose=1;
  int n_r = 25;
  int n_c = 50;
  float *A1_rc__=NULL,*A1_cr__=NULL;
  double *B1_rc__=NULL,*B1_cr__=NULL;
  double gamma = 0.01;
  int *out_xdrop_trn__=NULL;
  int n_iteration=0;
  double *trace_trn__=NULL;
  int nr=0,nc=0;
  float f=0;
  int out_xdrop_trn_ans_[150] = { 0 , -1 , -1 , 7 , -1 , 2 , -1 , 12 , 7 , -1 , -1 , 46 , -1 , 41 , -1 , 36 , 14 , -1 , -1 , 31 , -1 , 26 , -1 , 21 , 21 , -1 , -1 , 16 , -1 , 11 , -1 , 6 , 1 , -1 , -1 , 1 , -1 , 45 , -1 , 17 , 8 , -1 , -1 , 40 , -1 , 35 , -1 , 30 , 15 , -1 , -1 , 25 , -1 , 20 , 22 , -1 , -1 , 15 , -1 , 10 , 2 , -1 , -1 , 5 , -1 , 0 , 9 , -1 , -1 , 22 , -1 , 27 , 16 , -1 , -1 , 32 , -1 , 37 , 23 , -1 , -1 , 42 , -1 , 47 , 3 , -1 , -1 , 3 , -1 , 8 , 10 , -1 , -1 , 13 , -1 , 18 , 17 , -1 , -1 , 23 , -1 , 28 , 24 , -1 , -1 , 33 , -1 , 38 , 4 , -1 , -1 , 43 , -1 , 48 , 11 , -1 , -1 , 4 , -1 , 9 , 18 , -1 , -1 , 14 , -1 , 19 , 5 , -1 , -1 , 24 , 12 , -1 , -1 , 29 , 19 , -1 , -1 , 34 , 6 , -1 , -1 , 39 , 13 , -1 , 20 , -1 , -1 , 44 , -1 , 49 };
  double trace_trn_ans_[144] = { 1.0000000000000000 , 25.0000000000000000 , 50.0000000000000000 , 2483.6257279999995262 , 1241.8128639999999905 , 1.0000000000000000 , 2.0000000000000000 , 24.0000000000000000 , 47.0000000000000000 , 2376.4593843200004812 , 1188.2296921599997859 , 1.0000000000000000 , 3.0000000000000000 , 23.0000000000000000 , 44.0000000000000000 , 2297.7922294271997998 , 1148.8961147136003547 , 1.0000000000000000 , 4.0000000000000000 , 22.0000000000000000 , 41.0000000000000000 , 2245.6393145855995499 , 1122.8196572928000023 , 1.0000000000000000 , 5.0000000000000000 , 21.0000000000000000 , 38.0000000000000000 , 2218.7453260800002681 , 1109.3726630399996793 , 1.0000000000000000 , 6.0000000000000000 , 20.0000000000000000 , 35.0000000000000000 , 2138.7989967359999355 , 1069.3994983679999677 , 1.0000000000000000 , 7.0000000000000000 , 19.0000000000000000 , 32.0000000000000000 , 2069.1356022784002562 , 1034.5678011391999007 , 1.0000000000000000 , 8.0000000000000000 , 18.0000000000000000 , 30.0000000000000000 , 2021.5997698559999662 , 1010.7998849280000968 , 1.0000000000000000 , 9.0000000000000000 , 17.0000000000000000 , 28.0000000000000000 , 1988.3568463871997665 , 994.1784231936002243 , 1.0000000000000000 , 10.0000000000000000 , 16.0000000000000000 , 26.0000000000000000 , 1881.1977539583999715 , 940.5988769791998720 , 1.0000000000000000 , 11.0000000000000000 , 15.0000000000000000 , 24.0000000000000000 , 1761.0953908223996223 , 880.5476954112001522 , 1.0000000000000000 , 12.0000000000000000 , 14.0000000000000000 , 22.0000000000000000 , 1628.4384197631998177 , 814.2192098815997952 , 1.0000000000000000 , 13.0000000000000000 , 13.0000000000000000 , 20.0000000000000000 , 1484.1892474879996371 , 742.0946237439999322 , 1.0000000000000000 , 14.0000000000000000 , 12.0000000000000000 , 18.0000000000000000 , 1269.1680898559998241 , 634.5840449279999120 , 1.0000000000000000 , 15.0000000000000000 , 11.0000000000000000 , 16.0000000000000000 , 1061.6979671039998721 , 530.8489835519999360 , 1.0000000000000000 , 16.0000000000000000 , 10.0000000000000000 , 14.0000000000000000 , 864.8061715968000271 , 432.4030857984000136 , 1.0000000000000000 , 17.0000000000000000 , 9.0000000000000000 , 12.0000000000000000 , 681.5075844095999855 , 340.7537922048000496 , 1.0000000000000000 , 18.0000000000000000 , 8.0000000000000000 , 10.0000000000000000 , 487.4573020159998578 , 243.7286510080000141 , 1.0000000000000000 , 19.0000000000000000 , 7.0000000000000000 , 8.0000000000000000 , 325.7846609919999423 , 162.8923304960000280 , 1.0000000000000000 , 20.0000000000000000 , 6.0000000000000000 , 6.0000000000000000 , 197.5441136639999229 , 98.7720568320000183 , 1.0000000000000000 , 21.0000000000000000 , 5.0000000000000000 , 5.0000000000000000 , 116.1010931199999163 , 58.0505465600000363 , 1.0000000000000000 , 22.0000000000000000 , 4.0000000000000000 , 4.0000000000000000 , 59.9358410751999173 , 29.9679205376000404 , 1.0000000000000000 , 23.0000000000000000 , 3.0000000000000000 , 3.0000000000000000 , 24.9864972287998164 , 12.4932486144000450 , 1.0000000000000000 , 24.0000000000000000 , 2.0000000000000000 , 2.0000000000000000 , 5.6518589439998301 , 2.8259294719999986 , 1.0000000000000000 };
  if (verbose){ printf(" %% [entering dexcluster_nonbinary_test_error]\n");}
  A1_rc__ = (float *) malloc1(n_r*n_c*sizeof(float));
  A1_cr__ = (float *) malloc1(n_r*n_c*sizeof(float));
  B1_rc__ = (double *) malloc1(n_r*n_c*sizeof(double));
  B1_cr__ = (double *) malloc1(n_r*n_c*sizeof(double));
  for (nr=0;nr<n_r;nr++){ for (nc=0;nc<n_c;nc++){
      f = (float) (nr%7)-3 + (nc%5)-2 + (float)(nr+nc*n_r)/(float)(n_r*n_c);
      A1_rc__[nr+nc*n_r] = f;
      A1_cr__[nc+nr*n_c] = f;
      B1_rc__[nr+nc*n_r] = (double)f;
      B1_cr__[nc+nr*n_c] = (double)f;
      /* for (nr=0;nr<n_r;nr++){ for (nc=0;nc<n_c;nc++){ }} */}}
  /* %%%%%%%% */
  GLOBAL_tic(0);  
  dexcluster_nonbinary_f(n_r,n_c,A1_rc__,A1_cr__,gamma,&out_xdrop_trn__,&n_iteration,&trace_trn__);
  array_printf(trace_trn__,"double",6,minimum(8,n_iteration)," % trace_trn__: ");
  printf(" %% out_xdrop_trn_ans_ vs out_xdrop_trn__: relative error %0.16f (may not be small)\n",ifnormn(2*(n_r+n_c),out_xdrop_trn_ans_,out_xdrop_trn__));
  printf(" %% trace_trn_ans_ vs trace_trn__: relative error %0.16f (should be small)\n",dfnormn(n_iteration,trace_trn_ans_,trace_trn__));
  free1(&out_xdrop_trn__);
  free1(&trace_trn__);
  GLOBAL_toc(0,1," dexcluster_nonbinary_f: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);  
  dexcluster_nonbinary_d(n_r,n_c,B1_rc__,B1_cr__,gamma,&out_xdrop_trn__,&n_iteration,&trace_trn__);
  array_printf(trace_trn__,"double",6,minimum(8,n_iteration)," % trace_trn__: ");
  printf(" %% out_xdrop_trn_ans_ vs out_xdrop_trn__: relative error %0.16f (may not be small)\n",ifnormn(2*(n_r+n_c),out_xdrop_trn_ans_,out_xdrop_trn__));
  printf(" %% trace_trn_ans_ vs trace_trn__: relative error %0.16f (should be small)\n",dfnormn(n_iteration,trace_trn_ans_,trace_trn__));
  free1(&out_xdrop_trn__);
  free1(&trace_trn__);
  GLOBAL_toc(0,1," dexcluster_nonbinary_d: ");
  /* %%%%%%%% */
  free1(&A1_rc__);
  free1(&A1_cr__);
  free1(&B1_rc__);
  free1(&B1_cr__);
  /* %%%%%%%% */
  if (verbose){ printf(" %% [finished dexcluster_nonbinary_test_error]\n");}
}

void dexcluster_nonbinary_test_speed()
{
  int verbose=1;
  int n_r = 1000;
  int n_c = 20000;
  float *A1_rc__=NULL,*A1_cr__=NULL;
  double *B1_rc__=NULL,*B1_cr__=NULL;
  double gamma = 0.01;
  int *out_xdrop_trn__=NULL;
  int n_iteration=0;
  double *trace_trn__=NULL;
  int nr=0,nc=0;
  float f=0;
  if (verbose){ printf(" %% [entering dexcluster_nonbinary_test_speed]\n");}
  A1_rc__ = (float *) malloc1(n_r*n_c*sizeof(float));
  A1_cr__ = (float *) malloc1(n_r*n_c*sizeof(float));
  B1_rc__ = (double *) malloc1(n_r*n_c*sizeof(double));
  B1_cr__ = (double *) malloc1(n_r*n_c*sizeof(double));
  for (nr=0;nr<n_r;nr++){ for (nc=0;nc<n_c;nc++){
      f = (float) (nr%7)-3 + (nc%5)-2 + (float)(nr+nc*n_r)/(float)(n_r*n_c);
      A1_rc__[nr+nc*n_r] = f;
      A1_cr__[nc+nr*n_c] = f;
      B1_rc__[nr+nc*n_r] = (double)f;
      B1_cr__[nc+nr*n_c] = (double)f;
      /* for (nr=0;nr<n_r;nr++){ for (nc=0;nc<n_c;nc++){ }} */}}
  /* %%%%%%%% */
  GLOBAL_tic(0);  
  dexcluster_nonbinary_f(n_r,n_c,A1_rc__,A1_cr__,gamma,&out_xdrop_trn__,&n_iteration,&trace_trn__);
  array_printf(trace_trn__,"double",6,minimum(8,n_iteration)," % trace_trn__: ");
  free1(&out_xdrop_trn__);
  free1(&trace_trn__);
  GLOBAL_toc(0,1," dexcluster_nonbinary_f: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);  
  dexcluster_nonbinary_d(n_r,n_c,B1_rc__,B1_cr__,gamma,&out_xdrop_trn__,&n_iteration,&trace_trn__);
  array_printf(trace_trn__,"double",6,minimum(8,n_iteration)," % trace_trn__: ");
  free1(&out_xdrop_trn__);
  free1(&trace_trn__);
  GLOBAL_toc(0,1," dexcluster_nonbinary_d: ");
  /* %%%%%%%% */
  free1(&A1_rc__);
  free1(&A1_cr__);
  free1(&B1_rc__);
  free1(&B1_cr__);
  /* %%%%%%%% */
  if (verbose){ printf(" %% [finished dexcluster_nonbinary_test_speed]\n");}
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */  

void dexcluster_nonbinary_rdrop_f_bkp(int n_r,int n_c,float *A1_rc__,float *A1_cr__,double gamma,int **out_xdrop_p_,int *n_iteration_p,double **trace_p_)
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
  float *A1_rc_=NULL,*A1_cr_=NULL;
  float *er_=NULL,*ec_=NULL,*etA1n_=NULL,*etA2n_=NULL;
  float *QR_pre_=NULL,*QC_pre_=NULL,*QR_pos_=NULL,*QC_pos_=NULL;
  float *QR_pre_r_rem_=NULL,*QC_pre_c_rem_=NULL;
  float *QR_workspace_=NULL,*QC_workspace_=NULL;
  float *etx1n_=NULL,*etx2n_=NULL,*etx1n_stretch_=NULL;
  float *etw1n_=NULL,*etw2n_=NULL;
  float *ety1n_=NULL;
  float *etz1n_=NULL;
  float *ynyten_=NULL,*wnxten_=NULL;
  float *ynzten_=NULL,*y2ny2n_=NULL;
  int *tmp_index_r_=NULL,*tmp_index_c_=NULL;
  int n_r_rmv=0,n_r_rtn=0;
  int *r_rmv_=NULL,*r_rtn_=NULL;
  int n_c_rmv=0,n_c_rtn=0;
  int *c_rmv_=NULL,*c_rtn_=NULL;
  int nx=0,niteration=0;
  float tmp_f_0=0,tmp_f_1=0;
  if (verbose){ printf(" %% [entering dexcluster_nonbinary_rdrop_f] n_r %d n_c %d gamma %0.2f\n",n_r,n_c,gamma);}
  if (verbose>1){ farray_printf_margin(A1_rc__,n_r,n_c," % A1_rc__: ");printf("\n");}
  if (verbose>1){ farray_printf_margin(A1_cr__,n_c,n_r," % A1_cr__: ");printf("\n");}
  /* %%%%%%%%%%%%%%%% */  
  get_xdrop_logscale_array(n_r,n_c,gamma,&n_iteration,&rdrop_,&cdrop_,&rkeep_,&ckeep_);
  memset(cdrop_,0,niteration*sizeof(int));
  if (n_iteration_p!=NULL){ *n_iteration_p = n_iteration;}
  if (verbose>1){ printf(" %% n_iteration %d\n",n_iteration);}
  if (verbose>1){ array_printf(rdrop_,"int",1,n_iteration," % rdrop_: ");}
  if (verbose>1){ array_printf(cdrop_,"int",1,n_iteration," % cdrop_: ");}
  ulli_r1 = (unsigned long long int)n_r*(unsigned long long int)(1+cdrop_[0]);
  ulli_1c = (unsigned long long int)n_c*(unsigned long long int)(1+rdrop_[0]);
  /* %%%%%%%%%%%%%%%% */  
  if (*out_xdrop_p_==NULL){ (*out_xdrop_p_) = (int *) malloc1(2*n_r*sizeof(int));}
  if (*trace_p_==NULL){ (*trace_p_) = (double *) malloc1(6*n_iteration*sizeof(double));}
  out_xdrop_trn__ = *out_xdrop_p_; trace_trn__ = *trace_p_;
  /* %%%%%%%%%%%%%%%% */  
  n_index_r_rem = n_r; n_index_c_rem = n_c;
  index_r_rem_ = (int *) malloc1(n_index_r_rem*sizeof(int)); for (nr=0;nr<n_r;nr++){ index_r_rem_[nr]=nr;}
  index_c_rem_ = (int *) malloc1(n_index_c_rem*sizeof(int)); for (nc=0;nc<n_c;nc++){ index_c_rem_[nc]=nc;}
  etx1n_ = (float *) malloc1(n_c*sizeof(float));
  etx2n_ = (float *) malloc1(n_c*sizeof(float));
  etx1n_stretch_ = (float *) malloc1(n_c*sizeof(float));
  etw1n_ = (float *) malloc1(n_c*sizeof(float));
  etw2n_ = (float *) malloc1(n_c*sizeof(float));
  ety1n_ = (float *) malloc1(n_c*sizeof(float));
  etz1n_ = (float *) malloc1(n_c*sizeof(float));
  ynyten_ = (float *) malloc1(n_r*sizeof(float));
  wnxten_ = (float *) malloc1(n_r*sizeof(float));
  ynzten_ = (float *) malloc1(n_r*sizeof(float));
  y2ny2n_ = (float *) malloc1(n_r*sizeof(float));
  er_ = (float *) malloc1(n_r*sizeof(float));
  ec_ = (float *) malloc1(n_c*sizeof(float));
  etA1n_ = (float *) malloc1(n_c*sizeof(float));
  etA2n_ = (float *) malloc1(n_c*sizeof(float));
  QR_pre_ = (float *) malloc1(n_r*sizeof(float));
  QR_pre_r_rem_ = (float *) malloc1(n_r*sizeof(float));
  QC_pre_ = (float *) malloc1(n_c*sizeof(float));
  QC_pre_c_rem_ = (float *) malloc1(n_c*sizeof(float));
  QR_pos_ = (float *) malloc1(n_r*sizeof(float));
  QC_pos_ = (float *) malloc1(n_c*sizeof(float));
  QR_workspace_ = (float *) malloc1(n_r*sizeof(float));
  QC_workspace_ = (float *) malloc1(n_c*sizeof(float));
  tmp_index_r_ = (int *) malloc1(n_r*sizeof(int));
  tmp_index_c_ = (int *) malloc1(n_c*sizeof(int));
  r_rmv_ = (int *) malloc1(n_r*sizeof(int));
  r_rtn_ = (int *) malloc1(n_r*sizeof(int));
  c_rmv_ = (int *) malloc1(n_c*sizeof(int));
  c_rtn_ = (int *) malloc1(n_c*sizeof(int));
  /* %%%%%%%%%%%%%%%% */  
  for (nr=0;nr<n_r;nr++){ er_[nr] = (float)1.0;}
  for (nc=0;nc<n_c;nc++){ ec_[nc] = (float)1.0;}
  for (nc=0;nc<n_c;nc++){
    A1_rc_ = A1_rc__ + (unsigned long long int)nc*(unsigned long long int)n_r;
    dp_ps_immintrin_loadu(n_r,er_,A1_rc_,&(etA1n_[nc]));
    dp_ps_immintrin_loadu(n_r,A1_rc_,A1_rc_,&(etA2n_[nc]));
    QC_pre_[nc] = etA1n_[nc]*etA1n_[nc] - etA2n_[nc];
    /* for (nc=0;nc<n_c;nc++){ } */}
  for (nr=0;nr<n_r;nr++){
    A1_cr_ = A1_cr__ + (unsigned long long int)nr*(unsigned long long int)n_c;
    dp_ps_immintrin_loadu(n_c,etA1n_,A1_cr_,&(tmp_f_0));
    dp_ps_immintrin_loadu(n_c,A1_cr_,A1_cr_,&(tmp_f_1));
    QR_pre_[nr] = tmp_f_0 - tmp_f_1;
    /* for (nr=0;nr<n_r;nr++){ } */}
  if (verbose>2){ array_printf(etA1n_,"float",1,n_c," % etA1n_: ");}
  if (verbose>2){ array_printf(etA2n_,"float",1,n_c," % etA2n_: ");}
  if (verbose>2){ array_printf(QR_pre_,"float",1,n_r," % QR_pre_: ");}
  if (verbose>2){ array_printf(QC_pre_,"float",1,n_c," % QC_pre_: ");}
  if (verbose>2){ printf("\n");}
  /* %%%%%%%%%%%%%%%% */  
  nx=0;
  for (niteration=0;niteration<n_iteration;niteration++){
    /* %%%%%%%% */
    trace_trn__[0+6*niteration] = 1+niteration;
    trace_trn__[1+6*niteration] = n_index_r_rem;
    trace_trn__[2+6*niteration] = n_index_c_rem;
    dp_ps_immintrin_loadu(n_r,er_,QR_pre_,&(tmp_f_0));
    trace_trn__[3+6*niteration] = tmp_f_0/n_r;
    dp_ps_immintrin_loadu(n_c,ec_,QC_pre_,&(tmp_f_1));
    trace_trn__[4+6*niteration] = tmp_f_1/n_c;
    trace_trn__[5+6*niteration] = 1.0;
    if (verbose>2){ array_printf(&(trace_trn__[6*niteration]),"double",1,6," % trace_trn__: ");}
    /* %%%%%%%% */
    for (nr=0;nr<n_index_r_rem;nr++){ QR_pre_r_rem_[nr] = QR_pre_[index_r_rem_[nr]];}
    fquicksort_index_driver(n_index_r_rem,QR_pre_r_rem_,1,QR_workspace_,tmp_index_r_);
    if (verbose>2){ array_printf(tmp_index_r_,"int",1,n_index_r_rem," % tmp_index_r_: ");}
    n_r_rmv = rdrop_[niteration]; n_r_rtn = n_index_r_rem - n_r_rmv;
    for (nr=0;nr<n_r_rmv;nr++){ r_rmv_[nr] = index_r_rem_[tmp_index_r_[nr]];}
    for (nr=0;nr<n_r_rtn;nr++){ r_rtn_[nr] = index_r_rem_[tmp_index_r_[nr+n_r_rmv]];}
    if (verbose>2){ array_printf(r_rmv_,"int",1,n_r_rmv," % r_rmv_: ");}
    if (verbose>2){ array_printf(r_rtn_,"int",1,n_r_rtn," % r_rtn_: ");}
    /* %%%%%%%% */
    memset(etx1n_,0,n_c*sizeof(float));
    memset(etx2n_,0,n_c*sizeof(float));
    memset(etw1n_,0,n_c*sizeof(float));
    memset(etw2n_,0,n_c*sizeof(float));
    memset(QC_pos_,0,n_c*sizeof(float));
    for (nc=0;nc<n_index_c_rem;nc++){
      etx1n_[nc] = 0; etx2n_[nc] = 0; etw1n_[nc] = 0; etw2n_[nc] = 0;
      A1_rc_ = A1_rc__ + (unsigned long long int)(index_c_rem_[nc]) * (unsigned long long int)n_r;
      for (nr=0;nr<n_r_rmv;nr++){
	tmp_f_0 = A1_rc_[r_rmv_[nr]];
	etx1n_[nc] += tmp_f_0;
	etx2n_[nc] += tmp_f_0*tmp_f_0;
	/* for (nr=0;nr<n_r_rmv;nr++){ } */}
      etw1n_[nc] = etA1n_[index_c_rem_[nc]] - etx1n_[nc];
      etw2n_[nc] = etA2n_[index_c_rem_[nc]] - etx2n_[nc];
      etA1n_[index_c_rem_[nc]] -= etx1n_[nc];
      etA2n_[index_c_rem_[nc]] -= etx2n_[nc];
      QC_pos_[index_c_rem_[nc]] = QC_pre_[index_c_rem_[nc]] - (2*etw1n_[nc]*etx1n_[nc] + etx1n_[nc]*etx1n_[nc] - etx2n_[nc]);
      /* for (nc=0;nc<n_index_c_rem;nc++){ } */}
    memcpy(QC_pre_,QC_pos_,n_c*sizeof(float));
    if (verbose>2){ array_printf(etx1n_,"float",1,n_index_c_rem," % etx1n_: ");}
    if (verbose>2){ array_printf(etx2n_,"float",1,n_index_c_rem," % etx2n_: ");}
    if (verbose>2){ array_printf(etw1n_,"float",1,n_index_c_rem," % etw1n_: ");}
    if (verbose>2){ array_printf(etw2n_,"float",1,n_index_c_rem," % etw2n_: ");}
    if (verbose>2){ array_printf(QC_pos_,"float",1,n_c," % QC_pos_: ");}
    if (verbose>2){ array_printf(QC_pre_,"float",1,n_c," % QC_pre_: ");}
    /* %%%%%%%% */
    for (nc=0;nc<n_index_c_rem;nc++){ QC_pre_c_rem_[nc] = QC_pre_[index_c_rem_[nc]];}
    fquicksort_index_driver(n_index_c_rem,QC_pre_c_rem_,1,QC_workspace_,tmp_index_c_);
    if (verbose>2){ array_printf(tmp_index_c_,"int",1,n_index_c_rem," % tmp_index_c_: ");}
    n_c_rmv = 0; n_c_rtn = n_index_c_rem;
    memcpy(c_rtn_,index_c_rem_,n_c*sizeof(int));
    /* %%%%%%%% */
    memset(etx1n_,0,n_c*sizeof(float));
    memset(etx1n_stretch_,0,n_c*sizeof(float));
    for (nc=0;nc<n_c_rtn;nc++){
      etx1n_[nc] = 0;
      A1_rc_ = A1_rc__ + (unsigned long long int)(c_rtn_[nc]) * (unsigned long long int)n_r;
      for (nr=0;nr<n_r_rmv;nr++){
	tmp_f_0 = A1_rc_[r_rmv_[nr]];
	etx1n_[nc] += tmp_f_0;
	/* for (nr=0;nr<n_r_rmv;nr++){ } */}
      etx1n_stretch_[c_rtn_[nc]] = etx1n_[nc];
      /* for (nc=0;nc<n_c_rtn;nc++){ } */}
    memset(etz1n_,0,n_c*sizeof(float));
    for (nc=0;nc<n_c_rmv;nc++){
      etz1n_[nc] = 0;
      A1_rc_ = A1_rc__ + (unsigned long long int)(c_rmv_[nc]) * (unsigned long long int)n_r;
      for (nr=0;nr<n_r_rmv;nr++){
	tmp_f_0 = A1_rc_[r_rmv_[nr]];
	etz1n_[nc] += tmp_f_0;
	/* for (nr=0;nr<n_r_rmv;nr++){ } */}
      /* for (nc=0;nc<n_c_rmv;nc++){ } */}
    memset(ety1n_,0,n_c*sizeof(float));
    for (nc=0;nc<n_c_rmv;nc++){
      ety1n_[nc] = 0;
      A1_rc_ = A1_rc__ + (unsigned long long int)(c_rmv_[nc]) * (unsigned long long int)n_r;
      for (nr=0;nr<n_r_rtn;nr++){
	tmp_f_0 = A1_rc_[r_rtn_[nr]];
	ety1n_[nc] += tmp_f_0;
	/* for (nr=0;nr<n_r_rtn;nr++){ } */}
      /* for (nc=0;nc<n_c_rmv;nc++){ } */}
    memset(ynyten_,0,n_r*sizeof(float));
    memset(wnxten_,0,n_r*sizeof(float));
    memset(ynzten_,0,n_r*sizeof(float));
    memset(y2ny2n_,0,n_r*sizeof(float));
    for (nr=0;nr<n_r_rtn;nr++){
      ynyten_[nr] = 0; ynzten_[nr] = 0; y2ny2n_[nr] = 0;
      A1_cr_ = A1_cr__ + (unsigned long long int)(r_rtn_[nr]) * (unsigned long long int)n_c;
      for (nc=0;nc<n_c_rmv;nc++){
	tmp_f_0 = A1_cr_[c_rmv_[nc]];
	ynyten_[nr] += tmp_f_0*ety1n_[nc];
	ynzten_[nr] += tmp_f_0*etz1n_[nc];
	y2ny2n_[nr] += tmp_f_0*tmp_f_0;
	/* for (nc=0;nc<n_c_rmv;nc++){ } */}
      if (n_c_rtn>n_c/4){ dp_ps_immintrin_loadu(n_c,A1_cr_,etx1n_stretch_,&(wnxten_[nr]));}
      else{
	wnxten_[nr] = 0;
	for (nc=0;nc<n_c_rtn;nc++){
	  tmp_f_0 = A1_cr_[c_rtn_[nc]];
	  wnxten_[nr] += tmp_f_0*etx1n_[nc];
	  /* for (nc=0;nc<n_c_rtn;nc++){ } */}
	/* use dp_ps */}
      /* for (nr=0;nr<n_r_rtn;nr++){ } */}
    memset(QR_pos_,0,n_r*sizeof(float));
    for (nr=0;nr<n_r_rtn;nr++){
      QR_pos_[r_rtn_[nr]] = QR_pre_[r_rtn_[nr]] - (ynyten_[nr] + wnxten_[nr] + ynzten_[nr] - y2ny2n_[nr]);
      /* for (nr=0;nr<n_r_rtn;nr++){ } */}
    memcpy(QR_pre_,QR_pos_,n_r*sizeof(float));
    if (verbose>2){ array_printf(etx1n_,"float",1,n_c_rtn," % etx1n_: ");}
    if (verbose>2){ array_printf(etz1n_,"float",1,n_c_rmv," % etz1n_: ");}
    if (verbose>2){ array_printf(ety1n_,"float",1,n_c_rmv," % ety1n_: ");}
    if (verbose>2){ array_printf(ynyten_,"float",1,n_r_rtn," % ynyten_: ");}
    if (verbose>2){ array_printf(ynzten_,"float",1,n_r_rtn," % ynzten_: ");}
    if (verbose>2){ array_printf(y2ny2n_,"float",1,n_r_rtn," % y2ny2n_: ");}
    if (verbose>2){ array_printf(wnxten_,"float",1,n_r_rtn," % wnxten_: ");}
    if (verbose>2){ array_printf(QR_pos_,"float",1,n_r," % QR_pos_: ");}
    if (verbose>2){ array_printf(QR_pre_,"float",1,n_r," % QR_pre_: ");}
    /* %%%%%%%% */
    for (nc=0;nc<n_c_rmv;nc++){ etA1n_[c_rmv_[nc]] = 0; etA2n_[c_rmv_[nc]] = 0; QC_pre_[c_rmv_[nc]] = 0;}
    if (verbose>2){ array_printf(etA1n_,"float",1,n_c," % etA1n_: ");}
    if (verbose>2){ array_printf(etA2n_,"float",1,n_c," % etA2n_: ");}
    if (verbose>2){ array_printf(QC_pre_,"float",1,n_c," % QC_pre_: ");}
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
  free1(&index_r_rem_);
  free1(&index_c_rem_);
  free1(&etx1n_);
  free1(&etx2n_);
  free1(&etx1n_stretch_);
  free1(&etw1n_);
  free1(&etw2n_);
  free1(&ety1n_);
  free1(&etz1n_);
  free1(&ynyten_);
  free1(&wnxten_);
  free1(&ynzten_);
  free1(&y2ny2n_);
  free1(&er_);
  free1(&ec_);
  free1(&etA1n_);
  free1(&etA2n_);
  free1(&QR_pre_);
  free1(&QR_pre_r_rem_);
  free1(&QC_pre_);
  free1(&QC_pre_c_rem_);
  free1(&QR_pos_);
  free1(&QC_pos_);
  free1(&QR_workspace_);
  free1(&QC_workspace_);
  free1(&tmp_index_r_);
  free1(&tmp_index_c_);
  free1(&r_rmv_);
  free1(&r_rtn_);
  free1(&c_rmv_);
  free1(&c_rtn_);
  /* %%%% */
  free1(&rdrop_);
  free1(&cdrop_);
  free1(&rkeep_);
  free1(&ckeep_);
  /* %%%%%%%%%%%%%%%% */  
  if (verbose){ printf(" %% [finished dexcluster_nonbinary_rdrop_f] n_r %d n_c %d gamma %0.2f\n",n_r,n_c,gamma);}
}

void dexcluster_nonbinary_rdrop_f(int n_r,int n_c,float *A1_rc__,float *A1_cr__,double gamma,int **out_xdrop_p_,int *n_iteration_p,double **trace_p_)
{
  int verbose=0;
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
  float *A1_rc_=NULL,*A1_cr_=NULL;
  float *er_=NULL,*etA1n_=NULL;
  float *QR_pre_=NULL,*QR_pos_=NULL;
  float *QR_pre_r_rem_=NULL;
  float *QR_workspace_=NULL;
  float *etx1n_=NULL;
  float *wnxten_=NULL;
  int *tmp_index_r_=NULL,*tmp_index_c_=NULL;
  int n_r_rmv=0,n_r_rtn=0;
  int *r_rmv_=NULL,*r_rtn_=NULL;
  int nx=0,niteration=0;
  float tmp_f_0=0,tmp_f_1=0;
  if (verbose){ printf(" %% [entering dexcluster_nonbinary_rdrop_f] n_r %d n_c %d gamma %0.2f\n",n_r,n_c,gamma);}
  if (verbose>1){ farray_printf_margin(A1_rc__,n_r,n_c," % A1_rc__: ");printf("\n");}
  if (verbose>1){ farray_printf_margin(A1_cr__,n_c,n_r," % A1_cr__: ");printf("\n");}
  /* %%%%%%%%%%%%%%%% */  
  get_xdrop_logscale_array(n_r,n_c,gamma,&n_iteration,&rdrop_,&cdrop_,&rkeep_,&ckeep_);
  memset(cdrop_,0,niteration*sizeof(int));
  if (n_iteration_p!=NULL){ *n_iteration_p = n_iteration;}
  if (verbose>1){ printf(" %% n_iteration %d\n",n_iteration);}
  if (verbose>1){ array_printf(rdrop_,"int",1,n_iteration," % rdrop_: ");}
  if (verbose>1){ array_printf(cdrop_,"int",1,n_iteration," % cdrop_: ");}
  ulli_r1 = (unsigned long long int)n_r*(unsigned long long int)(1+cdrop_[0]);
  ulli_1c = (unsigned long long int)n_c*(unsigned long long int)(1+rdrop_[0]);
  /* %%%%%%%%%%%%%%%% */  
  if (*out_xdrop_p_==NULL){ (*out_xdrop_p_) = (int *) malloc1(2*n_r*sizeof(int));}
  if (*trace_p_==NULL){ (*trace_p_) = (double *) malloc1(6*n_iteration*sizeof(double));}
  out_xdrop_trn__ = *out_xdrop_p_; trace_trn__ = *trace_p_;
  /* %%%%%%%%%%%%%%%% */  
  n_index_r_rem = n_r;
  index_r_rem_ = (int *) malloc1(n_index_r_rem*sizeof(int)); for (nr=0;nr<n_r;nr++){ index_r_rem_[nr]=nr;}
  etx1n_ = (float *) malloc1(n_c*sizeof(float));
  wnxten_ = (float *) malloc1(n_r*sizeof(float));
  er_ = (float *) malloc1(n_r*sizeof(float));
  etA1n_ = (float *) malloc1(n_c*sizeof(float));
  QR_pre_ = (float *) malloc1(n_r*sizeof(float));
  QR_pre_r_rem_ = (float *) malloc1(n_r*sizeof(float));
  QR_pos_ = (float *) malloc1(n_r*sizeof(float));
  QR_workspace_ = (float *) malloc1(n_r*sizeof(float));
  tmp_index_r_ = (int *) malloc1(n_r*sizeof(int));
  r_rmv_ = (int *) malloc1(n_r*sizeof(int));
  r_rtn_ = (int *) malloc1(n_r*sizeof(int));
  /* %%%%%%%%%%%%%%%% */  
  for (nr=0;nr<n_r;nr++){ er_[nr] = (float)1.0;}
  for (nc=0;nc<n_c;nc++){
    A1_rc_ = A1_rc__ + (unsigned long long int)nc*(unsigned long long int)n_r;
    dp_ps_immintrin_loadu(n_r,er_,A1_rc_,&(etA1n_[nc]));
    /* for (nc=0;nc<n_c;nc++){ } */}
  for (nr=0;nr<n_r;nr++){
    A1_cr_ = A1_cr__ + (unsigned long long int)nr*(unsigned long long int)n_c;
    dp_ps_immintrin_loadu(n_c,etA1n_,A1_cr_,&(tmp_f_0));
    dp_ps_immintrin_loadu(n_c,A1_cr_,A1_cr_,&(tmp_f_1));
    QR_pre_[nr] = tmp_f_0 - tmp_f_1;
    /* for (nr=0;nr<n_r;nr++){ } */}
  if (verbose>2){ array_printf(etA1n_,"float",1,n_c," % etA1n_: ");}
  if (verbose>2){ array_printf(QR_pre_,"float",1,n_r," % QR_pre_: ");}
  if (verbose>2){ printf("\n");}
  /* %%%%%%%%%%%%%%%% */  
  nx=0;
  for (niteration=0;niteration<n_iteration;niteration++){
    /* %%%%%%%% */
    trace_trn__[0+6*niteration] = 1+niteration;
    trace_trn__[1+6*niteration] = n_index_r_rem;
    trace_trn__[2+6*niteration] = n_c;
    dp_ps_immintrin_loadu(n_r,er_,QR_pre_,&(tmp_f_0));
    trace_trn__[3+6*niteration] = tmp_f_0/n_r;
    trace_trn__[4+6*niteration] = 0.0;
    trace_trn__[5+6*niteration] = 1.0;
    if (verbose>2){ array_printf(&(trace_trn__[6*niteration]),"double",1,6," % trace_trn__: ");}
    /* %%%%%%%% */
    for (nr=0;nr<n_index_r_rem;nr++){ QR_pre_r_rem_[nr] = QR_pre_[index_r_rem_[nr]];}
    fquicksort_index_driver(n_index_r_rem,QR_pre_r_rem_,1,QR_workspace_,tmp_index_r_);
    if (verbose>2){ array_printf(tmp_index_r_,"int",1,n_index_r_rem," % tmp_index_r_: ");}
    n_r_rmv = rdrop_[niteration]; n_r_rtn = n_index_r_rem - n_r_rmv;
    for (nr=0;nr<n_r_rmv;nr++){ r_rmv_[nr] = index_r_rem_[tmp_index_r_[nr]];}
    for (nr=0;nr<n_r_rtn;nr++){ r_rtn_[nr] = index_r_rem_[tmp_index_r_[nr+n_r_rmv]];}
    if (verbose>2){ array_printf(r_rmv_,"int",1,n_r_rmv," % r_rmv_: ");}
    if (verbose>2){ array_printf(r_rtn_,"int",1,n_r_rtn," % r_rtn_: ");}
    /* %%%%%%%% */
    if (n_r_rmv==1){
      A1_cr_ = A1_cr__ + (unsigned long long int)(r_rmv_[0]) * (unsigned long long int)n_c;
      memcpy(etx1n_,A1_cr_,n_c*sizeof(float));
      /* if (n_r_rmv==1){ } */}
    else /* if (n_r_rmv!=1) */{
      memset(etx1n_,0,n_c*sizeof(float));
      for (nr=0;nr<n_r_rmv;nr++){
	A1_cr_ = A1_cr__ + (unsigned long long int)(r_rmv_[nr]) * (unsigned long long int)n_c;
	for (nc=0;nc<n_c;nc++){
	  etx1n_[nc] += A1_cr_[nc];
	  /* for (nc=0;nc<n_c;nc++){ } */}
	/* for (nr=0;nr<n_r_rmv;nr++){ } */}
      /* else if (n_r_rmv!=1){ } */}
    if (0){ /* old version */
      memset(etx1n_,0,n_c*sizeof(float));
      for (nc=0;nc<n_c;nc++){
	etx1n_[nc] = 0;
	A1_rc_ = A1_rc__ + (unsigned long long int)(nc) * (unsigned long long int)n_r;
	for (nr=0;nr<n_r_rmv;nr++){
	  tmp_f_0 = A1_rc_[r_rmv_[nr]];
	  etx1n_[nc] += tmp_f_0;
	  /* for (nr=0;nr<n_r_rmv;nr++){ } */}
	/* for (nc=0;nc<n_c;nc++){ } */}
      /* if (0){ } */}
    memset(wnxten_,0,n_r*sizeof(float));
    for (nr=0;nr<n_r_rtn;nr++){
      A1_cr_ = A1_cr__ + (unsigned long long int)(r_rtn_[nr]) * (unsigned long long int)n_c;
      dp_ps_immintrin_loadu(n_c,A1_cr_,etx1n_,&(wnxten_[nr]));
      /* for (nr=0;nr<n_r_rtn;nr++){ } */}
    memset(QR_pos_,0,n_r*sizeof(float));
    for (nr=0;nr<n_r_rtn;nr++){
      QR_pos_[r_rtn_[nr]] = QR_pre_[r_rtn_[nr]] - wnxten_[nr];
      /* for (nr=0;nr<n_r_rtn;nr++){ } */}
    memcpy(QR_pre_,QR_pos_,n_r*sizeof(float));
    if (verbose>2){ array_printf(etx1n_,"float",1,n_c," % etx1n_: ");}
    if (verbose>2){ array_printf(wnxten_,"float",1,n_r_rtn," % wnxten_: ");}
    if (verbose>2){ array_printf(QR_pos_,"float",1,n_r," % QR_pos_: ");}
    if (verbose>2){ array_printf(QR_pre_,"float",1,n_r," % QR_pre_: ");}
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
  free1(&index_r_rem_);
  free1(&etx1n_);
  free1(&wnxten_);
  free1(&er_);
  free1(&etA1n_);
  free1(&QR_pre_);
  free1(&QR_pre_r_rem_);
  free1(&QR_pos_);
  free1(&QR_workspace_);
  free1(&tmp_index_r_);
  free1(&r_rmv_);
  free1(&r_rtn_);
  /* %%%% */
  free1(&rdrop_);
  free1(&cdrop_);
  free1(&rkeep_);
  free1(&ckeep_);
  /* %%%%%%%%%%%%%%%% */  
  if (verbose){ printf(" %% [finished dexcluster_nonbinary_rdrop_f] n_r %d n_c %d gamma %0.2f\n",n_r,n_c,gamma);}
}

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
  /* %%%%%%%%%%%%%%%% */  
  if (*out_xdrop_p_==NULL){ (*out_xdrop_p_) = (int *) malloc1(2*n_r*sizeof(int));}
  if (*trace_p_==NULL){ (*trace_p_) = (double *) malloc1(6*n_iteration*sizeof(double));}
  out_xdrop_trn__ = *out_xdrop_p_; trace_trn__ = *trace_p_;
  /* %%%%%%%%%%%%%%%% */  
  n_index_r_rem = n_r; n_index_c_rem = n_c;
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
  free1(&index_r_rem_);
  free1(&index_c_rem_);
  free1(&etx1n_);
  free1(&etx2n_);
  free1(&etx1n_stretch_);
  free1(&etw1n_);
  free1(&etw2n_);
  free1(&ety1n_);
  free1(&etz1n_);
  free1(&ynyten_);
  free1(&wnxten_);
  free1(&ynzten_);
  free1(&y2ny2n_);
  free1(&er_);
  free1(&ec_);
  free1(&etA1n_);
  free1(&etA2n_);
  free1(&QR_pre_);
  free1(&QR_pre_r_rem_);
  free1(&QC_pre_);
  free1(&QC_pre_c_rem_);
  free1(&QR_pos_);
  free1(&QC_pos_);
  free1(&QR_workspace_);
  free1(&QC_workspace_);
  free1(&tmp_index_r_);
  free1(&tmp_index_c_);
  free1(&r_rmv_);
  free1(&r_rtn_);
  free1(&c_rmv_);
  free1(&c_rtn_);
  /* %%%% */
  free1(&rdrop_);
  free1(&cdrop_);
  free1(&rkeep_);
  free1(&ckeep_);
  /* %%%%%%%%%%%%%%%% */  
  if (verbose){ printf(" %% [finished dexcluster_nonbinary_rdrop_d] n_r %d n_c %d gamma %0.2f\n",n_r,n_c,gamma);}
}

void dexcluster_nonbinary_rdrop_d(int n_r,int n_c,double *A1_rc__,double *A1_cr__,double gamma,int **out_xdrop_p_,int *n_iteration_p,double **trace_p_)
{
  int verbose=0;
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
  /* %%%%%%%%%%%%%%%% */  
  if (*out_xdrop_p_==NULL){ (*out_xdrop_p_) = (int *) malloc1(2*n_r*sizeof(int));}
  if (*trace_p_==NULL){ (*trace_p_) = (double *) malloc1(6*n_iteration*sizeof(double));}
  out_xdrop_trn__ = *out_xdrop_p_; trace_trn__ = *trace_p_;
  /* %%%%%%%%%%%%%%%% */  
  n_index_r_rem = n_r;
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
	  etx1n_[nc] += A1_cr_[nc];
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
  free1(&index_r_rem_);
  free1(&etx1n_);
  free1(&wnxten_);
  free1(&er_);
  free1(&etA1n_);
  free1(&QR_pre_);
  free1(&QR_pre_r_rem_);
  free1(&QR_pos_);
  free1(&QR_workspace_);
  free1(&tmp_index_r_);
  free1(&r_rmv_);
  free1(&r_rtn_);
  /* %%%% */
  free1(&rdrop_);
  free1(&cdrop_);
  free1(&rkeep_);
  free1(&ckeep_);
  /* %%%%%%%%%%%%%%%% */  
  if (verbose){ printf(" %% [finished dexcluster_nonbinary_rdrop_d] n_r %d n_c %d gamma %0.2f\n",n_r,n_c,gamma);}
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%; */

void dexcluster_nonbinary_rdrop_test_error()
{
  int verbose=1;
  int n_r = 25;
  int n_c = 50;
  float *A1_rc__=NULL,*A1_cr__=NULL;
  double *B1_rc__=NULL,*B1_cr__=NULL;
  double gamma = 0.01;
  int *out_xdrop_trn__=NULL;
  int n_iteration=0;
  double *trace_trn__=NULL;
  int nr=0,nc=0;
  float f=0;
  int out_xdrop_trn_ans_[50] = { 0 , -1 , 7 , -1 , 14 , -1 , 21 , -1 , 1 , -1 , 8 , -1 , 15 , -1 , 22 , -1 , 2 , -1 , 9 , -1 , 16 , -1 , 23 , -1 , 3 , -1 , 10 , -1 , 17 , -1 , 24 , -1 , 4 , -1 , 11 , -1 , 18 , -1 , 5 , -1 , 12 , -1 , 19 , -1 , 6 , -1 , 13 , -1 , 20 , -1 };
  double trace_trn_ans_[144] = { 1.0000000000000000 , 25.0000000000000000 , 50.0000000000000000 , 2483.6257279999995262 , 1241.8128639999999905 , 1.0000000000000000 , 2.0000000000000000 , 24.0000000000000000 , 50.0000000000000000 , 2366.3089279999994687 , 1183.1544639999999617 , 1.0000000000000000 , 3.0000000000000000 , 23.0000000000000000 , 50.0000000000000000 , 2282.5318054400004257 , 1141.2659027200006676 , 1.0000000000000000 , 4.0000000000000000 , 22.0000000000000000 , 50.0000000000000000 , 2232.1261900800000149 , 1116.0630950400000074 , 1.0000000000000000 , 5.0000000000000000 , 21.0000000000000000 , 50.0000000000000000 , 2214.9242880000001605 , 1107.4621439999998529 , 1.0000000000000000 , 6.0000000000000000 , 20.0000000000000000 , 50.0000000000000000 , 2150.5564940799995384 , 1075.2782470399997692 , 1.0000000000000000 , 7.0000000000000000 , 19.0000000000000000 , 50.0000000000000000 , 2103.4818265600001723 , 1051.7409132799998588 , 1.0000000000000000 , 8.0000000000000000 , 18.0000000000000000 , 50.0000000000000000 , 2073.5993689599999925 , 1036.7996844800002236 , 1.0000000000000000 , 9.0000000000000000 , 17.0000000000000000 , 50.0000000000000000 , 2060.8085811199998716 , 1030.4042905600001632 , 1.0000000000000000 , 10.0000000000000000 , 16.0000000000000000 , 50.0000000000000000 , 1969.1464499200001228 , 984.5732249599998340 , 1.0000000000000000 , 11.0000000000000000 , 15.0000000000000000 , 50.0000000000000000 , 1886.6333708800000295 , 943.3166854400002421 , 1.0000000000000000 , 12.0000000000000000 , 14.0000000000000000 , 50.0000000000000000 , 1813.2356812800001080 , 906.6178406399998266 , 1.0000000000000000 , 13.0000000000000000 , 13.0000000000000000 , 50.0000000000000000 , 1748.9200947200004066 , 874.4600473600002033 , 1.0000000000000000 , 14.0000000000000000 , 12.0000000000000000 , 50.0000000000000000 , 1597.7971136000001025 , 798.8985567999999375 , 1.0000000000000000 , 15.0000000000000000 , 11.0000000000000000 , 50.0000000000000000 , 1455.7815872000001036 , 727.8907935999999381 , 1.0000000000000000 , 16.0000000000000000 , 10.0000000000000000 , 50.0000000000000000 , 1322.9071065600001020 , 661.4535532800000510 , 1.0000000000000000 , 17.0000000000000000 , 9.0000000000000000 , 50.0000000000000000 , 1199.2076390400002310 , 599.6038195200000018 , 1.0000000000000000 , 18.0000000000000000 , 8.0000000000000000 , 50.0000000000000000 , 1004.5341260800000782 , 502.2670630400000391 , 1.0000000000000000 , 19.0000000000000000 , 7.0000000000000000 , 50.0000000000000000 , 827.0289472000000615 , 413.5144736000000307 , 1.0000000000000000 , 20.0000000000000000 , 6.0000000000000000 , 50.0000000000000000 , 666.7929472000001851 , 333.3964736000000357 , 1.0000000000000000 , 21.0000000000000000 , 5.0000000000000000 , 50.0000000000000000 , 468.8390912000002686 , 234.4195455999999922 , 1.0000000000000000 , 22.0000000000000000 , 4.0000000000000000 , 50.0000000000000000 , 304.1831014400002005 , 152.0915507200000150 , 1.0000000000000000 , 23.0000000000000000 , 3.0000000000000000 , 50.0000000000000000 , 172.9930764800000418 , 86.4965382400000493 , 1.0000000000000000 , 24.0000000000000000 , 2.0000000000000000 , 50.0000000000000000 , 57.7428096000000366 , 28.8714048000000645 , 1.0000000000000000 };
  double trace0_trn_ans_[144] = { 1.0000000000000000 , 25.0000000000000000 , 50.0000000000000000 , 2483.6257279999995262 , 0.0000000000000000 , 1.0000000000000000 , 2.0000000000000000 , 24.0000000000000000 , 50.0000000000000000 , 2366.3089279999994687 , 0.0000000000000000 , 1.0000000000000000 , 3.0000000000000000 , 23.0000000000000000 , 50.0000000000000000 , 2282.5318054400004257 , 0.0000000000000000 , 1.0000000000000000 , 4.0000000000000000 , 22.0000000000000000 , 50.0000000000000000 , 2232.1261900800000149 , 0.0000000000000000 , 1.0000000000000000 , 5.0000000000000000 , 21.0000000000000000 , 50.0000000000000000 , 2214.9242880000001605 , 0.0000000000000000 , 1.0000000000000000 , 6.0000000000000000 , 20.0000000000000000 , 50.0000000000000000 , 2150.5564940799995384 , 0.0000000000000000 , 1.0000000000000000 , 7.0000000000000000 , 19.0000000000000000 , 50.0000000000000000 , 2103.4818265600001723 , 0.0000000000000000 , 1.0000000000000000 , 8.0000000000000000 , 18.0000000000000000 , 50.0000000000000000 , 2073.5993689599999925 , 0.0000000000000000 , 1.0000000000000000 , 9.0000000000000000 , 17.0000000000000000 , 50.0000000000000000 , 2060.8085811199998716 , 0.0000000000000000 , 1.0000000000000000 , 10.0000000000000000 , 16.0000000000000000 , 50.0000000000000000 , 1969.1464499200001228 , 0.0000000000000000 , 1.0000000000000000 , 11.0000000000000000 , 15.0000000000000000 , 50.0000000000000000 , 1886.6333708800000295 , 0.0000000000000000 , 1.0000000000000000 , 12.0000000000000000 , 14.0000000000000000 , 50.0000000000000000 , 1813.2356812800001080 , 0.0000000000000000 , 1.0000000000000000 , 13.0000000000000000 , 13.0000000000000000 , 50.0000000000000000 , 1748.9200947200004066 , 0.0000000000000000 , 1.0000000000000000 , 14.0000000000000000 , 12.0000000000000000 , 50.0000000000000000 , 1597.7971136000001025 , 0.0000000000000000 , 1.0000000000000000 , 15.0000000000000000 , 11.0000000000000000 , 50.0000000000000000 , 1455.7815872000001036 , 0.0000000000000000 , 1.0000000000000000 , 16.0000000000000000 , 10.0000000000000000 , 50.0000000000000000 , 1322.9071065600001020 , 0.0000000000000000 , 1.0000000000000000 , 17.0000000000000000 , 9.0000000000000000 , 50.0000000000000000 , 1199.2076390400002310 , 0.0000000000000000 , 1.0000000000000000 , 18.0000000000000000 , 8.0000000000000000 , 50.0000000000000000 , 1004.5341260800000782 , 0.0000000000000000 , 1.0000000000000000 , 19.0000000000000000 , 7.0000000000000000 , 50.0000000000000000 , 827.0289472000000615 , 0.0000000000000000 , 1.0000000000000000 , 20.0000000000000000 , 6.0000000000000000 , 50.0000000000000000 , 666.7929472000001851 , 0.0000000000000000 , 1.0000000000000000 , 21.0000000000000000 , 5.0000000000000000 , 50.0000000000000000 , 468.8390912000002686 , 0.0000000000000000 , 1.0000000000000000 , 22.0000000000000000 , 4.0000000000000000 , 50.0000000000000000 , 304.1831014400002005 , 0.0000000000000000 , 1.0000000000000000 , 23.0000000000000000 , 3.0000000000000000 , 50.0000000000000000 , 172.9930764800000418 , 0.0000000000000000 , 1.0000000000000000 , 24.0000000000000000 , 2.0000000000000000 , 50.0000000000000000 , 57.7428096000000366 , 0.0000000000000000 , 1.0000000000000000 };
  if (verbose){ printf(" %% [entering dexcluster_nonbinary_rdrop_test_error]\n");}
  A1_rc__ = (float *) malloc1(n_r*n_c*sizeof(float));
  A1_cr__ = (float *) malloc1(n_r*n_c*sizeof(float));
  B1_rc__ = (double *) malloc1(n_r*n_c*sizeof(double));
  B1_cr__ = (double *) malloc1(n_r*n_c*sizeof(double));
  for (nr=0;nr<n_r;nr++){ for (nc=0;nc<n_c;nc++){
      f = (float) (nr%7)-3 + (nc%5)-2 + (float)(nr+nc*n_r)/(float)(n_r*n_c);
      A1_rc__[nr+nc*n_r] = f;
      A1_cr__[nc+nr*n_c] = f;
      B1_rc__[nr+nc*n_r] = (double)f;
      B1_cr__[nc+nr*n_c] = (double)f;
      /* for (nr=0;nr<n_r;nr++){ for (nc=0;nc<n_c;nc++){ }} */}}
  /* %%%%%%%% */
  GLOBAL_tic(0);  
  dexcluster_nonbinary_rdrop_f_bkp(n_r,n_c,A1_rc__,A1_cr__,gamma,&out_xdrop_trn__,&n_iteration,&trace_trn__);
  array_printf(trace_trn__,"double",6,minimum(8,n_iteration)," % trace_trn__: ");
  printf(" %% out_xdrop_trn_ans_ vs out_xdrop_trn__: relative error %0.16f (may not be small)\n",ifnormn(2*n_r,out_xdrop_trn_ans_,out_xdrop_trn__));
  printf(" %% trace_trn_ans_ vs trace_trn__: relative error %0.16f (should be small)\n",dfnormn(n_iteration,trace_trn_ans_,trace_trn__));
  free1(&out_xdrop_trn__);
  free1(&trace_trn__);
  GLOBAL_toc(0,1," dexcluster_nonbinary_rdrop_f_bkp: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);  
  dexcluster_nonbinary_rdrop_f(n_r,n_c,A1_rc__,A1_cr__,gamma,&out_xdrop_trn__,&n_iteration,&trace_trn__);
  array_printf(trace_trn__,"double",6,minimum(8,n_iteration)," % trace_trn__: ");
  printf(" %% out_xdrop_trn_ans_ vs out_xdrop_trn__: relative error %0.16f (may not be small)\n",ifnormn(2*n_r,out_xdrop_trn_ans_,out_xdrop_trn__));
  printf(" %% trace0_trn_ans_ vs trace_trn__: relative error %0.16f (should be small)\n",dfnormn(n_iteration,trace0_trn_ans_,trace_trn__));
  free1(&out_xdrop_trn__);
  free1(&trace_trn__);
  GLOBAL_toc(0,1," dexcluster_nonbinary_rdrop_f: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);  
  dexcluster_nonbinary_rdrop_d_bkp(n_r,n_c,B1_rc__,B1_cr__,gamma,&out_xdrop_trn__,&n_iteration,&trace_trn__);
  array_printf(trace_trn__,"double",6,minimum(8,n_iteration)," % trace_trn__: ");
  printf(" %% out_xdrop_trn_ans_ vs out_xdrop_trn__: relative error %0.16f (may not be small)\n",ifnormn(2*n_r,out_xdrop_trn_ans_,out_xdrop_trn__));
  printf(" %% trace_trn_ans_ vs trace_trn__: relative error %0.16f (should be small)\n",dfnormn(n_iteration,trace_trn_ans_,trace_trn__));
  free1(&out_xdrop_trn__);
  free1(&trace_trn__);
  GLOBAL_toc(0,1," dexcluster_nonbinary_rdrop_d_bkp: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);  
  dexcluster_nonbinary_rdrop_d(n_r,n_c,B1_rc__,B1_cr__,gamma,&out_xdrop_trn__,&n_iteration,&trace_trn__);
  array_printf(trace_trn__,"double",6,minimum(8,n_iteration)," % trace_trn__: ");
  printf(" %% out_xdrop_trn_ans_ vs out_xdrop_trn__: relative error %0.16f (may not be small)\n",ifnormn(2*n_r,out_xdrop_trn_ans_,out_xdrop_trn__));
  printf(" %% trace0_trn_ans_ vs trace_trn__: relative error %0.16f (should be small)\n",dfnormn(n_iteration,trace0_trn_ans_,trace_trn__));
  free1(&out_xdrop_trn__);
  free1(&trace_trn__);
  GLOBAL_toc(0,1," dexcluster_nonbinary_rdrop_d: ");
  /* %%%%%%%% */
  free1(&A1_rc__);
  free1(&A1_cr__);
  free1(&B1_rc__);
  free1(&B1_cr__);
  /* %%%%%%%% */
  if (verbose){ printf(" %% [finished dexcluster_nonbinary_rdrop_test_error]\n");}
}

void dexcluster_nonbinary_rdrop_test_speed()
{
  int verbose=1;
  int n_r = 1000;
  int n_c = 20000;
  float *A1_rc__=NULL,*A1_cr__=NULL;
  double *B1_rc__=NULL,*B1_cr__=NULL;
  double gamma = 0.01;
  int *out_xdrop_trn__=NULL;
  int n_iteration=0;
  double *trace_trn__=NULL;
  int nr=0,nc=0;
  float f=0;
  if (verbose){ printf(" %% [entering dexcluster_nonbinary_rdrop_test_speed]\n");}
  A1_rc__ = (float *) malloc1(n_r*n_c*sizeof(float));
  A1_cr__ = (float *) malloc1(n_r*n_c*sizeof(float));
  B1_rc__ = (double *) malloc1(n_r*n_c*sizeof(double));
  B1_cr__ = (double *) malloc1(n_r*n_c*sizeof(double));
  for (nr=0;nr<n_r;nr++){ for (nc=0;nc<n_c;nc++){
      f = (float) (nr%7)-3 + (nc%5)-2 + (float)(nr+nc*n_r)/(float)(n_r*n_c);
      A1_rc__[nr+nc*n_r] = f;
      A1_cr__[nc+nr*n_c] = f;
      B1_rc__[nr+nc*n_r] = (double)f;
      B1_cr__[nc+nr*n_c] = (double)f;
      /* for (nr=0;nr<n_r;nr++){ for (nc=0;nc<n_c;nc++){ }} */}}
  /* %%%%%%%% */
  GLOBAL_tic(0);  
  dexcluster_nonbinary_rdrop_f_bkp(n_r,n_c,A1_rc__,A1_cr__,gamma,&out_xdrop_trn__,&n_iteration,&trace_trn__);
  array_printf(trace_trn__,"double",6,minimum(8,n_iteration)," % trace_trn__: ");
  free1(&out_xdrop_trn__);
  free1(&trace_trn__);
  GLOBAL_toc(0,1," dexcluster_nonbinary_rdrop_f_bkp: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);  
  dexcluster_nonbinary_rdrop_f(n_r,n_c,A1_rc__,A1_cr__,gamma,&out_xdrop_trn__,&n_iteration,&trace_trn__);
  array_printf(trace_trn__,"double",6,minimum(8,n_iteration)," % trace_trn__: ");
  free1(&out_xdrop_trn__);
  free1(&trace_trn__);
  GLOBAL_toc(0,1," dexcluster_nonbinary_rdrop_f: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);  
  dexcluster_nonbinary_rdrop_d_bkp(n_r,n_c,B1_rc__,B1_cr__,gamma,&out_xdrop_trn__,&n_iteration,&trace_trn__);
  array_printf(trace_trn__,"double",6,minimum(8,n_iteration)," % trace_trn__: ");
  free1(&out_xdrop_trn__);
  free1(&trace_trn__);
  GLOBAL_toc(0,1," dexcluster_nonbinary_rdrop_d_bkp: ");
  /* %%%%%%%% */
  GLOBAL_tic(0);  
  dexcluster_nonbinary_rdrop_d(n_r,n_c,B1_rc__,B1_cr__,gamma,&out_xdrop_trn__,&n_iteration,&trace_trn__);
  array_printf(trace_trn__,"double",6,minimum(8,n_iteration)," % trace_trn__: ");
  free1(&out_xdrop_trn__);
  free1(&trace_trn__);
  GLOBAL_toc(0,1," dexcluster_nonbinary_rdrop_d: ");
  /* %%%%%%%% */
  free1(&A1_rc__);
  free1(&A1_cr__);
  free1(&B1_rc__);
  free1(&B1_cr__);
  /* %%%%%%%% */
  if (verbose){ printf(" %% [finished dexcluster_nonbinary_rdrop_test_speed]\n");}
}
