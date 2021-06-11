#ifndef _MONOLITH
#include "playpark_header.h"
#endif /* _MONOLITH */

void dexcluster_nonbinary_f_recursive_helper_QR__
(
  int verbose
 ,int flag_rdrop_vs_rcdrop
 ,int flag_force_create
 ,int n_r
 ,int n_c
 ,float * E_base_rc__
 ,int n_r_index
 ,int *r_index_
 ,int n_c_index
 ,int *c_index_
 ,double gamma
 ,int n_shuffle
 ,char *fname_trace__
 ,char *fname_xdrop__
 ,char *fname_QR__
)
{
  struct stat stat_file = {0};
  int flag_not_exist=0;
  int trace_QR_index = 3;
  int n_iteration=0,tmp_n_iteration=0;
  int nshuffle=0;
  unsigned long int rseed=0;
  double *trace__=NULL;
  double *trace_shuffle__=NULL;
  double *QR_=NULL;
  double *QR__=NULL;
  int n_xdrop=0;
  int *xdrop__=NULL;
  int *xdrop_shuffle__=NULL;
  float *E_cr__=NULL;
  float *E_rc__=NULL;
  float *Q_rc__=NULL;
  float *QE_rc__=NULL;
  int MDA_n_dim=0;
  int *MDA_dim_=NULL;
  char MDA_fname[PNAMESIZE];
  if (verbose){ printf(" %% [entering dexcluster_nonbinary_f_recursive_helper_QR__]\n");}
  flag_not_exist = (stat(fname_trace__,&stat_file)==-1) || (stat(fname_QR__,&stat_file)==-1) || (stat(fname_xdrop__,&stat_file)==-1);
  if ( flag_force_create || flag_not_exist ){
    if (verbose){ printf(" %% %s not found, creating\n",fname_trace__);}
    MDA_dim_ = (int *) malloc1(2*sizeof(int));
    n_iteration = get_xdrop_logscale_length(n_r_index,n_c_index,gamma);
    trace__ = (double *) malloc1((unsigned long long int)6*(unsigned long long int)n_r*sizeof(double));
    trace_shuffle__ = (double *) malloc1((unsigned long long int)6*(unsigned long long int)n_r*sizeof(double));
    QR__ = (double *) malloc1((unsigned long long int)n_iteration*(unsigned long long int)(1+n_shuffle)*sizeof(double));
    if (flag_rdrop_vs_rcdrop==0){ n_xdrop = n_r_index + n_c_index;}
    if (flag_rdrop_vs_rcdrop==1){ n_xdrop = n_r_index;}
    xdrop__ = (int *) malloc1((unsigned long long int)2*(unsigned long long int)n_xdrop*sizeof(double));
    xdrop_shuffle__ = (int *) malloc1((unsigned long long int)2*(unsigned long long int)n_xdrop*sizeof(double));
    E_rc__ = (float *) malloc1((unsigned long long int)n_r_index*(unsigned long long int)n_c_index*sizeof(float));
    E_cr__ = (float *) malloc1((unsigned long long int)n_c_index*(unsigned long long int)n_r_index*sizeof(float));
    array_extract_f_from_f(n_r,n_c,E_base_rc__,n_r_index,r_index_,n_c_index,c_index_,&E_rc__,NULL);
    if (verbose>9){ array_printf_margin(E_rc__,"float",n_r_index,n_c_index," % E_rc__: "); printf(" %% %% %% %%\n");}
    array_mean_center_row(n_r_index,n_c_index,E_rc__,NULL,"float",&E_rc__,&E_cr__);
    if (verbose>9){ array_printf_margin(E_rc__,"float",n_r_index,n_c_index," % E_rc__: "); printf(" %% %% %% %%\n");}
    if (verbose>9){ array_printf_margin(E_cr__,"float",n_c_index,n_r_index," % E_cr__: "); printf(" %% %% %% %%\n");}
    if (flag_rdrop_vs_rcdrop==0){ dexcluster_nonbinary_f(n_r_index,n_c_index,E_rc__,E_cr__,gamma,&xdrop__,&tmp_n_iteration,&trace__);}
    if (flag_rdrop_vs_rcdrop==1){ dexcluster_nonbinary_rdrop_f(n_r_index,n_c_index,E_rc__,E_cr__,gamma,&xdrop__,&tmp_n_iteration,&trace__);}
    if (verbose>9){ array_printf_margin(xdrop__,"int",2,n_xdrop," % xdrop__: "); printf(" %% %% %% %%\n");}
    if (verbose>9){ array_printf_margin(trace__,"double",6,n_iteration," % trace__: "); printf(" %% %% %% %%\n");}
    if (tmp_n_iteration!=n_iteration){ printf(" %% Warning, tmp_n_iteration %d not equal n_iteration %d\n",tmp_n_iteration,n_iteration);}
    sprintf(MDA_fname,"%s",fname_trace__); MDA_n_dim = 2; MDA_dim_[0] = 6; MDA_dim_[1] = n_iteration;
    MDA_write_r8(MDA_n_dim,MDA_dim_,trace__,MDA_fname);
    if (verbose>9){ MDA_printf_r8_margin(MDA_fname);}
    sprintf(MDA_fname,"%s",fname_xdrop__); MDA_n_dim = 2; MDA_dim_[0] = 2; MDA_dim_[1] = n_xdrop;
    MDA_write_i4(MDA_n_dim,MDA_dim_,xdrop__,MDA_fname);
    if (verbose>9){ MDA_printf_i4_margin(MDA_fname);}
    QR_ = QR__ + (unsigned long long int)0*(unsigned long long int)n_iteration;
    array_extract_d_from_d(6,n_iteration,trace__,1,&trace_QR_index,0,NULL,&QR_,NULL);
    if (verbose>9){ array_printf_margin(QR_,"double",1,n_iteration," % QR_: ");}
    Q_rc__ = (float *) malloc1((unsigned long long int)n_r_index*(unsigned long long int)n_r_index*sizeof(float));
    QE_rc__ = (float *) malloc1((unsigned long long int)n_r_index*(unsigned long long int)n_c_index*sizeof(float));
    for (nshuffle=0;nshuffle<n_shuffle;nshuffle++){
      if (verbose>9){ printf(" %% nshuffle %d/%d\n",nshuffle,n_shuffle);}
      rseed = (1+nshuffle); RSEED_adv8(&rseed);
      array_extract_f_from_f(n_r,n_c,E_base_rc__,n_r_index,r_index_,n_c_index,c_index_,&E_rc__,NULL);
      array_mean_center_row(n_r_index,n_c_index,E_rc__,NULL,"float",&E_rc__,NULL);
      if (verbose>9){ array_printf_margin(E_rc__,"float",n_r_index,n_c_index," % E_rc__: "); printf(" %% %% %% %%\n");}
      array_orth_f(n_r_index,n_r_index,&Q_rc__,&rseed);
      if (verbose>9){ array_printf_margin(Q_rc__,"float",n_r_index,n_r_index," % Q_rc__: "); printf(" %% %% %% %%\n");}
      dp_ps_mult_immintrin_loadu(n_r_index,n_r_index,Q_rc__,n_c_index,E_rc__,&QE_rc__);
      if (verbose>9){ array_printf_margin(QE_rc__,"float",n_r_index,n_c_index," % QE_rc__: "); printf(" %% %% %% %%\n");}
      array_mean_center_row(n_r_index,n_c_index,QE_rc__,NULL,"float",&E_rc__,&E_cr__);
      if (verbose>9){ array_printf_margin(E_rc__,"float",n_r_index,n_c_index," % E_rc__: "); printf(" %% %% %% %%\n");}
      if (flag_rdrop_vs_rcdrop==0){ dexcluster_nonbinary_f(n_r_index,n_c_index,E_rc__,E_cr__,gamma,&xdrop_shuffle__,&tmp_n_iteration,&trace_shuffle__);}
      if (flag_rdrop_vs_rcdrop==1){ dexcluster_nonbinary_rdrop_f(n_r_index,n_c_index,E_rc__,E_cr__,gamma,&xdrop_shuffle__,&tmp_n_iteration,&trace_shuffle__);}
      if (verbose>9){ array_printf_margin(xdrop_shuffle__,"int",2,n_xdrop," % xdrop_shuffle__: "); printf(" %% %% %% %%\n");}
      if (verbose>9){ array_printf_margin(trace_shuffle__,"double",6,n_iteration," % trace_shuffle__: "); printf(" %% %% %% %%\n");}
      if (tmp_n_iteration!=n_iteration){ printf(" %% Warning, tmp_n_iteration %d not equal n_iteration %d\n",tmp_n_iteration,n_iteration);}
      QR_ = QR__ + (unsigned long long int)(1+nshuffle)*(unsigned long long int)n_iteration;
      array_extract_d_from_d(6,n_iteration,trace_shuffle__,1,&trace_QR_index,0,NULL,&QR_,NULL);
      if (verbose>9){ array_printf_margin(QR_,"double",1,n_iteration," % QR_: ");}
      /* for (nshuffle=0;nshuffle<n_shuffle;nshuffle++){ } */}
    if (verbose>9){ array_printf_margin(QR__,"double",n_iteration,(1+n_shuffle)," % QR__: ");}
    sprintf(MDA_fname,"%s",fname_QR__); MDA_n_dim = 2; MDA_dim_[0] = n_iteration; MDA_dim_[1] = (1+n_shuffle);
    MDA_write_r8(MDA_n_dim,MDA_dim_,QR__,MDA_fname);
    if (verbose>9){ MDA_printf_r8_margin(MDA_fname);}
    free1(&Q_rc__);
    free1(&QE_rc__);
    free1(&E_cr__);
    free1(&E_rc__);
    free1(&trace__);
    free1(&trace_shuffle__);
    free1(&QR__);
    free1(&xdrop__);
    free1(&xdrop_shuffle__);
    free1(&MDA_dim_);
    /* not found */}
  else{ if (verbose){ printf(" %% %s found, not creating\n",fname_trace__);}}
  if (verbose){ printf(" %% [finished dexcluster_nonbinary_f_recursive_helper_QR__]\n");}
}

void dexcluster_nonbinary_recursive_helper_ZR__
(
 int verbose
 ,int flag_rdrop_vs_rcdrop
 ,char *fname_trace__
 ,char *fname_xdrop__
 ,char *fname_QR__
 ,double p_use
 ,int n_member_lob
 ,double *nlp_ZR_max_p
 ,int *nlp_ZR_index_p
 ,int *n_r_rtn_index_p
 ,int **r_rtn_index_p_
 ,int *n_r_rmv_index_p
 ,int **r_rmv_index_p_
 ,int *n_c_rtn_index_p
 ,int **c_rtn_index_p_
 ,int *n_c_rmv_index_p
 ,int **c_rmv_index_p_
 ,double *nlp_gumb_opt_p
 ,double *nlp_gumb_emp_p
)
{
  int trace_r_rtn_index = 1;
  int trace_c_rtn_index = 2;
  int MDA_n_dim=0;
  int *MDA_dim_=NULL;
  int niteration=0,n_iteration=0;
  double *trace__=NULL;
  int nshuffle=0,n_shuffle=0;
  double QR_avg=0,QR_std=0,*QR_=NULL,*QR__=NULL;
  int n_xdrop=0,*xdrop_=NULL,*xdrop__=NULL;
  double ZR=0,nlp_ZR=0,*nlp_ZR_=NULL;
  int nsub=0,n_sub=0,*index_sub_=NULL;
  double *nlp_ZR_sub_=NULL;
  double nlp_ZR_max=0;
  int nlp_ZR_index;
  int r_rtn_max=0,*r_rtn_=NULL;
  int c_rtn_max=0,*c_rtn_=NULL;
  int n_r_rtn_index=0,*r_rtn_index_=NULL;
  int n_r_rmv_index=0,*r_rmv_index_=NULL;
  int n_c_rtn_index=0,*c_rtn_index_=NULL;
  int n_c_rmv_index=0,*c_rmv_index_=NULL;
  double *ZR__=NULL;
  double *ZR_sub_max_s_=NULL;
  double nlp_gumb_opt=0,nlp_gumb_emp=0,p_gumb_opt=0,p_gumb_emp=0;
  int nl=0,n_l=0;
  MDA_dim_ = (int *) malloc1(2*sizeof(int));
  if (verbose){ printf(" %% [entering dexcluster_nonbinary_recursive_helper_ZR__]\n");}
  /* %%%%%%%% */
  MDA_read_r8(&MDA_n_dim,&MDA_dim_,&trace__,fname_trace__);
  if (MDA_n_dim!=2){ printf(" %% Warning, improper n_dim %d in fanme_trace__\n",MDA_n_dim);}
  if (MDA_dim_[0]!=6){ printf(" %% Warning, improper dim_[0] %d in fname_trace__\n",MDA_dim_[0]);}
  n_iteration = MDA_dim_[1];
  if (verbose>2){ array_printf_margin(trace__,"double",6,n_iteration," % trace__: ");}
  r_rtn_ = (int *) malloc1((unsigned long long int)n_iteration*sizeof(int));
  for (niteration=0;niteration<n_iteration;niteration++){ r_rtn_[niteration] = (int)round(trace__[trace_r_rtn_index + niteration*6]);}
  r_rtn_max = r_rtn_[0];
  c_rtn_ = (int *) malloc1((unsigned long long int)n_iteration*sizeof(int));
  for (niteration=0;niteration<n_iteration;niteration++){ c_rtn_[niteration] = (int)round(trace__[trace_c_rtn_index + niteration*6]);}
  c_rtn_max = c_rtn_[0];
  if (flag_rdrop_vs_rcdrop==0){ n_xdrop = r_rtn_max + c_rtn_max;}
  if (flag_rdrop_vs_rcdrop==1){ n_xdrop = r_rtn_max;}
  MDA_read_r8(&MDA_n_dim,&MDA_dim_,&QR__,fname_QR__);
  if (MDA_n_dim!=2){ printf(" %% Warning, improper n_dim %d in fname_QR__\n",MDA_n_dim);}
  if (MDA_dim_[0]!=n_iteration){ printf(" %% Warning, improper dim_[0] %d in fname_QR__\n",MDA_dim_[0]);}
  n_shuffle = MDA_dim_[1]-1;
  if (verbose>2){ array_printf_margin(QR__,"double",n_iteration,1+n_shuffle," % QR__: ");}
  MDA_read_i4(&MDA_n_dim,&MDA_dim_,&xdrop__,fname_xdrop__);
  if (MDA_n_dim!=2){ printf(" %% Warning, improper n_dim %d in fname_xdrop__\n",MDA_n_dim);}
  if (MDA_dim_[0]!=2){ printf(" %% Warning, improper dim_[0] %d in fname_xdrop__\n",MDA_dim_[0]);}
  if (MDA_dim_[1]!=n_xdrop){ printf(" %% Warning, improper dim_[1] %d in fname_xdrop__\n",MDA_dim_[1]);}
  /* %%%%%%%% */
  QR_ = QR__;
  ZR__ = (double *) malloc1((unsigned long long int)n_iteration*(unsigned long long int)(1+n_shuffle)*sizeof(double));
  ZR_sub_max_s_ = (double *) malloc1((unsigned long long int)(1+n_shuffle)*sizeof(double));
  if (verbose>2){ array_printf_margin(QR_,"double",1,n_iteration," % QR_: ");}
  nlp_ZR_ = (double *) malloc1((unsigned long long int)n_iteration*sizeof(double));
  for (niteration=0;niteration<n_iteration;niteration++){
    QR_avg = 0; for (nshuffle=1;nshuffle<1+n_shuffle;nshuffle++){ QR_avg += QR__[niteration+nshuffle*n_iteration];}
    QR_avg /= maximum(1,n_shuffle);
    QR_std = 0; for (nshuffle=1;nshuffle<1+n_shuffle;nshuffle++){ QR_std += pow(QR__[niteration+nshuffle*n_iteration] - QR_avg,2);}
    QR_std /= maximum(1,n_shuffle); QR_std = sqrt(QR_std);
    ZR = (QR_[niteration] - QR_avg)/maximum(1e-12,QR_std);
    nlp_ZR = -z_to_lp_single_d(ZR);
    nlp_ZR_[niteration] = nlp_ZR;
    for (nshuffle=0;nshuffle<1+n_shuffle;nshuffle++){
      ZR = (QR__[niteration + nshuffle*n_iteration] - QR_avg)/maximum(1e-12,QR_std);
      ZR__[niteration + nshuffle*n_iteration] = ZR;
      /* for (nshuffle=0;nshuffle<1+n_shuffle;nshuffle++){ } */}
    /* for (niteration=0;niteration<n_iteration;niteration++){ } */}
  if (verbose>2){ array_printf_margin(nlp_ZR_,"double",1,n_iteration," % nlp_ZR_: ");}
  index_sub_ = (int *) malloc1((unsigned long long int)n_iteration*sizeof(int));
  nlp_ZR_sub_ = (double *) malloc1((unsigned long long int)n_iteration*sizeof(double));
  n_sub=0;
  for (niteration=0;niteration<n_iteration;niteration++){
    if ( (r_rtn_[niteration]>=n_member_lob) && (r_rtn_max - r_rtn_[niteration]>=n_member_lob) ){
      index_sub_[n_sub] = niteration; nlp_ZR_sub_[n_sub] = nlp_ZR_[niteration]; n_sub++;
      /* if sufficiently many members */}
    /* for (niteration=0;niteration<n_iteration;niteration++){ } */}
  if (n_sub<1){
    printf(" %% Warning, n_sub %d in dexcluster_nonbinary_recursive_helper_ZR__, using all of nlp_ZR_\n",n_sub);
    n_sub = n_iteration;
    for (niteration=0;niteration<n_iteration;niteration++){
      index_sub_[niteration] = niteration;
      nlp_ZR_sub_[niteration] = nlp_ZR_[niteration];
      /* for (niteration=0;niteration<n_iteration;niteration++){ } */}
    /* if (n_sub<1){ } */}
  for (nshuffle=0;nshuffle<1+n_shuffle;nshuffle++){
    ZR_sub_max_s_[nshuffle] = ZR__[index_sub_[0] + nshuffle*n_iteration];
    for (nsub=0;nsub<n_sub;nsub++){
      ZR = ZR__[index_sub_[nsub] + nshuffle*n_iteration];
      ZR_sub_max_s_[nshuffle] = maximum(ZR,ZR_sub_max_s_[nshuffle]);
      /* for (nsub=0;nsub<n_sub;nsub++){ } */}
    /* for (nshuffle=0;nshuffle<1+n_shuffle;nshuffle++){ } */}
  gumbel_fit(n_shuffle,ZR_sub_max_s_+1,ZR_sub_max_s_[0],NULL,&nlp_gumb_opt,&nlp_gumb_emp,&p_gumb_opt,&p_gumb_emp);
  if (nlp_gumb_opt_p!=NULL){ *nlp_gumb_opt_p = nlp_gumb_opt;}
  if (nlp_gumb_emp_p!=NULL){ *nlp_gumb_emp_p = nlp_gumb_emp;}
  if (verbose>2){ printf(" %% nlp_gumb_opt %f p_gumb_opt %f nlp_gumb_emp %f p_gumb_emp %f\n",nlp_gumb_opt,p_gumb_opt,nlp_gumb_emp,p_gumb_emp);}
  find_internal_maximum(0*verbose,n_sub,nlp_ZR_sub_,0,&nlp_ZR_max,&nlp_ZR_index);
  if (nlp_ZR_max>=-log(p_use)){ find_internal_maximum(0*verbose,n_sub,nlp_ZR_sub_,-log(p_use),&nlp_ZR_max,&nlp_ZR_index);}
  if (nlp_ZR_index<=-1){ printf(" %% Warning, no maximum found\n"); nlp_ZR_index = 0; nlp_ZR_max = nlp_ZR_sub_[nlp_ZR_index]; }
  if (verbose>2){ printf(" %% index_sub_[%d] = %d\n",nlp_ZR_index,index_sub_[nlp_ZR_index]);}
  nlp_ZR_index = index_sub_[nlp_ZR_index];
  n_r_rtn_index = trace__[trace_r_rtn_index + nlp_ZR_index*6]; n_r_rmv_index = r_rtn_max - n_r_rtn_index;
  n_c_rtn_index = trace__[trace_c_rtn_index + nlp_ZR_index*6]; n_c_rmv_index = c_rtn_max - n_c_rtn_index;
  if (verbose){ printf(" %% nlp_ZR_max %f nlp_ZR_index %d n_r_rtn_index %d/%d n_c_rtn_index %d/%d\n",nlp_ZR_max,nlp_ZR_index,n_r_rtn_index,n_r_rmv_index,n_c_rtn_index,n_c_rmv_index);}
  if (nlp_ZR_max_p!=NULL){ *nlp_ZR_max_p = nlp_ZR_max;}
  if (nlp_ZR_index_p!=NULL){ *nlp_ZR_index_p = nlp_ZR_index;}
  if (n_r_rtn_index_p!=NULL){ *n_r_rtn_index_p = n_r_rtn_index;}
  if (n_r_rmv_index_p!=NULL){ *n_r_rmv_index_p = n_r_rmv_index;}
  if (n_c_rtn_index_p!=NULL){ *n_c_rtn_index_p = n_c_rtn_index;}
  if (n_c_rmv_index_p!=NULL){ *n_c_rmv_index_p = n_c_rmv_index;}
  r_rtn_index_=NULL;
  if (r_rtn_index_p_!=NULL){
    if ( (*r_rtn_index_p_)==NULL ){ (*r_rtn_index_p_) = (int *) malloc1((unsigned long long int)n_r_rtn_index*sizeof(int));}
    r_rtn_index_ = *r_rtn_index_p_;
    /* if (r_rtn_index_p_!=NULL){ } */}
  if (r_rtn_index_!=NULL){
    xdrop_ = &(xdrop__[0 + (n_xdrop-1)*2]);
    nl = 0; while (nl<n_r_rtn_index){ if ((*xdrop_)>-1){ r_rtn_index_[nl++] = *xdrop_;} xdrop_-=2;}
    if (verbose>2){ array_printf_margin(r_rtn_index_,"int",1,n_r_rtn_index," % r_rtn_index_: ");}
    /* if (r_rtn_index_!=NULL){ } */}
  r_rmv_index_=NULL;
  if (r_rmv_index_p_!=NULL){
    if ( (*r_rmv_index_p_)==NULL ){ (*r_rmv_index_p_) = (int *) malloc1((unsigned long long int)n_r_rmv_index*sizeof(int));}
    r_rmv_index_ = *r_rmv_index_p_;
    /* if (r_rmv_index_p_!=NULL){ } */}
  if (r_rmv_index_!=NULL){
    xdrop_ = &(xdrop__[0 + (0)*2]);
    nl = 0; while (nl<n_r_rmv_index){ if ((*xdrop_)>-1){ r_rmv_index_[nl++] = *xdrop_;} xdrop_+=2;}
    if (verbose>2){ array_printf_margin(r_rmv_index_,"int",1,n_r_rmv_index," % r_rmv_index_: ");}
    /* if (r_rtn_index_!=NULL){ } */}
  if (flag_rdrop_vs_rcdrop==0){
    c_rtn_index_=NULL;
    if (c_rtn_index_p_!=NULL){
      if ( (*c_rtn_index_p_)==NULL ){ (*c_rtn_index_p_) = (int *) malloc1((unsigned long long int)n_c_rtn_index*sizeof(int));}
      c_rtn_index_ = *c_rtn_index_p_;
      /* if (c_rtn_index_p_!=NULL){ } */}
    if (c_rtn_index_!=NULL){
      xdrop_ = &(xdrop__[1 + (n_xdrop-1)*2]);
      nl = 0; while (nl<n_c_rtn_index){ if ((*xdrop_)>-1){ c_rtn_index_[nl++] = *xdrop_;} xdrop_-=2;}
      if (verbose>2){ array_printf_margin(c_rtn_index_,"int",1,n_c_rtn_index," % c_rtn_index_: ");}
      /* if (c_rtn_index_!=NULL){ } */}
    c_rmv_index_=NULL;
    if (c_rmv_index_p_!=NULL){
      if ( (*c_rmv_index_p_)==NULL ){ (*c_rmv_index_p_) = (int *) malloc1((unsigned long long int)n_c_rmv_index*sizeof(int));}
      c_rmv_index_ = *c_rmv_index_p_;
      /* if (c_rmv_index_p_!=NULL){ } */}
    if (c_rmv_index_!=NULL){
      xdrop_ = &(xdrop__[1 + (0)*2]);
      nl = 0; while (nl<n_c_rmv_index){ if ((*xdrop_)>-1){ c_rmv_index_[nl++] = *xdrop_;} xdrop_+=2;}
      if (verbose>2){ array_printf_margin(c_rmv_index_,"int",1,n_c_rmv_index," % c_rmv_index_: ");}
      /* if (c_rmv_index_!=NULL){ } */}    
    /* if (flag_rdrop_vs_rcdrop==0){ } */}
  free1(&index_sub_);
  free1(&nlp_ZR_sub_);
  free1(&nlp_ZR_);
  free1(&ZR__);
  free1(&ZR_sub_max_s_);
  free1(&r_rtn_);
  free1(&c_rtn_);
  free1(&trace__);
  free1(&QR__);
  free1(&xdrop__);
  free1(&MDA_dim_);
  /* %%%%%%%% */
  if (verbose){ printf(" %% [finished dexcluster_nonbinary_recursive_helper_ZR__]\n");}
}

void dexcluster_nonbinary_f_recursive
(
 char *dir_trunk_0in
 ,char *dir_out_0in
 ,char *prefix_base_0in
 ,int n_r
 ,int n_c
 ,float *E_base_rc__
 ,int n_r_index_0in
 ,int *r_index_0in_
 ,int n_c_index_0in
 ,int *c_index_0in_
 ,double gamma_0in
 ,int n_shuffle_0in
 ,double p_set_0in
 ,int n_member_lob_0in
 ,double p_prev_0in
 ,int flag_force_create_0in
 ,char ***output_label_p_
 ,char ***nlpbra_label_p_
 ,char ***nlpnex_label_p_
)
{
  int verbose=1;
  int flag_rcdrop = 0;
  int flag_rdrop = 1;
  char *cwd=NULL;
  int n_str_path=0;
  char dir_trunk[PNAMESIZE];
  char dir_out[PNAMESIZE],dir_A_out[PNAMESIZE],dir_B_out[PNAMESIZE];
  char prefix_base0[FNAMESIZE];
  char prefix_base1[FNAMESIZE];
  char dir_0in[PNAMESIZE];
  char prefix_n_member_lob[FNAMESIZE];
  char prefix_gamma[FNAMESIZE];
  struct stat stat_dir = {0};
  struct stat stat_file = {0};
  int nr=0,flag_free_r_index=0,nr_index=0,n_r_index=0,*r_index_=NULL;
  int nc=0,flag_free_c_index=0,n_c_index=0,*c_index_=NULL;
  double gamma=0;
  int nshuffle=0,n_shuffle=0;
  int n_member_lob=0;
  double p_set=0,p_use=0,p_prev=0;
  int flag_force_create=0;
  int MDA_n_dim=0;
  int *MDA_dim_=NULL;
  char MDA_fname[PNAMESIZE];
  char **output_label__=NULL;
  char **nlpbra_label__=NULL;
  char **nlpnex_label__=NULL;
  char **output_label_A__=NULL;
  char **nlpbra_label_A__=NULL;
  char **nlpnex_label_A__=NULL;
  char **output_label_B__=NULL;
  char **nlpbra_label_B__=NULL;
  char **nlpnex_label_B__=NULL;
  char fname_trace_E__[PNAMESIZE];
  char fname_QR_E__[PNAMESIZE];
  char fname_xdrop_E__[PNAMESIZE];
  char fname_trace_F__[PNAMESIZE];
  char fname_QR_F__[PNAMESIZE];
  char fname_xdrop_F__[PNAMESIZE];
  double nlp_ZR_E_max=0; int nlp_ZR_E_index=0;
  double nlp_ZR_F_max=0; int nlp_ZR_F_index=0;
  int n_c_rtn_index_E=0,*c_rtn_index_E_=NULL,*c_index_rtn_sub_=NULL;
  double nlp_gumb_opt=0,nlp_gumb_emp=0;
  int nr_rtn_index_F=0,n_r_rtn_index_F=0,*r_rtn_index_F_=NULL,*r_index_rtn_sub_=NULL;
  int nr_rmv_index_F=0,n_r_rmv_index_F=0,*r_rmv_index_F_=NULL,*r_index_rmv_sub_=NULL;
  double *nlp_p_split_=NULL;
  double p_branch=0,p_next=0;
  int flag_split=0;
  if (verbose>0){ printf(" %% [entering dexcluster_nonbinary_f_recursive]\n");}
  /* %%%%%%%%%%%%%%%% */
  n_str_path = pathconf(".",_PC_PATH_MAX);
  MDA_dim_ = (int *) malloc1(2*sizeof(int)); //%<-- should only ever need 2. ;
  nlp_p_split_ = (double *) malloc1(9*sizeof(double)); //%<-- gumb_emp, gumb_opt, F, E, p_branch, p_next, n_r_rtn_index_F, n_r_rmv_index_F, flag_split. ;
  cwd = (char *) malloc1((size_t)n_str_path); getcwd(cwd,n_str_path);
  if (n_str_path> FNAMESIZE){ printf(" %% Warning, n_str_path %d in dexcluster_nonbinary_f_recursive\n",n_str_path);}
  if ( (dir_trunk_0in!=NULL) && (strlen(dir_trunk_0in)> FNAMESIZE) ){ printf(" %% Warning, dir_trunk_0in %s too long in dexcluster_nonbinary_f_recursive\n",dir_trunk_0in);}
  if ( (dir_out_0in!=NULL) && (strlen(dir_out_0in)> FNAMESIZE) ){ printf(" %% Warning, dir_out_0in %s too long in dexcluster_nonbinary_f_recursive\n",dir_out_0in);}
  if ( (prefix_base_0in!=NULL) && (strlen(prefix_base_0in)> FNAMESIZE) ){ printf(" %% Warning, prefix_base_0in %s too long in dexcluster_nonbinary_f_recursive\n",prefix_base_0in);}
  if ( dir_trunk_0in==NULL ){ if (verbose>0){ printf(" %% cwd: [%s]\n",cwd);} sprintf(dir_trunk,"%s",cwd);} else{ sprintf(dir_trunk,"%s",dir_trunk_0in);}
  if ( prefix_base_0in==NULL ){ sprintf(prefix_base0,"%s","test");} else{ sprintf(prefix_base0,"%s",prefix_base_0in);}
  if (verbose>1){
    printf(" %% dir_trunk: %s\n",dir_trunk);
    printf(" %% prefix_base0: %s\n",prefix_base0);
    /* if (verbose>1){ } */}
  if ( (n_r_index_0in<=0) || (r_index_0in_==NULL) ){
    n_r_index = n_r;
    r_index_ = (int *) malloc1(n_r*sizeof(int));
    for (nr=0;nr<n_r;nr++){ r_index_[nr] = nr;}
    flag_free_r_index = 1;
    /* if ( (n_r_index_0in<=0) || (r_index_0in_==NULL) ){ } */}
  else /* if ( (n_r_index_0in>0) && (r_index_0in_!=NULL) ) */{
    n_r_index = n_r_index_0in;
    r_index_ = r_index_0in_;
    flag_free_r_index = 0;
    /* exists */}
  if ( (n_c_index_0in<=0) || (c_index_0in_==NULL) ){
    n_c_index = n_c;
    c_index_ = (int *) malloc1(n_c*sizeof(int));
    for (nc=0;nc<n_c;nc++){ c_index_[nc] = nc;}
    flag_free_c_index = 1;
    /* if ( (n_r_index_0in<=0) || (r_index_0in_==NULL) ){ } */}
  else /* if ( (n_c_index_0in>0) && (c_index_0in_!=NULL) ) */{
    n_c_index = n_c_index_0in;
    c_index_ = c_index_0in_;
    flag_free_c_index = 0;
    /* exists */}
  if (dir_out_0in!=NULL){
    sprintf(dir_out,"%s",dir_out_0in);
    /* if (dir_out_0in!=NULL){ } */}
  if (dir_out_0in==NULL){
    sprintf(dir_0in,"%s/dir_%s",dir_trunk,prefix_base0);
    if (stat(dir_0in,&stat_dir)==-1){ printf(" %% %s not found, creating\n",dir_0in); mkdir(dir_0in,0755);} else{ printf(" %% %s found, not creating\n",dir_0in);}
    if ( gamma_0in<=0 ){ gamma = 0.0; } else{ gamma = minimum(1,gamma_0in);}
    if (gamma> 0){ sprintf(prefix_gamma,"_g%.3d",floor(100*gamma));} else{ sprintf(prefix_gamma,"");}
    if (n_member_lob> 2){ sprintf(prefix_n_member_lob,"_n%.2d",n_member_lob);} else{ sprintf(prefix_n_member_lob,"");}
    sprintf(prefix_base1,"%s%s%s",prefix_base0,prefix_gamma,prefix_n_member_lob);
    sprintf(dir_out,"%s/dir_%s",dir_0in,prefix_base1);
    if (verbose>1){
      printf(" %% dir_0in: %s\n",dir_0in);
      printf(" %% prefix_base1: %s\n",prefix_base1);
      printf(" %% dir_out: %s\n",dir_out);
      /* if (verbose>1){ } */}
    /* if (dir_out_0in==NULL){ } */}
  if (stat(dir_out,&stat_dir)==-1){ printf(" %% %s not found, creating\n",dir_out); mkdir(dir_out,0755);} else{ printf(" %% %s found, not creating\n",dir_out);}
  if (n_shuffle_0in<=0){ n_shuffle = 64;} else{ n_shuffle = n_shuffle_0in;}
  if (p_set_0in<=0){ p_set = 0.05;} else{ p_set = p_set_0in;}
  p_use = (double)p_set / (double)(1+2*p_set);
  if (n_member_lob_0in<=0){ n_member_lob = 2;} else{ n_member_lob = n_member_lob_0in;}
  if (p_prev_0in<=0){ p_prev = 0.00;} else{ p_prev = p_prev_0in;}
  if (flag_force_create_0in<=0){ flag_force_create = 0;} else{ flag_force_create = flag_force_create_0in;}
  if (verbose>0){
    printf(" %% n_r %d n_c %d --> n_r_index %d n_c_index %d\n",n_r,n_c,n_r_index,n_c_index);
    printf(" %% gamma %0.3f n_shuffle %d p_set %0.3f (p_use %0.3f) n_member_lob %d p_prev %0.3f flag_force_create %d\n",gamma,n_shuffle,p_set,p_use,n_member_lob,p_prev,flag_force_create);
    /* if (verbose>0){ } */}
  sprintf(MDA_fname,"%s/r_index_.mda",dir_out); MDA_n_dim = 2; MDA_dim_[0] = n_r_index; MDA_dim_[1] = 1;
  MDA_write_i4(MDA_n_dim,MDA_dim_,r_index_,MDA_fname);
  if (verbose>2){ MDA_printf_i4_margin(MDA_fname);}
  //sprintf(MDA_fname,"%s/c_index_.mda",dir_out); MDA_n_dim = 2; MDA_dim_[0] = n_c_index; MDA_dim_[1] = 1;
  //MDA_write_i4(MDA_n_dim,MDA_dim_,c_index_,MDA_fname);
  //if (verbose>2){ MDA_printf_i4_margin(MDA_fname);}
  /* %%%%%%%%%%%%%%%% */
  /* initialize labels */
  /* %%%%%%%%%%%%%%%% */
  output_label__ = NULL;
  if (output_label_p_!=NULL){
    if ( (*output_label_p_)==NULL ){
      (*output_label_p_) = (char **) malloc1((unsigned long long int)n_r_index*sizeof(char *)); 
      for (nr_index=0;nr_index<n_r_index;nr_index++){ (*output_label_p_)[nr_index] = NULL;}
      /* if ( (*output_label_p_)==NULL ){ } */}
    output_label__ = *output_label_p_;
    /* if (output_label_p_!=NULL){ } */}
  if (output_label__!=NULL){
    for (nr_index=0;nr_index<n_r_index;nr_index++){
      if (output_label__[nr_index]==NULL){ output_label__[nr_index] = (char *) malloc1((unsigned long long int)FNAMESIZE*sizeof(char));}
      sprintf(output_label__[nr_index],"0");
      /* for (nr_index=0;nr_index<n_r_index;nr_index++){ } */}
    /* if (output_label__!=NULL){ } */}
  /* %%%% */
  nlpbra_label__ = NULL;
  if (nlpbra_label_p_!=NULL){
    if ( (*nlpbra_label_p_)==NULL ){ (*nlpbra_label_p_) = (char **) malloc1((unsigned long long int)n_r_index*sizeof(char *));
      for (nr_index=0;nr_index<n_r_index;nr_index++){ (*nlpbra_label_p_)[nr_index] = NULL;}
      /* if ( (*nlpbra_label_p_)==NULL ){ } */}
    nlpbra_label__ = *nlpbra_label_p_;
    /* if (nlpbra_label_p_!=NULL){ } */}
  if (nlpbra_label__!=NULL){
    for (nr_index=0;nr_index<n_r_index;nr_index++){
      if (nlpbra_label__[nr_index]==NULL){ nlpbra_label__[nr_index] = (char *) malloc1((unsigned long long int)FNAMESIZE*sizeof(char));}
      sprintf(nlpbra_label__[nr_index],"");
      /* for (nr_index=0;nr_index<n_r_index;nr_index++){ } */}
    /* if (nlpbra_label__!=NULL){ } */}
  /* %%%% */
  nlpnex_label__ = NULL;
  if (nlpnex_label_p_!=NULL){
    if ( (*nlpnex_label_p_)==NULL ){ (*nlpnex_label_p_) = (char **) malloc1((unsigned long long int)n_r_index*sizeof(char *));
      for (nr_index=0;nr_index<n_r_index;nr_index++){ (*nlpnex_label_p_)[nr_index] = NULL;}
      /* if ( (*nlpnex_label_p_)==NULL ){ } */}
    nlpnex_label__ = *nlpnex_label_p_;
    /* if (nlpnex_label_p_!=NULL){ } */}
  if (nlpnex_label__!=NULL){
    for (nr_index=0;nr_index<n_r_index;nr_index++){
      if (nlpnex_label__[nr_index]==NULL){ nlpnex_label__[nr_index] = (char *) malloc1((unsigned long long int)FNAMESIZE*sizeof(char));}
      sprintf(nlpnex_label__[nr_index],"");
      /* for (nr_index=0;nr_index<n_r_index;nr_index++){ } */}
    /* if (nlpnex_label__!=NULL){ } */}
  /* %%%%%%%%%%%%%%%% */
  /* perform calculation */
  /* %%%%%%%%%%%%%%%% */
  sprintf(fname_trace_E__,"%s/trace_E__.mda",dir_out);
  sprintf(fname_QR_E__,"%s/QR_E__.mda",dir_out);
  sprintf(fname_xdrop_E__,"%s/xdrop_E__.mda",dir_out);
  dexcluster_nonbinary_f_recursive_helper_QR__
    (
     verbose-1
     ,flag_rcdrop
     ,flag_force_create
     ,n_r
     ,n_c
     ,E_base_rc__
     ,n_r_index
     ,r_index_
     ,n_c_index
     ,c_index_
     ,gamma
     ,n_shuffle
     ,fname_trace_E__
     ,fname_xdrop_E__
     ,fname_QR_E__
     );
  dexcluster_nonbinary_recursive_helper_ZR__
    (
     verbose-1
     ,flag_rcdrop
     ,fname_trace_E__
     ,fname_xdrop_E__
     ,fname_QR_E__
     ,p_use
     ,0*n_member_lob
     ,&nlp_ZR_E_max
     ,&nlp_ZR_E_index
     ,NULL
     ,NULL
     ,NULL
     ,NULL
     ,&n_c_rtn_index_E
     ,&c_rtn_index_E_
     ,NULL
     ,NULL
     ,&nlp_gumb_opt
     ,&nlp_gumb_emp
     );
  array_extract_i_from_i(1,n_c_index,c_index_,0,NULL,n_c_rtn_index_E,c_rtn_index_E_,&c_index_rtn_sub_,NULL);
  iquicksort_index(0,c_index_rtn_sub_,1,NULL,0,n_c_rtn_index_E-1);
  if (verbose>2){ array_printf_margin(c_index_rtn_sub_,"int",1,n_c_rtn_index_E," % c_index_rtn_sub_: ");}
  /* %%%% */
  sprintf(fname_trace_F__,"%s/trace_F__.mda",dir_out);
  sprintf(fname_QR_F__,"%s/QR_F__.mda",dir_out);
  sprintf(fname_xdrop_F__,"%s/xdrop_F__.mda",dir_out);
  dexcluster_nonbinary_f_recursive_helper_QR__
    (
     verbose-1
     ,flag_rdrop
     ,flag_force_create
     ,n_r
     ,n_c
     ,E_base_rc__
     ,n_r_index
     ,r_index_
     ,n_c_rtn_index_E
     ,c_index_rtn_sub_
     ,0*gamma
     ,n_shuffle
     ,fname_trace_F__
     ,fname_xdrop_F__
     ,fname_QR_F__
     );
  dexcluster_nonbinary_recursive_helper_ZR__
    (
     verbose-1
     ,flag_rdrop
     ,fname_trace_F__
     ,fname_xdrop_F__
     ,fname_QR_F__
     ,p_use
     ,1*n_member_lob
     ,&nlp_ZR_F_max
     ,&nlp_ZR_F_index
     ,&n_r_rtn_index_F
     ,&r_rtn_index_F_
     ,&n_r_rmv_index_F
     ,&r_rmv_index_F_
     ,NULL
     ,NULL
     ,NULL
     ,NULL
     ,NULL
     ,NULL
     );
  array_extract_i_from_i(1,n_r_index,r_index_,0,NULL,n_r_rtn_index_F,r_rtn_index_F_,&r_index_rtn_sub_,NULL);
  iquicksort_index(0,r_index_rtn_sub_,1,NULL,0,n_r_rtn_index_F-1);
  if (verbose>2){ array_printf_margin(r_index_rtn_sub_,"int",1,n_r_rtn_index_F," % r_index_rtn_sub_: ");}
  array_extract_i_from_i(1,n_r_index,r_index_,0,NULL,n_r_rmv_index_F,r_rmv_index_F_,&r_index_rmv_sub_,NULL);
  iquicksort_index(0,r_index_rmv_sub_,1,NULL,0,n_r_rmv_index_F-1);
  if (verbose>2){ array_printf_margin(r_index_rmv_sub_,"int",1,n_r_rmv_index_F," % r_index_rmv_sub_: ");}
  /* %%%% */
  p_branch = exp(-nlp_gumb_emp); if (p_branch<(double)1.0/(double)n_shuffle){ p_branch = exp(-nlp_gumb_opt);}
  p_next = 1.0 - (1.0-p_prev)*(1-p_branch); //%<-- q_next = q_prev*q_branch;
  flag_split = (n_r_rtn_index_F>=n_member_lob) && (n_r_rmv_index_F>=n_member_lob) && (p_next<=p_use);
  if (verbose>0){ printf(" %% p_branch %f p_next %f n_r_rtn_index_F %d n_r_rmv_index_F %d flag_split %d\n",p_branch,p_next,n_r_rtn_index_F,n_r_rmv_index_F,flag_split);} 
  nlp_p_split_[0] = nlp_gumb_emp;
  nlp_p_split_[1] = nlp_gumb_opt;
  nlp_p_split_[2] = nlp_ZR_E_max;
  nlp_p_split_[3] = nlp_ZR_F_max;
  nlp_p_split_[4] = p_branch;
  nlp_p_split_[5] = p_next;
  nlp_p_split_[6] = (double)n_r_rtn_index_F;
  nlp_p_split_[7] = (double)n_r_rmv_index_F;
  nlp_p_split_[8] = (double)flag_split;
  sprintf(MDA_fname,"%s/nlp_p_split_.mda",dir_out); MDA_n_dim = 2; MDA_dim_[0] = 9; MDA_dim_[1] = 1;
  MDA_write_r8(MDA_n_dim,MDA_dim_,nlp_p_split_,MDA_fname);
  if (verbose>2){ MDA_printf_r8_margin(MDA_fname);}
  if (flag_split){
    /* %%%% */
    sprintf(dir_A_out,"%s/A",dir_out);
    if (stat(dir_A_out,&stat_dir)==-1){ printf(" %% %s not found, creating\n",dir_A_out); mkdir(dir_A_out,0755);} else{ printf(" %% %s found, not creating\n",dir_A_out);}
    dexcluster_nonbinary_f_recursive
      (
       dir_trunk
       ,dir_A_out
       ,prefix_base0
       ,n_r
       ,n_c
       ,E_base_rc__
       ,n_r_rtn_index_F
       ,r_index_rtn_sub_
       ,n_c_index
       ,c_index_
       ,gamma
       ,n_shuffle
       ,p_set
       ,n_member_lob
       ,p_next
       ,flag_force_create
       ,&output_label_A__
       ,&nlpbra_label_A__
       ,&nlpnex_label_A__
       );      
    /* %%%% */
    sprintf(dir_B_out,"%s/B",dir_out);
    if (stat(dir_B_out,&stat_dir)==-1){ printf(" %% %s not found, creating\n",dir_B_out); mkdir(dir_B_out,0755);} else{ printf(" %% %s found, not creating\n",dir_B_out);}
    dexcluster_nonbinary_f_recursive
      (
       dir_trunk
       ,dir_B_out
       ,prefix_base0
       ,n_r
       ,n_c
       ,E_base_rc__
       ,n_r_rmv_index_F
       ,r_index_rmv_sub_
       ,n_c_index
       ,c_index_
       ,gamma
       ,n_shuffle
       ,p_set
       ,n_member_lob
       ,p_next
       ,flag_force_create
       ,&output_label_B__
       ,&nlpbra_label_B__
       ,&nlpnex_label_B__
       );      
    /* %%%% */
    /* update labels */
    /* %%%% */
    for (nr_rtn_index_F=0;nr_rtn_index_F<n_r_rtn_index_F;nr_rtn_index_F++){
      nr_index = r_rtn_index_F_[nr_rtn_index_F];
      if (output_label__!=NULL){ if (strlen(output_label_A__[nr_rtn_index_F])>FNAMESIZE-16){ printf(" %% Warning, too many clusters, increase label length\n");} sprintf(output_label__[nr_index],"%s%s","A",output_label_A__[nr_rtn_index_F]);}
      if (nlpbra_label__!=NULL){ if (strlen(nlpbra_label_A__[nr_rtn_index_F])>FNAMESIZE-16){ printf(" %% Warning, too many clusters, increase label length\n");} sprintf(nlpbra_label__[nr_index],"%0.2f %s",-log(p_branch),nlpbra_label_A__[nr_rtn_index_F]);}
      if (nlpnex_label__!=NULL){ if (strlen(nlpnex_label_A__[nr_rtn_index_F])>FNAMESIZE-16){ printf(" %% Warning, too many clusters, increase label length\n");} sprintf(nlpnex_label__[nr_index],"%0.2f %s",-log(p_next),nlpnex_label_A__[nr_rtn_index_F]);}
      /* for (nr_rtn_index_F=0;nr_rtn_index_F<n_r_rtn_index_F;nr_rtn_index_F++){ } */}
    for (nr_rmv_index_F=0;nr_rmv_index_F<n_r_rmv_index_F;nr_rmv_index_F++){
      nr_index = r_rmv_index_F_[nr_rmv_index_F];
      if (output_label__!=NULL){ if (strlen(output_label_B__[nr_rmv_index_F])>FNAMESIZE-16){ printf(" %% Warning, too many clusters, increase label length\n");} sprintf(output_label__[nr_index],"%s%s","B",output_label_B__[nr_rmv_index_F]);}
      if (nlpbra_label__!=NULL){ if (strlen(nlpbra_label_B__[nr_rmv_index_F])>FNAMESIZE-16){ printf(" %% Warning, too many clusters, increase label length\n");} sprintf(nlpbra_label__[nr_index],"%0.2f %s",-log(p_branch),nlpbra_label_B__[nr_rmv_index_F]);}
      if (nlpnex_label__!=NULL){ if (strlen(nlpnex_label_B__[nr_rmv_index_F])>FNAMESIZE-16){ printf(" %% Warning, too many clusters, increase label length\n");} sprintf(nlpnex_label__[nr_index],"%0.2f %s",-log(p_next),nlpnex_label_B__[nr_rmv_index_F]);}
      /* for (nr_rmv_index_F=0;nr_rmv_index_F<n_r_rmv_index_F;nr_rmv_index_F++){ } */}
    /* %%%% */
    /* free labels */
    /* %%%% */
    for (nr_rtn_index_F=0;nr_rtn_index_F<n_r_rtn_index_F;nr_rtn_index_F++){
      free1(&output_label_A__[nr_rtn_index_F]); free1(&nlpbra_label_A__[nr_rtn_index_F]); free1(&nlpnex_label_A__[nr_rtn_index_F]);
      /* for (nr_rtn_index_F=0;nr_rtn_index_F<n_r_rtn_index_F;nr_rtn_index_F++){ } */}
    free1(&output_label_A__); free1(&nlpbra_label_A__); free1(&nlpnex_label_A__);
    for (nr_rmv_index_F=0;nr_rmv_index_F<n_r_rmv_index_F;nr_rmv_index_F++){
      free1(&output_label_B__[nr_rmv_index_F]); free1(&nlpbra_label_B__[nr_rmv_index_F]); free1(&nlpnex_label_B__[nr_rmv_index_F]);
      /* for (nr_rmv_index_F=0;nr_rmv_index_F<n_r_rmv_index_F;nr_rmv_index_F++){ } */}
    free1(&output_label_B__); free1(&nlpbra_label_B__); free1(&nlpnex_label_B__);
    /* %%%% */
    /* if (flag_split){ } */}
  /* %%%% */
  free1(&r_index_rmv_sub_);
  free1(&r_index_rtn_sub_);
  free1(&c_index_rtn_sub_);
  free1(&c_rtn_index_E_);
  free1(&r_rmv_index_F_);
  free1(&r_rtn_index_F_);
  /* %%%%%%%%%%%%%%%% */
  if (flag_free_r_index){ free1(&r_index_);}
  if (flag_free_c_index){ free1(&c_index_);}
  free1(&nlp_p_split_);
  free1(&cwd);
  free1(&MDA_dim_);
  /* %%%%%%%%%%%%%%%% */
  if (verbose>0){ printf(" %% [finished dexcluster_nonbinary_f_recursive]\n");}
}

void dexcluster_nonbinary_f_recursive_test()
{
  int flag_force_create=1;
  int nr=0,n_r = 128;
  int nc=0,n_c = 512;
  float x=0,y=0,z=0;
  unsigned long int rseed=0;
  unsigned long long int ulli=0;
  float *E_base_rc__=NULL;
  char **output_label__=NULL;
  char **nlpbra_label__=NULL;
  char **nlpnex_label__=NULL;
  E_base_rc__ = (float *) malloc1((unsigned long long int)n_r*(unsigned long long int)n_c*sizeof(float));
  rseed=1; RSEED_adv8(&rseed);
  ulli=0;
  for (nc=0;nc<n_c;nc++){ for (nr=0;nr<n_r;nr++){ 
      //E_base_rc__[ulli] = RNGET(&rseed);
      //if ( (nc<floor(sqrt((double)n_c))) && (nr<floor(sqrt((double)n_r))) ){ E_base_rc__[ulli] += 0.5; }
      x = 2.0*(float)nr/(float)(n_r-1) - 1.0;
      y = 2.0*(float)nc/(float)(n_c-1) - 1.0;
      z = (x+y)/4.0;
      E_base_rc__[ulli] = sin(2*PI*x) + cos(2*PI*2*y) + x*x + y*y*y + cos(2*PI*4*z);
      ulli++;
      /* for (nc=0;nc<n_c;nc++){ for (nr=0;nr<n_r;nr++){ }} */}}
  array_printf_margin(E_base_rc__,"float",n_r,n_c," % E_base_r__: ");
  GLOBAL_tic(0);
  dexcluster_nonbinary_f_recursive(
  NULL
 ,NULL
 ,NULL
 ,n_r
 ,n_c
 ,E_base_rc__
 ,0
 ,NULL
 ,0
 ,NULL
 ,0
 ,0
 ,0
 ,0
 ,0
 ,flag_force_create
 ,&output_label__
 ,&nlpbra_label__
 ,&nlpnex_label__
 );
  for (nr=0;nr<n_r;nr++){
    printf(" %% nr %.3d: %s --> %s --> %s \n",nr,output_label__[nr],nlpbra_label__[nr],nlpnex_label__[nr]);
    /* for (nr=0;nr<n_r;nr++){ } */}
  for (nr=0;nr<n_r;nr++){
    free1(&output_label__[nr]); free1(&nlpbra_label__[nr]); free1(&nlpnex_label__[nr]);
    /* for (nr=0;nr<n_r;nr++){ } */}
  free1(&output_label__); free1(&nlpbra_label__); free1(&nlpnex_label__);
  GLOBAL_toc(0,1," % dexcluster_nonbinary_recursive: ");
  free1(&E_base_rc__);
}
