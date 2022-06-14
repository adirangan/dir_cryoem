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
  int trace_QR_index = 3;
  int n_iteration=0,tmp_n_iteration=0;
  int nshuffle=0;
  unsigned long int rseed=0;
  double *trace__=NULL;
  double *trace_shuffle__=NULL;
  double *QR_=NULL;
  double *QR__=NULL;
  double *xdrop__=NULL;
  double *xdrop_shuffle__=NULL;
  float *E_cr__=NULL;
  float *E_rc__=NULL;
  float *Q_rc__=NULL;
  float *QE_rc__=NULL;
  int MDA_n_dim=0;
  int *MDA_dim_=NULL;
  char MDA_fname[PNAMESIZE];
  if (verbose){ printf(" %% [entering dexcluster_nonbinary_f_recursive_helper_QR__]\n");}
  if ( flag_force_create || (stat(fname_trace__,&stat_file)==-1) || (stat(fname_QR__,&stat_file)==-1) || (stat(fname_xdrop__,&stat_file)==-1) ){
    if (verbose){ printf(" %% %s not found, creating\n",fname_trace__);}
    MDA_dim_ = (int *) malloc1(2*sizeof(int));
    n_iteration = get_xdrop_logscale_length(n_r_index,n_c_index,gamma);
    trace__ = (double *) malloc1((unsigned long long int)6*(unsigned long long int)n_r*sizeof(double));
    trace_shuffle__ = (double *) malloc1((unsigned long long int)6*(unsigned long long int)n_r*sizeof(double));
    QR__ = (double *) malloc1((unsigned long long int)n_iteration*(unsigned long long int)(1+n_shuffle)*sizeof(double));
    xdrop__ = (int *) malloc1((unsigned long long int)2*(unsigned long long int)(n_r_index+n_c_index)*sizeof(double));
    xdrop_shuffle__ = (int *) malloc1((unsigned long long int)2*(unsigned long long int)(n_r_index+n_c_index)*sizeof(double));
    E_rc__ = (float *) malloc1((unsigned long long int)n_r_index*(unsigned long long int)n_c_index*sizeof(float));
    E_cr__ = (float *) malloc1((unsigned long long int)n_c_index*(unsigned long long int)n_r_index*sizeof(float));
    array_extract_f_from_f(n_r,n_c,E_base_rc__,n_r_index,r_index_,n_c_index,c_index_,&E_rc__,NULL);
    if (verbose>9){ array_printf_margin(E_rc__,"float",n_r_index,n_c_index," % E_rc__: "); printf(" %% %% %% %%\n");}
    array_mean_center_row(n_r_index,n_c_index,E_rc__,NULL,"float",&E_rc__,&E_cr__);
    if (verbose>9){ array_printf_margin(E_rc__,"float",n_r_index,n_c_index," % E_rc__: "); printf(" %% %% %% %%\n");}
    if (verbose>9){ array_printf_margin(E_cr__,"float",n_c_index,n_r_index," % E_cr__: "); printf(" %% %% %% %%\n");}
    if (flag_rdrop_vs_rcdrop==0){ dexcluster_nonbinary_f(n_r_index,n_c_index,E_rc__,E_cr__,gamma,&xdrop__,&tmp_n_iteration,&trace__);}
    if (flag_rdrop_vs_rcdrop==1){ dexcluster_nonbinary_rdrop_f(n_r_index,n_c_index,E_rc__,E_cr__,gamma,&xdrop__,&tmp_n_iteration,&trace__);}
    if (verbose>9){ array_printf_margin(xdrop__,"int",2,n_r_index+n_c_index," % xdrop__: "); printf(" %% %% %% %%\n");}
    if (verbose>9){ array_printf_margin(trace__,"double",6,n_iteration," % trace__: "); printf(" %% %% %% %%\n");}
    if (tmp_n_iteration!=n_iteration){ printf(" %% Warning, tmp_n_iteration %d not equal n_iteration %d\n",tmp_n_iteration,n_iteration);}
    sprintf(MDA_fname,"%s",fname_trace__); MDA_n_dim = 2; MDA_dim_[0] = 6; MDA_dim_[1] = n_iteration;
    MDA_write_r8(MDA_n_dim,MDA_dim_,trace__,MDA_fname);
    if (verbose>9){ MDA_printf_r8_margin(MDA_fname);}
    sprintf(MDA_fname,"%s",fname_xdrop__); MDA_n_dim = 2; MDA_dim_[0] = 2; MDA_dim_[1] = n_r_index + n_c_index;
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
      if (verbose>9){ array_printf_margin(xdrop_shuffle__,"int",2,n_r_index+n_c_index," % xdrop_shuffle__: "); printf(" %% %% %% %%\n");}
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

void dexcluster_nonbinary_f_recursive_helper_QR_E__
(
  int verbose
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
 ,char *fname_trace_E__
 ,char *fname_xdrop_E__
 ,char *fname_QR_E__
)
{
  struct stat stat_file = {0};
  int trace_QR_index = 3;
  int n_iteration_E=0,tmp_n_iteration_E=0;
  int nshuffle=0;
  unsigned long int rseed=0;
  double *trace_E__=NULL;
  double *trace_E_shuffle__=NULL;
  double *QR_E_=NULL;
  double *QR_E__=NULL;
  double *xdrop_E__=NULL;
  double *xdrop_E_shuffle__=NULL;
  float *E_cr__=NULL;
  float *E_rc__=NULL;
  float *Q_rc__=NULL;
  float *QE_rc__=NULL;
  int MDA_n_dim=0;
  int *MDA_dim_=NULL;
  char MDA_fname[PNAMESIZE];
  if (verbose){ printf(" %% [entering dexcluster_nonbinary_f_recursive_helper_QR_E__]\n");}
  if ( flag_force_create || (stat(fname_trace_E__,&stat_file)==-1) || (stat(fname_QR_E__,&stat_file)==-1) || (stat(fname_xdrop_E__,&stat_file)==-1) ){
    if (verbose){ printf(" %% %s not found, creating\n",fname_trace_E__);}
    MDA_dim_ = (int *) malloc1(2*sizeof(int));
    n_iteration_E = get_xdrop_logscale_length(n_r_index,n_c_index,gamma);
    trace_E__ = (double *) malloc1((unsigned long long int)6*(unsigned long long int)n_r*sizeof(double));
    trace_E_shuffle__ = (double *) malloc1((unsigned long long int)6*(unsigned long long int)n_r*sizeof(double));
    QR_E__ = (double *) malloc1((unsigned long long int)n_iteration_E*(unsigned long long int)(1+n_shuffle)*sizeof(double));
    xdrop_E__ = (int *) malloc1((unsigned long long int)2*(unsigned long long int)(n_r_index+n_c_index)*sizeof(double));
    xdrop_E_shuffle__ = (int *) malloc1((unsigned long long int)2*(unsigned long long int)(n_r_index+n_c_index)*sizeof(double));
    E_rc__ = (float *) malloc1((unsigned long long int)n_r_index*(unsigned long long int)n_c_index*sizeof(float));
    E_cr__ = (float *) malloc1((unsigned long long int)n_c_index*(unsigned long long int)n_r_index*sizeof(float));
    array_extract_f_from_f(n_r,n_c,E_base_rc__,n_r_index,r_index_,n_c_index,c_index_,&E_rc__,NULL);
    if (verbose>9){ array_printf_margin(E_rc__,"float",n_r_index,n_c_index," % E_rc__: "); printf(" %% %% %% %%\n");}
    array_mean_center_row(n_r_index,n_c_index,E_rc__,NULL,"float",&E_rc__,&E_cr__);
    if (verbose>9){ array_printf_margin(E_rc__,"float",n_r_index,n_c_index," % E_rc__: "); printf(" %% %% %% %%\n");}
    if (verbose>9){ array_printf_margin(E_cr__,"float",n_c_index,n_r_index," % E_cr__: "); printf(" %% %% %% %%\n");}
    dexcluster_nonbinary_f(n_r_index,n_c_index,E_rc__,E_cr__,gamma,&xdrop_E__,&tmp_n_iteration_E,&trace_E__);
    if (verbose>9){ array_printf_margin(xdrop_E__,"int",2,n_r_index+n_c_index," % xdrop_E__: "); printf(" %% %% %% %%\n");}
    if (verbose>9){ array_printf_margin(trace_E__,"double",6,n_iteration_E," % trace_E__: "); printf(" %% %% %% %%\n");}
    if (tmp_n_iteration_E!=n_iteration_E){ printf(" %% Warning, tmp_n_iteration_E %d not equal n_iteration_E %d\n",tmp_n_iteration_E,n_iteration_E);}
    sprintf(MDA_fname,"%s",fname_trace_E__); MDA_n_dim = 2; MDA_dim_[0] = 6; MDA_dim_[1] = n_iteration_E;
    MDA_write_r8(MDA_n_dim,MDA_dim_,trace_E__,MDA_fname);
    if (verbose>9){ MDA_printf_r8_margin(MDA_fname);}
    sprintf(MDA_fname,"%s",fname_xdrop_E__); MDA_n_dim = 2; MDA_dim_[0] = 2; MDA_dim_[1] = n_r_index + n_c_index;
    MDA_write_i4(MDA_n_dim,MDA_dim_,xdrop_E__,MDA_fname);
    if (verbose>9){ MDA_printf_i4_margin(MDA_fname);}
    QR_E_ = QR_E__ + (unsigned long long int)0*(unsigned long long int)n_iteration_E;
    array_extract_d_from_d(6,n_iteration_E,trace_E__,1,&trace_QR_index,0,NULL,&QR_E_,NULL);
    if (verbose>9){ array_printf_margin(QR_E_,"double",1,n_iteration_E," % QR_E_: ");}
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
      dexcluster_nonbinary_f(n_r_index,n_c_index,E_rc__,E_cr__,gamma,&xdrop_E_shuffle__,&tmp_n_iteration_E,&trace_E_shuffle__);
      if (verbose>9){ array_printf_margin(xdrop_E_shuffle__,"int",2,n_r_index+n_c_index," % xdrop_E_shuffle__: "); printf(" %% %% %% %%\n");}
      if (verbose>9){ array_printf_margin(trace_E_shuffle__,"double",6,n_iteration_E," % trace_E_shuffle__: "); printf(" %% %% %% %%\n");}
      if (tmp_n_iteration_E!=n_iteration_E){ printf(" %% Warning, tmp_n_iteration_E %d not equal n_iteration_E %d\n",tmp_n_iteration_E,n_iteration_E);}
      QR_E_ = QR_E__ + (unsigned long long int)(1+nshuffle)*(unsigned long long int)n_iteration_E;
      array_extract_d_from_d(6,n_iteration_E,trace_E_shuffle__,1,&trace_QR_index,0,NULL,&QR_E_,NULL);
      if (verbose>9){ array_printf_margin(QR_E_,"double",1,n_iteration_E," % QR_E_: ");}
      /* for (nshuffle=0;nshuffle<n_shuffle;nshuffle++){ } */}
    if (verbose>9){ array_printf_margin(QR_E__,"double",n_iteration_E,(1+n_shuffle)," % QR_E__: ");}
    sprintf(MDA_fname,"%s",fname_QR_E__); MDA_n_dim = 2; MDA_dim_[0] = n_iteration_E; MDA_dim_[1] = (1+n_shuffle);
    MDA_write_r8(MDA_n_dim,MDA_dim_,QR_E__,MDA_fname);
    if (verbose>9){ MDA_printf_r8_margin(MDA_fname);}
    free1(&Q_rc__);
    free1(&QE_rc__);
    free1(&E_cr__);
    free1(&E_rc__);
    free1(&trace_E__);
    free1(&trace_E_shuffle__);
    free1(&QR_E__);
    free1(&xdrop_E__);
    free1(&xdrop_E_shuffle__);
    free1(&MDA_dim_);
    /* not found */}
  else{ if (verbose){ printf(" %% %s found, not creating\n",fname_trace_E__);}}
  if (verbose){ printf(" %% [finished dexcluster_nonbinary_f_recursive_helper_QR_E__]\n");}
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
 ,int *n_c_rtn_index_p
 ,int **c_rtn_index_p_
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
  int n_r_rtn_index=0,n_c_rtn_index=0;
  int *r_rtn_index_=NULL;
  int *c_rtn_index_=NULL;
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
  find_internal_maximum(0*verbose,n_sub,nlp_ZR_sub_,0,&nlp_ZR_max,&nlp_ZR_index);
  if (nlp_ZR_max>=-log(p_use)){ find_internal_maximum(0*verbose,n_sub,nlp_ZR_sub_,-log(p_use),&nlp_ZR_max,&nlp_ZR_index);}
  if (nlp_ZR_index<=-1){ printf(" %% Warning, no maximum found\n"); nlp_ZR_index = 0; nlp_ZR_max = nlp_ZR_sub_[nlp_ZR_index]; }
  if (verbose>2){ printf(" %% index_sub_[%d] = %d\n",nlp_ZR_index,index_sub_[nlp_ZR_index]);}
  nlp_ZR_index = index_sub_[nlp_ZR_index];
  n_r_rtn_index = trace__[trace_r_rtn_index + nlp_ZR_index*6];
  n_c_rtn_index = trace__[trace_c_rtn_index + nlp_ZR_index*6];
  if (verbose){ printf(" %% nlp_ZR_max %f nlp_ZR_index %d n_r_rtn_index %d n_c_rtn_index %d\n",nlp_ZR_max,nlp_ZR_index,n_r_rtn_index,n_c_rtn_index);}
  if (nlp_ZR_max_p!=NULL){ *nlp_ZR_max_p = nlp_ZR_max;}
  if (nlp_ZR_index_p!=NULL){ *nlp_ZR_index_p = nlp_ZR_index;}
  if (n_r_rtn_index_p!=NULL){ *n_r_rtn_index_p = n_r_rtn_index;}
  if (n_c_rtn_index_p!=NULL){ *n_c_rtn_index_p = n_c_rtn_index;}
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
      /* if (r_rtn_index_!=NULL){ } */}
    /* if (flag_rdrop_vs_rcdrop==0){ } */}
  free1(&index_sub_);
  free1(&nlp_ZR_sub_);
  free1(&nlp_ZR_);
  free1(&r_rtn_);
  free1(&c_rtn_);
  free1(&trace__);
  free1(&QR__);
  free1(&xdrop__);
  free1(&MDA_dim_);
  /* %%%%%%%% */
  if (verbose){ printf(" %% [finished dexcluster_nonbinary_recursive_helper_ZR__]\n");}
}

void dexcluster_nonbinary_f_recursive(
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
 ,char ***lpFmax_label_p_
 ,char ***lpnext_label_p_
)
{
  int verbose=3;
  char *cwd=NULL;
  int n_str_path=0;
  char dir_trunk[PNAMESIZE];
  char dir_out[PNAMESIZE];
  char prefix_base0[FNAMESIZE];
  char prefix_base1[FNAMESIZE];
  char dir_0in[PNAMESIZE];
  char prefix_n_member_lob[FNAMESIZE];
  char prefix_gamma[FNAMESIZE];
  struct stat stat_dir = {0};
  struct stat stat_file = {0};
  int nr=0,flag_free_r_index=0,n_r_index=0,*r_index_=NULL;
  int nc=0,flag_free_c_index=0,n_c_index=0,*c_index_=NULL;
  double gamma=0;
  int nshuffle=0,n_shuffle=0;
  int n_member_lob=0;
  double p_set=0,p_use=0,p_prev=0;
  int flag_force_create=0;
  int MDA_n_dim=0;
  int *MDA_dim_=NULL;
  char MDA_fname[PNAMESIZE];
  int niteration_E=0,n_iteration_E=0,tmp_n_iteration_E=0;
  int trace_QR_index = 3;
  unsigned long int rseed=0;
  float *Q_rc__=NULL,*QE_rc__=NULL;
  char fname_trace_E__[PNAMESIZE];
  char fname_QR_E__[PNAMESIZE];
  char fname_xdrop_E__[PNAMESIZE];
  float *E_rc__=NULL;
  float *E_cr__=NULL;
  double *trace_E__=NULL;
  double *trace_E_shuffle__=NULL;
  double *QR_E_=NULL,*QR_E__=NULL;
  double *ZR_E__=NULL;
  int *xdrop_E__=NULL;
  int *xdrop_E_shuffle__=NULL;
  char fname_trace_F__[PNAMESIZE];
  char fname_QR_F__[PNAMESIZE];
  char fname_xdrop_F__[PNAMESIZE];
  double *trace_F__=NULL;
  double *trace_F_shuffle_=NULL;
  double *QR_F__=NULL;
  double *ZR_F__=NULL;
  int *xdrop_F__=NULL;
  int *xdrop_F_shuffle__=NULL;
  double nlp_ZR_E_max=0;
  int nlp_ZR_E_index=0,n_r_rtn_index_E=0,n_c_rtn_index_E=0;
  int *r_rtn_index_E_=NULL,*c_rtn_index_E_=NULL;
  int *c_index_sub_=NULL;
  if (verbose){ printf(" %% [entering dexcluster_nonbinary_f_recursive]\n");}
  /* %%%%%%%%%%%%%%%% */
  n_str_path = pathconf(".",_PC_PATH_MAX);
  MDA_dim_ = (int *) malloc1(2*sizeof(int)); //%<-- should only ever need 2. ;
  cwd = (char *) malloc1((size_t)n_str_path);
  getcwd(cwd,n_str_path);
  if (verbose){ printf(" %% cwd: [%s]\n",cwd);}
  if (n_str_path> FNAMESIZE){ printf(" %% Warning, n_str_path %d in dexcluster_nonbinary_f_recursive\n",n_str_path);}
  if ( (dir_trunk_0in!=NULL) && (strlen(dir_trunk_0in)> FNAMESIZE) ){ printf(" %% Warning, dir_trunk_0in %s too long in dexcluster_nonbinary_f_recursive\n",dir_trunk_0in);}
  if ( (dir_out_0in!=NULL) && (strlen(dir_out_0in)> FNAMESIZE) ){ printf(" %% Warning, dir_out_0in %s too long in dexcluster_nonbinary_f_recursive\n",dir_out_0in);}
  if ( (prefix_base_0in!=NULL) && (strlen(prefix_base_0in)> FNAMESIZE) ){ printf(" %% Warning, prefix_base_0in %s too long in dexcluster_nonbinary_f_recursive\n",prefix_base_0in);}
  if ( dir_trunk_0in==NULL ){ sprintf(dir_trunk,"%s",cwd);} else{ sprintf(dir_trunk,"%s",dir_trunk_0in);}
  if ( prefix_base_0in==NULL ){ sprintf(prefix_base0,"%s","test");} else{ sprintf(prefix_base0,"%s",prefix_base_0in);}
  sprintf(dir_0in,"%s/dir_%s",dir_trunk,prefix_base0);
  if (stat(dir_0in,&stat_dir)==-1){ printf(" %% %s not found, creating\n",dir_0in); mkdir(dir_0in,0755);} else{ printf(" %% %s found, not creating\n",dir_0in);}
  if ( (n_r_index_0in<=0) || (r_index_0in_==NULL) ){
    n_r_index = n_r;
    r_index_ = (int *) malloc1(n_r*sizeof(int));
    for (nr=0;nr<n_r;nr++){ r_index_[nr] = nr;}
    flag_free_r_index = 1;
    /* if ( (n_r_index_0in<=0) || (r_index_0in_==NULL) ){ } */}
  if ( (n_c_index_0in<=0) || (c_index_0in_==NULL) ){
    n_c_index = n_c;
    c_index_ = (int *) malloc1(n_c*sizeof(int));
    for (nc=0;nc<n_c;nc++){ c_index_[nc] = nc;}
    flag_free_c_index = 1;
    /* if ( (n_r_index_0in<=0) || (r_index_0in_==NULL) ){ } */}
  if ( gamma_0in<=0 ){ gamma = 0.0; } else{ gamma = minimum(1,gamma_0in);}
  if (gamma> 0){ sprintf(prefix_gamma,"_g%.3d",floor(100*gamma));} else{ sprintf(prefix_gamma,"");}
  if (n_member_lob> 2){ sprintf(prefix_n_member_lob,"_n%.2d",n_member_lob);} else{ sprintf(prefix_n_member_lob,"");}
  sprintf(prefix_base1,"%s%s%s",prefix_base0,prefix_gamma,prefix_n_member_lob);
  sprintf(dir_out,"%s/dir_%s",dir_0in,prefix_base1);
  if (stat(dir_out,&stat_dir)==-1){ printf(" %% %s not found, creating\n",dir_out); mkdir(dir_out,0755);} else{ printf(" %% %s found, not creating\n",dir_out);}
  if (n_shuffle_0in<=0){ n_shuffle = 64;} else{ n_shuffle = n_shuffle_0in;}
  if (p_set_0in<=0){ p_set = 0.05;} else{ p_set = p_set_0in;}
  p_use = (double)p_set / (double)(1+2*p_set);
  if (n_member_lob_0in<=0){ n_member_lob = 2;} else{ n_member_lob = n_member_lob_0in;}
  if (p_prev_0in<=0){ p_prev = 0.00;} else{ p_prev = p_prev_0in;}
  if (flag_force_create_0in<=0){ flag_force_create = 0;} else{ flag_force_create = flag_force_create_0in;}
  if (verbose){
    printf(" %% dir_trunk: %s\n",dir_trunk);
    printf(" %% prefix_base0: %s\n",prefix_base0);
    printf(" %% prefix_base1: %s\n",prefix_base1);
    printf(" %% dir_0in: %s\n",dir_0in);
    printf(" %% dir_out: %s\n",dir_out);
    printf(" %% n_r %d n_c %d --> n_r_index %d n_c_index %d\n",n_r,n_c,n_r_index,n_c_index);
    printf(" %% gamma %0.3f n_shuffle %d p_set %0.3f (p_use %0.3f) n_member_lob %d p_prev %0.3f flag_force_create %d\n",gamma,n_shuffle,p_set,p_use,n_member_lob,p_prev,flag_force_create);
    /* if (verbose){ } */}
  sprintf(MDA_fname,"%s/r_index_.mda",dir_out); MDA_n_dim = 2; MDA_dim_[0] = n_r_index; MDA_dim_[1] = 1;
  MDA_write_i4(MDA_n_dim,MDA_dim_,r_index_,MDA_fname);
  if (verbose>2){ MDA_printf_i4_margin(MDA_fname);}
  sprintf(MDA_fname,"%s/c_index_.mda",dir_out); MDA_n_dim = 2; MDA_dim_[0] = n_c_index; MDA_dim_[1] = 1;
  MDA_write_i4(MDA_n_dim,MDA_dim_,c_index_,MDA_fname);
  if (verbose>2){ MDA_printf_i4_margin(MDA_fname);}
  /* %%%%%%%%%%%%%%%% */
  if ( (output_label_p_!=NULL) || (lpFmax_label_p_!=NULL) || (lpnext_label_p_!=NULL) ){
    /* extract labels */
    /* if extract labels */}
  /* %%%%%%%%%%%%%%%% */
  if ( (output_label_p_==NULL) && (lpFmax_label_p_==NULL) && (lpnext_label_p_==NULL) ){
    /* calculation required */
    sprintf(fname_trace_E__,"%s/trace_E__.mda",dir_out);
    sprintf(fname_QR_E__,"%s/QR_E__.mda",dir_out);
    sprintf(fname_xdrop_E__,"%s/xdrop_E__.mda",dir_out);
    dexcluster_nonbinary_f_recursive_helper_QR__
      (
       verbose
       ,0
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
       verbose
       ,0
       ,fname_trace_E__
       ,fname_xdrop_E__
       ,fname_QR_E__
       ,p_use
       ,0*n_member_lob
       ,&nlp_ZR_E_max
       ,&nlp_ZR_E_index
       ,&n_r_rtn_index_E
       ,&r_rtn_index_E_
       ,&n_c_rtn_index_E
       ,&c_rtn_index_E_
       );
    array_extract_i_from_i(1,n_c_index,c_index_,0,NULL,n_c_rtn_index_E,c_rtn_index_E_,&c_index_sub_,NULL);
    if (verbose>2){ array_printf_margin(c_index_sub_,"int",1,n_c_rtn_index_E," % c_index_sub_: ");}
    sprintf(fname_trace_F__,"%s/trace_F__.mda",dir_out);
    sprintf(fname_QR_F__,"%s/QR_F__.mda",dir_out);
    sprintf(fname_xdrop_F__,"%s/xdrop_F__.mda",dir_out);
    dexcluster_nonbinary_f_recursive_helper_QR__
      (
       verbose
       ,1
       ,flag_force_create
       ,n_r
       ,n_c
       ,E_base_rc__
       ,n_r_index
       ,r_index_
       ,n_c_rtn_index_E
       ,c_index_sub_
       ,0*gamma
       ,n_shuffle
       ,fname_trace_F__
       ,fname_xdrop_F__
       ,fname_QR_F__
       );
    dexcluster_nonbinary_recursive_helper_ZR__
      (
       verbose
       ,1
       ,fname_trace_F__
       ,fname_xdrop_F__
       ,fname_QR_F__
       ,p_use
       ,1*n_member_lob
       ,&nlp_ZR_F_max
       ,&nlp_ZR_F_index
       ,&n_r_rtn_index_F
       ,&r_rtn_index_F_
       ,&n_c_rtn_index_F
       ,&c_rtn_index_F_
       );
    free1(&c_index_sub_);
    free1(&r_rtn_index_E_);
    free1(&c_rtn_index_E_);
    free1(&r_rtn_index_F_);
    free1(&c_rtn_index_F_);
    /* if calculation required */}  
  /* %%%%%%%%%%%%%%%% */
  if (flag_free_r_index){ free1(&r_index_);}
  if (flag_free_c_index){ free1(&c_index_);}
  free1(&cwd);
  free1(&MDA_dim_);
  /* %%%%%%%%%%%%%%%% */
  if (verbose){ printf(" %% [finished dexcluster_nonbinary_f_recursive]\n");}
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
 ,NULL
 ,NULL
 ,NULL
 );
  GLOBAL_toc(0,1," % dexcluster_nonbinary_recursive: ");
  free1(&E_base_rc__);
}
