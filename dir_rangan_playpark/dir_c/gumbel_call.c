#ifndef _MONOLITH
#include "playpark_header.h"
#endif /* _MONOLITH */

void gumbel_nll(int n_x,double *x_,double *g_,double tol,double **nll_p_,double *nll_sum_p)
{
  double mu=0,beta=0,z=0,nll=0;
  int nx=0;
  if (tol<=0){ tol=1e-6;}
  mu = g_[0]; beta = tol + maximum(0,g_[1])*maximum(0,g_[1]);
  if ((*nll_p_)==NULL){ (*nll_p_) = (double *) malloc1(n_x*sizeof(double)); }
  *nll_sum_p = 0;
  for (nx=0;nx<n_x;nx++){
    z = (x_[nx] - mu)/beta;
    nll = z + exp(-z) + log(beta);
    (*nll_p_)[nx] = nll;
    (*nll_sum_p) += nll;
    /* for (nx=0;nx<n_x;nx++){ } */}
}

void gumbel_nll_wrap(int n_g,double *g_,double *nll_sum_p,void *args)
{
  void **vp__ = args;
  int *n_x_p = (int *)(vp__[0]);
  int n_x = *n_x_p;
  double *x_ = (double *)(vp__[1]);
  double *nll_ = (double *)(vp__[2]); //%<-- workspace. ;
  gumbel_nll(n_x,x_,g_,0,&nll_,nll_sum_p);
}

void gumbel_nll_test()
{
  int n_x = 8,n_g = 2,nx=0;
  double *x_ = NULL;
  double *g_ = NULL;
  double tol = 1e-6;
  double *nll_ = NULL;
  double **nll_p_ = NULL;
  double nll_sum=0;
  double *nll_sum_p=NULL;
  void **vp__=NULL;
  double nll_ans_[8] = { 1.5395925088959164 , 1.7573290148568299 , 2.1420931669029204 , 2.6192865902943407 , 3.1476282541217899 , 3.7042741819732137 , 4.2765830408484504 , 4.8575594067837145 };
  double nll_sum_ans = 24.0443461646771723;
  x_ = (double *) malloc1(n_x*sizeof(double));
  g_ = (double *) malloc1(n_g*sizeof(double));
  nll_ = (double *) malloc1(n_x*sizeof(double));
  vp__ = (void **) malloc1(3*sizeof(void *));
  nll_p_ = &nll_;
  nll_sum_p = &nll_sum;
  g_[0] = 0.2; g_[1] = 1.3;
  for (nx=0;nx<n_x;nx++){ x_[nx] = nx+0.5;}
  memset(nll_,0,n_x*sizeof(double));
  vp__[0] = &n_x;
  vp__[1] = x_;
  vp__[2] = nll_;
  GLOBAL_tic(0);
  gumbel_nll(n_x,x_,g_,tol,nll_p_,nll_sum_p);
  printf(" %% n_x %d, g_ %0.2f %0.2f\n",n_x,g_[0],g_[1]);
  array_printf(x_,"double",1,n_x," % x_: ");  
  array_printf(nll_,"double",1,n_x," % nll_: ");
  printf(" %% nll_sum %0.2f\n",nll_sum);
  printf(" %% nll_ans_ vs nll_: relative error %0.16f\n",dfnormn(n_x,nll_,nll_ans_));
  GLOBAL_toc(0,1," gumbel_nll: ");
  GLOBAL_tic(0);
  memset(nll_,0,n_x*sizeof(double));
  gumbel_nll_wrap(n_g,g_,nll_sum_p,vp__);
  printf(" %% n_x %d, g_ %0.2f %0.2f\n",n_x,g_[0],g_[1]);
  array_printf(x_,"double",1,n_x," % x_: ");
  array_printf(nll_,"double",1,n_x," % nll_: ");
  printf(" %% nll_sum %0.2f\n",nll_sum);
  printf(" %% nll_ans_ vs nll_: relative error %0.16f\n",dfnormn(n_x,nll_,nll_ans_));
  GLOBAL_toc(0,1," gumbel_nll_wrap: ");
  free1(&x_);
  free1(&g_);
  free1(&nll_);
  free1(&vp__);
}

void gumbel_pdf(int n_x,double *x_,double *g_,double tol,double **pdf_p_)
{
  double mu=0,beta=0,z=0,pdf=0;
  int nx=0;
  if (tol<=0){ tol=1e-6;}
  mu = g_[0]; beta = tol + maximum(0,g_[1])*maximum(0,g_[1]);
  if ((*pdf_p_)==NULL){ (*pdf_p_) = (double *) malloc1(n_x*sizeof(double)); }
  for (nx=0;nx<n_x;nx++){
    z = (x_[nx] - mu)/beta;
    pdf = exp( -(z + exp(-z)) ) / beta ;
    (*pdf_p_)[nx] = pdf;
    /* for (nx=0;nx<n_x;nx++){ } */}
}

void gumbel_cdf(int n_x,double *x_,double *g_,double tol,double **cdf_p_)
{
  double mu=0,beta=0,z=0,cdf=0;
  int nx=0;
  if (tol<=0){ tol=1e-6;}
  mu = g_[0]; beta = tol + maximum(0,g_[1])*maximum(0,g_[1]);
  if ((*cdf_p_)==NULL){ (*cdf_p_) = (double *) malloc1(n_x*sizeof(double)); }
  for (nx=0;nx<n_x;nx++){
    z = (x_[nx] - mu)/beta;
    cdf = exp(-exp(-z));
    (*cdf_p_)[nx] = cdf;
    /* for (nx=0;nx<n_x;nx++){ } */}
}

double gumbel_cdf_single(double x,double *g_,double tol)
{
  double mu=0,beta=0,z=0,cdf=0;
  int nx=0;
  if (tol<=0){ tol=1e-6;}
  mu = g_[0]; beta = tol + maximum(0,g_[1])*maximum(0,g_[1]);
  z = (x - mu)/beta;
  return exp(-exp(-z));
}

void gumbel_fit(int n_x,double *x_,double x,double **g_opt_p_,double *nlp_opt_p,double *nlp_emp_p,double *p_opt_p,double *p_emp_p)
{
  int verbose=0;
  int nx=0;
  double nlp_opt=0,nlp_emp=0,p_opt=0,p_emp=0;
  double x_avg=0,x_std=0;
  double *y_=NULL,y=0;
  double *nll_=NULL;
  int n_g=2;
  double *g_start_=NULL;
  double *g_final_=NULL;
  void **gumbel_nll_wrap_args__=NULL;
  double option_tolx=0;
  double option_tolf=0;
  int option_maxiter=1e4;
  int option_maxfeval=1e4;
  double tmp_cdf=0;
  /* %%%%%%%%%%%%%%%% */
  if (verbose){ printf(" %% [entering gumbel_fit]\n");}
  /* %%%%%%%%%%%%%%%% */
  array_stats(x_,"double",n_x,NULL,NULL,&x_avg,&x_std);
  y_ = (double *) malloc1(n_x*sizeof(double));
  nll_ = (double *) malloc1(n_x*sizeof(double));
  for (nx=0;nx<n_x;nx++){ y_[nx] = (x_[nx]-x_avg)/maximum(1e-12,x_std);}
  if (verbose>1){ array_printf(y_,"double",1,n_x," % y_: ");}
  g_start_ = (double *) malloc1(n_g*sizeof(double));
  g_final_ = (double *) malloc1(n_g*sizeof(double));
  g_start_[0] = 0.0; g_start_[1] = 1.0;
  if (verbose>1){ array_printf(g_start_,"double",1,n_g," % g_start_: ");}
  gumbel_nll_wrap_args__ = (void **) malloc1(3*sizeof(void *));
  gumbel_nll_wrap_args__[0] = &n_x;
  gumbel_nll_wrap_args__[1] = y_;
  gumbel_nll_wrap_args__[2] = nll_;
  nelder_mead_optimization(n_g,g_start_,g_final_,gumbel_nll_wrap,gumbel_nll_wrap_args__,option_tolx,option_tolf,option_maxiter,option_maxfeval);
  if (verbose>1){ array_printf(g_final_,"double",1,n_g," % g_final_: ");}
  /* %%%%%%%%%%%%%%%% */
  nlp_opt=0;nlp_emp=0;
  p_emp=0; for (nx=0;nx<n_x;nx++){ p_emp += (x_[nx]> x ? 1.0 : (x_[nx]==x ? 0.5 : 0));} p_emp /= maximum(1,n_x);
  nlp_emp = -log(maximum(1e-12,p_emp));
  y = (x-x_avg)/maximum(1e-12,x_std);
  tmp_cdf = gumbel_cdf_single(y,g_final_,0);
  p_opt = (double)1.0 - tmp_cdf;
  nlp_opt = -log(maximum(1e-12,p_opt));
  /* %%%%%%%%%%%%%%%% */
  if (g_opt_p_!=NULL){ (*g_opt_p_)[0] = g_final_[0]; (*g_opt_p_)[1] = g_final_[1];}
  if (nlp_opt_p!=NULL){ *nlp_opt_p = nlp_opt;}
  if (nlp_emp_p!=NULL){ *nlp_emp_p = nlp_emp;}
  if (p_opt_p!=NULL){ *p_opt_p = p_opt;}
  if (p_emp_p!=NULL){ *p_emp_p = p_emp;}
  /* %%%%%%%%%%%%%%%% */
  free1(&gumbel_nll_wrap_args__);
  free1(&g_start_);
  free1(&g_final_);
  free1(&nll_);
  free1(&y_);
  /* %%%%%%%%%%%%%%%% */
  if (verbose){ printf(" %% [finished gumbel_fit]\n");}
}

void gumbel_fit_test()
{
  int n_x = 8,n_g = 2,nx=0;
  double *x_ = NULL;
  double x=0;
  double *g_opt_ = NULL;
  double nlp_opt=0,nlp_emp=0,p_opt=0,p_emp=0;
  double g_opt_ans_[2] = { -0.4975655217194898 , 0.9492539492472620 };
  double nlp_opt_ans = 1.6281473110900229;
  double nlp_emp_ans = 1.3862943611198906;
  double p_opt_ans = 0.1962929071436291;
  double p_emp_ans = 0.2500000000000000;
  x_ = (double *) malloc1(n_x*sizeof(double));
  g_opt_ = (double *) malloc1(n_g*sizeof(double));
  for (nx=0;nx<n_x;nx++){ x_[nx] = nx+0.5;}
  x = 6;
  GLOBAL_tic(0);
  gumbel_fit(n_x,x_,x,&g_opt_,&nlp_opt,&nlp_emp,&p_opt,&p_emp);
  printf(" %% g_opt_ %0.16f %0.16f nlp_opt %0.16f nlp_emp %0.16f p_opt %0.16f p_emp %0.16f\n",g_opt_[0],g_opt_[1],nlp_opt,nlp_emp,p_opt,p_emp);
  printf(" %% g_opt_ans_ vs g_opt_: relative error %0.16f\n",dfnormn(n_g,g_opt_ans_,g_opt_));
  printf(" %% nlp_opt_ans vs nlp_opt: error %0.16f\n",fabs(nlp_opt_ans - nlp_opt));
  printf(" %% nlp_emp_ans vs nlp_emp: error %0.16f\n",fabs(nlp_emp_ans - nlp_emp));
  printf(" %% p_opt_ans vs p_opt: error %0.16f\n",fabs(p_opt_ans - p_opt));
  printf(" %% p_emp_ans vs p_emp: error %0.16f\n",fabs(p_emp_ans - p_emp));
  GLOBAL_toc(0,1," gumbel_fit: ");
  free1(&x_);
  free1(&g_opt_);
}


