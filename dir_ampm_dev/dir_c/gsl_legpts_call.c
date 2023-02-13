#ifndef _MONOLITH
#include "ampm_header.h"
#endif /* _MONOLITH */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* testing gsl_legpts_call.c */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
void gsl_legpts(int n_node,double **x_p_,double **w_p_)
{
  int nnode=0;
  double *x_=NULL;
  double *w_=NULL;
  gsl_integration_glfixed_table *t=NULL;
  if (x_p_!=NULL){
    if ((*x_p_)==NULL){ (*x_p_) = (double *) malloc1(n_node*sizeof(double));}
    x_ = *x_p_;
    /* if (x_p_!=NULL){ } */}
  if (w_p_!=NULL){
    if ((*w_p_)==NULL){ (*w_p_) = (double *) malloc1(n_node*sizeof(double));}
    w_ = *w_p_;
    /* if (w_p_!=NULL){ } */}
  if ( (x_!=NULL) && (w_!=NULL)){
    t = gsl_integration_glfixed_table_alloc(n_node);
    for (nnode=0;nnode<n_node;nnode++){
      gsl_integration_glfixed_point(-1.0,+1.0,nnode,&(x_[nnode]),&(w_[nnode]),t);
      /* for (nnode=0;nnode<n_node;nnode++){ } */}
    gsl_integration_glfixed_table_free(t); t=NULL;
    /* if ( (x_!=NULL) && (w_!=NULL)){ } */}
}

void gsl_legpts_test()
{
  int n_node = 16;
  int nnode=0;
  double *x_ = NULL;
  double *w_ = NULL;
  double x_ans_[16] = { -0.9894009349916499 , -0.9445750230732326 , -0.8656312023878318 , -0.7554044083550030 , -0.6178762444026438 , -0.4580167776572274 , -0.2816035507792589 , -0.0950125098376374 , +0.0950125098376374 , +0.2816035507792589 , +0.4580167776572274 , +0.6178762444026438 , +0.7554044083550030 , +0.8656312023878318 , +0.9445750230732326 , +0.9894009349916499 };
  double w_ans_[16] = { +0.0271524594117541 , +0.0622535239386478 , +0.0951585116824928 , +0.1246289712555339 , +0.1495959888165768 , +0.1691565193950025 , +0.1826034150449235 , +0.1894506104550685 , +0.1894506104550685 , +0.1826034150449235 , +0.1691565193950025 , +0.1495959888165768 , +0.1246289712555339 , +0.0951585116824928 , +0.0622535239386478 , +0.0271524594117541 };
  /* %%%%%%%%; */
  GLOBAL_tic(0);
  gsl_legpts(n_node,&x_,&w_);
  array_printf(x_,"double",1,n_node," %% x_: ");
  array_printf(w_,"double",1,n_node," %% w_: ");
  printf(" %% x_ vs x_ans_: %0.16f\n",dfnormn((unsigned long long int)n_node,x_,x_ans_));
  printf(" %% w_ vs w_ans_: %0.16f\n",dfnormn((unsigned long long int)n_node,w_,w_ans_));
  free1(&x_);
  free1(&w_);
  GLOBAL_toc(0,1," gsl_legpts: ");
  /* %%%%%%%%; */
  GLOBAL_tic(0);
  double K_max = 7.0;
  double delta = 1.2;
  double Ix = 4*PI*(1/pow(delta,3))*( sin(K_max*delta) - (K_max*delta)*cos(K_max*delta) );
  n_node = 6;
  gsl_legpts(n_node,&x_,&w_);
  double *k_lx_ = NULL;
  double *k_lw_ = NULL;
  double *Ix_ = NULL;
  double k=0.0;
  double Il=0.0;
  k_lx_ = (double *) malloc1(n_node*sizeof(double));
  k_lw_ = (double *) malloc1(n_node*sizeof(double));
  Ix_ = (double *) malloc1(n_node*sizeof(double));
  Il=0.0;
  for (nnode=0;nnode<n_node;nnode++){
    k_lx_[nnode] = (x_[nnode] + 1.0)*K_max/2.0;
    k_lw_[nnode] = w_[nnode]*K_max/2.0;
    k = k_lx_[nnode];
    Ix_[nnode] = 4*PI*(1.0/delta)*sin(delta*k)*k;
    Il += k_lw_[nnode]*Ix_[nnode];
    /* for (nnode=0;nnode<n_node;nnode++){ } */}
  printf(" %% (Ix-Il)/Ix = %0.16f/%0.16f = %0.16f\n",Ix-Il,Ix,(Ix-Il)/Ix);
  free1(&k_lx_);
  free1(&k_lw_);
  free1(&Ix_);
  free1(&x_);
  free1(&w_);
  GLOBAL_toc(0,1," Test quadrature with r^2 weight: ");
  /* %%%%%%%%; */
}

