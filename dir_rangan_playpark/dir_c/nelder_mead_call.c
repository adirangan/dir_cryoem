#ifndef _MONOLITH
#include "playpark_header.h"
#endif /* _MONOLITH */

void nelder_mead_terminate(int n_d,double *simplex_point__,int *index_simplex_point_,double *simplex_cost_,double option_tolx,double option_tolf,int *flag_continue_p)
{
  int n_s=n_d+1,ns=0,nd=0;
  int ns_0best=0;
  double cost_difference=0;
  double dist_difference=0,tmp_dist=0,tmp_d=0;
  int flag_continue;
  double *simplex_point_=NULL,*simplex_0best_=NULL;
  ns_0best = index_simplex_point_[0];
  simplex_0best_ = simplex_point__ + ns_0best*n_d;
  cost_difference=-1; //%<-- approximates gradient across simplex. ;
  for (ns=1;ns<n_s;ns++){
    cost_difference = maximum(cost_difference,simplex_cost_[index_simplex_point_[ns]]-simplex_cost_[ns_0best]);
    /* for (ns=1;ns<n_s;ns++){ } */}
  dist_difference=-1; //%<-- approximates diameter of simplex. ;
  for (ns=1;ns<n_s;ns++){
    tmp_dist=0;
    simplex_point_ = simplex_point__ + index_simplex_point_[ns]*n_d;
    for (nd=0;nd<n_d;nd++){
      tmp_d = simplex_point_[nd] - simplex_0best_[nd];
      tmp_dist += tmp_d*tmp_d;
      /* for (nd=0;nd<n_d;nd++){ } */}
    tmp_dist = sqrt(tmp_dist);
    dist_difference = maximum(dist_difference,tmp_dist);
    /* for (ns=1;ns<n_s;ns++){ } */}
  flag_continue=1;
  if (cost_difference<=option_tolf){ flag_continue=0;}
  if (dist_difference<=option_tolx){ flag_continue=0;}
  *flag_continue_p = flag_continue;
}

void nelder_mead_simplex_centroid(int n_d,double *simplex_point__,int *index_simplex_point_,double *simplex_centroid_)
{
  /* calculate centroid of the best n_s-1 points in the simplex */
  int n_s=n_d+1,ns=0,nd=0;
  double *simplex_point_=NULL;
  memset(simplex_centroid_,0,n_d*sizeof(double));
  for (ns=0;ns<n_s-1;ns++){
    simplex_point_ = simplex_point__ + index_simplex_point_[ns]*n_d;
    for (nd=0;nd<n_d;nd++){ simplex_centroid_[nd] += simplex_point_[nd];}
    /* for (ns=0;ns<n_s-1;ns++){ } */}
  for (nd=0;nd<n_d;nd++){ simplex_centroid_[nd] /= maximum(1.0,(double)(n_s-1));}
}

void nelder_mead_update_point(int n_d,double *simplex_point_,double *simplex_centroid_,double lambda,double *point_)
{
  int nd=0;
  for (nd=0;nd<n_d;nd++){ point_[nd] = lambda*simplex_point_[nd] + (1.0-lambda)*simplex_centroid_[nd];}
}

void nelder_mead_optimization(int n_d,double *point_start_,double *point_final_, void cost_function(int,double *,double *,void *),void *cost_function_args,double option_tolx,double option_tolf,int option_maxiter,int option_maxfeval)
{
  int verbose=0;
  int nd=0;
  int n_iteration=0,niteration=0;
  int n_functeval=0,nfuncteval=0;
  int n_s = n_d+1; //%<-- number of points in simplex ;
  int ns=0,ns_worst=0,ns_2ndwr=0,ns_0best=0;
  double cost_start=0;
  double *point_0reflect_=NULL; double cost_0reflect=0;
  double *point_expanded_=NULL; double cost_expanded=0;
  double *point_contract_=NULL; double cost_contract=0;
  double *simplex_centroid_=NULL; double cost_centroid=0;
  double *simplex_point__ = NULL;
  double *simplex_point_ = NULL;
  double *simplex_0best_ = NULL;double cost_0best=0;
  double *simplex_worst_ = NULL;double cost_worst=0;
  double *simplex_2ndwr_ = NULL;double cost_2ndwr=0;
  double *simplex_cost_ = NULL;
  double *simplex_cost_workspace_=NULL;
  int *index_simplex_point_=NULL;
  int flag_continue=1;
  double tmp_d=0,tmp_dist=0,dist_difference=0,cost_difference=0;
  double nelder_mead_alpha = 1.0;
  double nelder_mead_gamma = 2.0;
  double nelder_mead_00rho = 0.5;
  double nelder_mead_sigma = 0.5;  
  if (verbose){ printf(" %% [entering nelder_mead_optimization]\n");}
  /* %%%%%%%%%%%%%%%% */
  if (option_tolx<=0){ option_tolx = 1e-6;}
  if (option_tolf<=0){ option_tolf = 1e-6;}
  if (option_maxiter<=0){ option_maxiter = 128;}
  if (option_maxfeval<=0){ option_maxfeval = 128*n_d;}
  n_iteration = option_maxiter;
  n_functeval = option_maxfeval;
  /* %%%%%%%%%%%%%%%% */
  if (verbose>1){
    printf(" %% n_d %d, option_tolx %f option_tolf %f option_maxiter %d option_maxfeval %d\n",n_d,option_tolx,option_tolf,option_maxiter,option_maxfeval);
    array_printf(point_start_,"double",1,n_d," % point_start_: ");
    cost_function(n_d,point_start_,&cost_start,cost_function_args);
    printf(" %% cost_start %f\n",cost_start);
    /* if (verbose){ } */}
  /* %%%%%%%%%%%%%%%% */
  point_0reflect_ = (double *) malloc1(n_d*sizeof(double));
  point_expanded_ = (double *) malloc1(n_d*sizeof(double));
  point_contract_ = (double *) malloc1(n_d*sizeof(double));
  simplex_centroid_ = (double *) malloc1(n_d*sizeof(double));
  /* %%%%%%%%%%%%%%%% */
  /* initialize simplex and sort simplex_point__ by simplex_cost_ */
  simplex_point__ = (double *) malloc1(n_d*n_s*sizeof(double));
  simplex_cost_ = (double *) malloc1(n_s*sizeof(double));
  simplex_cost_workspace_ = (double *) malloc1(n_s*sizeof(double));
  index_simplex_point_ = (int *) malloc1(n_s*sizeof(int));
  for (ns=0;ns<n_s;ns++){
    simplex_point_ = simplex_point__ + ns*n_d;
    memcpy(simplex_point_,point_start_,n_d*sizeof(double));
    if (ns>0){ nd=ns-1; simplex_point_[nd] *= 1.05; if (simplex_point_[nd]==0){ simplex_point_[nd]=0.00025;}}
    cost_function(n_d,simplex_point_,&(simplex_cost_[ns]),cost_function_args);
    if (verbose>2){ printf(" %% simplex_point %d, cost %f\n",ns,simplex_cost_[ns]); array_printf(simplex_point_,"double",1,n_d," % simplex_point_: ");}
    /* for (ns=0;ns<n_s;ns++){ } */}
  if (verbose>2){
    dquicksort_index_driver(n_s,simplex_cost_,1,simplex_cost_workspace_,index_simplex_point_);
    array_printf(simplex_cost_workspace_,"double",1,n_s," % simplex_cost_workspace_: ");
    array_printf(index_simplex_point_,"int",1,n_s," % index_simplex_point_: ");
    /* if (verbose){ } */}
  /* %%%%%%%%%%%%%%%% */
  do{
    /* %%%%%%%% */
    dquicksort_index_driver(n_s,simplex_cost_,1,simplex_cost_workspace_,index_simplex_point_);
    if (verbose>2){
      printf(" %% niteration %d/%d nfuncteval %d/%d\n",niteration,n_iteration,nfuncteval,n_functeval);
      array_printf(simplex_cost_workspace_,"double",1,n_s," % simplex_cost_workspace_: ");
      array_printf(index_simplex_point_,"int",1,n_s," % index_simplex_point_: ");
      ns_0best = index_simplex_point_[0];
      simplex_0best_ = simplex_point__ + ns_0best*n_d; cost_0best = simplex_cost_[ns_0best];
      array_printf(simplex_0best_,"double",1,n_d," simplex_0best_: ");
      /* if (verbose){ } */}
    nelder_mead_terminate(n_d,simplex_point__,index_simplex_point_,simplex_cost_,option_tolx,option_tolf,&flag_continue);
    if (niteration>=n_iteration){ flag_continue=0;} if (nfuncteval>=n_functeval){ flag_continue=0;}
    if (flag_continue){
      /* %%%%%%%% */    
      nelder_mead_simplex_centroid(n_d,simplex_point__,index_simplex_point_,simplex_centroid_);
      ns_0best = index_simplex_point_[0];
      simplex_0best_ = simplex_point__ + ns_0best*n_d; cost_0best = simplex_cost_[ns_0best];
      ns_worst = index_simplex_point_[n_s-1];
      simplex_worst_ = simplex_point__ + ns_worst*n_d; cost_worst = simplex_cost_[ns_worst];
      ns_2ndwr = index_simplex_point_[n_s-1-1];
      simplex_2ndwr_ = simplex_point__ + ns_2ndwr*n_d;  cost_2ndwr = simplex_cost_[ns_2ndwr];
      nelder_mead_update_point(n_d,simplex_worst_,simplex_centroid_,-nelder_mead_alpha,point_0reflect_);
      cost_function(n_d,point_0reflect_,&(cost_0reflect),cost_function_args); nfuncteval+=1;
      /* %%%%%%%% */    
      if (0){ /* do nothing */}
      else if ((cost_0best<=cost_0reflect) && (cost_0reflect< cost_2ndwr)){ /* reflected point better than 2ndwr, but not better than 0best */
	if (verbose>2){ printf(" %% reflection: replacing worst point. \n");}
	memcpy(simplex_worst_,point_0reflect_,n_d*sizeof(double)); simplex_cost_[ns_worst] = cost_0reflect;
	/* if ((cost_0best<=cost_0reflect) && (cost_0reflect< cost_2ndwr)){ } */}
      else if (cost_0reflect<cost_0best){ /* reflected point best so far */
	if (verbose>2){ printf(" %% expansion... \n");}
	nelder_mead_update_point(n_d,point_0reflect_,simplex_centroid_,+nelder_mead_gamma,point_expanded_);
	cost_function(n_d,point_expanded_,&(cost_expanded),cost_function_args); nfuncteval+=1;
	if (cost_expanded <  cost_0reflect){
	  if (verbose>2){ printf(" expansion: better than reflected point. \n");}
	  memcpy(simplex_worst_,point_expanded_,n_d*sizeof(double)); simplex_cost_[ns_worst] = cost_expanded;
	  /* if (cost_expanded <  cost_0reflect){ } */}
	if (cost_expanded >= cost_0reflect){
	  if (verbose>2){ printf(" %% expansion: worse than reflected point. \n");}
	  memcpy(simplex_worst_,point_0reflect_,n_d*sizeof(double)); simplex_cost_[ns_worst] = cost_0reflect;
	  /* if (cost_expanded >= cost_0reflect){ } */}
	/* if (cost_0reflect<cost_0best){ } */}
      else if (cost_0reflect>=cost_2ndwr){ /* reflected point at least as bad as second worst so far */
	if (verbose>2){ printf(" %% contraction... \n");}
	nelder_mead_update_point(n_d,simplex_worst_,simplex_centroid_,+nelder_mead_00rho,point_contract_);
	cost_function(n_d,point_contract_,&(cost_contract),cost_function_args); nfuncteval+=1;
	if (cost_contract <  cost_worst){
	  if (verbose>2){ printf(" contraction: better than worst point. \n");}
	  memcpy(simplex_worst_,point_contract_,n_d*sizeof(double)); simplex_cost_[ns_worst] = cost_contract;
	  /* if (cost_contract <  cost_worst){ } */}
	if (cost_contract >= cost_worst){ /* shrink */
	  if (verbose>2){ printf(" contraction: shrink simplex. \n");}
	  for (ns=1;ns<n_s;ns++){
	    simplex_point_ = simplex_point__ + index_simplex_point_[ns]*n_d;
	    nelder_mead_update_point(n_d,simplex_point_,simplex_0best_,+nelder_mead_sigma,simplex_point_);
	    cost_function(n_d,simplex_point_,&(simplex_cost_[index_simplex_point_[ns]]),cost_function_args); nfuncteval+=1;
	    /* for (ns=1;ns<n_s;ns++){ } */}
	  /* if (cost_contract >= cost_worst){ } */}
	/* if (cost_0reflect>=cost_2ndwr){ } */}    
      /* %%%%%%%% */
      /* if (flag_continue){ } */}
    niteration+=1;
  } while (flag_continue);
  if (verbose>1){ printf(" %% terminated with niteration %d/%d nfuncteval %d/%d cost_difference %f/%f dist_difference %f/%f\n",niteration,n_iteration,nfuncteval,n_functeval,cost_difference,option_tolf,dist_difference,option_tolx);}
  /* %%%%%%%%%%%%%%%% */
  ns_0best = index_simplex_point_[0];
  simplex_0best_ = simplex_point__ + ns_0best*n_d; cost_0best = simplex_cost_[ns_0best];
  memcpy(point_final_,simplex_0best_,n_d*sizeof(double));
  /* %%%%%%%%%%%%%%%% */
  free1(&point_0reflect_);
  free1(&point_expanded_);
  free1(&point_contract_);
  free1(&simplex_centroid_);
  free1(&simplex_point__);
  free1(&simplex_cost_);
  free1(&simplex_cost_workspace_);
  free1(&index_simplex_point_);
  if (verbose){ printf(" %% [finished nelder_mead_optimization]\n");}
}

void nelder_mead_test()
{
  int n_x = 8,n_g = 2,nx=0;
  double *x_ = NULL;
  double *g_start_ = NULL;
  double *g_final_ = NULL;
  double *nll_ = NULL;
  double nll_sum=0;
  void **gumbel_nll_wrap_args__=NULL;
  double option_tolx=0;
  double option_tolf=0;
  int option_maxiter=0;
  int option_maxfeval=0;
  double g_final_ans_[2] = { 2.8600292585746505 , 1.4369016493495146 };
  x_ = (double *) malloc1(n_x*sizeof(double));
  g_start_ = (double *) malloc1(n_g*sizeof(double));
  g_final_ = (double *) malloc1(n_g*sizeof(double));
  nll_ = (double *) malloc1(n_x*sizeof(double));
  gumbel_nll_wrap_args__ = (void **) malloc1(3*sizeof(void *));
  g_start_[0] = 0.2; g_start_[1] = 1.3;
  memset(g_final_,0,n_g*sizeof(double));
  for (nx=0;nx<n_x;nx++){ x_[nx] = nx+0.5;}
  memset(nll_,0,n_x*sizeof(double));
  gumbel_nll_wrap_args__[0] = &n_x;
  gumbel_nll_wrap_args__[1] = x_;
  gumbel_nll_wrap_args__[2] = nll_;
  GLOBAL_tic(0);
  nelder_mead_optimization(n_g,g_start_,g_final_,gumbel_nll_wrap,gumbel_nll_wrap_args__,option_tolx,option_tolf,option_maxiter,option_maxfeval);
  array_printf(g_final_,"double",1,n_g," % g_final_: ");
  printf(" %% g_final_ans_ vs g_final_: relative error %0.16f\n",dfnormn(n_g,g_final_ans_,g_final_));
  GLOBAL_toc(0,1," nelder_mead_optimization: ");
  free1(&x_);
  free1(&g_start_);
  free1(&g_final_);
  free1(&nll_);
  free1(&gumbel_nll_wrap_args__);
}
