#ifndef _MONOLITH
#include "playpark_header.h"
#endif /* _MONOLITH */

void halfloop_recursive_f(
 char *dir_trunk_0in
 ,char *dir_out_0in
 ,char *prefix_base_0in
 ,int n_r
 ,int n_c
 ,float *E_array_base_
 ,int n_E_array_r_index_0in
 ,int *E_array_r_index_0in_
 ,int n_E_array_c_index_0in
 ,int *E_array_c_index_0in_
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
  int verbose=1;
  char *cwd=NULL;
  int n_str_path=0;
  char dir_trunk[PNAMESIZE];
  char dir_out[PNAMESIZE];
  char prefix_base[FNAMESIZE];
  char prefix_string[FNAMESIZE];
  int nr=0,flag_free_E_array_r_index=0,n_E_array_r_index=0,*E_array_r_index_=NULL;
  int nc=0,flag_free_E_array_c_index=0,n_E_array_c_index=0,*E_array_c_index_=NULL;
  double gamma=0;
  if (verbose){ printf(" %% [entering halfloop_recursive_f]\n");}
  /* %%%%%%%%%%%%%%%% */
  n_str_path = pathconf(".",_PC_PATH_MAX);
  cwd = (char *) malloc1((size_t)n_str_path);
  getcwd(cwd,n_str_path);
  if (verbose){ printf(" %% cwd: [%s]\n",cwd);}
  if (n_str_path> FNAMESIZE){ printf(" %% Warning, n_str_path %d in halfloop_recursive_f\n",n_str_path);}
  if ( (dir_trunk_0in!=NULL) && (strlen(dir_trunk_0in)> FNAMESIZE) ){ printf(" %% Warning, dir_trunk_0in %s too long in halfloop_recursive_f\n",dir_trunk_0in);}
  if ( (dir_out_0in!=NULL) && (strlen(dir_out_0in)> FNAMESIZE) ){ printf(" %% Warning, dir_out_0in %s too long in halfloop_recursive_f\n",dir_out_0in);}
  if ( (prefix_base_0in!=NULL) && (strlen(prefix_base_0in)> FNAMESIZE) ){ printf(" %% Warning, prefix_base_0in %s too long in halfloop_recursive_f\n",prefix_base_0in);}
  if ( dir_trunk_0in==NULL ){ sprintf(dir_trunk,"%s",cwd);} else{ sprintf(dir_trunk,"%s",dir_trunk_0in);}
  if ( prefix_base_0in==NULL ){ sprintf(prefix_base,"%s","test");} else{ sprintf(prefix_base,"%s",prefix_base_0in);}
  if ( (n_E_array_r_index_0in<=0) || (E_array_r_index_0in_==NULL) ){
    n_E_array_r_index = n_r;
    E_array_r_index_ = (int *) malloc1(n_r*sizeof(int));
    for (nr=0;nr<n_r;nr++){ E_array_r_index_[nr] = nr;}
    flag_free_E_array_r_index = 1;
    /* if ( (n_E_array_r_index_0in<=0) || (E_array_r_index_0in_==NULL) ){ } */}
  if ( (n_E_array_c_index_0in<=0) || (E_array_c_index_0in_==NULL) ){
    n_E_array_c_index = n_c;
    E_array_c_index_ = (int *) malloc1(n_c*sizeof(int));
    for (nc=0;nc<n_c;nc++){ E_array_c_index_[nc] = nc;}
    flag_free_E_array_c_index = 1;
    /* if ( (n_E_array_r_index_0in<=0) || (E_array_r_index_0in_==NULL) ){ } */}
  if ( gamma_0in<=0 ){ gamma = 0.0; } else{ gamma = minimum(1,gamma_0in);}
  if (gamma> 0){ sprintf(prefix_string,"%s_g%.3d",prefix_base,floor(100*gamma));}
  if (gamma==0){ sprintf(prefix_string,"%s",prefix_base);}
  if (verbose){
    printf(" %% dir_trunk %s\n",dir_trunk);
    printf(" %% prefix_base %s\n",prefix_base);
    printf(" %% prefix_string %s\n",prefix_string);
    printf(" %% n_r %d n_c %d --> n_E_array_r_index %d n_E_array_c_index %d\n",n_r,n_c,n_E_array_r_index,n_E_array_c_index);
    /* if (verbose){ } */}
  /* %%%%%%%%%%%%%%%% */
  /* %%%%%%%%%%%%%%%% */
  if (flag_free_E_array_r_index){ free1(&E_array_r_index_);}
  if (flag_free_E_array_c_index){ free1(&E_array_c_index_);}
  free1(&cwd);
  /* %%%%%%%%%%%%%%%% */
  if (verbose){ printf(" %% [finished halfloop_recursive_f]\n");}
}

void halfloop_recursive_f_test()
{
  GLOBAL_tic(0);
  halfloop_recursive_f(
  NULL
 ,NULL
 ,NULL
 ,0
 ,0
 ,NULL
 ,0
 ,NULL
 ,0
 ,NULL
 ,0
 ,0
 ,0
 ,0
 ,0
 ,0
 ,NULL
 ,NULL
 ,NULL
 );
  GLOBAL_toc(0,1," % halfloop_recursive: ");
}
