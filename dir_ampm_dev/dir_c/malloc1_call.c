#ifndef _MONOLITH
#include "ampm_header.h"
#endif /* _MONOLITH */

void * malloc1(size_t size)
{ 
  void *vp=NULL;
  vp = malloc(1+size);
  if (vp==NULL){ printf(" %% Warning, out of memory: malloc failed with size %lld\n",(unsigned long long int)size); exit(EXIT_FAILURE);}
  memset(vp,0,1+size);
  if (GLOBAL_malloc1_notupdate==0){ GLOBAL_n_malloc1++;}
  return vp;
}

void free1(void **vp)
{
  if ( (vp!=NULL) && ((*vp)!=NULL) ){
    if (GLOBAL_malloc1_notupdate==0){ GLOBAL_n_malloc1--;} free(*vp); *vp=NULL; 
    /* if ( (vp!=NULL) && ((*vp)!=NULL) ){ } */}
}

void malloc1_char_FNAMESIZE__(int n_l,char ***str_p_)
{
  int nl=0;
  char **str__=NULL;
  str__=NULL;
  if (str_p_!=NULL){
    if ( (*str_p_)==NULL ){ (*str_p_) = (char **) malloc1((unsigned long long int)n_l*sizeof(char *));}
    str__ = *str_p_;
    /* if (str_p_!=NULL){ } */}
  if (str__!=NULL){
    for (nl=0;nl<n_l;nl++){
      if (str__[nl]==NULL){ str__[nl] = (char *) malloc1((unsigned long long int)FNAMESIZE*sizeof(char));}
      /* for (nl=0;nl<n_l;nl++){ } */}
    /* if (str__!=NULL){ } */}
}

void malloc1_char__(int n_l,int n_s,char ***str_p_)
{
  int nl=0;
  char **str__=NULL;
  str__=NULL;
  if (str_p_!=NULL){
    if ( (*str_p_)==NULL ){ (*str_p_) = (char **) malloc1((unsigned long long int)n_l*sizeof(char *));}
    str__ = *str_p_;
    /* if (str_p_!=NULL){ } */}
  if (str__!=NULL){
    for (nl=0;nl<n_l;nl++){
      if (str__[nl]==NULL){ str__[nl] = (char *) malloc1((unsigned long long int)n_s*sizeof(char));}
      /* for (nl=0;nl<n_l;nl++){ } */}
    /* if (str__!=NULL){ } */}
}

void free1_char__(int n_l,char ***str_p_)
{
  int nl=0;
  char **str__=NULL;
  str__=NULL;
  if (str_p_!=NULL){ if (*str_p_!=NULL){ str__ = *str_p_;}}
  if (str__!=NULL){
    for (nl=0;nl<n_l;nl++){ free1(&(str__[nl]));}
    free1(&str__);
    if (str_p_!=NULL){ *str_p_ = NULL;}
    /* if (str__!=NULL){ } */}
}
