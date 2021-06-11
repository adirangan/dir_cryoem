#ifndef _MONOLITH
#include "playpark_header.h"
#endif /* _MONOLITH */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

int iquicksort_partition_index(int *i_,int stride,int *index_,int l,int r) 
{
  int pivot=0,tmpl=0;
  int i=0,j=0,tmpi=0;
  pivot = i_[stride*l];
  i = l; j = r+1;
  do{
    do{ i++;} while( i<=r && i_[stride*i] <= pivot );
    do{ j--;} while( i_[stride*j] > pivot );
    if( i >= j ) break;
    tmpl = i_[stride*i]; i_[stride*i] = i_[stride*j]; i_[stride*j] = tmpl;
    if (index_!=NULL){ tmpi = index_[i]; index_[i] = index_[j]; index_[j] = tmpi;}
  }while(1);
  tmpl = i_[stride*l]; i_[stride*l] = i_[stride*j]; i_[stride*j] = tmpl;
  if (index_!=NULL){ tmpi = index_[l]; index_[l] = index_[j]; index_[j] = tmpi;}
  return j;
}

unsigned int iquicksort_index(unsigned int recursion_level,int *i_,int stride,int *index_,int l,int r)
{
  /* modified from http://www.comp.dit.ie/rlawlor/Alg_DS/sorting/quicksort.c */
  int j=0; unsigned int n1=recursion_level,n2=recursion_level;
  if( l < r ) { 
    j = iquicksort_partition_index(i_,stride,index_,l,r); 
    if (recursion_level<GLOBAL_recursion_limit){ 
      n1 = iquicksort_index(recursion_level+1,i_,stride,index_,l,j-1); 
      n2 = iquicksort_index(recursion_level+1,i_,stride,index_,j+1,r); 
      /* if (recursion_level<GLOBAL_recursion_limit){ } */}
    else /* (recursion_level>=GLOBAL_recursion_limit) */{ 
      printf(" %% Warning! recursion limit %d reached in iquicksort_index.\n",recursion_level);
      /* recursion_limit breached */}
    /* if( l < r ) { } */}
  return maximum(n1,n2);
}

void iquicksort_index_driver(int n_i,int *i_,int stride,int *i_workspace_,int *index_)
{
  /* finds index listing index_ so that i_[stride*index_[.]] is sorted. */
  int ni=0;
  for (ni=0;ni<n_i;ni++){ i_workspace_[ni] = i_[stride*ni]; index_[ni]=ni;}
  iquicksort_index(0,i_workspace_,1,index_,0,n_i-1);
}

void iquicksort_index_driver_test()
{
  int n_i = 5;
  int index_[n_i];
  int i_[2*n_i];
  int i_workspace_[1*n_i];
  memset(i_,0,2*n_i*sizeof(int));
  i_[0] = 15; i_[2] = 25; i_[4] = 13; i_[6] = 18; i_[8] = 20;
  array_printf(i_,"int",1,2*n_i," %% i_: ");
  iquicksort_index_driver(n_i,i_,2,i_workspace_,index_);
  array_printf(i_workspace_,"int",1,n_i," %% i_workspace_: ");
  array_printf(index_,"int",1,n_i," %% index_: ");
}

void iquicksort_index_index_driver(int n_i,int *i_,int stride,int *i_workspace_,int *index_orig_from_sort_,int *index_sort_from_orig_,int *index_workspace_)
{
  /* finds index listing index_orig_from_sort_ so that i_[stride*index_orig_from_sort_[.]] is sorted. */
  int ni=0;
  for (ni=0;ni<n_i;ni++){ i_workspace_[ni] = i_[stride*ni]; if (index_orig_from_sort_!=NULL){ index_orig_from_sort_[ni]=ni;}}
  iquicksort_index(0,i_workspace_,1,index_orig_from_sort_,0,n_i-1);
  if (index_sort_from_orig_!=NULL){
    for (ni=0;ni<n_i;ni++){ index_workspace_[ni] = index_orig_from_sort_[ni]; index_sort_from_orig_[ni]=ni;}
    iquicksort_index(0,index_workspace_,1,index_sort_from_orig_,0,n_i-1);
    /* if (index_sort_from_orig_!=NULL){ } */}
}

void iquicksort_index_index_driver_test()
{
  int n_i = 5;
  int index_orig_from_sort_[n_i];
  int index_sort_from_orig_[n_i];
  int i_[2*n_i];
  int i_workspace_[1*n_i];
  int index_workspace_[1*n_i];
  memset(i_,0,2*n_i*sizeof(int));
  i_[0] = 15; i_[2] = 25; i_[4] = 13; i_[6] = 18; i_[8] = 20;
  iquicksort_index_index_driver(n_i,i_,2,i_workspace_,index_orig_from_sort_,index_sort_from_orig_,index_workspace_);
  array_printf(i_,"int",1,2*n_i," %% i_: ");
  array_printf(i_workspace_,"int",1,n_i," %% i_workspace_: ");
  array_printf(index_orig_from_sort_,"int",1,n_i," %% index_orig_from_sort_: ");
  array_printf(index_sort_from_orig_,"int",1,n_i," %% index_sort_from_orig_: ");
}

void irandperm(int n_i,int **index_p_,unsigned long long int *rseed)
{
  double *d_=NULL;
  int nd=0;
  double *d_workspace_=NULL;
  int *index_;
  d_ = (double *) malloc1(n_i*sizeof(double));
  d_workspace_ = (double *) malloc1(n_i*sizeof(double));
  for (nd=0;nd<n_i;nd++){ d_[nd] = R01GET(rseed);}
  if (*index_p_==NULL){ (*index_p_) = (int *) malloc1(n_i*sizeof(int));}
  index_ = *index_p_;
  dquicksort_index_driver(n_i,d_,1,d_workspace_,index_);
  free1(&d_workspace_);
  free1(&d_);
}

void irandperm_test()
{
  int n_i = 10;
  int *index_ = NULL;
  unsigned long long int rseed = 67892;
  GLOBAL_tic(0);
  irandperm(n_i,&index_,&rseed); array_printf(index_,"int",1,n_i," irandperm index_: ");
  irandperm(n_i,&index_,&rseed); array_printf(index_,"int",1,n_i," irandperm index_: ");
  irandperm(n_i,&index_,&rseed); array_printf(index_,"int",1,n_i," irandperm index_: ");
  irandperm(n_i,&index_,&rseed); array_printf(index_,"int",1,n_i," irandperm index_: ");
  irandperm(n_i,&index_,&rseed); array_printf(index_,"int",1,n_i," irandperm index_: ");
  irandperm(n_i,&index_,&rseed); array_printf(index_,"int",1,n_i," irandperm index_: ");
  irandperm(n_i,&index_,&rseed); array_printf(index_,"int",1,n_i," irandperm index_: ");
  free1(&index_);
  GLOBAL_toc(0,1," irandperm_test: ");
}


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

int fquicksort_partition_index(float *f_,int stride,int *index_,int l,int r) 
{
  float pivot=0,tmpd=0;
  int i=0,j=0,tmpi=0;
  pivot = f_[stride*l];
  i = l; j = r+1;
  do{
    do{ i++;} while( i<=r && f_[stride*i] <= pivot );
    do{ j--;} while( f_[stride*j] > pivot );
    if( i >= j ) break;
    tmpd = f_[stride*i]; f_[stride*i] = f_[stride*j]; f_[stride*j] = tmpd;
    if (index_!=NULL){ tmpi = index_[i]; index_[i] = index_[j]; index_[j] = tmpi;}
  }while(1);
  tmpd = f_[stride*l]; f_[stride*l] = f_[stride*j]; f_[stride*j] = tmpd;
  if (index_!=NULL){ tmpi = index_[l]; index_[l] = index_[j]; index_[j] = tmpi;}
  return j;
}

unsigned int fquicksort_index(unsigned int recursion_level,float *f_,int stride,int *index_,int l,int r)
{
  /* modified from http://www.comp.dit.ie/rlawlor/Alg_DS/sorting/quicksort.c */
  int j=0; unsigned int n1=recursion_level,n2=recursion_level;
  if( l < r ) { 
    j = fquicksort_partition_index(f_,stride,index_,l,r); 
    if (recursion_level<GLOBAL_recursion_limit){ 
      n1 = fquicksort_index(recursion_level+1,f_,stride,index_,l,j-1); 
      n2 = fquicksort_index(recursion_level+1,f_,stride,index_,j+1,r); 
      /* if (recursion_level<GLOBAL_recursion_limit){ } */}
    else /* (recursion_level>=GLOBAL_recursion_limit) */{ 
      printf(" %% Warning! recursion limit %d reached in fquicksort_index.\n %% This means that many row- and col-scores are identical.\n %% It is likely that your data is mostly zero or constant.\n",recursion_level);
      /* recursion_limit breached */}
    /* if( l < r ) { } */}
  return maximum(n1,n2);
}

void fquicksort_index_driver(int n_f,float *f_,int stride,float *f_workspace_,int *index_)
{
  /* finds index listing index_ so that f_[stride*index_[.]] is sorted. */
  int nf=0;
  for (nf=0;nf<n_f;nf++){ f_workspace_[nf] = f_[stride*nf]; index_[nf]=nf;}
  fquicksort_index(0,f_workspace_,1,index_,0,n_f-1);
}

void fquicksort_index_driver_test()
{
  int n_f = 5;
  int index_[n_f];
  float f_[2*n_f];
  float f_workspace_[1*n_f];
  memset(f_,0,2*n_f*sizeof(float));
  f_[0] = 15; f_[2] = 25; f_[4] = 13; f_[6] = 18; f_[8] = 20;
  array_printf(f_,"float",1,2*n_f," %% f_: ");
  fquicksort_index_driver(n_f,f_,2,f_workspace_,index_);
  array_printf(f_workspace_,"float",1,n_f," %% f_workspace_: ");
  array_printf(index_,"int",1,n_f," %% index_: ");
}

void fquicksort_index_index_driver(int n_f,float *f_,int stride,float *f_workspace_,int *index_orig_from_sort_,int *index_sort_from_orig_,int *index_workspace_)
{
  /* finds index listing index_orig_from_sort_ so that f_[stride*index_orig_from_sort_[.]] is sorted. */
  int nf=0;
  for (nf=0;nf<n_f;nf++){ f_workspace_[nf] = f_[stride*nf]; if (index_orig_from_sort_!=NULL){ index_orig_from_sort_[nf]=nf;}}
  fquicksort_index(0,f_workspace_,1,index_orig_from_sort_,0,n_f-1);
  if (index_sort_from_orig_!=NULL){
    for (nf=0;nf<n_f;nf++){ index_workspace_[nf] = index_orig_from_sort_[nf]; index_sort_from_orig_[nf]=nf;}
    iquicksort_index(0,index_workspace_,1,index_sort_from_orig_,0,n_f-1);
    /* if (index_sort_from_orig_!=NULL){ } */}
}

void fquicksort_index_index_driver_test()
{
  int n_f = 5;
  int index_orig_from_sort_[n_f];
  int index_sort_from_orig_[n_f];
  float f_[2*n_f];
  float f_workspace_[1*n_f];
  int index_workspace_[1*n_f];
  memset(f_,0,2*n_f*sizeof(float));
  f_[0] = 15; f_[2] = 25; f_[4] = 13; f_[6] = 18; f_[8] = 20;
  fquicksort_index_index_driver(n_f,f_,2,f_workspace_,index_orig_from_sort_,index_sort_from_orig_,index_workspace_);
  array_printf(f_,"float",1,2*n_f," %% f_: ");
  array_printf(f_workspace_,"float",1,n_f," %% f_workspace_: ");
  array_printf(index_orig_from_sort_,"int",1,n_f," %% index_orig_from_sort_: ");
  array_printf(index_sort_from_orig_,"int",1,n_f," %% index_sort_from_orig_: ");
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

int dquicksort_partition_index(double *d_,int stride,int *index_,int l,int r) 
{
  double pivot=0,tmpd=0;
  int i=0,j=0,tmpi=0;
  pivot = d_[stride*l];
  i = l; j = r+1;
  do{
    do{ i++;} while( i<=r && d_[stride*i] <= pivot );
    do{ j--;} while( d_[stride*j] > pivot );
    if( i >= j ) break;
    tmpd = d_[stride*i]; d_[stride*i] = d_[stride*j]; d_[stride*j] = tmpd;
    if (index_!=NULL){ tmpi = index_[i]; index_[i] = index_[j]; index_[j] = tmpi;}
  }while(1);
  tmpd = d_[stride*l]; d_[stride*l] = d_[stride*j]; d_[stride*j] = tmpd;
  if (index_!=NULL){ tmpi = index_[l]; index_[l] = index_[j]; index_[j] = tmpi;}
  return j;
}

unsigned int dquicksort_index(unsigned int recursion_level,double *d_,int stride,int *index_,int l,int r)
{
  /* modified from http://www.comp.dit.ie/rlawlor/Alg_DS/sorting/quicksort.c */
  int j=0; unsigned int n1=recursion_level,n2=recursion_level;
  if( l < r ) { 
    j = dquicksort_partition_index(d_,stride,index_,l,r); 
    if (recursion_level<GLOBAL_recursion_limit){ 
      n1 = dquicksort_index(recursion_level+1,d_,stride,index_,l,j-1); 
      n2 = dquicksort_index(recursion_level+1,d_,stride,index_,j+1,r); 
      /* if (recursion_level<GLOBAL_recursion_limit){ } */}
    else /* (recursion_level>=GLOBAL_recursion_limit) */{ 
      printf(" %% Warning! recursion limit %d reached in dquicksort_index.\n %% This means that many row- and col-scores are identical.\n %% It is likely that your data is mostly zero or constant.\n",recursion_level);
      /* recursion_limit breached */}
    /* if( l < r ) { } */}
  return maximum(n1,n2);
}

void dquicksort_index_driver(int n_d,double *d_,int stride,double *d_workspace_,int *index_)
{
  /* finds index listing index_ so that d_[stride*index_[.]] is sorted. */
  int nd=0;
  for (nd=0;nd<n_d;nd++){ d_workspace_[nd] = d_[stride*nd]; index_[nd]=nd;}
  dquicksort_index(0,d_workspace_,1,index_,0,n_d-1);
}

void dquicksort_index_driver_test()
{
  int n_d = 5;
  int index_[n_d];
  double d_[2*n_d];
  double d_workspace_[1*n_d];
  memset(d_,0,2*n_d*sizeof(double));
  d_[0] = 15; d_[2] = 25; d_[4] = 13; d_[6] = 18; d_[8] = 20;
  array_printf(d_,"double",1,2*n_d," %% d_: ");
  dquicksort_index_driver(n_d,d_,2,d_workspace_,index_);
  array_printf(d_workspace_,"double",1,n_d," %% d_workspace_: ");
  array_printf(index_,"int",1,n_d," %% index_: ");
}

void dquicksort_index_index_driver(int n_d,double *d_,int stride,double *d_workspace_,int *index_orig_from_sort_,int *index_sort_from_orig_,int *index_workspace_)
{
  /* finds index listing index_orig_from_sort_ so that d_[stride*index_orig_from_sort_[.]] is sorted. */
  int nd=0;
  for (nd=0;nd<n_d;nd++){ d_workspace_[nd] = d_[stride*nd]; if (index_orig_from_sort_!=NULL){ index_orig_from_sort_[nd]=nd;}}
  dquicksort_index(0,d_workspace_,1,index_orig_from_sort_,0,n_d-1);
  if (index_sort_from_orig_!=NULL){
    for (nd=0;nd<n_d;nd++){ index_workspace_[nd] = index_orig_from_sort_[nd]; index_sort_from_orig_[nd]=nd;}
    iquicksort_index(0,index_workspace_,1,index_sort_from_orig_,0,n_d-1);
    /* if (index_sort_from_orig_!=NULL){ } */}
}

void dquicksort_index_index_driver_test()
{
  int n_d = 5;
  int index_orig_from_sort_[n_d];
  int index_sort_from_orig_[n_d];
  double d_[2*n_d];
  double d_workspace_[1*n_d];
  int index_workspace_[1*n_d];
  memset(d_,0,2*n_d*sizeof(double));
  d_[0] = 15; d_[2] = 25; d_[4] = 13; d_[6] = 18; d_[8] = 20;
  dquicksort_index_index_driver(n_d,d_,2,d_workspace_,index_orig_from_sort_,index_sort_from_orig_,index_workspace_);
  array_printf(d_,"double",1,2*n_d," %% d_: ");
  array_printf(d_workspace_,"double",1,n_d," %% d_workspace_: ");
  array_printf(index_orig_from_sort_,"int",1,n_d," %% index_orig_from_sort_: ");
  array_printf(index_sort_from_orig_,"int",1,n_d," %% index_sort_from_orig_: ");
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

