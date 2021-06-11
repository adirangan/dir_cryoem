#ifndef _MONOLITH
#include "playpark_header.h"
#endif /* _MONOLITH */

void MDA_write_i4(int n_dim,int *dim_,int *i4_,char *fname)
{
  int verbose=0;
  FILE *fp=NULL;
  size_t n_i=1,s=0;
  int ndim=0;
  if (verbose){ printf(" %% writing %s\n",fname);}
  if ((fp=fopen(fname,"w"))==NULL){ printf(" %% Warning! could not open %s.\n",fname); exit(RET_READ_FAIL);}
  s=fwrite(&n_dim,sizeof(int),1,fp); if (s!=1){ printf(" %% Warning, could not write to %s.\n",fname); exit(RET_READ_FAIL);}
  s=fwrite(dim_,sizeof(int),n_dim,fp); if (s!=n_dim){ printf(" %% Warning, could not write to %s.\n",fname); exit(RET_READ_FAIL);}
  for (ndim=0;ndim<n_dim;ndim++){ n_i*=dim_[ndim];}
  s=fwrite(i4_,sizeof(int),n_i,fp); if (s!=n_i){ printf(" %% Warning, could not write to %s.\n",fname); exit(RET_READ_FAIL);}
  fclose(fp);fp=NULL;
}

void MDA_read_i4(int *n_dim_p,int **dim_p_,int **i4_p_,char *fname)
{
  int verbose=0;
  FILE *fp=NULL;
  size_t n_i=1,s=0;
  int ndim=0;
  int n_dim=0;
  int *dim_=NULL;
  int *i4_=NULL;
  if (verbose){ printf(" %% reading %s\n",fname);}
  if ((fp=fopen(fname,"r"))==NULL){ printf(" %% Warning! could not open %s.\n",fname); exit(RET_READ_FAIL);}
  s=fread(&n_dim,sizeof(int),1,fp); if (s!=1){ printf(" %% Warning, could not read from%s.\n",fname); exit(RET_READ_FAIL);}
  if (n_dim_p!=NULL){ *n_dim_p = n_dim;}
  dim_=NULL;
  if (dim_p_!=NULL){
    if ( (*dim_p_)==NULL ){ (*dim_p_) = (int *) malloc1(n_dim*sizeof(int)); }
    dim_ = *dim_p_;
    /* if (dim_p_!=NULL){ } */}
  if (dim_!=NULL){
    s=fread(dim_,sizeof(int),n_dim,fp);
    if (s!=n_dim){ printf(" %% Warning, could not read from%s.\n",fname); exit(RET_READ_FAIL);}
    n_i = 1; for (ndim=0;ndim<n_dim;ndim++){ n_i *= (size_t)(dim_[ndim]);}
    i4_=NULL;
    if (i4_p_!=NULL){
      if ( (*i4_p_)==NULL ){ (*i4_p_) = (int *) malloc1(n_i*sizeof(int)); }
      i4_ = *i4_p_;
      /* if (i4_p_!=NULL){ } */}
    if (i4_!=NULL){
      s=fread(i4_,sizeof(int),n_i,fp);
      if (s!=n_i){ printf(" %% Warning, could not read from%s.\n",fname); exit(RET_READ_FAIL);}      
      /* if (i4_!=NULL){ } */}
    /* if (dim_!=NULL){ } */}
  fclose(fp);fp=NULL;
}

void MDA_printf_r4_margin(char *fname)
{
  int ndim=0,n_dim;
  int *dim_=NULL;
  float *r4_=NULL;
  size_t n_i=0;
  char prefix[FNAMESIZE];
  MDA_read_r4(&n_dim,&dim_,&r4_,fname);
  n_i=1; for (ndim=0;ndim<n_dim;ndim++){ n_i *= (size_t)dim_[ndim];}
  sprintf(prefix," %% %s: ",fname);
  if ( (n_dim==2) && (dim_[0]>1) && (dim_[1]>1) ){ array_printf_margin(r4_,"float",dim_[0],dim_[1],prefix);}
  else{ array_printf_margin(r4_,"float",1,n_i,prefix);}
  free1(&dim_);
  free1(&r4_);
}

void MDA_write_r4(int n_dim,int *dim_,float *r4_,char *fname)
{
  int verbose=0;
  FILE *fp=NULL;
  size_t n_i=1,s=0;
  int ndim=0;
  if (verbose){ printf(" %% writing %s\n",fname);}
  if ((fp=fopen(fname,"w"))==NULL){ printf(" %% Warning! could not open %s.\n",fname); exit(RET_READ_FAIL);}
  s=fwrite(&n_dim,sizeof(int),1,fp); if (s!=1){ printf(" %% Warning, could not write to %s.\n",fname); exit(RET_READ_FAIL);}
  s=fwrite(dim_,sizeof(int),n_dim,fp); if (s!=n_dim){ printf(" %% Warning, could not write to %s.\n",fname); exit(RET_READ_FAIL);}
  for (ndim=0;ndim<n_dim;ndim++){ n_i*=dim_[ndim];}
  s=fwrite(r4_,sizeof(float),n_i,fp); if (s!=n_i){ printf(" %% Warning, could not write to %s.\n",fname); exit(RET_READ_FAIL);}
  fclose(fp);fp=NULL;
}

void MDA_read_r4(int *n_dim_p,int **dim_p_,float **r4_p_,char *fname)
{
  int verbose=0;
  FILE *fp=NULL;
  size_t n_i=1,s=0;
  int ndim=0;
  int n_dim=0;
  int *dim_=NULL;
  float *r4_=NULL;
  if (verbose){ printf(" %% reading %s\n",fname);}
  if ((fp=fopen(fname,"r"))==NULL){ printf(" %% Warning! could not open %s.\n",fname); exit(RET_READ_FAIL);}
  s=fread(&n_dim,sizeof(int),1,fp); if (s!=1){ printf(" %% Warning, could not read from%s.\n",fname); exit(RET_READ_FAIL);}
  if (n_dim_p!=NULL){ *n_dim_p = n_dim;}
  dim_=NULL;
  if (dim_p_!=NULL){
    if ( (*dim_p_)==NULL ){ (*dim_p_) = (int *) malloc1(n_dim*sizeof(int)); }
    dim_ = *dim_p_;
    /* if (dim_p_!=NULL){ } */}
  if (dim_!=NULL){
    s=fread(dim_,sizeof(int),n_dim,fp);
    if (s!=n_dim){ printf(" %% Warning, could not read from%s.\n",fname); exit(RET_READ_FAIL);}
    n_i = 1; for (ndim=0;ndim<n_dim;ndim++){ n_i *= (size_t)(dim_[ndim]);}
    r4_=NULL;
    if (r4_p_!=NULL){
      if ( (*r4_p_)==NULL ){ (*r4_p_) = (float *) malloc1(n_i*sizeof(float)); }
      r4_ = *r4_p_;
      /* if (r4_p_!=NULL){ } */}
    if (r4_!=NULL){
      s=fread(r4_,sizeof(float),n_i,fp);
      if (s!=n_i){ printf(" %% Warning, could not read from%s.\n",fname); exit(RET_READ_FAIL);}      
      /* if (r4_!=NULL){ } */}
    /* if (dim_!=NULL){ } */}
  fclose(fp);fp=NULL;
}

void MDA_printf_i4_margin(char *fname)
{
  int ndim=0,n_dim;
  int *dim_=NULL;
  int *i4_=NULL;
  size_t n_i=0;
  char prefix[FNAMESIZE];
  MDA_read_i4(&n_dim,&dim_,&i4_,fname);
  n_i=1; for (ndim=0;ndim<n_dim;ndim++){ n_i *= (size_t)dim_[ndim];}
  sprintf(prefix," %% %s: ",fname);
  if ( (n_dim==2) && (dim_[0]>1) && (dim_[1]>1) ){ array_printf_margin(i4_,"int",dim_[0],dim_[1],prefix);}
  else{ array_printf_margin(i4_,"int",1,n_i,prefix);}
  free1(&dim_);
  free1(&i4_);
}

void MDA_write_r8(int n_dim,int *dim_,double *r8_,char *fname)
{
  int verbose=0;
  FILE *fp=NULL;
  size_t n_i=1,s=0;
  int ndim=0;
  if (verbose){ printf(" %% writing %s\n",fname);}
  if ((fp=fopen(fname,"w"))==NULL){ printf(" %% Warning! could not open %s.\n",fname); exit(RET_READ_FAIL);}
  s=fwrite(&n_dim,sizeof(int),1,fp); if (s!=1){ printf(" %% Warning, could not write to %s.\n",fname); exit(RET_READ_FAIL);}
  s=fwrite(dim_,sizeof(int),n_dim,fp); if (s!=n_dim){ printf(" %% Warning, could not write to %s.\n",fname); exit(RET_READ_FAIL);}
  for (ndim=0;ndim<n_dim;ndim++){ n_i*=dim_[ndim];}
  s=fwrite(r8_,sizeof(double),n_i,fp); if (s!=n_i){ printf(" %% Warning, could not write to %s.\n",fname); exit(RET_READ_FAIL);}
  fclose(fp);fp=NULL;
}

void MDA_read_r8(int *n_dim_p,int **dim_p_,double **r8_p_,char *fname)
{
  int verbose=0;
  FILE *fp=NULL;
  size_t n_i=1,s=0;
  int ndim=0;
  int n_dim=0;
  int *dim_=NULL;
  double *r8_=NULL;
  if (verbose){ printf(" %% reading %s\n",fname);}
  if ((fp=fopen(fname,"r"))==NULL){ printf(" %% Warning! could not open %s.\n",fname); exit(RET_READ_FAIL);}
  s=fread(&n_dim,sizeof(int),1,fp); if (s!=1){ printf(" %% Warning, could not read from%s.\n",fname); exit(RET_READ_FAIL);}
  if (n_dim_p!=NULL){ *n_dim_p = n_dim;}
  dim_=NULL;
  if (dim_p_!=NULL){
    if ( (*dim_p_)==NULL ){ (*dim_p_) = (int *) malloc1(n_dim*sizeof(int)); }
    dim_ = *dim_p_;
    /* if (dim_p_!=NULL){ } */}
  if (dim_!=NULL){
    s=fread(dim_,sizeof(int),n_dim,fp);
    if (s!=n_dim){ printf(" %% Warning, could not read from%s.\n",fname); exit(RET_READ_FAIL);}
    n_i = 1; for (ndim=0;ndim<n_dim;ndim++){ n_i *= (size_t)(dim_[ndim]);}
    r8_=NULL;
    if (r8_p_!=NULL){
      if ( (*r8_p_)==NULL ){ (*r8_p_) = (double *) malloc1(n_i*sizeof(double)); }
      r8_ = *r8_p_;
      /* if (r8_p_!=NULL){ } */}
    if (r8_!=NULL){
      s=fread(r8_,sizeof(double),n_i,fp);
      if (s!=n_i){ printf(" %% Warning, could not read from%s.\n",fname); exit(RET_READ_FAIL);}      
      /* if (r8_!=NULL){ } */}
    /* if (dim_!=NULL){ } */}
  fclose(fp);fp=NULL;
}

void MDA_printf_r8_margin(char *fname)
{
  int ndim=0,n_dim;
  int *dim_=NULL;
  double *r8_=NULL;
  size_t n_i=0;
  char prefix[FNAMESIZE];
  MDA_read_r8(&n_dim,&dim_,&r8_,fname);
  n_i=1; for (ndim=0;ndim<n_dim;ndim++){ n_i *= (size_t)dim_[ndim];}
  sprintf(prefix," %% %s: ",fname);
  if ( (n_dim==2) && (dim_[0]>1) && (dim_[1]>1) ){ array_printf_margin(r8_,"double",dim_[0],dim_[1],prefix);}
  else{ array_printf_margin(r8_,"double",1,n_i,prefix);}
  free1(&dim_);
  free1(&r8_);
}

void MDA_write_ulli(int n_dim,int *dim_,unsigned long long int *ulli_,char *fname)
{
  int verbose=0;
  FILE *fp=NULL;
  size_t n_i=1,s=0;
  int ndim=0;
  if (verbose){ printf(" %% writing %s\n",fname);}
  if ((fp=fopen(fname,"w"))==NULL){ printf(" %% Warning! could not open %s.\n",fname); exit(RET_READ_FAIL);}
  s=fwrite(&n_dim,sizeof(int),1,fp); if (s!=1){ printf(" %% Warning, could not write to %s.\n",fname); exit(RET_READ_FAIL);}
  s=fwrite(dim_,sizeof(int),n_dim,fp); if (s!=n_dim){ printf(" %% Warning, could not write to %s.\n",fname); exit(RET_READ_FAIL);}
  for (ndim=0;ndim<n_dim;ndim++){ n_i*=dim_[ndim];}
  s=fwrite(ulli_,sizeof(unsigned long long int),n_i,fp); if (s!=n_i){ printf(" %% Warning, could not write to %s.\n",fname); exit(RET_READ_FAIL);}
  fclose(fp);fp=NULL;
}

void MDA_read_ulli(int *n_dim_p,int **dim_p_,unsigned long long int **ulli_p_,char *fname)
{
  int verbose=0;
  FILE *fp=NULL;
  size_t n_i=1,s=0;
  int ndim=0;
  int n_dim=0;
  int *dim_=NULL;
  unsigned long long int *ulli_=NULL;
  if (verbose){ printf(" %% reading %s\n",fname);}
  if ((fp=fopen(fname,"r"))==NULL){ printf(" %% Warning! could not open %s.\n",fname); exit(RET_READ_FAIL);}
  s=fread(&n_dim,sizeof(int),1,fp); if (s!=1){ printf(" %% Warning, could not read from%s.\n",fname); exit(RET_READ_FAIL);}
  if (n_dim_p!=NULL){ *n_dim_p = n_dim;}
  dim_=NULL;
  if (dim_p_!=NULL){
    if ( (*dim_p_)==NULL ){ (*dim_p_) = (int *) malloc1(n_dim*sizeof(int)); }
    dim_ = *dim_p_;
    /* if (dim_p_!=NULL){ } */}
  if (dim_!=NULL){
    s=fread(dim_,sizeof(int),n_dim,fp);
    if (s!=n_dim){ printf(" %% Warning, could not read from%s.\n",fname); exit(RET_READ_FAIL);}
    n_i = 1; for (ndim=0;ndim<n_dim;ndim++){ n_i *= (size_t)(dim_[ndim]);}
    ulli_=NULL;
    if (ulli_p_!=NULL){
      if ( (*ulli_p_)==NULL ){ (*ulli_p_) = (unsigned long long int *) malloc1(n_i*sizeof(unsigned long long int)); }
      ulli_ = *ulli_p_;
      /* if (ulli_p_!=NULL){ } */}
    if (ulli_!=NULL){
      s=fread(ulli_,sizeof(unsigned long long int),n_i,fp);
      if (s!=n_i){ printf(" %% Warning, could not read from%s.\n",fname); exit(RET_READ_FAIL);}      
      /* if (ulli_!=NULL){ } */}
    /* if (dim_!=NULL){ } */}
  fclose(fp);fp=NULL;
}

void MDA_printf_ulli_margin(char *fname)
{
  int ndim=0,n_dim;
  int *dim_=NULL;
  unsigned long long int *ulli_=NULL;
  size_t n_i=0;
  char prefix[FNAMESIZE];
  MDA_read_ulli(&n_dim,&dim_,&ulli_,fname);
  n_i=1; for (ndim=0;ndim<n_dim;ndim++){ n_i *= (size_t)dim_[ndim];}
  sprintf(prefix," %% %s: ",fname);
  if ( (n_dim==2) && (dim_[0]>1) && (dim_[1]>1) ){ array_printf_margin(ulli_,"unsigned long long int",dim_[0],dim_[1],prefix);}
  else{ array_printf_margin(ulli_,"unsigned long long int",1,n_i,prefix);}
  free1(&dim_);
  free1(&ulli_);
}

void MDA_io_test()
{
  int verbose=1;
  int n_dim_0 = 2;
  int *dim_0_=NULL;
  int n_i=0;
  int ndim=0;
  int *i4_0_=NULL;
  int *i4_1_=NULL;
  int *dim_1_=NULL;
  int n_dim_1=NULL;
  char fname[32];
  double *r8_0_;
  int n_dim_2=NULL;
  int *r8_2_=NULL;
  int *dim_2_=NULL;
  float *r4_0_;
  int n_dim_3=NULL;
  int *r4_3_=NULL;
  int *dim_3_=NULL;
  unsigned long long int *ulli_0_;
  int n_dim_4=NULL;
  int *ulli_4_=NULL;
  int *dim_4_=NULL;
  GLOBAL_tic(0);
  dim_0_ = (int *) malloc1(n_dim_0*sizeof(int));
  dim_0_[0] = 2; dim_0_[1] = 3;
  n_i = 1; for (ndim=0;ndim<n_dim_0;ndim++){ n_i*=dim_0_[ndim];}
  /* %%%%%%%% */
  i4_0_ = (int *) malloc1(n_i*sizeof(int));
  i4_0_[0] = 0; i4_0_[2] = 10; i4_0_[4] = 100;
  i4_0_[1] = 1; i4_0_[3] = 20; i4_0_[5] = 200;
  sprintf(fname,"MDA_io.test");
  array_printf(i4_0_,"int",dim_0_[0],dim_0_[1]," %% i4_0_: ");
  MDA_write_i4(n_dim_0,dim_0_,i4_0_,fname);
  MDA_read_i4(&n_dim_1,&dim_1_,&i4_1_,fname);
  array_printf(i4_1_,"int",dim_1_[0],dim_1_[1]," %% i4_1_: ");
  /* %%%%%%%% */
  r8_0_ = (double *) malloc1(n_i*sizeof(double));
  r8_0_[0] = 0.1; r8_0_[2] = 10.1; r8_0_[4] = 100.1;
  r8_0_[1] = 1.1; r8_0_[3] = 20.1; r8_0_[5] = 200.1;
  array_printf(r8_0_,"double",dim_0_[0],dim_0_[1]," %% r8_0_: ");
  MDA_write_r8(n_dim_0,dim_0_,r8_0_,fname);
  MDA_read_r8(&n_dim_2,&dim_2_,&r8_2_,fname);
  array_printf(r8_2_,"double",dim_2_[0],dim_2_[1]," %% r8_2_: ");
  /* %%%%%%%% */
  r4_0_ = (float *) malloc1(n_i*sizeof(float));
  r4_0_[0] = 0.1; r4_0_[2] = 10.1; r4_0_[4] = 100.1;
  r4_0_[1] = 1.1; r4_0_[3] = 20.1; r4_0_[5] = 200.1;
  array_printf(r4_0_,"float",dim_0_[0],dim_0_[1]," %% r4_0_: ");
  MDA_write_r4(n_dim_0,dim_0_,r4_0_,fname);
  MDA_read_r4(&n_dim_3,&dim_3_,&r4_3_,fname);
  array_printf(r4_3_,"float",dim_3_[0],dim_3_[1]," %% r4_3_: ");
  /* %%%%%%%% */
  ulli_0_ = (unsigned long long int *) malloc1(n_i*sizeof(unsigned long long int));
  ulli_0_[0] = 100l; ulli_0_[2] = 1010l; ulli_0_[4] = 10100l;
  ulli_0_[1] = 101l; ulli_0_[3] = 1020l; ulli_0_[5] = 10200l;
  array_printf(ulli_0_,"unsigned long long int",dim_0_[0],dim_0_[1]," %% ulli_0_: ");
  MDA_write_ulli(n_dim_0,dim_0_,ulli_0_,fname);
  MDA_read_ulli(&n_dim_4,&dim_4_,&ulli_4_,fname);
  array_printf(ulli_4_,"unsigned long long int",dim_4_[0],dim_4_[1]," %% ulli_4_: ");
  /* %%%%%%%% */
  free1(&dim_0_);
  free1(&i4_0_);
  free1(&r8_0_);
  free1(&r4_0_);
  free1(&ulli_0_);
  free1(&ulli_4_);free1(&dim_4_);
  free1(&r4_3_);free1(&dim_3_);
  free1(&r8_2_);free1(&dim_2_);
  free1(&i4_1_);free1(&dim_1_);
  GLOBAL_toc(0,1," % MDA_io_test: ");
}
