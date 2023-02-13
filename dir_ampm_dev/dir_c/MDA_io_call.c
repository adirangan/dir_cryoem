#ifndef _MONOLITH
#include "ampm_header.h"
#endif /* _MONOLITH */

void MDA_write_i4(int n_dim,int *dim_,int *i4_,const char *fname)
{
  int verbose=0;
  FILE *fp=NULL;
  size_t n_i=1,s=0;
  int ndim=0;
  if (verbose){ printf(" %% writing %s\n",fname);}
  if ((fp=fopen(fname,"w"))==NULL){ printf(" %% Warning! could not open %s.\n",fname); exit(EXIT_FAILURE);}
  s=fwrite(&n_dim,sizeof(int),1,fp); if (s!=1){ printf(" %% Warning, could not write to %s.\n",fname); exit(EXIT_FAILURE);}
  s=fwrite(dim_,sizeof(int),n_dim,fp); if (s!=n_dim){ printf(" %% Warning, could not write to %s.\n",fname); exit(EXIT_FAILURE);}
  for (ndim=0;ndim<n_dim;ndim++){ n_i*=dim_[ndim];}
  s=fwrite(i4_,sizeof(int),n_i,fp); if (s!=n_i){ printf(" %% Warning, could not write to %s.\n",fname); exit(EXIT_FAILURE);}
  fclose(fp);fp=NULL;
}

void MDA_read_i4(int *n_dim_p,int **dim_p_,int **i4_p_,const char *fname)
{
  int verbose=0;
  FILE *fp=NULL;
  size_t n_i=1,s=0;
  int ndim=0;
  int n_dim=0;
  int *dim_=NULL;
  int *i4_=NULL;
  if (verbose){ printf(" %% reading %s\n",fname);}
  if ((fp=fopen(fname,"r"))==NULL){ printf(" %% Warning! could not open %s.\n",fname); exit(EXIT_FAILURE);}
  s=fread(&n_dim,sizeof(int),1,fp); if (s!=1){ printf(" %% Warning, could not read from%s.\n",fname); exit(EXIT_FAILURE);}
  if (n_dim_p!=NULL){ *n_dim_p = n_dim;}
  dim_=NULL;
  if (dim_p_!=NULL){
    if ( (*dim_p_)==NULL ){ (*dim_p_) = (int *) malloc1(n_dim*sizeof(int)); }
    dim_ = *dim_p_;
    /* if (dim_p_!=NULL){ } */}
  if (dim_!=NULL){
    s=fread(dim_,sizeof(int),n_dim,fp);
    if (s!=n_dim){ printf(" %% Warning, could not read from%s.\n",fname); exit(EXIT_FAILURE);}
    n_i = 1; for (ndim=0;ndim<n_dim;ndim++){ n_i *= (size_t)(dim_[ndim]);}
    i4_=NULL;
    if (i4_p_!=NULL){
      if ( (*i4_p_)==NULL ){ (*i4_p_) = (int *) malloc1(n_i*sizeof(int)); }
      i4_ = *i4_p_;
      /* if (i4_p_!=NULL){ } */}
    if (i4_!=NULL){
      s=fread(i4_,sizeof(int),n_i,fp);
      if (s!=n_i){ printf(" %% Warning, could not read from%s.\n",fname); exit(EXIT_FAILURE);}      
      /* if (i4_!=NULL){ } */}
    /* if (dim_!=NULL){ } */}
  fclose(fp);fp=NULL;
}

void MDA_write_r4(int n_dim,int *dim_,float *r4_,const char *fname)
{
  int verbose=0;
  FILE *fp=NULL;
  size_t n_i=1,s=0;
  int ndim=0;
  if (verbose){ printf(" %% writing %s\n",fname);}
  if ((fp=fopen(fname,"w"))==NULL){ printf(" %% Warning! could not open %s.\n",fname); exit(EXIT_FAILURE);}
  s=fwrite(&n_dim,sizeof(int),1,fp); if (s!=1){ printf(" %% Warning, could not write to %s.\n",fname); exit(EXIT_FAILURE);}
  s=fwrite(dim_,sizeof(int),n_dim,fp); if (s!=n_dim){ printf(" %% Warning, could not write to %s.\n",fname); exit(EXIT_FAILURE);}
  for (ndim=0;ndim<n_dim;ndim++){ n_i*=dim_[ndim];}
  s=fwrite(r4_,sizeof(float),n_i,fp); if (s!=n_i){ printf(" %% Warning, could not write to %s.\n",fname); exit(EXIT_FAILURE);}
  fclose(fp);fp=NULL;
}

void MDA_read_r4(int *n_dim_p,int **dim_p_,float **r4_p_,const char *fname)
{
  int verbose=0;
  FILE *fp=NULL;
  size_t n_i=1,s=0;
  int ndim=0;
  int n_dim=0;
  int *dim_=NULL;
  float *r4_=NULL;
  if (verbose){ printf(" %% reading %s\n",fname);}
  if ((fp=fopen(fname,"r"))==NULL){ printf(" %% Warning! could not open %s.\n",fname); exit(EXIT_FAILURE);}
  s=fread(&n_dim,sizeof(int),1,fp); if (s!=1){ printf(" %% Warning, could not read from%s.\n",fname); exit(EXIT_FAILURE);}
  if (n_dim_p!=NULL){ *n_dim_p = n_dim;}
  dim_=NULL;
  if (dim_p_!=NULL){
    if ( (*dim_p_)==NULL ){ (*dim_p_) = (int *) malloc1(n_dim*sizeof(int)); }
    dim_ = *dim_p_;
    /* if (dim_p_!=NULL){ } */}
  if (dim_!=NULL){
    s=fread(dim_,sizeof(int),n_dim,fp);
    if (s!=n_dim){ printf(" %% Warning, could not read from%s.\n",fname); exit(EXIT_FAILURE);}
    n_i = 1; for (ndim=0;ndim<n_dim;ndim++){ n_i *= (size_t)(dim_[ndim]);}
    r4_=NULL;
    if (r4_p_!=NULL){
      if ( (*r4_p_)==NULL ){ (*r4_p_) = (float *) malloc1(n_i*sizeof(float)); }
      r4_ = *r4_p_;
      /* if (r4_p_!=NULL){ } */}
    if (r4_!=NULL){
      s=fread(r4_,sizeof(float),n_i,fp);
      if (s!=n_i){ printf(" %% Warning, could not read from%s.\n",fname); exit(EXIT_FAILURE);}      
      /* if (r4_!=NULL){ } */}
    /* if (dim_!=NULL){ } */}
  fclose(fp);fp=NULL;
}

void MDA_write_r8(int n_dim,int *dim_,double *r8_,const char *fname)
{
  int verbose=0;
  FILE *fp=NULL;
  size_t n_i=1,s=0;
  int ndim=0;
  if (verbose){ printf(" %% writing %s\n",fname);}
  if ((fp=fopen(fname,"w"))==NULL){ printf(" %% Warning! could not open %s.\n",fname); exit(EXIT_FAILURE);}
  s=fwrite(&n_dim,sizeof(int),1,fp); if (s!=1){ printf(" %% Warning, could not write to %s.\n",fname); exit(EXIT_FAILURE);}
  s=fwrite(dim_,sizeof(int),n_dim,fp); if (s!=n_dim){ printf(" %% Warning, could not write to %s.\n",fname); exit(EXIT_FAILURE);}
  for (ndim=0;ndim<n_dim;ndim++){ n_i*=dim_[ndim];}
  s=fwrite(r8_,sizeof(double),n_i,fp); if (s!=n_i){ printf(" %% Warning, could not write to %s.\n",fname); exit(EXIT_FAILURE);}
  fclose(fp);fp=NULL;
}

void MDA_read_r8(int *n_dim_p,int **dim_p_,double **r8_p_,const char *fname)
{
  int verbose=0;
  FILE *fp=NULL;
  size_t n_i=1,s=0;
  int ndim=0;
  int n_dim=0;
  int *dim_=NULL;
  double *r8_=NULL;
  if (verbose){ printf(" %% reading %s\n",fname);}
  if ((fp=fopen(fname,"r"))==NULL){ printf(" %% Warning! could not open %s.\n",fname); exit(EXIT_FAILURE);}
  s=fread(&n_dim,sizeof(int),1,fp); if (s!=1){ printf(" %% Warning, could not read from%s.\n",fname); exit(EXIT_FAILURE);}
  if (n_dim_p!=NULL){ *n_dim_p = n_dim;}
  dim_=NULL;
  if (dim_p_!=NULL){
    if ( (*dim_p_)==NULL ){ (*dim_p_) = (int *) malloc1(n_dim*sizeof(int)); }
    dim_ = *dim_p_;
    /* if (dim_p_!=NULL){ } */}
  if (dim_!=NULL){
    s=fread(dim_,sizeof(int),n_dim,fp);
    if (s!=n_dim){ printf(" %% Warning, could not read from%s.\n",fname); exit(EXIT_FAILURE);}
    n_i = 1; for (ndim=0;ndim<n_dim;ndim++){ n_i *= (size_t)(dim_[ndim]);}
    r8_=NULL;
    if (r8_p_!=NULL){
      if ( (*r8_p_)==NULL ){ (*r8_p_) = (double *) malloc1(n_i*sizeof(double)); }
      r8_ = *r8_p_;
      /* if (r8_p_!=NULL){ } */}
    if (r8_!=NULL){
      s=fread(r8_,sizeof(double),n_i,fp);
      if (s!=n_i){ printf(" %% Warning, could not read from%s.\n",fname); exit(EXIT_FAILURE);}      
      /* if (r8_!=NULL){ } */}
    /* if (dim_!=NULL){ } */}
  fclose(fp);fp=NULL;
}

void MDA_write_c16(int n_dim,int *dim_,double complex *c16_,const char *fname)
{
  int verbose=0;
  FILE *fp=NULL;
  size_t n_i=1,s=0;
  int ndim=0;
  if (verbose){ printf(" %% writing %s\n",fname);}
  if ((fp=fopen(fname,"w"))==NULL){ printf(" %% Warning! could not open %s.\n",fname); exit(EXIT_FAILURE);}
  s=fwrite(&n_dim,sizeof(int),1,fp); if (s!=1){ printf(" %% Warning, could not write to %s.\n",fname); exit(EXIT_FAILURE);}
  s=fwrite(dim_,sizeof(int),n_dim,fp); if (s!=n_dim){ printf(" %% Warning, could not write to %s.\n",fname); exit(EXIT_FAILURE);}
  for (ndim=0;ndim<n_dim;ndim++){ n_i*=dim_[ndim];}
  s=fwrite(c16_,sizeof(double complex),n_i,fp); if (s!=n_i){ printf(" %% Warning, could not write to %s.\n",fname); exit(EXIT_FAILURE);}
  fclose(fp);fp=NULL;
}

void MDA_read_c16(int *n_dim_p,int **dim_p_,double complex **c16_p_,const char *fname)
{
  int verbose=0;
  FILE *fp=NULL;
  size_t n_i=1,s=0;
  int ndim=0;
  int n_dim=0;
  int *dim_=NULL;
  double complex *c16_=NULL;
  if (verbose){ printf(" %% reading %s\n",fname);}
  if ((fp=fopen(fname,"r"))==NULL){ printf(" %% Warning! could not open %s.\n",fname); exit(EXIT_FAILURE);}
  s=fread(&n_dim,sizeof(int),1,fp); if (s!=1){ printf(" %% Warning, could not read from%s.\n",fname); exit(EXIT_FAILURE);}
  if (n_dim_p!=NULL){ *n_dim_p = n_dim;}
  dim_=NULL;
  if (dim_p_!=NULL){
    if ( (*dim_p_)==NULL ){ (*dim_p_) = (int *) malloc(n_dim*sizeof(int)); }
    dim_ = *dim_p_;
    /* if (dim_p_!=NULL){ } */}
  if (dim_!=NULL){
    s=fread(dim_,sizeof(int),n_dim,fp);
    if (s!=n_dim){ printf(" %% Warning, could not read from%s.\n",fname); exit(EXIT_FAILURE);}
    n_i = 1; for (ndim=0;ndim<n_dim;ndim++){ n_i *= (size_t)(dim_[ndim]);}
    c16_=NULL;
    if (c16_p_!=NULL){
      if ( (*c16_p_)==NULL ){ (*c16_p_) = (double complex *) malloc(n_i*sizeof(double complex)); }
      c16_ = *c16_p_;
      /* if (c16_p_!=NULL){ } */}
    if (c16_!=NULL){
      s=fread(c16_,sizeof(double complex),n_i,fp);
      if (s!=n_i){ printf(" %% Warning, could not read from%s.\n",fname); exit(EXIT_FAILURE);}      
      /* if (c16_!=NULL){ } */}
    /* if (dim_!=NULL){ } */}
  fclose(fp);fp=NULL;
}

void MDA_write_ulli(int n_dim,int *dim_,unsigned long long int *ulli_,const char *fname)
{
  int verbose=0;
  FILE *fp=NULL;
  size_t n_i=1,s=0;
  int ndim=0;
  if (verbose){ printf(" %% writing %s\n",fname);}
  if ((fp=fopen(fname,"w"))==NULL){ printf(" %% Warning! could not open %s.\n",fname); exit(EXIT_FAILURE);}
  s=fwrite(&n_dim,sizeof(int),1,fp); if (s!=1){ printf(" %% Warning, could not write to %s.\n",fname); exit(EXIT_FAILURE);}
  s=fwrite(dim_,sizeof(int),n_dim,fp); if (s!=n_dim){ printf(" %% Warning, could not write to %s.\n",fname); exit(EXIT_FAILURE);}
  for (ndim=0;ndim<n_dim;ndim++){ n_i*=dim_[ndim];}
  s=fwrite(ulli_,sizeof(unsigned long long int),n_i,fp); if (s!=n_i){ printf(" %% Warning, could not write to %s.\n",fname); exit(EXIT_FAILURE);}
  fclose(fp);fp=NULL;
}

void MDA_read_ulli(int *n_dim_p,int **dim_p_,unsigned long long int **ulli_p_,const char *fname)
{
  int verbose=0;
  FILE *fp=NULL;
  size_t n_i=1,s=0;
  int ndim=0;
  int n_dim=0;
  int *dim_=NULL;
  unsigned long long int *ulli_=NULL;
  if (verbose){ printf(" %% reading %s\n",fname);}
  if ((fp=fopen(fname,"r"))==NULL){ printf(" %% Warning! could not open %s.\n",fname); exit(EXIT_FAILURE);}
  s=fread(&n_dim,sizeof(int),1,fp); if (s!=1){ printf(" %% Warning, could not read from%s.\n",fname); exit(EXIT_FAILURE);}
  if (n_dim_p!=NULL){ *n_dim_p = n_dim;}
  dim_=NULL;
  if (dim_p_!=NULL){
    if ( (*dim_p_)==NULL ){ (*dim_p_) = (int *) malloc1(n_dim*sizeof(int)); }
    dim_ = *dim_p_;
    /* if (dim_p_!=NULL){ } */}
  if (dim_!=NULL){
    s=fread(dim_,sizeof(int),n_dim,fp);
    if (s!=n_dim){ printf(" %% Warning, could not read from%s.\n",fname); exit(EXIT_FAILURE);}
    n_i = 1; for (ndim=0;ndim<n_dim;ndim++){ n_i *= (size_t)(dim_[ndim]);}
    ulli_=NULL;
    if (ulli_p_!=NULL){
      if ( (*ulli_p_)==NULL ){ (*ulli_p_) = (unsigned long long int *) malloc1(n_i*sizeof(unsigned long long int)); }
      ulli_ = *ulli_p_;
      /* if (ulli_p_!=NULL){ } */}
    if (ulli_!=NULL){
      s=fread(ulli_,sizeof(unsigned long long int),n_i,fp);
      if (s!=n_i){ printf(" %% Warning, could not read from%s.\n",fname); exit(EXIT_FAILURE);}      
      /* if (ulli_!=NULL){ } */}
    /* if (dim_!=NULL){ } */}
  fclose(fp);fp=NULL;
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
  double *r8_2_=NULL;
  int *dim_2_=NULL;
  float *r4_0_;
  int n_dim_3=NULL;
  float *r4_3_=NULL;
  int *dim_3_=NULL;
  unsigned long long int *ulli_0_;
  int n_dim_4=NULL;
  unsigned long long int *ulli_4_=NULL;
  int *dim_4_=NULL;
  double complex *c16_0_;
  int n_dim_5=NULL;
  double complex *c16_5_=NULL;
  int *dim_5_=NULL;
  if (verbose){ printf(" %% [entering MDA_io_test]\n");}
  /* %%%%%%%% */
  dim_0_ = (int *) malloc(n_dim_0*sizeof(int));
  dim_0_[0] = 2; dim_0_[1] = 3;
  n_i = 1; for (ndim=0;ndim<n_dim_0;ndim++){ n_i*=dim_0_[ndim];}
  /* %%%%%%%% */
  i4_0_ = (int *) malloc(n_i*sizeof(int));
  i4_0_[0] = 0; i4_0_[2] = 10; i4_0_[4] = 100;
  i4_0_[1] = 1; i4_0_[3] = 20; i4_0_[5] = 200;
  sprintf(fname,"MDA_io.test");
  array_printf(i4_0_,"int",dim_0_[0],dim_0_[1]," %% i4_0_: ");
  MDA_write_i4(n_dim_0,dim_0_,i4_0_,fname);
  MDA_read_i4(&n_dim_1,&dim_1_,&i4_1_,fname);
  array_printf(i4_1_,"int",dim_1_[0],dim_1_[1]," %% i4_1_: ");
  printf(" %% i4_0_ vs i4_1_: relative error %0.16f\n",ifnormn(n_i,i4_0_,i4_1_));
  /* %%%%%%%% */
  r8_0_ = (double *) malloc(n_i*sizeof(double));
  r8_0_[0] = 0.1; r8_0_[2] = 10.1; r8_0_[4] = 100.1;
  r8_0_[1] = 1.1; r8_0_[3] = 20.1; r8_0_[5] = 200.1;
  array_printf(r8_0_,"double",dim_0_[0],dim_0_[1]," %% r8_0_: ");
  MDA_write_r8(n_dim_0,dim_0_,r8_0_,fname);
  MDA_read_r8(&n_dim_2,&dim_2_,&r8_2_,fname);
  array_printf(r8_2_,"double",dim_2_[0],dim_2_[1]," %% r8_2_: ");
  printf(" %% r8_0_ vs r8_2_: relative error %0.16f\n",dfnormn(n_i,r8_0_,r8_2_));
  /* %%%%%%%% */
  r4_0_ = (float *) malloc(n_i*sizeof(float));
  r4_0_[0] = 0.1; r4_0_[2] = 10.1; r4_0_[4] = 100.1;
  r4_0_[1] = 1.1; r4_0_[3] = 20.1; r4_0_[5] = 200.1;
  array_printf(r4_0_,"float",dim_0_[0],dim_0_[1]," %% r4_0_: ");
  MDA_write_r4(n_dim_0,dim_0_,r4_0_,fname);
  MDA_read_r4(&n_dim_3,&dim_3_,&r4_3_,fname);
  array_printf(r4_3_,"float",dim_3_[0],dim_3_[1]," %% r4_3_: ");
  printf(" %% r4_0_ vs r4_3_: relative error %0.16f\n",ffnormn(n_i,r4_0_,r4_3_));
  /* %%%%%%%% */
  ulli_0_ = (unsigned long long int *) malloc(n_i*sizeof(unsigned long long int));
  ulli_0_[0] = 100l; ulli_0_[2] = 1010l; ulli_0_[4] = 10100l;
  ulli_0_[1] = 101l; ulli_0_[3] = 1020l; ulli_0_[5] = 10200l;
  array_printf(ulli_0_,"unsigned long long int",dim_0_[0],dim_0_[1]," %% ulli_0_: ");
  MDA_write_ulli(n_dim_0,dim_0_,ulli_0_,fname);
  MDA_read_ulli(&n_dim_4,&dim_4_,&ulli_4_,fname);
  array_printf(ulli_4_,"unsigned long long int",dim_4_[0],dim_4_[1]," %% ulli_4_: ");
  printf(" %% ulli_0_ vs ulli_4_: relative error %0.16f\n",ullifnormn(n_i,ulli_0_,ulli_4_));
  /* %%%%%%%% */
  c16_0_ = (double complex *) malloc(n_i*sizeof(double complex));
  c16_0_[0] = (double complex)100.0 + (double complex)10.0*_Complex_I;
  c16_0_[2] = (double complex)1010.0 + (double complex)12.0*_Complex_I;
  c16_0_[4] = (double complex)10100.0 + (double complex)14.0*_Complex_I;
  c16_0_[1] = (double complex)101.0 + (double complex)11.0*_Complex_I;
  c16_0_[3] = (double complex)1020.0 + (double complex)13.0*_Complex_I;
  c16_0_[5] = (double complex)10200.0 + (double complex)15.0*_Complex_I;
  array_printf(c16_0_,"double complex",dim_0_[0],dim_0_[1]," %% c16_0_: ");
  MDA_write_c16(n_dim_0,dim_0_,c16_0_,fname);
  MDA_read_c16(&n_dim_5,&dim_5_,&c16_5_,fname);
  array_printf(c16_5_,"double complex",dim_5_[0],dim_5_[1]," %% c16_5_: ");
  printf(" %% c16_0_ vs c16_5_: relative error %0.16f\n",zfnormn(n_i,c16_0_,c16_5_));
  /* %%%%%%%% */
  free(dim_0_); dim_0_=NULL;
  free(i4_0_); i4_0_=NULL;
  free(r8_0_); r8_0_=NULL;
  free(r4_0_); r4_0_=NULL;
  free(ulli_0_); ulli_0_=NULL;
  free(ulli_4_); ulli_4_=NULL;
  free(dim_5_); dim_5_=NULL;
  free(dim_4_); dim_4_=NULL;
  free(dim_3_); dim_3_=NULL;
  free(dim_2_); dim_2_=NULL;
  free(dim_1_); dim_1_=NULL;
  free(c16_0_); c16_0_=NULL;
  free(c16_5_); c16_5_=NULL;
  free(r4_3_); r4_3_=NULL;
  free(r8_2_); r8_2_=NULL;
  free(i4_1_); i4_1_=NULL;
  /* %%%%%%%% */
  if (verbose){ printf(" %% [finished MDA_io_test]\n");}
}

