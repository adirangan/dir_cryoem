#ifndef _MONOLITH
#include "ampm_header.h"
#endif /* _MONOLITH */

/* %%%%%%%% */
/* array_printf: */
/* %%%%%%%% */

void iarray_printf_margin(int *i_,int n_r,int n_c,const char *prefix)
{
  int nr=0,nc=0,margin=3;
  for (nr=0;nr<n_r;nr++){
    if ( (nr>margin-1) && (nr<n_r-margin) ){ printf("%s...............................................................\n",prefix); nr = maximum(0,n_r-margin);}
    if ( (nr<margin) || (nr>n_r-margin-1) ){
      printf("%s",prefix);
      for (nc=0;nc<n_c;nc++){
	if ( (nc>margin-1) && (nc<n_c-margin) ){ printf("... "); nc = maximum(0,n_c-margin);}
	if ( (nc<margin) || (nc>n_c-margin-1) ){
	  printf("%+.6d ",i_[nr+nc*n_r]);
	  /* margin */}
	/* for (nc=0;nc<n_c;nc++){ } */}
      printf("\n");
      /* margin */}
    /* for (nr=0;nr<n_r;nr++){ } */}
}

void farray_aligned_printf_margin(float *f_,int n_r,int n_c,const char *prefix)
{
  int nr=0,nc=0,margin=3;
  int n_r_rup = rup(n_r,8);
  for (nr=0;nr<n_r;nr++){
    if ( (nr>margin-1) && (nr<n_r-margin) ){ printf("%s...............................................................\n",prefix); nr = maximum(0,n_r-margin);}
    if ( (nr<margin) || (nr>n_r-margin-1) ){
      printf("%s",prefix);
      for (nc=0;nc<n_c;nc++){
	if ( (nc>margin-1) && (nc<n_c-margin) ){ printf("... "); nc = maximum(0,n_c-margin);}
	if ( (nc<margin) || (nc>n_c-margin-1) ){
	  printf("%+.6f ",f_[nr+nc*n_r_rup]);
	  /* margin */}
	/* for (nc=0;nc<n_c;nc++){ } */}
      printf("\n");
      /* margin */}
    /* for (nr=0;nr<n_r;nr++){ } */}
}

void farray_printf_margin(float *f_,int n_r,int n_c,const char *prefix)
{
  int nr=0,nc=0,margin=3;
  for (nr=0;nr<n_r;nr++){
    if ( (nr>margin-1) && (nr<n_r-margin) ){ printf("%s...............................................................\n",prefix); nr = maximum(0,n_r-margin);}
    if ( (nr<margin) || (nr>n_r-margin-1) ){
      printf("%s",prefix);
      for (nc=0;nc<n_c;nc++){
	if ( (nc>margin-1) && (nc<n_c-margin) ){ printf("... "); nc = maximum(0,n_c-margin);}
	if ( (nc<margin) || (nc>n_c-margin-1) ){
	  printf("%+.6f ",f_[nr+nc*n_r]);
	  /* margin */}
	/* for (nc=0;nc<n_c;nc++){ } */}
      printf("\n");
      /* margin */}
    /* for (nr=0;nr<n_r;nr++){ } */}
}

void darray_printf_margin(double *d_,int n_r,int n_c,const char *prefix)
{
  int nr=0,nc=0,margin=3;
  for (nr=0;nr<n_r;nr++){
    if ( (nr>margin-1) && (nr<n_r-margin) ){ printf("%s...............................................................\n",prefix); nr = maximum(0,n_r-margin);}
    if ( (nr<margin) || (nr>n_r-margin-1) ){
      printf("%s",prefix);
      for (nc=0;nc<n_c;nc++){
	if ( (nc>margin-1) && (nc<n_c-margin) ){ printf("... "); nc = maximum(0,n_c-margin);}
	if ( (nc<margin) || (nc>n_c-margin-1) ){
	  printf("%+.6f ",d_[nr+nc*n_r]);
	  /* margin */}
	/* for (nc=0;nc<n_c;nc++){ } */}
      printf("\n");
      /* margin */}
    /* for (nr=0;nr<n_r;nr++){ } */}
}

void array_printf_margin(void *v_,const char *type,int n_row,int n_col,const char *prefix)
{
  if (strcmp(type,"int")==0){ iarray_printf_margin((int *)v_,n_row,n_col,prefix);}
  if (strcmp(type,"float")==0){ farray_printf_margin((float *)v_,n_row,n_col,prefix);}
  if (strcmp(type,"double")==0){ darray_printf_margin((double *)v_,n_row,n_col,prefix);}
}

void array_sub_printf(
		      void *v_
		      ,const char *type
		      ,int n_row
		      ,int n_row_sub
		      ,int n_col
		      ,int n_col_sub
		      ,const char *prefix
		      )
{
  /* prints out sub-arrays of varying types */
  int nr=0,nc=0,tmp=0;
  float *f_=NULL;
  double *d_=NULL; int d_l_flag=0;
  int *i_=NULL;
  unsigned int *ui_=NULL;
  long *l_=NULL;
  long long *ll_=NULL;
  unsigned long int *ul_=NULL;
  unsigned long long int *ull_=NULL;
  char *char_=NULL;
  unsigned char *uchar_=NULL;
  float complex *c_=NULL;
  double complex *z_=NULL;
  double printftol=0.000000001;
  n_row_sub = minimum(n_row_sub,n_row);
  n_col_sub = minimum(n_col_sub,n_col);
  if (n_row_sub<=0){ n_row_sub=n_row;}
  if (n_col_sub<=0){ n_col_sub=n_col;}
  if (strcmp(type,"float")==0){
    f_ = (float *) v_;
    for (nr=0;nr<n_row_sub;nr++){ 
      printf("%s",prefix); 
      for (nc=0;nc<n_col_sub;nc++){ 
	if (fabs(f_[nr+nc*n_row]-(int)f_[nr+nc*n_row])<printftol){ printf(" %s%d",(int)f_[nr+nc*n_row]>0 ? "+" : ((int)f_[nr+nc*n_row]<0 ? "" : " "),(int)f_[nr+nc*n_row]);}
	else{ printf(" %s%f",f_[nr+nc*n_row]>0 ? "+" : (f_[nr+nc*n_row]<0 ? "" : " "),f_[nr+nc*n_row]);}}
      printf("\n");}}
  else if (strcmp(type,"double")==0){
    if (n_row>1 && n_col==1){ tmp=n_row; n_row=n_col; n_col=tmp;}
    d_ = (double *) v_;
    if (0){
      d_l_flag=1; for (nr=0;nr<n_row_sub;nr++){ for (nc=0;nc<n_col_sub;nc++){ if (fabs(d_[nr+nc*n_row]-(int)d_[nr+nc*n_row])>printftol){ d_l_flag=0;}}}
      for (nr=0;nr<n_row_sub;nr++){ 
	printf("%s",prefix); 
	for (nc=0;nc<n_col_sub;nc++){ 
	  if (d_l_flag && fabs(d_[nr+nc*n_row]-(int)d_[nr+nc*n_row])<printftol){ printf(" %s%d",(int)d_[nr+nc*n_row]>0 ? "+" : ((int)d_[nr+nc*n_row]<0 ? "" : " "),(int)d_[nr+nc*n_row]);}
	  else{ /* printf(" %s%.3f",d_[nr+nc*n_row]>0 ? "+" : (d_[nr+nc*n_row]<0 ? "" : " "),d_[nr+nc*n_row]); */ 
	    if (fabs(d_[nr+nc*n_row])<printftol){ printf("    .   ");} else{ printf(" %+07.3f",d_[nr+nc*n_row]);}
	    /* if !d_l_flag */}
	  /* for (nc=0;nc<n_col_sub;nc++){ } */}
	printf("\n");} /* if 0 */}
    if (1){
      for (nr=0;nr<n_row_sub;nr++){ printf("%s",prefix); 
	for (nc=0;nc<n_col_sub;nc++){ if (fabs(d_[nr+nc*n_row])<printftol){ printf("    .   ");} else{ printf(" %+07.3f",d_[nr+nc*n_row]);} /* for (nc=0;nc<n_col_sub;nc++){ } */}
	printf("\n"); /* for (nr=0;nr<n_row_sub;nr++){ } */} /* if 1 */}
    /* if double */}
  else if (strcmp(type,"double_trn")==0 && n_row>1){
    d_ = (double *) v_;
    for (nc=0;nc<n_col_sub;nc++){ printf("%s",prefix); 
      for (nr=0;nr<n_row_sub;nr++){ if (fabs(d_[nr+nc*n_row])<printftol){ printf("    .   ");} else{ printf(" %+07.3f",d_[nr+nc*n_row]);} /* for (nr=0;nr<n_row_sub;nr++){ } */}
      printf("\n"); /* for (nc=0;nc<n_col_sub;nc++){ } */}
    /* if double_trn */}
  else if (strcmp(type,"float complex")==0){
    if (n_row>1 && n_col==1){ tmp=n_row; n_row=n_col; n_col=tmp;}
    c_ = (float complex *) v_;
    for (nr=0;nr<n_row_sub;nr++){ printf("%s",prefix); 
      for (nc=0;nc<n_col_sub;nc++){
	if (cabsf(c_[nr+nc*n_row])<printftol){ printf("    .        .   ");} else{ printf(" [%+07.3f %+07.3fi]",crealf(c_[nr+nc*n_row]),cimagf(c_[nr+nc*n_row]));} 
	/* for (nc=0;nc<n_col_sub;nc++){ } */}
      printf("\n"); /* for (nr=0;nr<n_row_sub;nr++){ } */}
    /* else if (strcmp(type,"float complex")==0){ } */}
  else if (strcmp(type,"double complex")==0){
    if (n_row>1 && n_col==1){ tmp=n_row; n_row=n_col; n_col=tmp;}
    z_ = (double complex *) v_;
    for (nr=0;nr<n_row_sub;nr++){ printf("%s",prefix); 
      for (nc=0;nc<n_col_sub;nc++){
	if (cabs(z_[nr+nc*n_row])<printftol){ printf("    .        .   ");} else{ printf(" [%+07.3f %+07.3fi]",creal(z_[nr+nc*n_row]),cimag(z_[nr+nc*n_row]));} 
	/* for (nc=0;nc<n_col_sub;nc++){ } */}
      printf("\n"); /* for (nr=0;nr<n_row_sub;nr++){ } */}
    /* else if (strcmp(type,"double complex")==0){ } */}
  else if (strcmp(type,"long")==0){
    l_ = (long *) v_;
    for (nr=0;nr<n_row_sub;nr++){ printf("%s",prefix); for (nc=0;nc<n_col_sub;nc++){ printf(" %ld",l_[nr+nc*n_row]);} printf("\n");}}
  else if (strcmp(type,"long long int")==0){
    ll_ = (long long *) v_;
    for (nr=0;nr<n_row_sub;nr++){ printf("%s",prefix); for (nc=0;nc<n_col_sub;nc++){ printf(" %lld",ll_[nr+nc*n_row]);} printf("\n");}}
  else if (strcmp(type,"int")==0){
    i_ = (int *) v_;
    for (nr=0;nr<n_row_sub;nr++){ printf("%s",prefix); for (nc=0;nc<n_col_sub;nc++){ printf(" %d",i_[nr+nc*n_row]);} printf("\n");}}
  else if (strcmp(type,"unsigned int")==0){
    ui_ = (unsigned int *) v_;
    for (nr=0;nr<n_row_sub;nr++){ printf("%s",prefix); for (nc=0;nc<n_col_sub;nc++){ printf(" %d",(int)ui_[nr+nc*n_row]);} printf("\n");}}
  else if (strcmp(type,"unsigned long int")==0){
    ul_ = (unsigned long int *) v_;
    for (nr=0;nr<n_row_sub;nr++){ printf("%s",prefix); for (nc=0;nc<n_col_sub;nc++){ printf(" %ld",ul_[nr+nc*n_row]);} printf("\n");}}
  else if (strcmp(type,"unsigned long long int")==0){
    ull_ = (unsigned long long int *) v_;
    for (nr=0;nr<n_row_sub;nr++){ printf("%s",prefix); for (nc=0;nc<n_col_sub;nc++){ printf(" %lld",ull_[nr+nc*n_row]);} printf("\n");}}
  else if (strcmp(type,"alpha")==0){
    char_ = (char *) v_;
    for (nr=0;nr<n_row_sub;nr++){ printf("%s",prefix); for (nc=0;nc<n_col_sub;nc++){ printf(" %c",(int)char_[nr+nc*n_row]);} printf("\n");}}
  else if (strcmp(type,"char")==0){
    char_ = (char *) v_;
    for (nr=0;nr<n_row_sub;nr++){ printf("%s",prefix); for (nc=0;nc<n_col_sub;nc++){ printf(" %d",(int)char_[nr+nc*n_row]);} printf("\n");}}
  else if (strcmp(type,"unsigned char")==0){
    uchar_ = (unsigned char *) v_;
    for (nr=0;nr<n_row_sub;nr++){ printf("%s",prefix); for (nc=0;nc<n_col_sub;nc++){ printf(" %d",(int)uchar_[nr+nc*n_row]);} printf("\n");}}
  else{ printf(" warning, poor type %s in array_printf\n",type);}
}

void array_printf(
		  void *v_
		  ,const char *type
		  ,int n_row
		  ,int n_col
		  ,const char *prefix
		  )
{
  /* prints out arrays of varying types */
  int nr=0,nc=0,tmp=0;
  float *f_=NULL;
  double *d_=NULL; int d_l_flag=0;
  int *i_=NULL;
  unsigned int *ui_=NULL;
  long *l_=NULL;
  long long *ll_=NULL;
  unsigned long int *ul_=NULL;
  unsigned long long int *ull_=NULL;
  char *char_=NULL;
  unsigned char *uchar_=NULL;
  float complex *c_=NULL;
  double complex *z_=NULL;
  double printftol=0.000000001;
  if (strcmp(type,"float")==0){
    f_ = (float *) v_;
    for (nr=0;nr<n_row;nr++){ 
      printf("%s",prefix); 
      for (nc=0;nc<n_col;nc++){ 
	if (fabs(f_[nr+nc*n_row]-(int)f_[nr+nc*n_row])<printftol){ printf(" %s%d",(int)f_[nr+nc*n_row]>0 ? "+" : ((int)f_[nr+nc*n_row]<0 ? "" : " "),(int)f_[nr+nc*n_row]);}
	else{ printf(" %s%f",f_[nr+nc*n_row]>0 ? "+" : (f_[nr+nc*n_row]<0 ? "" : " "),f_[nr+nc*n_row]);}}
      printf("\n");}}
  else if (strcmp(type,"double")==0){
    if (n_row>1 && n_col==1){ tmp=n_row; n_row=n_col; n_col=tmp;}
    d_ = (double *) v_;
    if (0){
      d_l_flag=1; for (nr=0;nr<n_row;nr++){ for (nc=0;nc<n_col;nc++){ if (fabs(d_[nr+nc*n_row]-(int)d_[nr+nc*n_row])>printftol){ d_l_flag=0;}}}
      for (nr=0;nr<n_row;nr++){ 
	printf("%s",prefix); 
	for (nc=0;nc<n_col;nc++){ 
	  if (d_l_flag && fabs(d_[nr+nc*n_row]-(int)d_[nr+nc*n_row])<printftol){ printf(" %s%d",(int)d_[nr+nc*n_row]>0 ? "+" : ((int)d_[nr+nc*n_row]<0 ? "" : " "),(int)d_[nr+nc*n_row]);}
	  else{ /* printf(" %s%.3f",d_[nr+nc*n_row]>0 ? "+" : (d_[nr+nc*n_row]<0 ? "" : " "),d_[nr+nc*n_row]); */ 
	    if (fabs(d_[nr+nc*n_row])<printftol){ printf("    .   ");} else{ printf(" %+07.3f",d_[nr+nc*n_row]);}
	    /* if !d_l_flag */}
	  /* for (nc=0;nc<n_col;nc++){ } */}
	printf("\n");} /* if 0 */}
    if (1){
      for (nr=0;nr<n_row;nr++){ printf("%s",prefix); 
	for (nc=0;nc<n_col;nc++){ if (fabs(d_[nr+nc*n_row])<printftol){ printf("    .   ");} else{ printf(" %+07.3f",d_[nr+nc*n_row]);} /* for (nc=0;nc<n_col;nc++){ } */}
	printf("\n"); /* for (nr=0;nr<n_row;nr++){ } */} /* if 1 */}
    /* if double */}
  else if (strcmp(type,"double_trn")==0 && n_row>1){
    d_ = (double *) v_;
    for (nc=0;nc<n_col;nc++){ printf("%s",prefix); 
      for (nr=0;nr<n_row;nr++){ if (fabs(d_[nr+nc*n_row])<printftol){ printf("    .   ");} else{ printf(" %+07.3f",d_[nr+nc*n_row]);} /* for (nr=0;nr<n_row;nr++){ } */}
      printf("\n"); /* for (nc=0;nc<n_col;nc++){ } */}
    /* if double_trn */}
  else if (strcmp(type,"float complex")==0){
    if (n_row>1 && n_col==1){ tmp=n_row; n_row=n_col; n_col=tmp;}
    c_ = (float complex *) v_;
    for (nr=0;nr<n_row;nr++){ printf("%s",prefix); 
      for (nc=0;nc<n_col;nc++){
	if (cabsf(c_[nr+nc*n_row])<printftol){ printf("    .        .   ");} else{ printf(" [%+07.3f %+07.3fi]",crealf(c_[nr+nc*n_row]),cimagf(c_[nr+nc*n_row]));} 
	/* for (nc=0;nc<n_col;nc++){ } */}
      printf("\n"); /* for (nr=0;nr<n_row;nr++){ } */}
    /* else if (strcmp(type,"float complex")==0){ } */}
  else if (strcmp(type,"double complex")==0){
    if (n_row>1 && n_col==1){ tmp=n_row; n_row=n_col; n_col=tmp;}
    z_ = (double complex *) v_;
    for (nr=0;nr<n_row;nr++){ printf("%s",prefix); 
      for (nc=0;nc<n_col;nc++){
	if (cabs(z_[nr+nc*n_row])<printftol){ printf("    .        .   ");} else{ printf(" [%+07.3f %+07.3fi]",creal(z_[nr+nc*n_row]),cimag(z_[nr+nc*n_row]));} 
	/* for (nc=0;nc<n_col;nc++){ } */}
      printf("\n"); /* for (nr=0;nr<n_row;nr++){ } */}
    /* else if (strcmp(type,"double complex")==0){ } */}
  else if (strcmp(type,"long")==0){
    l_ = (long *) v_;
    for (nr=0;nr<n_row;nr++){ printf("%s",prefix); for (nc=0;nc<n_col;nc++){ printf(" %ld",l_[nr+nc*n_row]);} printf("\n");}}
  else if (strcmp(type,"long long int")==0){
    ll_ = (long long *) v_;
    for (nr=0;nr<n_row;nr++){ printf("%s",prefix); for (nc=0;nc<n_col;nc++){ printf(" %lld",ll_[nr+nc*n_row]);} printf("\n");}}
  else if (strcmp(type,"int")==0){
    i_ = (int *) v_;
    for (nr=0;nr<n_row;nr++){ printf("%s",prefix); for (nc=0;nc<n_col;nc++){ printf(" %d",i_[nr+nc*n_row]);} printf("\n");}}
  else if (strcmp(type,"unsigned int")==0){
    ui_ = (unsigned int *) v_;
    for (nr=0;nr<n_row;nr++){ printf("%s",prefix); for (nc=0;nc<n_col;nc++){ printf(" %d",(int)ui_[nr+nc*n_row]);} printf("\n");}}
  else if (strcmp(type,"unsigned long int")==0){
    ul_ = (unsigned long int *) v_;
    for (nr=0;nr<n_row;nr++){ printf("%s",prefix); for (nc=0;nc<n_col;nc++){ printf(" %ld",ul_[nr+nc*n_row]);} printf("\n");}}
  else if (strcmp(type,"unsigned long long int")==0){
    ull_ = (unsigned long long int *) v_;
    for (nr=0;nr<n_row;nr++){ printf("%s",prefix); for (nc=0;nc<n_col;nc++){ printf(" %lld",ull_[nr+nc*n_row]);} printf("\n");}}
  else if (strcmp(type,"alpha")==0){
    char_ = (char *) v_;
    for (nr=0;nr<n_row;nr++){ printf("%s",prefix); for (nc=0;nc<n_col;nc++){ printf(" %c",(int)char_[nr+nc*n_row]);} printf("\n");}}
  else if (strcmp(type,"char")==0){
    char_ = (char *) v_;
    for (nr=0;nr<n_row;nr++){ printf("%s",prefix); for (nc=0;nc<n_col;nc++){ printf(" %d",(int)char_[nr+nc*n_row]);} printf("\n");}}
  else if (strcmp(type,"unsigned char")==0){
    uchar_ = (unsigned char *) v_;
    for (nr=0;nr<n_row;nr++){ printf("%s",prefix); for (nc=0;nc<n_col;nc++){ printf(" %d",(int)uchar_[nr+nc*n_row]);} printf("\n");}}
  else{ printf(" warning, poor type %s in array_printf\n",type);}
}


void array_fprintf(
		   const char *fname
		   ,void *v_
		   ,const char *type
		   ,int n_row
		   ,int n_col
		   ,const char *prefix
		   )
{
  /* prints out arrays of varying types */
  FILE *fp=NULL;
  int nr=0,nc=0,tmp=0;
  float *f_=NULL;
  double *d_=NULL; int d_l_flag=1;
  int *i_=NULL;
  unsigned int *ui_=NULL;
  long *l_=NULL;
  long long *ll_=NULL;
  char *char_=NULL;
  unsigned char *uchar_=NULL;
  /* fftw_complex *ff_=NULL; */
  /* double dtmp=0; */
  double printftol=0.000000001;
  if ((fp=fopen(fname,"w"))==NULL){ printf(" %% Warning! cannot open %s in array_fprintf\n",fname); fp=stdout;}
  if (strcmp(type,"float")==0){
    f_ = (float *) v_;
    for (nr=0;nr<n_row;nr++){ 
      fprintf(fp,"%s",prefix); 
      for (nc=0;nc<n_col;nc++){ 
	if (fabs(f_[nr+nc*n_row]-(int)f_[nr+nc*n_row])<printftol){ fprintf(fp," %s%d",(int)f_[nr+nc*n_row]>0 ? "+" : ((int)f_[nr+nc*n_row]<0 ? "" : " "),(int)f_[nr+nc*n_row]);}
	else{ fprintf(fp," %s%f",f_[nr+nc*n_row]>0 ? "+" : (f_[nr+nc*n_row]<0 ? "" : " "),f_[nr+nc*n_row]);}}
      fprintf(fp,"\n");}}
  else if (strcmp(type,"double")==0){
    if (n_row>1 && n_col==1){ tmp=n_row; n_row=n_col; n_col=tmp;}
    d_ = (double *) v_;
    if (0){
      d_l_flag=1; for (nr=0;nr<n_row;nr++){ for (nc=0;nc<n_col;nc++){ if (fabs(d_[nr+nc*n_row]-(int)d_[nr+nc*n_row])>printftol){ d_l_flag=0;}}}
      for (nr=0;nr<n_row;nr++){ 
	fprintf(fp,"%s",prefix); 
	for (nc=0;nc<n_col;nc++){ 
	  if (d_l_flag && fabs(d_[nr+nc*n_row]-(int)d_[nr+nc*n_row])<printftol){ fprintf(fp," %s%d",(int)d_[nr+nc*n_row]>0 ? "+" : ((int)d_[nr+nc*n_row]<0 ? "" : " "),(int)d_[nr+nc*n_row]);}
	  else{ /* fprintf(fp," %s%.3f",d_[nr+nc*n_row]>0 ? "+" : (d_[nr+nc*n_row]<0 ? "" : " "),d_[nr+nc*n_row]); */ 
	    if (fabs(d_[nr+nc*n_row])<printftol){ fprintf(fp,"   0.0  ");} else{ fprintf(fp," %0.16f",d_[nr+nc*n_row]);}
	    /* if !d_l_flag */}
	  /* for (nc=0;nc<n_col;nc++){ } */}
	fprintf(fp,"\n");} /* if 0 */}
    if (1){
      for (nr=0;nr<n_row;nr++){ fprintf(fp,"%s",prefix); 
	for (nc=0;nc<n_col;nc++){ if (fabs(d_[nr+nc*n_row])<printftol){ fprintf(fp,"   0.0  ");} else{ fprintf(fp," %+07.3f",d_[nr+nc*n_row]);} /* for (nc=0;nc<n_col;nc++){ } */}
	fprintf(fp,"\n"); /* for (nr=0;nr<n_row;nr++){ } */} /* if 1 */}
    /* if double */}
  else if (strcmp(type,"double_trn")==0 && n_row>1){
    d_ = (double *) v_;
    for (nc=0;nc<n_col;nc++){ fprintf(fp,"%s",prefix); 
      for (nr=0;nr<n_row;nr++){ 
	if (d_l_flag && fabs(d_[nr+nc*n_row]-(int)d_[nr+nc*n_row])<printftol){ fprintf(fp," %s%d",(int)d_[nr+nc*n_row]>0 ? "+" : ((int)d_[nr+nc*n_row]<0 ? "" : " "),(int)d_[nr+nc*n_row]);}
	else{ if (fabs(d_[nr+nc*n_row])<printftol){ fprintf(fp,"   0.0  ");} else{ fprintf(fp," %0.16f",d_[nr+nc*n_row]);} /* if !d_l_flag */}
	/* for (nr=0;nr<n_row;nr++){ } */}
      fprintf(fp,"\n"); /* for (nc=0;nc<n_col;nc++){ } */}
    /* if double_trn */}
  else if (strcmp(type,"long")==0 || strcmp(type,"long int")==0){
    l_ = (long *) v_;
    for (nr=0;nr<n_row;nr++){ fprintf(fp,"%s",prefix); for (nc=0;nc<n_col;nc++){ fprintf(fp," %ld",l_[nr+nc*n_row]);} fprintf(fp,"\n");}}
  else if (strcmp(type,"long long int")==0){
    ll_ = (long long *) v_;
    for (nr=0;nr<n_row;nr++){ fprintf(fp,"%s",prefix); for (nc=0;nc<n_col;nc++){ fprintf(fp," %lld",ll_[nr+nc*n_row]);} fprintf(fp,"\n");}}
  else if (strcmp(type,"int")==0){
    i_ = (int *) v_;
    for (nr=0;nr<n_row;nr++){ fprintf(fp,"%s",prefix); for (nc=0;nc<n_col;nc++){ fprintf(fp," %d",i_[nr+nc*n_row]);} fprintf(fp,"\n");}}
  else if (strcmp(type,"unsigned int")==0){
    ui_ = (unsigned int *) v_;
    for (nr=0;nr<n_row;nr++){ fprintf(fp,"%s",prefix); for (nc=0;nc<n_col;nc++){ fprintf(fp," %d",(int)ui_[nr+nc*n_row]);} fprintf(fp,"\n");}}
  else if (strcmp(type,"alpha")==0){
    char_ = (char *) v_;
    for (nr=0;nr<n_row;nr++){ fprintf(fp,"%s",prefix); for (nc=0;nc<n_col;nc++){ fprintf(fp," %c",(int)char_[nr+nc*n_row]);} fprintf(fp,"\n");}}
  else if (strcmp(type,"char")==0){
    char_ = (char *) v_;
    for (nr=0;nr<n_row;nr++){ fprintf(fp,"%s",prefix); for (nc=0;nc<n_col;nc++){ fprintf(fp," %d",(int)char_[nr+nc*n_row]);} fprintf(fp,"\n");}}
  else if (strcmp(type,"unsigned char")==0){
    uchar_ = (unsigned char *) v_;
    for (nr=0;nr<n_row;nr++){ fprintf(fp,"%s",prefix); for (nc=0;nc<n_col;nc++){ fprintf(fp," %d",(int)uchar_[nr+nc*n_row]);} fprintf(fp,"\n");}}
  else{ fprintf(fp," warning, poor type %s in array_printf\n",type);}
  if (fp!=stdout){ fclose(fp); fp=NULL;}
}

