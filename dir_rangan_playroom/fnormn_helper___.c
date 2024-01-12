double fnormn_helper_z_vs_ff___(int n_0,int n_1,int n_2,int flag_rup_A,double complex *z_A_,int flag_rup_B,float *f_B_real_,float *f_B_imag_)
{
  int n_0_rup=0,n_1_rup=0,n_2_rup=0;
  int n_0_A_use=0,n_1_A_use=0,n_2_A_use=0;
  int n_0_B_use=0,n_1_B_use=0,n_2_B_use=0;
  int n0=0,n1=0,n2=0;
  unsigned long long int tabA=0,tabB=0;
  double d_fnorm=0.0,d_denom=0.0,d_fnormn=0.0;
  double complex z_A;
  float f_B_real=0,f_B_imag=0;
  n_0_rup = rup(n_0,8); n_1_rup = rup(n_1,8); n_2_rup = rup(n_2,8);
  n_0_A_use = n_0; if (flag_rup_A){ n_0_A_use = n_0_rup;}
  n_1_A_use = n_1; if (flag_rup_A){ n_1_A_use = n_1_rup;}
  n_2_A_use = n_2; if (flag_rup_A){ n_2_A_use = n_2_rup;}
  n_0_B_use = n_0; if (flag_rup_B){ n_0_B_use = n_0_rup;}
  n_1_B_use = n_1; if (flag_rup_B){ n_1_B_use = n_1_rup;}
  n_2_B_use = n_2; if (flag_rup_B){ n_2_B_use = n_2_rup;}
  d_fnorm=0.0;d_denom=0.0;
    for (n2=0;n2<n_2;n2++){
      for (n1=0;n1<n_1;n1++){
	for (n0=0;n0<n_0;n0++){
	  tabA =
	    (unsigned long long int)n0 +
	    ((unsigned long long int)n1 +
	     (unsigned long long int)n2
	     *(unsigned long long int)n_1_A_use)
	    *(unsigned long long int)n_0_A_use;
	  tabB =
	    (unsigned long long int)n0 +
	    ((unsigned long long int)n1 +
	     (unsigned long long int)n2
	     *(unsigned long long int)n_1_B_use)
	    *(unsigned long long int)n_0_B_use;
	  z_A = z_A_[tabA];
	  f_B_real = f_B_real_[tabB];
	  f_B_imag = f_B_imag_[tabB];
	  d_fnorm += (creal(z_A) - f_B_real)*(creal(z_A) - f_B_real) + (cimag(z_A) - f_B_imag)*(cimag(z_A) - f_B_imag);
	  d_denom += creal(z_A)*creal(z_A) + cimag(z_A)*cimag(z_A);
	  /* for (n0=0;n0<n_0;n0++){ } */}
	/* for (n1=0;n1<n_1;n1++){ } */}
      /* for (n2=0;n2<n_2;n2++){ } */}
  d_fnormn = d_fnorm / maximum(1e-12,d_denom);
  return d_fnormn;
}
double fnormn_helper_z_vs_c___(int n_0,int n_1,int n_2,int flag_rup_A,double complex *z_A_,int flag_rup_B,float complex *c_B_)
{
  int n_0_rup=0,n_1_rup=0,n_2_rup=0;
  int n_0_A_use=0,n_1_A_use=0,n_2_A_use=0;
  int n_0_B_use=0,n_1_B_use=0,n_2_B_use=0;
  int n0=0,n1=0,n2=0;
  unsigned long long int tabA=0,tabB=0;
  double d_fnorm=0.0,d_denom=0.0,d_fnormn=0.0;
  double complex z_A;
  float complex c_B;
  n_0_rup = rup(n_0,8); n_1_rup = rup(n_1,8); n_2_rup = rup(n_2,8);
  n_0_A_use = n_0; if (flag_rup_A){ n_0_A_use = n_0_rup;}
  n_1_A_use = n_1; if (flag_rup_A){ n_1_A_use = n_1_rup;}
  n_2_A_use = n_2; if (flag_rup_A){ n_2_A_use = n_2_rup;}
  n_0_B_use = n_0; if (flag_rup_B){ n_0_B_use = n_0_rup;}
  n_1_B_use = n_1; if (flag_rup_B){ n_1_B_use = n_1_rup;}
  n_2_B_use = n_2; if (flag_rup_B){ n_2_B_use = n_2_rup;}
  d_fnorm=0.0;d_denom=0.0;
    for (n2=0;n2<n_2;n2++){
      for (n1=0;n1<n_1;n1++){
	for (n0=0;n0<n_0;n0++){
	  tabA =
	    (unsigned long long int)n0 +
	    ((unsigned long long int)n1 +
	     (unsigned long long int)n2
	     *(unsigned long long int)n_1_A_use)
	    *(unsigned long long int)n_0_A_use;
	  tabB =
	    (unsigned long long int)n0 +
	    ((unsigned long long int)n1 +
	     (unsigned long long int)n2
	     *(unsigned long long int)n_1_B_use)
	    *(unsigned long long int)n_0_B_use;
	  z_A = z_A_[tabA];
	  c_B = c_B_[tabB];
	  d_fnorm += (creal(z_A) - crealf(c_B))*(creal(z_A) - crealf(c_B)) + (cimag(z_A) - cimagf(c_B))*(cimag(z_A) - cimagf(c_B));
	  d_denom += creal(z_A)*creal(z_A) + cimag(z_A)*cimag(z_A);
	  /* for (n0=0;n0<n_0;n0++){ } */}
	/* for (n1=0;n1<n_1;n1++){ } */}
      /* for (n2=0;n2<n_2;n2++){ } */}
  d_fnormn = d_fnorm / maximum(1e-12,d_denom);
  return d_fnormn;
}
double fnormn_helper_ff_vs_c___(int n_0,int n_1,int n_2,int flag_rup_A,float *f_A_real_,float *f_A_imag_,int flag_rup_B,float complex *c_B_)
{
  int n_0_rup=0,n_1_rup=0,n_2_rup=0;
  int n_0_A_use=0,n_1_A_use=0,n_2_A_use=0;
  int n_0_B_use=0,n_1_B_use=0,n_2_B_use=0;
  int n0=0,n1=0,n2=0;
  unsigned long long int tabA=0,tabB=0;
  double d_fnorm=0.0,d_denom=0.0,d_fnormn=0.0;
  float f_A_real=0,f_A_imag=0;
  float complex c_B;
  n_0_rup = rup(n_0,8); n_1_rup = rup(n_1,8); n_2_rup = rup(n_2,8);
  n_0_A_use = n_0; if (flag_rup_A){ n_0_A_use = n_0_rup;}
  n_1_A_use = n_1; if (flag_rup_A){ n_1_A_use = n_1_rup;}
  n_2_A_use = n_2; if (flag_rup_A){ n_2_A_use = n_2_rup;}
  n_0_B_use = n_0; if (flag_rup_B){ n_0_B_use = n_0_rup;}
  n_1_B_use = n_1; if (flag_rup_B){ n_1_B_use = n_1_rup;}
  n_2_B_use = n_2; if (flag_rup_B){ n_2_B_use = n_2_rup;}
  d_fnorm=0.0;d_denom=0.0;
    for (n2=0;n2<n_2;n2++){
      for (n1=0;n1<n_1;n1++){
	for (n0=0;n0<n_0;n0++){
	  tabA =
	    (unsigned long long int)n0 +
	    ((unsigned long long int)n1 +
	     (unsigned long long int)n2
	     *(unsigned long long int)n_1_A_use)
	    *(unsigned long long int)n_0_A_use;
	  tabB =
	    (unsigned long long int)n0 +
	    ((unsigned long long int)n1 +
	     (unsigned long long int)n2
	     *(unsigned long long int)n_1_B_use)
	    *(unsigned long long int)n_0_B_use;
	  f_A_real = f_A_real_[tabA];
	  f_A_imag = f_A_imag_[tabA];
	  c_B = c_B_[tabB];
	  d_fnorm += (crealf(c_B) - f_A_real)*(crealf(c_B) - f_A_real) + (cimagf(c_B) - f_A_imag)*(cimagf(c_B) - f_A_imag);
	  d_denom += f_A_real*f_A_real + f_A_imag*f_A_imag;
	  /* for (n0=0;n0<n_0;n0++){ } */}
	/* for (n1=0;n1<n_1;n1++){ } */}
      /* for (n2=0;n2<n_2;n2++){ } */}
  d_fnormn = d_fnorm / maximum(1e-12,d_denom);
  return d_fnormn;
}
