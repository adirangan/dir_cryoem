      integer *4 n_svd_r,n_svd_d,n_svd_l !used as n_svd_r_degree and n_svd_d_degree ;
      integer *4 nsvd_r,nsvd_d,nsvd_l !used as nsvd_r_degree and nsvd_d_degree ;
      real *8 svd_r_(0:0) !unused ;
      real *8 svd_d_(0:0) !unused ;
      integer *4 svd_l_(0:0) !l-index ;
      real *8 svd_U_d_(0:0) !real *8 array of size n_svd_d_degree -x- n_svd_l ; polynomial coefficients for evaluation of U_d_ ; note that these polynomials should be evaluated in the [-1,+1] interval (e.g., at (d_value - d_m) / d_c ) ;
      real *8 svd_s_(0:0) !real *8 array of size n_svd_l ; s-values ;
      real *8 svd_V_r_(0:0) !real *8 array of size n_svd_r_degree -x- n_svd_l ; polynomial coefficients for evaluation of V_r_ ; note that these polynomials should be evaluated in the [-1,+1] interval (e.g., at (r_value - r_m) / r_c ) ;
      integer *4 svd_unitnumber
      parameter (svd_unitnumber=143857)
      character(len=64) svd_fname,svd_format_string

