      integer *4 n_svd_r,n_svd_d,n_svd_l
      integer *4 nsvd_r,nsvd_d,nsvd_l
      real *8, allocatable :: svd_r_(:)
      real *8, allocatable :: svd_d_(:)
      integer *4, allocatable :: svd_l_(:)
      real *8, allocatable :: svd_U_d_(:)
      real *8, allocatable :: svd_s_(:)
      real *8, allocatable :: svd_V_r_(:)
      integer *4 svd_unitnumber
      parameter (svd_unitnumber=143857)
      character(len=64) svd_fname,svd_format_string
