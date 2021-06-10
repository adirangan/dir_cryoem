!> Doxygen comment: ;\n
!> Defines n_polar_a_, n_w_, n_Y_lm_sum, ld_M, n_polar_a_, n_w_, n_Y_lm_csum_ n_Y_l_ and n_Y_lm_. ;\n
!> Doxygen comment: ;\n
      subroutine get_n_polar_a_(
     $     n_k_p_max
     $     ,grid_k_p_r_
     $     ,n_polar_a_lowerbound
     $     ,n_polar_a_
     $     ,n_w_
     $     ,n_Y_l_
     $     ,n_Y_lm_
     $     ,n_Y_lm_csum_
     $     ,ld_M
     $     ,n_Y_lm_sum
     $     )
      implicit none
      integer n_k_p_max ! maximum number of k-values (in k-space polar coordinates). sometimes named ngridr. ;
      real *8 grid_k_p_r_(0:0) !real *8 array (size at least n_k_p_max): values for k on successive shells for k-space polar coordinates. sometimes called grid_p_ or xnodesr. ;
      integer *4 n_polar_a_lowerbound !integer *4 minimum value for n_polar_a_(nk). ;
      integer *4 n_polar_a_(0:0) !integer *4 array (size at least n_k_p_max): number of polar_a values on sphere of radius grid_k_p_r_(nk). also called nlats(nk). ;
      integer *4 n_w_(0:0) !integer *4 array (size at least n_k_p_max): n_w_(nk) is the number of points on the ring of radius grid_k_p_r_(nk) for a template or image in k-space polar coordinates. also called ngridc(nk). ;
      integer *4 n_Y_l_(0:0) !integer *4 array (size at least n_k_p_max): n_Y_l_(nk) is the order of the spherical harmonic expansion on sphere at radius grid_k_p_r_(nk). also called nterms_sph(nk). ;
      integer *4 n_Y_lm_(0:0) !integer *4 array (size at least n_k_p_max): n_Y_lm_(nk) is the number of terms (i.e., (l,m) pairs) of the spherical harmonic expansion on sphere at radius grid_k_p_r_(nk). also called nterms_sph(nk). ;
      integer *4 n_Y_lm_csum_(0:0) !integer *4 array (size at least n_k_p_max): n_Y_lm_csum_(nk) is the cumulative sum of n_Y_lm_(0:nk-1). ;
      integer ld_M !integer equal to the sum of n_w_. This is the total number of entries in each image (in k-space polar coordinates). i.e., leading-dimension of M array. sometimes called n_image_size or nimagesize. ;
      integer n_Y_lm_sum !integer equal to the sum of n_Y_lm_. This is the total number of basis functions used in spherical harmonic expansion (in k-space polar coordinates). sometimes called nsphstore. ;
      integer nk
      real *8 pi
      pi=4.0d0*datan(1.0d0)
      n_Y_lm_sum = 0
      ld_M = 0
      do nk = 0,n_k_p_max-1
         n_polar_a_(nk) = max(
     $        n_polar_a_lowerbound
     $        ,nint(pi*grid_k_p_r_(nk))
     $        )
         if (n_polar_a_(nk).lt.6) n_polar_a_(nk) = 6
         if (mod(n_polar_a_(nk),2).ne.0) n_polar_a_(nk) = n_polar_a_(nk)
     $        + 1
         n_w_(nk) = n_polar_a_(nk)*2
         n_Y_lm_csum_(nk) = n_Y_lm_sum
         n_Y_l_(nk) = nint(grid_k_p_r_(nk) + 2)
         n_Y_lm_(nk) = (n_Y_l_(nk)+1)**2
         n_Y_lm_sum = n_Y_lm_sum + n_Y_lm_(nk)
         ld_M = ld_M + n_w_(nk)
      enddo
      end !subroutine
