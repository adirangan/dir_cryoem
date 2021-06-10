!> Doxygen comment: ;\n
!> Defines grid_k_p_r_ and weight_k_p_r_. ;\n
!> Note that grid_k_p_r_ does not start at 0 (i.e., the zero-frequency is ommitted). ;\n
!> Doxygen comment: ;\n
      subroutine get_grid_k_p_r_(
     $     k_p_r_max
     $     ,n_k_p_r_max
     $     ,quadrature_type_radial
     $     ,grid_k_p_r_
     $     ,weight_k_p_r_
     $     )
      implicit none
      real *8 k_p_r_max ! maximum value of k in k-space polar coordinates. sometimes called rmax. ;
      integer n_k_p_r_max ! maximum number of k-values (in k-space polar coordinates). sometimes named ngridr. ;
      integer quadrature_type_radial ! quadrature type in the radial direction. sometimes named ityper. ;
      real *8 grid_k_p_r_(0:0) ! values for k on successive shells for k-space polar coordinates. sometimes called xnodesr(nk). ;
      real *8 weight_k_p_r_(0:0) ! weight associated with grid_k_p_r_(nk) in k-space polar coordinates. sometimes called wtsr(nk). ;
      real *8, allocatable :: tmp_(:) ! temporary array to hold quadrature output. ;
      integer nk
      if (quadrature_type_radial.eq.1) then
         allocate(tmp_(0:n_k_p_r_max-1))
         call gaussq(0,n_k_p_r_max,0,0,0,0,tmp_,grid_k_p_r_
     $        ,weight_k_p_r_)
         do nk=0,n_k_p_r_max-1
            grid_k_p_r_(nk) = (k_p_r_max/2.0d0)*(grid_k_p_r_(nk)+1.0d0)
            weight_k_p_r_(nk) = grid_k_p_r_(nk)*grid_k_p_r_(nk)
     $           *weight_k_p_r_(nk)*(k_p_r_max/2.0d0)
         enddo !do nk=0,n_k_p_r_max-1
      end if !if (quadrature_type_radial.eq.1) then
      if (quadrature_type_radial.eq.0) then
         do nk=0,n_k_p_r_max-1
            grid_k_p_r_(nk) = k_p_r_max*(nk+0.0d0)/n_k_p_r_max
            weight_k_p_r_(nk) = grid_k_p_r_(nk)*grid_k_p_r_(nk)
     $           *k_p_r_max/n_k_p_r_max
         enddo !do nk=0,n_k_p_r_max-1
      end if !if (quadrature_type_radial.eq.0) then
      end !subroutine
