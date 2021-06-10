!> Doxygen comment: ;\n
!> Evaluates polynomial representation of svd_V_r_. ;\n
!> This is the 'right' side of the displacement operator (involving k, images and templates). ;\n
      subroutine get_svd_polyval_V_r_(svd_r_max,n_svd_r,svd_r_,n_svd_l
     $     ,svd_l_,svd_V_r_,n_r,grid_p_,svd_polyval_V_r_)
      implicit none
      integer verbose
      data verbose / 0 /
      real *8 svd_r_max,svd_r_m,svd_r_c,svd_r(0:0)
      integer *4 n_svd_r,n_svd_l
      real *8 svd_r_(0:0)
      integer *4 svd_l_(0:0)
      real *8 svd_V_r_(0:0)
      real *8 svd_polyval_V_r_(0:0)
      logical flag_warning
      integer n_r,nr
      real *8 grid_p_(0:0)
      integer nl
      real *8 pi
      if (verbose.gt.0) then
         write(6,'(A)') '[entering get_svd_polyval_V_r_]'
      end if !if (verbose.gt.0) then
      pi = 4.0d0*datan(1.0d0)
      svd_r_m = svd_r_max / 2.0d0
      svd_r_c = svd_r_m
      do nr=0,n_r-1
         if (grid_p_(nr).gt.svd_r_max .and. flag_warning) then
            write(6,'(A,F6.3,A,F6.3,A,F6.3,A)')
     $           'Warning, grid_p_(nr) ',grid_p_(nr),'>',svd_r_max
     $           ,'; ratio = ',grid_p_(nr)/svd_r_max
     $           ,' in get_svd_polyval_V_r_0.f'
         end if
         svd_r(0) = (grid_p_(nr) - svd_r_m)/svd_r_c
         do nl=0,n_svd_l-1
            call polyval_r8_reverse_0(n_svd_r,svd_V_r_(0+nl*n_svd_r),1
     $           ,svd_r(0),svd_polyval_V_r_(nl+nr*n_svd_l))
            if (verbose.gt.1) then
               write (6,'(A,F6.3,A,F6.3,A,I0,A,F8.5)') ' % grid_p ' ,
     $              grid_p_(nr) , ' svd_r ' ,svd_r(0) , ' V_r_(' , nl ,
     $              ') ', svd_polyval_V_r_(nl+nr*n_svd_l)
            end if !if (verbose.gt.1) then
         enddo !do nl=0,n_svd_l-1         
      enddo !do nr=0,n_r-1
      if (verbose.gt.0) then
         write(6,'(A)') '[finished get_svd_polyval_V_r_]'
      end if !if (verbose.gt.0) then
      end
