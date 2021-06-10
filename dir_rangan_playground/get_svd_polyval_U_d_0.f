!> Doxygen comment: ;\n
!> Evaluates polynomial representation of svd_U_d_. ;\n
!> This is the 'left' side of the displacement operator (involving delta). ;\n
      subroutine get_svd_polyval_U_d_(svd_d_max,n_svd_d,svd_d_,n_svd_l
     $     ,svd_l_,svd_U_d_,n_delta_v,delta_x_,delta_y_
     $     ,svd_polyval_U_d_)
      implicit none
      integer verbose
      data verbose / 0 /
      real *8 svd_d_max,svd_d_m,svd_d_c,svd_d(0:0)
      integer *4 n_svd_d,n_svd_l
      real *8 svd_d_(0:0)
      integer *4 svd_l_(0:0)
      real *8 svd_U_d_(0:0)
      real *8 svd_polyval_U_d_(0:0)
      logical flag_warning
      integer n_delta_v,ndv
      real *8 delta_x_(0:0),delta_y_(0:0)
      real *8 delta_x,delta_y,delta,omega
      integer nl
      real *8 pi
      if (verbose.gt.0) then
         write(6,'(A)') '[entering get_svd_polyval_U_d_]'
      end if !if (verbose.gt.0) then
      pi = 4.0d0*datan(1.0d0)
      svd_d_m = svd_d_max / 2.0d0
      svd_d_c = svd_d_m
      do ndv=0,n_delta_v-1
         delta_x = delta_x_(ndv)
         delta_y = delta_y_(ndv)
         if (verbose.gt.1) then
            write (6,'(A,I0,1X,F6.3,1X,F6.3)')
     $           ' % ndv dx dy : ',ndv,delta_x,delta_y
         end if                 !if (verbose.gt.1) then
         delta = dsqrt(delta_x**2 + delta_y**2)
         if (delta.gt.svd_d_max .and. flag_warning) then
            write(6,'(A,F6.3,A,F6.3,A,F6.3,A)') 'Warning, delta '
     $           ,delta,'>',svd_d_max,'; ratio = ',delta/svd_d_max,
     $           ' in get_svd_polyval_U_d_'
         end if
         omega = datan2(delta_y,delta_x)
         if (verbose.gt.1) then
            write (6,'(A,F6.3,1X,F6.3)') ' % delta omega : ',delta
     $           ,omega
         end if
         svd_d(0) = (delta - svd_d_m)/svd_d_c
         do nl=0,n_svd_l-1
            call polyval_r8_reverse_0(n_svd_d,svd_U_d_(0+nl*n_svd_d),1
     $           ,svd_d(0),svd_polyval_U_d_(nl+ndv*n_svd_l))
            if (verbose.gt.1) then
               write (6,'(A,F6.3,A,F6.3,A,I0,A,F8.5)') ' % delta ' ,
     $              delta , ' svd_d ' ,svd_d(0) , ' U_d_(' , nl , ') ',
     $              svd_polyval_U_d_(nl+ndv*n_svd_l)
            end if !if (verbose.gt.1) then
         enddo !do nl=0,n_svd_l-1
      enddo !do ndv=0,n_delta_v-1
      if (verbose.gt.0) then
         write(6,'(A)') '[finished get_svd_polyval_U_d_]'
      end if !if (verbose.gt.0) then
      end
