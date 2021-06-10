      subroutine test_innerproduct_batch_excerpt_8(n_delta_x,delta_x_
     $     ,n_delta_y,delta_y_,n_gamma_z,gamma_z_,Z_p_,ndx_optimal
     $     ,ndy_optimal ,ngz_optimal ,C_Z_optimal)
      implicit none
      integer verbose
      data verbose / 0 /
      integer n_delta_x,n_delta_y,n_gamma_z,ndx_optimal,ndy_optimal
     $     ,ngz_optimal
      real *8 delta_x_(0:n_delta_x-1),delta_y_(0:n_delta_y-1)
      real *8 gamma_z_(0:n_gamma_z-1)
      complex *16 Z_p_(0:n_delta_x*n_delta_y*n_gamma_z-1)
      complex *16 C_Z_optimal
      integer ndx,ndy,ngz,nC
      if (verbose.gt.0) then
         write(6,'(A)') ' [entering test_innerproduct_batch_excerpt_7] '
      end if !if (verbose.gt.0) then
      if (verbose.gt.0) then
         write(6,'(A,I0)') ' n_delta_x: ' , n_delta_x
         write(6,'(A,I0)') ' n_delta_y: ' , n_delta_y
         write(6,'(A,I0)') ' n_gamma_z: ' , n_gamma_z
      end if !if (verbose.gt.0) then
      ndx_optimal = 0
      ndy_optimal = 0
      ngz_optimal = 0
      C_Z_optimal = (-1.0d0,0.0d0)
      nC=0
      do ngz=0,n_gamma_z-1
         do ndy=0,n_delta_y-1
            do ndx=0,n_delta_x-1
               nC = ndx + ndy*n_delta_x + ngz*n_delta_x*n_delta_y
               if (real(Z_p_(nC)).gt.real(C_Z_optimal)) then
                  C_Z_optimal = Z_p_(nC)
                  ndx_optimal = ndx
                  ndy_optimal = ndy
                  ngz_optimal = ngz
               end if           !if new value is bigger
            enddo ! n_delta_x
         enddo ! n_delta_y
      enddo ! n_gamma_z
      if (verbose.gt.0) then
         write(6,'(A)') ' [finished test_innerproduct_batch_excerpt_7] '
      end if !if (verbose.gt.0) then
      end
