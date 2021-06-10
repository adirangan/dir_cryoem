      subroutine test_innerproduct_batch_excerpt_5(delta_x_est
     $     ,delta_y_est,gamma_z_est,delta_x_,delta_y_,displacement_max
     $     ,n_delta_x ,n_delta_y,n_gamma_z,Z_q_,ndx_optimal,ndy_optimal
     $     ,ngz_optimal ,C_Z_optimal)
      implicit none
      integer verbose
      data verbose / 0 /
      real *8 delta_x_est,delta_y_est,gamma_z_est,displacement_max
      real *8 delta_x_(0:n_delta_x-1),delta_y_(0:n_delta_y-1)
      integer n_delta_x,n_delta_y,n_gamma_z,ndx_optimal,ndy_optimal
     $     ,ngz_optimal
      complex *16 Z_q_(0:n_delta_x*n_delta_y*n_gamma_z-1)
      complex *16 C_Z_optimal
      real *8 delta_x_tmp,delta_y_tmp
      real *8 delta_x_sort,delta_y_sort
      integer ndx,ndy,ngz,na
      complex *16 C_Z
      logical found_flag
      real *8 delta_tmp
      if (verbose.gt.0) then
         write(6,'(A)') ' [entering test_innerproduct_batch_excerpt_5] '
      end if !if (verbose.gt.0) then
      if (verbose.gt.0) then
         write(6,'(A,I0)') ' n_delta_x: ' , n_delta_x
         write(6,'(A,I0)') ' n_delta_y: ' , n_delta_y
         write(6,'(A,I0)') ' n_gamma_z: ' , n_gamma_z
         write(6,'(A,F8.4)') ' displacement_max: ' , displacement_max
      end if !if (verbose.gt.0) then
      found_flag = .false.
      ndx_optimal = 0
      ndy_optimal = 0
      ngz_optimal = 0
      na=0
      do ngz=0,n_gamma_z-1
         do ndy=0,n_delta_y-1
            do ndx=0,n_delta_x-1
               delta_x_sort = delta_x_(ndx)
               delta_y_sort = delta_y_(ndy)
               call get_interchange_delta(delta_x_sort,delta_y_sort
     $              ,gamma_z_est,delta_x_tmp,delta_y_tmp)
               delta_tmp = dsqrt((delta_x_est + delta_x_tmp)**2 +
     $              (delta_y_est + delta_y_tmp)**2)
               if (verbose.gt.1) then
                  write(6,'(A,F8.4)') ' delta_x_sort ' , delta_x_sort 
                  write(6,'(A,F8.4)') ' delta_y_sort ' , delta_y_sort 
                  write(6,'(A,F8.4)') ' delta_x_tmp ' , delta_x_tmp 
                  write(6,'(A,F8.4)') ' delta_y_tmp ' , delta_y_tmp 
                  write(6,'(A,F8.4)') ' delta_x_est ' , delta_x_est 
                  write(6,'(A,F8.4)') ' delta_y_est ' , delta_y_est 
                  write(6,'(A,F8.4)') ' delta_tmp ' , delta_tmp 
               end if !if (verbose.gt.1) then
               if (delta_tmp.le.displacement_max) then
                  if (found_flag.eqv..false.) then
                     C_Z = Z_q_(na)
                     C_Z_optimal = C_Z
                     ndx_optimal = ndx
                     ndy_optimal = ndy
                     ngz_optimal = ngz
                     found_flag = .true.
                  else !if (found_flag.eqv..true.) then
                     C_Z = Z_q_(na)
                     if (real(C_Z).gt.real(C_Z_optimal)) then
                        C_Z_optimal = C_Z
                        ndx_optimal = ndx
                        ndy_optimal = ndy
                        ngz_optimal = ngz
                     end if !if new value is bigger
                  end if !found_flag
               end if ! (delta_tmp.le.displacement_max)
               na = na+1
            enddo ! n_delta_x
         enddo ! n_delta_y
      enddo ! n_gamma_z
      if (verbose.gt.0) then
         write(6,'(A,I0)') ' ndx_optimal: ' , ndx_optimal
         write(6,'(A,I0)') ' ndy_optimal: ' , ndy_optimal
         write(6,'(A,I0)') ' ngz_optimal: ' , ngz_optimal
         write(6,'(A,2F8.4)') ' C_Z_optimal: ' , C_Z_optimal
      end if !if (verbose.gt.0) then
      if (found_flag.eqv..false.) then
         write(6,'(A)') 'Warning, found_flag.eqv..false. '
     $        ,'in test_innerproduct_batch_excerpt_5.f'
      end if !if (found_flag.eqv..false.) then
      if (verbose.gt.0) then
         write(6,'(A)') ' [finished test_innerproduct_batch_excerpt_5] '
      end if !if (verbose.gt.0) then
      end
