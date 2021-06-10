      subroutine test_innerproduct_batch_excerpt_4(delta_x_est
     $    ,delta_y_est,delta_x_,delta_y_,displacement_max,n_delta_x
     $    ,n_delta_y,n_gamma_z,Z_q_,ndx_max,ndy_max,ngz_max,C_Z_max)
      implicit none
      real *8 delta_x_est,delta_y_est,displacement_max
      real *8 delta_x_(0:n_delta_x-1),delta_y_(0:n_delta_y-1)
      integer n_delta_x,n_delta_y,n_gamma_z,ndx_max,ndy_max,ngz_max
      complex *16 Z_q_(0:n_delta_x*n_delta_y*n_gamma_z-1),C_Z_max
      integer ndx,ndy,ngz,na
      complex *16 C_Z
      logical found_flag
      real *8 delta_tmp
      found_flag = .false.
      ndx_max = 0
      ndy_max = 0
      ngz_max = 0
      na=0
      do ngz=0,n_gamma_z-1
         do ndy=0,n_delta_y-1
            do ndx=0,n_delta_x-1
               delta_tmp = dsqrt((delta_x_est + delta_x_(ndx))**2 +
     $              (delta_y_est + delta_y_(ndy))**2)
               if (delta_tmp.le.displacement_max) then
                  if (found_flag.eqv..false.) then
                     C_Z = Z_q_(na)
                     C_Z_max = C_Z
                     ndx_max = ndx
                     ndy_max = ndy
                     ngz_max = ngz
                     found_flag = .true.
                  else !if (found_flag.eqv..true.) then
                     C_Z = Z_q_(na)
                     if (real(C_Z).gt.real(C_Z_max)) then
                        C_Z_max = C_Z
                        ndx_max = ndx
                        ndy_max = ndy
                        ngz_max = ngz
                     end if !if new value is bigger
                  end if !found_flag
               end if ! (delta_tmp.le.displacement_max)
               na = na+1
            enddo ! n_delta_x
         enddo ! n_delta_y
      enddo ! n_gamma_z
      if (found_flag.eqv..false.) then
         write(6,'(A)') 'Warning, found_flag.eqv..false. '
     $        ,'in test_innerproduct_batch_excerpt_4.f'
      end if
      end
