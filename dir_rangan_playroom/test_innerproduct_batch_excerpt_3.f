      subroutine test_innerproduct_batch_excerpt_3(n_delta_x,n_delta_y
     $     ,n_gamma_z,Z_q_,ndx_max,ndy_max,ngz_max,C_Z_max)
      implicit none
      integer n_delta_x,n_delta_y,n_gamma_z,ndx_max,ndy_max,ngz_max
      complex *16 Z_q_(0:n_delta_x*n_delta_y*n_gamma_z-1),C_Z_max
      integer ndx,ndy,ngz,na
      complex *16 C_Z
      C_Z_max = Z_q_(0)
      ndx_max = 0
      ndy_max = 0
      ngz_max = 0
      na=0
      do ngz=0,n_gamma_z-1
         do ndy=0,n_delta_y-1
            do ndx=0,n_delta_x-1
               C_Z = Z_q_(na)
               if (real(C_Z).gt.real(C_Z_max)) then
                  C_Z_max = C_Z
                  ndx_max = ndx
                  ndy_max = ndy
                  ngz_max = ngz
               end if
               na = na+1
            enddo
         enddo
      enddo
      end
