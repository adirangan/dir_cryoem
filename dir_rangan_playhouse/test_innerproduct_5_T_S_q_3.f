      subroutine test_innerproduct_5_T_S_q_3(verbose,n_delta_x
     $     ,delta_x_,n_delta_y,delta_y_,fftw_plan_frwd_,fftw_in1_
     $     ,fftw_out_ ,n_r,grid_p_,n_w_,n_A,S_p_,S_q_,n_S_sub
     $     ,I_S_sample_,ld_S,S_p__ ,n_S_tot ,T_S_q__)
c$$$      %%%%%%%%%%%%%%%%%%%%%%%
c$$$      Calculate bessel coefficients of T(S).
c$$$      T = translation by delta_x,delta_y.
c$$$      S = template in k-space polar coord.
c$$$      T_S_q__(nr + n_r*(ndx + n_delta_x*(ndy + n_delta_y*(ns + n_S_tot*nw))))
c$$$      is equal to the bessel-coefficients
c$$$      conjg(S_q_(ic))
c$$$      for ic = nw + n_w_sum_(nr),
c$$$      where S_q_ = T(S), with
c$$$      T <-- delta_x_(ndx) and delta_y_(ndy)
c$$$      and S <-- S_p__(ns*ld_S).
c$$$      %%%%%%%%%%%%%%%%%%%%%%%
      implicit none
      include '/usr/include/fftw3.f'
      integer verbose
      integer n_delta_x,n_delta_y
      integer n_r,n_w_(0:n_r-1),n_A,n_S_sub,I_S_sample_(0:0)
     $     ,n_S_tot,ld_S
      real *8 delta_x_(0:n_delta_x-1)
      real *8 delta_y_(0:n_delta_y-1)
      real *8 grid_p_(0:n_r-1)
      integer *8 fftw_plan_frwd_(0:n_r-1)
      complex *16 fftw_in1_(0:0),fftw_out_(0:0)
      complex *16 Z_q
      complex *16 S_p_(0:0),S_q_(0:0),S_p__(0:0)
      complex *16 T_S_q__(0:0)
      complex *16 Z_tmp
      integer n_w_max,nr,nw,ns
      integer ndx,ndy,nx,ld_X
      logical flag_conjg
      real *8 delta_x,delta_y
      real *8 pi
      if (verbose.gt.0) then
         write(6,'(A)') '[entering test_innerproduct_5_T_S_q_3]'
      end if !if (verbose.gt.0) then
      pi = 4.0*atan(1.0)
      n_w_max = n_w_(n_r-1)

      do ns=0,n_S_sub-1
         if (verbose.gt.2) then
            call innerproduct_p(n_r,grid_p_,n_w_,n_A
     $           ,S_p__(I_S_sample_(ns)*ld_S)
     $           ,S_p__(I_S_sample_(ns)*ld_S),Z_tmp)
            Z_tmp = zsqrt(Z_tmp)/(n_r*n_r)
            write(6,'(A,I0,A,F8.4,1X,F8.4)') '|S_p_(',ns,')|: ' , Z_tmp
            call cp1_c16(n_A,S_p__(I_S_sample_(ns)*ld_S),S_p_)
            call interp_p_to_q_fftw(n_r,fftw_plan_frwd_,n_w_ ,n_A
     $           ,fftw_in1_,fftw_out_,S_p_,S_q_)
            call write_sub_c16(n_A,S_q_,5,'S_q_: ')
            call innerproduct_p(n_r,grid_p_,n_w_,n_A,S_q_,S_q_,Z_tmp)
            Z_tmp = zsqrt(Z_tmp)/(n_r*n_r)
            write(6,'(A,I0,A,F8.4,1X,F8.4)') '|S_q_(',ns,')|: ' , Z_tmp
         end if                 !if (verbose.gt.2) then
         do ndx=0,n_delta_x-1
            delta_x = delta_x_(ndx)
            do ndy=0,n_delta_y-1
               delta_y = delta_y_(ndy)
               call cp1_c16(n_A,S_p__(I_S_sample_(ns)*ld_S),S_p_)
               call transf_p_to_p(n_r,grid_p_,n_w_
     $              ,n_A,S_p_,+delta_x,+delta_y,S_p_)
               call interp_p_to_q_fftw(n_r,fftw_plan_frwd_,n_w_
     $              ,n_A,fftw_in1_,fftw_out_,S_p_,S_q_)
               flag_conjg = .true.
               nx = 0 + n_r*(ndx + n_delta_x*(ndy + n_delta_y*ns))
               ld_X = n_r*n_delta_x*n_delta_y*n_S_tot
               call innerproduct_q_stretch(n_r,grid_p_,n_w_,n_A,S_q_
     $              ,flag_conjg,ld_X,T_S_q__(nx))
               if (verbose.gt.2) then
                  write(6,'(A,I0,A,I0,A,I0)') ' ns: ' , ns , ' ndx: ' ,
     $                 ndx , ' ndy: ' , ndy
                  call write_sub_c16(n_A,S_q_,6,'S_q_: ')
                  call innerproduct_p(n_r,grid_p_,n_w_,n_A,S_q_,S_q_
     $                 ,Z_tmp)
                  Z_tmp = zsqrt(Z_tmp)/(n_r*n_r)
                  write(6,'(A,F8.4,1X,F8.4)') '|S_q_|: ' , Z_tmp
               end if ! if (verbose.gt.2) then
            enddo !do ndy=0,n_delta_y-1
         enddo ! do ndx=0,n_delta_x-1
      enddo ! do ns=0,n_S_sub-1

      if (verbose.gt.0) then
         write(6,'(A)') '[finished test_innerproduct_5_T_S_q_3]'
      end if !if (verbose.gt.0) then
      end
