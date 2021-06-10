      subroutine test_innerproduct_fast_CTF_p_R_S_p_0(n_gamma_z,gamma_z_
     $     ,fftw_plan_frwd_,fftw_in1_,fftw_out_,fftw_plan_back_last
     $     ,fftw_in1_last_,fftw_out_last_,n_r,grid_p_,n_w_,n_A,S_p_
     $     ,CTF_p_,Z_q_)
c$$$      %%%%%%%%%%%%%%%%%%%%%%%
c$$$      Calculates : 
c$$$      Z_q_(ngz) = \| CTF_p_ .* R ( S_p_ ) \|_{L^{2}} , 
c$$$      where:
c$$$      ".*" is a pointwise product, 
c$$$      CTF_p_ is a CTF-function given on nonuniform polar grid,
c$$$      S_p_ is a template given on a nonuniform polar grid,
c$$$      and R is rotation by +gamma_z_(ngz).
c$$$      (see test_innerproduct_timing_dr.f for examples).
c$$$      This calculation is carried out by considering
c$$$      the bessel-function expansions of: 
c$$$      conjg(CTF_p_).*(CTF_p_)
c$$$      and
c$$$      conjg(S_p_).*(S_p_),
c$$$      and appealing to the fft.
c$$$      %%%%%%%%%%%%%%%%%%%%%%%
      implicit none
      include '/usr/include/fftw3.f'
      integer n_gamma_z
      integer n_r,n_w_(0:n_r-1),n_A
      real *8 gamma_z_(0:n_gamma_z-1),grid_p_(0:n_r-1)
      integer *8 fftw_plan_frwd_(0:n_r-1),fftw_plan_back_last
      complex *16 fftw_in1_(0:0),fftw_out_(0:0)
      complex *16 fftw_in1_last_(0:0),fftw_out_last_(0:0)
      complex *16 Z_q
      complex *16 S_p_(0:0),CTF_p_(0:0),Z_q_(0:0)
      complex *16, allocatable :: S_q_(:)
      complex *16, allocatable :: CTF_q_(:)
      complex *16, allocatable :: Z_tmp_(:)
      integer n_w_max,nr,ngz,nw,l
      real *8 gamma_z
      real *8 pi
      pi = 4.0*atan(1.0)
      n_w_max = n_w_(n_r-1)
      allocate(S_q_(0:n_A-1))
      allocate(CTF_q_(0:n_A-1))
      allocate(Z_tmp_(0:n_w_max-1))
      call xc1_c16(n_A,S_p_,S_p_,S_q_)
      call interp_p_to_q_fftw(n_r,fftw_plan_frwd_,n_w_
     $     ,n_A,fftw_in1_,fftw_out_,S_q_,S_q_)
      call xc1_c16(n_A,CTF_p_,CTF_p_,CTF_q_)
      call interp_p_to_q_fftw(n_r,fftw_plan_frwd_,n_w_
     $     ,n_A,fftw_in1_,fftw_out_,CTF_q_,CTF_q_)
      call innerproduct_q_k_only(n_r,grid_p_,n_w_
     $     ,n_A,S_q_,CTF_q_,Z_tmp_)
      call cp1_c16(n_w_max,Z_tmp_,fftw_out_last_)
      call dfftw_execute_(fftw_plan_back_last)
      call cp1_c16(n_w_max,fftw_in1_last_,Z_tmp_)
      do ngz=0,n_gamma_z-1
         gamma_z = gamma_z_(ngz)
         call interp1_c16(n_w_max,0.0d0,2*pi,Z_tmp_,+gamma_z
     $        ,Z_q)
         Z_q_(ngz) = zsqrt(Z_q)/(n_r**2)
      enddo
      deallocate(Z_tmp_)
      deallocate(CTF_q_)
      deallocate(S_q_)
      end
