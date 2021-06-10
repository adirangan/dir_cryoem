!> Doxygen comment: ;\n
!> test_innerproduct_8 excerpt: ;\n
!> Calculates :  ;\n
!> CTF_R_S__(ngz + n_gamma_z*(nctf + n_CTF*ns)) = Z_q_(ngz), ;\n
!> where: ;\n
!> Z_q_(ngz) = \| CTF_p_(nctf) .* R(ngz) ( S_p__(ns*ld_S) ) \|_{L^{2}} ,  ;\n
!> where: ;\n
!> ".*" is a pointwise product,  ;\n
!> CTF_p_ is a CTF-function given on nonuniform polar grid, ;\n
!> S_p_ is a template given on a nonuniform polar grid, ;\n
!> and R is rotation by +gamma_z_(ngz). ;\n
!> (see test_transforms_dr.f for examples). ;\n
!> This calculation is carried out by considering ;\n
!> the bessel-function expansions of:  ;\n
!> dconjg(CTF_p_).*(CTF_p_) ;\n
!> and ;\n
!> dconjg(S_p_).*(S_p_), ;\n
!> and appealing to the fft. ;\n
!> Updated to use radial quadrature weight. ;
      subroutine ti8_build_CTF_R_S_quad_3(
     $     verbose
     $     ,n_gamma_z 
     $     ,gamma_z_
     $     ,fftw_plan_frwd_
     $     ,fftw_plan_back_
     $     ,fftw_0in_
     $     ,fftw_out_ 
     $     ,n_r
     $     ,grid_p_
     $     ,weight_p_
     $     ,n_w_
     $     ,n_w_sum 
     $     ,S_p_
     $     ,S_q_
     $     ,n_S
     $     ,I_S_sample_
     $     ,ld_S 
     $     ,S_p__ 
     $     ,CTF_p_ 
     $     ,CTF_q_ 
     $     ,n_CTF 
     $     ,ld_CTF
     $     ,CTF_p__
     $     ,CTF_R_S__ 
     $     ,Z_q_
     $     ,Z_p_
     $     )
c$$$      %%%%%%%%%%%%%%%%%%%%%%%
c$$$      Calculates : 
c$$$      CTF_R_S__(ngz + n_gamma_z*(nctf + n_CTF*ns)) = Z_q_(ngz),
c$$$      where:
c$$$      Z_q_(ngz) = \| CTF_p_(nctf) .* R(ngz) ( S_p__(I_S_sample_(ns)*ld_S) ) \|_{L^{2}} , 
c$$$      where:
c$$$      ".*" is a pointwise product, 
c$$$      CTF_p_ is a CTF-function given on nonuniform polar grid,
c$$$      S_p_ is a template given on a nonuniform polar grid,
c$$$      and R is rotation by +gamma_z_(ngz).
c$$$      (see test_transforms_dr.f for examples).
c$$$      This calculation is carried out by considering
c$$$      the bessel-function expansions of: 
c$$$      dconjg(CTF_p_).*(CTF_p_)
c$$$      and
c$$$      dconjg(S_p_).*(S_p_),
c$$$      and appealing to the fft.
c$$$      Updated to use radial quadrature weight. ;
c$$$      In this function weight_p_ refers to the radial weights only (size n_r), ;
c$$$      %%%%%%%%%%%%%%%%%%%%%%%
      implicit none
      include '/usr/include/fftw3.f'
      integer verbose
      integer n_gamma_z
      integer n_r,n_w_(0:n_r-1),n_w_sum
      integer n_S,I_S_sample_(0:n_S-1),ld_S,n_CTF,ld_CTF
      real *8 gamma_z_(0:n_gamma_z-1),grid_p_(0:n_r-1)
      real *8 weight_p_(0:n_r-1)
      integer *8 fftw_plan_frwd_(0:n_r-1)
      integer *8 fftw_plan_back_(0:n_r-1)
      complex *16 fftw_0in_(0:0),fftw_out_(0:0)
      pointer (p_fftw_plan_frwd_last,fftw_plan_frwd_last)
      pointer (p_fftw_plan_back_last,fftw_plan_back_last)
      pointer (p_fftw_0in_last_,fftw_0in_last_)
      pointer (p_fftw_out_last_,fftw_out_last_)
      integer *8 fftw_plan_frwd_last
      integer *8 fftw_plan_back_last
      complex *16 fftw_0in_last_(0:0),fftw_out_last_(0:0)
      complex *16 Z_p_(0:0)
      complex *16 Z_q_(0:0)
      complex *16 Z_q
      complex *16 S_p_(0:0),S_q_(0:0),S_p__(0:0)
      complex *16 CTF_p_(0:0),CTF_q_(0:0),CTF_p__(0:0)
      complex *16 CTF_R_S__(0:0)
      integer n_w_max,nr,ngz,nw,ns,nctf
      real *8 tmp_gamma_z,gamma_z
      logical flag_nufft
      real *8 pi
      real *8 dotw_c16_f !function output. ;
      if (verbose.gt.0) then
         write(6,'(A)') '[entering ti8_build_CTF_R_S_quad_3]'
      end if !if (verbose.gt.0) then

      if (n_w_sum.lt.n_gamma_z) then
         if (verbose.gt.1) then
            write(6,'(2(A,I0))') ' Note: n_w_sum: ' , n_w_sum ,
     $           ' n_gamma_z: ' ,n_gamma_z
         end if !if (verbose.gt.1) then
      end if !if (n_w_sum.lt.n_gamma_z) then
      call cl1_c16(n_w_sum,S_p_)
      call cl1_c16(n_w_sum,S_q_)
      call cl1_c16(n_gamma_z,Z_q_)
      if (verbose.gt.1) then
         write(6,'(A,I0)') ' verbose: ' , verbose
         write(6,'(A,I0)') ' n_gamma_z: ' , n_gamma_z
         call print_sub_r8(n_gamma_z,gamma_z_,' gamma_z_: ')
         call print_sub_c16(n_w_sum,fftw_0in_,' fftw_0in_: ')
         call print_sub_c16(n_w_sum,fftw_out_,' fftw_out_: ')
         write(6,'(A,I0)') ' n_r: ' , n_r
         call print_sub_r8(n_r,grid_p_,' grid_p_: ')
         call print_sub_r8(n_r,weight_p_,' weight_p_: ')
         call print_sub_i4(n_r,n_w_,' n_w_: ')
         write(6,'(A,I0)') ' n_w_sum: ' , n_w_sum
         call print_sub_c16(n_w_sum,S_p_,' S_p_: ')
         call print_sub_c16(n_w_sum,S_q_,' S_q_: ')
         write(6,'(A,I0)') ' n_S: ' , n_S
         call print_sub_i4(n_S,I_S_sample_,' I_S_sample_: ')
         write(6,'(A,I0)') ' ld_S: ' , ld_S
         call print_sub_c16(ld_S,S_p__(I_S_sample_(n_S-1)*ld_S)
     $        ,' S_p__: ')
         write(6,'(A,I0)') ' n_CTF: ' , n_CTF
         write(6,'(A,I0)') ' ld_CTF: ' , ld_CTF
         call print_sub_c16(n_CTF*ld_CTF,CTF_p__,' CTF_p__: ')
         call print_sub_c16(n_gamma_z*n_CTF*n_S,CTF_R_S__,'CTF_R_S__: ')
         call print_sub_c16(n_gamma_z,Z_q_,' Z_q_: ')
      end if !if (verbose.gt.1) then

      pi = 4.0d0*datan(1.0d0)
      n_w_max = n_w_(n_r-1)
      if (n_w_sum.lt.n_w_max) then
         if (verbose.gt.1) then
            write(6,'(2(A,I0))') ' Note: n_w_sum: ' , n_w_sum ,
     $           ' n_w_max: ' ,n_w_max
         end if !if (verbose.gt.1) then
      end if !if (n_w_sum.lt.n_w_max) then
      flag_nufft = .false.
      if (n_gamma_z.ne.n_w_max) then
         flag_nufft = .true.
      end if !if (n_gamma_z.ne.n_w_max) then
      if (n_gamma_z.eq.n_w_max) then
         do nw=0,n_w_max-1
            tmp_gamma_z = (2.0d0*pi*nw)/(1.0d0*n_gamma_z)
            ngz = nw
            gamma_z = gamma_z_(ngz)
            if (dabs(tmp_gamma_z-gamma_z).gt.1.0d-12) then
               flag_nufft = .true.
            end if !if (dabs(tmp_gamma_z-gamma_z).gt.1.0d-12) then
         enddo !do nw=0,n_w_max-1
      end if !if (n_gamma_z.eq.n_w_max) then
      if (verbose.gt.1) then
         write(6,'(A,I0)') ' n_w_max: ' , n_w_max
         write(6,'(A,I0)') ' n_gamma_z: ' , n_gamma_z
         write(6,'(A,L2)') ' flag_nufft: ' , flag_nufft
         write(6,'(A)') ' linking fftw_plan_frwd_last, etc. '
      end if !if (verbose.gt.1) then      
      p_fftw_plan_frwd_last = loc(fftw_plan_frwd_(n_r-1))
      p_fftw_plan_back_last = loc(fftw_plan_back_(n_r-1))
      p_fftw_0in_last_ = loc(fftw_0in_(n_w_sum-n_w_max))
      p_fftw_out_last_ = loc(fftw_out_(n_w_sum-n_w_max))
      if (verbose.gt.1) then
         write(6,'(A)') ' clearing Z_q_. '
      end if !if (verbose.gt.1) then      
      call cl1_c16(n_w_max,Z_q_)
      if (flag_nufft.eqv..true.) then
         call cs1_c16(n_gamma_z,Z_p_)
      end if !if (flag_nufft.eqv..true.) then
      do ns=0,n_S-1
         if (verbose.gt.2) then
            write(6,'(A,I0)') ' ns: ' , ns
         end if !if (verbose.gt.2) then      
         call cp1_c16(n_w_sum,S_p__(I_S_sample_(ns)*ld_S),S_p_)
         if (verbose.gt.2) then
            call innerproduct_p_quad(n_r,grid_p_,weight_p_,n_w_,n_w_sum
     $           ,S_p__(I_S_sample_(ns)*ld_S),S_p__(I_S_sample_(ns)
     $           *ld_S),Z_q)
            Z_q = zsqrt(Z_q)
            write(6,'(A,I0,A,F8.4,1X,F8.4)') ' |S_p__(',ns,'*ld_S)|: ' ,
     $           Z_q
         end if                 !if (verbose.gt.2) then
         do nctf=0,n_CTF-1
            if (verbose.gt.2) then
               write(6,'(A,I0)') ' nctf: ' , nctf
            end if !if (verbose.gt.2) then      
            call cp1_c16(n_w_sum,CTF_p__(nctf*ld_CTF),CTF_p_)
            call xc1_c16(n_w_sum,S_p_,S_p_,S_q_)
            call interp_p_to_q_fftw(n_r,fftw_plan_frwd_,n_w_
     $           ,n_w_sum,fftw_0in_,fftw_out_,S_q_,S_q_)
            call xc1_c16(n_w_sum,CTF_p_,CTF_p_,CTF_q_)
            call interp_p_to_q_fftw(n_r,fftw_plan_frwd_,n_w_
     $           ,n_w_sum,fftw_0in_,fftw_out_,CTF_q_,CTF_q_)
            call innerproduct_q_k_stretch_quad_0(n_r,grid_p_,weight_p_
     $           ,n_w_,n_w_sum,S_q_,CTF_q_,Z_q_)
            call cp1_c16(n_w_max,Z_q_,fftw_out_last_)
            call dfftw_execute_(fftw_plan_back_last)
            call cp1_c16(n_w_max,fftw_0in_last_,Z_q_)
            if (flag_nufft) then
               if (verbose.gt.2) then
                  write(6,'(A)') ' calling interp1_finufft_c16'
               end if !if (verbose.gt.2) then      
c$$$               Unfortunately, I could not get either finufft1d2 or nufft1d2 to work across multiple threads. ;
c$$$               call interp1_finufft_c16(
c$$$     $              n_w_max
c$$$     $              ,0.0d0
c$$$     $              ,2.0d0*pi
c$$$     $              ,Z_q_
c$$$     $              ,fftw_plan_frwd_last
c$$$     $              ,fftw_plan_back_last
c$$$     $              ,fftw_0in_last_
c$$$     $              ,fftw_out_last_
c$$$     $              ,n_gamma_z
c$$$     $              ,gamma_z_
c$$$     $              ,Z_p_
c$$$     $              )
               call interp1_lagrange_c16(
     $              n_w_max
     $              ,0.0d0
     $              ,2.0d0*pi
     $              ,Z_q_
     $              ,n_gamma_z
     $              ,gamma_z_
     $              ,Z_p_
     $              )
               call cp1_c16(n_gamma_z,Z_p_,Z_q_)
            end if !if (flag_nufft) then
            do ngz=0,n_gamma_z-1
               Z_q = Z_q_(ngz)
               CTF_R_S__(ngz + n_gamma_z*(nctf + n_CTF*ns)) = zsqrt(Z_q)
            enddo !do ngz=0,n_gamma_z-1
            if (verbose.gt.2) then
               call print_sub_r8(n_gamma_z,gamma_z_,' gamma_z_: ')
               call print_sub_c16(n_gamma_z,CTF_R_S__(0+n_gamma_z*(nctf
     $              + n_CTF*ns)),' CTF_R_S__: ')
            end if !if (verbose.gt.2) then
         enddo !do nctf=0,n_CTF-1
      enddo !do ns=0,n_S-1

      if (verbose.gt.0) then
         write(6,'(A)') '[finished ti8_build_CTF_R_S_quad_3]'
      end if !if (verbose.gt.0) then
      end
