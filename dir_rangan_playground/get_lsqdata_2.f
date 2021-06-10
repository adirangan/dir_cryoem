!> Doxygen comment: ;\n
!> set up right-hand-side for least-squares problem ;\n
      subroutine get_lsqdata_2(
     $     n_M
     $     ,I_M_sample_
     $     ,ld_M
     $     ,M_k_p__
     $     ,n_CTF
     $     ,ld_CTF
     $     ,CTF_k_p_
     $     ,k_p_r
     $     ,fftw_plan_frwd
     $     ,fftw_plan_back
     $     ,n_w_csum
     $     ,n_w
     $     ,fftw_0in_
     $     ,fftw_out_
     $     ,alpha2d__
     $     ,M_k_c_
     $     ,polar_a_
     $     ,azimu_b_
     $     ,weight_CTF_k_c_
     $     ,n_w_M
     $     )
      implicit none
      integer verbose
      data verbose / 0 /
      integer n_M,I_M_sample_(0:n_M-1),ld_M,n_CTF,ld_CTF,n_w_csum,n_w
     $     ,n_w_M
      integer *8 fftw_plan_frwd
      integer *8 fftw_plan_back
      complex *16 fftw_0in_(0:n_w-1),fftw_out_(0:n_w-1)
      include 'excerpt_define_nalpha.f'
      complex *16 M_k_p__(0:0)
      complex *16 CTF_k_p_(0:ld_CTF*n_CTF-1)
      real *8 k_p_r
      real *8 alpha2d__(0:n_alpha-1,0:n_M-1)
      complex *16 M_k_c_(0:n_w*n_M-1)
      real *8 polar_a_(0:n_w*n_M-1)
      real *8 azimu_b_(0:n_w*n_M-1)
      complex *16 weight_CTF_k_c_(0:n_w*n_M-1)
      real *8 pi,cos_polar_a,azimu_b
      real *8 k_c_0,k_c_1,k_c_2,k_p_r_1
      real *8 delta_x,delta_y,gamma_z,l2_norm
      real *8, allocatable :: k_c_0_(:)
      real *8, allocatable :: k_c_1_(:)
      real *8, allocatable :: k_c_2_(:)
      complex *16, allocatable ::  M_k_p_single_(:)
      complex *16, allocatable ::  CTF_k_p_single_(:)
      integer nM,nw,nctf
      logical flag_memory_checkset

c$$$         %%%%%%%%
      if (verbose.gt.0) then 
         write(6,'(A)') '[entering get_lsqdata_2]'
      end if !if (verbose.gt.0) then 
      allocate(k_c_0_(0:1+n_w-1))
      call cs1_r8(n_w,k_c_0_)
      allocate(k_c_1_(0:1+n_w-1))
      call cs1_r8(n_w,k_c_1_)
      allocate(k_c_2_(0:1+n_w-1))
      call cs1_r8(n_w,k_c_2_)
      allocate(M_k_p_single_(0:1+n_w-1))
      call cs1_c16(n_w,M_k_p_single_)
      allocate(CTF_k_p_single_(0:1+n_w-1))
      call cs1_c16(n_w,CTF_k_p_single_)
      pi=4.0d0*datan(1.0d0)
c$$$         %%%%%%%%

c$$$         %%%%%%%%
      n_w_M = 0
      do nM = 0,n_M-1
c$$$         %%%%%%%%
         nctf = nint(alpha2d__(nalpha_ctf_ind,nM))
         nctf = max(0,min(n_CTF-1,nctf))
         cos_polar_a = dcos(alpha2d__(nalpha_polar_a,nM))
         azimu_b = alpha2d__(nalpha_azimu_b,nM)
         gamma_z = 0.0d0
         k_p_r_1 = 1.0d0
         call mkonegreatcircle(
     $        k_p_r_1
     $        ,cos_polar_a
     $        ,azimu_b
     $        ,gamma_z
     $        ,n_w
     $        ,k_c_0_
     $        ,k_c_1_
     $        ,k_c_2_
     $        )
c$$$         %%%%%%%%
         if (n_alpha.ge.1+nalpha_l2_norm) then
            l2_norm = alpha2d__(nalpha_l2_norm,nM)
         else
            l2_norm = 1.0d0
         end if            
         if (dabs(l2_norm).lt.1.0d-15) then
            l2_norm=1.0d0
         end if
         do nw = 0,n_w-1
            M_k_p_single_(nw) = M_k_p__(n_w_csum-1+nw + I_M_sample_(nM)
     $           *ld_M)/l2_norm
            CTF_k_p_single_(nw) = CTF_k_p_(n_w_csum-1+nw + nctf*ld_CTF)
         enddo !do nw = 0,n_w-1
c$$$         %%%%%%%%
         gamma_z = -alpha2d__(nalpha_gamma_z,nM)
         call rotate_p_to_p_fftw_single(
     $        fftw_plan_frwd
     $        ,fftw_plan_back
     $        ,n_w
     $        ,fftw_0in_
     $        ,fftw_out_
     $        ,M_k_p_single_
     $        ,+gamma_z
     $        ,M_k_p_single_
     $        )
         call rotate_p_to_p_fftw_single(
     $        fftw_plan_frwd
     $        ,fftw_plan_back
     $        ,n_w
     $        ,fftw_0in_
     $        ,fftw_out_
     $        ,CTF_k_p_single_
     $        ,+gamma_z
     $        ,CTF_k_p_single_
     $        )
         delta_x = -alpha2d__(nalpha_delta_x,nM)
         delta_y = -alpha2d__(nalpha_delta_y,nM)
         call transf_p_to_p_single(
     $        k_p_r
     $        ,n_w
     $        ,M_k_p_single_
     $        ,+delta_x
     $        ,+delta_y
     $        ,M_k_p_single_
     $        )
c$$$         %%%%%%%%
         do nw = 0,n_w-1
            k_c_0 = k_c_0_(nw)
            k_c_1 = k_c_1_(nw)
            k_c_2 = k_c_2_(nw)
            cos_polar_a = k_c_2
            azimu_b = datan2(k_c_1,k_c_0)
            if (azimu_b.lt.0.0d0) azimu_b = azimu_b + 2.0d0*pi
            azimu_b_(n_w_M) = azimu_b
            polar_a_(n_w_M) = datan2(dsqrt(k_c_0**2+ k_c_1**2)
     $           ,cos_polar_a)
            M_k_c_(n_w_M) = M_k_p_single_(nw)
            weight_CTF_k_c_(n_w_M) = CTF_k_p_single_(nw)
            n_w_M = n_w_M+1;
         enddo !do nw = 0,n_w-1
c$$$         %%%%%%%%
      enddo !do nM = 0,n_M-1
      if (n_w_M.ne.(n_w*n_M)) then
         write(6,'(A,I0,A,I0)') 'Warning: n_w_M ',n_w_M,' n_w*n_M ', n_w
     $        *n_M
      end if !if (n_w_M.ne.(n_w*n_M)) then
c$$$         %%%%%%%%

c$$$         %%%%%%%%
      if (verbose.gt.1) then
         write(6,'(A)') ' checking memory allocation. '
      end if !if (verbose.gt.1) then
      flag_memory_checkset = .true.
      call cxs_r8(n_w,k_c_0_,'k_c_0_',flag_memory_checkset)
      call cxs_r8(n_w,k_c_1_,'k_c_1_',flag_memory_checkset)
      call cxs_r8(n_w,k_c_2_,'k_c_2_',flag_memory_checkset)
      call cxs_c16(n_w,M_k_p_single_,'M_k_p_single_'
     $     ,flag_memory_checkset)
      call cxs_c16(n_w,CTF_k_p_single_,'CTF_k_p_single_'
     $     ,flag_memory_checkset)
      if (flag_memory_checkset.eqv..true.) then
         if (verbose.gt.1) then 
            write(6,'(A)') '[checkset passed]'
         end if !if (verbose.gt.1) then 
      end if !if (flag_memory_checkset.eqv..true.) then
      if (flag_memory_checkset.eqv..false.) then
         write(6,'(A)') '[checkset failed] <-- WARNING'
         stop !exit program due to error.
      end if !if (flag_memory_checkset.eqv..false.) then
c$$$         %%%%%%%%
      
      if (verbose.gt.0) then 
         write(6,'(A)') '[finished get_lsqdata_2]'
      end if !if (verbose.gt.0) then 
      end !subroutine
