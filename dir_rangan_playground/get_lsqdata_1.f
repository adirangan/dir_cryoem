!> Doxygen comment: ;\n
!> Sets up data on the sphere for least-square solve. ;\n
!> This includes the values M_k_c_ and weights weight_CTF_k_c_. ;\n
      subroutine get_lsqdata_1(n_M,I_M_sample_,ld_M,M_k_p__,n_CTF
     $     ,ld_CTF,CTF_k_p_,k_p,n_w_csum,n_w,alpha2d__
     $     ,M_k_c_ ,polar_a_,azimu_b_,weight_CTF_k_c_,n_w_M)
      implicit none
      integer verbose
      data verbose / 0 /
      integer n_M,I_M_sample_(0:n_M-1),ld_M,n_CTF,ld_CTF,n_w_csum,n_w
     $     ,n_w_M
      include 'excerpt_define_nalpha.f'
      complex *16 M_k_p__(0:0)
      complex *16 CTF_k_p_(0:ld_CTF*n_CTF-1)
      real *8 k_p
      real *8 alpha2d__(0:n_alpha-1,0:n_M-1)
      complex *16 M_k_c_(0:n_w*n_M-1)
      real *8 polar_a_(0:n_w*n_M-1)
      real *8 azimu_b_(0:n_w*n_M-1)
      complex *16 weight_CTF_k_c_(0:n_w*n_M-1)
      real *8 pi,cos_polar_a,azimu_b
      real *8 k1_c,k2_c,k3_c,k_p_1
      real *8 delta_x,delta_y,gamma_z,l2_norm
      real *8, allocatable :: k1_c_(:)
      real *8, allocatable :: k2_c_(:)
      real *8, allocatable :: k3_c_(:)
      complex *16, allocatable ::  M_k_p_1_(:)
      complex *16, allocatable ::  CTF_k_p_1_(:)
      integer nM,nw,nctf
      if (verbose.gt.0) write(6,*)
     $     '[entering get_lsqdata_0]'
c
      allocate(k1_c_(0:n_w-1))
      allocate(k2_c_(0:n_w-1))
      allocate(k3_c_(0:n_w-1))
      allocate(M_k_p_1_(0:n_w-1))
      allocate(CTF_k_p_1_(0:n_w-1))
c
      pi=4.0d0*datan(1.0d0)
      n_w_M = 0
      do nM = 0,n_M-1
         nctf = nint(alpha2d__(nalpha_ctf_ind,nM))
         nctf = max(0,min(n_CTF-1,nctf))
         cos_polar_a = dcos(alpha2d__(nalpha_polar_a,nM))
         azimu_b = alpha2d__(nalpha_azimu_b,nM)
         gamma_z = 0.0d0
         k_p_1 = 1.0d0
         call mkonegreatcircle(k_p_1,cos_polar_a,azimu_b,gamma_z,
     1        n_w,k1_c_,k2_c_,k3_c_)
c
         if (n_alpha.ge.1+nalpha_l2_norm) then
            l2_norm = alpha2d__(nalpha_l2_norm,nM)
         else
            l2_norm = 1.0d0
         end if
            
         if (dabs(l2_norm).lt.1.0d-15) then
            l2_norm=1.0d0
         end if
         do nw = 0,n_w-1
            M_k_p_1_(nw) = M_k_p__(n_w_csum-1+nw + I_M_sample_(nM)*ld_M)
     $           /l2_norm
            CTF_k_p_1_(nw) = CTF_k_p_(n_w_csum-1+nw + nctf*ld_CTF)
         enddo !do nw = 0,n_w-1
c
         gamma_z = -alpha2d__(nalpha_gamma_z,nM)
         call rotate_p_to_p_1(n_w,M_k_p_1_,+gamma_z,M_k_p_1_)
         call rotate_p_to_p_1(n_w,CTF_k_p_1_,+gamma_z,CTF_k_p_1_)
         delta_x = -alpha2d__(nalpha_delta_x,nM)
         delta_y = -alpha2d__(nalpha_delta_y,nM)
         call transf_p_to_p_1(k_p,n_w,M_k_p_1_,+delta_x,+delta_y
     $        ,M_k_p_1_)
c
         do nw = 0,n_w-1
            k1_c = k1_c_(nw)
            k2_c = k2_c_(nw)
            k3_c = k3_c_(nw)
            cos_polar_a = k3_c
            azimu_b = datan2(k2_c,k1_c)
            if (azimu_b.lt.0.0d0) azimu_b = azimu_b + 2.0d0*pi
            azimu_b_(n_w_M) = azimu_b
            polar_a_(n_w_M) = datan2(dsqrt(k1_c**2+ k2_c**2)
     $           ,cos_polar_a)
            M_k_c_(n_w_M) = M_k_p_1_(nw)
            weight_CTF_k_c_(n_w_M) = CTF_k_p_1_(nw)
            n_w_M = n_w_M+1;
         enddo !do nw = 0,n_w-1
      enddo !do nM = 0,n_M-1
      
      if (n_w_M.ne.(n_w*n_M)) then
         write(6,'(A,I0,A,I0)') 'Warning: n_w_M ',n_w_M,' n_w*n_M ', n_w
     $        *n_M
      end if !if (n_w_M.ne.(n_w*n_M)) then
      
      if (verbose.gt.0) write(6,*)
     $     '[finished get_lsqdata_0]'
      return
      end
