      subroutine tmp_lsq_read(multaha)
      implicit none
      integer *4 nk
      complex *16, allocatable :: M_k_c_(:)
      real *8, allocatable :: polar_a_(:)
      real *8, allocatable :: azimu_b_(:)
      integer *4 n_w_M
      integer *4 n_Y_l,n_Y_lm
      integer *4 lsq_n_quad
      complex *16, allocatable :: weight_CTF_k_c_(:)
      external multaha
      integer *4 lsq_interpolation_order
      integer *4 lsq_n_iteration
      integer *4 lsq_niteration
      real *8 lsq_eps
      complex *16, allocatable :: Y_cur_(:)
      character(len=1024) tmp_dir
      character(len=1024) tmp_str
      integer *4 tmp_tab
      real *8, allocatable :: xnodesth(:)
      real *8, allocatable :: wtsth(:)
      tmp_dir = '/data/rangan/dir_cryoem/dir_rangan_playground/dir_tmp/'
c$$$  %%%%%%%%
      tmp_str = trim(adjustl(tmp_dir)) // 'nk'
      open(20,FILE=tmp_str)
      read(20,*) nk
      close(20)
c$$$  %%%%%%%%
      tmp_str = trim(adjustl(tmp_dir)) // 'n_w_M'
      open(20,FILE=tmp_str)
      read(20,*) n_w_M
      close(20)
c$$$  %%%%%%%%
      allocate(M_k_c_(0:1+n_w_M-1))
      call cs1_c16(n_w_M,M_k_c_)
      tmp_str = trim(adjustl(tmp_dir)) // 'M_k_c_'
      open(20,FILE=tmp_str)
      do tmp_tab=0,n_w_M-1
         read(20,*) M_k_c_(tmp_tab)
      enddo                     !do tmp_tab=0,n_w_M-1
      close(20)
c$$$  %%%%%%%%
      allocate(polar_a_(0:1+n_w_M-1))
      call cs1_r8(n_w_M,polar_a_)
      tmp_str = trim(adjustl(tmp_dir)) // 'polar_a_'
      open(20,FILE=tmp_str)
      do tmp_tab=0,n_w_M-1
         read(20,*) polar_a_(tmp_tab)
      enddo                     !do tmp_tab=0,n_w_M-1
      close(20)
c$$$  %%%%%%%%
      allocate(azimu_b_(0:1+n_w_M-1))
      call cs1_r8(n_w_M,azimu_b_)
      tmp_str = trim(adjustl(tmp_dir)) // 'azimu_b_'
      open(20,FILE=tmp_str)
      do tmp_tab=0,n_w_M-1
         read(20,*) azimu_b_(tmp_tab)
      enddo                     !do tmp_tab=0,n_w_M-1
      close(20)
c$$$  %%%%%%%%
      tmp_str = trim(adjustl(tmp_dir)) // 'n_Y_l'
      open(20,FILE=tmp_str)
      read(20,*) n_Y_l
      close(20)
c$$$  %%%%%%%%
      tmp_str = trim(adjustl(tmp_dir)) // 'lsq_n_quad'
      open(20,FILE=tmp_str)
      read(20,*) lsq_n_quad
      close(20)
c$$$  %%%%%%%%
      tmp_str = trim(adjustl(tmp_dir)) //
     $     'lsq_interpolation_order'
      open(20,FILE=tmp_str)
      read(20,*) lsq_interpolation_order
      close(20)
c$$$  %%%%%%%%
      allocate(weight_CTF_k_c_(0:1+n_w_M-1))
      call cs1_c16(n_w_M,weight_CTF_k_c_)
      tmp_str = trim(adjustl(tmp_dir)) // 'weight_CTF_k_c_'
      open(20,FILE=tmp_str)
      do tmp_tab=0,n_w_M-1
         read(20,*) weight_CTF_k_c_(tmp_tab)
      enddo                     !do tmp_tab=0,n_w_M-1
      close(20)
c$$$  %%%%%%%%
      tmp_str = trim(adjustl(tmp_dir)) // 'lsq_eps'
      open(20,FILE=tmp_str)
      read(20,*) lsq_eps
      close(20)
c$$$  %%%%%%%%
      tmp_str = trim(adjustl(tmp_dir)) // 'lsq_n_iteration'
      open(20,FILE=tmp_str)
      read(20,*) lsq_n_iteration
      close(20)
c$$$  %%%%%%%%
      n_Y_lm = (n_Y_l+1)**2
      allocate(Y_cur_(0:1+n_Y_lm-1))
      call cs1_c16(n_Y_lm,Y_cur_)
      write(6,'(A)') '[tmp_lsq]'
      write(6,'(7(A,I0))')
     $     ' nk: ' , nk
     $     , ' n_w_M: ' , n_w_M
     $     , ' n_Y_l: ' , n_Y_l
     $     , ' n_Y_lm: ' , n_Y_lm
     $     , ' lsq_n_quad: ' , lsq_n_quad
     $     , ' lsq_interpolation_order: ' 
     $     , lsq_interpolation_order
     $     , ' lsq_n_iteration: ' , lsq_n_iteration
      call print_sub_c16(n_w_M,M_k_c_,' M_k_c_: ')
      call print_sub_c16(n_w_M,weight_CTF_k_c_
     $     ,' C_k_c_: ')
      call print_sub_r8(n_w_M,polar_a_,' polar_a_: ')
      call print_sub_r8(n_w_M,azimu_b_,' azimu_b_: ')
      call lsqctfsolve(
     $     M_k_c_
     $     ,polar_a_
     $     ,azimu_b_
     $     ,n_w_M
     $     ,Y_cur_
     $     ,n_Y_l
     $     ,lsq_n_quad
     $     ,lsq_interpolation_order
     $     ,weight_CTF_k_c_
     $     ,multaha
     $     ,lsq_eps
     $     ,lsq_n_iteration
     $     ,lsq_niteration
     $     )
      call print_sub_c16(n_Y_lm,Y_cur_,' Y_cur_: ')
      allocate(xnodesth(0:1+lsq_n_quad-1))
      call cs1_r8(lsq_n_quad,xnodesth)
      allocate(wtsth(0:1+lsq_n_quad-1))
      call cs1_r8(lsq_n_quad,wtsth)
      call chebychev_quad(lsq_n_quad,xnodesth,wtsth)
      call lsqctfsolve_marina(
     $     M_k_c_
     $     ,polar_a_
     $     ,azimu_b_
     $     ,n_w_M
     $     ,Y_cur_
     $     ,n_Y_l
     $     ,lsq_n_quad
     $     ,2*lsq_n_quad
     $     ,lsq_interpolation_order
     $     ,xnodesth
     $     ,wtsth
     $     ,weight_CTF_k_c_
     $     ,multaha
     $     ,lsq_eps
     $     ,lsq_niteration
     $     )
      call print_sub_c16(n_Y_lm,Y_cur_,' Y_cur_: ')
      end !subroutine

c$$$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      

      subroutine chebychev_quad(ngrid,xnodes,wts)
      implicit none
      integer ngrid
      real *8 xnodes(ngrid),wts(ngrid)
      integer i
      integer itype
      real *8 t(2000),w(2000),b(2000)
      real *8 alpha, beta
      itype = 1
      alpha = 0.0d0
      beta = 0.0d0
      call chebexps(itype,ngrid,t,alpha,beta,w)
      do i = 1,ngrid
         xnodes(i) = t(ngrid+1-i)
         wts(i) = w(i)
      enddo
      return
      end
