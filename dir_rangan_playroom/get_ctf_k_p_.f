      subroutine get_ctf_k_p_(ngridr,CTF_k_p_,param_0,param_1,param_2
     $     ,param_3,param_4)
      implicit none
      integer verbose
      data verbose / 0 /
      integer ngridr
      complex *16 CTF_k_p_(0:0)
      real *8 param_0,param_1,param_2,param_3,param_4
c$$$      cryoem parameters
      integer ityper
      real *8 rmax
      real *8, allocatable :: xnodesr(:)
      real *8, allocatable :: wtsr(:)
      integer *4 ntemplatesize,ncur
      integer *4, allocatable :: nlats_(:)
      integer *4, allocatable :: ngridc_(:)
      integer *4, allocatable :: icstart_(:)
c$$$      indices
      integer *4 n_r,nr,n_w_max,n_A,na
      integer *4, allocatable :: n_w_(:)
c$$$      grids 
      integer *4 n_gridpoints,ng
      real *8 pi
      real *8 max_x_c,max_k_c,dtmp
      real *8, allocatable :: grid_x_c_(:)
      real *8, allocatable :: grid_k_c_(:)
      real *8, allocatable :: grid_k_p_(:)
c$$$      Template S
      complex *16, allocatable :: S_x_c_(:)
      complex *16, allocatable :: S_k_c_(:)
      complex *16, allocatable :: S_k_p_(:)
c$$$      innerproducts
      complex *16 C_S_x_c,C_S_k_c,C_S_k_p
      external get_ctf_x_c

      if (verbose.gt.0) then
         write(6,'(A,I0)') '[entering get_ctf_k_p_], ngridr: ',ngridr
      end if

c$$$         We define xnodesr via 'getgridr'
      rmax = 1.0d0*ngridr
      ityper = 0
      allocate(xnodesr(ngridr))
      allocate(wtsr(ngridr))
      call getgridr(rmax,ngridr,ityper,xnodesr,wtsr)

c$$$      Calculating template size using 'get_template_size'
      pi = 4*atan(1.0)
      ncur = ngridr
      allocate(nlats_(0:ncur-1));
      do ng=0,ncur-1
         nlats_(ng) = nint(pi*xnodesr(1+ng))
         if (nlats_(ng).lt.6) nlats_(ng) = 6
      enddo
      allocate(ngridc_(ncur))
      allocate(icstart_(ncur))
      call get_template_size(nlats_,ncur,ntemplatesize,ngridc_,icstart_)
      if (verbose.gt.0) then
         write(6,'(A,I0)') 'ntemplatesize = ',ntemplatesize
      end if
      n_gridpoints = ncur
      if (n_gridpoints.lt.2) then
         write(6,'(A,I0,A)') 'Error n_gridpoints',n_gridpoints,'<2'
      end if

c$$$      indices
      n_r = ncur
      if (n_r.lt.2) then
         write(6,'(A,I0,A)') 'Error n_r',n_r,'<2'
      end if
      allocate(n_w_(0:n_r-1))
      n_A = 0
      do nr=0,n_r-1
         n_w_(nr) = ngridc_(1+nr)
         n_A = n_A + n_w_(nr)
      enddo
      n_w_max = n_w_(nr-1)
      if (verbose.gt.1) then
         write(6,'(A,I0,A,I0)') 'n_w_max ',n_w_max,'; n_A ',n_A
      end if

c$$$      Calculating x-space template on regular cartesian-grid
      max_x_c = 1.0d0
      allocate(grid_x_c_(0:n_gridpoints-1))
      call linspace(0.0d0,max_x_c,n_gridpoints,grid_x_c_)
      allocate(S_x_c_(0:n_gridpoints*n_gridpoints-1))
      call get_ctf_x_c_(n_gridpoints,grid_x_c_,max_x_c,n_gridpoints
     $     ,grid_x_c_,max_x_c,S_x_c_,get_ctf_x_c,param_0,param_1,param_2
     $     ,param_3,param_4)
      call innerproduct_c(n_gridpoints,grid_x_c_,n_gridpoints ,grid_x_c_
     $     ,S_x_c_,S_x_c_,C_S_x_c)
      C_S_x_c = zsqrt(C_S_x_c)
      if (verbose.gt.0) write(6,'(A,2F16.3)') ' C_S_x_c: ',C_S_x_c
      
c$$$      Calculating k-space template on regular cartesian-grid      
      max_k_c = (1.0d0*n_gridpoints)/max_x_c
      allocate(grid_k_c_(0:n_gridpoints-1))
      call linspace(0.0d0,max_k_c,n_gridpoints,grid_k_c_)
      allocate(S_k_c_(0:n_gridpoints*n_gridpoints-1))
      call adi_fft2(-1,n_gridpoints,n_gridpoints,S_x_c_,S_k_c_)
      call innerproduct_c(n_gridpoints,grid_k_c_,n_gridpoints ,grid_k_c_
     $     ,S_k_c_,S_k_c_,C_S_k_c)
      C_S_k_c = zsqrt(C_S_k_c)/(n_gridpoints*n_gridpoints)
      if (verbose.gt.0) write(6,'(A,2F16.3)') ' C_S_k_c: ',C_S_k_c

c$$$      Calculating k-space template on quasi-uniform polar-grid
      allocate(S_k_p_(0:ntemplatesize-1))
      allocate(grid_k_p_(0:n_gridpoints-1))
c$$$      instead of:
c$$$      call linspace(0.0d0,max_k_c/2.0,n_gridpoints,grid_k_p_)
c$$$      we use:
      dtmp = max_k_c/2.0/n_gridpoints
      call linspace(dtmp,max_k_c/2.0+dtmp,n_gridpoints,grid_k_p_)
      call interp_c_to_p(n_gridpoints,max_k_c,n_gridpoints,max_k_c
     $     ,S_k_c_,n_gridpoints,grid_k_p_,ngridc_,ntemplatesize,S_k_p_)
      call innerproduct_p(n_gridpoints,grid_k_p_,ngridc_ ,ntemplatesize
     $     ,S_k_p_,S_k_p_,C_S_k_p)
      C_S_k_p = zsqrt(C_S_k_p)/(n_gridpoints*n_gridpoints)      
      if (verbose.gt.0) write(6,'(A,2F16.3)') ' C_S_k_p: ',C_S_k_p

      if (verbose.gt.0) write(6,'(A)') 'copying'
      call cp1_c16(ntemplatesize,S_k_p_,CTF_k_p_)

      if (verbose.gt.0) write(6,'(A)') 'deallocating'
      deallocate(xnodesr)
      deallocate(wtsr)
      deallocate(nlats_)
      deallocate(ngridc_)
      deallocate(icstart_)
      deallocate(n_w_)
      deallocate(grid_x_c_)
      deallocate(S_x_c_)
      deallocate(grid_k_c_)
      deallocate(S_k_c_)
      deallocate(S_k_p_)
      deallocate(grid_k_p_)

      if (verbose.gt.0) then
         write(6,'(A,I0)') '[finished get_ctf_k_p_], ngridr: ',ngridr
      end if

      end
