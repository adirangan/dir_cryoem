      subroutine get_f_x_c_from_file_0(n_x_c_max,str_dim,str_density
     $     ,f_x_c_)
      implicit none
      integer n_x_c_max !integer: maximum number of x-values (in x-space cartesian coordinates) used to generate initial molecule. sometimes named ngrid. ;
      character(len=*) str_dim !string array: name of dimension file. ;
      character(len=*) str_density !string array: name of density file. ;
      complex *16 f_x_c_(0:0) !complex *16 array (size at least n_x_c_max**3). ;
      integer n_x_c_from_file
      integer n_M_from_file
      integer nx,ny,nz,nt
      integer nx_from_file,ny_from_file,nz_from_file,nt_from_file
      real *8 value_from_file
      real *8 al2_c16_f

      call cl1_c16(n_x_c_max**3,f_x_c_)

      open(20,FILE=str_dim);
      read(20,*) n_x_c_from_file,n_M_from_file
      close(20)
      
      open(12,FILE=str_density)
      nt_from_file=0
      nz=0
      do nz_from_file=0,n_x_c_from_file-1
      ny=0
      do ny_from_file=0,n_x_c_from_file-1
      nx=0
      do nx_from_file=0,n_x_c_from_file-1
c$$$         nt_from_file = nx_from_file+ny_from_file*n_x_c_from_file+nz_from_file*n_x_c_from_file**2
         nt = nx+ny*n_x_c_max+nz*n_x_c_max**2
         read(12,*) value_from_file
         f_x_c_(nt) = f_x_c_(nt) + value_from_file
         nt_from_file=nt_from_file+1
      if ((nx_from_file*n_x_c_max)/n_x_c_from_file.gt.nx) nx=nx+1
      enddo !do nx_from_file=0,n_x_c_from_file-1
      if ((ny_from_file*n_x_c_max)/n_x_c_from_file.gt.ny) ny=ny+1
      enddo !do ny_from_file=0,n_x_c_from_file-1
      if ((nz_from_file*n_x_c_max)/n_x_c_from_file.gt.nz) nz=nz+1
      enddo !do nz_from_file=0,n_x_c_from_file-1
      close(12)

      end !subroutine
