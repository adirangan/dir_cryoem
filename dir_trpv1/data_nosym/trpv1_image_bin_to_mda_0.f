!> Doxygen comment: ;\n
!> converts n_image_sub images from image_bin into mda format. ;\n
c$$$      gfortran -o image_bin_to_mda_0.out image_bin_to_mda_0.f ; ./image_bin_to_mda_0.out ;
      implicit none
      integer n_image,n_image_sub
      integer n_x_c
      complex *16, allocatable :: M_x_c___(:,:,:)
      real *8, allocatable :: N_x_c___(:,:,:)
      character(len=1024) :: str_dir
      character(len=1024) :: str_fname
      integer nimage,nx_c,ny_c
      integer n_d
      integer, allocatable :: d_(:)
      
      str_dir = '/data/rangan/dir_cryoem/dir_trpv1/data_nosym'

      write(str_fname,'(A,A)') trim(str_dir) , '/dims'
      open(20,FILE=str_fname);
      read(20,*) n_x_c,n_image
      close(20)
      n_image_sub = min(1024,n_image)

      allocate(M_x_c___(0:n_x_c-1,0:n_x_c-1,0:n_image_sub-1))
      allocate(N_x_c___(0:n_x_c-1,0:n_x_c-1,0:n_image_sub-1))
      
      write(str_fname,'(A,A)') trim(str_dir) , '/images_bin'
      open(20,FILE=str_fname,FORM="unformatted");
      do nimage = 0,n_image_sub-1
      do nx_c = 0,n_x_c-1
      read(20) (M_x_c___(ny_c,nx_c,nimage),ny_c=0,n_x_c-1)
      enddo !do nx_c = 0,n_x_c-1
      enddo !do nimage = 0,n_image_sub-1
      close(20)

      do nimage = 0,n_image_sub-1
      do nx_c = 0,n_x_c-1
      do ny_c = 0,n_x_c-1
         N_x_c___(ny_c,nx_c,nimage) = real(M_x_c___(ny_c,nx_c,nimage))
      enddo !do ny_c = 0,n_x_c-1
      enddo !do nx_c = 0,n_x_c-1
      enddo !do nimage = 0,n_image_sub-1

      n_d = 3
      allocate(d_(0:n_d-1));
      
      write(str_fname,'(A,A)') trim(str_dir) , '/images_mda'
      d_(0) = n_x_c
      d_(1) = n_x_c
      d_(2) = n_image_sub
      call MDA_write_r8(n_d,d_,N_x_c___,str_fname)

      stop
      end

      include 'MDA_write_r8.f'
