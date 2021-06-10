      subroutine alpha_SM_write_0(n_SM_max,n_SM,alpha_SM_,n_p
     $     ,prefix_string)
      implicit none
      integer verbose
      data verbose / 0 /
c$$$      SM storage
      integer *4 n_SM_max ! total (maximum) number of templates to store per image. ;
      integer *4 n_SM ! the actual number of templates stored for this particular image. ;
      real *8 alpha_SM_(0:0) ! array of size n_alpha*n_SM_max storing the image-parameters for each stored template-image pair for this particular image. ;
      integer *4 n_p
      character(len=n_p) prefix_string
      integer nSM,nSM_sub
      include 'nalpha_define.f'
      character(len=1024) format_string

      write(format_string,'(A,I0,A)') '(A,A,I2,A,' , n_alpha ,
     $     '(F10.5,1X))'
      do nSM=0,n_SM-1
         write(6,format_string) prefix_string , 'nSM ' , nSM , ': ' ,
     $           (alpha_SM_(nSM_sub),nSM_sub=n_alpha*nSM,n_alpha*nSM
     $           +n_alpha-1)
      enddo                     !do nSM=0,n_SM_max-1

      end
