      subroutine convert_alpha2d_adi_to_marina_0(
     $     n_M
     $     ,alpha2d_adi__
     $     ,alpha2d_marina__
     $     )
      implicit none
      integer *4 n_M !integer: number of images. ;
      include 'excerpt_define_nalpha.f'
      real *8 alpha2d_adi__(0:n_alpha-1,0:n_M-1) !real *8 array (2d) ;
      real *8 alpha2d_marina__(0:n_alpha-1,0:n_M-1) !real *8 array (2d) ;
      real *8 delta_(0:1) !temporary: displacement. ;
      real *8 gamma_z !temporary: in-plane rotation. ;
      real *8 pi
      integer *4 nM
      pi = 4.0d0*datan(1.0d0)
      do nM=0,n_M-1
         call cp1_r8(n_alpha,alpha2d_adi__(0,nM),alpha2d_marina__(0,nM))
c$$$         delta_(0) = alpha2d_adi__(nalpha_delta_x,nM)
c$$$         delta_(1) = alpha2d_adi__(nalpha_delta_y,nM)
c$$$         gamma_z = alpha2d_adi__(nalpha_gamma_z,nM)
c$$$         call rotate_delta(+gamma_z,delta_,delta_)
c$$$         alpha2d_marina__(nalpha_delta_x,nM) = -2.0d0*pi*delta_(0)
c$$$         alpha2d_marina__(nalpha_delta_y,nM) = -2.0d0*pi*delta_(1)
c$$$         alpha2d_marina__(nalpha_gamma_z,nM) = -gamma_z
         call convert_gamma_delta_adi_to_marina_0(
     $        alpha2d_adi__(nalpha_gamma_z,nM)
     $        ,alpha2d_adi__(nalpha_delta_x,nM)
     $        ,alpha2d_adi__(nalpha_delta_y,nM)
     $        ,alpha2d_marina__(nalpha_gamma_z,nM)
     $        ,alpha2d_marina__(nalpha_delta_x,nM)
     $        ,alpha2d_marina__(nalpha_delta_y,nM)
     $        )
      enddo !do nM=0,n_M-1
      end !subroutine

      subroutine convert_alpha2d_marina_to_adi_0(
     $     n_M
     $     ,alpha2d_marina__
     $     ,alpha2d_adi__
     $     )
      implicit none
      integer *4 n_M !integer: number of images. ;
      include 'excerpt_define_nalpha.f'
      real *8 alpha2d_marina__(0:n_alpha-1,0:n_M-1) !real *8 array (2d) ;
      real *8 alpha2d_adi__(0:n_alpha-1,0:n_M-1) !real *8 array (2d) ;
      real *8 delta_(0:1) !temporary: displacement. ;
      real *8 gamma_z !temporary: in-plane rotation. ;
      real *8 pi
      integer *4 nM
      pi = 4.0d0*datan(1.0d0)
      do nM=0,n_M-1
         call cp1_r8(n_alpha,alpha2d_marina__(0,nM),alpha2d_adi__(0,nM))
c$$$         delta_(0) = -alpha2d_marina__(nalpha_delta_x,nM)/(2.0d0*pi)
c$$$         delta_(1) = -alpha2d_marina__(nalpha_delta_y,nM)/(2.0d0*pi)
c$$$         gamma_z = -alpha2d_marina__(nalpha_gamma_z,nM)
c$$$         call rotate_delta(-gamma_z,delta_,delta_)
c$$$         alpha2d_adi__(nalpha_delta_x,nM) = delta_(0)
c$$$         alpha2d_adi__(nalpha_delta_y,nM) = delta_(1)
c$$$         alpha2d_adi__(nalpha_gamma_z,nM) = gamma_z
         call convert_gamma_delta_marina_to_adi_0(
     $        alpha2d_marina__(nalpha_gamma_z,nM)
     $        ,alpha2d_marina__(nalpha_delta_x,nM)
     $        ,alpha2d_marina__(nalpha_delta_y,nM)
     $        ,alpha2d_adi__(nalpha_gamma_z,nM)
     $        ,alpha2d_adi__(nalpha_delta_x,nM)
     $        ,alpha2d_adi__(nalpha_delta_y,nM)
     $        )
      enddo !do nM=0,n_M-1
      end !subroutine

      subroutine convert_gamma_delta_adi_to_marina_0(
     $     gamma_z_adi
     $     ,delta_x_adi
     $     ,delta_y_adi
     $     ,gamma_z_marina
     $     ,delta_x_marina
     $     ,delta_y_marina
     $     )
      real *8 gamma_z_adi
      real *8 delta_x_adi
      real *8 delta_y_adi
      real *8 gamma_z_marina
      real *8 delta_x_marina
      real *8 delta_y_marina
      real *8 delta_adi_(0:1)
      real *8 delta_marina_(0:1)
      real *8 pi
      pi = 4.0d0*datan(1.0d0)
      delta_adi_(0) = delta_x_adi
      delta_adi_(1) = delta_y_adi
      call rotate_delta(+gamma_z_adi,delta_adi_,delta_marina_)
      delta_x_marina = -2.0d0*pi*delta_marina_(0)
      delta_y_marina = -2.0d0*pi*delta_marina_(1)
      gamma_z_marina = -gamma_z_adi
      call periodize_r8(gamma_z_marina,0.0d0,2.0d0*pi,gamma_z_marina)
      end !subroutine

      subroutine convert_gamma_delta_marina_to_adi_0(
     $     gamma_z_marina
     $     ,delta_x_marina
     $     ,delta_y_marina
     $     ,gamma_z_adi
     $     ,delta_x_adi
     $     ,delta_y_adi
     $     )
      real *8 gamma_z_marina
      real *8 delta_x_marina
      real *8 delta_y_marina
      real *8 gamma_z_adi
      real *8 delta_x_adi
      real *8 delta_y_adi
      real *8 delta_marina_(0:1)
      real *8 delta_adi_(0:1)
      real *8 pi
      pi = 4.0d0*datan(1.0d0)
      delta_marina_(0) = -delta_x_marina/(2.0d0*pi)
      delta_marina_(1) = -delta_y_marina/(2.0d0*pi)
      call rotate_delta(+gamma_z_marina,delta_marina_,delta_adi_)
      delta_x_adi = delta_adi_(0)
      delta_y_adi = delta_adi_(1)
      gamma_z_adi = -gamma_z_marina
      call periodize_r8(gamma_z_adi,0.0d0,2.0d0*pi,gamma_z_adi)
      end !subroutine
      
      
