!> Doxygen comment: ;\n
!> get an equispaced grid of gamma-values from [0,2.0d0*pi). ;\n
      subroutine get_gamma_0(n_gamma_z,gamma_z_)
      implicit none
      integer verbose
      data verbose / 0 /
      integer *4 n_gamma_z
      real *8 gamma_z_(0:n_gamma_z-1)
      integer ngz
      real *8 gamma_z
      real *8 pi
      if (verbose.gt.0) then
         write(6,'(A)') '[entering get_gamma_0]'
      end if !if (verbose.gt.0) then
      pi = 4.0d0*datan(1.0d0)
      do ngz=0,n_gamma_z-1
         if (n_gamma_z.gt.1) then
            gamma_z = (2.0d0*pi*ngz)/n_gamma_z
         else
            gamma_z = 0.0d0
         end if
         gamma_z_(ngz) = gamma_z
      enddo
      if (verbose.gt.0) then
         write(6,'(A)') '[finished get_gamma_0]'
      end if !if (verbose.gt.0) then
      end
      
