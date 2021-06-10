c     gfortran -o test_Jtaylor_dr.out test_Jtaylor_dr.f ; ./test_Jtaylor_dr.out ;
      program test_Jtaylor_dr
      implicit none
c$$$      arrays hold Jtaylor expansion
      integer l_max,n_J
      complex *16, allocatable :: f_(:)
      integer, allocatable :: b_(:)
      integer, allocatable :: p_(:)
      integer l,j,jc
      real *8 pi
      integer n_z,n_phi,nz,nphi
      real *8, allocatable :: z_(:)
      real *8, allocatable :: phi_(:)
      real *8 z,zl,phi,theta,omega
      complex *16, allocatable :: J_(:)
      complex *16, allocatable :: T_(:)
      real *8, allocatable :: E_(:)
      pi = 4.0*atan(1.0)

c$$$      Calculating Jtaylor expansion 
      write(6,'(A)') 'Enter expansion level: '
      read(5,*) l_max
      n_J = ((l_max+1) * (l_max+2)) / 2
      allocate(f_(0:n_J-1))
      allocate(b_(0:n_J-1))
      allocate(p_(0:n_J-1))
      call get_Jtaylor(l_max,n_J,f_,b_,p_)
      write(6,'(I0,A,I0,A)') l_max,'-level expansion: ',n_J
     $     ,'-terms total'
      write(6,'(A)') '--------------------------------'

      n_z = 16
      allocate(z_(0:n_z-1))
      allocate(E_(0:n_z-1))
      n_phi = 1024
      allocate(phi_(0:n_phi-1))
      allocate(J_(0:n_phi-1))
      allocate(T_(0:n_phi-1))

      do nz=0,n_z-1
         z_(nz) = 0.0 + nz*8.0/n_z
      enddo

      do nphi=0,n_phi-1
         phi_(nphi) = 0.0 + nphi*2*pi/n_phi
      enddo

c$$$      Testing against exp(-i*z*cos(phi))
      do nz=0,n_z-1
         E_(nz) = 0.0
         z = z_(nz)
         do nphi=0,n_phi-1
            phi = phi_(nphi)
            theta = -z*cos(phi)
            J_(nphi) = cmplx(cos(theta),sin(theta))
            T_(nphi) = cmplx( 0.0 , 0.0 )
            zl = 1.0
            jc=0
            do l=0,l_max
               if (l.gt.0) then
                  zl = zl*z
               end if
               do j=0,l
                  omega = p_(jc)*phi
                  T_(nphi) = T_(nphi) + f_(jc)*b_(jc)*cmplx(cos(omega)
     $                 ,sin(omega))*zl
                  jc = jc + 1
               enddo
            enddo
            E_(nz) = E_(nz) + abs(J_(nphi)-T_(nphi))**2
         enddo
         E_(nz) = E_(nz) / n_phi
         write(6,'(A,F3.1,A,F8.3)') 'z=',z,'; log10(error) = '
     $        ,log10(E_(nz))
      enddo

      stop
      end

      include 'get_Jtaylor.f'
      include 'linspace.f'

