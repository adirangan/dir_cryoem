c     gfortran -o quicksort_c16_dr.out quicksort_c16_dr.f ; ./quicksort_c16_dr.out ;
      program quicksort_c16_dr
      implicit none
      integer n_Z,nz
      parameter (n_Z=13)
      complex *16 Z_(0:n_Z-1)
      integer *4 I_(0:n_Z-1),J_(0:n_Z-1),K_(0:n_Z-1)
      real *8 D_(0:n_Z-1)
      external quicksort_c16
      external quicksort_i4
      do nz=0,n_Z-1
         Z_(nz) = dcmplx( 0.5d0*mod(1+7*nz,n_Z) , 0.0d0 )
         J_(nz) = nz
         D_(nz) = real(Z_(nz))
      enddo
      write(6,'(13I4)') (J_(nz),nz=0,n_Z-1)
      write(6,'(13F4.1)') (D_(nz),nz=0,n_Z-1)
      call quicksort_c16(0,n_Z-1,Z_,1,J_,1,quicksort_c16)
      do nz=0,n_Z-1
         I_(nz) = J_(nz)         
         K_(nz) = nz
      enddo
      call quicksort_i4(0,n_Z-1,I_,1,K_,1,quicksort_i4)
      do nz=0,n_Z-1
         D_(nz) = real(Z_(nz))
      enddo
      write(6,'(13I4)') (J_(nz),nz=0,n_Z-1)
      write(6,'(13I4)') (K_(nz),nz=0,n_Z-1)
      write(6,'(13F4.1)') (D_(nz),nz=0,n_Z-1)
      stop
      end

      include 'quicksort_c16.f'
      include 'quicksort_i4.f'
