!> Doxygen comment: ;\n
!> calculates pearson correlation between two c16 arrays. ;\n
      subroutine pearson_c16(n_z,A_,A_start,A_stride,B_,B_start,B_stride
     $     ,C)
      implicit none
      integer verbose
      data verbose / 0 /
      integer n_z,A_start,A_stride,B_start,B_stride
      complex *16 A_(0:(A_start + A_stride*n_z - 1))
      complex *16 B_(0:(A_start + A_stride*n_z - 1))
      complex *16 C,A_avg,A_std,B_avg,B_std
      real *8 A_std_use,B_std_use
      integer nz,na,nb
      if (verbose.gt.0) then
         write (6,'(A,I0)') ' % [entering pearson_c16] n_z ',n_z
      end if
      if (verbose.gt.0) then
         write (6,'(A)') ' A_: '
         write (6,'(4(2F8.3,1X))') (A_(na),na=A_start,A_start+A_stride
     $        *n_z,A_stride)
      end if
      if (verbose.gt.0) then
         write (6,'(A)') ' B_: '
         write (6,'(4(2F8.3,1X))') (B_(nb),nb=B_start,B_start+B_stride
     $        *n_z,B_stride)
      end if
      A_avg = dcmplx(0.0d0,0.0d0)
      A_std = dcmplx(0.0d0,0.0d0)
      B_avg = dcmplx(0.0d0,0.0d0)
      B_std = dcmplx(0.0d0,0.0d0)
      na = A_start
      nb = B_start
      do nz=0,n_z-1
         A_avg = A_avg + A_(na)
         A_std = A_std + dconjg(A_(na))*A_(na)
         B_avg = B_avg + B_(nb)
         B_std = B_std + dconjg(B_(nb))*B_(nb)
         na = na + A_stride
         nb = nb + B_stride
      enddo
      A_avg = A_avg/max(1,n_z)
      A_std = zsqrt(A_std/max(1,n_z) - dconjg(A_avg)*A_avg)
      if (real(A_std).lt.0.0d0) then
         A_std = dcmplx( 1.0 , 0.0)
      end if
      B_avg = B_avg/max(1,n_z)
      B_std = zsqrt(B_std/max(1,n_z) - dconjg(B_avg)*B_avg)
      if (real(B_std).lt.0.0d0) then
         B_std = dcmplx( 1.0 , 0.0)
      end if
      if (verbose.gt.0) then
         write (6,'(A,2F16.3,2F16.3)') ' % % A_avg A_std',A_avg,A_std
         write (6,'(A,2F16.3,2F16.3)') ' % % B_avg B_std',B_avg,B_std
      end if
      A_std_use = zabs(A_std)
      if (A_std_use.lt.0.000001d0) then
         A_std_use = 1.0d0
      end if
      B_std_use = zabs(B_std)
      if (B_std_use.lt.0.000001d0) then
         B_std_use = 1.0d0
      end if
      if (verbose.gt.0) then
         write (6,'(A,F16.3,1X,F16.3)') ' % % A_std_use B_std_use'
     $        ,A_std_use,B_std_use
      end if
      C = dcmplx(0.0d0,0.0d0)
      na = A_start
      nb = B_start
      do nz=0,n_z-1
         C = C + dconjg(A_(na) - A_avg)/A_std_use * (B_(nb) - B_avg)
     $        /B_std_use
         na = na + A_stride
         nb = nb + B_stride
      enddo
      C = C/max(1,n_z)
      if (verbose.gt.0) then
         write (6,'(A,2F16.3)') ' % % C: ',C
      end if
      if (verbose.gt.0) then
         write (6,'(A)') ' % [finished pearson_c16]'
      end if
      end
