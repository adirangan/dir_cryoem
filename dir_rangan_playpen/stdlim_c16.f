      subroutine stdlim_c16(n_A,real_or_imag,S_,avg,std)
      implicit none
      integer n_A
      logical real_or_imag
      complex *16 S_(0:n_A-1)
      real *8 avg,std,var
      integer na
      avg=0
      var=0
      std=0
      do na=0,n_A-1
         if (real_or_imag.eqv..true.) then
            avg = avg + real(S_(na))
            var = var + real(S_(na))**2
         else
            avg = avg + aimag(S_(na))
            var = var + aimag(S_(na))**2
         end if
      enddo
      avg = avg / max(1,n_A)
      var = var / max(1,n_A)
      std = dsqrt(var - avg*avg)
      end
      
