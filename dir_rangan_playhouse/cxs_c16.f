      subroutine cxs_c16(n_x,A_,prefix_string,flag_pass)
      implicit none
      integer n_x
      complex *16 A_(0:0)
      character(len=*) prefix_string
      logical flag_pass
      if (A_(n_x).ne.cmplx(-1.0d0,0.0d0)) then
         write(6,'(A,A,A,I0,A)') ' Warning: cxs_c16 failed for: ' ,
     $        prefix_string , '(' , n_x , ')'
         flag_pass = .false.
      end if !if (A_(n_x).ne.cmplx(-1.0d0,0.0d0)) then
      end
