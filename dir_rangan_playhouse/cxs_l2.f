      subroutine cxs_l2(n_x,A_,prefix_string,flag_pass)
      implicit none
      integer n_x
      logical A_(0:0)
      character(len=*) prefix_string
      logical flag_pass
      if (A_(n_x).neqv..true.) then
         write(6,'(A,A,A,I0,A)') ' Warning: cxs_l2 failed for: ' ,
     $        prefix_string , '(' , n_x , ')'
         flag_pass = .false.
      end if !if (A_(n_x).neqv..true.) then
      end
