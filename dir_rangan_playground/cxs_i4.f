!> Doxygen comment: ;\n
!> checks final+1 entry in integer *4 array ;\n
!> assumes array is set with cs1_i4 ;\n
      subroutine cxs_i4(n_x,A_,prefix_string,flag_pass)
      implicit none
      integer n_x
      integer *4 A_(0:0)
      character(len=*) prefix_string
      logical flag_pass
      if (A_(n_x).ne.-1) then         
         write(6,'(A,A,A,I0,A)') ' Warning: cxs_i4 failed for: ' ,
     $        prefix_string , '(' , n_x , ')'
         flag_pass = .false.
      end if !if (A_(n_x).ne.-1) then
      end
