!> Doxygen comment: ;\n
!> sum of logical array. ;\n
!> .true. counts as 1, and .false. counts as 0. ;\n
      integer *4 function sum_l2_f(n_x,x_)
      integer *4 n_x,nx
      logical x_(0:n_x-1)
      integer *4 s
      if (n_x.le.0) then
         s=0
      else
         s=0
         do nx=0,n_x-1
            if (x_(nx).eqv..true.) then
               s = s+1
            end if !if (x_(nx).eqv..true.) then
         enddo !do nx=0,n_x-1
      end if
      sum_l2_f = s
      return
      end
