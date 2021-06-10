!> Doxygen comment: ;\n
!> returns unique indices from witihin integer *4 list L_. ;\n
      subroutine unique_I4(n_L,L_,n_I_out,I_out)
      implicit none
      integer *4 n_L,n_I_out
      integer *4 L_(0:n_L-1)
      integer *4 I_out(0:0)
      integer *4 nL,n_I
      integer *4, allocatable :: J_(:)
      integer *4, allocatable :: K_(:)
      external quicksort_i4
      allocate(J_(0:n_L-1))
      allocate(K_(0:n_L-1))
      call cp1_i4(n_L,L_,J_)
      do nL=0,n_L-1
         K_(nL) = nL
      enddo !do nL=0,n_L-1
      call quicksort_i4(0,n_L-1,J_,1,K_,1,quicksort_i4)
      n_I=0
      do nL=0,n_L-1
         if (nL.eq.0) then
            I_out(n_I) = J_(nL)
            n_I = n_I+1
         end if !if (nL.eq.0) then
         if (nL.gt.0) then
            if (J_(nL).gt.J_(nL-1)) then
               I_out(n_I) = J_(nL)
               n_I = n_I+1
            end if !if (J_(nL).gt.J_(nL-1)) then
         end if !if (nL.gt.0) then
      enddo !do nL=0,n_L-1
      n_I_out = n_I
      deallocate(K_)
      deallocate(J_)
      end
      
      
