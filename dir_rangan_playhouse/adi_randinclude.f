      subroutine adi_randinclude(rseed,n_LT,LT_,n_S,S_used_,S_used_total
     $     ,n_add)
      implicit none
      integer *4 rseed
      integer *4 n_LT,nLT,n_S,nS,S_used_total,n_add
      integer *4 LT_(0:0)
      logical S_used_(0:0)
      integer *4, allocatable :: I_(:)
      integer *4 n_I,nI
      integer *4, allocatable :: J_(:)
      integer *4 n_J,nJ
      allocate(I_(0:n_S-1))
      allocate(J_(0:n_S-1))
      n_I = n_S - S_used_total
      n_J = min(n_add,n_I)
      do nJ=0,n_J-1
         J_(nJ) = nJ
      enddo
      call adi_randperm(rseed,n_I,I_)
      call quicksort_i4(0,n_J-1,I_,1,J_,1,quicksort_i4)
      nJ=0
      nI=0
      nS=0
      do while ((nS.lt.n_S) .and. (nI.lt.n_I) .and. (nJ.lt.n_J))
         if (S_used_(ns).eqv..true.) then
            ns = ns+1
         else if (S_used_(ns).eqv..false.) then
            if (nI.eq.I_(nJ)) then
               LT_(n_LT) = ns
               n_LT = n_LT+1
               S_used_(ns) = .true.
               S_used_total = S_used_total+1
               nJ = nJ+1
            end if !if not used but in I_
            nI = nI+1
            nS = nS+1
         end if !if
      enddo !do while
      deallocate(J_)
      deallocate(I_)
      end
