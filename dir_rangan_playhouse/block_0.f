      subroutine block_0(verbose,n_A_sub__in,n_A,n_A_per_,n_A_sum_
     $     ,n_A_sub_use,nA_per_min,nA_per_max)
      implicit none
      integer *4 verbose !verbosity level. ;
      integer *4 n_A_sub__in !requested (maximum) number of blocks to split n_A into. ;
      integer *4 n_A !total number of elements in array. ;
      integer *4 n_A_per_(0:0) !array (of length n_A_sub_use) storing the number of elements per block. ;
      integer *4 n_A_sum_(0:0) !array (of length n_A_sub_use) storing the cumulative sum of n_A_per_. ;
      integer *4 n_A_sub_use !actual number of blocks used when splitting n_A. ;
      integer *4 nA_per_min !minimum number of elements in a block. ;
      integer *4 nA_per_max !maximum number of elements in a block. ;
      integer *4 nA_sub !temporary index. ;
      integer *4 n_A_quo !temporary quotient. ;
      integer *4 n_A_rem !temporary remainder. ;
      integer *4 min_i4_f !function output. ;
      integer *4 max_i4_f !function output. ;

      if (verbose.gt.0) then
         write (6,'(A)') ' [entering block_0] '
      end if ! if (verbose.gt.0) then

      n_A_sub_use = min(n_A_sub__in,n_A)
      if (verbose.gt.0) then
         write (6,'(3(A,I0))') ' n_A_sub__in ' , n_A_sub__in , ' n_A ' ,
     $        n_A , ' n_A_sub_use = ',n_A_sub_use
      end if ! if (verbose.gt.0) then
      n_A_quo = n_A/n_A_sub_use
      n_A_rem = mod(n_A,n_A_sub_use)
      if (verbose.gt.0) then
         write(6,'(2(A,I0))') ' n_A_quo ' , n_A_quo , ' n_A_rem ' ,
     $        n_A_rem
      end if !if (verbose.gt.0) then
      do nA_sub=0,n_A_sub_use-1
         if (nA_sub.lt.n_A_rem) then
            n_A_per_(nA_sub) = n_A_quo + 1
         end if !if (nA_sub.lt.n_A_rem) then
         if (nA_sub.ge.n_A_rem) then
            n_A_per_(nA_sub) = n_A_quo + 0
         end if !if (nA_sub.ge.n_A_rem) then
      enddo ! do nA_sub=0,n_A_sub_use-1
      if (n_A_sub_use.ge.1) then
         n_A_sum_(0) = 0
      end if !if (n_A_sub_use.ge.1) then
      do nA_sub=1,n_A_sub_use-1
         n_A_sum_(nA_sub) = n_A_sum_(nA_sub-1) + n_A_per_(nA_sub-1)
      enddo ! do nA_sub=1,n_A_sub_use-1
      nA_per_min = min_i4_f(n_A_sub_use,n_A_per_)
      nA_per_max = max_i4_f(n_A_sub_use,n_A_per_)
      if (verbose.gt.0) then
         call write_all_i4(n_A_sub_use,n_A_per_,11,' n_A_per_: ')
         call write_all_i4(n_A_sub_use,n_A_sum_,11,' n_A_sum_: ')
         write(6,'(A,I0,A,I0)') ' nA_per_min: ' , nA_per_min ,
     $        ' nA_per_max: ' , nA_per_max
      end if ! if (verbose.gt.1) then

      if (verbose.gt.0) then
         write (6,'(A)') ' [finished block_0] '
      end if ! if (verbose.gt.0) then
      end 
