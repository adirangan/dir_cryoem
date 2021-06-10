      subroutine quicksort_i4(nl,nr,I_,I_stride,J_,J_stride
     $     ,p_quicksort_i4)
      implicit none
      integer verbose
      data verbose / 0 /
      integer nl,nr,I_stride,J_stride
      integer *4 I_(0:0)
      integer *4 J_(0:0)
      external p_quicksort_i4
      integer nj,nd
c$$$      character(len=20) format_string
      if (verbose.gt.0) then
         write(6,'(A,I0,1X,I0)') '[entering quicksort_i4] nl,nr:',nl,nr
c$$$         if (nl.le.nr) then
c$$$            write(format_string,'(A,I0,A)') '(',(nr-nl+1),'I0)'
c$$$            write(6,format_string) (I_(nd),nd=nl,nr)
c$$$         end if
      end if
      nj=0
      if (nl.lt.nr) then
         call quicksort_i4_excerpt(nl,nr,I_,I_stride,J_,J_stride,nj)
         call p_quicksort_i4(nl,nj-1,I_,I_stride,J_,J_stride
     $        ,p_quicksort_i4)
         call p_quicksort_i4(nj+1,nr,I_,I_stride,J_,J_stride
     $        ,p_quicksort_i4)
      end if
c$$$      if (verbose.gt.0) then
c$$$         if (nl.le.nr) then
c$$$            write(format_string,'(A,I0,A)') '(',(nr-nl+1),'I0)'
c$$$            write(6,format_string) (I_(nd),nd=nl,nr)
c$$$         end if
c$$$      end if
      if (verbose.gt.0) then
         write(6,'(A,I0,1X,I0)') '[finished quicksort_i4] nl,nr:',nl,nr
      end if
      end

      subroutine quicksort_i4_excerpt(nl,nr,I_,I_stride,J_,J_stride
     $     ,nj_out)
      implicit none
      integer verbose
      data verbose / 0 /
      integer nl,nr,I_stride,J_stride,nj_out
      integer *4 I_(0:0)
      integer *4 J_(0:0)
      integer *4 I_pivot,I_C
      integer ni,nj,nd,J_C
c$$$      character(len=20) format_string
      if (verbose.gt.0) then
         write(6,'(A,I0,1X,I0)')
     $        '[entering quicksort_i4_excerpt] nl,nr: ',nl,nr
c$$$         if (nl.le.nr) then
c$$$            write(format_string,'(A,I0,A)') '(',(nr-nl+1),'I0)'
c$$$            write(6,format_string) (I_(nd),nd=nl,nr)
c$$$         end if
      end if
      I_pivot = I_(nl*I_stride)
      ni = nl
      nj = nr+1
      if (verbose.gt.0) then
         write(6,'(A,I0,1X,I0,1X,I0)') 'ni,nj,I_pivot: '
     $        ,ni,nj,I_pivot
      end if
      do while (.true.)
         do while (.true.)
            ni = ni+1
            if ((I_(ni*I_stride).gt.I_pivot) .or.
     $           (ni.gt.nr)) then
               goto 10
            end if
         enddo
 10      continue
         if (verbose.gt.0) then
            write(6,'(A,I0,1X,I0)') 'stopping: ni,I_(ni): ',ni
     $           ,I_(ni*I_stride)
         end if
         do while (.true.)
            nj = nj-1
            if (I_(nj*I_stride).le.I_pivot) then
               goto 20
            end if
         enddo
 20      continue
         if (verbose.gt.0) then
            write(6,'(A,I0,1X,I0)') 'stopping: nj,I_(nj): ',nj
     $           ,I_(nj*I_stride)
         end if
         if (ni.ge.nj) then
            goto 30
         end if
         if (verbose.gt.0) then
            write(6,'(A,I0,1X,I0)') 'ni<nj: swapping ',I_(ni*I_stride)
     $           ,I_(nj*I_stride)
         end if
         I_C = I_(ni*I_stride)
         I_(ni*I_stride) = I_(nj*I_stride)
         I_(nj*I_stride) = I_C
         J_C = J_(ni*J_stride)
         J_(ni*J_stride) = J_(nj*J_stride)
         J_(nj*J_stride) = J_C
      enddo
 30   continue
      if (verbose.gt.0) then
         write(6,'(A,I0,1X,I0)') 'final: swapping ',I_(nl*I_stride)
     $        ,I_(nj*I_stride)
      end if
      I_C = I_(nl*I_stride)
      I_(nl*I_stride) = I_(nj*I_stride)
      I_(nj*I_stride) = I_C
      J_C = J_(nl*J_stride)
      J_(nl*J_stride) = J_(nj*J_stride)
      J_(nj*J_stride) = J_C
c$$$      if (verbose.gt.0) then
c$$$         if (nl.le.nr) then
c$$$            write(format_string,'(A,I0,A)') '(',(nr-nl+1),'I0)'
c$$$            write(6,format_string) (I_(nd),nd=nl,nr)
c$$$         end if
c$$$      end if
      nj_out = nj
      if (verbose.gt.0) then
         write(6,'(A,I0,1X,I0,1X,I0)')
     $        '[finished quicksort_i4_excerpt] nl,nr,nj_out: ' ,nl,nr
     $        ,nj_out
      end if
      end

