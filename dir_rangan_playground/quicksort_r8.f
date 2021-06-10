!> Doxygen comment: ;\n
!> Quicksort applied to real *8 array Z_. ;\n
!> Sorts in ascending order. ;\n
      subroutine quicksort_r8(nl,nr,Z_,Z_stride,J_,J_stride
     $     ,p_quicksort_r8)
      implicit none
      integer verbose
      data verbose / 0 /
      integer nl,nr,Z_stride,J_stride
      real *8 Z_(0:0)
      integer *4 J_(0:0)
c$$$      real *8 D_(0:nr-nl)
      external p_quicksort_r8
      integer nj,nd
c$$$      character(len=20) format_string
      if (verbose.gt.0) then
         write(6,'(A,I0,1X,I0)') '[entering quicksort_r8] nl,nr:',nl,nr
c$$$         if (nl.le.nr) then
c$$$            write(format_string,'(A,I0,A)') '(',(nr-nl+1),'F4.0)'
c$$$            do nd=0,nr-nl
c$$$               D_(nd) = Z_(nl+nd)
c$$$            enddo
c$$$            write(6,format_string) (D_(nd),nd=nl,nr)
c$$$         end if
      end if
      nj=0
      if (nl.lt.nr) then
         call quicksort_r8_excerpt(nl,nr,Z_,Z_stride,J_,J_stride,nj)
         call p_quicksort_r8(nl,nj-1,Z_,Z_stride,J_,J_stride
     $        ,p_quicksort_r8)
         call p_quicksort_r8(nj+1,nr,Z_,Z_stride,J_,J_stride
     $        ,p_quicksort_r8)
      end if
c$$$      if (verbose.gt.0) then
c$$$         if (nl.le.nr) then
c$$$            write(format_string,'(A,I0,A)') '(',(nr-nl+1),'F4.0)'
c$$$            do nd=0,nr-nl
c$$$               D_(nd) = Z_(nl+nd)
c$$$            enddo
c$$$            write(6,format_string) (D_(nd),nd=nl,nr)
c$$$         end if
c$$$      end if
      if (verbose.gt.0) then
         write(6,'(A,I0,1X,I0)') '[finished quicksort_r8] nl,nr:',nl,nr
      end if
      end

      subroutine quicksort_r8_excerpt(nl,nr,Z_,Z_stride,J_,J_stride
     $     ,nj_out)
      implicit none
      integer verbose
      data verbose / 0 /
      integer nl,nr,Z_stride,J_stride,nj_out
      real *8 Z_(0:0)
      integer *4 J_(0:0)
c$$$      real *8 D_(0:nr-nl)
      real *8 Z_pivot,Z_C
      integer ni,nj,nd,J_C
c$$$      character(len=20) format_string
      if (verbose.gt.0) then
         write(6,'(A,I0,1X,I0)')
     $        '[entering quicksort_r8_excerpt] nl,nr: ',nl,nr
c$$$         if (nl.le.nr) then
c$$$            write(format_string,'(A,I0,A)') '(',(nr-nl+1),'F4.0)'
c$$$            do nd=0,nr-nl
c$$$               D_(nd) = Z_(nl+nd)
c$$$            enddo
c$$$            write(6,format_string) (D_(nd),nd=nl,nr)
c$$$         end if
      end if
      Z_pivot = Z_(nl*Z_stride)
      ni = nl
      nj = nr+1
      if (verbose.gt.0) then
         write(6,'(A,I0,1X,I0,1X,F4.0)') 'ni,nj,Z_pivot: ' ,ni,nj
     $        ,Z_pivot
      end if
      do while (.true.)
         do while (.true.)
            ni = ni+1
            if ((Z_(ni*Z_stride).gt.Z_pivot) .or. (ni.gt.nr)) then
               goto 10
            end if
         enddo
 10      continue
         if (verbose.gt.0) then
            write(6,'(A,I0,1X,F4.0)') 'stopping: ni,Z_(ni): ',ni ,Z_(ni
     $           *Z_stride)
         end if
         do while (.true.)
            nj = nj-1
            if (Z_(nj*Z_stride).le.Z_pivot) then
               goto 20
            end if
         enddo
 20      continue
         if (verbose.gt.0) then
            write(6,'(A,I0,1X,F4.0)') 'stopping: nj,Z_(nj): ',nj ,Z_(nj
     $           *Z_stride)
         end if
         if (ni.ge.nj) then
            goto 30
         end if
         if (verbose.gt.0) then
            write(6,'(A,F4.0,1X,F4.0)') 'ni<nj: swapping ',Z_(ni
     $           *Z_stride),Z_(nj*Z_stride)
         end if
         Z_C = Z_(ni*Z_stride)
         Z_(ni*Z_stride) = Z_(nj*Z_stride)
         Z_(nj*Z_stride) = Z_C
         J_C = J_(ni*J_stride)
         J_(ni*J_stride) = J_(nj*J_stride)
         J_(nj*J_stride) = J_C
      enddo
 30   continue
      if (verbose.gt.0) then
         write(6,'(A,F4.0,1X,F4.0)') 'final: swapping ',Z_(nl *Z_stride)
     $        ,Z_(nj*Z_stride)
      end if
      Z_C = Z_(nl*Z_stride)
      Z_(nl*Z_stride) = Z_(nj*Z_stride)
      Z_(nj*Z_stride) = Z_C
      J_C = J_(nl*J_stride)
      J_(nl*J_stride) = J_(nj*J_stride)
      J_(nj*J_stride) = J_C
c$$$      if (verbose.gt.0) then
c$$$         if (nl.le.nr) then
c$$$            write(format_string,'(A,I0,A)') '(',(nr-nl+1),'F4.0)'
c$$$            do nd=0,nr-nl
c$$$               D_(nd) = Z_(nl+nd)
c$$$            enddo
c$$$            write(6,format_string) (D_(nd),nd=nl,nr)
c$$$         end if
c$$$      end if
      nj_out = nj
      if (verbose.gt.0) then
         write(6,'(A,I0,1X,I0,1X,I0)')
     $        '[finished quicksort_r8_excerpt] nl,nr,nj_out: ' ,nl,nr
     $        ,nj_out
      end if
      end

