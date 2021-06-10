!> Doxygen comment: ;\n
!> Quicksort applied to complex *16 array Z_. ;\n
!> Sorts in ascending order. ;\n
      subroutine quicksort_c16(nl,nr,Z_,Z_stride,J_,J_stride
     $     ,p_quicksort_c16)
      implicit none
      integer verbose
      data verbose / 0 /
      integer nl,nr,Z_stride,J_stride
      complex *16 Z_(0:0)
      integer *4 J_(0:0)
c$$$      real *8 D_(0:nr-nl)
      external p_quicksort_c16
      integer nj,nd
c$$$      character(len=20) format_string
      if (verbose.gt.0) then
         write(6,'(A,I0,1X,I0)') '[entering quicksort_c16] nl,nr:',nl,nr
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
         call quicksort_c16_excerpt(nl,nr,Z_,Z_stride,J_,J_stride,nj)
         call p_quicksort_c16(nl,nj-1,Z_,Z_stride,J_,J_stride
     $        ,p_quicksort_c16)
         call p_quicksort_c16(nj+1,nr,Z_,Z_stride,J_,J_stride
     $        ,p_quicksort_c16)
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
         write(6,'(A,I0,1X,I0)') '[finished quicksort_c16] nl,nr:',nl,nr
      end if
      end

      subroutine quicksort_c16_excerpt(nl,nr,Z_,Z_stride,J_,J_stride
     $     ,nj_out)
      implicit none
      integer verbose
      data verbose / 0 /
      integer nl,nr,Z_stride,J_stride,nj_out
      complex *16 Z_(0:0)
      integer *4 J_(0:0)
c$$$      real *8 D_(0:nr-nl)
      complex *16 Z_pivot,Z_C
      integer ni,nj,nd,J_C
c$$$      character(len=20) format_string
      if (verbose.gt.0) then
         write(6,'(A,I0,1X,I0)')
     $        '[entering quicksort_c16_excerpt] nl,nr: ',nl,nr
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
         write(6,'(A,I0,1X,I0,1X,F4.0)') 'ni,nj,Z_pivot: '
     $        ,ni,nj,real(Z_pivot)
      end if
      do while (.true.)
         do while (.true.)
            ni = ni+1
            if ((real(Z_(ni*Z_stride)).gt.real(Z_pivot)) .or.
     $           (ni.gt.nr)) then
               goto 10
            end if
         enddo
 10      continue
         if (verbose.gt.0) then
            write(6,'(A,I0,1X,F4.0)') 'stopping: ni,Z_(ni): ',ni
     $           ,real(Z_(ni*Z_stride))
         end if
         do while (.true.)
            nj = nj-1
            if (real(Z_(nj*Z_stride)).le.real(Z_pivot)) then
               goto 20
            end if
         enddo
 20      continue
         if (verbose.gt.0) then
            write(6,'(A,I0,1X,F4.0)') 'stopping: nj,Z_(nj): ',nj
     $           ,real(Z_(nj*Z_stride))
         end if
         if (ni.ge.nj) then
            goto 30
         end if
         if (verbose.gt.0) then
            write(6,'(A,F4.0,1X,F4.0)') 'ni<nj: swapping ',real(Z_(ni
     $           *Z_stride)),real(Z_(nj*Z_stride))
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
         write(6,'(A,F4.0,1X,F4.0)') 'final: swapping ',real(Z_(nl
     $        *Z_stride)),real(Z_(nj*Z_stride))
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
     $        '[finished quicksort_c16_excerpt] nl,nr,nj_out: ' ,nl,nr
     $        ,nj_out
      end if
      end

