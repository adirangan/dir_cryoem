!> Doxygen comment: ;\n
!> multiplies by (A) in the least-squares setup ;\n
!> calls 'multa' ;\n
      subroutine test_residual_multa(thetas,phis,nout,localp,nterms,
     1           nquad,k,ctfw,output)

      implicit real *8 (a-h,o-z)
      integer verbose
      data verbose / 0 /
      integer nterms,nquad,nquadm
      integer, allocatable :: icols(:,:)
      real *8, allocatable :: xnodes(:),wts(:),twork(:),sthetas(:)
      real *8 thetas(nout)
      real *8 phis(nout)
      real *8, allocatable :: values(:,:)
      complex *16 localp((nterms+1)*(nterms+1))
      complex *16 phiout(nout)
      complex *16 ctfw(nout)
      complex *16 output(nout)
      complex *16, allocatable :: ugrid(:)
      complex *16, allocatable :: phivals(:)
      complex *16 a,z1,z2,z3
      integer i
      character(len=80) str

      if (verbose.gt.0) then
         write(6,'(A,I0)') '[entering test_residual_multa], nout: ' , nout
      end if

      nquadm = 2*nquad
      allocate(xnodes(nquad))
      allocate(wts(nquad))
      allocate(sthetas(nquad))
      allocate(twork(nquad))
      call getgridth(nquad,twork,xnodes,sthetas,wts)

      np = (nterms+1)*(nterms+1)
      allocate(icols(nout,k*k))
      allocate(values(nout,k*k))
      allocate(ugrid(nquadm*nquad))
      allocate(phivals(nout))

      call mkshinterp(thetas,phis,nout,k,nquad,nquadm,values,icols)
      call multa(a,localp,output,np,nterms,nquad,nquadm,nout,k,icols,
     1    xnodes,wts,values,phivals,ugrid,ctfw)

      if (verbose.gt.0) then
         write(6,'(A)') '[finished test_residual_multa]'
      end if

      end
c
c

cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c    This does the A apply. 
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine multa(a,x,y,n,nterms,nquad,nquadm,nout,k,icols, xnodes
     $     ,wts,values,phiout,ugrid,ctfw)
      implicit real *8 (a-h,o-z)
c$$$      real *8 icols(nout,k*k)
      integer *4 icols(nout,k*k)
      real *8 xnodes(nquad)
      real *8 wts(nquad)
      real *8 values(nout,k*k)
      real *8 p3
      complex *16 ugrid(nquadm*nquad)
      complex *16 phiout(nout),z3
      complex *16 a,x(n),y(n),d
      complex *16 ctfw(nout)
c
ccc      call prinf(' calling shevalspherep *',n,1)
      call shevalspherep(x,ugrid,nterms,nquad,nquadm,xnodes)
ccc      call prinf(' calling shinterp *',n,1)
      call shinterp(nout,k,nquad,nquadm,values,icols,ugrid,phiout)
c
      do i = 1,nout
c$$$         phiout(i) = dconjg(ctfw(i))*ctfw(i)*phiout(i)    
         y(i) = ctfw(i)*phiout(i)    
      enddo
c
      return
      end
c
c
