!> Doxygen comment: ;\n
!> multiplies by transpose(A) in the least-squares setup ;\n
!> calls 'multah' ;\n
      subroutine test_residual_multah(thetas,phis,nout,residua,nterms,
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
      complex *16 residua(*)
      complex *16 phiout(nout)
      complex *16 ctfw(nout)
      complex *16 output(nout)
      complex *16, allocatable :: ugrid(:)
      complex *16, allocatable :: phivals(:)
      complex *16 a,z1,z2,z3
      integer i
      character(len=80) str

      if (verbose.gt.0) then
         write(6,'(A,I0)') '[entering test_residual_multah], nout: ' , nout
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
      call multah(a,residua,output,np,nterms,nquad,nquadm,nout,k,icols,
     1    xnodes,wts,values,phivals,ugrid,ctfw)

      if (verbose.gt.0) then
         write(6,'(A)') '[finished test_residual_multah]'
      end if

      end
c
c

cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c    This does the A* apply. 
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine multah(a,x,y,n,nterms,nquad,nquadm,nout,k,icols, xnodes
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
      do i = 1,nout
         phiout(i) = dconjg(ctfw(i))*x(i)    
      enddo

      call shinterp_adj(nout,k,nquad,nquadm,values,icols,phiout
     $     ,ugrid)
      call shevalspherep_adj(nterms,nquad,nquadm,xnodes,ugrid,y)

c
      return
      end
c
c
