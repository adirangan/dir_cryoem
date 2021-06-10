function Y_ = lsqctfsolve(f_out_,polar_a_,azimu_b_,n_out,l_max,n_polar_a,n_order,ctf_weight_,cg_eps);
n_azimu_b = 2*n_polar_a;
polar_a_ = zeros(n_polar_a,1);
weight_polar_a_ = zeros(n_polar_a,1);

      allocate(wts(nquad))
      allocate(sthetas(nquad))
      allocate(twork(nquad))
      call getgridth(nquad,twork,xnodes,sthetas,wts)

      np = (nterms+1)*(nterms+1)
      allocate(icols(nout,k*k))
      allocate(values(nout,k*k))
      allocate(y(np))
      allocate(ugrid(nquadm*nquad))
      allocate(phivals(nout))
      allocate(work(5*np+1000))
c
c     solving neq:   
c
      t1 = second()
ccc      call prinf(' calling mkshinterp nout is *',nout,1)
      call mkshinterp(thetas,phis,nout,k,nquad,nquadm,values,icols)
      t2 = second()
c      write(str,*) 'interp build time =',t2-t1,achar(15)
c      i = mexPrintf(str)
c
      do i = 1,nout
         phivals(i) = dconjg(ctfw(i))*phiout(i)    
      enddo

ccc      call prinf(' calling shinterp_adj nout is *',nout,1)
      call shinterp_adj(nout,k,nquad,nquadm,values,icols,phivals,ugrid)
ccc      call prinf(' calling shevalspherep_adj nout is *',nout,1)
      call shevalspherep_adj(nterms,nquad,nquadm,xnodes,ugrid,y)
ccc      call prin2(' ugrid is *',ugrid,2*nquad*nquadm)
ccc      call prin2(' y is *',y,2*np)
c
      call ccongr2(ier,np,a,multaha,nterms,nquad,nquadm,nout,k,icols,
     1    xnodes,wts,values,phivals,ugrid,ctfw,y,eps,numit,localp,niter,
     1    errs,work)
c      call prinf(' ier *',ier,1)
c$$$      call prinf(' niter *',niter,1)
c$$$      call prin2(' errs *',errs,niter)

c     alex printing to matlab:
c      write(str,*) 'niter =',niter,achar(5)
c      i = mexPrintf(str)
      return
      end

