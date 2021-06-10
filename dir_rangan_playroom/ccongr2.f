cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c          this is the end of the debugging code and the beginning 
c          of the actual conjugate residual routine 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
        subroutine ccongr2(ier,n,a,mult,p1,z1,y,eps,numit,x,
     1     niter,err,w)
        implicit real *8 (a-h,o-z)
        real *8 p1
        complex *16 z1
        complex *16 a(1),y(1),x(1),w(1)
        dimension err(1)
c
        external mult
c
c
c        this subroutine solves a  complex linear system with a self-adjoint
c        matrix by means of a conjugate gradient algorithm.
c        the algorithm has been put togethner without proper
c        analysis and with several cludges - user beware.
c
c                      input parameters:
c
c   n - the dimensionality of the system
c   a - the matrix of the system (or whatever parameter the user-supplied
c       matrix-vector multiplication routine mult (see below) requires.
c   mult - the user-supplied matrix-vector multiplication routine.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      description of the calling system of subroutine mult
c
c       the calling sequence of mult is as follows:
c
c          mult(a,x,y,n,p1,z1)
c
c             with the following parameters:
c      a - the matrix of the system (or whatever other parameter) (input)
c      x - the vectot to whichj a is to be applied (input)
c      y -=a(x) (output)
c      n - the dimensionality of a,x,y (input)
c
c      end of description of the calling system of subroutine mult
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   y - the right-hand side of the system to be solved
c   eps - thed accuracy to which the system is to be solved
c   numit - the maximum number of iterations permitted
c
c                       output parameters:
c
c   ier - error return code
c     ier=0 means successful conclusion
c     ier=4 means that the maximum pertmitted number of iterations
c           (numit) has been performed without the desired precision
c           being achieved
c     ier=8 means that the quadratic form failed to decrease before
c           numit iterations have been performed or the accuracy
c           eps has been achieved
c   x - the solution of the system
c   niter - the number of iterations actually performed
c   err - the array of errors produced by the algorithm on various steps,
c           that is err(i) = ||a(x_i)-y||, where x_i is the approximation
c           to the solution obtained on the i-th step of the conjugate
c           gradient process, and || * || is the l^2 norm
c
c                         work array :
c
c   w - must be at least 5*n+10 complex *16 elements long.
c
c        . . . allocate memory for the conjugate gradient scheme
c
        iaxk=1
        iek=iaxk+n
        iaek=iek+n
        iekm1=iaek+n
        iaekm1=iekm1+n
c
c       conduct the conjugate gradient iterations
c
        nrec=1
cccc        nrec=100
c
        call ccongr2a(ier,n,a,y,numit,x,w(iaxk),
     1    w(iek),w(iaek),w(iekm1),w(iaekm1),err,nrec,
     2    mult,p1,z1,eps,niter)
        return
        end
c
c
c
c
c
        subroutine ccongr2a(ier,n,a,y,numit,xk,axk,
     1     ek,aek,ekm1,aekm1,err,nrec,mult,p1,z1,eps,niter)
        implicit real *8 (a-h,o-z)
        complex *16 a(1),y(1),xk(1),axk(1),ek(1),aek(1),
     1     ekm1(1),aekm1(1)
        dimension err(1)
c
        external mult
c
c       initialize the conjugate gradient iterations
c
        done=1
        ifrec=0
        ier=0
        dold=1.0d60
cccc        call prin2('y=*',y,n*2)
c
        do 1400 i=1,n
        xk(i)=0
        axk(i)=0
cccc        xk(i)=y(i)
cccc        xk(i)=i
 1400 continue
c
        call mult(a,xk,axk,n,p1,z1)
c
c       ... find the first direction
c
        er=0
        do 1600 i=1,n
        ek(i)=y(i)-axk(i)
        ekm1(i)=ek(i)
        er=er+dconjg(ek(i))*ek(i)
 1600 continue
c
        er=dsqrt(er)
c
cccc        er=dznrm2(n,ek,1)
c
        err(1)=er
c
        if( er .lt. eps ) return
c
c
c      conduct conjugate gradient iterations
c
        do 4000 k=2,numit
        niter=k
cccc        call prinf('k=*',k,1)
c
        call mult(a,ekm1,aekm1,n,p1,z1)
c
        cd1=0
        cd2=0
        do 2200 i=1,n
        cd1=cd1+dconjg(ek(i))*ek(i)
        cd2=cd2+dconjg(ekm1(i))*aekm1(i)
 2200 continue
c
cccc        cd1=dznrm2(n,ek,1)**2
c
        gamma=cd1/cd2
c
c       ... update the solution
c
        do 2400 i=1,n
        xk(i)=xk(i)+gamma*ekm1(i)
 2400 continue
c
c       ... update the residual
c
        er=0
        do 2600 i=1,n
        ek(i)=ek(i)-gamma*aekm1(i)
        er=er+dconjg(ek(i))*ek(i)
 2600 continue
c
        if( 1 .eq. 2 ) then 
c
        call mult(a,xk,axk,n,p1,z1)
        er=0
        do 2800 i=1,n
        ek(i)=y(i)-axk(i)
        er=er+dconjg(ek(i))*ek(i)
 2800 continue
c
        endif
c
        er=dsqrt(er)
c
cccc        er=dznrm2(n,ek,1)
c
        err(niter)=er
c
c
c       evaluate the quadratic form being minimized
c
        d1=0
        d2=0
        do 3000 i=1,n
cccc        d1=d1+axk(i)*dconjg(xk(i))
        d1=d1+(y(i)-ek(i))*dconjg(xk(i))
        d2=d2+xk(i)*dconjg(y(i))
 3000 continue
cccc        call prin2('d1 as evaluated is*',d1,1)
cccc        call prin2('cd as evaluated is*',cd,2)
        d=d1/2-d2
cccc        call prin2('and quadratic form is*',d,1)
cccc        call prin2('and difference is*',d-dold,1)
c
c       ... if quadratic form failed to decrease, stop
c       
c        if( dold .lt. d ) then
c        ier=8
c        return
c        endif
c
        dold=d
c
c
        if( er .lt. eps ) return
c
c       ... update the search direction
c
        beta=err(niter)/err(niter-1)
        beta=beta**2
c
        do 3200 i=1,n
        ekm1(i)=ek(i)+beta*ekm1(i)
 3200 continue
c
 4000 continue
c
c        the maximum permitted number of iterations has been
c        performed without the desired accuracy being accuracy
c        achieved. inform the user and exit
c
        ier=4
        return
        end
c
c
c
c
c
