!> Doxygen comment: ;\n
!>      this subroutine solves a  complex linear system with a self-adjoint;\n
!>      matrix by means of a conjugate gradient algorithm.;\n
!>      the algorithm has been put togethner without proper;\n
!>      analysis and with several cludges - user beware.;\n
!>      ;\n
!>                    input parameters:;\n
!>      ;\n
!>  n - the dimensionality of the system;\n
!>  a - the matrix of the system (or whatever parameter the user-supplied;\n
!>      matrix-vector multiplication routine mult (see below) requires.;\n
!>  mult - the user-supplied matrix-vector multiplication routine.;\n
!>      ;\n
!>ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc;\n
!>     description of the calling system of subroutine mult;\n
!>     ;\n
!>      the calling sequence of mult is as follows:;\n
!>      ;\n
!>         mult(a,x,y,n,i1,i2,i3,i4,i5,i6,p1,p2,p3,z1,z2,z3);\n
!>      ;\n
!>            with the following parameters:;\n
!>     a - the matrix of the system (or whatever other parameter) (input);\n
!>     x - the vectot to whichj a is to be applied (input);\n
!>     y -=a(x) (output);\n
!>     n - the dimensionality of a,x,y (input);\n
!>      ;\n
!>     end of description of the calling system of subroutine mult;\n
!>cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc;\n
!>      ;\n
!>  y - the right-hand side of the system to be solved;\n
!>  eps - thed accuracy to which the system is to be solved;\n
!>  numit - the maximum number of iterations permitted;\n
!>     ;\n
!>                      output parameters:;\n
!>      ;\n
!>  ier - error return code;\n
!>    ier=0 means successful conclusion;\n
!>    ier=4 means that the maximum pertmitted number of iterations;\n
!>          (numit) has been performed without the desired precision;\n
!>          being achieved;\n
!>    ier=8 means that the quadrati!>form failed to decrease before;\n
!>          numit iterations have been performed or the accuracy;\n
!>          eps has been achieved;\n
!>  x - the solution of the system;\n
!>  niter - the number of iterations actually performed;\n
!>  err - the array of errors produced by the algorithm on various steps,;\n
!>          that is err(i) = ||a(x_i)-y||, where x_i is the approximation;\n
!>          to the solution obtained on the i-th step of the conjugate;\n
!>          gradient process, and || * || is the l^2 norm;\n
!>      ;\n
!>                        work array :;\n
!>      ;\n
!>  w - must be at least 5*n+10 complex *16 elements long.;\n
!>      ;\n
!>       . . . allocate memory for the conjugate gradient scheme;\n
!>      ;\n
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
        subroutine ccongr2(ier,n,a,mult,i1,i2,i3,i4,i5,i6,p1,p2,p3,
     1             z1,z2,z3,y,eps,numit,x,niter,err,w)
        implicit real *8 (a-h,o-z)
        integer i1,i2,i3,i4,i5,i6
        real *8 p1,p2,p3
        real *8 err(1)
        complex *16 z1,z2,z3
        complex *16 a(1),y(1),x(1),w(1)
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
c          mult(a,x,y,n,i1,i2,i3,i4,i5,i6,p1,p2,p3,z1,z2,z3)
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
c        call prinf(' calling ccongr2a *',n,1)
        call ccongr2a(ier,n,a,y,numit,x,w(iaxk),
     1    w(iek),w(iaek),w(iekm1),w(iaekm1),err,nrec,
     2    mult,i1,i2,i3,i4,i5,i6,p1,p2,p3,z1,z2,z3,eps,niter)
        return
        end
c
c
c
c
c
        subroutine ccongr2a(ier,n,a,y,numit,xk,axk,ek,aek,
     1             ekm1,aekm1,err,nrec,mult,i1,i2,i3,i4,i5,i6,
     2             p1,p2,p3,z1,z2,z3,eps,niter)
        implicit real *8 (a-h,o-z)
        integer i1,i2,i3,i4,i5,i6
        real *8 p1,p2,p3
        real *8 err(1)
        complex *16 a(1),y(1),xk(1),axk(1),ek(1),aek(1),
     1     ekm1(1),aekm1(1),z1,z2,z3
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
ccc        call prinf(' calling mult inside ccongr *',n,1)
        call mult(a,xk,axk,n,i1,i2,i3,i4,i5,i6,p1,p2,p3,z1,z2,z3)
ccc        call prin2(' axk *',axk,2*n)
ccc        call prin2(' y *',y,2*n)
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
        call mult(a,ekm1,aekm1,n,i1,i2,i3,i4,i5,i6,p1,p2,p3,z1,z2,z3)
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
        call mult(a,xk,axk,n,i1,i2,i3,i4,i5,i6,p1,p2,p3,z1,z2,z3)
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
