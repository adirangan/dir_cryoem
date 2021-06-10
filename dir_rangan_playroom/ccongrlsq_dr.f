
        implicit real *8 (a-h,o-z)
        complex *16 a(1000000),x(1000),y(1000),errs(1000),w(1000000)
c
        call prini(6,13)
C
C       SET ALL PARAMETERS
C
c
        PRINT *, 'ENTER n'
c        READ *,n
        n=12
        CALL PRINF('n=*',n,1)
c
        call congratest(n,a,x,y,errs,w)
c
        stop
        end
c
c
c
c
c     
        subroutine congratest(n,a,x,y,errs,w)
        implicit real *8 (a-h,o-z)
        integer i1,i2,i3,i4,i5,i6
        real *8 p1,p2,p3
        complex *16 a(n,n),x(1),y(1),errs(1),w(1)
        complex *16 ytemp(10 000)
        complex *16 ztemp(1000 000)
        complex *16 yout(10 000)
        complex *16 z1,z2,z3
c
        external multaha
c
        complex *16 ima,zsum
        data ima/(0.0d0,1.0d0)/
c
c
c     ... construct the matrix
c
        do 1400 i=1,n
        do 1200 j=1,n
           a(j,i)=dcmplx(dcos(i*j+0.0d0),dsin(i*j+0.0d0))
ccc           a(j,i)=1/dble(abs(i-j)+1)
           if( i .eq. j ) a(j,i)=a(i,j)+1
 1200 continue
 1400 continue
c
c
c     ... construct the right hand side
c
      do i=1,n
         y(i)=i+ima*i*2
      enddo
      call mult(a,y,ztemp,n)
      call multah(a,y,ytemp,n)
      zsum = 0.0d0
      do i=1,n
        zsum = zsum +dconjg(y(i))*ztemp(i)
      enddo
      write(6,*)' zsum is',zsum
c
c     ... perform conjugate gradient iterations
c
        eps=1d-12
        nitermax=100
c     
c        call prin2('a=*',a,n*n)
c        call corthom(a,n,w,cond)
c        call prin2('cond=*',cond,1)
c        stop
c
        call ccongr2(ier,n,a,multaha,i1,i2,i3,i4,i5,i6,p1,p2,p3,
     1       z1,z2,z3,ytemp,eps,nitermax,x,niter,errs,w)
c
        call prinf('ier=*',ier,1)
        call prinf('niter=*',niter,1)
        call prin2('errs=*',errs,niter)
c
        call prin2('x=*',x,2*n)
        call prin2('y=*',y,2*n)
c
        call mult(a,x,yout,n)
        call prin2('after mult, yout=*',yout,2*n)
c
        d=0
        do i=1,n
           d=d+abs(y(i)-yout(i))**2
        enddo
        d=sqrt(d)
        call prin2('and error=*',d,1)
c
        return
        end
c
c
c
c
c
        subroutine multaha(a,x,y,n,i1,i2,i3,i4,i5,i6,p1,p2,p3,z1,z2,z3)
        implicit real *8 (a-h,o-z)
        integer i1,i2,i3,i4,i5,i6
        real *8 p1,p2,p3
        complex *16 z1,z2,z3
        complex *16 a(n,n),x(1),y(1),d
        complex *16, allocatable :: xtemp(:)
c
        allocate(xtemp(n))
        call mult(a,x,xtemp,n)
        call multah(a,xtemp,y,n)
        return
        end
c
c
c
        subroutine multah(a,x,y,n)
        complex *16 a(n,n),x(1),y(1),d
        do 1400 i=1,n
        d=0
        do 1200 j=1,n
        d=d+dconjg(a(j,i))*x(j)
 1200 continue
        y(i)=d
 1400 continue
        return
        end
c
c
c
c
        subroutine mult(a,x,y,n)
        complex *16 a(n,n),x(1),y(1),d
        do 1400 i=1,n
        d=0
        do 1200 j=1,n
        d=d+a(i,j)*x(j)
 1200 continue
        y(i)=d
 1400 continue
        return
        end
c
c
c
c
c
        subroutine multa2(a,p1,p2,p3,p4,p5,p6,x,y,n)
        complex *16 a(n,n),x(1),y(1),d

        do i=1,n
        y(i)=0
        enddo
        do j=1,n
        do i=1,n
        y(i)=y(i)+a(i,j)*x(j)
        enddo
        enddo

        return
        end
c
c
c
c
c
        subroutine multa3(a,p1,p2,p3,p4,p5,p6,x,y,n)
        complex *16 a(n,n),x(1),y(1),d
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(d)
        do 1400 i=1,n
        d=0
        do 1200 j=1,n
        d=d+a(i,j)*x(j)
 1200 continue
        y(i)=d
 1400 continue
C$OMP END PARALLEL DO
        return
        end
c
c
c
c
c
