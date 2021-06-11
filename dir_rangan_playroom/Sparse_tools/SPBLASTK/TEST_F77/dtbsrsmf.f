
c---------------------------------------------------------------
c  Fortran tester for SPBLAS Toolkit
c---------------------------------------------------------------
c
c  Purpose:
c    Demonstrate the calling of the C Sparse BLAS Toolkit
c    routines from within a fortran program.
c---------------------------------------------------------------
c
      program tbsrsm
      implicit none
c
c--------------------------------------------------------------
c  Test the following functions from the Sparse BLAS:
c
c  C <-- alpha*D*inv(A)*B + beta*C
c  C <-- alpha*inv(A)*D*B + beta*C
c  C <-- alpha*D*inv(A')*B + beta*C
c  C <-- alpha*inv(A')*D*B + beta*C
c
c  This program reads two values from standard input,
c  specifying the scalar constants alpha and beta.
c  The underlying code takes different paths if:
c          alpha == 0.0 or 1.0
c      or  beta  == 0.0 or 1.0
c  Therefore, the testing program should be run with various
c  combinations of these parameters.
c
c--------------------------------------------------------------

      double precision dv(36)
      data dv /.5, 0, 0,  0, -.5, 0,  0, 0, .5,
     *       1, 0, 0,  0, 1, 0,  0, 0, 1,
     *       1, 0, 0,  0, 1, 0,  0, 0, 1,
     *       .25, 0, 0,  0, .25, 0,  0, 0, -.25/
      double precision diag(12)
      data diag /.5,-.5,.5,1,1,1,1,1,1,.25,.25,-.25/

c lower triangular part
      double precision a(63)
      data a /1,  0, 0,  0, 1, 0,  0, 0, 1,
     *             1,  4, 7,  2, 5, 8,  3, 6, 9,
     *             1,  0, 0,  0, 1, 0,  0, 0, 1,
     *             1,  0, 0,  0, 1, 0,  0, 0, 1,
     *            10,-13,16, 11,14,17, 12,15,18,
     *            -1,  4, 7,  2,-5, 8,  3, 6,-9,
     *             1,  0, 0,  0, 1, 0,  0, 0, 1/
      integer bindx(7)
      data bindx /1,1,2,3,1,3,4/
      integer bpntrb(4)
      data bpntrb /1,2,4,5/
      integer bpntre(4)
      data bpntre /2,4,5,8/

c  Upper triangular part of A:
      double precision a2(63)
      data a2 /1,  0, 0,  0, 1, 0,  0, 0, 1,
     *             1,  2, 3,  4, 5, 6,  7, 8, 9,
     *             10,11,12,-13,14,15, 16,17,18,
     *             1,  0, 0,  0, 1, 0,  0, 0, 1,
     *             1,  0, 0,  0, 1, 0,  0, 0, 1,
     *            -1,  2, 3,  4,-5, 6,  7, 8,-9,
     *             1,  0, 0,  0, 1, 0,  0, 0, 1 /
      integer bindx2(7)
      data bindx2 /1,2,4,2,3,4,4/
      integer bpntrb2(4)
      data bpntrb2 /1,4,5,7/
      integer bpntre2(4)
      data bpntre2 /4,5,7,8/

      double precision b(24)
      data b /1,2,3,4,5,6,7,8,9,10,11,12, 1,2,3,4,5,6,7,8,9,10,11,12/
      double precision c(24)
      data c /1,2,3,4,5,6,7,8,9,10,11,12, 1,2,3,4,5,6,7,8,9,10,11,12/
      double precision d(24)
      data d /1,2,3,4,5,6,7,8,9,10,11,12, 1,2,3,4,5,6,7,8,9,10,11,12/
      double precision check(24) 
      data check /1,2,3,4,5,6,7,8,9,10,11,12,1,2,3,4,5,6,7,8,9,10,11,12/
      integer mb, kb, m, lb, ldb, ldc
      integer i,j
      integer transa, n, unitd, lwork
      integer descra(9)
      integer errcount
      double precision alpha
      double precision invalpha
      double precision beta
      double precision dzero
      double precision error
      double precision tolerance
      double precision resid
      double precision work(36)
      lwork = 36
      mb=4
      kb=4
      lb=3
      m=12
      ldb=12
      ldc=12
      errcount=0
      dzero=0.0D0
      tolerance=.00001
      
c  Get input: alpha and beta 


      read(5,*) alpha, beta
      if ( alpha .ne. 0 ) then
        invalpha = 1/alpha
      endif
      
      descra(1) = 3
      descra(3) = 1
      descra(4) = 1
      descra(5) = 1
      
      
      print *,'-----------------------------------------------------'
      print *,'  alpha = ', alpha, ' beta = ', beta
      print *,'-----------------------------------------------------'
c  Loop on columns of C (test vector and matrix routines)
      do 10 n = 1,2
            print *,'*** n = ', n, ' ***'

c test non-transpose and transpose  
        do 20 transa=0,1     
          print *,'   << transa = ', transa, ' >>'
c test identity, left and right scaling 
          do 30 unitd=1,3
            print *,'      ++ unitd = ',unitd,' ++'
      
c lower triangular matrix 
            descra(2) = 1                      
            print *,'          -- lower triangular --'
      
c  Initialize C:
            do  40 i = 1,m
              do  50 j = 1,n
                c( (j-1)*m + i ) = i
  50          continue
  40        continue

      
c Call triangular solve with lower triangular matrix     
      
            call dbsrsm( transa, mb, n, unitd, dv, alpha, descra, a,
     *              bindx, bpntrb, bpntre, lb, b, ldb,
     *              beta, c, ldc, work, lwork)
      
      
c  Backtrack from solution using matrix multiply; after   
      
c  calculation, "check" should match "b"                  
      
            do 60 i=1,n*m
              d(i) = c(i) - beta * b(i)
 60         continue
        
            if ( alpha .ne. 0 ) then
              if ( unitd .eq. 2 ) then
                do 70 i=1,m
                  do 80 j=1,n
                    d((j-1)*m+i) =  d((j-1)*m+i) / diag(i)
 80               continue
 70             continue
              endif
      
              call dbsrmm( transa, mb, n, kb, invalpha, descra, a,
     *               bindx, bpntrb, bpntre, lb, d, ldb,
     *               dzero, check, ldc, work, lwork)
      
              if ( unitd .eq. 3 ) then
                do 90 i=1,m
                  do 100 j=1,n
                    check((j-1)*m+i) =  check((j-1)*m+i) / diag(i)
 100              continue
  90            continue
              endif
              error = resid(n*m, check, b)
            else
              error = 0
              do 110 i=1,n*m
                 check(i) = d(i)
                 error = error + dabs(d(i))
 110          continue
              error = error / n*m
            endif
            if ( error .ge. tolerance ) then
               errcount =  errcount + 1          
               print *,'Error for lower triangular solve with '
               print *,'n = ', n, ', transa = ', transa, ', unitd = ',
     *                unitd
               print *,'Residual: ', error
               do 120 i=1,n*m
                 print *, c(i), check(i)  
 120           continue
            endif
      
            descra(2) = 2
c  upper triangular matrix 
            print *,'          -- upper triangular --'
       
c  Initialize C:
            do  130 i = 1,m
              do  140 j = 1,n
                c( (j-1)*m + i ) = i
  140         continue
  130       continue
      
      
c  Call triangular solve with upper triangular matrix: 
      
            call dbsrsm( transa, mb, n, unitd, dv, alpha, descra, a2,
     *              bindx2, bpntrb2, bpntre2, lb, b, ldb,
     *              beta, c, ldc, work, lwork)
      
      
c  Backtrack from solution using matrix multiply; after   
      
c  calculation, "check" should match "b"                  
      
            do 160 i=1,n*m
              d(i) = c(i) - beta * b(i)
 160        continue
        
            if ( alpha .ne. 0 ) then
              if ( unitd .eq. 2 ) then
                do 170 i=1,m
                  do 180 j=1,n
                    d((j-1)*m+i)  =  d((j-1)*m+i) / diag(i)
 180              continue
 170            continue
              endif
        
      
              call dbsrmm( transa, mb, n, kb, invalpha, descra, a2,
     *                bindx2, bpntrb2, bpntre2, lb, d, ldb,
     *                dzero, check, ldc, work, lwork)
      
              if ( unitd .eq. 3 ) then
                do 190 i=1,m
                  do 200 j=1,n
                    check((j-1)*m+i) = check((j-1)*m+i) / diag(i)
 200              continue
 190            continue
              endif
              error = resid(n*m, check, b)
            else
              error = 0
              do 210 i = 1,n*m
                error = error + dabs(d(i))
                check(i) = d(i)
 210          continue
              error = error / n*m
            endif
            if ( error .ge. tolerance ) then
               errcount = errcount + 1
               print *,'Error for upper triangular solve with'
               print *,'n = ', n, ', transa = ', transa, ', unitd = ',
     *                unitd
               print *,'Residual: ', error
               do 220 i=1,n*m
                 print *, c(i), check(i)  
 220           continue
            endif
      
 30       continue
 20     continue                      
 10   continue
      
      if (errcount .gt. 0) then
         write(6,1000) errcount, alpha, beta
 1000    format(I4,' errors in dtbsrsm_f77 run for alpha = ',E10.4,
     *                                            ' beta = ',E10.4,/)
      endif


      call exit(errcount)
      stop
      end


      double precision function resid(m, x1, x2)
      integer m, i
      double precision x1(*), x2(*)
      double precision norm
      norm = 0.0

      do 10 i=1,m
         norm = norm + dabs(x1(i) - x2(i))
 10   continue

      if ( norm .lt. 0.0 .or. m .eq. 0)  then
       resid  = norm
      else
       resid = norm/m
      endif

      return
      end

      
