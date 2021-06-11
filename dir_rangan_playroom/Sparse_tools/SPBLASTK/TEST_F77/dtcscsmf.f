
c---------------------------------------------------------------
c  Fortran tester for SPBLAS Toolkit
c---------------------------------------------------------------
c
c  Purpose:
c    Demonstrate the calling of the C Sparse BLAS Toolkit
c    routines from within a fortran program.
c---------------------------------------------------------------
c
      program tcscsm
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

      double precision diag(11)
      data diag /.5,-.5,1,1,1,1,1,1,1,-.25,.25/
      double precision a(22)
      data a /1, 2, 8, -2, 1, 1, 4,  3, 1, 3, 1, 4, 1, 5,
     *             1, 4, 3, 1, 1, 1, 1, 1/
      integer pntrb(11)
      data pntrb /1,5, 9,11,13,15,18,19,20,21,22/
      integer pntre(11)
      data pntre /5,9,11,13,15,18,19,20,21,22,23/
      integer indx(22)
      data indx /1,6,10,11,2,6,10,11,3,6,4,6,5,6,6,7,8,7,8,9,10,11/
      double precision a2(22)
      data a2 /1, 1, 1, 1, 1, 2, 1, 3, 4, 5, 1,
     *            4, 1, 3, 1, 1, 8, 4, 1, -2, 3, 1/
      integer pntrb2(11)
      data pntrb2 /1,2,3,4,5, 6,12,14,16,17,20/
      integer pntre2(11)
      data pntre2 /2,3,4,5,6,12,14,16,17,20,23/
      integer indx2(22) 
      data indx2 /1,2,3,4,5,1,2,3,4,5,6,6,7,6,8,9,1,2,10,1,2,11/
      double precision b(22)
      data b /1,2,3,4,5,6,7,8,9,10,11, 1,2,3,4,5,6,7,8,9,10,11/
      double precision c(22)
      data c /1,2,3,4,5,6,7,8,9,10,11, 1,2,3,4,5,6,7,8,9,10,11/
      double precision d(22)
      data d /1,2,3,4,5,6,7,8,9,10,11, 1,2,3,4,5,6,7,8,9,10,11/
      double precision check(22) 
      data check /1,2,3,4,5,6,7,8,9,10,11, 1,2,3,4,5,6,7,8,9,10,11/
      integer m, k, ldb, ldc
      integer i,j
      integer transa, n, unitd, lwork
      integer descra(9)
      integer errcount
      double precision alpha
      double precision invalpha
      double precision beta
      double precision zero
      double precision error
      double precision tolerance
      double precision resid
      double precision work(24)
      lwork = 24
      m=11
      k=11
      ldb=11
      ldc=11
      errcount=0
      zero=0.0
      tolerance=.00001
      
c  Get input: alpha and beta 


      read(5,*) alpha, beta
      if ( alpha .ne. 0 ) then
        invalpha = 1/alpha
      endif

      
      descra(1) = 3
      descra(3) = 0
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
      
            call dcscsm( transa, m, n, unitd, diag, alpha, descra, a,
     *              indx, pntrb, pntre, b, ldb,
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
      
              call dcscmm( transa, m, n, k, invalpha, descra, a,
     *               indx, pntrb, pntre, d, ldb,
     *               zero, check, ldc, work, lwork)
      
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
      
            call dcscsm( transa, m, n, unitd, diag, alpha, descra, a2,
     *              indx2, pntrb2, pntre2, b, ldb,
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
        
      
              call dcscmm( transa, m, n, k, invalpha, descra, a2,
     *                indx2, pntrb2, pntre2, d, ldb,
     *                zero, check, ldc, work, lwork)
      
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
 1000    format(I4,' errors in dtcscsm_f77 run for alpha = ',E10.4,
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

      
