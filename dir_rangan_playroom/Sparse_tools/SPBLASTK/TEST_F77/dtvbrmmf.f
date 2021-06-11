c---------------------------------------------------------------
c  Fortran tester for SPBLAS Toolkit 
c---------------------------------------------------------------
c
c  Purpose:
c    Demonstrate the calling of the C Sparse BLAS Toolkit 
c    routines from within a fortran program.
c---------------------------------------------------------------
c
      program tbsrmm
      implicit none
c
c--------------------------------------------------------------
c  Test the following functions from the Sparse BLAS:         
c                                                             
c  C <-- alpha*A*B + beta*C                                   
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

c  All of matrix A:
      double precision a(51)
      data a / 1, 0,  0, 1, 1, 2, -1, 0,  1, -1,
     *      1, 0, 0,  0, 1, 0,  0, 0, 1, 3, 4, 5,
     *      1, 2, 3, 4, 5, 1, 6, 7, 8,
     *      6, 7, 8, 1, 0, 0,  0, 1, 0,  0, 0, 1,
     *      -1, 1,  0, -1, 1, 0,  0, 1/

c  All of skew matrix A: 
c     (A with diagonal zeroed out and upper triangular part negated)
      double precision ka(51)
      data ka / 0, 0,  0, 0, -1, -2, 1, 0,  -1, 1, 
     *      0, 0, 0,  0, 0, 0,  0, 0, 0,
     *      -3, -4, -5, 1, 2, 3, 4, 5, 0,
     *      -6, -7, -8, 6, 7, 8, 0, 0, 0,  0, 0, 0,  0, 0, 0,
     *      -1, 1,  0, -1, 0, 0,  0, 0 /
      integer indx(14)
      data indx /1,5,7,11,20,23,25,28,29,32,35,44,48,52/
      integer bindx(13)
      data bindx /1,3,5,2,3,1,2,3,4,3,4,1,5/
      integer rpntr(6)
      data rpntr /1,3,6,7,10,12/
      integer cpntr(6)
      data cpntr /1,3,6,7,10,12/
      integer bpntrb(5) 
      data bpntrb /1,4,6,10,12/
      integer bpntre(5)
      data bpntre /4,6,10,12,14/
      
c lower triangular part 
      double precision la(39)
      data la / 1, 0,  0, 1,
     *      1, 0, 0,  0, 1, 0,  0, 0, 1,
     *      1, 2, 3, 4, 5, 1, 6, 7, 8,
     *      1, 0, 0,  0, 1, 0,  0, 0, 1,
     *     -1, 1,  0, -1, 1, 0,  0, 1/
      integer lindx(10)
      data lindx /1,5,14,16,19,20,23,32,36,40/
      integer lbindx(9)
      data lbindx /1,2,1,2,3,3,4,1,5/
      integer lrpntr(6)
      data lrpntr /1,3,6,7,10,12/
      integer lcpntr(6)
      data lcpntr /1,3,6,7,10,12/
      integer lbpntrb(5)
      data lbpntrb /1,2,3,6,8/
      integer lbpntre(5)
      data lbpntre /2,3,6,8,10/

c  Upper triangular part of A:
      double precision ua(39)
      data ua / 1, 0,  0, 1, 1, 2, -1, 0,  1, -1,
     *       1, 0, 0,  0, 1, 0,  0, 0, 1,
     *       3, 4, 5, 1, 6, 7, 8,
     *       1, 0, 0,  0, 1, 0,  0, 0, 1, 1, 0, 0, 1/


      integer uindx(10)
      data uindx  /1,5,7,11,20,23,24,27,36,40/
      integer ubindx(9)
      data ubindx /1,3,5,2,3,3,4,4,5/
      integer urpntr(6)
      data urpntr /1,3,6,7,10,12/
      integer ucpntr(6)
      data ucpntr /1,3,6,7,10,12/
      integer ubpntrb(5)
      data ubpntrb /1,4,6,8,9/
      integer ubpntre(5)
      data ubpntre /4,6,8,9,10/

      
      
      double precision b(22)
      data b /1,2,3,4,5,6,7,8,9,10,11,
     *        1,2,3,4,5,6,7,8,9,10,11/
      double precision c(22)
      data c /1,2,3,4,5,6,7,8,9,10,11,
     *        1,2,3,4,5,6,7,8,9,10,11/
      double precision d(22)
      data d /1,2,3,4,5,6,7,8,9,10,11,
     *        1,2,3,4,5,6,7,8,9,10,11/
      double precision check(22)
      data check /1,2,3,4,5,6,7,8,9,10,11,
     *            1,2,3,4,5,6,7,8,9,10,11/
      integer mb, kb, m, ldb, ldc
c Begin description of rectangular matrix 
      double precision ra(62)
      data ra / 1, 0,  0, 1, 1, 2, -1, 0,  1, -1,
     *      1, 1, 1, 0, 0,  0, 1, 0,  0, 0, 1,
     *      3, 4, 5, 1, 1, 1, 1, 2, 3, 4, 5,
     *      1, 6, 7, 8, 1, 6, 7, 8,
     *      1, 0, 0,  0, 1, 0,  0, 0, 1,
     *      1, 1, 1, -1, 1,  0, -1, 1, 0,  0, 1, 1, 1/
      integer rindx(18)
      data rindx /1,5,7,11,13,22,25,28,30,33,34,37,38,41,50,53,57,61/
      integer rbindx(18)
      data rbindx /1,3,5,6,2,3,6,1,2,3,4,6,3,4,6,1,5,6/
      integer rrpntr(6)
      data rrpntr /1,3,6,7,10,12/
      integer rcpntr(7)
      data rcpntr  /1,3,6,7,10,12,13/
      integer rbpntrb(5)
      data rbpntrb /1,5,8,13,16/
      integer rbpntre(5)
      data rbpntre /5,8,13,16,19/

      double precision rb(24)
      data rb /1,2,3,4,5,6,7,8,9,10,11,1,
     *         1,2,3,4,5,6,7,8,9,10,11,1/
      double precision rc(24)
      data rc /1,2,3,4,5,6,7,8,9,10,11,1,
     *         1,2,3,4,5,6,7,8,9,10,11,1/
      double precision rd(24)
      data rd /1,2,3,4,5,6,7,8,9,10,11,1,
     *         1,2,3,4,5,6,7,8,9,10,11,1/
      double precision rcheck(24)
      data rcheck /1,2,3,4,5,6,7,8,9,10,11,1,
     *             1,2,3,4,5,6,7,8,9,10,11,1/
      double precision rsumb
      integer rmb, rkb, rlb, rm, rldb, rldc
      integer i,j
      integer transa, n, lwork
      integer descra(9)
      integer errcount
      double precision alpha, malpha
      double precision beta
      double precision dzero
      double precision error
      double precision tolerance
      double precision resid
      double precision work(24)
      mb=5
      kb=5
      m=11
      ldb=11
      ldc=11
      rmb=5
      rkb=6
      rlb=3
      rm=11
      rldb=12
      rldc=12
      rsumb = 66
      errcount = 0
      dzero = 0.0
      tolerance = .00001

c  Get input: alpha and beta */


      read(5,*) alpha, beta
      malpha = -1.0 * alpha
      
      descra(3) = 0
      descra(4) = 1
      descra(5) = 1
      
      
      print *,'-----------------------------------------------------'
      print *,'  alpha = ', alpha, ' beta = ', beta
      print *,'-----------------------------------------------------'

c  Loop on columns of C (test vector and matrix routines)
      do 10 n = 1,2
            print *,'*** n = ', n, ' ***' 
            print *, '   General matrices:'
           
c  Testing rectangular matrices 
            print *, '      rectangular'
            descra(1) = 0
            descra(2) = 1 
 
c  Initialize C:
            do  20 i = 1,m
              do  30 j = 1,n
                c( (j-1)*m + i ) = i
  30          continue
  20        continue     
      
            transa = 0
            call dvbrmm( transa, rmb, n, rkb, alpha, descra, ra,
     *        rindx, rbindx, rrpntr, rcpntr, rbpntrb,
     *        rbpntre, rb, rldb,
     *        beta, c, ldc, work, lwork)

            do 35 i= 1,n*m
              d(i) = c(i) - alpha
  35        continue

c  Initialize C:
            do  40 i = 1,m
              do  50 j = 1,n
                c( (j-1)*m + i ) = i
  50          continue
  40        continue     

c  Call mat-mult with explicit symmtric matrix           

            transa = 0
            call dvbrmm( transa, mb, n, kb, alpha, descra, a,
     *         indx, bindx, rpntr, cpntr, bpntrb, bpntre, b, ldb,
     *         beta, c, ldc, work, lwork)

            error = resid(n*m, d, c)
            if ( error .ge. tolerance ) then
               errcount = errcount + 1
               print *, 'Error for rectangular matmult (no transpose)'
               print *, 'n = ', n
               print *, 'Residual: ', error          
               do 55 i=1,n*m
                print *, d(i), c(i) 
 55            continue
             endif

c  Initialize rc:
             do 60 i = 1,m
              do 70 j = 1,n
                rc( (j-1)*(m+1) + i ) = i
  70          continue
  60        continue     
            do 66 j = 1,n
              rc( j*(m+1)) = 1
  66        continue
      
            transa = 1
            call dvbrmm( transa, rmb, n, rkb, alpha, descra, ra,
     *        rindx, rbindx, rrpntr, rcpntr, rbpntrb, rbpntre, b, ldb,
     *        beta, rc, rldc, work, lwork)

            error = resid(m, c, rc)
            error = error + alpha*rsumb + beta - rc(m+1)
            if ( error .ge. tolerance ) then
               errcount = errcount + 1
               print *, 'Error for rectangular matmult (transpose)'
               print *, 'n = ', n
               print *, 'Residual: ', error          
               do 75 i=1,m
                print *, c(i), rc(i) 
 75            continue
               print *, alpha*rsumb + beta, rc(m+1)
            endif
     
            descra(1) = 0
            descra(2) = 1
            print *, '      lower triangular'

c  Initialize C:
            do 100 i = 1,m
              do 110 j = 1,n
                c( (j-1)*m + i ) = i
 110          continue
 100        continue     
      
c  Call triangular mat-mult with lower triangular matrix:
            transa = 0
            call dvbrmm( transa, mb, n, kb, alpha, descra, la,
     *              lindx, lbindx, lrpntr, lcpntr, lbpntrb, 
     *              lbpntre, b, ldb,
     *              beta, c, ldc, work, lwork)
      
c  Save result in vector d:
            do 120 i = 1, n*m
              d(i) = c(i)
 120        continue
        
            descra(2) = 2 
            print *, '      upper triangular'
       
c  Initialize C:
            do 130 i = 1,m
              do 140 j = 1,n
                c( (j-1)*m + i ) = i
 140          continue
 130        continue     
      
c  Call triangular mat-mult with upper triangular matrix:
            transa = 1
            call dvbrmm( transa, mb, n, kb, alpha, descra, ua,
     *              uindx, ubindx, urpntr, ucpntr, ubpntrb, 
     *              ubpntre, b, ldb,
     *              beta, c, ldc, work, lwork)
      
c  Compare upper and lower triangular results:
            error = resid(n*m, d, c)
            if ( error .ge. tolerance ) then
               errcount = errcount +1
               print *, 'Error for upper(or lower) general matmult'
               print *, 'n = ', n
               print *, 'Residual: ', error          
               do 150 i=1,n*m
                print *, d(i), c(i) 
 150           continue
            endif

            print *, '   Symmetric matrices:'
      
c  First, calculate solution with explicit symmetric matrix:
            descra(1) = 0
      
c  Initialize C:
            do 160 i = 1,m
              do 170 j = 1,n
                c( (j-1)*m + i ) = i
 170          continue
 160        continue     
      
            transa = 0
            call dvbrmm( transa, mb, n, kb, alpha, descra, a,
     *              indx, bindx, rpntr, cpntr, bpntrb, 
     *              bpntre, b, ldb,
     *              beta, c, ldc, work, lwork)
      
c  Save result to vector d:
            do 180 i = 1, n*m
              d(i) = c(i)
 180        continue
        
            descra(1) = 1
            descra(2) = 1
            print *, '      lower triangular'
      
c  Initialize C:
            do 200 i = 1,m
              do 210 j = 1,n
                c( (j-1)*m + i ) = i
 210          continue
 200        continue     
      
c  Call symmetric mat-mult with lower triangular matrix:
            transa = 0
            call dvbrmm( transa, mb, n, kb, alpha, descra, la,
     *              lindx, lbindx, lrpntr, lcpntr, lbpntrb, 
     *              lbpntre, b, ldb,
     *              beta, c, ldc, work, lwork)
      
c  Compare explicit and implicit results:
            error = resid(n*m, d, c)
            if ( error .ge. tolerance ) then
               errcount = errcount + 1
               print *, 'Error for symmetric matmult (lower triangular)'
               print *, 'n = ', n
               print *, 'Residual: ', error           
               do 220 i=1,n*m
                print *, d(i), c(i) 
 220           continue
            endif
      
        
            descra(2) = 2
            print *, '      upper triangular'
      
c  Initialize C:
            do 230 i = 1,m
              do 240 j = 1,n
                c( (j-1)*m + i ) = i
 240          continue
 230        continue     
      
c  Call symmetric mat-mult with upper triangular matrix:
            transa = 0
            call dvbrmm( transa, mb, n, kb, alpha, descra, ua,
     *              uindx, ubindx, urpntr, ucpntr, ubpntrb, 
     *              ubpntre, b, ldb,
     *              beta, c, ldc, work, lwork)

c  Compare explicit and implicit results:
            error = resid(n*m, d, c)
            if ( error .ge. tolerance ) then
               errcount = errcount + 1
               print *, 'Error for symmetric matmult (upper triangular)'
               print *, 'n = ', n
               print *, 'Residual: ', error
               do 250 i=1,n*m
                print *, d(i), c(i) 
 250           continue
            endif
       
            print *, '   Skew-Symmetric matrices:'

c  First, calculate solution with explicit skew-symmetric matrix:
            descra(1) = 0
      
c  Initialize C:
            do 260 i = 1,m
              do 270 j = 1,n
                c( (j-1)*m + i ) = i
 270          continue
 260        continue     
      
            transa = 0
            call dvbrmm( transa, mb, n, kb, alpha, descra, ka,
     *              indx, bindx, rpntr, cpntr, bpntrb, 
     *              bpntre, b, ldb,
     *              beta, c, ldc, work, lwork)
      
            do 280 i = 1, n*m
              d(i) = c(i)
 280        continue
        
            descra(1) = 4
            descra(2) = 1
            print *, '      lower triangular (no transp)'
      
c  Initialize C:
            do 300 i = 1,m
              do 310 j = 1,n
                c( (j-1)*m + i ) = i
 310          continue
 300        continue     
      
c  Call skew-symmetric mat-mult with lower triangular matrix:
            transa = 0
            call dvbrmm( transa, mb, n, kb, alpha, descra, la,
     *              lindx, lbindx, lrpntr, lcpntr, lbpntrb, 
     *              lbpntre, b, ldb,
     *              beta, c, ldc, work, lwork)
      
c  Compare explicit and implicit results:
            error = resid(n*m, d, c)
            if ( error .ge. tolerance ) then
               errcount = errcount + 1
         print *, 'Error for skew-symmetric matmult (lower triangular)'
               print *, 'n = ', n
               print *, 'Residual: ', error
               do 320 i=1,n*m
                print *, d(i), c(i) 
 320           continue
            endif 
        
            descra(2) = 2
            print *, '      upper triangular (transp)'
      
c  Initialize C:
            do 330 i = 1,m
              do 340 j = 1,n
                c( (j-1)*m + i ) = i
 340          continue
 330        continue     
      
c  Call skew-symmetric mat-mult with upper triangular matrix:
            transa = 1
            call dvbrmm( transa, mb, n, kb, alpha, descra, ua,
     *              uindx, ubindx, urpntr, ucpntr, ubpntrb, 
     *              ubpntre, b, ldb,
     *              beta, c, ldc, work, lwork)
      
c  Compare explicit and implicit results:
            error = resid(n*m, d, c)
            if ( error .ge. tolerance ) then
               errcount = errcount + 1
         print *, 'Error for skew-symmetric matmult (upper triangular)'
               print *, 'n = ', n
               print *, 'Residual: ', error          
               do 350 i = 1, n*m
                print *, d(i), c(i) 
 350           continue
            endif
      
c  Now, work with transp of lower and upper triangular matrix,   
c  results should be negation of explicit matrix multiply    
c  check by taking alpha = -alpha  (malpha)

            descra(2) = 1
            print *, '      lower triangular (transp)'
      
c  Initialize C:
            do 360 i = 1,m
              do 370 j = 1,n
                c( (j-1)*m + i ) = i
 370          continue
 360        continue     
      
            transa = 1
            call dvbrmm( transa, mb, n, kb, malpha, descra, la,
     *              lindx, lbindx, lrpntr, lcpntr, lbpntrb, 
     *              lbpntre, b, ldb,
     *              beta, c, ldc, work, lwork)
      
c  Compare explicit and implicit results:
            error = resid(n*m, d, c)
            if ( error .ge. tolerance ) then
               errcount = errcount + 1
       print *, 'Error for skew-symmetric matmult (lower triangular)'
               print *, 'n = ', n 
               print *, 'Residual: ',error
               do 380 i = 1, n*m
                print *, d(i), c(i) 
 380           continue      
            endif
        
            descra(2) = 2
            print *, '      upper triangular (no transp)'
      
c  Initialize C:
            do 400 i = 1,m
              do 410 j = 1,n
                c( (j-1)*m + i ) = i
 410          continue
 400        continue     
      
            transa = 0

            call dvbrmm( transa, mb, n, kb, malpha, descra, ua,
     *              uindx, ubindx, urpntr, ucpntr, ubpntrb, 
     *              ubpntre, b, ldb,
     *              beta, c, ldc, work, lwork)
      
c  Compare explicit and implicit results:
            error = resid(n*m, d, c)
            if ( error .ge.  tolerance ) then
               errcount = errcount + 1
        print *, 'Error for skew-symmetric matmult (upper triangular)'
               print *, 'n = ', n
               print *, 'Residual: ', error          
               do 420 i = 1, n*m
                print *, d(i), c(i)
 420           continue
            endif
      
 10   continue
      
      if (errcount .gt. 0) then
         write(6,1000) errcount, alpha, beta
 1000    format(I4,' errors in dtvbrmm_f77 run for alpha = ',E10.4,
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

      
