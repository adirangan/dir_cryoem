      program test_chebexps
      implicit none
      integer n
      integer itype
      real *8 alpha
      real *8 beta
      real *8 x(8)
      real *8 w(8)
      integer j
      n = 8
      itype = 1
      alpha = 0.0d0
      beta = 0.0d0
      call chebexps(itype,n,x,alpha,beta,w)
      write(6,'(A,8(F8.4,1X))') ' nodes  : ' , (x(j),j=1,8)
      write(6,'(A,8(F8.4,1X))') ' weights: ' , (w(j),j=1,8)
      end !program
c$$$ type-1 I believe (i.e., chebpts(8,1)) ;
c$$$  nodes  :      -0.9807852804032304      -0.8314696123025453      -0.5555702330196020      -0.1950903220161282       0.1950903220161283       0.5555702330196023       0.8314696123025452       0.9807852804032304
c$$$  weights:       0.0669829456985898       0.2229879330145787       0.3241525190645244       0.3858766022223070       0.3858766022223071       0.3241525190645244       0.2229879330145788       0.0669829456985898
      stop

      include 'chebexps.f'
      include 'prini.f'
