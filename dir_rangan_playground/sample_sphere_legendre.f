      subroutine sample_sphere_legendre(
     $     verbose !integer: verbosity level. ;
     $     ,r !real *8: radius of sphere. ;
     $     ,d !real *8: distance between sampled points on equator. ;
     $     ,n_all_0in !integer *4: total number of points. if n_all_0in.le.0 then updated. ;
     $     ,polar_a_all_ !real *8 array (length n_all): polar_a of each point. only calculated if n_all_0in.gt.0 on entry. ;
     $     ,azimu_b_all_ !real *8 array (length n_all): azimu_b of each point. only calculated if n_all_0in.gt.0 on entry. ;
     $     ,weight_all_ !real *8 array (length n_all): quadrature weight of each point. only calculated if n_all_0in.gt.0 on entry. ;
     $     )
      implicit none
      integer verbose           !integer: verbosity level. ;
      real *8 r                 !real *8: radius. ;
      real *8 d                 !real *8: density. ;
      integer *4 n_all_0in      !integer *4: total number of points. ;
      integer *4 n_all_out      !integer *4: total number of points. ;
      real *8 polar_a_all_(0:0) !real *8 array (length n_all): polar_a of each point. ;
      real *8 azimu_b_all_(0:0) !real *8 array (length n_all): azimu_b of each point. ;
      real *8 weight_all_(0:0)  !real *8 array (length n_all): quadrature weight of each point. ;
      integer *4 n_equator !integer: number of sampled points on equator. ;
      integer *4 n_polar_a !integer: number of polar_a. ;
      real *8, allocatable :: legendre_work_(:) !real *8 array (length n_polar_a): scratch space. ;
      real *8, allocatable :: legendre_node_(:) !real *8 array (length n_polar_a): legendre nodes. ;
      real *8, allocatable :: polar_a_(:) !real *8 array (length n_polar_a): polar_a for each ring. ;
      real *8, allocatable :: legendre_weight_(:) !real *8 array (length n_polar_a): legendre weights for each ring. ;
      integer *4 n_azimu_b_ !integer: number of azimu_b for a particular polar_a. ;
      real *8 pi
      integer *4 npolar_a !integer: temporary: index. ;
      real *8 polar_a !real *8: temporary: angle. ;
      integer *4 n_azimu_b !integer: temporary: number of azimu_a for a single ring. ;
      integer *4 nazimu_b !integer: temporary: index. ;
      real *8 azimu_b,dazimu_b !real *8: temporary: angle. ;
      integer *4 nall !integer: temporary: index. ;
      real *8, allocatable :: f_all_(:) !real *8 array (length n_all): temporary: function. ;
      real *8, allocatable :: g_all_(:) !real *8 array (length n_all): temporary: function. ;
      real *8 dotw_r8_f !function output. ;
      real *8 ff_1_(0:0),ff_0_(0:0),gg_1_(0:0),gg_0_(0:0),xx_(0:0) !real *8 temporary: integral of function. ;
      real *8 Il,Ix !real *8 temporary: integrals. ;
      logical flag_memory_checkset
      pi = 4.0d0*atan(1.0d0)
      if (verbose.gt.0) then
         write(6,'(A)') ' [entering sample_sphere_legendre] '
      end if !if (verbose.gt.0) then
c$$$      %%%%%%%%
      if (verbose.gt.0) then
         write(6,'(A)') ' determining size n_all '
      end if !if (verbose.gt.0) then
      n_equator = 3 + nint(2.0d0*pi*r/max(1.0d-6,d))
      n_polar_a = 3 + nint(n_equator/2.0d0)
      allocate(legendre_work_(0:1+n_polar_a-1))
      call cs1_r8(n_polar_a,legendre_work_)
      allocate(legendre_node_(0:1+n_polar_a-1))
      call cs1_r8(n_polar_a,legendre_node_)
      allocate(polar_a_(0:1+n_polar_a-1))
      call cs1_r8(n_polar_a,polar_a_)
      allocate(legendre_weight_(0:1+n_polar_a-1))
      call cs1_r8(n_polar_a,legendre_weight_)
      call gaussq(
     $     1
     $     ,n_polar_a
     $     ,0.0d0
     $     ,0.0d0
     $     ,0
     $     ,0.0d0
     $     ,legendre_work_
     $     ,legendre_node_
     $     ,legendre_weight_
     $     )
      do npolar_a=0,n_polar_a-1
         polar_a_(npolar_a) = dacos(legendre_node_(npolar_a))
      enddo !do npolar_a=0,n_polar_a-1
      n_all_out = 0
      do npolar_a=0,n_polar_a-1
         polar_a = polar_a_(npolar_a)
         n_azimu_b = 3 + nint(2.0d0*pi*dsin(polar_a)*r/d)
         n_all_out = n_all_out + n_azimu_b
      enddo !do npolar_a=0,n_polar_a-1
      if (verbose.gt.0) then
         write(6,'(A,I0)') ' n_all_out: ' , n_all_out
      end if !if (verbose.gt.0) then
c$$$      %%%%%%%%
      flag_memory_checkset=.true.
      call cxs_r8(n_polar_a,legendre_work_,'legendre_work_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_polar_a,legendre_node_,'legendre_node_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_polar_a,polar_a_,'polar_a_',flag_memory_checkset)
      call cxs_r8(n_polar_a,legendre_weight_,'legendre_weight_'
     $     ,flag_memory_checkset)
      if (flag_memory_checkset.eqv..true.) then
         if (verbose.gt.0) then
            write(6,'(A)') '[checkset passed]'
         end if !if (verbose.gt.0) then
      end if !if (flag_memory_checkset.eqv..true.) then
      if (flag_memory_checkset.eqv..false.) then
         write(6,'(A)') '[checkset failed] <-- WARNING'
         stop !exit program due to error.
      end if !if (flag_memory_checkset.eqv..false.) then
c$$$      %%%%%%%%
      if (n_all_0in.le.0) then
         n_all_0in = n_all_out
         goto 10
      end if !if (n_all_0in.le.0) then
c$$$      %%%%%%%%
      if (n_all_0in.gt.0) then
         if (n_all_0in.lt.n_all_out) then
            write(6,'(2(A,I0))') ' Warning, n_all_0in: ' , n_all_0in ,
     $           ' .lt. n_all_out: ' , n_all_out
            stop !program error
         end if !if (n_all_0in.lt.n_all_out) then
         n_all_0in = n_all_out
      end if !if (n_all_0in.gt.0) then
c$$$      %%%%%%%%
      call cl1_r8(n_all_0in,polar_a_all_)
      call cl1_r8(n_all_0in,azimu_b_all_)
      call cl1_r8(n_all_0in,weight_all_)
      nall=0
      do npolar_a=0,n_polar_a-1
         polar_a = polar_a_(npolar_a)
         n_azimu_b = 3 + nint(2.0d0*pi*dsin(polar_a)*r/d)
         dazimu_b = (2.0d0*pi)/(1.0d0*n_azimu_b)
         do nazimu_b=0,n_azimu_b-1
            azimu_b = nazimu_b*dazimu_b
            polar_a_all_(nall) = polar_a
            azimu_b_all_(nall) = azimu_b
            weight_all_(nall) = r**2 * legendre_weight_(npolar_a) *
     $           dazimu_b
            nall = nall+1
         enddo !do nazimu_b=0,n_azimu_b-1
      enddo !do npolar_a=0,n_polar_a-1
      if (nall.ne.n_all_0in) then
         write(6,'(2(A,I0))') ' Warning, n_all_0in: ' , n_all_0in ,
     $        ' .ne. nall: ' , nall
         stop !program error         
      end if !if (nall.ne.n_all_0in) then
c$$$      %%%%%%%%
      if (verbose.gt.-1) then
         write(6,'(A)') ' testing: '
         allocate(f_all_(0:1+n_all_0in-1))
         call cs1_r8(n_all_0in,f_all_)
         allocate(g_all_(0:1+n_all_0in-1))
         call cs1_r8(n_all_0in,g_all_)
         call sample_sphere_f(n_all_0in,polar_a_all_,f_all_)
         call sample_sphere_g(n_all_0in,azimu_b_all_,g_all_)
         Il = dotw_r8_f(n_all_0in,f_all_,g_all_,weight_all_)
         xx_(0) = pi ; call sample_sphere_ff(1,xx_,ff_1_(0))
         xx_(0) = 0.0d0 ; call sample_sphere_ff(1,xx_,ff_0_(0))
         xx_(0) = 2.0d0*pi ; call sample_sphere_gg(1,xx_,gg_1_(0))
         xx_(0) = 0.0d0 ; call sample_sphere_gg(1,xx_,gg_0_(0))
         Ix = r**2 * (ff_1_(0) - ff_0_(0)) * (gg_1_(0) - gg_0_(0))
         write(6,'(2(A,F8.4),2(A,I0),2(A,F8.4),A,F24.16)')
     $         ' r: ' , r
     $        ,' d: ' , d
     $        ,' n_polar_a: ' , n_polar_a
     $        ,' n_all_0in: ' , n_all_0in
     $        ,' Ix: ' , Ix
     $        ,' Il: ' , Il
     $        ,' error: ' , dabs(Ix-Il)/dabs(Ix)
         flag_memory_checkset=.true.
         call cxs_r8(n_all_0in,f_all_,'f_all_',flag_memory_checkset)
         call cxs_r8(n_all_0in,g_all_,'g_all_',flag_memory_checkset)
         if (flag_memory_checkset.eqv..true.) then
            if (verbose.gt.0) then
               write(6,'(A)') '[checkset passed]'
            end if !if (verbose.gt.0) then
         end if !if (flag_memory_checkset.eqv..true.) then
         if (flag_memory_checkset.eqv..false.) then
            write(6,'(A)') '[checkset failed] <-- WARNING'
            stop !exit program due to error.
         end if !if (flag_memory_checkset.eqv..false.) then
      end if !if (verbose.gt.1) then
c$$$      %%%%%%%%
 10   continue
      if (verbose.gt.0) then
         write(6,'(A)') ' [finished sample_sphere_legendre] '
      end if !if (verbose.gt.0) then
      end !subroutine

      subroutine sample_sphere_f(n_polar_a,polar_a_,f_)
      implicit none
      integer *4 n_polar_a
      real *8 polar_a_(0:0)
      real *8 f_(0:0)
      integer *4 m
      integer *4 npolar_a
      real *8 polar_a
      m=6
      do npolar_a=0,n_polar_a-1
         polar_a = polar_a_(npolar_a)
         f_(npolar_a) = dexp(-(-m*dcos(polar_a))**2)*m
      enddo !do npolar_a=0,n_polar_a-1
      end !subroutine

      subroutine sample_sphere_ff(n_polar_a,polar_a_,ff_)
      implicit none
      integer *4 n_polar_a
      real *8 polar_a_(0:0)
      real *8 ff_(0:0)
      integer *4 m
      integer *4 npolar_a
      real *8 polar_a
      real *8 pi
      pi = 4.0d0*datan(1.0d0)
      m=6
      do npolar_a=0,n_polar_a-1
         polar_a = polar_a_(npolar_a)
         ff_(npolar_a) = dsqrt(pi)/2.0d0*derf(-m*dcos(polar_a))
      enddo !do npolar_a=0,n_polar_a-1
      end !subroutine

      subroutine sample_sphere_g(n_azimu_b,azimu_b_,g_)
      implicit none
      integer *4 n_azimu_b
      real *8 azimu_b_(0:0)
      real *8 g_(0:0)
      integer *4 w1,w2
      integer *4 nazimu_b
      real *8 azimu_b
      w1=3; w2=5
      do nazimu_b=0,n_azimu_b-1
         azimu_b = azimu_b_(nazimu_b)
         g_(nazimu_b) = (dcos(w1*azimu_b) + 1.0d0) + (dcos(w2*azimu_b) +
     $        1.0d0)
      enddo !do nazimu_b=0,n_azimu_b-1
      end !subroutine

      subroutine sample_sphere_gg(n_azimu_b,azimu_b_,gg_)
      implicit none
      integer *4 n_azimu_b
      real *8 azimu_b_(0:0)
      real *8 gg_(0:0)
      integer *4 w1,w2
      integer *4 nazimu_b
      real *8 azimu_b
      w1=3; w2=5
      do nazimu_b=0,n_azimu_b-1
         azimu_b = azimu_b_(nazimu_b)
         gg_(nazimu_b) = (dsin(w1*azimu_b)/w1 + azimu_b) + (dsin(w2
     $        *azimu_b)/w2 + azimu_b)
      enddo !do nazimu_b=0,n_azimu_b-1
      end !subroutine
      
