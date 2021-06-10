      program adi_rand_f_dr
c$$$      gfortran -w -o adi_rand_f_dr.out adi_rand_f_dr.f ; ./adi_rand_f_dr.out
      implicit none
      integer n_a,na
      include 'adi_rand_f_define.f'
      integer *4 rseed
      real *8 tmp_d,adi_rand_f
      integer *4, allocatable :: h_(:)
      integer n_h,nh
      
      n_h = 32
      allocate(h_(0:n_h-1))
      call cl1_i4(n_h,h_)
      rseed = 1
      n_a = 1024*32
      do na=0,n_a-1
      tmp_d = adi_rand_f(rseed)
      nh = max(0,min(n_h-1,floor(n_h*tmp_d)))
      h_(nh) = h_(nh) + 1
      enddo !      do na=0,n_a-1
      call write_all_i4(n_h,h_,5,' h_: ')
      stop
      end

      include 'cl1_i4.f'
      include 'write_xxx_xx.f'
      include 'adi_rand_f.f'
