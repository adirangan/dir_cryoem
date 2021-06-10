      subroutine additive_noise_k_p(rseed,n_r,grid_k_p_,n_w_,ld_M,M_k_p_
     $     ,sigma)
      implicit none
      integer *4 rseed
      real *8 adi_rand_f
      include 'adi_rand_f_define.f'
      integer n_r,ld_M
      real *8 grid_k_p_(0:0)
      integer n_w_(0:0)
      complex *16 M_k_p_(0:0)
      real *8 sigma
      integer nr,nw,na,nx,n_w_sum
      real *8 r_pre,r_pos,dA
      real *8 pi
      real *8 v1,v2,v3,v4
      complex *16 vE
      real *8 sigma_use
      pi = 4.0d0*atan(1.0d0)
      na=0
      do nr=0,n_r-1
         if (mod(n_w_(nr),2).ne.0) then
            write(6,'(A,I0,A,I0)') 'Warning! odd n_w_(' , nr , '): ' ,
     $           n_w_(nr)
         end if !if (mod(n_w_(nr),2).ne.0) then
         if (nr.eq.0) r_pre = 0.0d0
         if (nr.gt.0) r_pre = grid_k_p_(nr-1)
         r_pos = grid_k_p_(nr)
         dA = pi*(r_pos*r_pos - r_pre*r_pre)/n_w_(nr)
         sigma_use = sigma/dsqrt(dA)
         nx = n_w_(nr)/2
         do nw=0,nx-1
            v1 = adi_rand_f(rseed)
            v2 = adi_rand_f(rseed)
            v3 = sigma_use*dsqrt(-2.0d0*log(v1))*cos(2.0d0*pi*v2)
            v4 = sigma_use*dsqrt(-2.0d0*log(v1))*sin(2.0d0*pi*v2)
            vE = cmplx(v3,v4)
            M_k_p_(na) = M_k_p_(na) + vE
            M_k_p_(na+nx) = M_k_p_(na+nx) + conjg(ve)
            na = na+1
         enddo !do nw=0,nx-1
         na = na+nx
      enddo !do nr=0,n_r-1
      n_w_sum = 0
      do nr=0,n_r-1
         n_w_sum = n_w_sum + n_w_(nr)
      enddo !do nr=0,n_r-1
      if (na.ne.n_w_sum) then
         write(6,'(A,I0,A,I0)') 'Warning! na ' , na , ' n_w_sum ' ,
     $        n_w_sum
      end if !if (na.ne.n_w_sum) then
      if (na.gt.ld_M) then
         write(6,'(A,I0,A,I0)') 'Warning! na ' , na , ' ld_M ' ,
     $        ld_M
      end if !if (na.gt.ld_M) then
      end
