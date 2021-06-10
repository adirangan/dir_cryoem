      subroutine get_delta_1(N_pixels_in,n_r,half_diameter_x_c,n_delta_x
     $     ,n_delta_y,n_delta_v,delta_x_,delta_y_)
c$$$      This function outputs an array of displacements delta_x_, delta_y_ ;
c$$$      that are supported on a disc of radius ;
c$$$      Rmax = N_pixels*dsqrt(2)*half_diameter_x_c. ;
c$$$      This disc circumscribes a square of side-length ;
c$$$      N_pixels*half_diameter_x_c*2. ;
c$$$      The number of these displacements is stored in n_delta_v. ;
c$$$      We assume that delta_x_ and delta_y_ are preallocated. ;
c$$$      Note: n_delta_v should be close to ceiling(pi*n_delta_x*n_delta_y). ;
      implicit none
      integer verbose
      data verbose / 1 /
      real *8 N_pixels_in,half_diameter_x_c
      integer *4 n_r,n_delta_x,n_delta_y,n_delta_v
      real *8 delta_x_(0:0),delta_y_(0:0)
      integer ndx,ndy,ndv
      real *8 delta_x,delta_y,delta,delta_max
      
      if (verbose.gt.0) then
         write(6,'(A)') '[entering get_delta_1]'
      end if !if (verbose.gt.0) then

      delta_max = N_pixels_in/n_r*half_diameter_x_c*dsqrt(2.0d0)
      if (verbose.gt.0) then
         write(6,'(A,F16.8)') ' N_pixels_in: ' , N_pixels_in
         write(6,'(A,I0)') ' n_r: ' , n_r
         write(6,'(A,F16.8)') ' half_diameter_x_c: ' , half_diameter_x_c
         write(6,'(A,F16.8)') ' delta_max: ' , delta_max
      end if ! if (verbose.gt.0) then

      if ((n_delta_x.le.0) .or. (n_delta_y.le.0)) then
         delta_x_(0) = 0.0d0
         delta_y_(0) = 0.0d0
         n_delta_v = 1
         goto 10
      end if !if ((n_delta_x.le.0) .or. (n_delta_y.le.0)) then

      if ((n_delta_x.le.1) .and. (n_delta_y.le.1)) then
         delta_x_(0) = 0.0d0
         delta_y_(0) = 0.0d0
         n_delta_v = 1
      end if !if ((n_delta_x.le.1) .and. (n_delta_y.le.1)) then

      if ((n_delta_x.le.1) .and. (n_delta_y.gt.1)) then
         ndv=0
         do ndy=0,3*n_delta_y-3
            delta_y = (-3*N_pixels_in + ndy*6*N_pixels_in
     $           /(3*n_delta_y-3))/n_r*half_diameter_x_c
            if (dabs(delta_y).le.delta_max) then
               delta_x_(ndv) = 0.0d0
               delta_y_(ndv) = delta_y
               ndv = ndv+1
            end if !if (dabs(delta_y).le.delta_max) then
         enddo !do ndy=0,3*n_delta_y-1
         n_delta_v = ndv
      end if !if ((n_delta_x.le.1) .and. (n_delta_y.gt.1)) then

      if ((n_delta_x.gt.1) .and. (n_delta_y.le.0)) then
         ndv=0
         do ndx=0,3*n_delta_x-3
            delta_x = (-3*N_pixels_in + ndx*6*N_pixels_in
     $           /(3*n_delta_x-3))/n_r*half_diameter_x_c
            if (dabs(delta_x).le.delta_max) then
               delta_x_(ndv) = delta_x
               delta_y_(ndv) = 0.0d0
               ndv = ndv+1
            end if !if (dabs(delta_x).le.delta_max) then
         enddo !do ndx=0,3*n_delta_x-1
         n_delta_v = ndv
      end if !if ((n_delta_x.gt.1) .and. (n_delta_y.le.0)) then

      if ((n_delta_x.gt.1) .and. (n_delta_y.gt.1)) then
         ndv=0
         do ndy=0,3*n_delta_y-3
            delta_y = (-3*N_pixels_in + ndy*6*N_pixels_in
     $           /(3*n_delta_y-3))/n_r*half_diameter_x_c
            do ndx=0,3*n_delta_x-3
               delta_x = (-3*N_pixels_in + ndx*6*N_pixels_in
     $              /(3*n_delta_x-3))/n_r*half_diameter_x_c
               delta = dsqrt(delta_x**2 + delta_y**2)
               if (dabs(delta).le.delta_max) then
                  delta_x_(ndv) = delta_x
                  delta_y_(ndv) = delta_y
                  ndv = ndv+1
               end if           !if (dabs(delta_x).le.delta_max) then
            enddo               !do ndx=0,3*n_delta_x-1
         enddo                  !do ndy=0,3*n_delta_y-3
         n_delta_v = ndv
      end if !if ((n_delta_x.gt.1) .and. (n_delta_y.le.0)) then
      
      if (verbose.gt.1) then
         do ndv=0,n_delta_v-1
            write(6,'(F8.4,1X,F8.4)') delta_x_(ndv) , delta_y_(ndv)
         enddo !do ndv=0,n_delta_v-1
         write(6,'(A,I0)') ' n_delta_v: ' , n_delta_v
      end if !if (verbose.gt.1) then

      if (verbose.gt.0) then
         write(6,'(A)') '[finished get_delta_1]'
      end if !if (verbose.gt.0) then

 10   continue
      end
      
