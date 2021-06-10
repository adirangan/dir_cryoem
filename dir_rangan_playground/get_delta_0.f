!> Doxygen comment: ;\n
!> Sets up a square grid of displacements. ;\n
!> The displacements range from -T to +T in each direction. ;\n
!> T is equal to N_pixels_in/n_r * half_dimater_x_c. ;\n
!> We should update this to scale by an actual k-value (e.g., grid_k_p_(nr)). ;\n
      subroutine get_delta_0(N_pixels_in,n_r,half_diameter_x_c,n_delta_x
     $     ,delta_x_,n_delta_y,delta_y_)
      implicit none
      integer verbose
      data verbose / 0 /
      real *8 N_pixels_in,half_diameter_x_c
      integer *4 n_r,n_delta_x,n_delta_y
      real *8 delta_x_(0:n_delta_x-1),delta_y_(0:n_delta_y-1)
      integer ndx,ndy
      real *8 delta_x,delta_y
      if (verbose.gt.0) then
         write(6,'(A)') '[entering get_delta_0]'
      end if !if (verbose.gt.0) then
      do ndx=0,n_delta_x-1
         if (n_delta_x.gt.1) then
            delta_x = (-N_pixels_in + ndx*2*N_pixels_in/(n_delta_x-1))
     $           /n_r*half_diameter_x_c
         else
            delta_x = 0.0d0
         end if
         delta_x_(ndx) = delta_x
         do ndy=0,n_delta_y-1
            if (n_delta_y.gt.1) then
               delta_y = (-N_pixels_in + ndy*2*N_pixels_in/(n_delta_y
     $              -1))/n_r*half_diameter_x_c
            else
               delta_y = 0.0d0
            end if
            delta_y_(ndy) = delta_y
         enddo
      enddo      
      if (verbose.gt.0) then
         write(6,'(A)') '[finished get_delta_0]'
      end if !if (verbose.gt.0) then
      end
      
      
