!> Doxygen comment: ;\n
!> Decenters complex *16. ;\n
!> Note that this is only its own inverse when n_x and n_y are even! ;\n
!> Otherwise, use recenter_c16 to invert. ;\n
!> a+b = N
!> 0,1,2,...,a-1,a,a+1,...,N-1 --> a,a+1,a+2,...,a+a-1,0,1,...,n-a-1 ;\n
      subroutine decenter_c16(n_x,n_y,a_,b_)
      integer n_x,n_y,nx,ny,nxt,nyt
      complex *16 a_(0:n_x*n_y-1),b_(0:n_x*n_y-1)
      complex *16, allocatable :: c_(:)
      integer n_x_l,n_y_l
      integer n_x_u,n_y_u
      n_x_l = n_x/2
      n_y_l = n_y/2
      n_x_u = n_x - n_x_l
      n_y_u = n_y - n_y_l
      allocate(c_(0:n_x*n_y-1))
      do ny=0,n_y-1
         if (ny.lt.n_y_l) then
            nyt = ny + n_y_u
         else
            nyt = ny - n_y_u
         end if
         do nx=0,n_x-1
            if (nx.lt.n_x_l) then
               nxt = nx + n_x_u
            else
               nxt = nx - n_x_u
            end if
            c_(nxt+nyt*n_x) = a_(nx+ny*n_x)
         enddo
      enddo
      do ny=0,n_y-1
         do nx=0,n_x-1
            b_(nx+ny*n_x) = c_(nx+ny*n_x)
         enddo
      enddo
      deallocate(c_)
      end      
