!> Doxygen comment: ;\n
!> gets trivial identity ctf in k_p coordinates : ;\n
      subroutine get_ctf_ones_k_p_(n_r,grid_k_p_,n_w_,n_A,CTF_k_p_)
      implicit none
      integer verbose
      data verbose / 0 /
      integer n_r,n_w_(0:n_r-1),n_A
      real *8 grid_k_p_(0:n_r-1)
      complex *16 CTF_k_p_(0:0)
      integer nr,nw,na

      if (verbose.gt.0) then
         write(6,'(A,I0)') '[entering get_ctf_ones_k_p_]'
      end if
      
      na = 0
      do nr=0,n_r-1
         do nw=0,n_w_(nr)-1
            CTF_k_p_(na) = dcmplx( 1.0d0 , 0.0d0) 
            na = na+1
         enddo
      enddo
      
      if (verbose.gt.0) then
         write(6,'(A,I0)') '[finished get_ctf_ones_k_p_]'
      end if

      end
