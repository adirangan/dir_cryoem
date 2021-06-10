!> Doxygen comment: ;\n
!> gets star-shaped ctf in k_p coordinates: ;\n
      subroutine get_ctf_star_k_p_(n_r,grid_k_p_,n_w_,n_A,CTF_k_p_
     $     ,param_0,param_1,param_2)
      implicit none
      integer verbose
      data verbose / 0 /
      integer n_r,n_w_(0:n_r-1),n_A
      real *8 grid_k_p_(0:n_r-1)
      complex *16 CTF_k_p_(0:0)
      real *8 param_0,param_1,param_2
      real *8 pi,Rval,Wval,RR,WW,WS,FF
      integer nr,nw,na

      if (verbose.gt.0) then
         write(6,'(A,I0)') '[entering get_ctf_star_k_p_]'
      end if
      
      pi = 4*datan(1.0d0)
      na = 0
      do nr=0,n_r-1
         Rval = grid_k_p_(nr)
         if (param_0.gt.0.0d0) then
            RR = dexp(-Rval/param_0)
         else
            RR = 1.0d0
         end if
         do nw=0,n_w_(nr)-1
            Wval = (2.0d0*pi*nw)/n_w_(nr)
            WW = dcos(param_2*(Wval-param_1))
            if (WW.ge.0.0d0) then
               WS = 1.0d0
            else
               WS = -1.0d0
            end if
            FF = RR*WS
            CTF_k_p_(na) = dcmplx( FF , 0.0d0) 
            na = na+1
         enddo
      enddo
      
      if (verbose.gt.0) then
         write(6,'(A,I0)') '[finished get_ctf_star_k_p_]'
      end if

      end
