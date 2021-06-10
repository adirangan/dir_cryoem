      subroutine test_innerproduct_fast_ZZ_scan_1(flag_RTRT_vs_RTTR
     $     ,delta_x_est,delta_y_est,gamma_z_est,displacement_max
     $     ,n_delta_x,delta_x_,n_delta_y,delta_y_,n_gamma_z,gamma_z_,ZZ_
     $     ,ndx_optimal,ndy_optimal ,ngz_optimal ,C_Z_optimal)
      implicit none
      integer verbose
      data verbose / 0 /
      logical flag_RTRT_vs_RTTR ! flag indicating whether transformation operations are ordered as: R_{est}T_{est}R_{upd}T_{upd}(S) [flag_RTRT_vs_RTTR.eqv..true.] or R_{est}T_{est}T_{upd}R_{upd}(S) [flag_RTRT_vs_RTTR.eqv..false.] ;
      real *8 delta_x_est,delta_y_est,gamma_z_est,displacement_max
      real *8 delta_x_(0:n_delta_x-1),delta_y_(0:n_delta_y-1)
      real *8 gamma_z_(0:n_gamma_z-1)
      integer n_delta_x,n_delta_y,n_gamma_z,ndx_optimal,ndy_optimal
     $     ,ngz_optimal
      complex *16 ZZ_(0:n_delta_x*n_delta_y*n_gamma_z-1)
      complex *16 C_Z_optimal
      real *8 delta_x_pre,delta_y_pre
      real *8 delta_x_pos,delta_y_pos
      integer ndx,ndy,ngz,na
      complex *16 C_Z
      logical flag_found
      real *8 delta_pos
      if (verbose.gt.0) then
         write(6,'(A)') ' [entering test_innerproduct_fast_ZZ_scan_1] '
      end if !if (verbose.gt.0) then
      if (verbose.gt.0) then
         write(6,'(A,I0)') ' n_delta_x: ' , n_delta_x
         write(6,'(A,I0)') ' n_delta_y: ' , n_delta_y
         write(6,'(A,I0)') ' n_gamma_z: ' , n_gamma_z
         write(6,'(A,F8.4)') ' displacement_max: ' , displacement_max
      end if !if (verbose.gt.0) then
      flag_found = .false.
      ndx_optimal = 0
      ndy_optimal = 0
      ngz_optimal = 0
      na=0
      do ngz=0,n_gamma_z-1
         delta_x_pre = delta_x_est
         delta_y_pre = delta_y_est
         call get_interchange_delta(delta_x_pre,delta_y_pre
     $        ,gamma_z_(ngz),delta_x_pos,delta_y_pos)
         do ndy=0,n_delta_y-1
            do ndx=0,n_delta_x-1
               if (flag_RTRT_vs_RTTR.eqv..true.) then
                  delta_pos = dsqrt((delta_x_(ndx) + delta_x_pos)**2 +
     $                 (delta_y_(ndy) + delta_y_pos)**2)
               end if !if (flag_RTRT_vs_RTTR.eqv..true.) then
               if (flag_RTRT_vs_RTTR.eqv..false.) then
                  delta_pos = dsqrt((delta_x_(ndx) + delta_x_est)**2 +
     $                 (delta_y_(ndy) + delta_y_est)**2)
               end if !if (flag_RTRT_vs_RTTR.eqv..false.) then
               if (verbose.gt.1) then
                  write(6,'(A,F8.4)') ' delta_x_pre ' , delta_x_pre 
                  write(6,'(A,F8.4)') ' delta_y_pre ' , delta_y_pre 
                  write(6,'(A,F8.4)') ' delta_x_pos ' , delta_x_pos 
                  write(6,'(A,F8.4)') ' delta_y_pos ' , delta_y_pos 
                  write(6,'(A,F8.4)') ' delta_x_est ' , delta_x_est 
                  write(6,'(A,F8.4)') ' delta_y_est ' , delta_y_est 
                  write(6,'(A,F8.4)') ' delta_pos ' , delta_pos 
               end if !if (verbose.gt.1) then
               if (delta_pos.le.displacement_max) then
                  if (flag_found.eqv..false.) then
                     C_Z = ZZ_(na)
                     C_Z_optimal = C_Z
                     ndx_optimal = ndx
                     ndy_optimal = ndy
                     ngz_optimal = ngz
                     flag_found = .true.
                  else !if (flag_found.eqv..true.) then
                     C_Z = ZZ_(na)
                     if (real(C_Z).gt.real(C_Z_optimal)) then
                        C_Z_optimal = C_Z
                        ndx_optimal = ndx
                        ndy_optimal = ndy
                        ngz_optimal = ngz
                     end if !if new value is bigger
                  end if !flag_found
               end if ! (delta_pos.le.displacement_max)
               na = na+1
            enddo ! n_delta_x
         enddo ! n_delta_y
      enddo ! n_gamma_z
      if (verbose.gt.0) then
         write(6,'(A,I0)') ' ndx_optimal: ' , ndx_optimal
         write(6,'(A,I0)') ' ndy_optimal: ' , ndy_optimal
         write(6,'(A,I0)') ' ngz_optimal: ' , ngz_optimal
         write(6,'(A,2F8.4)') ' C_Z_optimal: ' , C_Z_optimal
      end if !if (verbose.gt.0) then
      if (flag_found.eqv..false.) then
         write(6,'(A)') 'Warning, flag_found.eqv..false. '
     $        ,'in test_innerproduct_fast_ZZ_scan_1.f'
      end if !if (flag_found.eqv..false.) then
      if (verbose.gt.0) then
         write(6,'(A)') ' [finished test_innerproduct_fast_ZZ_scan_1] '
      end if !if (verbose.gt.0) then
      end
