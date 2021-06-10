      subroutine test_innerproduct_5_ZZ_scan_2x(flag_RTRT_vs_RTTR
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
      integer ndx,ndy,ngz,na,na_tmp
      complex *16 C_Z
      logical flag_found
      real *8 delta_pos
      if (verbose.gt.0) then
         write(6,'(A)') ' [entering test_innerproduct_5_ZZ_scan_2x] '
      end if !if (verbose.gt.0) then
      if (verbose.gt.0) then
         write(6,'(A,I0)') ' n_delta_x: ' , n_delta_x
         write(6,'(A,I0)') ' n_delta_y: ' , n_delta_y
         write(6,'(A,I0)') ' n_gamma_z: ' , n_gamma_z
         write(6,'(A,F8.4)') ' displacement_max: ' , displacement_max
      end if !if (verbose.gt.0) then

      flag_found = .true.
      C_Z_optimal = cmplx( -1.0d0 , 0.0d0 )
      ndx_optimal = 0
      ndy_optimal = 0
      ngz_optimal = 0
      do na=0,n_gamma_z*n_delta_x*n_delta_y-1
         if (real(ZZ_(na)).gt.real(C_Z_optimal)) then
            C_Z_optimal = ZZ_(na)
            na_tmp = na
            ndx_optimal = mod(na_tmp,n_delta_x)
            na_tmp = na_tmp/n_delta_x
            ndy_optimal = mod(na_tmp,n_delta_y)
            na_tmp = na_tmp/n_delta_y
            ngz_optimal = na_tmp
         end if !if (real(ZZ_(na)).gt.real(C_Z_optimal)) then
      enddo !do na=0,n_gamma_z*n_delta_x*n_delta_y-1

c$$$      flag_found = .false.
c$$$      C_Z_optimal = cmplx( -1.0d0 , 0.0d0 )
c$$$      ndx_optimal = 0
c$$$      ndy_optimal = 0
c$$$      ngz_optimal = 0
c$$$      na=0
c$$$      do ngz=0,n_gamma_z-1
c$$$         delta_x_pre = delta_x_est
c$$$         delta_y_pre = delta_y_est
c$$$         call get_interchange_delta(delta_x_pre,delta_y_pre
c$$$     $        ,gamma_z_(ngz),delta_x_pos,delta_y_pos)
c$$$         do ndy=0,n_delta_y-1
c$$$            do ndx=0,n_delta_x-1
c$$$c$$$               if (flag_RTRT_vs_RTTR.eqv..true.) then
c$$$c$$$                  delta_pos = dsqrt((delta_x_(ndx) + delta_x_pos)**2 +
c$$$c$$$     $                 (delta_y_(ndy) + delta_y_pos)**2)
c$$$c$$$               end if !if (flag_RTRT_vs_RTTR.eqv..true.) then
c$$$c$$$               if (flag_RTRT_vs_RTTR.eqv..false.) then
c$$$c$$$                  delta_pos = dsqrt((delta_x_(ndx) + delta_x_est)**2 +
c$$$c$$$     $                 (delta_y_(ndy) + delta_y_est)**2)
c$$$c$$$               end if !if (flag_RTRT_vs_RTTR.eqv..false.) then
c$$$c$$$               if (verbose.gt.1) then
c$$$c$$$                  write(6,'(A,F8.4)') ' delta_x_pre ' , delta_x_pre 
c$$$c$$$                  write(6,'(A,F8.4)') ' delta_y_pre ' , delta_y_pre 
c$$$c$$$                  write(6,'(A,F8.4)') ' delta_x_pos ' , delta_x_pos 
c$$$c$$$                  write(6,'(A,F8.4)') ' delta_y_pos ' , delta_y_pos 
c$$$c$$$                  write(6,'(A,F8.4)') ' delta_x_est ' , delta_x_est 
c$$$c$$$                  write(6,'(A,F8.4)') ' delta_y_est ' , delta_y_est 
c$$$c$$$                  write(6,'(A,F8.4)') ' delta_pos ' , delta_pos 
c$$$c$$$               end if !if (verbose.gt.1) then
c$$$c$$$               if (delta_pos.le.displacement_max) then
c$$$                  if (flag_found.eqv..false.) then
c$$$                     C_Z = ZZ_(na)
c$$$                     C_Z_optimal = C_Z
c$$$                     ndx_optimal = ndx
c$$$                     ndy_optimal = ndy
c$$$                     ngz_optimal = ngz
c$$$                     flag_found = .true.
c$$$                  else !if (flag_found.eqv..true.) then
c$$$                     C_Z = ZZ_(na)
c$$$                     if (real(C_Z).gt.real(C_Z_optimal)) then
c$$$                        C_Z_optimal = C_Z
c$$$                        ndx_optimal = ndx
c$$$                        ndy_optimal = ndy
c$$$                        ngz_optimal = ngz
c$$$                     end if !if new value is bigger
c$$$                  end if !flag_found
c$$$c$$$               end if ! (delta_pos.le.displacement_max)
c$$$               na = na+1
c$$$            enddo ! n_delta_x
c$$$         enddo ! n_delta_y
c$$$      enddo ! n_gamma_z

      if (verbose.gt.0) then
         write(6,'(A,I0)') ' ndx_optimal: ' , ndx_optimal
         write(6,'(A,I0)') ' ndy_optimal: ' , ndy_optimal
         write(6,'(A,I0)') ' ngz_optimal: ' , ngz_optimal
         write(6,'(A,2F8.4)') ' C_Z_optimal: ' , C_Z_optimal
      end if !if (verbose.gt.0) then
      if (flag_found.eqv..false.) then
         write(6,'(A)') 'Warning, flag_found.eqv..false. '
     $        ,'in test_innerproduct_5_ZZ_scan_2x.f'
      end if !if (flag_found.eqv..false.) then
      if (verbose.gt.0) then
         write(6,'(A)') ' [finished test_innerproduct_5_ZZ_scan_2x] '
      end if !if (verbose.gt.0) then
      end
