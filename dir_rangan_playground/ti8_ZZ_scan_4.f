!> Doxygen comment: ;\n
!> test_innerproduct_8 excerpt: ;\n
!> Scans across valid displacements for best innerproduct. ;\n
!> Note that if flag_RTRT_vs_RTTR.eqv..true.: ;\n
!> the eventual displacement will be \delta <-- \delta_upd + R_{-\gamma_upd}*\delta_est. ;\n
!> However if flag_RTRT_vs_RTTR.eqv..false.: ;\n
!> the eventual displacement will be \delta <-- R_{-\gamma_upd}*(\delta_upd + \delta_est). ;\n
!> Therefore, when flag_RTRT_vs_RTTR.eqv..true. we define \delta_max using \delta_upd + R_{-\gamma_upd}*\delta_est ;\n
!> and when flag_RTRT_vs_RTTR.eqv..false. we define \delta_max using (\delta_upd + \delta_est). ;\n
      subroutine ti8_ZZ_scan_4(flag_RTRT_vs_RTTR
     $     ,delta_x_est,delta_y_est,gamma_z_est,displacement_max
     $     ,n_delta_v,delta_x_,delta_y_,n_gamma_z,gamma_z_,ZZ_
     $     ,ndv_optimal,ngz_optimal ,C_Z_optimal)
      implicit none
      integer verbose
      data verbose / 0 /
      logical flag_RTRT_vs_RTTR ! flag indicating whether transformation operations are ordered as: R_{est}T_{est}R_{upd}T_{upd}(S) [flag_RTRT_vs_RTTR.eqv..true.] or R_{est}T_{est}T_{upd}R_{upd}(S) [flag_RTRT_vs_RTTR.eqv..false.] ;
      real *8 delta_x_est,delta_y_est,gamma_z_est,displacement_max
      real *8 delta_x_(0:n_delta_v-1),delta_y_(0:n_delta_v-1)
      real *8 gamma_z_(0:n_gamma_z-1)
      integer n_delta_v,n_gamma_z,ndv_optimal
     $     ,ngz_optimal
      complex *16 ZZ_(0:n_delta_v*n_gamma_z-1)
      complex *16 C_Z_optimal
      real *8 delta_x_pre,delta_y_pre
      real *8 delta_x_pos,delta_y_pos
      integer ndv,ngz,na,na_tmp
      complex *16 C_Z
      logical flag_found
      real *8 delta_pos
      if (verbose.gt.0) then
         write(6,'(A)') ' [entering ti8_ZZ_scan_4] '
      end if !if (verbose.gt.0) then
      if (verbose.gt.0) then
         write(6,'(A,I0)') ' n_delta_v: ' , n_delta_v
         write(6,'(A,I0)') ' n_gamma_z: ' , n_gamma_z
         write(6,'(A,F8.4)') ' displacement_max: ' , displacement_max
      end if !if (verbose.gt.0) then
      flag_found = .false.
      C_Z_optimal = dcmplx( -1.0d0 , 0.0d0 )
      ndv_optimal = 0
      ngz_optimal = 0
      do na=0,n_gamma_z*n_delta_v-1
         if (verbose.gt.1) then
            write(6,'(A,I0,A,F8.4,F8.4)') ' ZZ_(' , na , ') ' , ZZ_(na)
         end if !if (verbose.gt.1) then
         if (real(ZZ_(na)).gt.real(C_Z_optimal)) then
            na_tmp = na
            ndv = mod(na_tmp,n_delta_v)
            na_tmp = na_tmp/n_delta_v
            ngz = na_tmp
            delta_x_pre = delta_x_est
            delta_y_pre = delta_y_est
            call get_interchange_delta(delta_x_pre,delta_y_pre
     $           ,gamma_z_(ngz),delta_x_pos,delta_y_pos)
            if (flag_RTRT_vs_RTTR.eqv..true.) then
               delta_pos = dsqrt((delta_x_(ndv) + delta_x_pos)**2 +
     $              (delta_y_(ndv) + delta_y_pos)**2)
            end if              !if (flag_RTRT_vs_RTTR.eqv..true.) then
            if (flag_RTRT_vs_RTTR.eqv..false.) then
               delta_pos = dsqrt((delta_x_(ndv) + delta_x_est)**2 +
     $              (delta_y_(ndv) + delta_y_est)**2)
            end if              !if (flag_RTRT_vs_RTTR.eqv..false.) then
            if (verbose.gt.1) then
               write(6,'(A,I0,A,I0)') ' ngz ' , ngz , ' ndv ' , ndv
               write(6,'(A,F8.4)') ' delta_x_pre ' , delta_x_pre 
               write(6,'(A,F8.4)') ' delta_y_pre ' , delta_y_pre 
               write(6,'(A,F8.4)') ' delta_x_pos ' , delta_x_pos 
               write(6,'(A,F8.4)') ' delta_y_pos ' , delta_y_pos 
               write(6,'(A,F8.4)') ' delta_x_est ' , delta_x_est 
               write(6,'(A,F8.4)') ' delta_y_est ' , delta_y_est 
               write(6,'(A,F8.4)') ' delta_pos ' , delta_pos 
            end if              !if (verbose.gt.1) then
            if (delta_pos.le.displacement_max) then 
               C_Z = ZZ_(na)
               C_Z_optimal = C_Z
               ndv_optimal = ndv
               ngz_optimal = ngz
               flag_found = .true.
            end if !if (delta_pos.le.displacement_max) then 
         end if !if (real(ZZ_(na)).gt.real(C_Z_optimal)) then
      enddo !do na=0,n_gamma_z*n_delta_v-1
      if (verbose.gt.0) then
         write(6,'(A,I0)') ' ndv_optimal: ' , ndv_optimal
         write(6,'(A,I0)') ' ngz_optimal: ' , ngz_optimal
         write(6,'(A,2F8.4)') ' C_Z_optimal: ' , C_Z_optimal
      end if !if (verbose.gt.0) then
      if (flag_found.eqv..false.) then
         write(6,'(A)') 'Warning, flag_found.eqv..false. '
     $        ,'in ti8_ZZ_scan_4.f'
         write(6,'(A,F8.4)') ' displacement_max ' ,
     $        displacement_max 
         write(6,'(A,L1)') ' flag_RTRT_vs_RTTR ' ,
     $        flag_RTRT_vs_RTTR 
         write(6,'(A,I0)') ' n_gamma_z ' , n_gamma_z 
         write(6,'(A,I0)') ' n_delta_v ' , n_delta_v 
         do na=0,n_gamma_z*n_delta_v-1
            if (verbose.gt.-1) then
               write(6,'(A,I0,A,F8.4,F8.4)') ' ZZ_(' , na , ') ' ,
     $              ZZ_(na)
            end if              !if (verbose.gt.-1) then
            if (real(ZZ_(na)).gt.real(C_Z_optimal)) then
               na_tmp = na
               ndv = mod(na_tmp,n_delta_v)
               na_tmp = na_tmp/n_delta_v
               ngz = na_tmp
               delta_x_pre = delta_x_est
               delta_y_pre = delta_y_est
               call get_interchange_delta(delta_x_pre,delta_y_pre
     $              ,gamma_z_(ngz),delta_x_pos,delta_y_pos)
               if (flag_RTRT_vs_RTTR.eqv..true.) then
                  delta_pos = dsqrt((delta_x_(ndv) + delta_x_pos)**2 +
     $                 (delta_y_(ndv) + delta_y_pos)**2)
               end if           !if (flag_RTRT_vs_RTTR.eqv..true.) then
               if (flag_RTRT_vs_RTTR.eqv..false.) then
                  delta_pos = dsqrt((delta_x_(ndv) + delta_x_est)**2 +
     $                 (delta_y_(ndv) + delta_y_est)**2)
               end if           !if (flag_RTRT_vs_RTTR.eqv..false.) then
               if (verbose.gt.-1) then
                  write(6,'(2(A,I0),10(A,F8.4))')
     $                 ' ngz ' , ngz ,
     $                 ' ndv ' , ndv, 
     $                 ' delta_x_(ndv) ' , delta_x_(ndv) ,
     $                 ' delta_y_(ndv) ' , delta_y_(ndv) ,
     $                 ' delta_x_pre ' , delta_x_pre ,
     $                 ' delta_y_pre ' , delta_y_pre ,
     $                 ' delta_x_est ' , delta_x_est ,
     $                 ' delta_y_est ' , delta_y_est ,
     $                 ' delta_x_pos ' , delta_x_pos ,
     $                 ' delta_y_pos ' , delta_y_pos ,
     $                 ' delta_pos ' , delta_pos ,
     $                 ' displacement_max ' , displacement_max 
               end if           !if (verbose.gt.-1) then
               if (delta_pos.le.displacement_max) then 
                  C_Z = ZZ_(na)
                  C_Z_optimal = C_Z
                  ndv_optimal = ndv
                  ngz_optimal = ngz
               end if           !if (delta_pos.le.displacement_max) then 
            end if              !if (real(ZZ_(na)).gt.real(C_Z_optimal)) then
         enddo                  !do na=0,n_gamma_z*n_delta_v-1
         stop !exit program due to error ;
      end if !if (flag_found.eqv..false.) then
      if (verbose.gt.0) then
         write(6,'(A)') ' [finished ti8_ZZ_scan_4] '
      end if !if (verbose.gt.0) then
      end
