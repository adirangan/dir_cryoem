
      nm0 = nM_9_sum + nm
      n_LT_add = 2

      if (verbose.gt.1) then
         write(6,'(A)') ' Set flag_S_use_ to false and clear LT_.'
      end if                    !if (verbose.gt.1) then
      call cl1_r8(3,vp_input_omp__(p_vp_input))
      call cl1_l2(nS_0_per,flag_S_use_omp__(p_flag_S_use))
      call cl1_i4(2*nS_0_per,LT_omp__(p_LT))
      n_S_use = 0

      if (0+verbose.gt.1) then
         write(6,'(A,I0,A)') ' Generating initial list for image ',
     $        nm0,' based on estimated '
         write(6,'(A,F8.4,A,F8.4,A)') ' polar_a ',polar_a_est_(nm0)
     $        ,' and azimu_b ',azimu_b_est_(nm0), '.'
         write(6,'(A)') ' '
         write(6,'(A)') ' This initial list comprises only the '
         write(6,'(A)') ' templates that have a polar_a and an '
         write(6,'(A)') ' azimu_b that are within distance  '
         write(6,'(A)') ' tesselation_distance_req of the image '
         write(6,'(A)') ' angles (as considered as vectors on '
         write(6,'(A)') ' the sphere). '
         write(6,'(A)') ' '
         write(6,'(A)') ' Note that if there are no templates '
         write(6,'(A)') ' within the requested distance, we  '
         write(6,'(A)') ' broaden our search until there are '
         write(6,'(A)') ' at least two templates in range. '
         write(6,'(A)') ' '
         write(6,'(A)') ' Note also that if you wish to perform a '
         write(6,'(A)') ' global search, simply set the required '
         write(6,'(A)') ' tesselation_distance_req to encompass '
         write(6,'(A)') ' the entire sphere (i.e., .ge. 2.0d0). '
         write(6,'(A)') ' '
         write(6,'(A)') ' Note also that we track which templates ' 
         write(6,'(A)') ' have been used for this particular image. '
         write(6,'(A)') ' (see temporary logical array flag_S_use_, '
         write(6,'(A)') ' as well as temporary indexing array '
         write(6,'(A)') ' LT_). '
         write(6,'(A)') ' '
      end if                    !if (0+verbose.gt.1) then
      call get_angle_to_vp_(polar_a_est_(nm0),azimu_b_est_(nm0)
     $     ,vp_input_omp__(p_vp_input))
      if (verbose.gt.1) then
         write(6,'(A,I0,2(A,F8.4),A,3(F8.4,1X))') ' nm0 ' , nm0 ,
     $        ' polar_a: ' ,polar_a_est_(nm0) , ' azimu_b: ' ,
     $        azimu_b_est_(nm0) ,' vp_input: ' , vp_input_omp__(0
     $        +p_vp_input) ,vp_input_omp__(1+p_vp_input) ,
     $        vp_input_omp__(2+p_vp_input)
      end if !if (verbose.gt.1) then
      n_LT = 0
      tesselation_distance_req_use = tesselation_distance_req
      do while (nS_0_per.gt.0 .and. n_LT.lt.min(nS_0_per,2))
         call tesselation_neighborhood_wrapper_0(nS_0_per,S_L_ ,nl_max
     $        ,nm_sum,ll_sum,T_nl_,T_vm_,T_tr_ ,T_ll_,T_lf_,T_c0_ ,T_c1_
     $        ,T_c2_,T_c3_,T_ls_,T_LT_,n_T_root_base ,T_root_base_
     $        ,vp_input_omp__(p_vp_input),tesselation_distance_req_use
     $        ,n_LT ,LT_omp__(0+p_LT))
         if (n_LT.lt.min(nS_0_per,2)) then
            if (verbose.gt.1) then
               include 'ti6_write_LT_0.f'
               write(6,'(2(A,F8.4))') ' increasing ' ,
     $              tesselation_distance_req_use , ' to ' ,
     $              tesselation_distance_req_use + 0.0125
            end if              !if (verbose.gt.1) then
            tesselation_distance_req_use =
     $           tesselation_distance_req_use + 0.0125
         end if                 !if (n_LT.lt.min(nS_0_per,2)) then
      enddo                     !do while (nS_0_per.gt.0 .and. n_LT.lt.min(nS_0_per,2))
      if (0+verbose.gt.1) then
         write(6,'(3(A,I0))') ' nm0: ' , nm0 , ' found ',n_LT
     $        ,' templates out of a possible nS_0_per: ' , nS_0_per
         include 'ti6_write_LT_0.f'
      end if                    !if (0+verbose.gt.1) then
      do nLT=0,n_LT-1
         flag_S_use_omp__(LT_omp__(nLT+p_LT)+p_flag_S_use)=.true.
         n_S_use = n_S_use+1
      enddo                     !do nLT=0,n_LT-1
      n_add = max(0,min(nS_0_per-n_S_use,n_LT_add))
      if (0+verbose.gt.1) then
         write(6,'(5(A,I0))') ' nS_0_per: ' , nS_0_per , ' n_LT: ' ,
     $        n_LT , ' n_S_use: ' , n_S_use, ' n_LT_add: ' , n_LT_add ,
     $        ' n_add: ' , n_add
      end if                    !if (0+verbose.gt.1) then
      if (0+verbose.gt.1) then
         write(6,'(A)') ' '
         write(6,'(A)') ' This step adds an additional number of'
         write(6,'(A)') ' randomly chosen (unique) templates from '
         write(6,'(A)') ' the list of templates. These additional '
         write(6,'(A)') ' templates will not necessarily  be close '
         write(6,'(A)') ' to the image in terms of polar_a and '
         write(6,'(A)') ' azimu_b. '
         write(6,'(A)') ' '
         write(6,'(A)') ' The number of additional templates '
         write(6,'(A)') ' is half the number originally considered. '
         write(6,'(A)') ' (i.e., half the number within '
         write(6,'(A)') ' tesselation_distance_req of the original '
         write(6,'(A)') ' estimated image-angle). '
         write(6,'(A)') ' If you wish to increase this number, '
         write(6,'(A)') ' modify the variable: n_LT_hlf above. '
         write(6,'(A)') ' '
         write(6,'(A)') ' Note that we always try to add at least '
         write(6,'(A)') ' a dozen templates if possible '
         write(6,'(A)') ' (see parameter n_LT_add above).'
         write(6,'(A)') ' '
         write(6,'(A)') ' Note that these extra templates are also '
         write(6,'(A)') ' accounted for in arrays flag_S_use_ and LT_. '
         write(6,'(A)') ' '
      end if                    !if (0+verbose.gt.1) then
      call adi_randinclude(rseed,n_LT,LT_omp__(0+p_LT),nS_0_per
     $     ,flag_S_use_omp__(0+p_flag_S_use),n_S_use,n_add)
      if (0+verbose.gt.1) then
         write(6,'(3(A,I0))') ' nS_0_per: ' , nS_0_per , ' n_S_use: '
     $        ,n_S_use,' n_LT: ' , n_LT
         include 'ti6_write_LT_0.f'
         call write_all_l2(nS_0_per,flag_S_use_omp__(p_flag_S_use),13
     $        ,' flag_S_use: ')
         call write_all_i4(n_LT,LT_omp__(p_LT),6,' LT_: ')
      end if                    !if (0+verbose.gt.1) then
      
      if (0+verbose.gt.1) then
         write(6,'(A)') ' '
         write(6,'(A)') ' At this point we step through the '
         write(6,'(A)') ' list of templates that has been '
         write(6,'(A)') ' generated for this particular '
         write(6,'(A)') ' image (i.e., including both nearby '
         write(6,'(A)') ' templates and a few randomly '
         write(6,'(A)') ' chosen templates from elsewhere on'
         write(6,'(A)') ' the sphere). '
         write(6,'(A)') ' '
      end if                    !if (0+verbose.gt.1) then

      include 'ti6_calculate_0.f'

c$$$      if (0+verbose.gt.1) then
c$$$         write (6,'(A)') ' After calculating the initial list of '
c$$$         write (6,'(A)') ' innerproducts for this particular image '
c$$$         write (6,'(A)') ' (see above), we pass through this list '
c$$$         write (6,'(A)') ' several times: '
c$$$         write (6,'(A)') ' '
c$$$         write (6,'(A)') ' Step 1: '
c$$$         write (6,'(A)') ' Each time we pass through the list we '
c$$$         write (6,'(A)') ' find the templates that correspond to '
c$$$         write (6,'(A)') ' the largest n_LT_hlf innerproducts. '
c$$$         write (6,'(A)') ' Note that this number can certainly be '
c$$$         write (6,'(A)') ' increased or decreased (just by '
c$$$         write (6,'(A)') ' modifying or replacing n_LT_hlf within '
c$$$         write (6,'(A)') ' the following code). '
c$$$         write (6,'(A)') ' '
c$$$         write (6,'(A)') ' Step 2: '
c$$$         write (6,'(A)') ' Once we find these "good" templates, '
c$$$         write (6,'(A)') ' we search for other templates nearby. '
c$$$         write (6,'(A)') ' As before, we use the distance '
c$$$         write (6,'(A)') ' tesselation_distance_req_use, which '
c$$$         write (6,'(A)') ' was previously used to generate our '
c$$$         write (6,'(A)') ' original list of templates. '
c$$$         write (6,'(A)') ' '
c$$$         write (6,'(A)') ' Now (hopefully) many of these nearby '
c$$$         write (6,'(A)') ' templates close to the "good" templates '
c$$$         write (6,'(A)') ' will have already been considered for '
c$$$         write (6,'(A)') ' this particular image. (i.e, many '
c$$$         write (6,'(A)') ' nearby templates will alread be listed '
c$$$         write (6,'(A)') ' in LT_ and tagged in flag_S_use_). '
c$$$         write (6,'(A)') ' However, some of these nearby templates '
c$$$         write (6,'(A)') ' will be new. '
c$$$         write (6,'(A)') ' '
c$$$         write (6,'(A)') ' Step 3: '
c$$$         write (6,'(A)') ' We add each of these new nearby '
c$$$         write (6,'(A)') ' templates to our list of templates '
c$$$         write (6,'(A)') ' (updating LT_ and flag_S_use_ as we go) '
c$$$         write (6,'(A)') ' and then compute the innerproducts for '
c$$$         write (6,'(A)') ' these new templates. '
c$$$         write (6,'(A)') ' '
c$$$         write (6,'(A)') ' After calculating these new '
c$$$         write (6,'(A)') ' innerproducts (for the new templates in '
c$$$         write (6,'(A)') ' our now-larger list) we go back to '
c$$$         write (6,'(A)') ' Step 1, and once again find the '
c$$$         write (6,'(A)') ' templates that correspond to the largest '
c$$$         write (6,'(A)') ' innerproducts. '
c$$$         write (6,'(A)') ' '
c$$$         write (6,'(A)') ' After passing through our slowly growing '
c$$$         write (6,'(A)') ' list multiple times, we will eventually '
c$$$         write (6,'(A)') ' reach a point where each of the nearby '
c$$$         write (6,'(A)') ' templates was already considered. '
c$$$         write (6,'(A)') ' (i.e., after Step-2 the list does not '
c$$$         write (6,'(A)') ' grow). At this point we consider our '
c$$$         write (6,'(A)') ' local-search complete.' 
c$$$         write (6,'(A)') ' '
c$$$         write (6,'(A)') ' Step 4: '
c$$$         write (6,'(A)') ' Now we run through our list of '
c$$$         write (6,'(A)') ' innerproducts. The largest few will be '
c$$$         write (6,'(A)') ' used to update the alpha1d_ parameters. '
c$$$         write (6,'(A)') ' '
c$$$      end if                    !if (0+verbose.gt.1) then
c$$$      
c$$$      n_pass = 0
c$$$      continue_flag = .true.
c$$$
c$$$      do while (continue_flag)
c$$$
c$$$         n_LT_srt = 0
c$$$         C_Z_optimal = (1.0d0,0.0d0)
c$$$         C_S_optimal = (0.0d0,0.0d0)
c$$$         do ns=0,n_S-1
c$$$            if (flag_S_use_(ns).eqv..true.) then
c$$$               C_Z_optimal = C_Z_SM_(ns)
c$$$               C_S_optimal = C_S_SM_(ns)
c$$$               C_Z_srt_(n_LT_srt) = C_Z_optimal
c$$$               if (real(C_Z_optimal).lt.real(C_Z_min)) then
c$$$                  C_Z_min = C_Z_optimal
c$$$               end if           !if (real(C_Z_optimal).lt.real(C_Z_min)) then
c$$$               if (real(C_S_optimal).lt.real(C_S_min)) then
c$$$                  C_S_min = C_S_optimal
c$$$               end if           !if (real(C_S_optimal).lt.real(C_S_min)) then
c$$$               LT_srt_(n_LT_srt) = ns
c$$$               n_LT_srt = n_LT_srt + 1
c$$$            end if              !if (flag_S_use_(ns).eqv..true.) then
c$$$         enddo
c$$$         
c$$$         if (0+verbose.gt.1) then
c$$$            write(6,'(A,I0)') 'n_pass: ' , n_pass
c$$$            write(6,'(A,I0)') 'found n_LT_srt: ' , n_LT_srt
c$$$            write(format_string,'(A,I0,A)') '(A,' , 2*n_LT_srt ,
c$$$     $           '(F8.4,1X))'
c$$$            write(6,format_string) 'pre: C_Z_srt_: ' , (C_Z_srt_(nLT)
c$$$     $           ,nLT=0,n_LT_srt-1)
c$$$            write(format_string,'(A,I0,A)') '(A,' , n_LT_srt ,
c$$$     $           '(I0,1X))'
c$$$            write(6,format_string) 'pre: LT_srt_: ' , (LT_srt_(nLT)
c$$$     $           ,nLT=0,n_LT_srt-1)
c$$$         end if                 !if (0+verbose.gt.1) then
c$$$         call quicksort_c16(0,n_LT_srt-1,C_Z_srt_,1,LT_srt_,1
c$$$     $        ,quicksort_c16)
c$$$         if (0+verbose.gt.1) then
c$$$            write(format_string,'(A,I0,A)') '(A,' , 2*n_LT_srt ,
c$$$     $           '(F8.4,1X))'
c$$$            write(6,format_string) 'pos: C_Z_srt_: ' , (C_Z_srt_(nLT)
c$$$     $           ,nLT=0,n_LT_srt-1)
c$$$            write(format_string,'(A,I0,A)') '(A,' , n_LT_srt ,
c$$$     $           '(I0,1X))'
c$$$            write(6,format_string) 'pos: LT_srt_: ' , (LT_srt_(nLT)
c$$$     $           ,nLT=0,n_LT_srt-1)
c$$$         end if                 !if (0+verbose.gt.1) then
c$$$
c$$$         n_LT = 0
c$$$         if (0+verbose.gt.1) then
c$$$            write(6,'(A,I0)') ' using last entries: ' , min(n_LT_srt
c$$$     $           ,n_LT_hlf)
c$$$         end if                 ! if (0+verbose.gt.1) then
c$$$         do nLT=0,min(n_LT_srt,n_LT_hlf)-1
c$$$            ns_use = LT_srt_(n_LT_srt-1-nLT)
c$$$            call get_angle_to_vp_(S_alpha_polar_a_(ns_use)
c$$$     $           ,S_alpha_azimu_b_(ns_use),vp_input_omp__(p_vp_input))
c$$$            if (0+verbose.gt.1) then
c$$$               write(6,'(A,I0,A,F8.4,F8.4)') ' accessing ns_use ' ,
c$$$     $              ns_use , ' with polar_a and azimu_b: ' ,
c$$$     $              S_alpha_polar_a_(ns_use) ,
c$$$     $              S_alpha_azimu_b_(ns_use)
c$$$            end if              !if (0+verbose.gt.1) then
c$$$            call tesselation_neighborhood_wrapper_0(n_S,S_L_ ,nl_max
c$$$     $           ,nm_sum,ll_sum,T_nl_,T_vm_,T_tr_ ,T_ll_,T_lf_,T_c0_
c$$$     $           ,T_c1_ ,T_c2_,T_c3_,T_ls_,T_LT_,n_T_root_base
c$$$     $           ,T_root_base_ ,vp_input_omp__(p_vp_input)
c$$$     $           ,tesselation_distance_req_use,n_LT_tmp ,LT_(n_LT))
c$$$            n_LT = n_LT + n_LT_tmp
c$$$            if (0+verbose.gt.1) then
c$$$               write(6,'(A,I0,A,I0)') ' pre: n_LT_tmp ' , n_LT_tmp ,
c$$$     $              ' n_LT ', n_LT
c$$$            end if              !if (0+verbose.gt.1) then
c$$$            call unique_i4(n_LT,LT_,n_LT,LT_)
c$$$            if (0+verbose.gt.1) then
c$$$               write(6,'(A,I0)') ' pos: n_LT ' , n_LT
c$$$            end if              !if (0+verbose.gt.1) then
c$$$         enddo                  !do nLT=0,max(n_LT_srt-1,n_LT_hlf-1)
c$$$
c$$$         n_LT_use=0
c$$$         do nLT=0,n_LT-1
c$$$            ns_use = LT_(nLT)
c$$$            if (flag_S_use_(ns_use).eqv..true.) then
c$$$               if (0+verbose.gt.1) then
c$$$                  write(6,'(A,I0,A,I0)') ' nLT ' , nLT ,
c$$$     $                 ', skipping ns_use ' , ns_use
c$$$               end if           !if (0+verbose.gt.1) then
c$$$            end if              !if (flag_S_use_(ns_use).eqv..true.) then
c$$$            if (flag_S_use_(ns_use).eqv..false.) then
c$$$
c$$$               if (0+verbose.gt.1) then
c$$$                  write(6,'(A,I0,A,I0,A,I0)') ' nLT ' , nLT ,
c$$$     $                 ' Calculating innerproduct for image ' , nm0 ,
c$$$     $                 ' and template ' , ns_use
c$$$               end if           !if (0+verbose.gt.1) then
c$$$               call test_innerproduct_batch_SM_stage_2a(verbose ,nm
c$$$     $              ,ns_use ,n_r,grid_p_,n_w_,n_S,ld_S,n_A,nctf
c$$$     $              ,n_w_max ,S_p_,S_p,S_q,M_p,M_q,Z_q_ ,C_M,C_S_
c$$$     $              ,C_S_optimal,C_Z_optimal ,svd_calculation_type
c$$$     $              ,n_svd_r ,svd_r_,n_svd_l ,svd_l_,svd_s_,svd_V_r_
c$$$     $              ,Z_svdd_ ,Z_svdr_,Z_tmpC_ ,displacement_max
c$$$     $              ,delta_x_est ,delta_y_est,gamma_z_est ,n_delta_x
c$$$     $              ,n_delta_y,n_gamma_z ,delta_x_,delta_y_ ,gamma_z_
c$$$     $              ,delta_x_est_ ,delta_y_est_ ,gamma_z_est_
c$$$     $              ,ndx_optimal ,ndy_optimal ,ngz_optimal
c$$$     $              ,delta_x_optimal ,delta_y_optimal
c$$$     $              ,gamma_z_optimal ,fftw_plan_frwd_
c$$$     $              ,fftw_plan_back_ ,fftw_in1_ ,fftw_out_
c$$$     $              ,timing_tic_1,timing_toc_1 ,timing_tot_1
c$$$     $              ,timing_tic_2,timing_toc_2 ,timing_tot_2
c$$$     $              ,timing_tic_3,timing_toc_3 ,timing_tot_3
c$$$     $              ,timing_tic_4,timing_toc_4 ,timing_tot_4)
c$$$               delta_x_SM_(ns_use) = delta_x_optimal
c$$$               delta_y_SM_(ns_use) = delta_y_optimal
c$$$               gamma_z_SM_(ns_use) = gamma_z_optimal
c$$$               C_S_SM_(ns_use) = C_S_optimal
c$$$               C_Z_SM_(ns_use) = C_Z_optimal
c$$$               flag_S_use_(ns_use) = .true.
c$$$               n_S_use = n_S_use+1
c$$$               n_LT_use = n_LT_use + 1
c$$$            end if              ! if (flag_S_use_(ns_use).eqv..true.) then
c$$$         enddo                  !do nLT=0,n_LT-1
c$$$
c$$$         if (0+verbose.gt.1) then
c$$$            write(6,'(A,I0)') ' n_LT_use ' , n_LT_use
c$$$            write(6,'(A,I0)') ' n_S_use ' , n_S_use
c$$$         end if                 !if (0+verbose.gt.1) then
c$$$         if (n_LT_use.gt.0) then
c$$$            if (0+verbose.gt.1) then
c$$$               write(6,'(A)') ' Our list grew: do another pass. '
c$$$            end if              !if (0+verbose.gt.1) then
c$$$            continue_flag = .true.
c$$$         else
c$$$            if (0+verbose.gt.1) then
c$$$               write(6,'(A)') ' Our list did not grow: stop. '
c$$$            end if              !if (0+verbose.gt.1) then
c$$$            continue_flag = .false.
c$$$         end if                 !if (n_LT_use.gt.0) then
c$$$
c$$$         n_pass = n_pass + 1
c$$$      enddo                     !continue_flag
