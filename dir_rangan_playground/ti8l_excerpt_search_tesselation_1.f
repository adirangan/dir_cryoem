!> Doxygen comment: ;\n
!> test_innerproduct_8 excerpt: ;\n
!> Searches tesselation-tree for local search. ;\n
      nm0 = nM_9_sum + nm

      if (verbose.gt.1) then
         write(6,'(A)') ' Set flag_S_use_ to false and clear LT_.'
      end if                    !if (verbose.gt.1) then
      call cl1_r8(3,vp_input_omp__(p_vp_input))
      call cl1_l2(nS_0_per,flag_S_use_omp__(p_flag_S_use))
      call cl1_i4(2*nS_0_per,LT_omp__(p_LT))
      n_S_use = 0


c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
c$$$      Initial list ;
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

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
         write(6,'(A)') ' (see temporary logical array '
         write(6,'(A)') ' flag_S_use_omp__, '
         write(6,'(A)') ' as well as temporary indexing array '
         write(6,'(A)') ' LT_omp__). '
         write(6,'(A)') ' '
         write(6,'(A)') ' Note that LT_omp__ stores indices j '
         write(6,'(A)') ' ranging from 0 to nS_0_per, '
         write(6,'(A)') ' corresponding to the order of the '
         write(6,'(A)') ' templates in the block associated with '
         write(6,'(A)') ' ns_0_sub. '
         write(6,'(A)') ' These indices are used to access '
         write(6,'(A)') ' elements of flag_S_use_omp__. '
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
               include 'ti8l_excerpt_write_LT_0.f'
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
         include 'ti8l_excerpt_write_LT_0.f'
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
         write(6,'(A)') ' is n_LT_add. '
         write(6,'(A)') ' If you wish to increase this number, '
         write(6,'(A)') ' modify the variable: n_LT_add above. '
         write(6,'(A)') ' '
         write(6,'(A)') ' Note that these extra templates are also '
         write(6,'(A)') ' accounted for in arrays '
         write(6,'(A)') ' flag_S_use_omp__ and LT_omp__. '
         write(6,'(A)') ' '
      end if                    !if (0+verbose.gt.1) then
      call adi_randinclude(rseed,n_LT,LT_omp__(0+p_LT),nS_0_per
     $     ,flag_S_use_omp__(0+p_flag_S_use),n_S_use,n_add)
      if (0+verbose.gt.1) then
         write(6,'(3(A,I0))') ' nS_0_per: ' , nS_0_per , ' n_S_use: '
     $        ,n_S_use,' n_LT: ' , n_LT
         include 'ti8l_excerpt_write_LT_0.f'
         call print_all_l2(nS_0_per,flag_S_use_omp__(p_flag_S_use)
     $        ,' flag_S_use: ')
         call print_all_i4(n_LT,LT_omp__(p_LT),' LT_: ')
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

      flag_tesselation_ref = .false.
      include 'ti8l_excerpt_calculate_1.f'

c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
c$$$      Updated list ;
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

      if (0+verbose.gt.1) then
         write (6,'(A)') ' After calculating the initial list of '
         write (6,'(A)') ' innerproducts for this particular image '
         write (6,'(A)') ' (see above), we pass through this list '
         write (6,'(A)') ' several times: '
         write (6,'(A)') ' '
         write (6,'(A)') ' Step 1: '
         write (6,'(A)') ' Each time we pass through the list we '
         write (6,'(A)') ' find the templates that correspond to '
         write (6,'(A)') ' the largest n_LT_ref innerproducts. '
         write (6,'(A)') ' Note that this number can certainly be '
         write (6,'(A)') ' increased or decreased (just by '
         write (6,'(A)') ' modifying or replacing n_LT_ref within '
         write (6,'(A)') ' the following code). '
         write (6,'(A)') ' However, n_LT_ref should be less than '
         write (6,'(A)') ' the total number of templates stored '
         write (6,'(A)') ' for each image (i.e., less than '
         write (6,'(A)') ' n_SM_max). '
         write (6,'(A)') ' '
         write (6,'(A)') ' Step 2: '
         write (6,'(A)') ' Once we find these "good" templates, '
         write (6,'(A)') ' we search for other templates nearby. '
         write (6,'(A)') ' As before, we use the distance '
         write (6,'(A)') ' tesselation_distance_req_use, which '
         write (6,'(A)') ' was previously used to generate our '
         write (6,'(A)') ' original list of templates. '
         write (6,'(A)') ' '
         write (6,'(A)') ' Now (hopefully) many of these nearby '
         write (6,'(A)') ' templates close to the "good" templates '
         write (6,'(A)') ' will have already been considered for '
         write (6,'(A)') ' this particular image. (i.e, many '
         write (6,'(A)') ' nearby templates will alread be listed '
         write (6,'(A)') ' in LT_omp__ and tagged in ' 
         write (6,'(A)') ' flag_S_use_omp__). '
         write (6,'(A)') ' However, some of these nearby templates '
         write (6,'(A)') ' will be new. '
         write (6,'(A)') ' '
         write (6,'(A)') ' Step 3: '
         write (6,'(A)') ' We add each of these new nearby '
         write (6,'(A)') ' templates to our list of templates '
         write (6,'(A)') ' (updating LT_omp__ and '
         write (6,'(A)') ' flag_S_use_omp__ as we go) '
         write (6,'(A)') ' and then compute the innerproducts for '
         write (6,'(A)') ' these new templates. '
         write (6,'(A)') ' '
         write (6,'(A)') ' After calculating these new '
         write (6,'(A)') ' innerproducts (for the new templates in '
         write (6,'(A)') ' our now-larger list) we go back to '
         write (6,'(A)') ' Step 1, and once again find the '
         write (6,'(A)') ' templates that correspond to the largest '
         write (6,'(A)') ' innerproducts. '
         write (6,'(A)') ' '
         write (6,'(A)') ' After passing through our slowly growing '
         write (6,'(A)') ' list multiple times, we will eventually '
         write (6,'(A)') ' reach a point where each of the nearby '
         write (6,'(A)') ' templates was already considered. '
         write (6,'(A)') ' (i.e., after Step-2 the list does not '
         write (6,'(A)') ' grow). At this point we consider our '
         write (6,'(A)') ' local-search complete.' 
         write (6,'(A)') ' '
         write (6,'(A)') ' Step 4: '
         write (6,'(A)') ' Now we run through our list of '
         write (6,'(A)') ' innerproducts. The largest few will be '
         write (6,'(A)') ' used to update the alpha1d_ parameters. '
         write (6,'(A)') ' '
      end if                    !if (0+verbose.gt.1) then
      
      n_pass = 0
      continue_flag = .true.

      do while (continue_flag)

         if (verbose.gt.2) then
            call alpha_SM_write_0(n_SM_max,n_SM_(nm0),alpha_SM__(n_alpha
     $           *n_SM_max*nm0),12,' alpha_SM_: ')
         end if !if (verbose.gt.2) then
         n_LT = 0
         n_ref = max(0,min(n_S,min(n_SM_(nm0),n_LT_ref)))
         do nref=0,n_ref-1
            nx0 = nalpha_S_index + n_alpha*(n_SM_(nm0)-1-nref + n_SM_max
     $           *nm0)
            nx1 = nint(alpha_SM__(nx0))
            ns = nx1
            if (verbose.gt.2) then
               write(6,'(5(A,I4),3(A,F8.4),A,I4)') ' n_ref: ' , n_ref ,
     $              ' nref: ', nref , ' nx0: ' , nx0 , ' nx1: ' , nx1 ,
     $              ' ns: ' , ns , ' S_index_local: ' ,
     $              S_alpha_S_index_(ns) ,' polar_a_local: ' ,
     $              S_alpha_polar_a_(ns) ,' azimu_b_local: ' ,
     $              S_alpha_azimu_b_(ns) ,' I_S_sample: ' ,
     $              I_S_sample_(ns)
            end if !if (verbose.gt.2) then

            call get_angle_to_vp_(S_alpha_polar_a_(ns)
     $           ,S_alpha_azimu_b_(ns),vp_input_omp__(p_vp_input))
            if (verbose.gt.2) then
               write(6,'(A,I0,2(A,F8.4),A,3(F8.4,1X))') ' ns ' , ns ,
     $              ' polar_a: ' , S_alpha_polar_a_(ns) , ' azimu_b: ' ,
     $              S_alpha_azimu_b_(ns) ,' vp_input: ' ,
     $              vp_input_omp__(0 +p_vp_input) ,vp_input_omp__(1
     $              +p_vp_input) , vp_input_omp__(2+p_vp_input)
            end if              !if (verbose.gt.2) then
            tesselation_distance_req_use = tesselation_distance_req
            call tesselation_neighborhood_wrapper_0(nS_0_per,S_L_
     $           ,nl_max,nm_sum,ll_sum,T_nl_,T_vm_,T_tr_ ,T_ll_,T_lf_
     $           ,T_c0_ ,T_c1_,T_c2_,T_c3_,T_ls_,T_LT_,n_T_root_base
     $           ,T_root_base_,vp_input_omp__(p_vp_input)
     $           ,tesselation_distance_req_use,n_LT ,LT_omp__(0+p_LT))
            if (verbose.gt.2) then
               write(6,'(A,I0)') ' n_LT: ' , n_LT
            end if              !if (verbose.gt.2) then

         enddo !do nref=0,n_ref-1

         flag_tesselation_ref = .true.
         include 'ti8l_excerpt_calculate_1.f'

         continue_flag = .false.
         do nLT=0,n_LT-1
            if (flag_S_use_omp__(LT_omp__(nLT+p_LT)
     $           +p_flag_S_use).eqv..true.) then
               if (verbose.gt.1) then
                  write(6,'(A,I0,A,L2,A)') ' nLT: ' , nLT ,
     $                 ' flag_S_use: ' ,flag_S_use_omp__(LT_omp__(nLT
     $                 +p_LT)+p_flag_S_use) , ' skipping... '
               end if !if (verbose.gt.1) then
            end if !if (flag_S_use_omp__(LT_omp__(nLT+p_LT)+p_flag_S_use).eqv..true.) then
            if (flag_S_use_omp__(LT_omp__(nLT+p_LT)
     $           +p_flag_S_use).eqv..false.) then
               if (verbose.gt.1) then
                  write(6,'(A,I0,A,L2,A)') ' nLT: ' , nLT ,
     $                 ' flag_S_use: ' ,flag_S_use_omp__(LT_omp__(nLT
     $                 +p_LT)+p_flag_S_use) , ' flipping... '
               end if !if (verbose.gt.1) then
               continue_flag = .true.
               flag_S_use_omp__(LT_omp__(nLT+p_LT)+p_flag_S_use)=.true.
               n_S_use = n_S_use+1
            end if              !if (flag_S_use_omp__(LT_omp__(nLT+p_LT)+p_flag_S_use).eqv..false.) then
         enddo                  !do nLT=0,n_LT-1
         
         if (0+verbose.gt.1) then
            write(6,'(3(A,I0))') 'n_pass: ' , n_pass ,
     $           '; found n_S_use: ' , n_S_use , ' out of n_S: ' ,
     $           n_S
         end if                 !if (0+verbose.gt.1) then

         n_pass = n_pass + 1

      enddo                     !continue_flag

      n_S_use_sum_(nM_9_sub) = n_S_use_sum_(nM_9_sub) + n_S_use
