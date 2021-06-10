c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      test_alpha_update_SM_0 <-- test_alpha_update_0
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      call rebuild_model_0 !rebuild model Y_est_ from alpha1d_est_
      call template_gen !generate templates S__ from Y_est_
      call getspheregrid !use spherical grid to define S_alpha_polar_a_all_, S_alpha_azimu_b_all_
      !define n_S_sample_
      call test_innerproduct_batch_SM_wrapper_0a !calculate innerproducts
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      test_innerproduct_batch_SM_wrapper_0a <-- test_innerproduct_batch_wrapper_0
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      do nS=0,n_S_sample_
         call get_angle_to_vp_ !create list of viewing-vectors L_ associated with templates
      enddo
      call tesselation_get_nl_nm_ll !allocate tesselation T_ using L_
      call tesselation_index_wrapper_0 !build tesselation T_ from L_
      call test_innerproduct_batch_SM_stage_0a ! batch innerproduct calculation for SM (passing T_)
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      test_innerproduct_batch_SM_stage_0a <-- test_innerproduct_batch_stage_0
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      call innerproduct_q_k_svdd !define Z_svdd_
      do ns=0,n_S-1 do nctf=0,n_CTF-1
         call test_innerproduct_batch_excerpt_0 !find C_S_ = \| CTF.*R(S) \| for each ctf-S pair
      enddo
      do nM_sub=0,n_M_sub-1 !omp parallel do across 8 blocks (i.e., many images in each)
         call test_innerproduct_batch_SM_stage_1a !calculate innerproducts for block nM_sub using T_
      enddo
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      test_innerproduct_batch_SM_stage_1a <-- test_innerproduct_batch_stage_1
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !allocate LT_ to hold list of unique template indices
      do nm=0,n_M-1 !step through images
         call tesselation_neighborhood_wrapper_0 !build LT_ using templates close to image nm
         call randinclude !add some other random templates
         do nLT=0,n_LT-1
            call test_innerproduct_batch_stage_2a ! <-- test_innerproduct_batch_stage_2 ! calculate innerproducts for image-template pair (nm,nLT)
         enddo !nLT=0,n_LT-1
         do while
            !determine which templates are 'good'
            !i.e., which templates are in the top n_LT_hlf of all calculated image-template pairs
            call tesselation_neighborhood_wrapper_0 !add to LT_ any templates nearby to these good templates
            call test_innerproduct_batch_stage_2a ! <-- test_innerproduct_batch_stage_2 ! calculate innerproducts for these newly added nearby templates
            !continue until no new nearby templates can be added
         enddo !while
         do ns=0,n_S-1
            if (S_use_(ns).eqv..false.) then
               C_Z_optimal_(ns) = C_Z_min !fill uncalculated image-template pairs with small values
            end if ! if (S_use_(ns).eqv..false.) then
         enddo !do ns=0,n_S-1
         call test_innerproduct_batch_sort_a !sort C_Z_optimal_
         call get_interchange_delta !determine update to delta_x and delta_y
         !update image parameters
      enddo !nm=0,n_M-1

c$$$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      test_ver9.f
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      call createfun_xspace !create f_x_c_A___ and f_x_c_B___
      call getgridr !create grid_k_p_
      call get_ntot !calculate n_azimu_b_polar_a_k_p_sum
      call fun_to_kspacegrid !create f_k_c_A_ and f_k_c_B_
      call kspacegrid_to_model !create Y_A_ and Y_B_ 
      call getspheregrid !create grid_cos_polar_a, grid_sin_polar_a
      call mk_simulated_slices_alpha2d !create Y_slice_A__, Y_slice_B__
      do nctf=0,n_CTF-1
         call get_ctf_star_k_p_ !define CTF_
      enddo
      do n_k_cur=n_k_bot,n_k_p_max
         !copy Y_slice_A__ and Y_slice_B__ to image stack M_sample_
         !copy alpha2d_est__ to alpha1d_est_
         call test_alpha_update_1 !update model
      enddo

c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      test_alpha_update_1.f
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      call rebuild_model_0 !rebuild model Y_est_ from alpha1d_est_
      call template_gen !generate templates S__ from Y_est_
      call getspheregrid !use spherical grid to define S_alpha_polar_a_all_, S_alpha_azimu_b_all_
      !define n_S_sample_
      call test_innerproduct_batch_wrapper_0 !calculate innerproducts
      call test_innerproduct_batch_wrapper_1 !sort innerproducts
      call test_innerproduct_batch_wrapper_2 !update alphas

c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      test_innerproduct_batch_wrapper_0
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      do nS=0,n_S_sample_
         call get_angle_to_vp_ !create list of viewing-vectors L_ associated with templates
      enddo
      call tesselation_get_nl_nm_ll !allocate tesselation T_ using L_
      call tesselation_index_wrapper_0 !build tesselation T_ from L_
      call test_innerproduct_batch_stage_0 !batch innerproduct calculation (passing T_)

c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      test_innerproduct_batch_stage_0
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      call innerproduct_q_k_svdd !define Z_svdd_
      do ns=0,n_S-1 do nctf=0,n_CTF-1
         call test_innerproduct_batch_excerpt_0 !find C_S_ = \| CTF.*R(S) \| for each ctf-S pair
      enddo
      do nM_sub=0,n_M_sub-1 !omp parallel do
         call test_innerproduct_batch_stage_1 !calculate innerproducts for block nM_sub using T_
      enddo

c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      test_innerproduct_batch_stage_1
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !allocate LT_ to hold list of unique template indices
      do nm=0,n_M-1 !step through images
         call tesselation_neighborhood_wrapper_0 !build LT_ using templates close to image nm
         call randinclude !add some other random templates
         do nLT=0,n_LT-1
            call test_innerproduct_batch_stage_2 ! calculate innerproducts for image-template pair (nm,nLT)
         enddo !nLT=0,n_LT-1
         do while
            !determine which templates are 'good'
            !i.e., which templates are in the top n_LT_hlf of all calculated image-template pairs
            call tesselation_neighborhood_wrapper_0 !add to LT_ any templates nearby to these good templates
            call test_innerproduct_batch_stage_2 ! calculate innerproducts for these newly added nearby templates
            !continue until no new nearby templates can be added
         enddo !while
         do ns=0,n_S-1
            if (S_use_(ns).eqv..false.) then
               C_Z_optimal_(ns + nm*n_S) = C_Z_min !fill uncalculated image-template pairs with small values
            end if ! if (S_use_(ns).eqv..false.) then
         enddo !do ns=0,n_S-1
      enddo !nm=0,n_M-1

c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      test_innerproduct_batch_stage_2
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !calculate single innerproduct for image-template pair
      !writes optimal values to C_S_optimal_(ns+nm*n_S), etc. 


