!> Doxygen comment: ;\n
!> test_innerproduct_8 excerpt: ;\n
!> Used for printing block strategy to screen. ;\n
      if (verbose_timing.gt.0) then
      write(6,'(16(A,I0))') 
     $     ' nS_0_sub: '
     $     ,nS_0_sub
     $     ,'/'
     $     ,n_S_0_sub_use
     $     , '; nS_0_per: '
     $     ,nS_0_per
     $     ,'; nS_0_sum: '
     $     ,nS_0_sum 
     $     ,'; nM_0_sub: '
     $     ,nM_0_sub
     $     ,'/'
     $     ,n_M_0_sub_use 
     $     ,'; nM_0_per: '
     $     ,nM_0_per
     $     ,'; nM_0_sum: '
     $     ,nM_0_sum
     $     ,' nS_1_sub: '
     $     ,nS_1_sub
     $     ,'/'
     $     ,n_S_1_sub_use
     $     , '; nS_1_per: '
     $     ,nS_1_per
     $     ,'; nS_1_sum: '
     $     ,nS_1_sum 
     $     ,'; nM_1_sub: '
     $     ,nM_1_sub
     $     ,'/'
     $     ,n_M_1_sub_use 
     $     ,'; nM_1_per: '
     $     ,nM_1_per
     $     ,'; nM_1_sum: '
     $     ,nM_1_sum
      end if !if (verbose_timing.gt.0) then
