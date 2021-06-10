!> Doxygen comment: ;\n
!> test_innerproduct_8 excerpt: ;\n
!> Displays indices within tesselation-tree. ;\n
      write(6,'(A,I0)') ' n_LT: ' , n_LT
      do nLT=0,n_LT-1
         npoint = LT_omp__(nLT+p_LT)
         write(6,'(3(A,I0),A,3(F8.4,1X),A,F8.4)') ' nLT: ' , nLT ,
     $        ' npoint: ',npoint , ' I_S_sample: ' ,
     $        I_S_sample_(nS_0_sum +npoint), ' vp: ' , S_L_(0 +3*npoint)
     $        , S_L_(1+3*npoint) , S_L_(2+3*npoint) , ' distance: ' ,
     $        dsqrt( (vp_input_omp__(0+p_vp_input) - S_L_(0 + 3*npoint))
     $        **2 + (vp_input_omp__(1+p_vp_input) - S_L_(1 + 3*npoint))
     $        **2 + (vp_input_omp__(2+p_vp_input) - S_L_(2 + 3*npoint))
     $        **2)
      enddo                     !do nLT=0,n_LT-1
