      tmp_dir = '/data/rangan/dir_cryoem/dir_rangan_playground/dir_tmp/'
c$$$  %%%%%%%%
      tmp_str = trim(adjustl(tmp_dir)) // 'nk'
      open(57,FILE=tmp_str)
      write(57,*) nk
      close(57)
c$$$  %%%%%%%%
      tmp_str = trim(adjustl(tmp_dir)) // 'n_w_M'
      open(57,FILE=tmp_str)
      write(57,*) n_w_M
      close(57)
c$$$  %%%%%%%%
      tmp_str = trim(adjustl(tmp_dir)) // 'M_k_c_'
      open(57,FILE=tmp_str)
      do tmp_tab=0,n_w_M-1
         write(57,*) M_k_c_(tmp_tab)
      enddo                     !do tmp_tab=0,n_w_M-1
      close(57)
c$$$  %%%%%%%%
      tmp_str = trim(adjustl(tmp_dir)) // 'polar_a_'
      open(57,FILE=tmp_str)
      do tmp_tab=0,n_w_M-1
         write(57,*) polar_a_(tmp_tab)
      enddo                     !do tmp_tab=0,n_w_M-1
      close(57)
c$$$  %%%%%%%%
      tmp_str = trim(adjustl(tmp_dir)) // 'azimu_b_'
      open(57,FILE=tmp_str)
      do tmp_tab=0,n_w_M-1
         write(57,*) azimu_b_(tmp_tab)
      enddo                     !do tmp_tab=0,n_w_M-1
      close(57)
c$$$  %%%%%%%%
      tmp_str = trim(adjustl(tmp_dir)) // 'n_Y_l'
      open(57,FILE=tmp_str)
      write(57,*) n_Y_l
      close(57)
c$$$  %%%%%%%%
      tmp_str = trim(adjustl(tmp_dir)) // 'lsq_n_quad'
      open(57,FILE=tmp_str)
      write(57,*) lsq_n_quad
      close(57)
c$$$  %%%%%%%%
      tmp_str = trim(adjustl(tmp_dir)) //
     $     'lsq_interpolation_order'
      open(57,FILE=tmp_str)
      write(57,*) lsq_interpolation_order
      close(57)
c$$$  %%%%%%%%
      tmp_str = trim(adjustl(tmp_dir)) // 'weight_CTF_k_c_'
      open(57,FILE=tmp_str)
      do tmp_tab=0,n_w_M-1
         write(57,*) weight_CTF_k_c_(tmp_tab)
      enddo                     !do tmp_tab=0,n_w_M-1
      close(57)
c$$$  %%%%%%%%
      tmp_str = trim(adjustl(tmp_dir)) // 'lsq_eps'
      open(57,FILE=tmp_str)
      write(57,*) lsq_eps
      close(57)
c$$$  %%%%%%%%
      tmp_str = trim(adjustl(tmp_dir)) // 'lsq_n_iteration'
      open(57,FILE=tmp_str)
      write(57,*) lsq_n_iteration
      close(57)
c$$$  %%%%%%%%
      
