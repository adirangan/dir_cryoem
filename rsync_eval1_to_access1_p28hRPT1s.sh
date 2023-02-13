rsync -avum \
      /home/rangan/dir_cryoem/dir_p28hRP \
      /home/rangan/dir_cryoem/dir_p28hRPT1 \
      /home/rangan/dir_cryoem/dir_p28hRPT2 \
      /home/rangan/dir_cryoem/dir_p28hRPT3 \
      /home/rangan/dir_cryoem/dir_p28hRPT1s \
      /home/rangan/dir_cryoem/dir_p28hRPT1r \
      rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/ \
      --exclude="*.mrcs" \
      --exclude="S_k_p__.mat" \
      --exclude="UX_Memp_N_k_p_wnM___.mat" \
      --exclude="X_2d_*.mat" \
      --exclude="X_2d_*.jpg" \
      --exclude="X_2d_*.eps" \
      --exclude="*.eps" \
      --exclude="*.tmp" \
      --exclude="*~" \
    ;
