rsync -avum \
      /home/rangan/dir_cryoem/dir_LSUbl17dep \
      rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/ \
      --include="L17Combine_weight_local_negated.mrcs" \
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
