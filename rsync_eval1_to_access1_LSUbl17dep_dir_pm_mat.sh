rsync -avum \
      /home/rangan/dir_cryoem/dir_LSUbl17dep \
      /home/rangan/dir_cryoem/dir_LSUbl17dep_x0 \
      rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/ \
      --include="L17Combine_weight_local_negated.mrcs" \
      --exclude="*.mrcs" \
      --exclude="*~" \
    ;
