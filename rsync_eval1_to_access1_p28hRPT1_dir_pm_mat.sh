rsync -avum \
      /home/rangan/dir_cryoem/dir_p28hRPT1_x0 \
      rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/ \
      --exclude="*.mrcs" \
      --exclude="*~" \
    ;
