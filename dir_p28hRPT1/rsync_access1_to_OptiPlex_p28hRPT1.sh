rsync -avum \
      rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_p28hRPT1/ \
      /home/rangan/dir_cryoem/dir_p28hRPT1/ \
      --include="*align*.mat" \
      --exclude="*X_2d*.mat" \
    ;
