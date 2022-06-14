rsync -avum \
      /home/rangan/dir_cryoem/ \
      rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/ \
      --include="*/" \
      --include="*X_relion_.*" \
      --exclude="*" \
    ;
