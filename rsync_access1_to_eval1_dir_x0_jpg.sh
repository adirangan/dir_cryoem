rsync -avum \
      rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_*_x0 \
      /home/rangan/dir_cryoem/ \
      --include="*/" \
      --include="*snapshot.jpg" \
      --exclude="*" \
      ;
