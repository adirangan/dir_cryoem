rsync \
    -avum \
    rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_LSUbl17dep \
    rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_LSUbl17dep_x0 \
    /home/rangan/dir_cryoem/ \
    --exclude "X_2d*" \
    ;
