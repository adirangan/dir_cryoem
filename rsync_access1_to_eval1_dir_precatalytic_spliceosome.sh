rsync \
    -avum \
    rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_precatalytic_spliceosome \
    rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_ps0 \
    /home/rangan/dir_cryoem/ \
    --exclude "X_2d*" \
    ;
