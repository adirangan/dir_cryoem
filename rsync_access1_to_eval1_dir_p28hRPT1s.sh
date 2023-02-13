rsync \
    -avum \
    rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_p28hRP \
    rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_p28hRPT1 \
    rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_p28hRPT2 \
    rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_p28hRPT3 \
    rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_p28hRPT1s \
    /home/rangan/dir_cryoem/ \
    --exclude "*.tmp" \
    ;
