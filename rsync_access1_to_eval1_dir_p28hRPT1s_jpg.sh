rsync \
    -avum \
    rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_p28hRPT1r \
    /home/rangan/dir_cryoem/ \
    --include="/*" \
    --include="*.jpg" \
    --exclude="*" \
    ;
