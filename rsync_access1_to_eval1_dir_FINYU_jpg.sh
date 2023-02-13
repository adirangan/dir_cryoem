rsync \
    -avum \
    rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_FINYU_nx64d0150s005nM1024nCTF20p40V300D20000r025A010 \
    /home/rangan/dir_cryoem/ \
    --include="/*" \
    --include="*.jpg" \
    --exclude="*.mrcs" \
    --exclude="*.mrc" \
    --exclude="*.mat" \
    --exclude="*.eps" \
    ;
