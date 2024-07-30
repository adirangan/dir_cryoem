rsync \
    -avum \
    rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_trpv1_x0/dir_ssnll* \
    /home/rangan/dir_cryoem/dir_trpv1_x0/ \
    --include "/*" \
    --include "/*/*" \
    --include "*.m" \
    --include "*.c" \
    --include "*.h" \
    --include "*.f" \
    --include "*.txt" \
    --include "*.sh" \
    --include "*.fig" \
    --exclude "*" \
    ;
