rsync \
    -avum \
    /home/rangan/dir_cryoem/dir_EMIODist2 \
    /home/rangan/dir_cryoem/dir_rangan_playpen \
    /home/rangan/dir_cryoem/dir_rangan_playroom \
    /home/rangan/dir_cryoem/dir_rangan_playhouse \
    /home/rangan/dir_cryoem/dir_rangan_playground \
    /home/rangan/dir_cryoem/dir_presentations \
    /home/rangan/dir_cryoem/dir_fig_xfig \
    rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/ \
    --include "/*" \
    --include "*.m" \
    --include "*.c" \
    --include "*.h" \
    --include "*.f" \
    --include "*.txt" \
    --include "*.sh" \
    --include "*.fig" \
    --exclude "*" \
    ;
