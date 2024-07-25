rsync \
    -avum \
    rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_EMIODist2 \
    rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_rangan_playpen \
    rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_rangan_playroom \
    rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_rangan_playhouse \
    rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_rangan_playground \
    rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_presentations \
    rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_fig_xfig \
    /home/rangan/dir_cryoem/ \
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

rsync \
    -avum \
    rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_rangan_playroom/dir_eig_ddssnll_lanczos_local \
    /home/rangan/dir_cryoem/dir_rangan_playroom/ \
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
