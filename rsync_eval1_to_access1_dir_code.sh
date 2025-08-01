rsync \
    -avum \
    /home/rangan/dir_cryoem/dir_EMIODist2 \
    /home/rangan/dir_cryoem/dir_CryoBIFE_MD \
    /home/rangan/dir_cryoem/dir_rangan_playpen \
    /home/rangan/dir_cryoem/dir_rangan_playroom \
    /home/rangan/dir_cryoem/dir_rangan_playhouse \
    /home/rangan/dir_cryoem/dir_rangan_playground \
    /home/rangan/dir_cryoem/dir_rangan_gpu \
    /home/rangan/dir_cryoem/dir_rangan_python \
    /home/rangan/dir_cryoem/dir_presentations \
    /home/rangan/dir_cryoem/dir_fig_xfig \
    avr209@access1.cims.nyu.edu:/data/rangan/dir_cryoem/ \
    --include "/*" \
    --exclude "EMPM/" \
    --include "*.py" \
    --include "*.m" \
    --include "*.c" \
    --include "*.h" \
    --include "*.f" \
    --include "*.txt" \
    --include "*.sh" \
    --include "*.fig" \
    --include "*.tex" \
    --include "*.bib" \
    --exclude "*" \
    ;

rsync \
    -avum \
    /home/rangan/dir_cryoem/dir_rangan_playroom/dir_eig_ddssnll_lanczos_local \
    avr209@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_rangan_playroom/ \
    --include "/*" \
    --exclude "EMPM/" \
    --include "*.py" \
    --include "*.m" \
    --include "*.c" \
    --include "*.h" \
    --include "*.f" \
    --include "*.txt" \
    --include "*.sh" \
    --include "*.fig" \
    --exclude "*" \
    ;
