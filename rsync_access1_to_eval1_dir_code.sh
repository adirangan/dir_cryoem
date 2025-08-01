rsync \
    -avum \
    avr209@access.cims.nyu.edu:/data/rangan/dir_cryoem/dir_EMIODist2 \
    avr209@access.cims.nyu.edu:/data/rangan/dir_cryoem/dir_CryoBIFE_MD \
    avr209@access.cims.nyu.edu:/data/rangan/dir_cryoem/dir_rangan_playpen \
    avr209@access.cims.nyu.edu:/data/rangan/dir_cryoem/dir_rangan_playroom \
    avr209@access.cims.nyu.edu:/data/rangan/dir_cryoem/dir_rangan_playhouse \
    avr209@access.cims.nyu.edu:/data/rangan/dir_cryoem/dir_rangan_playground \
    avr209@access.cims.nyu.edu:/data/rangan/dir_cryoem/dir_rangan_gpu \
    avr209@access.cims.nyu.edu:/data/rangan/dir_cryoem/dir_rangan_python \
    avr209@access.cims.nyu.edu:/data/rangan/dir_cryoem/dir_presentations \
    avr209@access.cims.nyu.edu:/data/rangan/dir_cryoem/dir_fig_xfig \
    /home/rangan/dir_cryoem/ \
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
    avr209@access.cims.nyu.edu:/data/rangan/dir_cryoem/dir_rangan_playroom/dir_eig_ddssnll_lanczos_local \
    /home/rangan/dir_cryoem/dir_rangan_playroom/ \
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
