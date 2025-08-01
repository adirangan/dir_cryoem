rsync \
    -avum \
    /home/rangan/dir_cryoem/dir_CryoBIFE_MD \
    avr209@access1.cims.nyu.edu:/data/rangan/dir_cryoem/ \
    --include "*/" \
    --include "*.tex" \
    --include "*.bib" \
    --include "*.pdf" \
    --include "*.fig" \
    --include "*.jpg" \
    --include "*.eps" \
    --include "*.m" \
    --include "*.c" \
    --include "*.h" \
    --include "*.f" \
    --include "*.txt" \
    --include "*.sh" \
    --include "*.csv" \
    --include "*.mrc" \
    --include "*.star" \
    --exclude "*" \
    ;
