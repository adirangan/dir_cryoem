#rsync -avum rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_nufftall-1.33 /home/rangan/dir_cryoem/ ;

rsync \
    -avum \
    rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_frtalign \
    /home/rangan/dir_cryoem/ \
    --include "/*" \
    --include "*.m" \
    --include "*.c" \
    --include "*.h" \
    --include "*.f" \
    --include "*.txt" \
    --include "*.tex" \
    --include "*.bib" \
    --include "*.sty" \
    --include "*.cls" \
    --include "*.clo" \
    --include "*.bst" \
    --include "*.sh" \
    --exclude "*" \
    ;
