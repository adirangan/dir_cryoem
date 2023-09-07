rsync -avum \
      rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_CryoBIFE_MD \
      /home/rangan/dir_cryoem/ \
    --include "/*" \
    --include "*.m" \
    --include "*.c" \
    --include "*.h" \
    --include "*.f" \
    --include "*.txt" \
    --include "*.sh" \
    --include "*.tex" \
    --include "*.fig" \
    --include "*.bib" \
    --include "*.pdf" \
    --exclude "*" \
    ;
