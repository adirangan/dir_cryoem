rsync -avum \
      ~/dir_cryoem/dir_ampm_dev \
      rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/ \
      --include "*/" \
      --include "*.sh" \
      --include "*.make" \
      --include "*.c" \
      --include "*.h" \
      --include "*.in" \
      --include "*.m" \
      --exclude "*" \
      ;
