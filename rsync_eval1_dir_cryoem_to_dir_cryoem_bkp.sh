rsync -avum \
      ~/dir_cryoem/* \
      ~/dir_cryoem_bkp/ \
      --include="*/" \
      --include="*.f" \
      --include="*.c" \
      --include="*.make" \
      --include "*.org" \
      --include="*.py" \
      --include="*.m" \
      --include="*.sh" \
      --include="*.tex" \
      --include="*.pdf" \
      --exclude="*" ;
