rsync -avum \
      ~/dir_cryoem/* \
      ~/dir_cryoem_bkp/ \
      --include="*/" \
      --include="*.f" \
      --include="*.c" \
      --include="*.make" \
      --include="*.m" \
      --include="*.sh" \
      --include="*.txt" \
      --include="*.tex" \
      --include="*.pdf" \
      --exclude="*" ;
