rsync -avum \
      rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_m_dependencies \
      ~/dir_cryoem/ \
      --include="*/" \
      --include="*.f" \
      --include="*.c" \
      --include="*.make" \
      --include="*.m" \
      --include="*.sh" \
      --include="*.tex" \
      --include="*.pdf" \
      --exclude="*" ;
rsync -avum \
      rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_rangan_play* \
      ~/dir_cryoem/ \
      --include="*/" \
      --include="*.f" \
      --include="*.c" \
      --include="*.make" \
      --include="*.m" \
      --include="*.sh" \
      --include="*.tex" \
      --include="*.pdf" \
      --exclude="*" ;
