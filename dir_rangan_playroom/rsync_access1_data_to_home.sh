rsync -avum \
      /data/rangan/dir_cryoem/dir_rangan_playroom \
      /home/rangan/dir_cryoem/ \
      --include="*/" \
      --include="*.sh" \
      --include="*.m" \
      --include="*.c" \
      --include="*.h" \
      --include="*.f" \
      --include="*.tex" \
      --include="*.bib" \
      --include="*.pdf" \
      --include="*.fig" \
      --exclude="*" \
      ;
