rsync -avum \
      rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_p28hRPT1/dir_pm_mat \
      /home/rangan/dir_cryoem/dir_p28hRPT1/ \
      --include="/*" \
      --include="*.jpg" \
      --exclude="*" \
;
