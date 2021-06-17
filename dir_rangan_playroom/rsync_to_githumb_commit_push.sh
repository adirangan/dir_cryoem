rsync -avum /data/rangan/dir_cryoem/dir_rangan_playroom /data/rangan/dir_cryoem/dir_cryoEM-SCDA/  --include="*/" --include="*.f" --include="*.c" --include="*.make" --include="*.m" --include="*.sh" --include="*.tex" --exclude="*" ;
cd /data/rangan/dir_cryoem/dir_cryoEM-SCDA/ ;
git add dir_rangan_playroom/*.f ;
git add dir_rangan_playroom/*.c ;
git add dir_rangan_playroom/*.make ;
git add dir_rangan_playroom/*.m ;
git add dir_rangan_playroom/*.tex ;
git commit -m "% finished get_loading_qbp_0. ;" ;
git push ;
git pull ;
cd /data/rangan/dir_cryoem/dir_rangan_playroom ;

