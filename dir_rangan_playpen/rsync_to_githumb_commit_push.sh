rsync -avum /data/rangan/dir_cryoem/dir_rangan_playpen /data/rangan/dir_cryoem/dir_cryoEM-SCDA/  --include="*/" --include="*.f" --include="*.c" --include="*.make" --include="*.m" --include="*.sh" --include="*.tex" --exclude="*" ;
cd /data/rangan/dir_cryoem/dir_cryoEM-SCDA/ ;
git add dir_rangan_playpen/*.f ;
git add dir_rangan_playpen/*.c ;
git add dir_rangan_playpen/*.make ;
git add dir_rangan_playpen/*.m ;
git add dir_rangan_playpen/*.tex ;
git commit -m "% see test_principal_marching_trpv1_12.m for translation testing. ; " ;
git push ;
git pull ;
cd /data/rangan/dir_cryoem/dir_rangan_playpen ;

