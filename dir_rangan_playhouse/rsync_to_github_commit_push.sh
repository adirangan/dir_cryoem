rsync -avum /data/rangan/dir_cryoem/dir_rangan_playhouse /data/rangan/dir_cryoem/dir_cryoEM-SCDA/  --include="*/" --include="*.f" --include="*.c" --include="*.make" --include="*.m" --include="*.sh" --exclude="*" ;
cd /data/rangan/dir_cryoem/dir_cryoEM-SCDA/ ;
git add dir_rangan_playhouse/*.f ;
git add dir_rangan_playhouse/*.c ;
git add dir_rangan_playhouse/*.make ;
git add dir_rangan_playhouse/*.m ;
git add dir_rangan_playhouse/dir_gen_Jsvd_6/*.txt ;
git add dir_rangan_playhouse/dir_gen_Jsvd_6/*.f ;
git commit -m "updating rangan_playhouse to test against ti8_check_SxZTRM_20200222.m " ;
git push ;
git pull ;
cd /data/rangan/dir_cryoem/dir_rangan_playhouse ;

