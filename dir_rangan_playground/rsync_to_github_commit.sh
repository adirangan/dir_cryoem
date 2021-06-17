#rsync -avum /data/rangan/dir_cryoem/dir_rangan_playpen /data/rangan/dir_cryoem/dir_cryoEM-SCDA/  --include="*/" --include="*.f" --include="*.c" --include="*.make" --include="*.m" --include="*.tex" --include="*.sh" --exclude="*" ;
#rsync -avum /data/rangan/dir_cryoem/dir_rangan_playroom /data/rangan/dir_cryoem/dir_cryoEM-SCDA/  --include="*/" --include="*.f" --include="*.c" --include="*.make" --include="*.m" --include="*.tex" --include="*.sh" --exclude="*" ;
#rsync -avum /data/rangan/dir_cryoem/dir_rangan_playhouse /data/rangan/dir_cryoem/dir_cryoEM-SCDA/  --include="*/" --include="*.f" --include="*.c" --include="*.make" --include="*.m" --include="*.tex" --include="*.sh" --exclude="*" ;
rsync -avum /data/rangan/dir_cryoem/dir_rangan_playground /data/rangan/dir_cryoem/dir_cryoEM-SCDA/  --include="*/" --include="*.f" --include="*.c" --include="*.make" --include="*.m" --include="*.txt" --include="*.sh" --exclude="*" ;
cd /data/rangan/dir_cryoem/dir_cryoEM-SCDA/ ;
git add dir_rangan_playground/*.f ;
git add dir_rangan_playground/*.c ;
git add dir_rangan_playground/*.make ;
git add dir_rangan_playground/*.m ;
git add dir_rangan_playground/dir_gen_Jsvd_6/*.txt ;
git add dir_rangan_playground/dir_gen_Jsvd_6/*.f ;
git commit -m "updating rangan_playground" ;
git pull ;
cd /data/rangan/dir_cryoem/dir_rangan_playground ;

