rsync -avum /data/rangan/dir_cryoem/dir_rangan_playground /data/rangan/dir_cryoem/dir_cryoEM-SCDA/  --include="*/" --include="*.f" --include="*.c" --include="*.make" --include="*.m" --include="*.txt" --include="*.sh" --include="*.tex" --exclude="*" ;
cd /data/rangan/dir_cryoem/dir_cryoEM-SCDA/ ;
git add dir_rangan_playground/*.f ;
git add dir_rangan_playground/*.c ;
git add dir_rangan_playground/*.make ;
git add dir_rangan_playground/*.m ;
git add dir_rangan_playground/*.sh ;
git add dir_rangan_playground/*.tex ;
git add dir_rangan_playground/dir_gen_Jsvd_6/*.txt ;
git add dir_rangan_playground/dir_gen_Jsvd_6/*.f ;
git commit -m " starting to match ti8_bestfitall_dr.f to recursive_lin.f. ;" ;
git push ;
git pull ;
cd /data/rangan/dir_cryoem/dir_rangan_playground ;

