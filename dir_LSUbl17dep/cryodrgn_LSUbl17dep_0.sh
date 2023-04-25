# stored in: ;
mkdir /home/rangan/cryodrgn/dir_LSUbl17dep ;
cd /home/rangan/cryodrgn/dir_LSUbl17dep ;
# conda activate ;
# conda activate cryodrgn ;

# copy over the appropriate mrcs files: ;
<<comment_download
%%%%%%%%;
head -1037 /home/rangan/dir_cryoem/dir_LSUbl17dep/Parameters_negated.star > /home/rangan/dir_cryoem/dir_LSUbl17dep/Parameters_negated_1024.star ;
%%%%%%%%;
comment_download
# This was suggested by https://github.com/zhonge/cryodrgn, ;
# Modified as the --relion31 flag is not recognized: ;
cryodrgn downsample -D 128 -o /home/rangan/cryodrgn/dir_LSUbl17dep/LSUbl17dep_cryodrgn_1024.mrcs --datadir /home/rangan/dir_cryoem/dir_LSUbl17dep/ /home/rangan/dir_cryoem/dir_LSUbl17dep/Parameters_negated_1024.star ;
# generated LSUbl17dep_cryodrgn_1024.mrcs (size 65M). ;

# Now trying to create ctf.pkl file. ;
# Adding phase-shift (manually). ;
# scp -p rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_LSUbl17dep/dir_pm_mat/M_k_p__.mat /home/rangan/dir_cryoem/dir_LSUbl17dep/dir_pm_mat/ ;
cryodrgn parse_ctf_star /home/rangan/dir_cryoem/dir_LSUbl17dep/Parameters_negated_1024.star --ps 0 -D 320 --Apix 1.31 -o /home/rangan/cryodrgn/dir_LSUbl17dep/LSUbl17dep_cryodrgn_ctf.pkl ;
# generated LSUbl17dep_cryodrgn_ctf.pkl

# Now trying crodrgn abinit_homo as suggested at: ;
# https://cryodrgn.notion.site/CryoDRGN2-quickstart-322823599fce4bd7a391d00bf749ab1f ;
cryodrgn abinit_homo /home/rangan/cryodrgn/dir_LSUbl17dep/LSUbl17dep_cryodrgn_1024.mrcs --ctf /home/rangan/cryodrgn/dir_LSUbl17dep/LSUbl17dep_cryodrgn_ctf.pkl --outdir /home/rangan/cryodrgn/dir_LSUbl17dep/dir_job_0 >> cryodrgn_job_0.log ;

cryodrgn abinit_homo /home/rangan/cryodrgn/dir_LSUbl17dep/LSUbl17dep_cryodrgn_1024.mrcs --ctf /home/rangan/cryodrgn/dir_LSUbl17dep/LSUbl17dep_cryodrgn_ctf.pkl --seed 1 --outdir /home/rangan/cryodrgn/dir_LSUbl17dep/dir_job_1 >> cryodrgn_job_1.log ;

cryodrgn abinit_homo /home/rangan/cryodrgn/dir_LSUbl17dep/LSUbl17dep_cryodrgn_1024.mrcs --ctf /home/rangan/cryodrgn/dir_LSUbl17dep/LSUbl17dep_cryodrgn_ctf.pkl --seed 2 --outdir /home/rangan/cryodrgn/dir_LSUbl17dep/dir_job_2 >> cryodrgn_job_2.log ;

cryodrgn abinit_homo /home/rangan/cryodrgn/dir_LSUbl17dep/LSUbl17dep_cryodrgn_1024.mrcs --ctf /home/rangan/cryodrgn/dir_LSUbl17dep/LSUbl17dep_cryodrgn_ctf.pkl --seed 3 --outdir /home/rangan/cryodrgn/dir_LSUbl17dep/dir_job_3 >> cryodrgn_job_3.log ;

scp -p /home/rangan/dir_cryoem/dir_LSUbl17dep/cryodrgn_LSUbl17dep_0.sh rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_LSUbl17dep/cryodrgn_LSUbl17dep_0.sh ;
