# stored in: ;
mkdir /home/rangan/cryodrgn/dir_LetB1 ;
cd /home/rangan/cryodrgn/dir_LetB1 ;
# conda activate ;
# conda activate cryodrgn ;

# copy over the appropriate mrcs files: ;
<<comment_download
%%%%%%%%;
scp -p rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_LetB1/job_569_model_1_350MB.star /home/rangan/dir_cryoem/dir_LetB1/ ;
head -1059 /home/rangan/dir_cryoem/dir_LetB1/job_569_model_1_350MB.star > /home/rangan/dir_cryoem/dir_LetB1/job_569_model_1_350MB_1024.star ;
trash /home/rangan/dir_cryoem/dir_LetB1/job_569_model_1_350MB.star;
%%%%%%%%;
scp -rp rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_LetB1/Particles_Stacks_Model_1 /home/rangan/dir_cryoem/dir_LetB1/ ;
%%%%%%%%;
comment_download
# This was suggested by https://github.com/zhonge/cryodrgn, ;
# Modified as the --relion31 flag is not recognized: ;
cryodrgn downsample -D 128 -o /home/rangan/cryodrgn/dir_LetB1/LetB1_cryodrgn_1024.mrcs --datadir /home/rangan/dir_cryoem/dir_LetB1/ /home/rangan/dir_cryoem/dir_LetB1/job_569_model_1_350MB_1024.star ;
# generated LetB1_cryodrgn_1024.mrcs (size 65M). ;

# Now trying to create ctf.pkl file. ;
# Adding phase-shift (manually). ; %<-- Note: phase-shift in star file is zero. ;
# scp -p rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_LetB1/dir_pm_mat/M_k_p__.mat /home/rangan/dir_cryoem/dir_LetB1/dir_pm_mat/ ;
cryodrgn parse_ctf_star /home/rangan/dir_cryoem/dir_LetB1/job_569_model_1_350MB_1024.star --ps 0 -D 280 --Apix 1.31 -o /home/rangan/cryodrgn/dir_LetB1/LetB1_cryodrgn_ctf.pkl ;
# generated LetB1_cryodrgn_ctf.pkl

# Now trying crodrgn abinit_homo as suggested at: ;
# https://cryodrgn.notion.site/CryoDRGN2-quickstart-322823599fce4bd7a391d00bf749ab1f ;
cryodrgn abinit_homo /home/rangan/cryodrgn/dir_LetB1/LetB1_cryodrgn_1024.mrcs --ctf /home/rangan/cryodrgn/dir_LetB1/LetB1_cryodrgn_ctf.pkl --seed 0 --outdir /home/rangan/cryodrgn/dir_LetB1/dir_job_0 >> cryodrgn_job_0.log ;

cryodrgn abinit_homo /home/rangan/cryodrgn/dir_LetB1/LetB1_cryodrgn_1024.mrcs --ctf /home/rangan/cryodrgn/dir_LetB1/LetB1_cryodrgn_ctf.pkl --seed 1 --outdir /home/rangan/cryodrgn/dir_LetB1/dir_job_1 >> cryodrgn_job_1.log ;

cryodrgn abinit_homo /home/rangan/cryodrgn/dir_LetB1/LetB1_cryodrgn_1024.mrcs --ctf /home/rangan/cryodrgn/dir_LetB1/LetB1_cryodrgn_ctf.pkl --seed 2 --outdir /home/rangan/cryodrgn/dir_LetB1/dir_job_2 >> cryodrgn_job_2.log ;

cryodrgn abinit_homo /home/rangan/cryodrgn/dir_LetB1/LetB1_cryodrgn_1024.mrcs --ctf /home/rangan/cryodrgn/dir_LetB1/LetB1_cryodrgn_ctf.pkl --seed 3 --outdir /home/rangan/cryodrgn/dir_LetB1/dir_job_3 >> cryodrgn_job_3.log ;

scp -p /home/rangan/dir_cryoem/dir_LetB1/cryodrgn_LetB1_0.sh rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_LetB1/cryodrgn_LetB1_0.sh ;
