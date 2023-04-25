# stored in: ;
cd /home/rangan/cryodrgn/dir_ISWINCP ;
# conda activate ;
# conda activate cryodrgn ;

# copy over the appropriate mrcs files: ;
<<comment_download
%%%%%%%%;
scp -p rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_ISWINCP/ADPBeF.star /home/rangan/dir_cryoem/dir_ISWINCP/ ;
head -1055 /home/rangan/dir_cryoem/dir_ISWINCP/ADPBeF.star > /home/rangan/dir_cryoem/dir_ISWINCP/ADPBeF_1024.star ;
trash /home/rangan/dir_cryoem/dir_ISWINCP/ADPBeF.star;
%%%%%%%%;
mkdir /home/rangan/dir_cryoem/dir_ISWINCP/Extract ;
mkdir /home/rangan/dir_cryoem/dir_ISWINCP/Extract/job204 ;
mkdir /home/rangan/dir_cryoem/dir_ISWINCP/Extract/job204/mics ;
%%%%%%%%;
scp -p rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_ISWINCP/Extract/job204/mics/stack_00[0-4]?_motion2_DW.mrcs /home/rangan/dir_cryoem/dir_ISWINCP/Extract/job204/mics/ ;
%%%%%%%%;
comment_download
# This was suggested by https://github.com/zhonge/cryodrgn, ;
# Modified as the --relion31 flag is not recognized: ;
cryodrgn downsample -D 128 -o /home/rangan/cryodrgn/dir_ISWINCP/ISWINCP_cryodrgn_1024.mrcs --datadir /home/rangan/dir_cryoem/dir_ISWINCP/ /home/rangan/dir_cryoem/dir_ISWINCP/ADPBeF_1024.star ;
# generated ISWINCP_cryodrgn_1024.mrcs (size 65M). ;

# Now trying to create ctf.pkl file. ;
# Adding phase-shift (manually). ;
# scp -p rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_ISWINCP/dir_pm_mat/M_k_p__.mat /home/rangan/dir_cryoem/dir_ISWINCP/dir_pm_mat/ ;
cryodrgn parse_ctf_star /home/rangan/dir_cryoem/dir_ISWINCP/ADPBeF_1024.star --ps 0 -D 240 --Apix 1.07 -o /home/rangan/cryodrgn/dir_ISWINCP/ISWINCP_cryodrgn_ctf.pkl ;
# generated ISWINCP_cryodrgn_ctf.pkl

# Now trying crodrgn abinit_homo as suggested at: ;
# https://cryodrgn.notion.site/CryoDRGN2-quickstart-322823599fce4bd7a391d00bf749ab1f ;
cryodrgn abinit_homo /home/rangan/cryodrgn/dir_ISWINCP/ISWINCP_cryodrgn_1024.mrcs --ctf /home/rangan/cryodrgn/dir_ISWINCP/ISWINCP_cryodrgn_ctf.pkl --outdir /home/rangan/cryodrgn/dir_ISWINCP/dir_job_0 >> cryodrgn_job_0.log ;

scp -p /home/rangan/dir_cryoem/dir_ISWINCP/cryodrgn_ISWINCP_0.sh rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_ISWINCP/cryodrgn_ISWINCP_0.sh ;
