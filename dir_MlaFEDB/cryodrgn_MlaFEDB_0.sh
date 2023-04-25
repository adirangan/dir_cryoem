# stored in: ;
mkdir /home/rangan/cryodrgn/dir_MlaFEDB ;
cd /home/rangan/cryodrgn/dir_MlaFEDB ;
# conda activate ;
# conda activate cryodrgn ;

# copy over the appropriate mrcs files: ;
<<comment_download
%%%%%%%%;
head -1053 /home/rangan/dir_cryoem/dir_MlaFEDB/Empiar_10536_00_to_23.star > /home/rangan/dir_cryoem/dir_MlaFEDB/Empiar_10536_00_to_23_1024.star ;
%%%%%%%%;
comment_download
# This was suggested by https://github.com/zhonge/cryodrgn, ;
# Modified as the --relion31 flag is not recognized: ;
cryodrgn downsample -D 128 -o /home/rangan/cryodrgn/dir_MlaFEDB/MlaFEDB_cryodrgn_1024.mrcs --datadir /home/rangan/dir_cryoem/dir_MlaFEDB/ /home/rangan/dir_cryoem/dir_MlaFEDB/Empiar_10536_00_to_23_1024.star ;
# generated MlaFEDB_cryodrgn_1024.mrcs (size 65M). ;

# Now trying to create ctf.pkl file. ;
# Adding phase-shift (manually). ;
# Adding voltage (manually). ;
# Adding spherical abberation (manually). ;
# Adding amplitude contrast (manually). ;
# scp -p rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_MlaFEDB/dir_pm_mat/M_k_p__.mat /home/rangan/dir_cryoem/dir_MlaFEDB/dir_pm_mat/ ;
cryodrgn parse_ctf_star /home/rangan/dir_cryoem/dir_MlaFEDB/Empiar_10536_00_to_23_1024.star --ps 0 --kv 300 --cs 0 -w 0 -D 500 --Apix 1.31 -o /home/rangan/cryodrgn/dir_MlaFEDB/MlaFEDB_cryodrgn_ctf.pkl ;
# generated MlaFEDB_cryodrgn_ctf.pkl

# Now trying crodrgn abinit_homo as suggested at: ;
# https://cryodrgn.notion.site/CryoDRGN2-quickstart-322823599fce4bd7a391d00bf749ab1f ;
cryodrgn abinit_homo /home/rangan/cryodrgn/dir_MlaFEDB/MlaFEDB_cryodrgn_1024.mrcs --ctf /home/rangan/cryodrgn/dir_MlaFEDB/MlaFEDB_cryodrgn_ctf.pkl --outdir /home/rangan/cryodrgn/dir_MlaFEDB/dir_job_0 >> cryodrgn_job_0.log ;

scp -p /home/rangan/dir_cryoem/dir_MlaFEDB/cryodrgn_MlaFEDB_0.sh rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_MlaFEDB/cryodrgn_MlaFEDB_0.sh ;
