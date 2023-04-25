# stored in: ;
cd /home/rangan/cryodrgn/dir_TMEM16F ;
# conda activate ;
# conda activate cryodrgn ;

# copy over the appropriate mrcs files: ;
<<comment_download
%%%%%%%%;
scp -p rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_TMEM16F/All_8192.star /home/rangan/dir_cryoem/dir_TMEM16F/ ;
head -1045 /home/rangan/dir_cryoem/dir_TMEM16F/All_8192.star > /home/rangan/dir_cryoem/dir_TMEM16F/All_8192_1024.star ;
%%%%%%%%;
mkdir /home/rangan/dir_cryoem/dir_TMEM16F/Extract ;
mkdir /home/rangan/dir_cryoem/dir_TMEM16F/Extract/job003 ;
mkdir /home/rangan/dir_cryoem/dir_TMEM16F/Extract/job003/Micrographs ;
%%%%%%%%;
scp -p rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_TMEM16F/Extract/job003/Micrographs/16F_Digitonin_Ca_000[0-4]_DW.mrcs /home/rangan/dir_cryoem/dir_TMEM16F/Extract/job003/Micrographs/ ;
%%%%%%%%;
comment_download
# This was suggested by https://github.com/zhonge/cryodrgn, ;
# Modified as the --relion31 flag is not recognized: ;
cryodrgn downsample -D 128 -o /home/rangan/cryodrgn/dir_TMEM16F/TMEM16F_cryodrgn_1024.mrcs --datadir /home/rangan/dir_cryoem/dir_TMEM16F/ /home/rangan/dir_cryoem/dir_TMEM16F/All_8192_1024.star ;
# generated TMEM16F_cryodrgn_1024.mrcs (size 65M). ;

# Now trying to create ctf.pkl file. ;
# Note, loading phase-shift. ;
# scp -p rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_TMEM16F/dir_pm_mat/M_k_p__.mat /home/rangan/dir_cryoem/dir_TMEM16F/dir_pm_mat/ ;
cryodrgn parse_ctf_star /home/rangan/dir_cryoem/dir_TMEM16F/All_8192_1024.star -D 256 --Apix 1.059 -o /home/rangan/cryodrgn/dir_TMEM16F/TMEM16F_cryodrgn_ctf.pkl ;
# generated TMEM16F_cryodrgn_ctf.pkl

# Now trying crodrgn abinit_homo as suggested at: ;
# https://cryodrgn.notion.site/CryoDRGN2-quickstart-322823599fce4bd7a391d00bf749ab1f ;
cryodrgn abinit_homo /home/rangan/cryodrgn/dir_TMEM16F/TMEM16F_cryodrgn_1024.mrcs --ctf /home/rangan/cryodrgn/dir_TMEM16F/TMEM16F_cryodrgn_ctf.pkl --outdir /home/rangan/cryodrgn/dir_TMEM16F/dir_job_0 >> cryodrgn_job_0.log ;

scp -p /home/rangan/dir_cryoem/dir_TMEM16F/cryodrgn_TMEM16F_0.sh rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_TMEM16F/cryodrgn_TMEM16F_0.sh ;
