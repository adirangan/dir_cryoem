# stored in: ;
cd /home/rangan/cryodrgn/dir_trpv1 ;
# conda activate ;
# conda activate cryodrgn ;

# This was suggested by https://github.com/zhonge/cryodrgn, ;
# but does not work (as the --relion31 flag is not recognized). ;
# cryodrgn downsample -D 128 -o /home/rangan/cryodrgn/dir_trpv1/tv1_cryodrgn_1024.mrcs --relion31 --datadir /home/rangan/dir_cryoem/dir_trpv1/ /home/rangan/dir_cryoem/dir_trpv1/dir_relion_mat/tv1_relion_job_1024.star ;
# trying this instead: ;
cryodrgn downsample -D 128 -o /home/rangan/cryodrgn/dir_trpv1/tv1_cryodrgn_1024.mrcs --datadir /home/rangan/dir_cryoem/dir_trpv1/ /home/rangan/dir_cryoem/dir_trpv1/dir_relion_mat/tv1_relion_job_1024.star ;
# generated tv1_cryodrgn_1024.mrcs (size 65M). ;

# Now trying to create ctf.pkl file. ;
cryodrgn parse_ctf_star /home/rangan/dir_cryoem/dir_trpv1/dir_relion_mat/tv1_relion_job_1024.star -D 256 --Apix 1.2156 -o /home/rangan/cryodrgn/dir_trpv1/tv1_cryodrgn_ctf.pkl ;
# This did not work, complaining about a _rlnPhaseShift variable. ;
# Trying again, this time with original star file. ;
cryodrgn parse_ctf_star /home/rangan/dir_cryoem/dir_trpv1/tv1_relion_data.star -D 256 --Apix 1.2156 -o /home/rangan/cryodrgn/dir_trpv1/tv1_cryodrgn_ctf.pkl ;
# Trying again, this time adding phase-shift (manually). ;
cryodrgn parse_ctf_star /home/rangan/dir_cryoem/dir_trpv1/dir_relion_mat/tv1_relion_job_1024.star --ps 0 -D 256 --Apix 1.2156 -o /home/rangan/cryodrgn/dir_trpv1/tv1_cryodrgn_ctf.pkl ;
# generated tv1_cryodrgn_ctf.pkl

# Now trying crodrgn abinit_homo as suggested at: ;
# https://cryodrgn.notion.site/CryoDRGN2-quickstart-322823599fce4bd7a391d00bf749ab1f ;
cryodrgn abinit_homo /home/rangan/cryodrgn/dir_trpv1/tv1_cryodrgn_1024.mrcs --ctf /home/rangan/cryodrgn/dir_trpv1/tv1_cryodrgn_ctf.pkl --outdir /home/rangan/cryodrgn/dir_trpv1/dir_job_0 >> cryodrgn_job_0.log ;

scp -p /home/rangan/dir_cryoem/dir_trpv1/cryodrgn_trpv1_0.sh rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_trpv1/cryodrgn_trpv1_0.sh ;
