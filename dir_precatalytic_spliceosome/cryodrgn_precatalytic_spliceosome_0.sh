# stored in: ;
mkdir /home/rangan/cryodrgn/dir_ps1 ;
cd /home/rangan/cryodrgn/dir_ps1 ;
# conda activate ;
# conda activate cryodrgn ;

# copy over the appropriate mrcs files: ;
<<comment_download
%%%%%%%%;
head -1065 /home/rangan/dir_cryoem/dir_precatalytic_spliceosome/consensus_data.star > /home/rangan/dir_cryoem/dir_precatalytic_spliceosome/consensus_data_1024.star ;
%%%%%%%%;
comment_download
# This was suggested by https://github.com/zhonge/cryodrgn, ;
# Modified as the --relion31 flag is not recognized: ;
cryodrgn downsample -D 128 -o /home/rangan/cryodrgn/dir_ps1/ps1_cryodrgn_1024.mrcs --datadir /home/rangan/dir_cryoem/dir_precatalytic_spliceosome/ /home/rangan/dir_cryoem/dir_precatalytic_spliceosome/consensus_data_1024.star ;
# generated ps1_cryodrgn_1024.mrcs (size 65M). ;

# Now trying to create ctf.pkl file. ;
cryodrgn parse_ctf_star /home/rangan/dir_cryoem/dir_precatalytic_spliceosome/consensus_data_1024.star -D 320 --Apix 1.699 -o /home/rangan/cryodrgn/dir_ps1/ps1_cryodrgn_ctf.pkl ;
# generated ps1_cryodrgn_ctf.pkl

# Now trying crodrgn abinit_homo as suggested at: ;
# https://cryodrgn.notion.site/CryoDRGN2-quickstart-322823599fce4bd7a391d00bf749ab1f ;
cryodrgn abinit_homo /home/rangan/cryodrgn/dir_ps1/ps1_cryodrgn_1024.mrcs --ctf /home/rangan/cryodrgn/dir_ps1/ps1_cryodrgn_ctf.pkl --outdir /home/rangan/cryodrgn/dir_ps1/dir_job_0 >> cryodrgn_job_0.log ;

cryodrgn abinit_homo /home/rangan/cryodrgn/dir_ps1/ps1_cryodrgn_1024.mrcs --ctf /home/rangan/cryodrgn/dir_ps1/ps1_cryodrgn_ctf.pkl --seed 1 --outdir /home/rangan/cryodrgn/dir_ps1/dir_job_1 >> cryodrgn_job_1.log ;

cryodrgn abinit_homo /home/rangan/cryodrgn/dir_ps1/ps1_cryodrgn_1024.mrcs --ctf /home/rangan/cryodrgn/dir_ps1/ps1_cryodrgn_ctf.pkl --seed 2 --outdir /home/rangan/cryodrgn/dir_ps1/dir_job_2 >> cryodrgn_job_2.log ;

cryodrgn abinit_homo /home/rangan/cryodrgn/dir_ps1/ps1_cryodrgn_1024.mrcs --ctf /home/rangan/cryodrgn/dir_ps1/ps1_cryodrgn_ctf.pkl --seed 3 --outdir /home/rangan/cryodrgn/dir_ps1/dir_job_3 >> cryodrgn_job_3.log ;


scp -p /home/rangan/dir_cryoem/dir_precatalytic_spliceosome/cryodrgn_precatalytic_spliceosome_0.sh rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_precatalytic_spliceosome/cryodrgn_precatalytic_spliceosome_0.sh ;
