# stored in: ;
cd /home/rangan/cryodrgn/dir_rib80s ;
# conda activate ;
# conda activate cryodrgn ;

# copy over the appropriate mrcs files: ;
<<comment_download
scp -p rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_rib80s/shiny_2sets.star /home/rangan/dir_cryoem/dir_rib80s/ ;
head -1044 /home/rangan/dir_cryoem/dir_rib80s/shiny_2sets.star > /home/rangan/dir_cryoem/dir_rib80s/shiny_2sets_1024.star ;
mkdir /home/rangan/dir_cryoem/dir_rib80s/Particles ;
mkdir /home/rangan/dir_cryoem/dir_rib80s/Particles/MRC_1901 ;
scp -p rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_rib80s/Particles/MRC_1901/001_particles_shiny_nb50_new.mrcs /home/rangan/dir_cryoem/dir_rib80s/Particles/MRC_1901/ ;
scp -p rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_rib80s/Particles/MRC_1901/002_particles_shiny_nb50_new.mrcs /home/rangan/dir_cryoem/dir_rib80s/Particles/MRC_1901/ ;
scp -p rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_rib80s/Particles/MRC_1901/003_particles_shiny_nb50_new.mrcs /home/rangan/dir_cryoem/dir_rib80s/Particles/MRC_1901/ ;
scp -p rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_rib80s/Particles/MRC_1901/004_particles_shiny_nb50_new.mrcs /home/rangan/dir_cryoem/dir_rib80s/Particles/MRC_1901/ ;
scp -p rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_rib80s/Particles/MRC_1901/005_particles_shiny_nb50_new.mrcs /home/rangan/dir_cryoem/dir_rib80s/Particles/MRC_1901/ ;
scp -p rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_rib80s/Particles/MRC_1901/006_particles_shiny_nb50_new.mrcs /home/rangan/dir_cryoem/dir_rib80s/Particles/MRC_1901/ ;
scp -p rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_rib80s/Particles/MRC_1901/007_particles_shiny_nb50_new.mrcs /home/rangan/dir_cryoem/dir_rib80s/Particles/MRC_1901/ ;
scp -p rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_rib80s/Particles/MRC_1901/008_particles_shiny_nb50_new.mrcs /home/rangan/dir_cryoem/dir_rib80s/Particles/MRC_1901/ ;
scp -p rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_rib80s/Particles/MRC_1901/009_particles_shiny_nb50_new.mrcs /home/rangan/dir_cryoem/dir_rib80s/Particles/MRC_1901/ ;
scp -p rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_rib80s/Particles/MRC_1901/010_particles_shiny_nb50_new.mrcs /home/rangan/dir_cryoem/dir_rib80s/Particles/MRC_1901/ ;
comment_download
# This was suggested by https://github.com/zhonge/cryodrgn, ;
# Modified as the --relion31 flag is not recognized: ;
cryodrgn downsample -D 128 -o /home/rangan/cryodrgn/dir_rib80s/rib80s_cryodrgn_1024.mrcs --datadir /home/rangan/dir_cryoem/dir_rib80s/ /home/rangan/dir_cryoem/dir_rib80s/shiny_2sets_1024.star ;
# generated rib80s_cryodrgn_1024.mrcs (size 65M). ;

# Now trying to create ctf.pkl file. ;
# Adding phase-shift (manually). ;
cryodrgn parse_ctf_star /home/rangan/dir_cryoem/dir_rib80s/shiny_2sets_1024.star --ps 0 -D 360 --Apix 1.34 -o /home/rangan/cryodrgn/dir_rib80s/rib80s_cryodrgn_ctf.pkl ;
# generated rib80s_cryodrgn_ctf.pkl

# Now trying crodrgn abinit_homo as suggested at: ;
# https://cryodrgn.notion.site/CryoDRGN2-quickstart-322823599fce4bd7a391d00bf749ab1f ;
cryodrgn abinit_homo /home/rangan/cryodrgn/dir_rib80s/rib80s_cryodrgn_1024.mrcs --ctf /home/rangan/cryodrgn/dir_rib80s/rib80s_cryodrgn_ctf.pkl --outdir /home/rangan/cryodrgn/dir_rib80s/dir_job_0 >> cryodrgn_job_0.log ;

scp -p /home/rangan/dir_cryoem/dir_rib80s/cryodrgn_rib80s_0.sh rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_rib80s/cryodrgn_rib80s_0.sh ;
