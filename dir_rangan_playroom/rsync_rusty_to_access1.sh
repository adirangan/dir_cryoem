#scp -p /mnt/home/rangan/dir_cryoem/dir_rangan_playroom/rib80s_bin_to_mda_0.f /mnt/home/rangan/dir_cryoem/dir_rangan_playroom/rsync_rusty_to_access1.sh rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_rib80s/new_joakim_relion_run/subs128/data/ ;
#scp -p /mnt/home/rangan/dir_cryoem/dir_rib80s/new_joakim_relion_run/subs128/data/images_mda rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_rib80s/new_joakim_relion_run/subs128/data/ ;
rsync -anvum ~/dir_cryoem/dir_rangan_playroom/ rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_rangan_playroom/ ;

