addpath('/mnt/home/rangan/ceph/dir_cryoem/chebfun');
addpath('/mnt/home/rangan/ceph/dir_bcc/dir_lakcluster_c_dev/dir_m');
addpath('/mnt/home/rangan/ceph/dir_bcc/dir_PGC_20180304/dir_m');
addpath('/mnt/home/rangan/ceph/dir_cryoem/dir_EMIODist2');
addpath('/mnt/home/rangan/ceph/dir_cryoem/dir_nufftall-1.33/');
addpath('/mnt/home/rangan/ceph/dir_cryoem/dir_rangan_playpen/');
addpath('/mnt/home/rangan/ceph/dir_cryoem/dir_rangan_playroom/');
addpath('/mnt/home/rangan/ceph/dir_cryoem/dir_rangan_playhouse/');
addpath('/mnt/home/rangan/ceph/dir_cryoem/dir_rangan_playhouse/dir_freq_march/');
addpath('/mnt/home/rangan/ceph/dir_cryoem/dir_rangan_playground/');
addpath('/mnt/home/rangan/ceph/dir_bcc/dir_ukb/dir_m/');
addpath('/mnt/home/rangan/ceph/dir_bcc/dir_jamison/dir_m/');
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig));
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'));
twopi = 2*pi;