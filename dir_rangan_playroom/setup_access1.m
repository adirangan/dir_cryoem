addpath('/data/rangan/dir_cryoem/chebfun');
addpath('/data/rangan/dir_bcc/dir_lakcluster_c_dev/dir_m');
addpath('/data/rangan/dir_bcc/dir_PGC_20180304/dir_m');
addpath('/data/rangan/dir_cryoem/dir_EMIODist2');
addpath('/data/rangan/dir_cryoem/dir_nufftall-1.33/');
addpath('/data/rangan/dir_cryoem/dir_rangan_playpen/');
addpath('/data/rangan/dir_cryoem/dir_rangan_playroom/');
addpath('/data/rangan/dir_cryoem/dir_rangan_playroom/dir_eig_ddssnll_lanczos_local/');
addpath('/data/rangan/dir_cryoem/dir_rangan_playhouse/');
addpath('/data/rangan/dir_cryoem/dir_rangan_playhouse/dir_freq_march/');
addpath('/data/rangan/dir_cryoem/dir_rangan_playground/');
addpath('/data/rangan/dir_bcc/dir_ukb/dir_m/');
addpath('/data/rangan/dir_bcc/dir_jamison/dir_m/');
addpath('/data/rangan/dir_bcc/dir_dolphin/dir_m/');
addpath('/home/rangan/dir_finufft/matlab/');
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig));
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'));
twopi = 2*pi;
setenv('NVIDIA_CUDNN','/usr/local/cuda/');
