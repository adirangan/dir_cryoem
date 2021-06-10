addpath('/data/rangan/dir_cryoem/chebfun');
addpath('/data/rangan/dir_bcc/dir_lakcluster_c_dev/dir_m');
addpath('/data/rangan/dir_bcc/dir_PGC_20180304/dir_m');
addpath('/data/rangan/dir_cryoem/dir_nufftall-1.33/');
addpath('/data/rangan/dir_cryoem/dir_rangan_playpen/');
addpath('/data/rangan/dir_cryoem/dir_rangan_playroom/');
addpath('/data/rangan/dir_cryoem/dir_rangan_playhouse/');
addpath('/data/rangan/dir_cryoem/dir_rangan_playhouse/dir_freq_march/');
addpath('/data/rangan/dir_cryoem/dir_rangan_playground/');
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig));
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'));
twopi = 2*pi;
