addpath('/data/rangan/dir_cryoem/chebfun');
addpath('/data/rangan/dir_cryoem/dir_rangan_playhouse');
addpath('/data/rangan/dir_bcc/dir_lakcluster_c_dev/dir_m');
addpath('/data/rangan/dir_bcc/dir_PGC_20180304/dir_m');
addpath('/data/rangan/dir_cryoem/nufftall-1.33/');
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig));
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'));
