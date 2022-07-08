function ...
[ ...
 global_parameter ...
] = ...
test_pm_24_combine_1( ...
 global_parameter ...
,fname_prefix ...
,dir_nopath_data_star ...
,Pixel_Spacing ...
,fname_nopath_volume ...
,fname_nopath_star ...
);
%%%%%%%%;
% instead of running each octile (as in test_pm_24_combine_0), ;
% we discard the bottom octiles (starting from the worst). ;
%%%%%%%%;

% try: ;
if (nargin<1);
nf=0;
nx=0;
global_parameter=struct('type','parameter');
global_parameter.flag_center_image = 1;
global_parameter.nM_start = 1024*nx;
global_parameter.nM_final = global_parameter.nM_start + 1024*8 - 1;
global_parameter.n_x_u_pack_0in = 64;
global_parameter.k_p_r_max_0in = 48/(2*pi);
global_parameter.delta_sigma_use_0in = 0.05;
global_parameter.dat_rseed_0in_ = [0:2];
global_parameter.rank_pm_0in_ = [16,18];
global_parameter.tolerance_pm_0in_ = 0.1.^[2:0.5:3];
fname_prefix=sprintf('ps1_x%dto%d_combine',floor(global_parameter.nM_start/1024),floor(global_parameter.nM_final/1024));
dir_nopath_data_star='precatalytic_spliceosome';
Pixel_Spacing=1.699;
fname_nopath_volume='consensus_half1_class001.mrc';
fname_nopath_star='consensus_data.star';
test_pm_24_combine_1( ...
	    global_parameter ...
	    ,fname_prefix ...
	    ,dir_nopath_data_star ...
	    ,Pixel_Spacing ...
	    ,fname_nopath_volume ...
	    ,fname_nopath_star ...
	    );
disp('returning');return;
end;%if (nargin<1);

%%%%%%%%;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
%%%%%%%%;

if isempty(global_parameter); global_parameter = struct('type','parameter'); end;
if (~isfield(global_parameter,'flag_verbose')); global_parameter.flag_verbose = 0; end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'flag_recalc')); global_parameter.flag_recalc = 0; end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'flag_replot')); global_parameter.flag_replot = 0; end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'flag_center_volume')); global_parameter.flag_center_volume = 0; end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'flag_center_image')); global_parameter.flag_center_image = 0; end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'flag_invert')); global_parameter.flag_invert = 0; end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'tolerance_master')); global_parameter.tolerance_master = 1e-2; end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'nM_start')); global_parameter.nM_start = 0; end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'nM_final')); global_parameter.nM_final = 1024; end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'n_x_u_pack_0in')); global_parameter.n_x_u_pack_0in = 64; end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'k_p_r_max_0in')); global_parameter.k_p_r_max_0in = 48/(2*pi); end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'date_diff_threshold_0in')); global_parameter.date_diff_threshold_0in = 0.25; end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'flag_force_create_mat_0in')); global_parameter.flag_force_create_mat_0in = 0; end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'flag_force_create_tmp_0in')); global_parameter.flag_force_create_tmp_0in = 0; end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'flag_alternate_MS_vs_SM_0in_')); global_parameter.flag_alternate_MS_vs_SM_0in_ = [1]; end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'delta_sigma_use_0in')); global_parameter.delta_sigma_use_0in = []; end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'dat_rseed_0in_')); global_parameter.dat_rseed_0in_ = [0]; end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'tolerance_pm_0in_')); global_parameter.tolerance_pm_0in_ = []; end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'rank_pm_0in_')); global_parameter.rank_pm_0in_ = [16]; end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'delta_r_max_factor_0in_')); global_parameter.delta_r_max_factor_0in_ = [1.00]; end; %<-- parameter_bookmark. ;
flag_verbose = global_parameter.flag_verbose;
flag_recalc = global_parameter.flag_recalc;
flag_replot = global_parameter.flag_replot;
flag_center_volume = global_parameter.flag_center_volume;
flag_center_image = global_parameter.flag_center_image;
flag_invert = global_parameter.flag_invert;
tolerance_master = global_parameter.tolerance_master;
nM_start = global_parameter.nM_start;
nM_final = global_parameter.nM_final;
n_x_u_pack_0in = global_parameter.n_x_u_pack_0in;
k_p_r_max_0in = global_parameter.k_p_r_max_0in;
date_diff_threshold_0in = global_parameter.date_diff_threshold_0in;
flag_force_create_mat_0in = global_parameter.flag_force_create_mat_0in;
flag_force_create_tmp_0in = global_parameter.flag_force_create_tmp_0in;
flag_alternate_MS_vs_SM_0in_ = global_parameter.flag_alternate_MS_vs_SM_0in_;
delta_sigma_use_0in = global_parameter.delta_sigma_use_0in;
dat_rseed_0in_ = global_parameter.dat_rseed_0in_;
tolerance_pm_0in_ = global_parameter.tolerance_pm_0in_;
rank_pm_0in_ = global_parameter.rank_pm_0in_;
delta_r_max_factor_0in_ = global_parameter.delta_r_max_factor_0in_;

nf=0;
if (flag_verbose); disp(sprintf(' %% [entering test_pm_24_combine_1]')); end;

%%%%%%%%;
fname_prefix_xfix = sprintf('%s',fname_prefix);
dir_pm = sprintf('/%s/rangan/dir_cryoem/dir_%s/dir_pm',string_root,fname_prefix_xfix);
if (~exist(sprintf('%s_mat',dir_pm),'dir')); disp(sprintf(' %% mkdir %s_mat',dir_pm)); mkdir(sprintf('%s_mat',dir_pm)); end;
if (~exist(sprintf('%s_jpg',dir_pm),'dir')); disp(sprintf(' %% mkdir %s_jpg',dir_pm)); mkdir(sprintf('%s_jpg',dir_pm)); end;
dir_relion = sprintf('/%s/rangan/dir_cryoem/dir_%s/dir_relion',string_root,fname_prefix_xfix);
if (~exist(sprintf('%s_mat',dir_relion),'dir')); disp(sprintf(' %% mkdir %s_mat',dir_relion)); mkdir(sprintf('%s_mat',dir_relion)); end;
if (~exist(sprintf('%s_jpg',dir_relion),'dir')); disp(sprintf(' %% mkdir %s_jpg',dir_relion)); mkdir(sprintf('%s_jpg',dir_relion)); end;
string_rusty_root = 'mnt/home';
dir_relion_bin = sprintf('/%s/rangan/relion/build/bin',string_rusty_root);
dir_data_star = sprintf('/%s/rangan/dir_cryoem/dir_%s',string_root,dir_nopath_data_star);
%%%%%%%%;
% all classes and subclasses. ;
%%%%%%%%;
fname_nopath_volume_ = ...
{ ...
 fname_nopath_volume ... %<-- class 0. ;
};
n_volume = numel(fname_nopath_volume_);
flag_het = 0; if (n_volume> 1); flag_het = 1; end;

%%%%%%%%;
% loading the individual runs. ;
%%%%%%%%;
narm_min = floor(nM_start/1024); %<-- assumes that narm_min==0;
if (narm_min~=0); disp(sprintf(' %% Warning, narm_min %d ~= 0 in test_pm_24_combine_1',narm_min)); end;
narm_max = floor(nM_final/1024);
narm_ = narm_min:narm_max; n_arm = numel(narm_);
fname_prefix_narm_ = cell(n_arm,1);
fname_prefix_sub = fname_prefix(1:strfind(fname_prefix,'_x0to')-1);
dir_pm_narm_ = cell(n_arm,1);
dir_pm_mat_narm_ = cell(n_arm,1);
for narm=0:n_arm-1;
fname_prefix_narm_{1+narm} = sprintf('%s_x%d',fname_prefix_sub,narm);
dir_pm_narm_{1+narm} = sprintf('/%s/rangan/dir_cryoem/dir_%s/dir_pm',string_root,fname_prefix_narm_{1+narm});
dir_pm_mat_narm_{1+narm} = sprintf('/%s/rangan/dir_cryoem/dir_%s/dir_pm_mat',string_root,fname_prefix_narm_{1+narm});
if (~exist(dir_pm_mat_narm_{1+narm},'dir')); disp(sprintf(' %% Warning, %s not found in test_pm_24_combine_1',dir_pm_mat_narm_{1+narm})); end;
end;%for narm=0:n_arm-1;

str_filter = 'X_2d_xcor_d0_a1t0122p*';
ampm__ = cell(n_arm,1);
for narm=0:n_arm-1;
if (flag_verbose); disp(sprintf(' %% narm %d/%d',narm,n_arm)); end;
global_parameter.flag_store_S_k_p__ = 1;
global_parameter.flag_store_M_k_p__ = 1;
[~,ampm__{1+narm}] ...
= ...
ampm_fromdisk_2( ...
 global_parameter ...
,dir_pm_narm_{1+narm} ...
,str_filter ...
);
end;%for narm=0:n_arm-1;

%%%%%%%%;
% Now just run a very simple version of ampmut_wrap_wrap_5 after segregating the images by octiles based on estimated quality. ;
% We will run ampmut_wrap_wrap_5 after discarding images in the bottom octiles. ;
%%%%%%%%;
n_w_max = ampm__{1+0}.n_w_max;
n_w_sum = ampm__{1+0}.n_w_sum;
n_k_p_r = ampm__{1+0}.n_k_p_r;
k_p_r_ = ampm__{1+0}.k_p_r_;
k_p_r_max = ampm__{1+0}.k_p_r_max;
weight_3d_k_p_r_ = ampm__{1+0}.weight_3d_k_p_r_;
weight_2d_k_p_r_ = ampm__{1+0}.weight_2d_k_p_r_;
n_w_max = ampm__{1+0}.n_w_max;
l_max_ = ampm__{1+0}.l_max_;
a_k_Y_true_ = ampm__{1+0}.a_k_Y_quad_;
%%%%%%%%;
n_M_narm_ = zeros(n_arm,1);
for narm=0:n_arm-1;
n_M_narm_(1+narm) = ampm__{1+narm}.n_M;
end;%for narm=0:n_arm-1;
n_M_csum_narm_ = cumsum([0;n_M_narm_]);
n_M_sum = sum(n_M_narm_);
N_k_p_wkMsum__ = zeros(n_w_sum,n_M_sum);
for narm=0:n_arm-1;
tmp_index_ = n_M_csum_narm_(1+narm) + [0:n_M_narm_(1+narm)-1];
N_k_p_wkMsum__(:,1+tmp_index_) = ampm__{1+narm}.N_k_p__;
end;%for narm=0:n_arm-1;
%%%%%%%%;
n_CTF_narm_ = zeros(n_arm,1);
for narm=0:n_arm-1;
n_CTF_narm_(1+narm) = ampm__{1+narm}.n_CTF;
end;%for narm=0:n_arm-1;
n_CTF_csum_narm_ = cumsum([0;n_CTF_narm_]);
n_CTF_sum = sum(n_CTF_narm_);
CTF_k_p_r_kCsum__ = zeros(n_k_p_r,n_CTF_sum);
index_nCTFsum_from_nMsum_ = zeros(n_M_sum,1);
for narm=0:n_arm-1;
tmp_index_ = n_CTF_csum_narm_(1+narm) + [0:n_CTF_narm_(1+narm)-1];
CTF_k_p_r_kCsum__(:,1+tmp_index_) = ampm__{1+narm}.CTF_k_p_r_kC__;
tmp_index_ = n_M_csum_narm_(1+narm) + [0:n_M_narm_(1+narm)-1];
index_nCTFsum_from_nMsum_(1+tmp_index_) = ampm__{1+narm}.index_nCTF_from_nM_ + n_CTF_csum_narm_(1+narm);
end;%for narm=0:n_arm-1;
%%%%%%%%;
image_X_value_est_Msum_ = zeros(n_M_sum,1);
image_R_value_est_Msum_ = zeros(n_M_sum,1);
image_X_value_tru_Msum_ = zeros(n_M_sum,1);
euler_polar_a_tru_Msum_ = zeros(n_M_sum,1);
euler_azimu_b_tru_Msum_ = zeros(n_M_sum,1);
euler_gamma_z_tru_Msum_ = zeros(n_M_sum,1);
image_delta_x_tru_Msum_ = zeros(n_M_sum,1);
image_delta_y_tru_Msum_ = zeros(n_M_sum,1);
for narm=0:n_arm-1;
%%%%%%%%;
disp(sprintf(' %% narm %d/%d',narm,n_arm));
ampm_ = ampm__{1+narm};
tmp_index_ = n_M_csum_narm_(1+narm) + [0:n_M_narm_(1+narm)-1];
if isfield(ampm_,'image_X_value_Mi__'); image_X_value_tru_Msum_(1+tmp_index_) = ampm_.image_X_value_Mi__(:,end-1); end;
tmp_f = @(s) ~isempty(strfind(s,'X_2d_xcor_d0_a1t0122p25r0.mat'));
tmp_index = efind(cellfun(tmp_f,ampm_.str_fname_mat_a_));
if numel(tmp_index)==1;
if (flag_verbose); disp(ampm_.str_fname_mat_a_{1+tmp_index}); end;
image_X_value_est_Msum_(1+tmp_index_) = ampm_.image_X_value_ampm_Ma__(:,1+tmp_index);
image_R_value_est_Msum_(1+tmp_index_) = ampm_.image_R_value_ampm_Ma__(:,1+tmp_index);
euler_polar_a_tru_Msum_(1+tmp_index_) = ampm_.euler_polar_a_Mi__(:,end);
euler_azimu_b_tru_Msum_(1+tmp_index_) = ampm_.euler_azimu_b_Mi__(:,end);
euler_gamma_z_tru_Msum_(1+tmp_index_) = ampm_.euler_gamma_z_Mi__(:,end);
image_delta_x_tru_Msum_(1+tmp_index_) = ampm_.image_delta_x_acc_Mi__(:,end) + ampm_.image_delta_x_upd_Mi__(:,end);
image_delta_y_tru_Msum_(1+tmp_index_) = ampm_.image_delta_y_acc_Mi__(:,end) + ampm_.image_delta_y_upd_Mi__(:,end);
end;%if numel(tmp_index)==1;
%%%%%%%%;
end;%for narm=0:n_arm-1;
%%%%%%%%;

%%%%%%%%;
n_M = ampm__{1+0}.n_M;
n_octile = 8;
[~,ij_image_X_value_est_Msum_from_sort_] = sort(image_X_value_est_Msum_,'ascend');
[~,ij_image_X_value_est_sort_from_Msum_] = sort(ij_image_X_value_est_Msum_from_sort_,'ascend');
ij_from_each_noctile__ = cell(n_octile,1);
ij_from_excl_noctile__ = cell(n_octile,1);
for noctile=0:n_octile-1;
ij_from_each_noctile__{1+noctile} = ij_image_X_value_est_Msum_from_sort_(1+[0:n_M-1] + floor((n_M_sum-n_M)*noctile/(n_octile-1)));
ij_from_excl_noctile__{1+noctile} = ij_image_X_value_est_Msum_from_sort_(1+[noctile*n_M:n_M_sum-1]);
end;%for noctile=0:n_octile-1;
%%%%%%%%;
flag_disp=0;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;
c_80s__ = colormap_80s; n_c_80s = size(c_80s__,1);
subplot(1,2,1);
hold on;
for noctile=0:n_octile-1;
ij_from_each_noctile_ = ij_from_each_noctile__{1+noctile};
nc_80s = max(0,min(n_c_80s-1,floor(n_c_80s*noctile/n_octile)));
plot(noctile*n_M+[0:n_M-1],image_X_value_est_Msum_(ij_from_each_noctile_),'.','Color',c_80s__(1+nc_80s,:));
end;%for noctile=0:n_octile-1;
hold off;
xlim([0,n_M_sum]);
xlabel('nMsum','Interpreter','none');
subplot(1,2,2);
hold on;
for noctile=0:n_octile-1;
ij_from_excl_noctile_ = ij_from_excl_noctile__{1+noctile};
nc_80s = max(0,min(n_c_80s-1,floor(n_c_80s*noctile/n_octile)));
plot(1+[noctile*n_M:n_M_sum-1],noctile*1e-2 + image_X_value_est_Msum_(ij_from_excl_noctile_),'.','Color',c_80s__(1+nc_80s,:));
end;%for noctile=0:n_octile-1;
hold off;
xlim([0,n_M_sum]);
xlabel('nMsum','Interpreter','none');
end;%if flag_disp;
%%%%%%%%;

n_CTF_noctile_ = zeros(n_octile,1);
CTF_k_p_r_kC_per_noctile___ = cell(n_octile,1);
for noctile=0:n_octile-1;
ij_from_excl_noctile_ = ij_from_excl_noctile__{1+noctile};
index_nCTFsum_from_nMsum_this_noctile_ = index_nCTFsum_from_nMsum_(ij_from_excl_noctile_);
[u_index_nCTFsum_from_nMsum_this_noctile_,ij_na_from_nu_,ij_nu_from_na_] = unique(index_nCTFsum_from_nMsum_this_noctile_);
CTF_k_p_r_kC_per_noctile___{1+noctile} = CTF_k_p_r_kCsum__(:,1+u_index_nCTFsum_from_nMsum_this_noctile_);
index_nCTF_from_nM_per_noctile__{1+noctile} = ij_nu_from_na_-1;
n_CTF_noctile_(1+noctile) = numel(u_index_nCTFsum_from_nMsum_this_noctile_);
end;%for noctile=0:n_octile-1;
tmp_error=0; n_val=0;
for noctile=0:n_octile-1;
ij_from_excl_noctile_ = ij_from_excl_noctile__{1+noctile};
assert(size(CTF_k_p_r_kC_per_noctile___{1+noctile},2)==n_CTF_noctile_(1+noctile));
for nl=0:numel(ij_from_excl_noctile_)-1;
nMsum = ij_from_excl_noctile_(1+nl)-1;
tmp_CTF_0_ = CTF_k_p_r_kCsum__(:,1+index_nCTFsum_from_nMsum_(1+nMsum));
tmp_CTF_1_ = CTF_k_p_r_kC_per_noctile___{1+noctile}(:,1+index_nCTF_from_nM_per_noctile__{1+noctile}(1+nl));
tmp_error = tmp_error + fnorm(tmp_CTF_0_ - tmp_CTF_1_);
n_val = n_val+1;
end;%for nl=0:numel(ij_from_excl_noctile_)-1;
end;%for noctile=0:n_octile-1;
if (flag_verbose); disp(sprintf(' %% CTF error: %0.16f (n_val %d)',tmp_error,n_val)); end;

date_diff_threshold = date_diff_threshold_0in; %<-- default: 0.25;
flag_force_create_mat = flag_force_create_mat_0in; %<-- default: 0;
flag_force_create_tmp = flag_force_create_tmp_0in; %<-- default: 0;
if ~isempty(flag_alternate_MS_vs_SM_0in_); flag_alternate_MS_vs_SM_ = flag_alternate_MS_vs_SM_0in_; else; flag_alternate_MS_vs_SM_ = [1]; end; %<-- default: [1];
if ~isempty(delta_sigma_use_0in); delta_sigma_use = delta_sigma_use_0in; else; delta_sigma_use = delta_sigma; end; %<-- default: delta_sigma;
dat_rseed_ = dat_rseed_0in_; %<-- default: [0];
n_dat_rseed = numel(dat_rseed_);
if ~isempty(tolerance_pm_0in_); tolerance_pm_ = tolerance_pm_0in_; else; tolerance_pm_ = []; end; %<-- default: [];
n_tolerance_pm = numel(tolerance_pm_);
if ~isempty(rank_pm_0in_); rank_pm_ = rank_pm_0in_; else; rank_pm_ = [16]; end; %<-- default: [16];
n_rank_pm = numel(rank_pm_);
if ~isempty(delta_r_max_factor_0in_); delta_r_max_factor_ = delta_r_max_factor_0in_; else; delta_r_max_factor_ = [1.00]; end; %<-- default: 1;
n_delta_r_max_factor = numel(delta_r_max_factor_);
%%%%%%%%;
for ndat_rseed=0:n_dat_rseed-1;
dat_rseed = dat_rseed_(1+ndat_rseed);
for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
delta_r_max_upb = 2.0 * delta_sigma_use*sqrt(log(20^2)); %<-- allow large accumulated translations. ;
if  isempty(delta_sigma_use_0in);
if ( exist(sprintf('%s_mat/delta_r_max_legacy.txt',dir_pm)));
delta_r_max_legacy = textread(sprintf('%s_mat/delta_r_max_legacy.txt',dir_pm));
disp(sprintf(' %% loading delta_r_max_legacy: %0.16f',delta_r_max_legacy));
delta_r_max_use = delta_r_max_factor * delta_r_max_legacy;
delta_r_max_upb = 2.0 * delta_r_max_legacy; %<-- allow large accumulated translations. ;
end;%if ( exist(sprintf('%s_mat/delta_r_max_legacy.txt',dir_pm)));
end;%if  isempty(delta_sigma_use_0in);
for flag_alternate_MS_vs_SM = flag_alternate_MS_vs_SM_;
%%%%;
for ntolerance_pm=0:n_tolerance_pm-1;
tolerance_pm = tolerance_pm_(1+ntolerance_pm);
for noctile=0:n_octile-1;
ij_from_excl_noctile_ = ij_from_excl_noctile__{1+noctile};
parameter = struct('type','parameter');
parameter.rseed = dat_rseed;
parameter.flag_rank_vs_tolerance = 0;
parameter.tolerance_pm = tolerance_pm;
parameter.delta_r_max = delta_r_max_use;
parameter.delta_r_upb = delta_r_max_upb;
parameter.dir_pm = sprintf('%sox%d',dir_pm,noctile); %<-- octiles to exclude. ;
if ~exist(sprintf('%s_mat',parameter.dir_pm),'dir'); mkdir(sprintf('%s_mat',parameter.dir_pm)); end;
if ~exist(sprintf('%s_jpg',parameter.dir_pm),'dir'); mkdir(sprintf('%s_jpg',parameter.dir_pm)); end;
parameter.flag_alternate_MS_vs_SM = flag_alternate_MS_vs_SM;
parameter = ...
ampmut_wrap_wrap_5( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_max ...
,l_max_ ...
,numel(ij_from_excl_noctile_) ...
,N_k_p_wkMsum__(:,ij_from_excl_noctile_) ...
,n_CTF_noctile_(1+noctile) ...
,index_nCTF_from_nM_per_noctile__{1+noctile} ...
,CTF_k_p_r_kC_per_noctile___{1+noctile} ...
,a_k_Y_true_ ...
,euler_polar_a_tru_Msum_(ij_from_excl_noctile_) ...
,euler_azimu_b_tru_Msum_(ij_from_excl_noctile_) ...
,euler_gamma_z_tru_Msum_(ij_from_excl_noctile_) ...
,image_delta_x_tru_Msum_(ij_from_excl_noctile_) ...
,image_delta_y_tru_Msum_(ij_from_excl_noctile_) ...
);
end;%for noctile=0:n_octile-1;
end;%for ntolerance_pm=0:n_tolerance_pm-1;
%%%%;
for nrank_pm=0:n_rank_pm-1;
rank_pm = rank_pm_(1+nrank_pm);
for noctile=0:n_octile-1;
ij_from_excl_noctile_ = ij_from_excl_noctile__{1+noctile};
parameter = struct('type','parameter');
parameter.rseed = dat_rseed;
parameter.flag_rank_vs_tolerance = 1;
parameter.rank_pm = rank_pm;
parameter.delta_r_max = delta_r_max_use;
parameter.delta_r_upb = delta_r_max_upb;
parameter.dir_pm = sprintf('%sox%d',dir_pm,noctile); %<-- octiles to exclude. ;
if ~exist(sprintf('%s_mat',parameter.dir_pm),'dir'); mkdir(sprintf('%s_mat',parameter.dir_pm)); end;
if ~exist(sprintf('%s_jpg',parameter.dir_pm),'dir'); mkdir(sprintf('%s_jpg',parameter.dir_pm)); end;
parameter.flag_alternate_MS_vs_SM = flag_alternate_MS_vs_SM;
parameter = ...
ampmut_wrap_wrap_5( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_max ...
,l_max_ ...
,numel(ij_from_excl_noctile_) ...
,N_k_p_wkMsum__(:,ij_from_excl_noctile_) ...
,n_CTF_noctile_(1+noctile) ...
,index_nCTF_from_nM_per_noctile__{1+noctile} ...
,CTF_k_p_r_kC_per_noctile___{1+noctile} ...
,a_k_Y_true_ ...
,euler_polar_a_tru_Msum_(ij_from_excl_noctile_) ...
,euler_azimu_b_tru_Msum_(ij_from_excl_noctile_) ...
,euler_gamma_z_tru_Msum_(ij_from_excl_noctile_) ...
,image_delta_x_tru_Msum_(ij_from_excl_noctile_) ...
,image_delta_y_tru_Msum_(ij_from_excl_noctile_) ...
);
end;%for noctile=0:n_octile-1;
end;%for nrank_pm=0:n_rank_pm-1;
%%%%;
if (parameter.n_complete_calculation> 1);
if strcmp(platform,'rusty');
disp(sprintf(' %% parameter.n_complete_calculation %d, returning',parameter.n_complete_calculation));
return; %<-- halt after one calculation on rusty. ;
end;%if strcmp(platform,'rusty');
end;%if (parameter.n_complete_calculation> 1);
%%%%;
end;%for flag_alternate_MS_vs_SM = flag_alternate_MS_vs_SM_;
end;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
end;%for ndat_rseed=0:n_dat_rseed-1;
%%%%%%%%;



verbose=1;
if (flag_verbose); disp(sprintf(' %% [finished test_pm_24_combine_1]')); end;

disp('returning'); return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
