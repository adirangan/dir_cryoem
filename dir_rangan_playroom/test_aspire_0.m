function ...
[ ...
 global_parameter ...
] = ...
test_aspire_0( ...
 global_parameter ...
,fname_prefix ...
,dir_nopath_data_star ...
,Pixel_Spacing ...
,fname_nopath_volume ...
,fname_nopath_star ...
);

if (nargin<1);
table_data__ = { ...
%'p28hRPT1_x0' , 'p28hRP' , 0.98 , 'emd_8674.map' , 'T1.star' ; ...
%'LetB1_x0' , 'LetB1' , 1.31 , 'emd_20993.map' , 'job_569_model_1_350MB.star' ; ...
'ISWINCP_x0' , 'ISWINCP' , 1.07 , 'emd_9718.map' , 'ADPBeF.star' ; ...
'trpv1_x0' , 'trpv1' , 1.2156 , 'emd_5778.mrc' , 'tv1_relion_data.star' ; ...
'rib80s_x0' , 'rib80s' , 1.34 , 'emd_2660.mrc' , 'shiny_2sets.star' ; ...
'MlaFEDB_x0' , 'MlaFEDB' , 0.832 , 'emd_22116.map' , 'Empiar_10536_00_to_23.star' ; ...
'LetB1_x0' , 'LetB1' , 1.31 , 'emd_20993.map' , 'job_569_model_1.star' ; ...
'TMEM16F_x0' , 'TMEM16F' , 1.059 , 'emd_20244.map' , 'All_8192.star' ; ...
'LSUbl17dep_x0' , 'LSUbl17dep' , 1.31 , 'emd_8434.map' , 'Parameters_negated.star' ; ...
'ps1_x0' , 'precatalytic_spliceosome' , 1.699 , 'consensus_half1_class001.mrc' , 'consensus_data.star' ; ...
'LSUbl17depE_x0' , 'LSUbl17dep' , 1.31 , 'emd_8450.map' , 'Parameters_negated.star' ; ...
};
n_experiment = size(table_data__,1);
%%%%%%%%;
tmp_ = clock;rng(tmp_(end));
for nexperiment=(randperm(n_experiment)-1);
na=0;
fname_prefix = table_data__{1+nexperiment,1+na}; na=na+1;
dir_nopath_data_star = table_data__{1+nexperiment,1+na}; na=na+1;
Pixel_Spacing = table_data__{1+nexperiment,1+na}; na=na+1;
fname_nopath_volume = table_data__{1+nexperiment,1+na}; na=na+1;
fname_nopath_star = table_data__{1+nexperiment,1+na}; na=na+1;
disp(sprintf(' %% nexperiment %d/%d: %16s %16s %0.3f %16s %32s' ...
	     ,nexperiment,n_experiment ...
	     ,fname_prefix,dir_nopath_data_star,Pixel_Spacing,fname_nopath_volume,fname_nopath_star ...
	     ));
global_parameter = struct('type','parameter');
if (strcmp(dir_nopath_data_star,'LSUbl17dep')); global_parameter.flag_invert = 0; end;
if (strcmp(dir_nopath_data_star,'LSUbl17dep')); global_parameter.flag_center_image = 1; end;
if (strcmp(dir_nopath_data_star,'precatalytic_spliceosome')); global_parameter.flag_center_image = 1; end;
global_parameter = ...
test_aspire_0( ...
 global_parameter ...
,fname_prefix ...
,dir_nopath_data_star ...
,Pixel_Spacing ...
,fname_nopath_volume ...
,fname_nopath_star ...
);
end;%for nexperiment=0:n_experiment-1;
%%%%%%%%;
disp('returning');return;
end;%if (nargin<1);

% try: ;
% global_parameter=[];fname_prefix='LetB1_x0';dir_nopath_data_star='LetB1';Pixel_Spacing=1.31;fname_nopath_volume='emd_20993.map';fname_nopath_star='job_569_model_1_350MB.star';
% global_parameter=[];fname_prefix='trpv1_x0';dir_nopath_data_star='trpv1';Pixel_Spacing=1.2156;fname_nopath_volume='emd_5778.mrc';fname_nopath_star='tv1_relion_data.star';
% global_parameter=[];fname_prefix='rib80s_x0';dir_nopath_data_star='rib80s';Pixel_Spacing=1.34;fname_nopath_volume='emd_2660.mrc';fname_nopath_star='shiny_2sets.star';
% global_parameter=[];fname_prefix='LSUbl17dep_x0';dir_nopath_data_star='LSUbl17dep';Pixel_Spacing=1.31;fname_nopath_volume='emd_8434.map';fname_nopath_star='Parameters_negated.star' ;
% global_parameter=[];fname_prefix='LetB1_x0';dir_nopath_data_star='LetB1';Pixel_Spacing=1.31;fname_nopath_volume='emd_20993.map';fname_nopath_star='job_569_model_1.star' ;

str_thisfunction = 'test_aspire_0';

verbose=1;
if (verbose); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

%%%%%%%%;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
%%%%%%%%;

if isempty(global_parameter); global_parameter = struct('type','parameter'); end;
if (~isfield(global_parameter,'flag_recalc')); global_parameter.flag_recalc = 0; end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'flag_replot')); global_parameter.flag_replot = 0; end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'flag_center_volume')); global_parameter.flag_center_volume = 0; end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'flag_center_image')); global_parameter.flag_center_image = 0; end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'flag_invert')); global_parameter.flag_invert = 0; end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'tolerance_master')); global_parameter.tolerance_master = 1e-2; end; %<-- parameter_bookmark. ;
flag_recalc = global_parameter.flag_recalc;
flag_replot = global_parameter.flag_replot;
flag_center_volume = global_parameter.flag_center_volume;
flag_center_image = global_parameter.flag_center_image;
flag_invert = global_parameter.flag_invert;
tolerance_master = global_parameter.tolerance_master;
nf=0;

%%%%%%%%;
fname_prefix_xfix = sprintf('%s',fname_prefix);
dir_data_pm = sprintf('/%s/rangan/dir_cryoem/dir_%s/dir_pm','data',fname_prefix_xfix);
dir_pm = sprintf('/%s/rangan/dir_cryoem/dir_%s/dir_pm',string_root,fname_prefix_xfix);
if (~exist(sprintf('%s_mat',dir_pm),'dir')); disp(sprintf(' %% mkdir %s_mat',dir_pm)); mkdir(sprintf('%s_mat',dir_pm)); end;
if (~exist(sprintf('%s_jpg',dir_pm),'dir')); disp(sprintf(' %% mkdir %s_jpg',dir_pm)); mkdir(sprintf('%s_jpg',dir_pm)); end;
dir_relion = sprintf('/%s/rangan/dir_cryoem/dir_%s/dir_relion',string_root,fname_prefix_xfix);
if (~exist(sprintf('%s_mat',dir_relion),'dir')); disp(sprintf(' %% mkdir %s_mat',dir_relion)); mkdir(sprintf('%s_mat',dir_relion)); end;
if (~exist(sprintf('%s_jpg',dir_relion),'dir')); disp(sprintf(' %% mkdir %s_jpg',dir_relion)); mkdir(sprintf('%s_jpg',dir_relion)); end;
string_rusty_root = 'mnt/home';
dir_relion_bin = sprintf('/%s/rangan/relion/build/bin',string_rusty_root);
dir_data_star = sprintf('/%s/rangan/dir_cryoem/dir_%s',string_root,dir_nopath_data_star);
dir_cryoaspi = sprintf('/%s/rangan/dir_aspire/dir_%s',string_root,dir_nopath_data_star);
if ~exist(dir_cryoaspi,'dir'); disp(sprintf(' %% mkdir %s',dir_cryoaspi)); mkdir(dir_cryoaspi); end;
dir_aspi_pm = sprintf('%s/dir_pm',dir_cryoaspi);
if (~exist(sprintf('%s_mat',dir_aspi_pm),'dir')); disp(sprintf(' %% mkdir %s_mat',dir_aspi_pm)); mkdir(sprintf('%s_mat',dir_aspi_pm)); end;
if (~exist(sprintf('%s_jpg',dir_aspi_pm),'dir')); disp(sprintf(' %% mkdir %s_jpg',dir_aspi_pm)); mkdir(sprintf('%s_jpg',dir_aspi_pm)); end;
%%%%%%%%;
% all classes and subclasses. ;
%%%%%%%%;
fname_nopath_volume_ = ...
{ ...
 fname_nopath_volume ... %<-- class 0. ;
};
n_volume = numel(fname_nopath_volume_);
flag_het = 0; if (n_volume> 1); flag_het = 1; end;

fname_data_mat = sprintf('%s_mat/a_x_u_base_.mat',dir_data_pm);
fname_mat = sprintf('%s_mat/a_x_u_base_.mat',dir_pm);
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, downloading',fname_mat));
str_command = sprintf('scp -p rangan@access1.cims.nyu.edu:%s %s ; ',fname_data_mat,fname_mat); system(str_command);
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));

fname_data_mat = sprintf('%s_mat/a_k_p_quad_.mat',dir_data_pm);
fname_mat = sprintf('%s_mat/a_k_p_quad_.mat',dir_pm);
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, downloading',fname_mat));
str_command = sprintf('scp -p rangan@access1.cims.nyu.edu:%s %s ; ',fname_data_mat,fname_mat); system(str_command);
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));

fname_data_mat = sprintf('%s_mat/a_k_Y_quad_.mat',dir_data_pm);
fname_mat = sprintf('%s_mat/a_k_Y_quad_.mat',dir_pm);
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, downloading',fname_mat));
str_command = sprintf('scp -p rangan@access1.cims.nyu.edu:%s %s ; ',fname_data_mat,fname_mat); system(str_command);
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
if ~exist('l_max_','var');
l_max_upb = round(2*pi*k_p_r_max); %<-- typically sufficient for 2-3 digits of precision. ;
l_max_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
l_max_(1+nk_p_r) = max(0,min(l_max_upb,1+ceil(2*pi*k_p_r_(1+nk_p_r))));
end;%for nk_p_r=0:n_k_p_r-1;
n_lm_ = (l_max_+1).^2;
n_lm_max = max(n_lm_);
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
l_max_max = max(l_max_); dWtdkd__l_max_max = 2*l_max_max;
m_max_ = -l_max_max : +l_max_max;
n_m_max = length(m_max_);
Y_l_val_ = zeros(n_lm_sum,1);
Y_m_val_ = zeros(n_lm_sum,1);
Y_k_val_ = zeros(n_lm_sum,1);
for nk_p_r=0:n_k_p_r-1;
l_max = l_max_(1+nk_p_r);
tmp_l_val_ = zeros(n_lm_(1+nk_p_r),1);
tmp_m_val_ = zeros(n_lm_(1+nk_p_r),1);
na=0; 
for l_val=0:l_max;
for m_val=-l_val:+l_val;
tmp_l_val_(1+na) = l_val;
tmp_m_val_(1+na) = m_val;
na=na+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
Y_l_val_(1+tmp_index_) = tmp_l_val_;
Y_m_val_(1+tmp_index_) = tmp_m_val_;
Y_k_val_(1+tmp_index_) = k_p_r_(1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
weight_Y_ = zeros(n_lm_sum,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
weight_Y_(1+tmp_index_) = weight_3d_k_p_r_(1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
end;%if ~exist('l_max_','var');

fname_data_mat = sprintf('%s_mat/S_k_p__.mat',dir_data_pm);
fname_mat = sprintf('%s_mat/S_k_p__.mat',dir_pm);
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, downloading',fname_mat));
str_command = sprintf('scp -p rangan@access1.cims.nyu.edu:%s %s ; ',fname_data_mat,fname_mat); system(str_command);
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));

fname_data_mat = sprintf('%s_mat/M_k_p__.mat',dir_data_pm);
fname_mat = sprintf('%s_mat/M_k_p__.mat',dir_pm);
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, downloading',fname_mat));
str_command = sprintf('scp -p rangan@access1.cims.nyu.edu:%s %s ; ',fname_data_mat,fname_mat); system(str_command);
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));

%%%%%%%%;
% set aspire paths. ;
%%%%%%%%;
fname_aspi_mat = sprintf('%s_mat/abinitio_0.mat',dir_aspi_pm);
if ~exist(fname_aspi_mat,'file');
run(sprintf('/%s/rangan/dir_aspire/initpath',string_root));
%%%%%%%%;
% create M_x_c__. ;
%%%%%%%%;
M_x_c_pack_01M___ = zeros(n_x_u_pack,n_x_u_pack,n_M);
for nM=0:n_M-1;
M_x_c_ = ...
interp_k_p_to_x_c_xxnufft( ...
 n_x_u_pack ...
,diameter_x_c ...
,n_x_u_pack ...
,diameter_x_c ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,M_k_p__(:,1+nM) ...
);
M_x_c_pack_01M___(:,:,1+nM) = reshape(real(M_x_c_),[n_x_u_pack,n_x_u_pack]);
clear M_x_c_;
end;%for nM=0:n_M-1;
%%%%%%%%;
% set parameters. ;
%%%%%%%%;
tmp_L = n_x_u_pack;         % Size of images. ;
tmp_n = n_M;                % Number of images. ;
tmp_SNR = 0.2;              % Signal-to-noise ratio of images. ;
tmp_n_r = ceil(tmp_L/2);    % Number of radial nodes in polar Fourier transform. ;
tmp_n_theta = 72;           % Number of angular nodes in polar Fourier transform. ;
%%%%%%%%;
% Mask projections. ;
%%%%;
tmp_mask_radius = floor(n_x_u_pack/2); %<-- reasonable default. ;
[tmp_np,~]=mask_fuzzy(M_x_c_pack_01M___,tmp_mask_radius);
%%%%;
[tmp_npf,~]=cryo_pft(tmp_np,tmp_n_r,tmp_n_theta);
tmp_clstack_est = cryo_clmatrix(tmp_npf);
tmp_S_est = cryo_syncmatrix_vote(tmp_clstack_est,tmp_n_theta);
tmp_inv_rots_est = cryo_syncrotations(tmp_S_est);
% Set up parameters for volume estimation.
tmp_params = struct();
tmp_params.rot_matrices = tmp_inv_rots_est;        % Estimated rotations. ;
tmp_params.ctf = ones(tmp_L*ones(1, 2));           % CTFs (none here). ;
tmp_params.ctf_idx = ones(1,tmp_n);                % CTF indices (all the same). ;
tmp_params.ampl = ones(1,tmp_n);                   % Amplitude multipliers (all one). ;
tmp_params.shifts = zeros(2,tmp_n);                % Shifts (none here). ;
tmp_basis = ffb_basis(tmp_L*ones(1, 3));
% Set up options for the volume estimation algorithm.
mean_est_opt = struct();
mean_est_opt.verbose = false;               % Do not output progress info. ;
mean_est_opt.max_iter = 10;                 % Maximum number of iterations. ;
mean_est_opt.rel_tolerance = 1e-3;          % Stopping tolerance. ;
mean_est_opt.half_pixel = true;             % Center volumes around half pixel. ;
mean_est_opt.verbose = 1;                   % Print progress information. ;
% Estimate volume using least squares.
tmp_vol_est = cryo_estimate_mean(tmp_np,tmp_params,tmp_basis,mean_est_opt);
a_x_u_aspi_base_ = tmp_vol_est;
tmp_fname=sprintf('%s_mat/example1.mrc',dir_aspi_pm);
WriteMRC(tmp_vol_est,1,tmp_fname); % Output density map reconstructed from projections. ;
%%%%;
save(fname_aspi_mat ...
,'tmp_L','tmp_n','tmp_SNR','tmp_n_r','tmp_n_theta','tmp_mask_radius' ...
,'tmp_clstack_est','tmp_S_est','tmp_inv_rots_est','a_x_u_aspi_base_','tmp_fname' ...
);
%%%%%%%%;
end;%if ~exist(fname_aspi_mat,'file');
%%%%%%%%;
if  exist(fname_aspi_mat,'file');
load(fname_aspi_mat);
end;%if  exist(fname_aspi_mat,'file');
if ~exist('a_x_u_aspi_base_','var');
tmp_fname = sprintf('%s_mat/example1.mrc',dir_aspi_pm);
a_x_u_aspi_base_ = cast(ReadMRC(tmp_fname),'double');
end;%if ~exist('a_x_u_aspi_base_','var');

%%%%%%%%;
% Now convert to a_k_p_aspi_ ;
%%%%%%%%;
verbose=0;
fname_aspi_mat = sprintf('%s_mat/a_k_p_aspi_quad_.mat',dir_aspi_pm);
if (flag_recalc | ~exist(fname_aspi_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_aspi_mat));
eta = pi/x_p_r_max; tmp_t = tic;
a_k_p_aspi_quad_ = xxnufft3d3(n_xxx_u,x_u_0___(:)*eta,x_u_1___(:)*eta,x_u_2___(:)*eta,a_x_u_aspi_base_(:).*xxx_u_weight_(:),-1,1e-12,n_k_all,2*pi*k_c_0_all_/eta,2*pi*k_c_1_all_/eta,2*pi*k_c_2_all_/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_k_p_aspi_quad_ time %0.2fs',tmp_t));
%%%%%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
a_x_u_aspi_reco_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,a_k_p_aspi_quad_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_x_u_aspi_reco_ time %0.2fs',tmp_t));
disp(sprintf(' %% xxnufft3d3: a_x_u_aspi_reco error: %0.16f',fnorm(a_x_u_aspi_base_(:)-a_x_u_aspi_reco_)/fnorm(a_x_u_aspi_base_(:))));
save(fname_aspi_mat ...
     ,'k_p_r_max','k_eq_d' ...
     ,'n_k_all','n_k_all_csum_' ...
     ,'k_p_r_all_','k_p_azimu_b_all_','k_p_polar_a_all_' ...
     ,'weight_3d_k_all_','weight_shell_k_' ...
     ,'n_k_p_r','k_p_r_' ...
     ,'weight_3d_k_p_r_' ...
     ,'k_c_0_all_','k_c_1_all_','k_c_2_all_' ...
     ,'J_node_','J_weight_','J_chebfun_','J_polyval_' ...
     ,'a_k_p_aspi_quad_' ...
     ,'a_x_u_aspi_reco_' ...
     );
end;%if (~exist(fname_aspi_mat,'file'));
%%%%%%%%;
if ( exist(fname_aspi_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_aspi_mat));
load(fname_aspi_mat);
end;%if ( exist(fname_aspi_mat,'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/a_k_p_aspi_quad_',dir_aspi_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;
clf;
plot(k_p_r_all_,log10(abs(a_k_p_aspi_quad_)),'.'); xlabel('k'); ylabel('log10(|a(k)|)');
figbig;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

%%%%%%%%;
% Now convert to a_k_Y_aspi_ ; 
%%%%%%%%;
verbose=0;
fname_aspi_mat = sprintf('%s_mat/a_k_Y_aspi_quad_.mat',dir_aspi_pm);
if (flag_recalc | ~exist(fname_aspi_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_aspi_mat));
%%%%%%%%;
tmp_t = tic;
[a_k_Y_aspi_quad_] = convert_k_p_to_spharm_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_p_aspi_quad_);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_aspi_quad_ time %0.2fs',tmp_t));
tmp_t = tic;
[a_k_p_aspi_reco_] = convert_spharm_to_k_p_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_Y_aspi_quad_);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_aspi_quad_ --> a_k_p_aspi_reco_ time %0.2fs',tmp_t));
disp(sprintf(' %% xxnufft3d3: a_k_p_aspi_reco error: %0.16f',fnorm(a_k_p_aspi_quad_-a_k_p_aspi_reco_)/fnorm(a_k_p_aspi_quad_))); %<-- this should be 2-3 digits. ;
%%%%%%%%;
a_k_Y_aspi_quad__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
a_k_Y_aspi_quad__(1:n_lm_(1+nk_p_r),1+nk_p_r) = a_k_Y_aspi_quad_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
save(fname_aspi_mat ...
     ,'l_max_' ...
     ,'n_lm_','n_lm_max','n_lm_sum','n_lm_csum_','l_max_max','m_max_','n_m_max' ...
     ,'Y_l_val_','Y_m_val_','Y_k_val_','weight_Y_' ...
     ,'a_k_Y_aspi_quad_' ...
     ,'a_k_p_aspi_reco_' ...
     ,'a_k_Y_aspi_quad__' ...
     );
end;%if (~exist(fname_aspi_mat,'file'));
%%%%%%%%;
if ( exist(fname_aspi_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_aspi_mat));
load(fname_aspi_mat);
end;%if ( exist(fname_aspi_mat,'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/a_k_Y_aspi_quad_A',dir_aspi_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;
clf;
subplot(1,3,1); plot(Y_l_val_,log10(abs(a_k_Y_aspi_quad_)),'.'); xlabel('Y_l_val_','Interpreter','none'); ylabel('log10(abs(a))');
subplot(1,3,2); plot(Y_m_val_,log10(abs(a_k_Y_aspi_quad_)),'.'); xlabel('Y_m_val_','Interpreter','none'); ylabel('log10(abs(a))');
subplot(1,3,3); plot(Y_k_val_,log10(abs(a_k_Y_aspi_quad_)),'.'); xlabel('Y_k_val_','Interpreter','none'); ylabel('log10(abs(a))');
figbig;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/a_k_Y_aspi_quad_',dir_aspi_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;
clf;
imagesc_Y(k_p_r_max,n_k_p_r,k_p_r_,l_max_,log10(abs(a_k_Y_aspi_quad_)),[-10,0],colormap_beach());
xlabel('l_val','Interpreter','none');
ylabel('m_val','Interpreter','none');
zlabel('k_p_r','Interpreter','none');
axis vis3d; view([-65,+20]);
title('a_k_Y_aspi_quad_','Interpreter','none');
%figbig;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

fname_aspi_mat = sprintf('%s_mat/X_aspi.mat',dir_aspi_pm);
%%%%%%%%;
if ~exist(fname_aspi_mat,'file');
%%%%%%%%;
[ ...
 X_aspi_best ...
] = ...
register_spharm_to_spharm_wigner_wrap_1( ...
 n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,0 ...
,l_max_ ...
,a_k_Y_quad_ ...
,a_k_Y_aspi_quad_ ...
);
disp(sprintf(' %% X_aspi_best = %0.6f',X_aspi_best));
%%%%;
x_u_r___ = sqrt(x_u_0___.^2 + x_u_1___.^2 + x_u_2___.^2);
a_x_u_aspi_mask_ = a_x_u_aspi_base_.*(x_u_r___<sqrt(0.5)); %<-- some kind of a mask. ;
%%%%;
eta = pi/x_p_r_max; tmp_t = tic;
a_k_p_aspi_mask_ = xxnufft3d3(n_xxx_u,x_u_0___(:)*eta,x_u_1___(:)*eta,x_u_2___(:)*eta,a_x_u_aspi_mask_(:).*xxx_u_weight_(:),-1,1e-12,n_k_all,2*pi*k_c_0_all_/eta,2*pi*k_c_1_all_/eta,2*pi*k_c_2_all_/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_k_p_aspi_mask_ time %0.2fs',tmp_t));
%%%%;
tmp_t = tic;
[a_k_Y_aspi_mask_] = convert_k_p_to_spharm_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_p_aspi_mask_);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_aspi_mask_ time %0.2fs',tmp_t));
tmp_t = tic;
%%%%;
[ ...
 X_aspi_mask ...
] = ...
register_spharm_to_spharm_wigner_wrap_1( ...
 n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,0 ...
,l_max_ ...
,a_k_Y_quad_ ...
,a_k_Y_aspi_mask_ ...
);
disp(sprintf(' %% X_aspi_mask = %0.6f',X_aspi_mask));
%%%%;
save(fname_aspi_mat,'X_aspi_best','X_aspi_mask');
%%%%%%%%;
end;%if ~exist(fname_aspi_mat,'file');
%%%%%%%%;
if  exist(fname_aspi_mat,'file');
load(fname_aspi_mat);
end;%if  exist(fname_aspi_mat,'file');

verbose=1;
if (verbose); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

disp('returning'); return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
