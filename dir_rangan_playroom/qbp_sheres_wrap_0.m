function ...
[ ...
  parameter ...
] = ...
qbp_sheres_wrap_0( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_max ...
,l_max_ ...
,n_M ...
,M_k_p_wkM__ ...
,image_delta_x_prev_M_ ...
,image_delta_y_prev_M_ ...
,n_CTF ...
,index_nCTF_from_nM_ ...
,CTF_k_p_r_kC__ ...
,a_k_Y_base_yk_ ...
,a_k_Y_true_yk_ ...
);

flag_replot = 0; nf=0;
str_thisfunction = 'qbp_sheres_wrap_0';

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); k_p_r_=[]; end; na=na+1;
if (nargin<1+na); k_p_r_max=[]; end; na=na+1;
if (nargin<1+na); weight_3d_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); weight_2d_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); n_w_max=[]; end; na=na+1;
if (nargin<1+na); l_max_=[]; end; na=na+1;
if (nargin<1+na); n_M=[]; end; na=na+1;
if (nargin<1+na); M_k_p_wkM__=[]; end; na=na+1;
if (nargin<1+na); image_delta_x_prev_M_=[]; end; na=na+1;
if (nargin<1+na); image_delta_y_prev_M_=[]; end; na=na+1;
if (nargin<1+na); n_CTF=[]; end; na=na+1;
if (nargin<1+na); index_nCTF_from_nM_=[]; end; na=na+1;
if (nargin<1+na); CTF_k_p_r_kC__=[]; end; na=na+1;
if (nargin<1+na); a_k_Y_base_yk_=[]; end; na=na+1;
if (nargin<1+na); a_k_Y_true_yk_=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;%if isempty(parameter);
%%%%%%%%;
if (~isfield(parameter,'tolerance_master')); parameter.tolerance_master = 1e-2; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_recalc')); parameter.flag_recalc = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_verbose')); parameter.flag_verbose = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_rank_vs_tolerance')); parameter.flag_rank_vs_tolerance = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'tolerance_pm')); parameter.tolerance_pm = parameter.tolerance_master; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'rank_pm')); parameter.rank_pm = 10; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'dir_pm')); parameter.dir_pm = pwd; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'rseed')); parameter.rseed = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'k_p_r_max')); parameter.k_p_r_max = k_p_r_max; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'delta_r_max')); parameter.delta_r_max = 0.1; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_xcor_vs_Memp')); parameter.flag_xcor_vs_Memp = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'sigma_sheres')); parameter.sigma_sheres = 0.1; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'str_strategy_prefix')); parameter.str_strategy_prefix = ''; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'str_base')); parameter.str_base = ''; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'n_iteration')); parameter.n_iteration = 16; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'half_diameter_x_c')); parameter.half_diameter_x_c = 1.0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'n_x_u_pack')); parameter.n_x_u_pack = 64; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'date_diff_threshold')); parameter.date_diff_threshold = 2.5; end; %<-- parameter_bookmark. ;
%%%%%%%%;
tolerance_master = parameter.tolerance_master;
flag_recalc = parameter.flag_recalc;
flag_verbose = parameter.flag_verbose;
flag_rank_vs_tolerance = parameter.flag_rank_vs_tolerance;
tolerance_pm = parameter.tolerance_pm;
rank_pm = parameter.rank_pm;
dir_pm = parameter.dir_pm;
string_root = dir_pm(2:strfind(dir_pm,'rangan')-2);
rseed = parameter.rseed;
k_p_r_max = parameter.k_p_r_max;
delta_r_max = parameter.delta_r_max;
flag_xcor_vs_Memp = parameter.flag_xcor_vs_Memp;
sigma_sheres = parameter.sigma_sheres;
str_strategy = parameter.str_strategy_prefix;
str_base = parameter.str_base;
n_iteration = parameter.n_iteration;
half_diameter_x_c = parameter.half_diameter_x_c;
n_x_u_pack = parameter.n_x_u_pack;
date_diff_threshold = parameter.date_diff_threshold;
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_verbose); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
str_xvM = sprintf('X_2d_Memp_d1'); if (flag_xcor_vs_Memp==1); str_xvM = sprintf('X_2d_xcor_d0'); end;
str_nls = sprintf('nls%.2d',round(10*-log(sigma_sheres)));
str_delta_r_max = sprintf('t%.4d',floor(1000*delta_r_max));
str_tolerance_pm = sprintf('p%.2d',floor(10*-log10(tolerance_pm)));
str_rseed = sprintf('r%d',rseed);
str_xfix = sprintf('%s%s%s%s',str_strategy,str_delta_r_max,str_tolerance_pm,str_rseed);
str_base_use = ''; if ~isempty(str_base); str_base_use = sprintf('_%s',str_base); end;
fname_sheres_wrap_pre = sprintf('%s_mat/%s_%s%s%s',dir_pm,str_xvM,str_nls,str_xfix,str_base_use);
if (flag_verbose); disp(sprintf(' %% fname_sheres_wrap_pre: %s',fname_sheres_wrap_pre)); end;
%%%%%%%%;

%%%%%%%%;
[flag_sheres_wrap_skip,fname_sheres_wrap_mat] = open_fname_tmp(fname_sheres_wrap_pre,date_diff_threshold);
if (flag_recalc | ~flag_sheres_wrap_skip);
%%%%%%%%;
l_max_max = n_w_max/2 - 1;
if (l_max_max~=max(l_max_));
disp(sprintf(' %% Warning, l_max_max %d max(l_max_) %d in %s',l_max_max,max(l_max_),str_thisfunction));
end;%if (l_max_max~=max(l_max_));
assert(l_max_max==max(l_max_));
n_w_ = n_w_max*ones(n_k_p_r,1);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);
n_lm_ = (l_max_+1).^2;
n_lm_max = max(n_lm_);
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
%%%%%%%%;
% construct CTF of same size as images. ;
%%%%%%%%
CTF_k_p_wkC__ = zeros(n_w_sum,n_CTF);
for nCTF=0:n_CTF-1;
for nk_p_r=0:n_k_p_r-1;
n_w = n_w_(1+nk_p_r);
tmp_index_ = n_w_csum_(1+nk_p_r) + [0:n_w-1];
CTF_k_p_wkC__(1+tmp_index_,1+nCTF) = CTF_k_p_r_kC__(1+nk_p_r,1+nCTF);
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nCTF=0:n_CTF-1;
%%%%%%%%;
% Initialize volume. ;
%%%%%%%%;
if isempty(a_k_Y_base_yk_);
if (flag_verbose); disp(sprintf(' %% creating a_k_Y_base_yk_')); end;
rng(rseed);
euler_polar_a_rand_M_ = 1*pi*rand(n_M,1);
euler_azimu_b_rand_M_ = 2*pi*rand(n_M,1);
euler_gamma_z_rand_M_ = 2*pi*rand(n_M,1);
image_delta_x_rand_M_ = zeros(n_M,1);
image_delta_y_rand_M_ = zeros(n_M,1);
image_I_value_rand_M_ = ones(n_M,1);
qbp_eps = tolerance_master;
[ ...
 a_k_Y_base_yk_ ...
] = ...
qbp_6( ...
 qbp_eps ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,n_M ...
,M_k_p_wkM__ ...
,index_nCTF_from_nM_ ...
,CTF_k_p_wkC__ ...
,euler_polar_a_rand_M_ ...
,euler_azimu_b_rand_M_ ...
,euler_gamma_z_rand_M_ ...
,image_delta_x_rand_M_ ...
,image_delta_y_rand_M_ ...
,image_I_value_rand_M_ ...
);
end;%if isempty(a_k_Y_base_yk_);
%%%%%%%%;
a_k_Y_zero_yki__ = zeros(n_lm_sum,n_iteration);
X_best_zero_alig_i_ = zeros(n_iteration,1);
a_k_Y_0qbp_yki__ = zeros(n_lm_sum,n_iteration);
X_best_0qbp_alig_i_ = zeros(n_iteration,1);
a_k_Y_sher_yki__ = zeros(n_lm_sum,n_iteration);
X_best_sher_alig_i_ = zeros(n_iteration,1);
for niteration=0:n_iteration-1;
%%%%;
index_ncluster_from_nCTF_ = [];
[ ...
 parameter ...
,a_k_Y_zero_yk_ ...
,a_k_Y_zero_yk__ ...
,a_k_Y_0qbp_yk_ ...
,a_k_Y_0qbp_yk__ ...
,a_k_Y_sher_yk_ ...
,a_k_Y_sher_yk__ ...
] = ...
qbp_sheres_0( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_max ...
,l_max_ ...
,n_M ...
,M_k_p_wkM__ ...
,image_delta_x_prev_M_ ...
,image_delta_y_prev_M_ ...
,n_CTF ...
,index_nCTF_from_nM_ ...
,CTF_k_p_r_kC__ ...
,index_ncluster_from_nCTF_ ...
,a_k_Y_base_yk_ ...
);
%%%%;
[ ... 
 a_k_Y_zero_alig_yk_ ...
,X_best_zero_alig ...
] = ...
spharm_register_and_rotate_2( ...
 n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_true_yk_ ...
,a_k_Y_zero_yk_ ...
);
[ ... 
 a_k_Y_0qbp_alig_yk_ ...
,X_best_0qbp_alig ...
] = ...
spharm_register_and_rotate_2( ...
 n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_true_yk_ ...
,a_k_Y_0qbp_yk_ ...
);
[ ... 
 a_k_Y_sher_alig_yk_ ...
,X_best_sher_alig ...
] = ...
spharm_register_and_rotate_2( ...
 n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_true_yk_ ...
,a_k_Y_sher_yk_ ...
);
%%%%;
a_k_Y_base_yk_ = a_k_Y_sher_yk_;
%%%%;
a_k_Y_zero_yki__(:,1+niteration) = a_k_Y_zero_yk_;
X_best_zero_alig_i_(1+niteration) = X_best_zero_alig;
a_k_Y_0qbp_yki__(:,1+niteration) = a_k_Y_0qbp_yk_;
X_best_0qbp_alig_i_(1+niteration) = X_best_0qbp_alig;
a_k_Y_sher_yki__(:,1+niteration) = a_k_Y_sher_yk_;
X_best_sher_alig_i_(1+niteration) = X_best_sher_alig;
%%%%%%%%;
end;%for niteration=0:n_iteration-1;
%%%%%%%%;
save(fname_sheres_wrap_mat ...
,'parameter' ...
,'a_k_Y_zero_yki__' ...
,'X_best_zero_alig_i_' ...
,'a_k_Y_0qbp_yki__' ...
,'X_best_0qbp_alig_i_' ...
,'a_k_Y_sher_yki__' ...
,'X_best_sher_alig_i_' ...
);     
close_fname_tmp(fname_sheres_wrap_pre);
end;%if (flag_recalc | ~flag_sheres_wrap_skip);
%%%%%%%%;

%%%%%%%%;
if  exist(fname_sheres_wrap_mat,'file');
fname_fig = sprintf('%s_mat/%s_%s%s%s_snapshot',dir_pm,str_xvM,str_nls,str_xfix,str_base_use);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
tmp_ = load(fname_sheres_wrap_mat); 
a_k_Y_sher_yk_ = tmp_.a_k_Y_sher_yki__(:,end);
X_best_sher_alig = tmp_.X_best_sher_alig_i_(end);
clear tmp_;
[ ... 
 a_k_Y_sher_alig_yk_ ...
,X_best_sher_alig ...
] = ...
spharm_register_and_rotate_2( ...
 n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_true_yk_ ...
,a_k_Y_sher_yk_ ...
);
%%%%;
sample_sphere_k_eq_d = 1/(2*pi);
[ ... 
  a_x_u_sher_alig_ ...
] = ...
convert_spharm_to_x_c_3( ...
 sample_sphere_k_eq_d ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_sher_alig_yk_ ...
,half_diameter_x_c ...
,n_x_u_pack ...
);
figure(1+nf);nf=nf+1;clf;figbig;
isosurface_f_x_u_0(reshape(real(a_x_u_sher_alig_),n_x_u_pack,n_x_u_pack,n_x_u_pack),[90,95,99]);
title(sprintf('%s: a_x_u_sher_alig_: corr %0.4f',fname_fig,X_best_sher_alig),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
end;%if  exist(fname_sheres_wrap_mat,'file');
%%%%%%%%;

%%%%%%%%;
flag_continue = ...
  1 ...
 & isfield(parameter,'Pixel_Spacing') ...
 & isfield(parameter,'n_x_M_u') ...
 & isfield(parameter,'n_x_u_pack') ...
 & isfield(parameter,'a_x_u_base_') ...
 & isfield(parameter,'a_x_u_reco_') ...
 & isfield(parameter,'n_k_all') ...
 & isfield(parameter,'n_k_all_csum_') ...
 & isfield(parameter,'weight_3d_k_all_') ...
 & isfield(parameter,'k_c_0_all_') ...
 & isfield(parameter,'k_c_1_all_') ...
 & isfield(parameter,'k_c_2_all_') ...
  ;
if flag_continue;
if  exist(fname_sheres_wrap_mat,'file');
fname_sheres_wrap_fsc_pre = sprintf('%s_mat/%s_%s%s%s_fsc',dir_pm,str_xvM,str_nls,str_xfix,str_base_use);
[flag_sheres_wrap_fsc_skip,fname_sheres_wrap_fsc_mat] = open_fname_tmp(fname_sheres_wrap_fsc_pre);
if (flag_recalc | ~flag_sheres_wrap_fsc_skip);
%%%%;
tmp_ = load(fname_sheres_wrap_mat); 
a_k_Y_sher_yk_ = tmp_.a_k_Y_sher_yki__(:,end);
X_best_sher_alig = tmp_.X_best_sher_alig_i_(end);
clear tmp_;
%%%%;
Pixel_Spacing = parameter.Pixel_Spacing;
n_x_M_u = parameter.n_x_M_u;
n_x_u_pack = parameter.n_x_u_pack;
a_x_u_base_ = parameter.a_x_u_base_;
a_x_u_reco_ = parameter.a_x_u_reco_;
n_k_all = parameter.n_k_all;
n_k_all_csum_ = parameter.n_k_all_csum_;
weight_3d_k_all_ = parameter.weight_3d_k_all_;
k_c_0_all_ = parameter.k_c_0_all_;
k_c_1_all_ = parameter.k_c_1_all_;
k_c_2_all_ = parameter.k_c_2_all_;
%%%%;
parameter_fsc = struct('type','parameter'); parameter_fsc.Pixel_Spacing = Pixel_Spacing;
[ ...
 parameter_fsc ...
,fsc_reco_stab_k_ ...
,corr_base_vs_reco_stab ...
,corr_reco_vs_reco_stab ...
,fsc_crop_reco_stab_kx__ ...
,corr_crop_base_vs_crop_reco_stab_x_ ...
,corr_full_base_vs_crop_reco_stab_x_ ...
,corr_crop_reco_vs_crop_reco_stab_x_ ...
,corr_full_reco_vs_crop_reco_stab_x_ ...
,k_Ainv_p_r_ ...
,k_Ainv_p_r_max ...
,kinv_A_p_r_ ...
] = ...
crop_from_x_c_0( ...
 parameter_fsc ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_true_yk_ ...
,a_k_Y_sher_yk_ ...
,[] ...
,half_diameter_x_c ...
,n_x_M_u ...
,n_x_u_pack ...
,a_x_u_base_ ...
,a_x_u_reco_ ...
,[] ...
,n_k_all ...
,n_k_all_csum_ ...
,weight_3d_k_all_ ...
,k_c_0_all_ ...
,k_c_1_all_ ...
,k_c_2_all_ ...
);
%%%%;
save(fname_sheres_wrap_fsc_mat ...
     ,'parameter_fsc' ...
,'fsc_reco_stab_k_' ...
,'corr_base_vs_reco_stab' ...
,'corr_reco_vs_reco_stab' ...
,'fsc_crop_reco_stab_kx__' ...
,'corr_crop_base_vs_crop_reco_stab_x_' ...
,'corr_full_base_vs_crop_reco_stab_x_' ...
,'corr_crop_reco_vs_crop_reco_stab_x_' ...
,'corr_full_reco_vs_crop_reco_stab_x_' ...
,'k_Ainv_p_r_' ...
,'k_Ainv_p_r_max' ...
,'kinv_A_p_r_' ...
     );
close_fname_tmp(fname_sheres_wrap_fsc_pre);
end;%if (flag_recalc | ~flag_sheres_wrap_fsc_skip);
end;%if  exist(fname_sheres_wrap_mat,'file');
end;%if flag_continue;
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_verbose); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
