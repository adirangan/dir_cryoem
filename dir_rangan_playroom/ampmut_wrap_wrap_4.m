function ...
parameter = ...
ampmut_wrap_wrap_4( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_ ...
,dat_n_UX_rank ...
,UX__ ...
,X_weight_r_ ...
,n_M ...
,M_k_p__ ...
,n_ctf ...
,n_CTF_rank ...
,CTF_index_ ...
,CTF_k_p__ ...
,l_max_ ...
,a_k_Y_true_ ...
,euler_polar_a_true_ ...
,euler_azimu_b_true_ ...
,euler_gamma_z_true_ ...
,image_delta_x_true_ ...
,image_delta_y_true_ ...
);

%{
%%%%%%%%;
% Should be unnecessary. ;
%%%%%%%%;
platform = 'rusty';%platform = 'access1';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
%}

if isempty(parameter);
parameter = struct('type','parameter');
end;%if isempty(parameter);
%%%%%%%%;
if (~isfield(parameter,'dir_pm')); parameter.dir_pm = pwd; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'rseed')); parameter.rseed = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'k_p_r_max')); parameter.k_p_r_max = k_p_r_max; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'sample_sphere_k_eq_d')); parameter.sample_sphere_k_eq_d = 1/(2*pi); end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'half_diameter_x_c')); parameter.half_diameter_x_c = 1.0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'n_x_u_pack')); parameter.n_x_u_pack = 64; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'delta_r_max')); parameter.delta_r_max = 0.1; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'delta_r_upb')); parameter.delta_r_upb = 0.2; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'cg_lsq_n_order')); parameter.cg_lsq_n_order = 5; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'date_diff_threshold')); parameter.date_diff_threshold = 0.25; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_force_create_mat')); parameter.flag_force_create_mat = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_force_create_tmp')); parameter.flag_force_create_tmp = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_alternate_MS_vs_SM')); parameter.flag_alternate_MS_vs_SM = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_local_exclusion')); parameter.flag_local_exclusion = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_qbp_vs_lsq')); parameter.flag_qbp_vs_lsq = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'str_strategy_prefix')); parameter.str_strategy_prefix = ''; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'n_complete_calculation')); parameter.n_complete_calculation = 0; end; %<-- parameter_bookmark. ;
%%%%%%%%;
flag_alternate_MS_vs_SM = parameter.flag_alternate_MS_vs_SM;
flag_local_exclusion = parameter.flag_local_exclusion;
flag_qbp_vs_lsq = parameter.flag_qbp_vs_lsq;
str_strategy = parameter.str_strategy_prefix;
if (flag_qbp_vs_lsq); str_strategy = sprintf('%sq1',str_strategy); end;
if (flag_alternate_MS_vs_SM); str_strategy = sprintf('%sa1',str_strategy); end;
if (flag_local_exclusion); str_strategy = sprintf('%se1',str_strategy); end;
dir_pm = parameter.dir_pm;
string_root = dir_pm(2:strfind(dir_pm,'rangan')-2);

rseed = parameter.rseed;
sample_sphere_k_eq_d = parameter.sample_sphere_k_eq_d;
half_diameter_x_c = parameter.half_diameter_x_c;
n_x_u_pack = parameter.n_x_u_pack;
delta_r_max = parameter.delta_r_max;
delta_r_upb = parameter.delta_r_upb;
cg_lsq_n_order = parameter.cg_lsq_n_order;
date_diff_threshold = parameter.date_diff_threshold;
flag_force_create_mat = parameter.flag_force_create_mat;
flag_force_create_tmp = parameter.flag_force_create_tmp;

XA_fname_mat = sprintf('%s_mat/X_2d_Memp_d1_%st%.4dn%.2dr%d.mat',dir_pm,str_strategy,floor(1000*delta_r_max),dat_n_UX_rank,rseed);
XA_fname_snapshot_jpg = sprintf('%s_mat/X_2d_Memp_d1_%st%.4dn%.2dr%d_snapshot.jpg',dir_pm,str_strategy,floor(1000*delta_r_max),dat_n_UX_rank,rseed);
XB_fname_mat = sprintf('%s_mat/X_2d_xcor_d0_%st%.4dn%.2dr%d.mat',dir_pm,str_strategy,floor(1000*delta_r_max),dat_n_UX_rank,rseed);
XB_fname_snapshot_jpg = sprintf('%s_mat/X_2d_xcor_d0_%st%.4dn%.2dr%d_snapshot.jpg',dir_pm,str_strategy,floor(1000*delta_r_max),dat_n_UX_rank,rseed);
flag_skip_all = 0;
if ( ...
        exist(XA_fname_mat,'file') ...
     &  exist(XB_fname_snapshot_jpg,'file') ...
     &  exist(XB_fname_mat,'file') ...
     &  exist(XB_fname_snapshot_jpg,'file') ...
   ); flag_skip_all = 1; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if ~flag_skip_all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

verbose=2;
if (verbose); disp(sprintf(' %% [entering ampmut_wrap_wrap_4]')); end;

n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);
l_max_max = max(l_max_);
n_lm_ = (l_max_+1).^2;
n_lm_max = max(n_lm_);
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
%%%%%%%%;
n_UX_rank = n_k_p_r-1; %<-- just to check dimensions. ;
CTF_k_p_r__ = zeros(n_k_p_r,n_ctf);
for nctf=0:n_ctf-1;
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
CTF_k_p_r__(1+nk_p_r,1+nctf) = mean(CTF_k_p__(1+tmp_index_,1+nctf));
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nctf=0:n_ctf-1;
CTF_avg_k_p_ = mean(CTF_k_p__,2);
%imagesc_p(n_k_p_r,k_p_r_,n_w_,sum(n_w_),real(CTF_avg_k_p_(:)),[-1,+1],colormap_beach());
CTF_avg_k_p_r_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
CTF_avg_k_p_r_(1+nk_p_r) = mean(CTF_avg_k_p_(1+tmp_index_));
end;%for nk_p_r=0:n_k_p_r-1;
CTF_k_p_r_xcor__ = CTF_k_p_r__(:,1+CTF_index_(1+(0:n_M-1))) * transpose(CTF_k_p_r__(:,1+CTF_index_(1+(0:n_M-1)))) / n_M;
%%%%%%%%;

%%%%%%%%;
% Now run the basic ampmut using the empirical-principal-modes. ;
%%%%%%%%;
XA_fname_pre = sprintf('%s_mat/X_2d_Memp_d1_%st%.4dn%.2dr%d',dir_pm,str_strategy,floor(1000*delta_r_max),dat_n_UX_rank,rseed);
[XA_flag_skip,XA_fname_mat] = open_fname_tmp(XA_fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~XA_flag_skip;
%%%%%%%%;
% First calculate empirical principal-modes. ;
%%%%%%%%;
[ ...
 X_2d_Memp_d1__ ...
,X_2d_Memp_d1_weight_r_ ...
] = ...
principled_marching_empirical_cost_matrix_0( ...
 n_k_p_r ...
,k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_ ...
,n_M ...
,M_k_p__ ...
);
%%%%%%%%;
[UX_2d_Memp_d1__,SX_2d_Memp_d1__,VX_2d_Memp_d1__] = svds(X_2d_Memp_d1__,n_UX_rank); SX_2d_Memp_d1_ = diag(SX_2d_Memp_d1__);
if (verbose); disp(sprintf(' %% cumsum(SX_2d_Memp_d1_): ')); disp(sprintf(' %% %0.4f',cumsum(SX_2d_Memp_d1_/sum(SX_2d_Memp_d1_)))); end;
%%%%%%%%;
%try;
disp(sprintf(' %% %s not found, creating',XA_fname_mat));
XA_parameter = struct('type','parameter');
XA_parameter = parameter; %<-- copy parameter and all its fields. ;
XA_parameter.delta_r_upb = delta_r_upb;
XA_parameter.delta_r_max = delta_r_max;
XA_parameter.rseed = rseed;
XA_parameter.fname_pre = XA_fname_pre;
XA_parameter.flag_alternate_MS_vs_SM = flag_alternate_MS_vs_SM;
XA_parameter.flag_local_exclusion = flag_local_exclusion;
XA_parameter.flag_qbp_vs_lsq = flag_qbp_vs_lsq;
ampmut_wrap_4( ...
 XA_parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_max ...
,dat_n_UX_rank ...
,UX_2d_Memp_d1__ ...
,X_2d_Memp_d1_weight_r_ ...
,n_M ...
,M_k_p__ ...
,n_CTF_rank ...
,CTF_index_ ...
,CTF_k_p_r__ ...
);
parameter.n_complete_calculation = parameter.n_complete_calculation + 1;
%catch; disp(sprintf(' %% WARNING: error generating %s',XA_fname_mat)); end;%try;
close_fname_tmp(XA_fname_pre);
end;%if ~XA_flag_skip;

%%%%%%%%;
if ( exist(XA_fname_mat,'file'));
%%%%%%%%;
% First collect statistics. ;
%%%%%%%%;
if (verbose); disp(sprintf(' %% %s found, aligning',XA_fname_mat)); end;
tmp_XA_ = load(XA_fname_mat);
if (~isfield(tmp_XA_.parameter,'fname_pre')); tmp_XA_.parameter.fname_pre = XA_fname_pre; end; %<-- parameter_bookmark. ;
%%%%%%%%;
tmp_XA_fname_pre = sprintf('%s_align_a_CTF_avg_UX_Y_',tmp_XA_.parameter.fname_pre);
tmp_XA_fname_pre = rootswitch(tmp_XA_fname_pre,string_root,'rangan');
tmp_XA_.parameter.fname_align_a_CTF_avg_UX_Y_pre = tmp_XA_fname_pre;
[tmp_XA_flag_skip,tmp_XA_fname_mat] = open_fname_tmp(tmp_XA_fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~tmp_XA_flag_skip;
%%%%%%%%;
% recalculate empirical principal-modes if necessary. ;
%%%%%%%%;
if ~exist('UX_2d_Memp_d1__','var');
[ ...
 X_2d_Memp_d1__ ...
,X_2d_Memp_d1_weight_r_ ...
] = ...
principled_marching_empirical_cost_matrix_0( ...
 n_k_p_r ...
,k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_ ...
,n_M ...
,M_k_p__ ...
);
[UX_2d_Memp_d1__,SX_2d_Memp_d1__,VX_2d_Memp_d1__] = svds(X_2d_Memp_d1__,n_UX_rank); SX_2d_Memp_d1_ = diag(SX_2d_Memp_d1__);
end;%if ~exist('UX_2d_Memp_d1__','var');
%%%%%%%%;
% First form a_CTF_avg_UX_2d_Memp_d1_Y_quad__ ;
%%%%%%%%;
a_CTF_avg_UX_2d_Memp_d1_Y_quad__ = zeros(n_lm_max,n_UX_rank);
for nUX_rank=0:n_UX_rank-1;
for nk_p_r=0:n_k_p_r-1;
tmp_l_max = l_max_(1+nk_p_r);
tmp_n_lm = (tmp_l_max+1).^2;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:tmp_n_lm-1);
a_CTF_avg_UX_2d_Memp_d1_Y_quad__(1:tmp_n_lm,1+nUX_rank) = ...
a_CTF_avg_UX_2d_Memp_d1_Y_quad__(1:tmp_n_lm,1+nUX_rank) + ...
UX_2d_Memp_d1__(1+nk_p_r,1+nUX_rank)*X_2d_Memp_d1_weight_r_(1+nk_p_r)*a_k_Y_true_(1+tmp_index_)*CTF_avg_k_p_r_(1+nk_p_r); ...
%<-- use average CTF here, under the assumption that a_CTF_UX_2d_Memp_d1_Y_quad_ will be used alone. ;
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nUX_rank=0:n_UX_rank-1;
%%%%%%%%;
try;
disp(sprintf(' %% %s not found, creating',tmp_XA_fname_mat));
%%%%;
ampmut_align_to_a_CTF_avg_UX_Y_0( ...
 tmp_XA_.parameter ...
,l_max_max ...
,dat_n_UX_rank ...
,reshape(a_CTF_avg_UX_2d_Memp_d1_Y_quad__(:,1:dat_n_UX_rank),[n_lm_max*dat_n_UX_rank,1]) ...
,n_M ...
,euler_polar_a_true_ ...
,euler_azimu_b_true_ ...
,euler_gamma_z_true_ ...
,image_delta_x_true_ ...
,image_delta_y_true_ ...
,[] ...
,tmp_XA_.a_CTF_avg_UX_Y__ ...
,tmp_XA_.euler_polar_a__ ...
,tmp_XA_.euler_azimu_b__ ...
,tmp_XA_.euler_gamma_z__ ...
,tmp_XA_.image_delta_x_acc__ + tmp_XA_.image_delta_x_upd__ ...
,tmp_XA_.image_delta_y_acc__ + tmp_XA_.image_delta_y_upd__ ...
);
catch; disp(sprintf(' %% WARNING: error generating %s',tmp_XA_fname_mat)); end;%try;
close_fname_tmp(tmp_XA_fname_pre);
end;%if ~tmp_XA_flag_skip;
%%%%%%%%;
tmp_XA_fname_pre = sprintf('%s_align_a_k_Y_',tmp_XA_.parameter.fname_pre);
tmp_XA_fname_pre = rootswitch(tmp_XA_fname_pre,string_root,'rangan');
tmp_XA_.parameter.fname_align_a_k_Y_pre = tmp_XA_fname_pre;
[tmp_XA_flag_skip,tmp_XA_fname_mat] = open_fname_tmp(tmp_XA_fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~tmp_XA_flag_skip;
try;
disp(sprintf(' %% %s not found, creating',tmp_XA_fname_mat));
ampmut_align_to_a_k_Y_0( ...
 tmp_XA_.parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_M ...
,M_k_p__ ...
,CTF_index_ ...
,CTF_k_p__ ...
,l_max_ ...
,a_k_Y_true_ ...
,[] ...
,tmp_XA_.euler_polar_a__ ...
,tmp_XA_.euler_azimu_b__ ...
,tmp_XA_.euler_gamma_z__ ...
,tmp_XA_.image_delta_x_acc__ + tmp_XA_.image_delta_x_upd__ ...
,tmp_XA_.image_delta_y_acc__ + tmp_XA_.image_delta_y_upd__ ...
);
catch; disp(sprintf(' %% WARNING: error generating %s',tmp_XA_fname_mat)); end;%try;
close_fname_tmp(tmp_XA_fname_pre);
end;%if ~tmp_XA_flag_skip;
%%%%%%%%;
% Now use the final step to reconstruct the molecule. ;
%%%%%%%%;
tmp_XA_euler_polar_a_ = tmp_XA_.euler_polar_a__(:,end);
tmp_XA_euler_azimu_b_ = tmp_XA_.euler_azimu_b__(:,end);
tmp_XA_euler_gamma_z_ = tmp_XA_.euler_gamma_z__(:,end);
tmp_XA_image_delta_x_ = tmp_XA_.image_delta_x_acc__(:,end) + tmp_XA_.image_delta_x_upd__(:,end);
tmp_XA_image_delta_y_ = tmp_XA_.image_delta_y_acc__(:,end) + tmp_XA_.image_delta_y_upd__(:,end);
fname_XA_k_Y_mat = sprintf('%s_a_k_Y_.mat',tmp_XA_.parameter.fname_pre);
fname_XA_k_Y_mat = rootswitch(fname_XA_k_Y_mat,string_root,'rangan');
if (~exist(fname_XA_k_Y_mat,'file'));
if (verbose); disp(sprintf(' %% %s not found, creating',fname_XA_k_Y_mat)); end;
a_k_Y_reco_ = ...
cg_lsq_4( ...
 cg_lsq_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,n_M ...
,M_k_p__ ...
,CTF_index_ ...
,CTF_k_p__ ...
,tmp_XA_euler_polar_a_ ...
,tmp_XA_euler_azimu_b_ ...
,tmp_XA_euler_gamma_z_ ...
,tmp_XA_image_delta_x_ ...
,tmp_XA_image_delta_y_ ...
);
save(fname_XA_k_Y_mat,'a_k_Y_reco_');
end;%if (~exist(fname_XA_k_Y_mat,'file'));

%%%%%%%%;
if ( exist(fname_XA_k_Y_mat,'file'));
%%%%%%%%;

%%%%%%%%;
% Now run ampmut once again, this time using the updated-principal-modes. ;
%%%%%%%%;
XB_fname_pre = sprintf('%s_mat/X_2d_xcor_d0_%st%.4dn%.2dr%d',dir_pm,str_strategy,floor(1000*delta_r_max),dat_n_UX_rank,rseed);
[XB_flag_skip,XB_fname_mat] = open_fname_tmp(XB_fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~XB_flag_skip;
%%%%%%%%;
% First use the final step from XA to recalculate the principal-modes. ;
%%%%%%%%;
tmp_ = load(fname_XA_k_Y_mat); XA_a_k_Y_reco_ = tmp_.a_k_Y_reco_; clear tmp_;
XA_a_k_Y_reco__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
XA_a_k_Y_reco__(1:n_lm_(1+nk_p_r),1+nk_p_r) = XA_a_k_Y_reco_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
[ ...
 X_2d_xcor_d0__ ...
,X_2d_xcor_d0_weight_r_ ...
] = ...
principled_marching_cost_matrix_3( ...
 n_k_p_r ...
,weight_2d_k_p_r_ ...
,l_max_max ...
,XA_a_k_Y_reco__ ...
,CTF_k_p_r_xcor__ ...
);
[UX_2d_xcor_d0__,SX_2d_xcor_d0__,VX_2d_xcor_d0__] = svds(X_2d_xcor_d0__,n_UX_rank); SX_2d_xcor_d0_ = diag(SX_2d_xcor_d0__);
if (verbose); disp(sprintf(' %% cumsum(SX_2d_xcor_d0_): ')); disp(sprintf(' %% %0.4f',cumsum(SX_2d_xcor_d0_/sum(SX_2d_xcor_d0_)))); end;
%%%%%%%%;
try;
disp(sprintf(' %% %s not found, creating',XB_fname_mat));
XB_parameter = struct('type','parameter');
XB_parameter = parameter; %<-- copy parameter and all its fields. ;
XB_parameter.delta_r_upb = delta_r_upb;
XB_parameter.delta_r_max = delta_r_max;
XB_parameter.rseed = rseed;
XB_parameter.fname_pre = XB_fname_pre;
XB_parameter.flag_alternate_MS_vs_SM = flag_alternate_MS_vs_SM;
XB_parameter.flag_local_exclusion = flag_local_exclusion;
XB_parameter.flag_qbp_vs_lsq = flag_qbp_vs_lsq;
ampmut_wrap_4( ...
 XB_parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_max ...
,dat_n_UX_rank ...
,UX_2d_xcor_d0__ ...
,X_2d_xcor_d0_weight_r_ ...
,n_M ...
,M_k_p__ ...
,n_CTF_rank ...
,CTF_index_ ...
,CTF_k_p_r__ ...
,tmp_XA_euler_polar_a_ ...
,tmp_XA_euler_azimu_b_ ...
,tmp_XA_euler_gamma_z_ ...
,tmp_XA_image_delta_x_ ...
,tmp_XA_image_delta_y_ ...
);
parameter.n_complete_calculation = parameter.n_complete_calculation + 1;
catch; disp(sprintf(' %% WARNING: error generating %s',XB_fname_mat)); end;%try;
close_fname_tmp(XB_fname_pre);
end;%if ~XB_flag_skip;

%%%%%%%%;
if ( exist(XB_fname_mat,'file'));
%%%%%%%%;
% First collect statistics. ;
%%%%%%%%;
if (verbose); disp(sprintf(' %% %s found, aligning',XB_fname_mat)); end;
tmp_XB_ = load(XB_fname_mat);
if (~isfield(tmp_XB_.parameter,'fname_pre')); tmp_XB_.parameter.fname_pre = XB_fname_pre; end; %<-- parameter_bookmark. ;
%%%%%%%%;
tmp_XB_fname_pre = sprintf('%s_align_a_CTF_avg_UX_Y_',tmp_XB_.parameter.fname_pre);
tmp_XB_fname_pre = rootswitch(tmp_XB_fname_pre,string_root,'rangan');
tmp_XB_.parameter.fname_align_a_CTF_avg_UX_Y_pre = tmp_XB_fname_pre;
[tmp_XB_flag_skip,tmp_XB_fname_mat] = open_fname_tmp(tmp_XB_fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~tmp_XB_flag_skip;
try;
%%%%%%%%;
% recalculate updated principal-modes if necessary. ;
%%%%%%%%;
if ~exist('UX_2d_xcor_d0__','var');
tmp_ = load(fname_XA_k_Y_mat); XA_a_k_Y_reco_ = tmp_.a_k_Y_reco_; clear tmp_;
XA_a_k_Y_reco__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
XA_a_k_Y_reco__(1:n_lm_(1+nk_p_r),1+nk_p_r) = XA_a_k_Y_reco_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
[ ...
 X_2d_xcor_d0__ ...
,X_2d_xcor_d0_weight_r_ ...
] = ...
principled_marching_cost_matrix_3( ...
 n_k_p_r ...
,weight_2d_k_p_r_ ...
,l_max_max ...
,XA_a_k_Y_reco__ ...
,CTF_k_p_r_xcor__ ...
);
[UX_2d_xcor_d0__,SX_2d_xcor_d0__,VX_2d_xcor_d0__] = svds(X_2d_xcor_d0__,n_UX_rank); SX_2d_xcor_d0_ = diag(SX_2d_xcor_d0__);
end;%if ~exist('UX_2d_xcor_d0__','var');
%%%%%%%%;
% First form a_CTF_avg_UX_2d_xcor_d0_Y_quad__ ;
%%%%%%%%;
a_CTF_avg_UX_2d_xcor_d0_Y_quad__ = zeros(n_lm_max,n_UX_rank);
for nUX_rank=0:n_UX_rank-1;
for nk_p_r=0:n_k_p_r-1;
tmp_l_max = l_max_(1+nk_p_r);
tmp_n_lm = (tmp_l_max+1).^2;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:tmp_n_lm-1);
a_CTF_avg_UX_2d_xcor_d0_Y_quad__(1:tmp_n_lm,1+nUX_rank) = ...
a_CTF_avg_UX_2d_xcor_d0_Y_quad__(1:tmp_n_lm,1+nUX_rank) + ...
UX_2d_xcor_d0__(1+nk_p_r,1+nUX_rank)*X_2d_xcor_d0_weight_r_(1+nk_p_r)*a_k_Y_true_(1+tmp_index_)*CTF_avg_k_p_r_(1+nk_p_r); ...
%<-- use average CTF here, under the assumption that a_CTF_UX_2d_xcor_d0_Y_quad_ will be used alone. ;
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nUX_rank=0:n_UX_rank-1;
%%%%%%%%;
disp(sprintf(' %% %s not found, creating',tmp_XB_fname_mat));
ampmut_align_to_a_CTF_avg_UX_Y_0( ...
 tmp_XB_.parameter ...
,l_max_max ...
,dat_n_UX_rank ...
,reshape(a_CTF_avg_UX_2d_xcor_d0_Y_quad__(:,1:dat_n_UX_rank),[n_lm_max*dat_n_UX_rank,1]) ...
,n_M ...
,euler_polar_a_true_ ...
,euler_azimu_b_true_ ...
,euler_gamma_z_true_ ...
,image_delta_x_true_ ...
,image_delta_y_true_ ...
,[] ...
,tmp_XB_.a_CTF_avg_UX_Y__ ...
,tmp_XB_.euler_polar_a__ ...
,tmp_XB_.euler_azimu_b__ ...
,tmp_XB_.euler_gamma_z__ ...
,tmp_XB_.image_delta_x_acc__ + tmp_XB_.image_delta_x_upd__ ...
,tmp_XB_.image_delta_y_acc__ + tmp_XB_.image_delta_y_upd__ ...
);
catch; disp(sprintf(' %% WARNING: error generating %s',tmp_XB_fname_mat)); end;%try;
close_fname_tmp(tmp_XB_fname_pre);
end;%if ~tmp_XB_flag_skip;
%%%%%%%%;
tmp_XB_fname_pre = sprintf('%s_align_a_k_Y_',tmp_XB_.parameter.fname_pre);
tmp_XB_fname_pre = rootswitch(tmp_XB_fname_pre,string_root,'rangan');
tmp_XB_.parameter.fname_align_a_k_Y_pre = tmp_XB_fname_pre;
[tmp_XB_flag_skip,tmp_XB_fname_mat] = open_fname_tmp(tmp_XB_fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~tmp_XB_flag_skip;
try;
disp(sprintf(' %% %s not found, creating',tmp_XB_fname_mat));
ampmut_align_to_a_k_Y_0( ...
 tmp_XB_.parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_M ...
,M_k_p__ ...
,CTF_index_ ...
,CTF_k_p__ ...
,l_max_ ...
,a_k_Y_true_ ...
,[] ...
,tmp_XB_.euler_polar_a__ ...
,tmp_XB_.euler_azimu_b__ ...
,tmp_XB_.euler_gamma_z__ ...
,tmp_XB_.image_delta_x_acc__ + tmp_XB_.image_delta_x_upd__ ...
,tmp_XB_.image_delta_y_acc__ + tmp_XB_.image_delta_y_upd__ ...
);
catch; disp(sprintf(' %% WARNING: error generating %s',tmp_XB_fname_mat)); end;%try;
close_fname_tmp(tmp_XB_fname_pre);
end;%if ~tmp_XB_flag_skip;
%%%%%%%%;
% Now use the final step to reconstruct the molecule. ;
%%%%%%%%;
tmp_XB_euler_polar_a_ = tmp_XB_.euler_polar_a__(:,end);
tmp_XB_euler_azimu_b_ = tmp_XB_.euler_azimu_b__(:,end);
tmp_XB_euler_gamma_z_ = tmp_XB_.euler_gamma_z__(:,end);
tmp_XB_image_delta_x_ = tmp_XB_.image_delta_x_acc__(:,end) + tmp_XB_.image_delta_x_upd__(:,end);
tmp_XB_image_delta_y_ = tmp_XB_.image_delta_y_acc__(:,end) + tmp_XB_.image_delta_y_upd__(:,end);
fname_XB_k_Y_mat = sprintf('%s_a_k_Y_.mat',tmp_XB_.parameter.fname_pre);
fname_XB_k_Y_mat = rootswitch(fname_XB_k_Y_mat,string_root,'rangan');
if (~exist(fname_XB_k_Y_mat,'file'));
if (verbose); disp(sprintf(' %% %s not found, creating',fname_XB_k_Y_mat)); end;
a_k_Y_reco_ = ...
cg_lsq_4( ...
 cg_lsq_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,n_M ...
,M_k_p__ ...
,CTF_index_ ...
,CTF_k_p__ ...
,tmp_XB_euler_polar_a_ ...
,tmp_XB_euler_azimu_b_ ...
,tmp_XB_euler_gamma_z_ ...
,tmp_XB_image_delta_x_ ...
,tmp_XB_image_delta_y_ ...
);
save(fname_XB_k_Y_mat,'a_k_Y_reco_');
end;%if (~exist(fname_XB_k_Y_mat,'file'));
%tmp_ = load(fname_XB_k_Y_mat); XB_a_k_Y_reco_ = tmp_.a_k_Y_reco_; clear tmp_; %<-- not used. ;
%%%%%%%%;
end;%if ( exist(XB_fname_mat,'file'));
%%%%%%%%;

%%%%%%%%;
end;%if ( exist(fname_XA_k_Y_mat,'file'));
%%%%%%%%;
%%%%%%%%;
end;%if ( exist(XA_fname_mat,'file'));
%%%%%%%%;

%%%%%%%%;
if ( exist(XA_fname_mat,'file') & ~exist(XA_fname_snapshot_jpg) );
%%%%%%%%;
if (verbose) disp(sprintf(' %% %s not found, creating',XA_fname_snapshot_jpg)); end;
tmp_XA_ = load(XA_fname_mat);
flag_found = 1; tmp_XA_align_a_k_Y_ = []; tmp_XA_a_k_Y_ = [];
%%%%;
tmp_XA_fname_pre = sprintf('%s_align_a_k_Y_',tmp_XA_.parameter.fname_pre);
tmp_XA_fname_pre = rootswitch(tmp_XA_fname_pre,string_root,'rangan');
tmp_XA_.parameter.fname_align_a_k_Y_pre = tmp_XA_fname_pre;
tmp_XA_fname_mat = sprintf('%s.mat',tmp_XA_fname_pre);
if (~exist(tmp_XA_fname_mat,'file'));
disp(sprintf(' %% Warning, %s not found in ampmut_wrap_wrap_4',tmp_XA_fname_mat));
flag_found = 0;
end;%if (~exist(tmp_XA_fname_mat,'file'));
if ( exist(tmp_XA_fname_mat,'file'));
tmp_XA_align_a_k_Y_ = load(tmp_XA_fname_mat);
end;%if ( exist(tmp_XA_fname_mat,'file'));
%%%%;
tmp_XA_fname_pre = sprintf('%s_a_k_Y_',tmp_XA_.parameter.fname_pre);
tmp_XA_fname_pre = rootswitch(tmp_XA_fname_pre,string_root,'rangan');
tmp_XA_.parameter.fname_a_k_Y_pre = tmp_XA_fname_pre;
tmp_XA_fname_mat = sprintf('%s.mat',tmp_XA_fname_pre);
if (~exist(tmp_XA_fname_mat,'file'));
disp(sprintf(' %% Warning, %s not found in ampmut_wrap_wrap_4',tmp_XA_fname_mat));
flag_found = 0;
end;%if (~exist(tmp_XA_fname_mat,'file'));
if ( exist(tmp_XA_fname_mat,'file'));
tmp_XA_a_k_Y_ = load(tmp_XA_fname_mat);
end;%if ( exist(tmp_XA_fname_mat,'file'));
%%%%;
if ( flag_found & ~isempty(tmp_XA_align_a_k_Y_) & ~isempty(tmp_XA_a_k_Y_) );
%%%%;
tmp_X_best = tmp_XA_align_a_k_Y_.X_best_(end);
tmp_flag_flip = tmp_XA_align_a_k_Y_.X_best_flag_flip_(end);
tmp_polar_a_best = tmp_XA_align_a_k_Y_.polar_a_best_(end);
tmp_azimu_b_best = tmp_XA_align_a_k_Y_.azimu_b_best_(end);
tmp_gamma_z_best = tmp_XA_align_a_k_Y_.gamma_z_best_(end);
tmp_delta_x_best = tmp_XA_align_a_k_Y_.delta_best__(1+0,end);
tmp_delta_y_best = tmp_XA_align_a_k_Y_.delta_best__(1+1,end);
tmp_a_k_Y_reco_ = tmp_XA_a_k_Y_.a_k_Y_reco_;
%%%%;
[ ... 
 tmp_b_k_Y_reco_ ...
] = ...
spharm_register_and_rotate_2( ...
 n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_true_ ...
,tmp_a_k_Y_reco_ ...
,0 ...
,tmp_X_best ...
,tmp_flag_flip ...
,tmp_polar_a_best ...
,tmp_azimu_b_best ...
,tmp_gamma_z_best ...
,tmp_delta_x_best ...
,tmp_delta_y_best ...
);
%%%%;
[ ... 
  tmp_b_x_u_reco_ ...
] = ...
convert_spharm_to_x_c_3( ...
 sample_sphere_k_eq_d ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,tmp_b_k_Y_reco_ ...
,half_diameter_x_c ...
,n_x_u_pack ...
);
%%%%;
figure(5);clf;figbig;
subplot(1,2,1);
isosurface_f_x_u_0(reshape(real(tmp_b_x_u_reco_),n_x_u_pack,n_x_u_pack,n_x_u_pack),[90,95,99]);
title(sprintf('tmp_b_x_u_: corr %0.4f',tmp_X_best),'Interpreter','none');
subplot(1,2,2);
isosurface_f_x_u_0(reshape(real(tmp_b_x_u_reco_),n_x_u_pack,n_x_u_pack,n_x_u_pack),[98.5]);
title(sprintf('tmp_b_x_u_: corr %0.4f',tmp_X_best),'Interpreter','none');
sgtitle(XA_fname_snapshot_jpg,'Interpreter','none');
disp(sprintf(' %% writing %s',XA_fname_snapshot_jpg));
print('-djpeg',XA_fname_snapshot_jpg);
close(gcf);
%%%%%%%%;
end;%if ( flag_found & ~isempty(tmp_XA_align_a_k_Y_) & ~isempty(tmp_XA_a_k_Y_) );
%%%%%%%%;
end;%if ( exist(XA_fname_mat,'file') & ~exist(XA_fname_snapshot_jpg) );
%%%%%%%%;

%%%%%%%%;
if ( exist(XB_fname_mat,'file') & ~exist(XB_fname_snapshot_jpg) );
%%%%%%%%;
if (verbose) disp(sprintf(' %% %s not found, creating',XB_fname_snapshot_jpg)); end;
tmp_XB_ = load(XB_fname_mat);
flag_found = 1; tmp_XB_align_a_k_Y_ = []; tmp_XB_a_k_Y_ = [];
%%%%;
tmp_XB_fname_pre = sprintf('%s_align_a_k_Y_',tmp_XB_.parameter.fname_pre);
tmp_XB_fname_pre = rootswitch(tmp_XB_fname_pre,string_root,'rangan');
tmp_XB_.parameter.fname_align_a_k_Y_pre = tmp_XB_fname_pre;
tmp_XB_fname_mat = sprintf('%s.mat',tmp_XB_fname_pre);
if (~exist(tmp_XB_fname_mat,'file'));
disp(sprintf(' %% Warning, %s not found in ampmut_wrap_wrap_4',tmp_XB_fname_mat));
flag_found = 0;
end;%if (~exist(tmp_XB_fname_mat,'file'));
if ( exist(tmp_XB_fname_mat,'file'));
tmp_XB_align_a_k_Y_ = load(tmp_XB_fname_mat);
end;%if ( exist(tmp_XB_fname_mat,'file'));
%%%%;
tmp_XB_fname_pre = sprintf('%s_a_k_Y_',tmp_XB_.parameter.fname_pre);
tmp_XB_fname_pre = rootswitch(tmp_XB_fname_pre,string_root,'rangan');
tmp_XB_.parameter.fname_a_k_Y_pre = tmp_XB_fname_pre;
tmp_XB_fname_mat = sprintf('%s.mat',tmp_XB_fname_pre);
if (~exist(tmp_XB_fname_mat,'file'));
disp(sprintf(' %% Warning, %s not found in ampmut_wrap_wrap_4',tmp_XB_fname_mat));
flag_found = 0;
end;%if (~exist(tmp_XB_fname_mat,'file'));
if ( exist(tmp_XB_fname_mat,'file'));
tmp_XB_a_k_Y_ = load(tmp_XB_fname_mat);
end;%if ( exist(tmp_XB_fname_mat,'file'));
%%%%;
if ( flag_found & ~isempty(tmp_XB_align_a_k_Y_) & ~isempty(tmp_XB_a_k_Y_) );
%%%%;
tmp_X_best = tmp_XB_align_a_k_Y_.X_best_(end);
tmp_flag_flip = tmp_XB_align_a_k_Y_.X_best_flag_flip_(end);
tmp_polar_a_best = tmp_XB_align_a_k_Y_.polar_a_best_(end);
tmp_azimu_b_best = tmp_XB_align_a_k_Y_.azimu_b_best_(end);
tmp_gamma_z_best = tmp_XB_align_a_k_Y_.gamma_z_best_(end);
tmp_delta_x_best = tmp_XB_align_a_k_Y_.delta_best__(1+0,end);
tmp_delta_y_best = tmp_XB_align_a_k_Y_.delta_best__(1+1,end);
tmp_a_k_Y_reco_ = tmp_XB_a_k_Y_.a_k_Y_reco_;
%%%%;
[ ... 
 tmp_b_k_Y_reco_ ...
] = ...
spharm_register_and_rotate_2( ...
 n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_true_ ...
,tmp_a_k_Y_reco_ ...
,0 ...
,tmp_X_best ...
,tmp_flag_flip ...
,tmp_polar_a_best ...
,tmp_azimu_b_best ...
,tmp_gamma_z_best ...
,tmp_delta_x_best ...
,tmp_delta_y_best ...
);
%%%%;
[ ... 
  tmp_b_x_u_reco_ ...
] = ...
convert_spharm_to_x_c_3( ...
 sample_sphere_k_eq_d ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,tmp_b_k_Y_reco_ ...
,half_diameter_x_c ...
,n_x_u_pack ...
);
%%%%;
figure(5);clf;figbig;
subplot(1,2,1);
isosurface_f_x_u_0(reshape(real(tmp_b_x_u_reco_),n_x_u_pack,n_x_u_pack,n_x_u_pack),[90,95,99]);
title(sprintf('tmp_b_x_u_: corr %0.4f',tmp_X_best),'Interpreter','none');
subplot(1,2,2);
isosurface_f_x_u_0(reshape(real(tmp_b_x_u_reco_),n_x_u_pack,n_x_u_pack,n_x_u_pack),[98.5]);
title(sprintf('tmp_b_x_u_: corr %0.4f',tmp_X_best),'Interpreter','none');
sgtitle(XB_fname_snapshot_jpg,'Interpreter','none');
disp(sprintf(' %% writing %s',XB_fname_snapshot_jpg));
print('-djpeg',XB_fname_snapshot_jpg);
close(gcf);
%%%%%%%%;
end;%if ( flag_found & ~isempty(tmp_XB_align_a_k_Y_) & ~isempty(tmp_XB_a_k_Y_) );
%%%%%%%%;
end;%if ( exist(XB_fname_mat,'file') & ~exist(XB_fname_snapshot_jpg) );
%%%%%%%%;

if (verbose); disp(sprintf(' %% [finished ampmut_wrap_wrap_4]')); end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if ~flag_skip_all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
