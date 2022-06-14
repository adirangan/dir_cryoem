%%%%%%%%;
% Attempt to plot density of viewing angles. ;
% Starting with ISWINCP. ;
%%%%%%%%;

test_pm_ISWINCP_10;

%%%%%%%%;
str_strategy_prefix = '';
date_diff_threshold = 0.25;
flag_force_create_mat=0;flag_force_create_tmp=0;
delta_sigma_use = delta_sigma;
dat_rseed_ = [0:2]; n_dat_rseed = numel(dat_rseed_);
dat_n_UX_rank_ = [2,4,6,8,10,12,14,16,18,20];  n_dat_n_UX_rank = numel(dat_n_UX_rank_);
delta_r_max_factor_ = [0.125,0.25,0.50,0.75]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%%%%%%%%;

%%%%%%%%;
% check to see if translations can be estimated (and appropriately updated) using only a few principal-modes. ;
%%%%%%%%;
a_k_Y_true_ = a_k_Y_quad_;
euler_polar_a_true_ = +euler_polar_a_true_(1+(0:n_M-1));
euler_azimu_b_true_ = +euler_azimu_b_true_(1+(0:n_M-1));
euler_gamma_z_true_ = +euler_gamma_z_true_(1+(0:n_M-1));
image_delta_x_true_ = +image_delta_x_true_(1+(0:n_M-1));
image_delta_y_true_ = +image_delta_y_true_(1+(0:n_M-1));
%%%%%%%%;
for ndat_rseed=0:n_dat_rseed-1;
dat_rseed = dat_rseed_(1+ndat_rseed);
for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
delta_r_max_upb = 2.0 * delta_sigma_use*sqrt(log(20^2)); %<-- allow large accumulated translations. ;
for ndat_n_UX_rank=0:n_dat_n_UX_rank-1;
dat_n_UX_rank = dat_n_UX_rank_(1+ndat_n_UX_rank);
for flag_alternate_MS_vs_SM = [0:1];
for flag_local_exclusion = [0];%for flag_local_exclusion = [0:1];
% test with: ;
% ndat_rseed=2; dat_rseed = dat_rseed_(1+ndat_rseed); ndelta_r_max_factor = 3; delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor); delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2)); delta_r_max_upb = 2.0 * delta_sigma_use*sqrt(log(20^2)); ndat_n_UX_rank = 9; dat_n_UX_rank = dat_n_UX_rank_(1+ndat_n_UX_rank); flag_alternate_MS_vs_SM=0; flag_local_exclusion=0;
parameter = struct('type','parameter');
parameter.rseed = dat_rseed;
parameter.delta_r_max = delta_r_max_use;
parameter.delta_r_upb = delta_r_max_upb;
parameter.dir_pm = dir_pm;
parameter.flag_alternate_MS_vs_SM = flag_alternate_MS_vs_SM;
parameter.flag_local_exclusion = flag_local_exclusion;
parameter.str_strategy_prefix = str_strategy_prefix;
parameter.flag_euler_polar_a_restrict = 0;
if (numel(unique(CTF_index_))>=n_M/8);
disp(sprintf(' %% Warning, setting flag_CTF_index_unused==1'));
parameter.flag_CTF_index_unused = 1;
end;%if (numel(unique(CTF_index_))>=n_M/8);
%%%%%%%%;
index_nCTF_from_nM_ = CTF_index_;
%%%%%%%%;
verbose=1;
if (verbose); disp(sprintf(' %% [entering ampmut_wrap_wrap_4]')); end;
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
%%%%%%%%;
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
%%%%%%%%;
for type_XX=0:1;
%%%%;
tmp_XX_str = 'Memp_d1'; if type_XX==1; tmp_XX_str = 'xcor_d0'; end;
%%%%%%%%;
n_w_max = n_w_max + mod(n_w_max,2); %<-- round up to nearest even number. ;
l_max_max = n_w_max/2 - 1; assert(l_max_max==max(l_max_));
n_w_ = n_w_max*ones(n_k_p_r,1);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);
n_lm_ = (l_max_+1).^2;
n_lm_max = max(n_lm_);
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
%%%%%%%%;
n_UX_rank = n_k_p_r-1; %<-- just to check dimensions. ;
%{
CTF_k_p_r__ = zeros(n_k_p_r,n_CTF);
for nCTF=0:n_CTF-1;
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
CTF_k_p_r__(1+nk_p_r,1+nCTF) = mean(CTF_k_p__(1+tmp_index_,1+nCTF));
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nCTF=0:n_CTF-1;
CTF_avg_k_p_ = mean(CTF_k_p__,2);
%imagesc_p(n_k_p_r,k_p_r_,n_w_,sum(n_w_),real(CTF_avg_k_p_(:)),[-1,+1],colormap_beach());
CTF_avg_k_p_r_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
CTF_avg_k_p_r_(1+nk_p_r) = mean(CTF_avg_k_p_(1+tmp_index_));
end;%for nk_p_r=0:n_k_p_r-1;
CTF_k_p_r_xcor__ = CTF_k_p_r__(:,1+index_nCTF_from_nM_(1+(0:n_M-1))) * transpose(CTF_k_p_r__(:,1+index_nCTF_from_nM_(1+(0:n_M-1)))) / n_M;
 %}
%%%%%%%%;
XX_fname_pre = sprintf('%s_mat/X_2d_%s_%st%.4dn%.2dr%d',dir_pm,tmp_XX_str,str_strategy,floor(1000*delta_r_max),dat_n_UX_rank,rseed);
XX_fname_mat = sprintf('%s.mat',XX_fname_pre);
if ( exist(XX_fname_mat,'file'));
%%%%%%%%;
if (verbose); disp(sprintf(' %% %s found, aligning',XX_fname_mat)); end;
tmp_XX_ = load(XX_fname_mat);
if (~isfield(tmp_XX_.parameter,'fname_pre')); tmp_XX_.parameter.fname_pre = XX_fname_pre; end; %<-- parameter_bookmark. ;
tmp_XX_fname_pre = sprintf('%s_align_a_CTF_avg_UX_Y_',tmp_XX_.parameter.fname_pre);
tmp_XX_fname_pre = rootswitch(tmp_XX_fname_pre,string_root,'rangan');
tmp_XX_.parameter.fname_align_a_CTF_avg_UX_Y_pre = tmp_XX_fname_pre;
tmp_XX_fname_mat = sprintf('%s.mat',tmp_XX_fname_pre);
%%%%%%%%;
a_CTF_avg_UX_Y_calc__ = tmp_XX_.a_CTF_avg_UX_Y__;
euler_polar_a_calc__ = tmp_XX_.euler_polar_a__;
euler_azimu_b_calc__ = tmp_XX_.euler_azimu_b__;
euler_gamma_z_calc__ = tmp_XX_.euler_gamma_z__;
image_delta_x_calc__ = tmp_XX_.image_delta_x_acc__ + tmp_XX_.image_delta_x_upd__;
image_delta_y_calc__ = tmp_XX_.image_delta_y_acc__ + tmp_XX_.image_delta_y_upd__;
n_iteration = size(euler_polar_a_calc__,2);
%%%%%%%%;
verbose=2;
if (verbose); disp(sprintf(' %% [entering ampmut_align_to_a_CTF_avg_UX_Y_0]')); end;
%%%%%%%%;
pm_n_UX_rank = dat_n_UX_rank;
pm_n_k_p_r = pm_n_UX_rank;
pm_k_p_r_ = ones(pm_n_k_p_r,1);
pm_k_p_r_max = 1;
pm_l_max_ = l_max_max*ones(pm_n_k_p_r,1);
pm_n_lm_ = (1+pm_l_max_).^2; pm_n_lm_sum = sum(pm_n_lm_);
pm_weight_3d_k_p_r_ = ones(pm_n_k_p_r,1);
pm_weight_2d_k_p_r_ = ones(pm_n_k_p_r,1);
%%%%%%%%;
if (~isempty(tmp_XX_.parameter.fname_align_a_CTF_avg_UX_Y_pre));
tmp_fname_mat = sprintf('%s.mat',tmp_XX_.parameter.fname_align_a_CTF_avg_UX_Y_pre);
load(tmp_fname_mat ...
     ,'X_best_','X_best_flag_flip_' ...
     ,'polar_a_best_' ...
     ,'azimu_b_best_' ...
     ,'gamma_z_best_' ...
     ,'delta_best__' ...
    );
end;%if (~isempty(tmp_XX_.parameter.fname_align_a_CTF_avg_UX_Y_pre));
%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%{
if (~isempty(tmp_XX_.parameter.fname_align_a_CTF_avg_UX_Y_pre));
tmp_fname_fig_jpg = sprintf('%s.jpg',tmp_XX_.parameter.fname_align_a_CTF_avg_UX_Y_pre);
if (~exist(tmp_fname_fig_jpg,'file'));
delta_sigma = 1.0 * std([image_delta_x_true_(1:n_M);image_delta_y_true_(1:n_M)],1); %<-- no reduction. ;
dscale = 2.5;
c2d_euler__ = colormap_polar_a_azimu_b_2d(+euler_polar_a_true_(1:n_M),+euler_azimu_b_true_(1:n_M),0.35);
markersize_euler = 25;
c2d_delta__ = colormap_gaussian_2d(+image_delta_x_true_(1:n_M),+image_delta_y_true_(1:n_M),dscale*delta_sigma,0.35);
markersize_delta = 25;
figure(3); clf; figbig; 
p_row = 4; p_col = 3*ceil((1+n_iteration)/p_row); np=0;
%%%%%%%%;
subplot(p_row,p_col,1+np+[0,1]); np=np+2;
hold on;
patch(2*pi*[+0;+1;+1;+0],1*pi*[+0;+0;+1;+1],0.65*[1,1,1],'EdgeColor','none');
scatter(0*pi+euler_azimu_b_true_(1:n_M),1*pi-euler_polar_a_true_(1:n_M),markersize_euler,c2d_euler__(1:n_M,:),'filled');
hold off;
axisnotick;
axis([0,2*pi,0,1*pi]); axis equal; xlabel('azimu_b'); ylabel('polar_a'); title('true view');
%%%%%%%%;
subplot(p_row,p_col,1+np); np=np+1;
hold on;
patch(dscale*delta_sigma*[-1;+1;+1;-1],dscale*delta_sigma*[-1;-1;+1;+1],'k');
scatter(+image_delta_x_true_(1:n_M),+image_delta_y_true_(1:n_M),markersize_delta,c2d_delta__(1:n_M,:),'filled');
hold off;
axisnotick;
axis(dscale*delta_sigma*[-1,+1,-1,+1]); axis square; xlabel('dx'); ylabel('dy'); title('true delta');
%%%%%%%%;
for niteration=0:n_iteration-1;
%%%%%%%%;
X_best = X_best_(1+niteration);
flag_flip = X_best_flag_flip_(1+niteration);
polar_a_best = polar_a_best_(1+niteration);
azimu_b_best = azimu_b_best_(1+niteration);
gamma_z_best = gamma_z_best_(1+niteration);
delta_best_ = delta_best__(:,1+niteration);
euler_best_ = [+azimu_b_best,+polar_a_best,+gamma_z_best];
R_best_ = euler_to_R(euler_best_);
%%%%%%%%;
tmp_euler_azimu_b_ = euler_azimu_b_calc__(:,1+niteration);
tmp_euler_polar_a_ = euler_polar_a_calc__(:,1+niteration);
tmp_k_c_0_ = sin(tmp_euler_polar_a_).*cos(tmp_euler_azimu_b_);
tmp_k_c_1_ = sin(tmp_euler_polar_a_).*sin(tmp_euler_azimu_b_);
tmp_k_c_2_ = cos(tmp_euler_polar_a_);
tmp__ = inv(R_best_) * transpose([tmp_k_c_0_(:) , tmp_k_c_1_(:) , tmp_k_c_2_(:)]);
tmp_k_c_0_(:) = transpose(tmp__(1+0,:));
tmp_k_c_1_(:) = transpose(tmp__(1+1,:));
tmp_k_c_2_(:) = transpose(tmp__(1+2,:));
tmp_k_c_01_ = sqrt(tmp_k_c_0_.^2 + tmp_k_c_1_.^2);
tmp_euler_azimu_b_ = periodize(atan2(tmp_k_c_1_,tmp_k_c_0_),0,2*pi);
tmp_euler_polar_a_ = atan2(tmp_k_c_01_,tmp_k_c_2_);
%%%%%%%%;
subplot(p_row,p_col,1+np+[0,1]); np=np+2; cla;
hold on;
patch(2*pi*[+0;+1;+1;+0],1*pi*[+0;+0;+1;+1],0.65*[1,1,1],'EdgeColor','none');
scatter(0*pi+tmp_euler_azimu_b_,1*pi-tmp_euler_polar_a_,markersize_euler,c2d_euler__(1:n_M,:),'filled');
hold off;
axisnotick;
axis([0,2*pi,0,1*pi]); axis equal; xlabel('azimu_b'); ylabel('polar_a');
title(sprintf('%d (flip %d) X %0.2f',niteration,flag_flip,X_best));
%%%%%%%%;
subplot(p_row,p_col,1+np); np=np+1;
hold on;
patch(dscale*delta_sigma*[-1;+1;+1;-1],dscale*delta_sigma*[-1;-1;+1;+1],'k');
scatter(+image_delta_x_calc__(1:n_M,1+niteration),+image_delta_y_calc__(1:n_M,1+niteration),markersize_delta,c2d_delta__(1:n_M,:),'filled');
hold off;
axisnotick;
axis(dscale*delta_sigma*[-1,+1,-1,+1]); axis square; xlabel('dx'); ylabel('dy'); 
title(sprintf('%d (delta)',niteration));
%%%%%%%%;
end;%for niteration=0:n_iteration-1;
%%%%%%%%;
sgtitle(sprintf('%s',tmp_fname_fig_jpg),'Interpreter','none');
print('-djpeg',tmp_fname_fig_jpg);
close(gcf);
end;%if (~exist(tmp_fname_fig_jpg,'file'));
end;%if (~isempty(tmp_XX_.parameter.fname_align_a_CTF_avg_UX_Y_pre));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
 %}
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (~isempty(tmp_XX_.parameter.fname_align_a_CTF_avg_UX_Y_pre));
tmp_fname_fig_jpg = sprintf('%s_sphere.jpg',tmp_XX_.parameter.fname_align_a_CTF_avg_UX_Y_pre);
if (~exist(tmp_fname_fig_jpg,'file'));
c2d_euler__ = colormap_polar_a_azimu_b_2d(+euler_polar_a_true_(1:n_M),+euler_azimu_b_true_(1:n_M),0.35);
markersize_euler = 25;
figure(3); clf; figbig; figbeach();
p_row = 4; p_col = 1*ceil((1+n_iteration)/p_row); np=0;
%%%%%%%%;
subplot(p_row,p_col,1+np+[0]); np=np+1;
tmp_x_ = cos(euler_azimu_b_true_(1:n_M)).*sin(euler_polar_a_true_(1:n_M));
tmp_y_ = sin(euler_azimu_b_true_(1:n_M)).*sin(euler_polar_a_true_(1:n_M));
tmp_z_ = cos(euler_polar_a_true_(1:n_M));
hold on;
scatter3(tmp_x_,tmp_y_,tmp_z_,markersize_euler,c2d_euler__,'filled','MarkerEdgeColor','k');
hold off;
axisnotick;
axis vis3d;
xlabel('x'); ylabel('x'); zlabel('x');
title('true view');
%%%%%%%%;
for niteration=0:n_iteration-1;
%%%%%%%%;
X_best = X_best_(1+niteration);
flag_flip = X_best_flag_flip_(1+niteration);
polar_a_best = polar_a_best_(1+niteration);
azimu_b_best = azimu_b_best_(1+niteration);
gamma_z_best = gamma_z_best_(1+niteration);
delta_best_ = delta_best__(:,1+niteration);
euler_best_ = [+azimu_b_best,+polar_a_best,+gamma_z_best];
R_best_ = euler_to_R(euler_best_);
%%%%%%%%;
tmp_euler_azimu_b_ = euler_azimu_b_calc__(:,1+niteration);
tmp_euler_polar_a_ = euler_polar_a_calc__(:,1+niteration);
tmp_k_c_0_ = sin(tmp_euler_polar_a_).*cos(tmp_euler_azimu_b_);
tmp_k_c_1_ = sin(tmp_euler_polar_a_).*sin(tmp_euler_azimu_b_);
tmp_k_c_2_ = cos(tmp_euler_polar_a_);
tmp__ = inv(R_best_) * transpose([tmp_k_c_0_(:) , tmp_k_c_1_(:) , tmp_k_c_2_(:)]);
tmp_k_c_0_(:) = transpose(tmp__(1+0,:));
tmp_k_c_1_(:) = transpose(tmp__(1+1,:));
tmp_k_c_2_(:) = transpose(tmp__(1+2,:));
tmp_k_c_01_ = sqrt(tmp_k_c_0_.^2 + tmp_k_c_1_.^2);
tmp_euler_azimu_b_ = periodize(atan2(tmp_k_c_1_,tmp_k_c_0_),0,2*pi);
tmp_euler_polar_a_ = atan2(tmp_k_c_01_,tmp_k_c_2_);
%%%%%%%%;
subplot(p_row,p_col,1+np+[0]); np=np+1; cla;
tmp_x_ = cos(tmp_euler_azimu_b_(1:n_M)).*sin(tmp_euler_polar_a_(1:n_M));
tmp_y_ = sin(tmp_euler_azimu_b_(1:n_M)).*sin(tmp_euler_polar_a_(1:n_M));
tmp_z_ = cos(tmp_euler_polar_a_(1:n_M));
hold on;
scatter3(tmp_x_,tmp_y_,tmp_z_,markersize_euler,c2d_euler__,'filled','MarkerEdgeColor','k');
hold off;
axisnotick;
axis vis3d;
xlabel('x'); ylabel('x'); zlabel('x');
title(sprintf('%d (flip %d) X %0.2f',niteration,flag_flip,X_best));
%%%%%%%%;
end;%for niteration=0:n_iteration-1;
%%%%%%%%;
sgtitle(sprintf('%s',tmp_fname_fig_jpg),'Interpreter','none');
%print('-djpeg',tmp_fname_fig_jpg);
%close(gcf);
end;%if (~exist(tmp_fname_fig_jpg,'file'));
end;%if (~isempty(tmp_XX_.parameter.fname_align_a_CTF_avg_UX_Y_pre));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
 %}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (~isempty(tmp_XX_.parameter.fname_align_a_CTF_avg_UX_Y_pre));
tmp_fname_fig_pre = sprintf('%s_h2',tmp_XX_.parameter.fname_align_a_CTF_avg_UX_Y_pre);
[tmp_flag_fig_skip,tmp_fname_fig_jpg] = open_fname_tmp(tmp_fname_fig_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp,'jpg');
if ~tmp_flag_fig_skip;
c2d_euler__ = colormap_polar_a_azimu_b_2d(+euler_polar_a_true_(1:n_M),+euler_azimu_b_true_(1:n_M),0.35);
markersize_euler = 25;
figure(3); clf; figbig; figbeach();
p_row = 4; p_col = 1*ceil((1+n_iteration)/p_row); np=0;
%%%%%%%%;
subplot(p_row,p_col,1+np+[0]); np=np+1;
hold on;
imagesc(log2(1+hist2d_0(euler_azimu_b_true_(1:n_M),euler_polar_a_true_(1:n_M),32,16,[0,2*pi],[0,pi])),[0,4]);
hold off;
axisnotick; axis image;
title('true view');
%%%%%%%%;
for niteration=0:n_iteration-1;
%%%%%%%%;
X_best = X_best_(1+niteration);
flag_flip = X_best_flag_flip_(1+niteration);
polar_a_best = polar_a_best_(1+niteration);
azimu_b_best = azimu_b_best_(1+niteration);
gamma_z_best = gamma_z_best_(1+niteration);
delta_best_ = delta_best__(:,1+niteration);
euler_best_ = [+azimu_b_best,+polar_a_best,+gamma_z_best];
R_best_ = euler_to_R(euler_best_);
%%%%%%%%;
tmp_euler_azimu_b_ = euler_azimu_b_calc__(:,1+niteration);
tmp_euler_polar_a_ = euler_polar_a_calc__(:,1+niteration);
tmp_k_c_0_ = sin(tmp_euler_polar_a_).*cos(tmp_euler_azimu_b_);
tmp_k_c_1_ = sin(tmp_euler_polar_a_).*sin(tmp_euler_azimu_b_);
tmp_k_c_2_ = cos(tmp_euler_polar_a_);
tmp__ = inv(R_best_) * transpose([tmp_k_c_0_(:) , tmp_k_c_1_(:) , tmp_k_c_2_(:)]);
tmp_k_c_0_(:) = transpose(tmp__(1+0,:));
tmp_k_c_1_(:) = transpose(tmp__(1+1,:));
tmp_k_c_2_(:) = transpose(tmp__(1+2,:));
tmp_k_c_01_ = sqrt(tmp_k_c_0_.^2 + tmp_k_c_1_.^2);
tmp_euler_azimu_b_ = periodize(atan2(tmp_k_c_1_,tmp_k_c_0_),0,2*pi);
tmp_euler_polar_a_ = atan2(tmp_k_c_01_,tmp_k_c_2_);
%%%%%%%%;
subplot(p_row,p_col,1+np+[0]); np=np+1; cla;
hold on;
imagesc(log2(1+hist2d_0(tmp_euler_azimu_b_(1:n_M),tmp_euler_polar_a_(1:n_M),32,16,[0,2*pi],[0,pi])),[0,4]);
hold off;
axisnotick; axis image;
title(sprintf('%d (flip %d) X %0.2f',niteration,flag_flip,X_best));
%%%%%%%%;
end;%for niteration=0:n_iteration-1;
%%%%%%%%;
sgtitle(sprintf('%s',tmp_fname_fig_jpg),'Interpreter','none');
print('-djpeg',tmp_fname_fig_jpg);
close(gcf);
close_fname_tmp(tmp_fname_fig_pre);
end;%if ~tmp_flag_fig_skip;
end;%if (~isempty(tmp_XX_.parameter.fname_align_a_CTF_avg_UX_Y_pre));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if ( exist(XX_fname_mat,'file'));
%%%%;
end;%for type_XX=0:1;
%%%%;
end;%for flag_local_exclusion = [0];%for flag_local_exclusion = [0:1];
end;%for flag_alternate_MS_vs_SM = [0:1];
end;%for ndat_n_UX_rank=0:n_dat_n_UX_rank-1;
end;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
end;%for ndat_rseed=0:n_dat_rseed-1;
%%%%%%%%;

disp('returning'); return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
