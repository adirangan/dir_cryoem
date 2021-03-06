% trying to set up wignert_ODE using lsq-interpolation to handle multiple translations. ;
clear; setup_OptiPlex;

%%%%%%%%;
% First define integral of <f,f>. ;
%%%%%%%%;
h_ = @(kd) 4*pi*( sin(kd) - (kd).*cos(kd) ) ./ kd.^3 ; % calculates <f_j,f_k>, normalized so that <f,f> = 4*pi/3, rather than 4*pi/3*R_max^3 ;
dh_ = @(kd) 12*pi*( (kd.^2/3 - 1) .* sin(kd) + (kd).*cos(kd) ) ./ kd.^4 ;

%%%%%%%%;
% Now set up and test k-quadrature on sphere. ;
%%%%%%%%;
verbose=0; k_p_r_max = 9.0d0; k_eq_d = 1.0/32*k_p_r_max*0.95; TorL = 'L';
[n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_k_p_r_,k_c_0_all_,k_c_1_all_,k_c_2_all_] = sample_sphere_6(verbose,k_p_r_max,k_eq_d,'L') ;
delta_a_c_ = [+0.15;-0.25;+0.35];
delta_a_p_r = sqrt(delta_a_c_(1+0)^2 + delta_a_c_(1+1)^2 + delta_a_c_(1+2)^2);
%disp(sprintf(' %% 2*pi*k_p_r_max*delta_a_p_r = %0.2f',2*pi*k_p_r_max*delta_a_p_r));
a_k_all_form_ = exp(+i*2*pi*(k_c_0_all_*delta_a_c_(1+0) + k_c_1_all_*delta_a_c_(1+1) + k_c_2_all_*delta_a_c_(1+2)));
I_quad = sum(a_k_all_form_.*weight_k_all_);
I_form = h_(2*pi*k_p_r_max*fnorm(delta_a_c_))*k_p_r_max^3;
disp(sprintf(' %% I_form vs I_quad %0.16f',fnorm(I_form-I_quad)/fnorm(I_form)));

%%%%%%%%;
% Now select a few (suboptimal) nodes for lsq-interpolation. ;
%%%%%%%%;
N_pixel = 1.0*sqrt(2)*2*pi; z_target = N_pixel/sqrt(2); delta_p_r_max = z_target/(2*pi*k_p_r_max);
% First define coarse nodes for lsq: ;
delta_eq_d = 1.0/2.0*delta_p_r_max; TorL = 'L';
[n_delta_all,n_delta_all_csum_,delta_p_r_all_,delta_p_azimu_b_all_,delta_p_polar_a_all_,weight_delta_all_,weight_shell_delta_,n_delta_p_r,delta_p_r_,weight_delta_p_r_,delta_c_0_all_,delta_c_1_all_,delta_c_2_all_] = sample_sphere_6(verbose,delta_p_r_max,delta_eq_d,'L') ;
n_delta_node = n_delta_all;
n_delta_node_csum_ = [n_delta_all_csum_;n_delta_node];
n_delta_node_p_r = n_delta_p_r;
delta_node_p_r_ = delta_p_r_;
weight_delta_node_all_ = weight_delta_all_;
delta_node_ = transpose([delta_c_0_all_,delta_c_1_all_,delta_c_2_all_]);
%delta_node_ = delta_node_ + 0.01*randn(3,n_delta_node)/sqrt(3)*delta_p_r_max ;
%delta_node_ = randn(3,n_delta_node)/sqrt(3); %delta_node_ = randn(3,n_delta_node)/sqrt(3)*delta_p_r_max;
disp(sprintf(' %% N_pixel %0.2f z_target (i.e., kd) %0.2f*2*pi n_delta_node: %d',N_pixel,z_target/(2*pi),n_delta_node));
flag_refine=0;
delta_eq_d = 1.0/32*delta_p_r_max; TorL = 'L';
[n_delta_all,n_delta_all_csum_,delta_p_r_all_,delta_p_azimu_b_all_,delta_p_polar_a_all_,weight_delta_all_,weight_shell_delta_,n_delta_p_r,delta_p_r_,weight_delta_p_r_,delta_c_0_all_,delta_c_1_all_,delta_c_2_all_] = sample_sphere_6(verbose,delta_p_r_max,delta_eq_d,'L') ;
svd_tolerance = 1e-6;
[~,E3_quad] = transf_3d_gradient(0,[delta_node_(:);0],svd_tolerance,k_p_r_max,delta_p_r_max,n_delta_all,delta_c_0_all_,delta_c_1_all_,delta_c_2_all_,weight_delta_all_);
if flag_refine;
% Then define finer samples for quadrature: ;
n_iteration = 8; dt = 0.01;
for niteration=1:n_iteration;
[tmp_gradient_,E3_quad] = transf_3d_gradient(0,[delta_node_(:);0],svd_tolerance,k_p_r_max,delta_p_r_max,n_delta_all,delta_c_0_all_,delta_c_1_all_,delta_c_2_all_,weight_delta_all_); tmp_gradient_ = reshape(tmp_gradient_(1:3*n_delta_node),3,n_delta_node);
disp(sprintf(' %% niteration %d/%d E3_quad: %0.16f',niteration,n_iteration,E3_quad));
delta_node_ = delta_node_ + dt*tmp_gradient_;
end;%for niteration=1:n_iteration;
end;%if flag_refine;
disp(sprintf(' %% using n_delta_node %d with E3_quad: %0.16f',n_delta_node,E3_quad));
% Then define moderate sampling for displacement-grid: ;
delta_eq_d = 1.0/4*delta_p_r_max; TorL = 'L';
[n_delta_all,n_delta_all_csum_,delta_p_r_all_,delta_p_azimu_b_all_,delta_p_polar_a_all_,weight_delta_all_,weight_shell_delta_,n_delta_p_r,delta_p_r_,weight_delta_p_r_,delta_c_0_all_,delta_c_1_all_,delta_c_2_all_] = sample_sphere_6(verbose,delta_p_r_max,delta_eq_d,'L') ;
delta_all_ = transpose([delta_c_0_all_,delta_c_1_all_,delta_c_2_all_]);
disp(sprintf(' %% displacement-grid has n_delta_all %d',n_delta_all));
% Now construct delta_lsq_interpolate_. ;
svd_tolerance = 1e-6;
[~,~,delta_lsq_interpolate_] = transf_3d_gradient(0,[delta_node_(:);0],svd_tolerance,k_p_r_max,delta_p_r_max,n_delta_all,delta_c_0_all_,delta_c_1_all_,delta_c_2_all_,weight_delta_all_); 

%%%%%%%%;
% Now try to reproduce a full array of shifts of simple plane-wave on a sphere. ;
%%%%%%%%;
delta_a_c_ = [+0.15;-0.25;+0.35];
delta_a_p_r = sqrt(delta_a_c_(1+0)^2 + delta_a_c_(1+1)^2 + delta_a_c_(1+2)^2);
%disp(sprintf(' %% 2*pi*k_p_r_max*delta_a_p_r = %0.2f',2*pi*k_p_r_max*delta_a_p_r));
a_k_all_form_ = exp(+i*2*pi*(k_c_0_all_*delta_a_c_(1+0) + k_c_1_all_*delta_a_c_(1+1) + k_c_2_all_*delta_a_c_(1+2)));
l_max_ = 1+ceil(2*pi*k_p_r_); n_lm_ = (l_max_+1).^2; n_lm_sum = sum(n_lm_);
l_max_max = max(l_max_); dWtdkd__l_max_max = 2*l_max_max;
dWtdkd__ = dwignertdkd__(dWtdkd__l_max_max);
tmp_t = tic;
[a_k_Y_quad_] = convert_k_p_to_spharm_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_k_p_r_,l_max_,a_k_all_form_);
tmp_t = toc(tmp_t); tmp_g = n_lm_sum/tmp_t/1e6;
disp(sprintf(' %% use convert_k_p_to_spharm_1 to form a_k_Y_quad_... time %0.6fs (%0.6f Mnumps)',tmp_t,tmp_g));
%%%%%%%%;

%%%%%%%%;
% Now translate by each of the delta_node_;
%%%%%%%%;
b_k_Y_node__ = zeros(n_lm_sum,n_delta_node);
sum_t1 = 0; sum_t2 = 0; sum_t3 = 0; sum_t4 = 0; sum_t5 = 0;
tmp_t = tic; dWtdkd__ = dwignertdkd__(dWtdkd__l_max_max);
tmp_t = toc(tmp_t); sum_t1 = sum_t1 + tmp_t; %disp(sprintf(' %% sphere: create dWtdkd__: time %0.6fs',tmp_t));
for ndelta_node_p_r=1:n_delta_node_p_r;
delta_node_p_r = delta_node_p_r_(ndelta_node_p_r);
tmp_t = tic; Wt___ = expm_dwignertdkd__(dWtdkd__,n_k_p_r,k_p_r_,l_max_,delta_node_p_r);
tmp_t = toc(tmp_t); sum_t2 = sum_t2 + tmp_t;  %disp(sprintf(' %% sphere: create Wt___: time %0.6fs',tmp_t));
tmp_t = tic; Wt_ = wignert_ODE_0(dWtdkd__,Wt___,n_k_p_r,k_p_r_,l_max_,delta_node_p_r);
tmp_t = toc(tmp_t); sum_t3 = sum_t3 + tmp_t;  %disp(sprintf(' %% sphere: create Wt_: time %0.6fs',tmp_t));
tmp_n_delta_node = numel(1+n_delta_node_csum_(ndelta_node_p_r):n_delta_node_csum_(1+ndelta_node_p_r));
for ndelta_node=1+n_delta_node_csum_(ndelta_node_p_r):n_delta_node_csum_(1+ndelta_node_p_r);
if (mod(ndelta_node,100)==0); disp(sprintf(' %% ndelta_node %d/%d',ndelta_node,n_delta_node)); end;
delta_z_c_ = transpose(delta_node_(:,ndelta_node));
delta_z_p_r = sqrt(delta_z_c_(1+0).^2 + delta_z_c_(1+1).^2 + delta_z_c_(1+2).^2);
delta_z_p_01 = sqrt(delta_z_c_(1+0).^2 + delta_z_c_(1+1).^2);
delta_z_p_azimu_b = atan2(delta_z_c_(1+1),delta_z_c_(1+0));
delta_z_p_polar_a = atan2(delta_z_p_01,delta_z_c_(1+2));
delta_z_p_euler_pos_ = [0,+delta_z_p_polar_a,+delta_z_p_azimu_b];
delta_z_p_euler_neg_ = [-delta_z_p_azimu_b,-delta_z_p_polar_a,0];
delta_z_c_ = [cos(delta_z_p_azimu_b)*sin(delta_z_p_polar_a);sin(delta_z_p_azimu_b)*sin(delta_z_p_polar_a);cos(delta_z_p_polar_a)]*delta_z_p_r;
%%%%%%%%;
tmp_t = tic;
W_beta_neg__ = wignerd_b(l_max_max,delta_z_p_euler_neg_(1+1));
W_beta_pos__ = wignerd_b(l_max_max,delta_z_p_euler_pos_(1+1));
tmp_t = toc(tmp_t);  sum_t4 = sum_t4 + tmp_t; %disp(sprintf(' %% sphere: create Wd_: time %0.6fs',tmp_t));
tmp_t = tic; 
tmp_Y_form_ = a_k_Y_quad_;
tmp_Y_form_ = rotate_spharm_to_spharm_2(0,W_beta_neg__,n_k_p_r,k_p_r_,l_max_,tmp_Y_form_,delta_z_p_euler_neg_);
tmp_Y_form_ = Wt_*tmp_Y_form_; 
tmp_Y_form_ = rotate_spharm_to_spharm_2(0,W_beta_pos__,n_k_p_r,k_p_r_,l_max_,tmp_Y_form_,delta_z_p_euler_pos_);
b_k_Y_form_ = tmp_Y_form_;
tmp_t = toc(tmp_t); sum_t5 = sum_t5 + tmp_t; 
flag_disp=1;
if flag_disp;
delta_b_c_ = delta_a_c_ - delta_z_c_; %<-- note sign of translation. ;
b_k_all_form_ = exp(+i*2*pi*(k_c_0_all_*delta_b_c_(1+0) + k_c_1_all_*delta_b_c_(1+1) + k_c_2_all_*delta_b_c_(1+2)));
tmp_t = tic;
[b_k_Y_quad_] = convert_k_p_to_spharm_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_k_p_r_,l_max_,b_k_all_form_);
tmp_t = toc(tmp_t); tmp_g = n_lm_sum/tmp_t/1e6;
disp(sprintf(' %% ndelta_node_p_r %d/%d ndelta_node %d/%d: convert_k_p_to_spharm_1: time %0.6fs: error vs b %0.16f <-- vs a %0.16f',ndelta_node_p_r,n_delta_node_p_r,ndelta_node,tmp_n_delta_node,tmp_t,fnorm(b_k_Y_quad_-b_k_Y_form_)/fnorm(b_k_Y_quad_),fnorm(b_k_Y_quad_-a_k_Y_quad_)/fnorm(b_k_Y_quad_)));
end;%if flag_disp;
b_k_Y_node__(:,ndelta_node) = b_k_Y_form_;
end;%for ndelta_node=1+n_delta_node_csum_(ndelta_node_p_r):n_delta_node_csum_(1+ndelta_node_p_r);
end;%for ndelta_node_p_r=1:n_delta_node_p_r;
disp(sprintf(' %% timing: sum_t1 %0.2fs sum_t2 %0.2fs sum_t3 %0.2fs sum_t4 %0.2fs sum_t5 %0.2fs',sum_t1,sum_t2,sum_t3,sum_t4,sum_t5));
b_k_Y_lsqi__ = b_k_Y_node__*delta_lsq_interpolate_;
flag_plot=0; if flag_plot; imagesc(abs(corr(b_k_Y_lsqi__)),[0,+1]); colormap(colormap_beach());axis image; end;%if flag_plot;
b_k_Y_lint__ = zeros(n_lm_sum,n_delta_all);
for nlm_sum=1:n_lm_sum;
if (mod(nlm_sum,1000)==0); disp(sprintf(' %% nlm_sum %d/%d',nlm_sum,n_lm_sum)); end;
fp_b_k_Y_lint__ = scatteredInterpolant(transpose(delta_node_(1,:)),transpose(delta_node_(2,:)),transpose(delta_node_(3,:)),transpose(b_k_Y_node__(nlm_sum,:)));
b_k_Y_lint__(nlm_sum,:) = transpose(fp_b_k_Y_lint__(transpose(delta_all_(1,:)),transpose(delta_all_(2,:)),transpose(delta_all_(3,:))));
end;%for nlm_sum=1:n_lm_sum;

flag_test=0;
if flag_test;
%%%%%%%%;
% Now test out b_k_Y_lsqi__ for one of the sampled delta_. ;
%%%%%%%%;
n_iteration = 16;
test_dist_ = zeros(n_iteration,1);
test_lsqi_error_ = zeros(n_iteration,1);
test_lint_error_ = zeros(n_iteration,1);
for niteration=1:n_iteration;
ndelta_all = 1+max(0,min(n_delta_all-1,floor(n_delta_all*rand())));
tmp_delta_all_ = delta_all_(:,ndelta_all);
tmp_dd_ = delta_node_ - repmat(tmp_delta_all_,1,n_delta_node); tmp_hypot_dd_ = sqrt(sum(tmp_dd_.^2,1));
[tmp_hypot_dd,tmp_ij] = min(tmp_hypot_dd_); 
disp(sprintf(' %% ndelta_all %d, closest to node %d, distance %0.6f',ndelta_all,tmp_ij,tmp_hypot_dd));
test_dist_(niteration) = tmp_hypot_dd;
delta_z_c_ = transpose(tmp_delta_all_);
delta_z_p_r = sqrt(delta_z_c_(1+0).^2 + delta_z_c_(1+1).^2 + delta_z_c_(1+2).^2);
delta_z_p_01 = sqrt(delta_z_c_(1+0).^2 + delta_z_c_(1+1).^2);
delta_z_p_azimu_b = atan2(delta_z_c_(1+1),delta_z_c_(1+0));
delta_z_p_polar_a = atan2(delta_z_p_01,delta_z_c_(1+2));
delta_z_p_euler_pos_ = [0,+delta_z_p_polar_a,+delta_z_p_azimu_b];
delta_z_p_euler_neg_ = [-delta_z_p_azimu_b,-delta_z_p_polar_a,0];
delta_z_c_ = [cos(delta_z_p_azimu_b)*sin(delta_z_p_polar_a);sin(delta_z_p_azimu_b)*sin(delta_z_p_polar_a);cos(delta_z_p_polar_a)]*delta_z_p_r;
delta_b_c_ = delta_a_c_ - delta_z_c_; %<-- note sign of translation. ;
b_k_all_form_ = exp(+i*2*pi*(k_c_0_all_*delta_b_c_(1+0) + k_c_1_all_*delta_b_c_(1+1) + k_c_2_all_*delta_b_c_(1+2)));
tmp_t = tic;
[b_k_Y_quad_] = convert_k_p_to_spharm_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_k_p_r_,l_max_,b_k_all_form_);
tmp_t = toc(tmp_t); tmp_g = n_lm_sum/tmp_t/1e6;
%disp(sprintf(' %% sphere: use convert_k_p_to_spharm_1 to form a_k_Y_quad_... time %0.6fs (%0.6f Mnumps)',tmp_t,tmp_g));
disp(sprintf(' %% b_k_Y_quad_ vs b_k_Y_lsqi__: %0.16f',fnorm(b_k_Y_quad_ - b_k_Y_lsqi__(:,ndelta_all))/fnorm(b_k_Y_quad_)));
test_lsqi_error_(niteration) = fnorm(b_k_Y_quad_ - b_k_Y_lsqi__(:,ndelta_all))/fnorm(b_k_Y_quad_);
disp(sprintf(' %% b_k_Y_quad_ vs b_k_Y_lint__: %0.16f',fnorm(b_k_Y_quad_ - b_k_Y_lint__(:,ndelta_all))/fnorm(b_k_Y_quad_)));
test_lint_error_(niteration) = fnorm(b_k_Y_quad_ - b_k_Y_lint__(:,ndelta_all))/fnorm(b_k_Y_quad_);
end;%for niteration=1:n_iteration;
figure();clf; hold on;
scatter(test_dist_,test_lsqi_error_,15,'ro','filled'); 
scatter(test_dist_,test_lint_error_,15,'bo','filled'); 
legend({'lsqi','lint'},'Location','NorthWest');
xlabel('distance to nearest node'); ylabel('error');
%%%%%%%%;
end;%if flag_test;

%%%%%%%%;
% Now calculate the dot-product between b_k_Y_lsqi__ and c_k_Y_quad_. ;
%%%%%%%%;
delta_c_c_ = [-0.65;+0.45;-0.15];
delta_c_p_r = sqrt(delta_c_c_(1+0)^2 + delta_c_c_(1+1)^2 + delta_c_c_(1+2)^2);
%disp(sprintf(' %% 2*pi*k_p_r_max*delta_c_p_r = %0.2f',2*pi*k_p_r_max*delta_c_p_r));
c_k_all_form_ = exp(+i*2*pi*(k_c_0_all_*delta_c_c_(1+0) + k_c_1_all_*delta_c_c_(1+1) + k_c_2_all_*delta_c_c_(1+2)));
tmp_t = tic;
[c_k_Y_quad_] = convert_k_p_to_spharm_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_k_p_r_,l_max_,c_k_all_form_);
tmp_t = toc(tmp_t); tmp_g = n_lm_sum/tmp_t/1e6;
disp(sprintf(' %% use convert_k_p_to_spharm_1 to form c_k_Y_quad_... time %0.6fs (%0.6f Mnumps)',tmp_t,tmp_g));
%%%%%%%%;
weight_Y_row_ = zeros(n_lm_sum,1);
weight_Y_col_ = zeros(n_lm_sum,1);
weight_Y_val_ = zeros(n_lm_sum,1);
na=0;
for nk_p_r=1:n_k_p_r;
tmp_ij_ = na + (1:n_lm_(nk_p_r));
weight_Y_row_(tmp_ij_) = tmp_ij_;
weight_Y_col_(tmp_ij_) = tmp_ij_;
weight_Y_val_(tmp_ij_) = weight_k_p_r_(nk_p_r);
na=na+n_lm_(nk_p_r);
end;%for nk_p_r=1:n_k_p_r;
weight_Y_ = sparse(weight_Y_row_,weight_Y_col_,weight_Y_val_,n_lm_sum,n_lm_sum);
X_lsqi_ = ctranspose(b_k_Y_lsqi__)*(weight_Y_*c_k_Y_quad_);
X_lint_ = ctranspose(b_k_Y_lint__)*(weight_Y_*c_k_Y_quad_);
X_form_ = zeros(n_delta_all,1);
for ndelta_all=1:n_delta_all;
delta_z_c_ = delta_all_(:,ndelta_all);
delta_d_c_ = delta_c_c_ - (delta_a_c_ - delta_z_c_);
delta_d_p_r = sqrt(delta_d_c_(1+0)^2 + delta_d_c_(1+1)^2 + delta_d_c_(1+2)^2);
X_form_(ndelta_all) = h_(2*pi*k_p_r_max*delta_d_p_r)*k_p_r_max^3;
end;%for ndelta_all=1:n_delta_all;
X_lsqi_diff_ = X_form_-X_lsqi_;
X_lsqi_l2error = sum(abs(X_lsqi_diff_).^2 .* weight_delta_all_);
disp(sprintf(' %% X_form_ vs X_lsqi_ %0.16f',X_lsqi_l2error));
X_lint_diff_ = X_form_-X_lint_;
X_lint_l2error = sum(abs(X_lint_diff_).^2 .* weight_delta_all_);
disp(sprintf(' %% X_form_ vs X_lint_ %0.16f',X_lint_l2error));
%%%%%%%%;
[tmp_idx_,tmp_dist_] = knnsearch(transpose(delta_node_),transpose(delta_all_));
flag_plot=1;
if flag_plot;
subplot(1,2,1);
hold on;
scatter(tmp_dist_,abs(X_lsqi_diff_),15,'ro','filled'); 
scatter(tmp_dist_,abs(X_lint_diff_),15,'bo','filled'); 
legend({'lsqi','lint'},'Location','NorthWest');
xlabel('distance to nearest node'); ylabel('abs(error)'); title('error vs distance');
subplot(1,2,2);
hold on;
scatter(abs(X_lsqi_diff_),abs(X_lint_diff_),15,'bo','filled'); 
plot(0:max(abs(X_lint_diff_)),0:max(abs(X_lint_diff_)),'k-');
axis equal;
xlabel('lsqi error');ylabel('lint error'); title('lsqi vs lint');
figbig;
end;%if flag_plot;

%%%%%%%%;
% Now calculate the dot-product between the b_k_Y_lsqi__ and c_k_Y_quad_ ;
% for multiple rotations across a single polar_a for uniform azimu_b and gamma_z. ;
%%%%%%%%;
polar_a = -1*pi/3;
l_max_max = l_max_(1+n_k_p_r-1);
n_m_max = 2*l_max_max+1;
Y_node___ = zeros(n_m_max,n_m_max,n_delta_node);
%%%%%%%%;
for ndelta_node=1:n_delta_node;
if (mod(ndelta_node,100)==0); disp(sprintf(' %% ndelta_node %d/%d',ndelta_node,n_delta_node)); end;
[Y_node_] = register_spharm_to_spharm_single_beta_2(verbose,n_k_p_r,k_p_r_,weight_k_p_r_,l_max_,c_k_Y_quad_,b_k_Y_node__(:,ndelta_node),polar_a,0,[],0,[],[],[],[],[]);
Y_node___(:,:,ndelta_node) = Y_node_;
end;%for ndelta_node=1:n_delta_node;
%%%%%%%%;
Y_lsqi___ = reshape(reshape(Y_node___,[n_m_max*n_m_max,n_delta_node])*delta_lsq_interpolate_,[n_m_max,n_m_max,n_delta_all]);
%%%%%%%%;
X_lsqi___ = zeros(n_m_max,n_m_max,n_delta_all);
X_lint___ = zeros(n_m_max,n_m_max,n_delta_all);
X_form___ = zeros(n_m_max,n_m_max,n_delta_all);
for ndelta_all=1:n_delta_all;
if (mod(ndelta_all,100)==0); disp(sprintf(' %% ndelta_all %d/%d',ndelta_all,n_delta_all)); end;
[X_lsqi_] = register_spharm_to_spharm_single_beta_2(verbose,n_k_p_r,k_p_r_,weight_k_p_r_,l_max_,c_k_Y_quad_,b_k_Y_lsqi__(:,ndelta_all),polar_a,0,[],0,[],[],[],[],[]);
[X_lint_] = register_spharm_to_spharm_single_beta_2(verbose,n_k_p_r,k_p_r_,weight_k_p_r_,l_max_,c_k_Y_quad_,b_k_Y_lint__(:,ndelta_all),polar_a,0,[],0,[],[],[],[],[]);
%%%%;
X_form_ = zeros(n_m_max,n_m_max);
for nazimu_b=0:n_m_max-1;
for ngamma_z=0:n_m_max-1;
azimu_b = 2*pi*nazimu_b/n_m_max;
gamma_z = 2*pi*ngamma_z/n_m_max;
ef_ = [-gamma_z;-polar_a;-azimu_b]; %<-- vector of euler angles: [azimu_b, polar_a, gamma_z];
R_ = euler_to_R(+ef_);
delta_z_c_ = R_*(delta_a_c_ - delta_all_(:,ndelta_all));
delta_d_c_ = delta_c_c_ - delta_z_c_;
delta_d_p_r = sqrt(delta_d_c_(1+0)^2 + delta_d_c_(1+1)^2 + delta_d_c_(1+2)^2);
X_form_(1+nazimu_b,1+ngamma_z) = h_(2*pi*k_p_r_max*delta_d_p_r)*k_p_r_max^3;
end;%for ngamma_z=0:n_m_max-1;
end;%for nazimu_b=0:n_m_max-1;
%%%%;
X_lsqi___(:,:,ndelta_all) = X_lsqi_;
X_lint___(:,:,ndelta_all) = X_lint_;
X_form___(:,:,ndelta_all) = X_form_;
end;%for ndelta_all=1:n_delta_all;
%%%%%%%%;
disp(sprintf(' %% Y_lsqi___ vs X_lsqi___: %0.16f',fnorm(Y_lsqi___-X_lsqi___)/fnorm(Y_lsqi___)));
disp(sprintf(' %% X_form___ vs X_lsqi___: %0.16f',fnorm(X_form___-X_lsqi___)/fnorm(X_form___)));
disp(sprintf(' %% X_form___ vs X_lint___: %0.16f',fnorm(X_form___-X_lint___)/fnorm(X_form___)));

l10X_lsqi_diff_ = log10( abs(X_form___(:) - X_lsqi___(:))./abs(X_form___(:)) );
l10X_lint_diff_ = log10( abs(X_form___(:) - X_lint___(:))./abs(X_form___(:)) );
flag_plot=0;
if flag_plot;
n_h = 24*7+1;
imagesc(log2(1+hist2d_0(l10X_lsqi_diff_,l10X_lint_diff_,n_h,n_h,[-5,1],[-5,1]))); 
colormap(colormap_beach()); 
set(gca,'XTick',round(linspace(1,n_h,7)),'XTickLabel',linspace(-5,1,7));
set(gca,'YTick',round(linspace(1,n_h,7)),'YTickLabel',linspace(-5,1,7));
xlabel('log10(lsqi rel error)'); ylabel('log10(lint rel error)');
title('scatterplot of lsqi vs lint');
axis image; 
figbig;
end;%if flag_plot;
