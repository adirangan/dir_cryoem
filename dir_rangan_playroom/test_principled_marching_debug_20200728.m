% Here we try and determine the difference between test_principled_marching_trpv1_8.m and test_principled_marching_trpv1_15.m. ;

% First load preliminary data. ;
clear; setup;

str_CTF_factor = 'original';

dir_trunk = sprintf('/data/rangan/dir_cryoem/dir_rangan_playroom/dir_principled_marching');

h2d_ = @(kd) (2*pi)^2*(besselj(0,kd) + besselj(2,kd)); % calculates <f_j,f_k>, normalized so that <f,f> = (2*pi)^2;
dh2d_ = @(kd) (2*pi)^2*0.5*(besselj(-1,kd) - besselj(+3,kd));
h3d_ = @(kd) (2*pi)^3*3*( sin(kd) - (kd).*cos(kd) ) ./ kd.^3 ; % calculates <f_j,f_k>, normalized so that <f,f> = (2*pi)^3;
dh3d_ = @(kd) (2*pi)^3*9*( (kd.^2/3 - 1) .* sin(kd) + (kd).*cos(kd) ) ./ kd.^4 ;

%%%%%%%%;
% First load trpv1 molecule on x_u grid. ;
%%%%%%%%;
dir_data = '/data/rangan/dir_cryoem/dir_trpv1/data_nosym';
fname_dims = sprintf('%s/dims',dir_data);
tmp_ = textread(fname_dims); n_x_u = tmp_(1); n_image = tmp_(2); clear tmp_;
fname_density = sprintf('%s/density_clean',dir_data);
a_x_u_load_ = textread(fname_density); a_x_u_load_ = reshape(a_x_u_load_,n_x_u,n_x_u,n_x_u);
%%%%%%%%;
half_diameter_x_c = 1.0d0;
diameter_x_c = 2.0d0*half_diameter_x_c;
x_p_r_max = 1.0;
x_u_res = 64;
n_pack = n_x_u/x_u_res;
pack_row_ij_ = zeros(x_u_res,1);
pack_col_ij_ = zeros(x_u_res,1);
pack_val_ij_ = zeros(x_u_res,1);
na=0;
for nx_u=0:n_x_u-1;
pack_row_ij_(1+na) = 1+nx_u;
pack_col_ij_(1+na) = 1+floor(nx_u/n_pack);
pack_val_ij_(1+na) = 1/n_pack;
na=na+1;
end;%for nx_u=0:n_x_u-1;
x_u_pack_ = sparse(pack_row_ij_,pack_col_ij_,pack_val_ij_,n_x_u,x_u_res);
a_x_u_pack_ = reshape(a_x_u_load_,[n_x_u*n_x_u,n_x_u])*x_u_pack_;
a_x_u_pack_ = reshape(permute(reshape(a_x_u_pack_,[n_x_u,n_x_u,x_u_res]),[3,1,2]),[n_x_u*x_u_res,n_x_u])*x_u_pack_;
a_x_u_pack_ = reshape(permute(reshape(a_x_u_pack_,[x_u_res,n_x_u,x_u_res]),[3,1,2]),[x_u_res*x_u_res,n_x_u])*x_u_pack_;
a_x_u_pack_ = permute(reshape(a_x_u_pack_,[x_u_res,x_u_res,x_u_res]),[3,1,2]);
flag_plot=0; if flag_plot; figure(1); subplot(1,2,1); isosurface_f_x_c_0(a_x_u_load_,98.5); subplot(1,2,2); isosurface_f_x_c_0(a_x_u_pack_,98.5); end;
%%%%%%%%;
x_u_0_ = linspace(-x_p_r_max,+x_p_r_max,x_u_res);
x_u_1_ = linspace(-x_p_r_max,+x_p_r_max,x_u_res);
x_u_2_ = linspace(-x_p_r_max,+x_p_r_max,x_u_res);
[X_u_0_,X_u_1_,X_u_2_] = ndgrid(x_u_0_,x_u_1_,x_u_2_); n_X_u = x_u_res^3;
X_u_weight_ = (2*x_p_r_max/x_u_res)^3;

%%%%%%%%;
% Now convert to a_k_p_ ;
%%%%%%%%;
verbose=0;
fname_mat = sprintf('%s_mat/test_principled_marching_trpv1_15_a_k_p_quad_.mat',dir_trunk);
load(fname_mat);

%%%%%%%%;
% Now convert to a_k_Y_ ; 
%%%%%%%%;
verbose=0;
fname_mat = sprintf('%s_mat/test_principled_marching_trpv1_15_a_k_Y_quad_.mat',dir_trunk);
load(fname_mat);

%%%%%%%%
% Now generate templates. ;
%%%%%%%%;
verbose=0;
fname_mat = sprintf('%s_mat/test_principled_marching_trpv1_15_S_k_p__.mat',dir_trunk);
load(fname_mat);

%%%%%%%%;
% For images we will use the same n_w_ as that used by templates. ;
% We will also use the same quadrature weights for integration in 2d. ;
%%%%%%%%;

%%%%%%%%;
% extract ctf function. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/test_principled_marching_trpv1_15_CTF_k_p__.mat',dir_trunk);
load(fname_mat);

fname_euler_angle = sprintf('%s/euler_angles',dir_data);
fp = fopen(fname_euler_angle,'r');
euler_angle_load_ = textscan(fp,'%f%f%f\n%f%f\n');
fclose(fp);
euler_angle_marina_ = zeros(3,n_image); %<-- [polar_a,azimu_b,gamma_z]. ;
for nimage=0:n_image-1;
euler_angle_marina_(:,1+nimage) = convert_euler_relion_to_marina([euler_angle_load_{1}(1+nimage),euler_angle_load_{2}(1+nimage),euler_angle_load_{3}(1+nimage)]);
end;%for nimage=0:n_image-1;
delta_read_x_ = euler_angle_load_{4}*(2/n_x_u);
delta_read_y_ = euler_angle_load_{5}*(2/n_x_u);
%%%%%%%%;
delta_sigma = 0.0 * std([delta_read_x_;delta_read_y_]);
str_cost = sprintf('2d_xcor_d%.4d',floor(1000*delta_sigma));
%%%%%%%%;
if ( strcmp(str_CTF_factor,'original'));
str_combine = sprintf('%s',str_cost);
end;%if ( strcmp(str_CTF_factor,'original'));
if (~strcmp(str_CTF_factor,'original'));
str_combine = sprintf('%s_%s',str_cost,str_CTF_factor);
end;%if (~strcmp(str_CTF_factor,'original'));

verbose=1;
%%%%%%%%;
% Now test out principled marching. ;
% First set up cost matrix. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/test_principled_marching_trpv1_15_X_%s.mat',dir_trunk,str_cost);
load(fname_mat);

%%%%%%%%;
% Now generate principled-[templates*CTF_avg]. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/test_principled_marching_trpv1_15_UX_%s_S_CTF_k_p___.mat',dir_trunk,str_cost);
load(fname_mat);

%%%%%%%%;
% Now look at some of the experimental images associated with trpv1. ;
% Note that these will *NOT* be corrected/centered according to the (presumed) displacements. ;
% Moreover, to save space, we will not store the M_x_c___. ;
%%%%%%%%;
fname_M_x_c_mat = sprintf('%s_mat/test_principled_marching_trpv1_15_M_x_c___.mat',dir_trunk);
load(fname_M_x_c_mat);
fname_mat = sprintf('%s_mat/test_principled_marching_trpv1_15_UX_%s_M_k_p___.mat',dir_trunk,str_combine);
load(fname_mat);

fname_mat = sprintf('%s_mat/test_principled_marching_trpv1_15_c_k_Y_.mat',dir_trunk);
load(fname_mat);

%%%%%%%%;
% load prior results. ;
%%%%%%%%;
tmp_B_X_ = load('./dir_principled_marching_mat/test_principled_marching_trpv1_15_X_2d_xcor_d0000.mat');
tmp_A_X_ = load('./dir_principled_marching_mat/test_principled_marching_trpv1_8_X_2d_xcor.mat');
disp(sprintf(' %% A vs B X__: %0.16f',fnorm(tmp_A_X_.X_ - tmp_B_X_.X__)/fnorm(tmp_B_X_.X__)));
disp(sprintf(' %% A vs B X_weight_r_: %0.16f',fnorm(tmp_A_X_.X_weight_r_ - tmp_B_X_.X_weight_r_)/fnorm(tmp_B_X_.X_weight_r_)));
tmp_B_MS_ = load('./dir_principled_marching_mat/tpmhtameux_2d_xcor_d0000_d0000le7v001_n1024_MS_nUX001rng000.mat');
disp(transpose(tmp_B_MS_.X_best_MS_));
tmp_B_SM_ = load('./dir_principled_marching_mat/tpmhtameux_2d_xcor_d0000_d0000le7v001_n1024_SM_nUX001rng000.mat');
disp(transpose(tmp_B_SM_.X_best_SM_));
tmp_A_MS_ = load('./dir_principled_marching_mat/tpmham_experimental_4_UX_2d_xcor_original_Y_n1024_allshell_rng000nUX001.mat');
disp(transpose(tmp_A_MS_.dat_X_best_allshell_));
tmp_A_SM_ = load('./dir_principled_marching_mat/tpmham_experimental_4_UX_2d_xcor_original_Y_n1024_allshell_SM_rng000nUX001.mat');
disp(transpose(tmp_A_SM_.dat_X_best_allshell_SM_));

%%%%%%%%;
% Now step through tpmhtam_experimental_10.m ;
%%%%%%%%;
delta_r_max = 0;
svd_eps = 1e-3;
dat_n_order = 5;
dat_n_M = n_image_sub;
dat_n_UX_rank = 5;
dat_n_iteration = 32;
dat_rseed = 0;
n_delta_v_requested=1;
dat_infix = 'delete_B';
dat_M_k_p__ = M_k_p__;
%%%%%%%%;
% enter tpmhtam_experimental_10.m ;
%%%%%%%%;
delta_r_p = 0.05;
delta_r_s = delta_r_max/sqrt(2*log(1/delta_r_p));
delta_r_N = delta_r_max * (2*pi*k_p_r_max) / (pi*sqrt(2));
tmp_t = tic();
FTK = ampmh_FTK_1(n_k_p_r,k_p_r_,k_p_r_max,delta_r_max,svd_eps,n_delta_v_requested);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% FTK: %0.3fs',tmp_t)); end;
disp(sprintf(' %% p-val %0.4f delta_r_max %0.6f sigma %0.4f N_pixel %0.4f --> FTK.n_svd_l %d',delta_r_p,delta_r_max,delta_r_s,delta_r_N,FTK.n_svd_l));
%%%%%%%%;
% prepare precomputation for ampm. ;
%%%%%%%%;
pm_n_UX_rank = dat_n_UX_rank;
tmp_t = tic();
dat_M_k_q__ = zeros(n_w_sum,dat_n_M);
for nM=0:dat_n_M-1;
dat_M_k_q__(:,1+nM) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,dat_M_k_p__(:,1+nM));
end;%for nM=0:dat_n_M-1;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% dat_M_k_q__: %0.3fs',tmp_t)); end;
tmp_t = tic();
svd_VUXM_lwnM____ = tpmh_VUXM_lwnM____0(FTK,n_k_p_r,n_w_,dat_n_M,dat_M_k_q__,pm_n_UX_rank,UX__,X_weight_r_);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% svd_VUXM_lwnM____: %0.3fs',tmp_t)); end;
tmp_t = tic();
UX_M_l2_dM__ = ampmh_UX_M_l2_dM__1(FTK,n_w_,dat_n_M,pm_n_UX_rank,svd_VUXM_lwnM____);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% UX_M_l2_dM__: %0.3fs',tmp_t)); end;
[dat_UX_M_k_q_wnM_d0___,dat_UX_M_k_p_wnM_d0___] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,pm_n_UX_rank,dat_n_M,svd_VUXM_lwnM____(:,:,1:pm_n_UX_rank,:),zeros(dat_n_M,1),zeros(dat_n_M,1));
tmp_t = tic();
dat_M_k_q__ = zeros(n_w_sum,dat_n_M);
%%%%%%%%;
for nM=0:dat_n_M-1;
dat_M_k_q__(:,1+nM) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,dat_M_k_p__(:,1+nM));
end;%for nM=0:dat_n_M-1;
%%%%%%%%;
pm_n_UX_rank = dat_n_UX_rank;
pm_n_k_p_r = pm_n_UX_rank;
pm_k_p_r_ = ones(pm_n_k_p_r,1);
pm_k_p_r_max = 1;
pm_n_w_ = n_w_max*ones(pm_n_k_p_r,1);
pm_n_w_max = n_w_max;
pm_n_w_sum = sum(pm_n_w_);
pm_n_w_csum_ = cumsum([0;pm_n_w_]);
pm_l_max_ = l_max_max*ones(pm_n_k_p_r,1);
pm_n_lm_ = (1+pm_l_max_).^2; pm_n_lm_sum = sum(pm_n_lm_);
pm_weight_3d_k_p_r_ = ones(pm_n_k_p_r,1);
pm_weight_2d_k_p_r_ = ones(pm_n_k_p_r,1);
pm_N_pixel = delta_r_max * (2*pi*k_p_r_max) / (pi*sqrt(2)) ; 
%%%%%%%%;
% Now set up alternating minimization for 'MS-phase' of successive combinations of principled-image-rings. ;
%%%%%%%%;
nUX_rank = 1;
tmp_rseed=dat_rseed;tmp_n_iteration=dat_n_iteration;tmp_n_iteration_register=1;
tmp_viewing_k_eq_d = 1/(2*pi)/sqrt(1);tmp_n_order = dat_n_order;
tmp_pm_n_UX_rank = 1+nUX_rank;
tmp_pm_n_k_p_r = tmp_pm_n_UX_rank;
tmp_pm_k_p_r_ = ones(tmp_pm_n_k_p_r,1);
tmp_pm_k_p_r_max = 1;
tmp_pm_n_w_ = n_w_max*ones(tmp_pm_n_k_p_r,1);
tmp_pm_n_w_max = n_w_max;
tmp_pm_n_w_sum = sum(tmp_pm_n_w_);
tmp_pm_n_w_csum_ = cumsum([0;tmp_pm_n_w_]);
tmp_pm_l_max_ = l_max_max*ones(tmp_pm_n_k_p_r,1);
tmp_pm_n_lm_ = (1+tmp_pm_l_max_).^2; tmp_pm_n_lm_sum = sum(tmp_pm_n_lm_);
tmp_pm_weight_3d_k_p_r_ = ones(tmp_pm_n_k_p_r,1);
tmp_pm_weight_2d_k_p_r_ = ones(tmp_pm_n_k_p_r,1);
flag_MS_vs_SM = 1;
f_rand = 0;
a_UX_Y_true_ = reshape(a_UX_Y_quad__(:,1:tmp_pm_n_UX_rank),[n_lm_max*tmp_pm_n_UX_rank,1]); %<-- use as ground truth for this set of principled-image-rings. ;
CTF_index_ = 1;
CTF_k_p__ = ones(n_w_sum,1);
tmp_t = tic();
[tmp_B_MS_X_best_MS_...
,tmp_B_MS_UX_Y_0lsq_MS__...
,tmp_B_MS_euler_polar_a_MS__...
,tmp_B_MS_euler_azimu_b_MS__...
,tmp_B_MS_euler_gamma_z_MS__...
,tmp_B_MS_image_delta_x_MS__...
,tmp_B_MS_image_delta_y_MS__...
,tmp_B_MS_image_X_value_MS__...
,tmp_B_MS_image_S_index_MS__...
] = ...
ampm_3(...
 tmp_rseed...
,tmp_n_iteration...
,tmp_n_iteration_register...
,tmp_viewing_k_eq_d...
,tmp_n_order...
,n_k_p_r...
,k_p_r_...
,k_p_r_max...
,weight_3d_k_p_r_...
,weight_2d_k_p_r_...
,delta_r_max...
,svd_eps...
,n_delta_v_requested...
,FTK...
,n_w_...
,tmp_pm_n_UX_rank...
,UX__(:,1:tmp_pm_n_UX_rank)...
,X_weight_r_...
,dat_n_M...
,dat_M_k_p__...
,dat_M_k_q__...
,CTF_index_...
,CTF_k_p__...
,svd_VUXM_lwnM____(:,:,1:tmp_pm_n_UX_rank,:)...
,[]...
,l_max_...
,a_UX_Y_true_...
,[]...
,[]...
,[]...
,[]...
,[]...
,flag_MS_vs_SM...
,f_rand...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% X_best_MS_: %0.3fs',tmp_t)); end;

%%%%%%%%;
% enter tpmham_experimental_4b.m ;
%%%%%%%%;
dat_UX_M_k_p_wMn_d0___ = permute(dat_UX_M_k_p_wnM_d0___,[1,3,2]);
dat_UX_M_k_q_wMn_d0___ = permute(dat_UX_M_k_q_wnM_d0___,[1,3,2]);
%%%%%%%%;
tmp_verbose=0;
tmp_n_k_p_r = 1+nUX_rank;
tmp_k_p_r_ = ones(1+nUX_rank,1); tmp_k_p_r_max = 1;
tmp_n_w_ = n_w_max*ones(1+nUX_rank,1); tmp_n_w_max = max(tmp_n_w_); tmp_n_w_sum = sum(tmp_n_w_); 
tmp_weight_k_p_r_ = ones(1+nUX_rank,1); tmp_weight_2d_k_p_r_ = ones(1+nUX_rank,1);
tmp_l_max_ = l_max_max*ones(1+nUX_rank,1);
tmp_n_lm_ = (1+tmp_l_max_).^2; tmp_n_lm_sum = sum(tmp_n_lm_);
tmp_a_k_Y_quad_ = reshape(a_UX_Y_quad__(:,1+(0:nUX_rank)),[tmp_n_lm_sum,1]); %<-- use as ground truth for this set of principled-image-rings. ;
tmp_M_k_p__ = reshape(permute(dat_UX_M_k_p_wMn_d0___(:,1+(0:dat_n_M-1),1+(0:nUX_rank)),[1,3,2]),[tmp_n_w_sum,dat_n_M]);
tmp_M_k_q__ = reshape(permute(dat_UX_M_k_q_wMn_d0___(:,1+(0:dat_n_M-1),1+(0:nUX_rank)),[1,3,2]),[tmp_n_w_sum,dat_n_M]);
flag_MS_vs_SM = 1;
%%%%%%%%;
tmp_rseed=dat_rseed;
tmp_n_iteration=dat_n_iteration;
tmp_n_iteration_register=1;
tmp_viewing_k_eq_d = 1/(2*pi);
tmp_n_order = dat_n_order;
[tmp_A_MS_dat_X_best_allshell_...
,~...
,tmp_A_MS_euler_polar_a__...
,tmp_A_MS_euler_azimu_b__...
,tmp_A_MS_euler_gamma_z__] = ...
am_2(...
 tmp_rseed...
,tmp_n_iteration...
,tmp_n_iteration_register...
,tmp_viewing_k_eq_d...
,tmp_n_order...
,tmp_n_k_p_r...
,tmp_k_p_r_...
,tmp_k_p_r_max...
,tmp_weight_k_p_r_...
,tmp_weight_2d_k_p_r_...
,tmp_n_w_...
,dat_n_M...
,tmp_M_k_p__...
,tmp_M_k_q__...
,1 ...
,ones(tmp_n_w_sum,1)...
,tmp_l_max_...
,tmp_a_k_Y_quad_...
,[]...
,[]...
,[]...
,flag_MS_vs_SM...
,[]...
);

%%%%%%%%;
% Note that this gives the same result as the tmp_A_ runs!. ;
% Note, however, that the actual UX__ and UX_M_k_p___ are the same for these two runs: ;
%%%%%%%%
tmp_A_M_ = load('./dir_principled_marching_mat/test_principled_marching_trpv1_8_UX_2d_xcor_original_M_k_p___.mat');
tmp_B_M_ = load('./dir_principled_marching_mat/test_principled_marching_trpv1_15_UX_2d_xcor_d0000_M_k_p___.mat');
tmp_sign_ = [1,1,-1,1,1,1,-1,-1,-1]; %<-- special case. ;
for tmp_nUX_rank=0:n_UX_rank-1; 
disp(sprintf(' %% A vs B UX__ mode %d error: %0.16f',tmp_nUX_rank,fnorm(tmp_sign_(1+tmp_nUX_rank)*tmp_B_X_.UX__(:,1+tmp_nUX_rank)-tmp_A_X_.UX_(:,1+tmp_nUX_rank)))); 
end;%for tmp_nUX_rank=0:n_UX_rank-1; 
for tmp_nUX_rank=0:n_UX_rank-1;
disp(sprintf(' %% A vs B UX_M_k_p___ mode %d error: %0.16f',tmp_nUX_rank,fnorm(tmp_sign_(1+tmp_nUX_rank)*tmp_B_M_.UX_M_k_p___(:,:,1+tmp_nUX_rank) - tmp_A_M_.UX_M_k_p___(:,:,1+tmp_nUX_rank))/fnorm(tmp_B_M_.UX_M_k_p___(:,:,1+tmp_nUX_rank))));
end;%for tmp_nUX_rank=0:n_UX_rank-1;
%%%%%%%%;
% However, for some reason the conversion from svd_VUXM_lwnM____ does not produce the same princpal images. ;
%%%%%%%%;
% See, e.g., the following: ;
%%%%%%%%;
tmp_M_k_p__ = reshape(permute(tmp_B_M_.UX_M_k_p___(:,1+(0:dat_n_M-1),1+(0:nUX_rank)),[1,3,2]),[tmp_n_w_sum,dat_n_M]);
tmp_M_k_q__ = reshape(permute(tmp_B_M_.UX_M_k_q___(:,1+(0:dat_n_M-1),1+(0:nUX_rank)),[1,3,2]),[tmp_n_w_sum,dat_n_M]);
[tmp_A_MS_dat_X_best_allshell_...
,~...
,tmp_A_MS_euler_polar_a__...
,tmp_A_MS_euler_azimu_b__...
,tmp_A_MS_euler_gamma_z__] = ...
am_2(...
 tmp_rseed...
,tmp_n_iteration...
,tmp_n_iteration_register...
,tmp_viewing_k_eq_d...
,tmp_n_order...
,tmp_n_k_p_r...
,tmp_k_p_r_...
,tmp_k_p_r_max...
,tmp_weight_k_p_r_...
,tmp_weight_2d_k_p_r_...
,tmp_n_w_...
,dat_n_M...
,tmp_M_k_p__...
,tmp_M_k_q__...
,1 ...
,ones(tmp_n_w_sum,1)...
,tmp_l_max_...
,tmp_a_k_Y_quad_...
,[]...
,[]...
,[]...
,flag_MS_vs_SM...
,[]...
);

%%%%%%%%;
% Note that the tmp_M_k_p__ used above is identical to tmp_B_M_.UX_M_k_p___. ;
%%%%%%%%;
tmp_M_k_p__ = reshape(permute(tmp_B_M_.UX_M_k_p___(:,1+(0:dat_n_M-1),1+(0:nUX_rank)),[1,3,2]),[tmp_n_w_sum,dat_n_M]);
tmp_M_k_q__ = reshape(permute(tmp_B_M_.UX_M_k_q___(:,1+(0:dat_n_M-1),1+(0:nUX_rank)),[1,3,2]),[tmp_n_w_sum,dat_n_M]);
disp(sprintf(' %% tmp_M_k_p__ - reshape(permute(tmp_B_M_.UX_M_k_p___)) %0.16f',fnorm(tmp_M_k_p__ - reshape(permute(tmp_B_M_.UX_M_k_p___(:,:,1:1+nUX_rank),[1,3,2]),tmp_pm_n_w_sum,dat_n_M))));
disp(sprintf(' %% tmp_M_k_q__ - reshape(permute(tmp_B_M_.UX_M_k_q___)) %0.16f',fnorm(tmp_M_k_q__ - reshape(permute(tmp_B_M_.UX_M_k_q___(:,:,1:1+nUX_rank),[1,3,2]),tmp_pm_n_w_sum,dat_n_M))));

%%%%%%%%;
% Now check the svd_VUXM calculation. ;
%%%%%%%%;
dat_n_M = n_image_sub;
M_k_q__ = zeros(n_w_sum,dat_n_M);
for nM=0:dat_n_M-1;
M_k_q__(:,1+nM) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_k_p__(:,1+nM));
end;%for nM=0:dat_n_M-1;
%%%%%%%%;
tmp_image_delta_j_ = max(0,min(FTK.n_delta_v-1,periodize(0:dat_n_M-1,0,FTK.n_delta_v)));
tmp_image_delta_x_ = FTK.delta_x_(1+tmp_image_delta_j_); tmp_image_delta_x_ = tmp_image_delta_x_(:);
tmp_image_delta_y_ = FTK.delta_y_(1+tmp_image_delta_j_); tmp_image_delta_y_ = tmp_image_delta_y_(:);
%%%%%%%%;
UX_M_k_q_wnM_0___ = zeros(n_w_max,pm_n_UX_rank,dat_n_M); %<-- brute. ;
UX_M_k_p_wnM_0___ = zeros(n_w_max,pm_n_UX_rank,dat_n_M); %<-- brute. ;
UX_M_k_q_wnM_1___ = zeros(n_w_max,pm_n_UX_rank,dat_n_M); %<-- blas. ;
UX_M_k_p_wnM_1___ = zeros(n_w_max,pm_n_UX_rank,dat_n_M); %<-- blas. ;
tmp_t = tic();
for nM=0:dat_n_M-1;
tmp_M_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p__(:,1+nM),tmp_image_delta_x_(1+nM),tmp_image_delta_y_(1+nM));
tmp_M_k_q_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,tmp_M_k_p_);
for tmp_nUX_rank=0:pm_n_UX_rank-1;
tmp_UX_M_k_q_ = zeros(n_w_max,1);
for nk_p_r=0:n_k_p_r-1;
tmp_n_w = n_w_(1+nk_p_r);
tmp_n_w_2 = round(tmp_n_w/2);
tmp_ij_set_ = (0:tmp_n_w_2-1); tmp_ij_get_ = n_w_csum_(1+nk_p_r)+(0:tmp_n_w_2-1);
tmp_UX_M_k_q_(1+tmp_ij_set_) = tmp_UX_M_k_q_(1+tmp_ij_set_) + UX__(1+nk_p_r,1+tmp_nUX_rank)*X_weight_r_(1+nk_p_r)*tmp_M_k_q_(1+tmp_ij_get_);
tmp_ij_set_ = (tmp_n_w_2+1:tmp_n_w-1)-tmp_n_w+n_w_max; tmp_ij_get_ = n_w_csum_(1+nk_p_r)+(tmp_n_w_2+1:tmp_n_w-1);
tmp_UX_M_k_q_(1+tmp_ij_set_) = tmp_UX_M_k_q_(1+tmp_ij_set_) + UX__(1+nk_p_r,1+tmp_nUX_rank)*X_weight_r_(1+nk_p_r)*tmp_M_k_q_(1+tmp_ij_get_);
end;%for nk_p_r=0:n_k_p_r-1;
UX_M_k_q_wnM_0___(:,1+tmp_nUX_rank,1+nM) = tmp_UX_M_k_q_;
UX_M_k_p_wnM_0___(:,1+tmp_nUX_rank,1+nM) = interp_q_to_p(1,n_w_max,n_w_max,tmp_UX_M_k_q_);
end;%for tmp_nUX_rank=0:pm_n_UX_rank-1;
end;%for nM=0:dat_n_M-1;
tmp_t = toc(tmp_t); disp(sprintf(' %% UX_M_k_p_wnM_0___: %0.6fs',tmp_t));
tmp_t = tic();
[UX_M_k_q_wnM_1___,UX_M_k_p_wnM_1___] = ampmh_M_k_p__to_UX_M_k_p_wnM___0(n_k_p_r,k_p_r_,n_w_,pm_n_UX_rank,UX__,X_weight_r_,dat_n_M,M_k_p__,M_k_q__,tmp_image_delta_x_,tmp_image_delta_y_);
tmp_t = toc(tmp_t); disp(sprintf(' %% UX_M_k_p_wnM_1___: %0.6fs',tmp_t));
for tmp_nUX_rank=0:pm_n_UX_rank-1;
tmp_0_ = UX_M_k_q_wnM_0___(:,1+tmp_nUX_rank,:);
tmp_1_ = permute(tmp_B_M_.UX_M_k_q___(:,:,1+tmp_nUX_rank),[1,3,2]);
subplot(2,pm_n_UX_rank,1+0*pm_n_UX_rank+tmp_nUX_rank); 
plot(real(tmp_0_(:)),real(tmp_1_(:)),'.'); title(sprintf('mode %d: real',tmp_nUX_rank)); axis equal; grid on;
subplot(2,pm_n_UX_rank,1+1*pm_n_UX_rank+tmp_nUX_rank); 
plot(imag(tmp_0_(:)),imag(tmp_1_(:)),'.'); title(sprintf('mode %d: imag',tmp_nUX_rank)); axis equal; grid on;
disp(sprintf(' %% mode %d: UX_M_k_q_wnM_0___ vs tmp_B_M_.UX_M_k_q___: %0.16f',tmp_nUX_rank,fnorm(tmp_0_-tmp_1_)/fnorm(tmp_0_)));
clear tmp_0_ tmp_1_;
end;%for tmp_nUX_rank=0:pm_n_UX_rank-1;
disp(sprintf(' %% UX_M_k_q_wnM_0___ vs UX_M_k_q_wnM_1___: %0.16f',fnorm(UX_M_k_q_wnM_0___-UX_M_k_q_wnM_1___)/fnorm(UX_M_k_q_wnM_0___)));
disp(sprintf(' %% UX_M_k_p_wnM_0___ vs UX_M_k_p_wnM_1___: %0.16f',fnorm(UX_M_k_p_wnM_0___-UX_M_k_p_wnM_1___)/fnorm(UX_M_k_p_wnM_0___)));
%%%%%%%%;
tmp_t = tic();
svd_VUXM_nMwl____ = tpmh_VUXM_nMwl____0(FTK,n_k_p_r,n_w_,dat_n_M,M_k_q__,pm_n_UX_rank,UX__,X_weight_r_);
tmp_t = toc(tmp_t); disp(sprintf(' %% svd_VUXM_nMwl____: %0.6fs',tmp_t));
UX_M_k_q_wnM_2___ = zeros(n_w_max,pm_n_UX_rank,dat_n_M); %<-- svd. ;
UX_M_k_p_wnM_2___ = zeros(n_w_max,pm_n_UX_rank,dat_n_M); %<-- svd. ;
tmp_t = tic();
for nM=0:dat_n_M-1;
for tmp_nUX_rank=0:pm_n_UX_rank-1;
for nl=0:FTK.n_svd_l-1;
UX_M_k_q_wnM_2___(:,1+tmp_nUX_rank,1+nM) = UX_M_k_q_wnM_2___(:,1+tmp_nUX_rank,1+nM) + FTK.svd_U_d_expiw_s__(1+tmp_image_delta_j_(1+nM),1+nl)*reshape(svd_VUXM_nMwl____(1+tmp_nUX_rank,1+nM,:,1+nl),[n_w_max,1]) * n_w_max;
end;%for nl=0:FTK.n_svd_l-1;
UX_M_k_p_wnM_2___(:,1+tmp_nUX_rank,1+nM) = interp_q_to_p(1,n_w_max,n_w_max,UX_M_k_q_wnM_2___(:,1+tmp_nUX_rank,1+nM));
end;%for tmp_nUX_rank=0:pm_n_UX_rank-1;
end;%for nM=0:dat_n_M-1;
tmp_t = toc(tmp_t); disp(sprintf(' %% UX_M_k_p_wnM_2___: %0.6fs',tmp_t));
disp(sprintf(' %% UX_M_k_q_wnM_0___ vs UX_M_k_q_wnM_2___: %0.16f',fnorm(UX_M_k_q_wnM_0___-UX_M_k_q_wnM_2___)/fnorm(UX_M_k_q_wnM_0___)));
disp(sprintf(' %% UX_M_k_p_wnM_0___ vs UX_M_k_p_wnM_2___: %0.16f',fnorm(UX_M_k_p_wnM_0___-UX_M_k_p_wnM_2___)/fnorm(UX_M_k_p_wnM_0___)));
%%%%%%%%;
UX_M_k_q_wnM_3___ = zeros(n_w_max,pm_n_UX_rank,dat_n_M); %<-- svd with blas. ;
UX_M_k_p_wnM_3___ = zeros(n_w_max,pm_n_UX_rank,dat_n_M); %<-- svd with blas. ;
tmp_t = tic();
%svd_VUXM_lwnM____ = permute(svd_VUXM_nMwl____,[4,3,1,2]);
svd_VUXM_lwnM____ = tpmh_VUXM_lwnM____0(FTK,n_k_p_r,n_w_,dat_n_M,M_k_q__,pm_n_UX_rank,UX__,X_weight_r_);
tmp_t = toc(tmp_t); disp(sprintf(' %% svd_VUXM_lwnM____: %0.6fs',tmp_t));
tmp_t = tic();
for nM=0:dat_n_M-1;
UX_M_k_q_wnM_3___(:,:,1+nM) = reshape(FTK.svd_U_d_expiw_s__(1+tmp_image_delta_j_(1+nM),:)*reshape(svd_VUXM_lwnM____(:,:,:,1+nM),[FTK.n_svd_l,n_w_max*pm_n_UX_rank]),[n_w_max,pm_n_UX_rank]) * n_w_max;
end;%for nM=0:dat_n_M-1;
tmp_t = toc(tmp_t); disp(sprintf(' %% UX_M_k_q_wnM_3___: %0.6fs',tmp_t));
tmp_t = tic();
for nM=0:dat_n_M-1;
for tmp_nUX_rank=0:pm_n_UX_rank-1;
UX_M_k_p_wnM_3___(:,1+tmp_nUX_rank,1+nM) = interp_q_to_p(1,n_w_max,n_w_max,UX_M_k_q_wnM_3___(:,1+tmp_nUX_rank,1+nM));
end;%for tmp_nUX_rank=0:pm_n_UX_rank-1;
end;%for nM=0:dat_n_M-1;
tmp_t = toc(tmp_t); disp(sprintf(' %% UX_M_k_p_wnM_3___: %0.6fs',tmp_t));
tmp_t = tic();
UX_M_k_p_wnM_3___ = ifft(UX_M_k_q_wnM_3___,[],1)*sqrt(n_w_max);
tmp_t = toc(tmp_t); disp(sprintf(' %% UX_M_k_p_wnM_3___: %0.6fs',tmp_t));
disp(sprintf(' %% UX_M_k_q_wnM_0___ vs UX_M_k_q_wnM_3___: %0.16f',fnorm(UX_M_k_q_wnM_0___-UX_M_k_q_wnM_3___)/fnorm(UX_M_k_q_wnM_0___)));
disp(sprintf(' %% UX_M_k_p_wnM_0___ vs UX_M_k_p_wnM_3___: %0.16f',fnorm(UX_M_k_p_wnM_0___-UX_M_k_p_wnM_3___)/fnorm(UX_M_k_p_wnM_0___)));
%%%%%%%%;
svd_eps = 1e-3;
FTK = ampmh_FTK_0(n_k_p_r,k_p_r_,k_p_r_max,delta_r_max,svd_eps);
[UX_M_k_q_wnM_5___,UX_M_k_p_wnM_5___] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,pm_n_UX_rank,dat_n_M,svd_VUXM_lwnM____,tmp_image_delta_x_,tmp_image_delta_y_);
disp(sprintf(' %% UX_M_k_q_wnM_3___ vs UX_M_k_q_wnM_5___: %0.16f',fnorm(UX_M_k_q_wnM_3___-UX_M_k_q_wnM_5___)/fnorm(UX_M_k_q_wnM_3___)));
disp(sprintf(' %% UX_M_k_p_wnM_3___ vs UX_M_k_p_wnM_5___: %0.16f',fnorm(UX_M_k_p_wnM_3___-UX_M_k_p_wnM_5___)/fnorm(UX_M_k_p_wnM_3___)));
%%%%%%%%;
% Note that the tmp_M_k_p__ used above is identical to UX_M_k_p_wnM_3___. ;
%%%%%%%%;
tmp_M_k_p__ = reshape(permute(dat_UX_M_k_p_wMn_d0___(:,1+(0:dat_n_M-1),1+(0:nUX_rank)),[1,3,2]),[tmp_n_w_sum,dat_n_M]);
tmp_M_k_q__ = reshape(permute(dat_UX_M_k_q_wMn_d0___(:,1+(0:dat_n_M-1),1+(0:nUX_rank)),[1,3,2]),[tmp_n_w_sum,dat_n_M]);
disp(sprintf(' %% tmp_M_k_p__ - UX_M_k_p_wnM_3___ %0.16f',fnorm(tmp_M_k_p__ - reshape(UX_M_k_p_wnM_3___(:,1:1+nUX_rank,:),tmp_pm_n_w_sum,dat_n_M))));
disp(sprintf(' %% tmp_M_k_q__ - UX_M_k_q_wnM_3___ %0.16f',fnorm(tmp_M_k_q__ - reshape(UX_M_k_q_wnM_3___(:,1:1+nUX_rank,:),tmp_pm_n_w_sum,dat_n_M))));
disp(sprintf(' %% tmp_M_k_p__ - UX_M_k_p_wnM_0___ %0.16f',fnorm(tmp_M_k_p__ - reshape(UX_M_k_p_wnM_0___(:,1:1+nUX_rank,:),tmp_pm_n_w_sum,dat_n_M))));
disp(sprintf(' %% tmp_M_k_q__ - UX_M_k_q_wnM_0___ %0.16f',fnorm(tmp_M_k_q__ - reshape(UX_M_k_q_wnM_0___(:,1:1+nUX_rank,:),tmp_pm_n_w_sum,dat_n_M))));

%%%%%%%%;
% Now try again but with uniformly sampled M_k_p_. ;
%%%%%%%%;
tmp_B_M_x_c_ = load('./dir_principled_marching_mat/test_principled_marching_trpv1_8_M_x_c___.mat');
M_k_p_uni__ = zeros(n_w_max*n_k_p_r,n_image_sub); %<-- holds yes centered images. ;
for nimage_sub=0:n_image_sub-1;
if (mod(nimage_sub,100)==0); disp(sprintf(' %% nimage_sub %d/%d',nimage_sub,n_image_sub)); end;
tmp_M_x_c_ = squeeze(tmp_B_M_x_c_.M_x_c___(:,:,1+nimage_sub));
tmp_O_k_p_ = interp_x_c_to_k_p_nufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,tmp_M_x_c_,n_k_p_r,k_p_r_,n_w_max*ones(n_k_p_r,1)) ;
tmp_M_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_max*ones(n_k_p_r,1),n_w_max*n_k_p_r,tmp_O_k_p_,+1*delta_read_x_(1+nimage_sub),+1*delta_read_y_(1+nimage_sub));
tmp_M_k_p__ = reshape(tmp_M_k_p_,n_w_max,n_k_p_r);
M_k_p_uni__(:,1+nimage_sub) = tmp_M_k_p__(:);
end;%for nimage_sub=0:n_image_sub-1;
%%%%%%%%;
n_w_uni_ = n_w_max*ones(n_k_p_r,1);
n_w_uni_sum = sum(n_w_uni_);
n_w_uni_max = n_w_max;
n_w_uni_csum_ = cumsum([0;n_w_uni_]);
dat_n_M = n_image_sub;
M_k_q_uni__ = zeros(n_w_uni_sum,dat_n_M);
for nM=0:dat_n_M-1;
M_k_q_uni__(:,1+nM) = interp_p_to_q(n_k_p_r,n_w_uni_,n_w_uni_sum,M_k_p_uni__(:,1+nM));
end;%for nM=0:dat_n_M-1;
%%%%%%%%;
tmp_image_delta_j_ = max(0,min(FTK.n_delta_v-1,periodize(0:dat_n_M-1,0,FTK.n_delta_v)));
tmp_image_delta_x_ = FTK.delta_x_(1+tmp_image_delta_j_); tmp_image_delta_x_ = tmp_image_delta_x_(:);
tmp_image_delta_y_ = FTK.delta_y_(1+tmp_image_delta_j_); tmp_image_delta_y_ = tmp_image_delta_y_(:);
%%%%%%%%;
UX_M_k_q_uni_wnM_0___ = zeros(n_w_uni_max,pm_n_UX_rank,dat_n_M); %<-- brute. ;
UX_M_k_p_uni_wnM_0___ = zeros(n_w_uni_max,pm_n_UX_rank,dat_n_M); %<-- brute. ;
UX_M_k_q_uni_wnM_1___ = zeros(n_w_uni_max,pm_n_UX_rank,dat_n_M); %<-- blas. ;
UX_M_k_p_uni_wnM_1___ = zeros(n_w_uni_max,pm_n_UX_rank,dat_n_M); %<-- blas. ;
tmp_t = tic();
for nM=0:dat_n_M-1;
tmp_M_k_p_uni_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_uni_,n_w_uni_sum,M_k_p_uni__(:,1+nM),tmp_image_delta_x_(1+nM),tmp_image_delta_y_(1+nM));
tmp_M_k_q_uni_ = interp_p_to_q(n_k_p_r,n_w_uni_,n_w_uni_sum,tmp_M_k_p_uni_);
for tmp_nUX_rank=0:pm_n_UX_rank-1;
tmp_UX_M_k_q_uni_ = zeros(n_w_uni_max,1);
for nk_p_r=0:n_k_p_r-1;
tmp_n_w_uni = n_w_uni_(1+nk_p_r);
tmp_n_w_uni_2 = round(tmp_n_w_uni/2);
tmp_ij_set_ = (0:tmp_n_w_uni_2-1); tmp_ij_get_ = n_w_uni_csum_(1+nk_p_r)+(0:tmp_n_w_uni_2-1);
tmp_UX_M_k_q_uni_(1+tmp_ij_set_) = tmp_UX_M_k_q_uni_(1+tmp_ij_set_) + UX__(1+nk_p_r,1+tmp_nUX_rank)*X_weight_r_(1+nk_p_r)*tmp_M_k_q_uni_(1+tmp_ij_get_);
tmp_ij_set_ = (tmp_n_w_uni_2+1:tmp_n_w_uni-1)-tmp_n_w_uni+n_w_uni_max; tmp_ij_get_ = n_w_uni_csum_(1+nk_p_r)+(tmp_n_w_uni_2+1:tmp_n_w_uni-1);
tmp_UX_M_k_q_uni_(1+tmp_ij_set_) = tmp_UX_M_k_q_uni_(1+tmp_ij_set_) + UX__(1+nk_p_r,1+tmp_nUX_rank)*X_weight_r_(1+nk_p_r)*tmp_M_k_q_uni_(1+tmp_ij_get_);
end;%for nk_p_r=0:n_k_p_r-1;
UX_M_k_q_uni_wnM_0___(:,1+tmp_nUX_rank,1+nM) = tmp_UX_M_k_q_uni_;
UX_M_k_p_uni_wnM_0___(:,1+tmp_nUX_rank,1+nM) = interp_q_to_p(1,n_w_uni_max,n_w_uni_max,tmp_UX_M_k_q_uni_);
end;%for tmp_nUX_rank=0:pm_n_UX_rank-1;
end;%for nM=0:dat_n_M-1;
tmp_t = toc(tmp_t); disp(sprintf(' %% UX_M_k_p_uni_wnM_0___: %0.6fs',tmp_t));
tmp_t = tic();
[UX_M_k_q_uni_wnM_1___,UX_M_k_p_uni_wnM_1___] = ampmh_M_k_p__to_UX_M_k_p_wnM___0(n_k_p_r,k_p_r_,n_w_uni_,pm_n_UX_rank,UX__,X_weight_r_,dat_n_M,M_k_p_uni__,M_k_q_uni__,tmp_image_delta_x_,tmp_image_delta_y_);
tmp_t = toc(tmp_t); disp(sprintf(' %% UX_M_k_p_uni_wnM_1___: %0.6fs',tmp_t));
for tmp_nUX_rank=0:pm_n_UX_rank-1;
tmp_0_ = UX_M_k_q_uni_wnM_0___(:,1+tmp_nUX_rank,:);
tmp_1_ = permute(tmp_B_M_.UX_M_k_q___(:,:,1+tmp_nUX_rank),[1,3,2]);
subplot(2,pm_n_UX_rank,1+0*pm_n_UX_rank+tmp_nUX_rank); 
plot(real(tmp_0_(:)),real(tmp_1_(:)),'.'); title(sprintf('mode %d: real',tmp_nUX_rank)); axis equal; grid on;
subplot(2,pm_n_UX_rank,1+1*pm_n_UX_rank+tmp_nUX_rank); 
plot(imag(tmp_0_(:)),imag(tmp_1_(:)),'.'); title(sprintf('mode %d: imag',tmp_nUX_rank)); axis equal; grid on;
disp(sprintf(' %% mode %d: UX_M_k_q_uni_wnM_0___ vs tmp_B_M_.UX_M_k_q___: %0.16f',tmp_nUX_rank,fnorm(tmp_0_-tmp_1_)/fnorm(tmp_0_)));
clear tmp_0_ tmp_1_;
end;%for tmp_nUX_rank=0:pm_n_UX_rank-1;
disp(sprintf(' %% UX_M_k_q_uni_wnM_0___ vs UX_M_k_q_uni_wnM_1___: %0.16f',fnorm(UX_M_k_q_uni_wnM_0___-UX_M_k_q_uni_wnM_1___)/fnorm(UX_M_k_q_uni_wnM_0___)));
disp(sprintf(' %% UX_M_k_p_uni_wnM_0___ vs UX_M_k_p_uni_wnM_1___: %0.16f',fnorm(UX_M_k_p_uni_wnM_0___-UX_M_k_p_uni_wnM_1___)/fnorm(UX_M_k_p_uni_wnM_0___)));

%%%%%%%%;
% Now rerun tmp_B_MS_ using uniformly sampled M_k_p__. ;
% Now set up alternating minimization for 'MS-phase' of successive combinations of principled-image-rings. ;
%%%%%%%%;
svd_VUXM_uni_lwnM____ = tpmh_VUXM_lwnM____0(FTK,n_k_p_r,n_w_uni_,dat_n_M,M_k_q_uni__,pm_n_UX_rank,UX__,X_weight_r_);
nUX_rank = 1;
tmp_rseed=dat_rseed;tmp_n_iteration=dat_n_iteration;tmp_n_iteration_register=1;
tmp_viewing_k_eq_d = 1/(2*pi)/sqrt(1);tmp_n_order = dat_n_order;
tmp_pm_n_UX_rank = 1+nUX_rank;
tmp_pm_n_k_p_r = tmp_pm_n_UX_rank;
tmp_pm_k_p_r_ = ones(tmp_pm_n_k_p_r,1);
tmp_pm_k_p_r_max = 1;
tmp_pm_n_w_ = n_w_max*ones(tmp_pm_n_k_p_r,1);
tmp_pm_n_w_max = n_w_max;
tmp_pm_n_w_sum = sum(tmp_pm_n_w_);
tmp_pm_n_w_csum_ = cumsum([0;tmp_pm_n_w_]);
tmp_pm_l_max_ = l_max_max*ones(tmp_pm_n_k_p_r,1);
tmp_pm_n_lm_ = (1+tmp_pm_l_max_).^2; tmp_pm_n_lm_sum = sum(tmp_pm_n_lm_);
tmp_pm_weight_3d_k_p_r_ = ones(tmp_pm_n_k_p_r,1);
tmp_pm_weight_2d_k_p_r_ = ones(tmp_pm_n_k_p_r,1);
flag_MS_vs_SM = 1;
f_rand = 0;
a_UX_Y_true_ = reshape(a_UX_Y_quad__(:,1:tmp_pm_n_UX_rank),[n_lm_max*tmp_pm_n_UX_rank,1]); %<-- use as ground truth for this set of principled-image-rings. ;
CTF_index_ = 1;
CTF_k_p_uni__ = ones(n_w_uni_sum,1);
tmp_t = tic();
[tmp_B_MS_X_best_MS_...
,tmp_B_MS_UX_Y_0lsq_MS__...
,tmp_B_MS_euler_polar_a_MS__...
,tmp_B_MS_euler_azimu_b_MS__...
,tmp_B_MS_euler_gamma_z_MS__...
,tmp_B_MS_image_delta_x_MS__...
,tmp_B_MS_image_delta_y_MS__...
,tmp_B_MS_image_X_value_MS__...
,tmp_B_MS_image_S_index_MS__...
] = ...
ampm_3(...
 tmp_rseed...
,tmp_n_iteration...
,tmp_n_iteration_register...
,tmp_viewing_k_eq_d...
,tmp_n_order...
,n_k_p_r...
,k_p_r_...
,k_p_r_max...
,weight_3d_k_p_r_...
,weight_2d_k_p_r_...
,delta_r_max...
,svd_eps...
,n_delta_v_requested...
,FTK...
,n_w_uni_...
,tmp_pm_n_UX_rank...
,UX__(:,1:tmp_pm_n_UX_rank)...
,X_weight_r_...
,dat_n_M...
,M_k_p_uni__...
,M_k_q_uni__...
,CTF_index_...
,CTF_k_p_uni__...
,svd_VUXM_uni_lwnM____(:,:,1:tmp_pm_n_UX_rank,:)...
,[]...
,l_max_...
,a_UX_Y_true_...
,[]...
,[]...
,[]...
,[]...
,[]...
,flag_MS_vs_SM...
,f_rand...
);

%%%%%%%%;
% Now try running tmp_B_SM_ using uniformly sampled M_k_p__. ;
% Now set up alternating minimization for 'MS-phase' of successive combinations of principled-image-rings. ;
%%%%%%%%;
nUX_rank = 1;
tmp_rseed=dat_rseed;tmp_n_iteration=dat_n_iteration;tmp_n_iteration_register=1;
tmp_viewing_k_eq_d = 1/(2*pi)/sqrt(1);tmp_n_order = dat_n_order;
tmp_pm_n_UX_rank = 1+nUX_rank;
tmp_pm_n_k_p_r = tmp_pm_n_UX_rank;
tmp_pm_k_p_r_ = ones(tmp_pm_n_k_p_r,1);
tmp_pm_k_p_r_max = 1;
tmp_pm_n_w_ = n_w_max*ones(tmp_pm_n_k_p_r,1);
tmp_pm_n_w_max = n_w_max;
tmp_pm_n_w_sum = sum(tmp_pm_n_w_);
tmp_pm_n_w_csum_ = cumsum([0;tmp_pm_n_w_]);
tmp_pm_l_max_ = l_max_max*ones(tmp_pm_n_k_p_r,1);
tmp_pm_n_lm_ = (1+tmp_pm_l_max_).^2; tmp_pm_n_lm_sum = sum(tmp_pm_n_lm_);
tmp_pm_weight_3d_k_p_r_ = ones(tmp_pm_n_k_p_r,1);
tmp_pm_weight_2d_k_p_r_ = ones(tmp_pm_n_k_p_r,1);
tmp_euler_polar_a_MS_ = tmp_A_MS_.euler_polar_a__(:,end);
tmp_euler_azimu_b_MS_ = tmp_A_MS_.euler_azimu_b__(:,end);
tmp_euler_gamma_z_MS_ = tmp_A_MS_.euler_gamma_z__(:,end);
tmp_image_delta_x_MS_ = zeros(dat_n_M,1);
tmp_image_delta_y_MS_ = zeros(dat_n_M,1);
flag_MS_vs_SM = 0;
a_UX_Y_true_ = reshape(a_UX_Y_quad__(:,1:tmp_pm_n_UX_rank),[n_lm_max*tmp_pm_n_UX_rank,1]); %<-- use as ground truth for this set of principled-image-rings. ;
CTF_index_ = 1;
CTF_k_p_uni__ = ones(n_w_uni_sum,1);
f_rand_ = [0.00,0.01,0.025,0.05,0.10,0.15]; n_f_rand = numel(f_rand_);
for nf_rand=0:n_f_rand-1;
f_rand = f_rand_(1+nf_rand);
tmp_t = tic();
[tmp_B_SM_X_best_SM_{1+nf_rand}...
,tmp_B_SM_UX_Y_0lsq_SM__{1+nf_rand}...
,tmp_B_SM_euler_polar_a_SM__{1+nf_rand}...
,tmp_B_SM_euler_azimu_b_SM__{1+nf_rand}...
,tmp_B_SM_euler_gamma_z_SM__{1+nf_rand}...
,tmp_B_SM_image_delta_x_SM__{1+nf_rand}...
,tmp_B_SM_image_delta_y_SM__{1+nf_rand}...
,tmp_B_SM_image_X_value_SM__{1+nf_rand}...
,tmp_B_SM_image_S_index_SM__{1+nf_rand}...
] = ...
ampm_3(...
 tmp_rseed...
,tmp_n_iteration...
,tmp_n_iteration_register...
,tmp_viewing_k_eq_d...
,tmp_n_order...
,n_k_p_r...
,k_p_r_...
,k_p_r_max...
,weight_3d_k_p_r_...
,weight_2d_k_p_r_...
,delta_r_max...
,svd_eps...
,n_delta_v_requested...
,FTK...
,n_w_uni_...
,tmp_pm_n_UX_rank...
,UX__(:,1:tmp_pm_n_UX_rank)...
,X_weight_r_...
,dat_n_M...
,M_k_p_uni__...
,M_k_q_uni__...
,CTF_index_...
,CTF_k_p_uni__...
,svd_VUXM_uni_lwnM____(:,:,1:tmp_pm_n_UX_rank,:)...
,[]...
,l_max_...
,a_UX_Y_true_...
,tmp_euler_polar_a_MS_...
,tmp_euler_azimu_b_MS_...
,tmp_euler_gamma_z_MS_...
,tmp_image_delta_x_MS_...
,tmp_image_delta_y_MS_...
,flag_MS_vs_SM...
,f_rand...
);
end;%for nf_rand=0:n_f_rand-1;
c_ = colormap_beach(); n_c = size(c_,1);
figure(1);clf;hold on;
for nf_rand=0:n_f_rand-1;
nc = max(0,min(n_c-1,floor(n_c*nf_rand/n_f_rand)));
plot(1:tmp_n_iteration,tmp_B_SM_X_best_SM_{1+nf_rand},'.-','LineWidth',2,'Color',c_(1+nc,:));
end;%for nf_rand=0:n_f_rand-1;
hold off;
figbig;
legend({'0','0.01','0.025','0.05','0.10','0.15'},'Location','NorthWest');
title('SM for various f_rand','Interpreter','none');
print('-depsc','./dir_principled_marching_jpg/f_rand_test.eps');
print('-djpeg','./dir_principled_marching_jpg/f_rand_test.jpg');
return;
