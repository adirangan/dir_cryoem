% See if we can construct an initial random function that has the following property: ;
% templates drawn from this function have bessel-coefficients with magnitude 1. ;
% Note that an iid-gaussian function should have bessel-coefficients iid-gaussian. ;
% We are trying to see if we can do better than that. ;
% Apparently not. ;

clear;

%platform = 'access1';
platform = 'OptiPlex';
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;

dir_trunk = sprintf('/%s/rangan/dir_cryoem/dir_rangan_playroom/dir_principled_marching',string_root);

verbose=0;

%%%%%%%%;
% generate spherical-shell. ;
%%%%%%%%;
tmp_t=tic();
k_p_r_max = 48/(2*pi);
k_eq_d = 1.0/(2*pi)*sqrt(0.5);
[ ...
n_k_all ...
,azimu_b_all_ ...
,polar_a_all_ ...
,weight_3d_all_ ...
,k_c_0_all_ ...
,k_c_1_all_ ...
,k_c_2_all_ ...
,n_polar_a ...
,polar_a_ ...
,n_azimu_b_ ...
] = ...
sample_shell_5( ...
 k_p_r_max ...
,k_eq_d ...
,'L' ...
);
tmp_t=toc(tmp_t); if (verbose>-1); disp(sprintf(' %% sample_shell_5: %0.6fs',tmp_t)); end;
n_k_p_r = 1;
k_p_r_ = k_p_r_max;
weight_3d_k_p_r_ = (1/3)*k_p_r_max^3; %<-- sum(weight_3d_k_p_r_)*(4*pi) = (4/3)*pi*k_p_r_max^3 --> sum(weight_3d_k_p_r_) = (1/3)*k_p_r_max^3 ;
% note that sum(weight_3d_all_) = (4*pi)*k_p_r_max^2 ;

%%%%%%%%;
% set up spherical-harmonic coefficients for an iid-gaussian function. ;
%%%%%%%%;
tmp_t=tic();
sigma_init = 1.0d0;
l_max_upb = 36;
%l_max_upb = 8;
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
tmp_ij_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
Y_l_val_(1+tmp_ij_) = tmp_l_val_;
Y_m_val_(1+tmp_ij_) = tmp_m_val_;
Y_k_val_(1+tmp_ij_) = k_p_r_(1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
weight_Y_ = zeros(n_lm_sum,1);
for nk_p_r=0:n_k_p_r-1;
tmp_ij_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
weight_Y_(1+tmp_ij_) = weight_3d_k_p_r_(1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;

fname_mat = sprintf('%s_mat/test_initial_random_function.mat',dir_trunk);
if ( exist(fname_mat,'file')); disp(sprintf(' %% %s found, not creating',fname_mat)); end;
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));

a_k_Y_rand_ = sigma_init*crandn(n_lm_sum,1);

%%%%%%%%;
% generate templates. ;
%%%%%%%%;
tmp_t=tic();
template_k_eq_d = 0.25*(2*pi*k_p_r_max) / (2*(l_max_max+1)) ; %<-- this only needs to be sufficiently large to sample each oscillation twice (i.e., nyquist). ;
viewing_k_eq_d = k_eq_d*8;
[ ...
 S_k_p_rand__ ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_k_all_ ...
,n_viewing_all ...
,viewing_azimu_b_all_ ...
,viewing_polar_a_all_ ...
,viewing_weight_all_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,template_k_c_0__ ...
,template_k_c_1__ ...
,template_k_c_2__ ...
] = ...
get_template_1( ...
 verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_rand_ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
);
n_S = n_viewing_all; n_w_max = max(n_w_); n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_]);
tmp_t=toc(tmp_t); if (verbose>-1); disp(sprintf(' %% get_template_0: %0.6fs',tmp_t)); end;

%%%%%%%%;
% Now reconstruct a_k_Y_0lsq_. ;
%%%%%%%%;
tmp_t=tic();
lsq_n_order = 5;
euler_polar_a_ = viewing_polar_a_all_ ;
euler_azimu_b_ = viewing_azimu_b_all_ ;
euler_gamma_z_ = zeros(n_viewing_all,1);
a_k_Y_0lsq_ = cg_lsq_1(lsq_n_order,n_k_p_r,l_max_,n_w_,n_S,S_k_p_rand__,euler_polar_a_,euler_azimu_b_,euler_gamma_z_);
tmp_t=toc(tmp_t); if (verbose>1); disp(sprintf(' %% cg_lsq_1: %0.6fs',tmp_t)); end;
disp(sprintf(' %% a_k_Y_rand_ vs a_k_Y_0lsq_: %0.16f',fnorm(a_k_Y_rand_-a_k_Y_0lsq_)/fnorm(a_k_Y_rand_)));

%%%%%%%%;
% generate templates with bessel-coefficients of uniform magnitude. ;
%%%%%%%%;
S_k_q_unif__ = zeros(n_w_max,n_S);
S_k_p_unif__ = zeros(n_w_max,n_S);
for nS=0:n_S-1;
S_k_q_unif__(:,1+nS) = exp(i*2*pi*rand(n_w_max,1));
S_k_p_unif__(:,1+nS) = interp_q_to_p(n_k_p_r,n_w_,n_w_sum,S_k_q_unif__(:,1+nS));
end;%for nS=0:n_S-1;

%%%%%%%%;
% Now reconstruct a_k_Y_ulsq_. ;
%%%%%%%%;
tmp_t=tic();
lsq_n_order = 5;
euler_polar_a_ = viewing_polar_a_all_ ;
euler_azimu_b_ = viewing_azimu_b_all_ ;
euler_gamma_z_ = zeros(n_viewing_all,1);
a_k_Y_ulsq_ = cg_lsq_1(lsq_n_order,n_k_p_r,l_max_,n_w_,n_S,S_k_p_unif__,euler_polar_a_,euler_azimu_b_,euler_gamma_z_);
tmp_t=toc(tmp_t); if (verbose>1); disp(sprintf(' %% cg_lsq_1: %0.6fs',tmp_t)); end;

%%%%%%%%;
% generate templates from a_k_Y_ulsq_. ;
%%%%%%%%;
tmp_t=tic();
[ ...
 S_k_p_ulsq__ ...
] = ...
get_template_1( ...
 verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_ulsq_ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
);
tmp_t=toc(tmp_t); if (verbose>-1); disp(sprintf(' %% get_template_0: %0.6fs',tmp_t)); end;

%%%%%%%%;
% try and apply alternating-minimization to the S_k_p_unif__. ;
%%%%%%%%;
rseed = 0;
n_iteration = 4;
n_iteration_register = n_iteration+1;
CTF_index_ = [];
CTF_k_p__ = [];
[ ...
 ~ ...
,a_k_Y_ulsq__ ...
,euler_polar_a__ ...
,euler_azimu_b__ ...
,euler_gamma_z__ ...
] = ...
am_0( ...
 rseed ...
,n_iteration ...
,n_iteration_register ...
,viewing_k_eq_d ...
,lsq_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_ ...
,n_S ...
,S_k_p_unif__ ...
,S_k_q_unif__ ...
,CTF_index_ ...
,CTF_k_p__ ...
,l_max_ ...
,zeros(n_lm_sum,1) ...
,[] ...
,[] ...
,[] ...
);

%%%%%%%%;
% generate templates from a_k_Y_ulsq__. ;
%%%%%%%%;
tmp_t=tic();
[ ...
 S_k_p_elsq__ ...
] = ...
get_template_1( ...
 verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_ulsq__(:,end) ...
,viewing_k_eq_d ...
,template_k_eq_d ...
);
tmp_t=toc(tmp_t); if (verbose>-1); disp(sprintf(' %% get_template_0: %0.6fs',tmp_t)); end;

%%%%%%%%;
% Now measure magnitude of the bessel-coefficients for the templates. ;
%%%%%%%%;
S_k_q_rand__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
S_k_q_rand__(:,1+nS) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,S_k_p_rand__(:,1+nS)); 
end;%for nS=0:n_S-1;
S_k_q_rand__ = S_k_q_rand__/fnorm(S_k_q_rand__);
q_ = periodize(0:n_w_sum-1,-n_w_sum/2,+n_w_sum/2);
S_k_q_rand_abs_avg_ = mean(abs(S_k_q_rand__),2);
S_k_q_rand_abs_std_ = std(abs(S_k_q_rand__),1,2);
S_k_q_rand_abs_p05_ = prctile(abs(S_k_q_rand__), 5,2);
S_k_q_rand_abs_p15_ = prctile(abs(S_k_q_rand__),15,2);
S_k_q_rand_abs_p50_ = prctile(abs(S_k_q_rand__),50,2);
S_k_q_rand_abs_p85_ = prctile(abs(S_k_q_rand__),85,2);
S_k_q_rand_abs_p95_ = prctile(abs(S_k_q_rand__),95,2);
%%%%%%%%;
S_k_q_ulsq__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
S_k_q_ulsq__(:,1+nS) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,S_k_p_ulsq__(:,1+nS)); 
end;%for nS=0:n_S-1;
S_k_q_ulsq__ = S_k_q_ulsq__/fnorm(S_k_q_ulsq__);
q_ = periodize(0:n_w_sum-1,-n_w_sum/2,+n_w_sum/2);
S_k_q_ulsq_abs_avg_ = mean(abs(S_k_q_ulsq__),2);
S_k_q_ulsq_abs_std_ = std(abs(S_k_q_ulsq__),1,2);
S_k_q_ulsq_abs_p05_ = prctile(abs(S_k_q_ulsq__), 5,2);
S_k_q_ulsq_abs_p15_ = prctile(abs(S_k_q_ulsq__),15,2);
S_k_q_ulsq_abs_p50_ = prctile(abs(S_k_q_ulsq__),50,2);
S_k_q_ulsq_abs_p85_ = prctile(abs(S_k_q_ulsq__),85,2);
S_k_q_ulsq_abs_p95_ = prctile(abs(S_k_q_ulsq__),95,2);
%%%%%%%%;
S_k_q_elsq__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
S_k_q_elsq__(:,1+nS) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,S_k_p_elsq__(:,1+nS)); 
end;%for nS=0:n_S-1;
S_k_q_elsq__ = S_k_q_elsq__/fnorm(S_k_q_elsq__);
q_ = periodize(0:n_w_sum-1,-n_w_sum/2,+n_w_sum/2);
S_k_q_elsq_abs_avg_ = mean(abs(S_k_q_elsq__),2);
S_k_q_elsq_abs_std_ = std(abs(S_k_q_elsq__),1,2);
S_k_q_elsq_abs_p05_ = prctile(abs(S_k_q_elsq__), 5,2);
S_k_q_elsq_abs_p15_ = prctile(abs(S_k_q_elsq__),15,2);
S_k_q_elsq_abs_p50_ = prctile(abs(S_k_q_elsq__),50,2);
S_k_q_elsq_abs_p85_ = prctile(abs(S_k_q_elsq__),85,2);
S_k_q_elsq_abs_p95_ = prctile(abs(S_k_q_elsq__),95,2);
%%%%%%%%;
figure(1);clf;
hold on;
plot(q_,S_k_q_rand_abs_p50_,'k-','LineWidth',2);
plot(q_,S_k_q_ulsq_abs_p50_,'r-','LineWidth',2);
plot(q_,S_k_q_elsq_abs_p50_,'g-','LineWidth',2);
%%%%%%%%;
plot(q_,S_k_q_rand_abs_p05_,'k-','LineWidth',0.5);
plot(q_,S_k_q_rand_abs_p15_,'k-','LineWidth',1);
plot(q_,S_k_q_rand_abs_p50_,'k-','LineWidth',2);
plot(q_,S_k_q_rand_abs_p85_,'k-','LineWidth',1);
plot(q_,S_k_q_rand_abs_p95_,'k-','LineWidth',0.5);
%%%%%%%%;
plot(q_,S_k_q_ulsq_abs_p05_,'r-','LineWidth',0.5);
plot(q_,S_k_q_ulsq_abs_p15_,'r-','LineWidth',1);
plot(q_,S_k_q_ulsq_abs_p50_,'r-','LineWidth',2);
plot(q_,S_k_q_ulsq_abs_p85_,'r-','LineWidth',1);
plot(q_,S_k_q_ulsq_abs_p95_,'r-','LineWidth',0.5);
%%%%%%%%;
plot(q_,S_k_q_elsq_abs_p05_,'g-','LineWidth',0.5);
plot(q_,S_k_q_elsq_abs_p15_,'g-','LineWidth',1);
plot(q_,S_k_q_elsq_abs_p50_,'g-','LineWidth',2);
plot(q_,S_k_q_elsq_abs_p85_,'g-','LineWidth',1);
plot(q_,S_k_q_elsq_abs_p95_,'g-','LineWidth',0.5);
%%%%%%%%;
xlim([-l_max_max-1,+l_max_max+1]); xlabel('q');
ylabel('bessel coord magnitude');
hold off;
legend({'rand','ulsq','elsq'});
set(gcf,'Position',1+[0,0,1536,512]);
print('-depsc',sprintf('%s_jpg/test_initial_random_function_FIGA.eps',dir_trunk));
print('-djpeg',sprintf('%s_jpg/test_initial_random_function_FIGA.jpg',dir_trunk));
 
%save(fname_mat);
end;%if (~exist(fname_mat,'file'));

