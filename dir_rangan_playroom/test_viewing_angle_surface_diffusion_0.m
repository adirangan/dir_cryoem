% test the relationship between viewing-angle dispersion and surface-diffusion. ;
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
k_eq_d = 1.0/(2*pi);
[n_k_all,azimu_b_all_,polar_a_all_,weight_3d_all_,k_c_0_all_,k_c_1_all_,k_c_2_all_,n_polar_a,polar_a_,n_azimu_b_] = sample_shell_5(k_p_r_max,k_eq_d,'L');
tmp_t=toc(tmp_t); if (verbose>-1); disp(sprintf(' %% sample_shell_5: %0.6fs',tmp_t)); end;
n_k_p_r = 1;
k_p_r_ = k_p_r_max;
weight_3d_k_p_r_ = (1/3)*k_p_r_max^3; %<-- sum(weight_3d_k_p_r_)*(4*pi) = (4/3)*pi*k_p_r_max^3 --> sum(weight_3d_k_p_r_) = (1/3)*k_p_r_max^3 ;
% note that sum(weight_3d_all_) = (4*pi)*k_p_r_max^2 ;

%%%%%%%%;
% set up spherical-harmonic coefficients. ;
%%%%%%%%;
tmp_t=tic();
l_max_upb = 36;%l_max_upb = 8;
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

fname_mat = sprintf('%s_mat/test_viewing_angle_surface_diffusion_0.mat',dir_trunk);
if ( exist(fname_mat,'file')); disp(sprintf(' %% %s found, not creating',fname_mat)); end;
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
figure(1); clf;

l_val_use_ = 0:6-1; n_l_val_use = numel(l_val_use_);
sigma_ = sqrt(linspace(0,0.1,32)); n_sigma = numel(sigma_);
F_diff___ = zeros(n_sigma,n_l_val_use,n_l_val_use);
F_0lsq___ = zeros(n_sigma,n_l_val_use,n_l_val_use);
L_diff___ = zeros(n_sigma,n_l_val_use,n_l_val_use);
L_0lsq___ = zeros(n_sigma,n_l_val_use,n_l_val_use);
a_k_Y_0lsq___ = cell(n_sigma,n_l_val_use,n_l_val_use);

%%%%%%%%;
% pick order and degree. ;
%%%%%%%%;
for nl_val_use = 0:n_l_val_use-1;;
l_val_use = l_val_use_(1+nl_val_use);

for m_val_use = -l_val_use:+l_val_use;
nm_val_use = m_val_use + l_val_use;
a_k_Y_true_ = zeros(n_lm_sum,1);
index_use = efind((Y_l_val_==l_val_use) & (Y_m_val_==m_val_use));
a_k_Y_true_(1+index_use) = 1;

%%%%%%%%;
% generate templates. ;
%%%%%%%%;
tmp_t=tic();
template_k_eq_d = k_eq_d*2;
viewing_k_eq_d = k_eq_d*8;
[ ...
 S_k_p__ ...
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
,a_k_Y_true_ ...
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
a_k_Y_0lsq_ = cg_lsq_1(lsq_n_order,n_k_p_r,l_max_,n_w_,n_S,S_k_p__,euler_polar_a_,euler_azimu_b_,euler_gamma_z_);
tmp_t=toc(tmp_t); if (verbose>1); disp(sprintf(' %% cg_lsq_1: %0.6fs',tmp_t)); end;
disp(sprintf(' %% a_k_Y_true_ vs a_k_Y_0lsq_: %0.16f',fnorm(a_k_Y_true_-a_k_Y_0lsq_)/fnorm(a_k_Y_true_)));

%%%%%%%%;
% Now perturb the viewing-angles by a gaussian-distribution with variance sigma^2 ;
%%%%%%%%;
for nsigma=0:n_sigma-1;
sigma = sigma_(1+nsigma);
euler_polar_a_ = viewing_polar_a_all_ ;
euler_azimu_b_ = viewing_azimu_b_all_ ;
euler_gamma_z_ = zeros(n_viewing_all,1);
for nviewing_all=0:n_viewing_all-1;
polar_a = viewing_polar_a_all_(1+nviewing_all);
azimu_b = viewing_azimu_b_all_(1+nviewing_all);
tmp_k0 = cos(azimu_b)*sin(polar_a);
tmp_k1 = sin(azimu_b)*sin(polar_a);
tmp_k2 = cos(polar_a);
tmp_k0 = tmp_k0 + sigma*randn();
tmp_k1 = tmp_k1 + sigma*randn();
tmp_k2 = tmp_k2 + sigma*randn();
tmp_k01 = sqrt(tmp_k0^2 + tmp_k1^2);
polar_a = atan2(tmp_k01,tmp_k2);
azimu_b = atan2(tmp_k1,tmp_k0);
gamma_z = euler_gamma_z_(1+nviewing_all);
gamma_z = gamma_z + sigma*randn();
euler_polar_a_(1+nviewing_all) = polar_a;
euler_azimu_b_(1+nviewing_all) = azimu_b;
euler_gamma_z_(1+nviewing_all) = gamma_z;
end;%for nviewing_all=0:n_viewing_all-1;
%%%%%%%%;
a_k_Y_diff_ = exp(-0.5*l_val_use*(l_val_use+1)*sigma^2)*a_k_Y_true_;
tmp_t=tic();
a_k_Y_0lsq_ = cg_lsq_1(lsq_n_order,n_k_p_r,l_max_,n_w_,n_S,S_k_p__,euler_polar_a_,euler_azimu_b_,euler_gamma_z_);
tmp_t=toc(tmp_t); if (verbose>1); disp(sprintf(' %% cg_lsq_1: %0.6fs',tmp_t)); end;
if (verbose>0); disp(sprintf(' %% sigma %0.2f: a_k_Y_diff_ vs a_k_Y_0lsq_: %0.16f',sigma,fnorm(a_k_Y_diff_-a_k_Y_0lsq_)/fnorm(a_k_Y_diff_))); end;
F_0lsq___(1+nsigma,1+nl_val_use,1+nm_val_use) = a_k_Y_0lsq_(1+index_use);
L_0lsq___(1+nsigma,1+nl_val_use,1+nm_val_use) = fnorm(a_k_Y_0lsq_);
a_k_Y_0lsq___{1+nsigma,1+nl_val_use,1+nm_val_use} = a_k_Y_0lsq_;
end;%for nsigma=0:n_sigma-1;
for nsigma=0:n_sigma-1;
sigma = sigma_(1+nsigma);
a_k_Y_diff_ = exp(-l_val_use*(l_val_use+1)*sigma^2)*a_k_Y_true_;
F_diff___(1+nsigma,1+nl_val_use,1+nm_val_use) = a_k_Y_diff_(1+index_use);
L_diff___(1+nsigma,1+nl_val_use,1+nm_val_use) = fnorm(a_k_Y_diff_);
end;%for nsigma=0:n_sigma-1;
subplot(2,3,1+l_val_use); hold on;
plot(sigma_.^2,F_diff___(:,1+nl_val_use,1+nm_val_use),'b.-');
plot(sigma_.^2,abs(F_0lsq___(:,1+nl_val_use,1+nm_val_use)),'ro-');
plot(sigma_.^2,abs(L_0lsq___(:,1+nl_val_use,1+nm_val_use)),'gx-');
xlim([0,sigma_(end).^2]); ylim([0,1]); grid on;
xlabel('variance'); ylabel('coefficient'); title(sprintf('l%+d m%+d',l_val_use,m_val_use)); legend({'Fdiff','F0lsq','L0lsq'});
hold off;
drawnow();
figbig;

end;%for m_val_use = -l_val_use:+l_val_use;
end;%for nl_val_use = 0:n_l_val_use-1;

save(fname_mat);
end;%if (~exist(fname_mat,'file'));
load(fname_mat);

figure(1);clf; hold on;
for nsigma=0:n_sigma-1;
sigma = sigma_(1+nsigma);
for nl_val_use=0:n_l_val_use-1;
l_val_use = l_val_use_(1+nl_val_use);
for nm_val_use=0:(1+2*l_val_use)-1;
m_val_use = -l_val_use + nm_val_use;
tmp_diff = exp(-l_val_use*(l_val_use+1)*sigma^2);
plot(abs(F_0lsq___(1+nsigma,1+nl_val_use,1+nm_val_use)),tmp_diff,'ro');
plot(abs(L_0lsq___(1+nsigma,1+nl_val_use,1+nm_val_use)),tmp_diff,'gx');
end;%for nm_val_use=0:(1+2*l_val_use)-1;
end;%for nl_val_use=0:n_l_val_use-1;
end;%for nsigma=0:n_sigma-1;
plot([0,1],[0,1],'k-');
hold off;
xlabel('0lsq'); ylabel('diff');
xlim([0,1]);ylim([0,1]);
axis square;

figure(1);clf;
var_ = sigma_.^2;
for nl_val_use=0:n_l_val_use-1;
l_val_use = l_val_use_(1+nl_val_use);
for nm_val_use=0:(1+2*l_val_use)-1;
m_val_use = -l_val_use + nm_val_use;
subplot(2,3,1+nl_val_use); hold on;
tmp_diff_ = exp(-l_val_use*(l_val_use+1)*sigma_.^2);
plot(var_,log(abs(tmp_diff_)),'b.-');
plot(var_,log(abs(F_0lsq___(:,1+nl_val_use,1+nm_val_use))),'ro-');
plot(var_,log(abs(L_0lsq___(:,1+nl_val_use,1+nm_val_use))),'gx-');
xlabel('variance'); ylabel('log(coefficient)');
xlim([0,var_(end)]);ylim([-5,0]);
end;%for nm_val_use=0:(1+2*l_val_use)-1;
end;%for nl_val_use=0:n_l_val_use-1;


