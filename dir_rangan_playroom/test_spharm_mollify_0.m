function test_spharm_mollify_0(l_val_use,m_val_use);
% tests mollification of spherical-harmonics. ;

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
% set up spherical-harmonic coefficients. ;
%%%%%%%%;
tmp_t=tic();
%l_max_upb = 36;
l_max_upb = 8;
l_max_(1) = l_val_use;
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

a_k_Y_ = zeros(n_lm_max,1);
index_use = efind(Y_l_val_==l_val_use & Y_m_val_==m_val_use);
a_k_Y_(1+index_use) = 1.0d0;

[a_k_all_] = ...
convert_spharm_to_k_p_1( ...
 verbose ...
,n_k_all ...
,0 ...
,k_p_r_max*ones(n_k_all,1) ...
,azimu_b_all_ ...
,polar_a_all_ ...
,weight_3d_all_ ...
,weight_3d_all_ ...
,1 ...
,k_p_r_max ...
,1 ...
,l_max_ ...
,a_k_Y_ ...
);

azimu_b = pi/6;
polar_a = pi/3;
gamma_z = pi/12;
euler_angle_ = [azimu_b,polar_a,gamma_z];
[b_k_Y_] = rotate_spharm_to_spharm_2(0,[],1,1,l_max_,a_k_Y_,euler_angle_);

[b_k_all_] = ...
convert_spharm_to_k_p_1( ...
 verbose ...
,n_k_all ...
,0 ...
,k_p_r_max*ones(n_k_all,1) ...
,azimu_b_all_ ...
,polar_a_all_ ...
,weight_3d_all_ ...
,weight_3d_all_ ...
,1 ...
,k_p_r_max ...
,1 ...
,l_max_ ...
,b_k_Y_ ...
);

flag_plot=1;
if flag_plot;
%%%%%%%%;
% plot a_k_Y_ and b_k_Y_ ;
%%%%%%%%;
subplot(2,2,1);
imagesc_polar_a_azimu_b_0(polar_a_all_,azimu_b_all_,real(a_k_all_)); xlim([0,2*pi]);ylim([0,pi]); axisnotick;
title(sprintf('l%d m%d real(a)',l_val_use,m_val_use));
subplot(2,2,2);
imagesc_polar_a_azimu_b_0(polar_a_all_,azimu_b_all_,imag(a_k_all_)); xlim([0,2*pi]);ylim([0,pi]); axisnotick;
title(sprintf('l%d m%d imag(a)',l_val_use,m_val_use));
subplot(2,2,3);
imagesc_polar_a_azimu_b_0(polar_a_all_,azimu_b_all_,real(b_k_all_)); xlim([0,2*pi]);ylim([0,pi]); axisnotick;
title(sprintf('l%d m%d real(b)',l_val_use,m_val_use));
subplot(2,2,4);
imagesc_polar_a_azimu_b_0(polar_a_all_,azimu_b_all_,imag(b_k_all_)); xlim([0,2*pi]);ylim([0,pi]); axisnotick;
title(sprintf('l%d m%d imag(b)',l_val_use,m_val_use));
end;%if flag_plot;

%%%%%%%%;
% sample a few points. ;
%%%%%%%%;
surface_area = 4*pi*k_p_r_max^2;
sigma = 0.025;
n_iteration = 12; flag_plot=0;
for niteration=0:n_iteration-1;
index_pt = round(n_k_all*niteration/n_iteration);
k_c_0_pt = k_c_0_all_(1+index_pt);
k_c_1_pt = k_c_1_all_(1+index_pt);
k_c_2_pt = k_c_2_all_(1+index_pt);
d2_on_sphere_ = (k_c_0_all_ - k_c_0_pt).^2 + (k_c_1_all_ - k_c_1_pt).^2 + (k_c_2_all_ - k_c_2_pt).^2;
gaussian_on_sphere_ = exp(-d2_on_sphere_/(2*(sigma*k_p_r_max)^2));
gaussian_on_sphere_ = gaussian_on_sphere_/sum(weight_3d_all_.*gaussian_on_sphere_) * surface_area;
if flag_plot;
subplot(3,4,1+niteration);
imagesc_polar_a_azimu_b_0(polar_a_all_,azimu_b_all_,gaussian_on_sphere_);
xlim([0,2*pi]);ylim([0,pi]); axisnotick;
end;%if flag_plot;
assert(abs(1-sum(weight_3d_all_.*gaussian_on_sphere_)/surface_area)<1e-6);
a_pt = a_k_all_(1+index_pt);
a_diffuse_pt = exp(-l_val_use*(l_val_use+1)*sigma^2)*a_pt;
a_mollify_pt = sum(weight_3d_all_.*a_k_all_.*gaussian_on_sphere_)/surface_area;
a_error_pt = fnorm(a_diffuse_pt - a_mollify_pt)/fnorm(a_diffuse_pt);
b_pt = b_k_all_(1+index_pt);
b_diffuse_pt = exp(-l_val_use*(l_val_use+1)*sigma^2)*b_pt;
b_mollify_pt = sum(weight_3d_all_.*b_k_all_.*gaussian_on_sphere_)/surface_area;
b_error_pt = fnorm(b_diffuse_pt - b_mollify_pt)/fnorm(b_diffuse_pt);
disp(sprintf(' %% a_error_pt %0.16f b_error_pt %0.16f',a_error_pt,b_error_pt));
end;%for niteration=0:n_iteration-1;
