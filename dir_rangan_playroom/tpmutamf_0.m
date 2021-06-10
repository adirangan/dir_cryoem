function ...
tpmutamf_0( ...
 dir_jpg ...
,fname_0 ...
,fname_1 ...
,fname_2 ...
,n_M ...
,M_k_p__ ...
,nUX_rank ...
,n_UX_rank ...
,n_iteration ...
,n_order ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,delta_r_max_use ...
,delta_r_max_upb ...
,svd_eps ...
,n_delta_v_requested ...
,n_CTF_rank ...
,CTF_idx_ ...
,CTF_k_p_r__ ...
,CTF_k_p__ ...
,l_max_ ...
,UX__ ...
,X_weight_r_ ...
,a_CTF_UX_Y_quad__ ...
,dat_f_rand ...
,a_k_Y_quad_ ...
,euler_angle_marina_ ...
,delta_read_x_ ...
,delta_read_y_ ...
,delta_sigma ...
);

% generates figures for data-file. ;
% test_principled_marching_updated_translation_alternating_minimization_figures.m ;

%%%%%%%%;
fname_mat = sprintf('%s_%s_%s.mat',fname_0,fname_1,fname_2);
tmp_ = load(fname_mat);
if strcmp(fname_1,'SM');
mat_X_best_ = tmp_.X_best_SM_;
mat_a_CTF_UX_Y_0lsq__ = tmp_.a_CTF_UX_Y_0lsq_SM__;
mat_euler_polar_a__ = tmp_.euler_polar_a_SM__;
mat_euler_azimu_b__ = tmp_.euler_azimu_b_SM__;
mat_euler_gamma_z__ = tmp_.euler_gamma_z_SM__;
mat_image_delta_x__ = tmp_.image_delta_x_SM__;
mat_image_delta_y__ = tmp_.image_delta_y_SM__;
mat_image_I_value__ = tmp_.image_I_value_SM__;
mat_image_X_value__ = tmp_.image_X_value_SM__;
mat_image_S_index__ = tmp_.image_S_index_SM__;
mat_X_best_time = tmp_.X_best_SM_time;
mat_dat_M_loading__ = tmp_.dat_M_loading_SM__;
mat_dat_M_loading_time = tmp_.dat_M_loading_SM_time;
end;%if strcmp(fname_1,'SM');
if strcmp(fname_1,'MS');
mat_X_best_ = tmp_.X_best_MS_;
mat_a_CTF_UX_Y_0lsq__ = tmp_.a_CTF_UX_Y_0lsq_MS__;
mat_euler_polar_a__ = tmp_.euler_polar_a_MS__;
mat_euler_azimu_b__ = tmp_.euler_azimu_b_MS__;
mat_euler_gamma_z__ = tmp_.euler_gamma_z_MS__;
mat_image_delta_x__ = tmp_.image_delta_x_MS__;
mat_image_delta_y__ = tmp_.image_delta_y_MS__;
mat_image_I_value__ = tmp_.image_I_value_MS__;
mat_image_X_value__ = tmp_.image_X_value_MS__;
mat_image_S_index__ = tmp_.image_S_index_MS__;
mat_X_best_time = tmp_.X_best_MS_time;
mat_dat_M_loading__ = tmp_.dat_M_loading_MS__;
mat_dat_M_loading_time = tmp_.dat_M_loading_MS_time;
end;%if strcmp(fname_1,'MS');
clear tmp_;
%%%%%%%%;
n_w_max = max(n_w_);
l_max_max = max(l_max_);
n_lm_max = (1+l_max_max).^2;
%%%%%%%%;
pm_n_UX_rank = nUX_rank+1;
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
a_CTF_UX_Y_true_ = reshape(a_CTF_UX_Y_quad__(:,1:pm_n_UX_rank),[n_lm_max*pm_n_UX_rank,1]);
%%%%%%%%;

c2d__ = colormap_polar_a_azimu_b_2d(+euler_angle_marina_(1,1:n_M),+euler_angle_marina_(2,1:n_M),0.35);
n_h = 64; lh_lim_ = [0,10];
figure(3);
figbeach();
figbig;
subplot(2,3,1);
hold on;
patch(2*pi*[+0;+1;+1;+0],1*pi*[+0;+0;+1;+1],0.65*[1,1,1],'EdgeColor','none');
scatter(0*pi+euler_angle_marina_(2,1:n_M),1*pi-euler_angle_marina_(1,1:n_M),64,c2d__(1:n_M,:),'filled');
hold off;
axisnotick;
axis([0,2*pi,0,1*pi]); axis equal; xlabel('azimu_b'); ylabel('polar_a'); title('true view');
subplot(2,3,4); cla;
imagesc(log2(1+hist2d_0(0*pi+euler_angle_marina_(2,1:n_M),1*pi-euler_angle_marina_(1,1:n_M),2*n_h,n_h,[0,2*pi],[0,1*pi])),lh_lim_);
set(gca,'Ydir','normal');
axisnotick; axis image;
for niteration=0:n_iteration-1;
%%%%%%%%;
N_wavelength = 0.0;
[ ...
 X_best_orig ...
,polar_a_best_orig ...
,azimu_b_best_orig ...
,gamma_z_best_orig ...
,delta_best_orig_ ...
,a_k_Y_best_orig_ ...
,b_k_Y_best_orig_ ...
] = ...
register_spharm_to_spharm_wigner_0( ...
 pm_n_k_p_r ...
,pm_k_p_r_ ...
,pm_k_p_r_max ...
,pm_weight_3d_k_p_r_ ...
,N_wavelength ...
,pm_l_max_ ...
,a_CTF_UX_Y_true_ ...
,mat_a_CTF_UX_Y_0lsq__(:,end) ...
);
[ ...
 X_best_flip ...
,polar_a_best_flip ...
,azimu_b_best_flip ...
,gamma_z_best_flip ...
,delta_best_flip_ ...
,a_k_Y_best_flip_ ...
,b_k_Y_best_flip_ ...
] = ...
register_spharm_to_spharm_wigner_0( ...
 pm_n_k_p_r ...
,pm_k_p_r_ ...
,pm_k_p_r_max ...
,pm_weight_3d_k_p_r_ ...
,N_wavelength ...
,pm_l_max_ ...
,a_CTF_UX_Y_true_ ...
,flipY(pm_n_k_p_r,pm_l_max_,mat_a_CTF_UX_Y_0lsq__(:,end)) ...
);
%%%%%%%%;
X_best = X_best_orig;
polar_a_best = polar_a_best_orig;
azimu_b_best = azimu_b_best_orig;
gamma_z_best = gamma_z_best_orig;
delta_best_ = delta_best_orig_;
a_k_Y_best_ = a_k_Y_best_orig_;
b_k_Y_best_ = b_k_Y_best_orig_;
if (X_best_flip>X_best_orig);
X_best = X_best_flip;
polar_a_best = polar_a_best_flip;
azimu_b_best = azimu_b_best_flip;
gamma_z_best = gamma_z_best_flip;
delta_best_ = delta_best_flip_;
a_k_Y_best_ = a_k_Y_best_flip_;
b_k_Y_best_ = b_k_Y_best_flip_;
end;%if (X_best_flip>X_best_orig);
euler_best_ = [+azimu_b_best,+polar_a_best,+gamma_z_best];
R_best_ = euler_to_R(euler_best_);
%%%%%%%%;
tmp_euler_azimu_b_ = mat_euler_azimu_b__(:,1+niteration);
tmp_euler_polar_a_ = mat_euler_polar_a__(:,1+niteration);
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
tmp_euler_azimu_b_dif_ = pi/1 + periodize((0*pi+transpose(euler_angle_marina_(2,1:n_M))) - (0*pi+tmp_euler_azimu_b_),-pi,+pi);
tmp_euler_polar_a_dif_ = pi/2 + (1*pi-transpose(euler_angle_marina_(1,1:n_M))) - (1*pi-tmp_euler_polar_a_);
%%%%%%%%;
figure(3);
subplot(2,3,2); cla;
hold on;
patch(2*pi*[+0;+1;+1;+0],1*pi*[+0;+0;+1;+1],0.65*[1,1,1],'EdgeColor','none');
scatter(0*pi+tmp_euler_azimu_b_,1*pi-tmp_euler_polar_a_,64,c2d__(1:n_M,:),'filled');
hold off;
axisnotick;
axis([0,2*pi,0,1*pi]); axis equal; xlabel('azimu_b'); ylabel('polar_a'); title('true view');
title(sprintf('%s %d',fname_1,niteration));
subplot(2,3,5); cla;
imagesc(log2(1+hist2d_0(0*pi+tmp_euler_azimu_b_,1*pi-tmp_euler_polar_a_,2*n_h,n_h,[0,2*pi],[0,1*pi])),lh_lim_);
set(gca,'Ydir','normal');
axisnotick; axis image;
subplot(2,3,3); cla;
hold on;
patch(2*pi*[+0;+1;+1;+0],1*pi*[+0;+0;+1;+1],0.65*[1,1,1],'EdgeColor','none');
scatter( ...
 tmp_euler_azimu_b_dif_ ...
,tmp_euler_polar_a_dif_ ...
,64 ...
,c2d__(1:n_M,:) ...
,'filled' ...
);
hold off;
axisnotick;
axis([0,2*pi,0,1*pi]); axis equal; xlabel('azimu_b'); ylabel('polar_a'); title('true view');
title(sprintf('true - %s %d',fname_1,niteration));
subplot(2,3,6); cla;
imagesc(log2(1+hist2d_0( ...
 tmp_euler_azimu_b_dif_ ...
,tmp_euler_polar_a_dif_ ...
,2*n_h,n_h,[0,2*pi],[0,1*pi])),lh_lim_);
set(gca,'Ydir','normal');
axisnotick; axis image;
drawnow();
end;%for niteration=0:n_iteration-1;

disp('return'); return;


%{
%%%%%%%%;
c_ = colormap_beach(); n_c = size(c_,1);
%hsv_ = colormap('hsv'); n_hsv = size(hsv_,1);
%nhsv_ = max(0,min(n_hsv-1,floor(n_hsv*(pi+atan2(delta_read_y_(1:n_M),delta_read_x_(1:n_M)))/(2*pi))));
%scatter(delta_read_x_(1:n_M),delta_read_y_(1:n_M),15,hsv_(1+nhsv_,:),'filled'); axis([-1,+1,-1,+1]); axis square;
figure(1); clf;
%%%%%%%%;
subplot(5,7,2);
c2d__ = colormap_gaussian_2d(delta_read_x_,delta_read_y_,delta_sigma,0.35);
hold on;
patch([-1;+1;+1;-1],[-1;-1;+1;+1],'k');
scatter(delta_read_x_(1:n_M),delta_read_y_(1:n_M),32,c2d__(1:n_M,:),'filled');
%scatter(delta_read_x_(1:end),delta_read_y_(1:end),25,c2d__(1:end,:),'filled');
hold off;
axisnotick;
axis([-1,+1,-1,+1]); axis square; xlabel('dx'); ylabel('dy'); title('true delta');
figbig;
%%%%%%%%;
subplot(5,7,1);
hold on;
plot(0:n_iteration-1,mat_X_best_,'o-','LineWidth',2);
xlim([0,n_iteration-1]);
ylim([-0.2,1.0]); grid on;
xlabel('iteration'); ylabel('correlation');
hold off;
for niteration=0:n_iteration-1;
subplot(5,7,3+niteration);
hold on;
patch([-1;+1;+1;-1],[-1;-1;+1;+1],'k');
scatter(mat_image_delta_x__(:,1+niteration),mat_image_delta_y__(:,1+niteration),32,c2d__(1:n_M,:),'filled');
hold off;
axisnotick;
axis([-1,+1,-1,+1]); axis square;
%xlabel('dx'); ylabel('dy');
title(sprintf('%s %d',fname_1,niteration));
end;%for niteration=0:n_iteration-1;
%%%%%%%%;
 %}

%{
c2d__ = colormap_gaussian_2d(delta_read_x_,delta_read_y_,delta_sigma,0.35);
n_h = 64; lh_lim_ = [0,10];
figure(2);
figbeach();
figbig;
subplot(2,3,1);
hold on;
patch([-1;+1;+1;-1],[-1;-1;+1;+1],'k');
scatter(delta_read_x_(1:n_M),delta_read_y_(1:n_M),64,c2d__(1:n_M,:),'filled');
hold off;
axisnotick;
axis([-1,+1,-1,+1]); axis square; xlabel('dx'); ylabel('dy'); title('true delta');
subplot(2,3,4); cla;
imagesc(log2(1+hist2d_0(delta_read_x_(1:n_M),delta_read_y_(1:n_M),n_h,n_h,[-1,+1],[-1,+1])),lh_lim_);
set(gca,'Ydir','normal');
axisnotick; axis image;
for niteration=0:n_iteration-1;
figure(2);
subplot(2,3,2); cla;
hold on;
patch([-1;+1;+1;-1],[-1;-1;+1;+1],'k');
scatter(mat_image_delta_x__(:,1+niteration),mat_image_delta_y__(:,1+niteration),64,c2d__(1:n_M,:),'filled');
hold off;
axisnotick;
axis([-1,+1,-1,+1]); axis square;
%xlabel('dx'); ylabel('dy');
title(sprintf('%s %d',fname_1,niteration));
subplot(2,3,5); cla;
imagesc(log2(1+hist2d_0(mat_image_delta_x__(1:n_M,1+niteration),mat_image_delta_y__(1:n_M,1+niteration),n_h,n_h,[-1,+1],[-1,+1])),lh_lim_);
set(gca,'Ydir','normal');
axisnotick; axis image;
subplot(2,3,3); cla;
hold on;
patch([-1;+1;+1;-1],[-1;-1;+1;+1],'k');
scatter(delta_read_x_(1:n_M) - mat_image_delta_x__(:,1+niteration),delta_read_y_(1:n_M) - mat_image_delta_y__(:,1+niteration),16,c2d__(1:n_M,:),'filled');
hold off;
axisnotick;
axis([-1,+1,-1,+1]); axis square;
%xlabel('dx'); ylabel('dy');
title(sprintf('true - %s %d',fname_1,niteration));
subplot(2,3,6); cla;
imagesc(log2(1+hist2d_0(delta_read_x_(1:n_M) - mat_image_delta_x__(1:n_M,1+niteration),delta_read_y_(1:n_M) - mat_image_delta_y__(1:n_M,1+niteration),n_h,n_h,[-1,+1],[-1,+1])),lh_lim_);
set(gca,'Ydir','normal');
axisnotick; axis image;
drawnow();
end;%for niteration=0:n_iteration-1;
 %}

%{
tmp_euler_polar_a_ = mat_euler_polar_a__(:,end);
tmp_euler_azimu_b_ = mat_euler_azimu_b__(:,end);
tmp_euler_gamma_z_ = mat_euler_gamma_z__(:,end);
tmp_image_delta_x_ = mat_image_delta_x__(:,end);
tmp_image_delta_y_ = mat_image_delta_y__(:,end);
%%%%%%%%;
tmp_t = tic;
c_k_Y_reco_ = ...
cg_lsq_4( ...
 n_order ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,n_M ...
,M_k_p__ ...
,CTF_idx_ ...
,CTF_k_p__ ...
,tmp_euler_polar_a_ ...
,tmp_euler_azimu_b_ ...
,tmp_euler_gamma_z_ ...
,tmp_image_delta_x_ ...
,tmp_image_delta_y_ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% M_k_p__ --> c_k_Y_reco_ time %0.2fs',tmp_t));
[tmp_X_best_orig,~,~,~,~,~,~] = register_spharm_to_spharm_wigner_0(n_k_p_r,k_p_r_,k_p_r_max,weight_3d_k_p_r_,0,l_max_,a_k_Y_quad_,c_k_Y_reco_);
[tmp_X_best_flip,~,~,~,~,~,~] = register_spharm_to_spharm_wigner_0(n_k_p_r,k_p_r_,k_p_r_max,weight_3d_k_p_r_,0,l_max_,a_k_Y_quad_,flipY(n_k_p_r,l_max_,c_k_Y_reco_));
tmp_X_best_reco = max(tmp_X_best_orig,tmp_X_best_flip);
disp(sprintf(' %% tmp_X_best_reco %0.3f',tmp_X_best_reco));
 %}
