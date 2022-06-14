function ...
[ ...
 parameter ...
] = ...
ampmut_align_to_a_CTF_avg_UX_Y_0( ...
 parameter ...
,l_max_max ...
,pm_n_UX_rank ...
,a_CTF_avg_UX_Y_true_ ...
,n_M ...
,euler_polar_a_true_ ...
,euler_azimu_b_true_ ...
,euler_gamma_z_true_ ...
,image_delta_x_true_ ...
,image_delta_y_true_ ...
,n_iteration ...
,a_CTF_avg_UX_Y_calc__ ...
,euler_polar_a_calc__ ...
,euler_azimu_b_calc__ ...
,euler_gamma_z_calc__ ...
,image_delta_x_calc__ ...
,image_delta_y_calc__ ...
);

verbose=2;
if (verbose); disp(sprintf(' %% [entering ampmut_align_to_a_CTF_avg_UX_Y_0]')); end;
if isempty(parameter);
parameter = struct('type','parameter');
end;%if isempty(parameter);

pm_n_k_p_r = pm_n_UX_rank;
pm_k_p_r_ = ones(pm_n_k_p_r,1);
pm_k_p_r_max = 1;
pm_l_max_ = l_max_max*ones(pm_n_k_p_r,1);
pm_n_lm_ = (1+pm_l_max_).^2; pm_n_lm_sum = sum(pm_n_lm_);
pm_weight_3d_k_p_r_ = ones(pm_n_k_p_r,1);
pm_weight_2d_k_p_r_ = ones(pm_n_k_p_r,1);

if ( isempty(n_iteration)); n_iteration=0; end;
if (n_iteration<=0); n_iteration = size(euler_polar_a_calc__,2); end;

flag_found=0;
if (~isempty(parameter.fname_align_a_CTF_avg_UX_Y_pre));
tmp_fname_mat = sprintf('%s.mat',parameter.fname_align_a_CTF_avg_UX_Y_pre);
if  exist(tmp_fname_mat,'file');
disp(sprintf(' %% %s found, loading',tmp_fname_mat));
load(tmp_fname_mat ...
     ,'X_best_','X_best_flag_flip_' ...
     ,'polar_a_best_' ...
     ,'azimu_b_best_' ...
     ,'gamma_z_best_' ...
     ,'delta_best__' ...
    );
flag_found=1;
end;%if  exist(tmp_fname_mat,'file');
end;%if (~isempty(parameter.fname_align_a_CTF_avg_UX_Y_pre));

if ~flag_found;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
X_best_ = zeros(n_iteration,1);
X_best_flag_flip_ = zeros(n_iteration,1);
polar_a_best_ = zeros(n_iteration,1);
azimu_b_best_ = zeros(n_iteration,1);
gamma_z_best_ = zeros(n_iteration,1);
delta_best__ = zeros(3,n_iteration);
%%%%%%%%;
for niteration=0:n_iteration-1;
%%%%%%%%;
N_wavelength = 0.0;
[ ...
 X_best ...
,flag_flip ...
,polar_a_best ...
,azimu_b_best ...
,gamma_z_best ...
,delta_best_ ...
,~ ...
,~ ...
] = ...
register_spharm_to_spharm_wigner_wrap_1( ...
 pm_n_k_p_r ...
,pm_k_p_r_ ...
,pm_k_p_r_max ...
,pm_weight_3d_k_p_r_ ...
,N_wavelength ...
,pm_l_max_ ...
,a_CTF_avg_UX_Y_true_ ...
,a_CTF_avg_UX_Y_calc__(:,1+niteration) ...
);
X_best_(1+niteration) = X_best;
X_best_flag_flip_(1+niteration) = flag_flip;
polar_a_best_(1+niteration) = polar_a_best;
azimu_b_best_(1+niteration) = azimu_b_best;
gamma_z_best_(1+niteration) = gamma_z_best;
delta_best__(:,1+niteration) = delta_best_;
if (verbose); disp(sprintf(' %% niteration %d/%d: X_best %0.4f X_best_flag_flip %d',niteration,n_iteration,X_best,flag_flip)); end;
end;%for niteration=0:n_iteration-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if ~flag_found;

flag_found=0;
if (~isempty(parameter.fname_align_a_CTF_avg_UX_Y_pre));
tmp_fname_mat = sprintf('%s.mat',parameter.fname_align_a_CTF_avg_UX_Y_pre);
if (~exist(tmp_fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',tmp_fname_mat));
save(tmp_fname_mat ...
     ,'X_best_','X_best_flag_flip_' ...
     ,'polar_a_best_' ...
     ,'azimu_b_best_' ...
     ,'gamma_z_best_' ...
     ,'delta_best__' ...
    );
end;%if (~exist(tmp_fname_mat,'file'));
flag_found=1;
end;%if (~isempty(parameter.fname_align_a_CTF_avg_UX_Y_pre));

if flag_found;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (~isempty(parameter.fname_align_a_CTF_avg_UX_Y_pre));
tmp_fname_fig_jpg = sprintf('%s.jpg',parameter.fname_align_a_CTF_avg_UX_Y_pre);
if (~exist(tmp_fname_fig_jpg,'file'));
if (verbose); disp(sprintf(' %% %s not found, creating',tmp_fname_fig_jpg)); end;
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
end;%if (~isempty(parameter.fname_align_a_CTF_avg_UX_Y_pre));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_found;

if (verbose); disp(sprintf(' %% [finished ampmut_align_to_a_CTF_avg_UX_Y_0]')); end;
