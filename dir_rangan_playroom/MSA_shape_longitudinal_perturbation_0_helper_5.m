%%%%%%%%;
% define dilation. ;
%%%%%%%%;
gamma_z_ = linspace(0,2*pi,n_w_max+1); gamma_z_ = transpose(gamma_z_(1:n_w_max));
g_dilate = @(gamma) sin(n_mode*gamma);
T_max = 2.0/n_mode; n_t = 16; dt = T_max/max(1,n_t);
f_dilate = @(gamma) gamma + dt*g_dilate(gamma);
%%%%%%%%;
gamma_z_pre_ = gamma_z_;
for nw=0:n_w_max-1;
gamma_z_pre_(1+nw) = fminsearch(@(gamma) abs(f_dilate(gamma)-gamma_z_(1+nw)),gamma_z_(1+nw),optimset('TolX',1e-6));
end;%for nw=0:n_w_max-1;
q_ = periodize(transpose(0:n_w_max-1),-n_w_max/2,+n_w_max/2);
%%%%;
tmp_p_from_q_wq__ = zeros(n_w_max,n_w_max);
for nq=0:n_w_max-1;
tmp_q = q_(1+nq);
tmp_p_from_q_wq__(:,1+nq) = exp(+i*gamma_z_*tmp_q)/sqrt(n_w_max);
end;%for nq=0:n_w_max-1;
tmp_N_k_q_form_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,tmp_N_k_p_form_);
tmp_N_k_p_reco_ = reshape(tmp_p_from_q_wq__*reshape(tmp_N_k_q_form_,[n_w_max,n_k_p_r]),[n_w_sum,1]);
disp(sprintf(' %% tmp_N_k_p_form_ vs tmp_N_k_p_reco_: %0.16f',fnorm(tmp_N_k_p_form_-tmp_N_k_p_reco_)/fnorm(tmp_N_k_p_form_)));
%%%%;
tmp_p_pos_from_q_wq__ = zeros(n_w_max,n_w_max);
for nq=0:n_w_max-1;
tmp_q = q_(1+nq);
tmp_p_pos_from_q_wq__(:,1+nq) = exp(+i*gamma_z_pre_*tmp_q)/sqrt(n_w_max);
end;%for nq=0:n_w_max-1;
%%%%;
%%%%%%%%;
% Now define tmp_O_k_p_pre_ using n_dt timesteps. ;
%%%%%%%%;
n_dt = n_dt_use;
tmp_O_k_p_pre_ = tmp_N_k_p_form_;
for nt=0:n_dt-1;
tmp_O_k_q_pre_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,tmp_O_k_p_pre_);
tmp_O_k_p_pos_ = reshape(tmp_p_pos_from_q_wq__*reshape(tmp_O_k_q_pre_,[n_w_max,n_k_p_r]),[n_w_sum,1]);
tmp_O_k_p_pre_ = tmp_O_k_p_pos_;
end;%for nt=0:n_t-1;
%%%%%%%%;
% Now define the final dilated picture. ;
%%%%%%%%;
tmp_R_x_p_pre_ = ...
radon_k_p_to_x_p_xxnufft( ...
 n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,tmp_O_k_p_pre_ ...
);
%%%%;
N_k_p_form_lim_ = prctile(abs(real(tmp_N_k_p_form_)), 95,'all')*[-1,+1];
N_x_c_reco_lim_ = prctile(abs(real(tmp_N_x_c_reco_)),100,'all')*[-1,+1];
R_x_p_form_lim_ = prctile(abs(real(tmp_R_x_p_form_)), 95,'all')*[-1,+1];
tmp_O_x_c_pre_ = interp_k_p_to_x_c_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_O_k_p_pre_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_c^2) * n_w_sum;
%%%%%%%%;

%%%%%%%%;
% Now collating figure: ;
%%%%%%%%;
tmp_flag_calc = 1;
%%%%%%%%%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figmed;
dx_1 = 3.5;
hold on;
%%%%%%%%;
tmp_x_0_offset = 0.0; tmp_x_1_offset = 0.0; dx_0 = 2.5;
n_omega_ = [3,5,7];
for n_omega = n_omega_;
tmp_omega_ = linspace(-pi/2,+pi/2,2*n_omega+1);
tmp_omega_ = transpose(tmp_omega_(2:2:end-1));
n_omega = numel(tmp_omega_);
tmp_omega_pos_ = tmp_omega_;
for ndt=0:n_dt-1; tmp_omega_pos_ = f_dilate(tmp_omega_pos_); end;%for ndt=0:n_dt-1;
%%%%%%%%;
tmp_parameter = struct('type','parameter');
tmp_parameter.strip_width = strip_width_use*min(n_omega_)/n_omega;
tmp_parameter.k_0_offset = tmp_x_0_offset;
tmp_parameter.k_1_offset = tmp_x_1_offset + 1.0*dx_1;
tmp_parameter.clim_ = R_x_p_form_lim_;
tmp_parameter.gamma_offset_ = tmp_omega_pos_ - tmp_omega_;
imagesc_p_strip_0( ...
 tmp_parameter ...
,n_x_p_r ...
,x_p_r_ ...
,x_p_r_max ...
,n_w_ ...
,n_w_sum ...
,tmp_R_x_p_pre_ ...
,n_omega ...
,tmp_omega_ ...
);
%%%%%%%%;
tmp_parameter = struct('type','parameter');
plot_arrow_0(tmp_parameter,[tmp_x_0_offset;tmp_x_1_offset + 0.5*dx_1],-pi/2,0.25);
%%%%%%%%;
tmp_S_x_p_pre_wko__ = zeros(n_w_sum,n_omega);
%%%%;
if tmp_flag_calc;
%%%%;
for nomega=0:n_omega-1;
tmp_omega_pos = tmp_omega_pos_(1+nomega);
tmp_n_0 = cos(tmp_omega_pos);
tmp_n_1 = sin(tmp_omega_pos);
tmp_x_c_xn_all_ = x_c_0_all_*tmp_n_0 + x_c_1_all_*tmp_n_1;
tmp_x_c_0_all_ = tmp_x_c_xn_all_*tmp_n_0;
tmp_x_c_1_all_ = tmp_x_c_xn_all_*tmp_n_1;
tmp_x_c_r_all_ = sqrt(tmp_x_c_0_all_.^2 + tmp_x_c_1_all_.^2);
tmp_x_c_w_all_ = atan2(tmp_x_c_1_all_,tmp_x_c_0_all_);
tmp_t = tic();
n_order = 7;
tmp_scatter_from_tensor_swk__ = ...
disk_k_p_scatter_from_tensor_interpolate_n_4( ...
 n_order ...
,n_w_max ...
,n_x_p_r ...
,x_p_r_ ...
,x_p_r_max ...
,n_w_sum ...
,tmp_x_c_w_all_ ...
,tmp_x_c_r_all_ ...
);
tmp_t = toc(tmp_t); if (flag_verbose> 0); disp(sprintf(' %% disk_k_p_scatter_from_tensor_interpolate_n_4: %0.2fs',tmp_t)); end;
tmp_S_x_p_pre_ = tmp_scatter_from_tensor_swk__*tmp_R_x_p_pre_;
tmp_S_x_p_pre_wko__(:,1+nomega) = tmp_S_x_p_pre_;
clear tmp_scatter_from_tensor_swk__ tmp_S_x_p_pre_;
end;%for nomega=0:n_omega-1;
%%%%;
end;%if tmp_flag_calc;
%%%%;
tmp_S_x_p_pre_sum_ = sum(tmp_S_x_p_pre_wko__,2);
S_x_p_form_lim_ = prctile(abs(real(tmp_S_x_p_pre_sum_)),100,'all')*[-1,+1];
tmp_parameter = struct('type','parameter');
tmp_parameter.strip_width = 1.0; %<-- not strip_width_use;
tmp_parameter.k_0_offset = tmp_x_0_offset;
tmp_parameter.k_1_offset = tmp_x_1_offset + 0*dx_1;
tmp_parameter.clim_ = S_x_p_form_lim_;
imagesc_p_strip_0( ...
 tmp_parameter ...
,n_x_p_r ...
,x_p_r_ ...
,x_p_r_max ...
,n_w_ ...
,n_w_sum ...
,tmp_S_x_p_pre_sum_ ...
);
%%%%%%%%;
tmp_parameter = struct('type','parameter');
plot_arrow_0(tmp_parameter,[tmp_x_0_offset;tmp_x_1_offset - 0.5*dx_1],-pi/2,0.25);
%%%%%%%%;
tmp_parameter = struct('type','parameter');
tmp_parameter.strip_width = strip_width_use*min(n_omega_)/n_omega;
tmp_parameter.k_0_offset = tmp_x_0_offset;
tmp_parameter.k_1_offset = tmp_x_1_offset - 1.0*dx_1;
tmp_parameter.clim_ = N_k_p_form_lim_;
tmp_parameter.c_use__ = colormap_80s;
tmp_parameter.gamma_offset_ = tmp_omega_pos_ - tmp_omega_;
imagesc_p_strip_0( ...
 tmp_parameter ...
,n_x_p_r ...
,x_p_r_ ...
,x_p_r_max ...
,n_w_ ...
,n_w_sum ...
,tmp_O_k_p_pre_ ...
,n_omega ...
,tmp_omega_ ...
);
%%%%%%%%;
tmp_x_0_offset = tmp_x_0_offset + dx_0;
%%%%%%%%;
end;%for n_omega = n_omega_;
%%%%%%%%;
tmp_parameter = struct('type','parameter');
tmp_parameter.strip_width = strip_width_use*min(n_omega_)/n_omega;
tmp_parameter.k_0_offset = tmp_x_0_offset;
tmp_parameter.k_1_offset = tmp_x_1_offset + 1.0*dx_1;
tmp_parameter.clim_ = R_x_p_form_lim_;
imagesc_p_strip_0( ...
 tmp_parameter ...
,n_x_p_r ...
,x_p_r_ ...
,x_p_r_max ...
,n_w_ ...
,n_w_sum ...
,tmp_R_x_p_pre_ ...
);
%%%%%%%%;
tmp_parameter = struct('type','parameter');
plot_arrow_0(tmp_parameter,[tmp_x_0_offset;tmp_x_1_offset + 0.5*dx_1],-pi/2,0.25);
%%%%%%%%;
imagesc_c(n_x_c,x_c_0_+tmp_x_0_offset,n_x_c,x_c_1_+tmp_x_1_offset,real(tmp_O_x_c_pre_),N_x_c_reco_lim_);
%%%%%%%%;
tmp_parameter = struct('type','parameter');
plot_arrow_0(tmp_parameter,[tmp_x_0_offset;tmp_x_1_offset - 0.5*dx_1],-pi/2,0.25);
%%%%%%%%;
tmp_parameter = struct('type','parameter');
tmp_parameter.strip_width = 0; %<-- not use strip_width_use. ;
tmp_parameter.k_0_offset = tmp_x_0_offset;
tmp_parameter.k_1_offset = tmp_x_1_offset - 1.0*dx_1;
tmp_parameter.clim_ = N_k_p_form_lim_;
tmp_parameter.c_use__ = colormap_80s;
imagesc_p_strip_0( ...
 tmp_parameter ...
,n_x_p_r ...
,x_p_r_ ...
,x_p_r_max ...
,n_w_ ...
,n_w_sum ...
,tmp_O_k_p_pre_ ...
);
%%%%%%%%;
set(gcf,'Position',1+[0,0,1024,1024]);
hold off;
xlim([-1.25,+8.75]);
ylim([-5.00,+5.00]);
axis equal;
axisnotick;
drawnow(); %error('stop');
fname_fig_pre = sprintf('%s/MSA_shape_demonstration_%s_ns%.2d_nt%.2d_ndt%.2d_reconstruction_FIGD',dir_manuscript_jpg,str_shape,n_source_gaussian,n_t,n_dt);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if ( flag_replot | ~exist(fname_fig_jpg,'file') );
sgtitle(fname_fig_pre,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
sgtitle('');
disp(sprintf(' %% writing %s',sprintf('%s_strip.jpg',fname_fig_pre)));
print('-djpeg',sprintf('%s_strip.jpg',fname_fig_pre));
close(gcf);
end;%if ( flag_replot | ~exist(fname_fig_jpg,'file') );

%{
figure(1);clf;figbig;
p_row = 2; p_col = 3; np=0;
subplot_{1+np} = subplot(p_row,p_col,1+np);np=np+1; imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(tmp_N_x_c_reco_),N_x_c_reco_lim_);axis image;axisnotick;
subplot_{1+np} = subplot(p_row,p_col,1+np);np=np+1; imagesc_p(n_k_p_r,k_p_r_/k_p_r_max,n_w_,n_w_sum,real(tmp_R_x_p_form_),R_x_p_form_lim_,colormap_beach);axis image;axisnotick;
subplot_{1+np} = subplot(p_row,p_col,1+np);np=np+1; imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(tmp_N_k_p_form_),N_k_p_form_lim_,colormap_80s);axis image;axisnotick;
subplot_{1+np} = subplot(p_row,p_col,1+np);np=np+1; imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(tmp_O_x_c_pre_),N_x_c_reco_lim_);axis image;axisnotick;
subplot_{1+np} = subplot(p_row,p_col,1+np);np=np+1; imagesc_p(n_k_p_r,k_p_r_/k_p_r_max,n_w_,n_w_sum,real(tmp_R_x_p_pre_),R_x_p_form_lim_,colormap_beach);axis image;axisnotick;
subplot_{1+np} = subplot(p_row,p_col,1+np);np=np+1; imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(tmp_O_k_p_pre_),N_k_p_form_lim_,colormap_80s);axis image;axisnotick;
sgtitle(sprintf('nt %d/%d',nt,n_t),'Interpreter','none');
np=0;
colormap(subplot_{1+np},colormap_beach()); np=np+1;
colormap(subplot_{1+np},colormap_beach()); np=np+1;
colormap(subplot_{1+np},colormap_80s()); np=np+1;
colormap(subplot_{1+np},colormap_beach()); np=np+1;
colormap(subplot_{1+np},colormap_beach()); np=np+1;
colormap(subplot_{1+np},colormap_80s()); np=np+1;
 %}
