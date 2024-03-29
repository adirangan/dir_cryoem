tmp_flag_calc = 1;
sigma_psi = pi/32;
n_psi = 32;
psi_p_ = transpose(linspace(-2*sigma_psi,+2*sigma_psi,n_psi));
psi_w_ = exp(-(psi_p_.^2)/(2*sigma_psi.^2));
psi_w_ = psi_w_/sum(psi_w_);
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
%%%%%%%%;
tmp_parameter = struct('type','parameter');
tmp_parameter.strip_width = strip_width_use*min(n_omega_)/n_omega;
tmp_parameter.strip_width = strip_width_use/8;
tmp_parameter.k_0_offset = tmp_x_0_offset;
tmp_parameter.k_1_offset = tmp_x_1_offset + 1.0*dx_1;
tmp_parameter.clim_ = R_x_p_form_lim_;
for npsi=0:n_psi-1;
tmp_parameter.gamma_offset = psi_p_(1+npsi);
imagesc_p_strip_0( ...
 tmp_parameter ...
,n_x_p_r ...
,x_p_r_ ...
,x_p_r_max ...
,n_w_ ...
,n_w_sum ...
,tmp_R_x_p_pre_*psi_w_(1+npsi)/max(psi_w_) ...
,n_omega ...
,tmp_omega_ ...
);
end;%for npsi=0:n_psi-1;
%%%%%%%%;
tmp_parameter = struct('type','parameter');
plot_arrow_0(tmp_parameter,[tmp_x_0_offset;tmp_x_1_offset + 0.5*dx_1],-pi/2,0.25);
%%%%%%%%;
tmp_S_x_p_pre_wko__ = zeros(n_w_sum,n_omega);
%%%%;
if tmp_flag_calc;
%%%%;
for nomega=0:n_omega-1;
tmp_omega = tmp_omega_(1+nomega);
tmp_n_0 = cos(tmp_omega);
tmp_n_1 = sin(tmp_omega);
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
for npsi=0:n_psi-1;
tmp_S_x_p_pre_wko__(:,1+nomega) = tmp_S_x_p_pre_wko__(:,1+nomega) + psi_w_(1+npsi)*rotate_p_to_p_fftw(n_x_p_r,n_w_,n_w_sum,tmp_S_x_p_pre_,psi_p_(1+npsi));
end;%for npsi=0:n_psi-1;
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
tmp_parameter.strip_width = strip_width_use/8;
tmp_parameter.k_0_offset = tmp_x_0_offset;
tmp_parameter.k_1_offset = tmp_x_1_offset - 1.0*dx_1;
tmp_parameter.clim_ = N_k_p_form_lim_;
tmp_parameter.c_use__ = colormap_80s;
for npsi=0:n_psi-1;
tmp_parameter.gamma_offset = psi_p_(1+npsi);
imagesc_p_strip_0( ...
 tmp_parameter ...
,n_x_p_r ...
,x_p_r_ ...
,x_p_r_max ...
,n_w_ ...
,n_w_sum ...
,tmp_N_k_p_form_*psi_w_(1+npsi)/max(psi_w_) ...
,n_omega ...
,tmp_omega_ ...
);
end;%for npsi=0:n_psi-1;
%%%%%%%%;
tmp_x_0_offset = tmp_x_0_offset + dx_0;
%%%%%%%%;
end;%for n_omega = n_omega_;
%%%%%%%%;
tmp_T_x_p_pre_ = zeros(n_w_sum,1);
for npsi=0:n_psi-1;
tmp_T_x_p_pre_ = tmp_T_x_p_pre_ + psi_w_(1+npsi)*rotate_p_to_p_fftw(n_x_p_r,n_w_,n_w_sum,tmp_R_x_p_pre_,psi_p_(1+npsi));
end;%for npsi=0:n_psi-1;
tmp_parameter = struct('type','parameter');
tmp_parameter.strip_width = strip_width_use*min(n_omega_)/n_omega;
tmp_parameter.strip_width = strip_width_use/8;
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
,tmp_T_x_p_pre_ ...
);
%%%%%%%%;
tmp_parameter = struct('type','parameter');
plot_arrow_0(tmp_parameter,[tmp_x_0_offset;tmp_x_1_offset + 0.5*dx_1],-pi/2,0.25);
%%%%%%%%;
tmp_T_x_c_reco_ = zeros(n_x_c,n_x_c);
for npsi=0:n_psi-1;
tmp_T_x_c_reco_ = tmp_T_x_c_reco_ + psi_w_(1+npsi)*fftshift(rotate_c_to_c(n_x_c,x_p_r_max,n_x_c,x_p_r_max,fftshift(tmp_N_x_c_reco_),psi_p_(1+npsi)));
end;%for npsi=0:n_psi-1;
T_x_c_reco_lim_ = prctile(abs(real(tmp_T_x_c_reco_)),100,'all')*[-1,+1];
imagesc_c(n_x_c,x_c_0_+tmp_x_0_offset,n_x_c,x_c_1_+tmp_x_1_offset,real(tmp_T_x_c_reco_),T_x_c_reco_lim_);
%%%%%%%%;
tmp_parameter = struct('type','parameter');
plot_arrow_0(tmp_parameter,[tmp_x_0_offset;tmp_x_1_offset - 0.5*dx_1],-pi/2,0.25);
%%%%%%%%;
tmp_T_k_p_form_ = zeros(n_w_sum,1);
for npsi=0:n_psi-1;
tmp_T_k_p_form_ = tmp_T_k_p_form_ + psi_w_(1+npsi)*rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,tmp_N_k_p_form_,psi_p_(1+npsi));
end;%for npsi=0:n_psi-1;
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
,tmp_T_k_p_form_ ...
);
%%%%%%%%;
set(gcf,'Position',1+[0,0,1024,1024]);
hold off;
xlim([-1.25,+8.75]);
ylim([-5.00,+5.00]);
axis equal;
axisnotick;
drawnow(); %error('stop');
fname_fig_pre = sprintf('%s/MSA_shape_demonstration_%s_ns%.2d_reconstruction_FIGC',dir_manuscript_jpg,str_shape,n_source_gaussian);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if ( flag_replot | ~exist(fname_fig_jpg,'file') );
sgtitle(fname_fig_pre,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
sgtitle('');
disp(sprintf(' %% writing %s',sprintf('%s_strip.jpg',fname_fig_pre)));
print('-djpeg',sprintf('%s_strip.jpg',fname_fig_pre));
%close(gcf);
end;%if ( flag_replot | ~exist(fname_fig_jpg,'file') );
