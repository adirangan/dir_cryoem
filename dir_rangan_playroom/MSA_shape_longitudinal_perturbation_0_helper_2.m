%%%%%%%%;
% visualize reconstruction. ;
%%%%%%%%;

gamma_z_ = linspace(0,2*pi,n_w_max+1); gamma_z_ = transpose(gamma_z_(1:n_w_max));
q_ = periodize(transpose(0:n_w_max-1),-n_w_max/2,+n_w_max/2);
%%%%%%%%;
tmp_p_from_q_wq__ = zeros(n_w_max,n_w_max);
for nq=0:n_w_max-1;
tmp_q = q_(1+nq);
tmp_p_from_q_wq__(:,1+nq) = exp(+i*gamma_z_*tmp_q)/sqrt(n_w_max);
end;%for nq=0:n_w_max-1;
tmp_N_k_q_form_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,tmp_N_k_p_form_);
tmp_N_k_p_reco_ = reshape(tmp_p_from_q_wq__*reshape(tmp_N_k_q_form_,[n_w_max,n_k_p_r]),[n_w_sum,1]);
disp(sprintf(' %% tmp_N_k_p_form_ vs tmp_N_k_p_reco_: %0.16f',fnorm(tmp_N_k_p_form_-tmp_N_k_p_reco_)/fnorm(tmp_N_k_p_form_)));
%%%%;
tmp_O_k_p_pre_ = tmp_N_k_p_form_;
tmp_O_k_q_pre_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,tmp_O_k_p_pre_);
%%%%;
n_x_p_r = n_k_p_r;
x_p_r_ = k_p_r_*x_p_r_max/k_p_r_max;
x_c_0_all_ = k_c_0_all_*x_p_r_max/k_p_r_max;
x_c_1_all_ = k_c_1_all_*x_p_r_max/k_p_r_max;
x_c_r_all_ = sqrt(x_c_0_all_.^2 + x_c_1_all_.^2);
x_c_w_all_ = atan2(x_c_1_all_,x_c_0_all_);
%%%%;
tmp_R_x_p_pre_ = ...
radon_k_p_to_x_p_xxnufft( ...
 n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,tmp_O_k_p_pre_ ...
);
%%%%%%%%;

%%%%%%%%%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figmed;
dx_1 = 2.0;
hold on;
%%%%%%%%;
tmp_x_0_offset = 0.0; tmp_x_1_offset = 0.0; dx_0 = 1.5;
n_omega=3; tmp_omega_ = transpose(linspace(-pi/3,+pi/3,n_omega));
%%%%%%%%;
tmp_parameter = struct('type','parameter');
tmp_parameter.strip_width = strip_width_use;
tmp_parameter.k_0_offset = tmp_x_0_offset - 2*dx_0;
tmp_parameter.k_1_offset = tmp_x_1_offset;
tmp_parameter.clim_ = R_x_p_form_lim_;
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
%%%%;
for nomega=0:n_omega-1;
tmp_omega = tmp_omega_(1+nomega);
tmp_parameter = struct('type','parameter');
plot_arrow_0(tmp_parameter,[tmp_x_0_offset - 1*dx_0;tmp_x_1_offset - dx_1 + dx_1*nomega],+tmp_omega/2,0.25);
end;%for nomega=0:n_omega-1;
%%%%%%%%;
tmp_S_x_p_pre_wko__ = zeros(n_w_sum,n_omega);
for nomega=0:n_omega-1;
tmp_omega = tmp_omega_(1+nomega);
tmp_n_0 = cos(tmp_omega);
tmp_n_1 = sin(tmp_omega);
tmp_x_c_xn_all_ = x_c_0_all_*tmp_n_0 + x_c_1_all_*tmp_n_1;
tmp_x_c_0_all_ = tmp_x_c_xn_all_*tmp_n_0;
tmp_x_c_1_all_ = tmp_x_c_xn_all_*tmp_n_1;
tmp_x_c_r_all_ = sqrt(tmp_x_c_0_all_.^2 + tmp_x_c_1_all_.^2);
tmp_x_c_w_all_ = atan2(tmp_x_c_1_all_,tmp_x_c_0_all_);
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
tmp_S_x_p_pre_ = tmp_scatter_from_tensor_swk__*tmp_R_x_p_pre_;
tmp_S_x_p_pre_wko__(:,1+nomega) = tmp_S_x_p_pre_;
clear tmp_scatter_from_tensor_swk__ tmp_S_x_p_pre_;
end;%for nomega=0:n_omega-1;
%%%%%%%%;
for nomega=0:n_omega-1;
tmp_omega = tmp_omega_(1+nomega);
tmp_parameter = struct('type','parameter');
tmp_parameter.strip_width = strip_width_use;
tmp_parameter.k_0_offset = tmp_x_0_offset;
tmp_parameter.k_1_offset = tmp_x_1_offset - dx_1 + dx_1*nomega;
tmp_parameter.clim_ = R_x_p_form_lim_;
imagesc_p_strip_0( ...
 tmp_parameter ...
,n_x_p_r ...
,x_p_r_ ...
,x_p_r_max ...
,n_w_ ...
,n_w_sum ...
,tmp_R_x_p_pre_ ...
,1 ...
,tmp_omega ...
);
end;%for nomega=0:n_omega-1;
%%%%%%%%;
tmp_x_0_offset = tmp_x_0_offset + dx_0;
%%%%%%%%;
for nomega=0:n_omega-1;
tmp_parameter = struct('type','parameter');
plot_arrow_0(tmp_parameter,[tmp_x_0_offset;tmp_x_1_offset - dx_1 + dx_1*nomega],0,0.25);
end;%for nomega=0:n_omega-1;
%%%%%%%%;
tmp_x_0_offset = tmp_x_0_offset + dx_0;
%%%%%%%%;
for nomega=0:n_omega-1;
tmp_parameter = struct('type','parameter');
tmp_parameter.strip_width = 1.0; %<-- not strip_width_use;
tmp_parameter.k_0_offset = tmp_x_0_offset;
tmp_parameter.k_1_offset = tmp_x_1_offset - dx_1 + dx_1*nomega;
tmp_parameter.clim_ = R_x_p_form_lim_;
imagesc_p_strip_0( ...
 tmp_parameter ...
,n_x_p_r ...
,x_p_r_ ...
,x_p_r_max ...
,n_w_ ...
,n_w_sum ...
,tmp_S_x_p_pre_wko__(:,1+nomega) ...
);
end;%for nomega=0:n_omega-1;
%%%%%%%%;
tmp_x_0_offset = tmp_x_0_offset + dx_0;
%%%%%%%%;
for nomega=0:n_omega-1;
tmp_omega = tmp_omega_(1+nomega);
tmp_parameter = struct('type','parameter');
plot_arrow_0(tmp_parameter,[tmp_x_0_offset - 0*dx_0;tmp_x_1_offset - dx_1 + dx_1*nomega],-tmp_omega/2,0.25);
end;%for nomega=0:n_omega-1;
%%%%%%%%;
tmp_x_0_offset = tmp_x_0_offset + dx_0;
%%%%%%%%;
tmp_S_x_p_pre_sum_ = sum(tmp_S_x_p_pre_wko__,2);
S_x_p_form_lim_ = prctile(abs(real(tmp_S_x_p_pre_sum_)),100,'all')*[-1,+1];
tmp_parameter = struct('type','parameter');
tmp_parameter.strip_width = 1.0; %<-- not strip_width_use;
tmp_parameter.k_0_offset = tmp_x_0_offset;
tmp_parameter.k_1_offset = tmp_x_1_offset;
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
hold off;
xlim([-4.25,+7.25]);
ylim([-1.00,+1.00]);
axis equal;
axisnotick;
drawnow(); %error('stop');
fname_fig_pre = sprintf('%s/MSA_shape_demonstration_%s_ns%.2d_reconstruction_FIGA',dir_manuscript_jpg,str_shape,n_source_gaussian);
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
%%%%%%%%%%%%%%%%;



