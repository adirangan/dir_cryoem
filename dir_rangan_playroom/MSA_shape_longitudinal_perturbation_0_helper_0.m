%%%%%%%%;
% visualize shadowsphere. ;
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for ntype=0:2-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if ntype==0;
str_ntype = 'shadowsphere';
end;%if ntype==0;
if ntype==1;
str_ntype = 'fourierspace';
end;%if ntype==1;
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
tmp_R_x_p_pre_ = ...
radon_k_p_to_x_p_xxnufft( ...
 n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,tmp_O_k_p_pre_ ...
);
flag_disp=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if flag_disp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
N_k_p_form_lim_ = prctile(abs(real(tmp_N_k_p_form_)), 95,'all')*[-1,+1];
N_x_c_reco_lim_ = prctile(abs(real(tmp_N_x_c_reco_)),100,'all')*[-1,+1];
R_x_p_form_lim_ = prctile(abs(real(tmp_R_x_p_form_)), 95,'all')*[-1,+1];
tmp_O_x_c_pre_ = interp_k_p_to_x_c_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_O_k_p_pre_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_c^2) * n_w_sum;
figure(1+nf);nf=nf+1;clf;figbig; strip_width_use = 1/16;
%%%%%%%%%%%%%%%%;
subplot(1,2,1); hold on;
%%%%%%%%;
imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(tmp_N_x_c_reco_),N_x_c_reco_lim_);
%%%%%%%%;
n_l = 16;
for nl=flip(0:n_l-1);
tmp_omega = (2*pi)*nl/n_l;
tmp_x_0_offset = sqrt(1.0+0.5*mod(nl,2))*diameter_x_c*cos(tmp_omega+pi/2);
tmp_x_1_offset = sqrt(1.0+0.5*mod(nl,2))*diameter_x_c*sin(tmp_omega+pi/2);
tmp_parameter = struct('type','parameter');
%tmp_parameter.arrow_head_length = 1.0+2.0*mod(nl,2);
%tmp_parameter.arrow_tail_length = 2.5+1.0*mod(nl,2);
tmp_parameter.arrow_head_length = 2.25+1.0*mod(nl,2);
tmp_parameter.arrow_tail_length = 1.0+2.0*mod(nl,2);
plot_arrow_0(tmp_parameter,[tmp_x_0_offset;tmp_x_1_offset],tmp_omega+pi/2,0.25+0.00*mod(nl,2));
tmp_x_0_offset = sqrt(2.0+1.0*mod(nl,2))*diameter_x_c*cos(tmp_omega+pi/2); tmp_x_1_offset = sqrt(2.0+1.0*mod(nl,2))*diameter_x_c*sin(tmp_omega+pi/2);
tmp_parameter = struct('type','parameter');
tmp_parameter.strip_width = strip_width_use;
tmp_parameter.k_0_offset = tmp_x_0_offset;
tmp_parameter.k_1_offset = tmp_x_1_offset;
tmp_parameter.clim_ = R_x_p_form_lim_;
imagesc_p_strip_0( ...
 tmp_parameter ...
,n_k_p_r ...
,k_p_r_*x_p_r_max/k_p_r_max ...
,k_p_r_max*x_p_r_max/k_p_r_max ...
,n_w_ ...
,n_w_sum ...
,tmp_R_x_p_pre_ ...
,1 ...
,tmp_omega ...
);
tmp_x_0_offset = (sqrt(2.0+1.0*mod(nl,2))*diameter_x_c + 2*strip_width_use)*cos(tmp_omega+pi/2); tmp_x_1_offset = (sqrt(2.0+1.0*mod(nl,2))*diameter_x_c + 2*strip_width_use)*sin(tmp_omega+pi/2);
tmp_parameter = struct('type','parameter');
tmp_parameter.strip_width = strip_width_use;
tmp_parameter.k_0_offset = tmp_x_0_offset;
tmp_parameter.k_1_offset = tmp_x_1_offset;
tmp_parameter.clim_ = N_k_p_form_lim_;
tmp_parameter.c_use__ = colormap_80s;
imagesc_p_strip_0( ...
 tmp_parameter ...
,n_k_p_r ...
,k_p_r_*x_p_r_max/k_p_r_max ...
,k_p_r_max*x_p_r_max/k_p_r_max ...
,n_w_ ...
,n_w_sum ...
,tmp_O_k_p_pre_ ...
,1 ...
,tmp_omega ...
);
end;%for nl=flip(0:n_l-1);
%%%%%%%%;
hold off;
xlim([-4,+4]);
ylim([-4,+4]);
axis square;
axisnotick;
%%%%%%%%%%%%%%%%;
subplot(1,2,2); hold on;
%%%%%%%%;
if ntype==0;
tmp_data_ = tmp_R_x_p_pre_;
tmp_clim_ = R_x_p_form_lim_;
tmp_c_use__ = colormap_beach();
end;%if ntype==0;
if ntype==1;
tmp_data_ = tmp_N_k_p_form_;
tmp_clim_ = N_k_p_form_lim_;
tmp_c_use__ = colormap_80s();
end;%if ntype==1;
n_l = 16;
for nl=flip(0:n_l-1);
tmp_omega = (2*pi)*nl/n_l;
tmp_x_0_offset = sqrt(1.0+1.0*mod(nl,2))*diameter_x_c*cos(tmp_omega+pi/2); tmp_x_1_offset = sqrt(1.0+1.0*mod(nl,2))*diameter_x_c*sin(tmp_omega+pi/2);
tmp_x_0_offset = sqrt(1.0+0.5*mod(nl,2))*diameter_x_c*cos(tmp_omega+pi/2);
tmp_x_1_offset = sqrt(1.0+0.5*mod(nl,2))*diameter_x_c*sin(tmp_omega+pi/2);
tmp_parameter = struct('type','parameter');
tmp_parameter.arrow_head_length = 1.0+2.0*mod(nl,2);
tmp_parameter.arrow_tail_length = 2.5+1.0*mod(nl,2);
plot_arrow_0(tmp_parameter,[tmp_x_0_offset;tmp_x_1_offset],tmp_omega+pi/2+pi,0.25+0.00*mod(nl,2));
tmp_x_0_offset = sqrt(2.0+1.0*mod(nl,2))*diameter_x_c*cos(tmp_omega+pi/2); tmp_x_1_offset = sqrt(2.0+1.0*mod(nl,2))*diameter_x_c*sin(tmp_omega+pi/2);
tmp_parameter = struct('type','parameter');
tmp_parameter.strip_width = strip_width_use;
tmp_parameter.k_0_offset = tmp_x_0_offset;
tmp_parameter.k_1_offset = tmp_x_1_offset;
tmp_parameter.clim_ = tmp_clim_;
tmp_parameter.c_use__ = tmp_c_use__;
imagesc_p_strip_0( ...
 tmp_parameter ...
,n_k_p_r ...
,k_p_r_*x_p_r_max/k_p_r_max ...
,k_p_r_max*x_p_r_max/k_p_r_max ...
,n_w_ ...
,n_w_sum ...
,tmp_data_ ...
,1 ...
,tmp_omega ...
);
end;%for nl=flip(0:n_l-1);
if ntype==0; tmp_gamma_z_ = linspace(-pi/2,+pi/2,1+n_l/2); end;
if ntype==1; tmp_gamma_z_ = linspace(0,pi,1+n_l/2); end;
tmp_parameter = struct('type','parameter');
tmp_parameter.strip_width = strip_width_use;
tmp_parameter.clim_ = tmp_clim_;
tmp_parameter.c_use__ = tmp_c_use__;
imagesc_p_strip_0( ...
 tmp_parameter ...
,n_k_p_r ...
,k_p_r_*x_p_r_max/k_p_r_max ...
,k_p_r_max*x_p_r_max/k_p_r_max ...
,n_w_ ...
,n_w_sum ...
,tmp_data_ ...
,1+n_l/2 ...
,tmp_gamma_z_ ...
);
%%%%%%%%;
hold off;
xlim([-4,+4]);
ylim([-4,+4]);
axis square;
axisnotick;
%%%%%%%%%%%%%%%%;
drawnow(); %error('stop');
fname_fig_pre = sprintf('%s/MSA_shape_demonstration_%s_ns%.2d_%s',dir_manuscript_jpg,str_shape,n_source_gaussian,str_ntype);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_disp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for ntype=0:2-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

figure(1+nf);nf=nf+1;clf;figmed;
hold on;
%%%%%%%%%%%%%%%%;
for ntype=0:2-1;
if ntype==0;
tmp_data_ = tmp_R_x_p_pre_;
tmp_clim_ = R_x_p_form_lim_;
tmp_c_use__ = colormap_beach();
end;%if ntype==0;
if ntype==1;
tmp_data_ = tmp_N_k_p_form_;
tmp_clim_ = N_k_p_form_lim_;
tmp_c_use__ = colormap_80s();
end;%if ntype==1;
%%%%%%%%;
tmp_x_0_offset = 0.0; dx_0 = sqrt(2.0);
if ntype==0; tmp_x_1_offset = +1.125; dx_1 = sqrt(0.0); end;
if ntype==1; tmp_x_1_offset = -1.125; dx_1 = sqrt(0.0); end;
%%%%;
for np=0:2;
n_l = 8*2.^(np);
if ntype==0; tmp_gamma_z_ = linspace(-pi/2,+pi/2,1+n_l/2); end;
if ntype==1; tmp_gamma_z_ = linspace(0,pi,1+n_l/2); end;
tmp_parameter = struct('type','parameter');
tmp_parameter.strip_width = strip_width_use/2.^(np);
tmp_parameter.k_0_offset = tmp_x_0_offset;
tmp_parameter.k_1_offset = tmp_x_1_offset;
tmp_parameter.clim_ = tmp_clim_;
tmp_parameter.c_use__ = tmp_c_use__;
imagesc_p_strip_0( ...
 tmp_parameter ...
,n_k_p_r ...
,k_p_r_*x_p_r_max/k_p_r_max ...
,k_p_r_max*x_p_r_max/k_p_r_max ...
,n_w_ ...
,n_w_sum ...
,tmp_data_ ...
,1+n_l/2 ...
,tmp_gamma_z_ ...
);
tmp_x_0_offset = tmp_x_0_offset + dx_0;
tmp_x_1_offset = tmp_x_1_offset + dx_1;
tmp_parameter = struct('type','parameter');
plot_arrow_0(tmp_parameter,[tmp_x_0_offset;tmp_x_1_offset],0,0.25);
tmp_x_0_offset = tmp_x_0_offset + dx_0;
tmp_x_1_offset = tmp_x_1_offset + dx_1;
end;%for np=0:3;
%%%%;
tmp_parameter = struct('type','parameter');
tmp_parameter.strip_width = strip_width_use;
tmp_parameter.k_0_offset = tmp_x_0_offset;
tmp_parameter.k_1_offset = tmp_x_1_offset;
tmp_parameter.clim_ = tmp_clim_;
tmp_parameter.c_use__ = tmp_c_use__;
imagesc_p_strip_0( ...
 tmp_parameter ...
,n_k_p_r ...
,k_p_r_*x_p_r_max/k_p_r_max ...
,k_p_r_max*x_p_r_max/k_p_r_max ...
,n_w_ ...
,n_w_sum ...
,tmp_data_ ...
);
%%%%%%%%;
end;%for ntype=0:2-1;
%%%%%%%%%%%%%%%%;
hold off;
xlim([-1.25,+9.75]);
ylim([-1.00,+1.00]);
axis equal;
axisnotick;
drawnow(); %error('stop');
fname_fig_pre = sprintf('%s/MSA_shape_demonstration_%s_ns%.2d_disc',dir_manuscript_jpg,str_shape,n_source_gaussian);
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




