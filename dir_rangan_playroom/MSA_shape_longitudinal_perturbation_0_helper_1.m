%%%%%%%%;
% longitudinal perturbation. ;
%%%%%%%%;
gamma_z_ = linspace(0,2*pi,n_w_max+1); gamma_z_ = transpose(gamma_z_(1:n_w_max));
%n_mode = 4;
g_dilate = @(gamma) sin(n_mode*gamma);
T_max = 2.0/n_mode; n_t = 16; dt = T_max/max(1,n_t);
f_dilate = @(gamma) gamma + dt*g_dilate(gamma);
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
tmp_O_k_p_pre_ = tmp_N_k_p_form_;
%%%%%%%%%%%%%%%%;
for nt=0:n_t-1;
%%%%%%%%%%%%%%%%;
tmp_O_k_q_pre_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,tmp_O_k_p_pre_);
tmp_O_k_p_pos_ = reshape(tmp_p_pos_from_q_wq__*reshape(tmp_O_k_q_pre_,[n_w_max,n_k_p_r]),[n_w_sum,1]);
tmp_O_k_p_pre_ = tmp_O_k_p_pos_;
tmp_R_x_p_pre_ = ...
radon_k_p_to_x_p_xxnufft( ...
 n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,tmp_O_k_p_pre_ ...
);
%%%%%%%%;
flag_disp=1;
%%%%%%%%;
if flag_disp;
%%%%%%%%;
n_contour = 32;
N_k_p_form_lim_ = prctile(abs(real(tmp_N_k_p_form_)), 95,'all')*[-1,+1];
N_x_c_reco_lim_ = prctile(abs(real(tmp_N_x_c_reco_)),100,'all')*[-1,+1];
R_x_p_form_lim_ = prctile(abs(real(tmp_R_x_p_form_)), 95,'all')*[-1,+1];
N_x_c_reco_val_ = linspace(min(N_x_c_reco_lim_),max(N_x_c_reco_lim_),n_contour);
tmp_O_x_c_pre_ = interp_k_p_to_x_c_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_O_k_p_pre_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_c^2) * n_w_sum;
figure(1);clf;figbig;
p_row = 2; p_col = 4; np=0;
subplot_{1+np} = subplot(p_row,p_col,1+np);np=np+1; imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(tmp_N_x_c_reco_),N_x_c_reco_lim_);axis image;axisnotick;
subplot_{1+np} = subplot(p_row,p_col,1+np);np=np+1; contour(real(transpose(tmp_N_x_c_reco_)),N_x_c_reco_val_);axis image;axisnotick;
subplot_{1+np} = subplot(p_row,p_col,1+np);np=np+1; imagesc_p(n_k_p_r,k_p_r_/k_p_r_max,n_w_,n_w_sum,real(tmp_R_x_p_form_),R_x_p_form_lim_,colormap_beach);axis image;axisnotick;
subplot_{1+np} = subplot(p_row,p_col,1+np);np=np+1; imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(tmp_N_k_p_form_),N_k_p_form_lim_,colormap_80s);axis image;axisnotick;
subplot_{1+np} = subplot(p_row,p_col,1+np);np=np+1; imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(tmp_O_x_c_pre_),N_x_c_reco_lim_);axis image;axisnotick;
subplot_{1+np} = subplot(p_row,p_col,1+np);np=np+1; contour(real(transpose(tmp_O_x_c_pre_)),N_x_c_reco_val_);axis image;axisnotick;
subplot_{1+np} = subplot(p_row,p_col,1+np);np=np+1; imagesc_p(n_k_p_r,k_p_r_/k_p_r_max,n_w_,n_w_sum,real(tmp_R_x_p_pre_),R_x_p_form_lim_,colormap_beach);axis image;axisnotick;
subplot_{1+np} = subplot(p_row,p_col,1+np);np=np+1; imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(tmp_O_k_p_pre_),N_k_p_form_lim_,colormap_80s);axis image;axisnotick;
sgtitle(sprintf('nt %d/%d',nt,n_t),'Interpreter','none');
np=0;
colormap(subplot_{1+np},colormap_beach()); np=np+1;
colormap(subplot_{1+np},1-colormap('gray')); np=np+1;
colormap(subplot_{1+np},colormap_beach()); np=np+1;
colormap(subplot_{1+np},colormap_80s()); np=np+1;
colormap(subplot_{1+np},colormap_beach()); np=np+1;
colormap(subplot_{1+np},1-colormap('gray')); np=np+1;
colormap(subplot_{1+np},colormap_beach()); np=np+1;
colormap(subplot_{1+np},colormap_80s()); np=np+1;
drawnow(); %error('stop');
fname_fig_pre = sprintf('%s/MSA_shape_longitudinal_perturbation_%s_ns%.2dnm%.2dnt%.2d',dir_manuscript_jpg,str_shape,n_source_gaussian,n_mode,nt);
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
%%%%%%%%;
end;%if flag_disp;
%%%%%%%%;
%%%%%%%%%%%%%%%%;
end;%for nt=0:n_t-1;
%%%%%%%%%%%%%%%%;
