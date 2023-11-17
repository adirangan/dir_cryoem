%%%%%%%%;
% Now rotate the volume. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/a_p312_x_u_base_.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
a_p312_x_u_base_ = permute(a_x_u_base_,[3,1,2]);
%%%%%%%%;
eta = pi/x_p_r_max; tmp_t = tic;
a_p312_k_p_quad_ = xxnufft3d3(n_xxx_u,x_u_0___(:)*eta,x_u_1___(:)*eta,x_u_2___(:)*eta,a_p312_x_u_base_(:).*xxx_u_weight_(:),-1,1e-12,n_k_all,2*pi*k_c_0_all_/eta,2*pi*k_c_1_all_/eta,2*pi*k_c_2_all_/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_p312_k_p_quad_ time %0.2fs',tmp_t));
%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
a_p312_x_u_reco_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,a_p312_k_p_quad_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_p312_x_u_reco_ time %0.2fs',tmp_t));
%%%%;
disp(sprintf(' %% xxnufft3d3: a_p312_x_u_reco error: %0.16f',fnorm(a_p312_x_u_base_(:)-a_p312_x_u_reco_)/fnorm(a_p312_x_u_base_(:))));
disp(sprintf(' %% at this point one should ensure that a_p312_k_p_quad_ on the outer shells (i.e., near k_p_r_max) has decayed to the desired precision.'));
%%%%%%%%;
tmp_t = tic;
[a_p312_k_Y_quad_] = convert_k_p_to_spharm_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_p312_k_p_quad_);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_p312_k_Y_quad_ time %0.2fs',tmp_t));
disp(sprintf(' %% Here we should ensure that a_p312_k_Y_ on the outer shells has decayed to the desired precision.'));
disp(sprintf(' %% Moreover, we should ensure that a_p312_k_Y_(l,m) has decayed for large l,m at each shell.'));
tmp_t = tic;
[a_p312_k_p_reco_] = convert_spharm_to_k_p_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_p312_k_Y_quad_);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_p312_k_Y_quad_ --> a_p312_k_p_reco_ time %0.2fs',tmp_t));
disp(sprintf(' %% xxnufft3d3: a_p312_k_p_reco error: %0.16f',fnorm(a_p312_k_p_quad_-a_p312_k_p_reco_)/fnorm(a_p312_k_p_quad_))); %<-- this should be 2-3 digits. ;
%%%%;
a_p312_k_Y_quad__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
a_p312_k_Y_quad__(1:n_lm_(1+nk_p_r),1+nk_p_r) = a_p312_k_Y_quad_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
tmp_t = tic();
[ ...
 S_p312_k_p__ ...
,~ ...
,tmp_n_S ...
,~ ...
,~ ...
] = ...
pm_template_2( ...
 verbose ...
,l_max_max...
,n_k_p_r ...
,reshape(a_p312_k_Y_quad__,[n_lm_max,n_k_p_r]) ...
,viewing_k_eq_d/k_p_r_max ...
,template_k_eq_d ...
,n_w_max ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% pm_template_2: %0.2fs',tmp_t));
S_p312_k_p__ = reshape(S_p312_k_p__,[n_w_max*n_k_p_r,tmp_n_S]);
%%%%%%%%;
save(fname_mat ...
     ,'a_p312_x_u_base_','a_p312_k_p_quad_','a_p312_x_u_reco_','a_p312_k_Y_quad_','a_p312_k_p_reco_','a_p312_k_Y_quad__','S_p312_k_p__' ...
     );
end;%if (flag_recalc | ~exist(fname_mat,'file'));
%%%%%%%%;
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;

%%%%%%%%;
% calculate norm of templates. ;
%%%%%%%%;
S_l2_S_ = zeros(n_S,1);
S_l2_S_ = sum(abs(S_k_p__).^2 .* weight_2d_k_all_,1);
[~,tmp_ij] = min(S_l2_S_); tmp_nS = tmp_ij-1;
S_k_p_ = S_k_p__(:,1+tmp_nS);
S_x_c_ = interp_k_p_to_x_c_xxnufft(n_x_u_pack,diameter_x_c,n_x_u_pack,diameter_x_c,n_k_p_r,k_p_r_,n_w_,S_k_p_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_u_pack^2) * n_w_sum;
%%%%%%%%;
S_p312_l2_S_ = zeros(n_S,1);
S_p312_l2_S_ = sum(abs(S_p312_k_p__).^2 .* weight_2d_k_all_,1);
[~,tmp_ij] = min(S_p312_l2_S_); tmp_nS = tmp_ij-1;
S_p312_k_p_ = S_p312_k_p__(:,1+tmp_nS);
S_p312_x_c_ = interp_k_p_to_x_c_xxnufft(n_x_u_pack,diameter_x_c,n_x_u_pack,diameter_x_c,n_k_p_r,k_p_r_,n_w_,S_p312_k_p_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_u_pack^2) * n_w_sum;
%%%%%%%%;
fname_fig_pre = sprintf('%s_jpg/a_p312_x_u_base_FIGB_',dir_pm);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_jpg_stripped = sprintf('%s_stripped.jpg',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_jpg));
%%%%;
figure(1+nf);nf=nf+1;clf;figmed;
fontsize_use = 12;
p_row = 2; p_col = 3; np=0;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
tmp_parameter = struct('type','parameter','percent_threshold_',99.5);
isosurface_f_x_u_1(tmp_parameter,a_x_u_base_); title('a_x_u_base_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
flag_2d_vs_3d = 0;
imagesc_polar_a_azimu_b_0( ...
 viewing_polar_a_all_ ... 
,viewing_azimu_b_all_ ... 
,S_l2_S_ ... 
,[] ... 
,colormap_81s ... 
,flag_2d_vs_3d ...
,k_p_r_max ...
);
set(gca,'XTick',[],'YTick',[],'ZTick',[]); xlabel('x0'); ylabel('x1'); zlabel('x2'); axis equal; axis vis3d;
title('S_l2_S_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
imagesc_c(n_x_u_pack,x_u_0_,n_x_u_pack,x_u_1_,fliplr(transpose(fliplr(real(S_x_c_)))),[],colormap_beach);
axisnotick; axis image;
title('S_x_c_ (min)','Interpreter','none');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
tmp_parameter = struct('type','parameter','percent_threshold_',99.5);
isosurface_f_x_u_1(tmp_parameter,a_p312_x_u_base_); title('a_p312_x_u_base_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
flag_2d_vs_3d = 0;
imagesc_polar_a_azimu_b_0( ...
 viewing_polar_a_all_ ... 
,viewing_azimu_b_all_ ... 
,S_p312_l2_S_ ... 
,[] ... 
,colormap_81s ... 
,flag_2d_vs_3d ...
,k_p_r_max ...
);
set(gca,'XTick',[],'YTick',[],'ZTick',[]); xlabel('x0'); ylabel('x1'); zlabel('x2'); axis equal; axis vis3d;
title('S_p312_l2_S_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
imagesc_c(n_x_u_pack,x_u_0_,n_x_u_pack,x_u_1_,fliplr(transpose(fliplr(real(S_p312_x_c_)))),[],colormap_beach);
axisnotick; axis image;
title('S_p312_x_c_ (min)','Interpreter','none');
%%%%;
disp(sprintf(' %% writing %s',fname_fig_jpg));
print('-djpeg',fname_fig_jpg_stripped);
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
close(gcf);
end;%if (~exist(fname_fig_jpg,'file'));
if ( exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s found, not creating',fname_fig_jpg));
end;%if ( exist(fname_fig_jpg,'file'));
%%%%%%%%;
