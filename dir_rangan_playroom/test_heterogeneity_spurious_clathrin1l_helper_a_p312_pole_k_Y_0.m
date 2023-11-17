%%%%%%%%;
% Now use the templates with polar viewing-angles to form a new volume. ;
%%%%%%%%;
equ_theta = pi*0.25 ; %<-- some band around the equator. ;
cap_theta = pi/2-pi/14;%cap_theta = equ_theta; %<-- cap_theta = pi*0.25;
cap_sigma = cap_theta/3; %<-- so now 3 stds captures cap_theta. ;
tmp_index_pole_ = efind(min(viewing_polar_a_all_,pi-viewing_polar_a_all_)<=cap_theta);
tmp_index_equa_ = setdiff(0:n_S-1,tmp_index_pole_);
%%%%%%%%;
fname_fig_pre = sprintf('%s_jpg/test_heterogeneity_spurious_clathrin1l_equ_vs_pol_sigma_%.2d_FIGA_',dir_pm,floor(100*cap_sigma));
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_jpg_stripped = sprintf('%s_stripped.jpg',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_jpg));
%%%%;
figure(1+nf);nf=nf+1;clf;figmed;
view_ab_ = [00,05];
%%%%;
%subplot(1,3,1);
%plot_sphere_grid_0(struct('linecolor_a',0.35*[1,1,1],'linecolor_b',0.65*[1,1,1]));
%xlabel('k0'); ylabel(''); zlabel('k2');
%xlim(1.5*[-1,+1]); ylim(1.5*[-1,+1]); zlim(1.5*[-1,+1]); axis vis3d;
%set(gca,'XTick',[],'YTick',[],'ZTick',[]);
%view(view_ab_);
%%%%;
subplot(1,3,3);
plot_sphere_grid_0(struct('linecolor_a',0.35*[1,1,1],'linecolor_b',0.65*[1,1,1]));
imagesc_S_k_p_3d_belt_3( ...
 struct('type','parameter','k_p_r_max',0.975,'flag_fill',0,'linewidth_use',2,'c_line_use_',[1,0,1]) ...
,n_S - numel(tmp_index_equa_) ...
,viewing_azimu_b_all_(1+tmp_index_pole_) ...
,viewing_polar_a_all_(1+tmp_index_pole_) ...
);
xlabel('k0'); ylabel(''); zlabel('k2');
xlim(1.05*[-1,+1]); ylim(1.05*[-1,+1]); zlim(1.25*[0.0,+1]); axis equal; %axis vis3d;
set(gca,'XTick',[],'YTick',[],'ZTick',[]);
view([00,0]);%view(view_ab_);
title('zoom in on polar cap');
%%%%;
subplot(1,3,2);
plot_sphere_grid_0(struct('linecolor_a',0.35*[1,1,1],'linecolor_b',0.65*[1,1,1]));
imagesc_S_k_p_3d_belt_3( ...
 struct('type','parameter','k_p_r_max',0.975,'flag_fill',0,'linewidth_use',2,'c_line_use_',[1,0,1]) ...
,n_S - numel(tmp_index_equa_) ...
,viewing_azimu_b_all_(1+tmp_index_pole_) ...
,viewing_polar_a_all_(1+tmp_index_pole_) ...
);
xlabel('k0'); ylabel(''); zlabel('k2');
xlim(1.5*[-1,+1]); ylim(1.5*[-1,+1]); zlim(1.5*[-1,+1]); axis vis3d;
set(gca,'XTick',[],'YTick',[],'ZTick',[]);
view(view_ab_);
title('non-equatorial templates');
%%%%;
subplot(1,3,1);
plot_sphere_grid_0(struct('linecolor_a',0.35*[1,1,1],'linecolor_b',0.65*[1,1,1]));
imagesc_S_k_p_3d_belt_3( ...
 struct('type','parameter','k_p_r_max',0.975,'flag_fill',0,'linewidth_use',2,'c_line_use_',[1,0,1]) ...
,numel(tmp_index_equa_) ...
,viewing_azimu_b_all_(1+tmp_index_equa_) ...
,viewing_polar_a_all_(1+tmp_index_equa_) ...
);
xlabel('k0'); ylabel(''); zlabel('k2');
xlim(1.5*[-1,+1]); ylim(1.5*[-1,+1]); zlim(1.5*[-1,+1]); axis vis3d;
set(gca,'XTick',[],'YTick',[],'ZTick',[]);
view(view_ab_);
title('equatorial templates');
%%%%;
sgtitle('template [$\alpha$,$\beta$] in cyan, data in magenta','Interpreter','latex');
set(gcf,'Position',1+[0,0,1024*1.5,512]);
%%%%;
disp(sprintf(' %% writing %s',fname_fig_jpg));
print('-djpeg',fname_fig_jpg_stripped);
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
end;%if (~exist(fname_fig_jpg,'file'));
if ( exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s found, not creating',fname_fig_jpg));
end;%if ( exist(fname_fig_jpg,'file'));
close(gcf);
%%%%%%%%;

%%%%%%%%;
% set cap_theta: ;
%%%%%%%%;
cap_theta = pi/2-pi/14;%cap_theta = pi/4;

fname_mat = sprintf('%s_mat/test_heterogeneity_spurious_clathrin1l_cap_theta_%.3d_a_p312_pole_k_Y.mat',dir_pm,floor(100*cap_theta));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if flag_recalc | ~exist(fname_mat,'file');
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Now build a_p312_pole_x_u_xxx_ using a more restrictive polar cap ;
% associated with a distribution of viewing-angles concentrated around the equator and poles ;
%%%%%%%%;
S_p312_k_p_wkS__ = S_p312_k_p__;
template_k_r01_wkS__ = sqrt(template_k_c_0__.^2 + template_k_c_1__.^2);
template_k_r012_wkS__ = sqrt(template_k_c_0__.^2 + template_k_c_1__.^2 + template_k_c_2__.^2);
template_azimu_b_wkS__ = atan2(template_k_c_1__,template_k_c_0__);
template_polar_a_wkS__ = atan2(template_k_r01_wkS__,template_k_c_2__);
tmp_template_polar_a_wkS__ = min(1,min(template_polar_a_wkS__,pi-template_polar_a_wkS__)/max(1e-12,cap_theta));
g_cap_wkS__ = 1-sin((pi/2)*tmp_template_polar_a_wkS__).^4;
S_amp_k_p_wkS__ = S_p312_k_p_wkS__;
ind_0in_wkS__ = (abs(template_polar_a_wkS__-pi/2)<=pi/2-cap_theta); %<-- near equator. ;
ind_out_wkS__ = (abs(template_polar_a_wkS__-pi/2)> pi/2-cap_theta); %<-- near polar cap. ;
S_cap_avg = real(mean(S_p312_k_p_wkS__(find(ind_out_wkS__)))); %<-- average of (unused) polar cap values. ;
S_amp_k_p_wkS__ = S_amp_k_p_wkS__.*ind_0in_wkS__ ....
  + ( 0*S_cap_avg.*g_cap_wkS__ + S_amp_k_p_wkS__.*(1-g_cap_wkS__) ).*ind_out_wkS__ ...
  ; %<-- modify the poles. ;
%%%%%%%%;
% Now we can actually use all the templates in S_amp_k_p_wkS__, since we have set the polar data to be consistent. ;
%%%%%%%%;
flag_disp=0; if flag_disp; plot(viewing_polar_a_all_(tmp_nllr_sort_ij_rem_),'.'); end;
%%%%%%%%;
tmp_t = tic();
tmp_n_order = 5;
a_p312_pole_k_Y_yk_ = ...
cg_lsq_6( ...
 tmp_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,n_S ...
,S_amp_k_p_wkS__ ...
,[] ...
,[] ...
,viewing_polar_a_all_ ...
,viewing_azimu_b_all_ ...
,zeros(n_S,1) ...
,[] ...
,[] ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% a_p312_pole_k_Y_yk_ time %0.2fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
[ ... 
  a_p312_pole_x_u_xxx_ ...
] = ...
convert_spharm_to_x_c_4( ...
 k_eq_d ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_p312_pole_k_Y_yk_ ...
,half_diameter_x_c ...
,n_x_u_pack ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% a_p312_pole_x_u_xxx_ time %0.2fs',tmp_t)); end;
%%%%%%%%;
flag_disp=1;
if flag_disp;
prct_ = [98,98.5,99.0,99.5]; n_prct = numel(prct_);
for nprct=0:n_prct-1;
tmp_prct = prct_(1+nprct);
%%%%;
figure(1+nf);nf=nf+1;clf;figmed;fontsize_use = 12;
p_row = 1; p_col = 2; pcol=0;
subplot(p_row,p_col,1+pcol);pcol=pcol+1;
tmp_parameter = struct('type','parameter','percent_threshold_',[tmp_prct]);
tmp_parameter = isosurface_f_x_u_1(tmp_parameter,a_p312_x_u_reco_    ); tmp_parameter = rmfield(tmp_parameter,'vval_');
title(sprintf('a_p312_x_u_reco_xxx_: %.1f%%',tmp_prct),'Interpreter','none');
%title(sprintf('A: %.1f%%',tmp_prct),'Interpreter','none');
set(gca,'FontSize',fontsize_use);
subplot(p_row,p_col,1+pcol);pcol=pcol+1;
tmp_parameter = struct('type','parameter','percent_threshold_',[tmp_prct]);
tmp_parameter = isosurface_f_x_u_1(tmp_parameter,a_p312_pole_x_u_xxx_); tmp_parameter = rmfield(tmp_parameter,'vval_');
title(sprintf('a_p312_pole_x_u_xxx_: %.1f%%',tmp_prct),'Interpreter','none');
%title(sprintf('B: %.1f%%',tmp_prct),'Interpreter','none');
set(gca,'FontSize',fontsize_use);
set(gcf,'Position',1+[0,0,1024*2,512]); drawnow();
%%%%;
fname_fig_pre = sprintf('%s_jpg/test_heterogeneity_spurious_clathrin1l_p312_pole_p%.4d_FIGD',dir_pm,round(100*tmp_prct));
fname_fig_pre_stripped = sprintf('%s_stripped',fname_fig_pre);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_jpg_stripped = sprintf('%s.jpg',fname_fig_pre_stripped);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_jpg));
print('-djpeg',fname_fig_jpg_stripped);
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
end;%if (flag_replot | ~exist(fname_fig_jpg,'file'));
close(gcf);
%%%%;
end;%for nprct=0:n_prct-1;
end;%if flag_disp;
%%%%%%%%
save(fname_mat ...
,'cap_theta' ...
,'cap_sigma' ...
,'tmp_index_pole_' ...
,'tmp_index_equa_' ...
,'a_p312_pole_k_Y_yk_' ...
);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if ~exist(fname_mat,'file');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if  exist(fname_mat,'file');
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if  exist(fname_mat,'file');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% Now generate all templates from restrictive polar cap volume. ;
%%%%%%%%;
a_p312_pole_k_Y_yk__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
a_p312_pole_k_Y_yk__(1:n_lm_(1+nk_p_r),1+nk_p_r) = a_p312_pole_k_Y_yk_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
template_k_eq_d = -1;
viewing_k_eq_d = 1.0; %<-- subsample just a few templates for visualization. ;
tmp_t = tic();
[ ...
 S_p312_pole_k_p_wkS__ ...
,~ ...
,n_S ...
,template_viewing_azimu_b_all_ ...
,template_viewing_polar_a_all_ ...
] = ...
pm_template_2( ...
 verbose ...
,l_max_max ...
,n_k_p_r ...
,reshape(a_p312_pole_k_Y_yk__,[n_lm_max,n_k_p_r]) ...
,viewing_k_eq_d/k_p_r_max ...
,template_k_eq_d ...
,n_w_max ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% pm_template_2: %0.6fs',tmp_t)); end;
S_p312_pole_k_p_wkS__ = reshape(S_p312_pole_k_p_wkS__,[n_w_max*n_k_p_r,n_S]);
S_p312_pole_k_q_wkS__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
S_p312_pole_k_q_wkS__(:,1+nS) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,S_p312_pole_k_p_wkS__(:,1+nS));
end;%for nS=0:n_S-1;
%%%%%%%%;
