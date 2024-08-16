%%%%%%%%;
% Testing RNA molecule from Pilar (20240813). ;
%%%%%%%%;

platform = 'rusty';%platform = 'access1';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data/rangan'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home/rangan'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home/rangan'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home/rangan'; end;
if (strcmp(platform,'ceph')); setup_ceph; string_root = 'mnt/home/rangan/ceph'; end;

flag_verbose = 1;
flag_recalc = 0;
flag_disp = 0; nf=0;
flag_replot = 0;
flag_center = 0;
tolerance_master = 1e-2;
nf=0;

dir_pm = sprintf('/%s/dir_cryoem/dir_RNA0/dir_pm',string_root);
if (~exist(sprintf('%s_mat',dir_pm),'dir')); disp(sprintf(' %% mkdir %s_mat',dir_pm)); mkdir(sprintf('%s_mat',dir_pm)); end;
if (~exist(sprintf('%s_jpg',dir_pm),'dir')); disp(sprintf(' %% mkdir %s_jpg',dir_pm)); mkdir(sprintf('%s_jpg',dir_pm)); end;
dir_data = sprintf('/%s/dir_cryoem/dir_RNA0',string_root);
Pixel_Spacing = 2.676; %<-- in angstroms. ;
n_x_u = 128; %<-- box diameter in pixels. ;
sigma_angstrom = 8.0; %<-- low resolution. ;
str_name_volume_ = {'WT_Folded','MUT_Open'};
fname_nopath_volume_ = {,'models/folded.pdb','models/open.pdb'};
n_volume = numel(fname_nopath_volume_);
fname_volume_ = cell(n_volume,1);
for nvolume=0:n_volume-1;
fname_volume_{1+nvolume} = sprintf('%s/%s',dir_data,fname_nopath_volume_{1+nvolume});
end;%for nvolume=0:n_volume-1;
fname_nopath_star_ = {'J738_passthrough_particles.star','J737_passthrough_particles.star'};
n_star = numel(fname_nopath_star_);
fname_star_ = cell(n_star,1);
for nstar=0:n_star-1;
fname_star_{1+nstar} = sprintf('%s/%s',dir_data,fname_nopath_star_{1+nstar});
end;%for nstar=0:n_star-1;

%%%%%%%%;
% First create volumes. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/a_x_u_base_xxxv__.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
half_diameter_x_c = 1.0d0;
diameter_x_c = 2.0d0*half_diameter_x_c;
x_p_r_max = 1.0;
n_x_u_pack = 64;
n_pack = n_x_u/n_x_u_pack;
pack_row_ij_ = zeros(n_x_u_pack,1);
pack_col_ij_ = zeros(n_x_u_pack,1);
pack_val_ij_ = zeros(n_x_u_pack,1);
na=0;
for nx_u=0:n_x_u-1;
pack_row_ij_(1+na) = 1+nx_u;
pack_col_ij_(1+na) = 1+floor(nx_u/n_pack);
pack_val_ij_(1+na) = 1/n_pack;
na=na+1;
end;%for nx_u=0:n_x_u-1;
x_u_pack_ = sparse(pack_row_ij_,pack_col_ij_,pack_val_ij_,n_x_u,n_x_u_pack);
x_u_0_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack);
x_u_1_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack);
x_u_2_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack);
[x_u_0___,x_u_1___,x_u_2___] = ndgrid(x_u_0_,x_u_1_,x_u_2_); n_xxx_u = n_x_u_pack^3;
xxx_u_weight_ = (2*x_p_r_max/n_x_u_pack)^3;
dx = diameter_x_c/max(1,n_x_u_pack);
%%%%%%%%;
a_x_u_base_xxxv__ = zeros(n_xxx_u,n_volume);
for nvolume=0:n_volume-1;
fname_pdb = fname_volume_{1+nvolume};
if (flag_verbose>0); disp(sprintf(' %% nvolume %d/%d: %s',nvolume,n_volume,fname_pdb)); end;
parameter = struct('type','parameter');
parameter.sigma_angstrom = sigma_angstrom;
[ ...
 parameter ...
,a_x_u_load_ ...
] = ...
a_x_u___from_pdb_0( ...
 parameter ...
,fname_pdb ...
,n_x_u ...
,Pixel_Spacing ...
);
%%%%%%%%;
a_x_u_pack_ = reshape(a_x_u_load_,[n_x_u*n_x_u,n_x_u])*x_u_pack_;
a_x_u_pack_ = reshape(permute(reshape(a_x_u_pack_,[n_x_u,n_x_u,n_x_u_pack]),[3,1,2]),[n_x_u*n_x_u_pack,n_x_u])*x_u_pack_;
a_x_u_pack_ = reshape(permute(reshape(a_x_u_pack_,[n_x_u_pack,n_x_u,n_x_u_pack]),[3,1,2]),[n_x_u_pack*n_x_u_pack,n_x_u])*x_u_pack_;
a_x_u_pack_ = permute(reshape(a_x_u_pack_,[n_x_u_pack,n_x_u_pack,n_x_u_pack]),[3,1,2]);
a_x_u_pack_ = a_x_u_pack_ - mean(a_x_u_pack_,'all');
a_x_u_pack_ = a_x_u_pack_/max(1e-12,sum(a_x_u_pack_,'all')*dx^3);
clear a_x_u_load_;
%%%%%%%%;
% Calculate moments. ;
%%%%%%%%;
a_rho_x_u_pack_ = a_x_u_pack_ + min(a_x_u_pack_,[],'all');
a_rho_x_c_0_avg = sum(x_u_0___.^1.*a_rho_x_u_pack_/sum(a_rho_x_u_pack_,'all'),'all');
a_rho_x_c_1_avg = sum(x_u_1___.^1.*a_rho_x_u_pack_/sum(a_rho_x_u_pack_,'all'),'all');
a_rho_x_c_2_avg = sum(x_u_2___.^1.*a_rho_x_u_pack_/sum(a_rho_x_u_pack_,'all'),'all');
a_rho_x_c_0_std = sum((x_u_0___ - a_rho_x_c_0_avg).^2.*a_rho_x_u_pack_/sum(a_rho_x_u_pack_,'all'),'all');
a_rho_x_c_1_std = sum((x_u_1___ - a_rho_x_c_1_avg).^2.*a_rho_x_u_pack_/sum(a_rho_x_u_pack_,'all'),'all');
a_rho_x_c_2_std = sum((x_u_2___ - a_rho_x_c_2_avg).^2.*a_rho_x_u_pack_/sum(a_rho_x_u_pack_,'all'),'all');
a_rho_x_c_avg_ = [a_rho_x_c_0_avg ; a_rho_x_c_1_avg ; a_rho_x_c_2_avg];
a_rho_x_c_std_ = [a_rho_x_c_0_std ; a_rho_x_c_1_std ; a_rho_x_c_2_std];
disp(sprintf(' %% a_rho_x_c_std_ vs a_rho_x_c_avg_: %0.2f',fnorm(a_rho_x_c_std_)/fnorm(a_rho_x_c_avg_)));
if (max(abs(a_rho_x_c_avg_))> diameter_x_c/n_x_u_pack);
disp(sprintf(' %% Warning, molecule may not be well centered. Consider recentering.'));
end;%if (max(abs(a_rho_x_c_avg_))> diameter_x_c/n_x_u_pack);
a_x_u_base_xxxv__(:,1+nvolume) = a_x_u_pack_(:);
clear a_x_u_pack_ a_rho_x_u_pack_ a_rho_x_c_0_avg a_rho_x_c_1_avg a_rho_x_c_2_avg a_rho_x_c_0_std a_rho_x_c_1_std a_rho_x_c_2_std a_rho_x_c_avg_ a_rho_x_c_std_ ;
%%%%%%%%;
end;%for nvolume=0:n_volume-1;
%%%%%%%%;
if flag_disp;
fname_fig_pre = sprintf('%s_jpg/a_x_u_base_xxxv__',dir_pm);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
disp(sprintf(' %% %s not found, creating',fname_fig_pre));
figure(1+nf);nf=nf+1;clf;figbig;
for nvolume=0:n_volume-1;
subplot(1,n_volume,1+nvolume);
isosurface_f_x_u_1([],a_x_u_base_xxxv__(:,1+nvolume));
title(str_name_volume_{1+nvolume},'Interpreter','none');
end;%for nvolume=0:n_volume-1;
sgtitle(fname_fig_pre,'Interpreter','none');
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
end;%if (~exist(fname_fig_jpg,'file'));
close(gcf);
end;%if flag_disp;
%%%%%%%%;
save(fname_mat ...
     ,'half_diameter_x_c' ...
     ,'diameter_x_c' ...
     ,'x_p_r_max' ...
     ,'n_x_u_pack' ...
     ,'n_pack' ...
     ,'pack_row_ij_' ...
     ,'pack_col_ij_' ...
     ,'pack_val_ij_' ...
     ,'x_u_pack_' ...
     ,'x_u_0_' ...
     ,'x_u_1_' ...
     ,'x_u_2_' ...
     ,'x_u_0___' ...
     ,'x_u_1___' ...
     ,'x_u_2___' ...
     ,'n_x_u' ...
     ,'n_xxx_u' ...
     ,'xxx_u_weight_' ...
     ,'a_x_u_base_xxxv__' ...
     );
end;%if (flag_recalc | ~exist(fname_mat,'file'));
%%%%%%%%;
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));

%%%%%%%%;
% Load class averages. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/O_k_p_wkOv___.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
for nvolume=0:n_volume-1;
fname_mrc = sprintf('%s/J834_templates_slected.mrc',dir_data); %<-- unfortunately only have one set of class-averages. ;
O_x_c_xxO___ = cast(ReadMRC_0(fname_mrc),'double');
n_O_sub = size(O_x_c_xxO___,3);
O_k_p_wkO__ = zeros(n_w_sum,n_O_sub);
for nO_sub=0:n_O_sub-1;
O_k_p_wkO__(:,1+nO_sub) = ...
interp_x_c_to_k_p_xxnufft( ...
 size(O_x_c_xxO___,1) ...
,diameter_x_c ...
,size(O_x_c_xxO___,2) ...
,diameter_x_c ...
,O_x_c_xxO___(:,:,1+nO_sub) ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
);
end;%for nO_sub=0:n_O_sub-1;
O_k_p_wkOv___(:,:,1+nvolume) = O_k_p_wkO__;
clear O_x_c_xxO___ O_k_p_wkO__ ;
end;%for nvolume=0:n_volume-1;
save(fname_mat ...
     ,'n_O_sub','O_k_p_wkOv___' ...
     );
end;%if (flag_recalc | ~exist(fname_mat,'file'));
%%%%%%%%;
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
if flag_disp;
fname_fig_pre = sprintf('%s_jpg/O_k_p_wkOv___',dir_pm);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 3; p_col = 6; np=0;
Olim_ = std(abs(O_k_p_wkOv___(:)))*2.5*[-1,+1];
for np=0:p_row*p_col-1;
subplot(p_row,p_col,1+np);
nOv = max(0,min(n_O_sub*n_volume-1,floor(n_O_sub*n_volume*np/(p_row*p_col))));
nO_sub = mod(nOv,n_O_sub);
nvolume = floor(nOv/n_O_sub);
O_k_p_wk_ = O_k_p_wkOv___(:,1+nO_sub,1+nvolume);
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(O_k_p_wk_),Olim_,colormap_80s);
axis image; axisnotick; title(sprintf('nO_sub %d nvolume %d',nO_sub,nvolume));
end;%for np=0:p_row*p_col-1;
sgtitle(fname_fig_pre,'Interpreter','none');
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_pre));
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
end;%if (~exist(fname_fig_jpg,'file'));
close(gcf);
end;%if flag_disp;
%%%%;
if flag_disp;
fname_fig_pre = sprintf('%s_jpg/O_x_c_xxOv___',dir_pm);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 3; p_col = 6; np=0;
for np=0:p_row*p_col-1;
subplot(p_row,p_col,1+np);
nOv = max(0,min(n_O_sub*n_volume-1,floor(n_O_sub*n_volume*np/(p_row*p_col))));
nO_sub = mod(nOv,n_O_sub);
nvolume = floor(nOv/n_O_sub);
O_k_p_wk_ = O_k_p_wkOv___(:,1+nO_sub,1+nvolume);
O_x_c_xx_ = ...
interp_k_p_to_x_c_xxnufft( ...
 n_x_u_pack ...
,diameter_x_c ...
,n_x_u_pack ...
,diameter_x_c ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,O_k_p_wk_ ...
);
Olim_ = prctile(real(O_x_c_xx_),[  0,100],'all'); Olim_  = mean(Olim_) + 1.25*0.5*diff(Olim_)*[-1,+1];
imagesc_c(n_x_u_pack,x_u_0_,n_x_u_pack,x_u_1_,real(O_x_c_xx_),Olim_,colormap_beach);
axis image; axisnotick; title(sprintf('nO_sub %d nvolume %d',nO_sub,nvolume));
end;%for np=0:p_row*p_col-1;
sgtitle(fname_fig_pre,'Interpreter','none');
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_pre));
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
end;%if (~exist(fname_fig_jpg,'file'));
close(gcf);
end;%if flag_disp;
%%%%%%%%;
n_O = n_O_sub*n_volume;
O_k_p_wkO__ = reshape(O_k_p_wkOv___,[n_w_sum,n_O]);
index_nCTF_from_nO_ = zeros(n_O,1);

if flag_disp>1;
%%%%%%%%;
% step through the class averages. ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 3; p_col = 6; np=0;
for nO=0:n_O-1;
subplot(p_row,p_col,1+np); cla;
O_k_p_wk_ = O_k_p_wkO__(:,1+nO);
O_x_c_xx_ = ...
interp_k_p_to_x_c_xxnufft( ...
 n_x_u_pack ...
,diameter_x_c ...
,n_x_u_pack ...
,diameter_x_c ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,O_k_p_wk_ ...
);
Olim_ = prctile(real(O_x_c_xx_),[  0,100],'all'); Olim_  = mean(Olim_) + 1.25*0.5*diff(Olim_)*[-1,+1];
imagesc_c(n_x_u_pack,x_u_0_,n_x_u_pack,x_u_1_,real(O_x_c_xx_),Olim_,colormap_beach);
axis image; axisnotick; title(sprintf('nO %d',nO));
drawnow();
np=np+1; if np>=p_row*p_col; np=0; end;
end;%for nO=0:n_O-1;
end;%if flag_disp>1;

%%%%%%%%;
if flag_disp;
fname_fig_pre = sprintf('%s_jpg/a_x_u_base_xxxv____vs_class_average',dir_pm);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
disp(sprintf(' %% %s not found, creating',fname_fig_pre));
fname_mrc = sprintf('%s/J834_templates_slected.mrc',dir_data);
tmp_S_x_c_xxS___ = cast(ReadMRC_0(fname_mrc),'double');
tmp_n_S = size(tmp_S_x_c_xxS___,3);
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 4; p_col = 6; np=0;
for nvolume=0:n_volume-1;
a_x_u_base_xxx___ = reshape(a_x_u_base_xxxv__(:,1+nvolume),[n_x_u_pack,n_x_u_pack,n_x_u_pack]);
S_x_c_0_xx__ = squeeze(mean(a_x_u_base_xxx___,1+0));
S_x_c_1_xx__ = squeeze(mean(a_x_u_base_xxx___,1+1));
S_x_c_2_xx__ = squeeze(mean(a_x_u_base_xxx___,1+2));
subplot(p_row,p_col,1+np);np=np+1;
imagesc_c(n_x_u_pack,x_u_0_,n_x_u_pack,x_u_1_,real(S_x_c_0_xx__),[],colormap_beach());
axis image; axisnotick; title(sprintf('nvolume %d',nvolume),'Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
imagesc_c(n_x_u_pack,x_u_0_,n_x_u_pack,x_u_1_,real(S_x_c_1_xx__),[],colormap_beach());
axis image; axisnotick; title(sprintf('nvolume %d',nvolume),'Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
imagesc_c(n_x_u_pack,x_u_0_,n_x_u_pack,x_u_1_,real(S_x_c_2_xx__),[],colormap_beach());
axis image; axisnotick; title(sprintf('nvolume %d',nvolume),'Interpreter','none');
clear a_x_u_base_xxx___ S_x_c_0_xx__ S_x_c_1_xx__ S_x_c_2_xx__ ;
end;%for nvolume=0:n_volume-1;
tmp_nS=0;
while np<p_row*p_col;
if tmp_nS<tmp_n_S;
subplot(p_row,p_col,1+np);np=np+1;
imagesc_c([],[],[],[],real(tmp_S_x_c_xxS___(:,:,1+tmp_nS)),[],colormap_beach());
axis image; axisnotick; title(sprintf('tmp_nS %d',tmp_nS),'Interpreter','none');
tmp_nS=tmp_nS+1;
end;%if tmp_nS<tmp_n_S;
end;%while np<p_row*p_col;
sgtitle(fname_fig_pre,'Interpreter','none');
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
end;%if (~exist(fname_fig_jpg,'file'));
close(gcf);
end;%if flag_disp;
%%%%%%%%;

%%%%%%%%;
% Now convert to a_k_p_ ;
%%%%%%%%;
verbose=0;
fname_mat = sprintf('%s_mat/a_k_p_quad_allv__.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
tmp_t = tic;
k_p_r_max = 48/(2*pi); k_eq_d = 1.0/(2*pi); str_T_vs_L = 'L';
flag_unif_vs_adap = 0; flag_tensor_vs_adap = 0; %<-- This is set to match test_ssnll_from_a_k_Y_12 ;
[ ...
 n_k_all ...
,n_k_all_csum_ ...
,k_p_r_all_ ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,weight_3d_k_all_ ...
,weight_shell_k_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,k_c_0_all_ ...
,k_c_1_all_ ...
,k_c_2_all_ ...
,J_node_ ...
,J_weight_ ...
,J_chebfun_ ...
,J_polyval_ ...
] = ...
sample_sphere_7( ...
 verbose ...
,k_p_r_max ...
,k_eq_d ...
,str_T_vs_L ...
,flag_unif_vs_adap ...
,flag_tensor_vs_adap ...
) ; %<-- sum(weight_3d_k_p_r_)*(4*pi) = (4/3)*pi*k_p_r_max^3 --> sum(weight_3d_k_p_r_) = (1/3)*k_p_r_max^3 ;
tmp_t = toc(tmp_t); disp(sprintf(' %% sample_sphere_7: time %0.2fs',tmp_t));
%%%%%%%%;
a_k_p_quad_allv__ = zeros(n_k_all,n_volume);
a_x_u_reco_v_xxxv__ = zeros(n_xxx_u,n_volume);
for nvolume=0:n_volume-1;
a_x_u_base_ = a_x_u_base_xxxv__(:,1+nvolume);
%%%%;
eta = pi/x_p_r_max; tmp_t = tic;
a_k_p_quad_ = xxnufft3d3(n_xxx_u,x_u_0___(:)*eta,x_u_1___(:)*eta,x_u_2___(:)*eta,a_x_u_base_.*xxx_u_weight_(:),-1,1e-12,n_k_all,2*pi*k_c_0_all_/eta,2*pi*k_c_1_all_/eta,2*pi*k_c_2_all_/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_k_p_quad_ time %0.2fs',tmp_t));
%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
a_x_u_reco_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,a_k_p_quad_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_x_u_reco_ time %0.2fs',tmp_t));
disp(sprintf(' %% xxnufft3d3: a_x_u_reco error: %0.16f',fnorm(a_x_u_base_-a_x_u_reco_)/fnorm(a_x_u_base_)));
%%%%;
a_k_p_quad_allv__(:,1+nvolume) = a_k_p_quad_;
a_x_u_reco_xxxv__(:,1+nvolume) = a_x_u_reco_(:);
clear a_x_u_base_ a_k_p_quad_ a_x_u_reco_;
end;%for nvolume=0:n_volume-1;
%%%%%%%%;
save(fname_mat ...
     ,'k_p_r_max','k_eq_d' ...
     ,'n_k_all','n_k_all_csum_' ...
     ,'k_p_r_all_','k_p_azimu_b_all_','k_p_polar_a_all_' ...
     ,'weight_3d_k_all_','weight_shell_k_' ...
     ,'n_k_p_r','k_p_r_' ...
     ,'weight_3d_k_p_r_' ...
     ,'k_c_0_all_','k_c_1_all_','k_c_2_all_' ...
     ,'J_node_','J_weight_','J_chebfun_','J_polyval_' ...
     ,'a_k_p_quad_allv__' ...
     ,'a_x_u_reco_xxxv__' ...
     );
end;%if (~exist(fname_mat,'file'));
%%%%%%%%;
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
if flag_disp;
fname_fig_pre = sprintf('%s_jpg/a_k_p_quad_allv__',dir_pm);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
figure(1+nf);nf=nf+1;clf;figbig;
for nvolume=0:n_volume-1;
subplot(1,n_volume,1+nvolume);
plot(k_p_r_all_,log10(abs(a_k_p_quad_allv__(:,1+nvolume))),'.');
xlabel('k'); ylabel('log10(|a(k)|)');
title(str_name_volume_{1+nvolume},'Interpreter','none');
end;%for nvolume=0:n_volume-1;
sgtitle(fname_fig_pre,'Interpreter','none');
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_pre));
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
end;%if (~exist(fname_fig_jpg,'file'));
close(gcf);
end;%if flag_disp;
%%%%%%%%;

if ~exist('Ylm_uklma___','var'); Ylm_uklma___ = []; end;
if ~exist('k_p_azimu_b_sub_uka__','var'); k_p_azimu_b_sub_uka__ = []; end;
if ~exist('k_p_polar_a_sub_uka__','var'); k_p_polar_a_sub_uka__ = []; end;
if ~exist('l_max_uk_','var'); l_max_uk_ = []; end;
if ~exist('index_nu_n_k_per_shell_from_nk_p_r_','var'); index_nu_n_k_per_shell_from_nk_p_r_ = []; end;
if ~exist('index_k_per_shell_uka__','var'); index_k_per_shell_uka__ = []; end;

%%%%%%%%;
% Now convert to a_k_Y_ ; 
%%%%%%%%;
verbose=0;
fname_mat = sprintf('%s_mat/a_k_Y_quad_ykv__.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
l_max_upb = round(2*pi*k_p_r_max); %<-- typically sufficient for 2-3 digits of precision. ;
l_max_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
l_max_(1+nk_p_r) = max(0,min(l_max_upb,1+ceil(2*pi*k_p_r_(1+nk_p_r))));
end;%for nk_p_r=0:n_k_p_r-1;
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
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
Y_l_val_(1+tmp_index_) = tmp_l_val_;
Y_m_val_(1+tmp_index_) = tmp_m_val_;
Y_k_val_(1+tmp_index_) = k_p_r_(1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
weight_Y_ = zeros(n_lm_sum,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
weight_Y_(1+tmp_index_) = weight_3d_k_p_r_(1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
a_k_Y_quad_ykv__ = zeros(n_lm_sum,n_volume);
a_k_Y_quad_ykv___ = zeros(n_lm_max,n_k_p_r,n_volume);
a_k_p_reco_allv__ = zeros(n_k_all,n_volume);
for nvolume=0:n_volume-1;
a_k_p_quad_ = a_k_p_quad_allv__(:,1+nvolume);
if (flag_verbose>0); disp(sprintf(' %% nvolume %d/%d',nvolume,n_volume)); end;
tmp_t = tic;
[ ...
 a_k_Y_quad_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
] = ...
convert_k_p_to_spharm_4( ...
 0*flag_verbose ...
,n_k_all ...
,n_k_all_csum_ ...
,k_p_r_all_ ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,weight_3d_k_all_ ...
,weight_shell_k_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_p_quad_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_quad_ time %0.2fs',tmp_t));
disp(sprintf(' %% Here we should ensure that a_k_Y_ on the outer shells has decayed to the desired precision.'));
disp(sprintf(' %% Moreover, we should ensure that a_k_Y_(l,m) has decayed for large l,m at each shell.'));
tmp_t = tic;
[...
 a_k_p_reco_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
] = ...
convert_spharm_to_k_p_4( ...
 0*flag_verbose ...
,n_k_all ...
,n_k_all_csum_ ...
,k_p_r_all_ ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,weight_3d_k_all_ ...
,weight_shell_k_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_quad_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_quad_ --> a_k_p_reco_ time %0.2fs',tmp_t));
disp(sprintf(' %% xxnufft3d3: a_k_p_reco error: %0.16f',fnorm(a_k_p_quad_-a_k_p_reco_)/fnorm(a_k_p_quad_))); %<-- this should be 2-3 digits. ;
a_k_Y_quad__ = local_yk__from_yk_(n_k_p_r,l_max_,a_k_Y_quad_);
a_k_Y_quad_ykv__(:,1+nvolume) = a_k_Y_quad_;
a_k_p_reco_allv__(:,1+nvolume) = a_k_p_reco_;
a_k_Y_quad_ykv___(:,:,1+nvolume) = a_k_Y_quad__;
clear a_k_p_quad_ a_k_Y_quad_ a_k_p_reco_ a_k_Y_quad__;
end;%for nvolume=0:n_volume-1;
%%%%%%%%;
save(fname_mat ...
     ,'l_max_' ...
     ,'n_lm_','n_lm_max','n_lm_sum','n_lm_csum_','l_max_max','m_max_','n_m_max' ...
     ,'Y_l_val_','Y_m_val_','Y_k_val_','weight_Y_' ...
     ,'a_k_Y_quad_ykv__' ...
     ,'a_k_p_reco_allv__' ...
     ,'a_k_Y_quad_ykv___' ...
     );
end;%if (~exist(fname_mat,'file'));
%%%%%%%%;
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
if flag_disp;
fname_fig_pre = sprintf('%s_jpg/a_k_Y_quad_A',dir_pm);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
figure(1+nf);nf=nf+1;clf;figbig;
p_row = n_volume; p_col = 3; np=0;
for nvolume=0:n_volume-1;
a_k_Y_quad_ = a_k_Y_quad_ykv__(:,1+nvolume);
subplot(p_row,p_col,1+np);np=np+1;
plot(Y_l_val_,log10(abs(a_k_Y_quad_)),'.'); xlabel('Y_l_val_','Interpreter','none'); ylabel('log10(abs(a))');
subplot(p_row,p_col,1+np);np=np+1;
plot(Y_m_val_,log10(abs(a_k_Y_quad_)),'.'); xlabel('Y_m_val_','Interpreter','none'); ylabel('log10(abs(a))');
subplot(p_row,p_col,1+np);np=np+1;
plot(Y_k_val_,log10(abs(a_k_Y_quad_)),'.'); xlabel('Y_k_val_','Interpreter','none'); ylabel('log10(abs(a))');
title(str_name_volume_{1+nvolume},'Interpreter','none');
clear a_k_Y_quad_;
end;%for nvolume=0:n_volume-1;
sgtitle(fname_fig_pre,'Interpreter','none');
if (flag_replot | ~exist(fname_fig_pre,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_pre));
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
end;%if (~exist(fname_fig_pre,'file'));
close(gcf);
end;%if flag_disp;
%%%%%%%%;

%%%%%%%%;
% generate templates S_k_p_ on k_p_ grid with uniform n_w_. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/S_k_p_wkSv___.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
n_w_max = 2*(l_max_max+1);
n_w_ = n_w_max*ones(n_k_p_r,1); n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_]);
template_k_eq_d = -1;
viewing_k_eq_d = 1.0; %<-- subsample just a few templates for visualization. ;
S_k_p_wkSv___ = [];
for nvolume=0:n_volume-1;
a_k_Y_quad_yk__ = a_k_Y_quad_ykv___(:,:,1+nvolume);
template_k_eq_d = -1;
tmp_t=tic();
[ ...
 S_k_p_wkS__ ...
,n_w ...
,n_viewing_all ...
,viewing_azimu_b_all_ ...
,viewing_polar_a_all_ ...
,viewing_weight_all_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
] = ...
pm_template_2( ...
 verbose ...
,l_max_max ...
,n_k_p_r ...
,a_k_Y_quad_yk__ ...
,viewing_k_eq_d/max(1e-12,k_p_r_max) ...
,template_k_eq_d ...
,n_w_max ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% pm_template_2: time %0.6fs',tmp_t)); end;
n_S = n_viewing_all;
S_k_p_wkS__ = reshape(S_k_p_wkS__,[n_w_max*n_k_p_r,n_S]);
if  isempty(S_k_p_wkSv___); S_k_p_wkSv___ = zeros(n_w_sum,n_S,n_volume); end;
S_k_p_wkSv___(:,:,1+nvolume) = S_k_p_wkS__;
clear a_k_Y_quad_yk__ S_k_p_wkS__;
end;%for nvolume=0:n_volume-1;
[ ...
 n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_wk_ ...
,k_p_r_wk_ ...
,k_p_w_wk_ ...
,k_c_0_wk_ ...
,k_c_1_wk_ ...
] = ...
get_weight_2d_2( ...
 verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,template_k_eq_d ...
,n_w_ ...
);
save(fname_mat ...
     ,'n_w_max','template_k_eq_d','viewing_k_eq_d' ...
     ,'S_k_p_wkSv___' ...
     ,'n_w_' ...
     ,'weight_2d_k_p_r_' ...
     ,'weight_2d_wk_' ...
     ,'n_viewing_all' ...
     ,'viewing_azimu_b_all_' ...
     ,'viewing_polar_a_all_' ...
     ,'viewing_weight_all_' ...
     ,'n_viewing_polar_a' ...
     ,'viewing_polar_a_' ...
     ,'n_viewing_azimu_b_' ...
     ,'n_S','n_w_max','n_w_sum','n_w_csum_' ...
     );
end;%if (flag_recalc | ~exist(fname_mat,'file'));
%%%%%%%%;
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
if flag_disp;
fname_fig_pre = sprintf('%s_jpg/S_k_p_wkSv___',dir_pm);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 3; p_col = 6; np=0;
Slim_ = std(abs(S_k_p_wkSv___(:)))*2.5*[-1,+1];
for np=0:p_row*p_col-1;
subplot(p_row,p_col,1+np);
nSv = max(0,min(n_S*n_volume-1,floor(n_S*n_volume*np/(p_row*p_col))));
nS = mod(nSv,n_S);
nvolume = floor(nSv/n_S);
S_k_p_wk_ = S_k_p_wkSv___(:,1+nS,1+nvolume);
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(S_k_p_wk_),Slim_,colormap_80s);
axis image; axisnotick; title(sprintf('nS %d nvolume %d',nS,nvolume));
end;%for np=0:p_row*p_col-1;
sgtitle(fname_fig_pre,'Interpreter','none');
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_pre));
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
end;%if (~exist(fname_fig_jpg,'file'));
close(gcf);
end;%if flag_disp;
%%%%;
if flag_disp;
fname_fig_pre = sprintf('%s_jpg/S_x_c_xxSv___',dir_pm);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 3; p_col = 6; np=0;
for np=0:p_row*p_col-1;
subplot(p_row,p_col,1+np);
nSv = max(0,min(n_S*n_volume-1,floor(n_S*n_volume*np/(p_row*p_col))));
nS = mod(nSv,n_S);
nvolume = floor(nSv/n_S);
S_k_p_wk_ = S_k_p_wkSv___(:,1+nS,1+nvolume);
S_x_c_xx_ = ...
interp_k_p_to_x_c_xxnufft( ...
 n_x_u_pack ...
,diameter_x_c ...
,n_x_u_pack ...
,diameter_x_c ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,S_k_p_wk_ ...
);
Slim_ = prctile(real(S_x_c_xx_),[  0,100],'all'); Slim_  = mean(Slim_) + 1.25*0.5*diff(Slim_)*[-1,+1];
imagesc_c(n_x_u_pack,x_u_0_,n_x_u_pack,x_u_1_,real(S_x_c_xx_),Slim_,colormap_beach);
axis image; axisnotick; title(sprintf('nS %d nvolume %d',nS,nvolume));
end;%for np=0:p_row*p_col-1;
sgtitle(fname_fig_pre,'Interpreter','none');
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_pre));
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
end;%if (~exist(fname_fig_jpg,'file'));
close(gcf);
end;%if flag_disp;
%%%%%%%%;

%%%%%%%%;
% Now load images and CTF parameters from the star-file. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/M_k_p_wkMs___.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
n_M_sub = 1024;
M_k_p_wkMs___=zeros(n_w_sum,n_M_sub,n_star);
index_nCTF_from_nM_Ms__ = zeros(n_M_sub,n_star);
index_nM_from_nCTF_s_ = cell(n_star,1);
Voltage_CTF_s_ = cell(n_star,1);
DefocusU_CTF_s_ = cell(n_star,1);
DefocusV_CTF_s_ = cell(n_star,1);
DefocusAngle_CTF_s_ = cell(n_star,1);
SphericalAberration_CTF_s_ = cell(n_star,1);
AmplitudeContrast_CTF_s_ = cell(n_star,1);
for nstar=0:n_star-1;
[ ...
 index_nCTF_from_nM_Ms__(:,1+nstar) ...
,index_nM_from_nCTF_s_{1+nstar} ...
,Voltage_CTF_s_{1+nstar} ...
,DefocusU_CTF_s_{1+nstar} ...
,DefocusV_CTF_s_{1+nstar} ...
,DefocusAngle_CTF_s_{1+nstar} ...
,SphericalAberration_CTF_s_{1+nstar} ...
,AmplitudeContrast_CTF_s_{1+nstar} ...
] = ...
rlnCTFParameters_from_star_2( ...
 dir_data ...
,fname_nopath_star_{1+nstar} ...
,n_M_sub ...
);
if (fnorm(Voltage_CTF_s_{1+nstar})< 1e-3); disp(sprintf(' %% Warning, Voltage not set, setting Voltage to 300kV')); Voltage_CTF_s_{1+nstar} = 300*ones(n_M_sub,1); end;
%%%%;
if strcmp(fname_nopath_star_{1+nstar},'J738_passthrough_particles.star'); str_infix = 'WT_Folded'; end;
if strcmp(fname_nopath_star_{1+nstar},'J737_passthrough_particles.star'); str_infix = 'MUT_Open'; end;
M_x_u_xxM___ = zeros(n_x_u,n_x_u,n_M_sub);
nM_sub=0; nbatch=0; flag_continue=1;
while flag_continue;
fname_mrc = sprintf('%s/%s/batch_%.6d_downsample.mrc',dir_data,str_infix,nbatch);
if ~exist(fname_mrc,'file'); disp(sprintf(' %% Warning, %s not found',fname_mrc)); flag_continue=0; end;
if  exist(fname_mrc,'file');
if (flag_verbose>0); disp(sprintf(' %% %s found',fname_mrc)); end;
tmp_M_x_u_xxM___ = cast(ReadMRC_0(fname_mrc),'double');
tmp_n_M_sub = size(tmp_M_x_u_xxM___,3);
tmp_n_M_sub = min(n_M_sub-nM_sub,tmp_n_M_sub);
if (flag_verbose>0); disp(sprintf(' %% %s found: using %d',fname_mrc,tmp_n_M_sub)); end;
M_x_u_xxM___(:,:,1+nM_sub+[0:tmp_n_M_sub-1]) = tmp_M_x_u_xxM___(:,:,1+[0:tmp_n_M_sub-1]);
nM_sub = nM_sub + tmp_n_M_sub;
nbatch = nbatch + 1;
if nM_sub>=n_M_sub; flag_continue=0; end;
clear tmp_M_x_u_xxM___ tmp_n_M_sub;
end;%if  exist(fname_mrc,'file');
end;%while flag_continue;
%%%%;
for nM_sub=0:n_M_sub-1;
if (flag_verbose>0); if (mod(nM_sub,32)==0); disp(sprintf(' %% nstar %d nM_sub %d/%d',nstar,nM_sub,n_M_sub)); end; end;
M_k_p_wkMs___(:,1+nM_sub,1+nstar) = ...
interp_x_c_to_k_p_xxnufft( ...
 n_x_u_pack ...
,diameter_x_c ...
,n_x_u_pack ...
,diameter_x_c ...
,M_x_u_xxM___(:,:,1+nM_sub) ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
);
end;%for nM_sub=0:n_M_sub-1;
clear M_x_u_xxM___ ;
end;%for nstar=0:n_star-1;
%%%%%%%%;
save(fname_mat ...
     ,'n_M_sub','M_k_p_wkMs___' ...
     ,'index_nCTF_from_nM_Ms__' ...
     ,'index_nM_from_nCTF_s_' ...
     ,'Voltage_CTF_s_' ...
     ,'DefocusU_CTF_s_' ...
     ,'DefocusV_CTF_s_' ...
     ,'DefocusAngle_CTF_s_' ...
     ,'SphericalAberration_CTF_s_' ...
     ,'AmplitudeContrast_CTF_s_' ...
     );
end;%if (flag_recalc | ~exist(fname_mat,'file'));
%%%%%%%%;
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));

%%%%%%%%;
if flag_disp;
fname_fig_pre = sprintf('%s_jpg/M_x_u_xxMs____',dir_pm);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
disp(sprintf(' %% %s not found, creating',fname_fig_pre));
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 4; p_col = 6; np=0;
nM=0;
while np<p_row*p_col;
for nstar=0:n_star-1;
M_k_p_wk_ = M_k_p_wkMs___(:,1+nM,1+nstar);
M_x_c_xx_ = ...
interp_k_p_to_x_c_xxnufft( ...
 n_x_u_pack ...
,diameter_x_c ...
,n_x_u_pack ...
,diameter_x_c ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,M_k_p_wk_ ...
);
subplot(p_row,p_col,1+np);np=np+1;
Mlim_ = prctile(real(M_x_c_xx_),[  0,100],'all'); Mlim_  = mean(Mlim_) + 1.25*0.5*diff(Mlim_)*[-1,+1];
imagesc_c(n_x_u_pack,x_u_0_,n_x_u_pack,x_u_1_,real(M_x_c_xx_),Mlim_,colormap_beach());
axis image; axisnotick; title(sprintf('nM %d nstar %d',nM,nstar),'Interpreter','none');
end;%for nstar=0:n_star-1;
nM=nM+1;
end;%while np<p_row*p_col;
sgtitle(fname_fig_pre,'Interpreter','none');
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
end;%if (~exist(fname_fig_jpg,'file'));
close(gcf);
end;%if flag_disp;
%%%%%%%%;

%%%%%%%%;
n_x_M_u = n_x_u;
n_M = n_M_sub*n_star;
M_k_p_wkM__ = reshape(M_k_p_wkMs___,[n_w_sum,n_M]);
index_nCTF_from_nM_ = zeros(n_M,1);
n_CTF=0;
for nstar=0:n_star-1;
index_nCTF_from_nM_(1+nstar*n_M_sub+[0:n_M_sub-1]) = n_CTF + index_nCTF_from_nM_Ms__(:,1+nstar);
n_CTF = n_CTF + numel(index_nM_from_nCTF_s_{1+nstar});
end;%for nstar=0:n_star-1;
Voltage_CTF_ = [];
DefocusU_CTF_ = [];
DefocusV_CTF_ = [];
DefocusAngle_CTF_ = [];
SphericalAberration_CTF_ = [];
AmplitudeContrast_CTF_ = [];
nCTF=0;
for nstar=0:n_star-1;
Voltage_CTF_ = cat(1,Voltage_CTF_,Voltage_CTF_s_{1+nstar});
DefocusU_CTF_ = cat(1,DefocusU_CTF_,DefocusU_CTF_s_{1+nstar});
DefocusV_CTF_ = cat(1,DefocusV_CTF_,DefocusV_CTF_s_{1+nstar});
DefocusAngle_CTF_ = cat(1,DefocusAngle_CTF_,DefocusAngle_CTF_s_{1+nstar});
SphericalAberration_CTF_ = cat(1,SphericalAberration_CTF_,SphericalAberration_CTF_s_{1+nstar});
AmplitudeContrast_CTF_ = cat(1,AmplitudeContrast_CTF_,AmplitudeContrast_CTF_s_{1+nstar});
end;%for nstar=0:n_star-1;
%%%%%%%%;

if flag_disp>1;
%%%%%%%%;
% step through the images. ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 3; p_col = 6; np=0;
for nM=0:n_M-1;
subplot(p_row,p_col,1+np); cla;
M_k_p_wk_ = M_k_p_wkM__(:,1+nM);
M_x_c_xx_ = ...
interp_k_p_to_x_c_xxnufft( ...
 n_x_u_pack ...
,diameter_x_c ...
,n_x_u_pack ...
,diameter_x_c ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,M_k_p_wk_ ...
);
Mlim_ = prctile(real(M_x_c_xx_),[  0,100],'all'); Mlim_  = mean(Mlim_) + 1.25*0.5*diff(Mlim_)*[-1,+1];
imagesc_c(n_x_u_pack,x_u_0_,n_x_u_pack,x_u_1_,real(M_x_c_xx_),Mlim_,colormap_beach);
axis image; axisnotick; title(sprintf('nM %d',nM));
drawnow();
np=np+1; if np>=p_row*p_col; np=0; end;
end;%for nM=0:n_M-1;
end;%if flag_disp>1;

%%%%%%%%;
fname_mat = sprintf('%s_mat/CTF_k_p_wkC__.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
CTF_k_p_wkC__ = zeros(n_w_sum,n_CTF);
for nCTF=0:n_CTF-1;
if (mod(nCTF,100)==0); disp(sprintf(' %% nCTF %d/%d',nCTF,n_CTF)); end;
CTF_Spherical_Aberration = SphericalAberration_CTF_(1+nCTF);% spherical aberration of the lens in mm ;
CTF_Spherical_Aberration=CTF_Spherical_Aberration*(10.0d0^7.0d0);% convert into Angstroms ;
CTF_Voltage_kV = Voltage_CTF_(1+nCTF);% voltage in kVolts ;
CTF_Voltage_1V=CTF_Voltage_kV*1000.0 ;% convert into Volts ;
CTF_lambda = 12.2643247/sqrt(CTF_Voltage_1V+CTF_Voltage_1V^2*0.978466d-6);% electron wavelength in Angstroms ;
CTF_Defocus_U = DefocusU_CTF_(1+nCTF);% defocus values (in Angstroms) ;
CTF_Defocus_V = DefocusV_CTF_(1+nCTF);% defocus values (in Angstroms) ;
CTF_Defocus_Angle = DefocusAngle_CTF_(1+nCTF);% angle of astigmatism ;
%CTF_Defocus_Angle = CTF_Defocus_Angle*pi/180.0d0;% convert into radians ; %<-- already in radians! make sure not to convert twice!;
CTF_Amplitude_Contrast = AmplitudeContrast_CTF_(1+nCTF);% CTF_Amplitude Contrast ;
tmp_w1=sqrt(1.0d0-CTF_Amplitude_Contrast^2);% weights for the amplitude and phase contrasts in CTF ;
tmp_w2=CTF_Amplitude_Contrast;% weights for the amplitude and phase contrasts in CTF ;
%  CTF_Object_Pixel_Size = CTF_Detector_Pixel_Size/CTF_Magnification;
CTF_Object_Pixel_Size = Pixel_Spacing;% pixel size of the scanner in physical space in Angstroms ;
CTF_lambda_per_box = CTF_lambda/(n_x_M_u*CTF_Object_Pixel_Size);% n_x_M_u*CTF_Object_Pixel_Size is the box size in Angstroms ;
%%%%;
na=0;
for nk = 0:n_k_p_r-1;
for nw=0:n_w_(1+nk)-1;
tmp_theta = (2.0d0*pi*nw)/n_w_(1+nk);
tmp_k_c_1 = (2.0d0*pi)*k_p_r_(1+nk)*cos(tmp_theta);
tmp_k_c_2 = (2.0d0*pi)*k_p_r_(1+nk)*sin(tmp_theta);
tmp_ctf_value = niko_ctf(CTF_Spherical_Aberration,CTF_lambda,tmp_w1,tmp_w2,CTF_Defocus_U,CTF_Defocus_V,CTF_Defocus_Angle,CTF_lambda_per_box/pi,tmp_k_c_1,tmp_k_c_2);
CTF_k_p_wkC__(1+na,1+nCTF) = -tmp_ctf_value;
na = na+1;
end;%for nw=0:n_w_(1+nk)-1;
end;%for nk = 0:n_k_p_r-1;
%%%%;
end;%for nCTF=0:n_CTF-1;
%%%%%%%%;
CTF_k_p_r_kC__ = zeros(n_k_p_r,n_CTF);
for nCTF=0:n_CTF-1;
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
CTF_k_p_r_kC__(1+nk_p_r,1+nCTF) = mean(CTF_k_p_wkC__(1+tmp_index_,1+nCTF));
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nCTF=0:n_CTF-1;
CTF_avg_k_p_wk_ = mean(CTF_k_p_wkC__,2);
%imagesc_p(n_k_p_r,k_p_r_,n_w_,sum(n_w_),real(CTF_avg_k_p_wk_(:)),[-1,+1],colormap_beach());
CTF_avg_k_p_r_k_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
CTF_avg_k_p_r_k_(1+nk_p_r) = mean(CTF_avg_k_p_wk_(1+tmp_index_));
end;%for nk_p_r=0:n_k_p_r-1;
CTF_k_p_r_kM__ = CTF_k_p_r_kC__(:,1+index_nCTF_from_nM_(1:n_M));
SCTF_ = svd(CTF_k_p_r_kM__);
n_CTF_rank = max(find(SCTF_/max(SCTF_)>tolerance_master));
[UCTF_kc__,SCTF_c__,VCTF_Mc__] = svds(CTF_k_p_r_kM__,n_CTF_rank);
VSCTF_Mc__ = VCTF_Mc__*SCTF_c__;
%%%%%%%%;
% Now plot out some of the CTF-functions for varying anisotropy. ;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/CTF_k_p_sample__',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;
[tmp_anisotropy_,CTF_anisotropy_index_] = sort(DefocusU_CTF_ - DefocusV_CTF_,'ascend'); CTF_anisotropy_index_ = CTF_anisotropy_index_ - 1;
for nl=0:15-1;
subplot(3,5,1+nl);
tmp_nCTF = max(0,min(n_CTF-1,floor(n_CTF*nl/(15-1)))); nCTF = CTF_anisotropy_index_(1+tmp_nCTF);
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,CTF_k_p_wkC__(:,1+nCTF),[-1,+1],colormap_beach());
axis image; axisnotick;
title(sprintf('nCTF %d anisotropy %0.2f',nCTF,tmp_anisotropy_(1+tmp_nCTF)));
end;%for nl=0:15-1;
sgtitle(sprintf('CTF_k_p_wkC__'),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
% Now determine the CTF cross correlation. ;
% This depends  on index_nCTF_from_nM_. ;
%%%%%%%%;
tmp_CTF_avg_k_p_wk_ = mean(CTF_k_p_wkC__(:,1+index_nCTF_from_nM_(1+(0:n_M-1))),2);
tmp_CTF_avg_k_p_r_k_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
tmp_CTF_avg_k_p_r_k_(1+nk_p_r) = mean(tmp_CTF_avg_k_p_wk_(1+tmp_index_));
end;%for nk_p_r=0:n_k_p_r-1;
CTF_k_p_r_xavg__ = tmp_CTF_avg_k_p_r_k_ * transpose(tmp_CTF_avg_k_p_r_k_);
CTF_k_p_r_xcor__ = CTF_k_p_r_kC__(:,1+index_nCTF_from_nM_(1+(0:n_M-1))) * transpose(CTF_k_p_r_kC__(:,1+index_nCTF_from_nM_(1+(0:n_M-1)))) / n_M;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/CTF_k_p_xcor__',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figmed;
figbeach();
subplot(1,2,1); imagesc(CTF_k_p_r_xavg__); axis image; axisnotick; title('CTF_k_p_r_xavg__','Interpreter','none');
subplot(1,2,2); imagesc(CTF_k_p_r_xcor__); axis image; axisnotick; title('CTF_k_p_r_xcor__','Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
save(fname_mat ...
     ,'n_CTF' ...
     ,'index_nCTF_from_nM_' ...
     ,'CTF_k_p_wkC__' ...
     ,'CTF_k_p_r_kC__' ...
     ,'CTF_k_p_r_kM__' ...
     ,'n_CTF_rank' ...
     ,'SCTF_','UCTF_kc__','VSCTF_Mc__' ...
     ,'CTF_avg_k_p_wk_' ...
     ,'CTF_avg_k_p_r_k_' ...
     ,'CTF_k_p_r_xavg__' ...
     ,'CTF_k_p_r_xcor__' ...
     );
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/CTF_k_p_r_xxxx__',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;
clf;
colormap(colormap_beach());
subplot(1,2,1);
imagesc(CTF_k_p_r_xavg__,[0,1]);
axis image;
set(gca,'XTick',1:n_k_p_r,'XTickLabel',k_p_r_); xtickangle(90);
set(gca,'YTick',1:n_k_p_r,'YTickLabel',k_p_r_);
title('CTF_k_p_r_xavg__ cross correlation','Interpreter','none');
xlabel('k_p_r','Interpreter','none');
ylabel('k_p_r','Interpreter','none');
subplot(1,2,2);
imagesc(CTF_k_p_r_xcor__,[0,1]);
axis image;
set(gca,'XTick',1:n_k_p_r,'XTickLabel',k_p_r_); xtickangle(90);
set(gca,'YTick',1:n_k_p_r,'YTickLabel',k_p_r_);
title('CTF_k_p_r_xcor__ cross correlation','Interpreter','none');
xlabel('k_p_r','Interpreter','none');
ylabel('k_p_r','Interpreter','none');
figbig;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

%%%%%%%%;
% Now use ideal principal-modes to align images to volumes. ;
%%%%%%%%;
tmp_delta_r_max = 0.030;
tmp_delta_r_upb = 0.250;
euler_polar_a_empi_Mv__ = zeros(n_M,n_volume);
euler_azimu_b_empi_Mv__ = zeros(n_M,n_volume);
euler_gamma_z_empi_Mv__ = zeros(n_M,n_volume);
image_delta_x_empi_Mv__ = zeros(n_M,n_volume);
image_delta_y_empi_Mv__ = zeros(n_M,n_volume);
corr_a_k_Y_empi_v_ = zeros(n_volume,1);
image_X_value_empi_Mv__ = zeros(n_M,n_volume);
for nvolume=0:n_volume-1;
fname_mat = sprintf('%s_mat/a_k_Y_reco_%s_from_M__.mat',dir_pm,str_name_volume_{1+nvolume});
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
parameter = struct('type','parameter');
parameter.tolerance_master = tolerance_master;
parameter.n_iteration = 8;
parameter.delta_r_max = tmp_delta_r_max;
parameter.n_delta_v_requested = 24;
parameter.delta_r_upb = tmp_delta_r_upb;
parameter.fname_align_a_k_Y_pre = sprintf('%s_mat/a_k_Y_reco_%s_from_M__',dir_pm,str_name_volume_{1+nvolume});
parameter.flag_recalc = flag_recalc;
n_cluster=[];
index_ncluster_from_nCTF_=[];
[ ...
 parameter ...
,a_k_Y_reco_yki__ ...
,corr_a_k_Y_i_ ...
,euler_polar_a_Mi__ ...
,euler_azimu_b_Mi__ ...
,euler_gamma_z_Mi__ ...
,image_delta_x_acc_Mi__ ...
,image_delta_y_acc_Mi__ ...
,image_delta_x_upd_Mi__ ...
,image_delta_y_upd_Mi__ ...
,flag_image_delta_upd_Mi__ ...
,image_I_value_Mi__ ...
,image_X_value_Mi__ ...
,image_S_index_Mi__ ...
] = ...
pm_align_M_k_p_to_a_k_Y_2( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_ ...
,n_M ...
,M_k_p_wkM__ ...
,n_cluster ...
,index_ncluster_from_nCTF_ ...
,n_CTF ...
,index_nCTF_from_nM_ ...
,CTF_k_p_r_kC__ ...
,l_max_ ...
,a_k_Y_quad_ykv__(:,1+nvolume) ...
);
if ~exist(fname_mat,'file');
save(fname_mat ...
     ,'parameter' ...
     ,'a_k_Y_reco_yki__' ...
     ,'corr_a_k_Y_i_' ...
     ,'euler_polar_a_Mi__' ...
     ,'euler_azimu_b_Mi__' ...
     ,'euler_gamma_z_Mi__' ...
     ,'image_delta_x_acc_Mi__' ...
     ,'image_delta_y_acc_Mi__' ...
     ,'image_delta_x_upd_Mi__' ...
     ,'image_delta_y_upd_Mi__' ...
     ,'flag_image_delta_upd_Mi__' ...
     ,'image_I_value_Mi__' ...
     ,'image_X_value_Mi__' ...
     ,'image_S_index_Mi__' ...
     );
end;%if ~exist(fname_mat,'file');
clear a_k_Y_reco_yki__ corr_a_k_Y_i_ euler_polar_a_Mi__ euler_azimu_b_Mi__ euler_gamma_z_Mi__ image_delta_x_acc_Mi__ image_delta_y_acc_Mi__ image_delta_x_upd_Mi__ image_delta_y_upd_Mi__ flag_image_delta_upd_Mi__ image_I_value_Mi__ image_X_value_Mi__ image_S_index_Mi__ ;
end;%if (~exist(fname_mat,'file'));
%%%%%%%%;
fname_mat = sprintf('%s_mat/a_k_Y_reco_%s_from_M__.mat',dir_pm,str_name_volume_{1+nvolume});
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
tmp_ = load(fname_mat);
euler_polar_a_empi_Mv__(:,1+nvolume) = tmp_.euler_polar_a_Mi__(:,end);
euler_azimu_b_empi_Mv__(:,1+nvolume) = tmp_.euler_azimu_b_Mi__(:,end);
euler_gamma_z_empi_Mv__(:,1+nvolume) = tmp_.euler_gamma_z_Mi__(:,end);
image_delta_x_empi_Mv__(:,1+nvolume) = tmp_.image_delta_x_acc_Mi__(:,end) + tmp_.image_delta_x_upd_Mi__(:,end);
image_delta_y_empi_Mv__(:,1+nvolume) = tmp_.image_delta_y_acc_Mi__(:,end) + tmp_.image_delta_y_upd_Mi__(:,end);
corr_a_k_Y_empi_v_(1+nvolume) = tmp_.corr_a_k_Y_i_(end);
image_X_value_empi_Mv__(:,1+nvolume) = tmp_.image_X_value_Mi__(:,end-1);
clear tmp_;
end;%if ( exist(fname_mat,'file'));
end;%for nvolume=0:n_volume-1;
%%%%%%%%;

%%%%%%%%;
% Now use ideal principal-modes to align class averages to volumes. ;
%%%%%%%%;
tmp_delta_r_max = 0.030;
tmp_delta_r_upb = 0.250;
euler_polar_a_clas_Ov__ = zeros(n_O,n_volume);
euler_azimu_b_clas_Ov__ = zeros(n_O,n_volume);
euler_gamma_z_clas_Ov__ = zeros(n_O,n_volume);
image_delta_x_clas_Ov__ = zeros(n_O,n_volume);
image_delta_y_clas_Ov__ = zeros(n_O,n_volume);
corr_a_k_Y_clas_v_ = zeros(n_volume,1);
image_X_value_clas_Ov__ = zeros(n_O,n_volume);
for nvolume=0:n_volume-1;
fname_mat = sprintf('%s_mat/a_k_Y_reco_%s_from_O__.mat',dir_pm,str_name_volume_{1+nvolume});
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
parameter = struct('type','parameter');
parameter.tolerance_master = tolerance_master;
parameter.n_iteration = 8;
parameter.delta_r_max = tmp_delta_r_max;
parameter.n_delta_v_requested = 24;
parameter.delta_r_upb = tmp_delta_r_upb;
parameter.fname_align_a_k_Y_pre = sprintf('%s_mat/a_k_Y_reco_%s_from_O__',dir_pm,str_name_volume_{1+nvolume});
parameter.flag_recalc = flag_recalc;
[ ...
 parameter ...
,a_k_Y_reco_yki__ ...
,corr_a_k_Y_i_ ...
,euler_polar_a_Mi__ ...
,euler_azimu_b_Mi__ ...
,euler_gamma_z_Mi__ ...
,image_delta_x_acc_Mi__ ...
,image_delta_y_acc_Mi__ ...
,image_delta_x_upd_Mi__ ...
,image_delta_y_upd_Mi__ ...
,flag_image_delta_upd_Mi__ ...
,image_I_value_Mi__ ...
,image_X_value_Mi__ ...
,image_S_index_Mi__ ...
] = ...
pm_align_M_k_p_to_a_k_Y_2( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_ ...
,n_O ...
,O_k_p_wkO__ ...
,[] ...
,[] ...
,1 ...
,zeros(n_O,1) ...
,ones(n_k_p_r,1) ...
,l_max_ ...
,a_k_Y_quad_ykv__(:,1+nvolume) ...
);
if ~exist(fname_mat,'file');
save(fname_mat ...
     ,'parameter' ...
     ,'a_k_Y_reco_yki__' ...
     ,'corr_a_k_Y_i_' ...
     ,'euler_polar_a_Mi__' ...
     ,'euler_azimu_b_Mi__' ...
     ,'euler_gamma_z_Mi__' ...
     ,'image_delta_x_acc_Mi__' ...
     ,'image_delta_y_acc_Mi__' ...
     ,'image_delta_x_upd_Mi__' ...
     ,'image_delta_y_upd_Mi__' ...
     ,'flag_image_delta_upd_Mi__' ...
     ,'image_I_value_Mi__' ...
     ,'image_X_value_Mi__' ...
     ,'image_S_index_Mi__' ...
     );
end;%if ~exist(fname_mat,'file');
clear a_k_Y_reco_yki__ corr_a_k_Y_i_ euler_polar_a_Mi__ euler_azimu_b_Mi__ euler_gamma_z_Mi__ image_delta_x_acc_Mi__ image_delta_y_acc_Mi__ image_delta_x_upd_Mi__ image_delta_y_upd_Mi__ flag_image_delta_upd_Mi__ image_I_value_Mi__ image_X_value_Mi__ image_S_index_Mi__ ;
end;%if (~exist(fname_mat,'file'));
%%%%%%%%;
fname_mat = sprintf('%s_mat/a_k_Y_reco_%s_from_O__.mat',dir_pm,str_name_volume_{1+nvolume});
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
tmp_ = load(fname_mat);
euler_polar_a_clas_Ov__(:,1+nvolume) = tmp_.euler_polar_a_Mi__(:,end);
euler_azimu_b_clas_Ov__(:,1+nvolume) = tmp_.euler_azimu_b_Mi__(:,end);
euler_gamma_z_clas_Ov__(:,1+nvolume) = tmp_.euler_gamma_z_Mi__(:,end);
image_delta_x_clas_Ov__(:,1+nvolume) = tmp_.image_delta_x_acc_Mi__(:,end) + tmp_.image_delta_x_upd_Mi__(:,end);
image_delta_y_clas_Ov__(:,1+nvolume) = tmp_.image_delta_y_acc_Mi__(:,end) + tmp_.image_delta_y_upd_Mi__(:,end);
corr_a_k_Y_clas_v_(1+nvolume) = tmp_.corr_a_k_Y_i_(end);
image_X_value_clas_Ov__(:,1+nvolume) = tmp_.image_X_value_Mi__(:,end-1);
clear tmp_;
end;%if ( exist(fname_mat,'file'));
end;%for nvolume=0:n_volume-1;
%%%%%%%%;

if flag_disp>1;
%%%%%%%%;
% now step through the class averages post alignment. ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 4; p_col = 4; np=0; np2=0;
for nO=0:min(n_O-1,4*32);%for nO=0:n_O-1;
O_k_p_wk_ = O_k_p_wkO__(:,1+nO);
for nvolume=0:n_volume-1;
tmp_euler_polar_a = euler_polar_a_clas_Ov__(1+nO,1+nvolume);
tmp_euler_azimu_b = euler_azimu_b_clas_Ov__(1+nO,1+nvolume);
tmp_image_X_value = image_X_value_clas_Ov__(1+nO,1+nvolume);
nS = knnsearch([viewing_polar_a_all_,viewing_azimu_b_all_],[tmp_euler_polar_a,tmp_euler_azimu_b]) - 1;
S_k_p_wk_ = S_k_p_wkSv___(:,1+nS,1+nvolume);
S_x_c_xx_ = ...
interp_k_p_to_x_c_xxnufft( ...
 n_x_u_pack ...
,diameter_x_c ...
,n_x_u_pack ...
,diameter_x_c ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,S_k_p_wk_ ...
);
Slim_ = prctile(real(S_x_c_xx_),[  0,100],'all'); Slim_  = mean(Slim_) + 1.25*0.5*diff(Slim_)*[-1,+1];
tmp_O_k_p_wk_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,O_k_p_wk_,+image_delta_x_clas_Ov__(1+nO,1+nvolume),+image_delta_y_clas_Ov__(1+nO,1+nvolume));
tmp_O_k_p_wk_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,tmp_O_k_p_wk_,-euler_gamma_z_clas_Ov__(1+nO,1+nvolume));
O_x_c_xx_ = ...
interp_k_p_to_x_c_xxnufft( ...
 n_x_u_pack ...
,diameter_x_c ...
,n_x_u_pack ...
,diameter_x_c ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,tmp_O_k_p_wk_ ...
);
Olim_ = prctile(real(O_x_c_xx_),[  0,100],'all'); Olim_  = mean(Olim_) + 1.25*0.5*diff(Olim_)*[-1,+1];
subplot(p_row,p_col,1+np); cla; np=np+1; if np>=p_row*p_col; np=0; end;
imagesc_c(n_x_u_pack,x_u_0_,n_x_u_pack,x_u_1_,real(O_x_c_xx_),Olim_,colormap_beach);
axis image; axisnotick; title(sprintf('nO %d nv %d',nO,nvolume));
subplot(p_row,p_col,1+np); cla; np=np+1; if np>=p_row*p_col; np=0; end;
imagesc_c(n_x_u_pack,x_u_0_,n_x_u_pack,x_u_1_,real(S_x_c_xx_),Slim_,colormap_beach);
axis image; axisnotick; title(sprintf('nS %d nv %d [X %+0.2f]',nS,nvolume,tmp_image_X_value));
clear tmp_O_k_p_wk_ O_x_c_xx_ ;
clear S_k_p_wk_ S_x_c_xx_ ;
drawnow();
if np==0;
fname_fig_pre = sprintf('%s_jpg/O_x_u_xxOv____aligned_%.3d',dir_pm,np2);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
disp(sprintf(' %% %s not found, creating',fname_fig_pre));
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
np2=np2+1;
end;%if np==0;
end;%for nvolume=0:n_volume-1;
end;%for nO=0:n_O-1;
%%%%%%%%;
end;%if flag_disp>1;

if flag_disp>1;
%%%%%%%%;
% now step through the images post alignment. ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 4; p_col = 4; np=0; np2=0;
for nM=0:min(n_M-1,4*32);%for nM=0:n_M-1;
M_k_p_wk_ = M_k_p_wkM__(:,1+nM);
for nvolume=0:n_volume-1;
tmp_euler_polar_a = euler_polar_a_empi_Mv__(1+nM,1+nvolume);
tmp_euler_azimu_b = euler_azimu_b_empi_Mv__(1+nM,1+nvolume);
tmp_image_X_value = image_X_value_empi_Mv__(1+nM,1+nvolume);
nS = knnsearch([viewing_polar_a_all_,viewing_azimu_b_all_],[tmp_euler_polar_a,tmp_euler_azimu_b]) - 1;
S_k_p_wk_ = S_k_p_wkSv___(:,1+nS,1+nvolume);
S_x_c_xx_ = ...
interp_k_p_to_x_c_xxnufft( ...
 n_x_u_pack ...
,diameter_x_c ...
,n_x_u_pack ...
,diameter_x_c ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,S_k_p_wk_ ...
);
Slim_ = prctile(real(S_x_c_xx_),[  0,100],'all'); Slim_  = mean(Slim_) + 1.25*0.5*diff(Slim_)*[-1,+1];
tmp_M_k_p_wk_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p_wk_,+image_delta_x_empi_Mv__(1+nM,1+nvolume),+image_delta_y_empi_Mv__(1+nM,1+nvolume));
tmp_M_k_p_wk_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,tmp_M_k_p_wk_,-euler_gamma_z_empi_Mv__(1+nM,1+nvolume));
M_x_c_xx_ = ...
interp_k_p_to_x_c_xxnufft( ...
 n_x_u_pack ...
,diameter_x_c ...
,n_x_u_pack ...
,diameter_x_c ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,tmp_M_k_p_wk_ ...
);
Mlim_ = prctile(real(M_x_c_xx_),[  0,100],'all'); Mlim_  = mean(Mlim_) + 1.25*0.5*diff(Mlim_)*[-1,+1];
subplot(p_row,p_col,1+np); cla; np=np+1; if np>=p_row*p_col; np=0; end;
imagesc_c(n_x_u_pack,x_u_0_,n_x_u_pack,x_u_1_,real(M_x_c_xx_),Mlim_,colormap_beach);
axis image; axisnotick; title(sprintf('nM %d nv %d',nM,nvolume));
subplot(p_row,p_col,1+np); cla; np=np+1; if np>=p_row*p_col; np=0; end;
imagesc_c(n_x_u_pack,x_u_0_,n_x_u_pack,x_u_1_,real(S_x_c_xx_),Slim_,colormap_beach);
axis image; axisnotick; title(sprintf('nS %d nv %d [X %+0.2f]',nS,nvolume,tmp_image_X_value));
clear M_tmp_M_k_p_wk_ M_x_c_xx_ ;
clear S_k_p_wk_ S_x_c_xx_ ;
drawnow();
if np==0;
fname_fig_pre = sprintf('%s_jpg/M_x_u_xxMs____aligned_%.3d',dir_pm,np2);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
disp(sprintf(' %% %s not found, creating',fname_fig_pre));
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
np2=np2+1;
end;%if np==0;
end;%for nvolume=0:n_volume-1;
end;%for nM=0:n_M-1;
%%%%%%%%;
end;%if flag_disp>1;

figure(1+nf);nf=nf+1;clf;figmed;
markersize_use = 12;
subplot(1,2,1);
hold on;
plot([0,1],[0,1],'k');
tmp_index_=0:n_O-1;
plot(image_X_value_clas_Ov__(1+tmp_index_,1),image_X_value_clas_Ov__(1+tmp_index_,2),'k.');
hold off;
xlabel('WT'); ylabel('MUT');
xlim([0,1]); ylim([0,1]); set(gca,'XTick',[0,1]); set(gca,'YTick',[0,1]); axis square;
title('class averages (best X)');
subplot(1,2,2);
hold on;
plot([0,1],[0,1],'k');
tmp_index_ = 0*n_M_sub + [0:n_M_sub-1];
plot(image_X_value_empi_Mv__(1+tmp_index_,1),image_X_value_empi_Mv__(1+tmp_index_,2),'k.','MarkerSize',markersize_use);
tmp_index_ = 1*n_M_sub + [0:n_M_sub-1];
plot(image_X_value_empi_Mv__(1+tmp_index_,1),image_X_value_empi_Mv__(1+tmp_index_,2),'r.','MarkerSize',markersize_use);
title('images (best X)');
xlabel('WT'); ylabel('MUT');
xlim([0,.4]); ylim([0,.4]); set(gca,'XTick',[0,.4]); set(gca,'YTick',[0,.4]); axis square;
fname_fig_pre = sprintf('%s_jpg/X_x_u_xxMs____aligned_scatterplot',dir_pm);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
disp(sprintf(' %% %s not found, creating',fname_fig_pre));
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);

return;

