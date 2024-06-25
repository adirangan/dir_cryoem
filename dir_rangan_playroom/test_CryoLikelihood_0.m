%%%%%%%%;
% look at the template-manifold for rib80s. ;
%%%%%%%%;
str_thisfunction = 'test_CryoLikelihood_0';

%%%%%%%%;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
%%%%%%%%;

global_parameter=[];fname_prefix='rib80s_x0';dir_nopath_data_star='rib80s';Pixel_Spacing=1.34;fname_nopath_volume='emd_2660.mrc';fname_nopath_star='shiny_2sets.star';

if isempty(global_parameter); global_parameter = struct('type','parameter'); end;
if (~isfield(global_parameter,'flag_verbose')); global_parameter.flag_verbose = 1; end; %<-- parameter_bookmark. ;
flag_verbose = global_parameter.flag_verbose;
if (~isfield(global_parameter,'flag_disp')); global_parameter.flag_disp = 1; end; %<-- parameter_bookmark. ;
flag_disp = global_parameter.flag_disp; nf=0;
if (~isfield(global_parameter,'tolerance_master')); global_parameter.tolerance_master = 1e-6; end; %<-- parameter_bookmark. ;
tolerance_master = global_parameter.tolerance_master;

if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

%%%%%%%%;
fname_prefix_xfix = sprintf('%s',fname_prefix);
dir_pm = sprintf('/%s/rangan/dir_cryoem/dir_%s/dir_pm',string_root,fname_prefix_xfix);
if (~exist(sprintf('%s_mat',dir_pm),'dir')); disp(sprintf(' %% mkdir %s_mat',dir_pm)); mkdir(sprintf('%s_mat',dir_pm)); end;
if (~exist(sprintf('%s_jpg',dir_pm),'dir')); disp(sprintf(' %% mkdir %s_jpg',dir_pm)); mkdir(sprintf('%s_jpg',dir_pm)); end;
dir_data_star = sprintf('/%s/rangan/dir_cryoem/dir_%s',string_root,dir_nopath_data_star);
%%%%%%%%;

%%%%%%%%;
% Load the volume and extract a few templates. ;
%%%%%%%%;
fname_emd = sprintf('%s/%s',dir_data_star,fname_nopath_volume);
a_x_u_load_ = cast(ReadMRC(fname_emd),'double');
n_x_u = size(a_x_u_load_,1);
n_S = 3;
S_x_c_xxS___ = zeros(n_x_u,n_x_u,n_S);
nS=0;
S_x_c_xxS___(:,:,1+nS) = squeeze(sum(a_x_u_load_,[1])); nS=nS+1;
S_x_c_xxS___(:,:,1+nS) = squeeze(sum(a_x_u_load_,[2])); nS=nS+1;
S_x_c_xxS___(:,:,1+nS) = squeeze(sum(a_x_u_load_,[3])); nS=nS+1;
n_x_S_u = size(S_x_c_xxS___,1);
assert(n_x_S_u==size(S_x_c_xxS___,2));
clear a_x_u_load_;
%%%%%%%%;
half_diameter_x_c = 1.0d0;
diameter_x_c = 2.0d0*half_diameter_x_c;
x_p_r_max = 1.0;
x_u_0_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u);
x_u_1_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u);
[x_u_0___,x_u_1___] = ndgrid(x_u_0_,x_u_1_); n_xx_u = n_x_u^2;
xx_u_weight_ = (2*x_p_r_max/n_x_u)^2;
%%%%%%%%;
if flag_disp>0;
figure(1+nf);nf=nf+1;clf;figmed;figbeach();
p_row = 1; p_col = n_S; np=0;
for nS=0:n_S-1;
subplot(p_row,p_col,1+np);np=np+1;
imagesc_c(n_x_S_u,x_u_0_,n_x_S_u,x_u_1_,S_x_c_xxS___(:,:,1+nS));
axis image; axisnotick; title(sprintf('nS %0.2d',nS),'Interpreter','none');
end;%for nS=0:n_S-1;
fname_fig_pre = sprintf('%s_jpg/test_CryoLikelihood_Templates',dir_pm);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if (flag_verbose>0); disp(sprintf(' %% writing %s',fname_fig_jpg)); end;
print('-djpeg',fname_fig_jpg);
end;%if flag_disp;
%%%%%%%%;

%%%%%%%%;
% Load a few images. ;
%%%%%%%%;
n_M = 4;
[ ...
 M_x_c_xxM___ ...
,index_nCTF_from_nM_ ...
,index_nM_from_nCTF_ ...
,Voltage_CTF_ ...
,DefocusU_CTF_ ...
,DefocusV_CTF_ ...
,DefocusAngle_CTF_ ...
,SphericalAberration_CTF_ ...
,AmplitudeContrast_CTF_ ...
] = ...
rlnImageName_from_star_1( ...
 dir_data_star ...
,fname_nopath_star ...
,n_M ...
);
%%%%%%%%;
% Remove any edge artefacts, mean center and normalize each image. ;
%%%%%%%%;
disp(sprintf(' %% Removing edge-artefacts'));
n_M_ext_ = zeros(n_M,1);
for nM=0:n_M-1;
if (mod(nM,128)==0); disp(sprintf(' %% nM %d/%d',nM,n_M)); end;
n_pixel = 4; edge_tolerance = 0.5; n_edge_overshoot = 8; rseed = 0;
[M_x_c_xxM___(:,:,1+nM),n_M_ext_(1+nM)] = image_replace_edge_artefact_0(M_x_c_xxM___(:,:,1+nM),4,0.5,2,0);
end;%for nM=0:n_M-1;
n_x_M_u = size(M_x_c_xxM___,1);
assert(n_x_M_u==size(M_x_c_xxM___,2));
disp(sprintf(' %% typical edge-artefact covers %0.6f = (%0.6f)^2 of image.',median(n_M_ext_/n_x_M_u^2),median(sqrt(n_M_ext_/n_x_M_u^2))));
disp(sprintf(' %% edge-artefacts detected in %d/%d images.',numel(find(n_M_ext_>0)),n_M));
%%%%%%%%;
if flag_disp>0;
figure(1+nf);nf=nf+1;clf;figmed;figbeach();
p_row = 1; p_col = n_M; np=0;
for nM=0:n_M-1;
subplot(p_row,p_col,1+np);np=np+1;
imagesc_c(n_x_M_u,x_u_0_,n_x_M_u,x_u_1_,M_x_c_xxM___(:,:,1+nM));
axis image; axisnotick; title(sprintf('nM %0.2d',nM),'Interpreter','none');
end;%for nM=0:n_M-1;
fname_fig_pre = sprintf('%s_jpg/test_CryoLikelihood_Images',dir_pm);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if (flag_verbose>0); disp(sprintf(' %% writing %s',fname_fig_jpg)); end;
print('-djpeg',fname_fig_jpg);
end;%if flag_disp;
%%%%%%%%;

%%%%%%%%;
% set grid. ;
%%%%%%%%;
k_p_r_max = 3*48/(2*pi); k_eq_d = 1.0/(2*pi); str_TorL = 'L';
%%%%;
[ ...
 n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
] = ...
get_weight_3d_1( ...
 flag_verbose ...
,k_p_r_max ...
,k_eq_d ...
,str_TorL ...
);
%%%%;
n_w_max = 2*ceil(2*pi*k_p_r_max + 1);
n_w_0in_ = n_w_max*ones(n_k_p_r,1);
template_k_eq_d = -1;
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
 flag_verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,template_k_eq_d ...
,n_w_0in_ ...
);
n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);
%%%%%%%%;

%%%%%%%%;
% Now convert templates to S_k_p_wkS__. ;
%%%%%%%%;
dx = diameter_x_c/n_x_S_u;
S_k_p_wkS__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
S_x_c_ = squeeze(S_x_c_xxS___(:,:,1+nS));
S_k_p_ = interp_x_c_to_k_p_xxnufft(n_x_S_u,diameter_x_c,n_x_S_u,diameter_x_c,S_x_c_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_S_u^2)*dx^2 ;
S_k_p_wkS__(:,1+nS) = S_k_p_;
end;%for nS=0:n_S-1;
%%%%%%%%;
% Now convert templates to S_k_q_wkS__. ;
%%%%%%%%;
tmp_t = tic();
S_k_q_wkS__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
S_k_q_wkS__(:,1+nS) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,S_k_p_wkS__(:,1+nS));
end;%for nS=0:n_S-1;
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% S_k_q_wkS__: %0.3fs',tmp_t)); end;
%%%%%%%%;

%%%%%%%%;
% Now convert images to M_k_p_wkM__. ;
%%%%%%%%;
dx = diameter_x_c/n_x_M_u;
M_k_p_wkM__ = zeros(n_w_sum,n_M);
for nM=0:n_M-1;
M_x_c_ = squeeze(M_x_c_xxM___(:,:,1+nM));
M_k_p_ = interp_x_c_to_k_p_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,M_x_c_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_M_u^2)*dx^2 ;
M_k_p_wkM__(:,1+nM) = M_k_p_;
end;%for nM=0:n_M-1;
%%%%%%%%;
% Now convert images to M_k_q_wkM__. ;
%%%%%%%%;
tmp_t = tic();
M_k_q_wkM__ = zeros(n_w_sum,n_M);
for nM=0:n_M-1;
M_k_q_wkM__(:,1+nM) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_k_p_wkM__(:,1+nM));
end;%for nM=0:n_M-1;
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% M_k_q_wkM__: %0.3fs',tmp_t)); end;

%%%%%%%%;
% Construct CTF. ;
%%%%%%%%;
n_CTF = numel(index_nM_from_nCTF_);
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
%%%%%%%%;
CTF_k_p_r_k_ = mean(CTF_k_p_r_kC__,2);
%%%%%%%%;

%%%%%%%%;
% Prepare displacement grid. ;
%%%%%%%%;
d_p_r_max = 0.15; d_eq_d = d_p_r_max/48.0; str_TorL = 'L';
%%%%;
[ ...
 n_d_p_r ...
,d_p_r_ ...
,weight_3d_d_p_r_ ...
] = ...
get_weight_3d_1( ...
 flag_verbose ...
,d_p_r_max ...
,d_eq_d ...
,str_TorL ...
);
%%%%;
n_o_max = 2*ceil(2*pi*d_p_r_max/d_eq_d + 1); %<-- here 'o' stands for 'omega', which is the angular component of delta in polar-coordinates. ;
template_d_eq_d = d_eq_d;
[ ...
 n_o_ ...
,weight_2d_d_p_r_ ...
,weight_2d_od_ ...
,d_p_r_od_ ...
,d_p_o_od_ ...
,d_c_0_od_ ...
,d_c_1_od_ ...
] = ...
get_weight_2d_2( ...
 flag_verbose ...
,n_d_p_r ...
,d_p_r_ ...
,d_p_r_max ...
,template_d_eq_d ...
);
n_o_max = max(n_o_);
n_o_sum = sum(n_o_);
n_o_csum_ = cumsum([0;n_o_]);
%%%%%%%%;

%%%%%%%%;
% Prepare FTK. ;
%%%%%%%%;
delta_p_r_max = d_p_r_max;
delta_r_p = delta_p_r_max;
delta_r_s = delta_p_r_max/sqrt(2*log(1/delta_r_p));
delta_r_N = delta_p_r_max * (2*pi*k_p_r_max) / (pi*sqrt(2));
svd_eps = 1e-4;
n_delta_v_requested = n_o_sum;
tmp_t = tic();
FTK = ampmh_FTK_1(n_k_p_r,k_p_r_,k_p_r_max,delta_p_r_max,svd_eps,n_delta_v_requested,d_c_0_od_,d_c_1_od_);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% FTK: %0.3fs',tmp_t)); end;
disp(sprintf(' %% p-val %0.4f delta_p_r_max %0.6f sigma %0.4f N_pixel %0.4f --> FTK.n_svd_l %d, n_delta_v_requested %d',delta_r_p,delta_p_r_max,delta_r_s,delta_r_N,FTK.n_svd_l,n_delta_v_requested));
%%%%%%%%;
% Prepare principal-modes. ;
%%%%%%%%;
tmp_t = tic();
n_UX_rank = n_k_p_r-1; %<-- just to check dimensions.; 
UX_kn__ = zeros(n_k_p_r,n_UX_rank);
X_weight_r_ = zeros(n_k_p_r,1);
[ ...
 X_kk__ ...
,X_weight_r_ ...
] = ...
principled_marching_empirical_cost_matrix_0( ...
 n_k_p_r ...
,k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_ ...
,n_M ...
,M_k_p_wkM__ ...
);
[tmp_UX_kn__,tmp_SX_k__,tmp_VX_kn__] = svds(X_kk__,n_UX_rank);
UX_kn__(:,1:n_UX_rank) = tmp_UX_kn__(:,1:n_UX_rank);
pm_n_UX_rank = n_UX_rank; %<-- just to check dimension. ;
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% principled_marching_empirical_cost_matrix_0: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Prepare quasi-images. ;
%%%%%%%%;
n_UX_rank = n_k_p_r;
tmp_t = tic();
svd_VUXM_lwnM____ = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,n_M,M_k_q_wkM__,pm_n_UX_rank,UX_kn__,X_weight_r_);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% svd_VUXM_lwnM____: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Now calculate norms of the non-translated images. ;
%%%%%%%%;
tmp_t = tic();
UX_M_l2_dM__ = ampmh_UX_M_l2_dM__1(FTK,n_w_,n_M,pm_n_UX_rank,svd_VUXM_lwnM____);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% UX_M_l2_dM__: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Prepare CTF_UX_S_k_q_wnS__. ;
%%%%%%%%;
tmp_t = tic();
S_k_q_wSk___ = permute(reshape(S_k_q_wkS__,[n_w_max,n_k_p_r,n_S]),[1,3,2]);
CTF_UX_S_k_q_wnS__ = reshape(permute(reshape(reshape(S_k_q_wSk___,[n_w_max*n_S,n_k_p_r])*diag(X_weight_r_.*CTF_k_p_r_k_)*UX_kn__,[n_w_max,n_S,pm_n_UX_rank]),[1,3,2]),[n_w_max*pm_n_UX_rank,n_S]);
CTF_UX_S_l2_S_ = sum(abs(CTF_UX_S_k_q_wnS__).^2,1)/n_w_max;
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% CTF_UX_S_l2_S_: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Calculate ampmh_X_wSM___8. ;
%%%%%%%%;
parameter = struct('type','parameter');
parameter.flag_compute_I_value = 0;
tmp_t = tic();
[ ...
 parameter ...
,X_dwSM____ ...
,Y_dwSM____ ...
,I_value_dwSM____ ...
] = ...
ampmh_X_dwSM____8( ...
 parameter ...
,FTK ...
,n_w_max ...
,pm_n_UX_rank ...
,n_S ...
,CTF_UX_S_k_q_wnS__ ...
,CTF_UX_S_l2_S_ ...
,n_M ...
,svd_VUXM_lwnM____ ...
,UX_M_l2_dM__ ...
);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% X_dwSM____: %0.3fs',tmp_t)); end;
X_dwSM____ = real(X_dwSM____);
%%%%%%%%;
% visualizing the landscape of innerproducts for one image-template pair. ;
%%%%%%%%;
if flag_disp>0;
lambda = 0.25;
n_nw = 10;
%nw_ = round(linspace(0.4*pi,0.7*pi,1+n_nw)*n_w_max/(2*pi)); nw_ = nw_(1:n_nw);
nw_ = round(linspace(0.0*pi,2.0*pi,1+n_nw)*n_w_max/(2*pi)); nw_ = nw_(1:n_nw);
for nM=0:n_M-1;
for nS=0:n_S-1;
%for nM=2;
%for nS=1;
for nw=nw_;
X_d_ = real(X_dwSM____(:,1+nw,1+nS,1+nM));
Y_d_ = real(Y_dwSM____(:,1+nw,1+nS,1+nM));
Z_d_ = exp(-(Y_d_ - min(Y_d_,[],'all'))./(2*lambda^2));
Xlim_ = prctile(X_d_,[ 5,95]); Xlim_ = mean(Xlim_) + 1.25*0.5*diff(Xlim_)*[-1,+1];
Ylim_ = prctile(Y_d_,[ 5,95]); Ylim_ = mean(Ylim_) + 1.25*0.5*diff(Ylim_)*[-1,+1];
%Zlim_ = prctile(Z_d_,[ 5,95]); Zlim_ = mean(Zlim_) + 1.25*0.5*diff(Zlim_)*[-1,+1];
Zlim_ = prctile(Z_d_,[ 0,100]);
figure(1+nf);nf=nf+1;clf;figmed;
p_row = 1; p_col = 3; np=0;
subplot(p_row,p_col,1+np);np=np+1; imagesc_p(n_d_p_r,d_p_r_,n_o_,n_o_sum,X_d_,Xlim_,colormap_beach());
axis image; axisnotick; title('$\chi(A,B)$','Interpreter','latex');% tmp_c_ = colorbar; set(tmp_c_,'Ticks',[],'TickLength',0);
subplot(p_row,p_col,1+np);np=np+1; imagesc_p(n_d_p_r,d_p_r_,n_o_,n_o_sum,Y_d_,Ylim_,colormap_80s());
axis image; axisnotick; title('$\|A-B\|^{2}$','Interpreter','latex');% tmp_c_ = colorbar; set(tmp_c_,'Ticks',[],'TickLength',0);
subplot(p_row,p_col,1+np);np=np+1; imagesc_p(n_d_p_r,d_p_r_,n_o_,n_o_sum,Z_d_,Zlim_,colormap_81s());
axis image; axisnotick; title(sprintf('$\\exp(-0.5\\times\\|A-B\\|^{2}/\\lambda^{2}): \\lambda=%0.3f$',lambda),'Interpreter','latex');% tmp_c_ = colorbar; set(tmp_c_,'Ticks',[],'TickLength',0);
sgtitle(sprintf('K %0.2f/(2$\\pi$), d %0.3f, nM %.2d nS %.2d nw %.4d $(\\gamma %0.2f\\pi)$',round(2*pi*k_p_r_max),d_p_r_max,nM,nS,nw,2*pi*nw/n_w_max),'Interpreter','latex');
drawnow();
fname_fig_pre = sprintf('%s_jpg/test_CryoLikelihood_nS%d_nM%d_nw%d_lambda%.3d',dir_pm,nS,nM,nw,floor(100*lambda));
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if (flag_verbose>0); disp(sprintf(' %% writing %s',fname_fig_jpg)); end;
print('-djpeg',fname_fig_jpg);
close(gcf);
end;%for nw=nw_;
end;%for nS=0:n_S-1;
end;%for nM=0:n_M-1;
end;%if flag_disp;
%%%%%%%%;

%%%%%%%%;
% Repeat, but this time on a cartesian displacement-grid. ;
%%%%%%%%;
n_d = 1+8*8; n_dd = n_d*n_d;
d_c_0_d_ = linspace(-d_p_r_max/sqrt(2),+d_p_r_max/sqrt(2),n_d);
d_c_1_d_ = linspace(-d_p_r_max/sqrt(2),+d_p_r_max/sqrt(2),n_d);
[d_c_0_dd__,d_c_1_dd__] = ndgrid(d_c_0_d_,d_c_1_d_);
%%%%%%%%;
% Prepare FTK. ;
%%%%%%%%;
delta_p_r_max = d_p_r_max;
delta_r_p = delta_p_r_max;
delta_r_s = delta_p_r_max/sqrt(2*log(1/delta_r_p));
delta_r_N = delta_p_r_max * (2*pi*k_p_r_max) / (pi*sqrt(2));
svd_eps = 1e-4;
n_delta_v_requested = n_dd;
tmp_t = tic();
FTK = ampmh_FTK_1(n_k_p_r,k_p_r_,k_p_r_max,delta_p_r_max,svd_eps,n_delta_v_requested,d_c_0_dd__,d_c_1_dd__);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% FTK: %0.3fs',tmp_t)); end;
disp(sprintf(' %% p-val %0.4f delta_p_r_max %0.6f sigma %0.4f N_pixel %0.4f --> FTK.n_svd_l %d, n_delta_v_requested %d',delta_r_p,delta_p_r_max,delta_r_s,delta_r_N,FTK.n_svd_l,n_delta_v_requested));
%%%%%%%%;
% Prepare principal-modes. ;
%%%%%%%%;
tmp_t = tic();
n_UX_rank = n_k_p_r-1; %<-- just to check dimensions.; 
UX_kn__ = zeros(n_k_p_r,n_UX_rank);
X_weight_r_ = zeros(n_k_p_r,1);
[ ...
 X_kk__ ...
,X_weight_r_ ...
] = ...
principled_marching_empirical_cost_matrix_0( ...
 n_k_p_r ...
,k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_ ...
,n_M ...
,M_k_p_wkM__ ...
);
[tmp_UX_kn__,tmp_SX_k__,tmp_VX_kn__] = svds(X_kk__,n_UX_rank);
UX_kn__(:,1:n_UX_rank) = tmp_UX_kn__(:,1:n_UX_rank);
pm_n_UX_rank = n_UX_rank; %<-- just to check dimension. ;
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% principled_marching_empirical_cost_matrix_0: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Prepare quasi-images. ;
%%%%%%%%;
n_UX_rank = n_k_p_r;
tmp_t = tic();
svd_VUXM_lwnM____ = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,n_M,M_k_q_wkM__,pm_n_UX_rank,UX_kn__,X_weight_r_);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% svd_VUXM_lwnM____: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Now calculate norms of the non-translated images. ;
%%%%%%%%;
tmp_t = tic();
UX_M_l2_dM__ = ampmh_UX_M_l2_dM__1(FTK,n_w_,n_M,pm_n_UX_rank,svd_VUXM_lwnM____);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% UX_M_l2_dM__: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Prepare CTF_UX_S_k_q_wnS__. ;
%%%%%%%%;
tmp_t = tic();
S_k_q_wSk___ = permute(reshape(S_k_q_wkS__,[n_w_max,n_k_p_r,n_S]),[1,3,2]);
CTF_UX_S_k_q_wnS__ = reshape(permute(reshape(reshape(S_k_q_wSk___,[n_w_max*n_S,n_k_p_r])*diag(X_weight_r_.*CTF_k_p_r_k_)*UX_kn__,[n_w_max,n_S,pm_n_UX_rank]),[1,3,2]),[n_w_max*pm_n_UX_rank,n_S]);
CTF_UX_S_l2_S_ = sum(abs(CTF_UX_S_k_q_wnS__).^2,1)/n_w_max;
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% CTF_UX_S_l2_S_: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Calculate ampmh_X_wSM___8. ;
%%%%%%%%;
parameter = struct('type','parameter');
parameter.flag_compute_I_value = 0;
tmp_t = tic();
[ ...
 parameter ...
,X_dwSM____ ...
,Y_dwSM____ ...
,I_value_dwSM____ ...
] = ...
ampmh_X_dwSM____8( ...
 parameter ...
,FTK ...
,n_w_max ...
,pm_n_UX_rank ...
,n_S ...
,CTF_UX_S_k_q_wnS__ ...
,CTF_UX_S_l2_S_ ...
,n_M ...
,svd_VUXM_lwnM____ ...
,UX_M_l2_dM__ ...
);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% X_dwSM____: %0.3fs',tmp_t)); end;
X_dwSM____ = real(X_dwSM____);
%%%%%%%%;
% visualizing the landscape of innerproducts for one image-template pair. ;
%%%%%%%%;
if flag_disp>0;
n_nw = 10;
%nw_ = round(linspace(0.4*pi,0.7*pi,1+n_nw)*n_w_max/(2*pi)); nw_ = nw_(1:n_nw);
nw_ = round(linspace(0.0*pi,2.0*pi,1+n_nw)*n_w_max/(2*pi)); nw_ = nw_(1:n_nw);
%for nM=0:n_M-1;
%for nS=0:n_S-1;
for nM=0;
for nS=1;

X_dw__ = real(X_dwSM____(:,:,1+nS,1+nM));
Xlim_ = prctile(X_dw__,[ 5,95],'all'); Xlim_ = mean(Xlim_) + 1.25*0.5*diff(Xlim_)*[-1,+1];
Y_dw__ = real(Y_dwSM____(:,:,1+nS,1+nM));
Ylim_ = prctile(Y_dw__,[ 5,95],'all'); Ylim_ = mean(Ylim_) + 1.25*0.5*diff(Ylim_)*[-1,+1];
lambda = 0.10;
Z_dw__ = exp(-(Y_dw__ - min(Y_dw__,[],'all'))./(2*lambda^2));
Zlim_ = prctile(Z_dw__,[ 0,100],'all');
figure(1+nf);nf=nf+1;clf;
p_row = 1; p_col = 1; np=0;
subplot(p_row,p_col,1+np);np=np+1; 
percent_threshold = [95]; isosurface_f_x_u_0(reshape(Z_dw__,[n_d,n_d,n_w_max]),percent_threshold); axis equal; axisnotick;
title(sprintf('$\\exp(-0.5\\times\\|A-B\\|^{2}/\\lambda^{2}): \\%%=%0.3f$',0.95),'Interpreter','latex');% tmp_c_ = colorbar; set(tmp_c_,'Ticks',[],'TickLength',0);
sgtitle(sprintf('K %0.2f/(2$\\pi$), d %0.3f, nM %.2d nS %.2d $(\\gamma %0.2f\\pi)$',round(2*pi*k_p_r_max),d_p_r_max,nM,nS,2*pi*nw/n_w_max),'Interpreter','latex');
set(gcf,'Position',1+[0,0,768,1024]);
drawnow();
fname_fig_pre = sprintf('%s_jpg/test_CryoLikelihood_nS%d_nM%d_nwx_lambda%.3d_d_c',dir_pm,nS,nM,floor(100*lambda));
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if (flag_verbose>0); disp(sprintf(' %% writing %s',fname_fig_jpg)); end;
print('-djpeg',fname_fig_jpg);
close(gcf);

%%%%;
for lambda = [0.10,0.15,0.20,0.25];
for nw=58;%for nw=nw_;
X_d_ = real(X_dwSM____(:,1+nw,1+nS,1+nM));
Y_d_ = real(Y_dwSM____(:,1+nw,1+nS,1+nM));
Z_d_ = exp(-(Y_d_ - min(Y_d_,[],'all'))./(2*lambda^2));
Xlim_ = prctile(X_d_,[ 5,95],'all'); Xlim_ = mean(Xlim_) + 1.25*0.5*diff(Xlim_)*[-1,+1];
Ylim_ = prctile(Y_d_,[ 5,95],'all'); Ylim_ = mean(Ylim_) + 1.25*0.5*diff(Ylim_)*[-1,+1];
%Zlim_ = prctile(Z_d_,[ 5,95],'all'); Zlim_ = mean(Zlim_) + 1.25*0.5*diff(Zlim_)*[-1,+1];
Zlim_ = prctile(Z_d_,[ 0,100],'all');
figure(1+nf);nf=nf+1;clf;figmed;
p_row = 1; p_col = 4; np=0;
subplot(p_row,p_col,1+np);np=np+1; imagesc_c(n_d,d_c_0_d_,n_d,d_c_1_d_,X_d_,Xlim_,colormap_beach());
axis image; axisnotick; title('$\chi(A,B)$','Interpreter','latex');% tmp_c_ = colorbar; set(tmp_c_,'Ticks',[],'TickLength',0);
subplot(p_row,p_col,1+np);np=np+1; imagesc_c(n_d,d_c_0_d_,n_d,d_c_1_d_,Y_d_,Ylim_,colormap_80s());
axis image; axisnotick; title('$\|A-B\|^{2}$','Interpreter','latex');% tmp_c_ = colorbar; set(tmp_c_,'Ticks',[],'TickLength',0);
subplot(p_row,p_col,1+np);np=np+1; imagesc_c(n_d,d_c_0_d_,n_d,d_c_1_d_,Z_d_,Zlim_,colormap_81s());
axis image; axisnotick; title(sprintf('$\\exp(-0.5\\times\\|A-B\\|^{2}/\\lambda^{2}): \\lambda=%0.3f$',lambda),'Interpreter','latex');% tmp_c_ = colorbar; set(tmp_c_,'Ticks',[],'TickLength',0);
subplot(p_row,p_col,1+np);np=np+1; contour(d_c_0_d_,d_c_1_d_,transpose(reshape(Z_d_,[n_d,n_d])),linspace(min(Zlim_),max(Zlim_),32));
axis image; axisnotick; title(sprintf('$\\exp(-0.5\\times\\|A-B\\|^{2}/\\lambda^{2}): \\lambda=%0.3f$',lambda),'Interpreter','latex');% tmp_c_ = colorbar; set(tmp_c_,'Ticks',[],'TickLength',0);
sgtitle(sprintf('K %0.2f/(2$\\pi$), d %0.3f, nM %.2d nS %.2d nw %.4d $(\\gamma %0.2f\\pi)$',round(2*pi*k_p_r_max),d_p_r_max,nM,nS,nw,2*pi*nw/n_w_max),'Interpreter','latex');
drawnow();
fname_fig_pre = sprintf('%s_jpg/test_CryoLikelihood_nS%d_nM%d_nw%d_lambda%.3d_d_c',dir_pm,nS,nM,nw,floor(100*lambda));
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if (flag_verbose>0); disp(sprintf(' %% writing %s',fname_fig_jpg)); end;
print('-djpeg',fname_fig_jpg);
close(gcf);
end;%for nw=nw_;
end;%for lambda = [0.10,0.15,0.20,0.25];
%%%%;
end;%for nS=0:n_S-1;
end;%for nM=0:n_M-1;
end;%if flag_disp;
%%%%%%%%;



if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

disp('returning');return;


