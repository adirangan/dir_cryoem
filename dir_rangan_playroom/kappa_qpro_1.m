function ...
[ ...
 parameter ...
,KAPPA ...
] = ...
kappa_qpro_1( ...
 parameter ...
,KAPPA ...
);

str_thisfunction = 'kappa_qpro_1';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if nargin<1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
disp(sprintf(' %% testing %s',str_thisfunction));
parameter = struct('type','parameter');
parameter.flag_verbose = 0; parameter.flag_disp = 0;
parameter.kernel_qpro_l_max_use = 49;
parameter.kernel_qpro_l_max_ext = ceil(1.25*49);
parameter.kernel_qpro_l_max_band = floor(parameter.kernel_qpro_l_max_use/2)-1;
parameter.flag_kernel_full = 0; %<-- if set to 1 then use full brute-force kernel. ;
parameter.kernel_qpro_polar_a_pole_north=5*pi/24;
parameter.kernel_qpro_polar_a_pole_south=3*pi/24;
parameter.kernel_qpro_qref_k_eq_d_double = 0.25;
[~,KAPPA] = kappa_qpro_1(parameter);
%%%%%%%%;
flag_verbose = 1; flag_disp = 1; nf=0;
l_max = KAPPA.l_max_use;
l_max_band = KAPPA.l_max_band;
Rz = KAPPA.Rz;
dRz = KAPPA.dRz;
Ry = KAPPA.Ry;
dRy = KAPPA.dRy;
kappa_norm_ = KAPPA.kappa_norm_;
flag_kernel_full = KAPPA.flag_kernel_full;
qref_k_eq_d = KAPPA.qref_k_eq_d;
n_ring_north = KAPPA.n_ring_north;
n_nearest_north = KAPPA.n_nearest_north;
n_ring_south = KAPPA.n_ring_south;
n_nearest_south = KAPPA.n_nearest_south;
n_nearest_total = KAPPA.n_nearest_total;
qref_n_shell = KAPPA.qref_n_shell;
qref_azimu_b_shell_ = KAPPA.qref_azimu_b_shell_;
qref_polar_a_shell_ = KAPPA.qref_polar_a_shell_;
qref_weight_shell_ = KAPPA.qref_weight_shell_;
qref_k_c_0_shell_ = KAPPA.qref_k_c_0_shell_;
qref_k_c_1_shell_ = KAPPA.qref_k_c_1_shell_;
qref_k_c_2_shell_ = KAPPA.qref_k_c_2_shell_;
qref_n_polar_a = KAPPA.qref_n_polar_a;
qref_polar_a_ = KAPPA.qref_polar_a_;
qref_n_azimu_b_ = KAPPA.qref_n_azimu_b_;
qref_k_c_qc__ = KAPPA.qref_k_c_qc__;
a_keep_node_ = KAPPA.a_keep_node_;
chebfun_kernel_norm_qpro_ = KAPPA.chebfun_kernel_norm_qpro_;
Y_l_val_ = KAPPA.Y_l_val_use_;
Y_m_val_ = KAPPA.Y_m_val_use_;
Ylm_weight_yq__ = KAPPA.Ylm_use_weight_yq__;
if (flag_verbose>0);
  disp(sprintf(' %% KAPPA.qref_k_eq_d %+0.6f',KAPPA.qref_k_eq_d));
  disp(sprintf(' %% KAPPA.n_ring_north %+0.6f',KAPPA.n_ring_north));
  disp(sprintf(' %% KAPPA.n_nearest_north %+0.6f',KAPPA.n_nearest_north));
  disp(sprintf(' %% KAPPA.n_ring_south %+0.6f',KAPPA.n_ring_south));
  disp(sprintf(' %% KAPPA.n_nearest_south %+0.6f',KAPPA.n_nearest_south));
  disp(sprintf(' %% KAPPA.n_nearest_total %+0.6f',KAPPA.n_nearest_total));
  disp(sprintf(' %% KAPPA.qref_n_shell %+0.6f',KAPPA.qref_n_shell));
  disp(sprintf(' %% KAPPA.qref_n_polar_a %+0.6f',KAPPA.qref_n_polar_a));
end;%if (flag_verbose>0);
%%%%;
% build grid on shell. ;
%%%%;
k_p_r_max = 48.0/(2*pi); k_eq_d = 0.5/(2*pi);
n_k_p_r = 1; k_p_r_1 = 1.0; k_p_r_ = k_p_r_1;
[ ...
 n_shell ...
,azimu_b_shell_ ...
,polar_a_shell_ ...
,weight_shell_ ...
,k_c_0_shell_ ...
,k_c_1_shell_ ...
,k_c_2_shell_ ...
] = ...
sample_shell_6( ...
 k_p_r_1 ...
,k_eq_d/k_p_r_max ...
) ;
k_p_r_shell_ = k_p_r_(1+0)*ones(n_shell,1);
%%%%;
% set up l_val and m_val. ;
%%%%;
n_lm = (l_max+1).^2;
l_val_ = transpose(0:l_max);
m_max_ = transpose(-l_max : +l_max);
n_m_max = length(m_max_);
Y_l_val_ = zeros(n_lm,1);
Y_m_val_ = zeros(n_lm,1);
tmp_l_val_ = zeros(n_lm,1);
tmp_m_val_ = zeros(n_lm,1);
na=0; 
for l_val=0:l_max;
for m_val=-l_val:+l_val;
tmp_l_val_(1+na) = l_val;
tmp_m_val_(1+na) = m_val;
na=na+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max;
tmp_index_ = 0:n_lm-1;
Y_l_val_(1+tmp_index_) = tmp_l_val_;
Y_m_val_(1+tmp_index_) = tmp_m_val_;
weight_Y_ = ones(n_lm,1);
weight_3d_k_p_r_ = 4*pi;
%%%%;
% build delta-function on shell. ;
%%%%;
a_deltafunct_band_I_I_k_Y_lm_ = sqrt(2*pi)*sqrt(4*pi)*sqrt(1+2*Y_l_val_).*(Y_m_val_==0).*(Y_l_val_<=l_max_band);
if ~exist('Ylm_uklma___','var'); Ylm_uklma___=[]; end;
if ~exist('k_p_azimu_b_sub_uka__','var'); k_p_azimu_b_sub_uka__ = []; end;
if ~exist('k_p_polar_a_sub_uka__','var'); k_p_polar_a_sub_uka__ = []; end;
if ~exist('l_max_uk_','var'); l_max_uk_ = []; end;
if ~exist('index_nu_n_k_per_shell_from_nk_p_r_','var'); index_nu_n_k_per_shell_from_nk_p_r_ = []; end;
if ~exist('index_k_per_shell_uka__','var'); index_k_per_shell_uka__ = []; end;
[ ...
 a_deltafunct_band_I_I_k_p_quad_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
] = ...
convert_spharm_to_k_p_4( ...
 flag_verbose ...
,n_shell ...
,[0,n_shell] ...
,k_p_r_shell_ ...
,azimu_b_shell_ ...
,polar_a_shell_ ...
,weight_3d_k_p_r_ ...
,weight_shell_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,l_max ...
,a_deltafunct_band_I_I_k_Y_lm_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
%%%%;
% build kernel. ;
%%%%;
kernel_qpro_k_Y_form_ = zeros(n_lm,1);
kernel_qpro_k_Y_form_(1+l_val_.^2 + l_val_) = kappa_norm_(1+l_val_);
deconvolve_lm_ = (sqrt(4*pi)*sqrt(1+2*Y_l_val_).*(Y_l_val_<=l_max_band))./max(1e-12,kappa_norm_(1+Y_l_val_));
if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;
linewidth_use = 3;
plot(Y_l_val_,log(deconvolve_lm_),'k.','LineWidth',linewidth_use);
xlabel('Y_l_val_','Interpreter','none');
ylabel('log(deconvolve_lm_)','Interpreter','none');
grid on;
end;%if flag_disp;
%%%%;
[ ...
 kernel_qpro_k_p_quad_ ...
] = ...
convert_spharm_to_k_p_4( ...
 flag_verbose ...
,n_shell ...
,[0,n_shell] ...
,k_p_r_shell_ ...
,azimu_b_shell_ ...
,polar_a_shell_ ...
,weight_3d_k_p_r_ ...
,weight_shell_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,l_max ...
,kernel_qpro_k_Y_form_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
if (flag_verbose>0); disp(sprintf(' %% sum(kernel_qpro_k_p_quad_.*weight_shell_): %0.6f',sum(kernel_qpro_k_p_quad_.*weight_shell_))); end;
if (flag_verbose>0); disp(sprintf(' %% chebfun_kernel_norm_qpro_(cos(polar_a_shell_)) vs kernel_qpro_k_p_quad_*sqrt(2*pi): %0.16f',fnorm(chebfun_kernel_norm_qpro_(cos(polar_a_shell_))-kernel_qpro_k_p_quad_*sqrt(2*pi))/max(1e-12,fnorm(chebfun_kernel_norm_qpro_(cos(polar_a_shell_)))))); end;
[ ...
 kernel_qpro_k_Y_quad_ ...
] = ...
convert_k_p_to_spharm_4( ...
 flag_verbose ...
,n_shell ...
,[0,n_shell] ...
,k_p_r_shell_ ...
,azimu_b_shell_ ...
,polar_a_shell_ ...
,weight_3d_k_p_r_ ...
,weight_shell_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,l_max ...
,kernel_qpro_k_p_quad_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
if (flag_verbose>0); disp(sprintf(' %% kernel_qpro_k_Y_form_ vs kernel_qpro_k_Y_quad_: %0.16f',fnorm(kernel_qpro_k_Y_form_-kernel_qpro_k_Y_quad_)/fnorm(kernel_qpro_k_Y_form_))); end;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;
klim_ = [0,prctile(chebfun_kernel_norm_qpro_(cos(a_keep_node_)),99.5)];
flag_2d_vs_3d = 0;
imagesc_polar_a_azimu_b_0( ...
 polar_a_shell_ ... 
,azimu_b_shell_ ... 
,real(kernel_qpro_k_p_quad_) ... 
,klim_ ... 
,colormap_81s ... 
,flag_2d_vs_3d ...
,k_p_r_1 ...
);
title('real(kernel_qpro_k_p_quad_)','Interpreter','none');
end;%if flag_disp;
%%%%;
% Now we place a point-source (i.e., delta-function) somewhere on the sphere, ;
% mollify it via the kernel_qpro_, ;
% then apply quadrature to calculate the associated spherical-harmonics, ;
% then compare the result with the rotated kernel. ;
%%%%;
% For this calculation we use a quadrature-reference grid (qref). ;
%%%%;
tmp_I__ = eye(3,3);
tmp_I_pole_ = tmp_I__*[0;0;1];
tmp_I_d2_ = (qref_k_c_0_shell_ - tmp_I_pole_(1+0)).^2 + (qref_k_c_1_shell_ - tmp_I_pole_(1+1)).^2 + (qref_k_c_2_shell_ - tmp_I_pole_(1+2)).^2 ;
tmp_I_arc_ = acos(1-tmp_I_d2_/2);
tmp_I_ker_ = chebfun_kernel_norm_qpro_(cos(tmp_I_arc_));
tmp_ij_north_ = []; if n_nearest_north> 0; [tmp_ij_north_] = knnsearch(qref_k_c_qc__,reshape(+tmp_I_pole_,[1,3]),'K',n_nearest_north); end;
tmp_ij_south_ = []; if n_nearest_south> 0; [tmp_ij_south_] = knnsearch(qref_k_c_qc__,reshape(-tmp_I_pole_,[1,3]),'K',n_nearest_south); end;
tmp_index_keep_ = union(tmp_ij_north_,tmp_ij_south_)-1; tmp_index_drop_ = setdiff(0:qref_n_shell-1,tmp_index_keep_);
tmp_I_ker_crop_ = tmp_I_ker_; tmp_I_ker_crop_(1+tmp_index_drop_)=0;
tmp_euler_a = -3*pi/8;
tmp_euler_b = pi/12;
tmp_euler_c = pi/4;
tmp_R__ = Rz(tmp_euler_b)*Ry(tmp_euler_a)*Rz(tmp_euler_c);
tmp_R_pole_ = tmp_R__*[0;0;1];
tmp_R_d2_ = (qref_k_c_0_shell_ - tmp_R_pole_(1+0)).^2 + (qref_k_c_1_shell_ - tmp_R_pole_(1+1)).^2 + (qref_k_c_2_shell_ - tmp_R_pole_(1+2)).^2 ;
tmp_R_arc_ = acos(1-tmp_R_d2_/2);
tmp_R_ker_ = chebfun_kernel_norm_qpro_(cos(tmp_R_arc_));
disp(sprintf(' %% tmp_R_ker_ arc_ vs 1-d2_/2: %0.16f',fnorm(tmp_R_ker_ - chebfun_kernel_norm_qpro_(1-tmp_R_d2_/2))/fnorm(tmp_R_ker_)));
tmp_ij_north_ = []; if n_nearest_north> 0; [tmp_ij_north_] = knnsearch(qref_k_c_qc__,reshape(+tmp_R_pole_,[1,3]),'K',n_nearest_north); end;
tmp_ij_south_ = []; if n_nearest_south> 0; [tmp_ij_south_] = knnsearch(qref_k_c_qc__,reshape(-tmp_R_pole_,[1,3]),'K',n_nearest_south); end;
tmp_index_keep_ = union(tmp_ij_north_,tmp_ij_south_)-1; tmp_index_drop_ = setdiff(0:qref_n_shell-1,tmp_index_keep_);
tmp_R_ker_crop_ = tmp_R_ker_; tmp_R_ker_crop_(1+tmp_index_drop_)=0;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 2; p_col = 3; np=0;
subplot(p_row,p_col,1+np);np=np+1;
alim_ = [0,prctile(tmp_R_arc_,99.5)];
imagesc_polar_a_azimu_b_0( ...
 qref_polar_a_shell_ ... 
,qref_azimu_b_shell_ ... 
,real(tmp_R_arc_) ... 
,alim_ ... 
,colormap_81s ... 
,flag_2d_vs_3d ...
,k_p_r_1 ...
);
hold on;
plot3(tmp_R_pole_(1+0),tmp_R_pole_(1+1),tmp_R_pole_(1+2),'wo','MarkerFaceColor','g');
hold off;
title('real(tmp_R_arc_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
klim_ = prctile(chebfun_kernel_norm_qpro_(cos(a_keep_node_)),75)*[-1,+1];
imagesc_polar_a_azimu_b_0( ...
 qref_polar_a_shell_ ... 
,qref_azimu_b_shell_ ... 
,real(tmp_I_ker_) ... 
,klim_ ... 
,colormap_80s ... 
,flag_2d_vs_3d ...
,k_p_r_1 ...
);
hold on;
plot3(tmp_I_pole_(1+0),tmp_I_pole_(1+1),tmp_I_pole_(1+2),'wo','MarkerFaceColor','g');
hold off;
title('real(tmp_I_ker_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
klim_ = prctile(chebfun_kernel_norm_qpro_(cos(a_keep_node_)),75)*[-1,+1];
imagesc_polar_a_azimu_b_0( ...
 qref_polar_a_shell_ ... 
,qref_azimu_b_shell_ ... 
,real(tmp_R_ker_) ...
,klim_ ... 
,colormap_80s ... 
,flag_2d_vs_3d ...
,k_p_r_1 ...
);
hold on;
plot3(tmp_R_pole_(1+0),tmp_R_pole_(1+1),tmp_R_pole_(1+2),'wo','MarkerFaceColor','g');
hold off;
title('real(tmp_R_ker_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
klim_ = prctile(chebfun_kernel_norm_qpro_(cos(a_keep_node_)),75)*[-1,+1];
imagesc_polar_a_azimu_b_0( ...
 qref_polar_a_shell_ ... 
,qref_azimu_b_shell_ ... 
,real(tmp_I_ker_crop_) ... 
,klim_ ... 
,colormap_80s ... 
,flag_2d_vs_3d ...
,k_p_r_1 ...
);
hold on;
plot3(tmp_I_pole_(1+0),tmp_I_pole_(1+1),tmp_I_pole_(1+2),'wo','MarkerFaceColor','g');
hold off;
title('real(tmp_I_ker_crop_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
klim_ = prctile(chebfun_kernel_norm_qpro_(cos(a_keep_node_)),75)*[-1,+1];
imagesc_polar_a_azimu_b_0( ...
 qref_polar_a_shell_ ... 
,qref_azimu_b_shell_ ... 
,real(tmp_R_ker_crop_) ...
,klim_ ... 
,colormap_80s ... 
,flag_2d_vs_3d ...
,k_p_r_1 ...
);
hold on;
plot3(tmp_R_pole_(1+0),tmp_R_pole_(1+1),tmp_R_pole_(1+2),'wo','MarkerFaceColor','g');
hold off;
title('real(tmp_R_ker_crop_)','Interpreter','none');
end;%if flag_disp;
%%%%;
a_I_I_k_Y_lm_ = conj(Ylm_weight_yq__)*tmp_I_ker_;
a_I_R_k_Y_lm_ = conj(Ylm_weight_yq__)*tmp_R_ker_;
tmp_euler_ = [tmp_euler_b,tmp_euler_a,tmp_euler_c];
tmp_euler_pos_ = [tmp_euler_c,tmp_euler_a,tmp_euler_b]; tmp_euler_neg_ = -flip(tmp_euler_pos_);
a_R_I_k_Y_lm_ = rotate_spharm_to_spharm_3(1,1,l_max,a_I_I_k_Y_lm_,tmp_euler_pos_);
a_R_R_k_Y_lm_ = rotate_spharm_to_spharm_3(1,1,l_max,a_I_R_k_Y_lm_,tmp_euler_neg_);
if (flag_verbose>0); disp(sprintf(' %% a_R_I_k_Y_lm_ vs a_I_R_k_Y_lm_: %0.16f',fnorm(a_R_I_k_Y_lm_ - a_I_R_k_Y_lm_)/fnorm(a_R_I_k_Y_lm_))); end;
if (flag_verbose>0); disp(sprintf(' %% a_I_I_k_Y_lm_ vs a_R_R_k_Y_lm_: %0.16f',fnorm(a_I_I_k_Y_lm_ - a_R_R_k_Y_lm_)/fnorm(a_I_I_k_Y_lm_))); end;
a_I_I_crop_k_Y_lm_ = conj(Ylm_weight_yq__)*tmp_I_ker_crop_;
a_I_R_crop_k_Y_lm_ = conj(Ylm_weight_yq__)*tmp_R_ker_crop_;
tmp_euler_ = [tmp_euler_b,tmp_euler_a,tmp_euler_c];
tmp_euler_pos_ = [tmp_euler_c,tmp_euler_a,tmp_euler_b]; tmp_euler_neg_ = -flip(tmp_euler_pos_);
a_R_I_crop_k_Y_lm_ = rotate_spharm_to_spharm_3(1,1,l_max,a_I_I_crop_k_Y_lm_,tmp_euler_pos_);
a_R_R_crop_k_Y_lm_ = rotate_spharm_to_spharm_3(1,1,l_max,a_I_R_crop_k_Y_lm_,tmp_euler_neg_);
if (flag_verbose>0); disp(sprintf(' %% a_I_I_crop_k_Y_lm_ vs a_I_I_k_Y_lm_: %0.16f',fnorm(a_I_I_crop_k_Y_lm_ - a_I_I_k_Y_lm_)/fnorm(a_I_I_k_Y_lm_))); end;
if (flag_verbose>0); disp(sprintf(' %% a_R_I_crop_k_Y_lm_ vs a_R_I_k_Y_lm_: %0.16f',fnorm(a_R_I_crop_k_Y_lm_ - a_R_I_k_Y_lm_)/fnorm(a_R_I_k_Y_lm_))); end;
if (flag_verbose>0); disp(sprintf(' %% a_I_R_crop_k_Y_lm_ vs a_I_R_k_Y_lm_: %0.16f',fnorm(a_I_R_crop_k_Y_lm_ - a_I_R_k_Y_lm_)/fnorm(a_I_R_k_Y_lm_))); end;
if (flag_verbose>0); disp(sprintf(' %% a_R_R_crop_k_Y_lm_ vs a_R_R_k_Y_lm_: %0.16f',fnorm(a_R_R_crop_k_Y_lm_ - a_R_R_k_Y_lm_)/fnorm(a_R_R_k_Y_lm_))); end;
if (flag_verbose>0); disp(sprintf(' %% a_R_I_crop_k_Y_lm_ vs a_I_R_crop_k_Y_lm_: %0.16f',fnorm(a_R_I_crop_k_Y_lm_ - a_I_R_crop_k_Y_lm_)/fnorm(a_R_I_crop_k_Y_lm_))); end;
if (flag_verbose>0); disp(sprintf(' %% a_I_I_crop_k_Y_lm_ vs a_R_R_crop_k_Y_lm_: %0.16f',fnorm(a_I_I_crop_k_Y_lm_ - a_R_R_crop_k_Y_lm_)/fnorm(a_I_I_crop_k_Y_lm_))); end;
%%%%;
% Now we continue the above experiment: ;
% place a point-source (i.e., delta-function) somewhere on the sphere, ;
% mollify it via the kernel_qpro_, ;
% then apply quadrature to calculate the associated spherical-harmonics, ;
% then deconvolve by the scaling associated with kernel_qpro_, ;
% then compare the result with the rotated kernel. ;
%%%%;
tmp_euler_ = [tmp_euler_b,tmp_euler_a,tmp_euler_c];
tmp_euler_pos_ = [tmp_euler_c,tmp_euler_a,tmp_euler_b]; tmp_euler_neg_ = -flip(tmp_euler_pos_);
a_deltafunct_band_R_I_k_Y_lm_ = rotate_spharm_to_spharm_3(1,1,l_max,a_deltafunct_band_I_I_k_Y_lm_,tmp_euler_pos_);
a_deltafunct_band_I_R_k_Y_lm_ = a_deltafunct_band_R_I_k_Y_lm_; %<-- use formula. ;
a_deltafunct_band_R_R_k_Y_lm_ = rotate_spharm_to_spharm_3(1,1,l_max,a_deltafunct_band_I_R_k_Y_lm_,tmp_euler_neg_);
if (flag_verbose>0); disp(sprintf(' %% a_deltafunct_band_I_I_k_Y_lm_ vs a_deltafunct_band_R_R_k_Y_lm_: %0.16f',fnorm(a_deltafunct_band_I_I_k_Y_lm_ - a_deltafunct_band_R_R_k_Y_lm_)/fnorm(a_deltafunct_band_I_I_k_Y_lm_))); end;
a_deconvolve_I_I_k_Y_lm_ = a_I_I_k_Y_lm_.*deconvolve_lm_;
a_deconvolve_I_R_k_Y_lm_ = a_I_R_k_Y_lm_.*deconvolve_lm_;
a_deconvolve_R_I_k_Y_lm_ = a_R_I_k_Y_lm_.*deconvolve_lm_;
a_deconvolve_R_R_k_Y_lm_ = a_R_R_k_Y_lm_.*deconvolve_lm_;
if (flag_verbose>0); disp(sprintf(' %% a_deltafunct_band_I_I_k_Y_lm_ vs a_deconvolve_I_I_k_Y_lm_: %0.16f',fnorm(a_deltafunct_band_I_I_k_Y_lm_ - a_deconvolve_I_I_k_Y_lm_)/fnorm(a_deltafunct_band_I_I_k_Y_lm_))); end;
if (flag_verbose>0); disp(sprintf(' %% a_deltafunct_band_I_R_k_Y_lm_ vs a_deconvolve_I_R_k_Y_lm_: %0.16f',fnorm(a_deltafunct_band_I_R_k_Y_lm_ - a_deconvolve_I_R_k_Y_lm_)/fnorm(a_deltafunct_band_I_R_k_Y_lm_))); end;
if (flag_verbose>0); disp(sprintf(' %% a_deltafunct_band_R_I_k_Y_lm_ vs a_deconvolve_R_I_k_Y_lm_: %0.16f',fnorm(a_deltafunct_band_R_I_k_Y_lm_ - a_deconvolve_R_I_k_Y_lm_)/fnorm(a_deltafunct_band_R_I_k_Y_lm_))); end;
if (flag_verbose>0); disp(sprintf(' %% a_deltafunct_band_R_R_k_Y_lm_ vs a_deconvolve_R_R_k_Y_lm_: %0.16f',fnorm(a_deltafunct_band_R_R_k_Y_lm_ - a_deconvolve_R_R_k_Y_lm_)/fnorm(a_deltafunct_band_R_R_k_Y_lm_))); end;
a_deconvolve_I_I_crop_k_Y_lm_ = a_I_I_crop_k_Y_lm_.*deconvolve_lm_;
a_deconvolve_I_R_crop_k_Y_lm_ = a_I_R_crop_k_Y_lm_.*deconvolve_lm_;
a_deconvolve_R_I_crop_k_Y_lm_ = a_R_I_crop_k_Y_lm_.*deconvolve_lm_;
a_deconvolve_R_R_crop_k_Y_lm_ = a_R_R_crop_k_Y_lm_.*deconvolve_lm_;
if (flag_verbose>0); disp(sprintf(' %% a_deltafunct_band_I_I_k_Y_lm_ vs a_deconvolve_I_I_crop_k_Y_lm_: %0.16f',fnorm(a_deltafunct_band_I_I_k_Y_lm_ - a_deconvolve_I_I_crop_k_Y_lm_)/fnorm(a_deltafunct_band_I_I_k_Y_lm_))); end;
if (flag_verbose>0); disp(sprintf(' %% a_deltafunct_band_I_R_k_Y_lm_ vs a_deconvolve_I_R_crop_k_Y_lm_: %0.16f',fnorm(a_deltafunct_band_I_R_k_Y_lm_ - a_deconvolve_I_R_crop_k_Y_lm_)/fnorm(a_deltafunct_band_I_R_k_Y_lm_))); end;
if (flag_verbose>0); disp(sprintf(' %% a_deltafunct_band_R_I_k_Y_lm_ vs a_deconvolve_R_I_crop_k_Y_lm_: %0.16f',fnorm(a_deltafunct_band_R_I_k_Y_lm_ - a_deconvolve_R_I_crop_k_Y_lm_)/fnorm(a_deltafunct_band_R_I_k_Y_lm_))); end;
if (flag_verbose>0); disp(sprintf(' %% a_deltafunct_band_R_R_k_Y_lm_ vs a_deconvolve_R_R_crop_k_Y_lm_: %0.16f',fnorm(a_deltafunct_band_R_R_k_Y_lm_ - a_deconvolve_R_R_crop_k_Y_lm_)/fnorm(a_deltafunct_band_R_R_k_Y_lm_))); end;
%%%%%%%%;
[ ...
 a_deconvolve_I_I_crop_k_p_quad_ ...
] = ...
convert_spharm_to_k_p_4( ...
 flag_verbose ...
,n_shell ...
,[0,n_shell] ...
,k_p_r_shell_ ...
,azimu_b_shell_ ...
,polar_a_shell_ ...
,weight_3d_k_p_r_ ...
,weight_shell_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,l_max ...
,a_deconvolve_I_I_crop_k_Y_lm_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
if (flag_disp);
figure(1+nf);nf=nf+1;clf;figmed;
p_row = 1; p_col = 3; np=0;
deltafunct_band_lim_ = [0,prctile(a_deltafunct_band_I_I_k_p_quad_,99.5)];
subplot(p_row,p_col,1+np);np=np+1;
flag_2d_vs_3d = 0;
k_p_r_max = 1;
imagesc_polar_a_azimu_b_0( ...
 polar_a_shell_ ... 
,azimu_b_shell_ ... 
,a_deltafunct_band_I_I_k_p_quad_ ... 
,deltafunct_band_lim_ ... 
,colormap_81s ... 
,flag_2d_vs_3d ...
,k_p_r_max ...
);
axisnotick3d; axis vis3d;
title('deltafunct_band_I_I_k_p_form_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
flag_2d_vs_3d = 0;
k_p_r_max = 1;
imagesc_polar_a_azimu_b_0( ...
 polar_a_shell_ ... 
,azimu_b_shell_ ... 
,real(a_deconvolve_I_I_crop_k_p_quad_) ... 
,deltafunct_band_lim_ ... 
,colormap_81s ... 
,flag_2d_vs_3d ...
,k_p_r_max ...
);
axisnotick3d; axis vis3d;
title('deconvolve_I_I_crop_k_p_form_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
plot(Y_l_val_,log10(abs(a_deltafunct_band_I_R_k_Y_lm_-a_deconvolve_I_R_crop_k_Y_lm_)./abs(a_deltafunct_band_I_R_k_Y_lm_)),'.');
xlabel('Y_l_val_','Interpreter','none'); xlim([0,l_max+1]);
ylabel('log10(abs(delta-decon)./abs(delta)))','Interpreter','none');
grid on;
end;%if (flag_disp);
%%%%;
if (flag_verbose>0); disp(sprintf(' %% testing first-derivative')); end;
%%%%;
tmp_I__ = eye(3,3);
tmp_I_pole_ = tmp_I__*[0;0;1];
tmp_I_d2_ = (qref_k_c_0_shell_ - tmp_I_pole_(1+0)).^2 + (qref_k_c_1_shell_ - tmp_I_pole_(1+1)).^2 + (qref_k_c_2_shell_ - tmp_I_pole_(1+2)).^2 ;
tmp_I_arc_ = acos(1-tmp_I_d2_/2);
tmp_I_ker_ = chebfun_kernel_norm_qpro_(cos(tmp_I_arc_));
tmp_ij_north_ = []; if n_nearest_north> 0; [tmp_ij_north_] = knnsearch(qref_k_c_qc__,reshape(+tmp_I_pole_,[1,3]),'K',n_nearest_north); end;
tmp_ij_south_ = []; if n_nearest_south> 0; [tmp_ij_south_] = knnsearch(qref_k_c_qc__,reshape(-tmp_I_pole_,[1,3]),'K',n_nearest_south); end;
tmp_index_keep_ = union(tmp_ij_north_,tmp_ij_south_)-1; tmp_index_drop_ = setdiff(0:qref_n_shell-1,tmp_index_keep_);
tmp_I_ker_crop_ = tmp_I_ker_; tmp_I_ker_crop_(1+tmp_index_drop_)=0;
%%%%;
tmp_euler_a = -3*pi/8;
tmp_euler_b = pi/12;
tmp_euler_c = pi/4;
tmp_R__ = Rz(tmp_euler_b)*Ry(tmp_euler_a)*Rz(tmp_euler_c);
tmp_R_pole_ = tmp_R__*[0;0;1];
tmp_R_d2_ = (qref_k_c_0_shell_ - tmp_R_pole_(1+0)).^2 + (qref_k_c_1_shell_ - tmp_R_pole_(1+1)).^2 + (qref_k_c_2_shell_ - tmp_R_pole_(1+2)).^2 ;
tmp_R_arc_ = acos(1-tmp_R_d2_/2);
tmp_R_ker_ = chebfun_kernel_norm_qpro_(cos(tmp_R_arc_));
tmp_ij_north_ = []; if n_nearest_north> 0; [tmp_ij_north_] = knnsearch(qref_k_c_qc__,reshape(+tmp_R_pole_,[1,3]),'K',n_nearest_north); end;
tmp_ij_south_ = []; if n_nearest_south> 0; [tmp_ij_south_] = knnsearch(qref_k_c_qc__,reshape(-tmp_R_pole_,[1,3]),'K',n_nearest_south); end;
tmp_index_keep_ = union(tmp_ij_north_,tmp_ij_south_)-1; tmp_index_drop_ = setdiff(0:qref_n_shell-1,tmp_index_keep_);
tmp_R_ker_crop_ = tmp_R_ker_; tmp_R_ker_crop_(1+tmp_index_drop_)=0;
%%%%;
da = +0.65*1e-4;
db = -1.00*1e-4;
dc = +0.25*1e-4;
tmp_dRda__ = Rz(tmp_euler_b)*dRy(tmp_euler_a)*Rz(tmp_euler_c);
tmp_dRdb__ = dRz(tmp_euler_b)*Ry(tmp_euler_a)*Rz(tmp_euler_c);
tmp_dRdc__ = Rz(tmp_euler_b)*Ry(tmp_euler_a)*dRz(tmp_euler_c);
tmp_dR_mid__ = (tmp_dRda__*da + tmp_dRdb__*db + tmp_dRdc__*dc)/max(1e-12,fnorm([da;db;dc]));
tmp_R_pos__ = Rz(tmp_euler_b+db)*Ry(tmp_euler_a+da)*Rz(tmp_euler_c+dc);
tmp_R_neg__ = Rz(tmp_euler_b-db)*Ry(tmp_euler_a-da)*Rz(tmp_euler_c-dc);
tmp_dR_dif__ = (tmp_R_pos__ - tmp_R_neg__)/max(1e-12,2*fnorm([da;db;dc]));;
disp(sprintf(' %% tmp_dR_dif__ vs tmp_dR_mid__: %0.16f',fnorm(tmp_dR_dif__ - tmp_dR_mid__)/fnorm(tmp_dR_dif__)));
tmp_R_pos_pole_ = tmp_R_pos__*[0;0;1];
tmp_R_neg_pole_ = tmp_R_neg__*[0;0;1];
tmp_dRda_pole_ = tmp_dRda__*[0;0;1];
tmp_dRdb_pole_ = tmp_dRdb__*[0;0;1];
tmp_dRdc_pole_ = tmp_dRdc__*[0;0;1];
tmp_dR_mid_d2_ = ...
 - 2*(qref_k_c_0_shell_ - tmp_R_pole_(1+0)).*(tmp_dRda_pole_(1+0)*da + tmp_dRdb_pole_(1+0)*db + tmp_dRdc_pole_(1+0)*dc) ...
 - 2*(qref_k_c_1_shell_ - tmp_R_pole_(1+1)).*(tmp_dRda_pole_(1+1)*da + tmp_dRdb_pole_(1+1)*db + tmp_dRdc_pole_(1+1)*dc) ...
 - 2*(qref_k_c_2_shell_ - tmp_R_pole_(1+2)).*(tmp_dRda_pole_(1+2)*da + tmp_dRdb_pole_(1+2)*db + tmp_dRdc_pole_(1+2)*dc) ...
;
tmp_R_pos_d2_ = (qref_k_c_0_shell_ - tmp_R_pos_pole_(1+0)).^2 + (qref_k_c_1_shell_ - tmp_R_pos_pole_(1+1)).^2 + (qref_k_c_2_shell_ - tmp_R_pos_pole_(1+2)).^2 ;
tmp_R_neg_d2_ = (qref_k_c_0_shell_ - tmp_R_neg_pole_(1+0)).^2 + (qref_k_c_1_shell_ - tmp_R_neg_pole_(1+1)).^2 + (qref_k_c_2_shell_ - tmp_R_neg_pole_(1+2)).^2 ;
tmp_dR_dif_d2_ = (tmp_R_pos_d2_ - tmp_R_neg_d2_)/2;
disp(sprintf(' %% tmp_dR_dif_d2_ vs tmp_dR_mid_d2_: %0.16f',fnorm(tmp_dR_dif_d2_ - tmp_dR_mid_d2_)/fnorm(tmp_dR_dif_d2_)));
dchebfun_kernel_norm_qpro_ = diff(chebfun_kernel_norm_qpro_);
tmp_dR_mid_ker_ = -dchebfun_kernel_norm_qpro_(1-tmp_R_d2_/2).*tmp_dR_mid_d2_/2;
tmp_R_pos_ker_ = chebfun_kernel_norm_qpro_(1-tmp_R_pos_d2_/2);
tmp_R_neg_ker_ = chebfun_kernel_norm_qpro_(1-tmp_R_neg_d2_/2);
tmp_dR_dif_ker_ = (tmp_R_pos_ker_ - tmp_R_neg_ker_)/2;
disp(sprintf(' %% tmp_dR_dif_ker_ vs tmp_dR_mid_ker_: %0.16f',fnorm(tmp_dR_dif_ker_ - tmp_dR_mid_ker_)/fnorm(tmp_dR_dif_ker_)));
tmp_dR_mid_ker_crop_ = tmp_dR_mid_ker_; tmp_dR_mid_ker_crop_(1+tmp_index_drop_)=0;
tmp_R_pos_ker_crop_ = tmp_R_pos_ker_; tmp_R_pos_ker_crop_(1+tmp_index_drop_)=0;
tmp_R_neg_ker_crop_ = tmp_R_neg_ker_; tmp_R_neg_ker_crop_(1+tmp_index_drop_)=0;
%%%%;
a_I_dR_mid_k_Y_lm_ = conj(Ylm_weight_yq__)*tmp_dR_mid_ker_;
a_I_R_pos_k_Y_lm_ = conj(Ylm_weight_yq__)*tmp_R_pos_ker_;
a_I_R_neg_k_Y_lm_ = conj(Ylm_weight_yq__)*tmp_R_neg_ker_;
a_I_dR_dif_k_Y_lm_ = (a_I_R_pos_k_Y_lm_ - a_I_R_neg_k_Y_lm_)/2;
disp(sprintf(' %% a_I_dR_dif_k_Y_lm_ vs a_I_dR_mid_k_Y_lm_: %0.16f',fnorm(a_I_dR_dif_k_Y_lm_ - a_I_dR_mid_k_Y_lm_)/fnorm(a_I_dR_dif_k_Y_lm_)));
%%%%;
a_I_dR_mid_crop_k_Y_lm_ = conj(Ylm_weight_yq__)*tmp_dR_mid_ker_crop_;
a_I_R_pos_crop_k_Y_lm_ = conj(Ylm_weight_yq__)*tmp_R_pos_ker_crop_;
a_I_R_neg_crop_k_Y_lm_ = conj(Ylm_weight_yq__)*tmp_R_neg_ker_crop_;
a_I_dR_dif_crop_k_Y_lm_ = (a_I_R_pos_crop_k_Y_lm_ - a_I_R_neg_crop_k_Y_lm_)/2;
disp(sprintf(' %% a_I_dR_dif_crop_k_Y_lm_ vs a_I_dR_mid_crop_k_Y_lm_: %0.16f',fnorm(a_I_dR_dif_crop_k_Y_lm_ - a_I_dR_mid_crop_k_Y_lm_)/fnorm(a_I_dR_dif_crop_k_Y_lm_)));
%%%%;
a_deconvolve_I_dR_mid_k_Y_lm_ = a_I_dR_mid_k_Y_lm_.*deconvolve_lm_;
a_deconvolve_I_dR_mid_crop_k_Y_lm_ = a_I_dR_mid_crop_k_Y_lm_.*deconvolve_lm_;
a_deconvolve_I_R_pos_k_Y_lm_ = a_I_R_pos_k_Y_lm_.*deconvolve_lm_;
a_deconvolve_I_R_pos_crop_k_Y_lm_ = a_I_R_pos_crop_k_Y_lm_.*deconvolve_lm_;
a_deconvolve_I_R_neg_k_Y_lm_ = a_I_R_neg_k_Y_lm_.*deconvolve_lm_;
a_deconvolve_I_R_neg_crop_k_Y_lm_ = a_I_R_neg_crop_k_Y_lm_.*deconvolve_lm_;
a_deconvolve_I_dR_dif_k_Y_lm_ = (a_deconvolve_I_R_pos_k_Y_lm_ - a_deconvolve_I_R_neg_k_Y_lm_)/2;
a_deconvolve_I_dR_dif_crop_k_Y_lm_ = (a_deconvolve_I_R_pos_crop_k_Y_lm_ - a_deconvolve_I_R_neg_crop_k_Y_lm_)/2;
disp(sprintf(' %% a_deconvolve_I_dR_dif_k_Y_lm_ vs a_deconvolve_I_dR_mid_k_Y_lm_: %0.16f',fnorm(a_deconvolve_I_dR_dif_k_Y_lm_ - a_deconvolve_I_dR_mid_k_Y_lm_)/fnorm(a_deconvolve_I_dR_dif_k_Y_lm_)));
disp(sprintf(' %% a_deconvolve_I_dR_dif_crop_k_Y_lm_ vs a_deconvolve_I_dR_mid_crop_k_Y_lm_: %0.16f',fnorm(a_deconvolve_I_dR_dif_crop_k_Y_lm_ - a_deconvolve_I_dR_mid_crop_k_Y_lm_)/fnorm(a_deconvolve_I_dR_dif_crop_k_Y_lm_)));
%%%%;
% Now testing with image-data. ;
%%%%;
n_M = 4;
n_w = 98;
rng(0);
viewing_polar_a_M_ = 2*pi*rand(n_M,1);
viewing_azimu_b_M_ = 2*pi*rand(n_M,1);
viewing_gamma_z_M_ = 2*pi*rand(n_M,1);
[ ...
 k_p_polar_a_wM__ ...
,k_p_azimu_b_wM__ ...
,k_c_0_wM__ ...
,k_c_1_wM__ ...
,k_c_2_wM__ ...
,k_p_r01_wM__ ...
,dtau_k_p_polar_a_wM3___ ...
,dtau_k_p_azimu_b_wM3___ ...
,dtau_k_c_0_wM3___ ...
,dtau_k_c_1_wM3___ ...
,dtau_k_c_2_wM3___ ...
,dtau_k_p_r01_wM3___ ...
,dtau_dtau_k_p_polar_a_wM33____ ...
,dtau_dtau_k_p_azimu_b_wM33____ ...
,dtau_dtau_k_c_0_wM33____ ...
,dtau_dtau_k_c_1_wM33____ ...
,dtau_dtau_k_c_2_wM33____ ...
,dtau_dtau_k_p_r01_wM33____ ...
] = ...
cg_rhs_2( ...
 n_M ...
,n_w ...
,viewing_polar_a_M_ ...
,viewing_azimu_b_M_ ...
,viewing_gamma_z_M_ ...
);
data_k_c_wMd__ = [ k_c_0_wM__(:) , k_c_1_wM__(:) , k_c_2_wM__(:) ];
%%%%;
tmp_error = 0;
for nM=0:n_M-1; for nw=0:n_w-1;
tmp_euler_a = +viewing_polar_a_M_(1+nM); tmp_euler_b = +viewing_azimu_b_M_(1+nM); tmp_euler_c = -viewing_gamma_z_M_(1+nM) + (2*pi*nw)/n_w ;
tmp_R__ = Rz(tmp_euler_b)*Ry(tmp_euler_a)*Rz(tmp_euler_c);
tmp_R_point_k_c_ = tmp_R__*[1;0;0];
tmp_diff_ = data_k_c_wMd__(1+nw+nM*n_w,:) - reshape(tmp_R_point_k_c_,[1,3]);
tmp_error = tmp_error + fnorm(tmp_diff_);
end;end;%for nw=0:n_w-1; for nM=0:n_M-1;
if (flag_verbose>0); disp(sprintf(' %% data_k_c_wMd__ error: %0.16f',tmp_error)); end;
%%%%%%%%;
tmp_ij_north_wMn__ = zeros(n_w*n_M,n_nearest_north);
tmp_d1_north_wMn__ = zeros(n_w*n_M,n_nearest_north);
if n_nearest_north> 0;
[tmp_ij_north_wMn__,tmp_d1_north_wMn__] = knnsearch(qref_k_c_qc__,+data_k_c_wMd__,'K',n_nearest_north); %<-- obtain tmp_d1_north_wMn__. ;
end;%if n_nearest_north> 0;
tmp_d2_north_wMn__ = sum(bsxfun(@minus,reshape(qref_k_c_qc__(tmp_ij_north_wMn__(:),:),[n_w*n_M,n_nearest_north,3]),reshape(+data_k_c_wMd__,[n_w*n_M,1,3])).^2,3);
tmp_ij_south_wMn__ = zeros(n_w*n_M,n_nearest_south);
if n_nearest_south> 0;
[tmp_ij_south_wMn__] = knnsearch(qref_k_c_qc__,-data_k_c_wMd__,'K',n_nearest_south); %<-- cannot obtain tmp_d1_south_wMn__. ;
end;%if n_nearest_south> 0;
tmp_d1_south_wMn__ = sqrt(sum(bsxfun(@minus,reshape(qref_k_c_qc__(tmp_ij_south_wMn__(:),:),[n_w*n_M,n_nearest_south,3]),reshape(+data_k_c_wMd__,[n_w*n_M,1,3])).^2,3));
tmp_d2_south_wMn__ = sum(bsxfun(@minus,reshape(qref_k_c_qc__(tmp_ij_south_wMn__(:),:),[n_w*n_M,n_nearest_south,3]),reshape(+data_k_c_wMd__,[n_w*n_M,1,3])).^2,3);
index_keep_wMn__ = [ tmp_ij_north_wMn__ , tmp_ij_south_wMn__ ] - 1;
index_qref_from_data_wMn__ = index_keep_wMn__;
d1_keep_wMn__ = [ tmp_d1_north_wMn__ , tmp_d1_south_wMn__ ] ;
d2_keep_wMn__ = [ tmp_d2_north_wMn__ , tmp_d2_south_wMn__ ] ;
mollify_qref_from_data_wMn__ = chebfun_kernel_norm_qpro_(1-d2_keep_wMn__/2);
index_data_wMn__ = repmat(transpose(0:n_w*n_M-1),[1,n_nearest_total]);
qref_from_data_qwM__ = sparse(1+index_qref_from_data_wMn__,1+index_data_wMn__(:),mollify_qref_from_data_wMn__,qref_n_shell,n_w*n_M);
%%%%%%%%;
rng(0);
T_k_p_wM__ = randn(n_w,n_M); T_k_p_wM_ = reshape(T_k_p_wM__,[n_w*n_M,1]);
T_k_p_q_ = qref_from_data_qwM__*T_k_p_wM_;
T_k_Y_lm_ = conj(Ylm_weight_yq__)*T_k_p_q_;
T_restore_k_Y_lm_ = T_k_Y_lm_.*deconvolve_lm_;
%%%%;
R_k_p_q_ = zeros(qref_n_shell,1);
for nM=0:n_M-1; for nw=0:n_w-1;
tmp_euler_a = +viewing_polar_a_M_(1+nM); tmp_euler_b = +viewing_azimu_b_M_(1+nM); tmp_euler_c = -viewing_gamma_z_M_(1+nM) + (2*pi*nw)/n_w ;
tmp_R__ = Rz(tmp_euler_b)*Ry(tmp_euler_a)*Rz(tmp_euler_c); tmp_R_point_k_c_ = tmp_R__*[1;0;0];
tmp_R_d2_ = (qref_k_c_0_shell_ - tmp_R_point_k_c_(1+0)).^2 + (qref_k_c_1_shell_ - tmp_R_point_k_c_(1+1)).^2 + (qref_k_c_2_shell_ - tmp_R_point_k_c_(1+2)).^2 ;
R_k_p_q_ = R_k_p_q_ + T_k_p_wM_(1+nw+nM*n_w)*chebfun_kernel_norm_qpro_(1-tmp_R_d2_/2);
end;end;%for nw=0:n_w-1; for nM=0:n_M-1;
R_k_Y_lm_ = conj(Ylm_weight_yq__)*R_k_p_q_;
R_restore_k_Y_lm_ = R_k_Y_lm_.*deconvolve_lm_;
if (flag_verbose>0); disp(sprintf(' %% T_k_p_q_ vs R_k_p_q_: %0.16f',fnorm(T_k_p_q_-R_k_p_q_)/fnorm(T_k_p_q_))); end;
if (flag_verbose>0); disp(sprintf(' %% T_k_Y_lm_ vs R_k_Y_lm_: %0.16f',fnorm(T_k_Y_lm_-R_k_Y_lm_)/fnorm(T_k_Y_lm_))); end;
if (flag_verbose>0); disp(sprintf(' %% T_restore_k_Y_lm_ vs R_restore_k_Y_lm_: %0.16f',fnorm(T_restore_k_Y_lm_-R_restore_k_Y_lm_)/fnorm(T_restore_k_Y_lm_))); end;
%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,2,1);
plot(Y_l_val_,log10(abs(T_k_Y_lm_-R_k_Y_lm_)),'.');
xlim([0,l_max]); xlabel('Y_l_val_','Interpreter','none');
ylabel('log10(abs(T_k_Y_lm_-R_k_Y_lm_))','Interpreter','none');
title('log10(abs(T_crop - R_full))','Interpreter','none');
subplot(1,2,2);
plot(Y_l_val_,log10(abs(T_restore_k_Y_lm_-R_restore_k_Y_lm_)),'.');
xlim([0,l_max]); xlabel('Y_l_val_','Interpreter','none');
ylabel('log10(abs(T_restore_k_Y_lm_-R_restore_k_Y_lm_))','Interpreter','none');
title('log10(abs(T_restore_crop - R_restore_full))','Interpreter','none');
end;%if flag_disp;
%%%%;
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% testing first-derivative (directional derivative only)')); end;
%%%%%%%%;
n_M = 4;
n_w = 98;
rng(0);
viewing_polar_a_M_ = 2*pi*rand(n_M,1);
viewing_azimu_b_M_ = 2*pi*rand(n_M,1);
viewing_gamma_z_M_ = 2*pi*rand(n_M,1);
rng(1);
dtau_viewing_polar_a_M_ = 2*pi*rand(n_M,1);
dtau_viewing_azimu_b_M_ = 2*pi*rand(n_M,1);
dtau_viewing_gamma_z_M_ = 2*pi*rand(n_M,1);
%dtau_fnorm = fnorm([dtau_viewing_polar_a_M_,dtau_viewing_azimu_b_M_,dtau_viewing_gamma_z_M_]);
%dtau_viewing_polar_a_M_ = dtau_viewing_polar_a_M_/max(1e-12,dtau_fnorm);
%dtau_viewing_azimu_b_M_ = dtau_viewing_azimu_b_M_/max(1e-12,dtau_fnorm);
%dtau_viewing_gamma_z_M_ = dtau_viewing_gamma_z_M_/max(1e-12,dtau_fnorm);
rng(2);
T_k_p_wM__ = randn(n_w,n_M);
parameter_apply = parameter;
parameter_apply.flag_verbose = 0;
parameter_apply.flag_check = 0;
parameter_apply.flag_disp = 0;
[ ...
 parameter ...
,a_restore_mid_k_Y_lm_ ...
,dtau_a_restore_mid_k_Y_lm_ ...
,dtau_dtau_a_restore_mid_k_Y_lm_ ...
] = ...
kappa_qpro_apply_0( ...
 parameter_apply ...
,qref_n_shell ...
,qref_k_c_qc__ ...
,l_max ...
,Y_l_val_ ...
,Y_m_val_ ...
,Ylm_weight_yq__ ...
,chebfun_kernel_norm_qpro_ ...
,n_nearest_north ...
,n_nearest_south ...
,deconvolve_lm_ ...
,n_M ...
,n_w ...
,T_k_p_wM__ ...
,viewing_polar_a_M_ ...
,viewing_azimu_b_M_ ...
,viewing_gamma_z_M_ ...
,dtau_viewing_polar_a_M_ ...
,dtau_viewing_azimu_b_M_ ...
,dtau_viewing_gamma_z_M_ ...
);
dtau = 1e-4;
[ ...
 parameter ...
,a_restore_pos_k_Y_lm_ ...
] = ...
kappa_qpro_apply_0( ...
 parameter_apply ...
,qref_n_shell ...
,qref_k_c_qc__ ...
,l_max ...
,Y_l_val_ ...
,Y_m_val_ ...
,Ylm_weight_yq__ ...
,chebfun_kernel_norm_qpro_ ...
,n_nearest_north ...
,n_nearest_south ...
,deconvolve_lm_ ...
,n_M ...
,n_w ...
,T_k_p_wM__ ...
,viewing_polar_a_M_ + dtau*dtau_viewing_polar_a_M_ ...
,viewing_azimu_b_M_ + dtau*dtau_viewing_azimu_b_M_ ...
,viewing_gamma_z_M_ + dtau*dtau_viewing_gamma_z_M_ ...
);
[ ...
 parameter ...
,a_restore_neg_k_Y_lm_ ...
] = ...
kappa_qpro_apply_0( ...
 parameter_apply ...
,qref_n_shell ...
,qref_k_c_qc__ ...
,l_max ...
,Y_l_val_ ...
,Y_m_val_ ...
,Ylm_weight_yq__ ...
,chebfun_kernel_norm_qpro_ ...
,n_nearest_north ...
,n_nearest_south ...
,deconvolve_lm_ ...
,n_M ...
,n_w ...
,T_k_p_wM__ ...
,viewing_polar_a_M_ - dtau*dtau_viewing_polar_a_M_ ...
,viewing_azimu_b_M_ - dtau*dtau_viewing_azimu_b_M_ ...
,viewing_gamma_z_M_ - dtau*dtau_viewing_gamma_z_M_ ...
);
dtau_a_restore_dif_k_Y_lm_ = (a_restore_pos_k_Y_lm_ - a_restore_neg_k_Y_lm_)/max(1e-12,2*dtau);
disp(sprintf(' %% dtau_a_restore_dif_k_Y_lm_ vs dtau_a_restore_mid_k_Y_lm_: %0.16f',fnorm(dtau_a_restore_dif_k_Y_lm_ - dtau_a_restore_mid_k_Y_lm_)/max(1e-12,fnorm(dtau_a_restore_dif_k_Y_lm_))));
%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,2,1);
plot(Y_l_val_,log10(abs(dtau_a_restore_dif_k_Y_lm_-dtau_a_restore_mid_k_Y_lm_)),'.');
xlim([0,l_max]); xlabel('Y_l_val_','Interpreter','none');
ylabel('log10(abs(error))','Interpreter','none');
title('log10(abs(dtau error))','Interpreter','none');
subplot(1,2,2);
plot(Y_l_val_,log10(abs(dtau_a_restore_dif_k_Y_lm_-dtau_a_restore_mid_k_Y_lm_)./abs(dtau_a_restore_dif_k_Y_lm_)),'.');
xlim([0,l_max]); xlabel('Y_l_val_','Interpreter','none');
ylabel('log10(abs(relative error))','Interpreter','none');
title('log10(abs(dtau relative error))','Interpreter','none');
end;%if flag_disp;
%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 7; p_col = ceil((1+l_max)/p_row); np=0;
for l_val=0:l_max;
subplot(p_row,p_col,1+np);np=np+1;cla;
tmp_index_ = efind(Y_l_val_==l_val);
hold on;
rlim_dif_ = prctile(real(dtau_a_restore_dif_k_Y_lm_(1+tmp_index_)),[  0,100]); rlim_dif_ = mean(rlim_dif_) + 0.5*1.25*diff(rlim_dif_)*[-1,+1]; if diff(rlim_dif_)==0; rlim_dif_=[-1,+1]; end;
rlim_mid_ = prctile(real(dtau_a_restore_mid_k_Y_lm_(1+tmp_index_)),[  0,100]); rlim_mid_ = mean(rlim_mid_) + 0.5*1.25*diff(rlim_mid_)*[-1,+1]; if diff(rlim_mid_)==0; rlim_mid_=[-1,+1]; end;
rlim_ = [min([rlim_dif_,rlim_mid_]),max([rlim_dif_,rlim_mid_])];
ilim_dif_ = prctile(imag(dtau_a_restore_dif_k_Y_lm_(1+tmp_index_)),[  0,100]); ilim_dif_ = mean(ilim_dif_) + 0.5*1.25*diff(ilim_dif_)*[-1,+1]; if diff(ilim_dif_)==0; ilim_dif_=[-1,+1]; end;
ilim_mid_ = prctile(imag(dtau_a_restore_mid_k_Y_lm_(1+tmp_index_)),[  0,100]); ilim_mid_ = mean(ilim_mid_) + 0.5*1.25*diff(ilim_mid_)*[-1,+1]; if diff(ilim_mid_)==0; ilim_mid_=[-1,+1]; end;
ilim_ = [min([ilim_dif_,ilim_mid_]),max([ilim_dif_,ilim_mid_])];
dlim_ = [min([rlim_,ilim_]),max([rlim_,ilim_])];
plot(dlim_,dlim_,'-','Color',0.85*[1,1,1]);
plot(real(dtau_a_restore_dif_k_Y_lm_(1+tmp_index_)),real(dtau_a_restore_mid_k_Y_lm_(1+tmp_index_)),'ro');
plot(imag(dtau_a_restore_dif_k_Y_lm_(1+tmp_index_)),imag(dtau_a_restore_mid_k_Y_lm_(1+tmp_index_)),'bx');
hold off;
xlim(dlim_); ylim(dlim_);
axisnotick;
xlabel('dif'); ylabel('mid')
title(sprintf('l_val %.2d',l_val),'Interpreter','none');
end;%for l_val=0:l_max;
sgtitle(sprintf('dtau_scatterplot'),'Interpreter','none');
end;%if flag_disp;
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% testing second-derivative (directional derivative only)')); end;
%%%%%%%%;
dtau_dtau_a_restore_dif_k_Y_lm_ = (a_restore_pos_k_Y_lm_ - 2*a_restore_mid_k_Y_lm_ + a_restore_neg_k_Y_lm_)/max(1e-12,dtau*dtau);
disp(sprintf(' %% dtau_dtau_a_restore_dif_k_Y_lm_ vs dtau_dtau_a_restore_mid_k_Y_lm_: %0.16f',fnorm(dtau_dtau_a_restore_dif_k_Y_lm_ - dtau_dtau_a_restore_mid_k_Y_lm_)/max(1e-12,fnorm(dtau_dtau_a_restore_dif_k_Y_lm_))));
%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,2,1);
plot(Y_l_val_,log10(abs(dtau_dtau_a_restore_dif_k_Y_lm_-dtau_dtau_a_restore_mid_k_Y_lm_)),'.');
xlim([0,l_max]); xlabel('Y_l_val_','Interpreter','none');
ylabel('log10(abs(error))','Interpreter','none');
title('log10(abs(dtau_dtau error))','Interpreter','none');
subplot(1,2,2);
plot(Y_l_val_,log10(abs(dtau_dtau_a_restore_dif_k_Y_lm_-dtau_dtau_a_restore_mid_k_Y_lm_)./abs(dtau_dtau_a_restore_dif_k_Y_lm_)),'.');
xlim([0,l_max]); xlabel('Y_l_val_','Interpreter','none');
ylabel('log10(abs(relative error))','Interpreter','none');
title('log10(abs(dtau_dtau relative error))','Interpreter','none');
end;%if flag_disp;
%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 7; p_col = ceil((1+l_max)/p_row); np=0;
for l_val=0:l_max;
subplot(p_row,p_col,1+np);np=np+1;cla;
tmp_index_ = efind(Y_l_val_==l_val);
hold on;
rlim_dif_ = prctile(real(dtau_dtau_a_restore_dif_k_Y_lm_(1+tmp_index_)),[  0,100]); rlim_dif_ = mean(rlim_dif_) + 0.5*1.25*diff(rlim_dif_)*[-1,+1]; if diff(rlim_dif_)==0; rlim_dif_=[-1,+1]; end;
rlim_mid_ = prctile(real(dtau_dtau_a_restore_mid_k_Y_lm_(1+tmp_index_)),[  0,100]); rlim_mid_ = mean(rlim_mid_) + 0.5*1.25*diff(rlim_mid_)*[-1,+1]; if diff(rlim_mid_)==0; rlim_mid_=[-1,+1]; end;
rlim_ = [min([rlim_dif_,rlim_mid_]),max([rlim_dif_,rlim_mid_])];
ilim_dif_ = prctile(imag(dtau_dtau_a_restore_dif_k_Y_lm_(1+tmp_index_)),[  0,100]); ilim_dif_ = mean(ilim_dif_) + 0.5*1.25*diff(ilim_dif_)*[-1,+1]; if diff(ilim_dif_)==0; ilim_dif_=[-1,+1]; end;
ilim_mid_ = prctile(imag(dtau_dtau_a_restore_mid_k_Y_lm_(1+tmp_index_)),[  0,100]); ilim_mid_ = mean(ilim_mid_) + 0.5*1.25*diff(ilim_mid_)*[-1,+1]; if diff(ilim_mid_)==0; ilim_mid_=[-1,+1]; end;
ilim_ = [min([ilim_dif_,ilim_mid_]),max([ilim_dif_,ilim_mid_])];
dlim_ = [min([rlim_,ilim_]),max([rlim_,ilim_])];
plot(dlim_,dlim_,'-','Color',0.85*[1,1,1]);
plot(real(dtau_dtau_a_restore_dif_k_Y_lm_(1+tmp_index_)),real(dtau_dtau_a_restore_mid_k_Y_lm_(1+tmp_index_)),'ro');
plot(imag(dtau_dtau_a_restore_dif_k_Y_lm_(1+tmp_index_)),imag(dtau_dtau_a_restore_mid_k_Y_lm_(1+tmp_index_)),'bx');
hold off;
xlim(dlim_); ylim(dlim_);
axisnotick;
xlabel('dif'); ylabel('mid')
title(sprintf('l_val %.2d',l_val),'Interpreter','none');
end;%for l_val=0:l_max;
sgtitle(sprintf('dtau_dtau_scatterplot'),'Interpreter','none');
end;%if flag_disp;
%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
disp('returning'); return;
end;%if nargin<1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); KAPPA=[]; end; na=na+1;

%%%%;
if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
flag_verbose=parameter.flag_verbose;
if ~isfield(parameter,'flag_disp'); parameter.flag_disp=0; end;
flag_disp=parameter.flag_disp; nf=0;
if ~isfield(parameter,'tolerance_master'); parameter.tolerance_master=1e-2; end;
tolerance_master=parameter.tolerance_master;
if ~isfield(parameter,'kernel_qpro_l_max_use'); parameter.kernel_qpro_l_max_use=49; end;
kernel_qpro_l_max_use=parameter.kernel_qpro_l_max_use;
if ~isfield(parameter,'kernel_qpro_l_max_ext'); parameter.kernel_qpro_l_max_ext=ceil(1.25*kernel_qpro_l_max_use); end;
kernel_qpro_l_max_ext=parameter.kernel_qpro_l_max_ext;
if ~isfield(parameter,'kernel_qpro_l_max_band'); parameter.kernel_qpro_l_max_band=+Inf; end;
kernel_qpro_l_max_band=parameter.kernel_qpro_l_max_band;
if ~isfield(parameter,'flag_kernel_qpro_d0'); parameter.flag_kernel_qpro_d0=1; end;
flag_kernel_qpro_d0=parameter.flag_kernel_qpro_d0;
if ~isfield(parameter,'flag_kernel_qpro_d1'); parameter.flag_kernel_qpro_d1=0; end;
flag_kernel_qpro_d1=parameter.flag_kernel_qpro_d1;
if ~isfield(parameter,'flag_kernel_qpro_d2'); parameter.flag_kernel_qpro_d2=0; end;
flag_kernel_qpro_d2=parameter.flag_kernel_qpro_d2;
if ~isfield(parameter,'kernel_qpro_polar_a_pole_north'); parameter.kernel_qpro_polar_a_pole_north=2.5*pi/12; end;
kernel_qpro_polar_a_pole_north=min(pi/2,parameter.kernel_qpro_polar_a_pole_north);
parameter.kernel_qpro_polar_a_pole_north = kernel_qpro_polar_a_pole_north;
if ~isfield(parameter,'kernel_qpro_polar_a_pole_south'); parameter.kernel_qpro_polar_a_pole_south=1.5*pi/12; end;
kernel_qpro_polar_a_pole_south=min(pi/2,parameter.kernel_qpro_polar_a_pole_south);
parameter.kernel_qpro_polar_a_pole_south = kernel_qpro_polar_a_pole_south;
if ~isfield(parameter,'kernel_qpro_deconvolution_factor_max'); parameter.kernel_qpro_deconvolution_factor_max=1024; end;
kernel_qpro_deconvolution_factor_max=parameter.kernel_qpro_deconvolution_factor_max;
if ~isfield(parameter,'kernel_qpro_qref_k_eq_d_double'); parameter.kernel_qpro_qref_k_eq_d_double=0.5; end;
kernel_qpro_qref_k_eq_d_double=parameter.kernel_qpro_qref_k_eq_d_double;
if ~isfield(parameter,'kernel_qpro_MaxIterations'); parameter.kernel_qpro_MaxIterations=1024; end;
kernel_qpro_MaxIterations=parameter.kernel_qpro_MaxIterations;
if ~isfield(parameter,'flag_kernel_k_Y_use'); parameter.flag_kernel_k_Y_use=1; end;
flag_kernel_k_Y_use=parameter.flag_kernel_k_Y_use;
%%%%;
if ~isfield(parameter,'flag_kernel_full'); parameter.flag_kernel_full=0; end;
flag_kernel_full=parameter.flag_kernel_full;
if (kernel_qpro_polar_a_pole_north + kernel_qpro_polar_a_pole_south > pi-1e-12);
flag_kernel_full = 1;
parameter.flag_kernel_full = flag_kernel_full;
end;%if (kernel_qpro_polar_a_pole_north + kernel_qpro_polar_a_pole_south > pi-1e-12);
%%%%;

if isempty(KAPPA); KAPPA = struct('type','KAPPA'); end;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

%%%%;
if ~isfield(KAPPA,'Rz');
Rz = @(azimu_b) ...
[ +cos(azimu_b) -sin(azimu_b) 0 ; ...
  +sin(azimu_b) +cos(azimu_b) 0 ; ...
   0             0            1 ; ...
] ;
KAPPA.Rz = Rz;
end;%if ~isfield(KAPPA,'Rz');
Rz = KAPPA.Rz;
%%%%;
if ~isfield(KAPPA,'dRz');
dRz = @(azimu_b) ...
[ -sin(azimu_b) -cos(azimu_b) 0 ; ...
  +cos(azimu_b) -sin(azimu_b) 0 ; ...
   0             0            0 ; ...
] ;
KAPPA.dRz = dRz;
end;%if ~isfield(KAPPA,'dRz');
dRz = KAPPA.dRz;
%%%%;
if ~isfield(KAPPA,'Ry');
Ry = @(polar_a) ...
[ +cos(polar_a) 0 +sin(polar_a) ; ...
   0            1  0            ; ...
  -sin(polar_a) 0 +cos(polar_a) ; ...
];
KAPPA.Ry = Ry;
end;%if ~isfield(KAPPA,'Ry');
Ry = KAPPA.Ry;
%%%%;
if ~isfield(KAPPA,'dRy');
dRy = @(polar_a) ...
[ -sin(polar_a) 0 +cos(polar_a) ; ...
   0            0  0            ; ...
  -cos(polar_a) 0 -sin(polar_a) ; ...
];
KAPPA.dRy = dRy;
end;%if ~isfield(KAPPA,'dRy');
dRy = KAPPA.dRy;
%%%%;
if flag_disp;
b = 2*pi*rand(); db = 1e-4; 
R_mid = Rz(b+0*db);
R_pos = Rz(b+1*db);
R_neg = Rz(b-1*db);
db_R_mid = dRz(b+0*db);
db_R_dif = (R_pos - R_neg)/max(1e-12,2*db);
if (flag_verbose>0); disp(sprintf(' %% db_R_dif vs db_R_mid: %0.16f',fnorm(db_R_dif-db_R_mid)/fnorm(db_R_dif))); end;
end;%if flag_disp;
%%%%;
if flag_disp;
a = 1*pi*rand(); da = 1e-3; 
R_mid = Ry(a+0*da);
R_pos = Ry(a+1*da);
R_neg = Ry(a-1*da);
da_R_mid = dRy(a+0*da);
da_R_dif = (R_pos - R_neg)/max(1e-12,2*da);
if (flag_verbose>0); disp(sprintf(' %% da_R_dif vs da_R_mid: %0.16f',fnorm(da_R_dif-da_R_mid)/fnorm(da_R_dif))); end;
end;%if flag_disp;
%%%%;

%%%%%%%%;
if ~isfield(KAPPA,'l_max_use');
l_max_use = kernel_qpro_l_max_use;
KAPPA.l_max_use = l_max_use;
end;%if ~isfield(KAPPA,'l_max_use');
l_max_use = KAPPA.l_max_use;
%%%%;
if ~isfield(KAPPA,'l_max_band');
l_max_band = kernel_qpro_l_max_band;
KAPPA.l_max_band = l_max_band;
end;%if ~isfield(KAPPA,'l_max_band');
l_max_band = KAPPA.l_max_band;
%%%%;
if ~isfield(KAPPA,'l_val_use_');
l_val_use_ = transpose([0:l_max_use]);
KAPPA.l_val_use_ = l_val_use_;
KAPPA.l_val_use_ = l_val_use_;
end;%if ~isfield(KAPPA,'l_val_use_');
l_val_use_ = KAPPA.l_val_use_;
%%%%;
if ~isfield(KAPPA,'l_max_ext');
l_max_ext = kernel_qpro_l_max_ext;
KAPPA.l_max_ext = l_max_ext;
end;%if ~isfield(KAPPA,'l_max_ext');
l_max_ext = KAPPA.l_max_ext;
%%%%;
if ~isfield(KAPPA,'l_val_ext_');
l_val_ext_ = transpose([0:l_max_ext]);
KAPPA.l_val_ext_ = l_val_ext_;
KAPPA.l_val_ext_ = l_val_ext_;
end;%if ~isfield(KAPPA,'l_val_ext_');
l_val_ext_ = KAPPA.l_val_ext_;
%%%%%%%%;
if ~isfield(KAPPA,'chebleg_ext_d_');
chebleg_ext_d_ = cell(1+l_max_ext,1);
for l_val_ext=0:l_max_ext;
tmp_c_ = [zeros(l_val_ext,1);1];
chebleg_ext_d_{1+l_val_ext} = chebfun(leg2cheb(tmp_c_,'norm'),'coeffs');
end;%for l_val_ext=0:l_max_ext;
KAPPA.chebleg_ext_d_ = chebleg_ext_d_;
end;%if ~isfield(KAPPA,'chebleg_ext_d_');
chebleg_ext_d_ = KAPPA.chebleg_ext_d_;
%%%%%%%%;

%%%%%%%%;
if ~isfield(KAPPA,'chebfun_kernel_norm_qpro_');
[ ...
 parameter ...
,kappa_norm_ ...
,chebfun_kernel_norm_qpro_ ...
,deconvolve_ext_l_ ...
,kappa_sparse_f ...
,relative_error_crop_ ...
,relative_error_full_ ...
,l_max_ext ...
,chebleg_ext_d_ ...
] = ...
kappa_qpro_0( ...
 parameter ...
,l_max_ext ...
,chebleg_ext_d_ ...
);
KAPPA.kappa_norm_ = kappa_norm_;
KAPPA.chebfun_kernel_norm_qpro_ = chebfun_kernel_norm_qpro_;
KAPPA.deconvolve_ext_l_ = deconvolve_ext_l_;
KAPPA.kappa_sparse_f = kappa_sparse_f;
KAPPA.relative_error_crop_ = relative_error_crop_;
KAPPA.relative_error_full_ = relative_error_full_;
end;%if ~isfield(KAPPA,'chebfun_kernel_norm_qpro_');
kappa_norm_ = KAPPA.kappa_norm_;
chebfun_kernel_norm_qpro_ = KAPPA.chebfun_kernel_norm_qpro_;
deconvolve_ext_l_ = KAPPA.deconvolve_ext_l_;
kappa_sparse_f = KAPPA.kappa_sparse_f;
relative_error_crop_ = KAPPA.relative_error_crop_;
relative_error_full_ = KAPPA.relative_error_full_;
%%%%%%%%;

%%%%%%%%;
if ~isfield(KAPPA,'a_full_weight_');
n_a_use = 1+2*l_max_ext + 16; %<-- Need to integrate polynomials of degree l_max_ext^2. ;
[a_drop_node_,a_drop_weight_] = legpts(n_a_use,[0+kernel_qpro_polar_a_pole_north,pi-kernel_qpro_polar_a_pole_south]);
[a_keep_node_north_,a_keep_weight_north_] = legpts(n_a_use,[0 ,kernel_qpro_polar_a_pole_north]);
[a_keep_node_south_,a_keep_weight_south_] = legpts(n_a_use,[pi-kernel_qpro_polar_a_pole_south,pi]);
a_keep_node_ = [a_keep_node_north_;a_keep_node_south_]; %<-- col vector. ;
a_keep_weight_ = [a_keep_weight_north_,a_keep_weight_south_]; %<-- row vector. ;
a_full_node_ = [a_keep_node_;a_drop_node_]; %<-- col vector. ;
a_full_weight_ = [a_keep_weight_,a_drop_weight_]; %<-- row vector. ;
KAPPA.n_a_use = n_a_use;
KAPPA.a_drop_node_ = a_drop_node_;
KAPPA.a_drop_weight_ = a_drop_weight_;
KAPPA.a_keep_node_north_ = a_keep_node_north_;
KAPPA.a_keep_weight_north_ = a_keep_weight_north_;
KAPPA.a_keep_node_south_ = a_keep_node_south_;
KAPPA.a_keep_weight_south_ = a_keep_weight_south_;
KAPPA.a_keep_node_ = a_keep_node_;
KAPPA.a_keep_weight_ = a_keep_weight_;
KAPPA.a_full_node_ = a_full_node_;
KAPPA.a_full_weight_ = a_full_weight_;
end;%if ~isfield(KAPPA,'a_full_weight_');
n_a_use = KAPPA.n_a_use;
a_drop_node_ = KAPPA.a_drop_node_;
a_drop_weight_ = KAPPA.a_drop_weight_;
a_keep_node_north_ = KAPPA.a_keep_node_north_;
a_keep_weight_north_ = KAPPA.a_keep_weight_north_;
a_keep_node_south_ = KAPPA.a_keep_node_south_;
a_keep_weight_south_ = KAPPA.a_keep_weight_south_;
a_keep_node_ = KAPPA.a_keep_node_;
a_keep_weight_ = KAPPA.a_keep_weight_;
a_full_node_ = KAPPA.a_full_node_;
a_full_weight_ = KAPPA.a_full_weight_;
%%%%%%%%;

%%%%%%%%;
if ~isfield(KAPPA,'qref_k_c_qc__');
n_lm_use = (l_max_use+1).^2;
n_lm_ext = (l_max_ext+1).^2;
qref_k_eq_d = kernel_qpro_qref_k_eq_d_double*sqrt(4*pi./max(1,n_lm_ext)); %<-- one half the older setting from qbp_6. increased density of quadrature points. ;
%%%%;
tmp_t = tic();
[ ...
 qref_n_shell ...
,qref_azimu_b_shell_ ...
,qref_polar_a_shell_ ...
,qref_weight_shell_ ...
,qref_k_c_0_shell_ ...
,qref_k_c_1_shell_ ...
,qref_k_c_2_shell_ ...
,qref_n_polar_a ...
,qref_polar_a_ ...
,qref_n_azimu_b_ ...
] = ...
sample_shell_5( ...
 1.0 ...
,qref_k_eq_d ...
,'L' ...
) ;
qref_k_c_qc__ = [ qref_k_c_0_shell_ , qref_k_c_1_shell_ , qref_k_c_2_shell_ ];
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% sample_shell_5: %0.2fs',tmp_t)); end;
%%%%;
n_ring_north = ceil(kernel_qpro_polar_a_pole_north/max(1e-12,qref_k_eq_d)); %<-- number of nearest neighbor-rings requested for each point. ;
n_nearest_north = 1+6*n_ring_north*(n_ring_north+1)/2; %<-- rough number of neighbors on hexagonal grid (i.e., 1+6+12+18+...). ;
n_nearest_north = min(floor(qref_n_shell/2),n_nearest_north);
if (flag_verbose>0); disp(sprintf(' %% qref_k_eq_d %0.6f n_nearest_north %d',qref_k_eq_d,n_nearest_north)); end;
n_ring_south = ceil(kernel_qpro_polar_a_pole_south/max(1e-12,qref_k_eq_d)); %<-- number of nearest neighbor-rings requested for each point. ;
n_nearest_south = 1+6*n_ring_south*(n_ring_south+1)/2; %<-- rough number of neighbors on hexagonal grid (i.e., 1+6+12+18+...). ;
n_nearest_south = min(floor(qref_n_shell/2),n_nearest_south);
if (flag_verbose>0); disp(sprintf(' %% qref_k_eq_d %0.6f n_nearest_south %d',qref_k_eq_d,n_nearest_south)); end;
n_nearest_total = min(qref_n_shell,n_nearest_north + n_nearest_south);
%%%%;
if flag_kernel_full;
n_nearest_total = qref_n_shell;
n_nearest_north = qref_n_shell;
n_nearest_south = 0;
end;%if flag_kernel_full;
%%%%;
KAPPA.n_lm_use = n_lm_use;
KAPPA.n_lm_ext = n_lm_ext;
KAPPA.flag_kernel_full = flag_kernel_full;
KAPPA.qref_k_eq_d = qref_k_eq_d;
KAPPA.n_ring_north = n_ring_north;
KAPPA.n_nearest_north = n_nearest_north;
KAPPA.n_ring_south = n_ring_south;
KAPPA.n_nearest_south = n_nearest_south;
KAPPA.n_nearest_total = n_nearest_total;
KAPPA.qref_n_shell = qref_n_shell;
KAPPA.qref_azimu_b_shell_ = qref_azimu_b_shell_;
KAPPA.qref_polar_a_shell_ = qref_polar_a_shell_;
KAPPA.qref_weight_shell_ = qref_weight_shell_;
KAPPA.qref_k_c_0_shell_ = qref_k_c_0_shell_;
KAPPA.qref_k_c_1_shell_ = qref_k_c_1_shell_;
KAPPA.qref_k_c_2_shell_ = qref_k_c_2_shell_;
KAPPA.qref_n_polar_a = qref_n_polar_a;
KAPPA.qref_polar_a_ = qref_polar_a_;
KAPPA.qref_n_azimu_b_ = qref_n_azimu_b_;
KAPPA.qref_k_c_qc__ = qref_k_c_qc__;
end;%if ~isfield(KAPPA,'qref_k_c_qc__');
n_lm_use = KAPPA.n_lm_use;
n_lm_ext = KAPPA.n_lm_ext;
qref_k_eq_d = KAPPA.qref_k_eq_d;
n_ring_north = KAPPA.n_ring_north;
n_nearest_north = KAPPA.n_nearest_north;
n_ring_south = KAPPA.n_ring_south;
n_nearest_south = KAPPA.n_nearest_south;
n_nearest_total = KAPPA.n_nearest_total;
qref_n_shell = KAPPA.qref_n_shell;
qref_azimu_b_shell_ = KAPPA.qref_azimu_b_shell_;
qref_polar_a_shell_ = KAPPA.qref_polar_a_shell_;
qref_weight_shell_ = KAPPA.qref_weight_shell_;
qref_k_c_0_shell_ = KAPPA.qref_k_c_0_shell_;
qref_k_c_1_shell_ = KAPPA.qref_k_c_1_shell_;
qref_k_c_2_shell_ = KAPPA.qref_k_c_2_shell_;
qref_n_polar_a = KAPPA.qref_n_polar_a;
qref_polar_a_ = KAPPA.qref_polar_a_;
qref_n_azimu_b_ = KAPPA.qref_n_azimu_b_;
qref_k_c_qc__ = KAPPA.qref_k_c_qc__;
%%%%%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;
markersize_use = 1;
plot_sphere_grid_0;
hold on;
plot3(qref_k_c_0_shell_,qref_k_c_1_shell_,qref_k_c_2_shell_,'k.','MarkerSize',markersize_use);
hold off;
axis equal; axisnotick3d; axis vis3d;
title('qref_k_c_?_shell_','Interpreter','none');
end;%if flag_disp;
%%%%%%%%;

%%%%%%%%;
if ~isfield(KAPPA,'Ylm_ext_weight_yq__');
%%%%;
if flag_kernel_k_Y_use;
tmp_t = tic();
Ylm_ext__ = get_Ylm__(1+l_max_ext,0:l_max_ext,qref_n_shell,qref_azimu_b_shell_,qref_polar_a_shell_);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% get_Ylm__: %0.2fs',tmp_t)); end;
end;%if flag_kernel_k_Y_use;
tmp_t = tic();
if flag_kernel_k_Y_use; Ylm_ext_yq__ = zeros(n_lm_ext,qref_n_shell); end;
Y_l_val_ext_ = zeros(n_lm_ext,1);
Y_m_val_ext_ = zeros(n_lm_ext,1);
nml=0;
for l_val_ext=0:l_max_ext;
for m_val_ext=-l_val_ext:+l_val_ext;
Y_l_val_ext_(1+nml) = l_val_ext;
Y_m_val_ext_(1+nml) = m_val_ext;
if flag_kernel_k_Y_use; Ylm_ext_yq__(1+nml,:) = Ylm_ext__{1+l_val_ext}(1+l_val_ext+m_val_ext,:); end;
nml=nml+1;
end;%for m_val_ext=-l_val_ext:+l_val_ext;
end;%for l_val_ext=0:l_max_ext;
if flag_kernel_k_Y_use; Ylm_ext_weight_yq__ = Ylm_ext_yq__ * sparse(1:qref_n_shell,1:qref_n_shell,qref_weight_shell_,qref_n_shell,qref_n_shell); end;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% Ylm_ext_weight_yq__: %0.2fs',tmp_t)); end;
%%%%;
KAPPA.Y_l_val_ext_ = Y_l_val_ext_;
KAPPA.Y_m_val_ext_ = Y_m_val_ext_;
KAPPA.Ylm_ext_weight_yq__ = []; if flag_kernel_k_Y_use; KAPPA.Ylm_ext_weight_yq__ = Ylm_ext_weight_yq__; end;
end;%if ~isfield(KAPPA,'Ylm_ext_weight_yq__');
%%%%%%%%;
Y_l_val_ext_ = KAPPA.Y_l_val_ext_;
Y_m_val_ext_ = KAPPA.Y_m_val_ext_;
Ylm_ext_weight_yq__ = KAPPA.Ylm_ext_weight_yq__;
%%%%%%%%;

%%%%%%%%;
if ~isfield(KAPPA,'Ylm_use_weight_yq__');
KAPPA.Y_l_val_use_ = KAPPA.Y_l_val_ext_(1:n_lm_use);
KAPPA.Y_m_val_use_ = KAPPA.Y_m_val_ext_(1:n_lm_use);
KAPPA.Ylm_use_weight_yq__ = []; if flag_kernel_k_Y_use; KAPPA.Ylm_use_weight_yq__ = KAPPA.Ylm_ext_weight_yq__(1:n_lm_use,:); end;
end;%if ~isfield(KAPPA,'Ylm_use_weight_yq__');
%%%%%%%%;
Y_l_val_use_ = KAPPA.Y_l_val_use_;
Y_m_val_use_ = KAPPA.Y_m_val_use_;
Ylm_use_weight_yq__ = KAPPA.Ylm_use_weight_yq__;
%%%%%%%%;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

