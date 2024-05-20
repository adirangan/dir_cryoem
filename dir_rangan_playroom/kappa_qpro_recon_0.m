function ...
[ ...
 parameter ...
,l_max ...
,kappa_norm_ ...
,chebfun_kernel_norm_qpro_ ...
,deconvolve_l_ ...
,kappa_sparse_f ...
,relative_error_crop_ ...
,relative_error_full_ ...
,chebleg_d_ ...
,k_p_r_max ...
,k_eq_d ...
,n_k_p_r ...
,k_p_r_1 ...
,k_p_r_ ...
,n_shell ...
,azimu_b_shell_ ...
,polar_a_shell_ ...
,weight_shell_ ...
,k_c_0_shell_ ...
,k_c_1_shell_ ...
,k_c_2_shell_ ...
,k_p_r_shell_ ...
,n_lm ...
,m_max_ ...
,n_m_max ...
,Y_l_val_ ...
,Y_m_val_ ...
,tmp_l_val_ ...
,tmp_m_val_ ...
,weight_Y_ ...
,weight_3d_k_p_r_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
] = ...
kappa_qpro_recon_0( ...
 parameter ...
,l_max ...
,kappa_norm_ ...
,chebfun_kernel_norm_qpro_ ...
,deconvolve_l_ ...
,kappa_sparse_f ...
,relative_error_crop_ ...
,relative_error_full_ ...
,chebleg_d_ ...
,k_p_r_max ...
,k_eq_d ...
,n_k_p_r ...
,k_p_r_1 ...
,k_p_r_ ...
,n_shell ...
,azimu_b_shell_ ...
,polar_a_shell_ ...
,weight_shell_ ...
,k_c_0_shell_ ...
,k_c_1_shell_ ...
,k_c_2_shell_ ...
,k_p_r_shell_ ...
,n_lm ...
,m_max_ ...
,n_m_max ...
,Y_l_val_ ...
,Y_m_val_ ...
,tmp_l_val_ ...
,tmp_m_val_ ...
,weight_Y_ ...
,weight_3d_k_p_r_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);

str_thisfunction = 'kappa_qpro_recon_0';

%%%%%%%%;
if nargin<1;
disp(sprintf(' %% testing %s',str_thisfunction));
disp('returning'); return;
end;%if nargin<1;
%%%%%%%%;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); l_max=[]; end; na=na+1;
if (nargin<1+na); kappa_norm_=[]; end; na=na+1;
if (nargin<1+na); chebfun_kernel_norm_qpro_=[]; end; na=na+1;
if (nargin<1+na); deconvolve_l_=[]; end; na=na+1;
if (nargin<1+na); kappa_sparse_f=[]; end; na=na+1;
if (nargin<1+na); relative_error_crop_=[]; end; na=na+1;
if (nargin<1+na); relative_error_full_=[]; end; na=na+1;
if (nargin<1+na); chebleg_d_=[]; end; na=na+1;
if (nargin<1+na); k_p_r_max=[]; end; na=na+1;
if (nargin<1+na); k_eq_d=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); k_p_r_1=[]; end; na=na+1;
if (nargin<1+na); k_p_r_=[]; end; na=na+1;
if (nargin<1+na); n_shell=[]; end; na=na+1;
if (nargin<1+na); azimu_b_shell_=[]; end; na=na+1;
if (nargin<1+na); polar_a_shell_=[]; end; na=na+1;
if (nargin<1+na); weight_shell_=[]; end; na=na+1;
if (nargin<1+na); k_c_0_shell_=[]; end; na=na+1;
if (nargin<1+na); k_c_1_shell_=[]; end; na=na+1;
if (nargin<1+na); k_c_2_shell_=[]; end; na=na+1;
if (nargin<1+na); k_p_r_shell_=[]; end; na=na+1;
if (nargin<1+na); n_lm=[]; end; na=na+1;
if (nargin<1+na); m_max_=[]; end; na=na+1;
if (nargin<1+na); n_m_max=[]; end; na=na+1;
if (nargin<1+na); Y_l_val_=[]; end; na=na+1;
if (nargin<1+na); Y_m_val_=[]; end; na=na+1;
if (nargin<1+na); tmp_l_val_=[]; end; na=na+1;
if (nargin<1+na); tmp_m_val_=[]; end; na=na+1;
if (nargin<1+na); weight_Y_=[]; end; na=na+1;
if (nargin<1+na); weight_3d_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); Ylm_uklma___ =[]; end; na=na+1;
if (nargin<1+na); k_p_azimu_b_sub_uka__ =[]; end; na=na+1;
if (nargin<1+na); k_p_polar_a_sub_uka__ =[]; end; na=na+1;
if (nargin<1+na); l_max_uk_ =[]; end; na=na+1;
if (nargin<1+na); index_nu_n_k_per_shell_from_nk_p_r_ =[]; end; na=na+1;
if (nargin<1+na); index_k_per_shell_uka__ =[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
flag_verbose=parameter.flag_verbose;
if ~isfield(parameter,'flag_disp'); parameter.flag_disp=0; end;
flag_disp=parameter.flag_disp; nf=0;
if ~isfield(parameter,'tolerance_master'); parameter.tolerance_master=1e-2; end;
tolerance_master=parameter.tolerance_master;
if ~isfield(parameter,'kernel_qpro_polar_a_pole_north'); parameter.kernel_qpro_polar_a_pole_north=1.0*pi/12; end;
kernel_qpro_polar_a_pole_north=parameter.kernel_qpro_polar_a_pole_north;
if ~isfield(parameter,'kernel_qpro_polar_a_pole_south'); parameter.kernel_qpro_polar_a_pole_south=0.0*pi/12; end;
kernel_qpro_polar_a_pole_south=parameter.kernel_qpro_polar_a_pole_south;
if ~isfield(parameter,'kernel_qpro_deconvolution_factor_max'); parameter.kernel_qpro_deconvolution_factor_max=1024; end;
kernel_qpro_deconvolution_factor_max=parameter.kernel_qpro_deconvolution_factor_max;
if ~isfield(parameter,'kernel_qpro_MaxIterations'); parameter.kernel_qpro_MaxIterations=1024; end;
kernel_qpro_MaxIterations=parameter.kernel_qpro_MaxIterations;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

%%%%;
Rz = @(azimu_b) ...
[ +cos(azimu_b) -sin(azimu_b) 0 ; ...
  +sin(azimu_b) +cos(azimu_b) 0 ; ...
   0             0            1 ; ...
] ;
dRz = @(azimu_b) ...
[ -sin(azimu_b) -cos(azimu_b) 0 ; ...
  +cos(azimu_b) -sin(azimu_b) 0 ; ...
   0             0            0 ; ...
] ;
Ry = @(polar_a) ...
[ +cos(polar_a) 0 +sin(polar_a) ; ...
   0            1  0            ; ...
  -sin(polar_a) 0 +cos(polar_a) ; ...
];
dRy = @(polar_a) ...
[ -sin(polar_a) 0 +cos(polar_a) ; ...
   0            0  0            ; ...
  -cos(polar_a) 0 -sin(polar_a) ; ...
];
%%%%;
b = 2*pi*rand(); db = 1e-4; 
R_mid = Rz(b+0*db);
R_pos = Rz(b+1*db);
R_pre = Rz(b-1*db);
DR_mid = dRz(b+0*db);
DR_dif = (R_pos - R_pre)/max(1e-12,2*db);
disp(sprintf(' %% DR_dif vs DR_mid: %0.16f',fnorm(DR_dif-DR_mid)/fnorm(DR_dif)));
%%%%;
a = 1*pi*rand(); da = 1e-3; 
R_mid = Ry(a+0*da);
R_pos = Ry(a+1*da);
R_pre = Ry(a-1*da);
DR_mid = dRy(a+0*da);
DR_dif = (R_pos - R_pre)/max(1e-12,2*da);
disp(sprintf(' %% DR_dif vs DR_mid: %0.16f',fnorm(DR_dif-DR_mid)/fnorm(DR_dif)));
%%%%;

%%%%%%%%;
if isempty(l_max); l_max = 49; end;
l_val_ = transpose([0:l_max]);
%%%%%%%%;
if isempty(chebleg_d_);
chebleg_d_ = cell(1+l_max,1);
for l_val=0:l_max;
tmp_c_ = [zeros(l_val,1);1];
chebleg_d_{1+l_val} = chebfun(leg2cheb(tmp_c_,'norm'),'coeffs');
end;%for l_val=0:l_max;
end;%if isempty(chebleg_d_);
%%%%%%%%;

%%%%%%%%;
if isempty(kappa_norm_);
[ ...
 parameter ...
,kappa_norm_ ...
,chebfun_kernel_norm_qpro_ ...
,deconvolve_l_ ...
,kappa_sparse_f ...
,relative_error_crop_ ...
,relative_error_full_ ...
,l_max ...
,chebleg_d_ ...
] = ...
kappa_qpro_0( ...
 parameter ...
,l_max ...
,chebleg_d_ ...
);
end;%if isempty(kappa_norm_);
%%%%%%%%;

n_a_use = 1+2*l_max + 16; %<-- Need to integrate polynomials of degree l_max^2. ;
[a_drop_node_,a_drop_weight_] = legpts(n_a_use,[0+kernel_qpro_polar_a_pole_north,pi-kernel_qpro_polar_a_pole_south]);
[a_keep_node_north_,a_keep_weight_north_] = legpts(n_a_use,[0 ,kernel_qpro_polar_a_pole_north]);
[a_keep_node_south_,a_keep_weight_south_] = legpts(n_a_use,[pi-kernel_qpro_polar_a_pole_south,pi]);
a_keep_node_ = [a_keep_node_north_;a_keep_node_south_]; %<-- col vector. ;
a_keep_weight_ = [a_keep_weight_north_,a_keep_weight_south_]; %<-- row vector. ;
a_full_node_ = [a_keep_node_;a_drop_node_]; %<-- col vector. ;
a_full_weight_ = [a_keep_weight_,a_drop_weight_]; %<-- row vector. ;

if isempty(k_p_r_max);
%%%%;
% Now we reconstruct the kernel in k_p_ ;
% and see how well resolved it is on polar_a_shell_ and azimu_b_shell_. ;
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
%%%%%%%%;
n_lm = (l_max+1).^2;
m_max_ = -l_max : +l_max;
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
%%%%%%%%;
end;%if isempty(k_p_r_max);

%%%%%%%%;
a_deltafunct_I_I_k_Y_lm_ = sqrt(2*pi)*sqrt(4*pi)*sqrt(1+2*Y_l_val_).*(Y_m_val_==0);
if ~exist('Ylm_uklma___','var'); Ylm_uklma___=[]; end;
if ~exist('k_p_azimu_b_sub_uka__','var'); k_p_azimu_b_sub_uka__ = []; end;
if ~exist('k_p_polar_a_sub_uka__','var'); k_p_polar_a_sub_uka__ = []; end;
if ~exist('l_max_uk_','var'); l_max_uk_ = []; end;
if ~exist('index_nu_n_k_per_shell_from_nk_p_r_','var'); index_nu_n_k_per_shell_from_nk_p_r_ = []; end;
if ~exist('index_k_per_shell_uka__','var'); index_k_per_shell_uka__ = []; end;
[ ...
 a_deltafunct_I_I_k_p_quad_ ...
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
,a_deltafunct_I_I_k_Y_lm_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
%%%%%%%%;

%%%%%%%%;
kernel_qpro_k_Y_form_ = zeros(n_lm,1);
kernel_qpro_k_Y_form_(1+l_val_.^2 + l_val_) = kappa_norm_;
%%%%;
deconvolve_lm_ = (sqrt(4*pi)*sqrt(1+2*Y_l_val_))./max(1e-12,kappa_norm_(1+Y_l_val_));
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
%qref_k_eq_d = 1.00*sqrt(4*pi./max(1,n_lm)); %<-- older setting from qbp_6. consider changing. ;
 qref_k_eq_d = 0.50*sqrt(4*pi./max(1,n_lm)); %<-- one half the older setting from qbp_6. increased density of quadrature points. ;
%qref_k_eq_d = 0.25*sqrt(4*pi./max(1,n_lm)); %<-- one quarter the older setting from qbp_6. increased density of quadrature points. ;
n_ring_north = ceil(kernel_qpro_polar_a_pole_north/max(1e-12,qref_k_eq_d)); %<-- number of nearest neighbor-rings requested for each point. ;
n_nearest_north = 1+6*n_ring_north*(n_ring_north+1)/2; %<-- rough number of neighbors on hexagonal grid (i.e., 1+6+12+18+...). ;
if (flag_verbose>0); disp(sprintf(' %% qref_k_eq_d %0.6f n_nearest_north %d',qref_k_eq_d,n_nearest_north)); end;
n_ring_south = ceil(kernel_qpro_polar_a_pole_south/max(1e-12,qref_k_eq_d)); %<-- number of nearest neighbor-rings requested for each point. ;
n_nearest_south = 1+6*n_ring_south*(n_ring_south+1)/2; %<-- rough number of neighbors on hexagonal grid (i.e., 1+6+12+18+...). ;
if (flag_verbose>0); disp(sprintf(' %% qref_k_eq_d %0.6f n_nearest_south %d',qref_k_eq_d,n_nearest_south)); end;
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
qref_k_c_qd__ = [ qref_k_c_0_shell_ , qref_k_c_1_shell_ , qref_k_c_2_shell_ ];
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% sample_shell_5 (should be a precomputation): %0.2fs',tmp_t)); end;
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
%%%%;
tmp_t = tic();
Ylm__ = get_Ylm__(1+l_max,0:l_max,qref_n_shell,qref_azimu_b_shell_,qref_polar_a_shell_);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% get_Ylm__ (should be a precomputation): %0.2fs',tmp_t)); end;
tmp_t = tic();
Ylm_yq__ = zeros(n_lm,qref_n_shell);
Y_l_val_ = zeros(n_lm,1);
Y_m_val_ = zeros(n_lm,1);
nml=0;
for l_val=0:l_max;
for m_val=-l_val:+l_val;
Y_l_val_(1+nml) = l_val;
Y_m_val_(1+nml) = m_val;
Ylm_yq__(1+nml,:) = Ylm__{1+l_val}(1+l_val+m_val,:);
nml=nml+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max;
Ylm_weight_yq__ = Ylm_yq__ * sparse(1:qref_n_shell,1:qref_n_shell,qref_weight_shell_,qref_n_shell,qref_n_shell);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% Ylm_weight_yq__ (should be a precomputation): %0.2fs',tmp_t)); end;
%%%%;
tmp_I__ = eye(3,3);
tmp_I_pole_ = tmp_I__*[0;0;1];
tmp_I_d2_ = (qref_k_c_0_shell_ - tmp_I_pole_(1+0)).^2 + (qref_k_c_1_shell_ - tmp_I_pole_(1+1)).^2 + (qref_k_c_2_shell_ - tmp_I_pole_(1+2)).^2 ;
tmp_I_arc_ = acos(1-tmp_I_d2_/2);
tmp_I_ker_ = chebfun_kernel_norm_qpro_(cos(tmp_I_arc_));
[tmp_ij_north_] = knnsearch(qref_k_c_qd__,reshape(+tmp_I_pole_,[1,3]),'K',n_nearest_north);
[tmp_ij_south_] = knnsearch(qref_k_c_qd__,reshape(-tmp_I_pole_,[1,3]),'K',n_nearest_south);
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
[tmp_ij_north_] = knnsearch(qref_k_c_qd__,reshape(+tmp_R_pole_,[1,3]),'K',n_nearest_north);
[tmp_ij_south_] = knnsearch(qref_k_c_qd__,reshape(-tmp_R_pole_,[1,3]),'K',n_nearest_south);
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
a_deltafunct_R_I_k_Y_lm_ = rotate_spharm_to_spharm_3(1,1,l_max,a_deltafunct_I_I_k_Y_lm_,tmp_euler_pos_);
a_deltafunct_I_R_k_Y_lm_ = a_deltafunct_R_I_k_Y_lm_; %<-- use formula. ;
a_deltafunct_R_R_k_Y_lm_ = rotate_spharm_to_spharm_3(1,1,l_max,a_deltafunct_I_R_k_Y_lm_,tmp_euler_neg_);
if (flag_verbose>0); disp(sprintf(' %% a_deltafunct_I_I_k_Y_lm_ vs a_deltafunct_R_R_k_Y_lm_: %0.16f',fnorm(a_deltafunct_I_I_k_Y_lm_ - a_deltafunct_R_R_k_Y_lm_)/fnorm(a_deltafunct_I_I_k_Y_lm_))); end;
a_deconvolve_I_I_k_Y_lm_ = a_I_I_k_Y_lm_.*deconvolve_lm_;
a_deconvolve_I_R_k_Y_lm_ = a_I_R_k_Y_lm_.*deconvolve_lm_;
a_deconvolve_R_I_k_Y_lm_ = a_R_I_k_Y_lm_.*deconvolve_lm_;
a_deconvolve_R_R_k_Y_lm_ = a_R_R_k_Y_lm_.*deconvolve_lm_;
if (flag_verbose>0); disp(sprintf(' %% a_deltafunct_I_I_k_Y_lm_ vs a_deconvolve_I_I_k_Y_lm_: %0.16f',fnorm(a_deltafunct_I_I_k_Y_lm_ - a_deconvolve_I_I_k_Y_lm_)/fnorm(a_deltafunct_I_I_k_Y_lm_))); end;
if (flag_verbose>0); disp(sprintf(' %% a_deltafunct_I_R_k_Y_lm_ vs a_deconvolve_I_R_k_Y_lm_: %0.16f',fnorm(a_deltafunct_I_R_k_Y_lm_ - a_deconvolve_I_R_k_Y_lm_)/fnorm(a_deltafunct_I_R_k_Y_lm_))); end;
if (flag_verbose>0); disp(sprintf(' %% a_deltafunct_R_I_k_Y_lm_ vs a_deconvolve_R_I_k_Y_lm_: %0.16f',fnorm(a_deltafunct_R_I_k_Y_lm_ - a_deconvolve_R_I_k_Y_lm_)/fnorm(a_deltafunct_R_I_k_Y_lm_))); end;
if (flag_verbose>0); disp(sprintf(' %% a_deltafunct_R_R_k_Y_lm_ vs a_deconvolve_R_R_k_Y_lm_: %0.16f',fnorm(a_deltafunct_R_R_k_Y_lm_ - a_deconvolve_R_R_k_Y_lm_)/fnorm(a_deltafunct_R_R_k_Y_lm_))); end;
a_deconvolve_I_I_crop_k_Y_lm_ = a_I_I_crop_k_Y_lm_.*deconvolve_lm_;
a_deconvolve_I_R_crop_k_Y_lm_ = a_I_R_crop_k_Y_lm_.*deconvolve_lm_;
a_deconvolve_R_I_crop_k_Y_lm_ = a_R_I_crop_k_Y_lm_.*deconvolve_lm_;
a_deconvolve_R_R_crop_k_Y_lm_ = a_R_R_crop_k_Y_lm_.*deconvolve_lm_;
if (flag_verbose>0); disp(sprintf(' %% a_deltafunct_I_I_k_Y_lm_ vs a_deconvolve_I_I_crop_k_Y_lm_: %0.16f',fnorm(a_deltafunct_I_I_k_Y_lm_ - a_deconvolve_I_I_crop_k_Y_lm_)/fnorm(a_deltafunct_I_I_k_Y_lm_))); end;
if (flag_verbose>0); disp(sprintf(' %% a_deltafunct_I_R_k_Y_lm_ vs a_deconvolve_I_R_crop_k_Y_lm_: %0.16f',fnorm(a_deltafunct_I_R_k_Y_lm_ - a_deconvolve_I_R_crop_k_Y_lm_)/fnorm(a_deltafunct_I_R_k_Y_lm_))); end;
if (flag_verbose>0); disp(sprintf(' %% a_deltafunct_R_I_k_Y_lm_ vs a_deconvolve_R_I_crop_k_Y_lm_: %0.16f',fnorm(a_deltafunct_R_I_k_Y_lm_ - a_deconvolve_R_I_crop_k_Y_lm_)/fnorm(a_deltafunct_R_I_k_Y_lm_))); end;
if (flag_verbose>0); disp(sprintf(' %% a_deltafunct_R_R_k_Y_lm_ vs a_deconvolve_R_R_crop_k_Y_lm_: %0.16f',fnorm(a_deltafunct_R_R_k_Y_lm_ - a_deconvolve_R_R_crop_k_Y_lm_)/fnorm(a_deltafunct_R_R_k_Y_lm_))); end;
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
deltafunct_lim_ = [0,prctile(a_deltafunct_I_I_k_p_quad_,99.5)];
subplot(p_row,p_col,1+np);np=np+1;
flag_2d_vs_3d = 0;
k_p_r_max = 1;
imagesc_polar_a_azimu_b_0( ...
 polar_a_shell_ ... 
,azimu_b_shell_ ... 
,a_deltafunct_I_I_k_p_quad_ ... 
,deltafunct_lim_ ... 
,colormap_81s ... 
,flag_2d_vs_3d ...
,k_p_r_max ...
);
axisnotick3d; axis vis3d;
title('deltafunct_I_I_k_p_form_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
flag_2d_vs_3d = 0;
k_p_r_max = 1;
imagesc_polar_a_azimu_b_0( ...
 polar_a_shell_ ... 
,azimu_b_shell_ ... 
,real(a_deconvolve_I_I_crop_k_p_quad_) ... 
,deltafunct_lim_ ... 
,colormap_81s ... 
,flag_2d_vs_3d ...
,k_p_r_max ...
);
axisnotick3d; axis vis3d;
title('deconvolve_I_I_crop_k_p_form_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
plot(Y_l_val_,log10(abs(a_deltafunct_I_R_k_Y_lm_-a_deconvolve_I_R_crop_k_Y_lm_)./abs(a_deltafunct_I_R_k_Y_lm_)),'.');
xlabel('Y_l_val_','Interpreter','none'); xlim([0,l_max+1]);
ylabel('log10(abs(delta-decon)./abs(delta)))','Interpreter','none');
grid on;
end;%if (flag_disp);
%%%%%%%%;

%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% testing first-derivative')); end;
%%%%%%%%;
tmp_I__ = eye(3,3);
tmp_I_pole_ = tmp_I__*[0;0;1];
tmp_I_d2_ = (qref_k_c_0_shell_ - tmp_I_pole_(1+0)).^2 + (qref_k_c_1_shell_ - tmp_I_pole_(1+1)).^2 + (qref_k_c_2_shell_ - tmp_I_pole_(1+2)).^2 ;
tmp_I_arc_ = acos(1-tmp_I_d2_/2);
tmp_I_ker_ = chebfun_kernel_norm_qpro_(cos(tmp_I_arc_));
[tmp_ij_north_] = knnsearch(qref_k_c_qd__,reshape(+tmp_I_pole_,[1,3]),'K',n_nearest_north);
[tmp_ij_south_] = knnsearch(qref_k_c_qd__,reshape(-tmp_I_pole_,[1,3]),'K',n_nearest_south);
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
[tmp_ij_north_] = knnsearch(qref_k_c_qd__,reshape(+tmp_R_pole_,[1,3]),'K',n_nearest_north);
[tmp_ij_south_] = knnsearch(qref_k_c_qd__,reshape(-tmp_R_pole_,[1,3]),'K',n_nearest_south);
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
tmp_R_pre__ = Rz(tmp_euler_b-db)*Ry(tmp_euler_a-da)*Rz(tmp_euler_c-dc);
tmp_dR_dif__ = (tmp_R_pos__ - tmp_R_pre__)/max(1e-12,2*fnorm([da;db;dc]));;
disp(sprintf(' %% tmp_dR_dif__ vs tmp_dR_mid__: %0.16f',fnorm(tmp_dR_dif__ - tmp_dR_mid__)/fnorm(tmp_dR_dif__)));
tmp_R_pos_pole_ = tmp_R_pos__*[0;0;1];
tmp_R_pre_pole_ = tmp_R_pre__*[0;0;1];
tmp_dRda_pole_ = tmp_dRda__*[0;0;1];
tmp_dRdb_pole_ = tmp_dRdb__*[0;0;1];
tmp_dRdc_pole_ = tmp_dRdc__*[0;0;1];
tmp_dR_mid_d2_ = ...
 - 2*(qref_k_c_0_shell_ - tmp_R_pole_(1+0)).*(tmp_dRda_pole_(1+0)*da + tmp_dRdb_pole_(1+0)*db + tmp_dRdc_pole_(1+0)*dc) ...
 - 2*(qref_k_c_1_shell_ - tmp_R_pole_(1+1)).*(tmp_dRda_pole_(1+1)*da + tmp_dRdb_pole_(1+1)*db + tmp_dRdc_pole_(1+1)*dc) ...
 - 2*(qref_k_c_2_shell_ - tmp_R_pole_(1+2)).*(tmp_dRda_pole_(1+2)*da + tmp_dRdb_pole_(1+2)*db + tmp_dRdc_pole_(1+2)*dc) ...
;
tmp_R_pos_d2_ = (qref_k_c_0_shell_ - tmp_R_pos_pole_(1+0)).^2 + (qref_k_c_1_shell_ - tmp_R_pos_pole_(1+1)).^2 + (qref_k_c_2_shell_ - tmp_R_pos_pole_(1+2)).^2 ;
tmp_R_pre_d2_ = (qref_k_c_0_shell_ - tmp_R_pre_pole_(1+0)).^2 + (qref_k_c_1_shell_ - tmp_R_pre_pole_(1+1)).^2 + (qref_k_c_2_shell_ - tmp_R_pre_pole_(1+2)).^2 ;
tmp_dR_dif_d2_ = (tmp_R_pos_d2_ - tmp_R_pre_d2_)/2;
disp(sprintf(' %% tmp_dR_dif_d2_ vs tmp_dR_mid_d2_: %0.16f',fnorm(tmp_dR_dif_d2_ - tmp_dR_mid_d2_)/fnorm(tmp_dR_dif_d2_)));
dchebfun_kernel_norm_qpro_ = diff(chebfun_kernel_norm_qpro_);
tmp_dR_mid_ker_ = -dchebfun_kernel_norm_qpro_(1-tmp_R_d2_/2).*tmp_dR_mid_d2_/2;
tmp_R_pos_ker_ = chebfun_kernel_norm_qpro_(1-tmp_R_pos_d2_/2);
tmp_R_pre_ker_ = chebfun_kernel_norm_qpro_(1-tmp_R_pre_d2_/2);
tmp_dR_dif_ker_ = (tmp_R_pos_ker_ - tmp_R_pre_ker_)/2;
disp(sprintf(' %% tmp_dR_dif_ker_ vs tmp_dR_mid_ker_: %0.16f',fnorm(tmp_dR_dif_ker_ - tmp_dR_mid_ker_)/fnorm(tmp_dR_dif_ker_)));
tmp_dR_mid_ker_crop_ = tmp_dR_mid_ker_; tmp_dR_mid_ker_crop_(1+tmp_index_drop_)=0;
tmp_R_pos_ker_crop_ = tmp_R_pos_ker_; tmp_R_pos_ker_crop_(1+tmp_index_drop_)=0;
tmp_R_pre_ker_crop_ = tmp_R_pre_ker_; tmp_R_pre_ker_crop_(1+tmp_index_drop_)=0;
%%%%;
a_I_dR_mid_k_Y_lm_ = conj(Ylm_weight_yq__)*tmp_dR_mid_ker_;
a_I_R_pos_k_Y_lm_ = conj(Ylm_weight_yq__)*tmp_R_pos_ker_;
a_I_R_pre_k_Y_lm_ = conj(Ylm_weight_yq__)*tmp_R_pre_ker_;
a_I_dR_dif_k_Y_lm_ = (a_I_R_pos_k_Y_lm_ - a_I_R_pre_k_Y_lm_)/2;
disp(sprintf(' %% a_I_dR_dif_k_Y_lm_ vs a_I_dR_mid_k_Y_lm_: %0.16f',fnorm(a_I_dR_dif_k_Y_lm_ - a_I_dR_mid_k_Y_lm_)/fnorm(a_I_dR_dif_k_Y_lm_)));
%%%%;
a_I_dR_mid_crop_k_Y_lm_ = conj(Ylm_weight_yq__)*tmp_dR_mid_ker_crop_;
a_I_R_pos_crop_k_Y_lm_ = conj(Ylm_weight_yq__)*tmp_R_pos_ker_crop_;
a_I_R_pre_crop_k_Y_lm_ = conj(Ylm_weight_yq__)*tmp_R_pre_ker_crop_;
a_I_dR_dif_crop_k_Y_lm_ = (a_I_R_pos_crop_k_Y_lm_ - a_I_R_pre_crop_k_Y_lm_)/2;
disp(sprintf(' %% a_I_dR_dif_crop_k_Y_lm_ vs a_I_dR_mid_crop_k_Y_lm_: %0.16f',fnorm(a_I_dR_dif_crop_k_Y_lm_ - a_I_dR_mid_crop_k_Y_lm_)/fnorm(a_I_dR_dif_crop_k_Y_lm_)));
%%%%;
a_deconvolve_I_dR_mid_k_Y_lm_ = a_I_dR_mid_k_Y_lm_.*deconvolve_lm_;
a_deconvolve_I_dR_mid_crop_k_Y_lm_ = a_I_dR_mid_crop_k_Y_lm_.*deconvolve_lm_;
a_deconvolve_I_R_pos_k_Y_lm_ = a_I_R_pos_k_Y_lm_.*deconvolve_lm_;
a_deconvolve_I_R_pos_crop_k_Y_lm_ = a_I_R_pos_crop_k_Y_lm_.*deconvolve_lm_;
a_deconvolve_I_R_pre_k_Y_lm_ = a_I_R_pre_k_Y_lm_.*deconvolve_lm_;
a_deconvolve_I_R_pre_crop_k_Y_lm_ = a_I_R_pre_crop_k_Y_lm_.*deconvolve_lm_;
a_deconvolve_I_dR_dif_k_Y_lm_ = (a_deconvolve_I_R_pos_k_Y_lm_ - a_deconvolve_I_R_pre_k_Y_lm_)/2;
a_deconvolve_I_dR_dif_crop_k_Y_lm_ = (a_deconvolve_I_R_pos_crop_k_Y_lm_ - a_deconvolve_I_R_pre_crop_k_Y_lm_)/2;
disp(sprintf(' %% a_deconvolve_I_dR_dif_k_Y_lm_ vs a_deconvolve_I_dR_mid_k_Y_lm_: %0.16f',fnorm(a_deconvolve_I_dR_dif_k_Y_lm_ - a_deconvolve_I_dR_mid_k_Y_lm_)/fnorm(a_deconvolve_I_dR_dif_k_Y_lm_)));
disp(sprintf(' %% a_deconvolve_I_dR_dif_crop_k_Y_lm_ vs a_deconvolve_I_dR_mid_crop_k_Y_lm_: %0.16f',fnorm(a_deconvolve_I_dR_dif_crop_k_Y_lm_ - a_deconvolve_I_dR_mid_crop_k_Y_lm_)/fnorm(a_deconvolve_I_dR_dif_crop_k_Y_lm_)));
%%%%%%%%;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

