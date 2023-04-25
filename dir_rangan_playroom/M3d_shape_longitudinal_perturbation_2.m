function ...
[ ...
 parameter ...
] = ...
M3d_shape_longitudinal_perturbation_2( ...
 parameter ...
);

str_thisfunction = 'M3d_shape_longitudinal_perturbation_2';

if (nargin<1);
disp(sprintf(' %% testing %s',str_thisfunction));
parameter = struct('type','parameter');
parameter.flag_verbose = 1;
for flag_sign = [-1,+1];
for d_source = [0,0.125*pi,0.250*pi];
for sigma_x_c = [0.15,0.10,0.05,0.025];
parameter.flag_sign = flag_sign;
parameter.d_source = d_source;
parameter.sigma_x_c = sigma_x_c;
M3d_shape_longitudinal_perturbation_2(parameter);
end;%for sigma_x_c = [0.10,0.05,0.025];
end;%for d_source = [0,0.125*pi,0.250*pi];
end;%for flag_sign = [-1,+1];
disp(sprintf(' %% returning')); return;
end;%if (nargin<1);

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;

%clear; parameter = []; str_thisfunction = 'M3d_shape_longitudinal_perturbation_2';

if isempty(parameter); parameter = struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
flag_verbose = parameter.flag_verbose;
if ~isfield(parameter,'flag_disp'); parameter.flag_disp = 0; end;
flag_disp = parameter.flag_disp; nf=0;
if ~isfield(parameter,'flag_replot'); parameter.flag_replot = 0; end;
flag_replot = parameter.flag_replot;
if ~isfield(parameter,'str_shape'); parameter.str_shape = 'rand'; end;
str_shape = parameter.str_shape;
if ~isfield(parameter,'flag_sign'); parameter.flag_sign = -1; end;
flag_sign = parameter.flag_sign;
if ~isfield(parameter,'n_mode'); parameter.n_mode = 2; end;
n_mode = parameter.n_mode;
if ~isfield(parameter,'d_source'); parameter.d_source = 0.25*pi; end;
d_source = parameter.d_source;
if ~isfield(parameter,'sigma_x_c'); parameter.sigma_x_c = 0.05; end;
sigma_x_c = parameter.sigma_x_c;
if ~isfield(parameter,'str_shape'); parameter.str_shape = 'rand'; end;
str_shape = parameter.str_shape;
if ~isfield(parameter,'rseed'); parameter.rseed = 0; end;
rseed = parameter.rseed;
if ~isfield(parameter,'n_x_c'); parameter.n_x_c = 256; end;
n_x_c = parameter.n_x_c;
if ~isfield(parameter,'k_p_r_max'); parameter.k_p_r_max = 2.0*48.0/(2*pi); end;
k_p_r_max = parameter.k_p_r_max;
if ~isfield(parameter,'k_eq_d'); parameter.k_eq_d = sqrt(1.0)/(2*pi); end;
k_eq_d = parameter.k_eq_d;

if flag_verbose; disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

%%%%%%%%;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
%%%%%%%%;
dir_manuscript = sprintf('/%s/rangan/dir_cryoem/dir_spurious_heterogeneity_manuscript',string_root);
if ~exist(dir_manuscript,'dir'); disp(sprintf(' %% mkdir %s',dir_manuscript)); mkdir(dir_manuscript); end;
dir_manuscript_jpg = sprintf('%s/dir_M3d_shape_longitudinal_perturbation_jpg',dir_manuscript);
if ~exist(dir_manuscript_jpg,'dir'); disp(sprintf(' %% mkdir %s',dir_manuscript_jpg)); mkdir(dir_manuscript_jpg); end;

if strcmp(str_shape,'rand'); rng(rseed); end;

%%%%%%%%;
flag_unif_vs_adap = 1;
flag_tensor_vs_adap = 1;
half_diameter_x_c = 1.0d0;
diameter_x_c = 2.0d0*half_diameter_x_c;
x_p_r_max = 1.0;
x_c_0_ = -half_diameter_x_c + transpose([0:n_x_c-1]/n_x_c)*diameter_x_c;
x_c_1_ = -half_diameter_x_c + transpose([0:n_x_c-1]/n_x_c)*diameter_x_c;
x_c_2_ = -half_diameter_x_c + transpose([0:n_x_c-1]/n_x_c)*diameter_x_c;
dx = diameter_x_c/n_x_c;
[x_c_0_012___,x_c_1_012___,x_c_2_012___] = ndgrid(x_c_0_,x_c_1_,x_c_2_); n_xxx_c = n_x_c^3; xxx_c_weight = (2*x_p_r_max/n_x_c)^3;
k_c_0_ = periodize(0:n_x_c-1,-n_x_c/2,+n_x_c/2)/2; %<-- box has diameter 2. ;
k_c_1_ = periodize(0:n_x_c-1,-n_x_c/2,+n_x_c/2)/2; %<-- box has diameter 2. ;
k_c_2_ = periodize(0:n_x_c-1,-n_x_c/2,+n_x_c/2)/2; %<-- box has diameter 2. ;
dk = 1.0/diameter_x_c;
[k_c_0_012___,k_c_1_012___,k_c_2_012___] = ndgrid(k_c_0_,k_c_1_,k_c_2_); n_kkk_c = n_x_c^3;
%%%%%%%%;
tmp_t = tic;
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
,n_polar_a_k_ ...
,polar_a_ka__ ...
,n_azimu_b_ka__ ...
] = ...
sample_sphere_7( ...
 flag_verbose ...
,k_p_r_max ...
,k_eq_d ...
,'L' ...
,flag_unif_vs_adap ...
,flag_tensor_vs_adap ...
) ; %<-- sum(weight_3d_k_p_r_)*(4*pi) = (4/3)*pi*k_p_r_max^3 --> sum(weight_3d_k_p_r_) = (1/3)*k_p_r_max^3 ;
tmp_t = toc(tmp_t); disp(sprintf(' %% sample_sphere_7: time %0.2fs',tmp_t));
%%%%%%%%;
n_w_max = 2*(n_k_p_r);
template_k_eq_d = -1;%template_k_eq_d = 1.0/(2*pi);
n_w_0in_ = n_w_max*ones(n_k_p_r,1);
[ ...
 n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_k_all_ ...
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
n_w_max = max(n_w_); n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_]);
%%%%%%%%;

sigma_k_p = 1/sigma_x_c;
%%%%%%%%;
% check formula for fourier-transform of gaussian. ;
%%%%%%%%;
flag_check=0;
if flag_check;
b_x_c_0_ = 1/(sqrt(2*pi)*sigma_x_c) * exp(-(x_c_0_).^2/(2*sigma_x_c^2));
b_k_c_0_ = exp(-(2*pi*k_c_0_).^2/(2*sigma_k_p^2));
b_x_c_l2 = sum(abs(b_x_c_0_).^2,'all')*dx;
disp(sprintf(' %% b_x_c_l2 = %0.16f',b_x_c_l2));
b_k_c_l2 = sum(abs(b_k_c_0_).^2,'all')*dk;
disp(sprintf(' %% b_k_c_l2 = %0.16f',b_k_c_l2));
b_x_c_form_ = 1/(sqrt(2*pi)*sigma_x_c)^3 * exp( -( ...
						    + (x_c_0_012___).^2 ...
						    + (x_c_1_012___).^2 ...
						    + (x_c_2_012___).^2 ...
						    ) / (2*sigma_x_c^2) );
b_k_p_form_ = exp( -(2*pi*k_p_r_all_).^2 / (2*sigma_k_p^2) ) ;
b_x_c_l2 = sum(abs(b_x_c_form_).^2,'all')*dx^3;
disp(sprintf(' %% b_x_c_l2 = %0.16f',b_x_c_l2));
b_k_p_l2 = sum(abs(b_k_p_form_).^2.*weight_3d_k_all_,'all');
disp(sprintf(' %% b_k_p_l2 = %0.16f',b_k_p_l2));
end;%if flag_check;

%%%%%%%%;
% Sum of gaussians. ;
%%%%%%%%;
if d_source==0;
n_source = 2; x_c_0_source_ = [-0.2;+0.4]; x_c_1_source_ = [-0.7;+0.3]; x_c_2_source_ = [-0.4;+0.6];
end;%if d_source==0;
if d_source> 0;
[ ...
 n_source ...
,azimu_b_source_ ...
,polar_a_source_ ...
,weight_source_ ...
,x_c_0_source_ ...
,x_c_1_source_ ...
,x_c_2_source_ ...
] = ...
sample_shell_5( ...
 0.75*half_diameter_x_c ...
,0.75*half_diameter_x_c*d_source ...
) ;
end;%if d_source> 0;
%%%%%%%%;
a_x_c_form_ = zeros(n_x_c,n_x_c,n_x_c);
a_k_p_form_ = zeros(n_k_all,1);
for nsource=0:n_source-1;
b_x_c_form_ = 1/(sqrt(2*pi)*sigma_x_c)^3 * exp( -( ...
						    + (x_c_0_012___ - x_c_0_source_(1+nsource)).^2 ...
						    + (x_c_1_012___ - x_c_1_source_(1+nsource)).^2 ...
						    + (x_c_2_012___ - x_c_2_source_(1+nsource)).^2 ...
						    ) / (2*sigma_x_c^2) );
a_x_c_form_ = a_x_c_form_ + b_x_c_form_;
b_k_p_form_ = exp( -(2*pi*k_p_r_all_).^2 / (2*sigma_k_p^2) ) .* ...
  exp(-i*2*pi*( ...
		+k_c_0_all_*x_c_0_source_(1+nsource) ...
		+k_c_1_all_*x_c_1_source_(1+nsource) ...
		+k_c_2_all_*x_c_2_source_(1+nsource) ...
		));
a_k_p_form_ = a_k_p_form_ + b_k_p_form_;
clear b_x_c_form_ b_k_p_form_;
end;%for nsource=0:n_source-1;
%%%%%%%%;
a_x_c_l2 = sum(abs(a_x_c_form_).^2,'all')*dx^3;
disp(sprintf(' %% a_x_c_l2 = %0.16f',a_x_c_l2));
a_k_p_l2 = sum(abs(a_k_p_form_).^2.*weight_3d_k_all_,'all');
disp(sprintf(' %% a_k_p_l2 = %0.16f',a_k_p_l2));
%%%%%%%%;
flag_disp=0;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;
isosurface_f_x_c_0(a_x_c_form_,98.5);
end;%if flag_disp;
%%%%%%%%;

%%%%%%%%;
% Now convert to a_k_p_ ;
%%%%%%%%;
%%%%%%%%;
eta = pi/x_p_r_max; tmp_t = tic;
a_k_p_quad_ = xxnufft3d3(n_xxx_c,x_c_0_012___(:)*eta,x_c_1_012___(:)*eta,x_c_2_012___(:)*eta,a_x_c_form_(:).*xxx_c_weight,-1,1e-12,n_k_all,2*pi*k_c_0_all_/eta,2*pi*k_c_1_all_/eta,2*pi*k_c_2_all_/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi) *  sqrt(2*pi)^3;
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_k_p_quad_ time %0.2fs',tmp_t));
disp(sprintf(' %% xxnufft3d3: a_k_p_l2 ratio: (sum(abs(a_k_p_form_).^2.*weight_3d_k_all_)/sum(abs(a_k_p_quad_).^2.*weight_3d_k_all_)): %0.16f',(sum(abs(a_k_p_form_).^2.*weight_3d_k_all_)/sum(abs(a_k_p_quad_).^2.*weight_3d_k_all_))));
disp(sprintf(' %% xxnufft3d3: a_k_p_quad error: %0.16f',fnorm(a_k_p_form_-a_k_p_quad_)/fnorm(a_k_p_form_)));
%%%%%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
a_x_c_reco_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,a_k_p_quad_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_xxx_c,x_c_0_012___(:)/eta,x_c_1_012___(:)/eta,x_c_2_012___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi) /  sqrt(2*pi)^3;
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_x_c_reco_ time %0.2fs',tmp_t));
disp(sprintf(' %% xxnufft3d3: a_x_c_reco error: %0.16f',fnorm(a_x_c_form_(:)-a_x_c_reco_)/fnorm(a_x_c_form_(:))));
%%%%%%%%;
flag_disp=1;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,2,1); isosurface_f_x_u_1([],a_x_c_form_); title('a_x_c_form_','Interpreter','none');
subplot(1,2,2); isosurface_f_x_u_1([],a_x_c_reco_); title('a_x_c_reco_','Interpreter','none');
end;%if flag_disp;
%%%%%%%%;

%%%%%%%%;
% Define deformation via longitudinal perturbation. ;
%%%%%%%%;
% Below we assume that the integral of a function on a ring can be approximated by the laplacian: ;
% B_{r} = ring of radius r. ;
% I(r) = \frac{1}{|\partial B_{r}|} \int_{B_{r}} u(\vx) d\Gamma. ;
% \partial_{r} I(r) = \frac{1}{|\partial B_{r}|} \int_{B_{r}} \Delta u(\vx) d\vx. ;
% \partial_{r} I(r) \approx \frac{1}{|\partial B_{r}|} \pi r^{2} \Delta u(0) \approx 0.5*r*\Delta u(0). ;
% I(r) = I(0) + \int_{0}^{r} \partial_{r} I(s)ds \approx u(0) + 0.25*r^{2}*\Delta u(0). ;
%%%%%%%%;
% Similarly, the integral of a function on a collapsed ring (e.g., only in the azimu_b direction) is: ;
% B_{r} = ring of radius r projected onto the x_{0} coordinate. ;
% I(r) = \frac{1}{|\partial B_{r}|} \int_{B_{r}} u(x_{0}) d\Gamma \approx u(0) + 0.25*r^{2}*\partial^{2}_{x_{0}} u(0). ;
% This is what one would expect if u did not vary orthogonal to the x_{0} direction. ;
%%%%%%%%;

%%%%%%%%;
% Set up interpolation matrix. ;
% This strongly assumes that each shell and ring is discretized identically. ;
% (i.e., flag_unif_vs_adap==1 & flag_tensor_vs_adap==1). ;
%%%%%%%%;
nk_p_r = 0; %<-- any individual shell will do. ;
t_max = 1.0; n_t = 16;
dt = t_max/max(1,n_t);
g_dilation = @(polar_a) dt*sin(2*polar_a); %<-- first mode. ;
n_order = 5;
%%%%;
n_k_all_csum_pre = n_k_all_csum_(1+nk_p_r+0);
n_k_all_csum_pos = n_k_all_csum_(1+nk_p_r+1);
shell_n_k_all = n_k_all_csum_pos-n_k_all_csum_pre;
shell_n_polar_a = n_polar_a_k_(1+nk_p_r);
tmp_n_azimu_b_a_ = n_azimu_b_ka__{1+nk_p_r};
if (mean(diff(tmp_n_azimu_b_a_))> 0); disp(sprintf(' %% Warning, nonuniform tmp_n_azimu_b_a_')); end;
shell_n_azimu_b = tmp_n_azimu_b_a_(1+0);
assert(shell_n_k_all==shell_n_polar_a*shell_n_azimu_b);
shell_k_p_azimu_b_ = k_p_azimu_b_all_(1+n_k_all_csum_pre+[0:shell_n_k_all-1]);
shell_k_p_polar_a_ = k_p_polar_a_all_(1+n_k_all_csum_pre+[0:shell_n_k_all-1]);
shell_longitudinal_r_ = flag_sign*0.5*g_dilation(shell_k_p_polar_a_);
%%%%;
tmp_t = tic();
shell_n_scatter = shell_n_k_all; %<-- single point for each. ;
shell_azimu_b_scatter__ = shell_k_p_azimu_b_;
shell_polar_a_scatter_ = shell_k_p_polar_a_ + shell_longitudinal_r_;
[ ...
 shell_scatter_from_tensor_sba__ ...
,shell_scatter_from_tensor_db1da0_sba__ ...
,shell_scatter_from_tensor_db0da1_sba__ ...
,shell_scatter_from_tensor_db1da1_sba__ ...
,shell_scatter_from_tensor_db2da0_sba__ ...
,shell_scatter_from_tensor_db0da2_sba__ ...
] = ...
cg_interpolate_n_3( ...
 n_order ...
,shell_n_azimu_b ...
,shell_n_polar_a ...
,shell_n_scatter ...
,shell_azimu_b_scatter__ ...
,shell_polar_a_scatter_ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% cg_interpolate_n_3: time %0.2fs',tmp_t));
%%%%%%%%;

%%%%%%%%;
% Define slices to visualize during deformation. ;
%%%%%%%%;
n_viewing_sub = 1;
template_viewing_sub_azimu_b_ori_ = pi*[linspace(-1.0/3.0,-1.0/3.0,n_viewing_sub)];
template_viewing_sub_polar_a_ori_ = pi*[linspace(0,+0.5,n_viewing_sub+2)];
template_viewing_sub_polar_a_ori_ = template_viewing_sub_polar_a_ori_(2:end-1);
template_viewing_sub_azimu_b_ = template_viewing_sub_azimu_b_ori_;
template_viewing_sub_polar_a_ = template_viewing_sub_polar_a_ori_;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
flag_disp=1;
if flag_disp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
S_k_p_wkn__ = zeros(n_w_sum,n_viewing_sub);
%%%%%%%%;
for nviewing_sub=0:n_viewing_sub-1;
template_viewing_sub_azimu_b = template_viewing_sub_azimu_b_(1+nviewing_sub);
template_viewing_sub_polar_a = template_viewing_sub_polar_a_(1+nviewing_sub);
[ ...
 template_ring_k_c_0__ ...
,template_ring_k_c_1__ ...
,template_ring_k_c_2__ ...
,template_ring_azimu_b__ ...
,template_ring_polar_a__ ...
] = ...
get_template_single_ring_k_c_0( ...
 flag_verbose ...
,1 ...
,template_viewing_sub_azimu_b ...
,template_viewing_sub_polar_a ...
,n_w_max ...
);
%%%%%%%%;
n_order = 5;
ring_azimu_b_scatter_ = template_ring_azimu_b__(:,1+0);
ring_polar_a_scatter_ = template_ring_polar_a__(:,1+0);
[ ...
 ring_scatter_from_tensor_sba__ ...
] = ...
cg_interpolate_n_3( ...
 n_order ...
,shell_n_azimu_b ...
,shell_n_polar_a ...
,n_w_max ...
,ring_azimu_b_scatter_ ...
,ring_polar_a_scatter_ ...
);
%%%%%%%%;
S_k_p_wk_ = zeros(n_w_sum,1);
for nk_p_r=0:n_k_p_r-1;
n_k_all_csum_pre = n_k_all_csum_(1+nk_p_r+0);
n_k_all_csum_pos = n_k_all_csum_(1+nk_p_r+1);
assert(shell_n_k_all==n_k_all_csum_pos-n_k_all_csum_pre);
a_k_p_sub_ = a_k_p_quad_(1+n_k_all_csum_pre+[0:shell_n_k_all-1]);
S_k_p_sub_ = ring_scatter_from_tensor_sba__*a_k_p_sub_;
S_k_p_wk_(1+n_w_csum_(1+nk_p_r)+[0:n_w_max-1]) = S_k_p_sub_;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
S_k_p_wkn__(:,1+nviewing_sub) = S_k_p_wk_;
clear template_ring* ring_azimu_b_scatter_ ring_polar_a_scatter_ ring_scatter_from_tensor_sba__ S_k_p_wk_;
end;%for nviewing_sub=0:n_viewing_sub-1;
%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figmed;
%%%%%%%%;
subplot(1,2,1);
tmp_parameter = struct('type','parameter');
tmp_parameter.c_use__ = colormap_beach();
tmp_parameter.vlim_ = prctile(mean(reshape(a_x_c_form_,[n_x_c,n_x_c,n_x_c]),3),[ 5,95]);
isosurface_f_x_u_1(tmp_parameter,a_x_c_form_);
xlabel('x0'); ylabel('x1'); zlabel('x2'); axisnotick3d;
%%%%%%%%;
subplot(1,2,2);
%%%%;
nk_p_r = floor(0.50*n_k_p_r);
k_p_r = k_p_r_(1+nk_p_r);
n_k_all_csum_pre = n_k_all_csum_(1+nk_p_r+0);
n_k_all_csum_pos = n_k_all_csum_(1+nk_p_r+1);
assert(shell_n_k_all==n_k_all_csum_pos-n_k_all_csum_pre);
a_k_p_sub_fin_ = a_k_p_quad_(1+n_k_all_csum_pre+[0:shell_n_k_all-1]);
a_k_p_sub_fin_lim_ = prctile(abs(a_k_p_sub_fin_),95)*[-1,+1];
%%%%;
nk_p_r = floor(0.33*n_k_p_r);%nk_p_r = n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
n_k_all_csum_pre = n_k_all_csum_(1+nk_p_r+0);
n_k_all_csum_pos = n_k_all_csum_(1+nk_p_r+1);
assert(shell_n_k_all==n_k_all_csum_pos-n_k_all_csum_pre);
a_k_p_sub_ = a_k_p_quad_(1+n_k_all_csum_pre+[0:shell_n_k_all-1]);
a_k_p_sub_lim_ = prctile(abs(a_k_p_sub_),95)*[-1,+1];
%%%%;
hold on;
imagesc_polar_a_azimu_b_0(shell_k_p_polar_a_,shell_k_p_azimu_b_,real(a_k_p_sub_),a_k_p_sub_fin_lim_,colormap_80s,0,k_p_r);
n_contour = 16;
tmp_parameter = struct('type','parameter');
tmp_parameter.flag_percent_use=0;
tmp_parameter.vlim_ = 1.0*a_k_p_sub_fin_lim_;
tmp_parameter.vval_ = transpose(linspace(min(tmp_parameter.vlim_),max(tmp_parameter.vlim_),n_contour));
imagesc_S_k_p_3d_2( ...
 tmp_parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,weight_2d_k_all_ ...
,n_viewing_sub ...
,real(S_k_p_wkn__) ...
,template_viewing_sub_azimu_b_ ...
,template_viewing_sub_polar_a_ ...
);
hold off;
xlabel('k0'); ylabel('k1'); zlabel('k2'); axisnotick3d;
%%%%%%%%;
drawnow;
str_d_source = sprintf('ds%.2d',floor(100*d_source));
str_sigma_x_c = sprintf('sg%.2d',floor(100*sigma_x_c));
if flag_sign> 0; str_flag_sign = sprintf('dp'); end;
if flag_sign< 0; str_flag_sign = sprintf('dn'); end;
str_nt = sprintf('nt00');
fname_fig_pre = sprintf('%s/M3d_shape_longitudinal_perturbation_%s%s%s%s_FIGA',dir_manuscript_jpg,str_d_source,str_sigma_x_c,str_flag_sign,str_nt);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if ( flag_replot | ~exist(fname_fig_jpg,'file'));
if (flag_verbose> 0); disp(sprintf(' %% writing %s',fname_fig_pre)); end;
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
end;%if ( flag_replot | ~exist(fname_fig_jpg,'file'));
%%%%%%%%;
clear S_k_p_wkn__;
%close(gcf);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_disp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% Now perform deformation. ;
%%%%%%%%;
a_k_p_pre_ = a_k_p_form_;
for nt=0:5;%for nt=0:n_t-1;
tmp_t = tic;
a_k_p_pos_ = a_k_p_pre_;
for nk_p_r=0:n_k_p_r-1;
n_k_all_csum_pre = n_k_all_csum_(1+nk_p_r+0);
n_k_all_csum_pos = n_k_all_csum_(1+nk_p_r+1);
assert(shell_n_k_all==n_k_all_csum_pos-n_k_all_csum_pre);
a_k_p_sub_ = a_k_p_pos_(1+n_k_all_csum_pre+[0:shell_n_k_all-1]);
a_k_p_sub_ = shell_scatter_from_tensor_sba__*a_k_p_sub_ + 0.25 * shell_longitudinal_r_.^2 .* (shell_scatter_from_tensor_db2da0_sba__*a_k_p_sub_ + shell_scatter_from_tensor_db0da2_sba__*a_k_p_sub_);
a_k_p_pos_(1+n_k_all_csum_pre+[0:shell_n_k_all-1]) = a_k_p_sub_;
end;%for nk_p_r=0:n_k_p_r-1;
a_k_p_pre_ = a_k_p_pos_;
tmp_t = toc(tmp_t); disp(sprintf(' %% shell_scatter_from_tensor_sba__: time %0.2fs',tmp_t));
%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
a_x_c_pos_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,a_k_p_pos_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_xxx_c,x_c_0_012___(:)/eta,x_c_1_012___(:)/eta,x_c_2_012___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi) /  sqrt(2*pi)^3;
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_x_c_pos_ time %0.2fs',tmp_t));
%%%%;
template_viewing_sub_polar_a_ = template_viewing_sub_polar_a_ + flag_sign*0.5*g_dilation(template_viewing_sub_polar_a_);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
flag_disp=1;
if flag_disp;
%%%%%%%%;
S_k_p_pos_wkn__ = zeros(n_w_sum,n_viewing_sub);
%%%%%%%%;
for nviewing_sub=0:n_viewing_sub-1;
template_viewing_sub_azimu_b = template_viewing_sub_azimu_b_(1+nviewing_sub);
template_viewing_sub_polar_a = template_viewing_sub_polar_a_(1+nviewing_sub);
[ ...
 template_ring_k_c_0__ ...
,template_ring_k_c_1__ ...
,template_ring_k_c_2__ ...
,template_ring_azimu_b__ ...
,template_ring_polar_a__ ...
] = ...
get_template_single_ring_k_c_0( ...
 flag_verbose ...
,1 ...
,template_viewing_sub_azimu_b ...
,template_viewing_sub_polar_a ...
,n_w_max ...
);
%%%%%%%%;
n_order = 5;
ring_azimu_b_scatter_ = template_ring_azimu_b__(:,1+0);
ring_polar_a_scatter_ = template_ring_polar_a__(:,1+0);
[ ...
 ring_scatter_from_tensor_sba__ ...
] = ...
cg_interpolate_n_3( ...
 n_order ...
,shell_n_azimu_b ...
,shell_n_polar_a ...
,n_w_max ...
,ring_azimu_b_scatter_ ...
,ring_polar_a_scatter_ ...
);
%%%%%%%%;
S_k_p_pos_wk_ = zeros(n_w_sum,1);
for nk_p_r=0:n_k_p_r-1;
n_k_all_csum_pre = n_k_all_csum_(1+nk_p_r+0);
n_k_all_csum_pos = n_k_all_csum_(1+nk_p_r+1);
assert(shell_n_k_all==n_k_all_csum_pos-n_k_all_csum_pre);
a_k_p_pos_sub_ = a_k_p_pos_(1+n_k_all_csum_pre+[0:shell_n_k_all-1]);
S_k_p_pos_sub_ = ring_scatter_from_tensor_sba__*a_k_p_pos_sub_;
S_k_p_pos_wk_(1+n_w_csum_(1+nk_p_r)+[0:n_w_max-1]) = S_k_p_pos_sub_;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
S_k_p_pos_wkn__(:,1+nviewing_sub) = S_k_p_pos_wk_;
clear template_ring* ring_azimu_b_scatter_ ring_polar_a_scatter_ ring_scatter_from_tensor_sba__ S_k_p_pos_wk_;
end;%for nviewing_sub=0:n_viewing_sub-1;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figmed;
%%%%%%%%;
subplot(1,2,1);
tmp_parameter = struct('type','parameter');
tmp_parameter.c_use__ = colormap_beach();
tmp_parameter.vlim_ = prctile(mean(reshape(real(a_x_c_pos_),[n_x_c,n_x_c,n_x_c]),3),[ 5,95]);
isosurface_f_x_u_1(tmp_parameter,real(a_x_c_pos_));
xlabel('x0'); ylabel('x1'); zlabel('x2'); axisnotick3d;
%%%%%%%%;
subplot(1,2,2);
%%%%;
nk_p_r = floor(0.50*n_k_p_r);
k_p_r = k_p_r_(1+nk_p_r);
n_k_all_csum_pre = n_k_all_csum_(1+nk_p_r+0);
n_k_all_csum_pos = n_k_all_csum_(1+nk_p_r+1);
assert(shell_n_k_all==n_k_all_csum_pos-n_k_all_csum_pre);
a_k_p_pos_sub_fin_ = a_k_p_pos_(1+n_k_all_csum_pre+[0:shell_n_k_all-1]);
a_k_p_pos_sub_fin_lim_ = prctile(abs(a_k_p_pos_sub_fin_),95)*[-1,+1];
%%%%;
nk_p_r = floor(0.33*n_k_p_r);%nk_p_r = n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
n_k_all_csum_pre = n_k_all_csum_(1+nk_p_r+0);
n_k_all_csum_pos = n_k_all_csum_(1+nk_p_r+1);
assert(shell_n_k_all==n_k_all_csum_pos-n_k_all_csum_pre);
a_k_p_pos_sub_ = a_k_p_pos_(1+n_k_all_csum_pre+[0:shell_n_k_all-1]);
a_k_p_pos_sub_lim_ = prctile(abs(a_k_p_pos_sub_),95)*[-1,+1];
%%%%;
hold on;
imagesc_polar_a_azimu_b_0(shell_k_p_polar_a_,shell_k_p_azimu_b_,real(a_k_p_pos_sub_),a_k_p_pos_sub_fin_lim_,colormap_80s,0,k_p_r);
n_contour = 16;
tmp_parameter = struct('type','parameter');
tmp_parameter.flag_percent_use=0;
tmp_parameter.vlim_ = 1.0*a_k_p_pos_sub_fin_lim_;
tmp_parameter.vval_ = transpose(linspace(min(tmp_parameter.vlim_),max(tmp_parameter.vlim_),n_contour));
imagesc_S_k_p_3d_2( ...
 tmp_parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,weight_2d_k_all_ ...
,n_viewing_sub ...
,real(S_k_p_pos_wkn__) ...
,template_viewing_sub_azimu_b_ ...
,template_viewing_sub_polar_a_ ...
);
hold off;
xlabel('k0'); ylabel('k1'); zlabel('k2'); axisnotick3d;
%%%%%%%%;
drawnow;
str_d_source = sprintf('ds%.2d',floor(100*d_source));
str_sigma_x_c = sprintf('sg%.2d',floor(100*sigma_x_c));
str_flag_sign = sprintf('d%.1d',flag_sign>0);
str_nt = sprintf('nt%.2d',1+nt);
fname_fig_pre = sprintf('%s/M3d_shape_longitudinal_perturbation_%s%s%s%s_FIGA',dir_manuscript_jpg,str_d_source,str_sigma_x_c,str_flag_sign,str_nt);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if ( flag_replot | ~exist(fname_fig_jpg,'file'));
if (flag_verbose> 0); disp(sprintf(' %% writing %s',fname_fig_pre)); end;
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
end;%if ( flag_replot | ~exist(fname_fig_jpg,'file'));
%%%%%%%%;
clear S_k_p_pos_wkn__;
close(gcf);
end;%if flag_disp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for nt=0:n_t-1;

if flag_verbose; disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
