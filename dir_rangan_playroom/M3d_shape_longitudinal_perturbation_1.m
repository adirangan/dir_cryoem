function ...
[ ...
 parameter ...
] = ...
M3d_shape_longitudinal_perturbation_1( ...
 parameter ...
);

if (nargin<1);
disp(sprintf(' %% testing %s',str_thisfunction));
parameter = struct('type','parameter');
parameter.flag_verbose = 1;
parameter.str_shape = 'rand';
M3d_shape_longitudinal_perturbation_1(parameter);
disp(sprintf(' %% returning')); return;
end;%if (nargin<1);

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;

%clear; parameter = [];
str_thisfunction = 'M3d_shape_longitudinal_perturbation_1';

if isempty(parameter); parameter = struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
flag_verbose = parameter.flag_verbose;
if ~isfield(parameter,'flag_disp'); parameter.flag_disp = 0; end;
flag_disp = parameter.flag_disp; nf=0;
if ~isfield(parameter,'flag_replot'); parameter.flag_replot = 0; end;
flag_replot = parameter.flag_replot;
if ~isfield(parameter,'str_shape'); parameter.str_shape = 'rand'; end;
str_shape = parameter.str_shape;
if ~isfield(parameter,'n_mode'); parameter.n_mode = 2; end;
n_mode = parameter.n_mode;
if ~isfield(parameter,'d_source'); parameter.d_source = pi/2; end;
d_source = parameter.d_source;
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

%%%%%%%%;
% check formula for fourier-transform of gaussian. ;
%%%%%%%%;
flag_check=0;
if flag_check;
sigma_x_c = 0.05;
sigma_k_p = 1/sigma_x_c;
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
%%%%%%%%;
sigma_x_c = 0.05;
sigma_k_p = 1/sigma_x_c;
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
subplot(1,2,1); isosurface_f_x_u_0(a_x_c_form_,98.5); title('a_x_c_form_','Interpreter','none');
subplot(1,2,2); isosurface_f_x_u_0(a_x_c_reco_,98.5); title('a_x_c_reco_','Interpreter','none');
end;%if flag_disp;
%%%%%%%%;

%%%%%%%%;
% Now deform via longitudinal perturbation. ;
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
shell_longitudinal_r_ = 0.5*g_dilation(shell_k_p_polar_a_);
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
a_k_p_pre_ = a_k_p_form_;
for nt=0:n_t-1;
tmp_t = tic;
a_k_p_pos_ = a_k_p_pre_;
for nk_p_r=0:n_k_p_r-1;
n_k_all_csum_pre = n_k_all_csum_(1+nk_p_r+0);
n_k_all_csum_pos = n_k_all_csum_(1+nk_p_r+1);
shell_n_k_all = n_k_all_csum_pos-n_k_all_csum_pre;
a_k_p_sub_ = a_k_p_pos_(1+n_k_all_csum_pre+[0:shell_n_k_all-1]);
a_k_p_sub_ = shell_scatter_from_tensor_sba__*a_k_p_sub_ + 0.25 * shell_longitudinal_r_.^2 .* (shell_scatter_from_tensor_db2da0_sba__*a_k_p_sub_ + shell_scatter_from_tensor_db0da2_sba__*a_k_p_sub_);
a_k_p_pos_(1+n_k_all_csum_pre+[0:shell_n_k_all-1]) = a_k_p_sub_;
end;%for nk_p_r=0:n_k_p_r-1;
a_k_p_pre_ = a_k_p_pos_;
tmp_t = toc(tmp_t); disp(sprintf(' %% shell_scatter_from_tensor_sba__: time %0.2fs',tmp_t));
%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
a_x_c_pre_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,a_k_p_pre_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_xxx_c,x_c_0_012___(:)/eta,x_c_1_012___(:)/eta,x_c_2_012___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi) /  sqrt(2*pi)^3;
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_x_c_pre_ time %0.2fs',tmp_t));
%%%%;
flag_disp=1;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,2,1); isosurface_f_x_u_0(a_x_c_pre_,98.5); title(sprintf('a_x_c_pre_ nt %d/%d',nt,n_t),'Interpreter','none');
nk_p_r = floor(n_k_p_r/3);%nk_p_r = n_k_p_r-1;
n_k_all_csum_pre = n_k_all_csum_(1+nk_p_r+0);
n_k_all_csum_pos = n_k_all_csum_(1+nk_p_r+1);
shell_n_k_all = n_k_all_csum_pos-n_k_all_csum_pre;
a_k_p_sub_ = a_k_p_pos_(1+n_k_all_csum_pre+[0:shell_n_k_all-1]);
a_k_p_sub_lim_ = prctile(abs(a_k_p_sub_),95)*[-1,+1];
subplot(1,2,2); imagesc_polar_a_azimu_b_0(shell_k_p_polar_a_,shell_k_p_azimu_b_,real(a_k_p_sub_),a_k_p_sub_lim_,colormap_80s,0); axisnotick;
drawnow();
end;%if flag_disp;
%%%%;
end;%for nt=0:n_t-1;

if flag_verbose; disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
