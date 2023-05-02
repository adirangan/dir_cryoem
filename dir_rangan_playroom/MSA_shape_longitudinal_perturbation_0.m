function ...
[ ...
 parameter ...
] = ...
MSA_shape_longitudinal_perturbation_0( ...
 parameter ...
);

str_thisfunction = 'MSA_shape_longitudinal_perturbation_0';

if (nargin<1);
disp(sprintf(' %% testing %s',str_thisfunction));
%%%%%%%%;
for n_source_gaussian=[1,8];
for n_mode=[2,4];
parameter = struct('type','parameter');
parameter.flag_verbose = 1;
parameter.flag_replot = 1;
parameter.str_shape = 'rand';
parameter.n_mode = n_mode;
parameter.n_source_gaussian = n_source_gaussian;
MSA_shape_longitudinal_perturbation_0(parameter);
end;%for n_mode=[2,4];
end;%for n_source_gaussian=[1,8];
%%%%%%%%;
disp(sprintf(' %% returning')); return;
end;%if (nargin<1);

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;

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
if ~isfield(parameter,'n_source_gaussian'); parameter.n_source_gaussian = 1; end;
n_source_gaussian = parameter.n_source_gaussian;
if ~isfield(parameter,'rseed'); parameter.rseed = 0; end;
rseed = parameter.rseed;
if ~isfield(parameter,'n_x_c'); parameter.n_x_c = 256; end;
n_x_c = parameter.n_x_c;
if ~isfield(parameter,'k_p_r_max'); parameter.k_p_r_max = 3.0*48.0/(2*pi); end;
k_p_r_max = parameter.k_p_r_max;
if ~isfield(parameter,'k_eq_d'); parameter.k_eq_d = 0.5/(2*pi); end;
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
dir_manuscript_jpg = sprintf('%s/dir_MSA_shape_longitudinal_perturbation_jpg',dir_manuscript);
if ~exist(dir_manuscript_jpg,'dir'); disp(sprintf(' %% mkdir %s',dir_manuscript_jpg)); mkdir(dir_manuscript_jpg); end;

if strcmp(str_shape,'rand'); rng(rseed); end;

%%%%%%%%;
half_diameter_x_c = 1.0d0;
diameter_x_c = 2.0d0*half_diameter_x_c;
x_p_r_max = 1.0;
x_c_0_ = -half_diameter_x_c + transpose([0:n_x_c-1]/n_x_c)*diameter_x_c;
x_c_1_ = -half_diameter_x_c + transpose([0:n_x_c-1]/n_x_c)*diameter_x_c;
dx = diameter_x_c/n_x_c;
[x_c_0_01__,x_c_1_01__] = ndgrid(x_c_0_,x_c_1_); n_xx_c = n_x_c^2; xx_c_weight = (2*x_p_r_max/n_x_c)^2;
k_c_0_ = periodize(0:n_x_c-1,-n_x_c/2,+n_x_c/2)/2; %<-- box has diameter 2. ;
k_c_1_ = periodize(0:n_x_c-1,-n_x_c/2,+n_x_c/2)/2; %<-- box has diameter 2. ;
[k_c_0_01__,k_c_1_01__] = ndgrid(k_c_0_,k_c_1_); n_kk_c = n_x_c^2;
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
] = ...
sample_sphere_7( ...
 flag_verbose ...
,k_p_r_max ...
,k_eq_d ...
,'L' ...
) ; %<-- sum(weight_3d_k_p_r_)*(4*pi) = (4/3)*pi*k_p_r_max^3 --> sum(weight_3d_k_p_r_) = (1/3)*k_p_r_max^3 ;
tmp_t = toc(tmp_t); disp(sprintf(' %% sample_sphere_7: time %0.2fs',tmp_t));
%%%%%%%%;
n_w_max = 2*(n_k_p_r);
n_w_0in_ = n_w_max*ones(n_k_p_r,1);
%{
[ ...
 ~ ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_k_all_ ...
] = ...
get_template_1( ...
 flag_verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,zeros(n_k_p_r,1) ...
,zeros(n_k_p_r,1) ...
,1 ...
,-1 ...
,n_w_0in_ ...
);
n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_]);
%}

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

k_c_0_all_ = zeros(n_w_sum,1);
k_c_1_all_ = zeros(n_w_sum,1);
na=0;
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
n_w = n_w_(1+nk_p_r);
for nw=0:n_w-1;
k_c_0_all = k_p_r*cos(2*pi*nw/n_w);
k_c_1_all = k_p_r*sin(2*pi*nw/n_w);
k_c_0_all_(1+na) = k_c_0_all;
k_c_1_all_(1+na) = k_c_1_all;
na=na+1;
end;%for nw=0:n_w-1;
end;%for nk_p_r=0:n_k_p_r-1;
assert(na==n_w_sum);
%%%%%%%%;
h2d_ = @(kd) 4*pi^2*(besselj(0,kd) + besselj(2,kd)); % calculates <f_j,f_k>, normalized so that <f,f> = (4*pi^2);
dh2d_ = @(kd) 4*pi^3*(besselj(-1,kd) - besselj(+3,kd));
h3d_ = @(kd) 4*pi*( sin(kd) - (kd).*cos(kd) ) ./ kd.^3 ; % calculates <f_j,f_k>, normalized so that <f,f> = 4*pi/3;
dh3d_ = @(kd) 12*pi*( (kd.^2/3 - 1) .* sin(kd) + (kd).*cos(kd) ) ./ kd.^4 ;
%%%%%%%%;
n_source_planewave = 4;
rng(0);
delta_a_c__ = zeros(2,n_source_planewave);
delta_b_c__ = zeros(2,n_source_planewave);
for nsource=0:n_source_planewave-1;
delta_a_c_ = 0.25*(2*rand(2,1)-1);
delta_a_c__(:,1+nsource) = delta_a_c_;
delta_b_c_ = 0.25*(2*rand(2,1)-1);
delta_b_c__(:,1+nsource) = delta_b_c_;
end;%for nsource=0:n_source_planewave-1;
a_k_p_form_ = zeros(n_w_sum,1);
b_k_p_form_ = zeros(n_w_sum,1);
for nsource=0:n_source_planewave-1;
delta_a_c_ = delta_a_c__(:,1+nsource);
a_k_p_form_ = a_k_p_form_ + exp(+i*2*pi*(k_c_0_all_*delta_a_c_(1+0) + k_c_1_all_*delta_a_c_(1+1)));
delta_b_c_ = delta_b_c__(:,1+nsource);
b_k_p_form_ = b_k_p_form_ + exp(+i*2*pi*(k_c_0_all_*delta_b_c_(1+0) + k_c_1_all_*delta_b_c_(1+1)));
end;%for nsource=0:n_source_planewave-1;
I_a_quad = sum(a_k_p_form_.*weight_2d_k_all_)*(2*pi)^2;
I_b_quad = sum(b_k_p_form_.*weight_2d_k_all_)*(2*pi)^2;
I_a_form = 0;
I_b_form = 0;
for nsource=0:n_source_planewave-1;
delta_a_c_ = delta_a_c__(:,1+nsource);
delta_b_c_ = delta_b_c__(:,1+nsource);
tmp_kd = 2*pi*k_p_r_max*fnorm(delta_a_c_);
I_a_form = I_a_form + h2d_(tmp_kd)*k_p_r_max^2/(4*pi);
tmp_kd = 2*pi*k_p_r_max*fnorm(delta_b_c_);
I_b_form = I_b_form + h2d_(tmp_kd)*k_p_r_max^2/(4*pi);
end;%for nsource=0:n_source_planewave-1;
disp(sprintf(' %% I_a_form vs I_a_quad %0.16f',fnorm(I_a_form-I_a_quad)/fnorm(I_a_form)));
disp(sprintf(' %% I_b_form vs I_b_quad %0.16f',fnorm(I_b_form-I_b_quad)/fnorm(I_b_form)));
%%%%%%%%;
a_k_p_form_l2 = sum(abs(a_k_p_form_).^2.*weight_2d_k_all_)*(2*pi)^2;
a_x_c_quad_ = interp_k_p_to_x_c_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,n_k_p_r,k_p_r_,n_w_,a_k_p_form_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_c^2) * n_w_sum;
a_x_c_quad_l2 = sum(abs(a_x_c_quad_).^2,'all')*dx^2;
b_k_p_form_l2 = sum(abs(b_k_p_form_).^2.*weight_2d_k_all_)*(2*pi)^2;
b_x_c_quad_ = interp_k_p_to_x_c_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,n_k_p_r,k_p_r_,n_w_,b_k_p_form_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_c^2) * n_w_sum;
b_x_c_quad_l2 = sum(abs(b_x_c_quad_).^2,'all')*dx^2;
disp(sprintf(' %% a_k_p_form_l2 vs a_x_c_quad_l2: %0.16f',fnorm(a_k_p_form_l2 - a_x_c_quad_l2)/fnorm(a_k_p_form_l2)));
disp(sprintf(' %% b_k_p_form_l2 vs b_x_c_quad_l2: %0.16f',fnorm(b_k_p_form_l2 - b_x_c_quad_l2)/fnorm(b_k_p_form_l2)));
%%%%%%%%;
a_k_p_form_l2 = sum(abs(a_k_p_form_).^2.*weight_2d_k_all_)*(2*pi)^2;
a_x_c_quad_ = interp_k_p_to_x_c_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,n_k_p_r,k_p_r_,n_w_,a_k_p_form_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_c^2) * n_w_sum;
a_x_c_quad_l2 = sum(abs(a_x_c_quad_).^2,'all')*dx^2;
b_k_p_form_l2 = sum(abs(b_k_p_form_).^2.*weight_2d_k_all_)*(2*pi)^2;
b_x_c_quad_ = interp_k_p_to_x_c_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,n_k_p_r,k_p_r_,n_w_,b_k_p_form_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_c^2) * n_w_sum;
b_x_c_quad_l2 = sum(abs(b_x_c_quad_).^2,'all')*dx^2;
disp(sprintf(' %% a_k_p_form_l2 vs a_x_c_quad_l2: %0.16f',fnorm(a_k_p_form_l2 - a_x_c_quad_l2)/fnorm(a_k_p_form_l2)));
disp(sprintf(' %% b_k_p_form_l2 vs b_x_c_quad_l2: %0.16f',fnorm(b_k_p_form_l2 - b_x_c_quad_l2)/fnorm(b_k_p_form_l2)));
%%%%%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figbig;
subplot(2,2,1); imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(a_x_c_quad_)); axis image; axisnotick;
subplot(2,2,2); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(a_k_p_form_)); axis image; axisnotick;
subplot(2,2,3); imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(b_x_c_quad_)); axis image; axisnotick;
subplot(2,2,4); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(b_k_p_form_)); axis image; axisnotick;
end;%if flag_disp;
%%%%%%%%;
a_x_c_form_ = zeros(n_x_c,n_x_c);
b_x_c_form_ = zeros(n_x_c,n_x_c);
for nsource=0:n_source_planewave-1;
delta_a_c_ = delta_a_c__(:,1+nsource);
delta_b_c_ = delta_b_c__(:,1+nsource);
tmp_kd_01__ = 2*pi * k_p_r_max * sqrt( (x_c_0_01__ + delta_a_c_(1+0)).^2 + (x_c_1_01__ + delta_a_c_(1+1)).^2 ) ;
a_x_c_form_ = a_x_c_form_ + h2d_(tmp_kd_01__)*k_p_r_max^2/(4*pi);
tmp_kd_01__ = 2*pi * k_p_r_max * sqrt( (x_c_0_01__ + delta_b_c_(1+0)).^2 + (x_c_1_01__ + delta_b_c_(1+1)).^2 ) ;
b_x_c_form_ = b_x_c_form_ + h2d_(tmp_kd_01__)*k_p_r_max^2/(4*pi);
end;%for nsource=0:n_source_planewave-1;
disp(sprintf(' %% a_x_c_form_(:) vs a_x_c_quad_(:): %0.16f',fnorm(a_x_c_form_(:) - a_x_c_quad_(:))/fnorm(a_x_c_form_(:))));
disp(sprintf(' %% b_x_c_form_(:) vs b_x_c_quad_(:): %0.16f',fnorm(b_x_c_form_(:) - b_x_c_quad_(:))/fnorm(b_x_c_form_(:))));
a_x_c_form_l2 = sum(abs(a_x_c_form_(:)).^2)*dx^2;
b_x_c_form_l2 = sum(abs(b_x_c_form_(:)).^2)*dx^2;
a_x_c_quad_l2 = sum(abs(a_x_c_quad_(:)).^2)*dx^2;
b_x_c_quad_l2 = sum(abs(b_x_c_quad_(:)).^2)*dx^2;
disp(sprintf(' %% a_x_c_form_l2 vs a_x_c_quad_l2: %0.16f',fnorm(a_x_c_form_l2 - a_x_c_quad_l2)/fnorm(a_x_c_form_l2)));
disp(sprintf(' %% b_x_c_form_l2 vs b_x_c_quad_l2: %0.16f',fnorm(b_x_c_form_l2 - b_x_c_quad_l2)/fnorm(b_x_c_form_l2)));
%%%%%%%%;
a_k_p_quad_ = interp_x_c_to_k_p_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,a_x_c_form_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_c^2)*dx^2;
a_k_p_quad_l2 = sum(abs(a_k_p_quad_).^2.*weight_2d_k_all_)*(2*pi)^2;
b_k_p_quad_ = interp_x_c_to_k_p_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,b_x_c_form_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_c^2)*dx^2;
b_k_p_quad_l2 = sum(abs(b_k_p_quad_).^2.*weight_2d_k_all_)*(2*pi)^2;
disp(sprintf(' %% a_x_c_form_l2 vs a_k_p_quad_l2: %0.16f',fnorm(a_x_c_form_l2 - a_k_p_quad_l2)/fnorm(a_x_c_form_l2)));
disp(sprintf(' %% b_x_c_form_l2 vs b_k_p_quad_l2: %0.16f',fnorm(b_x_c_form_l2 - b_k_p_quad_l2)/fnorm(b_x_c_form_l2)));
disp(sprintf(' %% a_k_p_form_l2 vs a_k_p_quad_l2: %0.16f',fnorm(a_k_p_form_l2 - a_k_p_quad_l2)/fnorm(a_k_p_form_l2)));
disp(sprintf(' %% b_k_p_form_l2 vs b_k_p_quad_l2: %0.16f',fnorm(b_k_p_form_l2 - b_k_p_quad_l2)/fnorm(b_k_p_form_l2)));
%%%%%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figbig;
subplot(2,2,1); imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(a_x_c_form_)); axis image; axisnotick;
subplot(2,2,2); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(a_k_p_quad_)); axis image; axisnotick;
subplot(2,2,3); imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(b_x_c_form_)); axis image; axisnotick;
subplot(2,2,4); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(b_k_p_quad_)); axis image; axisnotick;
end;%if flag_disp;
%%%%%%%%;
% Sum of gaussians. ;
%%%%%%%%;
%n_source_gaussian = 8;
tmp_sigma_x_c = 0.025;
tmp_sigma_k_p = 1/tmp_sigma_x_c;
tmp_delta_2s__ = 1.25/sqrt(2)*[cos(linspace(0,2*pi,1+n_source_gaussian)) ; sin(linspace(0,2*pi,1+n_source_gaussian))]; tmp_delta_2s__ = tmp_delta_2s__(:,1:n_source_gaussian);
tmp_N_x_c_ = zeros(n_x_c,n_x_c);
for nsource_gaussian=0:n_source_gaussian-1;
tmp_delta_ = tmp_delta_2s__(:,1+nsource_gaussian);
tmp_M_x_c_ = 1/(sqrt(2*pi)*tmp_sigma_x_c)^2 * exp( -( (x_c_0_01__-tmp_delta_(1+0)).^2 + (x_c_1_01__-tmp_delta_(1+1)).^2 ) / (2*tmp_sigma_x_c^2) );
tmp_N_x_c_ = tmp_N_x_c_ + tmp_M_x_c_;
clear tmp_delta_;
end;%for nsource_gaussian=0:n_source_gaussian-1;
tmp_N_x_c_l2 = sum(tmp_N_x_c_.^2,'all')*dx^2;
disp(sprintf(' %% sum(tmp_N_x_c_*dx^2,''all'') = %0.16f',sum(tmp_N_x_c_*dx^2,'all')));
disp(sprintf(' %% tmp_N_x_c_l2 = %0.16f',tmp_N_x_c_l2));
tmp_N_k_p_ = interp_x_c_to_k_p_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,tmp_N_x_c_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_c^2)*dx^2;
tmp_N_k_p_l2 = sum(abs(tmp_N_k_p_).^2 .* weight_2d_k_all_) * (2*pi)^2;
disp(sprintf(' %% tmp_N_k_p_l2 = %0.16f',tmp_N_k_p_l2));
tmp_N_k_p_form_ = zeros(n_w_sum,1);
for nsource_gaussian=0:n_source_gaussian-1;
tmp_delta_ = tmp_delta_2s__(:,1+nsource_gaussian);
tmp_M_k_p_form_ = zeros(n_w_sum,1);
na=0;
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
n_w = n_w_(1+nk_p_r);
for nw=0:n_w-1;
k_x_c_0 = k_p_r*cos(2*pi*nw/n_w);
k_x_c_1 = k_p_r*sin(2*pi*nw/n_w);
tmp_M_k_p_form_(1+na) = exp( -( (2*pi*k_x_c_0).^2 + (2*pi*k_x_c_1).^2 ) / (2/tmp_sigma_x_c^2) ) .* exp( - 2*pi*i*( k_x_c_0*tmp_delta_(1+0) + k_x_c_1*tmp_delta_(1+1) ) );
na=na+1;
end;%for nw=0:n_w-1;
end;%for nk_p_r=0:n_k_p_r-1;
tmp_N_k_p_form_ = tmp_N_k_p_form_ + tmp_M_k_p_form_;
clear tmp_delta_;
end;%for nsource_gaussian=0:n_source_gaussian-1;
%%%%%%%%;
n_x_p_r = n_k_p_r;
x_p_r_ = k_p_r_*x_p_r_max/k_p_r_max;
x_c_0_all_ = k_c_0_all_*x_p_r_max/k_p_r_max;
x_c_1_all_ = k_c_1_all_*x_p_r_max/k_p_r_max;
%%%%%%%%;
tmp_R_x_p_form_ = ...
radon_k_p_to_x_p_xxnufft( ...
 n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,tmp_N_k_p_form_ ...
);
tmp_N_k_p_form_l2 = sum(abs(tmp_N_k_p_form_).^2 .* weight_2d_k_all_) * (2*pi)^2;
disp(sprintf(' %% tmp_N_k_p_form_l2 = %0.16f',tmp_N_k_p_form_l2));
disp(sprintf(' %% tmp_N_k_p_ vs tmp_N_k_p_form: %0.16f',fnorm(tmp_N_k_p_ - tmp_N_k_p_form_)/fnorm(tmp_N_k_p_)));
tmp_N_x_c_reco_ = interp_k_p_to_x_c_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_N_k_p_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_c^2) * n_w_sum;
tmp_N_x_c_reco_l2 = sum(abs(tmp_N_x_c_reco_).^2,'all')*dx^2;
disp(sprintf(' %% tmp_N_x_c_reco_l2 = %0.16f',tmp_N_x_c_reco_l2));
disp(sprintf(' %% tmp_N_x_c_ vs tmp_N_x_c_reco: %0.16f',fnorm(tmp_N_x_c_ - tmp_N_x_c_reco_)/fnorm(tmp_N_x_c_)));
%%%%%%%%;
flag_disp=0;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figbig;
subplot(2,2,1);imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,tmp_N_x_c_);axis image;axisnotick;
subplot(2,2,2);imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(tmp_N_k_p_));axis image;axisnotick;
subplot(2,2,3);imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(tmp_N_k_p_form_));axis image;axisnotick;
subplot(2,2,4);imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(tmp_N_x_c_reco_));axis image;axisnotick;
error('stopping');
end;%if flag_disp;
%%%%%%%%;

%%%%%%%%;
% visualize shadowsphere. ;
%%%%%%%%;
MSA_shape_longitudinal_perturbation_0_helper_0;

%%%%%%%%;
% visualize reconstruction. ;
%%%%%%%%;
MSA_shape_longitudinal_perturbation_0_helper_2;
MSA_shape_longitudinal_perturbation_0_helper_3;

%%%%%%%%;
% longitudinal perturbation. ;
%%%%%%%%;
MSA_shape_longitudinal_perturbation_0_helper_1;

if flag_verbose; disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
