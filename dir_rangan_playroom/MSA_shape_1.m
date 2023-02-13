function ...
[ ...
 parameter ...
,q_ ...
,gamma_w_ ...
,F_wq__ ...
,F_inv_qw__ ...
,A_q_ ...
,A_w_ ...
,x_ ...
,x_0__ ...
,x_1__ ...
,Ag__ ...
] = ...
MSA_shape_1( ...
 parameter ...
,q_max ...
,n_gamma ...
);

str_thisfunction = 'MSA_shape_1';

if (nargin<1);
disp(sprintf(' %% testing %s',str_thisfunction));
parameter = struct('type','parameter');
parameter.flag_verbose = 1;
parameter.flag_disp = 1;
parameter.str_shape = 'figure8'; parameter.x_max = 2.5; parameter.t_diffuse = 1e-4; q_max = 128; n_gamma = 1024*2;
MSA_shape_1(parameter,q_max,n_gamma);
disp(sprintf(' %% returning')); return;
end;%if (nargin<1);

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); q_max=[]; end; na=na+1;
if (nargin<1+na); n_gamma=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
if ~isfield(parameter,'flag_disp'); parameter.flag_disp = 0; end;
if ~isfield(parameter,'str_shape'); parameter.str_shape = 'butterfly'; end;
if ~isfield(parameter,'t_diffuse'); parameter.t_diffuse = 0.0; end;
if ~isfield(parameter,'x_max'); parameter.x_max = 1.35; end;
if ~isfield(parameter,'n_x'); parameter.n_x = 128; end;
if ~isfield(parameter,'sigma'); parameter.sigma = 0.075; end;
flag_verbose = parameter.flag_verbose;
flag_disp = parameter.flag_disp; nf=0;
str_shape = parameter.str_shape;
t_diffuse = parameter.t_diffuse;
x_max = parameter.x_max;
n_x = parameter.n_x;
sigma = parameter.sigma;

if flag_verbose; disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

if isempty(q_max); q_max = 64; end;
if isempty(n_gamma); n_gamma = max(1024,q_max*8); end;

if strcmp(str_shape,'figure8'); n_gamma = (14*4)*ceil(n_gamma/(14*4)); end;
n_gamma = n_gamma + mod(n_gamma,2); %<-- ensure that n_gamma is even. ;

%%%%%%%%;
q_ = transpose(-q_max:1:+q_max); n_q = numel(q_);
gamma_w_ = linspace(0,2*pi,1+n_gamma);
gamma_w_ = reshape(gamma_w_(1:n_gamma),[n_gamma,1]);
dgamma = mean(diff(gamma_w_));
n_w = n_gamma;
F_wq__ = exp(+i*gamma_w_*transpose(q_));
F_inv_qw__ = ctranspose(F_wq__)/n_gamma;
%%%%%%%%;

n_type = 3;
A_q_ = zeros(n_q,1);
B_q_ = zeros(n_q,1);
C_q_ = zeros(n_q,1);
%%%%%%%%;
if strcmp(str_shape,'figure8');
%%%%;
% 14 segments to the path. ;
% 00: (+0,+0) -> (+1,+1) ; <-- speed 1. ;
% 01: (+1,+1) -> (+0,+2) ; <-- speed 4. ;
% 02: (+0,+2) -> (+1,+1) ; <-- speed 4. ;
% 03: (+1,+1) -> (-1,+1) ; <-- speed 2. ;
% 04: (-1,+1) -> (+0,+2) ; <-- speed 4. ;
% 05: (+0,+2) -> (-1,+1) ; <-- speed 4. ;
% 06: (-1,+1) -> (+0,+0) ; <-- speed 1. ;
% 07: (+0,+0) -> (+1,-1) ; <-- speed 1. ;
% 08: (+1,-1) -> (+0,-2) ; <-- speed 4. ;
% 09: (+0,-2) -> (+1,-1) ; <-- speed 4. ;
% 10: (+1,-1) -> (-1,-1) ; <-- speed 2. ;
% 11: (-1,-1) -> (+0,-2) ; <-- speed 4. ;
% 12: (+0,-2) -> (-1,-1) ; <-- speed 4. ;
% 13: (-1,-1) -> (+0,+0) ; <-- speed 1. ;
%%%%;
assert(n_gamma==56*floor(n_gamma/56));
n_gamma_d07 = n_gamma/ 7; index_d07_ = transpose([0:n_gamma_d07-1]); gamma_d07_ = transpose(linspace(0,pi/2,n_gamma_d07)); %<-- note both endpoints included. ;
n_gamma_d14 = n_gamma/14; index_d14_ = transpose([0:n_gamma_d14-1]); gamma_d14_ = transpose(linspace(0,pi/2,n_gamma_d14)); %<-- note both endpoints included. ;
n_gamma_d28 = n_gamma/28; index_d28_ = transpose([0:n_gamma_d28-1]); gamma_d28_ = transpose(linspace(0,pi/2,n_gamma_d28)); %<-- note both endpoints included. ;
ng=0;
index_segment_00_ = ng + index_d07_; ng = ng + n_gamma_d07; A_00_w_ = (+0 + cos(gamma_d07_ - 1*pi/2)) + i*(+1 + sin(gamma_d07_ - 1*pi/2));
index_segment_01_ = ng + index_d28_; ng = ng + n_gamma_d28; A_01_w_ = (+0 + cos(gamma_d28_ - 0*pi/2)) + i*(+1 + sin(gamma_d28_ - 0*pi/2));
index_segment_02_ = ng + index_d28_; ng = ng + n_gamma_d28; A_02_w_ = flip(A_01_w_);
index_segment_03_ = ng + index_d14_; ng = ng + n_gamma_d14; A_03_w_ = transpose(linspace(+1 + i*1,-1+i*1,n_gamma_d14));
index_segment_04_ = ng + index_d28_; ng = ng + n_gamma_d28; A_04_w_ = (+0 - cos(gamma_d28_ - 0*pi/2)) + i*(+1 + sin(gamma_d28_ - 0*pi/2));
index_segment_05_ = ng + index_d28_; ng = ng + n_gamma_d28; A_05_w_ = flip(A_04_w_);
index_segment_06_ = ng + index_d07_; ng = ng + n_gamma_d07; A_06_w_ = (+0 + cos(gamma_d07_ + 2*pi/2)) + i*(+1 + sin(gamma_d07_ + 2*pi/2));
index_segment_07_ = ng + index_d07_; ng = ng + n_gamma_d07; A_07_w_ = conj(A_00_w_);
index_segment_08_ = ng + index_d28_; ng = ng + n_gamma_d28; A_08_w_ = conj(A_01_w_);
index_segment_09_ = ng + index_d28_; ng = ng + n_gamma_d28; A_09_w_ = conj(A_02_w_);
index_segment_10_ = ng + index_d14_; ng = ng + n_gamma_d14; A_10_w_ = conj(A_03_w_);
index_segment_11_ = ng + index_d28_; ng = ng + n_gamma_d28; A_11_w_ = conj(A_04_w_);
index_segment_12_ = ng + index_d28_; ng = ng + n_gamma_d28; A_12_w_ = conj(A_05_w_);
index_segment_13_ = ng + index_d07_; ng = ng + n_gamma_d07; A_13_w_ = conj(A_06_w_);
assert(ng==n_gamma);
A_w_pre_ =  [ ...
A_00_w_ ; ...
A_01_w_ ; ...
A_02_w_ ; ...
A_03_w_ ; ...
A_04_w_ ; ...
A_05_w_ ; ...
A_06_w_ ; ...
A_07_w_ ; ...
A_08_w_ ; ...
A_09_w_ ; ...
A_10_w_ ; ...
A_11_w_ ; ...
A_12_w_ ; ...
A_13_w_ ; ...
];
if t_diffuse>0;
A_w_pre_ = MSA_shape_diffuse_0(n_gamma,gamma_w_,A_w_pre_,t_diffuse);
end;%if t_diffuse>0;
A_q_ = F_inv_qw__*A_w_pre_;
%%%%%%%%;
n_gamma_d02 = n_gamma/ 2; index_d02_ = transpose([0:n_gamma_d02-1]); gamma_d02_ = transpose(linspace(0,2*pi,n_gamma_d02));
B_w_pre_ = zeros(n_gamma,1);
B_w_pre_(1 + 0*n_gamma_d02 + index_d02_) = cos(gamma_d02_ - 1*pi/2) + i*(1 + sin(gamma_d02_ - 1*pi/2));
B_w_pre_(1 + 1*n_gamma_d02 + index_d02_) = conj(B_w_pre_(1 + 0*n_gamma_d02 + index_d02_));
if t_diffuse>0;
B_w_pre_ = MSA_shape_diffuse_0(n_gamma,gamma_w_,B_w_pre_,t_diffuse);
end;%if t_diffuse>0;
B_q_ = F_inv_qw__*B_w_pre_;
%%%%%%%%;
C_w_pre_ = real(B_w_pre_) + i*max(-1,min(+1,imag(B_w_pre_)));
if t_diffuse>0;
C_w_pre_ = MSA_shape_diffuse_0(n_gamma,gamma_w_,C_w_pre_,t_diffuse);
end;%if t_diffuse>0;
C_q_ = F_inv_qw__*C_w_pre_;
end;%if strcmp(str_shape,'figure8');
%%%%%%%%;

A_w_ = F_wq__*A_q_;
%%%%;
lim_ = x_max*[-1,+1]; dx = diff(lim_)/n_x; dx2 = dx*dx;
x_ = transpose(linspace(lim_(1+0),lim_(1+1),n_x));
[x_0__,x_1__] = ndgrid(x_,x_);
x_r2__ = x_0__.^2 + x_1__.^2;
%%%%;
index_r_ = max(0,min(n_x-1,round(n_x*(real(A_w_)-min(lim_))/diff(lim_))));
index_c_ = max(0,min(n_x-1,round(n_x*(imag(A_w_)-min(lim_))/diff(lim_))));
A__ = sparse(1+index_r_,1+index_c_,1,n_x,n_x);
Ilim_ = [0,max(A__,[],'all')];
g__ = 1/(2*pi) * 1/(sigma).^2 * exp(-x_r2__/(2*sigma.^2));
g__ = g__/sum(g__,'all');
Ag__ = conv2(full(A__),g__,'same')*dx2;

B_w_ = F_wq__*B_q_;
index_r_ = max(0,min(n_x-1,round(n_x*(real(B_w_)-min(lim_))/diff(lim_))));
index_c_ = max(0,min(n_x-1,round(n_x*(imag(B_w_)-min(lim_))/diff(lim_))));
B__ = sparse(1+index_r_,1+index_c_,1,n_x,n_x);
Bg__ = conv2(full(B__),g__,'same')*dx2;

C_w_ = F_wq__*C_q_;
index_r_ = max(0,min(n_x-1,round(n_x*(real(C_w_)-min(lim_))/diff(lim_))));
index_c_ = max(0,min(n_x-1,round(n_x*(imag(C_w_)-min(lim_))/diff(lim_))));
C__ = sparse(1+index_r_,1+index_c_,1,n_x,n_x);
Cg__ = conv2(full(C__),g__,'same')*dx2;

if flag_disp;
%%%%%%%%;
% illustrate curve A_w_pre_ in the complex plane. ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
%%%%;
for ntype=0:n_type-1;
if ntype==0; tmp_w_pre_ = A_w_pre_; tmp_str = 'A'; end;
if ntype==1; tmp_w_pre_ = B_w_pre_; tmp_str = 'D';  end;
if ntype==2; tmp_w_pre_ = C_w_pre_; tmp_str = 'C';  end;
subplot_{1+ntype} = subplot(1,n_type,1+ntype);
linewidth_use = 8;
s = surfline_0(real(tmp_w_pre_),imag(tmp_w_pre_),gamma_w_); set(s,'LineWidth',linewidth_use); colormap('hsv');
xlabel(sprintf('real(%s_w_pre_)',tmp_str),'Interpreter','none');
ylabel(sprintf('imag(%s_w_pre_)',tmp_str),'Interpreter','none');
xlim(x_max*[-1,+1]);
ylim(x_max*[-1,+1]);
set(gca,'XTick',-x_max:0.1:+x_max);
set(gca,'YTick',-x_max:0.1:+x_max);
xtickangle(90); grid on;
axis square;
title('image ring');
%%%%;
colormap(subplot_{1+ntype},colormap('hsv'));
%%%%;
sgtitle(str_shape,'Interpreter','none');
close(gcf);
end;%for ntype=0:n_type-1;
end;%if flag_disp;

if flag_disp;
%%%%%%%%;
% illustrate curve A_w_ in the complex plane. ;
%%%%%%%%;
for ntype=0:n_type-1;
if ntype==0; tmp_w_ = A_w_; tmp_str = 'A'; tmp__ = A__; end;
if ntype==1; tmp_w_ = B_w_; tmp_str = 'D'; tmp__ = B__;   end;
if ntype==2; tmp_w_ = C_w_; tmp_str = 'C'; tmp__ = C__;   end;
figure(1+nf);nf=nf+1;clf;figbig; ns=0;
%%%%;
subplot_{1+ns} = subplot(2,3,1+ns); ns=ns+1;
linewidth_use = 8;
s = surfline_0(real(tmp_w_),imag(tmp_w_),gamma_w_); set(s,'LineWidth',linewidth_use); colormap('hsv');
xlabel(sprintf('real(%s_w_)',tmp_str),'Interpreter','none');
ylabel(sprintf('imag(%s_w_)',tmp_str),'Interpreter','none');
xlim(x_max*[-1,+1]);
ylim(x_max*[-1,+1]);
set(gca,'XTick',-x_max:0.1:+x_max);
set(gca,'YTick',-x_max:0.1:+x_max);
xtickangle(90); grid on;
axis square;
title('image ring');
%%%%;
sigma_ = [0.001;0.050;0.075;0.100;0.125]; n_sigma = numel(sigma_);
for nsigma=0:n_sigma-1;
tmp_sigma = sigma_(1+nsigma);
g__ = 1/(2*pi) * 1/(tmp_sigma).^2 * exp(-x_r2__/(2*tmp_sigma.^2));
g__ = g__/sum(g__,'all');
tmp_g__ = conv2(full(tmp__),g__,'same')*dx2;
tmp_Ilim_ = [0,max(tmp_g__,[],'all')];
subplot_{1+ns} = subplot(2,3,1+ns); ns=ns+1;
imagesc(tmp_g__,tmp_Ilim_); set(gca,'ydir','normal'); axis image; axisnotick;
xlabel(sprintf('imag(%s_w_)',tmp_str),'Interpreter','none');
ylabel(sprintf('real(%s_w_)',tmp_str),'Interpreter','none');
title(sprintf('sigma=%0.3f',tmp_sigma),'Interpreter','none');
end;%for nsigma=0:n_sigma-1;
%%%%;
ns=0;
colormap(subplot_{1+ns},colormap('hsv')); ns=ns+1;
for nsigma=0:n_sigma-1;
colormap(subplot_{1+ns},colormap_80s); ns=ns+1;
end;%for nsigma=0:n_sigma-1;
%%%%;
sgtitle(str_shape,'Interpreter','none');
close(gcf);
end;%for ntype=0:n_type-1;
%%%%%%%%;
end;%if flag_disp;

if flag_disp;
%%%%%%%%;
% Illustrate the corresponding real-space function. ;
%%%%%%%%;
k_p_r_max = 48/(2*pi); k_eq_d = 1.0/(2*pi);
l_max_upb = round(2*pi*k_p_r_max); %<-- typically sufficient for 2-3 digits of precision. ;
l_max_max = l_max_upb;
n_w_max = 2*(l_max_max+1); %<-- overwrite with provided gamma_w_. ;
n_w_max = n_gamma;
[ ...
 n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
] = ...
get_weight_3d_1( ...
 flag_verbose ...
,k_p_r_max ...
,k_eq_d ...
);
%%%%;
[ ...
 n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_k_all_ ...
] = ...
get_weight_2d_1( ...
 flag_verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,0 ...
,n_w_max*ones(n_k_p_r,1) ...
);
n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);
k_c_0_all_ = reshape(cos(gamma_w_)*transpose(k_p_r_),[n_w_sum,1]);
k_c_1_all_ = reshape(sin(gamma_w_)*transpose(k_p_r_),[n_w_sum,1]);
%%%%;
gr_ = exp(-(k_p_r_-(k_p_r_max*0.50)).^2/(2*(k_p_r_max/12)^2));
%gr_ = k_p_r_.*exp(-(k_p_r_.^2)/(k_p_r_max/2));
gr_ = gr_/sum(gr_.*weight_2d_k_p_r_);
k_shift = 0.5;
a_k_p_wk__ = zeros(n_w_sum,n_k_p_r);
a_k_p_wk__ = (A_w_ + k_shift)*transpose(gr_);
b_k_p_wk__ = zeros(n_w_sum,n_k_p_r);
b_k_p_wk__ = (B_w_ + k_shift)*transpose(gr_);
c_k_p_wk__ = zeros(n_w_sum,n_k_p_r);
c_k_p_wk__ = (C_w_ + k_shift)*transpose(gr_);
%%%%;
x_p_r_max = 1.0;
n_x_c = 128;
x_c_0_ = transpose(linspace(-x_p_r_max,+x_p_r_max,n_x_c));
x_c_1_ = transpose(linspace(-x_p_r_max,+x_p_r_max,n_x_c));
[x_c_0__,x_c_1__] = ndgrid(x_c_0_,x_c_1_); n_xx_c = n_x_c^2;
xx_c_weight = (2*x_p_r_max/n_x_c)^2;
%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
a_x_c_xx__ = reshape(xxnufft2d3(n_w_sum,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,a_k_p_wk__(:).*(2*pi)^2.*weight_2d_k_all_,+1,1e-12,n_xx_c,x_c_0__(:)/eta,x_c_1__(:)/eta)/sqrt(2*pi)/sqrt(2*pi),[n_x_c,n_x_c]);
b_x_c_xx__ = reshape(xxnufft2d3(n_w_sum,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,b_k_p_wk__(:).*(2*pi)^2.*weight_2d_k_all_,+1,1e-12,n_xx_c,x_c_0__(:)/eta,x_c_1__(:)/eta)/sqrt(2*pi)/sqrt(2*pi),[n_x_c,n_x_c]);
c_x_c_xx__ = reshape(xxnufft2d3(n_w_sum,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,c_k_p_wk__(:).*(2*pi)^2.*weight_2d_k_all_,+1,1e-12,n_xx_c,x_c_0__(:)/eta,x_c_1__(:)/eta)/sqrt(2*pi)/sqrt(2*pi),[n_x_c,n_x_c]);
%%%%%%%%;
for ntype=0:n_type-1;
if ntype==0; tmp_k_p_wk__ = a_k_p_wk__; tmp_str = 'A'; tmp_x_c_xx__ = a_x_c_xx__; end;
if ntype==1; tmp_k_p_wk__ = b_k_p_wk__; tmp_str = 'B'; tmp_x_c_xx__ = b_x_c_xx__; end;
if ntype==2; tmp_k_p_wk__ = c_k_p_wk__; tmp_str = 'C'; tmp_x_c_xx__ = c_x_c_xx__; end;
figure(1+nf);nf=nf+1;clf;figbig; ns=0;
p_row = 2; p_col = 2;
%%%%;
tmp_k_p_lim_ = prctile(abs(tmp_k_p_wk__(:)),[95])*[-1,+1];
tmp_x_c_lim_ = prctile(abs(tmp_x_c_xx__(:)),[95])*[-1,+1];
subplot(p_row,p_col,1+0+0*p_col); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(tmp_k_p_wk__(:)),tmp_k_p_lim_,colormap_beach); axis image;
title(sprintf('real(%s_k_p_wk__(:))',tmp_str),'Interpreter','none');
subplot(p_row,p_col,1+0+1*p_col); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(tmp_k_p_wk__(:)),tmp_k_p_lim_,colormap_beach); axis image;
title(sprintf('imag(%s_k_p_wk__(:))',tmp_str),'Interpreter','none');
subplot(p_row,p_col,1+1+0*p_col); imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(tmp_x_c_xx__(:)),tmp_x_c_lim_,colormap_beach); axis image;
title(sprintf('real(%s_x_c_xx__(:))',tmp_str),'Interpreter','none');
subplot(p_row,p_col,1+1+1*p_col); imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,imag(tmp_x_c_xx__(:)),tmp_x_c_lim_,colormap_beach); axis image;
title(sprintf('imag(%s_x_c_xx__(:))',tmp_str),'Interpreter','none');
%%%%;
sgtitle(str_shape,'Interpreter','none');
close(gcf);
end;%for ntype=0:n_type-1;
%%%%%%%%;
end;%if flag_disp;

if flag_disp;
%%%%%%%%;
% show all three cross-sections. ;
%%%%%%%%;
dir_jpg = sprintf('%s/dir_jpg',pwd); if ~exist(dir_jpg,'dir'); mkdir(dir_jpg); end;
fname_fig_pre = sprintf('%s/MSA_shape_ABC_crosssection_FIGA',dir_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
figure(1+nf);nf=nf+1;clf;figmed;
p_row=1;p_col=3;
for ntype=0:n_type-1;
if ntype==0; tmp_k_p_wk__ = a_k_p_wk__; tmp_str = 'A'; tmp_x_c_xx__ = a_x_c_xx__; end;
if ntype==1; tmp_k_p_wk__ = b_k_p_wk__; tmp_str = 'B'; tmp_x_c_xx__ = b_x_c_xx__; end;
if ntype==2; tmp_k_p_wk__ = c_k_p_wk__; tmp_str = 'C'; tmp_x_c_xx__ = c_x_c_xx__; end;
tmp_x_c_lim_ = prctile(abs(tmp_x_c_xx__(:)),[95])*[-1,+1];
subplot(p_row,p_col,1+ntype);cla;
imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(tmp_x_c_xx__(:)),tmp_x_c_lim_,colormap_beach); axis image; axisnotick;
title(sprintf('%s',tmp_str),'Interpreter','none');
end;%for ntype=0:n_type-1;
disp(sprintf(' %% writing %s',fname_fig_jpg));
print('-djpeg',fname_fig_jpg);
%close(gcf);
end;%if flag_disp;

flag_disp=0;
if flag_disp;
%%%%%%%%;
% illustrate the associated molecule. ;
%%%%%%%%;
dir_jpg = sprintf('%s/dir_jpg',pwd); if ~exist(dir_jpg,'dir'); mkdir(dir_jpg); end;
x_c_2_ = linspace(-x_p_r_max,+x_p_r_max,n_x_c);
[x_c_0___,x_c_1___,x_c_2___] = ndgrid(x_c_0_,x_c_1_,x_c_2_); n_xxx_c = n_x_c^3;
xxx_c_weight_ = (2*x_p_r_max/n_x_c)^3;
sigma_z = 0.50;
gz_ = 1/sqrt(2*pi)/sigma_z * exp(-(x_c_2_.^2)/(2*sigma_z^2));
x_shift = 0.0;
a_x_c_xxx___ = bsxfun(@times,repmat((a_x_c_xx__ + x_shift),[1,1,n_x_c]),reshape(gz_,[1,1,n_x_c]));
b_x_c_xxx___ = bsxfun(@times,repmat((b_x_c_xx__ + x_shift),[1,1,n_x_c]),reshape(gz_,[1,1,n_x_c]));
c_x_c_xxx___ = bsxfun(@times,repmat((c_x_c_xx__ + x_shift),[1,1,n_x_c]),reshape(gz_,[1,1,n_x_c]));
%%%%%%%%;
for ntype=0:n_type-1;
if ntype==0; tmp_x_c_xxx___ = a_x_c_xxx___; tmp_str = 'A'; end;
if ntype==1; tmp_x_c_xxx___ = b_x_c_xxx___; tmp_str = 'B'; end;
if ntype==2; tmp_x_c_xxx___ = c_x_c_xxx___; tmp_str = 'C'; end;
fname_fig_pre = sprintf('%s/MSA_shape_%s_FIGA',dir_jpg,tmp_str);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
figure(1+nf);nf=nf+1;clf;figbig;
prctile_ = [98.5];%prctile_ = [97.5,98.5,99.5];
subplot(1,1,1); isosurface_f_x_u_0(tmp_x_c_xxx___,prctile_); title(sprintf('%s_x_c_xxx___',tmp_str),'Interpreter','none');
xlabel('x'); ylabel('y'); zlabel('z'); axis vis3d;
disp(sprintf(' %% writing %s',fname_fig_jpg));
print('-djpeg',fname_fig_jpg);
close(gcf);
end;%for ntype=0:n_type-1;
%%%%%%%%;
end;%if flag_disp;

if flag_verbose; disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
