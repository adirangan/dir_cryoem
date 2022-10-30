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
MSA_shape_0( ...
 parameter ...
,q_max ...
,n_gamma ...
);

str_thisfunction = 'MSA_shape_0';

if (nargin<1);
disp(sprintf(' %% testing %s',str_thisfunction));
parameter = struct('type','parameter');
parameter.flag_verbose = 1;
parameter.flag_disp = 1;
parameter.str_shape = 'horseshoe';
q_max = 12; n_gamma = 1024;
MSA_shape_0(parameter,q_max,n_gamma);
disp(sprintf(' %% returning')); return;
end;%if (nargin<1);

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); q_max=[]; end; na=na+1;
if (nargin<1+na); n_gamma=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
if ~isfield(parameter,'flag_disp'); parameter.flag_disp = 0; end;
if ~isfield(parameter,'str_shape'); parameter.str_shape = 'circle'; end;
if ~isfield(parameter,'x_max'); parameter.x_max = 1.35; end;
if ~isfield(parameter,'n_x'); parameter.n_x = 128; end;
if ~isfield(parameter,'sigma'); parameter.sigma = 0.075; end;
flag_verbose = parameter.flag_verbose;
flag_disp = parameter.flag_disp;
str_shape = parameter.str_shape;
x_max = parameter.x_max;
n_x = parameter.n_x;
sigma = parameter.sigma;

if flag_verbose; disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

if isempty(q_max); q_max = 8; end;
if isempty(n_gamma); n_gamma = max(1024,q_max*8); end;

%%%%%%%%;
q_ = transpose(-q_max:1:+q_max); n_q = numel(q_);
gamma_w_ = linspace(0,2*pi,1+n_gamma);
gamma_w_ = reshape(gamma_w_(1:n_gamma),[n_gamma,1]);
dgamma = mean(diff(gamma_w_));
n_w = n_gamma;
F_wq__ = exp(+i*gamma_w_*transpose(q_));
F_inv_qw__ = ctranspose(F_wq__)/n_gamma;
%%%%%%%%;

A_q_ = zeros(n_q,1);
%%%%%%%%;
if strcmp(str_shape,'circle');
A_q_(1+q_max+1) = 1;
end;%if strcmp(str_shape,'circle');
%%%%%%%%;
if strcmp(str_shape,'horseshoe');
pw_ = periodize(gamma_w_,-pi,+pi);
%pw_ = pi*(abs(pw_)/pi).^(1.5) .* sign(pw_);
g0_w_ = exp(-pw_.^2/(2*(pi/64)^2)); g0_w_ = g0_w_/max(abs(g0_w_));
g1_w_ = exp(-pw_.^2/(2*(pi/16)^2)); g1_w_ = g1_w_/max(abs(g1_w_));
h0_w_ = pw_.*exp(-pw_.^2/(2*(pi/128)^2)); h0_w_ = h0_w_/max(abs(h0_w_));
A_w_ = (cos(pw_) - 0.60*g0_w_ + 0.20*g1_w_) + i*(sin(pw_) + 0.15*h0_w_);
A_w_ = MSA_shape_speed1_0(n_w,gamma_w_,A_w_);
A_q_ = F_inv_qw__*A_w_;
end;%if strcmp(str_shape,'horseshoe');
%%%%%%%%;

A_w_ = F_wq__*A_q_;
lim_ = x_max*[-1,+1]; dx = diff(lim_)/n_x; dx2 = dx*dx;
x_ = transpose(linspace(lim_(1+0),lim_(1+1),n_x));
[x_0__,x_1__] = ndgrid(x_,x_);
x_r2__ = x_0__.^2 + x_1__.^2;
index_r_ = max(0,min(n_x-1,round(n_x*(real(A_w_)-min(lim_))/diff(lim_))));
index_c_ = max(0,min(n_x-1,round(n_x*(imag(A_w_)-min(lim_))/diff(lim_))));
A__ = sparse(1+index_r_,1+index_c_,1,n_x,n_x);
Ilim_ = [0,max(A__,[],'all')];
g__ = 1/(2*pi) * 1/(sigma).^2 * exp(-x_r2__/(2*sigma.^2));
g__ = g__/sum(g__,'all');
Ag__ = conv2(full(A__),g__,'same')*dx2;

if flag_disp;
figure(1);clf;figbig; ns=0;
%%%%;
subplot_{1+ns} = subplot(2,3,1+ns); ns=ns+1;
linewidth_use = 8;
s = surfline_0(real(A_w_),imag(A_w_),gamma_w_); set(s,'LineWidth',linewidth_use); colormap('hsv');
xlabel('real(A_w_)','Interpreter','none');
ylabel('imag(A_w_)','Interpreter','none');
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
tmp_Ag__ = conv2(full(A__),g__,'same')*dx2;
tmp_Ilim_ = [0,max(tmp_Ag__,[],'all')];
subplot_{1+ns} = subplot(2,3,1+ns); ns=ns+1;
imagesc(tmp_Ag__,tmp_Ilim_); set(gca,'ydir','normal'); axis image; axisnotick;
xlabel('imag(A_w_)','Interpreter','none');
ylabel('real(A_w_)','Interpreter','none');
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
end;%if flag_disp;

if flag_verbose; disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
