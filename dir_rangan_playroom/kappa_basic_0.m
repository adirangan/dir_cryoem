function ...
[ ...
 parameter ...
,kappa_norm_ ...
,chebfun_kernel_norm_ ...
,deconvolve_q ...
,relative_error_full_ ...
,l_max ...
,chebleg_d_ ...
] = ...
kappa_basic_0( ...
 parameter ...
,l_max ...
,chebleg_d_ ...
);
%%%%%%%%;
% Defines basic kappa (with no deconvolution) with no quadratic-program. ;
%%%%%%%%;

str_thisfunction = 'kappa_basic_0';

%%%%%%%%;
if nargin<1;
disp(sprintf(' %% testing %s',str_thisfunction));
nf=0;
l_max = 49; l_val_ = transpose([0:l_max]);
chebleg_d_ = [];
parameter = struct('type','parameter');
parameter.flag_verbose = 1;
parameter.flag_disp = 1;
[ ...
 parameter ...
,kappa_norm_ ...
,chebfun_kernel_norm_ ...
,deconvolve_q ...
,relative_error_full_ ...
,l_max ...
,chebleg_d_ ...
] = ...
kappa_basic_0( ...
 parameter ...
,l_max ...
,chebleg_d_ ...
);
%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
markersize_use = 12;
subplot(1,1,1);
hold on;
plot(l_val_,log10(relative_error_full_),'k.','MarkerSize',markersize_use);
hold off;
xlabel('l_val_','Interpreter','none');
ylabel('log10(re)','Interpreter','none');
xlim([-1,1+l_max]); set(gca,'XTick',l_val_(1:7:end));
ylim([-16,+1]); set(gca,'YTick',-15:0);
grid on;
title(sprintf('basic'),'Interpreter','none');
%%%%%%%%;
disp('returning'); return;
end;%if nargin<1;
%%%%%%%%;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); l_max=[]; end; na=na+1;
if (nargin<1+na); chebleg_d_=[]; end; na=na+1;

%%%%;
if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
flag_verbose=parameter.flag_verbose;
if ~isfield(parameter,'flag_disp'); parameter.flag_disp=0; end;
flag_disp=parameter.flag_disp; nf=0;
if ~isfield(parameter,'tolerance_master'); parameter.tolerance_master=1e-2; end;
tolerance_master=parameter.tolerance_master;
if ~isfield(parameter,'kernel_basic_l_max_band'); parameter.kernel_basic_l_max_band=+Inf; end;
kernel_basic_l_max_band=parameter.kernel_basic_l_max_band;
%%%%;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

%%%%%%%%;
if isempty(l_max); l_max = 49; end;
l_val_ = transpose([0:l_max]);
%%%%%%%%;
l_max_band = min(l_max,kernel_basic_l_max_band);
%%%%%%%%;
if isempty(chebleg_d_);
chebleg_d_ = cell(1+l_max,1);
for l_val=0:l_max;
tmp_c_ = [zeros(l_val,1);1];
chebleg_d_{1+l_val} = chebfun(leg2cheb(tmp_c_,'norm'),'coeffs');
end;%for l_val=0:l_max;
end;%if isempty(chebleg_d_);
%%%%%%%%;
deltafunc_formula_ = sqrt(4*pi)*sqrt(1+2*l_val_);
deltafunc_band_formula_ = deltafunc_formula_;
if isfinite(l_max_band); tmp_index_ = efind(l_val_> l_max_band); deltafunc_band_formula_(1+tmp_index_) = 0; end;

%%%%%%%%;
n_a_use = 1+2*l_max + 16; %<-- Need to integrate polynomials of degree l_max^2. ;
[a_full_node_,a_full_weight_] = legpts(n_a_use,[0,pi]);
%%%%%%%%;
% Note, this is normalized so that the integral of: ;
% chebleg_d_{1+l_val_}^{2} ;
% over the surface of the sphere is (2*pi), and not 1.0d0. ;
% This means that the integral of: ;
% chebleg_d_{1+0}^{1} ;
% over the surface of the sphere will be (2*pi)*sqrt(2). ;
% Thus, the surface-integral of chebfun_kernel_norm_ (below) ;
% will be (2*pi)*sqrt(2). ;
% Finally, after deconvolving by deconvolve_l_(1+0) = sqrt(4*pi), ;
% the surface-integral: ;
% sum(2*pi*a_full_weight_*(chebfun_kernel_norm_(cos(a_full_node_)).*chebleg_d_{1+l_val}(cos(a_full_node_)).*sin(a_full_node_))) * deconvolve_l_(1+l_val) / (4*pi) ;
% will equal sqrt(pi), ;
% and the surface-integral: ;
% sum(2*pi*a_full_weight_*(chebfun_kernel_norm_(cos(a_full_node_)).*sin(a_full_node_))) * deconvolve_l_(1+l_val) / (4*pi) ;
% will equal sqrt(2*pi). ;
%%%%%%%%;
if (flag_disp>1);
I_dd__ = zeros(1+l_max,1+l_max);
for l_val_0=0:l_max;
for l_val_1=0:l_max;
tmp_I = sum(2*pi*a_full_weight_*(chebleg_d_{1+l_val_0}(cos(a_full_node_)).*chebleg_d_{1+l_val_1}(cos(a_full_node_)).*sin(a_full_node_)));
I_dd__(1+l_val_0,1+l_val_1) = tmp_I;
end;%for l_val_1=0:l_max;
end;%for l_val_0=0:l_max;
figure(1+nf);nf=nf+1;clf;figsml; fig80s;
imagesc(log10(abs(I_dd__-2*pi*eye(1+l_max))),[-15,0]); colorbar;
axis image; axisnotick;
title('log10(abs(I_dd__ - eye))','Interpreter','none');
end;%if (flag_disp>1);
%%%%%%%%;

%%%%%%%%%%%%%%%%;
flag_kernel_full = 1; if flag_kernel_full==1; kappa_ = deltafunc_band_formula_; end; %<-- basic. ;
%%%%%%%%%%%%%%%%;

%%%%%%%%
chebfun_kernel_basic_ = chebfun(0);
for l_val=0:l_max;
chebfun_kernel_basic_ = chebfun_kernel_basic_ + kappa_(1+l_val)*chebleg_d_{1+l_val};
end;%for l_val=0:l_max;
const_avg = sum(2*pi*a_full_weight_*(chebleg_d_{1+0}(cos(a_full_node_)).^1.*sin(a_full_node_)));
if (flag_verbose>0); disp(sprintf(' %% const_avg vs 2*sqrt(2)*pi %0.16f',fnorm(const_avg - 2*sqrt(2)*pi))); end;
chebfun_kernel_basic_avg = sum(2*pi*a_full_weight_*(chebfun_kernel_basic_(cos(a_full_node_)).^1.*sin(a_full_node_)));
chebfun_kernel_basic_avg = chebfun_kernel_basic_avg / max(1e-12,abs(const_avg));
if (flag_verbose>0); disp(sprintf(' %% chebfun_kernel_basic_avg: %0.16f',chebfun_kernel_basic_avg)); end;
chebfun_kernel_norm_ = chebfun_kernel_basic_ / max(1e-12,chebfun_kernel_basic_avg) ;
kappa_norm_ = kappa_ / max(1e-12,chebfun_kernel_basic_avg) ;
%%%%%%%%;

%%%%%%%%;
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figmed;
linewidth_use = 2;
c_use__ = colormap_81s; n_c_use = size(c_use__,1);
subplot(1,3,[1,2]);
n_z = 1+1024; z_ = transpose(linspace(-1,+1,n_z));
hold on;
plot(log(1-z_),chebfun_kernel_norm_(z_),'-','Color',0.85*[1,1,1],'LineWidth',linewidth_use);
hold off;
xlabel('log(1-z)'); ylabel('P(z)'); title('chebfun_kernel_norm_','Interpreter','none'); grid on;
subplot(1,3,3);
hold on;
plot(l_val_,kappa_norm_,'x-','Color',0.85*[1,1,1],'LineWidth',3);
hold off;
xlabel('l_val','Interpreter','none');
ylabel('kappa_norm_','Interpreter','none');
title('scaling factor','Interpreter','none');
end;%if (flag_disp>1);
sum_l1_full = sum(2*pi*a_full_weight_*(chebfun_kernel_norm_(cos(a_full_node_)).^1.*sin(a_full_node_)));
sum_l2_full = sum(2*pi*a_full_weight_*(chebfun_kernel_norm_(cos(a_full_node_)).^2.*sin(a_full_node_)));
if (flag_verbose>0); disp(sprintf(' %% sum_l1_full: %0.16f',sum_l1_full)); end;
if (flag_verbose>0); disp(sprintf(' %% sum_l2_full: %0.16f',sum_l2_full)); end;
%%%%%%%%;

%%%%%%%%;
deconvolve_q = sqrt(4*pi);
deconvolve_l_ = deltafunc_band_formula_./max(1e-12,kappa_norm_(1+l_val_));
deltafunc_band_from_chebfun_mollify_ = zeros(1+l_max,1);
for l_val=0:l_max;
tmp_I = sum(2*pi*a_full_weight_*(chebfun_kernel_norm_(cos(a_full_node_)).*chebleg_d_{1+l_val}(cos(a_full_node_)).*sin(a_full_node_)));
deltafunc_band_from_chebfun_mollify_(1+l_val) = tmp_I;
end;%for l_val=0:l_max;
deltafunc_band_from_chebfun_restore_ = deltafunc_band_from_chebfun_mollify_ .* deconvolve_l_ ;
relative_error_full_ = abs(deltafunc_band_formula_ - deltafunc_band_from_chebfun_restore_/(2*pi))./abs(deltafunc_band_formula_);
%%%%%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;
markersize_use = 16;
hold on;
plot(l_val_,log10(relative_error_full_),'k.','MarkerSize',markersize_use);
hold off;
xlabel('l_val_','Interpreter','none');
ylabel('log10(abs(form-reco)./abs(form))','Interpreter','none');
xlim([-1,1+l_max]); set(gca,'XTick',l_val_(1:7:end));
ylim([-16,+1]); set(gca,'YTick',-15:0);
grid on;
end;%if flag_disp;
%%%%%%%%;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

