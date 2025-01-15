function ...
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
%%%%%%%%;
% Defines kappa (with option for deconvolution) using quadratic-program. ;
%%%%%%%%;

str_thisfunction = 'kappa_qpro_0';

%%%%%%%%;
if nargin<1;
disp(sprintf(' %% testing %s',str_thisfunction));
nf=0;
l_max = 49; l_val_ = transpose([0:l_max]);
n_pole_north = 4;
pole_north_ = transpose(linspace(1,12,n_pole_north)*pi/24);
n_pole_south = 3;
pole_south_ = transpose(linspace(0,12,n_pole_south)*pi/24);
chebleg_d_ = [];
[pole_north_ns__,pole_south_ns__] = ndgrid(pole_north_,pole_south_);
kappa_sparse_f_ns__ = zeros(n_pole_north,n_pole_south);
relative_error_full_lns___ = zeros(1+l_max,n_pole_north,n_pole_south);
relative_error_crop_lns___ = zeros(1+l_max,n_pole_north,n_pole_south);
for npole_north=0:n_pole_north-1;
disp(sprintf(' %% npole_north %d/%d',npole_north,n_pole_north));
for npole_south=0:n_pole_south-1;
pole_north = pole_north_(1+npole_north);
pole_south = pole_south_(1+npole_south);
parameter = struct('type','parameter');
parameter.flag_verbose = 1;
parameter.flag_disp = 1;
parameter.kernel_qpro_polar_a_pole_north = pole_north;
parameter.kernel_qpro_polar_a_pole_south = pole_south;
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
kappa_sparse_f_ns__(1+npole_north,1+npole_south) = kappa_sparse_f;
relative_error_full_lns___(:,1+npole_north,1+npole_south) = relative_error_full_;
relative_error_crop_lns___(:,1+npole_north,1+npole_south) = relative_error_crop_;
end;%for npole_south=0:n_pole_south-1;
end;%for npole_north=0:n_pole_north-1;
%%%%;
figure(1+nf);nf=nf+1;clf;figsml;fig80s;
imagesc(kappa_sparse_f_ns__); colorbar;
axis image; axisnotick; set(gca,'Ydir','normal');
xlabel('pole_south_','Interpreter','none');
ylabel('pole_north_','Interpreter','none');
%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
markersize_use = 12;
p_row = n_pole_north; p_col = n_pole_south; np=0;
for npole_north=0:n_pole_north-1;
for npole_south=0:n_pole_south-1;
subplot(p_row,p_col,1+np);np=np+1;
pole_north = pole_north_(1+npole_north); str_pole_north = strtrim(rats(pole_north/(pi/24)));
pole_south = pole_south_(1+npole_south); str_pole_south = strtrim(rats(pole_south/(pi/24)));
kappa_sparse_f = kappa_sparse_f_ns__(1+npole_north,1+npole_south);
relative_error_full_ = relative_error_full_lns___(:,1+npole_north,1+npole_south);
relative_error_crop_ = relative_error_crop_lns___(:,1+npole_north,1+npole_south);
hold on;
plot(l_val_,log10(relative_error_full_),'k.','MarkerSize',markersize_use);
plot(l_val_,log10(relative_error_crop_),'r.','MarkerSize',markersize_use);
hold off;
xlabel('l_val_','Interpreter','none');
ylabel('log10(re)','Interpreter','none');
xlim([-1,1+l_max]); set(gca,'XTick',l_val_(1:7:end));
ylim([-16,+1]); set(gca,'YTick',-15:0);
grid on;
title(sprintf('[%s,%s]$\\times\\pi/24$ f %0.2f',str_pole_north,str_pole_south,kappa_sparse_f),'Interpreter','latex');
end;%for npole_south=0:n_pole_south-1;
end;%for npole_north=0:n_pole_north-1;
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
if ~isfield(parameter,'flag_kernel_qpro_d0'); parameter.flag_kernel_qpro_d0=1; end;
flag_kernel_qpro_d0=parameter.flag_kernel_qpro_d0;
if ~isfield(parameter,'flag_kernel_qpro_d1'); parameter.flag_kernel_qpro_d1=0; end;
flag_kernel_qpro_d1=parameter.flag_kernel_qpro_d1;
if ~isfield(parameter,'flag_kernel_qpro_d2'); parameter.flag_kernel_qpro_d2=0; end;
flag_kernel_qpro_d2=parameter.flag_kernel_qpro_d2;
if ~isfield(parameter,'kernel_qpro_l_max_band'); parameter.kernel_qpro_l_max_band=+Inf; end;
kernel_qpro_l_max_band=parameter.kernel_qpro_l_max_band;
if ~isfield(parameter,'kernel_qpro_polar_a_pole_north'); parameter.kernel_qpro_polar_a_pole_north=1.0*pi/12; end;
kernel_qpro_polar_a_pole_north=min(pi/2,parameter.kernel_qpro_polar_a_pole_north);
parameter.kernel_qpro_polar_a_pole_north = kernel_qpro_polar_a_pole_north;
if ~isfield(parameter,'kernel_qpro_polar_a_pole_south'); parameter.kernel_qpro_polar_a_pole_south=0.5*pi/12; end;
kernel_qpro_polar_a_pole_south=min(pi/2,parameter.kernel_qpro_polar_a_pole_south);
parameter.kernel_qpro_polar_a_pole_south = kernel_qpro_polar_a_pole_south;
if ~isfield(parameter,'kernel_qpro_deconvolution_factor_max'); parameter.kernel_qpro_deconvolution_factor_max=1024; end;
kernel_qpro_deconvolution_factor_max=parameter.kernel_qpro_deconvolution_factor_max;
if ~isfield(parameter,'kernel_qpro_MaxIterations'); parameter.kernel_qpro_MaxIterations=1024; end;
kernel_qpro_MaxIterations=parameter.kernel_qpro_MaxIterations;
%%%%;
if ~isfield(parameter,'flag_kernel_full'); parameter.flag_kernel_full=0; end;
flag_kernel_full=parameter.flag_kernel_full;
if (kernel_qpro_polar_a_pole_north + kernel_qpro_polar_a_pole_south > pi-1e-12);
flag_kernel_full = 1;
parameter.flag_kernel_full = flag_kernel_full;
end;%if (kernel_qpro_polar_a_pole_north + kernel_qpro_polar_a_pole_south > pi-1e-12);
%%%%;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

%%%%%%%%;
if isempty(l_max); l_max = 49; end;
l_val_ = transpose([0:l_max]);
%%%%%%%%;
l_max_band = min(l_max,kernel_qpro_l_max_band);
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
[a_drop_node_,a_drop_weight_] = legpts(n_a_use,[0+kernel_qpro_polar_a_pole_north,pi-kernel_qpro_polar_a_pole_south]);
[a_keep_node_north_,a_keep_weight_north_] = legpts(n_a_use,[0 ,kernel_qpro_polar_a_pole_north]);
[a_keep_node_south_,a_keep_weight_south_] = legpts(n_a_use,[pi-kernel_qpro_polar_a_pole_south,pi]);
a_keep_node_ = [a_keep_node_north_;a_keep_node_south_]; %<-- col vector. ;
a_keep_weight_ = [a_keep_weight_north_,a_keep_weight_south_]; %<-- row vector. ;
a_full_node_ = [a_keep_node_;a_drop_node_]; %<-- col vector. ;
a_full_weight_ = [a_keep_weight_,a_drop_weight_]; %<-- row vector. ;
kappa_sparse_f = sum(2*pi*a_keep_weight_*sin(a_keep_node_))./max(1e-12,sum(2*pi*a_full_weight_*sin(a_full_node_)));
%%%%%%%%;
% Note, this is normalized so that the integral of: ;
% chebleg_d_{1+l_val_}^{2} ;
% over the surface of the sphere is (2*pi), and not 1.0d0. ;
% This means that the integral of: ;
% chebleg_d_{1+0}^{1} ;
% over the surface of the sphere will be (2*pi)*sqrt(2). ;
% Thus, the surface-integral of chebfun_kernel_norm_qpro_ (below) ;
% will be (2*pi)*sqrt(2). ;
% Finally, after deconvolving by deconvolve_l_(1+0) = sqrt(4*pi), ;
% the surface-integral: ;
% sum(2*pi*a_full_weight_*(chebfun_kernel_norm_qpro_(cos(a_full_node_)).*chebleg_d_{1+l_val}(cos(a_full_node_)).*sin(a_full_node_))) * deconvolve_l_(1+l_val) / (4*pi) ;
% will equal sqrt(pi), ;
% and the surface-integral: ;
% sum(2*pi*a_full_weight_*(chebfun_kernel_norm_qpro_(cos(a_full_node_)).*sin(a_full_node_))) * deconvolve_l_(1+l_val) / (4*pi) ;
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
if flag_kernel_full==1; kappa_ = deltafunc_band_formula_; end;
%%%%%%%%%%%%%%%%;
if flag_kernel_full==0;
%%%%%%%%%%%%%%%%;
leg_drop_da__ = zeros(1+l_max,n_a_use);
for l_val=0:l_max_band; leg_drop_da__(1+l_val,:) = chebleg_d_{1+l_val}(cos(a_drop_node_)); end;%for l_val=0:l_max_band;
%%%%;
C_cc__ = zeros(1+l_max,1+l_max);
for l_val_0=0:l_max_band;
for l_val_1=0:l_max_band;
C_cc__(1+l_val_0,1+l_val_1) = ...
  2*pi ...
* a_drop_weight_ ...
* ( ...
      chebleg_d_{1+l_val_0}(cos(a_drop_node_)) ...
    .*chebleg_d_{1+l_val_1}(cos(a_drop_node_)) ...
    .*sin(a_drop_node_) ...
   ) ...
;
end;%for l_val_1=0:l_max_band;
end;%for l_val_0=0:l_max_band;
%%%%;
if flag_kernel_qpro_d1;
D_cc__ = zeros(1+l_max,1+l_max);
for l_val_0=0:l_max_band;
tmp_chebleg_0_ = chebleg_d_{1+l_val_0};
tmp_dchebleg_0_ = diff(tmp_chebleg_0_);
for l_val_1=0:l_max_band;
tmp_chebleg_1_ = chebleg_d_{1+l_val_1};
tmp_dchebleg_1_ = diff(tmp_chebleg_1_);
D_cc__(1+l_val_0,1+l_val_1) = ...
  2*pi ...
* a_drop_weight_ ...
* ( ...
      tmp_dchebleg_0_(cos(a_drop_node_)) ...
    .*tmp_dchebleg_1_(cos(a_drop_node_)) ...
    .*sin(a_drop_node_) ...
   ) ...
;
end;%for l_val_1=0:l_max_band;
end;%for l_val_0=0:l_max_band;
end;%if flag_kernel_qpro_d1;
%%%%;
if flag_kernel_qpro_d1;
E_cc__ = zeros(1+l_max,1+l_max);
for l_val_0=0:l_max_band;
tmp_chebleg_0_ = chebleg_d_{1+l_val_0};
tmp_dchebleg_0_ = diff(tmp_chebleg_0_);
tmp_ddchebleg_0_ = diff(tmp_dchebleg_0_);
for l_val_1=0:l_max_band;
tmp_chebleg_1_ = chebleg_d_{1+l_val_1};
tmp_dchebleg_1_ = diff(tmp_chebleg_1_);
tmp_ddchebleg_1_ = diff(tmp_dchebleg_1_);
E_cc__(1+l_val_0,1+l_val_1) = ...
  2*pi ...
* a_drop_weight_ ...
* ( ...
      tmp_ddchebleg_0_(cos(a_drop_node_)) ...
    .*tmp_ddchebleg_1_(cos(a_drop_node_)) ...
    .*sin(a_drop_node_) ...
   ) ...
;
end;%for l_val_1=0:l_max_band;
end;%for l_val_0=0:l_max_band;
end;%if flag_kernel_qpro_d1;
%%%%;
H_cc__ = C_cc__;
if flag_kernel_qpro_d1; H_cc__ = H_cc__ + D_cc__ ; end;
if flag_kernel_qpro_d2; H_cc__ = H_cc__ + E_cc__ ; end;
%%%%;
lb_ = zeros(1+l_max,1); %<-- lower-bound. ;
lb_(1+0) = sqrt(1+2*0)*sqrt(4*pi); %<-- ensure average is constant. ;
lb_(1+[1:l_max]) = deltafunc_band_formula_(1+[1:l_max])./max(1e-12,kernel_qpro_deconvolution_factor_max); %<-- ensure kernel_qpro_deconvolution_factor_max. ;
ub_ = zeros(1+l_max,1); %<-- lower-bound. ;
ub_(1+0) = sqrt(1+2*0)*sqrt(4*pi); %<-- ensure average is constant. ;
ub_(1+[1:l_max]) = deltafunc_band_formula_(1+[1:l_max]).*max(1e-12,kernel_qpro_deconvolution_factor_max); %<-- ensure kernel_qpro_deconvolution_factor_max. ;
tmp_H = H_cc__;
tmp_f = []; tmp_A = []; tmp_b = []; tmp_Aeq = []; tmp_beq = []; tmp_lb = lb_; tmp_ub = ub_; tmp_x0 = deltafunc_band_formula_;
if isempty(tmp_H) | fnorm(tmp_H)==0;
kappa_ = tmp_x0;
else;
kernel_qpro_optimoptions = optimoptions('quadprog');
kernel_qpro_optimoptions.MaxIterations = kernel_qpro_MaxIterations;
kernel_qpro_optimoptions.Display = 'off';
kappa_ = quadprog( ...
 tmp_H ...
,tmp_f ...
,tmp_A ...
,tmp_b ...
,tmp_Aeq ...
,tmp_beq ...
,tmp_lb ...
,tmp_ub ...
,tmp_x0 ...
,kernel_qpro_optimoptions ...
);
end;%if isempty(tmp_H) | fnorm(tmp_H)==0;
%%%%%%%%%%%%%%%%;
end;%if flag_kernel_full==0;
%%%%%%%%%%%%%%%%;

%%%%%%%%
chebfun_kernel_qpro_ = chebfun(0);
for l_val=0:l_max;
chebfun_kernel_qpro_ = chebfun_kernel_qpro_ + kappa_(1+l_val)*chebleg_d_{1+l_val};
end;%for l_val=0:l_max;
const_avg = sum(2*pi*a_full_weight_*(chebleg_d_{1+0}(cos(a_full_node_)).^1.*sin(a_full_node_)));
if (flag_verbose>0); disp(sprintf(' %% const_avg vs 2*sqrt(2)*pi %0.16f',fnorm(const_avg - 2*sqrt(2)*pi))); end;
chebfun_kernel_qpro_avg = sum(2*pi*a_full_weight_*(chebfun_kernel_qpro_(cos(a_full_node_)).^1.*sin(a_full_node_)));
chebfun_kernel_qpro_avg = chebfun_kernel_qpro_avg / max(1e-12,abs(const_avg));
if (flag_verbose>0); disp(sprintf(' %% chebfun_kernel_qpro_avg: %0.16f',chebfun_kernel_qpro_avg)); end;
chebfun_kernel_norm_qpro_ = chebfun_kernel_qpro_ / max(1e-12,chebfun_kernel_qpro_avg) ;
kappa_norm_ = kappa_ / max(1e-12,chebfun_kernel_qpro_avg) ;
%%%%%%%%;

%%%%%%%%;
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figmed;
linewidth_use = 2;
c_use__ = colormap_81s; n_c_use = size(c_use__,1);
subplot(1,3,[1,2]);
n_z = 1+1024; z_ = transpose(linspace(-1,+1,n_z));
hold on;
plot(log(1-z_),chebfun_kernel_norm_qpro_(z_),'-','Color',0.85*[1,1,1],'LineWidth',linewidth_use);
hold off;
xlabel('log(1-z)'); ylabel('P(z)'); title('chebfun_kernel_norm_qpro_','Interpreter','none'); grid on;
subplot(1,3,3);
hold on;
plot(l_val_,kappa_norm_,'x-','Color',0.85*[1,1,1],'LineWidth',3);
hold off;
xlabel('l_val','Interpreter','none');
ylabel('kappa_norm_','Interpreter','none');
title('scaling factor','Interpreter','none');
end;%if (flag_disp>1);
sum_l1_drop = sum(2*pi*a_drop_weight_*(chebfun_kernel_norm_qpro_(cos(a_drop_node_)).^1.*sin(a_drop_node_)));
sum_l1_keep = sum(2*pi*a_keep_weight_*(chebfun_kernel_norm_qpro_(cos(a_keep_node_)).^1.*sin(a_keep_node_)));
sum_l2_drop = sum(2*pi*a_drop_weight_*(chebfun_kernel_norm_qpro_(cos(a_drop_node_)).^2.*sin(a_drop_node_)));
sum_l2_keep = sum(2*pi*a_keep_weight_*(chebfun_kernel_norm_qpro_(cos(a_keep_node_)).^2.*sin(a_keep_node_)));
if (flag_verbose>0); disp(sprintf(' %% sum_l1_drop: %0.16f',sum_l1_drop)); end;
if (flag_verbose>0); disp(sprintf(' %% sum_l1_keep: %0.16f',sum_l1_keep)); end;
if (flag_verbose>0); disp(sprintf(' %% sum_l2_drop: %0.16f',sum_l2_drop)); end;
if (flag_verbose>0); disp(sprintf(' %% sum_l2_keep: %0.16f',sum_l2_keep)); end;
%%%%%%%%;

%%%%%%%%;
deconvolve_l_ = deltafunc_band_formula_./max(1e-12,kappa_norm_(1+l_val_));
deltafunc_band_from_chebfun_mollify_ = zeros(1+l_max,1);
deltafunc_band_from_chebfun_crop_mollify_ = zeros(1+l_max,1);
for l_val=0:l_max;
tmp_I = sum(2*pi*a_full_weight_*(chebfun_kernel_norm_qpro_(cos(a_full_node_)).*chebleg_d_{1+l_val}(cos(a_full_node_)).*sin(a_full_node_)));
deltafunc_band_from_chebfun_mollify_(1+l_val) = tmp_I;
tmp_I = sum(2*pi*a_keep_weight_*(chebfun_kernel_norm_qpro_(cos(a_keep_node_)).*chebleg_d_{1+l_val}(cos(a_keep_node_)).*sin(a_keep_node_)));
deltafunc_band_from_chebfun_crop_mollify_(1+l_val) = tmp_I;
end;%for l_val=0:l_max;
deltafunc_band_from_chebfun_restore_ = deltafunc_band_from_chebfun_mollify_ .* deconvolve_l_ ;
%figure(1);clf;plot(deltafunc_band_from_chebfun_restore_- 2*pi*deltafunc_band_formula_ , '.');return;
deltafunc_band_from_chebfun_crop_restore_ = deltafunc_band_from_chebfun_crop_mollify_ .* deconvolve_l_ ;
relative_error_full_ = abs(deltafunc_band_formula_ - deltafunc_band_from_chebfun_restore_/(2*pi))./abs(deltafunc_band_formula_);
relative_error_crop_ = abs(deltafunc_band_formula_ - deltafunc_band_from_chebfun_crop_restore_/(2*pi))./abs(deltafunc_band_formula_);
%%%%%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;
markersize_use = 16;
hold on;
plot(l_val_,log10(relative_error_full_),'k.','MarkerSize',markersize_use);
plot(l_val_,log10(relative_error_crop_),'r.','MarkerSize',markersize_use);
hold off;
xlabel('l_val_','Interpreter','none');
ylabel('log10(abs(form-reco)./abs(form))','Interpreter','none');
xlim([-1,1+l_max]); set(gca,'XTick',l_val_(1:7:end));
ylim([-16,+1]); set(gca,'YTick',-15:0);
grid on;
title(sprintf('f %0.6f',kappa_sparse_f),'Interpreter','none');
end;%if flag_disp;
%%%%%%%%;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

