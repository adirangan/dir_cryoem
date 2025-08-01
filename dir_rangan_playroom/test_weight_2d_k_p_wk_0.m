function ...
[ ...
 parameter ...
] = ...
test_weight_2d_k_p_wk_0( ...
 parameter ...
);

str_thisfunction = 'test_weight_2d_k_p_wk_0';

if nargin<1;
disp(sprintf(' %% testing %s',str_thisfunction));
%%%%;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
%%%%;
flag_verbose = 1; %<-- verbosity level. ;
flag_disp=1; %<-- display level. ;
flag_replot = 1;
%%%%;
k_int = 48; %<-- highest frequency (2*pi*k_p_r_max), watch out for aliasing! ;
k_eq_d_double = 0.5; %<-- prefactor for k_eq_d, determines density of sampling in frequency-space. ;
n_w_int = 2.0; %<-- prefactor for n_w_max, determines the number of distinct angles (i.e., n_gamma_z) used in frequency-space 2d-polar-grid. ;
%%%%;
dir_base = sprintf('/%s/rangan/dir_cryoem/dir_rangan_playroom',string_root);
dir_jpg = sprintf('%s/dir_test_transforms_jpg',dir_base);
if (~exist(dir_jpg,'dir')); disp(sprintf(' %% mkdir %s',dir_jpg)); mkdir(dir_jpg); end;
%%%%;
parameter = struct('type','parameter');
parameter.flag_verbose=flag_verbose;
parameter.flag_disp=flag_disp;
parameter.flag_replot=flag_replot;
parameter.k_int=k_int;
parameter.k_eq_d_double=k_eq_d_double;
parameter.n_w_int=n_w_int;
parameter.dir_base=dir_base;
parameter.dir_jpg=dir_jpg;
[ ...
 parameter ...
] = ...
test_weight_2d_k_p_wk_0( ...
 parameter ...
);
%%%%;
disp('returning'); return;
end;%if nargin<1;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
flag_verbose=parameter.flag_verbose;
if ~isfield(parameter,'flag_disp'); parameter.flag_disp=0; end;
flag_disp=parameter.flag_disp; nf=0;
if ~isfield(parameter,'flag_replot'); parameter.flag_replot=0; end;
flag_replot=parameter.flag_replot;
if ~isfield(parameter,'k_int'); parameter.k_int=32; end;
k_int=parameter.k_int;
if ~isfield(parameter,'k_eq_d_double'); parameter.k_eq_d_double=0.5; end;
k_eq_d_double=parameter.k_eq_d_double;
if ~isfield(parameter,'n_w_int'); parameter.n_w_int=1.0; end;
n_w_int=parameter.n_w_int;
if ~isfield(parameter,'dir_base'); parameter.dir_base=pwd; end;
dir_base=parameter.dir_base;
if ~isfield(parameter,'dir_jpg'); parameter.dir_jpg=pwd; end;
dir_jpg=parameter.dir_jpg;
if ~isfield(parameter,'frequency_max'); parameter.frequency_max=3.0; end;
frequency_max=parameter.frequency_max;
if ~isfield(parameter,'n_frequency'); parameter.n_frequency=1+3*8; end;
n_frequency=parameter.n_frequency;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

%%%%%%%%;
% tests 2d quadrature. ;
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% First define integral of <f,f>. ;
%%%%%%%%;
h2d_ = @(kd) 4*pi^2*(besselj(0,kd) + besselj(2,kd)); % calculates <f_j,f_k>, normalized so that <f,f> = (4*pi^2);
dh2d_ = @(kd) 4*pi^3*(besselj(-1,kd) - besselj(+3,kd));

%%%%%%%%;
% Now set up k-quadrature on sphere. ;
%%%%%%%%;
k_p_r_max = k_int/(2*pi); k_eq_d = k_eq_d_double/(2*pi); str_L = 'L';
[ ...
 n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
] = ...
get_weight_3d_1( ...
 0*flag_verbose ...
,k_p_r_max ...
,k_eq_d ...
,str_L ...
);
%%%%%%%%;
% Now set up and test polar-quadrature-weights on disk. ;
%%%%%%%%;
l_max_upb = round(2*pi*k_p_r_max);
l_max_max = min(l_max_upb,1+ceil(2*pi*k_p_r_(end)));
n_w_max = n_w_int*2*(l_max_max+1); n_w_0in_ = n_w_max*ones(n_k_p_r,1);
[ ...
 n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_k_p_wk_ ...
,k_p_r_wk_ ...
,k_p_w_wk_ ...
,k_c_0_wk_ ...
,k_c_1_wk_ ...
] = ...
get_weight_2d_2( ...
 0*flag_verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,-1 ...
,n_w_0in_ ...
,weight_3d_k_p_r_ ...
);
n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_]);
n_gamma_z = n_w_max;
gamma_z_ = transpose(linspace(0,2*pi,n_gamma_z+1)); gamma_z_ = gamma_z_(1:n_gamma_z);
%%%%%%%%;
if (flag_verbose>0);
disp(sprintf(' %% 2*pi*k_p_r_max %+0.2f <-- k_p_r_max +%0.6f',2*pi*k_p_r_max,k_p_r_max));
disp(sprintf(' %% n_k_p_r %d',n_k_p_r));
disp(sprintf(' %% n_w_max %d',n_w_max));
end;%if (flag_verbose>0);
%%%%%%%%;
frequency_ = transpose(linspace(0,frequency_max,n_frequency));
errrel_f_ = zeros(n_frequency,1);
if flag_disp>0; figure(1+nf);nf=nf+1;clf;figbig;p_row=4;p_col=ceil(n_frequency/p_row);np=0; end;
for nfrequency=0:n_frequency-1;
frequency = frequency_(1+nfrequency);
disp(sprintf(' %% nfrequency %.2d/%.2d %0.6f',nfrequency,n_frequency,frequency));
delta_T = frequency;
omega_T = 2*pi*rand();
delta_T_0 = delta_T * cos(omega_T);
delta_T_1 = delta_T * sin(omega_T);
T_k_p_wk_ = exp(+2*pi*i*(k_c_0_wk_*delta_T_0 + k_c_1_wk_*delta_T_1));
I_quad = sum(T_k_p_wk_.*weight_2d_k_p_wk_)*(2*pi)^2;
T_kd = 2*pi*k_p_r_max*delta_T;
I_form = h2d_(T_kd) / (2*pi)^2 * pi*k_p_r_max^2;
tmp_errel = fnorm_disp(flag_verbose,'I_form',I_form,'I_quad',I_quad);
errel_f_(1+nfrequency) = tmp_errel;
if flag_disp>0; 
subplot(p_row,p_col,1+np);np=np+1;cla;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(T_k_p_wk_)); axis image; axisnotick;
title(sprintf('f %+0.1f <-- log10(errrel) %+0.1f',frequency,log10(tmp_errel)));
drawnow();
end;%if flag_disp>0; 
end;%for nfrequency=0:n_frequency-1;
%%%%%%%%;
if flag_disp>0;
fname_fig_pre = sprintf('%s/test_transforms_test_weight_2d_k_p_wk_FIGB',dir_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
%close(gcf);
end;%if (flag_disp>0);
%%%%%%%%;

%%%%%%%%;
if (flag_disp>0);
figure(1+nf);nf=nf+1;clf;figsml;
plot(frequency_,min(+1,log10(errel_f_)),'ko','MarkerSize',16,'MarkerFaceColor','r');
xlim([0,frequency_max]); xlabel('frequency');
ylim([-16,+1]); ylabel('errrel');
%%%%;
fname_fig_pre = sprintf('%s/test_transforms_test_weight_2d_k_p_wk_FIGA',dir_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
%close(gcf);
end;%if (flag_disp>0);
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

return;





