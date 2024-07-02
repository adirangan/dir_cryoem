function ...
[ ...
 parameter ...
,polar_a_St__ ...
,azimu_b_St__ ...
] = ...
sphere_point_cloud_dynamics_1( ...
 parameter ...
,n_S ...
,polar_a_S_ ...
,azimu_b_S_ ...
);

str_thisfunction = 'sphere_point_cloud_dynamics_1';

%%%%%%%%;
if nargin<1;
disp(sprintf(' %% testing %s',str_thisfunction));
flag_disp = 1; nf=0;
%%%%;
hist2dab_k_eq_d = 0.5/(48/(2*pi));
[ ...
 n_hist2dab_S ...
,hist2dab_azimu_b_S_ ...
,hist2dab_polar_a_S_ ...
,hist2dab_weight_S_ ...
,hist2dab_k_c_0_S_ ...
,hist2dab_k_c_1_S_ ...
,hist2dab_k_c_2_S_ ...
] = ...
sample_shell_6( ...
 1.0 ...
,hist2dab_k_eq_d ...
,'L' ... %<-- exclude pole. ;
,0 ... %<-- adaptive grid. ;
) ;
%%%%;
n_S = 1024*32;
s_ = linspace(0,1,1+n_S); s_ = transpose(s_(1:n_S));
polar_a_S_ = periodize(7*pi*s_,0,1*pi);
azimu_b_S_ = periodize(31*pi*(s_-0.5),0,2*pi);
[polar_a_S_,azimu_b_S_] = periodize_polar_a_azimu_b_0(polar_a_S_,azimu_b_S_);
[k_0_S_,k_1_S_,k_2_S_] = local_sphere_012_from_ab(polar_a_S_,azimu_b_S_);
%%%%;
T_MAX = 0.50;
n_T = 64;
parameter = struct('type','parameter');
parameter.T_MAX = T_MAX;
parameter.n_T = n_T;
[ ...
 parameter ...
,polar_a_St__ ...
,azimu_b_St__ ...
] = ...
sphere_point_cloud_dynamics_1( ...
 parameter ...
,n_S ...
,polar_a_S_ ...
,azimu_b_S_ ...
);
[k_0_St__,k_1_St__,k_2_St__] = local_sphere_012_from_ab(polar_a_St__,azimu_b_St__);
%%%%;
if flag_disp>1;
figure(1+nf);nf=nf+1;clf;figbig;
markersize_use = 8;
c_use__ = colormap_81s(); n_c_use = size(c_use__,1);
p_row = 3; p_col = 3; np=0;
nT_ = round(linspace(0,n_T-1,p_row*p_col));
for np=0:p_row*p_col-1;
subplot(p_row,p_col,1+np);
nT = nT_(1+np);
nc_use = max(0,min(n_c_use-1,floor(n_c_use*nT/max(1,n_T-1))));
plot_sphere_grid_0;
hold on;
plot3(k_0_St__(:,1+nT),k_1_St__(:,1+nT),k_2_St__(:,1+nT),'ko','MarkerSize',markersize_use,'MarkerFaceColor',c_use__(1+nc_use,:));
hold off;
axis(1.25*[-1,+1,-1,+1,-1,+1]);
axisnotick3d;
axis vis3d;
title(sprintf('nT %d/%d',nT,n_T));
drawnow();
end;%for np=0:p_row*p_col-1;
end;% if flag_disp>1;
%%%%;
if flag_disp>0;
figure(1+nf);nf=nf+1;clf;figbig;
markersize_use = 8;
p_row = 3; p_col = 3; np=0;
nT_ = round(linspace(0,n_T-1,p_row*p_col));
for np=0:p_row*p_col-1;
subplot(p_row,p_col,1+np);
nT = nT_(1+np);
hlim_ = [];
c_use__ = colormap_81s;
flag_2d_vs_3d = 0;
flag_loghist_vs_hist = 0;
[ ...
 h_raw_ab_ ...
 h_w3d_ab_ ...
] = ...
hist2d_polar_a_azimu_b_0( ...
 hist2dab_polar_a_S_ ...
,hist2dab_azimu_b_S_ ...
,hist2dab_weight_S_ ...
,polar_a_St__(:,1+nT) ...
,azimu_b_St__(:,1+nT) ...
,hlim_ ...
,c_use__ ...
,flag_2d_vs_3d ...
,flag_loghist_vs_hist ...
);
if flag_2d_vs_3d==0; axis(1.25*[-1,+1,-1,+1,-1,+1]); axisnotick3d; end;
if flag_2d_vs_3d==1; axis image; axisnotick; end;
title(sprintf('nT %d/%d',nT,n_T));
drawnow();
end;%for np=0:p_row*p_col-1;
end;%if flag_disp>0;

%%%%;
disp('returning'); return;
end;%if nargin<1;
%%%%%%%%;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_S=[]; end; na=na+1;
if (nargin<1+na); polar_a_S_=[]; end; na=na+1;
if (nargin<1+na); azimu_b_S_=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
flag_verbose=parameter.flag_verbose;
if ~isfield(parameter,'T_MAX'); parameter.T_MAX=1.0; end;
T_MAX=parameter.T_MAX;
if ~isfield(parameter,'n_T'); parameter.n_T=32; end;
n_T=parameter.n_T;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;
sigma_t = sqrt(T_MAX/max(1,n_T));
polar_a_St__ = zeros(n_S,n_T);
azimu_b_St__ = zeros(n_S,n_T);
nT=0;
polar_a_St__(:,1+nT) = polar_a_S_;
azimu_b_St__(:,1+nT) = azimu_b_S_;
for nT=1:n_T-1;
[k_0_S_,k_1_S_,k_2_S_] = local_sphere_012_from_ab(polar_a_St__(:,1+nT-1),azimu_b_St__(:,1+nT-1));
W_0_S_ = sigma_t*randn(n_S,1);
W_1_S_ = sigma_t*randn(n_S,1);
W_2_S_ = sigma_t*randn(n_S,1);
k_0_S_ = k_0_S_ + W_0_S_;
k_1_S_ = k_1_S_ + W_1_S_;
k_2_S_ = k_2_S_ + W_2_S_;
[polar_a_S_,azimu_b_S_] = local_sphere_ab_from_012(k_0_S_,k_1_S_,k_2_S_);
[polar_a_S_,azimu_b_S_] = periodize_polar_a_azimu_b_0(polar_a_S_,azimu_b_S_);
polar_a_St__(:,1+nT) = polar_a_S_;
azimu_b_St__(:,1+nT) = azimu_b_S_;
end;%for nT=1:n_T-1;
if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

function [k_0_,k_1_,k_2_] = local_sphere_012_from_ab(polar_a_,azimu_b_);
k_0_ = sin(polar_a_).*cos(azimu_b_);
k_1_ = sin(polar_a_).*sin(azimu_b_);
k_2_ = cos(polar_a_);

function [polar_a_,azimu_b_] = local_sphere_ab_from_012(k_0_,k_1_,k_2_);
k_r_ = sqrt(k_0_.^2 + k_1_.^2 + k_2_.^2);
k_01_r_ = sqrt(k_0_.^2 + k_1_.^2);
azimu_b_ = atan2(k_1_,k_0_);
polar_a_ = atan2(k_01_r_,k_2_);

