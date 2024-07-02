function ...
[ ...
 parameter ...
,polar_a_St__ ...
,azimu_b_St__ ...
] = ...
sphere_point_cloud_dynamics_0( ...
 parameter ...
,n_S ...
,polar_a_S_ ...
,azimu_b_S_ ...
);

str_thisfunction = 'sphere_point_cloud_dynamics_0';

%%%%%%%%;
if nargin<1;
disp(sprintf(' %% testing %s',str_thisfunction));
flag_disp = 1; nf=0;
n_S = 16;
s_ = linspace(0,1,1+n_S); s_ = transpose(s_(1:n_S));
polar_a_S_ = periodize(4*pi*s_.^2,0,1*pi);
azimu_b_S_ = periodize(16*pi*(s_-0.5).^3,0,2*pi);
[k_0_S_,k_1_S_,k_2_S_] = local_sphere_012_from_ab(polar_a_S_,azimu_b_S_);
%%%%;
figure(1+nf);nf=nf+1;clf;figsml;
markersize_use = 12;
plot_sphere_grid_0;
hold on;
plot3(k_0_S_,k_1_S_,k_2_S_,'ko','MarkerSize',markersize_use,'MarkerFaceColor','g');
hold off;
axis(1.25*[-1,+1,-1,+1,-1,+1]);
axisnotick3d;
axis vis3d;
drawnow();
%%%%;
parameter = struct('type','parameter');
parameter.T_MAX = 0.10;
[ ...
 parameter ...
,polar_a_St__ ...
,azimu_b_St__ ...
] = ...
sphere_point_cloud_dynamics_0( ...
 parameter ...
,n_S ...
,polar_a_S_ ...
,azimu_b_S_ ...
);
[k_0_St__,k_1_St__,k_2_St__] = local_sphere_012_from_ab(polar_a_St__,azimu_b_St__);
%%%%;
n_t = size(polar_a_St__,2);
disp(sprintf(' %% n_t: %0.3d',n_t));
%%%%;
figure(1+nf);nf=nf+1;clf;figsml;
markersize_start = 16;
markersize_final = 4;
markersize_use_ = linspace(markersize_start,markersize_final,n_t);
plot_sphere_grid_0;
c_use__ = colormap_81s(); n_c_use = size(c_use__,1);
hold on;
for nt=0:n_t-1;
nc_use = max(0,min(n_c_use-1,floor(n_c_use*nt/max(1,n_t-1))));
markersize_use = markersize_use_(1+nt);
plot3(k_0_St__(:,1+nt),k_1_St__(:,1+nt),k_2_St__(:,1+nt),'ko','MarkerSize',markersize_use,'MarkerFaceColor',c_use__(1+nc_use,:));
end;%for nt=0:n_t-1;
hold off;
axis(1.25*[-1,+1,-1,+1,-1,+1]);
axisnotick3d;
axis vis3d;
drawnow();
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

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;
n_Y = 2*n_S;
Y_0_ = reshape([polar_a_S_,azimu_b_S_],[n_Y,1]);
[t_,Y_t__] = ode15s(@local_ODEFUN,[0,T_MAX],Y_0_);
n_t = numel(t_);
polar_a_St__ = periodize(reshape(Y_t__(1+0:2:end),[n_S,n_t]),0,1*pi);
azimu_b_St__ = periodize(reshape(Y_t__(1+1:2:end),[n_S,n_t]),0,2*pi);
if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

function dY_ = local_ODEFUN(t,Y_);
n_Y = size(Y_,1);
n_S = floor(n_Y/2);
polar_a_S_ = Y_(1+0:2:end);
azimu_b_S_ = Y_(1+1:2:end);
[k_0_S_,k_1_S_,k_2_S_] = local_sphere_012_from_ab(polar_a_S_,azimu_b_S_);
k_S3__ = [k_0_S_,k_1_S_,k_2_S_];
K_0_SS__ = bsxfun(@minus,reshape(k_0_S_,[n_S,1]),reshape(k_0_S_,[1,n_S]));
K_1_SS__ = bsxfun(@minus,reshape(k_1_S_,[n_S,1]),reshape(k_1_S_,[1,n_S]));
K_2_SS__ = bsxfun(@minus,reshape(k_2_S_,[n_S,1]),reshape(k_2_S_,[1,n_S]));
KK_SS__ = K_0_SS__.^2 + K_1_SS__.^2 + K_2_SS__.^2;
dK_0_S_ = sum(+K_0_SS__./max(1e-12,KK_SS__),2);
dK_1_S_ = sum(+K_1_SS__./max(1e-12,KK_SS__),2);
dK_2_S_ = sum(+K_2_SS__./max(1e-12,KK_SS__),2);
dK_S3__ = [dK_0_S_,dK_1_S_,dK_2_S_];
dK_S3__ = dK_S3__ - k_S3__.*(sum(bsxfun(@times,k_S3__,dK_S3__)));
dK_0_S_ = dK_S3__(:,1+0);
dK_1_S_ = dK_S3__(:,1+1);
dK_2_S_ = dK_S3__(:,1+2);
k_01_S_ = sqrt(k_0_S_.^2 + k_1_S_.^2);
dK_01_S_ = (k_0_S_.*dK_0_S_ + k_1_S_.*dK_1_S_)./max(1e-12,sqrt(k_0_S_.^2 + k_1_S_.^2));
%% Note that: datan(y/x) = 1/(1+y^2/x^2) * (y'/x - yx'/x^2) = (x^2)/(x^2+y^2)*(y'x - yx')/(x^2) = (y'x-yx')/(x^2+y^2) ;
dazimu_b_S_ = (dK_1_S_.*k_0_S_ - k_1_S_.*dK_0_S_)./max(1e-12,k_0_S_.^2 + k_1_S_.^2);
dpolar_a_S_ = (dK_01_S_.*k_2_S_ - k_01_S_.*dK_2_S_)./max(1e-12,k_2_S_.^2 + k_01_S_.^2);
dY_ = reshape([dpolar_a_S_,dazimu_b_S_],[n_Y,1]);

function distance = local_sphere_distance_0(polar_a_0,azimu_b_0,polar_a_1,azimu_b_1);
[k_0_0_,k_1_0_,k_2_0_] = local_sphere_012_from_ab(polar_a_0_,azimu_b_0_);
[k_0_1_,k_1_1_,k_2_1_] = local_sphere_012_from_ab(polar_a_1_,azimu_b_1_);
distance = sqrt( (k_0_0-k_0_1).^2 + (k_1_0-k_1_1).^2 + (k_2_0-k_2_1).^2 );

function [k_0_,k_1_,k_2_] = local_sphere_012_from_ab(polar_a_,azimu_b_);
k_0_ = sin(polar_a_).*cos(azimu_b_);
k_1_ = sin(polar_a_).*sin(azimu_b_);
k_2_ = cos(polar_a_);

function [polar_a_,azimu_b_] = local_sphere_ab_from_012(k_0_,k_1_,k_2_);
k_r_ = sqrt(k_0_.^2 + k_1_.^2 + k_2_.^2);
k_01_r_ = sqrt(k_0_.^2 + k_1_.^2);
azimu_b_ = atan2(k_1_,k_0_);
polar_a_ = atan2(k_01_r_,k_2_);

