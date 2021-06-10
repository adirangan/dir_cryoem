function ...
[ ...
 image_gamma_z_out_ ...
,image_delta_out_dM__ ...
,image_sigma_gamma_z_out_ ...
,image_sigma_delta_out_ ...
,image_flag_update_out_ ...
,fnorm_d_image_gamma_z_out_ ...
,fnorm_d_image_delta_out_ ...
] = ...
mp_init_align_gamma_z_delta_0( ...
 n_M ...
,k_c_Mk__ ...
,d_gamma_z_MM__ ...
,d_delta_dMM___ ...
,image_flag_update_0in_ ...
,image_gamma_z_0in_ ...
,image_delta_0in_dM__ ...
,n_iteration ...
,eta_step ...
);

if nargin<1;
verbose=1;
if (verbose); disp(sprintf(' %% [testing mp_init_align_gamma_z_delta_0]')); end;
k_p_r_max = 1; k_eq_d = (2*pi)/64; TorL = 'L';
[ ...
 n_all ...
,azimu_b_all_ ...
,polar_a_all_ ...
,weight_all_ ...
,k_c_0_all_ ...
,k_c_1_all_ ...
,k_c_2_all_ ...
,n_polar_a ...
,polar_a_ ...
,n_azimu_b_ ...
] = ...
sample_shell_5( ...
 k_p_r_max ...
,k_eq_d ...
,TorL ...
) ;
n_S = n_all; n_M = n_S;
k_c_Mk__ = [ k_c_0_all_ , k_c_1_all_ , k_c_2_all_ ];

%%%%%%%%;
% Assume true molecule is a sum of plane-waves: ;
% F = \sum_{j} exp(+2*pi*i*dot(k_c_,delta_true__(:,1+nj))) ;
% Now consider a template taken from viewing_polar_a and viewing_azimu_b has values: ;
% let sa and ca be sin(polar_a) and cos(polar_a), respectively. ;
% let sb and cb be sin(azimu_b) and cos(azimu_b), respectively. ;
% let sc and cc be sin(gamma_z) and cos(gamma_z), respectively. ;
% And rotation by azimu_b about the +z-axis is represented as: ;
% Rz__(azimu_b) = ;
% [ +cb -sb 0 ] ;
% [ +sb +cb 0 ] ;
% [  0   0  1 ] ;
% And rotation by polar_a about the +y-axis is represented as: ;
% Ry__(polar_a) = ;
% [ +ca 0 +sa ] ;
% [  0  1  0  ] ;
% [ -sa 0 +ca ] ;
% And rotation by gamma_z about the +z-axis is represented as: ;
% Rz__(gamma_z) = ;
% [ +cc -sc 0 ] ;
% [ +sc +cc 0 ] ;
% [  0   0  1 ] ;
%%%%%%%% ;
% Now the template S_ =: S_(vewing_polar_a,viewing_azimu_b) will have values given by: ;
% S_(inplane_gamma_z) = \sum_{j} exp(+2*pi*i*dot( Rz__(viewing_azimu_b) * Ry__(viewing_polar_a) * Rz__(inplane_gamma_z) * [1;0;0] , delta_true__(:,1+nj))) ;
% S_(inplane_gamma_z) = \sum_{j} exp(+2*pi*i*dot( Rz__(inplane_gamma_z) * [1;0;0] , Ry_(-viewing_polar_a) * Rz_(-viewing_azimu_b) * delta_true__(:,1+nj))) ;
% corresponding to a sum of 2d-plane-waves with frequencies: ;
% Ry_(-viewing_polar_a) * Rz_(-viewing_azimu_b) * delta_true__(:,1+nj) ;
% projected onto the 01-plane: ;
% [ +ca 0 -sa ] [ +cb +sb 0 ] [d0]   [ +ca*cb  +ca*sb -sa ] [d0]   [ +ca*cb*d0 + ca*sb*d1 - sa*d2 ] ;
% [  0  1  0  ] [ -sb +cb 0 ] [d1] = [ -sb     -cb     0  ] [d1] = [ -sb*d0 -cb*d1 + 0*d2         ] ;
% [ +sa 0 +ca ] [  0   0  1 ] [d2]   [ +sa*cb  +sa*sb +ca ] [d2]   [ +sa*cb*d0 + sa*sb*d1 + ca*d2 ] ;
% producing a 2d-frequency: ;
% [ +ca*cb*d0 + ca*sb*d1 - sa*d2 ] ;
% [ -sb*d0 -cb*d1 + 0*d2         ] ;
% directed at angle: ;
% omega = atan2( -sb*d0 -cb*d1 , +ca*cb*d0 + ca*sb*d1 - sa*d2 );
%%%%%%%%;
% The correlation matrix X__(1+nS0,1+nS1) with parameters d_gamma_z and d_delta_ corresponds to: ;
% dot( Rz__(+d_gamma_z) * S_0 , T__(+d_delta_) * S_1 ) ;
% Thus, because there is no displacement to the templates, ;
% two nearby templates S_0 and S_1 built from the same single plane-wave will have an alignment given by: ;
% d_delta_ = [0;0] and ;
% d_gamma_z = +omega_1 - omega_0 ;
%%%%%%%%;
% An image M with image_gamma_z and image_delta_ will be constructed as: ;
% M_ = T__(-image_delta_) * Rz__(+image_gamma_z) * S_ = Rz__(+image_gamma_z) * T__(Rz__(-image_gamma_z)*-image_delta_) * S_ . ;
% with T__(delta_) := exp(-2*pi*i*dot(k_c_ , delta_)). ;
%%%%%%%%;
% Thus, two nearby images M_0 and M_1 with similar viewing-angles will have an alignment given by: ;
% X = dot( Rz__(+d_gamma_z) * T__(-image_delta_0_) * Rz__(+image_gamma_z_0) * S_0 ,    ;
%          T__(+d_delta_) * T__(-image_delta_1_) * Rz__(+image_gamma_z_1) * S_1     ). ;
% X = dot( R(+dg) * T(-d0_) * R(+g0) * S_0 ,    ;
%          T(+dd_) * T(-d1_) * R(+g1) * S_1     ). ;
% X = dot( R(+dg) * R(+g0) * T(R(-g0)*-d0_) * S_0 ,    ;
%          T(+dd_) * T(-d1_) * R(+g1) * S_1     ). ;
% X = dot( R(+dg+g0) * T(R(-g0)*-d0_) * S_0 ,    ;
%          T(+dd_-d1_) * R(+g1) * S_1     ). ;
% X = dot( T( R(+dg+g0) * R(-g0) * -d0_ ) * R(+dg+g0) * S_0 ,    ;
%          T(+dd_-d1_) * R(+g1) * S_1     ). ;
% So for the translation: ;
% R(+dg+g0) * R(-g0) * -d0_ = +dd_ - d1_ ;
% dd_ = d1_ - R(+dg) * d0_ ;
% And for the rotation: ;
% dg + g0 + omega_0 = g1 + omega_1 ;
% dg = +g1 + omega_1 - g0 - omega_0 ;
%%%%%%%%;
% A consequence of these formulae is: ;
% Given two nearby images, as well as the dg and dd_ between them, ;
% we can determine d1_ and g1_ from d0_ and g0_ as follows: ;
% g1 \approx g0 + dg ; %<-- assuming that omega_1 - omega_0 is negligible, ;
% d1_ = dd_ + R(+dg) * d0_ ; %<-- no explicit error. ;
%%%%%%%%;

%%%%%%%%;
% Now set up links between nodes. ;
%%%%%%%%;
rng(0);
sigma_delta = 0.05; %<-- not too large. ;
%delta_true_ = [0;0;1]; if (verbose); disp(sprintf(' %% omega == 0')); end; %<-- omega_ will be identically 0 for all templates. ;
delta_true_ = [0;1;0]; if (verbose); disp(sprintf(' %% omega ~= 0')); end; %<-- omega_ will not be identically 0 for all templates. ;
ca_ = cos(polar_a_all_); sa_ = sin(polar_a_all_);
cb_ = cos(azimu_b_all_); sb_ = sin(azimu_b_all_);
omega_ = atan2( -sb_.*delta_true_(1+0) - cb_.*delta_true_(1+1) , +ca_.*cb_.*delta_true_(1+0) + ca_.*sb_.*delta_true_(1+1) - sa_.*delta_true_(1+2) );
omega_ = periodize(omega_,-pi,+pi);
%omega_ = zeros(n_M,1); %<-- no shift in frequency. ;
image_gamma_z_true_ = 2*pi*rand(n_M,1);
image_delta_true_dM__ = sigma_delta * randn(2,n_M);
d_gamma_z_MM__ = zeros(n_M,n_M);
d_delta_dMM___ = zeros(2,n_M,n_M);
for nM0=0:n_M-1;
omega_0 = omega_(1+nM0);
image_gamma_z_true_0 = image_gamma_z_true_(1+nM0);
image_delta_true_0_ = image_delta_true_dM__(:,1+nM0);
for nM1=0:n_M-1;
omega_1 = omega_(1+nM1);
image_gamma_z_true_1 = image_gamma_z_true_(1+nM1);
image_delta_true_1_ = image_delta_true_dM__(:,1+nM1);
d_gamma_z = periodize( image_gamma_z_true_1 + omega_1 - image_gamma_z_true_0 - omega_0 , -pi , +pi );
Rdg__ = [ +cos(d_gamma_z) , -sin(d_gamma_z) ; +sin(d_gamma_z) , +cos(d_gamma_z) ];
d_delta_ = image_delta_true_1_ - Rdg__*image_delta_true_0_;
d_gamma_z_MM__(1+nM0,1+nM1) = d_gamma_z;
d_delta_dMM___(:,1+nM0,1+nM1) = d_delta_;
end;%for nM1=0:n_M-1;
end;%for nM0=0:n_M-1;

%%%%%%%%;
% Now implement simple message passing, starting with one correctly aligned image. ;
%%%%%%%%;
image_flag_update_ = zeros(n_M,1);
image_gamma_z_0est_ = zeros(n_M,1);
image_delta_0est_dM__ = zeros(2,n_M);
%%%%;
nM0 = floor(3*n_M/4);
image_gamma_z_0est_(1+nM0) = image_gamma_z_true_(1+nM0);
image_delta_0est_dM__(:,1+nM0) = image_delta_true_dM__(:,1+nM0);
image_flag_update_(1+nM0) = 1;
%%%%;
flag_check=1;
if flag_check;
if (verbose); disp(sprintf(' %% introducing frustration')); end;
nM0 = floor(1*n_M/4);
image_gamma_z_0est_(1+nM0) = image_gamma_z_true_(1+nM0) + 1.0*pi ; %<-- opposite orientation. ;
image_delta_0est_dM__(:,1+nM0) = image_delta_true_dM__(:,1+nM0) + sigma_delta * [+1;-1]/sqrt(2); %<-- different translation. ;
image_flag_update_(1+nM0) = 1;
end;%if flag_check;
%%%%;
[ ...
 image_gamma_z_0est_ ...
,image_delta_0est_dM__ ...
,image_sigma_gamma_z_0est_ ...
,image_sigma_delta_0est_ ...
,image_flag_update_ ...
] = ...
mp_init_align_gamma_z_delta_0( ...
 n_M ...
,k_c_Mk__ ...
,d_gamma_z_MM__ ...
,d_delta_dMM___ ...
,image_flag_update_ ...
,image_gamma_z_0est_ ...
,image_delta_0est_dM__ ...
);

flag_plot=1;
if flag_plot;
figure(1);clf;
markersize_use = 32;
if (verbose); disp(sprintf(' %% n_all %d',n_all)); end;
%{
%%%%%%%%;
subplot(2,2,1);
c3d__ = colorpsphere_8points_gamma(polar_a_all_,azimu_b_all_);
scatter3(k_c_0_all_,k_c_1_all_,k_c_2_all_,markersize_use,c3d__,'filled'); axis vis3d;
%%%%%%%%;
subplot(2,2,2);
c2d__ = colormap_polar_a_azimu_b_2d(polar_a_all_,azimu_b_all_,0.35);
hold on;
patch(2*pi*[+0;+1;+1;+0],1*pi*[+0;+0;+1;+1],0.65*[1,1,1],'EdgeColor','none');
scatter(0*pi+azimu_b_all_,1*pi-polar_a_all_,markersize_use,c2d__,'filled');
hold off;
axisnotick;
axis([0,2*pi,0,1*pi]); axis equal; xlabel('azimu_b'); ylabel('polar_a'); title('true view');
%%%%%%%%;
%}
subplot(2,2,1);
hold on;
c__ = colormap_beach(); n_c = size(c__,1);
nc_ = max(0,min(n_c-1,floor(n_c*abs(periodize(image_gamma_z_0est_ - image_gamma_z_true_,-pi,+pi))/pi)));
patch(2*pi*[+0;+1;+1;+0],1*pi*[+0;+0;+1;+1],0.65*[1,1,1],'EdgeColor','none');
scatter(0*pi+azimu_b_all_,1*pi-polar_a_all_,markersize_use,c__(1+nc_,:),'filled');
scatter(0*pi+azimu_b_all_(1+nM0),1*pi-polar_a_all_(1+nM0),markersize_use,[0,0,0]);
hold off;
axisnotick;
axis([0,2*pi,0,1*pi]); axis equal; xlabel('azimu_b'); ylabel('polar_a'); title('d_gamma_z_true','Interpreter','none');
%%%%%%%%;
subplot(2,2,2);
hold on;
c__ = colormap_beach(); n_c = size(c__,1);
nc_ = max(0,min(n_c-1,floor(n_c*sqrt(sum((image_delta_0est_dM__ - image_delta_true_dM__).^2,1))/sigma_delta)));
patch(2*pi*[+0;+1;+1;+0],1*pi*[+0;+0;+1;+1],0.65*[1,1,1],'EdgeColor','none');
scatter(0*pi+azimu_b_all_,1*pi-polar_a_all_,markersize_use,c__(1+nc_,:),'filled');
scatter(0*pi+azimu_b_all_(1+nM0),1*pi-polar_a_all_(1+nM0),markersize_use,[0,0,0]);
hold off;
axisnotick;
axis([0,2*pi,0,1*pi]); axis equal; xlabel('azimu_b'); ylabel('polar_a'); title('d_delta_true','Interpreter','none');
%%%%%%%%;
subplot(2,2,3);
hold on;
c__ = colormap_beach(); n_c = size(c__,1);
nc_ = max(0,min(n_c-1,floor(n_c*image_sigma_gamma_z_0est_/pi)));
patch(2*pi*[+0;+1;+1;+0],1*pi*[+0;+0;+1;+1],0.65*[1,1,1],'EdgeColor','none');
scatter(0*pi+azimu_b_all_,1*pi-polar_a_all_,markersize_use,c__(1+nc_,:),'filled');
scatter(0*pi+azimu_b_all_(1+nM0),1*pi-polar_a_all_(1+nM0),markersize_use,[0,0,0]);
hold off;
axisnotick;
axis([0,2*pi,0,1*pi]); axis equal; xlabel('azimu_b'); ylabel('polar_a'); title('sigma_gamma_z_0est','Interpreter','none');
%%%%%%%%;
subplot(2,2,4);
hold on;
c__ = colormap_beach(); n_c = size(c__,1);
nc_ = max(0,min(n_c-1,floor(n_c*image_sigma_delta_0est_/sigma_delta)));
patch(2*pi*[+0;+1;+1;+0],1*pi*[+0;+0;+1;+1],0.65*[1,1,1],'EdgeColor','none');
scatter(0*pi+azimu_b_all_,1*pi-polar_a_all_,markersize_use,c__(1+nc_,:),'filled');
scatter(0*pi+azimu_b_all_(1+nM0),1*pi-polar_a_all_(1+nM0),markersize_use,[0,0,0]);
hold off;
axisnotick;
axis([0,2*pi,0,1*pi]); axis equal; xlabel('azimu_b'); ylabel('polar_a'); title('sigma_delta_0est','Interpreter','none');
%%%%%%%%;
figbig;
end;%

disp('returning'); return;
end;%if nargin<1;

if (nargin<5); image_flag_update_0in_ = []; end;
if (nargin<6); image_gamma_z_0in_ = []; end;
if (nargin<7); image_delta_0in_dM__ = []; end;
if (nargin<8); n_iteration = []; end;
if (nargin<9); eta_step = []; end;
if (isempty(image_flag_update_0in_)); image_flag_update_0in_ = zeros(n_M,1); nM0 = max(0,min(n_M-1,floor(n_M*rand()))); image_flag_update_0in_(1+nM0) = 0; end;
if (isempty(image_gamma_z_0in_)); image_gamma_z_0in_ = zeros(n_M,1); end;
if (isempty(image_delta_0in_dM__)); image_delta_0in_dM__ = zeros(2,n_M); end;
if (isempty(n_iteration)); n_iteration = 64; end;
if (isempty(eta_step)); eta_step = 1.0; end;

verbose=0;
if (verbose); disp(sprintf(' %% [entering mp_init_align_gamma_z_delta_0]')); end;

n_neighbor = 6;
k_knn_index__ = knnsearch(k_c_Mk__,k_c_Mk__,'K',1+n_neighbor) - 1;
image_gamma_z_out_ = image_gamma_z_0in_;
image_delta_out_dM__ = image_delta_0in_dM__;
image_sigma_gamma_z_out_ = zeros(n_M,1);
image_sigma_delta_out_ = zeros(n_M,1);
image_flag_update_out_ = image_flag_update_0in_;
cos_d_gamma_z_MM__ = cos(d_gamma_z_MM__);
sin_d_gamma_z_MM__ = sin(d_gamma_z_MM__);
%%%%%%%%;
% Given two nearby images, as well as the dg and dd_ between them, ;
% we can determine d1_ and g1_ from d0_ and g0_ as follows: ;
% g1 \approx g0 + dg ; %<-- assuming that omega_1 - omega_0 is negligible, ;
% d1_ = dd_ + R(+dg) * d0_ ; %<-- no explicit error. ;
%%%%%%%%;
fnorm_d_image_gamma_z_out_ = zeros(n_iteration,1);
fnorm_d_image_delta_out_ = zeros(n_iteration,1);
niteration=0;
while (niteration<n_iteration);
image_gamma_z_tmp_ = image_gamma_z_out_;
image_delta_tmp_dM__ = image_delta_out_dM__;
image_flag_update_tmp_ = image_flag_update_out_;
for nM0=0:n_M-1;
k_knn_index_ = k_knn_index__(1+nM0,:);
nM1_ = k_knn_index_(1+1+[0:n_neighbor-1]);
%%%%%%%%;
if (image_flag_update_out_(1+nM0)==0);
if (sum(image_flag_update_out_(1+nM1_))==0);
% do nothing. ;
else;
tmp_index_ = efind(image_flag_update_out_(1+nM1_)> 0);
nM1_sub_ = nM1_(1+tmp_index_); n_neighbor_sub = numel(nM1_sub_);
dg_ = reshape(d_gamma_z_MM__(1+nM1_sub_,1+nM0),[n_neighbor_sub,1]); %<-- note transposition of index 0 and 1. ;
cdg_ = transpose(cos_d_gamma_z_MM__(1+nM1_sub_,1+nM0)); sdg_ = transpose(sin_d_gamma_z_MM__(1+nM1_sub_,1+nM0));
dd__ = reshape(d_delta_dMM___(:,1+nM1_sub_,1+nM0),[2,n_neighbor_sub]); %<-- note transposition of index 0 and 1. ;
g1 = image_gamma_z_out_(1+nM1_sub_);
g0 = circle_mean(g1 + dg_) ;
sigma_g0 = sqrt(mean(periodize((g1 + dg_) - g0,-pi,+pi).^2));
d1x_ = image_delta_out_dM__(1+0,1+nM1_sub_);
d1y_ = image_delta_out_dM__(1+1,1+nM1_sub_);
d0_ = mean( dd__ +  [ +cdg_.*d1x_ - sdg_.*d1y_ ; +sdg_.*d1x_ + cdg_.*d1y_ ] , 2 );
sigma_d0 = sqrt(mean((sum((( dd__ +  [ +cdg_.*d1x_ - sdg_.*d1y_ ; +sdg_.*d1x_ + cdg_.*d1y_ ] ) - repmat(d0_,[1,n_neighbor_sub])).^2,1)),2));
image_gamma_z_tmp_(1+nM0) = periodize(g0,-pi,+pi);
image_delta_tmp_dM__(:,1+nM0) = d0_;
image_sigma_gamma_z_out_(1+nM0) = sigma_g0;
image_sigma_delta_out_(1+nM0) = sigma_d0;
image_flag_update_tmp_(1+nM0) = 1;
end;%if (sum(image_flag_update_out_(1+nM1_))==0);
end;%if (image_flag_update_out_(1+nM0)==0);
%%%%%%%%;
if (image_flag_update_out_(1+nM0)> 0);
if (sum(image_flag_update_out_(1+nM1_))<=image_flag_update_out_(1+nM0));
% do nothing. ;
else;
tmp_index_ = efind(image_flag_update_out_(1+nM1_)> 0);
nM1_sub_ = setdiff(nM1_(1+tmp_index_),nM0); n_neighbor_sub = numel(nM1_sub_);
dg_ = reshape(d_gamma_z_MM__(1+nM1_sub_,1+nM0),[n_neighbor_sub,1]); %<-- note transposition of index 0 and 1. ;
cdg_ = transpose(cos_d_gamma_z_MM__(1+nM1_sub_,1+nM0)); sdg_ = transpose(sin_d_gamma_z_MM__(1+nM1_sub_,1+nM0));
dd__ = reshape(d_delta_dMM___(:,1+nM1_sub_,1+nM0),[2,n_neighbor_sub]); %<-- note transposition of index 0 and 1. ;
g1 = image_gamma_z_out_(1+nM1_sub_);
g0 = circle_mean(g1 + dg_) ;
sigma_g0 = sqrt(mean(periodize((g1 + dg_) - g0,-pi,+pi).^2));
d1x_ = image_delta_out_dM__(1+0,1+nM1_sub_);
d1y_ = image_delta_out_dM__(1+1,1+nM1_sub_);
d0_ = mean( dd__ +  [ +cdg_.*d1x_ - sdg_.*d1y_ ; +sdg_.*d1x_ + cdg_.*d1y_ ] , 2 ) ;
sigma_d0 = sqrt(mean((sum((( dd__ +  [ +cdg_.*d1x_ - sdg_.*d1y_ ; +sdg_.*d1x_ + cdg_.*d1y_ ] ) - repmat(d0_,[1,n_neighbor_sub])).^2,1)),2));
image_gamma_z_tmp_(1+nM0) = image_gamma_z_tmp_(1+nM0) + eta_step * periodize(g0 - image_gamma_z_tmp_(1+nM0),-pi/2,+pi/2) ; %<-- eta-step. ;
image_gamma_z_tmp_(1+nM0) = periodize(image_gamma_z_tmp_(1+nM0),-pi,+pi);
image_sigma_gamma_z_out_(1+nM0) = sigma_g0;
image_sigma_delta_out_(1+nM0) = sigma_d0;
image_delta_tmp_dM__(:,1+nM0) = image_delta_tmp_dM__(:,1+nM0) + eta_step * (d0_ - image_delta_tmp_dM__(:,1+nM0)) ; %<-- eta-step. ;
end;%if (sum(image_flag_update_out_(1+nM1_))<=image_flag_update_out_(1+nM0));
end;%if (image_flag_update_out_(1+nM0)> 0);
%%%%%%%%;
end;%for nM0=0:n_M-1;
d_gamma_z = periodize( image_gamma_z_out_ - image_gamma_z_tmp_ , -pi , + pi );
fnorm_d_image_gamma_z_out_(1+niteration) = fnorm(d_gamma_z);
fnorm_d_image_delta_out_(1+niteration) = fnorm(image_delta_out_dM__ - image_delta_tmp_dM__);
image_gamma_z_out_ = image_gamma_z_tmp_;
image_delta_out_dM__ = image_delta_tmp_dM__;
image_flag_update_out_ = image_flag_update_tmp_;
%%%%%%%%;
if ( (fnorm_d_image_gamma_z_out_(1+niteration)<1e-3) & (fnorm_d_image_delta_out_(1+niteration)<1e-3) );
if (verbose); disp(sprintf(' %% niteration %d; terminating. ',niteration)); end;
n_iteration = 1+niteration;
fnorm_d_image_gamma_z_out_ = fnorm_d_image_gamma_z_out_(1:n_iteration);
fnorm_d_image_delta_out_ = fnorm_d_image_delta_out_(1:n_iteration);
end;%if ( (fnorm_d_image_gamma_z_out_(1+niteration)<1e-3) & (fnorm_d_image_delta_out_(1+niteration)<1e-3) );
%%%%%%%%;
niteration=niteration+1;
end;%while (niteration<n_iteration);
if (verbose); disp( [ fnorm_d_image_gamma_z_out_ , fnorm_d_image_delta_out_ ] ); end;

if (verbose); disp(sprintf(' %% [finished mp_init_align_gamma_z_delta_0]')); end;
