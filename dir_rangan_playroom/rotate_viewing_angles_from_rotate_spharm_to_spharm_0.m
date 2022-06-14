function ...
[ ...
 R_ ...
,viewing_euler_polar_a_pos_ ...
,viewing_euler_azimu_b_pos_ ...
,viewing_euler_gamma_z_pos_ ...
] = ...
rotate_viewing_angles_from_rotate_spharm_to_spharm_0( ...
 flag_verbose ...
,rotate_euler_polar_a ...
,rotate_euler_azimu_b ...
,rotate_euler_gamma_z ...
,viewing_euler_polar_a_pre_ ...
,viewing_euler_azimu_b_pre_ ...
,viewing_euler_gamma_z_pre_ ...
);

if (nargin<1);
rng(0);
n_M  = 16;
viewing_euler_polar_a_pre_ = 1*pi*rand(n_M,1);
viewing_euler_azimu_b_pre_ = 2*pi*rand(n_M,1);
viewing_euler_gamma_z_pre_ = 2*pi*rand(n_M,1);
rotate_euler_polar_a = +pi/3;
rotate_euler_azimu_b = -pi/4;
rotate_euler_gamma_z = +6*pi/7;
flag_verbose=0;
[ ...
 R_ ...
,viewing_euler_polar_a_pos_ ...
,viewing_euler_azimu_b_pos_ ...
,viewing_euler_gamma_z_pos_ ...
] = ...
rotate_viewing_angles_from_rotate_spharm_to_spharm_0( ...
 flag_verbose ...
,rotate_euler_polar_a ...
,rotate_euler_azimu_b ...
,rotate_euler_gamma_z ...
,viewing_euler_polar_a_pre_ ...
,viewing_euler_azimu_b_pre_ ...
,viewing_euler_gamma_z_pre_ ...
);
n_w = 1024;
[ ...
 k_p_polar_a_pre_wM__ ...
,k_p_azimu_b_pre_wM__ ...
,k_c_0_pre_wM__ ...
,k_c_1_pre_wM__ ...
,k_c_2_pre_wM__ ...
] = ...
cg_rhs_1( ...
 n_M ...
,n_w ...
,viewing_euler_polar_a_pre_ ...
,viewing_euler_azimu_b_pre_ ...
,viewing_euler_gamma_z_pre_ ...
);
[ ...
 k_p_polar_a_pos_wM__ ...
,k_p_azimu_b_pos_wM__ ...
,k_c_0_pos_wM__ ...
,k_c_1_pos_wM__ ...
,k_c_2_pos_wM__ ...
] = ...
cg_rhs_1( ...
 n_M ...
,n_w ...
,viewing_euler_polar_a_pos_ ...
,viewing_euler_azimu_b_pos_ ...
,viewing_euler_gamma_z_pos_ ...
);
%%%%%%%%;
% Now we expect: ;
% R * k_c_pos__ = k_c_pre__. ;
%%%%%%%%;
R_k_c_pos_wMd__ = transpose( R_ * transpose( [ k_c_0_pos_wM__(:) , k_c_1_pos_wM__(:) , k_c_2_pos_wM__(:) ] ) );
R_k_c_0_pos_wM__ = reshape(R_k_c_pos_wMd__(:,1+0),[n_w,n_M]);
R_k_c_1_pos_wM__ = reshape(R_k_c_pos_wMd__(:,1+1),[n_w,n_M]);
R_k_c_2_pos_wM__ = reshape(R_k_c_pos_wMd__(:,1+2),[n_w,n_M]);
tmp_error_0 = fnorm(R_k_c_0_pos_wM__ - k_c_0_pre_wM__)/fnorm(k_c_0_pre_wM__);
tmp_error_1 = fnorm(R_k_c_1_pos_wM__ - k_c_1_pre_wM__)/fnorm(k_c_1_pre_wM__);
tmp_error_2 = fnorm(R_k_c_2_pos_wM__ - k_c_2_pre_wM__)/fnorm(k_c_2_pre_wM__);
disp(sprintf(' %% tmp_error_0: %0.16f',tmp_error_0));
disp(sprintf(' %% tmp_error_1: %0.16f',tmp_error_1));
disp(sprintf(' %% tmp_error_2: %0.16f',tmp_error_2));
%%%%%%%%;
figure(1);clf;figmed;
c_80s__ = colormap_80s; n_c_80s = size(c_80s__,1);
linewidth_use = 3;
markersize_use = 8;
%%%%%%%%;
subplot(1,2,1);
hold on;
for nM=0:n_M-1;
nc_80s = max(0,min(n_c_80s-1,floor(n_c_80s*nM/n_M)));
if (n_M> 2);
plot3(k_c_0_pre_wM__(:,1+nM),k_c_1_pre_wM__(:,1+nM),k_c_2_pre_wM__(:,1+nM),'-','Color',c_80s__(1+nc_80s,:),'LineWidth',linewidth_use);
end;%if (n_M> 2);
if (n_M<=2);
for nw=0:n_w-1;
nc_80s = max(0,min(n_c_80s-1,floor(n_c_80s*nw/n_w)));
plot3(k_c_0_pre_wM__(1+nw,1+nM),k_c_1_pre_wM__(1+nw,1+nM),k_c_2_pre_wM__(1+nw,1+nM),'o','Color',c_80s__(1+nc_80s,:),'MarkerSize',markersize_use);
end;%for nw=0:n_w-1;
end;%if (n_M<=2);
end;%for nM=0:n_M-1;
hold off;
xlabel('x'); ylabel('y'); zlabel('z');
axis vis3d;
title('pre');
%%%%%%%%;
subplot(1,2,2);
hold on;
for nM=0:n_M-1;
nc_80s = max(0,min(n_c_80s-1,floor(n_c_80s*nM/n_M)));
if (n_M> 2);
plot3(R_k_c_0_pos_wM__(:,1+nM),R_k_c_1_pos_wM__(:,1+nM),R_k_c_2_pos_wM__(:,1+nM),'-','Color',c_80s__(1+nc_80s,:),'LineWidth',linewidth_use);
end;%if (n_M> 2);
if (n_M<=2);
for nw=0:n_w-1;
nc_80s = max(0,min(n_c_80s-1,floor(n_c_80s*nw/n_w)));
plot3(R_k_c_0_pos_wM__(1+nw,1+nM),R_k_c_1_pos_wM__(1+nw,1+nM),R_k_c_2_pos_wM__(1+nw,1+nM),'o','Color',c_80s__(1+nc_80s,:),'MarkerSize',markersize_use);
end;%for nw=0:n_w-1;
end;%if (n_M<=2);
end;%for nM=0:n_M-1;
hold off;
xlabel('x'); ylabel('y'); zlabel('z');
axis vis3d;
title('R*pos');
%%%%%%%%;
disp('returning');return;
end;%if (nargin<1);

str_thisfunction = 'rotate_viewing_angles_from_rotate_spharm_to_spharm_0';

if (flag_verbose); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

n_M = numel(viewing_euler_polar_a_pre_);
assert(n_M==numel(viewing_euler_azimu_b_pre_));
assert(n_M==numel(viewing_euler_gamma_z_pre_));
viewing_euler_polar_a_pre_ = reshape(viewing_euler_polar_a_pre_,[n_M,1]);
viewing_euler_azimu_b_pre_ = reshape(viewing_euler_azimu_b_pre_,[n_M,1]);
viewing_euler_gamma_z_pre_ = reshape(viewing_euler_gamma_z_pre_,[n_M,1]);

%%%%%%%%;
% Here we assume that the rotation given by: ;
% rotate_euler_ = [+rotate_euler_azimu_b,+rotate_euler_polar_a,+rotate_euler_gamma_z];
% is applied to the latter of the two volumes(say b_k_Y_) passed into ;
% the function register_spharm_to_spharm_wigner_wrap_1. ;
%%%%%%%%;
% For this particular function, ;
% we align the latter volume b_k_Y_ to the former by calling: ;
% rotate_spharm_to_spharm_2 using rotate_euler_inv_ given by: ;
% rotate_euler_inv_ = [-rotate_euler_gamma_z,-rotate_euler_polar_a,-rotate_euler_azimu_b];
% (see spharm_register_and_rotate_2). ;
%%%%%%%%;
% Consequently, we expect the new volume c_k_Y_ (i.e., the rotated version of b_k_Y_) ;
% to be given by: ;
% c_k_Y_(k_) = Rz[-rotate_euler_azimu_b]*Ry[-rotate_euler_polar_a]*Rz[-rotate_euler_gamma_z]*b_k_Y_(k_), ;
% or: ;
% c_k_Y_(k_) = b_k_Y_( Rz[+rotate_euler_gamma_z]*Ry[+rotate_euler_polar_a]*Rz[+rotate_euler_azimu_b]*k_ ). ;
%            = b_k_Y_( R_*k_ ). ;
%%%%%%%%;
rotate_euler_ = [+rotate_euler_azimu_b,+rotate_euler_polar_a,+rotate_euler_gamma_z];
R_ = euler_to_R(rotate_euler_);
R_inv_ = inv(R_);
%%%%%%%%;
% Now for a particular template associated with the viewing_euler_pre_: ;
% viewing_euler_pre_ = [+viewing_euler_azimu_b_pre,+viewing_euler_polar_a_pre,+viewing_euler_gamma_z_pre];
% we expect the value of psi_(1+nw) = 2*pi*nw/n_w to correspond to the point: ;
% k_pre_(1+nw) = Rz[+viewing_euler_azimu_b_pre]*Ry[+viewing_euler_polar_a_pre]*Rz[psi_(1+nw) - viewing_euler_gamma_z_pre]*[1;0;0]. ;
% Consequently, we expect the new viewing_euler_pos to have the property: ;
% c_k_Y_( Rz[+viewing_euler_azimu_b_pos]*Ry[+viewing_euler_polar_a_pos]*Rz[psi_(1+nw) - viewing_euler_gamma_z_pos]*[1;0;0] ) ;
% = ;
% b_k_Y_( Rz[+viewing_euler_azimu_b_pre]*Ry[+viewing_euler_polar_a_pre]*Rz[psi_(1+nw) - viewing_euler_gamma_z_pre]*[1;0;0] ) ;
% = ;
% b_k_Y_( R_ * Rz[+viewing_euler_azimu_b_pos]*Ry[+viewing_euler_polar_a_pos]*Rz[psi_(1+nw) - viewing_euler_gamma_z_pos]*[1;0;0] ) ;
% implying that: ;
% Rz[+viewing_euler_azimu_b_pos]*Ry[+viewing_euler_polar_a_pos]*Rz[psi_(1+nw) - viewing_euler_gamma_z_pos]*[1;0;0] ;
% = ;
% R_inv_ * Rz[+viewing_euler_azimu_b_pre]*Ry[+viewing_euler_polar_a_pre]*Rz[psi_(1+nw) - viewing_euler_gamma_z_pre]*[1;0;0]. ;
%%%%%%%%;
% With this in mind, we can define: ;
% R0_ = R_inv_ * Rz[+viewing_euler_azimu_b_pre]*Ry[+viewing_euler_polar_a_pre], ;
% such that R0_ rotates the pole associated with the viewing_euler_pre_ (i.e., at [0;0;1]) ;
% to a new location, which must be the pole for viewing_euler_pos_. ;
% This new pole defines viewing_euler_polar_a_pos_ and viewing_euler_azimu_b_pos_. ;
% This process can be simplified by defining: ;
% Rz[+viewing_euler_azimu_b_pre]*Ry[+viewing_euler_polar_a_pre]*[0;0;1] = [ +sa_pre * cb_pre ; +sa_pre * sb_pre ; +ca_pre ];
%%%%%%%%;
sa_pre_ = sin(viewing_euler_polar_a_pre_); ca_pre_ = cos(viewing_euler_polar_a_pre_);
sb_pre_ = sin(viewing_euler_azimu_b_pre_); cb_pre_ = cos(viewing_euler_azimu_b_pre_);
k0_pre_ = sa_pre_.*cb_pre_;
k1_pre_ = sa_pre_.*sb_pre_;
k2_pre_ = ca_pre_;
k_pos__ = transpose( R_inv_ * transpose( [k0_pre_ , k1_pre_ , k2_pre_ ] ) );
k0_pos_ = k_pos__(:,1+0);
k1_pos_ = k_pos__(:,1+1);
k2_pos_ = k_pos__(:,1+2);
k01_pos_ = sqrt(k0_pos_.^2 + k1_pos_.^2);
viewing_euler_azimu_b_pos_ = periodize(atan2(k1_pos_,k0_pos_),0,2*pi);
viewing_euler_polar_a_pos_ = atan2(k01_pos_,k2_pos_);
%%%%%%%%;
% Once we have defined viewing_euler_polar_a_pos and viewing_euler_azimu_b_pos, we can define: ;
% R1_ = Ry[-viewing_euler_polar_a_pos]*Rz[-viewing_euler_azimu_b_pos]*R0_, ;
% such that R1_ rotates the xy-plane by some angle nu. ;
% If R1_ corresponds to Rz[nu] = [ +cos(nu) , -sin(nu) , 0 ; +sin(nu) , +cos(nu) , 0 ; 0 , 0 , 1 ], then: ;
% -viewing_euler_gamma_z_pos = nu - viewing_euler_gamma_z_pre, ;
% implying that viewing_euler_gamma_z_pos = viewing_euler_gamma_z_pre - nu. ;
%%%%%%%%;
% This calculation can be further simplified by noting that: ;
% R0_ * [1;0;0] = R_inv_ * Rz[+viewing_euler_azimu_b_pre] * Ry[+viewing_euler_polar_a_pre] * [1;0;0] ;
%               = R_inv_ * Rz[+viewing_euler_azimu_b_pre] * Ry[+viewing_euler_polar_a_pre] * [1;0;0] ;
% So R1_ * [1;0;0] = ;
% [ +ca_pos , +0 , -sa_pos ]   [ +cb_pos , +sb_pos , +0 ]            [ +ca_pre*cb_pre ] ;
% [ +0      , +1 , +0      ] * [ -sb_pos , +cb_pos , +0 ] * R_inv_ * [ +ca_pre*sb_pre ] ;
% [ +sa_pos , +0 , +ca_pos ]   [ +0      , +0      , +1 ]            [ -sa_pre        ] ;
% = ;
% [ +cos(nu) ] ;
% [ +sin(nu) ] ;
% [ ~        ] ;
%%%%%%%%;
sa_pos_ = sin(viewing_euler_polar_a_pos_); ca_pos_ = cos(viewing_euler_polar_a_pos_);
sb_pos_ = sin(viewing_euler_azimu_b_pos_); cb_pos_ = cos(viewing_euler_azimu_b_pos_);
viewing_euler_gamma_z_pos_ = zeros(n_M,1);
for nM=0:n_M-1;
viewing_euler_gamma_z_pre = viewing_euler_gamma_z_pre_(1+nM);
sa_pre = sa_pre_(1+nM); ca_pre = ca_pre_(1+nM);
sb_pre = sb_pre_(1+nM); cb_pre = cb_pre_(1+nM);
sa_pos = sa_pos_(1+nM); ca_pos = ca_pos_(1+nM);
sb_pos = sb_pos_(1+nM); cb_pos = cb_pos_(1+nM);
Ry_ = [ +ca_pos , +0      , -sa_pos ; +0      , +1      , +0 ; +sa_pos , +0 , +ca_pos ];
Rz_ = [ +cb_pos , +sb_pos , +0      ; -sb_pos , +cb_pos , +0 ; +0      , +0 , +1      ];
tmp_ = Ry_ * Rz_ * R_inv_ * [ +ca_pre*cb_pre ; +ca_pre*sb_pre ; -sa_pre ];
nu = atan2(tmp_(1+1),tmp_(1+0));
viewing_euler_gamma_z_pos = viewing_euler_gamma_z_pre - nu;
viewing_euler_gamma_z_pos_(1+nM) = viewing_euler_gamma_z_pos;
end;%for nM=0:n_M-1;

if (flag_verbose); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
