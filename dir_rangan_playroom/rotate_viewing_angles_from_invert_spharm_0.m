function ...
[ ...
 viewing_euler_gamma_z_pos_ ...
] = ...
rotate_viewing_angles_from_invert_spharm_0( ...
 flag_verbose ...
,viewing_euler_gamma_z_pre_ ...
);

if (nargin<1);
rng(0);
n_M  = 16;
viewing_euler_polar_a_pre_ = 1*pi*rand(n_M,1);
viewing_euler_azimu_b_pre_ = 2*pi*rand(n_M,1);
viewing_euler_gamma_z_pre_ = 2*pi*rand(n_M,1);
flag_verbose=0;
viewing_euler_polar_a_pos_ = viewing_euler_polar_a_pre_;
viewing_euler_azimu_b_pos_ = viewing_euler_azimu_b_pre_;
[ ...
 viewing_euler_gamma_z_pos_ ...
] = ...
rotate_viewing_angles_from_invert_spharm_0( ...
 flag_verbose ...
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
% k_c_pos__ = -k_c_pre__. ;
%%%%%%%%;
I_k_c_pos_wMd__ = transpose( (-1) * transpose( [ k_c_0_pos_wM__(:) , k_c_1_pos_wM__(:) , k_c_2_pos_wM__(:) ] ) );
I_k_c_0_pos_wM__ = reshape(I_k_c_pos_wMd__(:,1+0),[n_w,n_M]);
I_k_c_1_pos_wM__ = reshape(I_k_c_pos_wMd__(:,1+1),[n_w,n_M]);
I_k_c_2_pos_wM__ = reshape(I_k_c_pos_wMd__(:,1+2),[n_w,n_M]);
tmp_error_0 = fnorm(I_k_c_0_pos_wM__ - k_c_0_pre_wM__)/fnorm(k_c_0_pre_wM__);
tmp_error_1 = fnorm(I_k_c_1_pos_wM__ - k_c_1_pre_wM__)/fnorm(k_c_1_pre_wM__);
tmp_error_2 = fnorm(I_k_c_2_pos_wM__ - k_c_2_pre_wM__)/fnorm(k_c_2_pre_wM__);
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
plot3(I_k_c_0_pos_wM__(:,1+nM),I_k_c_1_pos_wM__(:,1+nM),I_k_c_2_pos_wM__(:,1+nM),'-','Color',c_80s__(1+nc_80s,:),'LineWidth',linewidth_use);
end;%if (n_M> 2);
if (n_M<=2);
for nw=0:n_w-1;
nc_80s = max(0,min(n_c_80s-1,floor(n_c_80s*nw/n_w)));
plot3(I_k_c_0_pos_wM__(1+nw,1+nM),I_k_c_1_pos_wM__(1+nw,1+nM),I_k_c_2_pos_wM__(1+nw,1+nM),'o','Color',c_80s__(1+nc_80s,:),'MarkerSize',markersize_use);
end;%for nw=0:n_w-1;
end;%if (n_M<=2);
end;%for nM=0:n_M-1;
hold off;
xlabel('x'); ylabel('y'); zlabel('z');
axis vis3d;
title('(-1)*pos');
%%%%%%%%;
disp('returning');return;
end;%if (nargin<1);

str_thisfunction = 'rotate_viewing_angles_from_invert_spharm_0';

if (flag_verbose); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

n_M = numel(viewing_euler_gamma_z_pre_);
viewing_euler_gamma_z_pre_ = reshape(viewing_euler_gamma_z_pre_,[n_M,1]);

%%%%%%%%;
% We assume the inversion (i.e., flipY) sends +x_c_ (a vector in 3-space) to -x_c_. ;
% Correspondingly, this sends +k_c_ (a vector in 3-space) to -k_c_. ;
% and the euler-angles polar_a and azimu_b are sent, respectively, to: ;
% (pi - polar_a) and (pi + azimu_b). ;
% For a viewing-angle the equivalent is simply: ;
% rotate_euler_gamma_z sent to rotate_euler_gamma_z + pi. ;
%%%%%%%%;
viewing_euler_gamma_z_pos_ = periodize(viewing_euler_gamma_z_pre_ + pi,0,2*pi);

if (flag_verbose); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
