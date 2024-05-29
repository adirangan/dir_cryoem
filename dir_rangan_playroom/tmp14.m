%%%%%%%%;
% look at the template-manifold for rib80s. ;
% start with test_pm_23. ;
%%%%%%%%;
flag_disp=1;
dx_u_pack = diameter_x_c/max(1,n_x_u_pack);
a_k_Y_quad_norm_ = spharm_normalize_2(n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_Y_quad_);
a_k_Y_quad_norm__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
a_k_Y_quad_norm__(1:n_lm_(1+nk_p_r),1+nk_p_r) = a_k_Y_quad_norm_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;


[ ...
 n_net ...
,azimu_b_net_ ...
,polar_a_net_ ...
,weight_net_ ...
,k_c_0_net_ ...
,k_c_1_net_ ...
,k_c_2_net_ ...
,n_polar_a_net_per ...
,polar_a_net_per_ ...
,n_azimu_b_net_per_ ...
] = ...
sample_shell_6( ...
 1.0d0 ...
,1.0/(2*pi) ...
,'L' ...
,1 ...
) ;
n_azimu_b_net_per = n_azimu_b_net_per_(1+0);
azimu_b_net_per_ = (2*pi)*transpose(0:n_azimu_b_net_per-1)/max(1,n_azimu_b_net_per);
k_c_0_net_ba__ = reshape(k_c_0_net_,[n_azimu_b_net_per,n_polar_a_net_per]); k_c_0_net_ba__ = [k_c_0_net_ba__;k_c_0_net_ba__(1,:)];
k_c_1_net_ba__ = reshape(k_c_1_net_,[n_azimu_b_net_per,n_polar_a_net_per]); k_c_1_net_ba__ = [k_c_1_net_ba__;k_c_1_net_ba__(1,:)];
k_c_2_net_ba__ = reshape(k_c_2_net_,[n_azimu_b_net_per,n_polar_a_net_per]); k_c_2_net_ba__ = [k_c_2_net_ba__;k_c_2_net_ba__(1,:)];
k_c_0_net_ab__ = transpose(k_c_0_net_ba__);
k_c_1_net_ab__ = transpose(k_c_1_net_ba__);
k_c_2_net_ab__ = transpose(k_c_2_net_ba__);
%%%%%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;
hold on;
line(k_c_0_net_ba__,k_c_1_net_ba__,k_c_2_net_ba__,'Color',0.85*[1,1,1]);
line(k_c_0_net_ab__,k_c_1_net_ab__,k_c_2_net_ab__,'Color',0.85*[1,1,1]);
hold off;
end;%if flag_disp;
%%%%%%%%;
n_R = n_net;
R_k_p__ = ...
pm_template_2( ...
 verbose ...
,l_max_max ...
,n_k_p_r ...
,a_k_Y_quad_norm__ ...
,[] ...
,-1 ...
,n_w_max ...
,n_net ...
,azimu_b_net_ ...
,polar_a_net_ ...
,weight_net_ ...
,n_polar_a_net_per ...
,polar_a_net_per_ ...
,n_azimu_b_net_per_ ...
);
R_k_p__ = reshape(R_k_p__,[n_w_max*n_k_p_r,n_net]);
%%%%%%%%;

%%%%%%%%;
n_gamma_z_net_per = 3;
gamma_z_net_per_ = (2*pi)*transpose(0:n_gamma_z_net_per-1)/max(1,n_gamma_z_net_per);
R_x_c_xxzR____ = zeros(n_x_u_pack,n_x_u_pack,n_gamma_z_net_per,n_R);
for nR=0:n_R-1;
R_k_p_ = R_k_p__(:,1+nR);
for ngamma_z_net_per=0:n_gamma_z_net_per-1;
gamma_z_net_per = gamma_z_net_per_(1+ngamma_z_net_per);
R_x_c__ = ...
interp_k_p_to_x_c_xxnufft( ...
 n_x_u_pack ...
,diameter_x_c ...
,n_x_u_pack ...
,diameter_x_c ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,R_k_p_,gamma_z_net_per).*CTF_avg_k_p_wk_ ...
);
R_x_c_xxzR____(:,:,1+ngamma_z_net_per,1+nR) = real(R_x_c__);
end;%for ngamma_z_net_per=0:n_gamma_z_net_per-1;
end;%for nR=0:n_R-1;
%%%%%%%%;

scale_R = 2.0e4;
%%%%;
[U_R__,S_R__,V_R__] = svds(reshape(scale_R*R_x_c_xxzR____,[n_x_u_pack.^2,n_gamma_z_net_per*n_R]),3); S_R_ = diag(S_R__);
SV_R__ = transpose(transpose(U_R__)*reshape(scale_R*R_x_c_xxzR____,[n_x_u_pack.^2,n_gamma_z_net_per*n_R]));
SV_R_0_net_ = SV_R__(:,1+0);
SV_R_1_net_ = SV_R__(:,1+1);
SV_R_2_net_ = SV_R__(:,1+2);
SV_R_0_net_baz___ = permute(reshape(SV_R_0_net_,[n_gamma_z_net_per,n_azimu_b_net_per,n_polar_a_net_per]),[2,3,1]);
SV_R_1_net_baz___ = permute(reshape(SV_R_1_net_,[n_gamma_z_net_per,n_azimu_b_net_per,n_polar_a_net_per]),[2,3,1]);
SV_R_2_net_baz___ = permute(reshape(SV_R_2_net_,[n_gamma_z_net_per,n_azimu_b_net_per,n_polar_a_net_per]),[2,3,1]);
SV_R_0_net_baz___ = cat(1,SV_R_0_net_baz___,SV_R_0_net_baz___(1,:,:));
SV_R_1_net_baz___ = cat(1,SV_R_1_net_baz___,SV_R_1_net_baz___(1,:,:));
SV_R_2_net_baz___ = cat(1,SV_R_2_net_baz___,SV_R_2_net_baz___(1,:,:));
SV_R_0_net_abz___ = permute(SV_R_0_net_baz___,[2,1,3]);
SV_R_1_net_abz___ = permute(SV_R_1_net_baz___,[2,1,3]);
SV_R_2_net_abz___ = permute(SV_R_2_net_baz___,[2,1,3]);
%%%%%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;
c_0_ = 0.85*[1.00,0.65,0.65];
c_1_ = 0.85*[0.65,1.00,0.65];
c_2_ = 0.85*[0.65,0.65,1.00];
hold on;
line(SV_R_0_net_baz___(:,:,1+0),SV_R_1_net_baz___(:,:,1+0),SV_R_2_net_baz___(:,:,1+0),'Color',c_0_);
line(SV_R_0_net_abz___(:,:,1+0),SV_R_1_net_abz___(:,:,1+0),SV_R_2_net_abz___(:,:,1+0),'Color',c_0_);
line(SV_R_0_net_baz___(:,:,1+1),SV_R_1_net_baz___(:,:,1+1),SV_R_2_net_baz___(:,:,1+1),'Color',c_1_);
line(SV_R_0_net_abz___(:,:,1+1),SV_R_1_net_abz___(:,:,1+1),SV_R_2_net_abz___(:,:,1+1),'Color',c_1_);
line(SV_R_0_net_baz___(:,:,1+2),SV_R_1_net_baz___(:,:,1+2),SV_R_2_net_baz___(:,:,1+2),'Color',c_2_);
line(SV_R_0_net_abz___(:,:,1+2),SV_R_1_net_abz___(:,:,1+2),SV_R_2_net_abz___(:,:,1+2),'Color',c_2_);
hold off;
axis equal; axis vis3d;
end;%if flag_disp;

%%%%;
M_x_c___ = zeros(n_x_u_pack,n_x_u_pack,n_M);
for nM=0:n_M-1;
M_k_p_ = M_k_p__(:,1+nM);
M_x_c__ = ...
interp_k_p_to_x_c_xxnufft( ...
 n_x_u_pack ...
,diameter_x_c ...
,n_x_u_pack ...
,diameter_x_c ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,M_k_p_ ...
);
M_x_c___(:,:,1+nM) = real(M_x_c__)/max(1e-12,sqrt(sum(real(M_x_c__).^2,'all')*dx_u_pack^2));
end;%for nM=0:n_M-1;
%%%%;
V_M__ = transpose(inv(S_R__)*transpose(U_R__)*reshape(M_x_c___,[n_x_u_pack.^2,n_M]));
SV_M__ = transpose(transpose(U_R__)*reshape(M_x_c___,[n_x_u_pack.^2,n_M]));
%%%%;
%[tmp_ij_,tmp_d2_] = knnsearch(SV_R__,SV_M__,'K',1); [~,tmp_ij_min] = min(tmp_d2_); 
%nN_from_nM = tmp_ij_min-1;
%%%%%%%%;

%%%%%%%%;
nN_from_nM = 3;
n_shift = 8;
n_N = (1+2*n_shift)^2;
N_x_c__ = M_x_c___(:,:,1+nN_from_nM);
N_x_c_xxN___ = zeros(n_x_u_pack,n_x_u_pack,n_N);
nN=0;
for nshift0=-n_shift:+n_shift;
for nshift1=-n_shift:+n_shift;
N_x_c_xxN___(:,:,1+nN) = circshift(N_x_c__,[nshift0,nshift1]);
nN = nN+1;
end;%for nshift1=-n_shift:+n_shift;
end;%for nshift0=-n_shift:+n_shift;
%%%%;
V_N__ = transpose(inv(S_R__)*transpose(U_R__)*reshape(N_x_c_xxN___,[n_x_u_pack.^2,n_N]));
nN_mid = (n_shift) + (n_shift)*(1+2*n_shift);
SV_N__ = transpose(transpose(U_R__)*reshape(N_x_c_xxN___,[n_x_u_pack.^2,n_N]));
SV_N_0_ = SV_N__(:,1+0);
SV_N_1_ = SV_N__(:,1+1);
SV_N_2_ = SV_N__(:,1+2);
SV_N_0_01__ = reshape(SV_N_0_,[1+2*n_shift,1+2*n_shift]);
SV_N_1_01__ = reshape(SV_N_1_,[1+2*n_shift,1+2*n_shift]);
SV_N_2_01__ = reshape(SV_N_2_,[1+2*n_shift,1+2*n_shift]);
SV_N_0_10__ = transpose(SV_N_0_01__);
SV_N_1_10__ = transpose(SV_N_1_01__);
SV_N_2_10__ = transpose(SV_N_2_01__);
%%%%%%%%;

%%%%%%%%;
if flag_disp;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
markersize_use = 12;
linewidth_use = 2;
subplot(1,1,1);
c_0_ = 0.85*[1.00,0.65,0.65];
c_1_ = 0.85*[0.65,1.00,0.65];
c_2_ = 0.85*[0.65,0.65,1.00];
c_N_ = 0.85*[1,1,1];
hold on;
line(SV_R_0_net_baz___(:,:,1+0),SV_R_1_net_baz___(:,:,1+0),SV_R_2_net_baz___(:,:,1+0),'Color',c_0_);
line(SV_R_0_net_abz___(:,:,1+0),SV_R_1_net_abz___(:,:,1+0),SV_R_2_net_abz___(:,:,1+0),'Color',c_0_);
line(SV_R_0_net_baz___(:,:,1+1),SV_R_1_net_baz___(:,:,1+1),SV_R_2_net_baz___(:,:,1+1),'Color',c_1_);
line(SV_R_0_net_abz___(:,:,1+1),SV_R_1_net_abz___(:,:,1+1),SV_R_2_net_abz___(:,:,1+1),'Color',c_1_);
line(SV_R_0_net_baz___(:,:,1+2),SV_R_1_net_baz___(:,:,1+2),SV_R_2_net_baz___(:,:,1+2),'Color',c_2_);
line(SV_R_0_net_abz___(:,:,1+2),SV_R_1_net_abz___(:,:,1+2),SV_R_2_net_abz___(:,:,1+2),'Color',c_2_);
line(SV_N_0_10__,SV_N_1_10__,SV_N_2_10__,'Color',c_N_);
line(SV_N_0_01__,SV_N_1_01__,SV_N_2_01__,'Color',c_N_);
plot3(SV_N_0_(1+nN_mid),SV_N_1_(1+nN_mid),SV_N_2_(1+nN_mid),'o','MarkerSize',markersize_use,'MarkerFaceColor',c_N_);
hold off;
axis equal; axis vis3d;
%%%%%%%%;
end;%if flag_disp;
%%%%%%%%;

disp('returning');return;

%%%%%%%%;
