verbose=1;
n_w_max = 2*(l_max_max+1);
%%%%%%%%;
% Load the images. ;
%%%%%%%%;
n_M = 1024;
[ ...
 M_x_c___ ...
,index_nCTF_from_nM_ ...
,index_nM_from_nCTF_ ...
,Voltage_CTF_ ...
,DefocusU_CTF_ ...
,DefocusV_CTF_ ...
,DefocusAngle_CTF_ ...
,SphericalAberration_CTF_ ...
,AmplitudeContrast_CTF_ ...
] = ...
rlnImageName_from_star_0( ...
 dir_data_star ...
,fname_nopath_star ...
,n_M ...
);
if (fnorm(Voltage_CTF_)< 1e-3); disp(sprintf(' %% Warning, Voltage not set, setting Voltage to 300kV')); Voltage_CTF_ = 300*ones(n_M,1); end;
%%%%%%%%;
% Remove any edge artefacts, mean center and normalize each image. ;
%%%%%%%%;
disp(sprintf(' %% Removing edge-artefacts'));
n_M_ext_ = zeros(n_M,1);
for nM=0:n_M-1;
if (mod(nM,128)==0); disp(sprintf(' %% nM %d/%d',nM,n_M)); end;
n_pixel = 4; edge_tolerance = 0.5; n_edge_overshoot = 8; rseed = 0;
[M_x_c___(:,:,1+nM),n_M_ext_(1+nM)] = image_replace_edge_artefact_0(M_x_c___(:,:,1+nM),4,0.5,8,0);
end;%for nM=0:n_M-1;
disp(sprintf(' %% edge-artefacts detected in %d/%d images.',numel(find(n_M_ext_>0)),n_M));
%%%%%%%%;
% Now examine image-centroids. ;
%%%%%%%%;
n_x_M_u = size(M_x_c___,1);
assert(n_x_M_u==size(M_x_c___,2));
x_c_0_ = -half_diameter_x_c + transpose([0:n_x_M_u-1]/n_x_M_u)*diameter_x_c;
x_c_1_ = -half_diameter_x_c + transpose([0:n_x_M_u-1]/n_x_M_u)*diameter_x_c;
[x_c_0__,x_c_1__] = ndgrid(x_c_0_,x_c_1_);
x_p_r__ = sqrt(x_c_0__.^2 + x_c_1__.^2);
x_c_mask__ = x_p_r__<=half_diameter_x_c;
M_abs_x_c_0_avg_ = zeros(n_M,1);
M_abs_x_c_1_avg_ = zeros(n_M,1);
M_mask_abs_x_c_0_avg_ = zeros(n_M,1);
M_mask_abs_x_c_1_avg_ = zeros(n_M,1);
for nM=0:n_M-1;
M_abs_x_c_ = abs(squeeze(M_x_c___(:,:,1+nM))); %<-- no mask. ;
M_abs_avg = mean(M_abs_x_c_,'all');
M_abs_x_c_0_avg = mean(M_abs_x_c_/M_abs_avg.*x_c_0__,'all');
M_abs_x_c_1_avg = mean(M_abs_x_c_/M_abs_avg.*x_c_1__,'all');
M_abs_x_c_0_avg_(1+nM) = M_abs_x_c_0_avg;
M_abs_x_c_1_avg_(1+nM) = M_abs_x_c_1_avg;
clear M_abs_x_c_;
M_mask_abs_x_c_ = abs(squeeze(M_x_c___(:,:,1+nM)).*x_c_mask__); %<-- radial mask. ;
M_mask_abs_avg = mean(M_mask_abs_x_c_,'all');
M_mask_abs_x_c_0_avg = mean(M_mask_abs_x_c_/M_mask_abs_avg.*x_c_0__,'all');
M_mask_abs_x_c_1_avg = mean(M_mask_abs_x_c_/M_mask_abs_avg.*x_c_1__,'all');
M_mask_abs_x_c_0_avg_(1+nM) = M_mask_abs_x_c_0_avg;
M_mask_abs_x_c_1_avg_(1+nM) = M_mask_abs_x_c_1_avg;
clear M_mask_abs_x_c_;
end;%for nM=0:n_M-1;
%%%%%%%%;
% Now convert images to M_k_p__. ;
%%%%%%%%;
dx = diameter_x_c/n_x_M_u;
%%%%;
M_k_p__ = zeros(n_w_sum,n_M);
M_k_q__ = zeros(n_w_sum,n_M);
M_k_p_l2_ = zeros(n_M,1);
M_k_q_l2_ = zeros(n_M,1);
for nM=0:n_M-1;
if (mod(nM,128)==0); disp(sprintf(' %% nM %d/%d',nM,n_M)); end;
M_x_c_ = squeeze(M_x_c___(:,:,1+nM));
%M_k_p_ = interp_x_c_to_k_p_nufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,M_x_c_,n_k_p_r,k_p_r_,n_w_);
%this time try: ;
M_k_p_ = interp_x_c_to_k_p_nufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,M_x_c_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_M_u^2)*dx^2 ;
% Now *DO* translate according to M_abs_x_c_0_avg_ and M_abs_x_c_1_avg_. ;
M_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p_,-1*M_abs_x_c_0_avg_(1+nM),-1*M_abs_x_c_1_avg_(1+nM));
M_k_p__(:,1+nM) = M_k_p_;
M_k_q_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_k_p_);
M_k_q__(:,1+nM) = M_k_q_;
M_k_p_l2_(1+nM) = sum(abs(M_k_p_).^2 .* weight_2d_k_all_) * (2*pi)^2;
M_k_q_l2_(1+nM) = sum(abs(M_k_q_).^2 .* weight_2d_k_all_) * (2*pi)^2;
end;%for nM=0:n_M-1;
wUX__ = diag(X_weight_r_)*UX__;
UX_M_k_p_wnM___ = zeros(n_w_max,n_UX_rank,n_M);
UX_M_k_q_wnM___ = zeros(n_w_max,n_UX_rank,n_M);
UX_M_k_p_l2_ = zeros(n_M,1);
UX_M_k_q_l2_ = zeros(n_M,1);
for nM=0:n_M-1;
tmp_M_k_p__ = reshape(M_k_p__(:,1+nM),[n_w_max,n_k_p_r]);
tmp_M_k_q__ = reshape(M_k_q__(:,1+nM),[n_w_max,n_k_p_r]);
UX_M_k_p_wnM___(:,:,1+nM) = reshape(M_k_p__(:,1+nM),[n_w_max,n_k_p_r])*wUX__;
UX_M_k_q_wnM___(:,:,1+nM) = reshape(M_k_q__(:,1+nM),[n_w_max,n_k_p_r])*wUX__;
UX_M_k_p_l2_(1+nM) = sum(abs(UX_M_k_p_wnM___(:,:,1+nM).^2),'all') * (2*pi)^2 / n_w_max / (2*pi)^2 ;
UX_M_k_q_l2_(1+nM) = sum(abs(UX_M_k_q_wnM___(:,:,1+nM).^2),'all') * (2*pi)^2 / n_w_max / (2*pi)^2 ;
end;%for nM=0:n_M-1;
%%%%%%%%;
% do the same for a set of random images. ;
%%%%%%%%;
x_p_r__ = sqrt(x_c_0__.^2 + x_c_1__.^2);
edge_M_x_c_ = x_p_r__> 1/sqrt(2);
cent_M_x_c_ = x_p_r__<=1/sqrt(2);
index_edge_M_x_c_ = efind(edge_M_x_c_);
index_cent_M_x_c_ = efind(cent_M_x_c_);
R_x_c___ = zeros(n_x_M_u,n_x_M_u,n_M);
for nM=0:n_M-1;
tmp_M_x_c__ = M_x_c___(:,:,1+nM);
rng(nM); R_x_c___(:,:,1+nM) = std(tmp_M_x_c__(1+index_edge_M_x_c_),1)*randn(n_x_M_u,n_x_M_u);
end;%for nM=0:n_M-1;
R_k_p__ = zeros(n_w_sum,n_M);
R_k_q__ = zeros(n_w_sum,n_M);
R_k_p_l2_ = zeros(n_M,1);
R_k_q_l2_ = zeros(n_M,1);
for nM=0:n_M-1;
if (mod(nM,128)==0); disp(sprintf(' %% nM %d/%d',nM,n_M)); end;
R_x_c_ = squeeze(R_x_c___(:,:,1+nM));
R_k_p_ = interp_x_c_to_k_p_nufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,R_x_c_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_M_u^2)*dx^2 ;
R_k_p__(:,1+nM) = R_k_p_;
R_k_q_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,R_k_p_);
R_k_q__(:,1+nM) = R_k_q_;
R_k_p_l2_(1+nM) = sum(abs(R_k_p_).^2 .* weight_2d_k_all_) * (2*pi)^2;
R_k_q_l2_(1+nM) = sum(abs(R_k_q_).^2 .* weight_2d_k_all_) * (2*pi)^2;
end;%for nM=0:n_M-1;
wUX__ = diag(X_weight_r_)*UX__;
UX_R_k_p_wnM___ = zeros(n_w_max,n_UX_rank,n_M);
UX_R_k_q_wnM___ = zeros(n_w_max,n_UX_rank,n_M);
UX_R_k_p_l2_ = zeros(n_M,1);
UX_R_k_q_l2_ = zeros(n_M,1);
for nM=0:n_M-1;
tmp_R_k_p__ = reshape(R_k_p__(:,1+nM),[n_w_max,n_k_p_r]);
tmp_R_k_q__ = reshape(R_k_q__(:,1+nM),[n_w_max,n_k_p_r]);
UX_R_k_p_wnM___(:,:,1+nM) = reshape(R_k_p__(:,1+nM),[n_w_max,n_k_p_r])*wUX__;
UX_R_k_q_wnM___(:,:,1+nM) = reshape(R_k_q__(:,1+nM),[n_w_max,n_k_p_r])*wUX__;
UX_R_k_p_l2_(1+nM) = sum(abs(UX_R_k_p_wnM___(:,:,1+nM).^2),'all') * (2*pi)^2 / n_w_max / (2*pi)^2 ;
UX_R_k_q_l2_(1+nM) = sum(abs(UX_R_k_q_wnM___(:,:,1+nM).^2),'all') * (2*pi)^2 / n_w_max / (2*pi)^2 ;
end;%for nM=0:n_M-1;
%%%%%%%%;
% and for the templates. ;
%%%%%%%%;
template_k_eq_d = -1;
viewing_k_eq_d = 1.0; %<-- subsample just a few templates for visualization. ;
tmp_t = tic();
[ ...
 S_k_p__ ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_k_all_ ...
,n_viewing_all ...
,viewing_azimu_b_all_ ...
,viewing_polar_a_all_ ...
,viewing_weight_all_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,template_k_c_0__ ...
,template_k_c_1__ ...
,template_k_c_2__ ...
] = ...
get_template_1( ...
 0*verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_quad_ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,n_w_max*ones(n_k_p_r,1) ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% get_template_1: %0.2fs',tmp_t)); end;
if (verbose); disp(sprintf(' %% n_viewing_all %d n_viewing_polar_a %d n_w_max %d',n_viewing_all,n_viewing_polar_a,max(n_w_))); end;
n_S = n_viewing_all; n_w_max = max(n_w_); n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_]);
%%%%;
CTF_S_k_p_l2_avg = mean(M_k_p_l2_) - mean(R_k_p_l2_);
CTF_S_k_p__ = zeros(n_w_sum,n_S);
CTF_S_k_q__ = zeros(n_w_sum,n_S);
CTF_S_k_p_l2_ = zeros(n_S,1);
CTF_S_k_q_l2_ = zeros(n_S,1);
for nS=0:n_S-1;
if (mod(nS,128)==0); disp(sprintf(' %% nS %d/%d',nS,n_S)); end;
CTF_S_k_p_ = reshape(reshape(S_k_p__(:,1+nS),[n_w_max,n_k_p_r])*diag(CTF_avg_k_p_r_),[n_w_sum,1]);
tmp_CTF_S_k_p_l2 = sum(abs(CTF_S_k_p_).^2 .* weight_2d_k_all_) * (2*pi)^2;
CTF_S_k_p_ = CTF_S_k_p_/sqrt(tmp_CTF_S_k_p_l2)*sqrt(CTF_S_k_p_l2_avg);
CTF_S_k_p__(:,1+nS) = CTF_S_k_p_;
CTF_S_k_q_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,CTF_S_k_p_);
CTF_S_k_q__(:,1+nS) = CTF_S_k_q_;
CTF_S_k_p_l2_(1+nS) = sum(abs(CTF_S_k_p_).^2 .* weight_2d_k_all_) * (2*pi)^2;
CTF_S_k_q_l2_(1+nS) = sum(abs(CTF_S_k_q_).^2 .* weight_2d_k_all_) * (2*pi)^2;
end;%for nS=0:n_S-1;
wUX__ = diag(X_weight_r_)*UX__;
UX_CTF_S_k_p_wnS___ = zeros(n_w_max,n_UX_rank,n_S);
UX_CTF_S_k_q_wnS___ = zeros(n_w_max,n_UX_rank,n_S);
UX_CTF_S_k_p_l2_ = zeros(n_S,1);
UX_CTF_S_k_q_l2_ = zeros(n_S,1);
for nS=0:n_S-1;
tmp_CTF_S_k_p__ = reshape(CTF_S_k_p__(:,1+nS),[n_w_max,n_k_p_r]);
tmp_CTF_S_k_q__ = reshape(CTF_S_k_q__(:,1+nS),[n_w_max,n_k_p_r]);
UX_CTF_S_k_p_wnS___(:,:,1+nS) = reshape(CTF_S_k_p__(:,1+nS),[n_w_max,n_k_p_r])*wUX__;
UX_CTF_S_k_q_wnS___(:,:,1+nS) = reshape(CTF_S_k_q__(:,1+nS),[n_w_max,n_k_p_r])*wUX__;
UX_CTF_S_k_p_l2_(1+nS) = sum(abs(UX_CTF_S_k_p_wnS___(:,:,1+nS).^2),'all') * (2*pi)^2 / n_w_max / (2*pi)^2 ;
UX_CTF_S_k_q_l2_(1+nS) = sum(abs(UX_CTF_S_k_q_wnS___(:,:,1+nS).^2),'all') * (2*pi)^2 / n_w_max / (2*pi)^2 ;
end;%for nS=0:n_S-1;
%%%%%%%%;

flag_check = 1;
if flag_check;
%%%%;
% assess random image. ;
%%%%;
tmp_R_x_c_ = R_x_c___(:,:,1+nM);
tmp_R_k_p_ = interp_x_c_to_k_p_nufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,tmp_R_x_c_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_M_u^2)*dx^2;
tmp_R_k_q_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,tmp_R_k_p_);
tmp_R_x_c_rec0_ = interp_k_p_to_x_c_nufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_R_k_p_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_M_u^2) * n_w_sum;
tmp_R_k_p_rec0_ = interp_x_c_to_k_p_nufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,tmp_R_x_c_rec0_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_M_u^2)*dx^2;
tmp_R_x_c_rec1_ = interp_k_p_to_x_c_nufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_R_k_p_rec0_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_M_u^2) * n_w_sum;
tmp_R_k_p_rec1_ = interp_x_c_to_k_p_nufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,tmp_R_x_c_rec1_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_M_u^2)*dx^2;
tmp_R_x_c_rec2_ = interp_k_p_to_x_c_nufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_R_k_p_rec1_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_M_u^2) * n_w_sum;
tmp_R_k_p_wk__ = reshape(tmp_R_k_p_,[n_w_max,n_k_p_r]);
tmp_R_k_q_wk__ = reshape(tmp_R_k_q_,[n_w_max,n_k_p_r]);
tmp_UX_R_k_p_wn__ = tmp_R_k_p_wk__*diag(X_weight_r_)*UX__;
tmp_UX_R_k_q_wn__ = tmp_R_k_q_wk__*diag(X_weight_r_)*UX__;
tmp_UX_R_k_p_ = tmp_UX_R_k_p_wn__(:);
tmp_UX_R_k_q_ = tmp_UX_R_k_q_wn__(:);
tmp_R_x_c_l2 = sum(abs(tmp_R_x_c_).^2,'all')*dx^2;
tmp_R_x_c_rec0_l2 = sum(abs(tmp_R_x_c_rec0_).^2,'all')*dx^2;
tmp_R_x_c_rec1_l2 = sum(abs(tmp_R_x_c_rec1_).^2,'all')*dx^2;
tmp_R_x_c_rec2_l2 = sum(abs(tmp_R_x_c_rec2_).^2,'all')*dx^2;
tmp_R_k_p_l2 = sum(abs(tmp_R_k_p_).^2 .* weight_2d_k_all_) * (2*pi)^2;
tmp_R_k_q_l2 = sum(abs(tmp_R_k_q_).^2 .* weight_2d_k_all_) * (2*pi)^2;
tmp_R_k_p_rec0_l2 = sum(abs(tmp_R_k_p_rec0_).^2 .* weight_2d_k_all_) * (2*pi)^2;
tmp_R_k_p_rec1_l2 = sum(abs(tmp_R_k_p_rec1_).^2 .* weight_2d_k_all_) * (2*pi)^2;
tmp_UX_R_k_p_l2 = sum(abs(tmp_UX_R_k_p_).^2) * (2*pi)^2 / n_w_max / (2*pi)^2 ;
tmp_UX_R_k_q_l2 = sum(abs(tmp_UX_R_k_q_).^2) * (2*pi)^2 / n_w_max / (2*pi)^2 ;
%%%%;
disp(sprintf(' %% tmp_R_k_q_l2 = %0.16f',tmp_R_k_q_l2));
disp(sprintf(' %% tmp_R_k_p_l2 = %0.16f',tmp_R_k_p_l2));
disp(sprintf(' %% tmp_R_k_p_rec0_l2 = %0.16f',tmp_R_k_p_rec0_l2));
disp(sprintf(' %% tmp_R_k_p_rec1_l2 = %0.16f',tmp_R_k_p_rec1_l2));
disp(sprintf(' %% tmp_UX_R_k_q_l2 = %0.16f',tmp_UX_R_k_q_l2));
disp(sprintf(' %% tmp_UX_R_k_p_l2 = %0.16f',tmp_UX_R_k_p_l2));
disp(sprintf(' %% tmp_R_x_c_l2 = %0.16f',tmp_R_x_c_l2));
disp(sprintf(' %% tmp_R_x_c_rec0_l2 = %0.16f',tmp_R_x_c_rec0_l2));
disp(sprintf(' %% tmp_R_x_c_rec1_l2 = %0.16f',tmp_R_x_c_rec1_l2));
disp(sprintf(' %% tmp_R_x_c_rec2_l2 = %0.16f',tmp_R_x_c_rec2_l2));
disp(sprintf(' %% tmp_R_k_p_ vs tmp_R_k_p_rec0_: %0.16f',fnorm(tmp_R_k_p_ - tmp_R_k_p_rec0_)/fnorm(tmp_R_k_p_)));
disp(sprintf(' %% tmp_R_x_c_ vs tmp_R_x_c_rec0_: %0.16f',fnorm(tmp_R_x_c_ - tmp_R_x_c_rec0_)/fnorm(tmp_R_x_c_)));
disp(sprintf(' %% tmp_R_k_p_rec0_ vs tmp_R_k_p_rec1_: %0.16f',fnorm(tmp_R_k_p_rec0_ - tmp_R_k_p_rec1_)/fnorm(tmp_R_k_p_rec0_)));
disp(sprintf(' %% tmp_R_x_c_rec0_ vs tmp_R_x_c_rec1_: %0.16f',fnorm(tmp_R_x_c_rec0_ - tmp_R_x_c_rec1_)/fnorm(tmp_R_x_c_rec0_)));
disp(sprintf(' %% tmp_R_x_c_rec1_ vs tmp_R_x_c_rec2_: %0.16f',fnorm(tmp_R_x_c_rec1_ - tmp_R_x_c_rec2_)/fnorm(tmp_R_x_c_rec1_)));
%%%%;
% and one actual image. ;
%%%%;
nM=0;
tmp_M_x_c_ = M_x_c___(:,:,1+nM);
tmp_M_k_p_ = interp_x_c_to_k_p_nufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,tmp_M_x_c_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_M_u^2)*dx^2;
tmp_M_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,tmp_M_k_p_,-1*M_abs_x_c_0_avg_(1+nM),-1*M_abs_x_c_1_avg_(1+nM));
flag_test=0; if flag_test; tmp_M_k_p_ = tmp_M_k_p_form_; end;
tmp_M_k_q_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,tmp_M_k_p_);
tmp_M_x_c_rec0_ = interp_k_p_to_x_c_nufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_M_k_p_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_M_u^2) * n_w_sum;
tmp_M_k_p_rec0_ = interp_x_c_to_k_p_nufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,tmp_M_x_c_rec0_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_M_u^2)*dx^2;
tmp_M_x_c_rec1_ = interp_k_p_to_x_c_nufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_M_k_p_rec0_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_M_u^2) * n_w_sum;
tmp_M_k_p_rec1_ = interp_x_c_to_k_p_nufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,tmp_M_x_c_rec1_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_M_u^2)*dx^2;
tmp_M_x_c_rec2_ = interp_k_p_to_x_c_nufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_M_k_p_rec1_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_M_u^2) * n_w_sum;
tmp_M_k_p_wk__ = reshape(tmp_M_k_p_,[n_w_max,n_k_p_r]);
tmp_M_k_q_wk__ = reshape(tmp_M_k_q_,[n_w_max,n_k_p_r]);
tmp_UX_M_k_p_wn__ = tmp_M_k_p_wk__*diag(X_weight_r_)*UX__;
tmp_UX_M_k_q_wn__ = tmp_M_k_q_wk__*diag(X_weight_r_)*UX__;
tmp_UX_M_k_p_ = tmp_UX_M_k_p_wn__(:);
tmp_UX_M_k_q_ = tmp_UX_M_k_q_wn__(:);
tmp_M_x_c_l2 = sum(abs(tmp_M_x_c_).^2,'all')*dx^2;
tmp_M_x_c_rec0_l2 = sum(abs(tmp_M_x_c_rec0_).^2,'all')*dx^2;
tmp_M_x_c_rec1_l2 = sum(abs(tmp_M_x_c_rec1_).^2,'all')*dx^2;
tmp_M_x_c_rec2_l2 = sum(abs(tmp_M_x_c_rec2_).^2,'all')*dx^2;
tmp_M_k_p_l2 = sum(abs(tmp_M_k_p_).^2 .* weight_2d_k_all_) * (2*pi)^2;
tmp_M_k_q_l2 = sum(abs(tmp_M_k_q_).^2 .* weight_2d_k_all_) * (2*pi)^2;
tmp_M_k_p_rec0_l2 = sum(abs(tmp_M_k_p_rec0_).^2 .* weight_2d_k_all_) * (2*pi)^2;
tmp_M_k_p_rec1_l2 = sum(abs(tmp_M_k_p_rec1_).^2 .* weight_2d_k_all_) * (2*pi)^2;
tmp_UX_M_k_p_l2 = sum(abs(tmp_UX_M_k_p_).^2) * (2*pi)^2 / n_w_max / (2*pi)^2 ;
tmp_UX_M_k_q_l2 = sum(abs(tmp_UX_M_k_q_).^2) * (2*pi)^2 / n_w_max / (2*pi)^2 ;
%%%%;
disp(sprintf(' %% tmp_M_k_q_l2 = %0.16f',tmp_M_k_q_l2));
disp(sprintf(' %% tmp_M_k_p_l2 = %0.16f',tmp_M_k_p_l2));
disp(sprintf(' %% tmp_M_k_p_rec0_l2 = %0.16f',tmp_M_k_p_rec0_l2));
disp(sprintf(' %% tmp_M_k_p_rec1_l2 = %0.16f',tmp_M_k_p_rec1_l2));
disp(sprintf(' %% tmp_UX_M_k_q_l2 = %0.16f',tmp_UX_M_k_q_l2));
disp(sprintf(' %% tmp_UX_M_k_p_l2 = %0.16f',tmp_UX_M_k_p_l2));
disp(sprintf(' %% tmp_M_x_c_l2 = %0.16f',tmp_M_x_c_l2));
disp(sprintf(' %% tmp_M_x_c_rec0_l2 = %0.16f',tmp_M_x_c_rec0_l2));
disp(sprintf(' %% tmp_M_x_c_rec1_l2 = %0.16f',tmp_M_x_c_rec1_l2));
disp(sprintf(' %% tmp_M_x_c_rec2_l2 = %0.16f',tmp_M_x_c_rec2_l2));
disp(sprintf(' %% tmp_M_k_p_ vs tmp_M_k_p_rec0_: %0.16f',fnorm(tmp_M_k_p_ - tmp_M_k_p_rec0_)/fnorm(tmp_M_k_p_)));
disp(sprintf(' %% tmp_M_x_c_ vs tmp_M_x_c_rec0_: %0.16f',fnorm(tmp_M_x_c_ - tmp_M_x_c_rec0_)/fnorm(tmp_M_x_c_)));
disp(sprintf(' %% tmp_M_k_p_rec0_ vs tmp_M_k_p_rec1_: %0.16f',fnorm(tmp_M_k_p_rec0_ - tmp_M_k_p_rec1_)/fnorm(tmp_M_k_p_rec0_)));
disp(sprintf(' %% tmp_M_x_c_rec0_ vs tmp_M_x_c_rec1_: %0.16f',fnorm(tmp_M_x_c_rec0_ - tmp_M_x_c_rec1_)/fnorm(tmp_M_x_c_rec0_)));
disp(sprintf(' %% tmp_M_x_c_rec1_ vs tmp_M_x_c_rec2_: %0.16f',fnorm(tmp_M_x_c_rec1_ - tmp_M_x_c_rec2_)/fnorm(tmp_M_x_c_rec1_)));
%%%%;
% and one actual template. ;
%%%%;
nS=0;
tmp_CTF_S_k_p_ = CTF_S_k_p__(:,1+nS);
tmp_CTF_S_k_q_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,tmp_CTF_S_k_p_);
tmp_CTF_S_x_c_rec0_ = interp_k_p_to_x_c_nufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_CTF_S_k_p_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_M_u^2) * n_w_sum;
tmp_CTF_S_k_p_rec0_ = interp_x_c_to_k_p_nufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,tmp_CTF_S_x_c_rec0_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_M_u^2)*dx^2;
tmp_CTF_S_x_c_rec1_ = interp_k_p_to_x_c_nufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_CTF_S_k_p_rec0_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_M_u^2) * n_w_sum;
tmp_CTF_S_k_p_rec1_ = interp_x_c_to_k_p_nufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,tmp_CTF_S_x_c_rec1_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_M_u^2)*dx^2;
tmp_CTF_S_x_c_rec2_ = interp_k_p_to_x_c_nufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_CTF_S_k_p_rec1_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_M_u^2) * n_w_sum;
tmp_CTF_S_k_p_wk__ = reshape(tmp_CTF_S_k_p_,[n_w_max,n_k_p_r]);
tmp_CTF_S_k_q_wk__ = reshape(tmp_CTF_S_k_q_,[n_w_max,n_k_p_r]);
tmp_UX_CTF_S_k_p_wn__ = tmp_CTF_S_k_p_wk__*diag(X_weight_r_)*UX__;
tmp_UX_CTF_S_k_q_wn__ = tmp_CTF_S_k_q_wk__*diag(X_weight_r_)*UX__;
tmp_UX_CTF_S_k_p_ = tmp_UX_CTF_S_k_p_wn__(:);
tmp_UX_CTF_S_k_q_ = tmp_UX_CTF_S_k_q_wn__(:);
tmp_CTF_S_x_c_rec0_l2 = sum(abs(tmp_CTF_S_x_c_rec0_).^2,'all')*dx^2;
tmp_CTF_S_x_c_rec1_l2 = sum(abs(tmp_CTF_S_x_c_rec1_).^2,'all')*dx^2;
tmp_CTF_S_x_c_rec2_l2 = sum(abs(tmp_CTF_S_x_c_rec2_).^2,'all')*dx^2;
tmp_CTF_S_k_p_l2 = sum(abs(tmp_CTF_S_k_p_).^2 .* weight_2d_k_all_) * (2*pi)^2;
tmp_CTF_S_k_q_l2 = sum(abs(tmp_CTF_S_k_q_).^2 .* weight_2d_k_all_) * (2*pi)^2;
tmp_CTF_S_k_p_rec0_l2 = sum(abs(tmp_CTF_S_k_p_rec0_).^2 .* weight_2d_k_all_) * (2*pi)^2;
tmp_CTF_S_k_p_rec1_l2 = sum(abs(tmp_CTF_S_k_p_rec1_).^2 .* weight_2d_k_all_) * (2*pi)^2;
tmp_UX_CTF_S_k_p_l2 = sum(abs(tmp_UX_CTF_S_k_p_).^2) * (2*pi)^2 / n_w_max / (2*pi)^2 ;
tmp_UX_CTF_S_k_q_l2 = sum(abs(tmp_UX_CTF_S_k_q_).^2) * (2*pi)^2 / n_w_max / (2*pi)^2 ;
%%%%;
disp(sprintf(' %% tmp_CTF_S_k_q_l2 = %0.16f',tmp_CTF_S_k_q_l2));
disp(sprintf(' %% tmp_CTF_S_k_p_l2 = %0.16f',tmp_CTF_S_k_p_l2));
disp(sprintf(' %% tmp_CTF_S_k_p_rec0_l2 = %0.16f',tmp_CTF_S_k_p_rec0_l2));
disp(sprintf(' %% tmp_CTF_S_k_p_rec1_l2 = %0.16f',tmp_CTF_S_k_p_rec1_l2));
disp(sprintf(' %% tmp_UX_CTF_S_k_q_l2 = %0.16f',tmp_UX_CTF_S_k_q_l2));
disp(sprintf(' %% tmp_UX_CTF_S_k_p_l2 = %0.16f',tmp_UX_CTF_S_k_p_l2));
disp(sprintf(' %% tmp_CTF_S_x_c_rec0_l2 = %0.16f',tmp_CTF_S_x_c_rec0_l2));
disp(sprintf(' %% tmp_CTF_S_x_c_rec1_l2 = %0.16f',tmp_CTF_S_x_c_rec1_l2));
disp(sprintf(' %% tmp_CTF_S_x_c_rec2_l2 = %0.16f',tmp_CTF_S_x_c_rec2_l2));
disp(sprintf(' %% tmp_CTF_S_k_p_ vs tmp_CTF_S_k_p_rec0_: %0.16f',fnorm(tmp_CTF_S_k_p_ - tmp_CTF_S_k_p_rec0_)/fnorm(tmp_CTF_S_k_p_)));
disp(sprintf(' %% tmp_CTF_S_k_p_rec0_ vs tmp_CTF_S_k_p_rec1_: %0.16f',fnorm(tmp_CTF_S_k_p_rec0_ - tmp_CTF_S_k_p_rec1_)/fnorm(tmp_CTF_S_k_p_rec0_)));
disp(sprintf(' %% tmp_CTF_S_x_c_rec0_ vs tmp_CTF_S_x_c_rec1_: %0.16f',fnorm(tmp_CTF_S_x_c_rec0_ - tmp_CTF_S_x_c_rec1_)/fnorm(tmp_CTF_S_x_c_rec0_)));
disp(sprintf(' %% tmp_CTF_S_x_c_rec1_ vs tmp_CTF_S_x_c_rec2_: %0.16f',fnorm(tmp_CTF_S_x_c_rec1_ - tmp_CTF_S_x_c_rec2_)/fnorm(tmp_CTF_S_x_c_rec1_)));
%%%%;
flag_plot=1;
if flag_plot;
figure(1);figbig;figbeach();
ns=0;
subplot(3,5,1+ns);ns=ns+1;imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,real(tmp_R_x_c_)); axis image;axisnotick; title('tmp_R_x_c_','Interpreter','none');
subplot(3,5,1+ns);ns=ns+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(tmp_R_k_p_)); axis image;axisnotick; title('tmp_R_k_p_','Interpreter','none');
subplot(3,5,1+ns);ns=ns+1;imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(tmp_R_k_q_)); axis tight; axisnotick; title('tmp_R_k_q_','Interpreter','none');
subplot(3,5,1+ns);ns=ns+1;imagesc(real(tmp_UX_R_k_p_wn__)); axis tight; axisnotick; title('tmp_UX_R_k_p_','Interpreter','none');
subplot(3,5,1+ns);ns=ns+1;imagesc_q(n_UX_rank,1:n_UX_rank,n_w_,n_w_sum,real(tmp_UX_R_k_q_)); axis tight; axisnotick; title('tmp_UX_R_k_q_','Interpreter','none');
subplot(3,5,1+ns);ns=ns+1;imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,real(tmp_M_x_c_)); axis image;axisnotick; title('tmp_M_x_c_','Interpreter','none');
subplot(3,5,1+ns);ns=ns+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(tmp_M_k_p_)); axis image;axisnotick; title('tmp_M_k_p_','Interpreter','none');
subplot(3,5,1+ns);ns=ns+1;imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(tmp_M_k_q_)); axis tight; axisnotick; title('tmp_M_k_q_','Interpreter','none');
subplot(3,5,1+ns);ns=ns+1;imagesc(real(tmp_UX_M_k_p_wn__)); axis tight; axisnotick; title('tmp_UX_M_k_p_','Interpreter','none');
subplot(3,5,1+ns);ns=ns+1;imagesc_q(n_UX_rank,1:n_UX_rank,n_w_,n_w_sum,real(tmp_UX_M_k_q_)); axis tight; axisnotick; title('tmp_UX_M_k_q_','Interpreter','none');
subplot(3,5,1+ns);ns=ns+1;imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,real(tmp_CTF_S_x_c_rec0_)); axis image;axisnotick; title('tmp_CTF_S_x_c_rec0_','Interpreter','none');
subplot(3,5,1+ns);ns=ns+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(tmp_CTF_S_k_p_)); axis image;axisnotick; title('tmp_CTF_S_k_p_','Interpreter','none');
subplot(3,5,1+ns);ns=ns+1;imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(tmp_CTF_S_k_q_)); axis tight; axisnotick; title('tmp_CTF_S_k_q_','Interpreter','none');
subplot(3,5,1+ns);ns=ns+1;imagesc(real(tmp_UX_CTF_S_k_p_wn__)); axis tight; axisnotick; title('tmp_UX_CTF_S_k_p_','Interpreter','none');
subplot(3,5,1+ns);ns=ns+1;imagesc_q(n_UX_rank,1:n_UX_rank,n_w_,n_w_sum,real(tmp_UX_CTF_S_k_q_)); axis tight; axisnotick; title('tmp_UX_CTF_S_k_q_','Interpreter','none');
sgtitle(sprintf('random (top) vs image %d (middle) and CTF-x-template %d (bottom)',nM,nS));
end;%if flag_plot;
%%%%%%%%;
end;%if flag_check;

%%%%%%%%;
% Note that the actual images are slightly different from noise. ;
%%%%%%%%;
h_v_ = 8*linspace(-1,1,128);
h_M_x_c_ = hist(M_x_c___(:),h_v_);
h_R_x_c_ = hist(R_x_c___(:),h_v_);
if flag_plot;
figure(1);clf;figsml;
hold on;
stairs(h_v_,h_M_x_c_,'k');
stairs(h_v_,h_R_x_c_,'r');
hold off;
xlabel('value'); ylabel('number');
end;%if flag_plot;

if flag_plot;
%%%%%%%%;
% now look at the average variance (per degree-of-freedom) for random images. ;
% repeat for image-stack. ;
% repeat again for template-stack. ;
%%%%%%%%;
figure(1);figbig;figbeach();
ns=0;
subplot(3,5,1+ns);ns=ns+1;imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,var(R_x_c___,1,3),[0.85,1.15]); axis image;axisnotick; title('R_x_c_','Interpreter','none');
subplot(3,5,1+ns);ns=ns+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,var(R_k_p__,1,2),[0,2e-3]); axis image;axisnotick; title('R_k_p_','Interpreter','none');
subplot(3,5,1+ns);ns=ns+1;imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,var(R_k_q__,1,2),[0,2e-3]); axis tight; axisnotick; title('R_k_q_','Interpreter','none');
subplot(3,5,1+ns);ns=ns+1;imagesc(var(UX_R_k_p_wnM___,1,3),[0,5e-3]); axis tight; axisnotick; title('UX_R_k_p_','Interpreter','none');
subplot(3,5,1+ns);ns=ns+1;imagesc_q(n_UX_rank,1:n_UX_rank,n_w_,n_w_sum,var(UX_R_k_q_wnM___,1,3),[0,5e-3]); axis tight; axisnotick; title('UX_R_k_q_','Interpreter','none');
subplot(3,5,1+ns);ns=ns+1;imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,var(M_x_c___,1,3),[0.85,1.15]); axis image;axisnotick; title('M_x_c_','Interpreter','none');
subplot(3,5,1+ns);ns=ns+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,var(M_k_p__,1,2),[0,2e-3]); axis image;axisnotick; title('M_k_p_','Interpreter','none');
subplot(3,5,1+ns);ns=ns+1;imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,var(M_k_q__,1,2),[0,2e-3]); axis tight; axisnotick; title('M_k_q_','Interpreter','none');
subplot(3,5,1+ns);ns=ns+1;imagesc(var(UX_M_k_p_wnM___,1,3),[0,5e-3]); axis tight; axisnotick; title('UX_M_k_p_','Interpreter','none');
subplot(3,5,1+ns);ns=ns+1;imagesc_q(n_UX_rank,1:n_UX_rank,n_w_,n_w_sum,var(UX_M_k_q_wnM___,1,3),[0,5e-3]); axis tight; axisnotick; title('UX_M_k_q_','Interpreter','none');
%subplot(3,5,1+ns);ns=ns+1;imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,var(CTF_S_x_c___,1,3),[0.85,1.15]); axis image;axisnotick; title('CTF_S_x_c_','Interpreter','none');
ns=ns+1;
subplot(3,5,1+ns);ns=ns+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,var(CTF_S_k_p__,1,2),[0,2e-3]); axis image;axisnotick; title('CTF_S_k_p_','Interpreter','none');
subplot(3,5,1+ns);ns=ns+1;imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,var(CTF_S_k_q__,1,2),[0,2e-3]); axis tight; axisnotick; title('CTF_S_k_q_','Interpreter','none');
subplot(3,5,1+ns);ns=ns+1;imagesc(var(UX_CTF_S_k_p_wnS___,1,3),[0,5e-3]); axis tight; axisnotick; title('UX_CTF_S_k_p_','Interpreter','none');
subplot(3,5,1+ns);ns=ns+1;imagesc_q(n_UX_rank,1:n_UX_rank,n_w_,n_w_sum,var(UX_CTF_S_k_q_wnS___,1,3),[0,5e-3]); axis tight; axisnotick; title('UX_CTF_S_k_q_','Interpreter','none');
sgtitle(sprintf('random (top) vs image (middle) and CTF-x-template (bottom)'));
%%%%%%%%;
end;%if flag_plot;

%%%%%%%%;
% Now generate synthetic-images N_k_p_ (using the templates CTF_S_k_p_) with the same noise-level as the original images M_k_p_. ;
%%%%%%%%;
N_k_p__ = zeros(n_w_sum,n_M);
N_k_q__ = zeros(n_w_sum,n_M);
N_k_p_l2_ = zeros(n_M,1);
N_k_q_l2_ = zeros(n_M,1);
euler_polar_a_from_nM_ = zeros(n_M,1);
euler_azimu_b_from_nM_ = zeros(n_M,1);
euler_gamma_z_from_nM_ = zeros(n_M,1);
image_delta_x_from_nM_ = zeros(n_M,1);
image_delta_y_from_nM_ = zeros(n_M,1);
nS=0;
for nM=0:n_M-1;
N_k_p_ = CTF_S_k_p__(:,1+nS);
N_k_p_ = N_k_p_ + R_k_p__(:,1+nM);
N_k_p__(:,1+nM) = N_k_p_;
N_k_q_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,N_k_p_);
N_k_q__(:,1+nM) = N_k_q_;
N_k_p_l2_(1+nM) = sum(abs(N_k_p_).^2 .* weight_2d_k_all_) * (2*pi)^2;
N_k_q_l2_(1+nM) = sum(abs(N_k_q_).^2 .* weight_2d_k_all_) * (2*pi)^2;
euler_polar_a_from_nM_(1+nM) = viewing_polar_a_all_(1+nS);
euler_azimu_b_from_nM_(1+nM) = viewing_azimu_b_all_(1+nS);
euler_gamma_z_from_nM_(1+nM) = 0;
image_delta_x_from_nM_(1+nM) = 0;
image_delta_y_from_nM_(1+nM) = 0;
nS = nS+1; if (nS==n_S); nS=0; end;
end;%for nM=0:n_M-1;
wUX__ = diag(X_weight_r_)*UX__;
UX_N_k_p_wnM___ = zeros(n_w_max,n_UX_rank,n_M);
UX_N_k_q_wnM___ = zeros(n_w_max,n_UX_rank,n_M);
UX_N_k_p_l2_ = zeros(n_M,1);
UX_N_k_q_l2_ = zeros(n_M,1);
for nM=0:n_M-1;
tmp_N_k_p__ = reshape(N_k_p__(:,1+nM),[n_w_max,n_k_p_r]);
tmp_N_k_q__ = reshape(N_k_q__(:,1+nM),[n_w_max,n_k_p_r]);
UX_N_k_p_wnM___(:,:,1+nM) = reshape(N_k_p__(:,1+nM),[n_w_max,n_k_p_r])*wUX__;
UX_N_k_q_wnM___(:,:,1+nM) = reshape(N_k_q__(:,1+nM),[n_w_max,n_k_p_r])*wUX__;
UX_N_k_p_l2_(1+nM) = sum(abs(UX_N_k_p_wnM___(:,:,1+nM).^2),'all') * (2*pi)^2 / n_w_max / (2*pi)^2 ;
UX_N_k_q_l2_(1+nM) = sum(abs(UX_N_k_q_wnM___(:,:,1+nM).^2),'all') * (2*pi)^2 / n_w_max / (2*pi)^2 ;
end;%for nM=0:n_M-1;
%%%%%%%%;
% ensure that variance roughly matches. ;
%%%%%%%%;
if flag_plot;
figure(1);figbig;figbeach();
ns=0;
subplot(2,5,1+ns);ns=ns+1;imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,var(M_x_c___,1,3),[0.85,1.15]); axis image;axisnotick; title('M_x_c_','Interpreter','none');
subplot(2,5,1+ns);ns=ns+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,var(M_k_p__,1,2),[0,2e-3]); axis image;axisnotick; title('M_k_p_','Interpreter','none');
subplot(2,5,1+ns);ns=ns+1;imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,var(M_k_q__,1,2),[0,2e-3]); axis tight; axisnotick; title('M_k_q_','Interpreter','none');
subplot(2,5,1+ns);ns=ns+1;imagesc(var(UX_M_k_p_wnM___,1,3),[0,5e-3]); axis tight; axisnotick; title('UX_M_k_p_','Interpreter','none');
subplot(2,5,1+ns);ns=ns+1;imagesc_q(n_UX_rank,1:n_UX_rank,n_w_,n_w_sum,var(UX_M_k_q_wnM___,1,3),[0,5e-3]); axis tight; axisnotick; title('UX_M_k_q_','Interpreter','none');
%subplot(2,5,1+ns);ns=ns+1;imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,var(N_x_c___,1,3),[0.85,1.15]); axis image;axisnotick; title('N_x_c_','Interpreter','none');
ns=ns+1;
subplot(2,5,1+ns);ns=ns+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,var(N_k_p__,1,2),[0,2e-3]); axis image;axisnotick; title('N_k_p_','Interpreter','none');
subplot(2,5,1+ns);ns=ns+1;imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,var(N_k_q__,1,2),[0,2e-3]); axis tight; axisnotick; title('N_k_q_','Interpreter','none');
subplot(2,5,1+ns);ns=ns+1;imagesc(var(UX_N_k_p_wnM___,1,3),[0,5e-3]); axis tight; axisnotick; title('UX_N_k_p_','Interpreter','none');
subplot(2,5,1+ns);ns=ns+1;imagesc_q(n_UX_rank,1:n_UX_rank,n_w_,n_w_sum,var(UX_N_k_q_wnM___,1,3),[0,5e-3]); axis tight; axisnotick; title('UX_N_k_q_','Interpreter','none');
end;%if flag_plot;
%%%%%%%%;
% Now try reconstruction using all information and exact angles. ;
%%%%%%%%;
tmp_t = tic();
a_k_Y_reco_ ...
= ...
qbp_5(...
 1e-2 ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,n_M ...
,N_k_p__ ...
,zeros(n_M,1) ...
,CTF_avg_k_p_ ...
,euler_polar_a_from_nM_ ...
,euler_azimu_b_from_nM_ ...
,euler_gamma_z_from_nM_ ...
,image_delta_x_from_nM_ ...
,image_delta_y_from_nM_ ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% qbp_5: %0.6fs',tmp_t)); end;
disp(sprintf(' %% corr(a_k_Y_quad_,a_k_Y_reco_): %0.16f',corr(a_k_Y_quad_,a_k_Y_reco_)));
%%%%%%%%;
% Now form a_CTF_avg_UX_Y_quad__ ;
%%%%%%%%;
a_CTF_avg_UX_Y_quad__ = zeros(n_lm_max,n_UX_rank);
for nUX_rank=0:n_UX_rank-1;
for nk_p_r=0:n_k_p_r-1;
tmp_l_max = l_max_(1+nk_p_r);
tmp_n_lm = (tmp_l_max+1).^2;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:tmp_n_lm-1);
a_CTF_avg_UX_Y_quad__(1:tmp_n_lm,1+nUX_rank) = a_CTF_avg_UX_Y_quad__(1:tmp_n_lm,1+nUX_rank) + UX__(1+nk_p_r,1+nUX_rank)*X_weight_r_(1+nk_p_r)*a_k_Y_quad_(1+tmp_index_)*CTF_avg_k_p_r_(1+nk_p_r); %<-- use average CTF here, under the assumption that a_CTF_UX_Y_quad_ will be used alone. ;
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nUX_rank=0:n_UX_rank-1;
%%%%%%%%;
% now try pm reconstruction using exact angles. ;
%%%%%%%%;
FTK_0 = ampmh_FTK_1(n_k_p_r,k_p_r_,k_p_r_max,0,1e-3,1);
tmp_n_CTF_rank = 1;
%[tmp_UCTF_kc__,tmp_SCTF_c__,tmp_VCTF_Mc__] = svds(repmat(CTF_avg_k_p_r_,[1,n_M]),tmp_n_CTF_rank);
%tmp_VSCTF_Mc__ = tmp_VCTF_Mc__*tmp_SCTF_c__;
CTF_avg_k_p_r_norm = fnorm(CTF_avg_k_p_r_);
tmp_UCTF_kc__ = CTF_avg_k_p_r_/CTF_avg_k_p_r_norm;
tmp_VSCTF_Mc__ = CTF_avg_k_p_r_norm*ones(n_M,1);


flag_X = 0;
sigma_diff_ = 2.^[-2:2]; n_sigma_diff = numel(sigma_diff_)
tmp_corr_a_k_Y_diff_ns__ = zeros(n_UX_rank,n_sigma_diff);
tmp_corr_a_k_Y_reco_ns__ = zeros(n_UX_rank,n_sigma_diff);
tmp_corr_a_UCTF_UX_Y_diff_ns__ = zeros(n_UX_rank,n_sigma_diff);
tmp_corr_a_UCTF_UX_Y_reco_ns__ = zeros(n_UX_rank,n_sigma_diff);
tmp_X_a_UCTF_UX_Y_diff_ns__ = zeros(n_UX_rank,n_sigma_diff);
tmp_X_a_UCTF_UX_Y_reco_ns__ = zeros(n_UX_rank,n_sigma_diff);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for nsigma_diff=n_sigma_diff-1:-1:0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
sigma_diff = sigma_diff_(1+nsigma_diff); %<-- sqrt( diffusion time ). ;
tmp_corr_a_k_Y_diff_ = zeros(n_UX_rank,1);
tmp_corr_a_k_Y_reco_ = zeros(n_UX_rank,1);
tmp_corr_a_UCTF_UX_Y_diff_ = zeros(n_UX_rank,1);
tmp_corr_a_UCTF_UX_Y_reco_ = zeros(n_UX_rank,1);
tmp_X_a_UCTF_UX_Y_diff_ = zeros(n_UX_rank,1);
tmp_X_a_UCTF_UX_Y_reco_ = zeros(n_UX_rank,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for npm=0:n_UX_rank-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (verbose); disp(sprintf(' %% npm %d/%d',npm,n_UX_rank)); end;
tmp_pm_n_UX_rank = 1+npm;
tmp_pm_n_k_p_r = tmp_pm_n_UX_rank;
tmp_pm_k_p_r_ = ones(tmp_pm_n_k_p_r,1);
tmp_pm_n_w_ = n_w_max*ones(tmp_pm_n_k_p_r,1);
tmp_pm_n_w_sum = sum(n_w_max*ones(tmp_pm_n_k_p_r,1));
tmp_pm_l_max_ = l_max_max*ones(tmp_pm_n_UX_rank,1);
tmp_pm_n_lm_ = (1+tmp_pm_l_max_).^2;
tmp_pm_n_lm_max = max(tmp_pm_n_lm_);
tmp_pm_n_lm_sum = sum(tmp_pm_n_lm_);
tmp_svd_VUXN_lwnM____ = tpmh_VUXM_lwnM____3(FTK_0,n_k_p_r,n_w_,n_M,N_k_q__,tmp_pm_n_UX_rank,UX__,X_weight_r_);
tmp_UX_N_l2_dM__ = ampmh_UX_M_l2_dM__1(FTK_0,n_w_,n_M,tmp_pm_n_UX_rank,tmp_svd_VUXN_lwnM____);
[tmp_UX_N_k_q_wnM___,tmp_UX_N_k_p_wnM___] = ampmh_UX_M_k_p_wnM___0(FTK_0,n_w_,tmp_pm_n_UX_rank,n_M,tmp_svd_VUXN_lwnM____);
tmp_t = tic();
tmp_parameter = struct('type','parameter');
tmp_parameter.flag_qbp_vs_lsq = 1;
[ ...
 tmp_parameter ...
,tmp_a_UCTF_UX_Y_ync__ ... 
] = ...
a_UCTF_UX_Y_wrap_ync__0( ...
 tmp_parameter ...
,tmp_pm_n_k_p_r ...
,tmp_pm_l_max_ ...
,tmp_pm_n_w_ ...
,n_M ...
,reshape(tmp_UX_N_k_p_wnM___,[n_w_max*tmp_pm_n_k_p_r,n_M]) ...
,tmp_n_CTF_rank ...
,tmp_VSCTF_Mc__ ...
,euler_polar_a_from_nM_ ...
,euler_azimu_b_from_nM_ ...
,euler_gamma_z_from_nM_ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% tmp_a_UCTF_UX_Y_ync__: %0.3fs',tmp_t)); end;
tmp_parameter = parameter_timing_update(tmp_parameter,'a_UCTF_UX_Y_wrap_ync__0',tmp_t);
disp(sprintf(' %% corr(reshape(a_CTF_avg_UX_Y_quad__(:,1:tmp_pm_n_UX_rank),[tmp_pm_n_lm_sum,1]),tmp_a_UCTF_UX_Y_ync__): %0.16f',corr(reshape(a_CTF_avg_UX_Y_quad__(:,1:tmp_pm_n_UX_rank),[tmp_pm_n_lm_sum,1]),tmp_a_UCTF_UX_Y_ync__)));
%%%%%%%%;
% Now attempt similar with scattered angles. ;
% Because we do not want to correlated sets of angles with one another (which will be confounded by symmetries in the molecule), ;
% we measure success in terms of the correlation with a_k_Y_quad. ;
%%%%%%%%;
tmp_pm_Y_l_val_ = zeros(tmp_pm_n_lm_sum,1);
na=0;
for tmp_pm_nk_p_r=0:tmp_pm_n_k_p_r-1;
for l_val=0:l_max_max;
for m_val=-l_val:+l_val;
tmp_pm_Y_l_val_(1+na) = l_val;
na=na+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max_max;
end;%for tmp_pm_nk_p_r=0:tmp_pm_n_k_p_r-1;
assert(na==tmp_pm_n_lm_sum);
%%%%%%%%;
tmp_a_k_Y_diff_ = exp(-Y_l_val_.*(1+Y_l_val_)*sigma_diff^2).*a_k_Y_quad_;
tmp_corr_a_k_Y_diff = corr(tmp_a_k_Y_diff_,a_k_Y_quad_);
tmp_a_UCTF_UX_Y_diff_ync__ = exp(-tmp_pm_Y_l_val_.*(1+tmp_pm_Y_l_val_)*sigma_diff^2).*tmp_a_UCTF_UX_Y_ync__;
tmp_corr_a_UCTF_UX_Y_diff = corr(tmp_a_UCTF_UX_Y_ync__,tmp_a_UCTF_UX_Y_diff_ync__);
tmp_X_a_UCTF_UX_Y_diff = 0;
if flag_X;
[ ...
 tmp_X_a_UCTF_UX_Y_diff ...
,tmp_flag_flip ...
] = ...
register_spharm_to_spharm_wigner_wrap_1( ...
 tmp_pm_n_k_p_r ...
,ones(tmp_pm_n_k_p_r,1) ...
,1 ...
,ones(tmp_pm_n_k_p_r,1) ...
,0 ...
,tmp_pm_l_max_ ...
,tmp_a_UCTF_UX_Y_ync__ ...
,tmp_a_UCTF_UX_Y_diff_ync__ ...
);
end;%if flag_X;
disp(sprintf(' %% sigma_diff %0.2f, tmp_corr_a_k_Y_diff: %0.16f, tmp_corr_a_UCTF_UX_Y_diff: %0.16f, tmp_X_a_UCTF_UX_Y_diff %0.16f',sigma_diff,tmp_corr_a_k_Y_diff,tmp_corr_a_UCTF_UX_Y_diff,tmp_X_a_UCTF_UX_Y_diff));
tmp_corr_a_k_Y_diff_(1+npm) = tmp_corr_a_k_Y_diff;
tmp_corr_a_UCTF_UX_Y_diff_(1+npm) = tmp_corr_a_UCTF_UX_Y_diff;
tmp_X_a_UCTF_UX_Y_diff_(1+npm) = tmp_X_a_UCTF_UX_Y_diff;
%%%%%%%%;
% Use diffused principal-model to align principal-images. ;
%%%%%%%%;
tmp_t = tic();
tmp_parameter.template_viewing_k_eq_d = 1/k_p_r_max;
[ ...
 tmp_parameter ...
,tmp_n_S ...
,tmp_template_viewing_polar_a_all_ ...
,tmp_template_viewing_azimu_b_all_ ...
,tmp_X_SM__ ...
,tmp_delta_x_SM__ ...
,tmp_delta_y_SM__ ...
,tmp_gamma_z_SM__ ...
,tmp_I_value_SM__ ...
] = ...
ampmh_X_wrap_wrap_SM__8( ...
 tmp_parameter ...
,FTK_0 ...
,n_w_max ...
,l_max_max ...
,tmp_pm_n_UX_rank ...
,tmp_n_CTF_rank ...
,tmp_a_UCTF_UX_Y_diff_ync__ ...
,n_M ...
,zeros(n_M,1) ...
,tmp_VSCTF_Mc__ ...
,tmp_svd_VUXN_lwnM____ ...
,tmp_UX_N_l2_dM__ ...
,[] ...
,[] ...
,[] ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_SM__: %0.3fs',tmp_t)); end;
tmp_parameter = parameter_timing_update(tmp_parameter,'ampmh_X_wrap_wrap_SM__8',tmp_t);
%%%%%%%%;
% Use correlations to udate euler-angles. ;
%%%%%%%%;
tmp_t = tic();
tmp_parameter.flag_MS_vs_SM = 1;
[ ...
 tmp_parameter ...
,tmp_euler_polar_a_ ...
,tmp_euler_azimu_b_ ...
,tmp_euler_gamma_z_ ...
,tmp_image_delta_x_bit_ ...
,tmp_image_delta_y_bit_ ...
,tmp_image_I_value_ ...
,tmp_image_X_value_ ...
,tmp_image_S_index_ ...
] = ...
ampmh_MS_vs_SM_2( ...
 tmp_parameter ...
,n_w_max ...
,tmp_n_S ...
,tmp_template_viewing_polar_a_all_ ...
,tmp_template_viewing_azimu_b_all_ ...
,n_M ...
,tmp_X_SM__ ...
,tmp_delta_x_SM__ ...
,tmp_delta_y_SM__ ...
,tmp_gamma_z_SM__ ...
,tmp_I_value_SM__ ...
);
tmp_str = 'SM'; if (tmp_parameter.flag_MS_vs_SM); tmp_str = 'MS'; end;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% %s: update euler_polar_a_ euler_azimu_b_ euler_gamma_z_ : %0.3fs',tmp_str,tmp_t)); end;
tmp_parameter = parameter_timing_update(tmp_parameter,'ampmh_MS_vs_SM_2',tmp_t);
%%%%%%%%;
% now realign volume. ;
%%%%%%%%;
tmp_t = tic();
tmp_a_k_Y_reco_ ...
= ...
qbp_5(...
 1e-2 ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,n_M ...
,N_k_p__ ...
,zeros(n_M,1) ...
,CTF_avg_k_p_ ...
,tmp_euler_polar_a_ ...
,tmp_euler_azimu_b_ ...
,tmp_euler_gamma_z_ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% tmp_a_k_Y_reco_: %0.3fs',tmp_t)); end;
tmp_parameter = parameter_timing_update(tmp_parameter,'qbp_5',tmp_t);
%%%%%%%%;
tmp_corr_a_k_Y_reco = real(corr(tmp_a_k_Y_reco_,a_k_Y_quad_));
%%%%%%%%;
% now realign principal-volume. ;
%%%%%%%%;
tmp_t = tic();
tmp_parameter = struct('type','parameter');
tmp_parameter.flag_qbp_vs_lsq = 1;
[ ...
 tmp_parameter ...
,tmp_a_UCTF_UX_Y_reco_ync__ ... 
] = ...
a_UCTF_UX_Y_wrap_ync__0( ...
 tmp_parameter ...
,tmp_pm_n_k_p_r ...
,tmp_pm_l_max_ ...
,tmp_pm_n_w_ ...
,n_M ...
,reshape(tmp_UX_N_k_p_wnM___,[n_w_max*tmp_pm_n_k_p_r,n_M]) ...
,tmp_n_CTF_rank ...
,tmp_VSCTF_Mc__ ...
,tmp_euler_polar_a_ ...
,tmp_euler_azimu_b_ ...
,tmp_euler_gamma_z_ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% tmp_a_UCTF_UX_Y_reco_ync__: %0.3fs',tmp_t)); end;
tmp_parameter = parameter_timing_update(tmp_parameter,'a_UCTF_UX_Y_wrap_ync__0',tmp_t);
%%%%%%%%;
tmp_corr_a_UCTF_UX_Y_reco = real(corr(tmp_a_UCTF_UX_Y_reco_ync__,tmp_a_UCTF_UX_Y_ync__));
tmp_X_a_UCTF_UX_Y_reco = 0;
if flag_X;
[ ...
 tmp_X_a_UCTF_UX_Y_reco ...
,tmp_flag_flip ...
] = ...
register_spharm_to_spharm_wigner_wrap_1( ...
 tmp_pm_n_k_p_r ...
,ones(tmp_pm_n_k_p_r,1) ...
,1 ...
,ones(tmp_pm_n_k_p_r,1) ...
,0 ...
,tmp_pm_l_max_ ...
,tmp_a_UCTF_UX_Y_ync__ ...
,tmp_a_UCTF_UX_Y_reco_ync__ ...
);
end;%if flag_X;
disp(sprintf(' %% sigma_diff %0.2f, tmp_corr_a_k_Y_reco: %0.16f, tmp_corr_a_UCTF_UX_Y_reco: %0.16f, tmp_X_a_UCTF_UX_Y_reco %0.16f',sigma_diff,tmp_corr_a_k_Y_reco,tmp_corr_a_UCTF_UX_Y_reco,tmp_X_a_UCTF_UX_Y_reco));
tmp_corr_a_k_Y_reco_(1+npm) = tmp_corr_a_k_Y_reco;
tmp_corr_a_UCTF_UX_Y_reco_(1+npm) = tmp_corr_a_UCTF_UX_Y_reco;
tmp_X_a_UCTF_UX_Y_reco_(1+npm) = tmp_X_a_UCTF_UX_Y_reco;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for npm=0:n_UX_rank-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
tmp_corr_a_k_Y_diff_ns__(:,1+nsigma_diff) = tmp_corr_a_k_Y_diff_;
tmp_corr_a_k_Y_reco_ns__(:,1+nsigma_diff) = tmp_corr_a_k_Y_reco_;
tmp_corr_a_UCTF_UX_Y_diff_ns__(:,1+nsigma_diff) = tmp_corr_a_UCTF_UX_Y_diff_;
tmp_corr_a_UCTF_UX_Y_reco_ns__(:,1+nsigma_diff) = tmp_corr_a_UCTF_UX_Y_reco_;
tmp_X_a_UCTF_UX_Y_diff_ns__(:,1+nsigma_diff) = tmp_X_a_UCTF_UX_Y_diff_;
tmp_X_a_UCTF_UX_Y_reco_ns__(:,1+nsigma_diff) = tmp_X_a_UCTF_UX_Y_reco_;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for nsigma_diff=0:n_sigma_diff-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% verdict: not very clear cut. ;
%%%%%%%%;

