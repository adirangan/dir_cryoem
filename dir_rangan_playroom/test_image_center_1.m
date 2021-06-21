% test_pm_LetB1r_8;
verbose=1; nf=0;
n_w_max = 2*(l_max_max+1);
%%%%%%%%;
% Load the images. ;
%%%%%%%%;
disp(sprintf(' %% Warning, running rlnImageName_from_star_LetB1_0 rather than rlnImageName_from_star_0'));
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
rlnImageName_from_star_LetB1_0( ...
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
[M_x_c___(:,:,1+nM),n_M_ext_(1+nM)] = image_replace_edge_artefact_0(M_x_c___(:,:,1+nM),4,0.5,2,0);
end;%for nM=0:n_M-1;
disp(sprintf(' %% edge-artefacts detected in %d/%d images.',numel(find(n_M_ext_>0)),n_M));
%%%%%%%%;
% Now examine image-centroids. ;
%%%%%%%%;
n_x_M_u = size(M_x_c___,1);
assert(n_x_M_u==size(M_x_c___,2));
disp(sprintf(' %% typical edge-artefact covers %0.6f = (%0.6f)^2 of image.',median(n_M_ext_/n_x_M_u^2),median(sqrt(n_M_ext_/n_x_M_u^2))));
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
nM=0;
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
flag_plot=0;
if flag_plot;
figure(1+nf);nf=nf+1;figbig;figbeach();
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
flag_plot=0;
if flag_plot;
figure(1+nf);nf=nf+1;clf;figsml;
hold on;
stairs(h_v_,h_M_x_c_,'k');
stairs(h_v_,h_R_x_c_,'r');
hold off;
xlabel('value'); ylabel('number');
end;%if flag_plot;

flag_plot=0;
if flag_plot;
%%%%%%%%;
% now look at the average variance (per degree-of-freedom) for random images. ;
% repeat for image-stack. ;
% repeat again for template-stack. ;
%%%%%%%%;
figure(1+nf);nf=nf+1;figbig;figbeach();
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
% Now try and determine center from x_c_reco coordinates (filtered). ;
%%%%%%%%;
for nM=0:8-1;%for nM=0:n_M-1;
tmp_M_x_c_ = M_x_c___(:,:,1+nM);
tmp_M_k_p_ = interp_x_c_to_k_p_nufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,tmp_M_x_c_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_M_u^2)*dx^2;
tmp_M_x_c_rec0_ = interp_k_p_to_x_c_nufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_M_k_p_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_M_u^2) * n_w_sum;
%%%%%%%%;
flag_plot=0;
if flag_plot;
figure(1+nf);nf=nf+1;clf;figbig;figbeach();
ns=0;
subplot(2,2,1+ns);ns=ns+1;imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,real(tmp_M_x_c_),0.5*[-1,1]); axis image;axisnotick; title('M_x_c_','Interpreter','none');
subplot(2,2,1+ns);ns=ns+1;imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,real(tmp_M_x_c_rec0_),0.5*[-1,1]); axis image;axisnotick; title('M_x_c_rec0_','Interpreter','none');
subplot(2,2,1+ns);ns=ns+1;imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,real(tmp_M_x_c_rec1_),0.5*[-1,1]); axis image;axisnotick; title('M_x_c_rec1_','Interpreter','none');
end;%if flag_plot;
%%%%%%%%;
p_cut_ = 5:5:95; n_p_cut = numel(p_cut_);
for np_cut=0:n_p_cut-1;
p_cut = p_cut_(1+np_cut);
tmp_cut = prctile(real(tmp_M_x_c_rec0_(1+index_edge_M_x_c_)),p_cut);
M_rect_x_c_ = max(0,real(tmp_M_x_c_rec0_)-tmp_cut); %<-- no mask. ;
M_rect_avg = mean(M_rect_x_c_,'all');
M_rect_x_c_0_avg = mean(M_rect_x_c_/M_rect_avg.*x_c_0__,'all');
M_rect_x_c_1_avg = mean(M_rect_x_c_/M_rect_avg.*x_c_1__,'all');
M_rect_x_c_0_avg_(1+np_cut) = M_rect_x_c_0_avg;
M_rect_x_c_1_avg_(1+np_cut) = M_rect_x_c_1_avg;
M_mask_rect_x_c_ = max(0,real(tmp_M_x_c_rec0_.*x_c_mask__)-tmp_cut); %<-- radial mask. ;
M_mask_rect_avg = mean(M_mask_rect_x_c_,'all');
M_mask_rect_x_c_0_avg = mean(M_mask_rect_x_c_/M_mask_rect_avg.*x_c_0__,'all');
M_mask_rect_x_c_1_avg = mean(M_mask_rect_x_c_/M_mask_rect_avg.*x_c_1__,'all');
M_mask_rect_x_c_0_avg_(1+np_cut) = M_mask_rect_x_c_0_avg;
M_mask_rect_x_c_1_avg_(1+np_cut) = M_mask_rect_x_c_1_avg;
end;%for np_cut=0:n_p_cut-1;
np_cut = n_p_cut-1; %<-- use 95th percentile. ;
tmp_M_k_p_cent_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,tmp_M_k_p_,-1*M_mask_rect_x_c_0_avg_(1+np_cut),-1*M_mask_rect_x_c_1_avg_(1+np_cut));
tmp_M_x_c_rec1_ = interp_k_p_to_x_c_nufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_M_k_p_cent_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_M_u^2) * n_w_sum;
p_cut_ = 5:5:95; n_p_cut = numel(p_cut_);
for np_cut=0:n_p_cut-1;
p_cut = p_cut_(1+np_cut);
tmp_cut = prctile(real(tmp_M_x_c_rec1_(1+index_edge_M_x_c_)),p_cut);
M_mask_rec1_rect_x_c_ = max(0,real(tmp_M_x_c_rec1_.*x_c_mask__)-tmp_cut); %<-- radial mask. ;
M_mask_rec1_rect_avg = mean(M_mask_rec1_rect_x_c_,'all');
M_mask_rec1_rect_x_c_0_avg = mean(M_mask_rec1_rect_x_c_/M_mask_rec1_rect_avg.*x_c_0__,'all');
M_mask_rec1_rect_x_c_1_avg = mean(M_mask_rec1_rect_x_c_/M_mask_rec1_rect_avg.*x_c_1__,'all');
M_mask_rec1_rect_x_c_0_avg_(1+np_cut) = M_mask_rec1_rect_x_c_0_avg;
M_mask_rec1_rect_x_c_1_avg_(1+np_cut) = M_mask_rec1_rect_x_c_1_avg;
end;%for np_cut=0:n_p_cut-1;
%%%%%%%%;
flag_plot=1;
if flag_plot;
%figure(1+nf);nf=nf+1;clf;figmed;figbeach();
figure(1);clf;figmed;figbeach();
c_80s__ = colormap_80s(); n_c_80s = size(c_80s__,1);
ns=0;
subplot(1,3,1+ns);ns=ns+1;hold on;
imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,real(tmp_M_x_c_rec0_),0.5*[-1,1]);
for np_cut=0:n_p_cut-1;
nc_80s = max(0,min(n_c_80s-1,floor(n_c_80s*np_cut/n_p_cut)));
plot(M_rect_x_c_0_avg_(1+np_cut),M_rect_x_c_1_avg_(1+np_cut),'o','MarkerEdgeColor','w','MarkerFaceColor',c_80s__(1+nc_80s,:));
end;%for np_cut=0:n_p_cut-1;
axis image;axisnotick; title('M_rect_x_c_avg_','Interpreter','none');
subplot(1,3,1+ns);ns=ns+1;hold on;
imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,real(tmp_M_x_c_rec0_),0.5*[-1,1]);
for np_cut=0:n_p_cut-1;
nc_80s = max(0,min(n_c_80s-1,floor(n_c_80s*np_cut/n_p_cut)));
plot(M_mask_rect_x_c_0_avg_(1+np_cut),M_mask_rect_x_c_1_avg_(1+np_cut),'o','MarkerEdgeColor','w','MarkerFaceColor',c_80s__(1+nc_80s,:));
end;%for np_cut=0:n_p_cut-1;
axis image;axisnotick; title('M_mask_rect_x_c_avg_','Interpreter','none');
subplot(1,3,1+ns);ns=ns+1;hold on;
imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,real(tmp_M_x_c_rec1_),0.5*[-1,1]);
for np_cut=0:n_p_cut-1;
nc_80s = max(0,min(n_c_80s-1,floor(n_c_80s*np_cut/n_p_cut)));
plot(M_mask_rec1_rect_x_c_0_avg_(1+np_cut),M_mask_rec1_rect_x_c_1_avg_(1+np_cut),'o','MarkerEdgeColor','w','MarkerFaceColor',c_80s__(1+nc_80s,:));
end;%for np_cut=0:n_p_cut-1;
axis image;axisnotick; title('M_mask_rec1_rect_x_c_avg_','Interpreter','none');
sgtitle(sprintf(' %% centering filtered image %d',nM),'Interpreter','none');
end;%if flag_plot;
%%%%%%%%;
end;%for nM=0:n_M-1;

%%%%%%%%;
% Now systematize the centering above. ;
% main steps are: filter -> mask -> define noise -> rectify -> center -> repeat. ;
%%%%%%%%;
flag_plot=1;
if flag_plot;
%%%%%%%%;
figure(1+nf);nf=nf+1;figbig;figbeach();ns=0;
tmp_n_x = n_x_M_u/4;
tmp_x_c_0_ = -half_diameter_x_c + transpose([0:tmp_n_x-1]/tmp_n_x)*diameter_x_c;
tmp_x_c_1_ = -half_diameter_x_c + transpose([0:tmp_n_x-1]/tmp_n_x)*diameter_x_c;
for nM=12:24-1;%for nM=0:n_M-1;
tmp_M_k_p_ = M_k_p__(:,1+nM);
tmp_M_x_c_fil_ = interp_k_p_to_x_c_nufft(tmp_n_x,diameter_x_c,tmp_n_x,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_M_k_p_.*weight_2d_k_all_*(2*pi)^2)*sqrt(tmp_n_x^2) * n_w_sum;
tmp_M_x_c_lim_ = std(tmp_M_x_c_fil_,1,'all')*3.5*[-1,+1];
tmp_parameter = struct('type','parameter');
[ ...
 tmp_parameter ...
,tmp_M_k_p_out_ ...
,tmp_M_x_c_out_ ...
,tmp_delta_x_c_0 ...
,tmp_delta_x_c_1 ...
] = ...
image_center_0( ...
 tmp_parameter ...
,tmp_n_x ...
,diameter_x_c ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,tmp_M_k_p_ ...
,weight_2d_k_all_ ...
);
%%%%%%%%;
subplot(4,6,1+ns);ns=ns+1;
imagesc_c(tmp_n_x,tmp_x_c_0_,tmp_n_x,tmp_x_c_1_,real(tmp_M_x_c_fil_),tmp_M_x_c_lim_);
axis image;axisnotick; title(sprintf('M_x_c_fil_ %d',nM),'Interpreter','none');
subplot(4,6,1+ns);ns=ns+1;
imagesc_c(tmp_n_x,tmp_x_c_0_,tmp_n_x,tmp_x_c_1_,real(tmp_M_x_c_out_),tmp_M_x_c_lim_);
axis image;axisnotick; title(sprintf('M_x_c_out_ %d',nM),'Interpreter','none');
drawnow();
if (ns==4*6); ns=0; end;
%%%%%%%%;
end;%for nM=0:n_M-1;
%%%%%%%%;
end;%if flag_plot;


