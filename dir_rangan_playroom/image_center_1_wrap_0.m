%%%%%%%%;
n_w_max = 2*(l_max_max+1);
%%%%%%%%;
% Now load images and CTF parameters from the star-file. ;
%%%%%%%%;
disp(sprintf(' %% %s not found, creating',fname_mat));
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
if flag_invert; M_x_c___ = -M_x_c___; end;
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
fname_fig = sprintf('%s_jpg/M_abs_x_c_avg_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figmed;
%%%%;
subplot(1,2,1);
plot(M_abs_x_c_0_avg_,M_abs_x_c_1_avg_,'.');
xlabel('M_abs_x_c_0_avg_','Interpreter','none');
ylabel('M_abs_x_c_1_avg_','Interpreter','none');
axis equal; grid on;
title('M_abs_x_c_ (no radial mask)','Interpreter','none');
%%%%;
subplot(1,2,2);
plot(M_mask_abs_x_c_0_avg_,M_mask_abs_x_c_1_avg_,'.');
xlabel('M_mask_abs_x_c_0_avg_','Interpreter','none');
ylabel('M_mask_abs_x_c_1_avg_','Interpreter','none');
axis equal; grid on;
title('M_mask_abs_x_c_ (max radial mask)','Interpreter','none');
%%%%;
tmp_corr_ = corr([M_abs_x_c_0_avg_,M_abs_x_c_1_avg_],[M_mask_abs_x_c_0_avg_,M_mask_abs_x_c_1_avg_]);
sgtitle(sprintf('correlation x,y = ( %0.4f , %0.4f )',tmp_corr_(1,1),tmp_corr_(2,2)));
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));

%%%%%%%%;
% Now convert images to M_k_p__. ;
%%%%%%%%;
dx = diameter_x_c/n_x_M_u;
%%%%;
M_x_c_rec0___ = zeros(n_x_u,n_x_u,n_M);
M_k_p__ = zeros(n_w_sum,n_M);
M_k_q__ = zeros(n_w_sum,n_M);
M_k_p_l2_ = zeros(n_M,1);
M_k_q_l2_ = zeros(n_M,1);
for nM=0:n_M-1;
if (mod(nM,128)==0); disp(sprintf(' %% nM %d/%d',nM,n_M)); end;
M_x_c_ = squeeze(M_x_c___(:,:,1+nM));
M_k_p_ = interp_x_c_to_k_p_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,M_x_c_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_M_u^2)*dx^2 ;
% Now *DO* translate according to M_abs_x_c_0_avg_ and M_abs_x_c_1_avg_. ;
M_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p_,-1*M_abs_x_c_0_avg_(1+nM),-1*M_abs_x_c_1_avg_(1+nM));
M_k_p__(:,1+nM) = M_k_p_;
M_x_c_rec0_ = interp_k_p_to_x_c_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,M_k_p_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_M_u^2) * n_w_sum;
M_x_c_rec0___(:,:,1+nM) = reshape(M_x_c_rec0_,[n_x_u,n_x_u]);
M_k_q_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_k_p_);
M_k_q__(:,1+nM) = M_k_q_;
M_k_p_l2_(1+nM) = sum(abs(M_k_p_).^2 .* weight_2d_k_all_) * (2*pi)^2;
M_k_q_l2_(1+nM) = sum(abs(M_k_q_).^2 .* weight_2d_k_all_) * (2*pi)^2;
end;%for nM=0:n_M-1;
%%%%%%%%;
% do the same for a set of random images. ;
%%%%%%%%;
x_p_r__ = sqrt(x_c_0__.^2 + x_c_1__.^2);
edge_M_x_c_ = x_p_r__> 1/sqrt(2);
cent_M_x_c_ = x_p_r__<=1/sqrt(2);
index_edge_M_x_c_ = efind(edge_M_x_c_);
index_cent_M_x_c_ = efind(cent_M_x_c_);
%%%%;
% visualize with: ;
% tmp_ = zeros(n_x_M_u); tmp_(1+index_edge_M_x_c_) = 1; imagesc(tmp_);axis image; axisnotick;
%%%%;
R_x_c___ = zeros(n_x_M_u,n_x_M_u,n_M);
R_x_c_rec0___ = zeros(n_x_M_u,n_x_M_u,n_M);
for nM=0:n_M-1;
tmp_M_x_c__ = M_x_c___(:,:,1+nM);
tmp_M_avg = mean(tmp_M_x_c__(1+index_edge_M_x_c_));
tmp_M_std = std(tmp_M_x_c__(1+index_edge_M_x_c_),1);
rng(nM); R_x_c___(:,:,1+nM) = tmp_M_avg + tmp_M_std*randn(n_x_M_u,n_x_M_u);
end;%for nM=0:n_M-1;
R_k_p__ = zeros(n_w_sum,n_M);
R_k_q__ = zeros(n_w_sum,n_M);
R_k_p_l2_ = zeros(n_M,1);
R_k_q_l2_ = zeros(n_M,1);
for nM=0:n_M-1;
if (mod(nM,128)==0); disp(sprintf(' %% nM %d/%d',nM,n_M)); end;
R_x_c_ = squeeze(R_x_c___(:,:,1+nM));
R_k_p_ = interp_x_c_to_k_p_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,R_x_c_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_M_u^2)*dx^2 ;
R_k_p__(:,1+nM) = R_k_p_;
R_x_c_rec0_ = interp_k_p_to_x_c_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,R_k_p_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_M_u^2) * n_w_sum;
R_x_c_rec0___(:,:,1+nM) = reshape(R_x_c_rec0_,[n_x_u,n_x_u]);
R_k_q_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,R_k_p_);
R_k_q__(:,1+nM) = R_k_q_;
R_k_p_l2_(1+nM) = sum(abs(R_k_p_).^2 .* weight_2d_k_all_) * (2*pi)^2;
R_k_q_l2_(1+nM) = sum(abs(R_k_q_).^2 .* weight_2d_k_all_) * (2*pi)^2;
end;%for nM=0:n_M-1;

flag_check = 1;
if flag_check;
%%%%;
% assess random image. ;
%%%%;
nM=0;
tmp_R_x_c_ = R_x_c___(:,:,1+nM);
tmp_R_k_p_ = interp_x_c_to_k_p_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,tmp_R_x_c_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_M_u^2)*dx^2;
tmp_R_k_q_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,tmp_R_k_p_);
tmp_R_x_c_rec0_ = interp_k_p_to_x_c_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_R_k_p_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_M_u^2) * n_w_sum;
tmp_R_k_p_rec0_ = interp_x_c_to_k_p_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,tmp_R_x_c_rec0_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_M_u^2)*dx^2;
tmp_R_x_c_rec1_ = interp_k_p_to_x_c_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_R_k_p_rec0_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_M_u^2) * n_w_sum;
tmp_R_k_p_rec1_ = interp_x_c_to_k_p_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,tmp_R_x_c_rec1_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_M_u^2)*dx^2;
tmp_R_x_c_rec2_ = interp_k_p_to_x_c_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_R_k_p_rec1_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_M_u^2) * n_w_sum;
tmp_R_k_p_wk__ = reshape(tmp_R_k_p_,[n_w_max,n_k_p_r]);
tmp_R_k_q_wk__ = reshape(tmp_R_k_q_,[n_w_max,n_k_p_r]);
tmp_R_x_c_l2 = sum(abs(tmp_R_x_c_).^2,'all')*dx^2;
tmp_R_x_c_rec0_l2 = sum(abs(tmp_R_x_c_rec0_).^2,'all')*dx^2;
tmp_R_x_c_rec1_l2 = sum(abs(tmp_R_x_c_rec1_).^2,'all')*dx^2;
tmp_R_x_c_rec2_l2 = sum(abs(tmp_R_x_c_rec2_).^2,'all')*dx^2;
tmp_R_k_p_l2 = sum(abs(tmp_R_k_p_).^2 .* weight_2d_k_all_) * (2*pi)^2;
tmp_R_k_q_l2 = sum(abs(tmp_R_k_q_).^2 .* weight_2d_k_all_) * (2*pi)^2;
tmp_R_k_p_rec0_l2 = sum(abs(tmp_R_k_p_rec0_).^2 .* weight_2d_k_all_) * (2*pi)^2;
tmp_R_k_p_rec1_l2 = sum(abs(tmp_R_k_p_rec1_).^2 .* weight_2d_k_all_) * (2*pi)^2;
%%%%;
disp(sprintf(' %% tmp_R_k_q_l2 = %0.16f',tmp_R_k_q_l2));
disp(sprintf(' %% tmp_R_k_p_l2 = %0.16f',tmp_R_k_p_l2));
disp(sprintf(' %% tmp_R_k_p_rec0_l2 = %0.16f',tmp_R_k_p_rec0_l2));
disp(sprintf(' %% tmp_R_k_p_rec1_l2 = %0.16f',tmp_R_k_p_rec1_l2));
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
tmp_M_k_p_ = interp_x_c_to_k_p_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,tmp_M_x_c_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_M_u^2)*dx^2;
tmp_M_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,tmp_M_k_p_,-1*M_abs_x_c_0_avg_(1+nM),-1*M_abs_x_c_1_avg_(1+nM));
flag_test=0; if flag_test; tmp_M_k_p_ = tmp_M_k_p_form_; end;
tmp_M_k_q_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,tmp_M_k_p_);
tmp_M_x_c_rec0_ = interp_k_p_to_x_c_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_M_k_p_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_M_u^2) * n_w_sum;
tmp_M_k_p_rec0_ = interp_x_c_to_k_p_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,tmp_M_x_c_rec0_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_M_u^2)*dx^2;
tmp_M_x_c_rec1_ = interp_k_p_to_x_c_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_M_k_p_rec0_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_M_u^2) * n_w_sum;
tmp_M_k_p_rec1_ = interp_x_c_to_k_p_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,tmp_M_x_c_rec1_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_M_u^2)*dx^2;
tmp_M_x_c_rec2_ = interp_k_p_to_x_c_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_M_k_p_rec1_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_M_u^2) * n_w_sum;
tmp_M_k_p_wk__ = reshape(tmp_M_k_p_,[n_w_max,n_k_p_r]);
tmp_M_k_q_wk__ = reshape(tmp_M_k_q_,[n_w_max,n_k_p_r]);
tmp_M_x_c_l2 = sum(abs(tmp_M_x_c_).^2,'all')*dx^2;
tmp_M_x_c_rec0_l2 = sum(abs(tmp_M_x_c_rec0_).^2,'all')*dx^2;
tmp_M_x_c_rec1_l2 = sum(abs(tmp_M_x_c_rec1_).^2,'all')*dx^2;
tmp_M_x_c_rec2_l2 = sum(abs(tmp_M_x_c_rec2_).^2,'all')*dx^2;
tmp_M_k_p_l2 = sum(abs(tmp_M_k_p_).^2 .* weight_2d_k_all_) * (2*pi)^2;
tmp_M_k_q_l2 = sum(abs(tmp_M_k_q_).^2 .* weight_2d_k_all_) * (2*pi)^2;
tmp_M_k_p_rec0_l2 = sum(abs(tmp_M_k_p_rec0_).^2 .* weight_2d_k_all_) * (2*pi)^2;
tmp_M_k_p_rec1_l2 = sum(abs(tmp_M_k_p_rec1_).^2 .* weight_2d_k_all_) * (2*pi)^2;
%%%%;
disp(sprintf(' %% tmp_M_k_q_l2 = %0.16f',tmp_M_k_q_l2));
disp(sprintf(' %% tmp_M_k_p_l2 = %0.16f',tmp_M_k_p_l2));
disp(sprintf(' %% tmp_M_k_p_rec0_l2 = %0.16f',tmp_M_k_p_rec0_l2));
disp(sprintf(' %% tmp_M_k_p_rec1_l2 = %0.16f',tmp_M_k_p_rec1_l2));
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
fname_fig = sprintf('%s_jpg/R_vs_M_single_FIGA',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;figbig;figbeach();
ns=0;
subplot(2,4,1+ns);ns=ns+1;imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,real(tmp_R_x_c_)); axis image;axisnotick; title('tmp_R_x_c_','Interpreter','none');
subplot(2,4,1+ns);ns=ns+1;imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,real(tmp_R_x_c_rec0_)); axis image;axisnotick; title('tmp_R_x_c_rec0_','Interpreter','none');
subplot(2,4,1+ns);ns=ns+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(tmp_R_k_p_)); axis image;axisnotick; title('tmp_R_k_p_','Interpreter','none');
subplot(2,4,1+ns);ns=ns+1;imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(tmp_R_k_q_)); axis tight; axisnotick; title('tmp_R_k_q_','Interpreter','none');
subplot(2,4,1+ns);ns=ns+1;imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,real(tmp_M_x_c_)); axis image;axisnotick; title('tmp_M_x_c_','Interpreter','none');
subplot(2,4,1+ns);ns=ns+1;imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,real(tmp_M_x_c_rec0_)); axis image;axisnotick; title('tmp_M_x_c_rec0_','Interpreter','none');
subplot(2,4,1+ns);ns=ns+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(tmp_M_k_p_)); axis image;axisnotick; title('tmp_M_k_p_','Interpreter','none');
subplot(2,4,1+ns);ns=ns+1;imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(tmp_M_k_q_)); axis tight; axisnotick; title('tmp_M_k_q_','Interpreter','none');
sgtitle(sprintf('random (top) vs image %d (bottom)',nM));
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
end;%if flag_check;

%%%%%%%%;
% Note that the actual images are slightly different from noise. ;
%%%%%%%%;
h_v_ = 8*linspace(-1,1,128);
h_M_x_c_ = hist(M_x_c___(:),h_v_);
h_R_x_c_ = hist(R_x_c___(:),h_v_);
%%%%;
fname_fig = sprintf('%s_jpg/R_vs_M_hist_FIGA',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,2,1);
hold on;
stairs(h_v_,h_M_x_c_,'k');
stairs(h_v_,h_R_x_c_,'r');
hold off;
xlabel('value'); ylabel('number');
title('h_M_x_c_ (black) vs h_R_x_c_ (red)','Interpreter','none');
subplot(1,2,2);
hold on;
stairs(h_v_,h_M_x_c_ - h_R_x_c_,'k');
hold off;
xlabel('value'); ylabel('number'); title('difference');
sgtitle(sprintf('histogram of image-values'));
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));

%%%%%%%%;
% now look at the average variance (per degree-of-freedom) for random images. ;
% repeat for image-stack. ;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/R_vs_M_variance_FIGA',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
M_x_c_var_lim_ = prctile(var(M_x_c___,1,3),[ 5,95],'all');
M_x_c_rec0_var_lim_ = prctile(var(M_x_c_rec0___,1,3),[ 5,95],'all');
M_k_p_var_lim_ = prctile(var(M_k_p__,1,2),[ 5,95],'all');
M_k_q_var_lim_ = prctile(var(M_k_q__,1,2),[ 5,95],'all');
figure(1+nf);nf=nf+1;figbig;figbeach();
ns=0;
subplot(2,4,1+ns);ns=ns+1;imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,var(R_x_c___,1,3),M_x_c_var_lim_); axis image;axisnotick; title('R_x_c_','Interpreter','none');
subplot(2,4,1+ns);ns=ns+1;imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,var(R_x_c_rec0___,1,3),M_x_c_rec0_var_lim_); axis image;axisnotick; title('R_x_c_rec0_','Interpreter','none');
subplot(2,4,1+ns);ns=ns+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,var(R_k_p__,1,2),M_k_p_var_lim_); axis image;axisnotick; title('R_k_p_','Interpreter','none');
subplot(2,4,1+ns);ns=ns+1;imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,var(R_k_q__,1,2),M_k_q_var_lim_); axis tight; axisnotick; title('R_k_q_','Interpreter','none');
subplot(2,4,1+ns);ns=ns+1;imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,var(M_x_c___,1,3),M_x_c_var_lim_); axis image;axisnotick; title('M_x_c_','Interpreter','none');
subplot(2,4,1+ns);ns=ns+1;imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,var(M_x_c_rec0___,1,3),M_x_c_rec0_var_lim_); axis image;axisnotick; title('M_x_c_rec0_','Interpreter','none');
subplot(2,4,1+ns);ns=ns+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,var(M_k_p__,1,2),M_k_p_var_lim_); axis image;axisnotick; title('M_k_p_','Interpreter','none');
subplot(2,4,1+ns);ns=ns+1;imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,var(M_k_q__,1,2),M_k_q_var_lim_); axis tight; axisnotick; title('M_k_q_','Interpreter','none');
sgtitle(sprintf('random (top) vs image (bottom)'));
%%%%%%%%;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));

%%%%%%%%;
% Now try and determine center from x_c_reco coordinates (filtered). ;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/M_sample_center_FIGA',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;figbig;figbeach();
%%%%%%%%;
for nM=0:6-1;%for nM=0:n_M-1;
tmp_M_x_c_ = M_x_c___(:,:,1+nM);
tmp_M_k_p_ = interp_x_c_to_k_p_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,tmp_M_x_c_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_M_u^2)*dx^2;
tmp_M_x_c_rec0_ = interp_k_p_to_x_c_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_M_k_p_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_M_u^2) * n_w_sum;
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
tmp_M_x_c_rec1_ = interp_k_p_to_x_c_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_M_k_p_cent_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_M_u^2) * n_w_sum;
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
tmp_M_x_c_lim_ = mean(real(tmp_M_x_c_rec0_),'all') + std(real(tmp_M_x_c_rec0_),1,'all')*2.5*[-1,+1];
c_80s__ = colormap_80s(); n_c_80s = size(c_80s__,1);
ns=0;
subplot(3,6,1+ns*6+nM);ns=ns+1;
hold on;
imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,real(tmp_M_x_c_rec0_),tmp_M_x_c_lim_);
for np_cut=0:n_p_cut-1;
nc_80s = max(0,min(n_c_80s-1,floor(n_c_80s*np_cut/n_p_cut)));
plot(M_rect_x_c_0_avg_(1+np_cut),M_rect_x_c_1_avg_(1+np_cut),'o','MarkerEdgeColor','w','MarkerFaceColor',c_80s__(1+nc_80s,:));
end;%for np_cut=0:n_p_cut-1;
axis image;axisnotick; title(sprintf('M_rect_x_c_avg_ %d',nM),'Interpreter','none');
subplot(3,6,1+ns*6+nM);ns=ns+1;
hold on;
imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,real(tmp_M_x_c_rec0_),tmp_M_x_c_lim_);
for np_cut=0:n_p_cut-1;
nc_80s = max(0,min(n_c_80s-1,floor(n_c_80s*np_cut/n_p_cut)));
plot(M_mask_rect_x_c_0_avg_(1+np_cut),M_mask_rect_x_c_1_avg_(1+np_cut),'o','MarkerEdgeColor','w','MarkerFaceColor',c_80s__(1+nc_80s,:));
end;%for np_cut=0:n_p_cut-1;
axis image;axisnotick; title(sprintf('M_mask_rect_x_c_avg_ %d',nM),'Interpreter','none');
subplot(3,6,1+ns*6+nM);ns=ns+1;
hold on;
imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,real(tmp_M_x_c_rec1_),tmp_M_x_c_lim_);
for np_cut=0:n_p_cut-1;
nc_80s = max(0,min(n_c_80s-1,floor(n_c_80s*np_cut/n_p_cut)));
plot(M_mask_rec1_rect_x_c_0_avg_(1+np_cut),M_mask_rec1_rect_x_c_1_avg_(1+np_cut),'o','MarkerEdgeColor','w','MarkerFaceColor',c_80s__(1+nc_80s,:));
end;%for np_cut=0:n_p_cut-1;
axis image;axisnotick; title(sprintf('M_mask_rec1_rect_x_c_avg_ %d',nM),'Interpreter','none');
drawnow();
%%%%%%%%;
end;%for nM=0:n_M-1;
%%%%%%%%;
sgtitle(sprintf(' %% centering filtered images'),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));

%%%%%%%%;
% Now systematize the centering above. ;
% main steps are: filter -> mask -> define noise -> rectify -> center -> repeat. ;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/M_sample_center_FIGB',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;figbig;figbeach();ns=0;
%%%%%%%%;
tmp_n_x = n_x_M_u/4;
tmp_x_c_0_ = -half_diameter_x_c + transpose([0:tmp_n_x-1]/tmp_n_x)*diameter_x_c;
tmp_x_c_1_ = -half_diameter_x_c + transpose([0:tmp_n_x-1]/tmp_n_x)*diameter_x_c;
for nM=12:24-1;%for nM=0:n_M-1;
tmp_M_k_p_ = M_k_p__(:,1+nM);
tmp_M_x_c_fil_ = interp_k_p_to_x_c_xxnufft(tmp_n_x,diameter_x_c,tmp_n_x,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_M_k_p_.*weight_2d_k_all_*(2*pi)^2)*sqrt(tmp_n_x^2) * n_w_sum;
tmp_M_x_c_lim_ = std(tmp_M_x_c_fil_,1,'all')*3.5*[-1,+1];
tmp_parameter = struct('type','parameter');
[ ...
 tmp_parameter ...
,tmp_M_k_p_out_ ...
,tmp_M_x_c_0in_ ...
,tmp_M_x_c_out_ ...
,tmp_delta_x_c_0 ...
,tmp_delta_x_c_1 ...
] = ...
image_center_1( ...
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
sgtitle(sprintf(' %% centering filtered images'),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));

%{
%%%%%%%%;
% Now recenter all the images. ;
% M_k_p_ stores images centered on M_abs_x_c_?_avg_, ;
% N_k_p_ stores images centered on M_abs_x_c_?_avg_ as well as image_center_1. ;
%%%%%%%%;
dx = diameter_x_c/n_x_M_u;
M_k_p__ = zeros(n_w_sum,n_M);
N_k_p__ = zeros(n_w_sum,n_M);
image_center_delta_x_c_0_ = zeros(n_M,1);
image_center_delta_x_c_1_ = zeros(n_M,1);
tmp_n_x = max(n_w_max,n_x_M_u/4);
if (verbose); disp(sprintf(' %% tmp_n_x %d',tmp_n_x)); end;
M_x_c_avg__ = zeros(tmp_n_x,tmp_n_x);
M_x_c_std__ = zeros(tmp_n_x,tmp_n_x);
N_x_c_avg__ = zeros(tmp_n_x,tmp_n_x);
N_x_c_std__ = zeros(tmp_n_x,tmp_n_x);
for nM=0:n_M-1;
if (mod(nM,128)==0); if (verbose); disp(sprintf(' %% nM %d/%d',nM,n_M)); end; end;
M_x_c_ = squeeze(M_x_c___(:,:,1+nM));
M_k_p_ = interp_x_c_to_k_p_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,M_x_c_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_M_u^2)*dx^2 ;
% Now *DO* translate according to M_abs_x_c_0_avg_ and M_abs_x_c_1_avg_. ;
M_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p_,-1*M_abs_x_c_0_avg_(1+nM),-1*M_abs_x_c_1_avg_(1+nM));
% Now *ALSO* center images after filtering/masking/thresholding: ;
tmp_parameter = struct('type','parameter');
[ ...
 tmp_parameter ...
,N_k_p_ ...
,M_x_c_0in_ ...
,M_x_c_out_ ...
,tmp_delta_x_c_0 ...
,tmp_delta_x_c_1 ...
] = ...
image_center_1( ...
 tmp_parameter ...
,tmp_n_x ...
,diameter_x_c ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,M_k_p_ ...
,weight_2d_k_all_ ...
);
M_k_p__(:,1+nM) = M_k_p_;
N_k_p__(:,1+nM) = N_k_p_;
image_center_delta_x_c_0_(1+nM) = tmp_delta_x_c_0;
image_center_delta_x_c_1_(1+nM) = tmp_delta_x_c_1;
M_x_c_avg__ = M_x_c_avg__ + reshape(real(M_x_c_0in_),[tmp_n_x,tmp_n_x]);
M_x_c_std__ = M_x_c_std__ + reshape(real(M_x_c_0in_).^2,[tmp_n_x,tmp_n_x]);
N_x_c_avg__ = N_x_c_avg__ + reshape(real(M_x_c_out_),[tmp_n_x,tmp_n_x]);
N_x_c_std__ = N_x_c_std__ + reshape(real(M_x_c_out_).^2,[tmp_n_x,tmp_n_x]);
end;%for nM=0:n_M-1;
M_x_c_avg__ = M_x_c_avg__/n_M;
M_x_c_std__ = M_x_c_std__/n_M; M_x_c_std__ = sqrt(M_x_c_std__ - M_x_c_avg__.^2);
N_x_c_avg__ = N_x_c_avg__/n_M;
N_x_c_std__ = N_x_c_std__/n_M; N_x_c_std__ = sqrt(N_x_c_std__ - N_x_c_avg__.^2);
%%%%%%%%;
fname_fig = sprintf('%s_jpg/M_all_center_FIGB',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;figbig;figbeach();
subplot(1,2,1); imagesc(M_x_c_std__); axis image; axisnotick;
title('M_x_c_std__','Interpreter','none');
subplot(1,2,2); imagesc(N_x_c_std__); axis image; axisnotick;
title('N_x_c_std__','Interpreter','none');
sgtitle(sprintf(' %% centering filtered images'),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%}

