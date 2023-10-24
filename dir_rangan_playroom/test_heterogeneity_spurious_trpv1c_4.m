%%%%%%%%;
% tests heterogeneity_spurious on trpv1c. ;
% Adding calculation of bayesian likelihood. ;
%%%%%%%%;

clear;
%%%%%%%%;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
%%%%%%%%;
test_pm_trpv1c_9b;

%%%%%%%%;
% save with tmp9.m. ;
%%%%%%%%;
dir_base = '/data/rangan/dir_cryoem/dir_rangan_playroom';
fname_mat = sprintf('%s/test_heterogeneity_spurious_trpv1c_3_workspace_20231019.mat',dir_base);
load(fname_mat);

%%%%%%%%;
% Now load a larger set of images. ;
%%%%%%%%;
n_big_M = 1024*32;
[ ...
 big_M_x_c___ ...
,index_nbig_CTF_from_nbig_M_ ...
,index_nbig_M_from_nbig_CTF_ ...
,Voltage_big_CTF_ ...
,DefocusU_big_CTF_ ...
,DefocusV_big_CTF_ ...
,DefocusAngle_big_CTF_ ...
,SphericalAberration_big_CTF_ ...
,AmplitudeContrast_big_CTF_ ...
] = ...
rlnImageName_from_star_0( ...
 dir_data_star ...
,fname_nopath_star ...
,n_big_M ...
);
n_x_M_u = size(big_M_x_c___,1);
assert(n_x_M_u==size(big_M_x_c___,2));
disp(sprintf(' %% Removing edge-artefacts'));
n_big_M_ext_ = zeros(n_big_M,1);
for nbig_M=0:n_big_M-1;
if (mod(nbig_M,128)==0); disp(sprintf(' %% nbig_M %d/%d',nbig_M,n_big_M)); end;
n_pixel = 4; edge_tolerance = 0.5; n_edge_overshoot = 8; rseed = 0;
[big_M_x_c___(:,:,1+nbig_M),n_big_M_ext_(1+nbig_M)] = image_replace_edge_artefact_0(big_M_x_c___(:,:,1+nbig_M),4,0.5,2,0);
end;%for nbig_M=0:n_big_M-1;
disp(sprintf(' %% edge-artefacts detected in %d/%d images.',numel(find(n_big_M_ext_>0)),n_big_M));
disp(sprintf(' %% typical edge-artefact covers %0.6f = (%0.6f)^2 of image.',median(n_big_M_ext_/n_x_M_u^2),median(sqrt(n_big_M_ext_/n_x_M_u^2))));
%%%%%%%%;
% Now examine image-centroids. ;
%%%%%%%%;
x_c_0_ = -half_diameter_x_c + transpose([0:n_x_M_u-1]/n_x_M_u)*diameter_x_c;
x_c_1_ = -half_diameter_x_c + transpose([0:n_x_M_u-1]/n_x_M_u)*diameter_x_c;
[x_c_0__,x_c_1__] = ndgrid(x_c_0_,x_c_1_);
x_p_r__ = sqrt(x_c_0__.^2 + x_c_1__.^2);
x_c_mask__ = x_p_r__<=half_diameter_x_c;
big_M_abs_x_c_0_avg_ = zeros(n_big_M,1);
big_M_abs_x_c_1_avg_ = zeros(n_big_M,1);
big_M_mask_abs_x_c_0_avg_ = zeros(n_big_M,1);
big_M_mask_abs_x_c_1_avg_ = zeros(n_big_M,1);
for nbig_M=0:n_big_M-1;
big_M_abs_x_c_ = abs(squeeze(big_M_x_c___(:,:,1+nbig_M))); %<-- no mask. ;
big_M_abs_avg = mean(big_M_abs_x_c_,'all');
big_M_abs_x_c_0_avg = mean(big_M_abs_x_c_/big_M_abs_avg.*x_c_0__,'all');
big_M_abs_x_c_1_avg = mean(big_M_abs_x_c_/big_M_abs_avg.*x_c_1__,'all');
big_M_abs_x_c_0_avg_(1+nbig_M) = big_M_abs_x_c_0_avg;
big_M_abs_x_c_1_avg_(1+nbig_M) = big_M_abs_x_c_1_avg;
clear big_M_abs_x_c_;
big_M_mask_abs_x_c_ = abs(squeeze(big_M_x_c___(:,:,1+nbig_M)).*x_c_mask__); %<-- radial mask. ;
big_M_mask_abs_avg = mean(big_M_mask_abs_x_c_,'all');
big_M_mask_abs_x_c_0_avg = mean(big_M_mask_abs_x_c_/big_M_mask_abs_avg.*x_c_0__,'all');
big_M_mask_abs_x_c_1_avg = mean(big_M_mask_abs_x_c_/big_M_mask_abs_avg.*x_c_1__,'all');
big_M_mask_abs_x_c_0_avg_(1+nbig_M) = big_M_mask_abs_x_c_0_avg;
big_M_mask_abs_x_c_1_avg_(1+nbig_M) = big_M_mask_abs_x_c_1_avg;
clear big_M_mask_abs_x_c_;
end;%for nbig_M=0:n_big_M-1;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/big_M_abs_x_c_avg_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figmed;
%%%%;
subplot(1,2,1);
plot(big_M_abs_x_c_0_avg_,big_M_abs_x_c_1_avg_,'.');
xlabel('big_M_abs_x_c_0_avg_','Interpreter','none');
ylabel('big_M_abs_x_c_1_avg_','Interpreter','none');
axis equal; grid on;
title('big_M_abs_x_c_ (no radial mask)','Interpreter','none');
%%%%;
subplot(1,2,2);
plot(big_M_mask_abs_x_c_0_avg_,big_M_mask_abs_x_c_1_avg_,'.');
xlabel('big_M_mask_abs_x_c_0_avg_','Interpreter','none');
ylabel('big_M_mask_abs_x_c_1_avg_','Interpreter','none');
axis equal; grid on;
title('big_M_mask_abs_x_c_ (max radial mask)','Interpreter','none');
%%%%;
tmp_corr_ = corr([big_M_abs_x_c_0_avg_,big_M_abs_x_c_1_avg_],[big_M_mask_abs_x_c_0_avg_,big_M_mask_abs_x_c_1_avg_]);
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
big_M_k_p__ = zeros(n_w_sum,n_big_M);
image_center_delta_x_c_0_ = zeros(n_big_M,1);
image_center_delta_x_c_1_ = zeros(n_big_M,1);
for nbig_M=0:n_big_M-1;
if (mod(nbig_M,128)==0); disp(sprintf(' %% nbig_M %d/%d',nbig_M,n_big_M)); end;
big_M_x_c_ = squeeze(big_M_x_c___(:,:,1+nbig_M));
big_M_k_p_ = interp_x_c_to_k_p_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,big_M_x_c_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_M_u^2)*dx^2 ;
% Now *DO* translate according to big_M_abs_x_c_0_avg_ and big_M_abs_x_c_1_avg_. ;
big_M_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,big_M_k_p_,-1*big_M_abs_x_c_0_avg_(1+nbig_M),-1*big_M_abs_x_c_1_avg_(1+nbig_M));
big_delta_x_c_0 = 0; big_delta_x_c_1 = 0;
if flag_center;
% if flag_center, then *ALSO* center images after filtering/masking/thresholding: ;
tmp_parameter = struct('type','parameter');
[ ...
 tmp_parameter ...
,big_M_k_p_ ...
,~ ...
,big_delta_x_c_0 ...
,big_delta_x_c_1 ...
] = ...
image_center_0( ...
 tmp_parameter ...
,max(n_w_max,n_x_M_u/4) ...
,diameter_x_c ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,big_M_k_p_ ...
,weight_2d_k_all_ ...
);
end;%if flag_center;
big_M_k_p__(:,1+nbig_M) = big_M_k_p_;
big_image_center_delta_x_c_0_(1+nbig_M) = big_delta_x_c_0;
big_image_center_delta_x_c_1_(1+nbig_M) = big_delta_x_c_1;
end;%for nbig_M=0:n_big_M-1;
%%%%%%%%;
% Now convert images to big_M_k_q__. ;
%%%%%%%%;
big_M_k_q__ = zeros(n_w_sum,n_big_M);
for nbig_M=0:n_big_M-1;
big_M_k_q__(:,1+nbig_M) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,big_M_k_p__(:,1+nbig_M));
end;%for nbig_M=0:n_big_M-1;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/big_M_k_q__',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;
big_Mlim_ = std(abs(big_M_k_p__(:)))*2.5*[-1,+1];
for nl=0:15-1;
subplot(3,5,1+nl); nbig_M = max(0,min(n_big_M-1,floor(n_big_M*nl/(15-1))));
%imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,abs(big_M_k_q__(:,1+nbig_M)),big_Mlim_,colormap_80s);
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,log10(abs(big_M_k_q__(:,1+nbig_M))/max(abs(big_M_k_q__(:)))),[-4,0],colormap_80s);
title(sprintf('nbig_M %d',nbig_M));
end;%for nl=0:15-1;
sgtitle(sprintf('big_M_k_q__'),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/big_M_k_q_norm_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;
subplot(2,2,1);
S_k_q_norm_ = sqrt(sum(abs(S_k_q__).^2,2));
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,log10(S_k_q_norm_/max(S_k_q_norm_)),[-4,0],colormap_80s);
title('S_k_q_norm_','Interpreter','none');
subplot(2,2,2);
big_M_k_q_norm_ = sqrt(sum(abs(big_M_k_q__).^2,2));
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,log10(big_M_k_q_norm_/max(big_M_k_q_norm_)),[-4,0],colormap_80s);
title('big_M_k_q_norm_','Interpreter','none');
subplot(2,2,3);
S_k_q_norm_ = sqrt(sum(abs(S_k_q__).^2,2)); tmp_eps = S_k_q_norm_(1+n_w_csum_(1+n_k_p_r-1)+n_w_(1+n_k_p_r-1)/2);
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,log10(max(0,S_k_q_norm_-tmp_eps)/max(S_k_q_norm_-tmp_eps)),[-4,0],colormap_80s);
title('S_k_q_norm_ - tmp_eps','Interpreter','none');
subplot(2,2,4);
big_M_k_q_norm_ = sqrt(sum(abs(big_M_k_q__).^2,2)); tmp_eps = big_M_k_q_norm_(1+n_w_csum_(1+n_k_p_r-1)+n_w_(1+n_k_p_r-1)/2);
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,log10(max(0,big_M_k_q_norm_-tmp_eps)/max(big_M_k_q_norm_-tmp_eps)),[-4,0],colormap_80s);
title('big_M_k_q_norm_ - tmp_eps','Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

%%%%%%%%;
% Now calculate big_CTF functions. ;
%%%%%%%%;
n_big_CTF = numel(index_nbig_M_from_nbig_CTF_);
big_CTF_index_ = index_nbig_CTF_from_nbig_M_;
big_CTF_k_p__ = zeros(n_w_sum,n_big_CTF);
for nbig_CTF=0:n_big_CTF-1;
if (mod(nbig_CTF,100)==0); disp(sprintf(' %% nbig_CTF %d/%d',nbig_CTF,n_big_CTF)); end;
big_CTF_Spherical_Aberration = SphericalAberration_big_CTF_(1+nbig_CTF);% spherical aberation of the lens in mm ;
big_CTF_Spherical_Aberration=big_CTF_Spherical_Aberration*(10.0d0^7.0d0);% convert into Angstroms ;
big_CTF_Voltage_kV = Voltage_big_CTF_(1+nbig_CTF);% voltage in kVolts ;
big_CTF_Voltage_1V=big_CTF_Voltage_kV*1000.0 ;% convert into Volts ;
big_CTF_lambda = 12.2643247/sqrt(big_CTF_Voltage_1V+big_CTF_Voltage_1V^2*0.978466d-6);% electron wavelength in Angstroms ;
big_CTF_Defocus_U = DefocusU_big_CTF_(1+nbig_CTF);% defocus values (in Angstroms) ;
big_CTF_Defocus_V = DefocusV_big_CTF_(1+nbig_CTF);% defocus values (in Angstroms) ;
big_CTF_Defocus_Angle = DefocusAngle_big_CTF_(1+nbig_CTF);% angle of astigmatism ;
%big_CTF_Defocus_Angle = big_CTF_Defocus_Angle*pi/180.0d0;% convert into radians ; %<-- already in radians! make sure not to convert twice!;
big_CTF_Amplitude_Contrast = AmplitudeContrast_big_CTF_(1+nbig_CTF);% big_CTF_Amplitude Contrast ;
tmp_w1=sqrt(1.0d0-big_CTF_Amplitude_Contrast^2);% weights for the amplitude and phase contrasts in big_CTF ;
tmp_w2=big_CTF_Amplitude_Contrast;% weights for the amplitude and phase contrasts in big_CTF ;
%  big_CTF_Object_Pixel_Size = big_CTF_Detector_Pixel_Size/big_CTF_Magnification;
big_CTF_Object_Pixel_Size = Pixel_Spacing;% pixel size of the scanner in physical space in Angstroms ;
big_CTF_lambda_per_box = big_CTF_lambda/(n_x_M_u*big_CTF_Object_Pixel_Size);% n_x_M_u*big_CTF_Object_Pixel_Size is the box size in Angstroms ;
%%%%;
na=0;
for nk = 0:n_k_p_r-1;
for nw=0:n_w_(1+nk)-1;
tmp_theta = (2.0d0*pi*nw)/n_w_(1+nk);
tmp_k_c_1 = (2.0d0*pi)*k_p_r_(1+nk)*cos(tmp_theta);
tmp_k_c_2 = (2.0d0*pi)*k_p_r_(1+nk)*sin(tmp_theta);
tmp_ctf_value = niko_ctf(big_CTF_Spherical_Aberration,big_CTF_lambda,tmp_w1,tmp_w2,big_CTF_Defocus_U,big_CTF_Defocus_V,big_CTF_Defocus_Angle,big_CTF_lambda_per_box/pi,tmp_k_c_1,tmp_k_c_2);
big_CTF_k_p__(1+na,1+nbig_CTF) = -tmp_ctf_value;
na = na+1;
end;%for nw=0:n_w_(1+nk)-1;
end;%for nk = 0:n_k_p_r-1;
%%%%;
end;%for nbig_CTF=0:n_big_CTF-1;
%%%%%%%%;
big_CTF_k_p_r__ = zeros(n_k_p_r,n_big_CTF);
for nbig_CTF=0:n_big_CTF-1;
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
big_CTF_k_p_r__(1+nk_p_r,1+nbig_CTF) = mean(big_CTF_k_p__(1+tmp_index_,1+nbig_CTF));
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nbig_CTF=0:n_big_CTF-1;
big_CTF_avg_k_p_ = mean(big_CTF_k_p__,2);
%imagesc_p(n_k_p_r,k_p_r_,n_w_,sum(n_w_),real(big_CTF_avg_k_p_(:)),[-1,+1],colormap_beach());
big_CTF_avg_k_p_r_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
big_CTF_avg_k_p_r_(1+nk_p_r) = mean(big_CTF_avg_k_p_(1+tmp_index_));
end;%for nk_p_r=0:n_k_p_r-1;
Sbig_CTF_ = svd(big_CTF_k_p_r__(:,1+big_CTF_index_(1:n_big_M)));
n_big_CTF_rank = max(find(Sbig_CTF_/max(Sbig_CTF_)>tolerance_master));
[Ubig_CTF_kc__,Sbig_CTF_c__,Vbig_CTF_Mc__] = svds(big_CTF_k_p_r__(:,1+big_CTF_index_(1:n_big_M)),n_big_CTF_rank);
VSbig_CTF_Mc__ = Vbig_CTF_Mc__*Sbig_CTF_c__;
%%%%%%%%;
% Now plot out some of the big_CTF-functions for varying anisotropy. ;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/big_CTF_k_p_sample__',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;
[tmp_anisotropy_,big_CTF_anisotropy_index_] = sort(DefocusU_big_CTF_ - DefocusV_big_CTF_,'ascend'); big_CTF_anisotropy_index_ = big_CTF_anisotropy_index_ - 1;
for nl=0:15-1;
subplot(3,5,1+nl);
tmp_nbig_CTF = max(0,min(n_big_CTF-1,floor(n_big_CTF*nl/(15-1)))); nbig_CTF = big_CTF_anisotropy_index_(1+tmp_nbig_CTF);
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,big_CTF_k_p__(:,1+nbig_CTF),[-1,+1],colormap_beach());
axis image; axisnotick;
title(sprintf('nbig_CTF %d anisotropy %0.2f',nbig_CTF,tmp_anisotropy_(1+tmp_nbig_CTF)));
end;%for nl=0:15-1;
sgtitle(sprintf('big_CTF_k_p__'),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
% Now determine the big_CTF cross correlation. ;
% This depends  on big_CTF_index_. ;
%%%%%%%%;
tmp_big_CTF_avg_k_p_ = mean(big_CTF_k_p__(:,1+big_CTF_index_(1+(0:n_big_M-1))),2);
tmp_big_CTF_avg_k_p_r_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
tmp_big_CTF_avg_k_p_r_(1+nk_p_r) = mean(tmp_big_CTF_avg_k_p_(1+tmp_index_));
end;%for nk_p_r=0:n_k_p_r-1;
big_CTF_k_p_r_xavg__ = tmp_big_CTF_avg_k_p_r_ * transpose(tmp_big_CTF_avg_k_p_r_);
big_CTF_k_p_r_xcor__ = big_CTF_k_p_r__(:,1+big_CTF_index_(1+(0:n_big_M-1))) * transpose(big_CTF_k_p_r__(:,1+big_CTF_index_(1+(0:n_big_M-1)))) / n_big_M;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/big_CTF_k_p_xcor__',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figmed;
figbeach();
subplot(1,2,1); imagesc(big_CTF_k_p_r_xavg__); axis image; axisnotick; title('big_CTF_k_p_r_xavg__','Interpreter','none');
subplot(1,2,2); imagesc(big_CTF_k_p_r_xcor__); axis image; axisnotick; title('big_CTF_k_p_r_xcor__','Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/big_CTF_k_p_r_xxxx__',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;
clf;
colormap(colormap_beach());
subplot(1,2,1);
imagesc(big_CTF_k_p_r_xavg__,[0,1]);
axis image;
set(gca,'XTick',1:n_k_p_r,'XTickLabel',k_p_r_); xtickangle(90);
set(gca,'YTick',1:n_k_p_r,'YTickLabel',k_p_r_);
title('big_CTF_k_p_r_xavg__ cross correlation','Interpreter','none');
xlabel('k_p_r','Interpreter','none');
ylabel('k_p_r','Interpreter','none');
subplot(1,2,2);
imagesc(big_CTF_k_p_r_xcor__,[0,1]);
axis image;
set(gca,'XTick',1:n_k_p_r,'XTickLabel',k_p_r_); xtickangle(90);
set(gca,'YTick',1:n_k_p_r,'YTickLabel',k_p_r_);
title('big_CTF_k_p_r_xcor__ cross correlation','Interpreter','none');
xlabel('k_p_r','Interpreter','none');
ylabel('k_p_r','Interpreter','none');
figbig;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

%%%%%%%%;
% Now calculate the bayesian likelihood. ;
% Note that this is limited to FTK_0. ;
%%%%%%%%;
%%%%%%%%%%%%%%%%;
%%%%%%%%;
% First estimate distances between templates and images. ;
%%%%%%%%;
ncluster = 0; %<-- pick one of the CTF clusters. ;
UX_kn__ = UX_knc___(:,:,1+ncluster);
n_UX_rank = n_k_p_r-1; %<-- just to check dimensions. ;
CTF_k_p_r_xavg_k_ = CTF_k_p_r_xavg_kc__(:,1+ncluster);
T_k_p_wkS__ = reshape(bsxfun(@times,reshape(S_k_p_wkS__,[n_w_max,n_k_p_r,n_S]),reshape(CTF_k_p_r_xavg_k_,[1,n_k_p_r,1])),[n_w_sum,n_S]);
T_k_q_wkS__ = reshape(bsxfun(@times,reshape(S_k_q_wkS__,[n_w_max,n_k_p_r,n_S]),reshape(CTF_k_p_r_xavg_k_,[1,n_k_p_r,1])),[n_w_sum,n_S]);
T_half_k_p_wkS__ = reshape(bsxfun(@times,reshape(S_half_k_p_wkS__,[n_w_max,n_k_p_r,n_S]),reshape(CTF_k_p_r_xavg_k_,[1,n_k_p_r,1])),[n_w_sum,n_S]);
T_half_k_q_wkS__ = reshape(bsxfun(@times,reshape(S_half_k_q_wkS__,[n_w_max,n_k_p_r,n_S]),reshape(CTF_k_p_r_xavg_k_,[1,n_k_p_r,1])),[n_w_sum,n_S]);
T_rem2_k_p_wkS__ = reshape(bsxfun(@times,reshape(S_rem2_k_p_wkS__,[n_w_max,n_k_p_r,n_S]),reshape(CTF_k_p_r_xavg_k_,[1,n_k_p_r,1])),[n_w_sum,n_S]);
T_rem2_k_q_wkS__ = reshape(bsxfun(@times,reshape(S_rem2_k_q_wkS__,[n_w_max,n_k_p_r,n_S]),reshape(CTF_k_p_r_xavg_k_,[1,n_k_p_r,1])),[n_w_sum,n_S]);
tmp_t = tic();
svd_0_VUXT_lwnS____ = tpmh_VUXM_lwnM____3(FTK_0,n_k_p_r,n_w_,n_S,T_k_q_wkS__,n_UX_rank,UX_kn__,sqrt(weight_2d_k_p_r_));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_0_VUXT_lwnS____: %0.3fs',tmp_t)); end;
tmp_t = tic();
svd_0_VUXT_half_lwnS____ = tpmh_VUXM_lwnM____3(FTK_0,n_k_p_r,n_w_,n_S,T_half_k_q_wkS__,n_UX_rank,UX_kn__,sqrt(weight_2d_k_p_r_));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_0_VUXT_half_lwnS____: %0.3fs',tmp_t)); end;
tmp_t = tic();
svd_0_VUXT_rem2_lwnS____ = tpmh_VUXM_lwnM____3(FTK_0,n_k_p_r,n_w_,n_S,T_rem2_k_q_wkS__,n_UX_rank,UX_kn__,sqrt(weight_2d_k_p_r_));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_0_VUXT_rem2_lwnS____: %0.3fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
[UX_0_T_k_q_wnS___,UX_0_T_k_p_wnS___] = ampmh_UX_M_k_p_wnM___0(FTK_0,n_w_,n_UX_rank,n_S,svd_0_VUXT_lwnS____,zeros(n_S,1),zeros(n_S,1));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_0_T_k_q_wnS___: %0.6fs',tmp_t)); end;
UX_0_T_k_p_wnS__ = reshape(UX_0_T_k_p_wnS___(:,1:n_UX_rank,:),[n_w_max*n_UX_rank,n_S]);
UX_0_T_k_q_wnS__ = reshape(UX_0_T_k_q_wnS___(:,1:n_UX_rank,:),[n_w_max*n_UX_rank,n_S]);
%%%%%%%%;
tmp_t = tic();
[UX_0_T_half_k_q_wnS___,UX_0_T_half_k_p_wnS___] = ampmh_UX_M_k_p_wnM___0(FTK_0,n_w_,n_UX_rank,n_S,svd_0_VUXT_half_lwnS____,zeros(n_S,1),zeros(n_S,1));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_0_T_half_k_q_wnS___: %0.6fs',tmp_t)); end;
UX_0_T_half_k_p_wnS__ = reshape(UX_0_T_half_k_p_wnS___(:,1:n_UX_rank,:),[n_w_max*n_UX_rank,n_S]);
UX_0_T_half_k_q_wnS__ = reshape(UX_0_T_half_k_q_wnS___(:,1:n_UX_rank,:),[n_w_max*n_UX_rank,n_S]);
%%%%%%%%;
tmp_t = tic();
[UX_0_T_rem2_k_q_wnS___,UX_0_T_rem2_k_p_wnS___] = ampmh_UX_M_k_p_wnM___0(FTK_0,n_w_,n_UX_rank,n_S,svd_0_VUXT_rem2_lwnS____,zeros(n_S,1),zeros(n_S,1));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_0_T_rem2_k_q_wnS___: %0.6fs',tmp_t)); end;
UX_0_T_rem2_k_p_wnS__ = reshape(UX_0_T_rem2_k_p_wnS___(:,1:n_UX_rank,:),[n_w_max*n_UX_rank,n_S]);
UX_0_T_rem2_k_q_wnS__ = reshape(UX_0_T_rem2_k_q_wnS___(:,1:n_UX_rank,:),[n_w_max*n_UX_rank,n_S]);
%%%%%%%%;
% Now calculate norms of the templates. ;
%%%%%%%%;
tmp_TT_S_ = zeros(n_S,1);
tmp_half_TT_S_ = zeros(n_S,1);
tmp_rem2_TT_S_ = zeros(n_S,1);
for nS=0:n_S-1;
tmp_TT = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,T_k_p_wkS__(:,1+nS),T_k_p_wkS__(:,1+nS))/(2*pi);
tmp_TT_S_(1+nS) = tmp_TT;
tmp_half_TT = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,T_half_k_p_wkS__(:,1+nS),T_half_k_p_wkS__(:,1+nS))/(2*pi);
tmp_half_TT_S_(1+nS) = tmp_half_TT;
tmp_rem2_TT = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,T_rem2_k_p_wkS__(:,1+nS),T_rem2_k_p_wkS__(:,1+nS))/(2*pi);
tmp_rem2_TT_S_(1+nS) = tmp_rem2_TT;
end;%for nS=0:n_S-1;
UX_0_T_l2_S_ = tmp_TT_S_;
UX_0_T_half_l2_S_ = tmp_half_TT_S_;
UX_0_T_rem2_l2_S_ = tmp_rem2_TT_S_;
%%%%%%%%;
% Now generate big_svd_0_VUXM_lwnM____. ;
%%%%%%%%;
tmp_t = tic();
big_svd_0_VUXM_lwnM____ = tpmh_VUXM_lwnM____3(FTK_0,n_k_p_r,n_w_,n_big_M,big_M_k_q__,n_UX_rank,UX_kn__,sqrt(weight_2d_k_p_r_));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_0_VUXM_lwnM____: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Now calculate norms of the images. ;
%%%%%%%%;
tmp_t = tic();
UX_0_big_M_l2_dM__ = ampmh_UX_M_l2_dM__1(FTK_0,n_w_,n_big_M,n_UX_rank,big_svd_0_VUXM_lwnM____);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_0_big_M_l2_dM__: %0.3fs',tmp_t)); end;
if (verbose); disp(sprintf(' %% average l2-norm of images: %0.16f',mean(UX_0_big_M_l2_dM__(:))/(pi*k_p_r_max^2))); end;
tmp_big_MM_M_ = zeros(n_big_M,1);
for nbig_M=0:n_big_M-1;
tmp_big_MM = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,big_M_k_p__(:,1+nbig_M),big_M_k_p__(:,1+nbig_M))/(2*pi);
tmp_big_MM_M_(1+nbig_M) = tmp_big_MM;
end;%for nbig_M=0:n_big_M-1;
tmp_index = efind((FTK_0.delta_x_.^2 + FTK_0.delta_y_.^2)<1e-6);
UX_0_big_M_l2_M_ = transpose(UX_0_big_M_l2_dM__(1+tmp_index,:));
if (verbose); disp(sprintf(' %% tmp_big_MM_M_ vs UX_0_big_M_l2_M_: %0.16f',fnorm(tmp_big_MM_M_ - UX_0_big_M_l2_M_)/fnorm(tmp_big_MM_M_))); end;
flag_plot=0;
if flag_plot;
tmp_index = efind((FTK_0.delta_x_.^2 + FTK_0.delta_y_.^2)<1e-6);
subplot(1,2,1); hold on; 
plot(0:n_big_M-1,UX_0_big_M_l2_M_/(pi*k_p_r_max^2),'rx'); xlabel('nbig_M'); ylabel('l2');
plot(0:n_big_M-1,tmp_big_MM_M_/(pi*k_p_r_max^2),'bo'); xlabel('nbig_M'); ylabel('l2');
hold off;
subplot(1,2,2); plot(UX_0_big_M_l2_M_/(pi*k_p_r_max^2),tmp_big_MM_M_/(pi*k_p_r_max^2),'g.'); xlabel('l2_A'); ylabel('l2_B');
end;%if flag_plot;
%%%%%%%%;
tmp_t = tic();
[UX_0_big_M_k_q_wnM___,UX_0_big_M_k_p_wnM___] = ampmh_UX_M_k_p_wnM___0(FTK_0,n_w_,n_UX_rank,n_big_M,big_svd_0_VUXM_lwnM____,zeros(n_big_M,1),zeros(n_big_M,1));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_0_big_M_k_q_wnM___: %0.6fs',tmp_t)); end;
UX_0_big_M_k_p_wnM__ = reshape(UX_0_big_M_k_p_wnM___(:,1:n_UX_rank,:),[n_w_max*n_UX_rank,n_big_M]);
UX_0_big_M_k_q_wnM__ = reshape(UX_0_big_M_k_q_wnM___(:,1:n_UX_rank,:),[n_w_max*n_UX_rank,n_big_M]);
%%%%%%%%;
% Visualize: ;
%%%%%%%%;
flag_plot=0;
if flag_plot;
figure(3);clf;figbig;fig80s;
subplot(2,2,1);imagesc(reshape(permute(log10(abs(big_svd_0_VUXM_lwnM____)),[1,3,4,2]),[FTK_0.n_svd_l*n_UX_rank*n_big_M,n_w_max]));axisnotick; colorbar;
subplot(2,2,2);plot(0:n_w_max-1,log10(sqrt(mean(reshape(permute(abs(big_svd_0_VUXM_lwnM____).^2,[1,3,4,2]),[FTK_0.n_svd_l*n_UX_rank*n_big_M,n_w_max]),1))));
subplot(2,2,3);imagesc(reshape(permute(reshape(log10(abs(UX_0_big_M_k_q_wnM__)),[n_w_max,n_UX_rank,n_big_M]),[2,3,1]),[n_UX_rank*n_big_M,n_w_max]));axisnotick;colorbar;
subplot(2,2,4);plot(0:n_w_max-1,log10(sqrt(mean(reshape(permute(reshape(abs(UX_0_big_M_k_q_wnM__).^2,[n_w_max,n_UX_rank,n_big_M]),[2,3,1]),[n_UX_rank*n_big_M,n_w_max]),1))));
end;%if flag_plot;

%%%%%%%%;
% Calculate ampmh_X_wSM___8. ;
%%%%%%%%;
parameter = struct('type','parameter');
parameter.flag_optimize_over_gamma_z = 1;
parameter.flag_compute_I_value = 0;
parameter.tolerance_master = tolerance_master;
parameter.pm_n_UX_rank_use = n_UX_rank;
tmp_t = tic();
[ ...
 parameter ...
,big_X_0_TM__ ...
] = ...
ampmh_X_wSM___8( ...
 parameter ...
,FTK_0 ...
,n_w_max ...
,n_UX_rank ...
,n_S ...
,UX_0_T_k_q_wnS__ ...
,UX_0_T_l2_S_ ...
,n_big_M ...
,big_svd_0_VUXM_lwnM____ ...
,UX_0_big_M_l2_dM__ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% big_X_0_TM__: %0.3fs',tmp_t)); end;
%%%%%%%%;
parameter = struct('type','parameter');
parameter.flag_optimize_over_gamma_z = 1;
parameter.flag_compute_I_value = 0;
parameter.tolerance_master = tolerance_master;
parameter.pm_n_UX_rank_use = n_UX_rank;
tmp_t = tic();
[ ...
 parameter ...
,big_X_0_half_TM__ ...
] = ...
ampmh_X_wSM___8( ...
 parameter ...
,FTK_0 ...
,n_w_max ...
,n_UX_rank ...
,n_S ...
,UX_0_T_half_k_q_wnS__ ...
,UX_0_T_half_l2_S_ ...
,n_big_M ...
,big_svd_0_VUXM_lwnM____ ...
,UX_0_big_M_l2_dM__ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% big_X_0_half_TM__: %0.3fs',tmp_t)); end;
%%%%%%%%;
parameter = struct('type','parameter');
parameter.flag_optimize_over_gamma_z = 1;
parameter.flag_compute_I_value = 0;
parameter.tolerance_master = tolerance_master;
parameter.pm_n_UX_rank_use = n_UX_rank;
tmp_t = tic();
[ ...
 parameter ...
,big_X_0_rem2_TM__ ...
] = ...
ampmh_X_wSM___8( ...
 parameter ...
,FTK_0 ...
,n_w_max ...
,n_UX_rank ...
,n_S ...
,UX_0_T_rem2_k_q_wnS__ ...
,UX_0_T_rem2_l2_S_ ...
,n_big_M ...
,big_svd_0_VUXM_lwnM____ ...
,UX_0_big_M_l2_dM__ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% big_X_0_rem2_TM__: %0.3fs',tmp_t)); end;
%%%%%%%%;
%%%%%%%%%%%%%%%%;
% imagesc(hist2d_0(reshape(big_X_0_TM__,[n_S*n_big_M,1]),reshape(max(big_X_0_half_TM__,big_X_0_rem2_TM__),[n_S*n_big_M,1]),129,129));set(gca,'Ydir','normal'); %<-- histogram. ;
%%%%%%%%;
flag_disp=0;
if flag_disp;
%%%%;
% visualize distance. ;
%%%%;
big_D2_0_TM__ = bsxfun(@plus,reshape(UX_0_T_l2_S_,[n_S,1]),reshape(UX_0_big_M_l2_M_,[1,n_big_M])) ...
  - 2*big_X_0_TM__ .* (reshape(sqrt(UX_0_T_l2_S_),[n_S,1])*reshape(sqrt(UX_0_big_M_l2_M_),[1,n_big_M]));
big_D2_0_half_TM__ = bsxfun(@plus,reshape(UX_0_T_half_l2_S_,[n_S,1]),reshape(UX_0_big_M_l2_M_,[1,n_big_M])) ...
  - 2*big_X_0_half_TM__ .* (reshape(sqrt(UX_0_T_half_l2_S_),[n_S,1])*reshape(sqrt(UX_0_big_M_l2_M_),[1,n_big_M]));
big_D2_0_rem2_TM__ = bsxfun(@plus,reshape(UX_0_T_rem2_l2_S_,[n_S,1]),reshape(UX_0_big_M_l2_M_,[1,n_big_M])) ...
  - 2*big_X_0_rem2_TM__ .* (reshape(sqrt(UX_0_T_rem2_l2_S_),[n_S,1])*reshape(sqrt(UX_0_big_M_l2_M_),[1,n_big_M]));
hlim_ = [0.050,0.125];
tmp_h2__ = hist2d_0(reshape(big_D2_0_TM__,[n_S*n_big_M,1]),reshape(min(big_D2_0_half_TM__,big_D2_0_rem2_TM__),[n_S*n_big_M,1]),129,129,hlim_,hlim_);
figure(1+nf);nf=nf+1;clf;figsml;fig81s;
imagesc(tmp_h2__); xlabel('big_D2_0_TM__','Interpreter','none'); ylabel('min(big_D2_0_half_TM__,big_D2_0_rem2_TM__)','Interpreter','none');
axis image; axisnotick;
set(gca,'Ydir','normal');
end%;if flag_disp;

%%%%%%%%;
% Now combine these correlations into a likelihood. ;
%%%%%%%%;
viewing_weight_emp_ = viewing_weight_all_;
u_viewing_polar_a_ = sort(unique(viewing_polar_a_all_),'ascend');
assert(fnorm(u_viewing_polar_a_ - sort(viewing_polar_a_,'ascend'))<1e-12); %<-- viewing_polar_a_ defined in test_pm_trpv1c_9b.m ;
tmp_h_polar_a_all_ = hist(viewing_polar_a_all_,u_viewing_polar_a_);
tmp_h_polar_a_emp_ = hist(euler_polar_a_true_,u_viewing_polar_a_);
tmp_h_polar_a_rel_ = tmp_h_polar_a_emp_./max(1,tmp_h_polar_a_all_);
tmp_ij_ = knnsearch(u_viewing_polar_a_,viewing_polar_a_all_);
viewing_weight_emp_ = viewing_weight_emp_.*reshape(tmp_h_polar_a_rel_(tmp_ij_),[n_S,1]);
viewing_weight_emp_ = viewing_weight_emp_.*sum(viewing_weight_all_)./max(1e-12,sum(viewing_weight_emp_));
%%%%;
% visualize: ;
% subplot(1,2,1); imagesc_polar_a_azimu_b_0(viewing_polar_a_all_,viewing_azimu_b_all_,viewing_weight_all_,[0,1]); title('uni');
% subplot(1,2,2); imagesc_polar_a_azimu_b_0(viewing_polar_a_all_,viewing_azimu_b_all_,viewing_weight_emp_,[0,1]); title('emp');
%%%%;
sigma_bayesian_ = transpose([0,2.^[-12:0.125:0]]); n_sigma_bayesian = numel(sigma_bayesian_);
%%%%;
[ ...
 ~ ...
,big_ssnll_0_uni_Ms__ ...
,big_D2_0_uni_SM__ ...
] = ...
ssnll_from_X_0( ...
 [] ...
,n_S ...
,n_big_M ...
,big_X_0_TM__ ...
,UX_0_T_l2_S_ ...
,UX_0_big_M_l2_M_ ...
,viewing_weight_all_ ...
,n_sigma_bayesian ...
,sigma_bayesian_ ...
);
%%%%;
[ ...
 ~ ...
,big_ssnll_0_uni_half_Ms__ ...
,big_D2_0_uni_half_SM__ ...
] = ...
ssnll_from_X_0( ...
 [] ...
,n_S ...
,n_big_M ...
,big_X_0_half_TM__ ...
,UX_0_T_half_l2_S_ ...
,UX_0_big_M_l2_M_ ...
,viewing_weight_all_ ...
,n_sigma_bayesian ...
,sigma_bayesian_ ...
);
%%%%;
[ ...
 ~ ...
,big_ssnll_0_uni_rem2_Ms__ ...
,big_D2_0_uni_rem2_SM__ ...
] = ...
ssnll_from_X_0( ...
 [] ...
,n_S ...
,n_big_M ...
,big_X_0_rem2_TM__ ...
,UX_0_T_rem2_l2_S_ ...
,UX_0_big_M_l2_M_ ...
,viewing_weight_all_ ...
,n_sigma_bayesian ...
,sigma_bayesian_ ...
);
%%%%%%%%;
big_ssnllr_0_uni_Ms__ = big_ssnll_0_uni_Ms__ - min(big_ssnll_0_uni_half_Ms__,big_ssnll_0_uni_rem2_Ms__); %<-- assumes hard assignment of images to volumes in 2-volume-model (see min). ;
big_ssnllr_0_uni_s_ = sum(big_ssnllr_0_uni_Ms__,1); %<-- assumes hard assignment of images to volumes in 2-volume-model (see min). ;
%%%%%%%%;
[ ...
 ~ ...
,big_ssnll_0_emp_Ms__ ...
,big_D2_0_emp_SM__ ...
] = ...
ssnll_from_X_0( ...
 [] ...
,n_S ...
,n_big_M ...
,big_X_0_TM__ ...
,UX_0_T_l2_S_ ...
,UX_0_big_M_l2_M_ ...
,viewing_weight_emp_ ...
,n_sigma_bayesian ...
,sigma_bayesian_ ...
);
%%%%;
[ ...
 ~ ...
,big_ssnll_0_emp_half_Ms__ ...
,big_D2_0_emp_half_SM__ ...
] = ...
ssnll_from_X_0( ...
 [] ...
,n_S ...
,n_big_M ...
,big_X_0_half_TM__ ...
,UX_0_T_half_l2_S_ ...
,UX_0_big_M_l2_M_ ...
,viewing_weight_emp_ ...
,n_sigma_bayesian ...
,sigma_bayesian_ ...
);
%%%%;
[ ...
 ~ ...
,big_ssnll_0_emp_rem2_Ms__ ...
,big_D2_0_emp_rem2_SM__ ...
] = ...
ssnll_from_X_0( ...
 [] ...
,n_S ...
,n_big_M ...
,big_X_0_rem2_TM__ ...
,UX_0_T_rem2_l2_S_ ...
,UX_0_big_M_l2_M_ ...
,viewing_weight_emp_ ...
,n_sigma_bayesian ...
,sigma_bayesian_ ...
);
%%%%%%%%;
big_ssnllr_0_emp_Ms__ = big_ssnll_0_emp_Ms__ - min(big_ssnll_0_emp_half_Ms__,big_ssnll_0_emp_rem2_Ms__); %<-- assumes hard assignment of images to volumes in 2-volume-model (see min). ;
big_ssnllr_0_emp_s_ = sum(big_ssnllr_0_emp_Ms__,1); %<-- assumes hard assignment of images to volumes in 2-volume-model (see min). ;
%%%%%%%%;

%%%%%%%%;
% Plot the log-unlikelihood-ratio vs temperature: ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figmed;
linewidth_use = 2;
markersize_use = 12;
fontsize_use = 12;
hold on;
plot(big_ssnllr_0_uni_s_(2:end),'ko-','LineWidth',linewidth_use,'MarkerSize',markersize_use,'MarkerFaceColor','g');
plot(big_ssnllr_0_emp_s_(2:end),'ko-','LineWidth',linewidth_use,'MarkerSize',markersize_use,'MarkerFaceColor','r');
hold off;
xlim([0.5,n_sigma_bayesian-1]);
xlabel('$-\log_{2}$(temperature)','Interpreter','latex');
set(gca,'XTick',1:4:n_sigma_bayesian-1,'XTickLabel',num2str(-log2(sigma_bayesian_(2:4:end)),'%0.2f'));xtickangle(90);
ylabel('$\sigma^{2}\cdot$ nllr','Interpreter','latex');
legend({'uniform','empirical'},'Location','SouthWest');
set(gca,'FontSize',fontsize_use);
%%%%;
fname_fig_pre = sprintf('%s_jpg/test_heterogeneity_spurious_trpv1c_big_ssnllr_0_FIGA_stripped',dir_pm);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_jpg));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
end;%if (flag_replot | ~exist(fname_fig_jpg,'file'));
%%%%;
fname_fig_pre = sprintf('%s_jpg/test_heterogeneity_spurious_trpv1c_big_ssnllr_0_FIGA',dir_pm);
sgtitle(fname_fig_pre,'Interpreter','none');
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_jpg));
print('-djpeg',fname_fig_jpg);
end;%if (flag_replot | ~exist(fname_fig_jpg,'file'));
close(gcf);

%%%%%%%%;
% Plot the optimal correlation, optimal distance and template-magnitudes. ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figmed;
linewidth_use = 2;
markersize_use = 6;
fontsize_use = 12;
p_row = 2; p_col = 3; np=0;
tmp_X_lim_ = prctile(max(big_X_0_TM__,[],1), [ 1,99]); tmp_X_lim_ = mean(tmp_X_lim_) + 0.5*1.75*diff(tmp_X_lim_)*[-1,+1];
tmp_D2_lim_ = prctile(min(big_D2_0_uni_SM__,[],1), [ 1,99]); tmp_D2_lim_ = mean(tmp_D2_lim_) + 0.5*1.75*diff(tmp_D2_lim_)*[-1,+1];
tmp_l2_lim_ = prctile(UX_T_l2_S_, [ 1,99]); tmp_l2_lim_ = mean(tmp_l2_lim_) + 0.5*1.75*diff(tmp_l2_lim_)*[-1,+1];
%%%%;
subplot(p_row,p_col,1+np);np=np+1; cla;
plot(tmp_X_lim_,tmp_X_lim_,'k-',max(big_X_0_TM__,[],1),max(big_X_0_half_TM__,[],1),'r.','MarkerSize',markersize_use); axis square;
xlim(tmp_X_lim_);ylim(tmp_X_lim_);
xlabel('big_X_0_TM__','Interpreter','none'); ylabel('big_X_0_half_TM__','Interpreter','none'); axisnotick; set(gca,'XTick',tmp_X_lim_,'YTick',tmp_X_lim_);
set(gca,'FontSize',fontsize_use);
subplot(p_row,p_col,1+np);np=np+1; cla;
plot(tmp_D2_lim_,tmp_D2_lim_,'k-',min(big_D2_0_uni_SM__,[],1),min(big_D2_0_uni_half_SM__,[],1),'r.','MarkerSize',markersize_use); axis square;
xlim(tmp_D2_lim_);ylim(tmp_D2_lim_);
xlabel('big_D2_0_uni_SM__','Interpreter','none'); ylabel('big_D2_0_uni_half_SM__','Interpreter','none'); axisnotick; set(gca,'XTick',tmp_D2_lim_,'YTick',tmp_D2_lim_);
set(gca,'FontSize',fontsize_use);
subplot(p_row,p_col,1+np);np=np+1; cla;
plot(tmp_l2_lim_,tmp_l2_lim_,'k-',UX_T_l2_S_,UX_T_half_l2_S_,'r.','MarkerSize',markersize_use); axis square;
xlim(tmp_l2_lim_);ylim(tmp_l2_lim_);
xlabel('UX_T_l2_S_','Interpreter','none'); ylabel('UX_T_half_l2_S_','Interpreter','none'); axisnotick; set(gca,'XTick',tmp_l2_lim_,'YTick',tmp_l2_lim_);
set(gca,'FontSize',fontsize_use);
%%%%;
subplot(p_row,p_col,1+np);np=np+1; cla;
plot(tmp_X_lim_,tmp_X_lim_,'k-',max(big_X_0_TM__,[],1),max(big_X_0_rem2_TM__,[],1),'r.','MarkerSize',markersize_use); axis square;
xlabel('big_X_0_TM__','Interpreter','none'); ylabel('big_X_0_rem2_TM__','Interpreter','none'); axisnotick; set(gca,'XTick',tmp_X_lim_,'YTick',tmp_X_lim_);
xlim(tmp_X_lim_);ylim(tmp_X_lim_);
set(gca,'FontSize',fontsize_use);
subplot(p_row,p_col,1+np);np=np+1; cla;
plot(tmp_D2_lim_,tmp_D2_lim_,'k-',min(big_D2_0_uni_SM__,[],1),min(big_D2_0_uni_rem2_SM__,[],1),'r.','MarkerSize',markersize_use); axis square;
xlim(tmp_D2_lim_);ylim(tmp_D2_lim_);
xlabel('big_D2_0_uni_SM__','Interpreter','none'); ylabel('big_D2_0_uni_rem2_SM__','Interpreter','none'); axisnotick; set(gca,'XTick',tmp_D2_lim_,'YTick',tmp_D2_lim_);
set(gca,'FontSize',fontsize_use);
subplot(p_row,p_col,1+np);np=np+1; cla;
plot(tmp_l2_lim_,tmp_l2_lim_,'k-',UX_T_l2_S_,UX_T_rem2_l2_S_,'r.','MarkerSize',markersize_use); axis square;
xlim(tmp_l2_lim_);ylim(tmp_l2_lim_);
xlabel('UX_T_l2_S_','Interpreter','none'); ylabel('UX_T_rem2_l2_S_','Interpreter','none'); axisnotick; set(gca,'XTick',tmp_l2_lim_,'YTick',tmp_l2_lim_);
set(gca,'FontSize',fontsize_use);
%%%%;
fname_fig_pre = sprintf('%s_jpg/test_heterogeneity_spurious_trpv1c_big_half_vs_rem2_0_FIGE',dir_pm);
fname_fig_str = sprintf('%s_stripped',fname_fig_pre);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_str);
fname_fig_eps = sprintf('%s.eps',fname_fig_str);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_jpg));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
end;%if (flag_replot | ~exist(fname_fig_jpg,'file'));
%%%%;
sgtitle(fname_fig_pre,'Interpreter','none');
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_jpg));
print('-djpeg',fname_fig_jpg);
end;%if (flag_replot | ~exist(fname_fig_jpg,'file'));
close(gcf);
%%%%%%%%;
