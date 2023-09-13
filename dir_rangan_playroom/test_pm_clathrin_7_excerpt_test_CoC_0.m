n_x_M_u = n_x_u_pack;
Voltage_CTF = 300; %<-- set to 300kV: see https://www.ebi.ac.uk/empiar/EMPIAR-10998/ ; 
DefocusU_CTF = 5000; %<-- see kexin email: defocus 1 = 5000A ;
DefocusV_CTF = 5000; %<-- see kexin email: defocus 1 = 5000A ;
DefocusAngle_CTF = 18.39; %<-- see kexin email: defocus angle = 18.39A. ;
SphericalAberration_CTF = 2.7; %<-- Spherical aberration is usually set to 2.7mm. ;
AmplitudeContrast_CTF = 0.07; %<-- and amplitude contrast to 0.07. ;
%%%%%%%%;
% determine original CTF. ;
%%%%%%%%;
CTF_Spherical_Aberration = SphericalAberration_CTF;% spherical aberration of the lens in mm ;
CTF_Spherical_Aberration=CTF_Spherical_Aberration*(10.0d0^7.0d0);% convert into Angstroms ;
CTF_Voltage_kV = Voltage_CTF;% voltage in kVolts ;
CTF_Voltage_1V=CTF_Voltage_kV*1000.0 ;% convert into Volts ;
CTF_lambda = 12.2643247/sqrt(CTF_Voltage_1V+CTF_Voltage_1V^2*0.978466d-6);% electron wavelength in Angstroms ;
CTF_Defocus_U = DefocusU_CTF;% defocus values (in Angstroms) ;
CTF_Defocus_V = DefocusV_CTF;% defocus values (in Angstroms) ;
CTF_Defocus_Angle = DefocusAngle_CTF;% angle of astigmatism ;
CTF_Defocus_Angle = CTF_Defocus_Angle*pi/180.0d0;% convert into radians ; %<-- already in radians! make sure not to convert twice!;
CTF_Amplitude_Contrast = AmplitudeContrast_CTF;% CTF_Amplitude Contrast ;
tmp_w1=sqrt(1.0d0-CTF_Amplitude_Contrast^2);% weights for the amplitude and phase contrasts in CTF ;
tmp_w2=CTF_Amplitude_Contrast;% weights for the amplitude and phase contrasts in CTF ;
%  CTF_Object_Pixel_Size = CTF_Detector_Pixel_Size/CTF_Magnification;
CTF_Object_Pixel_Size = Pixel_Spacing;% pixel size of the scanner in physical space in Angstroms ;
CTF_lambda_per_box = CTF_lambda/(n_x_M_u*CTF_Object_Pixel_Size);% n_x_M_u*CTF_Object_Pixel_Size is the box size in Angstroms ;
%%%%;
CTF_k_ori_p_wk_ = zeros(n_w_sum,1);
na=0;
for nk = 0:n_k_p_r-1;
for nw=0:n_w_(1+nk)-1;
tmp_theta = (2.0d0*pi*nw)/max(1,n_w_(1+nk));
tmp_k_c_1 = (2.0d0*pi)*k_p_r_(1+nk)*cos(tmp_theta);
tmp_k_c_2 = (2.0d0*pi)*k_p_r_(1+nk)*sin(tmp_theta);
tmp_ctf_value = ...
niko_ctf( ...
 CTF_Spherical_Aberration ... 
,CTF_lambda ... 
,tmp_w1 ... 
,tmp_w2 ... 
,CTF_Defocus_U ... 
,CTF_Defocus_V ... 
,CTF_Defocus_Angle ... 
,CTF_lambda_per_box/pi ... 
,tmp_k_c_1 ... 
,tmp_k_c_2 ...
);
clear tmp_k_c_1 tmp_k_c_2 ;
CTF_k_ori_p_wk_(1+na) = -tmp_ctf_value;
na = na+1;
end;%for nw=0:n_w_(1+nk)-1;
end;%for nk = 0:n_k_p_r-1;
%%%%%%%%;
% determine refined CTF. ;
%%%%%%%%;
CTF_k_ref_p_wk_ = zeros(n_w_ref_sum,1);
na=0;
for nk_ref = 0:n_k_ref_p_r-1;
for nw_ref=0:n_w_ref_(1+nk_ref)-1;
tmp_theta = (2.0d0*pi*nw_ref)/max(1,n_w_ref_(1+nk_ref));
tmp_k_ref_c_1 = (2.0d0*pi)*k_ref_p_r_(1+nk_ref)*cos(tmp_theta);
tmp_k_ref_c_2 = (2.0d0*pi)*k_ref_p_r_(1+nk_ref)*sin(tmp_theta);
tmp_ctf_value = ...
niko_ctf( ...
 CTF_Spherical_Aberration ... 
,CTF_lambda ... 
,tmp_w1 ... 
,tmp_w2 ... 
,CTF_Defocus_U ... 
,CTF_Defocus_V ... 
,CTF_Defocus_Angle ... 
,CTF_lambda_per_box/pi ... 
,tmp_k_ref_c_1 ... 
,tmp_k_ref_c_2 ...
);
clear tmp_k_ref_c_1 tmp_k_ref_c_2 ;
CTF_k_ref_p_wk_(1+na) = -tmp_ctf_value;
na = na+1;
end;%for nw_ref=0:n_w_ref_(1+nk_ref)-1;
end;%for nk_ref = 0:n_k_ref_p_r-1;
%%%%;
% normalize refined CTF. ;
%%%%;
CTF_k_ori_p_std = sqrt( sum(abs(CTF_k_ori_p_wk_).^2.*weight_2d_wk_)/(pi*k_p_r_max^2)*(4*pi^2) );
CTF_k_ref_p_wk_ = CTF_k_ref_p_wk_/max(1e-12,CTF_k_ori_p_std);
%%%%;
CTF_k_ref_p_r_k_ = zeros(n_k_ref_p_r,1);
for nk_ref_p_r=0:n_k_ref_p_r-1;
tmp_index_ = n_w_ref_csum_(1+nk_ref_p_r) + (0:n_w_ref_(1+nk_ref_p_r)-1);
CTF_k_ref_p_r_k_(1+nk_ref_p_r) = mean(CTF_k_ref_p_wk_(1+tmp_index_));
end;%for nk_ref_p_r=0:n_k_ref_p_r-1;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figsml;
plot(k_p_r_,CTF_k_p_r_k_,'rx-',k_ref_p_r_,CTF_k_ref_p_r_k_,'go-');
xlabel('k'); ylabel('CTF');
%%%%%%%%;

%%%%%%%%;
% simple 2d-spatial grid ;
%%%%%%%%;
x_u_0_ = linspace(-half_diameter_x_c,+half_diameter_x_c,n_x_u_pack);
x_u_1_ = linspace(-half_diameter_x_c,+half_diameter_x_c,n_x_u_pack);
[x_u_0__,x_u_1__] = ndgrid(x_u_0_,x_u_1_); n_xx = n_x_u_pack^2;
tmp_d2 = @(delta_) (x_u_0__-delta_(1+0)).^2 + (x_u_1__-delta_(1+1)).^2 ;
tmp_g2 = @(delta_,sigma) exp(-tmp_d2(delta_)./2./sigma.^2) ./ sigma.^2 ./ (2*pi) ;

nS = min(128,n_S-1); %<-- pick a single nS. ;
ngamma_z = min(floor(n_w_ref_max/12),n_w_ref_max-1); %<-- pick a single nw. ;
gamma_z = (2*pi*ngamma_z)/max(1,n_gamma_z);
%%%%;
% define rotated template. ;
%%%%;
S_x_u_xx__ = S_x_u_xxS___(:,:,1+nS);
%%%%;
% redefine template to be a sum of gaussians. ;
%%%%;
flag_source=0;
if flag_source;
n_source = 32*4; rng(0); sigma_source = (1.0/8.0);
S_x_u_xx__ = zeros(n_x_u_pack,n_x_u_pack);
for nsource=0:n_source-1;
delta_source_ = 0.25*(2*rand(2,1)-1);
S_x_u_xx__ = S_x_u_xx__ + tmp_g2(delta_source_,sigma_source);
end;%for nsource=0:n_source-1;
end;%if flag_source;
%%%%;
figure(1+nf);nf=nf+1;clf;figsml;
imagesc_c(n_x_u_pack,x_u_0_,n_x_u_pack,x_u_1_,S_x_u_xx__,[],colormap_beach);
axis image; axisnotick; title('S_x_u_xx__','Interpreter','none');
%%%%;
n_S_x_0 = n_x_u_pack; n_S_x_1 = n_x_u_pack; dx_0 = dx_pack; dx_1 = dx_pack;
S_k_p_wk_ = interp_x_c_to_k_p_xxnufft(n_S_x_0,diameter_x_c,n_S_x_1,diameter_x_c,S_x_u_xx__,n_k_ref_p_r,k_ref_p_r_,n_w_ref_)*sqrt(n_S_x_0*n_S_x_1)*dx_0*dx_1;
S_k_p_wk_ = rotate_p_to_p_fftw(n_k_ref_p_r,n_w_ref_,n_w_ref_sum,S_k_p_wk_,+gamma_z);
S_x_u_xx__ = real( interp_k_p_to_x_c_xxnufft(n_x_u_pack,diameter_x_c,n_x_u_pack,diameter_x_c,n_k_ref_p_r,k_ref_p_r_,n_w_ref_,S_k_p_wk_.*weight_2d_ref_wk_*(2*pi)^2)*sqrt(n_x_u_pack*n_x_u_pack) * n_w_ref_sum );
%%%%;
% define rotated mask. ;
%%%%;
M_x_u_xx__ = M_x_u_xxS___(:,:,1+nS);
M_k_p_wk_ = interp_x_c_to_k_p_xxnufft(n_S_x_0,diameter_x_c,n_S_x_1,diameter_x_c,M_x_u_xx__,n_k_ref_p_r,k_ref_p_r_,n_w_ref_)*sqrt(n_S_x_0*n_S_x_1)*dx_0*dx_1;
M_k_p_wk_ = rotate_p_to_p_fftw(n_k_ref_p_r,n_w_ref_,n_w_ref_sum,M_k_p_wk_,+gamma_z);
M_x_u_xx__ = real( interp_k_p_to_x_c_xxnufft(n_x_u_pack,diameter_x_c,n_x_u_pack,diameter_x_c,n_k_ref_p_r,k_ref_p_r_,n_w_ref_,M_k_p_wk_.*weight_2d_ref_wk_*(2*pi)^2)*sqrt(n_x_u_pack*n_x_u_pack) * n_w_ref_sum );
%%%%;
CTF_S_k_p_wk_ = zeros(n_w_ref_sum,1);
na=0;
for nk_ref_p_r=0:n_k_ref_p_r-1;
n_w_ref = n_w_ref_(1+nk_ref_p_r);
CTF_k_ref_p_r = CTF_k_ref_p_r_k_(1+nk_ref_p_r);
for nw_ref=0:n_w_ref-1;
CTF_S_k_p_wk_(1+na) = S_k_p_wk_(1+na)*CTF_k_ref_p_r;
na=na+1;
end;%for nw_ref=0:n_w_ref-1;
end;%for nk_ref_p_r=0:n_k_ref_p_r-1;
CTF_S_x_c_xx__ = real( interp_k_p_to_x_c_xxnufft(n_x_u_pack,diameter_x_c,n_x_u_pack,diameter_x_c,n_k_ref_p_r,k_ref_p_r_,n_w_ref_,CTF_S_k_p_wk_.*weight_2d_ref_wk_*(2*pi)^2)*sqrt(n_x_u_pack*n_x_u_pack) * n_w_ref_sum );
CTF_S_B_x_c_xx__ = conv2(S_x_u_xx__,PSF_avg_x_u_xx__,'same')*dx_pack*dx_pack;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figmed;np=0;
subplot(1,3,1+np);np=np+1;
imagesc_c(n_x_u_pack,x_u_0_,n_x_u_pack,x_u_1_,S_x_u_xx__,[],colormap_beach);
axis image; axis square; axisnotick; title('S_x_u_xx__','Interpreter','none');
subplot(1,3,1+np);np=np+1;
imagesc_c(n_x_u_pack,x_u_0_,n_x_u_pack,x_u_1_,CTF_S_x_c_xx__,[],colormap_beach);
axis image; axis square; axisnotick; title('CTF_S_x_c_xx__','Interpreter','none');
subplot(1,3,1+np);np=np+1;
imagesc_c(n_x_u_pack,x_u_0_,n_x_u_pack,x_u_1_,CTF_S_B_x_c_xx__,[],colormap_beach);
axis image; axis square; axisnotick; title('CTF_S_B_x_c_xx__','Interpreter','none');
clear CTF_S_B_x_c_xx__ ;
%%%%%%%%;

%%%%%%%%;
% R1S1_A_xxz___ <-- rotation of < Q*CTF , S > ;
% R1S1_B_xxz___ <-- < Q*CTF , rotation of S > ;
% R1S1_C_xxz___ <-- < Q , rotation of CTF*S > ;
%%%%%%%%;
R1S1_A_xxz___ = circshift(R1S1_xxzS____(:,:,:,1+nS),-ngamma_z,3); %<-- recall rotation by tmp_gamma_z. ;
%%%%;
parameter_innerproduct = struct('type','parameter','nw_stride',2,'flag_conv2_vs_svt',1);
tmp_t = tic();
[ ...
 ~ ...
,R1S1_B_xxz___ ...
,scatter_from_tensor_order ...
,scatter_from_tensor_zswk___ ...
] = ...
innerproduct_O_x_u_vs_S_x_u_2( ...
 parameter_innerproduct ...
,n_R_sub_x_u_pack_(1+0) ...
,n_R_sub_x_u_pack_(1+1) ...
,R_sub_x_u_pack_.^1 ...
,n_x_u_pack ...
,n_x_u_pack ...
,S_x_u_xx__ ...
,n_k_ref_p_r ...
,k_ref_p_r_ ...
,k_ref_p_r_max ...
,n_w_ref_ ...
,[] ...
,weight_2d_ref_wk_ ...
,[] ...
,scatter_from_tensor_order ...
,scatter_from_tensor_zswk___ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% innerproduct_O_x_u_vs_S_x_u_2: %0.2fs',tmp_t)); end;
%%%%;
nx0 = 7; nx1 = 3;
R1S1_D = sum(R_sub_x_u_pack_(1+nx0+[0:n_x_u_pack-1],1+nx1+[0:n_x_u_pack-1]).*fliplr(flipud(S_x_u_xx__)),'all');
R1S1_B = R1S1_B_xxz___(1+nx0+floor(n_x_u_pack/2)-1,1+nx1+floor(n_x_u_pack/2)-1,1+0);
disp(sprintf(' %% R1S1_D vs R1S1_B: %0.16f',fnorm(R1S1_D-R1S1_B)/fnorm(R1S1_D)));
%%%%;
parameter_innerproduct = struct('type','parameter','nw_stride',2,'flag_conv2_vs_svt',1);
tmp_t = tic();
[ ...
 ~ ...
,R1S1_C_xxz___ ...
,scatter_from_tensor_order ...
,scatter_from_tensor_zswk___ ...
] = ...
innerproduct_O_x_u_vs_S_x_u_2( ...
 parameter_innerproduct ...
,n_R_sub_x_u_pack_(1+0) ...
,n_R_sub_x_u_pack_(1+1) ...
,Q_sub_x_u_pack_.^1 ...
,n_x_u_pack ...
,n_x_u_pack ...
,CTF_S_x_c_xx__ ...
,n_k_ref_p_r ...
,k_ref_p_r_ ...
,k_ref_p_r_max ...
,n_w_ref_ ...
,[] ...
,weight_2d_ref_wk_ ...
,[] ...
,scatter_from_tensor_order ...
,scatter_from_tensor_zswk___ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% innerproduct_O_x_u_vs_S_x_u_2: %0.2fs',tmp_t)); end;
%%%%;
n_h = 128; np=0;
hlim_A_ = [ min(R1S1_A_xxz___,[],'all'),max(R1S1_A_xxz___,[],'all') ];
h2_AB__ = hist2d_0(R1S1_A_xxz___,R1S1_B_xxz___,n_h,n_h,hlim_A_,hlim_A_);
h2_AC__ = hist2d_0(R1S1_A_xxz___,R1S1_C_xxz___,n_h,n_h,hlim_A_,hlim_A_);
h2_BC__ = hist2d_0(R1S1_B_xxz___,R1S1_C_xxz___,n_h,n_h,hlim_A_,hlim_A_);
figure(1+nf);nf=nf+1;clf;figmed;np=0;
subplot(1,3,1+np);np=np+1;
imagesc(log2(1+h2_AB__)); axis image; axisnotick; colorbar; set(gca,'ydir','normal');
title('R1S1_A_xxz___ vs R1S1_B_xxz___','Interpreter','none');
subplot(1,3,1+np);np=np+1;
imagesc(log2(1+h2_AC__)); axis image; axisnotick; colorbar; set(gca,'ydir','normal');
title('R1S1_A_xxz___ vs R1S1_C_xxz___','Interpreter','none');
subplot(1,3,1+np);np=np+1;
imagesc(log2(1+h2_BC__)); axis image; axisnotick; colorbar; set(gca,'ydir','normal');
title('R1S1_B_xxz___ vs R1S1_C_xxz___','Interpreter','none');
%%%%;
nx0 = 7; nx1 = 3;
R1S1_D = sum(Q_sub_x_u_pack_(1+nx0+[0:n_x_u_pack-1],1+nx1+[0:n_x_u_pack-1]).*fliplr(flipud(CTF_S_x_c_xx__)),'all');
R1S1_C = R1S1_C_xxz___(1+nx0+floor(n_x_u_pack/2)-1,1+nx1+floor(n_x_u_pack/2)-1,1+0);
disp(sprintf(' %% R1S1_D vs R1S1_C: %0.16f',fnorm(R1S1_D-R1S1_C)/fnorm(R1S1_D)));
%%%%%%%%;

%%%%%%%%;
% SCTFCTFS_xxz___ <-- < rotation of CTF*S , rotation of CTF*S > ;
%%%%%%%%;
%%%%;
CTF_S_flip_x_c_xx__ = fliplr(flipud(CTF_S_x_c_xx__));
zero_xx__ = zeros(n_x_u_pack,n_x_u_pack);
CTF_S_flip_pad_x_c_xx__ = repmat(zero_xx__,[3,3]); 
CTF_S_flip_pad_x_c_xx__(n_x_u_pack+[1:n_x_u_pack],n_x_u_pack+[1:n_x_u_pack]) = CTF_S_flip_x_c_xx__;
parameter_innerproduct = struct('type','parameter','nw_stride',2,'flag_conv2_vs_svt',1);
tmp_t = tic();
[ ...
 ~ ...
,SCTFCTFS_xxz___ ...
] = ...
innerproduct_O_x_u_vs_S_x_u_2( ...
 parameter_innerproduct ...
,n_x_u_pack ...
,n_x_u_pack ...
,CTF_S_flip_x_c_xx__ ...
,n_x_u_pack ...
,n_x_u_pack ...
,CTF_S_x_c_xx__ ...
,n_k_ref_p_r ...
,k_ref_p_r_ ...
,k_ref_p_r_max ...
,n_w_ref_ ...
,[] ...
,weight_2d_ref_wk_ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% innerproduct_O_x_u_vs_S_x_u_2: %0.2fs',tmp_t)); end;
%%%%;
nx0 = 0; nx1 = 0;
SCTFCTFS_D = sum(CTF_S_flip_x_c_xx__(1+nx0+[0:n_x_u_pack-1],1+nx1+[0:n_x_u_pack-1]).*CTF_S_flip_x_c_xx__,'all');
SCTFCTFS_C = SCTFCTFS_xxz___(1+nx0+floor(n_x_u_pack/2)-1,1+nx1+floor(n_x_u_pack/2)-1,1+0);
disp(sprintf(' %% SCTFCTFS_D vs SCTFCTFS_C: %0.16f',fnorm(SCTFCTFS_D-SCTFCTFS_C)/fnorm(SCTFCTFS_D)));
%%%%;
nx0 = 3; nx1 = -7;
tmp_n2 = floor(n_x_u_pack/2);
SCTFCTFS_D = sum(CTF_S_flip_pad_x_c_xx__(n_x_u_pack+1+nx0+[0:n_x_u_pack-1],n_x_u_pack+1+nx1+[0:n_x_u_pack-1]).*CTF_S_flip_x_c_xx__,'all');
SCTFCTFS_C = SCTFCTFS_xxz___(1+nx0+floor(n_x_u_pack/2)-1,1+nx1+floor(n_x_u_pack/2)-1,1+0);
disp(sprintf(' %% SCTFCTFS_D vs SCTFCTFS_C: %0.16f',fnorm(SCTFCTFS_D-SCTFCTFS_C)/fnorm(SCTFCTFS_D)));
%%%%;
nx0 = -13; nx1 = +17;
tmp_n2 = floor(n_x_u_pack/2);
SCTFCTFS_D = sum(CTF_S_flip_pad_x_c_xx__(n_x_u_pack+1+nx0+[0:n_x_u_pack-1],n_x_u_pack+1+nx1+[0:n_x_u_pack-1]).*CTF_S_flip_x_c_xx__,'all');
SCTFCTFS_C = SCTFCTFS_xxz___(1+nx0+floor(n_x_u_pack/2)-1,1+nx1+floor(n_x_u_pack/2)-1,1+0);
disp(sprintf(' %% SCTFCTFS_D vs SCTFCTFS_C: %0.16f',fnorm(SCTFCTFS_D-SCTFCTFS_C)/fnorm(SCTFCTFS_D)));
%%%%%%%%;

%%%%%%%%;
% A direct calculation of CoC2 and CoC3. ;
% Used to compare R1S1_B_xxz___ to R1S1_C_xxz___ via CoC2 and CoC3. ;
%%%%%%%%;
n_R_0 = n_R_sub_x_u_pack_(1+0);
n_R_1 = n_R_sub_x_u_pack_(1+1);
n_R = prod(n_R_sub_x_u_pack_);
R_x_0_ = transpose(linspace(-1,+1,n_R_sub_x_u_pack_(1+0)));
R_x_1_ = transpose(linspace(-1,+1,n_R_sub_x_u_pack_(1+1)));
[R_x_0__,R_x_1__] = ndgrid(R_x_0_,R_x_1_);
R1S1_pad_B_xxz___ = zeros(2*n_x_u_pack+n_R_0,2*n_x_u_pack+n_R_1,n_gamma_z);
R1S1_pad_B_xxz___(n_x_u_pack+[1:n_R_0],n_x_u_pack+[1:n_R_1],:) = R1S1_B_xxz___;
R1S1_pad_C_xxz___ = zeros(2*n_x_u_pack+n_R_0,2*n_x_u_pack+n_R_1,n_gamma_z);
R1S1_pad_C_xxz___(n_x_u_pack+[1:n_R_0],n_x_u_pack+[1:n_R_1],:) = R1S1_C_xxz___;
CoC2_B_xx__ = zeros(n_R_0,n_R_1);
CoC3_B_xx__ = zeros(n_R_0,n_R_1);
CoC2_C_xx__ = zeros(n_R_0,n_R_1);
CoC3_C_xx__ = zeros(n_R_0,n_R_1);
tmp_t = tic();
nR=0;
for nx0=0:n_R_0-1; for nx1=0:n_R_1-1;
if (mod(nR,1024)==0); disp(sprintf(' %% nR %d/%d',nR,n_R)); end;
% nx0 = min( 2,n_R_0-1); nx1 = min(13,n_R_1-1);
tmp_n2 = floor(n_x_u_pack/2); dnx = 1;
OCTFS_xxz___ = R1S1_pad_B_xxz___(n_x_u_pack + nx0-(tmp_n2-1-dnx)+[1:n_x_u_pack],n_x_u_pack + nx1-(tmp_n2-1-dnx)+[1:n_x_u_pack],:);
%CoC2 = corr( reshape(SCTFCTFS_xxz___(:,:,1+0),[n_x_u_pack^2,1]) , reshape(OCTFS_xxz___(:,:,1+0),[n_x_u_pack^2,1]) );
%CoC3 = corr( reshape(SCTFCTFS_xxz___,[n_x_u_pack^2*n_gamma_z,1]) , reshape(OCTFS_xxz___,[n_x_u_pack^2*n_gamma_z,1]) );
CoC2 = dot( reshape(SCTFCTFS_xxz___(:,:,1+0),[n_x_u_pack^2,1]) , reshape(OCTFS_xxz___(:,:,1+0),[n_x_u_pack^2,1]) );
CoC3 = dot( reshape(SCTFCTFS_xxz___,[n_x_u_pack^2*n_gamma_z,1]) , reshape(OCTFS_xxz___,[n_x_u_pack^2*n_gamma_z,1]) );
CoC2_B_xx__(1+nx0,1+nx1) = CoC2;
CoC3_B_xx__(1+nx0,1+nx1) = CoC3;
%%%%;
OCTFS_xxz___ = R1S1_pad_C_xxz___(n_x_u_pack + nx0-(tmp_n2-1-dnx)+[1:n_x_u_pack],n_x_u_pack + nx1-(tmp_n2-1-dnx)+[1:n_x_u_pack],:);
%CoC2 = corr( reshape(SCTFCTFS_xxz___(:,:,1+0),[n_x_u_pack^2,1]) , reshape(OCTFS_xxz___(:,:,1+0),[n_x_u_pack^2,1]) );
%CoC3 = corr( reshape(SCTFCTFS_xxz___,[n_x_u_pack^2*n_gamma_z,1]) , reshape(OCTFS_xxz___,[n_x_u_pack^2*n_gamma_z,1]) );
CoC2 = dot( reshape(SCTFCTFS_xxz___(:,:,1+0),[n_x_u_pack^2,1]) , reshape(OCTFS_xxz___(:,:,1+0),[n_x_u_pack^2,1]) );
CoC3 = dot( reshape(SCTFCTFS_xxz___,[n_x_u_pack^2*n_gamma_z,1]) , reshape(OCTFS_xxz___,[n_x_u_pack^2*n_gamma_z,1]) );
CoC2_C_xx__(1+nx0,1+nx1) = CoC2;
CoC3_C_xx__(1+nx0,1+nx1) = CoC3;
%%%%;
clear OCTFS_xxz___ CoC2 CoC3 ;
nR=nR+1;
end;end;%for nx0=0:n_R_0-1; for nx1=0:n_R_1-1;
tmp_t = toc(tmp_t);
if (verbose>-1); disp(sprintf(' %% CoC2 and CoC3 direct calculation: %0.2fs',tmp_t/2)); end;
%%%%;
%imagesc((1+ij0_final_xx__-ij0_start_xx__).*(1+ij1_final_xx__-ij1_start_xx__)); %<-- area of lhs ij-array. ;
%imagesc((1+kl0_final_xx__-kl0_start_xx__).*(1+kl1_final_xx__-kl1_start_xx__)); %<-- area of rhs ij-array. ;
%%%%%%%%;
n_h = 128; np=0;
hlim_CoC2_ = [ min(CoC2_B_xx__,[],'all'),max(CoC2_B_xx__,[],'all') ];
hlim_CoC3_ = [ min(CoC3_B_xx__,[],'all'),max(CoC3_B_xx__,[],'all') ];
h2_CoC2_BC__ = hist2d_0(CoC2_B_xx__,CoC2_C_xx__,n_h,n_h,hlim_CoC2_,hlim_CoC2_);
h2_CoC3_BC__ = hist2d_0(CoC3_B_xx__,CoC3_C_xx__,n_h,n_h,hlim_CoC3_,hlim_CoC3_);
figure(1+nf);nf=nf+1;clf;figmed;p_row=2;p_col=3;np=0;
subplot(p_row,p_col,1+np);np=np+1;
imagesc_c(n_R_0,R_x_0_,n_R_1,R_x_1_,CoC2_B_xx__,[],colormap_beach);
axis image; axis square; axisnotick; title('CoC2_B_xx__','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
imagesc_c(n_R_0,R_x_0_,n_R_1,R_x_1_,CoC2_C_xx__,[],colormap_beach);
axis image; axis square; axisnotick; title('CoC2_C_xx__','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
imagesc(log2(1+h2_CoC2_BC__)); axis image; axisnotick; colorbar; set(gca,'ydir','normal');
title('CoC2_B_xx__ vs CoC2_C_xx__','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
imagesc_c(n_R_0,R_x_0_,n_R_1,R_x_1_,CoC3_B_xx__,[],colormap_beach);
axis image; axis square; axisnotick; title('CoC3_B_xx__','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
imagesc_c(n_R_0,R_x_0_,n_R_1,R_x_1_,CoC3_C_xx__,[],colormap_beach);
axis image; axis square; axisnotick; title('CoC3_C_xx__','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
imagesc(log2(1+h2_CoC3_BC__)); axis image; axisnotick; colorbar; set(gca,'ydir','normal');
title('CoC3_B_xx__ vs CoC3_C_xx__','Interpreter','none');
%%%%%%%%;

%%%%%%%%;
% Now construct kernel K2 and K3. ;
% Use K2 to construct CoC2_D and use K3 to construct CoC3_D. ;
%%%%%%%%;
K2_xx__ = SCTFCTFS_xxz___(:,:,1+0)*n_x_u_pack*n_x_u_pack/(n_x_u_pack*n_x_u_pack*dx_0*dx_1);
K2_k_p_wk_ = interp_x_c_to_k_p_xxnufft(n_S_x_0,diameter_x_c,n_S_x_1,diameter_x_c,K2_xx__,n_k_ref_p_r,k_ref_p_r_,n_w_ref_)*sqrt(n_S_x_0*n_S_x_1)*dx_0*dx_1;
K2_k_q_wk_ = interp_p_to_q(n_k_ref_p_r,n_w_ref_,n_w_ref_sum,K2_k_p_wk_);
K3_k_q_wk_ = zeros(n_w_ref_sum,1);
na=0;
for nk_ref_p_r=0:n_k_ref_p_r-1;
n_w_ref = n_w_ref_(1+nk_ref_p_r);
for nw_ref=0:n_w_ref-1;
if nw_ref==0; K3_k_q_wk_(1+na) = K2_k_q_wk_(1+na); end;
na=na+1;
end;%for nw_ref=0:n_w_ref-1;
end;%for nk_ref_p_r=0:n_k_ref_p_r-1;
K3_k_p_wk_ = interp_q_to_p(n_k_ref_p_r,n_w_ref_,n_w_ref_sum,K3_k_q_wk_);
K3_xx__ = real( interp_k_p_to_x_c_xxnufft(n_x_u_pack,diameter_x_c,n_x_u_pack,diameter_x_c,n_k_ref_p_r,k_ref_p_r_,n_w_ref_,K3_k_p_wk_.*weight_2d_ref_wk_*(2*pi)^2)*sqrt(n_x_u_pack*n_x_u_pack) * n_w_ref_sum );
%%%%;
if flag_check>1;
K3_B_xx__ = zeros(n_x_u_pack,n_x_u_pack);
for ngamma_z=0:n_gamma_z-1;
gamma_z = (2*pi*ngamma_z)/max(1,n_gamma_z);
tmp_K2_k_p_wk_ = rotate_p_to_p_fftw(n_k_ref_p_r,n_w_ref_,n_w_ref_sum,K2_k_p_wk_,+gamma_z);
tmp_K2_xx__ = real( interp_k_p_to_x_c_xxnufft(n_x_u_pack,diameter_x_c,n_x_u_pack,diameter_x_c,n_k_ref_p_r,k_ref_p_r_,n_w_ref_,tmp_K2_k_p_wk_.*weight_2d_ref_wk_*(2*pi)^2)*sqrt(n_x_u_pack*n_x_u_pack) * n_w_ref_sum );
K3_B_xx__ = K3_B_xx__ + tmp_K2_xx__;
clear tmp_K2_k_p_wk_ tmp_K2_xx__ ;
end;%for ngamma_z=0:n_gamma_z-1;
K3_B_xx__ = K3_B_xx__ / max(1,n_gamma_z);
figure(1+nf);nf=nf+1;clf;figsml;
plot(1:n_x_u_pack^2,K3_xx__(:),'ro',1:n_x_u_pack^2,K3_B_xx__(:),'gx');
xlabel('pixel'); ylabel('K3'); axis tight; grid on;
end;%if flag_check>1;
%%%%;
K3_xx__ = K3_xx__*n_gamma_z;
%%%%%%%%;
parameter_innerproduct = struct('type','parameter','nw_stride',2,'flag_conv2_vs_svt',1);
tmp_t = tic();
[ ...
 ~ ...
,CoC2_D_xxz___ ...
,scatter_from_tensor_order ...
,scatter_from_tensor_zswk___ ...
] = ...
innerproduct_O_x_u_vs_S_x_u_2( ...
 parameter_innerproduct ...
,n_R_sub_x_u_pack_(1+0) ...
,n_R_sub_x_u_pack_(1+1) ...
,Q_sub_x_u_pack_.^1 ...
,n_x_u_pack ...
,n_x_u_pack ...
,conv2(CTF_S_x_c_xx__,K2_xx__,'same')*dx_pack*dx_pack ...
,n_k_ref_p_r ...
,k_ref_p_r_ ...
,k_ref_p_r_max ...
,n_w_ref_ ...
,[] ...
,weight_2d_ref_wk_ ...
,[] ...
,scatter_from_tensor_order ...
,scatter_from_tensor_zswk___ ...
);
tmp_t = toc(tmp_t); if (verbose>-1); disp(sprintf(' %% CoC2 via innerproduct_O_x_u_vs_S_x_u_2: %0.2fs',tmp_t)); end;
%%%%%%%%;
parameter_innerproduct = struct('type','parameter','nw_stride',2,'flag_conv2_vs_svt',1);
tmp_t = tic();
[ ...
 ~ ...
,CoC3_D_xxz___ ...
,scatter_from_tensor_order ...
,scatter_from_tensor_zswk___ ...
] = ...
innerproduct_O_x_u_vs_S_x_u_2( ...
 parameter_innerproduct ...
,n_R_sub_x_u_pack_(1+0) ...
,n_R_sub_x_u_pack_(1+1) ...
,Q_sub_x_u_pack_.^1 ...
,n_x_u_pack ...
,n_x_u_pack ...
,conv2(CTF_S_x_c_xx__,K3_xx__,'same')*dx_pack*dx_pack ...
,n_k_ref_p_r ...
,k_ref_p_r_ ...
,k_ref_p_r_max ...
,n_w_ref_ ...
,[] ...
,weight_2d_ref_wk_ ...
,[] ...
,scatter_from_tensor_order ...
,scatter_from_tensor_zswk___ ...
);
tmp_t = toc(tmp_t); if (verbose>-1); disp(sprintf(' %% CoC3 via innerproduct_O_x_u_vs_S_x_u_2: %0.2fs',tmp_t)); end;
%%%%%%%%;
CoC2_D_xx__ = CoC2_D_xxz___(:,:,1+0);
CoC3_D_xx__ = CoC3_D_xxz___(:,:,1+0);
n_h = 128; np=0;
hlim_CoC2_ = [ min(CoC2_C_xx__,[],'all'),max(CoC2_C_xx__,[],'all') ];
hlim_CoC3_ = [ min(CoC3_C_xx__,[],'all'),max(CoC3_C_xx__,[],'all') ];
h2_CoC2_CD__ = hist2d_0(CoC2_C_xx__,CoC2_D_xx__,n_h,n_h,hlim_CoC2_,hlim_CoC2_);
h2_CoC3_CD__ = hist2d_0(CoC3_C_xx__,CoC3_D_xx__,n_h,n_h,hlim_CoC3_,hlim_CoC3_);
figure(1+nf);nf=nf+1;clf;figmed;p_row=2;p_col=3;np=0;
subplot(p_row,p_col,1+np);np=np+1;
imagesc_c(n_R_0,R_x_0_,n_R_1,R_x_1_,CoC2_C_xx__,[],colormap_beach);
axis image; axis square; axisnotick; title('CoC2_C_xx__','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
imagesc_c(n_R_0,R_x_0_,n_R_1,R_x_1_,CoC2_D_xx__,[],colormap_beach);
axis image; axis square; axisnotick; title('CoC2_D_xx__','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
imagesc(log2(1+h2_CoC2_CD__)); axis image; axisnotick; colorbar; set(gca,'ydir','normal');
title('CoC2_C_xx__ vs CoC2_D_xx__','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
imagesc_c(n_R_0,R_x_0_,n_R_1,R_x_1_,CoC3_C_xx__,[],colormap_beach);
axis image; axis square; axisnotick; title('CoC3_C_xx__','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
imagesc_c(n_R_0,R_x_0_,n_R_1,R_x_1_,CoC3_D_xx__,[],colormap_beach);
axis image; axis square; axisnotick; title('CoC3_D_xx__','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
imagesc(log2(1+h2_CoC3_CD__)); axis image; axisnotick; colorbar; set(gca,'ydir','normal');
title('CoC3_C_xx__ vs CoC3_D_xx__','Interpreter','none');
%%%%%%%%;

%clear S_x_u_xx__ S_k_p_wk_ M_x_u_xx__ M_k_p_wk_ ;
%clear CTF_S_k_p_wk_ CTF_S_x_c_xx__ ;
%clear SCTFCTFS_x_u_xxz___ OCTFS_xxz___ ;
