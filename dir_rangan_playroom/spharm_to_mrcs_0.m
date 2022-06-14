function ...
[ ...
 parameter ...
,fname_mrcs_nopath ...
,fname_star_nopath ...
,fname_xfix ...
,M_x_c___ ...
,euler_polar_a_true_ ...
,euler_azimu_b_true_ ...
,euler_gamma_z_true_ ...
,image_delta_x_true_ ...
,image_delta_y_true_ ...
] = ...
spharm_to_mrcs_0( ...
 parameter ...
,k_p_r_max ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_quad_ ...
,half_diameter_x_c ...
,n_x_u ...
,delta_sigma_true ...
,snr_per_AA ...
,n_M ...
,n_CTF ...
,Pixel_A ...
,Voltage ...
,Defocus ...
,Defocus_relative_spread ...
,Amplitude_Contrast ...
,dir_data ...
,str_prefix ...
);

verbose=1;
if (verbose); disp(sprintf(' %% [entering spharm_to_mrcs_0]')); end;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); k_p_r_max=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); k_p_r_=[]; end; na=na+1;
if (nargin<1+na); weight_3d_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); l_max_=[]; end; na=na+1;
if (nargin<1+na); a_k_Y_quad_=[]; end; na=na+1;
if (nargin<1+na); half_diameter_x_c=[]; end; na=na+1;
if (nargin<1+na); n_x_u=[]; end; na=na+1;
if (nargin<1+na); delta_sigma_true=[]; end; na=na+1;
if (nargin<1+na); snr_per_AA=[]; end; na=na+1;
if (nargin<1+na); n_M=[]; end; na=na+1;
if (nargin<1+na); n_CTF=[]; end; na=na+1;
if (nargin<1+na); Pixel_A=[]; end; na=na+1;
if (nargin<1+na); Voltage=[]; end; na=na+1;
if (nargin<1+na); Defocus=[]; end; na=na+1;
if (nargin<1+na); Defocus_relative_spread=[]; end; na=na+1;
if (nargin<1+na); Amplitude_Contrast=[]; end; na=na+1;
if (nargin<1+na); dir_data=[]; end; na=na+1;
if (nargin<1+na); str_prefix=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;

if isempty(half_diameter_x_c); half_diameter_x_c = 1.0; end; %<-- box diameter. ;
if isempty(n_x_u); n_x_u = 64; end;
if isempty(delta_sigma_true); delta_sigma_true = 0; end;
if isempty(snr_per_AA); snr_per_AA=0; end; %<-- signal-to-noise per angstrom-squared. ;
if isempty(n_M); n_M = 1024; end;
if isempty(n_CTF); n_CTF = max(1,floor(n_M/50)); end;
if isempty(Pixel_A); Pixel_A = 4.0; end; %<-- in angstroms. ;
if isempty(Voltage); Voltage = 300.0; end; %<-- in mV. ;
if isempty(Defocus); Defocus = 2e4; end; %<-- for niko. ;
if isempty(Defocus_relative_spread); Defocus_relative_spread = 0.25; end; %<-- for niko. ;
if isempty(Amplitude_Contrast); Amplitude_Contrast = 0.1; end; %<-- for niko. ;
if isempty(dir_data); dir_data = pwd(); end;
if isempty(str_prefix); str_prefix = 'test'; end;

% test with: ;
% parameter = struct('type','parameter'); n_x_u = 64; delta_sigma_true = 0.1; snr_per_AA = 0.15; n_M = 1024; n_CTF = max(1,floor(n_M/50)); Pixel_A = 4.0; Voltage = 300; Defocus = 2e4; Defocus_relative_spread = 0.25; Amplitude_Contrast = 0.1; dir_data = '/home/rangan/dir_cryoem/dir_FINYU0'; str_prefix = 'FINYU';

n_lm_ = (l_max_+1).^2;
n_lm_max = max(n_lm_);
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
l_max_max = max(l_max_); dWtdkd__l_max_max = 2*l_max_max;
n_w_max = 2*(l_max_max+1);
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

index_nS_from_nM_=0;
euler_polar_a_true_ = zeros(n_M,1);
euler_azimu_b_true_ = zeros(n_M,1);
euler_gamma_z_true_ = zeros(n_M,1);
image_delta_x_true_ = zeros(n_M,1);
image_delta_y_true_ = zeros(n_M,1);
rng(0);
nS=0;
for nM=0:n_M-1;
index_nS_from_nM_(1+nM) = nS;
euler_polar_a_true_(1+nM) = viewing_polar_a_all_(1+nS);
euler_azimu_b_true_(1+nM) = viewing_azimu_b_all_(1+nS);
euler_gamma_z_true_(1+nM) = 2*pi*rand(1);
image_delta_x_true_(1+nM) = delta_sigma_true*randn(1);
image_delta_y_true_(1+nM) = delta_sigma_true*randn(1);
nS=nS+1; if (nS>=n_S); nS=0; end;
end;%for nM=0:n_M-1;

%%%%%%%%;
% construct CTF, no spherical aberration. ;
%%%%%%%%;
SphericalAberration_CTF_ = 0.0*ones(n_CTF,1);
Voltage_CTF_ = Voltage*ones(n_CTF,1);
DefocusU_CTF_ = Defocus*ones(n_CTF,1).*(1.0 + Defocus_relative_spread*(2*rand(n_CTF,1)-1));
DefocusV_CTF_ = DefocusU_CTF_;
DefocusAngle_CTF_ = 0.0*ones(n_CTF,1);
AmplitudeContrast_CTF_ = Amplitude_Contrast*ones(n_CTF,1);
n_x_M_u = n_x_u;
%%%%%%%%;
index_nCTF_from_nM_ = max(0,min(n_CTF-1,floor(n_CTF*(0:n_M-1)/n_M)));
CTF_index_ = index_nCTF_from_nM_;
CTF_k_p__ = zeros(n_w_sum,n_CTF);
for nCTF=0:n_CTF-1;
if (mod(nCTF,100)==0); disp(sprintf(' %% nCTF %d/%d',nCTF,n_CTF)); end;
CTF_Spherical_Aberration = SphericalAberration_CTF_(1+nCTF);% spherical aberration of the lens in mm ;
CTF_Spherical_Aberration=CTF_Spherical_Aberration*(10.0d0^7.0d0);% convert into Angstroms ;
CTF_Voltage_kV = Voltage_CTF_(1+nCTF);% voltage in kVolts ;
CTF_Voltage_1V=CTF_Voltage_kV*1000.0 ;% convert into Volts ;
CTF_lambda = 12.2643247/sqrt(CTF_Voltage_1V+CTF_Voltage_1V^2*0.978466d-6);% electron wavelength in Angstroms ;
CTF_Defocus_U = DefocusU_CTF_(1+nCTF);% defocus values (in Angstroms) ;
CTF_Defocus_V = DefocusV_CTF_(1+nCTF);% defocus values (in Angstroms) ;
CTF_Defocus_Angle = DefocusAngle_CTF_(1+nCTF);% angle of astigmatism ;
CTF_Amplitude_Contrast = AmplitudeContrast_CTF_(1+nCTF);% CTF_Amplitude Contrast ;
tmp_w1=sqrt(1.0d0-CTF_Amplitude_Contrast^2);% weights for the amplitude and phase contrasts in CTF ;
tmp_w2=CTF_Amplitude_Contrast;% weights for the amplitude and phase contrasts in CTF ;
%  CTF_Object_Pixel_Size = CTF_Detector_Pixel_Size/CTF_Magnification;
CTF_Object_Pixel_Size = Pixel_A;% pixel size of the scanner in physical space in Angstroms ;
CTF_lambda_per_box = CTF_lambda/(n_x_M_u*CTF_Object_Pixel_Size);% n_x_M_u*CTF_Object_Pixel_Size is the box size in Angstroms ;
%%%%;
na=0;
for nk = 0:n_k_p_r-1;
for nw=0:n_w_(1+nk)-1;
tmp_theta = (2.0d0*pi*nw)/n_w_(1+nk);
tmp_k_c_1 = (2.0d0*pi)*k_p_r_(1+nk)*cos(tmp_theta);
tmp_k_c_2 = (2.0d0*pi)*k_p_r_(1+nk)*sin(tmp_theta);
tmp_ctf_value = niko_ctf(CTF_Spherical_Aberration,CTF_lambda,tmp_w1,tmp_w2,CTF_Defocus_U,CTF_Defocus_V,CTF_Defocus_Angle,CTF_lambda_per_box/pi,tmp_k_c_1,tmp_k_c_2);
CTF_k_p__(1+na,1+nCTF) = -tmp_ctf_value;
na = na+1;
end;%for nw=0:n_w_(1+nk)-1;
end;%for nk = 0:n_k_p_r-1;
%%%%;
end;%for nCTF=0:n_CTF-1;
%%%%%%%%;
CTF_k_p_r__ = zeros(n_k_p_r,n_CTF);
for nCTF=0:n_CTF-1;
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
CTF_k_p_r__(1+nk_p_r,1+nCTF) = mean(CTF_k_p__(1+tmp_index_,1+nCTF));
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nCTF=0:n_CTF-1;
CTF_avg_k_p_ = mean(CTF_k_p__,2);
%imagesc_p(n_k_p_r,k_p_r_,n_w_,sum(n_w_),real(CTF_avg_k_p_(:)),[-1,+1],colormap_beach());
CTF_avg_k_p_r_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
CTF_avg_k_p_r_(1+nk_p_r) = mean(CTF_avg_k_p_(1+tmp_index_));
end;%for nk_p_r=0:n_k_p_r-1;

%%%%%%%%;
% Now generate clean images. ;
%%%%%%%%;
tmp_t = tic();
T_k_p__ = zeros(n_w_sum,n_M);
T_x_c___ = zeros(n_x_u,n_x_u,n_M);
S_x_c___ = zeros(n_x_u,n_x_u,n_M);
for nM=0:n_M-1;
nS = index_nS_from_nM_(1+nM);
tmp_euler_gamma_z = +euler_gamma_z_true_(1+nM);
tmp_image_delta_x = +1.0*image_delta_x_true_(1+nM);
tmp_image_delta_y = +1.0*image_delta_y_true_(1+nM);
S_k_p_ = S_k_p__(:,1+nS);
S_x_c__ = interp_k_p_to_x_c_xxnufft(n_x_u,2*half_diameter_x_c,n_x_u,2*half_diameter_x_c,n_k_p_r,k_p_r_,n_w_,S_k_p_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_M_u^2) * n_w_sum;
S_x_c__ = real(S_x_c__);
S_x_c__ = S_x_c__/fnorm(S_x_c__);
S_x_c___(:,:,1+nM) = S_x_c__;
T_k_p_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_,+tmp_euler_gamma_z);
T_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,T_k_p_,-tmp_image_delta_x,-tmp_image_delta_y);
T_k_p_ = sparse(1:n_w_sum,1:n_w_sum,CTF_k_p__(:,1+CTF_index_(1+nM)),n_w_sum,n_w_sum)*T_k_p_;
T_k_p__(:,1+nM) = T_k_p_;
T_x_c__ = interp_k_p_to_x_c_xxnufft(n_x_u,2*half_diameter_x_c,n_x_u,2*half_diameter_x_c,n_k_p_r,k_p_r_,n_w_,T_k_p_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_M_u^2) * n_w_sum;
T_x_c__ = real(T_x_c__);
T_x_c__ = T_x_c__/fnorm(T_x_c__);
T_x_c___(:,:,1+nM) = T_x_c__;
end;%for nM=0:n_M-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% T_k_p__: %0.2fs',tmp_t)); end;
%%%%%%%%;
% add noise. ;
%%%%%%%%;
T_x_c_l2 = mean(sqrt(sum(sum(abs(T_x_c___).^2,1),2)));
sigma_noise_per_Pixel = 0;
if (snr_per_AA>0); sigma_noise_per_Pixel = T_x_c_l2/max(1e-12,snr_per_AA)/sqrt(n_x_u^2)/sqrt(Pixel_A^2); end;
M_x_c___ = T_x_c___ + sigma_noise_per_Pixel*rand(n_x_u,n_x_u,n_M);
%%%%%%%%;
flag_plot=0;
if flag_plot;
figure(1);clf;figmed;figbeach();
nM = 500; 
Mlim_ = [min(real(M_x_c___),[],'all') , max(real(M_x_c___),[],'all')];
subplot(1,2,1); imagesc(T_x_c___(:,:,1+nM),Mlim_); axisnotick;
subplot(1,2,2); imagesc(M_x_c___(:,:,1+nM),Mlim_); axisnotick;
end;%if flag_plot;
%%%%%%%%;

%%%%%%%%;
% Now write mrcs. ;
%%%%%%%%;
[ ...
 parameter ...
,fname_xfix ...
] = ...
spharm_to_mrcs_xfix_0( ...
 parameter ...
,n_x_u ...
,delta_sigma_true ...
,snr_per_AA ...
,n_M ...
,n_CTF ...
,Pixel_A ...
,Voltage ...
,Defocus ...
,Defocus_relative_spread ...
,Amplitude_Contrast ...
);
fname_mrcs_prefix = sprintf('%s_%s',str_prefix,fname_xfix);
fname_mrcs_nopath = sprintf('%s.mrcs',fname_mrcs_prefix);
fname_mrcs = sprintf('%s/%s',dir_data,fname_mrcs_nopath);
fp_mrcs = WriteMRCHeader(real(M_x_c___),Pixel_A,fname_mrcs,n_M);
fwrite(fp_mrcs,real(M_x_c___),'float32');
fclose(fp_mrcs);

%%%%%%%%;
% Now write primitive star file. ;
%%%%%%%%;
fname_star_prefix = sprintf('%s_%s',str_prefix,fname_xfix);
fname_star_nopath = sprintf('%s.star',fname_star_prefix);
fname_star = sprintf('%s/%s',dir_data,fname_star_nopath);
fp_star = fopen(fname_star,'w');
tmp_ImageName_list_ = cell(n_M,1);
tmp_OpticsGroup_list_ = cell(n_M,1);
tmp_DefocusU_list_ = cell(n_M,1);
tmp_DefocusV_list_ = cell(n_M,1);
tmp_DefocusAngle_list_ = cell(n_M,1);
for nM=0:n_M-1;
nCTF = index_nCTF_from_nM_(1+nM);
tmp_rlnImageName = sprintf('%d@%s',1+nM,fname_mrcs_nopath);
tmp_ImageName_list_{1+nM} = tmp_rlnImageName;
tmp_OpticsGroup_list_{1+nM} = sprintf('1');
tmp_DefocusU_list_{1+nM} = num2str(DefocusU_CTF_(1+nCTF));
tmp_DefocusV_list_{1+nM} = num2str(DefocusV_CTF_(1+nCTF));
tmp_DefocusAngle_list_{1+nM} = num2str(DefocusAngle_CTF_(1+nCTF));
end;%for nM=0:n_M-1;
tmp_Voltage = Voltage;
tmp_SphericalAberration = SphericalAberration_CTF_(1+nCTF);;
tmp_AmplitudeContrast = Amplitude_Contrast;
fprintf(fp_star,'# version 30001\n');
fprintf(fp_star,'\n');
fprintf(fp_star,'data_optics\n');
fprintf(fp_star,'\n');
fprintf(fp_star,'loop_\n');
fprintf(fp_star,'_rlnOpticsGroupName #1\n');
fprintf(fp_star,'_rlnOpticsGroup #2\n');
fprintf(fp_star,'_rlnMtfFileName #3\n');
fprintf(fp_star,'_rlnVoltage #4\n');
fprintf(fp_star,'_rlnSphericalAberration #5\n');
fprintf(fp_star,'_rlnAmplitudeContrast #6\n');
fprintf(fp_star,'_rlnImagePixelSize #7\n');
fprintf(fp_star,'_rlnImageSize #8\n');
fprintf(fp_star,'_rlnImageDimensionality #9\n');
fprintf(fp_star,'opticsGroup1 1 mtf_k2_%dkV.star %f %f %f %f %d 2\n',tmp_Voltage,tmp_Voltage,tmp_SphericalAberration,tmp_AmplitudeContrast,Pixel_A,n_x_u);
fprintf(fp_star,'\n');
fprintf(fp_star,'# version 30001\n');
fprintf(fp_star,'\n');
fprintf(fp_star,'data_particles\n');
fprintf(fp_star,'\n');
fprintf(fp_star,'loop_\n');
fprintf(fp_star,'_rlnImageName #1\n');
fprintf(fp_star,'_rlnOpticsGroup #2\n');
fprintf(fp_star,'_rlnDefocusU #3\n');
fprintf(fp_star,'_rlnDefocusV #4\n');
fprintf(fp_star,'_rlnDefocusAngle #5\n');
for nM=0:n_M-1;
fprintf(fp_star,'%s %s %s %s %s\n',tmp_ImageName_list_{1+nM},tmp_OpticsGroup_list_{1+nM},tmp_DefocusU_list_{1+nM},tmp_DefocusV_list_{1+nM},tmp_DefocusAngle_list_{1+nM});
end;%for nM=0:n_M-1;
fclose(fp_star);
%%%%%%%%;

flag_check=1;
if flag_check;
%%%%%%%%;
% check result. ;
%%%%%%%%;
[ ...
 tmp_M_x_c___ ...
,tmp_index_nCTF_from_nM_ ...
,tmp_index_nM_from_nCTF_ ...
,tmp_Voltage_CTF_ ...
,tmp_DefocusU_CTF_ ...
,tmp_DefocusV_CTF_ ...
,tmp_DefocusAngle_CTF_ ...
,tmp_SphericalAberration_CTF_ ...
,tmp_AmplitudeContrast_CTF_ ...
] = ...
rlnImageName_from_star_1( ...
 dir_data ...
,fname_star_nopath ...
,n_M ...
);
disp(sprintf(' %% M_x_c___ vs tmp_M_x_c___: %0.16f',fnorm(M_x_c___ - tmp_M_x_c___)/fnorm(M_x_c___)));
%%%%;
error_Voltage=0;
for nM=0:n_M-1;
error_Voltage = error_Voltage + abs(Voltage_CTF_(1+index_nCTF_from_nM_(1+nM)) - tmp_Voltage_CTF_(1+tmp_index_nCTF_from_nM_(1+nM)));
end;%for nM=0:n_M-1;
error_Voltage = error_Voltage / max(1e-12,fnorm(Voltage_CTF_(1+index_nCTF_from_nM_)));
disp(sprintf(' %% error_Voltage: %0.16f',error_Voltage));
%%%%;
error_DefocusU=0;
for nM=0:n_M-1;
error_DefocusU = error_DefocusU + abs(DefocusU_CTF_(1+index_nCTF_from_nM_(1+nM)) - tmp_DefocusU_CTF_(1+tmp_index_nCTF_from_nM_(1+nM)));
end;%for nM=0:n_M-1;
error_DefocusU = error_DefocusU / max(1e-12,fnorm(DefocusU_CTF_(1+index_nCTF_from_nM_)));
disp(sprintf(' %% error_DefocusU: %0.16f',error_DefocusU));
%%%%;
error_DefocusV=0;
for nM=0:n_M-1;
error_DefocusV = error_DefocusV + abs(DefocusV_CTF_(1+index_nCTF_from_nM_(1+nM)) - tmp_DefocusV_CTF_(1+tmp_index_nCTF_from_nM_(1+nM)));
end;%for nM=0:n_M-1;
error_DefocusV = error_DefocusV / max(1e-12,fnorm(DefocusV_CTF_(1+index_nCTF_from_nM_)));
disp(sprintf(' %% error_DefocusV: %0.16f',error_DefocusV));
%%%%;
error_SphericalAberration=0;
for nM=0:n_M-1;
error_SphericalAberration = error_SphericalAberration + abs(SphericalAberration_CTF_(1+index_nCTF_from_nM_(1+nM)) - tmp_SphericalAberration_CTF_(1+tmp_index_nCTF_from_nM_(1+nM)));
end;%for nM=0:n_M-1;
error_SphericalAberration = error_SphericalAberration / max(1e-12,fnorm(SphericalAberration_CTF_(1+index_nCTF_from_nM_)));
disp(sprintf(' %% error_SphericalAberration: %0.16f',error_SphericalAberration));
%%%%;
error_AmplitudeContrast=0;
for nM=0:n_M-1;
error_AmplitudeContrast = error_AmplitudeContrast + abs(AmplitudeContrast_CTF_(1+index_nCTF_from_nM_(1+nM)) - tmp_AmplitudeContrast_CTF_(1+tmp_index_nCTF_from_nM_(1+nM)));
end;%for nM=0:n_M-1;
error_AmplitudeContrast = error_AmplitudeContrast / max(1e-12,fnorm(AmplitudeContrast_CTF_(1+index_nCTF_from_nM_)));
disp(sprintf(' %% error_AmplitudeContrast: %0.16f',error_AmplitudeContrast));
%%%%;
end;%if flag_check;

if (verbose); disp(sprintf(' %% [finished spharm_to_mrcs_0]')); end;










