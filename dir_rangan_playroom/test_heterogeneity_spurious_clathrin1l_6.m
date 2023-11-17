%%%%%%%%;
% Tests heterogeneity_spurious on 6sct_one_leg (clathrin). ;
% Constructs images and various spurious volumes. ;
% Adding calculation of bayesian likelihood. ;
%{
% contents of dir_Kexin_mip_vs_CoC: ;
result_search_5_2_radius_10.csv
result_search_5_2_radius_50.csv
result_search_5_3_radius_10.csv
result_search_5_3_radius_50.csv
result_search_9_2_radius_10.csv
result_search_9_2_radius_50.csv
result_search_9_3_radius_10.csv
result_search_9_3_radius_50.csv
simulator_default_clathrin_exp2_t1000A_d5000A_1.06Apix_45e_montage_boxsize_512_coc_9_2.mrc
simulator_default_clathrin_exp2_t1000A_d5000A_1.06Apix_45e_montage_boxsize_512_coc_9_3.mrc
simulator_default_clathrin_exp2_t1000A_d5000A_1.06Apix_45e_montage_boxsize_512_scaled_mip_9_2.mrc
simulator_default_clathrin_exp2_t1000A_d5000A_1.06Apix_45e_montage_boxsize_512_scaled_mip_9_3.mrc
simulator_default_mature60S_exp2_t1000A_d5000A_1.06Apix_45e_box512_montage_coc_5_2.mrc
simulator_default_mature60S_exp2_t1000A_d5000A_1.06Apix_45e_box512_montage_coc_5_3.mrc
simulator_default_mature60S_exp2_t1000A_d5000A_1.06Apix_45e_box512_montage_scaled_mip_5_2.mrc
simulator_default_mature60S_exp2_t1000A_d5000A_1.06Apix_45e_box512_montage_scaled_mip_5_3.mrc
  %}
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
flag_recalc = 0;
flag_realign = 0;
flag_replot = 0;
flag_center = 0;
tolerance_master = 1e-2;
nf=0;
%%%%;
dir_base = sprintf('/%s/rangan/dir_cryoem',string_root);
dir_Kexin = sprintf('%s/dir_Kexin_mip_vs_CoC',dir_base);
dir_Kexin_jpg = sprintf('%s/dir_Kexin_jpg',dir_Kexin);
if ~exist(dir_Kexin_jpg,'dir'); disp(sprintf(' %% mkdir %s',dir_Kexin_jpg)); mkdir(dir_Kexin_jpg); end;
dir_Kexin_star = sprintf('%s/dir_Ground_Truth_star',dir_Kexin);
%%%%;
dir_data = sprintf('%s/dir_clathrin_montage_Micrograph',dir_base);
dir_pm = sprintf('%s/dir_pm',dir_data);
if (~exist(sprintf('%s_mat',dir_pm),'dir')); disp(sprintf(' %% mkdir %s_mat',dir_pm)); mkdir(sprintf('%s_mat',dir_pm)); end;
if (~exist(sprintf('%s_jpg',dir_pm),'dir')); disp(sprintf(' %% mkdir %s_jpg',dir_pm)); mkdir(sprintf('%s_jpg',dir_pm)); end;
Pixel_Spacing = 1.00; %<-- in angstroms; I am unsure about this, but 1.00 is unlikely to be very wrong. ;
fname_nopath_volume = '6sct_one_leg_bf60_center.mrc';
flag_center_volume = 0;
%%%%%%%%;

%%%%%%%%;
% load the original volume, as well as rotated versions of the original volume. ;
%%%%%%%%;
test_heterogeneity_spurious_clathrin1l_helper_a_k_Y_0;
x_u_0_ = transpose(linspace(-half_diameter_x_c,+half_diameter_x_c,n_x_u));
x_u_1_ = transpose(linspace(-half_diameter_x_c,+half_diameter_x_c,n_x_u));
x_u_pack_0_ = transpose(linspace(-half_diameter_x_c,+half_diameter_x_c,n_x_u_pack));
x_u_pack_1_ = transpose(linspace(-half_diameter_x_c,+half_diameter_x_c,n_x_u_pack));
dx_u = diameter_x_c/max(1,n_x_u);
dx_u_pack = diameter_x_c/max(1,n_x_u_pack);
a_p123_x_u_base_ = a_x_u_base_;
a_p123_k_p_quad_ = a_k_p_quad_;
a_p123_x_u_reco_ = a_x_u_reco_;
a_p123_k_Y_quad_ = a_k_Y_quad_;
a_p123_k_p_reco_ = a_k_p_reco_;
a_p123_k_Y_quad__ = a_k_Y_quad__;
S_p123_k_p__ = S_k_p__;
test_heterogeneity_spurious_clathrin1l_helper_a_p231_k_Y_0;
test_heterogeneity_spurious_clathrin1l_helper_a_p312_k_Y_0;
S_p123_l2_S_ = S_l2_S_;
test_heterogeneity_spurious_clathrin1l_helper_a_pxxx_k_Y_1; %<-- normalize volumes. ;
%%%%%%%%;
% Now load CTFs and set up principal-modes. ;
%%%%%%%%;
test_heterogeneity_spurious_clathrin1l_helper_CTF_0; %<-- create CTFs. ;
test_heterogeneity_spurious_clathrin1l_helper_CTF_1; %<-- sets up FTK, FTK_0 ;
%%%%%%%%;
% Now load the montage micrograph. ;
%%%%%%%%;
test_heterogeneity_spurious_clathrin1l_helper_O_x_u_0; %<-- load with n_pack_6. ;
test_heterogeneity_spurious_clathrin1l_helper_O_x_u_1; %<-- find residuals and form simulated images. ;
test_heterogeneity_spurious_clathrin1l_helper_readtable_0; %<-- load Kexin table of assigned viewing-angles. ;
%%%%%%%%;
% Now create the spurious volumes. ;
%%%%%%%%;
test_heterogeneity_spurious_clathrin1l_helper_a_p123_pole_k_Y_0;
test_heterogeneity_spurious_clathrin1l_helper_a_p231_pole_k_Y_0;
test_heterogeneity_spurious_clathrin1l_helper_a_p312_pole_k_Y_0;
test_heterogeneity_spurious_clathrin1l_helper_viewing_weight_0; %<-- define empirical viewing weights. ;
%%%%%%%%;
flag_replot = 1;
flag_recalc = 0;
flag_realign = 0; 
n_R_ = [0,32,64,128,256,512,1024]; n_n_R = numel(n_R_);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for nn_R=0:n_n_R-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
n_R = n_R_(1+nn_R);
M2_k_p__ = [M_sim_k_p_wkM__ , R_mon_k_p_wkM__(:,1:n_R)]; %<-- define the expanded set of images. ;
M2_k_q__ = [M_sim_k_q_wkM__ , R_mon_k_q_wkM__(:,1:n_R)]; %<-- define the expanded set of images. ;
n_M2 = n_M + n_R;
if (verbose); disp(sprintf(' %% n_R %d --> n_M2 %d',n_R,n_M2)); end;
%%%%%%%%%%%%%%%%;
fname_ssnll_pole_pre = sprintf('%s_mat/test_heterogeneity_spurious_clathrin1l_helper_ssnll_pole_R%.4d_',dir_pm,n_R);
[flag_ssnll_pole_skip,fname_ssnll_pole_mat] = open_fname_tmp(fname_ssnll_pole_pre);
fname_ssnll_pol3_pre = sprintf('%s_mat/test_heterogeneity_spurious_clathrin1l_helper_ssnll_pol3_R%.4d_',dir_pm,n_R);
[flag_ssnll_pol3_skip,fname_ssnll_pol3_mat] = open_fname_tmp(fname_ssnll_pol3_pre);
if flag_recalc | ~flag_ssnll_pole_skip | ~flag_ssnll_pol3_skip ; 
test_heterogeneity_spurious_clathrin1l_helper_p123_pol3_0; %<-- align noisy-images to polar-cap. ;
test_heterogeneity_spurious_clathrin1l_helper_p231_pol3_0; %<-- align noisy-images to polar-cap. ;
test_heterogeneity_spurious_clathrin1l_helper_p312_pol3_0; %<-- align noisy-images to polar-cap. ;
test_heterogeneity_spurious_clathrin1l_helper_ssnll_pole_0; %<-- calculate ssnll for polar-cap. ;
test_heterogeneity_spurious_clathrin1l_helper_ssnll_pol3_0; %<-- calculate ssnll for noise-aligned polar-cap. ;
end;%if flag_recalc | ~flag_ssnll_pole_skip | ~flag_ssnll_pol3_skip ; 
%%%%%%%%;
if ~exist(fname_ssnll_pole_mat,'file'); disp(sprintf(' %% Warning, %s not found',fname_ssnll_pole_mat)); end;
if  exist(fname_ssnll_pole_mat,'file'); load(fname_ssnll_pole_mat);
test_heterogeneity_spurious_clathrin1l_helper_ssnll_pole_1; %<-- plot graph. ;
end;%if  exist(fname_ssnll_pol3_mat,'file'); 
%%%%;
if ~exist(fname_ssnll_pol3_mat,'file'); disp(sprintf(' %% Warning, %s not found',fname_ssnll_pol3_mat)); end;
if  exist(fname_ssnll_pol3_mat,'file'); load(fname_ssnll_pol3_mat);
test_heterogeneity_spurious_clathrin1l_helper_ssnll_pol3_1; %<-- plot graph. ;
end;%if  exist(fname_ssnll_pol3_mat,'file'); 
%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for nn_R=0:n_n_R-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

disp('returning'); return;



