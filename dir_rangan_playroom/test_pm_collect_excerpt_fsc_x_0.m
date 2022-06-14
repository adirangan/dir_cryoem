function ...
[ ...
 parameter ...
,fsc_reco_k_ ...
,fsc_crop_reco_kx__ ...
,k_Ainv_p_r_ ...
,k_Ainv_p_r_max ...
,kinv_A_p_r_ ...
] = ...
test_pm_collect_excerpt_fsc_x_0( ...
 parameter ...
,dir_pm ...
);

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); dir_pm=[]; end; na=na+1;
if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_recalc'); parameter.flag_recalc=0; end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
if ~isfield(parameter,'Pixel_Spacing'); parameter.Pixel_Spacing=1; end;
flag_recalc = parameter.flag_recalc;
flag_verbose = parameter.flag_verbose;
Pixel_Spacing = parameter.Pixel_Spacing;
if isempty(dir_pm); dir_pm = pwd; end;

dir_pm_mat = sprintf('%s_mat',dir_pm);

%%%%%%%%;
% find cropped (i.e., masked) fsc between ground-truth reconstruction and published molecule. ;
%%%%%%%%;
fsc_reco_k_ = [];
corr_reco = [];
fsc_crop_reco_kx__ = [];
corr_crop_reco_x_ = [];
k_Ainv_p_r_ = [];
k_Ainv_p_r_max = [];
kinv_A_p_r_ = [];
%%%%%%%%;
fname_pre = sprintf('%s/test_pm_collect_fsc_crop_reco_',dir_pm_mat);
[fname_skip,fname_mat] = open_fname_tmp(fname_pre);
if ( flag_recalc | ~fname_skip );
disp(sprintf(' %% %s not found, creating',fname_mat));
tmp_fname = sprintf('%s/M_k_p__.mat',dir_pm_mat);
if ~exist(tmp_fname,'file'); disp(sprintf(' %% Warning, %s not found in test_pm_collect_excerpt_fsc_x_0',tmp_fname)); return; end;
tmp_ = load(tmp_fname,'n_x_M_u','x_c_0_');
tmp_flag_exist=0;
if (isfield(tmp_,'n_x_M_u')); tmp_flag_exist=1; n_x_M_u = tmp_.n_x_M_u; end;
if (isfield(tmp_,'x_c_0_')); tmp_flag_exist=1; n_x_M_u = numel(tmp_.x_c_0_); end;
if ~tmp_flag_exist; disp(sprintf(' %% Warning, %s not found in test_pm_collect_excerpt_fsc_x_0','n_x_M_u')); return; end;
clear tmp_;
%%%%%%%%;
tmp_fname = sprintf('%s/a_k_p_quad_.mat',dir_pm_mat);
if ~exist(tmp_fname,'file'); disp(sprintf(' %% Warning, %s not found in test_pm_collect_excerpt_fsc_x_0',tmp_fname)); return; end;
tmp_ = load(tmp_fname);
k_p_r_max = tmp_.k_p_r_max;
k_eq_d = tmp_.k_eq_d;
n_k_all = tmp_.n_k_all;
n_k_all_csum_ = tmp_.n_k_all_csum_;
k_p_r_all_ = tmp_.k_p_r_all_;
k_p_azimu_b_all_ = tmp_.k_p_azimu_b_all_;
k_p_polar_a_all_ = tmp_.k_p_polar_a_all_;
weight_3d_k_all_ = tmp_.weight_3d_k_all_;
weight_shell_k_ = tmp_.weight_shell_k_;
n_k_p_r = tmp_.n_k_p_r;
k_p_r_ = tmp_.k_p_r_;
weight_3d_k_p_r_ = tmp_.weight_3d_k_p_r_;
k_c_0_all_ = tmp_.k_c_0_all_;
k_c_1_all_ = tmp_.k_c_1_all_;
k_c_2_all_ = tmp_.k_c_2_all_;
a_k_p_quad_ = tmp_.a_k_p_quad_;
a_x_u_reco_ = tmp_.a_x_u_reco_;
clear tmp_;
%%%%%%%%;
tmp_fname = sprintf('%s/a_k_Y_quad_.mat',dir_pm_mat);
if ~exist(tmp_fname,'file'); disp(sprintf(' %% Warning, %s not found in test_pm_collect_excerpt_fsc_x_0',tmp_fname)); return; end;
tmp_ = load(tmp_fname,'l_max_','a_k_Y_quad_');
l_max_ = tmp_.l_max_;
a_k_Y_quad_ = tmp_.a_k_Y_quad_;
clear tmp_;
%%%%%%%%;
tmp_flag_exist=0;
tmp_fname = sprintf('%s/a_x_u_base_.mat',dir_pm_mat);
if exist(tmp_fname,'file');
tmp_flag_exist=1;
tmp_ = load(tmp_fname,'half_diameter_x_c','x_u_0_','x_u_1_','x_u_2_','n_x_u_pack','a_x_u_base_');
a_x_u_base_ = tmp_.a_x_u_base_;
end;%if exist(tmp_fname,'file');
tmp_fname = sprintf('%s/a_x_u_pack_.mat',dir_pm_mat);
if exist(tmp_fname,'file');
tmp_flag_exist=1;
tmp_ = load(tmp_fname,'half_diameter_x_c','x_u_0_','x_u_1_','x_u_2_','n_x_u_pack','a_x_u_pack_');
a_x_u_base_ = tmp_.a_x_u_pack_;
end;%if exist(tmp_fname,'file');
if ~tmp_flag_exist; disp(sprintf(' %% Warning, %s not found in test_pm_collect_excerpt_fsc_x_0',tmp_fname)); return; end;
half_diameter_x_c = tmp_.half_diameter_x_c;
x_u_0_ = tmp_.x_u_0_; x_u_1_ = tmp_.x_u_1_; x_u_2_ = tmp_.x_u_2_; n_x_u_pack = tmp_.n_x_u_pack;
[x_u_0___,x_u_1___,x_u_2___] = ndgrid(x_u_0_,x_u_1_,x_u_2_); n_xxx_u = n_x_u_pack^3;
x_u_r___ = sqrt(x_u_0___.^2 + x_u_1___.^2 + x_u_2___.^2);
clear tmp_;
%%%%%%%%;
tmp_flag_exist = 0;
tmp_fname = sprintf('%s/a_k_Y_0lsq_reco_.mat',dir_pm_mat);
if exist(tmp_fname,'file');
tmp_flag_exist = 1;
tmp_ = load(tmp_fname);
b_k_Y_reco_ = tmp_.a_k_Y_0lsq_reco_;
if isfield(tmp_,'a_k_p_0lsq_reco_'); b_k_p_reco_ = tmp_.a_k_p_0lsq_reco_; end;
if isfield(tmp_,'a_k_p_0lsq_reco'); b_k_p_reco_ = tmp_.a_k_p_0lsq_reco; end; %<-- typo in earlier version. ;
b_x_u_reco_ = tmp_.a_x_u_0lsq_reco_;
end;%if exist(tmp_fname,'file');
tmp_fname = sprintf('%s/c_k_Y_.mat',dir_pm_mat);
if exist(tmp_fname,'file');
tmp_flag_exist = 1;
tmp_ = load(tmp_fname,'c_k_Y_reco_','X_best_reco','c_k_p_reco_','c_x_u_reco_');
b_k_Y_reco_ = tmp_.c_k_Y_reco_;
b_k_p_reco_ = tmp_.c_k_p_reco_;
b_x_u_reco_ = tmp_.c_x_u_reco_;
end;%if exist(tmp_fname,'file');
if ~tmp_flag_exist; disp(sprintf(' %% Warning, %s not found in test_pm_collect_excerpt_fsc_x_0',tmp_fname)); return; end;
clear tmp_;
%%%%;
fsc_reco_k_ = zeros(n_k_p_r,1);
corr_reco = zeros(1,1);
fsc_crop_reco_kx__ = zeros(n_k_p_r,n_x_u_pack);
corr_crop_reco_x_ = zeros(n_x_u_pack,1);
k_Ainv_p_r_ = (2*k_p_r_)/(n_x_M_u * Pixel_Spacing);
k_Ainv_p_r_max = (2*k_p_r_max)/(n_x_M_u * Pixel_Spacing);
kinv_A_p_r_ = 1./max(1e-12,k_Ainv_p_r_);
%%%%%%%%;
parameter_fsc = struct('type','parameter');
parameter_fsc.flag_register = 0;
[ ...
 parameter_fsc ...
,fsc_reco_k_ ...
] = ...
fsc_from_a_k_Y_0( ...
 parameter_fsc ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,0 ...
,l_max_ ...
,a_k_Y_quad_ ...
,b_k_Y_reco_ ...
);
%%%%;
corr_reco = real(corr( a_x_u_base_(:) , b_x_u_reco_(:) ));
n_crop = n_x_u_pack;
for ncrop=0:n_crop-1;
if (flag_verbose>-1); disp(sprintf(' %% ncrop %d/%d',ncrop,n_crop)); end;
r_crop = 1.0*ncrop/(n_crop-1);
tmp_m_x_u_ = reshape(x_u_r___<=r_crop,[n_xxx_u,1]);
corr_crop_reco_x_(1+ncrop) = real(corr( a_x_u_base_(:).*tmp_m_x_u_ , b_x_u_reco_(:).*tmp_m_x_u_ ));
[ ...
 c_k_p_base_ ...
] = ...
convert_x_c_to_k_p_1( ...
 flag_verbose ...
,n_k_all ...
,weight_3d_k_all_ ...
,k_c_0_all_ ...
,k_c_1_all_ ...
,k_c_2_all_ ...
,k_p_r_max ...
,n_xxx_u ...
,[] ...
,x_u_0___ ...
,x_u_1___ ...
,x_u_2___ ...
,half_diameter_x_c ...
,a_x_u_base_(:).*tmp_m_x_u_ ...
);
[ ...
 d_k_p_reco_ ...
] = ...
convert_x_c_to_k_p_1( ...
 flag_verbose ...
,n_k_all ...
,weight_3d_k_all_ ...
,k_c_0_all_ ...
,k_c_1_all_ ...
,k_c_2_all_ ...
,k_p_r_max ...
,n_xxx_u ...
,[] ...
,x_u_0___ ...
,x_u_1___ ...
,x_u_2___ ...
,half_diameter_x_c ...
,b_x_u_reco_(:).*tmp_m_x_u_ ...
);
parameter_fsc = struct('type','parameter');
[ ...
 parameter_fsc ...
,fsc_crop_reco_k_ ...
] = ...
fsc_from_a_k_p_0( ...
 parameter_fsc ...
,n_k_all ...
,n_k_p_r ...
,n_k_all_csum_ ...
,weight_3d_k_all_ ...
,c_k_p_base_ ...
,d_k_p_reco_ ...
);
fsc_crop_reco_kx__(:,1+ncrop) = fsc_crop_reco_k_;
end;%for ncrop=0:n_crop-1;
%%%%;
save( ...
 fname_mat ...
,'fsc_reco_k_' ...
,'corr_reco' ...
,'fsc_crop_reco_kx__' ...
,'corr_crop_reco_x_' ...
,'k_Ainv_p_r_' ...
,'k_Ainv_p_r_max' ...
,'kinv_A_p_r_' ...
);
%%%%%%%%;
close_fname_tmp(fname_pre);
end;%if ( flag_recalc | ~fname_skip );
%%%%%%%%;
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
