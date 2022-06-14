function ...
[ ...
 parameter ...
,ampm_ ...
] = ...
ampm_fromdisk_2( ...
 parameter ...
,dir_pm ...
,str_filter ...
,dir_relion ...
);

str_thisfunction = 'ampm_fromdisk_2';

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); dir_pm=[]; end; na=na+1;
if (nargin<1+na); str_filter=[]; end; na=na+1;
if (nargin<1+na); dir_relion=[]; end; na=na+1;
if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
if ~isfield(parameter,'flag_center_image'); parameter.flag_center_image= ~isempty(strfind(dir_pm,'ps')) | ~isempty(strfind(dir_pm,'LSUbl17dep')); end;
if ~isfield(parameter,'flag_store_S_k_p__'); parameter.flag_store_S_k_p__=1; end;
if ~isfield(parameter,'flag_store_M_k_p__'); parameter.flag_store_M_k_p__=1; end;
flag_verbose = parameter.flag_verbose;
flag_center_image = parameter.flag_center_image;
flag_store_S_k_p__ = parameter.flag_store_S_k_p__;
flag_store_M_k_p__ = parameter.flag_store_M_k_p__;
if isempty(dir_pm); dir_pm = pwd; end;
if isempty(str_filter); str_filter = 'X_2d_xcor_d0_a1t*'; end;
if isempty(dir_relion); 
dir_relion = dir_pm; tmp_ij = max(strfind(dir_pm,'dir_pm')); 
if ~isempty(tmp_ij); dir_relion = sprintf('%s/dir_relion',dir_pm(1:tmp_ij-2)); end;
end;%if isempty(dir_relion); 
%%%%%%%%;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); string_root = 'home'; end;
if (strcmp(platform,'eval1')); string_root = 'home'; end;
if (strcmp(platform,'rusty')); string_root = 'mnt/home'; end;
%%%%%%%%;
dir_base = sprintf('/%s/rangan/dir_cryoem',string_root);

dir_pm_mat = sprintf('%s_mat',dir_pm);
ampm_ = struct('type','ampm');
ampm_.str_filter = str_filter;
ampm_.dir_pm = dir_pm;
ampm_.dir_pm_mat = dir_pm_mat;
dir_relion_mat = sprintf('%s_mat',dir_relion);
ampm_.dir_relion = dir_relion;
ampm_.dir_relion_mat = dir_relion_mat;
%%%%%%%%;
fname_X_relion_mat = sprintf('%s/X_relion_.mat',dir_relion_mat);
if (~exist(fname_X_relion_mat,'file'));
fname_X_relion_mat = sprintf('%s/job_1024/X_relion_.mat',dir_relion_mat);
end;%if (~exist(fname_X_relion_mat,'file'));
if (~exist(fname_X_relion_mat,'file'));
if (flag_verbose); disp(sprintf(' %% %s not found',fname_X_relion_mat)); end;
end;%if (~exist(fname_X_relion_mat,'file'));
ampm_.X_relion_ = [0];
if ( exist(fname_X_relion_mat,'file'));
tmp_ = load(fname_X_relion_mat);
if ( isfield(tmp_,'X_relion_')); ampm_.X_relion_ = tmp_.X_relion_; end;
clear tmp_;
end;%if ( exist(fname_X_relion_mat,'file'));
ampm_.X_relion = ampm_.X_relion_(end);
if (flag_verbose); disp(sprintf(' %% ampm_.X_relion %0.2f',ampm_.X_relion)); end;
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% load the mat-files. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
str_a_k_Y_reco_from_x__ = 'a_k_Y_reco_from_M__.mat';
if flag_center_image; str_a_k_Y_reco_from_x__ = 'a_k_Y_reco_from_N__.mat'; end;
str_a_k_Y_reco_stab_alig_from_x__ = 'a_k_Y_reco_stab_alig_from_M__.mat';
if flag_center_image; str_a_k_Y_reco_stab_alig_from_x__ = 'a_k_Y_reco_stab_alig_from_N__.mat'; end;
str_fsc_crop_reco_stab_alig_from_x__ = 'fsc_crop_reco_stab_alig_from_M__.mat';
if flag_center_image; str_fsc_crop_reco_stab_alig_from_x__ = 'fsc_crop_reco_stab_alig_from_N__.mat'; end;
tmp_fname_mat__ = { ...
,{'a_x_u_pack_.mat','a_x_u_base_.mat'} ...
,{'a_k_p_quad_.mat'} ...
,{'a_k_Y_quad_.mat'} ...
,{'S_k_p__.mat'} ...
,{'M_k_p__.mat'} ...
,{'CTF_k_p__.mat','CTF_k_p_wkC__.mat'} ...
,{'a_CTF_avg_UX_Y_reco__.mat',str_a_k_Y_reco_from_x__} ...
,{str_a_k_Y_reco_stab_alig_from_x__} ...
,{str_fsc_crop_reco_stab_alig_from_x__} ...
,{'X_TM_.mat'} ...
,{'a_k_Y_0lsq_reco_.mat'} ...
,{'a_k_Y_0qbp_reco_.mat'} ...
};
for nl0=0:numel(tmp_fname_mat__)-1;
tmp_fname_mat_ = tmp_fname_mat__{1+nl0};
%%%%;
tmp_flag_exist = 0; tmp_ = {};
for nl1=0:numel(tmp_fname_mat_)-1;
tmp_fname_mat = sprintf('%s/%s',dir_pm_mat,tmp_fname_mat_{1+nl1});
if (~tmp_flag_exist &  exist(tmp_fname_mat,'file'));
tmp_flag_exist = 1; tmp_ = load(tmp_fname_mat);
end;%if (~tmp_flag_exist &  exist(tmp_fname_mat,'file'));
end;%for nl1=0:numel(tmp_fname_mat_)-1;
%%%%;
if (~tmp_flag_exist); disp(sprintf(' %% Warning, could not find %s in %s',tmp_fname_mat,str_thisfunction)); end;%if (~tmp_flag_exist);
%%%%;
if ( tmp_flag_exist & ~isempty(tmp_));
tmp_fieldname_ = fieldnames(tmp_);
for nf=0:numel(tmp_fieldname_)-1; 
tmp_fieldname = tmp_fieldname_{1+nf};
ampm_ = setfield(ampm_,tmp_fieldname,getfield(tmp_,tmp_fieldname));
end;%for nf=0:numel(tmp_fieldname_)-1; 
clear tmp_; tmp_ = {};
end;%if ( tmp_flag_exist & ~isempty(tmp_));
%%%%;
end;%for nl0=0:numel(tmp_fname_mat__)-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% clean up discrepancies in naming convention. ;
%%%%%%%%;
if isfield(ampm_,'X_u_0_'); ampm_.x_u_0___ = ampm_.X_u_0_; end;
if isfield(ampm_,'X_u_1_'); ampm_.x_u_1___ = ampm_.X_u_1_; end;
if isfield(ampm_,'X_u_2_'); ampm_.x_u_2___ = ampm_.X_u_2_; end;
if isfield(ampm_,'n_X_u'); ampm_.n_xxx_u = ampm_.n_X_u; end;
if isfield(ampm_,'X_u_weight_'); ampm_.xxx_u_weight_ = ampm_.X_u_weight_; end;
if isfield(ampm_,'a_x_u_pack_'); ampm_.a_x_u_base_ = ampm_.a_x_u_pack_; end;
%%%%%%%%;
if isfield(ampm_,'CTF_index_'); ampm_.index_nCTF_from_nM_ = ampm_.CTF_index_; end;
if isfield(ampm_,'CTF_k_p__'); ampm_.CTF_k_p_kwC__ = ampm_.CTF_k_p__; end;
if isfield(ampm_,'CTF_k_p_r__'); ampm_.CTF_k_p_kC__ = ampm_.CTF_k_p_r__; end;
%%%%%%%%;
if ~flag_store_S_k_p__;
if isfield(ampm_,'S_k_p__'); ampm_.S_k_p__ = []; end;
if isfield(ampm_,'template_k_c_0__'); ampm_.template_k_c_0__ = []; end;
if isfield(ampm_,'template_k_c_1__'); ampm_.template_k_c_1__ = []; end;
if isfield(ampm_,'template_k_c_2__'); ampm_.template_k_c_2__ = []; end;
end;%if ~flag_store_S_k_p__;
if ~flag_store_M_k_p__;
if isfield(ampm_,'M_k_p__'); ampm_.M_k_p__ = []; end;
if isfield(ampm_,'N_k_p__'); ampm_.N_k_p__ = []; end;
end;%if ~flag_store_M_k_p__;

%%%%%%%%;
% load final iterations of the various ampm runs. ;
%%%%%%%%;
tmp_str_filter_ = ls(sprintf('%s/%s.mat',dir_pm_mat,str_filter));
tmp_index_start_ = strfind(tmp_str_filter_,sprintf('/%s',string_root))-1;
tmp_index_final_ = strfind(tmp_str_filter_,'.mat')+4-1;
n_str_filter = min(numel(tmp_index_start_),numel(tmp_index_final_));
tmp_str_filter__ = cell(n_str_filter,1);
%%%%%%%%;
% exclude secondary analyses. ;
%%%%%%%%;
str_exclude_ = {'_a_k_Y_.mat','_align_a_k_Y_.mat','_compare_image_rank.mat'}; n_str_exclude = numel(str_exclude_);
%%%%%%%%;
% also only include files of the form *t????[n,p]??r?.mat. ;
%%%%%%%%;
nb=0;
for na=0:n_str_filter-1;
tmp_fname_mat = tmp_str_filter_(1+tmp_index_start_(1+na)+0:1+tmp_index_final_(1+na)-1);
if (tmp_fname_mat(1)==sprintf('\n')); tmp_fname_mat = tmp_fname_mat(2:end); end;
if  exist(tmp_fname_mat,'file');
tmp_f = @(s) ~isempty(strfind(tmp_fname_mat,s));
flag_exclude = sum(cellfun(tmp_f,str_exclude_));
flag_include = ...
  (tmp_fname_mat(end-3-2)=='r') ...
& ((tmp_fname_mat(end-3-2-3)=='n') | (tmp_fname_mat(end-3-2-3)=='p')) ...
& (tmp_fname_mat(end-3-2-3-5)=='t') ...
  ;
if ~flag_exclude &  flag_include; tmp_str_filter__{1+nb} = tmp_fname_mat; nb=nb+1; end;
end;%if  exist(tmp_fname_mat,'file');
end;%for na=0:n_str_filter-1;
n_str_filter = nb;
tmp_str_filter__ = tmp_str_filter__(1:n_str_filter);
%%%%%%%%;
ampm_.n_str_filter = n_str_filter;
ampm_.str_fname_mat_a_ = cell(n_str_filter,1);
ampm_.str_fname_a_k_Y_mat_a_ = cell(n_str_filter,1);
ampm_.str_fname_compare_image_rank_mat_a_ = cell(n_str_filter,1);
ampm_.a_k_Y_ampm_yka__ = zeros(ampm_.n_lm_sum,n_str_filter);
ampm_.euler_polar_a_ampm_Ma__ = zeros(ampm_.n_M,n_str_filter);
ampm_.euler_azimu_b_ampm_Ma__ = zeros(ampm_.n_M,n_str_filter);
ampm_.euler_gamma_z_ampm_Ma__ = zeros(ampm_.n_M,n_str_filter);
ampm_.image_delta_x_ampm_Ma__ = zeros(ampm_.n_M,n_str_filter);
ampm_.image_delta_y_ampm_Ma__ = zeros(ampm_.n_M,n_str_filter);
ampm_.image_I_value_ampm_Ma__ = zeros(ampm_.n_M,n_str_filter);
ampm_.image_X_value_ampm_Ma__ = zeros(ampm_.n_M,n_str_filter);
ampm_.image_S_index_ampm_Ma__ = zeros(ampm_.n_M,n_str_filter);
ampm_.image_R_value_ampm_Ma__ = zeros(ampm_.n_M,n_str_filter);
ampm_.X_best_ampm_ia__ = zeros(0,0);
ampm_.X_best_flag_flip_ampm_ia__ = zeros(0,0);
ampm_.polar_a_best_ampm_ia__ = zeros(0,0);
ampm_.azimu_b_best_ampm_ia__ = zeros(0,0);
ampm_.gamma_z_best_ampm_ia__ = zeros(0,0);
ampm_.delta_best_ampm_dia__ = zeros(0,0,0);
%%%%;
for na=0:n_str_filter-1;
%%%%%%%%;
tmp_fname_mat = tmp_str_filter__{1+na};
ampm_.str_fname_mat_a_{1+na} = tmp_fname_mat;
if ( exist(tmp_fname_mat,'file'));
tmp_fname_mat;
tmp_ = load(tmp_fname_mat);
%ampm_.a_k_Y_ampm_yka__(:,1+na) = tmp_.a_k_Y_reco_yki__(:,end); %<-- should be empty. ;
if ( isfield(tmp_,'euler_polar_a_Mi__')); ampm_.euler_polar_a_ampm_Ma__(:,1+na) = tmp_.euler_polar_a_Mi__(:,end); end;
if ( isfield(tmp_,'euler_azimu_b_Mi__')); ampm_.euler_azimu_b_ampm_Ma__(:,1+na) = tmp_.euler_azimu_b_Mi__(:,end); end;
if ( isfield(tmp_,'euler_gamma_z_Mi__')); ampm_.euler_gamma_z_ampm_Ma__(:,1+na) = tmp_.euler_gamma_z_Mi__(:,end); end;
if ( isfield(tmp_,'image_delta_x_acc_Mi__')); ampm_.image_delta_x_ampm_Ma__(:,1+na) = tmp_.image_delta_x_acc_Mi__(:,end) + tmp_.image_delta_x_upd_Mi__(:,end); end;
if ( isfield(tmp_,'image_delta_y_acc_Mi__')); ampm_.image_delta_y_ampm_Ma__(:,1+na) = tmp_.image_delta_y_acc_Mi__(:,end) + tmp_.image_delta_y_upd_Mi__(:,end); end;
if ( isfield(tmp_,'image_I_value_Mi__')); ampm_.image_I_value_ampm_Ma__(:,1+na) = tmp_.image_I_value_Mi__(:,end); end;
if ( isfield(tmp_,'image_X_value_Mi__')); ampm_.image_X_value_ampm_Ma__(:,1+na) = tmp_.image_X_value_Mi__(:,end-1); end;
if ( isfield(tmp_,'image_S_index_Mi__')); ampm_.image_S_index_ampm_Ma__(:,1+na) = tmp_.image_S_index_Mi__(:,end); end;
clear tmp_;
end;%if ( exist(tmp_fname_mat,'file'));
%%%%%%%%;
tmp_ij = strfind(tmp_fname_mat,'.mat');
tmp_fname_a_k_Y_mat = sprintf('%s_a_k_Y_.mat',tmp_fname_mat([1:tmp_ij-1]));
ampm_.str_fname_mat_a_k_Y_a_{1+na} = tmp_fname_a_k_Y_mat;
if ( exist(tmp_fname_a_k_Y_mat,'file'));
tmp_ = load(tmp_fname_a_k_Y_mat);
if ( isfield(tmp_,'a_k_Y_reco_')); ampm_.a_k_Y_ampm_yka__(:,1+na) = tmp_.a_k_Y_reco_; end;
clear tmp_;
end;%if ( exist(tmp_fname_a_k_Y_mat,'file'));
%%%%%%%%;
tmp_ij = strfind(tmp_fname_mat,'.mat');
tmp_fname_align_a_k_Y_mat = sprintf('%s_align_a_k_Y_.mat',tmp_fname_mat([1:tmp_ij-1]));
ampm_.str_fname_mat_align_a_k_Y_a_{1+na} = tmp_fname_align_a_k_Y_mat;
if ( exist(tmp_fname_align_a_k_Y_mat,'file'));
tmp_ = load(tmp_fname_align_a_k_Y_mat);
if ( isfield(tmp_,'X_best_')); ampm_.X_best_ampm_ia__(:,1+na) = tmp_.X_best_; end;
if ( isfield(tmp_,'X_best_flag_flip_')); ampm_.X_best_flag_flip_ampm_ia__(:,1+na) = tmp_.X_best_flag_flip_; end;
if ( isfield(tmp_,'polar_a_best_')); ampm_.polar_a_best_ampm_ia__(:,1+na) = tmp_.polar_a_best_; end;
if ( isfield(tmp_,'azimu_b_best_')); ampm_.azimu_b_best_ampm_ia__(:,1+na) = tmp_.azimu_b_best_; end;
if ( isfield(tmp_,'gamma_z_best_')); ampm_.gamma_z_best_ampm_ia__(:,1+na) = tmp_.gamma_z_best_; end;
if ( isfield(tmp_,'delta_best__')); ampm_.delta_best_ampm_dia__(:,:,1+na) = tmp_.delta_best__; end;
clear tmp_;
end;%if ( exist(tmp_fname_align_a_k_Y_mat,'file'));
%%%%%%%%;
tmp_ij = strfind(tmp_fname_mat,'.mat');
tmp_fname_compare_image_rank_mat = sprintf('%s_compare_image_rank.mat',tmp_fname_mat([1:tmp_ij-1]));
ampm_.str_fname_mat_compare_image_rank_a_{1+na} = tmp_fname_compare_image_rank_mat;
if ( exist(tmp_fname_compare_image_rank_mat,'file'));
tmp_ = load(tmp_fname_compare_image_rank_mat);
if ( isfield(tmp_,'tmp_image_X_value_')); ampm_.image_X_value_ampm_Ma__(:,1+na) = tmp_.tmp_image_X_value_; end;
if ( isfield(tmp_,'tmp_R_k_p_l2_')); ampm_.image_R_value_ampm_Ma__(:,1+na) = tmp_.tmp_R_k_p_l2_; end;
clear tmp_;
end;%if ( exist(tmp_fname_compare_image_rank_mat,'file'));
%%%%%%%%;
tmp_ij = strfind(tmp_fname_mat,'.mat');
tmp_fname_align_crop_a_k_Y_mat = sprintf('%s_align_crop_a_k_Y_.mat',tmp_fname_mat([1:tmp_ij-1]));
ampm_.str_fname_mat_align_crop_a_k_Y_a_{1+na} = tmp_fname_align_crop_a_k_Y_mat;
if ( exist(tmp_fname_align_crop_a_k_Y_mat,'file'));
tmp_ = load(tmp_fname_align_crop_a_k_Y_mat);
if ( isfield(tmp_,'corr_crop_')); ampm_.corr_crop_ampm_xa__(:,1+na) = tmp_.corr_crop_; end;
clear tmp_;
end;%if ( exist(tmp_fname_align_crop_a_k_Y_mat,'file'));
%%%%%%%%;
tmp_ij = strfind(tmp_fname_mat,'.mat');
tmp_fname_fsc_crop_ampm_kx__mat = sprintf('%s_fsc_crop_ampm_kx__.mat',tmp_fname_mat([1:tmp_ij-1]));
ampm_.str_fname_mat_fsc_crop_ampm_kx__a_{1+na} = tmp_fname_fsc_crop_ampm_kx__mat;
if ( exist(tmp_fname_fsc_crop_ampm_kx__mat,'file'));
tmp_ = load(tmp_fname_fsc_crop_ampm_kx__mat);

if ( isfield(tmp_,'corr_base_vs_ampm')); ampm_.corr_base_vs_ampm_a_(1+na) = tmp_.corr_base_vs_ampm; end;
if ( isfield(tmp_,'corr_reco_vs_ampm')); ampm_.corr_reco_vs_ampm_a_(1+na) = tmp_.corr_reco_vs_ampm; end;
if ( isfield(tmp_,'corr_full_base_vs_crop_ampm_x_')); ampm_.corr_full_base_vs_crop_ampm_xa__(:,1+na) = tmp_.corr_full_base_vs_crop_ampm_x_; end;
if ( isfield(tmp_,'corr_crop_base_vs_crop_ampm_x_')); ampm_.corr_crop_base_vs_crop_ampm_xa__(:,1+na) = tmp_.corr_crop_base_vs_crop_ampm_x_; end;
if ( isfield(tmp_,'corr_full_reco_vs_crop_ampm_x_')); ampm_.corr_full_reco_vs_crop_ampm_xa__(:,1+na) = tmp_.corr_full_reco_vs_crop_ampm_x_; end;
if ( isfield(tmp_,'corr_crop_reco_vs_crop_ampm_x_')); ampm_.corr_crop_reco_vs_crop_ampm_xa__(:,1+na) = tmp_.corr_crop_reco_vs_crop_ampm_x_; end;


if ( isfield(tmp_,'fsc_crop_ampm_kx__')); ampm_.fsc_crop_ampm_kxa___(:,:,1+na) = tmp_.fsc_crop_ampm_kx__; end;
if ( isfield(tmp_,'fsc_ampm_k_')); ampm_.fsc_ampm_ka__(:,1+na) = tmp_.fsc_ampm_k_; end;
clear tmp_;
end;%if ( exist(tmp_fname_fsc_crop_ampm_kx__mat,'file'));
%%%%%%%%;
end;%for na=0:n_str_filter-1;

ampm_.t_sval_a_ = zeros(ampm_.n_str_filter,1);
ampm_.p_vs_n_a_ = zeros(ampm_.n_str_filter,1);
ampm_.n_sval_a_ = zeros(ampm_.n_str_filter,1);
ampm_.p_sval_a_ = zeros(ampm_.n_str_filter,1);
ampm_.r_sval_a_ = zeros(ampm_.n_str_filter,1);
for nstr_filter=0:ampm_.n_str_filter-1;
str_fname_mat = ampm_.str_fname_mat_a_{1+nstr_filter};
tmp_ij = strfind(str_fname_mat,'.mat');
tmp_ij = tmp_ij-2; assert(str_fname_mat(tmp_ij)=='r');
r_sval = str2num(str_fname_mat(tmp_ij+1));
tmp_ij = tmp_ij-3; assert((str_fname_mat(tmp_ij)=='n') | (str_fname_mat(tmp_ij)=='p'));
p_vs_n = (str_fname_mat(tmp_ij)=='p');
if (p_vs_n==1); p_sval = str2num(str_fname_mat(tmp_ij+[1:2])); n_sval = 0; end;
if (p_vs_n==0); n_sval = str2num(str_fname_mat(tmp_ij+[1:2])); p_sval = 0; end;
tmp_ij = tmp_ij-5; assert(str_fname_mat(tmp_ij)=='t');
t_sval = str2num(str_fname_mat(tmp_ij+[1:4]));
ampm_.t_sval_a_(1+nstr_filter) = t_sval;
ampm_.p_vs_n_a_(1+nstr_filter) = p_vs_n;
ampm_.n_sval_a_(1+nstr_filter) = n_sval;
ampm_.p_sval_a_(1+nstr_filter) = p_sval;
ampm_.r_sval_a_(1+nstr_filter) = r_sval;
end;%for nstr_filter=0:ampm_.n_str_filter-1;
