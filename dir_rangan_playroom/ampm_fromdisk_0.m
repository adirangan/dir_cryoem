function ...
[ ...
 parameter ...
,ampm_ ...
] = ...
ampm_fromdisk_0( ...
 parameter ...
,dir_pm ...
,str_filter ...
);

str_thisfunction = 'ampm_fromdisk_0';

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); dir_pm=[]; end; na=na+1;
if (nargin<1+na); str_filter=[]; end; na=na+1;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% load the mat-files. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
str_a_k_Y_reco_from_x__ = 'a_k_Y_reco_from_M__.mat'; if flag_center_image; str_a_k_Y_reco_from_x__ = 'a_k_Y_reco_from_N__.mat'; end;
tmp_fname_mat__ = { ...
,{'a_x_u_pack_.mat','a_x_u_base_.mat'} ...
,{'a_k_p_quad_.mat'} ...
,{'a_k_Y_quad_.mat'} ...
,{'S_k_p__.mat'} ...
,{'M_k_p__.mat'} ...
,{'CTF_k_p__.mat','CTF_k_p_wkC__.mat'} ...
,{'a_CTF_avg_UX_Y_reco__.mat',str_a_k_Y_reco_from_x__} ...
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
% load final iterations of the various ampm X_2d_xcor_d0_a1t runs. ;
%%%%%%%%;
tmp_X_2d_xcor_d0_a1t_ = ls(sprintf('%s/%s.mat',dir_pm_mat,str_filter));
tmp_index_start_ = strfind(tmp_X_2d_xcor_d0_a1t_,sprintf('/%s',string_root))-1;
tmp_index_final_ = strfind(tmp_X_2d_xcor_d0_a1t_,'.mat')+4-1;
n_X_2d_xcor_d0_a1t = min(numel(tmp_index_start_),numel(tmp_index_final_));
tmp_X_2d_xcor_d0_a1t__ = cell(n_X_2d_xcor_d0_a1t,1);
nb=0;
for na=0:n_X_2d_xcor_d0_a1t-1;
tmp_fname_mat = tmp_X_2d_xcor_d0_a1t_(1+tmp_index_start_(1+na)+0:1+tmp_index_final_(1+na)-1);
if (tmp_fname_mat(1)==sprintf('\n')); tmp_fname_mat = tmp_fname_mat(2:end); end;
%%%%%%%%;
% ensure format is 'X_2d_xcor_d0_a1t????p??r?.mat'. ;
%%%%%%%%;
tmp_ij = strfind(tmp_fname_mat,'.mat');
flag_include = tmp_ij>25;
if flag_include;
tmp_str = tmp_fname_mat(tmp_ij-25:tmp_ij-1);
flag_include = strcmp(tmp_str(1:16),'X_2d_xcor_d0_a1t');
end;%if flag_include;
if flag_include;
tmp_X_2d_xcor_d0_a1t__{1+nb} = tmp_fname_mat;
nb=nb+1;
end;%if flag_include;
end;%for na=0:n_X_2d_xcor_d0_a1t-1;
n_X_2d_xcor_d0_a1t = nb;
tmp_X_2d_xcor_d0_a1t__ = tmp_X_2d_xcor_d0_a1t__(1:n_X_2d_xcor_d0_a1t);
%%%%%%%%;
ampm_.n_X_2d_xcor_d0_a1t = n_X_2d_xcor_d0_a1t;
ampm_.str_fname_mat_a_ = cell(n_X_2d_xcor_d0_a1t,1);
ampm_.str_fname_a_k_Y_mat_a_ = cell(n_X_2d_xcor_d0_a1t,1);
ampm_.str_fname_compare_image_rank_mat_a_ = cell(n_X_2d_xcor_d0_a1t,1);
ampm_.a_k_Y_ampm_yka__ = zeros(ampm_.n_lm_sum,n_X_2d_xcor_d0_a1t);
ampm_.euler_polar_a_ampm_Ma__ = zeros(ampm_.n_M,n_X_2d_xcor_d0_a1t);
ampm_.euler_azimu_b_ampm_Ma__ = zeros(ampm_.n_M,n_X_2d_xcor_d0_a1t);
ampm_.euler_gamma_z_ampm_Ma__ = zeros(ampm_.n_M,n_X_2d_xcor_d0_a1t);
ampm_.image_delta_x_ampm_Ma__ = zeros(ampm_.n_M,n_X_2d_xcor_d0_a1t);
ampm_.image_delta_y_ampm_Ma__ = zeros(ampm_.n_M,n_X_2d_xcor_d0_a1t);
ampm_.image_I_value_ampm_Ma__ = zeros(ampm_.n_M,n_X_2d_xcor_d0_a1t);
ampm_.image_X_value_ampm_Ma__ = zeros(ampm_.n_M,n_X_2d_xcor_d0_a1t);
ampm_.image_S_index_ampm_Ma__ = zeros(ampm_.n_M,n_X_2d_xcor_d0_a1t);
ampm_.image_R_value_ampm_Ma__ = zeros(ampm_.n_M,n_X_2d_xcor_d0_a1t);
%%%%;
for na=0:n_X_2d_xcor_d0_a1t-1;
%%%%%%%%;
tmp_fname_mat = tmp_X_2d_xcor_d0_a1t__{1+na};
ampm_.str_fname_mat_a_{1+na} = tmp_fname_mat;
if ( exist(tmp_fname_mat,'file'));
tmp_ = load(tmp_fname_mat);
%ampm_.a_k_Y_ampm_yka__(:,1+na) = tmp_.a_k_Y_reco_yki__(:,end); %<-- should be empty. ;
ampm_.euler_polar_a_ampm_Ma__(:,1+na) = tmp_.euler_polar_a_Mi__(:,end);
ampm_.euler_azimu_b_ampm_Ma__(:,1+na) = tmp_.euler_azimu_b_Mi__(:,end);
ampm_.euler_gamma_z_ampm_Ma__(:,1+na) = tmp_.euler_gamma_z_Mi__(:,end);
ampm_.image_delta_x_ampm_Ma__(:,1+na) = tmp_.image_delta_x_acc_Mi__(:,end) + tmp_.image_delta_x_upd_Mi__(:,end);
ampm_.image_delta_y_ampm_Ma__(:,1+na) = tmp_.image_delta_y_acc_Mi__(:,end) + tmp_.image_delta_y_upd_Mi__(:,end);
ampm_.image_I_value_ampm_Ma__(:,1+na) = tmp_.image_I_value_Mi__(:,end);
ampm_.image_X_value_ampm_Ma__(:,1+na) = tmp_.image_X_value_Mi__(:,end-1);
ampm_.image_S_index_ampm_Ma__(:,1+na) = tmp_.image_S_index_Mi__(:,end);
clear tmp_;
end;%if ( exist(tmp_fname_mat,'file'));
%%%%%%%%;
tmp_ij = strfind(tmp_fname_mat,'.mat');
tmp_fname_a_k_Y_mat = sprintf('%s_a_k_Y_.mat',tmp_fname_mat([1:tmp_ij-1]));
ampm_.str_fname_mat_a_k_Y_a_{1+na} = tmp_fname_a_k_Y_mat;
if ( exist(tmp_fname_a_k_Y_mat,'file'));
tmp_ = load(tmp_fname_a_k_Y_mat);
ampm_.a_k_Y_ampm_yka__(:,1+na) = tmp_.a_k_Y_reco_;
clear tmp_;
end;%if ( exist(tmp_fname_a_k_Y_mat,'file'));
%%%%%%%%;
tmp_ij = strfind(tmp_fname_mat,'.mat');
tmp_fname_compare_image_rank_mat = sprintf('%s_compare_image_rank.mat',tmp_fname_mat([1:tmp_ij-1]));
ampm_.str_fname_mat_compare_image_rank_a_{1+na} = tmp_fname_compare_image_rank_mat;
if ( exist(tmp_fname_compare_image_rank_mat,'file'));
tmp_ = load(tmp_fname_compare_image_rank_mat);
ampm_.image_X_value_ampm_Ma__(:,1+na) = tmp_.tmp_image_X_value_; %<-- should be identical to the above. ;
ampm_.image_R_value_ampm_Ma__(:,1+na) = tmp_.tmp_R_k_p_l2_;
clear tmp_;
end;%if ( exist(tmp_fname_compare_image_rank_mat,'file'));
%%%%%%%%;
end;%for na=0:n_X_2d_xcor_d0_a1t-1;
