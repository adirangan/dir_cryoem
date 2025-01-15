%{
if ~exist('wS_from_single_shell_sba__','var'); wS_from_single_shell_sba__=[]; end;
if ~exist('dwSda_from_single_shell_sba__','var'); dwSda_from_single_shell_sba__=[]; end;
if ~exist('dwSdb_from_single_shell_sba__','var'); dwSdb_from_single_shell_sba__=[]; end;
if ~exist('ddwSdaa_from_single_shell_sba__','var'); ddwSdaa_from_single_shell_sba__=[]; end;
if ~exist('ddwSdab_from_single_shell_sba__','var'); ddwSdab_from_single_shell_sba__=[]; end;
if ~exist('ddwSdbb_from_single_shell_sba__','var'); ddwSdbb_from_single_shell_sba__=[]; end;
if ~exist('R_use__','var'); R_use__ = []; end;
if ~exist('ba_from_single_shell_baba__','var'); ba_from_single_shell_baba__=[]; end;
if ~exist('a_R_k_p_qk_','var'); a_R_k_p_qk_=[]; end;
if ~exist('dvol_a_R_k_p_qk_','var'); dvol_a_R_k_p_qk_=[]; end;
if ~exist('ba_from_single_shell_baba__','var'); ba_from_single_shell_baba__=[]; end;
if ~exist('wS_from_R_single_shell_sba__','var'); wS_from_R_single_shell_sba__=[]; end;
if ~exist('dwSda_from_R_single_shell_sba__','var'); dwSda_from_R_single_shell_sba__=[]; end;
if ~exist('dwSdb_from_R_single_shell_sba__','var'); dwSdb_from_R_single_shell_sba__=[]; end;
if ~exist('ddwSdaa_from_R_single_shell_sba__','var'); ddwSdaa_from_R_single_shell_sba__=[]; end;
if ~exist('ddwSdab_from_R_single_shell_sba__','var'); ddwSdab_from_R_single_shell_sba__=[]; end;
if ~exist('ddwSdbb_from_R_single_shell_sba__','var'); ddwSdbb_from_R_single_shell_sba__=[]; end;
%}
function ...
[ ...
 parameter ...
,Hvt_qkabc_ ...
,Hv_q3d_k_p_qk_ ...
,Ht_q2d_M3__ ...
,a_restore_C2M0_k_p_qk__ ...
,Hvv_q3d_k_p_qk_ ...
,Hvt_q3d_k_p_qk_ ...
,Htv_q2d_M3__ ...
,Htt_q2d_M3__ ...
,dtau_euler_polar_a_M_ ...
,dtau_euler_azimu_b_M_ ...
,dtau_euler_gamma_z_M_ ...
,n_dvt ... 
,dvt_ ... 
,dvt ... 
,ssnll_tmp_q2d_dvt_ ... 
,dssnll_mid_q2d ... 
,dssnll_dif_q2d ... 
,dssnll_lsq_q2d ... 
,ddssnll_mid_q2d ... 
,ddssnll_dif_q2d ... 
,ddssnll_lsq_q2d ... 
,wS_from_single_shell_sba__ ...
,dwSda_from_single_shell_sba__ ...
,dwSdb_from_single_shell_sba__ ...
,ddwSdaa_from_single_shell_sba__ ...
,ddwSdab_from_single_shell_sba__ ...
,ddwSdbb_from_single_shell_sba__ ...
,R_use__ ...
,a_R_k_p_qk_ ...
,dvol_a_R_k_p_qk_ ...
,ba_from_single_shell_baba__ ...
,wS_from_R_single_shell_sba__ ...
,dwSda_from_R_single_shell_sba__ ...
,dwSdb_from_R_single_shell_sba__ ...
,ddwSdaa_from_R_single_shell_sba__ ...
,ddwSdab_from_R_single_shell_sba__ ...
,ddwSdbb_from_R_single_shell_sba__ ...
] = ...
ddssnll_from_a_k_p_4( ...
 parameter ...
,n_qk ...
,n_qk_csum_ ...
,k_p_r_qk_ ...
,k_p_azimu_b_qk_ ...
,k_p_polar_a_qk_ ...
,weight_3d_k_p_qk_ ...
,weight_shell_qk_ ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,k_c_0_qk_ ...
,k_c_1_qk_ ...
,k_c_2_qk_ ...
,n_polar_a_k_ ...
,polar_a_ka__ ...
,n_azimu_b_ka__ ...
,qref_k_eq_d ...
,a_k_p_qk_ ...
,dvol_a_k_p_qk_ ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_wk_ ...
,n_viewing_S ...
,S_k_p_q2d_wkS__ ...
,viewing_polar_a_S_ ...
,viewing_azimu_b_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,viewing_gamma_z_S_ ...
,n_M ...
,weight_imagecount_M_ ...
,M_k_p_wkM__ ...
,n_CTF ...
,index_nCTF_from_nM_ ...
,CTF_k_p_r_kC__ ...
,CTF_k_p_wkC__ ...
,n_eta ...
,index_neta_from_nM_ ...
,eta_k_p_r_ke__ ...
,eta_k_p_wke__ ...
,euler_polar_a_M_ ...
,euler_azimu_b_M_ ...
,euler_gamma_z_M_ ...
,dtau_euler_polar_a_M_ ...
,dtau_euler_azimu_b_M_ ...
,dtau_euler_gamma_z_M_ ...
,KAPPA ...
,wS_from_single_shell_sba__ ...
,dwSda_from_single_shell_sba__ ...
,dwSdb_from_single_shell_sba__ ...
,ddwSdaa_from_single_shell_sba__ ...
,ddwSdab_from_single_shell_sba__ ...
,ddwSdbb_from_single_shell_sba__ ...
,R_use__ ...
,a_R_k_p_qk_ ...
,dvol_a_R_k_p_qk_ ...
,ba_from_single_shell_baba__ ...
,wS_from_R_single_shell_sba__ ...
,dwSda_from_R_single_shell_sba__ ...
,dwSdb_from_R_single_shell_sba__ ...
,ddwSdaa_from_R_single_shell_sba__ ...
,ddwSdab_from_R_single_shell_sba__ ...
,ddwSdbb_from_R_single_shell_sba__ ...
);

str_thisfunction = 'ddssnll_from_a_k_p_4';

%%%%%%%%;
if (nargin<1);
%%%%%%%%;
test_slice_vs_volume_integral_from_a_k_p_6;
%%%%%%%%;
disp(sprintf(' %% returning')); return;
%%%%%%%%;
end;%if (nargin<1);
%%%%%%%%;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_qk=[]; end; na=na+1;
if (nargin<1+na); n_qk_csum_=[]; end; na=na+1;
if (nargin<1+na); k_p_r_qk_=[]; end; na=na+1;
if (nargin<1+na); k_p_azimu_b_qk_=[]; end; na=na+1;
if (nargin<1+na); k_p_polar_a_qk_=[]; end; na=na+1;
if (nargin<1+na); weight_3d_k_p_qk_=[]; end; na=na+1;
if (nargin<1+na); weight_shell_qk_=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); k_p_r_=[]; end; na=na+1;
if (nargin<1+na); k_p_r_max=[]; end; na=na+1;
if (nargin<1+na); weight_3d_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); k_c_0_qk_=[]; end; na=na+1;
if (nargin<1+na); k_c_1_qk_=[]; end; na=na+1;
if (nargin<1+na); k_c_2_qk_=[]; end; na=na+1;
if (nargin<1+na); n_polar_a_k_=[]; end; na=na+1;
if (nargin<1+na); polar_a_ka__=[]; end; na=na+1;
if (nargin<1+na); n_azimu_b_ka__=[]; end; na=na+1;
if (nargin<1+na); qref_k_eq_d=[]; end; na=na+1;
if (nargin<1+na); a_k_p_qk_=[]; end; na=na+1;
if (nargin<1+na); dvol_a_k_p_qk_=[]; end; na=na+1;
if (nargin<1+na); n_w_=[]; end; na=na+1;
if (nargin<1+na); weight_2d_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); weight_2d_wk_=[]; end; na=na+1;
if (nargin<1+na); n_viewing_S=[]; end; na=na+1;
if (nargin<1+na); S_k_p_q2d_wkS__=[]; end; na=na+1;
if (nargin<1+na); viewing_polar_a_S_=[]; end; na=na+1;
if (nargin<1+na); viewing_azimu_b_S_=[]; end; na=na+1;
if (nargin<1+na); viewing_weight_S_=[]; end; na=na+1;
if (nargin<1+na); n_viewing_polar_a=[]; end; na=na+1;
if (nargin<1+na); viewing_polar_a_=[]; end; na=na+1;
if (nargin<1+na); n_viewing_azimu_b_=[]; end; na=na+1;
if (nargin<1+na); viewing_gamma_z_S_=[]; end; na=na+1;
if (nargin<1+na); n_M=[]; end; na=na+1;
if (nargin<1+na); weight_imagecount_M_=[]; end; na=na+1;
if (nargin<1+na); M_k_p_wkM__=[]; end; na=na+1;
if (nargin<1+na); n_CTF=[]; end; na=na+1;
if (nargin<1+na); index_nCTF_from_nM_=[]; end; na=na+1;
if (nargin<1+na); CTF_k_p_r_kC__=[]; end; na=na+1;
if (nargin<1+na); CTF_k_p_wkC__=[]; end; na=na+1;
if (nargin<1+na); n_eta=[]; end; na=na+1;
if (nargin<1+na); index_neta_from_nM_=[]; end; na=na+1;
if (nargin<1+na); eta_k_p_r_ke__=[]; end; na=na+1;
if (nargin<1+na); eta_k_p_wke__=[]; end; na=na+1;
if (nargin<1+na); euler_polar_a_M_=[]; end; na=na+1;
if (nargin<1+na); euler_azimu_b_M_=[]; end; na=na+1;
if (nargin<1+na); euler_gamma_z_M_=[]; end; na=na+1;
if (nargin<1+na); dtau_euler_polar_a_M_=[]; end; na=na+1;
if (nargin<1+na); dtau_euler_azimu_b_M_=[]; end; na=na+1;
if (nargin<1+na); dtau_euler_gamma_z_M_=[]; end; na=na+1;
if (nargin<1+na); KAPPA=[]; end; na=na+1;
if (nargin<1+na); wS_from_single_shell_sba__=[]; end; na=na+1;
if (nargin<1+na); dwSda_from_single_shell_sba__=[]; end; na=na+1;
if (nargin<1+na); dwSdb_from_single_shell_sba__=[]; end; na=na+1;
if (nargin<1+na); ddwSdaa_from_single_shell_sba__=[]; end; na=na+1;
if (nargin<1+na); ddwSdab_from_single_shell_sba__=[]; end; na=na+1;
if (nargin<1+na); ddwSdbb_from_single_shell_sba__=[]; end; na=na+1;
if (nargin<1+na); R_use__=[]; end; na=na+1;
if (nargin<1+na); a_R_k_p_qk_=[]; end; na=na+1;
if (nargin<1+na); dvol_a_R_k_p_qk_=[]; end; na=na+1;
if (nargin<1+na); ba_from_single_shell_baba__=[]; end; na=na+1;
if (nargin<1+na); wS_from_R_single_shell_sba__=[]; end; na=na+1;
if (nargin<1+na); dwSda_from_R_single_shell_sba__=[]; end; na=na+1;
if (nargin<1+na); dwSdb_from_R_single_shell_sba__=[]; end; na=na+1;
if (nargin<1+na); ddwSdaa_from_R_single_shell_sba__=[]; end; na=na+1;
if (nargin<1+na); ddwSdab_from_R_single_shell_sba__=[]; end; na=na+1;
if (nargin<1+na); ddwSdbb_from_R_single_shell_sba__=[]; end; na=na+1;

if isempty(parameter);
parameter = struct('type','parameter');
end;%if isempty(parameter);
%%%%%%%%;
if (~isfield(parameter,'flag_verbose')); parameter.flag_verbose = 0; end; %<-- parameter_bookmark. ;
flag_verbose = parameter.flag_verbose;
if (~isfield(parameter,'flag_check')); parameter.flag_check = 0; end; %<-- parameter_bookmark. ;
flag_check = parameter.flag_check;
if (~isfield(parameter,'dvt_check')); parameter.dvt_check = 1e-3; end; %<-- parameter_bookmark. ;
dvt_check = parameter.dvt_check;
if (~isfield(parameter,'flag_disp')); parameter.flag_disp = 0; end; %<-- parameter_bookmark. ;
flag_disp = parameter.flag_disp; nf=0;
if (~isfield(parameter,'tolerance_master')); parameter.tolerance_master = 1e-2; end; %<-- parameter_bookmark. ;
tolerance_master = parameter.tolerance_master;
if (~isfield(parameter,'viewing_k_eq_d')); parameter.viewing_k_eq_d = []; end; %<-- parameter_bookmark. ;
viewing_k_eq_d = parameter.viewing_k_eq_d;
if (~isfield(parameter,'template_k_eq_d')); parameter.template_k_eq_d = -1; end; %<-- parameter_bookmark. ;
template_k_eq_d = parameter.template_k_eq_d;
if (~isfield(parameter,'n_order')); parameter.n_order = 5; end; %<-- parameter_bookmark. ;
n_order = parameter.n_order;
if ~isfield(parameter,'kernel_basic_qref_k_eq_d_double'); parameter.kernel_basic_qref_k_eq_d_double=[]; end;
kernel_basic_qref_k_eq_d_double=parameter.kernel_basic_qref_k_eq_d_double;
if (~isfield(parameter,'kernel_basic_l_max_use')); parameter.kernel_basic_l_max_use = round(2*pi*k_p_r_max); end; %<-- parameter_bookmark. ;
kernel_basic_l_max_use = parameter.kernel_basic_l_max_use;
if (~isfield(parameter,'kernel_basic_l_max_ext')); parameter.kernel_basic_l_max_ext = ceil(1.25*kernel_basic_l_max_use); end; %<-- parameter_bookmark. ;
kernel_basic_l_max_ext = parameter.kernel_basic_l_max_ext;
if (~isfield(parameter,'kernel_basic_l_max_band')); parameter.kernel_basic_l_max_band = +Inf; end; %<-- parameter_bookmark. ;
kernel_basic_l_max_band = parameter.kernel_basic_l_max_band;
if (~isfield(parameter,'tolerance_pinv')); parameter.tolerance_pinv = 1e-6; end; %<-- parameter_bookmark. ;
tolerance_pinv = parameter.tolerance_pinv;
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

a_restore_C2M0_k_p_qk_ = []; a_restore_C2M0_k_p_qk__ = [];
a_restore_C1M1_k_p_qk_ = []; a_restore_C1M1_k_p_qk__ = [];
a_restore_C0M2_k_p_qk_ = []; a_restore_C0M2_k_p_qk__ = [];
dtau_a_restore_C2M0_k_p_qk_ = []; dtau_a_restore_C2M0_k_p_qk__ = [];
dtau_a_restore_C1M1_k_p_qk_ = []; dtau_a_restore_C1M1_k_p_qk__ = [];
dtau_a_restore_C0M2_k_p_qk_ = []; dtau_a_restore_C0M2_k_p_qk__ = [];
dtau_dtau_a_restore_C2M0_k_p_qk_ = []; dtau_dtau_a_restore_C2M0_k_p_qk__ = [];
dtau_dtau_a_restore_C1M1_k_p_qk_ = []; dtau_dtau_a_restore_C1M1_k_p_qk__ = [];
dtau_dtau_a_restore_C0M2_k_p_qk_ = []; dtau_dtau_a_restore_C0M2_k_p_qk__ = [];
dtau_ssnll_q2d_M3__ = []; dtau_ssnll_q2d_3M__ = [];
dtau_dtau_ssnll_q2d_M33___ = []; dtau_dtau_ssnll_q2d_33M___ = [];
dtau_dvol_ssnll_q2d_M3__ = []; dtau_dvol_ssnll_q2d_3M__ = [];

if isempty(weight_imagecount_M_); weight_imagecount_M_ = ones(n_M,1); end;
flag_dtau = ~isempty(dtau_euler_polar_a_M_) | ~isempty(dtau_euler_azimu_b_M_) | ~isempty(dtau_euler_gamma_z_M_) ;
flag_dvol = ~isempty(dvol_a_k_p_qk_) ;

ddssnll_from_a_k_p_helper_n_w_4; %<-- set indices. ;
ddssnll_from_a_k_p_helper_qref_4; %<-- set indices and check for consistency. ;
ddssnll_from_a_k_p_helper_n_qk_4; %<-- set indices (some redundant definitions). ;
ddssnll_from_a_k_p_helper_weight_3d_riesz_4; %<-- set weights. ;
ddssnll_from_a_k_p_helper_reconfigure_a_4; %<-- reconfigure a. ;

if  flag_dtau &  flag_dvol;
if (flag_verbose>0); disp(sprintf(' %%  flag_dtau %d &  flag_dvol %d = %d', flag_dtau, flag_dvol, flag_dtau &  flag_dvol)); end;
ddssnll_from_a_k_p_helper_reconfigure_dvol_a_4; %<-- reconfigure dvol_a. ;
ddssnll_from_a_k_p_helper_q2d_4; %<-- determine derivatives using dvol_a. ;
ddssnll_helper_reconfigure_dtau_3; %<-- reconfigure dtau. ;
ddssnll_from_a_k_p_helper_a_restore_from_kappa_4; %<-- use kappa_basic_apply to construct a_restore. ;
% Note: a_restore_C2M0_k_p_qk_ etc. needed for Hvv_q3d. ;
end;%if  flag_dtau &  flag_dvol;

if ~flag_dtau &  flag_dvol;
if (flag_verbose>0); disp(sprintf(' %% ~flag_dtau %d &  flag_dvol %d = %d',~flag_dtau, flag_dvol,~flag_dtau &  flag_dvol)); end;
ddssnll_from_a_k_p_helper_reconfigure_dvol_a_4; %<-- reconfigure dvol_a. ;
ddssnll_helper_reconfigure_dtau_3; %<-- reconfigure dtau. ;
ddssnll_from_a_k_p_helper_q2d_4; %<-- determine derivatives using dvol_a. ;
ddssnll_helper_reconfigure_dtau_3; %<-- reconfigure dtau. ;
ddssnll_helper_dtau_firstorder_3; %<-- find critical dtau, given dvol_a. ;
ddssnll_helper_dtau_from_dtau_firstorder_3; %<-- if dtau is empty. ;
ddssnll_from_a_k_p_helper_q2d_4; %<-- determine derivatives using dvol_a. ;
ddssnll_from_a_k_p_helper_a_restore_from_kappa_4; %<-- use kappa_qpro_apply to construct a_restore. ;
% Note: a_restore_C2M0_k_p_qk_ etc. needed for Hvv_q3d. ;
end;%if ~flag_dtau &  flag_dvol;

if  flag_dtau & ~flag_dvol;
if (flag_verbose>0); disp(sprintf(' %%  flag_dtau %d & ~flag_dvol %d = %d', flag_dtau,~flag_dvol, flag_dtau & ~flag_dvol)); end;
ddssnll_helper_reconfigure_dtau_3; %<-- reconfigure dtau. ;
ddssnll_from_a_k_p_helper_a_restore_from_kappa_4; %<-- use kappa_qpro_apply to construct a_restore. ;
% Note: a_restore_C2M0_k_p_qk_ etc. needed for Hvv_q3d. ;
ddssnll_from_a_k_p_helper_convert_dvol_firstorder_from_dtau_4; %<-- find critical dvol_a, given dtau. ;
ddssnll_from_a_k_p_helper_dvol_from_dvol_firstorder_4; %<-- if dvol is empty. ;
ddssnll_from_a_k_p_helper_reconfigure_dvol_a_4; %<-- reconfigure dvol_a. ;
ddssnll_helper_reconfigure_dtau_3; %<-- reconfigure dtau. ;
ddssnll_from_a_k_p_helper_q2d_4; %<-- determine derivatives using dvol_a. ;
ddssnll_helper_reconfigure_dtau_3; %<-- reconfigure dtau. ;
end;%if  flag_dtau & ~flag_dvol;

ddssnll_from_a_k_p_helper_convert_a_times_dtau_a_restore_C2M0_4;
ddssnll_from_a_k_p_helper_convert_dvol_a_times_a_restore_C2M0_4;

%%%%%%%%;
dtau_M3__ = [ dtau_euler_polar_a_M_ , dtau_euler_azimu_b_M_ , dtau_euler_gamma_z_M_ ] ; dtau_3M__ = transpose(dtau_M3__);
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if flag_check;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

ddssnll_from_a_k_p_helper_convert_a_times_a_4;
ddssnll_from_a_k_p_helper_convert_dvol_a_times_dvol_a_4;
ddssnll_from_a_k_p_helper_convert_a_times_a_restore_C2M0_4;

ddssnll_helper_dtau_criticality_check_3;
ddssnll_from_a_k_p_helper_dvol_criticality_check_4;
ddssnll_from_a_k_p_helper_ssnll_check_4;
ddssnll_from_a_k_p_helper_dtau_ssnll_check_4;
ddssnll_from_a_k_p_helper_dtau_dtau_ssnll_check_4;
ddssnll_helper_note_hessian_3;
ddssnll_from_a_k_p_helper_Hvv_check_4;
ddssnll_from_a_k_p_helper_Htt_check_4;
ddssnll_from_a_k_p_helper_Htv_check_4;
ddssnll_from_a_k_p_helper_Hvt_check_4;
ddssnll_from_a_k_p_helper_imagesc_shell_4;
ddssnll_from_a_k_p_helper_Hv_q3d_check_4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_check;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% Now we calculate the left-hand-side: ;
%%%%%%%%;
Hvv_q3d_k_p_qk_ = dvol_a_k_p_qk_.*a_restore_C2M0_k_p_qk_;
Htv_q2d_M3__ = dtau_dvol_ssnll_q2d_M3__;
Hvt_q3d_k_p_qk_ = (a_k_p_qk_.*dtau_a_restore_C2M0_k_p_qk_ - dtau_a_restore_C1M1_k_p_qk_);
Htt_q2d_M3__ = sum(bsxfun(@times,dtau_dtau_ssnll_q2d_M33___,reshape(dtau_M3__,[n_M,1,3])),[3]);
%%%%%%%%;
Hv_q3d_k_p_qk_ = Hvv_q3d_k_p_qk_ + Hvt_q3d_k_p_qk_;
Ht_q2d_M3__ = Htt_q2d_M3__ + Htv_q2d_M3__;
%%%%%%%%;
Hvt_qkabc_ = cat(1,Hv_q3d_k_p_qk_,Ht_q2d_M3__(:,1+0),Ht_q2d_M3__(:,1+1),Ht_q2d_M3__(:,1+2));

dvt = [];
dvt_ = [];
n_dvt = [];
ssnll_tmp_q2d_dvt_ = [];
dssnll_mid_q2d = [];
dssnll_dif_q2d = [];
dssnll_lsq_q2d = [];
ddssnll_mid_q2d = [];
ddssnll_dif_q2d = [];
ddssnll_lsq_q2d = [];
ddssnll_from_a_k_p_helper_ddssnll_check_4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
