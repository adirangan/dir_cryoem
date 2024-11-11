function ...
[ ...
 parameter ...
,Hvt_ykabc_ ...
,Hv_q3d_k_Y_quad_yk_ ...
,Hv_q3d_k_Y_quad_yk__ ...
,Hv_q3d_k_p_quad_ ...
,Ht_q2d_M3__ ...
,a_restore_C2M0_k_Y_quad_yk_ ...
,a_restore_C2M0_k_p_quad_ ...
,Hvv_q3d_k_Y_quad_yk_ ...
,Hvt_q3d_k_Y_quad_yk_ ...
,Htv_q2d_M3__ ...
,Htt_q2d_M3__ ...
,dvol_a_k_Y_quad_yk_ ...
,dvol_a_k_Y_quad_yk__ ...
,dvol_a_k_p_quad_ ...
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
,KAPPA ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
,V_lmm___ ...
,L_lm__ ...
,d0W_betazeta_mlma____ ...
,d1W_betazeta_mlma____ ...
,d2W_betazeta_mlma____ ...
] = ...
ddssnll_3( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,l_max_ ...
,a_k_Y_quad_yk_ ...
,a_k_Y_quad_yk__ ...
,dvol_a_k_Y_quad_yk_ ...
,dvol_a_k_Y_quad_yk__ ...
,n_k_all ...
,n_k_all_csum_ ...
,k_p_r_all_ ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,weight_3d_k_all_ ...
,weight_shell_k_ ...
,weight_3d_k_p_r_ ...
,a_k_p_quad_ ...
,dvol_a_k_p_quad_ ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_wk_ ...
,n_S ...
,S_k_p_q2d_wkS__ ...
,viewing_polar_a_S_ ...
,viewing_azimu_b_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
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
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
,V_lmm___ ...
,L_lm__ ...
,d0W_betazeta_mlma____ ...
,d1W_betazeta_mlma____ ...
,d2W_betazeta_mlma____ ...
);

str_thisfunction = 'ddssnll_3';

%%%%%%%%;
if (nargin<1);
%%%%%%%%;
test_slice_vs_volume_integral_5;
%%%%%%%%;
disp(sprintf(' %% returning')); return;
%%%%%%%%;
end;%if (nargin<1);
%%%%%%%%;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); k_p_r_=[]; end; na=na+1;
if (nargin<1+na); k_p_r_max=[]; end; na=na+1;
if (nargin<1+na); l_max_=[]; end; na=na+1;
if (nargin<1+na); a_k_Y_quad_yk_=[]; end; na=na+1;
if (nargin<1+na); a_k_Y_quad_yk__=[]; end; na=na+1;
if (nargin<1+na); dvol_a_k_Y_quad_yk_=[]; end; na=na+1;
if (nargin<1+na); dvol_a_k_Y_quad_yk__=[]; end; na=na+1;
if (nargin<1+na); n_k_all=[]; end; na=na+1;
if (nargin<1+na); n_k_all_csum_=[]; end; na=na+1;
if (nargin<1+na); k_p_r_all_=[]; end; na=na+1;
if (nargin<1+na); k_p_azimu_b_all_=[]; end; na=na+1;
if (nargin<1+na); k_p_polar_a_all_=[]; end; na=na+1;
if (nargin<1+na); weight_3d_k_all_=[]; end; na=na+1;
if (nargin<1+na); weight_shell_k_=[]; end; na=na+1;
if (nargin<1+na); weight_3d_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); a_k_p_quad_=[]; end; na=na+1;
if (nargin<1+na); dvol_a_k_p_quad_=[]; end; na=na+1;
if (nargin<1+na); n_w_=[]; end; na=na+1;
if (nargin<1+na); weight_2d_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); weight_2d_wk_=[]; end; na=na+1;
if (nargin<1+na); n_S=[]; end; na=na+1;
if (nargin<1+na); S_k_p_q2d_wkS__=[]; end; na=na+1;
if (nargin<1+na); viewing_polar_a_S_=[]; end; na=na+1;
if (nargin<1+na); viewing_azimu_b_S_=[]; end; na=na+1;
if (nargin<1+na); viewing_weight_S_=[]; end; na=na+1;
if (nargin<1+na); n_viewing_polar_a=[]; end; na=na+1;
if (nargin<1+na); viewing_polar_a_=[]; end; na=na+1;
if (nargin<1+na); n_viewing_azimu_b_=[]; end; na=na+1;
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
if (nargin<1+na); Ylm_uklma___=[]; end; na=na+1;
if (nargin<1+na); k_p_azimu_b_sub_uka__=[]; end; na=na+1;
if (nargin<1+na); k_p_polar_a_sub_uka__=[]; end; na=na+1;
if (nargin<1+na); l_max_uk_=[]; end; na=na+1;
if (nargin<1+na); index_nu_n_k_per_shell_from_nk_p_r_=[]; end; na=na+1;
if (nargin<1+na); index_k_per_shell_uka__=[]; end; na=na+1;
if (nargin<1+na); V_lmm___=[]; end; na=na+1;
if (nargin<1+na); L_lm__=[]; end; na=na+1;
if (nargin<1+na); d0W_betazeta_mlma____=[]; end; na=na+1;
if (nargin<1+na); d1W_betazeta_mlma____=[]; end; na=na+1;
if (nargin<1+na); d2W_betazeta_mlma____=[]; end; na=na+1;

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
if (~isfield(parameter,'flag_kernel_qpro_d0')); parameter.flag_kernel_qpro_d0 = 1; end; %<-- parameter_bookmark. ;
flag_kernel_qpro_d0 = parameter.flag_kernel_qpro_d0;
if (~isfield(parameter,'flag_kernel_qpro_d1')); parameter.flag_kernel_qpro_d1 = 1; end; %<-- parameter_bookmark. ;
flag_kernel_qpro_d1 = parameter.flag_kernel_qpro_d1;
if (~isfield(parameter,'kernel_qpro_polar_a_pole_north')); parameter.kernel_qpro_polar_a_pole_north = 4.5*pi/24; end; %<-- parameter_bookmark. ;
kernel_qpro_polar_a_pole_north = min(pi/2,parameter.kernel_qpro_polar_a_pole_north);
parameter.kernel_qpro_polar_a_pole_north = kernel_qpro_polar_a_pole_north;
if (~isfield(parameter,'kernel_qpro_polar_a_pole_south')); parameter.kernel_qpro_polar_a_pole_south = 3.5*pi/24; end; %<-- parameter_bookmark. ;
kernel_qpro_polar_a_pole_south = min(pi/2,parameter.kernel_qpro_polar_a_pole_south);
parameter.kernel_qpro_polar_a_pole_south = kernel_qpro_polar_a_pole_south;
if ~isfield(parameter,'kernel_qpro_qref_k_eq_d_double'); parameter.kernel_qpro_qref_k_eq_d_double=0.5; end;
kernel_qpro_qref_k_eq_d_double=parameter.kernel_qpro_qref_k_eq_d_double;
if (~isfield(parameter,'kernel_qpro_l_max_use')); parameter.kernel_qpro_l_max_use = max(l_max_); end; %<-- parameter_bookmark. ;
kernel_qpro_l_max_use = parameter.kernel_qpro_l_max_use;
if (~isfield(parameter,'tolerance_pinv')); parameter.tolerance_pinv = 1e-6; end; %<-- parameter_bookmark. ;
tolerance_pinv = parameter.tolerance_pinv;
%%%%;
if ~isfield(parameter,'flag_kernel_full'); parameter.flag_kernel_full=0; end;
flag_kernel_full=parameter.flag_kernel_full;
if (kernel_qpro_polar_a_pole_north + kernel_qpro_polar_a_pole_south > pi-1e-12);
flag_kernel_full = 1;
parameter.flag_kernel_full = flag_kernel_full;
end;%if (kernel_qpro_polar_a_pole_north + kernel_qpro_polar_a_pole_south > pi-1e-12);
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

a_restore_C2M0_k_Y_quad_yk_ = []; a_restore_C2M0_k_Y_quad_yk__ = []; a_restore_C2M0_k_p_quad_ = [];
a_restore_C1M1_k_Y_quad_yk_ = []; a_restore_C1M1_k_Y_quad_yk__ = []; a_restore_C1M1_k_p_quad_ = [];
a_restore_C0M2_k_Y_quad_yk_ = []; a_restore_C0M2_k_Y_quad_yk__ = []; a_restore_C0M2_k_p_quad_ = [];
dtau_a_restore_C2M0_k_Y_quad_yk_ = []; dtau_a_restore_C2M0_k_Y_quad_yk__ = []; dtau_a_restore_C2M0_k_p_quad_ = [];
dtau_a_restore_C1M1_k_Y_quad_yk_ = []; dtau_a_restore_C1M1_k_Y_quad_yk__ = []; dtau_a_restore_C1M1_k_p_quad_ = [];
dtau_a_restore_C0M2_k_Y_quad_yk_ = []; dtau_a_restore_C0M2_k_Y_quad_yk__ = []; dtau_a_restore_C0M2_k_p_quad_ = [];
dtau_dtau_a_restore_C2M0_k_Y_quad_yk_ = []; dtau_dtau_a_restore_C2M0_k_Y_quad_yk__ = []; dtau_dtau_a_restore_C2M0_k_p_quad_ = [];
dtau_dtau_a_restore_C1M1_k_Y_quad_yk_ = []; dtau_dtau_a_restore_C1M1_k_Y_quad_yk__ = []; dtau_dtau_a_restore_C1M1_k_p_quad_ = [];
dtau_dtau_a_restore_C0M2_k_Y_quad_yk_ = []; dtau_dtau_a_restore_C0M2_k_Y_quad_yk__ = []; dtau_dtau_a_restore_C0M2_k_p_quad_ = [];
dtau_ssnll_q2d_M3__ = []; dtau_ssnll_q2d_3M__ = [];
dtau_dtau_ssnll_q2d_M33___ = []; dtau_dtau_ssnll_q2d_33M___ = [];
dtau_dvol_ssnll_q2d_M3__ = []; dtau_dvol_ssnll_q2d_3M__ = [];

if isempty(weight_imagecount_M_); weight_imagecount_M_ = ones(n_M,1); end;

flag_dtau = ~isempty(dtau_euler_polar_a_M_) | ~isempty(dtau_euler_azimu_b_M_) | ~isempty(dtau_euler_gamma_z_M_) ;
flag_dvol = ~isempty(dvol_a_k_Y_quad_yk_) | ~isempty(dvol_a_k_Y_quad_yk__) | ~isempty(dvol_a_k_p_quad_) ;

n_3 = 3;
ddssnll_helper_n_w_3; %<-- set indices. ;
ddssnll_helper_Y_lmk_3; %<-- set indices. ;
ddssnll_helper_weight_3d_riesz_3; %<-- set weights. ;
%ddssnll_helper_hires_3; %<-- set up k_hires_p_quad grid. ;
%disp('returning'); return;

ddssnll_helper_reconfigure_a_3; %<-- reconfigure a. ;

if  flag_dtau &  flag_dvol;
if (flag_verbose>0); disp(sprintf(' %%  flag_dtau %d &  flag_dvol %d = %d', flag_dtau, flag_dvol, flag_dtau &  flag_dvol)); end;
ddssnll_helper_reconfigure_dvol_a_3; %<-- reconfigure dvol_a. ;
ddssnll_helper_q2d_3; %<-- determine derivatives using dvol_a. ;
ddssnll_helper_reconfigure_dtau_3; %<-- reconfigure dtau. ;
ddssnll_helper_a_restore_from_kappa_3; %<-- use kappa_qpro_apply to construct a_restore. ;
ddssnll_helper_convert_a_restore_C2M0_3; %<-- Needed for Hvv_q3d. ;
end;%if  flag_dtau &  flag_dvol;

if ~flag_dtau &  flag_dvol;
if (flag_verbose>0); disp(sprintf(' %% ~flag_dtau %d &  flag_dvol %d = %d',~flag_dtau, flag_dvol,~flag_dtau &  flag_dvol)); end;
ddssnll_helper_reconfigure_dvol_a_3; %<-- reconfigure dvol_a. ;
ddssnll_helper_reconfigure_dtau_3; %<-- reconfigure dtau. ;
ddssnll_helper_q2d_3; %<-- determine derivatives using dvol_a. ;
ddssnll_helper_reconfigure_dtau_3; %<-- reconfigure dtau. ;
ddssnll_helper_dtau_firstorder_3; %<-- find critical dtau, given dvol_a. ;
ddssnll_helper_dtau_from_dtau_firstorder_3; %<-- if dtau is empty. ;
ddssnll_helper_q2d_3; %<-- determine derivatives using dvol_a. ;
ddssnll_helper_a_restore_from_kappa_3; %<-- use kappa_qpro_apply to construct a_restore. ;
ddssnll_helper_convert_a_restore_C2M0_3; %<-- Needed for Hvv_q3d. ;
end;%if ~flag_dtau &  flag_dvol;

if  flag_dtau & ~flag_dvol;
if (flag_verbose>0); disp(sprintf(' %%  flag_dtau %d & ~flag_dvol %d = %d', flag_dtau,~flag_dvol, flag_dtau & ~flag_dvol)); end;
ddssnll_helper_reconfigure_dtau_3; %<-- reconfigure dtau. ;
ddssnll_helper_a_restore_from_kappa_3; %<-- use kappa_qpro_apply to construct a_restore. ;
ddssnll_helper_convert_a_restore_C2M0_3; %<-- Needed for Hvv_q3d. ;
ddssnll_helper_convert_a_restore_C1M1_3;
ddssnll_helper_convert_a_restore_C0M2_3;
ddssnll_helper_convert_dtau_a_restore_C2M0_3;
ddssnll_helper_convert_dtau_a_restore_C1M1_3;
ddssnll_helper_convert_dtau_a_restore_C0M2_3;
ddssnll_helper_convert_dvol_firstorder_from_dtau_3; %<-- find critical dvol_a, given dtau. ;
ddssnll_helper_dvol_from_dvol_firstorder_3; %<-- if dvol is empty. ;
ddssnll_helper_reconfigure_dvol_a_3; %<-- reconfigure dvol_a. ;
ddssnll_helper_reconfigure_dtau_3; %<-- reconfigure dtau. ;
ddssnll_helper_q2d_3; %<-- determine derivatives using dvol_a. ;
ddssnll_helper_reconfigure_dtau_3; %<-- reconfigure dtau. ;
end;%if  flag_dtau & ~flag_dvol;

ddssnll_helper_convert_dtau_a_restore_C2M0_3; %<-- Needed for Hvv_q3d. ;
ddssnll_helper_convert_dtau_a_restore_C1M1_3; %<-- Needed for Hvv_q3d. ;

ddssnll_helper_convert_a_times_dtau_a_restore_C2M0_3;
ddssnll_helper_convert_dvol_a_times_a_restore_C2M0_3;

%%%%%%%%;
dtau_M3__ = [ ...
,dtau_euler_polar_a_M_ ...
,dtau_euler_azimu_b_M_ ...
,dtau_euler_gamma_z_M_ ...
] ;
dtau_3M__ = transpose(dtau_M3__);
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if flag_check;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

ddssnll_helper_convert_a_restore_C1M1_3;
ddssnll_helper_convert_a_restore_C0M2_3;
ddssnll_helper_convert_dtau_a_restore_C0M2_3;
ddssnll_helper_convert_dtau_dtau_a_restore_C2M0_3;
ddssnll_helper_convert_dtau_dtau_a_restore_C1M1_3;
ddssnll_helper_convert_dtau_dtau_a_restore_C0M2_3;
ddssnll_helper_convert_a_times_a_3;
ddssnll_helper_convert_dvol_a_times_dvol_a_3;
ddssnll_helper_convert_a_times_a_restore_C2M0_3;

ddssnll_helper_dtau_criticality_check_3;
ddssnll_helper_dvol_criticality_check_3;
ddssnll_helper_reconstruct_check_3;
ddssnll_helper_ssnll_check_3;
ddssnll_helper_dtau_ssnll_check_3;
ddssnll_helper_dtau_dtau_ssnll_check_3;
ddssnll_helper_note_hessian_3;
ddssnll_helper_Hvv_check_3;
ddssnll_helper_Htt_check_3;
ddssnll_helper_Htv_check_3;
ddssnll_helper_Hvt_check_3;
ddssnll_helper_imagesc_shell_3;
ddssnll_helper_Hv_q3d_check_3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_check;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% Now we calculate the left-hand-side: ;
%%%%%%%%;
Hvv_q3d_k_p_quad_ = dvol_a_k_p_quad_.*a_restore_C2M0_k_p_quad_;
Hvv_q3d_k_Y_quad_yk__ = dvol_a_times_a_restore_C2M0_k_Y_quad_yk__;
Htv_q2d_M3__ = dtau_dvol_ssnll_q2d_M3__;
Hvt_q3d_k_p_quad_ = (a_k_p_quad_.*dtau_a_restore_C2M0_k_p_quad_ - dtau_a_restore_C1M1_k_p_quad_);
Hvt_q3d_k_Y_quad_yk__ = (a_times_dtau_a_restore_C2M0_k_Y_quad_yk__ - dtau_a_restore_C1M1_k_Y_quad_yk__);
Htt_q2d_M3__ = sum(bsxfun(@times,dtau_dtau_ssnll_q2d_M33___,reshape(dtau_M3__,[n_M,1,3])),[3]);
%%%%%%%%;
Hv_q3d_k_p_quad_ = Hvv_q3d_k_p_quad_ + Hvt_q3d_k_p_quad_;
Hv_q3d_k_Y_quad_yk__ = Hvv_q3d_k_Y_quad_yk__ + Hvt_q3d_k_Y_quad_yk__;
Ht_q2d_M3__ = Htt_q2d_M3__ + Htv_q2d_M3__;
%%%%%%%%;
Hvv_q3d_k_Y_quad_yk_ = local_yk_from_yk__(n_k_p_r,l_max_,Hvv_q3d_k_Y_quad_yk__);
Hvt_q3d_k_Y_quad_yk_ = local_yk_from_yk__(n_k_p_r,l_max_,Hvt_q3d_k_Y_quad_yk__);
Hv_q3d_k_Y_quad_yk_ = local_yk_from_yk__(n_k_p_r,l_max_,Hv_q3d_k_Y_quad_yk__);

Hvt_ykabc_ = cat(1,Hv_q3d_k_Y_quad_yk_,Ht_q2d_M3__(:,1+0),Ht_q2d_M3__(:,1+1),Ht_q2d_M3__(:,1+2));

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
ddssnll_helper_ddssnll_check_3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
