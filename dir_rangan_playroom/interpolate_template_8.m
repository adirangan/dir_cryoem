%{
if ~exist('wS_from_single_shell_sba__','var'); wS_from_single_shell_sba__=[]; end;
if ~exist('dwSda_from_single_shell_sba__','var'); dwSda_from_single_shell_sba__=[]; end;
if ~exist('dwSdb_from_single_shell_sba__','var'); dwSdb_from_single_shell_sba__=[]; end;
if ~exist('ddwSdaa_from_single_shell_sba__','var'); ddwSdaa_from_single_shell_sba__=[]; end;
if ~exist('ddwSdab_from_single_shell_sba__','var'); ddwSdab_from_single_shell_sba__=[]; end;
if ~exist('ddwSdbb_from_single_shell_sba__','var'); ddwSdbb_from_single_shell_sba__=[]; end;
if ~exist('n_R','var'); n_R = []; end;
if ~exist('R_use_R___','var'); R_use_R___ = []; end;
if ~exist('a_R_k_p_Rqk__','var'); a_R_k_p_Rqk__=[]; end;
if ~exist('ba_from_single_shell_Rbaba___','var'); ba_from_single_shell_Rbaba___=[]; end;
if ~exist('wT_from_R_single_shell_Rsba___','var'); wT_from_R_single_shell_Rsba___=[]; end;
if ~exist('dwTda_from_R_single_shell_Rsba___','var'); dwTda_from_R_single_shell_Rsba___=[]; end;
if ~exist('dwTdb_from_R_single_shell_Rsba___','var'); dwTdb_from_R_single_shell_Rsba___=[]; end;
if ~exist('ddwTdaa_from_R_single_shell_Rsba___','var'); ddwTdaa_from_R_single_shell_Rsba___=[]; end;
if ~exist('ddwTdab_from_R_single_shell_Rsba___','var'); ddwTdab_from_R_single_shell_Rsba___=[]; end;
if ~exist('ddwTdbb_from_R_single_shell_Rsba___','var'); ddwTdbb_from_R_single_shell_Rsba___=[]; end;
%}
function ...
[ ...
 parameter ...
,template_ori_wkS__ ...
,n_w ...
,n_viewing_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,viewing_gamma_z_S_ ...
,wS_from_single_shell_sba__ ...
,n_T ...
,index_nS_from_nT_ ...
,n_R ...
,R_use_R___ ...
,viewing_R_polar_a_mod_TR__ ...
,viewing_R_azimu_b_mod_TR__ ...
,nR_from_nT_ ...
,a_R_k_p_Rqk__ ...
,ba_from_single_shell_Rbaba___ ...
,wT_from_R_single_shell_Rsba___ ...
,dtemplateda_ori_wkS__ ...
,dtemplatedb_ori_wkS__ ...
,dtemplatedc_ori_wkS__ ...
,dwSda_from_single_shell_sba__ ...
,dwSdb_from_single_shell_sba__ ...
,dwTda_from_R_single_shell_Rsba___ ...
,dwTdb_from_R_single_shell_Rsba___ ...
,ddtemplatedaa_ori_wkS__ ...
,ddtemplatedab_ori_wkS__ ...
,ddtemplatedac_ori_wkS__ ...
,ddtemplatedbb_ori_wkS__ ...
,ddtemplatedbc_ori_wkS__ ...
,ddtemplatedcc_ori_wkS__ ...
,ddwSdaa_from_single_shell_sba__ ...
,ddwSdab_from_single_shell_sba__ ...
,ddwSdbb_from_single_shell_sba__ ...
,ddwTdaa_from_R_single_shell_Rsba___ ...
,ddwTdab_from_R_single_shell_Rsba___ ...
,ddwTdbb_from_R_single_shell_Rsba___ ...
] = ...
interpolate_template_8( ...
 parameter ...
,n_order ...
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
,a_k_p_qk_ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,n_w_0in ...
,n_viewing_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,viewing_gamma_z_S_ ...
,wS_from_single_shell_sba__ ...
,dwSda_from_single_shell_sba__ ...
,dwSdb_from_single_shell_sba__ ...
,ddwSdaa_from_single_shell_sba__ ...
,ddwSdab_from_single_shell_sba__ ...
,ddwSdbb_from_single_shell_sba__ ...
,n_R ...
,R_use_R___ ...
,a_R_k_p_Rqk__ ...
,ba_from_single_shell_Rbaba___ ...
,wT_from_R_single_shell_Rsba___ ...
,dwTda_from_R_single_shell_Rsba___ ...
,dwTdb_from_R_single_shell_Rsba___ ...
,ddwTdaa_from_R_single_shell_Rsba___ ...
,ddwTdab_from_R_single_shell_Rsba___ ...
,ddwTdbb_from_R_single_shell_Rsba___ ...
);
%%%%;
% Uses a_k_p_qk_ to evaluate templates on a collection of points on spherical shells. ;
% Each spherical-shell has the same resolution, determined by viewing_k_eq_d and template_k_eq_d and/or n_w_max. ;
% We assume that each spherical-shell is discretized with a (potentially adaptive) discretization with: ;
% a collection of equispaced longitudes (i.e., azimuthal-points) for each latitude (i.e., polar-point), with: ;
% polar-points equispaced from [0,pi]. ;
% Furthermore, we assume that the number of longitudes on each pole is 1, ;
% and that the number of longitudes at each other latitude is even. ;
%%%%;
% This code is largely similar to interpolate_template_5.m ;
% As an additional feature we calculate each quantity n_R times, ;
% once for the original quadrature-grid, and once for each rotated version of the grid. ;
%%%%;
% ;
% inputs: ;
% ;
% flag_verbose = integer verbosity_level. ;
% n_qk = integer total number of points in spherical discretization. ;
% n_qk_csum_ = integer array of size (n_k_p_r). n_qk_csum_(1+nk_p_r) is the number of points prior to shell nk_p_r. ;
% k_p_r_qk_ = double array of size (n_qk). k_p_r_qk_(1+na) = radius of point na. ;
% k_p_azimu_b_qk_ = double array of size (n_qk). k_p_azimu_b_qk_(1+na) = azimu_b of point na. ;
% k_p_polar_a_qk_ = double array of size (n_qk). k_p_polar_a_qk_(1+na) = azimu_b of point na. ;
% weight_3d_k_p_qk_ = double array of size (n_qk). weight_3d_k_p_qk_(1+na) = quadrature-weight (3d) for point na. ;
% weight_shell_qk_ = double array of size (n_qk). weight_shell_qk_(1+na) = quadrature-weight (shell) for point na. ;
% n_k_p_r = integer number of shells. ;
% k_p_r_ = double array of size (n_k_p_r). k_p_r_(1+nk_p_r) = radius of shell nk_p_r.
% k_p_r_max = double maximum k_value intended for radial integration. ;
% weight_3d_k_p_r_ = double array of size (n_k_p_r). weight_3d_k_p_r_(1+nk_p_r) = quadrature-weight for shell nk_p_r. ;
% k_c_0_qk_ = double array of size (n_qk). k_c_0_qk_(1+na) = k_c_0 of point na. ;
% k_c_1_qk_ = double array of size (n_qk). k_c_1_qk_(1+na) = k_c_1 of point na. ;
% k_c_2_qk_ = double array of size (n_qk). k_c_2_qk_(1+na) = k_c_2 of point na. ;
% n_polar_a_k_ = integer array of size (n_k_p_r). n_polar_a_k_(1+nk_p_r) = number of latitudinal-lines for shell nk_p_r. ;
% polar_a_ka__ = cell-array of size (n_k_p_r). polar_a_ka__{1+nk_p_r} = double array of size n_polar_a_k_(1+nk_p_r) storing latitudes. ;
% n_azimu_b_ka__ = cell-array of size (n_k_p_r). n_azimu_b_ka__{1+nk_p_r} = integer array of size n_polar_a_k_(1+nk_p_r) storing n_azimu_b per latitude. ;
% a_k_p_qk_ = complex array of size (n_qk). a_k_p_qk_(1+na) = function-value a for point na. ;
% viewing_k_eq_d = real equatorial-distance used for sampling viewing angles and templates. ;
% template_k_eq_d = real equatorial-distance used for sampling inplane-shifts along each template. ;
% n_w_0in = integer. used if template_k_eq_d <=0; desired n_w for templates. ;
% n_viewing_S = integer. number of viewing angles (i.e., number of templates) .;
% viewing_azimu_b_S_ = real array of size (n_viewing_S,1). ;
%                        azimu_b values for each template. ;
% viewing_polar_a_S_ = real array of size (n_viewing_S,1). ;
%                        polar_a values for each template. ;
% viewing_weight_S_ = real array of size (n_viewing_S,1). ;
%                       integration weight (on shell of radius 1) for each template. ;
% n_viewing_polar_a = integer. number of distinct polar_a across the viewing angles. ;
% viewing_polar_a_ = real array of size (n_viewing_polar_a,1). ;
%                    polar_a values for each viewing_polar_a_. ;
% n_viewing_azimu_b_ = integer array of size (n_viewing_polar_a,1). ;
%                      number of azimu_b values for each polar_a. ;
%                      These azimu_b values are assumed to be equispaced on [0,2*pi). ;
% viewing_gamma_z_S_ = real gamma_z value for each template. (typically 0.0). ;
% ;
% outputs: ;
% ;
% template_wkS__ = complex array of templates. ;
%                  template_wkS__(1+nw+nk_p_r*n_w_max,1+nS) ;
%                  stores template value for angle-index nw, radial-index nk_p_r, ;
%                  and viewing_azimu_b = viewing_azimu_b_all_(1+nS). ;
%                  and viewing_polar_a = viewing_polar_a_all_(1+nS). ;
% dtemplatedx_wkS__ = complex array analogous to template_wkS__. ;
%                     stores first-derivative of template with respect to: ;
%                     x==a: polar_a ; %<-- note that the first-derivative with respect to polar_a has a different sign than wignerd_c produces. ;
%                     x==b: azimu_b ;
%                     x==c: gamma_z ;
% ddtemplatedxy_wkS__ = complex array analogous to template_wkS__. ;
%                       stores second-derivative of template with respect to x and y (see above). ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

str_thisfunction = 'interpolate_template_8';

if nargin<1;
rng(0);
flag_verbose = 2; nf=1;
if (flag_verbose); disp(sprintf(' %% testing %s',str_thisfunction)); end;
test_interpolate_template_8;
disp(sprintf(' %% returning')); return;
end;% if nargin<7;
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_order=[]; end; na=na+1;
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
if (nargin<1+na); a_k_p_qk_=[]; end; na=na+1;
if (nargin<1+na); viewing_k_eq_d=[]; end; na=na+1;
if (nargin<1+na); template_k_eq_d=[]; end; na=na+1;
if (nargin<1+na); n_w_0in=[]; end; na=na+1;
if (nargin<1+na); n_viewing_S=[]; end; na=na+1;
if (nargin<1+na); viewing_azimu_b_S_=[]; end; na=na+1;
if (nargin<1+na); viewing_polar_a_S_=[]; end; na=na+1;
if (nargin<1+na); viewing_weight_S_=[]; end; na=na+1;
if (nargin<1+na); n_viewing_polar_a=[]; end; na=na+1;
if (nargin<1+na); viewing_polar_a_=[]; end; na=na+1;
if (nargin<1+na); n_viewing_azimu_b_=[]; end; na=na+1;
if (nargin<1+na); viewing_gamma_z_S_=[]; end; na=na+1;
if (nargin<1+na); wS_from_single_shell_sba__=[]; end; na=na+1;
if (nargin<1+na); dwSda_from_single_shell_sba__=[]; end; na=na+1;
if (nargin<1+na); dwSdb_from_single_shell_sba__=[]; end; na=na+1;
if (nargin<1+na); ddwSdaa_from_single_shell_sba__=[]; end; na=na+1;
if (nargin<1+na); ddwSdab_from_single_shell_sba__=[]; end; na=na+1;
if (nargin<1+na); ddwSdbb_from_single_shell_sba__=[]; end; na=na+1;
if (nargin<1+na); n_R=[]; end; na=na+1;
if (nargin<1+na); R_use_R___=[]; end; na=na+1;
if (nargin<1+na); a_R_k_p_Rqk__=[]; end; na=na+1;
if (nargin<1+na); ba_from_single_shell_Rbaba___=[]; end; na=na+1;
if (nargin<1+na); wT_from_R_single_shell_Rsba___=[]; end; na=na+1;
if (nargin<1+na); dwTda_from_R_single_shell_Rsba___=[]; end; na=na+1;
if (nargin<1+na); dwTdb_from_R_single_shell_Rsba___=[]; end; na=na+1;
if (nargin<1+na); ddwTdaa_from_R_single_shell_Rsba___=[]; end; na=na+1;
if (nargin<1+na); ddwTdab_from_R_single_shell_Rsba___=[]; end; na=na+1;
if (nargin<1+na); ddwTdbb_from_R_single_shell_Rsba___=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
flag_verbose=parameter.flag_verbose;
if ~isfield(parameter,'flag_check'); parameter.flag_check=0; end;
flag_check=parameter.flag_check; nf=0;
if ~isfield(parameter,'tolerance_pinv'); parameter.tolerance_pinv=1e-6; end;
tolerance_pinv=parameter.tolerance_pinv;
if ~isfield(parameter,'tolerance_w'); parameter.tolerance_w = 0.5*(2*pi)/max(1,n_w_0in); end;
tolerance_w=parameter.tolerance_w;
if ~isfield(parameter,'flag_parsimonious'); parameter.flag_parsimonious=1; end;
flag_parsimonious=parameter.flag_parsimonious; nf=0;

if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

if isempty(n_order); n_order = 5; end;

flag_1 = 1;
flag_d = (nargout>=17);
flag_dd = (nargout>=27);
if (flag_verbose>0); disp(sprintf(' %% flag_1 %d flag_d %d flag_dd %d',flag_1,flag_d,flag_dd)); end;

%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% Set n_w_.')); end;
%%%%%%%%;
n_w = n_w_0in;
n_w_max = n_w;
n_w_ = n_w_max*ones(n_k_p_r,1);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);
n_S = n_viewing_S;

%%%%%%%%;
% Construct templates under original frame. ;
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% First establish baseline calculation under original frame.')); end;
%%%%%%%%;
if ~exist('wS_from_single_shell_sba__','var'); wS_from_single_shell_sba__=[]; end;
if ~exist('dwSda_from_single_shell_sba__','var'); dwSda_from_single_shell_sba__=[]; end;
if ~exist('dwSdb_from_single_shell_sba__','var'); dwSdb_from_single_shell_sba__=[]; end;
if ~exist('ddwSdaa_from_single_shell_sba__','var'); ddwSdaa_from_single_shell_sba__=[]; end;
if ~exist('ddwSdab_from_single_shell_sba__','var'); ddwSdab_from_single_shell_sba__=[]; end;
if ~exist('ddwSdbb_from_single_shell_sba__','var'); ddwSdbb_from_single_shell_sba__=[]; end;
tmp_parameter = parameter;
tmp_parameter.flag_verbose = max(0,flag_verbose-1);
tmp_parameter.flag_attend = 0; %<-- Here we bypass any calls to cartesian_from_shell. ;
%%%%%%%%;
if  flag_1  & ~flag_d  & ~flag_dd;
tmp_t = tic();
[ ...
 tmp_parameter ...
,template_ori_wkS__ ...
,n_w ...
,n_viewing_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,viewing_gamma_z_S_ ...
,wS_from_single_shell_sba__ ...
] = ...
interpolate_template_5( ...
 tmp_parameter ...
,n_order ...
,n_qk ...
,n_qk_csum_ ...
,k_p_r_qk_ ...
,k_p_azimu_b_qk_ ...
,k_p_polar_a_qk_ ...
,weight_3d_k_p_qk_ ...
,weight_shell_qk_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,k_c_0_qk_ ...
,k_c_1_qk_ ...
,k_c_2_qk_ ...
,n_polar_a_k_ ...
,polar_a_ka__ ...
,n_azimu_b_ka__ ...
,a_k_p_qk_ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,n_w_0in ...
,n_viewing_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,viewing_gamma_z_S_ ...
,wS_from_single_shell_sba__ ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% interpolate_template_5: time %0.6fs',tmp_t)); end;
end;%if  flag_1  & ~flag_d  & ~flag_dd;
%%%%%%%%;
if  flag_1  &  flag_d  & ~flag_dd;
tmp_t = tic();
[ ...
 tmp_parameter ...
,template_ori_wkS__ ...
,n_w ...
,n_viewing_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,viewing_gamma_z_S_ ...
,wS_from_single_shell_sba__ ...
,dtemplateda_ori_wkS__ ...
,dtemplatedb_ori_wkS__ ...
,dtemplatedc_ori_wkS__ ...
,dwSda_from_single_shell_sba__ ...
,dwSdb_from_single_shell_sba__ ...
] = ...
interpolate_template_5( ...
 tmp_parameter ...
,n_order ...
,n_qk ...
,n_qk_csum_ ...
,k_p_r_qk_ ...
,k_p_azimu_b_qk_ ...
,k_p_polar_a_qk_ ...
,weight_3d_k_p_qk_ ...
,weight_shell_qk_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,k_c_0_qk_ ...
,k_c_1_qk_ ...
,k_c_2_qk_ ...
,n_polar_a_k_ ...
,polar_a_ka__ ...
,n_azimu_b_ka__ ...
,a_k_p_qk_ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,n_w_0in ...
,n_viewing_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,viewing_gamma_z_S_ ...
,wS_from_single_shell_sba__ ...
,dwSda_from_single_shell_sba__ ...
,dwSdb_from_single_shell_sba__ ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% interpolate_template_5: time %0.6fs',tmp_t)); end;
end;%if  flag_1  &  flag_d  & ~flag_dd;
%%%%%%%%;
if  flag_1  &  flag_d  &  flag_dd;
tmp_t = tic();
[ ...
 tmp_parameter ...
,template_ori_wkS__ ...
,n_w ...
,n_viewing_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,viewing_gamma_z_S_ ...
,wS_from_single_shell_sba__ ...
,dtemplateda_ori_wkS__ ...
,dtemplatedb_ori_wkS__ ...
,dtemplatedc_ori_wkS__ ...
,dwSda_from_single_shell_sba__ ...
,dwSdb_from_single_shell_sba__ ...
,ddtemplatedaa_ori_wkS__ ...
,ddtemplatedab_ori_wkS__ ...
,ddtemplatedac_ori_wkS__ ...
,ddtemplatedbb_ori_wkS__ ...
,ddtemplatedbc_ori_wkS__ ...
,ddtemplatedcc_ori_wkS__ ...
,ddwSdaa_from_single_shell_sba__ ...
,ddwSdab_from_single_shell_sba__ ...
,ddwSdbb_from_single_shell_sba__ ...
] = ...
interpolate_template_5( ...
 tmp_parameter ...
,n_order ...
,n_qk ...
,n_qk_csum_ ...
,k_p_r_qk_ ...
,k_p_azimu_b_qk_ ...
,k_p_polar_a_qk_ ...
,weight_3d_k_p_qk_ ...
,weight_shell_qk_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,k_c_0_qk_ ...
,k_c_1_qk_ ...
,k_c_2_qk_ ...
,n_polar_a_k_ ...
,polar_a_ka__ ...
,n_azimu_b_ka__ ...
,a_k_p_qk_ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,n_w_0in ...
,n_viewing_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,viewing_gamma_z_S_ ...
,wS_from_single_shell_sba__ ...
,dwSda_from_single_shell_sba__ ...
,dwSdb_from_single_shell_sba__ ...
,ddwSdaa_from_single_shell_sba__ ...
,ddwSdab_from_single_shell_sba__ ...
,ddwSdbb_from_single_shell_sba__ ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% interpolate_template_5: time %0.6fs',tmp_t)); end;
end;%if  flag_1  &  flag_d  &  flag_dd;
%%%%%%%%;
if flag_dd;
%ddtemplatedba_ori_wkS__ = ddtemplatedab_ori_wkS__; %<-- not explicitly referenced. ;
%ddtemplatedca_ori_wkS__ = ddtemplatedac_ori_wkS__; %<-- not explicitly referenced. ;
%ddtemplatedcb_ori_wkS__ = ddtemplatedbc_ori_wkS__; %<-- not explicitly referenced. ;
end;%if flag_dd;
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Now step through each of the additional rotations. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if isempty(n_R); n_R = 2; end;
if isempty(R_use_R___); R_use_R___ = cell(n_R,1); end;
if isempty(a_R_k_p_Rqk__); a_R_k_p_Rqk__ = cell(n_R,1); end;
if isempty(ba_from_single_shell_Rbaba___); ba_from_single_shell_Rbaba___ = cell(n_R,1); end;
if isempty(wT_from_R_single_shell_Rsba___); wT_from_R_single_shell_Rsba___ = cell(n_R,1); end;
if isempty(dwTda_from_R_single_shell_Rsba___); dwTda_from_R_single_shell_Rsba___ = cell(n_R,1); end;
if isempty(dwTdb_from_R_single_shell_Rsba___); dwTdb_from_R_single_shell_Rsba___ = cell(n_R,1); end;
if isempty(ddwTdaa_from_R_single_shell_Rsba___); ddwTdaa_from_R_single_shell_Rsba___ = cell(n_R,1); end;
if isempty(ddwTdab_from_R_single_shell_Rsba___); ddwTdab_from_R_single_shell_Rsba___ = cell(n_R,1); end;
if isempty(ddwTdbb_from_R_single_shell_Rsba___); ddwTdbb_from_R_single_shell_Rsba___ = cell(n_R,1); end;
viewing_R_polar_a_SR__ = zeros(n_S,n_R);
viewing_R_azimu_b_SR__ = zeros(n_S,n_R);
viewing_R_gamma_z_SR__ = zeros(n_S,n_R);
template_rec_RwkT___ = cell(n_R,1);
if flag_d;
dtemplateda_rec_RwkT___ = cell(n_R,1);
dtemplatedb_rec_RwkT___ = cell(n_R,1);
dtemplatedc_rec_RwkT___ = cell(n_R,1);
end;%if flag_d;
if flag_dd;
ddtemplatedaa_rec_RwkT___ = cell(n_R,1);
ddtemplatedab_rec_RwkT___ = cell(n_R,1);
ddtemplatedac_rec_RwkT___ = cell(n_R,1);
ddtemplatedbb_rec_RwkT___ = cell(n_R,1);
ddtemplatedbc_rec_RwkT___ = cell(n_R,1);
ddtemplatedcc_rec_RwkT___ = cell(n_R,1);
end;%if flag_dd;
%%%%%%%%%%%%%%%%;
for nR=0:n_R-1;
%%%%%%%%%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% nR %d/%d calling interpolate_template_R_1.',nR,n_R)); end;
R_use__ = R_use_R___{1+nR};
a_R_k_p_qk_ = a_R_k_p_Rqk__{1+nR};
ba_from_single_shell_baba__ = ba_from_single_shell_Rbaba___{1+nR};
wT_from_R_single_shell_sba__ = wT_from_R_single_shell_Rsba___{1+nR};
if flag_d;
dwTda_from_R_single_shell_sba__ = dwTda_from_R_single_shell_Rsba___{1+nR};
dwTdb_from_R_single_shell_sba__ = dwTdb_from_R_single_shell_Rsba___{1+nR};
end;%if flag_d;
if flag_dd;
ddwTdaa_from_R_single_shell_sba__ = ddwTdaa_from_R_single_shell_Rsba___{1+nR};
ddwTdab_from_R_single_shell_sba__ = ddwTdab_from_R_single_shell_Rsba___{1+nR};
ddwTdbb_from_R_single_shell_sba__ = ddwTdbb_from_R_single_shell_Rsba___{1+nR};
end;%if flag_dd;
%%%%%%%%;
if  flag_1 & ~flag_d & ~flag_dd;
tmp_t = tic();
[ ...
 parameter ...
,n_T ...
,index_nS_from_nT_ ...
,R_use__ ...
,a_R_k_p_qk_ ...
,ba_from_single_shell_baba__ ...
,viewing_R_polar_a_T_ ...
,viewing_R_azimu_b_T_ ...
,viewing_R_gamma_z_T_ ...
,template_rec_wkT__ ...
,wT_from_R_single_shell_sba__ ...
] = ...
interpolate_template_R_1( ...
 parameter ...
,flag_1 ...
,flag_d ...
,flag_dd ...
,n_order ...
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
,a_k_p_qk_ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,n_w_0in ...
,n_viewing_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,viewing_gamma_z_S_ ...
,wS_from_single_shell_sba__ ...
,[] ...
,[] ...
,[] ...
,[] ...
,[] ...
,template_ori_wkS__ ...
,[] ...
,[] ...
,[] ...
,[] ...
,[] ...
,[] ...
,[] ...
,[] ...
,[] ...
,nR ...
,R_use__ ...
,a_R_k_p_qk_ ...
,ba_from_single_shell_baba__ ...
,wT_from_R_single_shell_sba__ ...
,[] ...
,[] ...
,[] ...
,[] ...
,[] ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% interpolate_template_R_1: time %0.6fs',tmp_t)); end;
end;%if  flag_1 & ~flag_d & ~flag_dd;
%%%%%%%%;
if  flag_1 &  flag_d & ~flag_dd;
tmp_t = tic();
[ ...
 parameter ...
,n_T ...
,index_nS_from_nT_ ...
,R_use__ ...
,a_R_k_p_qk_ ...
,ba_from_single_shell_baba__ ...
,viewing_R_polar_a_T_ ...
,viewing_R_azimu_b_T_ ...
,viewing_R_gamma_z_T_ ...
,template_rec_wkT__ ...
,wT_from_R_single_shell_sba__ ...
,dtemplateda_rec_wkT__ ...
,dtemplatedb_rec_wkT__ ...
,dtemplatedc_rec_wkT__ ...
,dwTda_from_R_single_shell_sba__ ...
,dwTdb_from_R_single_shell_sba__ ...
] = ...
interpolate_template_R_1( ...
 parameter ...
,flag_1 ...
,flag_d ...
,flag_dd ...
,n_order ...
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
,a_k_p_qk_ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,n_w_0in ...
,n_viewing_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,viewing_gamma_z_S_ ...
,wS_from_single_shell_sba__ ...
,dwSda_from_single_shell_sba__ ...
,dwSdb_from_single_shell_sba__ ...
,[] ...
,[] ...
,[] ...
,template_ori_wkS__ ...
,dtemplateda_ori_wkS__ ...
,dtemplatedb_ori_wkS__ ...
,dtemplatedc_ori_wkS__ ...
,[] ...
,[] ...
,[] ...
,[] ...
,[] ...
,[] ...
,nR ...
,R_use__ ...
,a_R_k_p_qk_ ...
,ba_from_single_shell_baba__ ...
,wT_from_R_single_shell_sba__ ...
,dwTda_from_R_single_shell_sba__ ...
,dwTdb_from_R_single_shell_sba__ ...
,[] ...
,[] ...
,[] ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% interpolate_template_R_1: time %0.6fs',tmp_t)); end;
end;%if  flag_1 &  flag_d & ~flag_dd;
%%%%%%%%;
if  flag_1 &  flag_d &  flag_dd;
tmp_t = tic();
[ ...
 parameter ...
,n_T ...
,index_nS_from_nT_ ...
,R_use__ ...
,a_R_k_p_qk_ ...
,ba_from_single_shell_baba__ ...
,viewing_R_polar_a_T_ ...
,viewing_R_azimu_b_T_ ...
,viewing_R_gamma_z_T_ ...
,template_rec_wkT__ ...
,wT_from_R_single_shell_sba__ ...
,dtemplateda_rec_wkT__ ...
,dtemplatedb_rec_wkT__ ...
,dtemplatedc_rec_wkT__ ...
,dwTda_from_R_single_shell_sba__ ...
,dwTdb_from_R_single_shell_sba__ ...
,ddtemplatedaa_rec_wkT__ ...
,ddtemplatedab_rec_wkT__ ...
,ddtemplatedac_rec_wkT__ ...
,ddtemplatedbb_rec_wkT__ ...
,ddtemplatedbc_rec_wkT__ ...
,ddtemplatedcc_rec_wkT__ ...
,ddwTdaa_from_R_single_shell_sba__ ...
,ddwTdab_from_R_single_shell_sba__ ...
,ddwTdbb_from_R_single_shell_sba__ ...
] = ...
interpolate_template_R_1( ...
 parameter ...
,flag_1 ...
,flag_d ...
,flag_dd ...
,n_order ...
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
,a_k_p_qk_ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,n_w_0in ...
,n_viewing_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,viewing_gamma_z_S_ ...
,wS_from_single_shell_sba__ ...
,dwSda_from_single_shell_sba__ ...
,dwSdb_from_single_shell_sba__ ...
,ddwSdaa_from_single_shell_sba__ ...
,ddwSdab_from_single_shell_sba__ ...
,ddwSdbb_from_single_shell_sba__ ...
,template_ori_wkS__ ...
,dtemplateda_ori_wkS__ ...
,dtemplatedb_ori_wkS__ ...
,dtemplatedc_ori_wkS__ ...
,ddtemplatedaa_ori_wkS__ ...
,ddtemplatedab_ori_wkS__ ...
,ddtemplatedac_ori_wkS__ ...
,ddtemplatedbb_ori_wkS__ ...
,ddtemplatedbc_ori_wkS__ ...
,ddtemplatedcc_ori_wkS__ ...
,nR ...
,R_use__ ...
,a_R_k_p_qk_ ...
,ba_from_single_shell_baba__ ...
,wT_from_R_single_shell_sba__ ...
,dwTda_from_R_single_shell_sba__ ...
,dwTdb_from_R_single_shell_sba__ ...
,ddwTdaa_from_R_single_shell_sba__ ...
,ddwTdab_from_R_single_shell_sba__ ...
,ddwTdbb_from_R_single_shell_sba__ ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% interpolate_template_R_1: time %0.6fs',tmp_t)); end;
end;%if  flag_1 &  flag_d &  flag_dd;
%%%%%%%%;
R_use_R___{1+nR} = R_use__;
a_R_k_p_Rqk__{1+nR} = a_R_k_p_qk_;
ba_from_single_shell_Rbaba___{1+nR} = ba_from_single_shell_baba__;
viewing_R_polar_a_SR__(1+index_nS_from_nT_,1+nR) = viewing_R_polar_a_T_;
viewing_R_azimu_b_SR__(1+index_nS_from_nT_,1+nR) = viewing_R_azimu_b_T_;
viewing_R_gamma_z_SR__(1+index_nS_from_nT_,1+nR) = viewing_R_gamma_z_T_;
wT_from_R_single_shell_Rsba___{1+nR} = wT_from_R_single_shell_sba__;
if flag_d;
dwTda_from_R_single_shell_Rsba___{1+nR} = dwTda_from_R_single_shell_sba__;
dwTdb_from_R_single_shell_Rsba___{1+nR} = dwTdb_from_R_single_shell_sba__;
end;%if flag_d;
if flag_dd;
ddwTdaa_from_R_single_shell_Rsba___{1+nR} = ddwTdaa_from_R_single_shell_sba__;
ddwTdab_from_R_single_shell_Rsba___{1+nR} = ddwTdab_from_R_single_shell_sba__;
ddwTdbb_from_R_single_shell_Rsba___{1+nR} = ddwTdbb_from_R_single_shell_sba__;
end;%if flag_dd;
template_rec_RwkT___{1+nR} = template_rec_wkT__;
if flag_d;
dtemplateda_rec_RwkT___{1+nR} = dtemplateda_rec_wkT__;
dtemplatedb_rec_RwkT___{1+nR} = dtemplatedb_rec_wkT__;
dtemplatedc_rec_RwkT___{1+nR} = dtemplatedc_rec_wkT__;
end;%if flag_d;
if flag_dd;
ddtemplatedaa_rec_RwkT___{1+nR} = ddtemplatedaa_rec_wkT__;
ddtemplatedab_rec_RwkT___{1+nR} = ddtemplatedab_rec_wkT__;
ddtemplatedac_rec_RwkT___{1+nR} = ddtemplatedac_rec_wkT__;
ddtemplatedbb_rec_RwkT___{1+nR} = ddtemplatedbb_rec_wkT__;
ddtemplatedbc_rec_RwkT___{1+nR} = ddtemplatedbc_rec_wkT__;
ddtemplatedcc_rec_RwkT___{1+nR} = ddtemplatedcc_rec_wkT__;
end;%if flag_dd;
%%%%%%%%%%%%%%%%;
end;%for nR=0:n_R-1;
%%%%%%%%%%%%%%%%;

%%%%%%%%;
% Now select which representative to use. ;
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% Measure quality of each representative for each nR.')); end;
viewing_R_polar_a_mod_TR__ = zeros(n_T,n_R);
viewing_R_azimu_b_mod_TR__ = zeros(n_T,n_R);
for nT=0:n_T-1;
nS = index_nS_from_nT_(1+nT);
for nR=0:n_R-1;
viewing_R_polar_a_mod_TR__(1+nT,1+nR) = abs(periodize(viewing_R_polar_a_SR__(1+nS,1+nR),-pi/4,+pi/4));
viewing_R_azimu_b_mod_TR__(1+nT,1+nR) = abs(periodize(viewing_R_azimu_b_SR__(1+nS,1+nR),-pi/4,+pi/4));
end;%for nR=0:n_R-1;
end;%for nT=0:n_T-1;
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% Select representative for each nR.')); end;
nR_from_nT_ = zeros(n_T,1);
for nT=0:n_T-1;
[~,ijR] = max(min(viewing_R_polar_a_mod_TR__(1+nT,:),viewing_R_azimu_b_mod_TR__(1+nT,:))); nR = ijR-1;
nR_from_nT_(1+nT) = nR;
end;%for nT=0:n_T-1;
%%%%%%%%;

%%%%%%%%;
% Now overwrite original. ;
%%%%%%%%;
if flag_1;
if (flag_verbose>0); disp(sprintf(' %% Overwriting template_ori_wkS__.')); end;
for nT=0:n_T-1;
nS = index_nS_from_nT_(1+nT);
nR = nR_from_nT_(1+nT);
template_ori_wkS__(:,1+nS) = template_rec_RwkT___{1+nR}(:,1+nT);
end;%for nT=0:n_T-1;
end;%if flag_1;
%%%%%%%%;
if flag_d;
if (flag_verbose>0); disp(sprintf(' %% Overwriting dtemplateda_ori_wkS__.')); end;
for nT=0:n_T-1;
nS = index_nS_from_nT_(1+nT);
nR = nR_from_nT_(1+nT);
dtemplateda_ori_wkS__(:,1+nS) = dtemplateda_rec_RwkT___{1+nR}(:,1+nT);
dtemplatedb_ori_wkS__(:,1+nS) = dtemplatedb_rec_RwkT___{1+nR}(:,1+nT);
dtemplatedc_ori_wkS__(:,1+nS) = dtemplatedc_rec_RwkT___{1+nR}(:,1+nT);
end;%for nT=0:n_T-1;
end;%if flag_d;
%%%%%%%%;
if flag_dd;
if (flag_verbose>0); disp(sprintf(' %% Overwriting ddtemplatedaa_ori_wkS__.')); end;
for nT=0:n_T-1;
nS = index_nS_from_nT_(1+nT);
nR = nR_from_nT_(1+nT);
ddtemplatedaa_ori_wkS__(:,1+nS) = ddtemplatedaa_rec_RwkT___{1+nR}(:,1+nT);
ddtemplatedab_ori_wkS__(:,1+nS) = ddtemplatedab_rec_RwkT___{1+nR}(:,1+nT);
ddtemplatedac_ori_wkS__(:,1+nS) = ddtemplatedac_rec_RwkT___{1+nR}(:,1+nT);
ddtemplatedbb_ori_wkS__(:,1+nS) = ddtemplatedbb_rec_RwkT___{1+nR}(:,1+nT);
ddtemplatedbc_ori_wkS__(:,1+nS) = ddtemplatedbc_rec_RwkT___{1+nR}(:,1+nT);
ddtemplatedcc_ori_wkS__(:,1+nS) = ddtemplatedcc_rec_RwkT___{1+nR}(:,1+nT);
end;%for nT=0:n_T-1;
end;%if flag_dd;
%%%%%%%%;

if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;



