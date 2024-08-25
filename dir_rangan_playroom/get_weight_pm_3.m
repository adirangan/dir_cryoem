function ...
[ ...
 parameter ...
,pm_n_k_p_r ...
,pm_k_p_r_ ...
,pm_k_p_r_max ...
,pm_l_max_ ...
,pm_n_k_all ...
,pm_n_k_all_csum_ ...
,pm_k_p_r_all_ ...
,pm_k_p_azimu_b_all_ ...
,pm_k_p_polar_a_all_ ...
,pm_weight_3d_k_all_ ...
,pm_weight_shell_k_ ...
,pm_weight_3d_k_p_r_ ...
,pm_n_w_ ...
,pm_weight_2d_k_p_r_ ...
,pm_weight_2d_wk_ ...
] = ...
get_weight_pm_3( ...
 parameter ...
,pm_n_UX_rank ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,l_max_ ...
,n_k_all ...
,n_k_all_csum_ ...
,k_p_r_all_ ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,weight_3d_k_all_ ...
,weight_shell_k_ ...
,weight_3d_k_p_r_ ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_wk_ ...
);

str_thisfunction = 'get_weight_pm_3';

if nargin<1;
disp(sprintf(' %% testing %s',str_thisfunction));

disp('returning'); return;
end;%if nargin<1;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); pm_n_UX_rank=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); k_p_r_=[]; end; na=na+1;
if (nargin<1+na); k_p_r_max=[]; end; na=na+1;
if (nargin<1+na); l_max_=[]; end; na=na+1;
if (nargin<1+na); n_k_all=[]; end; na=na+1;
if (nargin<1+na); n_k_all_csum_=[]; end; na=na+1;
if (nargin<1+na); k_p_r_all_=[]; end; na=na+1;
if (nargin<1+na); k_p_azimu_b_all_=[]; end; na=na+1;
if (nargin<1+na); k_p_polar_a_all_=[]; end; na=na+1;
if (nargin<1+na); weight_3d_k_all_=[]; end; na=na+1;
if (nargin<1+na); weight_shell_k_=[]; end; na=na+1;
if (nargin<1+na); weight_3d_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); n_w_=[]; end; na=na+1;
if (nargin<1+na); weight_2d_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); weight_2d_wk_=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
flag_verbose=parameter.flag_verbose;
if ~isfield(parameter,'flag_disp'); parameter.flag_disp=0; end;
flag_disp=parameter.flag_disp; nf=0;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);
%%%%%%%%;
pm_n_k_p_r = pm_n_UX_rank;
pm_k_p_r_ = ones(pm_n_k_p_r,1);
pm_k_p_r_max = 1.0;
pm_l_max_max = max(l_max_);
pm_l_max_ = pm_l_max_max*ones(pm_n_k_p_r,1);
pm_n_lm_ = (1+pm_l_max_).^2; pm_n_lm_sum = sum(pm_n_lm_);
pm_weight_3d_k_p_r_ = ones(pm_n_k_p_r,1);
pm_weight_2d_k_p_r_ = ones(pm_n_k_p_r,1);
pm_n_w_ = n_w_max*ones(pm_n_k_p_r,1);
pm_n_w_max = max(pm_n_w_);
pm_n_w_sum = sum(pm_n_w_);
pm_n_w_csum_ = cumsum([0;pm_n_w_]);
pm_weight_2d_wk_ = ones(pm_n_w_sum,1);
if (flag_verbose>0); disp(sprintf(' %% pm_n_k_p_r %d',pm_n_k_p_r)); end;
if (flag_verbose>0); disp(sprintf(' %% pm_l_max_max %d',pm_l_max_max)); end;
if (flag_verbose>0); disp(sprintf(' %% pm_n_w_max %d',pm_n_w_max)); end;
%%%%;
nk_p_r = n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
weight_3d_k_p_r = weight_3d_k_p_r_(1+nk_p_r);
weight_2d_k_p_r = weight_2d_k_p_r_(1+nk_p_r);
pm_k_p_r = 1.0;
pm_weight_3d_k_p_r = 1.0;
pm_weight_2d_k_p_r = 1.0;
%%%%;
tmp_index_3d_ = [n_k_all_csum_(1+nk_p_r):n_k_all_csum_(1+nk_p_r+1)-1];
pm_n_shell = numel(tmp_index_3d_);
pm_n_k_all = pm_n_k_p_r*pm_n_shell;
pm_n_k_all_csum_ = pm_n_shell*transpose([0:pm_n_k_p_r]);
assert(pm_n_shell==numel(tmp_index_3d_));
pm_k_p_r_all_ = ones(pm_n_k_all,1);
pm_k_p_azimu_b_all_ = repmat(k_p_azimu_b_all_(1+tmp_index_3d_),[pm_n_k_p_r,1]);
pm_k_p_polar_a_all_ = repmat(k_p_polar_a_all_(1+tmp_index_3d_),[pm_n_k_p_r,1]);
if (flag_verbose>0); disp(sprintf(' %% pm_n_shell %d',pm_n_shell)); end;
if (flag_disp>0);
pm_k_c_0_ = cos(pm_k_p_azimu_b_all_).*sin(pm_k_p_polar_a_all_);
pm_k_c_1_ = sin(pm_k_p_azimu_b_all_).*sin(pm_k_p_polar_a_all_);
pm_k_c_2_ = cos(pm_k_p_polar_a_all_);
figure(1+nf);nf=nf+1;clf;figsml;
plot3(pm_k_c_0_,pm_k_c_1_,pm_k_c_2_,'.'); axis equal; axis vis3d; axisnotick3d;
title('pm_k_c_x_','Interpreter','none');
end;%if (flag_disp>0);
%%%%;
weight_shell_single_ = weight_shell_k_(1+tmp_index_3d_);
weight_shell_single_sum = sum(weight_shell_single_);
if (flag_verbose>0); disp(sprintf(' %% weight_shell_single_sum %+0.6f',weight_shell_single_sum)); end;
if (flag_verbose>0); disp(sprintf(' %% 4*pi*k_p_r^2 %+0.6f',4*pi*k_p_r^2)); end;
pm_weight_shell_single_ = 4*pi*1.0^2*weight_shell_single_/max(1e-12,weight_shell_single_sum);
pm_weight_shell_single_sum = sum(pm_weight_shell_single_);
if (flag_verbose>0); disp(sprintf(' %% pm_weight_shell_single_sum %+0.6f',pm_weight_shell_single_sum)); end;
if (flag_verbose>0); disp(sprintf(' %% 4*pi*pm_k_p_r^2 %+0.6f',4*pi*pm_k_p_r^2)); end;
pm_weight_shell_k_ = repmat(pm_weight_shell_single_,[pm_n_k_p_r,1]);
%%%%;
weight_3d_single_ = weight_3d_k_all_(1+tmp_index_3d_);
weight_3d_single_sum = sum(weight_3d_single_);
if (flag_verbose>0); disp(sprintf(' %% weight_3d_single_sum %+0.6f',weight_3d_single_sum)); end;
if (flag_verbose>0); disp(sprintf(' %% weight_shell_single_sum*weight_3d_k_p_r/max(1e-12,k_p_r^2) %+0.6f',weight_shell_single_sum*weight_3d_k_p_r/max(1e-12,k_p_r^2))); end;
pm_weight_3d_single_ = pm_weight_shell_single_*1.0*1.0^2;
pm_weight_3d_single_sum = sum(pm_weight_3d_single_);
if (flag_verbose>0); disp(sprintf(' %% pm_weight_3d_single_sum %+0.6f',pm_weight_3d_single_sum)); end;
if (flag_verbose>0); disp(sprintf(' %% pm_weight_shell_single_sum*pm_weight_3d_k_p_r/max(1e-12,pm_k_p_r^2) %+0.6f',pm_weight_shell_single_sum*pm_weight_3d_k_p_r/max(1e-12,pm_k_p_r^2))); end;
pm_weight_3d_k_all_ = repmat(pm_weight_3d_single_,[pm_n_k_p_r,1]);
%%%%;
tmp_index_2d_ = [n_w_csum_(1+nk_p_r):n_w_csum_(1+nk_p_r+1)-1];
weight_2d_single_ = weight_2d_wk_(1+tmp_index_2d_);
weight_2d_single_sum = sum(weight_2d_single_);
if (flag_verbose>0); disp(sprintf(' %% weight_2d_single_sum %+0.6f',weight_2d_single_sum)); end;
if (flag_verbose>0); disp(sprintf(' %% weight_2d_k_p_r/(4*pi^2) %+0.6f',weight_2d_k_p_r/(4*pi^2))); end;
pm_weight_2d_single_ = 1.0/(4*pi^2) * weight_2d_single_/max(1e-12,weight_2d_single_sum);
pm_weight_2d_single_sum = sum(pm_weight_2d_single_);
if (flag_verbose>0); disp(sprintf(' %% pm_weight_2d_single_sum %+0.6f',pm_weight_2d_single_sum)); end;
if (flag_verbose>0); disp(sprintf(' %% pm_weight_2d_k_p_r/(4*pi^2) %+0.6f',pm_weight_2d_k_p_r/(4*pi^2))); end;
pm_weight_2d_wk_ = repmat(pm_weight_2d_single_,[pm_n_k_p_r,1]);
%%%%%%%%;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

