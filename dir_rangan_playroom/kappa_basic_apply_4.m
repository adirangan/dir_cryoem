function ...
[ ...
 parameter ...
,KAPPA ...
,a_restore_C2M0_k_p_qk__ ...
,a_restore_C1M1_k_p_qk__ ...
,a_restore_C0M2_k_p_qk__ ...
,dtau_a_restore_C2M0_k_p_qk__ ...
,dtau_a_restore_C1M1_k_p_qk__ ...
,dtau_a_restore_C0M2_k_p_qk__ ...
,dtau_dtau_a_restore_C2M0_k_p_qk__ ...
,dtau_dtau_a_restore_C1M1_k_p_qk__ ...
,dtau_dtau_a_restore_C0M2_k_p_qk__ ...
,qref_k_eq_d ...
,qref_n_shell ...
,qref_azimu_b_shell_ ...
,qref_polar_a_shell_ ...
,qref_weight_shell_ ...
,qref_k_c_0_shell_ ...
,qref_k_c_1_shell_ ...
,qref_k_c_2_shell_ ...
,qref_n_polar_a ...
,qref_polar_a_ ...
,qref_n_azimu_b_ ...
] = ...
kappa_basic_apply_4( ...
 parameter ...
,KAPPA ...
,n_w_max ...
,n_M ...
,weight_imagecount_M_ ...
,viewing_polar_a_M_ ...
,viewing_azimu_b_M_ ...
,viewing_gamma_z_M_ ...
,dtau_viewing_polar_a_M_ ...
,dtau_viewing_azimu_b_M_ ...
,dtau_viewing_gamma_z_M_ ...
,n_k_p_r ...
,M_k_p_wkM__ ...
,n_CTF ...
,index_nCTF_from_nM_ ...
,CTF_k_p_wkC__ ...
,qref_k_eq_d ...
,qref_n_shell ...
,qref_azimu_b_shell_ ...
,qref_polar_a_shell_ ...
,qref_weight_shell_ ...
,qref_k_c_0_shell_ ...
,qref_k_c_1_shell_ ...
,qref_k_c_2_shell_ ...
,qref_n_polar_a ...
,qref_polar_a_ ...
,qref_n_azimu_b_ ...
);

str_thisfunction = 'kappa_basic_apply_4';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if nargin<1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
disp(sprintf(' %% testing %s',str_thisfunction));
disp(sprintf(' %% see kappa_basic_recon_4'));
disp('returning'); return;
end;%if nargin<1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); KAPPA=[]; end; na=na+1;
if (nargin<1+na); n_w_max=[]; end; na=na+1;
if (nargin<1+na); n_M=[]; end; na=na+1;
if (nargin<1+na); weight_imagecount_M_=[]; end; na=na+1;
if (nargin<1+na); viewing_polar_a_M_=[]; end; na=na+1;
if (nargin<1+na); viewing_azimu_b_M_=[]; end; na=na+1;
if (nargin<1+na); viewing_gamma_z_M_=[]; end; na=na+1;
if (nargin<1+na); dtau_viewing_polar_a_M_=[]; end; na=na+1;
if (nargin<1+na); dtau_viewing_azimu_b_M_=[]; end; na=na+1;
if (nargin<1+na); dtau_viewing_gamma_z_M_=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); M_k_p_wkM__=[]; end; na=na+1;
if (nargin<1+na); n_CTF=[]; end; na=na+1;
if (nargin<1+na); index_nCTF_from_nM_=[]; end; na=na+1;
if (nargin<1+na); CTF_k_p_wkC__=[]; end; na=na+1;
if (nargin<1+na); qref_k_eq_d=[]; end; na=na+1;
if (nargin<1+na); qref_n_shell=[]; end; na=na+1;
if (nargin<1+na); qref_azimu_b_shell_=[]; end; na=na+1;
if (nargin<1+na); qref_polar_a_shell_=[]; end; na=na+1;
if (nargin<1+na); qref_weight_shell_=[]; end; na=na+1;
if (nargin<1+na); qref_k_c_0_shell_=[]; end; na=na+1;
if (nargin<1+na); qref_k_c_1_shell_=[]; end; na=na+1;
if (nargin<1+na); qref_k_c_2_shell_=[]; end; na=na+1;
if (nargin<1+na); qref_n_polar_a=[]; end; na=na+1;
if (nargin<1+na); qref_polar_a_=[]; end; na=na+1;
if (nargin<1+na); qref_n_azimu_b_=[]; end; na=na+1;

%%%%;
if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
flag_verbose=parameter.flag_verbose;
if ~isfield(parameter,'flag_disp'); parameter.flag_disp=0; end;
flag_disp=parameter.flag_disp; nf=0;
if ~isfield(parameter,'flag_check'); parameter.flag_check=0; end;
flag_check=parameter.flag_check;
if ~isfield(parameter,'tolerance_master'); parameter.tolerance_master=1e-2; end;
tolerance_master=parameter.tolerance_master;
if ~isfield(parameter,'kernel_basic_l_max_use'); parameter.kernel_basic_l_max_use=49; end;
kernel_basic_l_max_use=parameter.kernel_basic_l_max_use;
if ~isfield(parameter,'kernel_basic_l_max_ext'); parameter.kernel_basic_l_max_ext=ceil(1.25*kernel_basic_l_max_use); end;
kernel_basic_l_max_ext=parameter.kernel_basic_l_max_ext;
if ~isfield(parameter,'kernel_basic_l_max_band'); parameter.kernel_basic_l_max_band=+Inf; end;
kernel_basic_l_max_band=parameter.kernel_basic_l_max_band;
if ~isfield(parameter,'kernel_basic_qref_k_eq_d_double'); parameter.kernel_basic_qref_k_eq_d_double=[]; end;
kernel_basic_qref_k_eq_d_double=parameter.kernel_basic_qref_k_eq_d_double;
%%%%;
if ~isfield(parameter,'flag_recalc_qref_from_data'); parameter.flag_recalc_qref_from_data=1; end;
flag_recalc_qref_from_data=parameter.flag_recalc_qref_from_data;
if ~isfield(parameter,'flag_recalc_dtau_qref_from_data'); parameter.flag_recalc_dtau_qref_from_data=1; end;
flag_recalc_dtau_qref_from_data=parameter.flag_recalc_dtau_qref_from_data;
if ~isfield(parameter,'flag_recalc_dtau_dtau_qref_from_data'); parameter.flag_recalc_dtau_dtau_qref_from_data=1; end;
flag_recalc_dtau_dtau_qref_from_data=parameter.flag_recalc_dtau_dtau_qref_from_data;
%%%%;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

flag_kernel_full = 1;
if flag_kernel_full~=1; disp(sprintf(' %% Warning, flag_kernel_full %d in %s',flag_kernel_full,str_thisfunction)); end;

tmp_t=tic();
n_w = n_w_max ;
n_w_sum = n_w_max*n_k_p_r;
if isempty(weight_imagecount_M_); weight_imagecount_M_ = ones(n_M,1); end;
CTF_k_p_wkM__ = [];
if ~isempty(n_M) & ~isempty(n_k_p_r);
if isempty(n_CTF); n_CTF = 1; end;
if isempty(index_nCTF_from_nM_); index_nCTF_from_nM_ = zeros(n_M,1); end;
if isempty(CTF_k_p_wkC__); CTF_k_p_wkC__ = ones(n_w_sum,n_CTF); end;
if (size(CTF_k_p_wkC__,1)==n_k_p_r); CTF_k_p_r_kC__ = CTF_k_p_wkC__; CTF_k_p_wkC__ = reshape(repmat(reshape(CTF_k_p_r_kC__,[1,n_k_p_r,n_CTF]),[n_w_max,1,1]),[n_w_sum,n_CTF]); end;
CTF_k_p_wkM__ = CTF_k_p_wkC__(:,1+index_nCTF_from_nM_);
end;%if ~isempty(n_M) & ~isempty(n_k_p_r);
tmp_t=toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% CTF_k_p_wkM__: %0.6fs',tmp_t)); end;

flag_dtau = ~isempty(dtau_viewing_polar_a_M_) | ~isempty(dtau_viewing_azimu_b_M_) | ~isempty(dtau_viewing_gamma_z_M_) ;
if isempty(dtau_viewing_polar_a_M_); dtau_viewing_polar_a_M_ = zeros(n_M,1); end;
if isempty(dtau_viewing_azimu_b_M_); dtau_viewing_azimu_b_M_ = zeros(n_M,1); end;
if isempty(dtau_viewing_gamma_z_M_); dtau_viewing_gamma_z_M_ = zeros(n_M,1); end;

if isempty(KAPPA);
tmp_t=tic();
parameter.flag_verbose = 0; parameter.flag_disp = 0; parameter.flag_check = 0;
[ ...
 parameter ...
,KAPPA ...
,qref_k_eq_d ...
,qref_n_shell ...
,qref_azimu_b_shell_ ...
,qref_polar_a_shell_ ...
,qref_weight_shell_ ...
,qref_k_c_0_shell_ ...
,qref_k_c_1_shell_ ...
,qref_k_c_2_shell_ ...
,qref_n_polar_a ...
,qref_polar_a_ ...
,qref_n_azimu_b_ ...
] = ...
kappa_basic_1( ...
 parameter ...
,KAPPA ...
,qref_k_eq_d ...
,qref_n_shell ...
,qref_azimu_b_shell_ ...
,qref_polar_a_shell_ ...
,qref_weight_shell_ ...
,qref_k_c_0_shell_ ...
,qref_k_c_1_shell_ ...
,qref_k_c_2_shell_ ...
,qref_n_polar_a ...
,qref_polar_a_ ...
,qref_n_azimu_b_ ...
);
parameter.flag_verbose = flag_verbose; parameter.flag_disp = flag_disp; parameter.flag_check = flag_check;
tmp_t=toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% kappa_basic_1: %0.6fs',tmp_t)); end;
end;%if isempty(KAPPA);
%%%%%%%%;

tmp_t=tic();
l_max = KAPPA.l_max_use;
Rz = KAPPA.Rz;
dRz = KAPPA.dRz;
Ry = KAPPA.Ry;
dRy = KAPPA.dRy;
kappa_norm_ = KAPPA.kappa_norm_;
qref_k_eq_d = KAPPA.qref_k_eq_d;
qref_n_shell = KAPPA.qref_n_shell;
n_q = qref_n_shell;
n_3 = 3;
qref_k_c_qc__ = KAPPA.qref_k_c_qc__;
a_full_node_ = KAPPA.a_full_node_;
chebfun_kernel_norm_ = KAPPA.chebfun_kernel_norm_;
deconvolve_q = sqrt(4*pi); %<-- basic. ;
tmp_t=toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% deconvolve_q: %0.6fs',tmp_t)); end;

if ~isfield(KAPPA,'n_w_max'); KAPPA.n_w_max = n_w_max; end;
if ~isfield(KAPPA,'n_M'); KAPPA.n_M = n_M; end;
if ~isfield(KAPPA,'viewing_polar_a_M_'); KAPPA.viewing_polar_a_M_ = viewing_polar_a_M_; end;
if ~isfield(KAPPA,'viewing_azimu_b_M_'); KAPPA.viewing_azimu_b_M_ = viewing_azimu_b_M_; end;
if ~isfield(KAPPA,'viewing_gamma_z_M_'); KAPPA.viewing_gamma_z_M_ = viewing_gamma_z_M_; end;
if ~isfield(KAPPA,'dtau_viewing_polar_a_M_'); KAPPA.dtau_viewing_polar_a_M_ = dtau_viewing_polar_a_M_; end;
if ~isfield(KAPPA,'dtau_viewing_azimu_b_M_'); KAPPA.dtau_viewing_azimu_b_M_ = dtau_viewing_azimu_b_M_; end;
if ~isfield(KAPPA,'dtau_viewing_gamma_z_M_'); KAPPA.dtau_viewing_gamma_z_M_ = dtau_viewing_gamma_z_M_; end;

flag_d0 = nargout>=1+2;
flag_d1 = (flag_dtau & nargout>=1+2+3);
flag_d2 = (flag_dtau & nargout>=1+2+3+3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if flag_recalc_qref_from_data | flag_recalc_dtau_qref_from_data | flag_recalc_dtau_dtau_qref_from_data;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%;
if flag_d0 & ~flag_d1 & ~flag_d2 ;
tmp_t=tic();
[ ...
 data_k_p_polar_a_wM__ ...
,data_k_p_azimu_b_wM__ ...
,data_k_c_0_wM__ ...
,data_k_c_1_wM__ ...
,data_k_c_2_wM__ ...
,data_k_p_r01_wM__ ...
] = ...
cg_rhs_2( ...
 n_M ...
,n_w ...
,viewing_polar_a_M_ ...
,viewing_azimu_b_M_ ...
,viewing_gamma_z_M_ ...
);
tmp_t=toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% flag_d0 cg_rhs_2: %0.6fs',tmp_t)); end;
end;%if flag_d0 & ~flag_d1 & ~flag_d2 ;
%%%%;
if  flag_d0 &  flag_d1 & ~flag_d2 ;
tmp_t=tic();
[ ...
 data_k_p_polar_a_wM__ ...
,data_k_p_azimu_b_wM__ ...
,data_k_c_0_wM__ ...
,data_k_c_1_wM__ ...
,data_k_c_2_wM__ ...
,data_k_p_r01_wM__ ...
,dtau_data_k_p_polar_a_wM3___ ...
,dtau_data_k_p_azimu_b_wM3___ ...
,dtau_data_k_c_0_wM3___ ...
,dtau_data_k_c_1_wM3___ ...
,dtau_data_k_c_2_wM3___ ...
,dtau_data_k_p_r01_wM3___ ...
] = ...
cg_rhs_2( ...
 n_M ...
,n_w ...
,viewing_polar_a_M_ ...
,viewing_azimu_b_M_ ...
,viewing_gamma_z_M_ ...
);
tmp_t=toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% flag_d1 cg_rhs_2: %0.6fs',tmp_t)); end;
end;%if  flag_d0 &  flag_d1 & ~flag_d2 ;
%%%%;
if  flag_d0 &  flag_d1 &  flag_d2 ;
tmp_t=tic();
[ ...
 data_k_p_polar_a_wM__ ...
,data_k_p_azimu_b_wM__ ...
,data_k_c_0_wM__ ...
,data_k_c_1_wM__ ...
,data_k_c_2_wM__ ...
,data_k_p_r01_wM__ ...
,dtau_data_k_p_polar_a_wM3___ ...
,dtau_data_k_p_azimu_b_wM3___ ...
,dtau_data_k_c_0_wM3___ ...
,dtau_data_k_c_1_wM3___ ...
,dtau_data_k_c_2_wM3___ ...
,dtau_data_k_p_r01_wM3___ ...
,dtau_dtau_data_k_p_polar_a_wM33____ ...
,dtau_dtau_data_k_p_azimu_b_wM33____ ...
,dtau_dtau_data_k_c_0_wM33____ ...
,dtau_dtau_data_k_c_1_wM33____ ...
,dtau_dtau_data_k_c_2_wM33____ ...
,dtau_dtau_data_k_p_r01_wM33____ ...
] = ...
cg_rhs_2( ...
 n_M ...
,n_w ...
,viewing_polar_a_M_ ...
,viewing_azimu_b_M_ ...
,viewing_gamma_z_M_ ...
);
tmp_t=toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% flag_d2 cg_rhs_2: %0.6fs',tmp_t)); end;
end;%if  flag_d0 &  flag_d1 &  flag_d2 ;
%%%%;
data_k_c_wMc__ = [ data_k_c_0_wM__(:) , data_k_c_1_wM__(:) , data_k_c_2_wM__(:) ];
%%%%;
if flag_check;
tmp_error = 0;
for nM=0:n_M-1; for nw=0:n_w-1;
tmp_euler_a = +viewing_polar_a_M_(1+nM); tmp_euler_b = +viewing_azimu_b_M_(1+nM); tmp_euler_c = -viewing_gamma_z_M_(1+nM) + (2*pi*nw)/max(1,n_w) ;
tmp_R__ = Rz(tmp_euler_b)*Ry(tmp_euler_a)*Rz(tmp_euler_c);
tmp_R_point_k_c_ = tmp_R__*[1;0;0];
tmp_diff_ = data_k_c_wMc__(1+nw+nM*n_w,:) - reshape(tmp_R_point_k_c_,[1,3]);
tmp_error = tmp_error + fnorm(tmp_diff_);
end;end;%for nw=0:n_w-1; for nM=0:n_M-1;
if (flag_verbose>0); disp(sprintf(' %% data_k_c_wMc__ error: %0.16f',tmp_error)); end;
end;%if flag_check;
%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_recalc_qref_from_data | flag_recalc_dtau_qref_from_data | flag_recalc_dtau_dtau_qref_from_data;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

if ~isfield(KAPPA,'qref_from_data_qwM__'); KAPPA.qref_from_data_qwM__ = []; end;
%%%%%%%%;
if flag_recalc_qref_from_data;
d2_full_wMq__ = sum(bsxfun(@minus,reshape(qref_k_c_qc__,[1,n_q,3]),reshape(+data_k_c_wMc__,[n_w*n_M,1,3])).^2,3);
d1_full_wMq__ = sqrt(d2_full_wMq__);
mollify_qref_from_data_wMq__ = chebfun_kernel_norm_(1-d2_full_wMq__/2);
weight_mollify_qref_from_data_wMq__ = reshape(bsxfun(@times,reshape(mollify_qref_from_data_wMq__,[n_w,n_M,n_q]),reshape(weight_imagecount_M_,[1,n_M,1])),[n_w*n_M,n_q]);
index_qref_wMq__ = repmat([0:n_q-1],[n_w*n_M,1]);
index_data_wMq__ = repmat(transpose(0:n_w*n_M-1),[1,n_q]);
qref_from_data_qwM__ = sparse(1+index_qref_wMq__(:),1+index_data_wMq__(:),weight_mollify_qref_from_data_wMq__,n_q,n_w*n_M);
KAPPA.qref_from_data_qwM__ = qref_from_data_qwM__;
end%if flag_recalc_qref_from_data;
%%%%%%%%;
qref_from_data_qwM__ = KAPPA.qref_from_data_qwM__;

if ~isfield(KAPPA,'dtau_qref_from_data_qwM__'); KAPPA.dtau_qref_from_data_qwM__ = []; end;
%%%%%%%%;
if flag_recalc_dtau_qref_from_data;
chebfun_kernel_ = chebfun_kernel_norm_;
qref_c11q____ = reshape(permute(qref_k_c_qc__,[2,1]),[n_3,1,1,n_q]);
k_c_cwM___ = permute(reshape(data_k_c_wMc__,[n_w,n_M,n_3]),[3,1,2]);
da_M_ = dtau_viewing_polar_a_M_;
db_M_ = dtau_viewing_azimu_b_M_;
dc_M_ = dtau_viewing_gamma_z_M_;
dtau_data_k_c_wM3c____ = cat(4,dtau_data_k_c_0_wM3___,dtau_data_k_c_1_wM3___,dtau_data_k_c_2_wM3___); %<-- 3c is abc,012 ;
dtau_data_k_c_cwM3____ = permute(dtau_data_k_c_wM3c____,[4,1,2,3]);
dtau_data_k_c_cwM___ = ...
  + bsxfun(@times,dtau_data_k_c_cwM3____(:,:,:,1+0),reshape(da_M_,[1,1,n_M,1])) ...
  + bsxfun(@times,dtau_data_k_c_cwM3____(:,:,:,1+1),reshape(db_M_,[1,1,n_M,1])) ...
  + bsxfun(@times,dtau_data_k_c_cwM3____(:,:,:,1+2),reshape(dc_M_,[1,1,n_M,1])) ...
  ;
distsquared_1wMq___ = reshape(sum(bsxfun(@minus,qref_c11q____,k_c_cwM___).^2,1),[n_w,n_M,n_q]);
d_distsquared_cwMq____ = -2*bsxfun(@minus,qref_c11q____,k_c_cwM___);
dtau_distsquared_1wMq___ = reshape(sum(bsxfun(@times,d_distsquared_cwMq____,dtau_data_k_c_cwM___),1),[n_w,n_M,n_q]);
dchebfun_kernel_ = diff(chebfun_kernel_);
ddchebfun_kernel_ = diff(dchebfun_kernel_);
dtau_data_kappa_1wMq___ = dchebfun_kernel_(1-distsquared_1wMq___/2).*(-0.5).*dtau_distsquared_1wMq___;
weight_dtau_data_kappa_1wMq___ = reshape(bsxfun(@times,reshape(dtau_data_kappa_1wMq___,[1,n_w,n_M,n_q]),reshape(weight_imagecount_M_,[1,1,n_M,1])),[1,n_w,n_M,n_q]);
dtau_qref_from_data_qwM__ = sparse(1+index_qref_wMq__(:),1+index_data_wMq__(:),weight_dtau_data_kappa_1wMq___(:),n_q,n_w*n_M);
KAPPA.dtau_qref_from_data_qwM__ = dtau_qref_from_data_qwM__;
end;%if flag_recalc_dtau_qref_from_data;
%%%%%%%%;
dtau_qref_from_data_qwM__ = KAPPA.dtau_qref_from_data_qwM__;

if ~isfield(KAPPA,'dtau_dtau_qref_from_data_qwM__'); KAPPA.dtau_dtau_qref_from_data_qwM__ = []; end;
%%%%%%%%;
if flag_recalc_dtau_dtau_qref_from_data;
dtau_dtau_data_k_c_wM33c_____ = cat(5,dtau_dtau_data_k_c_0_wM33____,dtau_dtau_data_k_c_1_wM33____,dtau_dtau_data_k_c_2_wM33____); %<-- 33c is abc,abc,012 ;
dtau_dtau_data_k_c_cwM33_____ = permute(dtau_dtau_data_k_c_wM33c_____,[5,1,2,3,4]);
dtau_dtau_data_k_c_cwM___ = ...
  + 1.0*bsxfun(@times,dtau_dtau_data_k_c_cwM33_____(:,:,:,1+0,1+0),reshape(da_M_,[1,1,n_M,1,1]).*reshape(da_M_,[1,1,n_M,1,1])) ...
  + 1.0*bsxfun(@times,dtau_dtau_data_k_c_cwM33_____(:,:,:,1+0,1+1),reshape(da_M_,[1,1,n_M,1,1]).*reshape(db_M_,[1,1,n_M,1,1])) ...
  + 1.0*bsxfun(@times,dtau_dtau_data_k_c_cwM33_____(:,:,:,1+0,1+2),reshape(da_M_,[1,1,n_M,1,1]).*reshape(dc_M_,[1,1,n_M,1,1])) ...
  + 1.0*bsxfun(@times,dtau_dtau_data_k_c_cwM33_____(:,:,:,1+1,1+0),reshape(db_M_,[1,1,n_M,1,1]).*reshape(da_M_,[1,1,n_M,1,1])) ...
  + 1.0*bsxfun(@times,dtau_dtau_data_k_c_cwM33_____(:,:,:,1+1,1+1),reshape(db_M_,[1,1,n_M,1,1]).*reshape(db_M_,[1,1,n_M,1,1])) ...
  + 1.0*bsxfun(@times,dtau_dtau_data_k_c_cwM33_____(:,:,:,1+1,1+2),reshape(db_M_,[1,1,n_M,1,1]).*reshape(dc_M_,[1,1,n_M,1,1])) ...
  + 1.0*bsxfun(@times,dtau_dtau_data_k_c_cwM33_____(:,:,:,1+2,1+0),reshape(dc_M_,[1,1,n_M,1,1]).*reshape(da_M_,[1,1,n_M,1,1])) ...
  + 1.0*bsxfun(@times,dtau_dtau_data_k_c_cwM33_____(:,:,:,1+2,1+1),reshape(dc_M_,[1,1,n_M,1,1]).*reshape(db_M_,[1,1,n_M,1,1])) ...
  + 1.0*bsxfun(@times,dtau_dtau_data_k_c_cwM33_____(:,:,:,1+2,1+2),reshape(dc_M_,[1,1,n_M,1,1]).*reshape(dc_M_,[1,1,n_M,1,1])) ...
  ;
dtau_dtau_distsquared_1wMq___ = ...
 + reshape(sum(bsxfun(@times,2*ones(1,1,1,n_q),bsxfun(@times,dtau_data_k_c_cwM___,dtau_data_k_c_cwM___)),1),[n_w,n_M,n_q]) ...
 + reshape(sum(bsxfun(@times,d_distsquared_cwMq____,dtau_dtau_data_k_c_cwM___),1),[n_w,n_M,n_q]) ...
;
dtau_dtau_data_kappa_1wMq___ = ...
 + dchebfun_kernel_(1-distsquared_1wMq___/2).*(-0.5).*dtau_dtau_distsquared_1wMq___ ...
 + ddchebfun_kernel_(1-distsquared_1wMq___/2).*(-0.5).*dtau_distsquared_1wMq___.*(-0.5).*dtau_distsquared_1wMq___ ...
;
weight_dtau_dtau_data_kappa_1wMq___ = reshape(bsxfun(@times,reshape(dtau_dtau_data_kappa_1wMq___,[1,n_w,n_M,n_q]),reshape(weight_imagecount_M_,[1,1,n_M,1])),[1,n_w,n_M,n_q]);
dtau_dtau_qref_from_data_qwM__ = sparse(1+index_qref_wMq__(:),1+index_data_wMq__(:),weight_dtau_dtau_data_kappa_1wMq___(:),n_q,n_w*n_M);
KAPPA.dtau_dtau_qref_from_data_qwM__ = dtau_dtau_qref_from_data_qwM__;
end;%if flag_recalc_dtau_dtau_qref_from_data;
%%%%%%%%;
dtau_dtau_qref_from_data_qwM__ = KAPPA.dtau_dtau_qref_from_data_qwM__;

tmp_t=tic();
M_k_p_wMk__ = [];
if ~isempty(M_k_p_wkM__);
M_k_p_wMk__ = reshape(permute(reshape(M_k_p_wkM__,[n_w,n_k_p_r,n_M]),[1,3,2]),[n_w*n_M,n_k_p_r]);
end;%if ~isempty(M_k_p_wkM__);
CTF_k_p_wMk__ = [];
if ~isempty(CTF_k_p_wkM__);
CTF_k_p_wMk__ = reshape(permute(reshape(CTF_k_p_wkM__,[n_w,n_k_p_r,n_M]),[1,3,2]),[n_w*n_M,n_k_p_r]);
end;%if ~isempty(CTF_k_p_wkM__);
tmp_t=toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% M_k_p_wMk__: %0.6fs',tmp_t)); end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
a_restore_C2M0_k_p_qk__ = [];
a_restore_C1M1_k_p_qk__ = [];
a_restore_C0M2_k_p_qk__ = [];
dtau_a_restore_C2M0_k_p_qk__ = [];
dtau_a_restore_C1M1_k_p_qk__ = [];
dtau_a_restore_C0M2_k_p_qk__ = [];
dtau_dtau_a_restore_C2M0_k_p_qk__ = [];
dtau_dtau_a_restore_C1M1_k_p_qk__ = [];
dtau_dtau_a_restore_C0M2_k_p_qk__ = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for index_calc=0:3-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if index_calc==0; flag_calc = (~2|~isempty(CTF_k_p_wMk__)) & (~0|~isempty(M_k_p_wMk__)); if flag_calc; N_k_p_wMk__ = abs(CTF_k_p_wMk__).^2                        ; end; end;
if index_calc==1; flag_calc = (~1|~isempty(CTF_k_p_wMk__)) & (~1|~isempty(M_k_p_wMk__)); if flag_calc; N_k_p_wMk__ =     CTF_k_p_wMk__ .^1 .*     M_k_p_wMk__ .^1 ; end; end;
if index_calc==2; flag_calc = (~0|~isempty(CTF_k_p_wMk__)) & (~2|~isempty(M_k_p_wMk__)); if flag_calc; N_k_p_wMk__ =                          abs(M_k_p_wMk__).^2 ; end; end;
if (flag_verbose>1); disp(sprintf(' %% index_calc %d <-- N_k_p_wMk__ (%d,%d)',index_calc,size(N_k_p_wMk__))); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if flag_calc;
N_k_p_wMk__  = N_k_p_wMk__ * (2*pi) / max(1,n_w_max) ; %<-- scale by number of points on image-ring. ;
%%%%%%%%%%%%%%%%;
if flag_d0;
%%%%%%%%;
a_k_p_qk__ = qref_from_data_qwM__*N_k_p_wMk__;
%%%%;
if flag_check;
%%%%;
b_k_p_qk__ = zeros(qref_n_shell,n_k_p_r);
for nM=0:n_M-1; for nw=0:n_w-1;
tmp_euler_a = +viewing_polar_a_M_(1+nM); tmp_euler_b = +viewing_azimu_b_M_(1+nM); tmp_euler_c = -viewing_gamma_z_M_(1+nM) + (2*pi*nw)/max(1,n_w) ;
tmp_R__ = Rz(tmp_euler_b)*Ry(tmp_euler_a)*Rz(tmp_euler_c); tmp_R_point_k_c_ = tmp_R__*[1;0;0];
tmp_R_distsquared_ = (qref_k_c_qc__(:,1+0) - tmp_R_point_k_c_(1+0)).^2 + (qref_k_c_qc__(:,1+1) - tmp_R_point_k_c_(1+1)).^2 + (qref_k_c_qc__(:,1+2) - tmp_R_point_k_c_(1+2)).^2 ;
b_k_p_qk__ = b_k_p_qk__ + chebfun_kernel_norm_(1-tmp_R_distsquared_/2)*N_k_p_wMk__(1+nw+nM*n_w,:);
end;end;%for nM=0:n_M-1; for nw=0:n_w-1;
if (flag_verbose>0); disp(sprintf(' %% b_k_p_qk__ vs a_k_p_qk__: %0.16f',fnorm(b_k_p_qk__-a_k_p_qk__)/max(1e-12,fnorm(b_k_p_qk__)))); end;
end;%if flag_check;
%%%%%%%%;
end;%if flag_d0;
%%%%;
if flag_d1;
dtau_a_k_p_qk__ = dtau_qref_from_data_qwM__*N_k_p_wMk__;
end;%if flag_d1;
if flag_d2;
dtau_dtau_a_k_p_qk__ = dtau_dtau_qref_from_data_qwM__*N_k_p_wMk__;
end;%if flag_d2;
%%%%%%%%%%%%%%%%;
if index_calc==0;
if flag_d0; a_restore_C2M0_k_p_qk__ = a_k_p_qk__*deconvolve_q; end;
if flag_d1; dtau_a_restore_C2M0_k_p_qk__ = dtau_a_k_p_qk__*deconvolve_q; end;
if flag_d2; dtau_dtau_a_restore_C2M0_k_p_qk__ = dtau_dtau_a_k_p_qk__*deconvolve_q; end; 
end;%if index_calc==0;
%%%%;
if index_calc==1;
if flag_d0; a_restore_C1M1_k_p_qk__ = a_k_p_qk__*deconvolve_q; end;
if flag_d1; dtau_a_restore_C1M1_k_p_qk__ = dtau_a_k_p_qk__*deconvolve_q; end;
if flag_d2; dtau_dtau_a_restore_C1M1_k_p_qk__ = dtau_dtau_a_k_p_qk__*deconvolve_q; end; 
end;%if index_calc==1;
%%%%;
if index_calc==2;
if flag_d0; a_restore_C0M2_k_p_qk__ = a_k_p_qk__*deconvolve_q; end;
if flag_d1; dtau_a_restore_C0M2_k_p_qk__ = dtau_a_k_p_qk__*deconvolve_q; end;
if flag_d2; dtau_dtau_a_restore_C0M2_k_p_qk__ = dtau_dtau_a_k_p_qk__*deconvolve_q; end; 
end;%if index_calc==2;
%%%%%%%%%%%%%%%%;
end;%if flag_calc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for index_calc=0:3-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
