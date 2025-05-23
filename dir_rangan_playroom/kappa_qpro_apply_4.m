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
,a_restore_C2M0_k_Y_yk__ ...
,a_restore_C1M1_k_Y_yk__ ...
,a_restore_C0M2_k_Y_yk__ ...
,dtau_a_restore_C2M0_k_Y_yk__ ...
,dtau_a_restore_C1M1_k_Y_yk__ ...
,dtau_a_restore_C0M2_k_Y_yk__ ...
,dtau_dtau_a_restore_C2M0_k_Y_yk__ ...
,dtau_dtau_a_restore_C1M1_k_Y_yk__ ...
,dtau_dtau_a_restore_C0M2_k_Y_yk__ ...
] = ...
kappa_qpro_apply_4( ...
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
,CTF_k_p_r_kC__ ...
);

str_thisfunction = 'kappa_qpro_apply_4';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if nargin<1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
disp(sprintf(' %% testing %s',str_thisfunction));
disp(sprintf(' %% see kappa_qpro_apply_3 and kappa_qpro_recon_4'));
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
if (nargin<1+na); CTF_k_p_r_kC__=[]; end; na=na+1;

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
if ~isfield(parameter,'kernel_qpro_l_max_use'); parameter.kernel_qpro_l_max_use=49; end;
kernel_qpro_l_max_use=parameter.kernel_qpro_l_max_use;
if ~isfield(parameter,'kernel_qpro_l_max_ext'); parameter.kernel_qpro_l_max_ext=ceil(1.25*kernel_qpro_l_max_use); end;
kernel_qpro_l_max_ext=parameter.kernel_qpro_l_max_ext;
if ~isfield(parameter,'kernel_qpro_l_max_band'); parameter.kernel_qpro_l_max_band=+Inf; end;
kernel_qpro_l_max_band=parameter.kernel_qpro_l_max_band;
if ~isfield(parameter,'kernel_qpro_polar_a_pole_north'); parameter.kernel_qpro_polar_a_pole_north=2.5*pi/24; end;
kernel_qpro_polar_a_pole_north=min(pi/2,parameter.kernel_qpro_polar_a_pole_north);
parameter.kernel_qpro_polar_a_pole_north = kernel_qpro_polar_a_pole_north;
if ~isfield(parameter,'kernel_qpro_polar_a_pole_south'); parameter.kernel_qpro_polar_a_pole_south=1.5*pi/24; end;
kernel_qpro_polar_a_pole_south=min(pi/2,parameter.kernel_qpro_polar_a_pole_south);
parameter.kernel_qpro_polar_a_pole_south = kernel_qpro_polar_a_pole_south;
if ~isfield(parameter,'kernel_qpro_deconvolution_factor_max'); parameter.kernel_qpro_deconvolution_factor_max=1024; end;
kernel_qpro_deconvolution_factor_max=parameter.kernel_qpro_deconvolution_factor_max;
if ~isfield(parameter,'kernel_qpro_qref_k_eq_d_double'); parameter.kernel_qpro_qref_k_eq_d_double=0.5; end;
kernel_qpro_qref_k_eq_d_double=parameter.kernel_qpro_qref_k_eq_d_double;
if ~isfield(parameter,'kernel_qpro_MaxIterations'); parameter.kernel_qpro_MaxIterations=1024; end;
kernel_qpro_MaxIterations=parameter.kernel_qpro_MaxIterations;
%%%%;
if ~isfield(parameter,'flag_kernel_full'); parameter.flag_kernel_full=0; end;
flag_kernel_full=parameter.flag_kernel_full;
if (kernel_qpro_polar_a_pole_north + kernel_qpro_polar_a_pole_south > pi-1e-12);
flag_kernel_full = 1;
parameter.flag_kernel_full = flag_kernel_full;
end;%if (kernel_qpro_polar_a_pole_north + kernel_qpro_polar_a_pole_south > pi-1e-12);
%%%%;
if ~isfield(parameter,'flag_recalc_knn'); parameter.flag_recalc_knn=1; end;
flag_recalc_knn=parameter.flag_recalc_knn;
if ~isfield(parameter,'flag_recalc_qref_from_data'); parameter.flag_recalc_qref_from_data=1; end;
flag_recalc_qref_from_data=parameter.flag_recalc_qref_from_data;
if ~isfield(parameter,'flag_recalc_dtau_qref_from_data'); parameter.flag_recalc_dtau_qref_from_data=1; end;
flag_recalc_dtau_qref_from_data=parameter.flag_recalc_dtau_qref_from_data;
if ~isfield(parameter,'flag_recalc_dtau_dtau_qref_from_data'); parameter.flag_recalc_dtau_dtau_qref_from_data=1; end;
flag_recalc_dtau_dtau_qref_from_data=parameter.flag_recalc_dtau_dtau_qref_from_data;
%%%%;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

if flag_kernel_full~=1; disp(sprintf(' %% Warning, flag_kernel_full %d in %s',flag_kernel_full,str_thisfunction)); end;
if (kernel_qpro_polar_a_pole_north+kernel_qpro_polar_a_pole_south<pi-1e-12); disp(sprintf(' %% Warning, kernel_qpro_polar_a_pole_north+kernel_qpro_polar_a_pole_south %0.6f in %s',kernel_qpro_polar_a_pole_north+kernel_qpro_polar_a_pole_south,str_thisfunction)); end;

flag_kernel_k_p_use = 1;
flag_kernel_k_Y_use = (nargout>=12);

tmp_t=tic();
n_w = n_w_max ;
n_w_sum = n_w_max*n_k_p_r;
if isempty(weight_imagecount_M_); weight_imagecount_M_ = ones(n_M,1); end;
CTF_k_p_wkM__ = [];
if ~isempty(n_M) & ~isempty(n_k_p_r);
if isempty(n_CTF); n_CTF = 1; end;
if isempty(index_nCTF_from_nM_); index_nCTF_from_nM_ = zeros(n_M,1); end;
if isempty(CTF_k_p_r_kC__); CTF_k_p_r_kC__ = ones(n_k_p_r,n_CTF); end;
CTF_k_p_wkM__ = reshape(permute(repmat(CTF_k_p_r_kC__(:,1+index_nCTF_from_nM_),[1,1,n_w]),[3,1,2]),[n_w_sum,n_M]);
end;%if ~isempty(n_M) & ~isempty(n_k_p_r);
tmp_t=toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% CTF_k_p_wkM__: %0.6fs',tmp_t)); end;

flag_dtau = ~isempty(dtau_viewing_polar_a_M_) | ~isempty(dtau_viewing_azimu_b_M_) | ~isempty(dtau_viewing_gamma_z_M_) ;
if isempty(dtau_viewing_polar_a_M_); dtau_viewing_polar_a_M_ = zeros(n_M,1); end;
if isempty(dtau_viewing_azimu_b_M_); dtau_viewing_azimu_b_M_ = zeros(n_M,1); end;
if isempty(dtau_viewing_gamma_z_M_); dtau_viewing_gamma_z_M_ = zeros(n_M,1); end;

if isempty(KAPPA);
tmp_t=tic();
parameter.flag_verbose = 0; parameter.flag_disp = 0; parameter.flag_check = 0;
parameter.flag_kernel_k_Y_use = flag_kernel_k_Y_use;
[parameter,KAPPA] = kappa_qpro_1(parameter);
parameter.flag_verbose = flag_verbose; parameter.flag_disp = flag_disp; parameter.flag_check = flag_check;
tmp_t=toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% kappa_qpro_1: %0.6fs',tmp_t)); end;
end;%if isempty(KAPPA);

tmp_t=tic();
l_max = KAPPA.l_max_use;
Rz = KAPPA.Rz;
dRz = KAPPA.dRz;
Ry = KAPPA.Ry;
dRy = KAPPA.dRy;
kappa_norm_ = KAPPA.kappa_norm_;
qref_k_eq_d = KAPPA.qref_k_eq_d;
n_ring_north = KAPPA.n_ring_north;
n_nearest_north = KAPPA.n_nearest_north;
n_ring_south = KAPPA.n_ring_south;
n_nearest_south = KAPPA.n_nearest_south;
n_nearest_total = KAPPA.n_nearest_total;
n_n = n_nearest_total;
qref_n_shell = KAPPA.qref_n_shell;
qref_k_c_qc__ = KAPPA.qref_k_c_qc__;
a_keep_node_ = KAPPA.a_keep_node_;
chebfun_kernel_norm_qpro_ = KAPPA.chebfun_kernel_norm_qpro_;
Y_l_val_ = KAPPA.Y_l_val_use_;
Y_m_val_ = KAPPA.Y_m_val_use_;
if flag_kernel_k_Y_use; Ylm_weight_yq__ = KAPPA.Ylm_use_weight_yq__; end;
deconvolve_lm_ = (sqrt(4*pi)*sqrt(1+2*Y_l_val_).*(Y_l_val_<=KAPPA.l_max_band))./max(1e-12,kappa_norm_(1+Y_l_val_));
tmp_t=toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% deconvolve_lm_: %0.6fs',tmp_t)); end;

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
if flag_recalc_knn | flag_recalc_qref_from_data | flag_recalc_dtau_qref_from_data | flag_recalc_dtau_dtau_qref_from_data;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%;
if ( flag_d0 & ~flag_d1 & ~flag_d2) | (~flag_d0 & ~flag_d1 & ~flag_d2 &  flag_recalc_knn);
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
end;%if ( flag_d0 & ~flag_d1 & ~flag_d2) | (~flag_d0 & ~flag_d1 & ~flag_d2 &  flag_recalc_knn);
%%%%;
if  flag_d0 &  flag_d1 & ~flag_d2;
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
end;%if  flag_d0 &  flag_d1 & ~flag_d2;
%%%%;
if  flag_d0 &  flag_d1 &  flag_d2;
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
end;%if  flag_d0 &  flag_d1 &  flag_d2;
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
if flag_recalc_knn;
tmp_t=tic();
ij_north_wMn__ = zeros(n_w*n_M,n_nearest_north);
if n_nearest_north> 0;
[ij_north_wMn__] = knnsearch(qref_k_c_qc__,+data_k_c_wMc__,'K',n_nearest_north);
end;%if n_nearest_north> 0;
ij_south_wMn__ = zeros(n_w*n_M,n_nearest_south);
if n_nearest_south> 0;
[ij_south_wMn__] = knnsearch(qref_k_c_qc__,-data_k_c_wMc__,'K',n_nearest_south);
end;%if n_nearest_south> 0;
index_qref_from_data_wMn__ = [ ij_north_wMn__ , ij_south_wMn__ ] - 1;
index_data_wMn__ = repmat(transpose(0:n_w*n_M-1),[1,n_nearest_total]);
qref_cwMn____ = cat( ...
		     4 ...
		     ,permute(reshape(qref_k_c_qc__(ij_north_wMn__(:),:),[n_w,n_M,n_nearest_north,3]),[4,1,2,3]) ...
		     ,permute(reshape(qref_k_c_qc__(ij_south_wMn__(:),:),[n_w,n_M,n_nearest_south,3]),[4,1,2,3]) ...
		     );
data_k_c_cwM___ = permute(reshape(data_k_c_wMc__,[n_w,n_M,3]),[3,1,2]);
distsquared_1wMn___ = reshape(sum(bsxfun(@minus,qref_cwMn____,data_k_c_cwM___).^2,1),[n_w,n_M,n_n]);
tmp_t=toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% qref_cwMn____: %0.6fs',tmp_t)); end;
KAPPA.ij_north_wMn__ = ij_north_wMn__;
KAPPA.ij_south_wMn__ = ij_south_wMn__;
KAPPA.index_qref_from_data_wMn__ = index_qref_from_data_wMn__;
KAPPA.index_data_wMn__ = index_data_wMn__;
KAPPA.qref_cwMn____ = qref_cwMn____;
KAPPA.data_k_c_cwM___ = data_k_c_cwM___;
KAPPA.distsquared_1wMn___ = distsquared_1wMn___;
end;%if flag_recalc_knn;
ij_north_wMn__  = KAPPA.ij_north_wMn__;
ij_south_wMn__  = KAPPA.ij_south_wMn__;
index_qref_from_data_wMn__ = KAPPA.index_qref_from_data_wMn__;
index_data_wMn__ = KAPPA.index_data_wMn__;
qref_cwMn____  = KAPPA.qref_cwMn____;
data_k_c_cwM___  = KAPPA.data_k_c_cwM___;
distsquared_1wMn___  = KAPPA.distsquared_1wMn___;
%%%%;
if flag_d0;
%%%%;
if flag_check;
tmp_distsquared_north_wMn__ = sum(bsxfun(@minus,reshape(qref_k_c_qc__(ij_north_wMn__(:),:),[n_w*n_M,n_nearest_north,3]),reshape(+data_k_c_wMc__,[n_w*n_M,1,3])).^2,3);
tmp_distsquared_south_wMn__ = sum(bsxfun(@minus,reshape(qref_k_c_qc__(ij_south_wMn__(:),:),[n_w*n_M,n_nearest_south,3]),reshape(+data_k_c_wMc__,[n_w*n_M,1,3])).^2,3);
distsquared_keep_wMn__ = [ tmp_distsquared_north_wMn__ , tmp_distsquared_south_wMn__ ] ;
disp(sprintf(' %% distsquared_1wMn___ vs distsquared_keep_wMn__: %0.16f',fnorm(distsquared_1wMn___(:)-distsquared_keep_wMn__(:))));
end;%if flag_check;
tmp_t=tic();
distsquared_keep_wMn__ = reshape(distsquared_1wMn___,[n_w*n_M,n_n]);
mollify_qref_from_data_wMn__ = chebfun_kernel_norm_qpro_(1-distsquared_keep_wMn__/2);
weight_mollify_qref_from_data_wMn__ = reshape(bsxfun(@times,reshape(mollify_qref_from_data_wMn__,[n_w,n_M,n_nearest_total]),reshape(weight_imagecount_M_,[1,n_M,1])),[n_w*n_M,n_nearest_total]);
qref_from_data_qwM__ = ...
sparse( ...
 1+index_qref_from_data_wMn__ ...
,1+index_data_wMn__(:) ...
,weight_mollify_qref_from_data_wMn__ ...
,qref_n_shell ...
,n_w*n_M ...
);
tmp_t=toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% qref_from_data_qwM__: %0.6fs',tmp_t)); end;
%%%%;
end;%if flag_d0;
%%%%;
if flag_d1; %<-- directional derivative only: accumulate gradient over dtau_viewing_M3__. ;
%%%%;
tmp_t=tic();
chebfun_kernel_ = chebfun_kernel_norm_qpro_;
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
d_distsquared_cwMn____ = -2*bsxfun(@minus,qref_cwMn____,data_k_c_cwM___);
dtau_distsquared_1wMn___ = reshape(sum(bsxfun(@times,d_distsquared_cwMn____,dtau_data_k_c_cwM___),1),[n_w,n_M,n_n]);
dchebfun_kernel_ = diff(chebfun_kernel_);
dtau_kappa_1wMn___ = dchebfun_kernel_(1-distsquared_1wMn___/2).*(-0.5).*dtau_distsquared_1wMn___;
weight_dtau_kappa_1wMn___ = reshape(bsxfun(@times,reshape(dtau_kappa_1wMn___,[n_w,n_M,n_nearest_total]),reshape(weight_imagecount_M_,[1,n_M,1])),[n_w,n_M,n_nearest_total]);
dtau_qref_from_data_qwM__ = ...
sparse( ...
 1+index_qref_from_data_wMn__ ...
,1+index_data_wMn__(:) ...
,weight_dtau_kappa_1wMn___(:) ...
,qref_n_shell ...
,n_w*n_M ...
);
tmp_t=toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% dtau_qref_from_data_qwM__: %0.6fs',tmp_t)); end;
%%%%;
end;%if flag_d1; %<-- directional derivative only: accumulate gradient over dtau_viewing_M3__. ;
%%%%;
if flag_d2; %<-- directional derivative only: accumulate hessian over dtau_viewing_M3__. ;
%%%%;
tmp_t=tic();
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
dtau_dtau_distsquared_1wMn___ = ...
 + reshape(sum(bsxfun(@times,2*ones(1,1,1,n_n),bsxfun(@times,dtau_data_k_c_cwM___,dtau_data_k_c_cwM___)),1),[n_w,n_M,n_n]) ...
 + reshape(sum(bsxfun(@times,d_distsquared_cwMn____,dtau_dtau_data_k_c_cwM___),1),[n_w,n_M,n_n]) ...
;
ddchebfun_kernel_ = diff(dchebfun_kernel_);
dtau_dtau_kappa_1wMn___ = ...
 + dchebfun_kernel_(1-distsquared_1wMn___/2).*(-0.5).*dtau_dtau_distsquared_1wMn___ ...
 + ddchebfun_kernel_(1-distsquared_1wMn___/2).*(-0.5).*dtau_distsquared_1wMn___.*(-0.5).*dtau_distsquared_1wMn___ ...
;
weight_dtau_dtau_kappa_1wMn___ = reshape(bsxfun(@times,reshape(dtau_dtau_kappa_1wMn___,[n_w,n_M,n_nearest_total]),reshape(weight_imagecount_M_,[1,n_M,1])),[n_w,n_M,n_nearest_total]);
dtau_dtau_qref_from_data_qwM__ = ...
sparse( ...
 1+index_qref_from_data_wMn__ ...
,1+index_data_wMn__(:) ...
,weight_dtau_dtau_kappa_1wMn___(:) ...
,qref_n_shell ...
,n_w*n_M ...
);
tmp_t=toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% dtau_dtau_qref_from_data_qwM__: %0.6fs',tmp_t)); end;
%%%%;
end;%if flag_d2; %<-- directional derivative only: accumulate hessian over dtau_viewing_M3__. ;
%%%%;
if flag_d0; KAPPA.qref_from_data_qwM__ = qref_from_data_qwM__; end;
if flag_d1; KAPPA.dtau_qref_from_data_qwM__ = dtau_qref_from_data_qwM__; end;
if flag_d2; KAPPA.dtau_dtau_qref_from_data_qwM__ = dtau_dtau_qref_from_data_qwM__; end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_recalc_knn | flag_recalc_qref_from_data | flag_recalc_dtau_qref_from_data | flag_recalc_dtau_dtau_qref_from_data;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if flag_d0; qref_from_data_qwM__ = KAPPA.qref_from_data_qwM__; end;
if flag_d1; dtau_qref_from_data_qwM__ = KAPPA.dtau_qref_from_data_qwM__; end;
if flag_d2; dtau_dtau_qref_from_data_qwM__ = KAPPA.dtau_dtau_qref_from_data_qwM__; end;

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
a_restore_C2M0_k_Y_yk__ = [];
a_restore_C1M1_k_Y_yk__ = [];
a_restore_C0M2_k_Y_yk__ = [];
dtau_a_restore_C2M0_k_Y_yk__ = [];
dtau_a_restore_C1M1_k_Y_yk__ = [];
dtau_a_restore_C0M2_k_Y_yk__ = [];
dtau_dtau_a_restore_C2M0_k_Y_yk__ = [];
dtau_dtau_a_restore_C1M1_k_Y_yk__ = [];
dtau_dtau_a_restore_C0M2_k_Y_yk__ = [];
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
if flag_kernel_k_Y_use;
tmp_t=tic();
a_k_Y_yk__ = conj(Ylm_weight_yq__)*a_k_p_qk__;
a_restore_k_Y_yk__ = bsxfun(@times,a_k_Y_yk__,deconvolve_lm_);
tmp_t=toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% a_restore_k_Y_yk__: %0.6fs',tmp_t)); end;
end;%if flag_kernel_k_Y_use;
%%%%;
if 1 | flag_check;
%%%%;
b_k_p_qk__ = zeros(qref_n_shell,n_k_p_r);
for nM=0:n_M-1; for nw=0:n_w-1;
tmp_euler_a = +viewing_polar_a_M_(1+nM); tmp_euler_b = +viewing_azimu_b_M_(1+nM); tmp_euler_c = -viewing_gamma_z_M_(1+nM) + (2*pi*nw)/max(1,n_w) ;
tmp_R__ = Rz(tmp_euler_b)*Ry(tmp_euler_a)*Rz(tmp_euler_c); tmp_R_point_k_c_ = tmp_R__*[1;0;0];
tmp_R_distsquared_ = (qref_k_c_qc__(:,1+0) - tmp_R_point_k_c_(1+0)).^2 + (qref_k_c_qc__(:,1+1) - tmp_R_point_k_c_(1+1)).^2 + (qref_k_c_qc__(:,1+2) - tmp_R_point_k_c_(1+2)).^2 ;
b_k_p_qk__ = b_k_p_qk__ + chebfun_kernel_norm_qpro_(1-tmp_R_distsquared_/2)*N_k_p_wMk__(1+nw+nM*n_w,:);
end;end;%for nM=0:n_M-1; for nw=0:n_w-1;
if (flag_verbose>0); disp(sprintf(' %% b_k_p_qk__ vs a_k_p_qk__: %0.16f',fnorm(b_k_p_qk__-a_k_p_qk__)/max(1e-12,fnorm(b_k_p_qk__)))); end;
%%%%;
if flag_kernel_k_Y_use;
b_k_Y_yk__ = conj(Ylm_weight_yq__)*b_k_p_qk__;
b_restore_k_Y_yk__ = bsxfun(@times,b_k_Y_yk__,deconvolve_lm_);
if (flag_verbose>0); disp(sprintf(' %% b_k_Y_yk__ vs a_k_Y_yk__: %0.16f',fnorm(b_k_Y_yk__-a_k_Y_yk__)/max(1e-12,fnorm(b_k_Y_yk__)))); end;
if (flag_verbose>0); disp(sprintf(' %% b_restore_k_Y_yk__ vs a_restore_k_Y_yk__: %0.16f',fnorm(b_restore_k_Y_yk__-a_restore_k_Y_yk__)/max(1e-12,fnorm(b_restore_k_Y_yk__)))); end;
end;%if flag_kernel_k_Y_use;
%%%%;
if flag_kernel_k_Y_use;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figmed;
nk_p_r = n_k_p_r-1;
subplot(1,2,1);
plot(Y_l_val_,log10(abs(a_k_Y_yk__(:,1+nk_p_r)-b_k_Y_yk__(:,1+nk_p_r))),'.');
xlim([0,l_max]); xlabel('Y_l_val_','Interpreter','none');
ylabel('log10(abs(a_k_Y_yk__-b_k_Y_yk__))','Interpreter','none');
title('log10(abs(a_crop - b_full))','Interpreter','none');
subplot(1,2,2);
plot(Y_l_val_,log10(abs(a_restore_k_Y_yk__(:,1+nk_p_r)-b_restore_k_Y_yk__(:,1+nk_p_r))),'.');
xlim([0,l_max]); xlabel('Y_l_val_','Interpreter','none');
ylabel('log10(abs(a_restore_k_Y_yk__-b_restore_k_Y_yk__))','Interpreter','none');
title('log10(abs(a_restore_crop - b_restore_full))','Interpreter','none');
end;%if flag_disp;
end;%if flag_kernel_k_Y_use;
%%%%;
if flag_kernel_k_Y_use;
if ~exist('qref_Ylm_uklma___','var'); qref_Ylm_uklma___=[]; end;
if ~exist('qref_k_p_azimu_b_sub_uka__','var'); qref_k_p_azimu_b_sub_uka__=[]; end;
if ~exist('qref_k_p_polar_a_sub_uka__','var'); qref_k_p_polar_a_sub_uka__=[]; end;
if ~exist('qref_l_max_uk_','var'); qref_l_max_uk_=[]; end;
if ~exist('qref_index_nu_n_k_per_shell_from_nk_p_r_','var'); qref_index_nu_n_k_per_shell_from_nk_p_r_=[]; end;
if ~exist('qref_index_k_per_shell_uka__','var'); qref_index_k_per_shell_uka__=[]; end;
tmp_t = tic();
[ ...
 c_k_p_qk__ ...
,qref_Ylm_uklma___ ...
,qref_k_p_azimu_b_sub_uka__ ...
,qref_k_p_polar_a_sub_uka__ ...
,qref_l_max_uk_ ...
,qref_index_nu_n_k_per_shell_from_nk_p_r_ ...
,qref_index_k_per_shell_uka__ ...
] = ...
convert_spharm_to_k_p_4( ...
 0*flag_verbose ...
,KAPPA.qref_n_shell ...
,[0,KAPPA.qref_n_shell] ...
,ones(KAPPA.qref_n_shell,1) ...
,KAPPA.qref_azimu_b_shell_ ...
,KAPPA.qref_polar_a_shell_ ...
,4*pi ...
,KAPPA.qref_weight_shell_ ...
,1 ...
,1.0 ...
,4*pi ...
,l_max ...
,a_restore_k_Y_yk__ ...
,qref_Ylm_uklma___ ...
,qref_k_p_azimu_b_sub_uka__ ...
,qref_k_p_polar_a_sub_uka__ ...
,qref_l_max_uk_ ...
,qref_index_nu_n_k_per_shell_from_nk_p_r_ ...
,qref_index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% convert_spharm_to_k_p_4: %0.6fs',tmp_t)); end;
fnorm_disp(flag_verbose,'a_k_p_qk__*sqrt(4*pi)',a_k_p_qk__*sqrt(4*pi),'c_k_p_qk__',c_k_p_qk__);
end;%if flag_kernel_k_Y_use;
%%%%;
end;%if flag_check;
%%%%%%%%;
end;%if flag_d0;
%%%%;
if flag_d1;
dtau_a_k_p_qk__ = dtau_qref_from_data_qwM__*N_k_p_wMk__;
if flag_kernel_k_Y_use;
tmp_t=tic();
dtau_a_k_Y_yk__ = conj(Ylm_weight_yq__)*dtau_a_k_p_qk__;
dtau_a_restore_k_Y_yk__ = bsxfun(@times,dtau_a_k_Y_yk__,deconvolve_lm_);
tmp_t=toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% dtau_a_restore_k_Y_yk__: %0.6fs',tmp_t)); end;
end;%if flag_kernel_k_Y_use;
end;%if flag_d1;
if flag_d2;
dtau_dtau_a_k_p_qk__ = dtau_dtau_qref_from_data_qwM__*N_k_p_wMk__;
if flag_kernel_k_Y_use;
tmp_t=tic();
dtau_dtau_a_k_Y_yk__ = conj(Ylm_weight_yq__)*dtau_dtau_a_k_p_qk__;
dtau_dtau_a_restore_k_Y_yk__ = bsxfun(@times,dtau_dtau_a_k_Y_yk__,deconvolve_lm_);
tmp_t=toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% dtau_dtau_a_restore_k_Y_yk__: %0.6fs',tmp_t)); end;
end;%if flag_kernel_k_Y_use;
end;%if flag_d2;
%%%%%%%%%%%%%%%%;
if index_calc==0;
if flag_d0; a_restore_C2M0_k_p_qk__ = a_k_p_qk__*sqrt(4*pi); end;
if flag_d1; dtau_a_restore_C2M0_k_p_qk__ = dtau_a_k_p_qk__*sqrt(4*pi); end;
if flag_d2; dtau_dtau_a_restore_C2M0_k_p_qk__ = dtau_dtau_a_k_p_qk__*sqrt(4*pi); end; 
end;%if index_calc==2;
%%%%;
if index_calc==1;
if flag_d0; a_restore_C1M1_k_p_qk__ = a_k_p_qk__*sqrt(4*pi); end;
if flag_d1; dtau_a_restore_C1M1_k_p_qk__ = dtau_a_k_p_qk__*sqrt(4*pi); end;
if flag_d2; dtau_dtau_a_restore_C1M1_k_p_qk__ = dtau_dtau_a_k_p_qk__*sqrt(4*pi); end; 
end;%if index_calc==3;
%%%%;
if index_calc==2;
if flag_d0; a_restore_C0M2_k_p_qk__ = a_k_p_qk__*sqrt(4*pi); end;
if flag_d1; dtau_a_restore_C0M2_k_p_qk__ = dtau_a_k_p_qk__*sqrt(4*pi); end;
if flag_d2; dtau_dtau_a_restore_C0M2_k_p_qk__ = dtau_dtau_a_k_p_qk__*sqrt(4*pi); end; 
end;%if index_calc==4;
%%%%%%%%;
if flag_kernel_k_Y_use;
%%%%;
if index_calc==0;
if flag_d0; a_restore_C2M0_k_Y_yk__ = a_restore_k_Y_yk__; end;
if flag_d1; dtau_a_restore_C2M0_k_Y_yk__ = dtau_a_restore_k_Y_yk__; end;
if flag_d2; dtau_dtau_a_restore_C2M0_k_Y_yk__ = dtau_dtau_a_restore_k_Y_yk__; end; 
end;%if index_calc==2;
%%%%;
if index_calc==1;
if flag_d0; a_restore_C1M1_k_Y_yk__ = a_restore_k_Y_yk__; end;
if flag_d1; dtau_a_restore_C1M1_k_Y_yk__ = dtau_a_restore_k_Y_yk__; end;
if flag_d2; dtau_dtau_a_restore_C1M1_k_Y_yk__ = dtau_dtau_a_restore_k_Y_yk__; end; 
end;%if index_calc==3;
%%%%;
if index_calc==2;
if flag_d0; a_restore_C0M2_k_Y_yk__ = a_restore_k_Y_yk__; end;
if flag_d1; dtau_a_restore_C0M2_k_Y_yk__ = dtau_a_restore_k_Y_yk__; end;
if flag_d2; dtau_dtau_a_restore_C0M2_k_Y_yk__ = dtau_dtau_a_restore_k_Y_yk__; end; 
end;%if index_calc==4;
%%%%;
end;%if flag_kernel_k_Y_use;
%%%%%%%%%%%%%%%%;
end;%if flag_calc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for index_calc=0:3-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
