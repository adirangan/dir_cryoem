function ...
[ ...
 parameter ...
,ssnll_M_ ...
,ssnll ...
,S_k_p_wkS__ ...
,dvol_ssnll_M_ ...
,dvol_ssnll ...
,dvol_S_k_p_wkS__ ...
,dvol_dvol_ssnll ...
,dtau_ssnll_M3__ ...
,dtau_ssnll ...
,dtau_S_k_p_wkS3___ ...
,dtau_dvol_ssnll_M3__ ...
,dtau_dvol_ssnll ...
,dtau_dvol_S_k_p_wkS3___ ...
,dtau_dtau_ssnll_M33___ ...
,dtau_dtau_ssnll ...
,dtau_dtau_S_k_p_wkS33____ ...
] = ...
ssnll_from_plane_wave_expansion_2( ...
 parameter ...
,n_source_a ...
,v_source_a_ ...
,delta_a_c__ ...
,n_source_dvol_a ...
,v_source_dvol_a_ ...
,delta_dvol_a_c__ ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_max ...
,weight_2d_wk_ ...
,n_CTF ...
,index_nCTF_from_nM_ ...
,CTF_phi_C_ ...
,n_eta ...
,index_neta_from_nM_ ...
,eta_k_p_wke__ ...
,n_M ...
,weight_imagecount_M_ ...
,euler_polar_a_M_ ...
,euler_azimu_b_M_ ...
,euler_gamma_z_M_ ...
,dtau_euler_polar_a_M_ ...
,dtau_euler_azimu_b_M_ ...
,dtau_euler_gamma_z_M_ ...
,n_source_b ...
,v_source_b_ ...
,delta_b_c__ ...
,fromb_polar_a_M_ ...
,fromb_azimu_b_M_ ...
,fromb_gamma_z_M_ ...
,M_phi_M_ ...
,n_S ...
,viewing_polar_a_S_ ...
,viewing_azimu_b_S_ ...
,viewing_gamma_z_S_ ...
);

str_thisfunction = 'ssnll_from_plane_wave_expansion_2';

if nargin<1;
disp(sprintf(' %% testing %s',str_thisfunction));
test_ssnll_from_plane_wave_expansion_2;
disp('returning'); return;
end;%if nargin<1;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_source_a=[]; end; na=na+1;
if (nargin<1+na); v_source_a_=[]; end; na=na+1;
if (nargin<1+na); delta_a_c__=[]; end; na=na+1;
if (nargin<1+na); n_source_dvol_a=[]; end; na=na+1;
if (nargin<1+na); v_source_dvol_a_=[]; end; na=na+1;
if (nargin<1+na); delta_dvol_a_c__=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); k_p_r_=[]; end; na=na+1;
if (nargin<1+na); k_p_r_max=[]; end; na=na+1;
if (nargin<1+na); n_w_max=[]; end; na=na+1;
if (nargin<1+na); weight_2d_wk_=[]; end; na=na+1;
if (nargin<1+na); n_CTF=[]; end; na=na+1;
if (nargin<1+na); index_nCTF_from_nM_=[]; end; na=na+1;
if (nargin<1+na); CTF_phi_C_=[]; end; na=na+1;
if (nargin<1+na); n_eta=[]; end; na=na+1;
if (nargin<1+na); index_neta_from_nM_=[]; end; na=na+1;
if (nargin<1+na); eta_k_p_wke__=[]; end; na=na+1;
if (nargin<1+na); n_M=[]; end; na=na+1;
if (nargin<1+na); weight_imagecount_M_=[]; end; na=na+1;
if (nargin<1+na); euler_polar_a_M_=[]; end; na=na+1;
if (nargin<1+na); euler_azimu_b_M_=[]; end; na=na+1;
if (nargin<1+na); euler_gamma_z_M_=[]; end; na=na+1;
if (nargin<1+na); dtau_euler_polar_a_M_=[]; end; na=na+1;
if (nargin<1+na); dtau_euler_azimu_b_M_=[]; end; na=na+1;
if (nargin<1+na); dtau_euler_gamma_z_M_=[]; end; na=na+1;
if (nargin<1+na); n_source_b=[]; end; na=na+1;
if (nargin<1+na); v_source_b_=[]; end; na=na+1;
if (nargin<1+na); delta_b_c__=[]; end; na=na+1;
if (nargin<1+na); fromb_polar_a_M_=[]; end; na=na+1;
if (nargin<1+na); fromb_azimu_b_M_=[]; end; na=na+1;
if (nargin<1+na); fromb_gamma_z_M_=[]; end; na=na+1;
if (nargin<1+na); M_phi_M_=[]; end; na=na+1;
if (nargin<1+na); n_S=[]; end; na=na+1;
if (nargin<1+na); viewing_polar_a_S_=[]; end; na=na+1;
if (nargin<1+na); viewing_azimu_b_S_=[]; end; na=na+1;
if (nargin<1+na); viewing_gamma_z_S_=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
flag_verbose=parameter.flag_verbose;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

if isempty(weight_imagecount_M_); weight_imagecount_M_ = ones(n_M,1); end;

P__ = ...
[1 , 0 , 0 ; ...
 0 , 1 , 0 ; ...
];
%%%%%%%%;
% define rotations. ;
%%%%%%%%;
Rz = @(azimu_b) ...
[ +cos(azimu_b) -sin(azimu_b) 0 ; ...
  +sin(azimu_b) +cos(azimu_b) 0 ; ...
   0             0            1 ; ...
] ;
%%%%%%%%;
Ry = @(polar_a) ...
[ +cos(polar_a) 0 +sin(polar_a) ; ...
   0            1  0            ; ...
  -sin(polar_a) 0 +cos(polar_a) ; ...
];
%%%%%%%%;
dRz = @(azimu_b) ...
[ -sin(azimu_b) -cos(azimu_b) 0 ; ...
  +cos(azimu_b) -sin(azimu_b) 0 ; ...
   0             0            0 ; ...
] ;
%%%%%%%%;
dRy = @(polar_a) ...
[ -sin(polar_a) 0 +cos(polar_a) ; ...
   0            0  0            ; ...
  -cos(polar_a) 0 -sin(polar_a) ; ...
];
%%%%%%%%;
ddRz = @(azimu_b) ...
[ -cos(azimu_b) +sin(azimu_b) 0 ; ...
  -sin(azimu_b) -cos(azimu_b) 0 ; ...
   0             0            0 ; ...
] ;
%%%%%%%%;
ddRy = @(polar_a) ...
[ -cos(polar_a) 0 -sin(polar_a) ; ...
   0            0  0            ; ...
  +sin(polar_a) 0 -cos(polar_a) ; ...
];
%%%%%%%%;
flag_check=0;
if flag_check;
tmp_da = 1e-4;
tmp_a = 1*pi*rand();
Ry_mid__ = Ry(tmp_a); Ry_pos__ = Ry(tmp_a + tmp_da); Ry_neg__ = Ry(tmp_a - tmp_da);
dRy_dif__ = (Ry_pos__ - Ry_neg__)/max(1e-12,2*tmp_da);
dRy_mid__ = dRy(tmp_a);
fnorm_disp(flag_verbose,'dRy_mid__',dRy_mid__,'dRy_dif__',dRy_dif__);
ddRy_dif__ = (Ry_pos__ -2*Ry_mid__ + Ry_neg__)/max(1e-12,tmp_da.^2);
ddRy_mid__ = ddRy(tmp_a);
fnorm_disp(flag_verbose,'ddRy_mid__',ddRy_mid__,'ddRy_dif__',ddRy_dif__,' %<-- should be <1e-6');
tmp_db = 1e-5;
tmp_b = 2*pi*rand();
Rz_mid__ = Rz(tmp_b); Rz_pos__ = Rz(tmp_b + tmp_db); Rz_neg__ = Rz(tmp_b - tmp_db);
dRz_dif__ = (Rz_pos__ - Rz_neg__)/max(1e-12,2*tmp_db);
dRz_mid__ = dRz(tmp_b);
fnorm_disp(flag_verbose,'dRz_mid__',dRz_mid__,'dRz_dif__',dRz_dif__);
ddRz_dif__ = (Rz_pos__ -2*Rz_mid__ + Rz_neg__)/max(1e-12,tmp_db.^2);
ddRz_mid__ = ddRz(tmp_b);
fnorm_disp(flag_verbose,'ddRz_mid__',ddRz_mid__,'ddRz_dif__',ddRz_dif__,' %<-- should be <1e-6');
end;%if flag_check;
%%%%%%%%;

n_w_sum = n_w_max*n_k_p_r;
psi_z_ = transpose((2*pi)*[0:n_w_max-1]/max(1,n_w_max));
k_c_0_wk_ = reshape(cos(psi_z_)*reshape(k_p_r_,[1,n_k_p_r]),[n_w_sum,1]);
k_c_1_wk_ = reshape(sin(psi_z_)*reshape(k_p_r_,[1,n_k_p_r]),[n_w_sum,1]);

if  isempty(v_source_a_); v_source_a_ = ones(n_source_a,1); end;
if  isempty(v_source_b_); v_source_b_ = ones(n_source_b,1); end;
if  isempty(v_source_dvol_a_); v_source_dvol_a_ = ones(n_source_dvol_a,1); end;
if  isempty(euler_gamma_z_M_); euler_gamma_z_M_ = zeros(n_M,1); end;
if  isempty(viewing_gamma_z_S_); viewing_gamma_z_S_ = zeros(n_S,1); end;
if  isempty(dtau_euler_polar_a_M_); dtau_euler_polar_a_M_ = zeros(n_M,1); end;
if  isempty(dtau_euler_azimu_b_M_); dtau_euler_azimu_b_M_ = zeros(n_M,1); end;
if  isempty(dtau_euler_gamma_z_M_); dtau_euler_gamma_z_M_ = zeros(n_M,1); end;

flag_ssnll = 1;
flag_dvol_ssnll = 0; if (nargout>=1+4); flag_dvol_ssnll = 1; end;
flag_dtau_ssnll = 0; if (nargout>=1+8); flag_dtau_ssnll = 1; end;
flag_dvol_dvol_ssnll = 0; if (nargout>=1+7); flag_dvol_dvol_ssnll = 1; end;
flag_dtau_dvol_ssnll = 0; if (nargout>=1+11); flag_dtau_dvol_ssnll = 1; end;
flag_dtau_dtau_ssnll = 0; if (nargout>=1+14); flag_dtau_dtau_ssnll = 1; end;

flag_S = 0; if (~isempty(n_S) & nargout>=1+3); flag_S = 1; end;
flag_dvol_S = 0; if (~isempty(n_S) & nargout>=1+6); flag_dvol_S = 1; end;
flag_dtau_S = 0; if (~isempty(n_S) & nargout>=1+10); flag_dtau_S = 1; end;
flag_dtau_dvol_S = 0; if (~isempty(n_S) & nargout>=1+13); flag_dtau_dvol_S = 1; end;
flag_dtau_dtau_S = 0; if (~isempty(n_S) & nargout>=1+16); flag_dtau_dtau_S = 1; end;

if ~flag_dvol_S;
n_source_dvol_a = 1;
v_source_dvol_a_ = zeros(n_source_dvol_a,1);
delta_dvol_a_c__ = ones(3,n_source_dvol_a);
end;%if ~flag_dvol_S;

if flag_S;
tmp_t = tic();
S_k_p_wkS__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
tmp_viewing_polar_a = +viewing_polar_a_S_(1+nS);
tmp_viewing_azimu_b = +viewing_azimu_b_S_(1+nS);
tmp_viewing_gamma_z = -viewing_gamma_z_S_(1+nS);
tmp_R_a__ = Rz(-tmp_viewing_gamma_z)*Ry(-tmp_viewing_polar_a)*Rz(-tmp_viewing_azimu_b);
U_k_p_wk_ = zeros(n_w_sum,1);
for nsource_a=0:n_source_a-1;
tmp_v_source_a = v_source_a_(1+nsource_a);
tmp_delta_U_ = tmp_R_a__*delta_a_c__(:,1+nsource_a);
U_k_p_wk_ = U_k_p_wk_ + tmp_v_source_a*exp(+i*2*pi*(k_c_0_wk_*tmp_delta_U_(1+0) + k_c_1_wk_*tmp_delta_U_(1+1)));
end;%for nsource_a=0:n_source_a-1;
S_k_p_wkS__(:,1+nS) = U_k_p_wk_;
end;%for nS=0:n_S-1;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% S_k_p_wkS__: %0.6fs',tmp_t)); end;
end;%if flag_S;

if flag_dtau_S;
tmp_t = tic();
dtau_S_k_p_wkS3___ = zeros(n_w_sum,n_S,3);
for nS=0:n_S-1;
tmp_viewing_polar_a = +viewing_polar_a_S_(1+nS);
tmp_viewing_azimu_b = +viewing_azimu_b_S_(1+nS);
tmp_viewing_gamma_z = -viewing_gamma_z_S_(1+nS);
tmp_R_a__ = Rz(-tmp_viewing_gamma_z)*Ry(-tmp_viewing_polar_a)*Rz(-tmp_viewing_azimu_b);
tmp_da_R_a__ = -Rz(-tmp_viewing_gamma_z)*dRy(-tmp_viewing_polar_a)*Rz(-tmp_viewing_azimu_b);
tmp_db_R_a__ = -Rz(-tmp_viewing_gamma_z)*Ry(-tmp_viewing_polar_a)*dRz(-tmp_viewing_azimu_b);
tmp_dc_R_a__ = +dRz(-tmp_viewing_gamma_z)*Ry(-tmp_viewing_polar_a)*Rz(-tmp_viewing_azimu_b);
U_k_p_wk_ = zeros(n_w_sum,1);
da_U_k_p_wk_ = zeros(n_w_sum,1);
db_U_k_p_wk_ = zeros(n_w_sum,1);
dc_U_k_p_wk_ = zeros(n_w_sum,1);
for nsource_a=0:n_source_a-1;
tmp_v_source_a = v_source_a_(1+nsource_a);
tmp_delta_U_ = tmp_R_a__*delta_a_c__(:,1+nsource_a);
tmp_da_delta_U_ = tmp_da_R_a__*delta_a_c__(:,1+nsource_a);
tmp_db_delta_U_ = tmp_db_R_a__*delta_a_c__(:,1+nsource_a);
tmp_dc_delta_U_ = tmp_dc_R_a__*delta_a_c__(:,1+nsource_a);
planewave_wk_ = exp(+i*2*pi*(k_c_0_wk_*tmp_delta_U_(1+0) + k_c_1_wk_*tmp_delta_U_(1+1)));
U_k_p_wk_ = U_k_p_wk_ + tmp_v_source_a*planewave_wk_;
da_U_k_p_wk_ = da_U_k_p_wk_ + tmp_v_source_a*planewave_wk_.*(+i*2*pi*(k_c_0_wk_*tmp_da_delta_U_(1+0) + k_c_1_wk_*tmp_da_delta_U_(1+1)));
db_U_k_p_wk_ = db_U_k_p_wk_ + tmp_v_source_a*planewave_wk_.*(+i*2*pi*(k_c_0_wk_*tmp_db_delta_U_(1+0) + k_c_1_wk_*tmp_db_delta_U_(1+1)));
dc_U_k_p_wk_ = dc_U_k_p_wk_ + tmp_v_source_a*planewave_wk_.*(+i*2*pi*(k_c_0_wk_*tmp_dc_delta_U_(1+0) + k_c_1_wk_*tmp_dc_delta_U_(1+1)));
end;%for nsource_a=0:n_source_a-1;
dtau_S_k_p_wkS3___(:,1+nS,1+0) = da_U_k_p_wk_;
dtau_S_k_p_wkS3___(:,1+nS,1+1) = db_U_k_p_wk_;
dtau_S_k_p_wkS3___(:,1+nS,1+2) = dc_U_k_p_wk_;
end;%for nS=0:n_S-1;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% dtau_S_k_p_wkS3___: %0.6fs',tmp_t)); end;
end;%if flag_dtau_S;

if flag_dtau_dtau_S;
tmp_t = tic();
dtau_dtau_S_k_p_wkS33____ = zeros(n_w_sum,n_S,3,3);
for nS=0:n_S-1;
tmp_viewing_polar_a = +viewing_polar_a_S_(1+nS);
tmp_viewing_azimu_b = +viewing_azimu_b_S_(1+nS);
tmp_viewing_gamma_z = -viewing_gamma_z_S_(1+nS);
tmp_R_a__ = Rz(-tmp_viewing_gamma_z)*Ry(-tmp_viewing_polar_a)*Rz(-tmp_viewing_azimu_b);
tmp_da_R_a__ = -Rz(-tmp_viewing_gamma_z)*dRy(-tmp_viewing_polar_a)*Rz(-tmp_viewing_azimu_b);
tmp_db_R_a__ = -Rz(-tmp_viewing_gamma_z)*Ry(-tmp_viewing_polar_a)*dRz(-tmp_viewing_azimu_b);
tmp_dc_R_a__ = +dRz(-tmp_viewing_gamma_z)*Ry(-tmp_viewing_polar_a)*Rz(-tmp_viewing_azimu_b);
tmp_dada_R_a__ = +Rz(-tmp_viewing_gamma_z)*ddRy(-tmp_viewing_polar_a)*Rz(-tmp_viewing_azimu_b);
tmp_dadb_R_a__ = +Rz(-tmp_viewing_gamma_z)*dRy(-tmp_viewing_polar_a)*dRz(-tmp_viewing_azimu_b);
tmp_dadc_R_a__ = -dRz(-tmp_viewing_gamma_z)*dRy(-tmp_viewing_polar_a)*Rz(-tmp_viewing_azimu_b);
tmp_dbda_R_a__ = +Rz(-tmp_viewing_gamma_z)*dRy(-tmp_viewing_polar_a)*dRz(-tmp_viewing_azimu_b);
tmp_dbdb_R_a__ = +Rz(-tmp_viewing_gamma_z)*Ry(-tmp_viewing_polar_a)*ddRz(-tmp_viewing_azimu_b);
tmp_dbdc_R_a__ = -dRz(-tmp_viewing_gamma_z)*Ry(-tmp_viewing_polar_a)*dRz(-tmp_viewing_azimu_b);
tmp_dcda_R_a__ = -dRz(-tmp_viewing_gamma_z)*dRy(-tmp_viewing_polar_a)*Rz(-tmp_viewing_azimu_b);
tmp_dcdb_R_a__ = -dRz(-tmp_viewing_gamma_z)*Ry(-tmp_viewing_polar_a)*dRz(-tmp_viewing_azimu_b);
tmp_dcdc_R_a__ = +ddRz(-tmp_viewing_gamma_z)*Ry(-tmp_viewing_polar_a)*Rz(-tmp_viewing_azimu_b);
U_k_p_wk_ = zeros(n_w_sum,1);
da_U_k_p_wk_ = zeros(n_w_sum,1);
db_U_k_p_wk_ = zeros(n_w_sum,1);
dc_U_k_p_wk_ = zeros(n_w_sum,1);
dada_U_k_p_wk_ = zeros(n_w_sum,1);
dadb_U_k_p_wk_ = zeros(n_w_sum,1);
dadc_U_k_p_wk_ = zeros(n_w_sum,1);
dbda_U_k_p_wk_ = zeros(n_w_sum,1);
dbdb_U_k_p_wk_ = zeros(n_w_sum,1);
dbdc_U_k_p_wk_ = zeros(n_w_sum,1);
dcda_U_k_p_wk_ = zeros(n_w_sum,1);
dcdb_U_k_p_wk_ = zeros(n_w_sum,1);
dcdc_U_k_p_wk_ = zeros(n_w_sum,1);
for nsource_a=0:n_source_a-1;
tmp_v_source_a = v_source_a_(1+nsource_a);
tmp_delta_U_ = tmp_R_a__*delta_a_c__(:,1+nsource_a);
tmp_da_delta_U_ = tmp_da_R_a__*delta_a_c__(:,1+nsource_a);
tmp_db_delta_U_ = tmp_db_R_a__*delta_a_c__(:,1+nsource_a);
tmp_dc_delta_U_ = tmp_dc_R_a__*delta_a_c__(:,1+nsource_a);
tmp_dada_delta_U_ = tmp_dada_R_a__*delta_a_c__(:,1+nsource_a);
tmp_dadb_delta_U_ = tmp_dadb_R_a__*delta_a_c__(:,1+nsource_a);
tmp_dadc_delta_U_ = tmp_dadc_R_a__*delta_a_c__(:,1+nsource_a);
tmp_dbda_delta_U_ = tmp_dbda_R_a__*delta_a_c__(:,1+nsource_a);
tmp_dbdb_delta_U_ = tmp_dbdb_R_a__*delta_a_c__(:,1+nsource_a);
tmp_dbdc_delta_U_ = tmp_dbdc_R_a__*delta_a_c__(:,1+nsource_a);
tmp_dcda_delta_U_ = tmp_dcda_R_a__*delta_a_c__(:,1+nsource_a);
tmp_dcdb_delta_U_ = tmp_dcdb_R_a__*delta_a_c__(:,1+nsource_a);
tmp_dcdc_delta_U_ = tmp_dcdc_R_a__*delta_a_c__(:,1+nsource_a);
planewave_wk_ = exp(+i*2*pi*(k_c_0_wk_*tmp_delta_U_(1+0) + k_c_1_wk_*tmp_delta_U_(1+1)));
da_factor_wk_ = (+i*2*pi*(k_c_0_wk_*tmp_da_delta_U_(1+0) + k_c_1_wk_*tmp_da_delta_U_(1+1)));
db_factor_wk_ = (+i*2*pi*(k_c_0_wk_*tmp_db_delta_U_(1+0) + k_c_1_wk_*tmp_db_delta_U_(1+1)));
dc_factor_wk_ = (+i*2*pi*(k_c_0_wk_*tmp_dc_delta_U_(1+0) + k_c_1_wk_*tmp_dc_delta_U_(1+1)));
dada_factor_wk_ = (+i*2*pi*(k_c_0_wk_*tmp_dada_delta_U_(1+0) + k_c_1_wk_*tmp_dada_delta_U_(1+1)));
dadb_factor_wk_ = (+i*2*pi*(k_c_0_wk_*tmp_dadb_delta_U_(1+0) + k_c_1_wk_*tmp_dadb_delta_U_(1+1)));
dadc_factor_wk_ = (+i*2*pi*(k_c_0_wk_*tmp_dadc_delta_U_(1+0) + k_c_1_wk_*tmp_dadc_delta_U_(1+1)));
dbda_factor_wk_ = (+i*2*pi*(k_c_0_wk_*tmp_dbda_delta_U_(1+0) + k_c_1_wk_*tmp_dbda_delta_U_(1+1)));
dbdb_factor_wk_ = (+i*2*pi*(k_c_0_wk_*tmp_dbdb_delta_U_(1+0) + k_c_1_wk_*tmp_dbdb_delta_U_(1+1)));
dbdc_factor_wk_ = (+i*2*pi*(k_c_0_wk_*tmp_dbdc_delta_U_(1+0) + k_c_1_wk_*tmp_dbdc_delta_U_(1+1)));
dcda_factor_wk_ = (+i*2*pi*(k_c_0_wk_*tmp_dcda_delta_U_(1+0) + k_c_1_wk_*tmp_dcda_delta_U_(1+1)));
dcdb_factor_wk_ = (+i*2*pi*(k_c_0_wk_*tmp_dcdb_delta_U_(1+0) + k_c_1_wk_*tmp_dcdb_delta_U_(1+1)));
dcdc_factor_wk_ = (+i*2*pi*(k_c_0_wk_*tmp_dcdc_delta_U_(1+0) + k_c_1_wk_*tmp_dcdc_delta_U_(1+1)));
dada_U_k_p_wk_ = dada_U_k_p_wk_ + tmp_v_source_a*planewave_wk_.*(da_factor_wk_.*da_factor_wk_ + dada_factor_wk_);
dadb_U_k_p_wk_ = dadb_U_k_p_wk_ + tmp_v_source_a*planewave_wk_.*(da_factor_wk_.*db_factor_wk_ + dadb_factor_wk_);
dadc_U_k_p_wk_ = dadc_U_k_p_wk_ + tmp_v_source_a*planewave_wk_.*(da_factor_wk_.*dc_factor_wk_ + dadc_factor_wk_);
dbda_U_k_p_wk_ = dbda_U_k_p_wk_ + tmp_v_source_a*planewave_wk_.*(db_factor_wk_.*da_factor_wk_ + dbda_factor_wk_);
dbdb_U_k_p_wk_ = dbdb_U_k_p_wk_ + tmp_v_source_a*planewave_wk_.*(db_factor_wk_.*db_factor_wk_ + dbdb_factor_wk_);
dbdc_U_k_p_wk_ = dbdc_U_k_p_wk_ + tmp_v_source_a*planewave_wk_.*(db_factor_wk_.*dc_factor_wk_ + dbdc_factor_wk_);
dcda_U_k_p_wk_ = dcda_U_k_p_wk_ + tmp_v_source_a*planewave_wk_.*(dc_factor_wk_.*da_factor_wk_ + dcda_factor_wk_);
dcdb_U_k_p_wk_ = dcdb_U_k_p_wk_ + tmp_v_source_a*planewave_wk_.*(dc_factor_wk_.*db_factor_wk_ + dcdb_factor_wk_);
dcdc_U_k_p_wk_ = dcdc_U_k_p_wk_ + tmp_v_source_a*planewave_wk_.*(dc_factor_wk_.*dc_factor_wk_ + dcdc_factor_wk_);
end;%for nsource_a=0:n_source_a-1;
dtau_dtau_S_k_p_wkS33____(:,1+nS,1+0,1+0) = dada_U_k_p_wk_;
dtau_dtau_S_k_p_wkS33____(:,1+nS,1+0,1+1) = dadb_U_k_p_wk_;
dtau_dtau_S_k_p_wkS33____(:,1+nS,1+0,1+2) = dadc_U_k_p_wk_;
dtau_dtau_S_k_p_wkS33____(:,1+nS,1+1,1+0) = dbda_U_k_p_wk_;
dtau_dtau_S_k_p_wkS33____(:,1+nS,1+1,1+1) = dbdb_U_k_p_wk_;
dtau_dtau_S_k_p_wkS33____(:,1+nS,1+1,1+2) = dbdc_U_k_p_wk_;
dtau_dtau_S_k_p_wkS33____(:,1+nS,1+2,1+0) = dcda_U_k_p_wk_;
dtau_dtau_S_k_p_wkS33____(:,1+nS,1+2,1+1) = dcdb_U_k_p_wk_;
dtau_dtau_S_k_p_wkS33____(:,1+nS,1+2,1+2) = dcdc_U_k_p_wk_;
end;%for nS=0:n_S-1;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% dtau_dtau_S_k_p_wkS33____: %0.6fs',tmp_t)); end;
end;%if flag_dtau_dtau_S;

if flag_dvol_S;
tmp_t = tic();
dvol_S_k_p_wkS__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
tmp_viewing_polar_a = +viewing_polar_a_S_(1+nS);
tmp_viewing_azimu_b = +viewing_azimu_b_S_(1+nS);
tmp_viewing_gamma_z = -viewing_gamma_z_S_(1+nS);
tmp_R_a__ = Rz(-tmp_viewing_gamma_z)*Ry(-tmp_viewing_polar_a)*Rz(-tmp_viewing_azimu_b);
dvol_U_k_p_wk_ = zeros(n_w_sum,1);
for nsource_dvol_a=0:n_source_dvol_a-1;
tmp_v_source_dvol_a = v_source_dvol_a_(1+nsource_dvol_a);
tmp_delta_dvol_U_ = tmp_R_a__*delta_dvol_a_c__(:,1+nsource_dvol_a);
dvol_U_k_p_wk_ = dvol_U_k_p_wk_ + tmp_v_source_dvol_a*exp(+i*2*pi*(k_c_0_wk_*tmp_delta_dvol_U_(1+0) + k_c_1_wk_*tmp_delta_dvol_U_(1+1)));
end;%for nsource_dvol_a=0:n_source_dvol_a-1;
dvol_S_k_p_wkS__(:,1+nS) = dvol_U_k_p_wk_;
end;%for nS=0:n_S-1;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% dvol_S_k_p_wkS__: %0.6fs',tmp_t)); end;
end;%if flag_dvol_S;

if flag_dtau_dvol_S;
tmp_t = tic();
dtau_dvol_S_k_p_wkS3___ = zeros(n_w_sum,n_S,3);
for nS=0:n_S-1;
tmp_viewing_polar_a = +viewing_polar_a_S_(1+nS);
tmp_viewing_azimu_b = +viewing_azimu_b_S_(1+nS);
tmp_viewing_gamma_z = -viewing_gamma_z_S_(1+nS);
tmp_R_a__ = Rz(-tmp_viewing_gamma_z)*Ry(-tmp_viewing_polar_a)*Rz(-tmp_viewing_azimu_b);
tmp_da_R_a__ = -Rz(-tmp_viewing_gamma_z)*dRy(-tmp_viewing_polar_a)*Rz(-tmp_viewing_azimu_b);
tmp_db_R_a__ = -Rz(-tmp_viewing_gamma_z)*Ry(-tmp_viewing_polar_a)*dRz(-tmp_viewing_azimu_b);
tmp_dc_R_a__ = +dRz(-tmp_viewing_gamma_z)*Ry(-tmp_viewing_polar_a)*Rz(-tmp_viewing_azimu_b);
dvol_U_k_p_wk_ = zeros(n_w_sum,1);
da_dvol_U_k_p_wk_ = zeros(n_w_sum,1);
db_dvol_U_k_p_wk_ = zeros(n_w_sum,1);
dc_dvol_U_k_p_wk_ = zeros(n_w_sum,1);
for nsource_dvol_a=0:n_source_dvol_a-1;
tmp_v_source_dvol_a = v_source_dvol_a_(1+nsource_dvol_a);
tmp_delta_dvol_U_ = tmp_R_a__*delta_dvol_a_c__(:,1+nsource_dvol_a);
tmp_da_delta_dvol_U_ = tmp_da_R_a__*delta_dvol_a_c__(:,1+nsource_dvol_a);
tmp_db_delta_dvol_U_ = tmp_db_R_a__*delta_dvol_a_c__(:,1+nsource_dvol_a);
tmp_dc_delta_dvol_U_ = tmp_dc_R_a__*delta_dvol_a_c__(:,1+nsource_dvol_a);
dvol_planewave_wk_ = exp(+i*2*pi*(k_c_0_wk_*tmp_delta_dvol_U_(1+0) + k_c_1_wk_*tmp_delta_dvol_U_(1+1)));
dvol_U_k_p_wk_ = dvol_U_k_p_wk_ + tmp_v_source_dvol_a*dvol_planewave_wk_;
da_dvol_U_k_p_wk_ = da_dvol_U_k_p_wk_ + tmp_v_source_dvol_a*dvol_planewave_wk_.*(+i*2*pi*(k_c_0_wk_*tmp_da_delta_dvol_U_(1+0) + k_c_1_wk_*tmp_da_delta_dvol_U_(1+1)));
db_dvol_U_k_p_wk_ = db_dvol_U_k_p_wk_ + tmp_v_source_dvol_a*dvol_planewave_wk_.*(+i*2*pi*(k_c_0_wk_*tmp_db_delta_dvol_U_(1+0) + k_c_1_wk_*tmp_db_delta_dvol_U_(1+1)));
dc_dvol_U_k_p_wk_ = dc_dvol_U_k_p_wk_ + tmp_v_source_dvol_a*dvol_planewave_wk_.*(+i*2*pi*(k_c_0_wk_*tmp_dc_delta_dvol_U_(1+0) + k_c_1_wk_*tmp_dc_delta_dvol_U_(1+1)));
end;%for nsource_dvol_a=0:n_source_dvol_a-1;
dtau_dvol_S_k_p_wkS3___(:,1+nS,1+0) = da_dvol_U_k_p_wk_;
dtau_dvol_S_k_p_wkS3___(:,1+nS,1+1) = db_dvol_U_k_p_wk_;
dtau_dvol_S_k_p_wkS3___(:,1+nS,1+2) = dc_dvol_U_k_p_wk_;
end;%for nS=0:n_S-1;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% dtau_dvol_S_k_p_wkS3___: %0.6fs',tmp_t)); end;
end;%if flag_dtau_dvol_S;

%%%%%%%%%%%%%%%%;
tmp_t = tic();
ssnll_M_ = zeros(n_M,1);
ssnll = 0.0d0;
dvol_ssnll_M_ = zeros(n_M,1);
dvol_ssnll = 0.0d0;
dvol_dvol_ssnll = 0.0d0;
dtau_ssnll_M3__ = zeros(n_M,3);
dtau_ssnll = 0.0d0;
dtau_dvol_ssnll_M3__ = zeros(n_M,3);
dtau_dvol_ssnll = 0.0d0;
dtau_dtau_ssnll_M33___ = zeros(n_M,3,3);
dtau_dtau_ssnll = 0.0d0;
%%%%;
da_ssnll_M_ = zeros(n_M,1);
db_ssnll_M_ = zeros(n_M,1);
dc_ssnll_M_ = zeros(n_M,1);
dada_ssnll_M_ = zeros(n_M,1);
dadb_ssnll_M_ = zeros(n_M,1);
dadc_ssnll_M_ = zeros(n_M,1);
dbdb_ssnll_M_ = zeros(n_M,1);
dbdc_ssnll_M_ = zeros(n_M,1);
dcdc_ssnll_M_ = zeros(n_M,1);
dvol_ssnll_M_ = zeros(n_M,1);
da_dvol_ssnll_M_ = zeros(n_M,1);
db_dvol_ssnll_M_ = zeros(n_M,1);
dc_dvol_ssnll_M_ = zeros(n_M,1);
dvol_dvol_ssnll_M_ = zeros(n_M,1);
%%%%%%%%%%%%%%%%;
for nM=0:n_M-1;
%%%%%%%%%%%%%%%%;
M_phi = M_phi_M_(1+nM);
nCTF = index_nCTF_from_nM_(1+nM);
CTF_phi = CTF_phi_C_(1+nCTF);
neta = index_neta_from_nM_(1+nM);
eta_k_p_wk_ = eta_k_p_wke__(:,1+neta);
%%%%;
% the template U_k_p_wk_ (from volume a) is associated with [ euler_azimu_b_M_(1+nM) , euler_polar_a_M_(1+nM) , euler_gamma_z_M_(1+nM) ] ;
% the template V_k_p_wk_ (from volume b) is associated with [ fromb_azimu_b_M_(1+nM) , fromb_polar_a_M_(1+nM) , fromb_gamma_z_M_(1+nM) ] ;
% We assume the CTF*S is 2*k_p_r_wk_.*cos(k_p_w_wk_ - CTF_phi_M_(1+nM)) .* U_k_p_wk_. ;
% We assume the image is 2*k_p_r_wk_.*cos(k_p_w_wk_ -   M_phi_M_(1+nM)) .* V_k_p_wk_. ;
%%%%;
tmp_euler_polar_a = +euler_polar_a_M_(1+nM);
tmp_euler_azimu_b = +euler_azimu_b_M_(1+nM);
tmp_euler_gamma_z = -euler_gamma_z_M_(1+nM);
tmp_R_a__ = Rz(-tmp_euler_gamma_z)*Ry(-tmp_euler_polar_a)*Rz(-tmp_euler_azimu_b);
tmp_da_R_a__ = -Rz(-tmp_euler_gamma_z)*dRy(-tmp_euler_polar_a)*Rz(-tmp_euler_azimu_b);
tmp_db_R_a__ = -Rz(-tmp_euler_gamma_z)*Ry(-tmp_euler_polar_a)*dRz(-tmp_euler_azimu_b);
tmp_dc_R_a__ = +dRz(-tmp_euler_gamma_z)*Ry(-tmp_euler_polar_a)*Rz(-tmp_euler_azimu_b);
tmp_dada_R_a__ = +Rz(-tmp_euler_gamma_z)*ddRy(-tmp_euler_polar_a)*Rz(-tmp_euler_azimu_b);
tmp_dadb_R_a__ = +Rz(-tmp_euler_gamma_z)*dRy(-tmp_euler_polar_a)*dRz(-tmp_euler_azimu_b);
tmp_dadc_R_a__ = -dRz(-tmp_euler_gamma_z)*dRy(-tmp_euler_polar_a)*Rz(-tmp_euler_azimu_b);
%tmp_dbda_R_a__ = +Rz(-tmp_euler_gamma_z)*dRy(-tmp_euler_polar_a)*dRz(-tmp_euler_azimu_b);
tmp_dbdb_R_a__ = +Rz(-tmp_euler_gamma_z)*Ry(-tmp_euler_polar_a)*ddRz(-tmp_euler_azimu_b);
tmp_dbdc_R_a__ = -dRz(-tmp_euler_gamma_z)*Ry(-tmp_euler_polar_a)*dRz(-tmp_euler_azimu_b);
%tmp_dcda_R_a__ = -dRz(-tmp_euler_gamma_z)*dRy(-tmp_euler_polar_a)*Rz(-tmp_euler_azimu_b);
%tmp_dcdb_R_a__ = -dRz(-tmp_euler_gamma_z)*Ry(-tmp_euler_polar_a)*dRz(-tmp_euler_azimu_b);
tmp_dcdc_R_a__ = +ddRz(-tmp_euler_gamma_z)*Ry(-tmp_euler_polar_a)*Rz(-tmp_euler_azimu_b);
%%%%;
tmp_delta_U_2a__ = zeros(2,n_source_a);
tmp_da_delta_U_2a__ = zeros(2,n_source_a);
tmp_db_delta_U_2a__ = zeros(2,n_source_a);
tmp_dc_delta_U_2a__ = zeros(2,n_source_a);
tmp_dada_delta_U_2a__ = zeros(2,n_source_a);
tmp_dadb_delta_U_2a__ = zeros(2,n_source_a);
tmp_dadc_delta_U_2a__ = zeros(2,n_source_a);
%tmp_dbda_delta_U_2a__ = zeros(2,n_source_a);
tmp_dbdb_delta_U_2a__ = zeros(2,n_source_a);
tmp_dbdc_delta_U_2a__ = zeros(2,n_source_a);
%tmp_dcda_delta_U_2a__ = zeros(2,n_source_a);
%tmp_dcdb_delta_U_2a__ = zeros(2,n_source_a);
tmp_dcdc_delta_U_2a__ = zeros(2,n_source_a);
for nsource_a=0:n_source_a-1;
tmp_v_source_a = v_source_a_(1+nsource_a);
tmp_delta_U_2a__(:,1+nsource_a) = P__ * tmp_R_a__*delta_a_c__(:,1+nsource_a);
tmp_da_delta_U_2a__(:,1+nsource_a) = P__ * tmp_da_R_a__*delta_a_c__(:,1+nsource_a);
tmp_db_delta_U_2a__(:,1+nsource_a) = P__ * tmp_db_R_a__*delta_a_c__(:,1+nsource_a);
tmp_dc_delta_U_2a__(:,1+nsource_a) = P__ * tmp_dc_R_a__*delta_a_c__(:,1+nsource_a);
tmp_dada_delta_U_2a__(:,1+nsource_a) = P__ * tmp_dada_R_a__*delta_a_c__(:,1+nsource_a);
tmp_dadb_delta_U_2a__(:,1+nsource_a) = P__ * tmp_dadb_R_a__*delta_a_c__(:,1+nsource_a);
tmp_dadc_delta_U_2a__(:,1+nsource_a) = P__ * tmp_dadc_R_a__*delta_a_c__(:,1+nsource_a);
%tmp_dbda_delta_U_2a__(:,1+nsource_a) = P__ * tmp_dbda_R_a__*delta_a_c__(:,1+nsource_a);
tmp_dbdb_delta_U_2a__(:,1+nsource_a) = P__ * tmp_dbdb_R_a__*delta_a_c__(:,1+nsource_a);
tmp_dbdc_delta_U_2a__(:,1+nsource_a) = P__ * tmp_dbdc_R_a__*delta_a_c__(:,1+nsource_a);
%tmp_dcda_delta_U_2a__(:,1+nsource_a) = P__ * tmp_dcda_R_a__*delta_a_c__(:,1+nsource_a);
%tmp_dcdb_delta_U_2a__(:,1+nsource_a) = P__ * tmp_dcdb_R_a__*delta_a_c__(:,1+nsource_a);
tmp_dcdc_delta_U_2a__(:,1+nsource_a) = P__ * tmp_dcdc_R_a__*delta_a_c__(:,1+nsource_a);
end;%for nsource_a=0:n_source_a-1;
%%%%;
% repeat for the perturbed volume. ;
%%%%;
if flag_dvol_ssnll;
tmp_euler_polar_a = +euler_polar_a_M_(1+nM);
tmp_euler_azimu_b = +euler_azimu_b_M_(1+nM);
tmp_euler_gamma_z = -euler_gamma_z_M_(1+nM);
tmp_R_a__ = Rz(-tmp_euler_gamma_z)*Ry(-tmp_euler_polar_a)*Rz(-tmp_euler_azimu_b);
tmp_da_R_a__ = -Rz(-tmp_euler_gamma_z)*dRy(-tmp_euler_polar_a)*Rz(-tmp_euler_azimu_b);
tmp_db_R_a__ = -Rz(-tmp_euler_gamma_z)*Ry(-tmp_euler_polar_a)*dRz(-tmp_euler_azimu_b);
tmp_dc_R_a__ = +dRz(-tmp_euler_gamma_z)*Ry(-tmp_euler_polar_a)*Rz(-tmp_euler_azimu_b);
%%;
tmp_delta_dvol_U_2a__ = zeros(2,n_source_dvol_a);
tmp_da_delta_dvol_U_2a__ = zeros(2,n_source_dvol_a);
tmp_db_delta_dvol_U_2a__ = zeros(2,n_source_dvol_a);
tmp_dc_delta_dvol_U_2a__ = zeros(2,n_source_dvol_a);
for nsource_dvol_a=0:n_source_dvol_a-1;
tmp_v_source_dvol_a = v_source_dvol_a_(1+nsource_dvol_a);
tmp_delta_dvol_U_2a__(:,1+nsource_dvol_a) = P__ * tmp_R_a__*delta_dvol_a_c__(:,1+nsource_dvol_a);
tmp_da_delta_dvol_U_2a__(:,1+nsource_dvol_a) = P__ * tmp_da_R_a__*delta_dvol_a_c__(:,1+nsource_dvol_a);
tmp_db_delta_dvol_U_2a__(:,1+nsource_dvol_a) = P__ * tmp_db_R_a__*delta_dvol_a_c__(:,1+nsource_dvol_a);
tmp_dc_delta_dvol_U_2a__(:,1+nsource_dvol_a) = P__ * tmp_dc_R_a__*delta_dvol_a_c__(:,1+nsource_dvol_a);
end;%for nsource_dvol_a=0:n_source_dvol_a-1;
end;%if flag_dvol_ssnll;
%%%%;
% the nM image M_k_p_wk_ = 2*k_p_r_wk_.*cos(k_p_w_wk_ - M_phi).*V_k_p_wk_ is assumed to come from volume b. ;
%%%%;
tmp_fromb_polar_a = +fromb_polar_a_M_(1+nM);
tmp_fromb_azimu_b = +fromb_azimu_b_M_(1+nM);
tmp_fromb_gamma_z = -fromb_gamma_z_M_(1+nM);
tmp_R_b__ = Rz(-tmp_fromb_gamma_z)*Ry(-tmp_fromb_polar_a)*Rz(-tmp_fromb_azimu_b);
tmp_delta_V_2b__ = zeros(2,n_source_b);
for nsource_b=0:n_source_b-1;
tmp_v_source_b = v_source_b_(1+nsource_b);
tmp_delta_V_2b__(:,1+nsource_b) = P__ * tmp_R_b__*delta_b_c__(:,1+nsource_b);
end;%for nsource_b=0:n_source_b-1;
%%%%%%%%;
if ~exist('tmp_delta_dvol_U_2a__','var'); tmp_delta_dvol_U_2a__=[]; end;
if ~exist('tmp_da_delta_dvol_U_2a__','var'); tmp_da_delta_dvol_U_2a__=[]; end;
if ~exist('tmp_db_delta_dvol_U_2a__','var'); tmp_db_delta_dvol_U_2a__=[]; end;
if ~exist('tmp_dc_delta_dvol_U_2a__','var'); tmp_dc_delta_dvol_U_2a__=[]; end;
%%%%;
%%%%%%%%;
if  flag_ssnll & ~flag_dtau_ssnll & ~flag_dtau_dtau_ssnll ;
[ ...
 ssnll_M ...
,dvol_ssnll_M ...
,dvol_dvol_ssnll_M ...
] = ...
I_xP__plus_xP__vs_xP__0( ...
 k_p_r_max ...
,CTF_phi ...
,n_source_a ...
,v_source_a_ ...
,tmp_delta_U_2a__ ...
,n_source_dvol_a ...
,v_source_dvol_a_ ...
,tmp_delta_dvol_U_2a__ ...
,M_phi ...
,n_source_b ...
,v_source_b_ ...
,tmp_delta_V_2b__ ...
);
end;%if  flag_ssnll & ~flag_dtau_ssnll & ~flag_dtau_dtau_ssnll ;
%%%%%%%%;
if  flag_ssnll &  flag_dtau_ssnll & ~flag_dtau_dtau_ssnll ;
[ ...
 ssnll_M ...
,dvol_ssnll_M ...
,dvol_dvol_ssnll_M ...
,da_ssnll_M ...
,db_ssnll_M ...
,dc_ssnll_M ...
,da_dvol_ssnll_M ...
,db_dvol_ssnll_M ...
,dc_dvol_ssnll_M ...
] = ...
I_xP__plus_xP__vs_xP__0( ...
 k_p_r_max ...
,CTF_phi ...
,n_source_a ...
,v_source_a_ ...
,tmp_delta_U_2a__ ...
,n_source_dvol_a ...
,v_source_dvol_a_ ...
,tmp_delta_dvol_U_2a__ ...
,M_phi ...
,n_source_b ...
,v_source_b_ ...
,tmp_delta_V_2b__ ...
,tmp_da_delta_U_2a__ ...
,tmp_db_delta_U_2a__ ...
,tmp_dc_delta_U_2a__ ...
);
end;%if  flag_ssnll &  flag_dtau_ssnll & ~flag_dtau_dtau_ssnll ;
%%%%%%%%;
if  flag_ssnll &  flag_dtau_ssnll &  flag_dtau_dtau_ssnll ;
[ ...
 ssnll_M ...
,dvol_ssnll_M ...
,dvol_dvol_ssnll_M ...
,da_ssnll_M ...
,db_ssnll_M ...
,dc_ssnll_M ...
,da_dvol_ssnll_M ...
,db_dvol_ssnll_M ...
,dc_dvol_ssnll_M ...
,dada_ssnll_M ...
,dadb_ssnll_M ...
,dadc_ssnll_M ...
,dbdb_ssnll_M ...
,dbdc_ssnll_M ...
,dcdc_ssnll_M ...
] = ...
I_xP__plus_xP__vs_xP__0( ...
 k_p_r_max ...
,CTF_phi ...
,n_source_a ...
,v_source_a_ ...
,tmp_delta_U_2a__ ...
,n_source_dvol_a ...
,v_source_dvol_a_ ...
,tmp_delta_dvol_U_2a__ ...
,M_phi ...
,n_source_b ...
,v_source_b_ ...
,tmp_delta_V_2b__ ...
,tmp_da_delta_U_2a__ ...
,tmp_db_delta_U_2a__ ...
,tmp_dc_delta_U_2a__ ...
,tmp_da_delta_dvol_U_2a__ ...
,tmp_db_delta_dvol_U_2a__ ...
,tmp_dc_delta_dvol_U_2a__ ...
,tmp_dada_delta_U_2a__ ...
,tmp_dadb_delta_U_2a__ ...
,tmp_dadc_delta_U_2a__ ...
,tmp_dbdb_delta_U_2a__ ...
,tmp_dbdc_delta_U_2a__ ...
,tmp_dcdc_delta_U_2a__ ...
);
end;%if  flag_ssnll &  flag_dtau_ssnll &  flag_dtau_dtau_ssnll ;
%%%%%%%%;
%%%%;
%%%%%%%%;
if flag_ssnll;
ssnll_M_(1+nM) = ssnll_M;
end;%if flag_ssnll;
if flag_dtau_ssnll;
da_ssnll_M_(1+nM) = da_ssnll_M;
db_ssnll_M_(1+nM) = db_ssnll_M;
dc_ssnll_M_(1+nM) = dc_ssnll_M;
dtau_ssnll_M3__(1+nM,1+0) = da_ssnll_M;
dtau_ssnll_M3__(1+nM,1+1) = db_ssnll_M;
dtau_ssnll_M3__(1+nM,1+2) = dc_ssnll_M;
end;%if flag_dtau_ssnll;
if flag_dtau_dtau_ssnll;
dada_ssnll_M_(1+nM) = dada_ssnll_M;
dadb_ssnll_M_(1+nM) = dadb_ssnll_M;
dadc_ssnll_M_(1+nM) = dadc_ssnll_M;
dbdb_ssnll_M_(1+nM) = dbdb_ssnll_M;
dbdc_ssnll_M_(1+nM) = dbdc_ssnll_M;
dcdc_ssnll_M_(1+nM) = dcdc_ssnll_M;
dtau_dtau_ssnll_M33___(1+nM,1+0,1+0) = dada_ssnll_M;
dtau_dtau_ssnll_M33___(1+nM,1+0,1+1) = dadb_ssnll_M;
dtau_dtau_ssnll_M33___(1+nM,1+0,1+2) = dadc_ssnll_M;
dtau_dtau_ssnll_M33___(1+nM,1+1,1+0) = dadb_ssnll_M;
dtau_dtau_ssnll_M33___(1+nM,1+1,1+1) = dbdb_ssnll_M;
dtau_dtau_ssnll_M33___(1+nM,1+1,1+2) = dbdc_ssnll_M;
dtau_dtau_ssnll_M33___(1+nM,1+2,1+0) = dadc_ssnll_M;
dtau_dtau_ssnll_M33___(1+nM,1+2,1+1) = dbdc_ssnll_M;
dtau_dtau_ssnll_M33___(1+nM,1+2,1+2) = dcdc_ssnll_M;
end;%if flag_dtau_dtau_ssnll;
if flag_dvol_ssnll;
dvol_ssnll_M_(1+nM) = dvol_ssnll_M;
end;%if flag_dvol_ssnll;
if flag_dtau_dvol_ssnll;
da_dvol_ssnll_M_(1+nM) = da_dvol_ssnll_M;
db_dvol_ssnll_M_(1+nM) = db_dvol_ssnll_M;
dc_dvol_ssnll_M_(1+nM) = dc_dvol_ssnll_M;
dtau_dvol_ssnll_M3__(1+nM,1+0) = da_dvol_ssnll_M;
dtau_dvol_ssnll_M3__(1+nM,1+1) = db_dvol_ssnll_M;
dtau_dvol_ssnll_M3__(1+nM,1+2) = dc_dvol_ssnll_M;
end;%if flag_dtau_dvol_ssnll;
if flag_dvol_dvol_ssnll;
dvol_dvol_ssnll_M_(1+nM) = dvol_dvol_ssnll_M;
end;%if flag_dvol_dvol_ssnll;
%%%%%%%%%%%%%%%%;
end;%for nM=0:n_M-1;
%%%%%%%%%%%%%%%%;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% ssnll: %0.6fs',tmp_t)); end;
%%%%%%%%%%%%%%%%;
ssnll = sum(ssnll_M_.*weight_imagecount_M_,[1]);
dvol_ssnll = sum(dvol_ssnll_M_.*weight_imagecount_M_,[1]);
dvol_dvol_ssnll = sum(dvol_dvol_ssnll_M_.*weight_imagecount_M_,[1]);
dtau_M3__ = [dtau_euler_polar_a_M_,dtau_euler_azimu_b_M_,dtau_euler_gamma_z_M_];
dtau_ssnll = sum(bsxfun(@times,bsxfun(@times,dtau_ssnll_M3__,dtau_M3__),weight_imagecount_M_),[1,2]);
dtau_dvol_ssnll = sum(bsxfun(@times,bsxfun(@times,dtau_dvol_ssnll_M3__,dtau_M3__),weight_imagecount_M_),[1,2]);
dtau_M33___ = bsxfun(@times,reshape(dtau_M3__,[n_M,3,1]),reshape(dtau_M3__,[n_M,1,3]));
dtau_dtau_ssnll = sum(bsxfun(@times,bsxfun(@times,dtau_dtau_ssnll_M33___,dtau_M33___),weight_imagecount_M_),[1,2,3]);
%%%%%%%%%%%%%%%%;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

