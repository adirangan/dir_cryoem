function ...
[ ...
 parameter ...
,dSd0_k_p_wS_ ...
,dSd1_k_p_wS_ ...
,dSd2_k_p_wS_ ...
,ddSd00_k_p_wS_ ...
,ddSd01_k_p_wS_ ...
,ddSd02_k_p_wS_ ...
,ddSd11_k_p_wS_ ...
,ddSd12_k_p_wS_ ...
,ddSd22_k_p_wS_ ...
] = ...
cartesian_from_shell_10( ...
 parameter ...
,n_w ...
,n_S ...
,k_p_polar_a_wS_ ...
,k_p_azimu_b_wS_ ...
,S_k_p_wS_ ...
,dSda_k_p_wS_ ...
,dSdb_k_p_wS_ ...
,ddSdaa_k_p_wS_ ...
,ddSdab_k_p_wS_ ...
,ddSdbb_k_p_wS_ ...
);
%%%%%%%%;
% Here we assume that: ;
% k_p_polar_a_wS_ ...
% k_p_azimu_b_wS_ ...
% have been passed into ;
% shell_k_p_scatter_from_adaptive_interpolate_n_9, ;
% with outputs: ;
% S_k_p_wS_ = wS_from_single_shell_sba__*a_k_p_single_shell_;
% dSda_k_p_wS_ = dwSda_from_single_shell_sba__*a_k_p_single_shell_;
% dSdb_k_p_wS_ = dwSdb_from_single_shell_sba__*a_k_p_single_shell_;
% ddSdaa_k_p_wS_ = ddwSdaa_from_single_shell_sba__*a_k_p_single_shell_;
% ddSdab_k_p_wS_ = ddwSdab_from_single_shell_sba__*a_k_p_single_shell_;
% ddSdbb_k_p_wS_ = ddwSdbb_from_single_shell_sba__*a_k_p_single_shell_;
% from which we compute the corresponding: ;
% dSd0_k_p_wS_ ;
% dSd1_k_p_wS_ ;
% dSd2_k_p_wS_ ;
% ddSd00_k_p_wS_ ;
% ddSd01_k_p_wS_ ;
% ddSd02_k_p_wS_ ;
% ddSd11_k_p_wS_ ;
% ddSd12_k_p_wS_ ;
% ddSd22_k_p_wS_ ;
%%%%%%%%;

str_thisfunction = 'cartesian_from_shell_10';

if nargin<1;
disp(sprintf(' %% testing %s',str_thisfunction));
disp(sprintf(' %% see test_slice_vs_volume_integral_from_a_k_p_6.m'));
disp('returning'); return;
end;%if nargin<1;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_w=[]; end; na=na+1;
if (nargin<1+na); n_S=[]; end; na=na+1;
if (nargin<1+na); k_p_polar_a_wS_=[]; end; na=na+1;
if (nargin<1+na); k_p_azimu_b_wS_=[]; end; na=na+1;
if (nargin<1+na); S_k_p_wS_=[]; end; na=na+1;
if (nargin<1+na); dSda_k_p_wS_=[]; end; na=na+1;
if (nargin<1+na); dSdb_k_p_wS_=[]; end; na=na+1;
if (nargin<1+na); ddSdaa_k_p_wS_=[]; end; na=na+1;
if (nargin<1+na); ddSdab_k_p_wS_=[]; end; na=na+1;
if (nargin<1+na); ddSdbb_k_p_wS_=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
flag_verbose=parameter.flag_verbose;
if ~isfield(parameter,'tolerance_pinv'); parameter.tolerance_pinv=1e-6; end;
tolerance_pinv=parameter.tolerance_pinv;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

flag_d = ~isempty(dSda_k_p_wS_) & ~isempty(dSdb_k_p_wS_) ;
flag_dd = ~isempty(ddSdaa_k_p_wS_) & ~isempty(ddSdab_k_p_wS_) & ~isempty(ddSdbb_k_p_wS_) ;

n_wS = n_w*n_S;
%%%%%%%%;
k_c_0_wS_ = cos(k_p_azimu_b_wS_).*sin(k_p_polar_a_wS_);
k_c_1_wS_ = sin(k_p_azimu_b_wS_).*sin(k_p_polar_a_wS_);
k_c_2_wS_ = cos(k_p_polar_a_wS_);
%%%%%%%%;
if flag_d;
da_k_c_0_wS_ = +cos(k_p_azimu_b_wS_).*cos(k_p_polar_a_wS_);
db_k_c_0_wS_ = -sin(k_p_azimu_b_wS_).*sin(k_p_polar_a_wS_);
da_k_c_1_wS_ = +sin(k_p_azimu_b_wS_).*cos(k_p_polar_a_wS_);
db_k_c_1_wS_ = +cos(k_p_azimu_b_wS_).*sin(k_p_polar_a_wS_);
da_k_c_2_wS_ = -sin(k_p_polar_a_wS_);
db_k_c_2_wS_ = 0*k_p_polar_a_wS_;
end;%if flag_d;
%%%%%%%%;
if flag_dd;
daa_k_c_0_wS_ = -cos(k_p_azimu_b_wS_).*sin(k_p_polar_a_wS_);
dab_k_c_0_wS_ = -sin(k_p_azimu_b_wS_).*cos(k_p_polar_a_wS_);
dba_k_c_0_wS_ = -sin(k_p_azimu_b_wS_).*cos(k_p_polar_a_wS_);
dbb_k_c_0_wS_ = -cos(k_p_azimu_b_wS_).*sin(k_p_polar_a_wS_);
daa_k_c_1_wS_ = -sin(k_p_azimu_b_wS_).*sin(k_p_polar_a_wS_);
dab_k_c_1_wS_ = +cos(k_p_azimu_b_wS_).*cos(k_p_polar_a_wS_);
dba_k_c_1_wS_ = +cos(k_p_azimu_b_wS_).*cos(k_p_polar_a_wS_);
dbb_k_c_1_wS_ = -sin(k_p_azimu_b_wS_).*sin(k_p_polar_a_wS_);
daa_k_c_2_wS_ = -cos(k_p_polar_a_wS_);
dab_k_c_2_wS_ = 0*k_p_polar_a_wS_;
dba_k_c_2_wS_ = 0*k_p_polar_a_wS_;
dbb_k_c_2_wS_ = 0*k_p_polar_a_wS_;
end;%if flag_dd;
%%%%%%%%;

dSd0_k_p_wS_=[];
dSd1_k_p_wS_=[];
dSd2_k_p_wS_=[];
%%%%%%%%;
if flag_d;
dSd0_k_p_wS_ = zeros(n_wS,1);
dSd1_k_p_wS_ = zeros(n_wS,1);
dSd2_k_p_wS_ = zeros(n_wS,1);
for nwS=0:n_wS-1;
%%%%;
dSda = dSda_k_p_wS_(1+nwS);
dSdb = dSdb_k_p_wS_(1+nwS);
%%%%;
d0da = da_k_c_0_wS_(1+nwS);
d0db = db_k_c_0_wS_(1+nwS);
d1da = da_k_c_1_wS_(1+nwS);
d1db = db_k_c_1_wS_(1+nwS);
d2da = da_k_c_2_wS_(1+nwS);
d2db = db_k_c_2_wS_(1+nwS);
%%%%;
RHS_ = ...
[ ...
  dSda ...
; dSdb ...
];
J__ = ...
[ ...
  d0da , d1da , d2da ...
; d0db , d1db , d2db ...
];
LHS_ = pinv(J__,tolerance_pinv)*RHS_;
dSd0 = LHS_(1+0);
dSd1 = LHS_(1+1);
dSd2 = LHS_(1+2);
dSd0_k_p_wS_(1+nwS) = dSd0;
dSd1_k_p_wS_(1+nwS) = dSd1;
dSd2_k_p_wS_(1+nwS) = dSd2;
end;%for nwS=0:n_wS-1;
end;%if flag_d;
%%%%%%%%;
%{
%%%%%%%%;
if flag_d;
J_wS32___ = ...
cat(3 ...
    ,cat(2 ...
	 ,da_k_c_0_wS_(:) ...
	 ,db_k_c_0_wS_(:) ...
	) ...
    ,cat(2 ...
	 ,da_k_c_1_wS_(:) ...
	 ,db_k_c_1_wS_(:) ...
	) ...
    ,cat(2 ...
	 ,da_k_c_2_wS_(:) ...
	 ,db_k_c_2_wS_(:) ...
	) ...
   ) ...
;
J_32wS___ = permute(J_wS32___,[2,3,1]);
dSd0_k_p_wS_ = zeros(n_wS,1);
dSd1_k_p_wS_ = zeros(n_wS,1);
dSd2_k_p_wS_ = zeros(n_wS,1);
for nwS=0:n_wS-1;
J_32__ = J_32wS___(:,:,1+nwS);
inv_J_23__ = pinv(J_32__,tolerance_pinv);
tmp_ = inv_J_23__*[dSda_k_p_wS_(1+nwS) ; dSdb_k_p_wS_(1+nwS)];
dSd0_k_p_wS_(1+nwS) = tmp_(1+0);
dSd1_k_p_wS_(1+nwS) = tmp_(1+1);
dSd2_k_p_wS_(1+nwS) = tmp_(1+2);
end;%for nwS=0:n_wS-1;
end;%if flag_d;
%%%%%%%%;
%}

ddSd00_k_p_wS_=[];
ddSd01_k_p_wS_=[];
ddSd02_k_p_wS_=[];
ddSd11_k_p_wS_=[];
ddSd12_k_p_wS_=[];
ddSd22_k_p_wS_=[];
%%%%%%%%;
if flag_dd;
ddSd00_k_p_wS_ = zeros(n_wS,1);
ddSd01_k_p_wS_ = zeros(n_wS,1);
ddSd02_k_p_wS_ = zeros(n_wS,1);
ddSd11_k_p_wS_ = zeros(n_wS,1);
ddSd12_k_p_wS_ = zeros(n_wS,1);
ddSd22_k_p_wS_ = zeros(n_wS,1);
for nwS=0:n_wS-1;
%%%%;
dSda = dSda_k_p_wS_(1+nwS);
dSdb = dSdb_k_p_wS_(1+nwS);
ddSdaa = ddSdaa_k_p_wS_(1+nwS);
ddSdab = ddSdab_k_p_wS_(1+nwS);
ddSdba = ddSdab;
ddSdbb = ddSdbb_k_p_wS_(1+nwS);
%%%%;
d0da = da_k_c_0_wS_(1+nwS);
d0db = db_k_c_0_wS_(1+nwS);
d1da = da_k_c_1_wS_(1+nwS);
d1db = db_k_c_1_wS_(1+nwS);
d2da = da_k_c_2_wS_(1+nwS);
d2db = db_k_c_2_wS_(1+nwS);
dd0daa = daa_k_c_0_wS_(1+nwS);
dd0dab = dab_k_c_0_wS_(1+nwS);
dd0dba = dba_k_c_0_wS_(1+nwS);
dd0dbb = dbb_k_c_0_wS_(1+nwS);
dd1daa = daa_k_c_1_wS_(1+nwS);
dd1dab = dab_k_c_1_wS_(1+nwS);
dd1dba = dba_k_c_1_wS_(1+nwS);
dd1dbb = dbb_k_c_1_wS_(1+nwS);
dd2daa = daa_k_c_2_wS_(1+nwS);
dd2dab = dab_k_c_2_wS_(1+nwS);
dd2dba = dba_k_c_2_wS_(1+nwS);
dd2dbb = dbb_k_c_2_wS_(1+nwS);
%%%%;
dSd0 = dSd0_k_p_wS_(1+nwS);
dSd1 = dSd1_k_p_wS_(1+nwS);
dSd2 = dSd2_k_p_wS_(1+nwS);
RHS_ = ...
[ ...
  ddSdaa - dSd0*dd0daa - dSd1*dd1daa - dSd2*dd2daa ...
; ddSdab - dSd0*dd0dab - dSd1*dd1dab - dSd2*dd2dab ...
; ddSdba - dSd0*dd0dba - dSd1*dd1dba - dSd2*dd2dba ...
; ddSdbb - dSd0*dd0dbb - dSd1*dd1dbb - dSd2*dd2dbb ...
];
H__ = ...
[ ...
  d0da*d0da , d0da*d1da , d0da*d2da , d1da*d0da , d1da*d1da , d1da*d2da , d2da*d0da , d2da*d1da , d2da*d2da ...
; d0da*d0db , d0da*d1db , d0da*d2db , d1da*d0db , d1da*d1db , d1da*d2db , d2da*d0db , d2da*d1db , d2da*d2db ...
; d0db*d0da , d0db*d1da , d0db*d2da , d1db*d0da , d1db*d1da , d1db*d2da , d2db*d0da , d2db*d1da , d2db*d2da ...
; d0db*d0db , d0db*d1db , d0db*d2db , d1db*d0db , d1db*d1db , d1db*d2db , d2db*d0db , d2db*d1db , d2db*d2db ...
];
LHS_ = pinv(H__,tolerance_pinv)*RHS_;
LHS__ = reshape(LHS_,3,3);
LHS__ = 0.5*(LHS__ + transpose(LHS__));
ddSd00 = LHS__(1+0,1+0);
ddSd01 = LHS__(1+0,1+1);
ddSd02 = LHS__(1+0,1+2);
ddSd11 = LHS__(1+1,1+1);
ddSd12 = LHS__(1+1,1+2);
ddSd22 = LHS__(1+2,1+2);
ddSd00_k_p_wS_(1+nwS) = ddSd00;
ddSd01_k_p_wS_(1+nwS) = ddSd01;
ddSd02_k_p_wS_(1+nwS) = ddSd02;
ddSd11_k_p_wS_(1+nwS) = ddSd11;
ddSd12_k_p_wS_(1+nwS) = ddSd12;
ddSd22_k_p_wS_(1+nwS) = ddSd22;
end;%for nwS=0:n_wS-1;
end;%if flag_dd;
%%%%%%%%;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;



