function ...
[ ...
 a_k_p_qk_ ...
,sqrt_2lp1_ ...
,sqrt_2mp1_ ...
,sqrt_rat0_m_ ...
,sqrt_rat3_lm__ ...
,sqrt_rat4_lm__ ...
] = ...
convert_spharm_to_k_p_uniform_over_n_k_p_r_5( ...
 flag_verbose ...
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
,l_max_ ...
,a_k_Y_yk_ ...
,sqrt_2lp1_ ...
,sqrt_2mp1_ ...
,sqrt_rat0_m_ ...
,sqrt_rat3_lm__ ...
,sqrt_rat4_lm__ ...
);
% uses spherical-harmonic-expansion a_k_Y_ to evaluate a_k_p_ on a collection of points on spherical shells determined by k_p_r_. ;
% We assume that the polar-representation and quadrature weights associated with these points have been previously calculated. ; 
% We also assume that flag_uniform_over_n_k_p_r==1 when generating the spherical grid. ;
% ;
% inputs: ;
% ;
% flag_verbose = integer verbosity_level. ;
% n_qk = integer total number of points. ;
% n_qk_csum_ = integer array of starting indices associated with each k-value. ;
% k_p_r_qk_ = real array of k-values for each point. ;
% k_p_azimu_b_qk_ = real array of azimu_b-values for each point. ;
% k_p_polar_a_qk_ = real array of polar_a-values for each point. ;
% weight_3d_k_p_qk_ = real array of quadrature weights for volume integral (for each point) (unused). ;
% weight_shell_qk_ = real array of quadrature weights for shell integral (for each point). ;
% n_k_p_r = integer maximum k. ;
% k_p_r_ = real array of length n_k_p_r; k_p_r_(nk_p_r) = k_p_r_value for shell nk_p_r. ;
% weight_3d_k_p_r_ = real array of length n_k_p_r; radial quadrature weights (already assumed to be a factor of weight_3d_k_p_qk_) (unused). ;
% l_max_ = integer array of length n_k_p_r; l_max_(nk_p_r) = spherical harmonic order on shell nk_p_r; l_max_(nk_p_r) corresponds to n_y_(nk_p_r) = (l_max_(nk_p_r)+1)^2 coefficients. ;
% a_k_Y_ = complex array of length \sum_{nk_p_r} (n_y_(nk_p_r)+1)^2 ; coefficients are ordered in a row, with m varying quickly, l varying slowly and k varying most slowly. ;
% ;
% outputs: ;
% ;
% a_k_p_qk_ = complex array of a-values for each point. ;

str_thisfunction = 'convert_spharm_to_k_p_uniform_over_n_k_p_r_5';

na=0;
if (nargin<1+na); flag_verbose=[]; end; na=na+1;
if (nargin<1+na); n_qk=[]; end; na=na+1;
if (nargin<1+na); n_qk_csum_=[]; end; na=na+1;
if (nargin<1+na); k_p_r_qk_=[]; end; na=na+1;
if (nargin<1+na); k_p_azimu_b_qk_=[]; end; na=na+1;
if (nargin<1+na); k_p_polar_a_qk_=[]; end; na=na+1;
if (nargin<1+na); weight_3d_k_p_qk_=[]; end; na=na+1;
if (nargin<1+na); weight_shell_qk_=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); k_p_r_=[]; end; na=na+1;
if (nargin<1+na); weight_3d_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); l_max_=[]; end; na=na+1;
if (nargin<1+na); a_k_Y_yk_=[]; end; na=na+1;
if (nargin<1+na); sqrt_2lp1_=[]; end; na=na+1;
if (nargin<1+na); sqrt_2mp1_=[]; end; na=na+1;
if (nargin<1+na); sqrt_rat0_m_=[]; end; na=na+1;
if (nargin<1+na); sqrt_rat3_lm__=[]; end; na=na+1;
if (nargin<1+na); sqrt_rat4_lm__=[]; end; na=na+1;

n_y_ = (1+l_max_(:)).^2; n_y_sum = sum(n_y_); n_y_max = max(n_y_); n_y_csum_ = cumsum([0;n_y_]);
if (flag_verbose>0); disp(sprintf(' %% [entering %s] n_qk %d, n_y_sum %d',str_thisfunction,n_qk,sum(n_y_))); end;

if (std(diff(n_qk_csum_),1)>1e-6); disp(sprintf(' %% Warning, flag_uniform_over_n_k_p_r==1 expected in %s',str_thisfunction)); end;
n_q = n_qk_csum_(1+1);
assert(n_qk==n_q*n_k_p_r);
k_p_r_qk__ = reshape(k_p_r_qk_,[n_q,n_k_p_r]);
if (max(std(k_p_r_qk__,1,1))>1e-6); disp(sprintf(' %% Warning, flag_uniform_over_n_k_p_r==1 expected in %s',str_thisfunction)); end;
k_p_azimu_b_qk__ = reshape(k_p_azimu_b_qk_,[n_q,n_k_p_r]);
if (max(std(k_p_azimu_b_qk__,1,2))>1e-6); disp(sprintf(' %% Warning, flag_uniform_over_n_k_p_r==1 expected in %s',str_thisfunction)); end;
k_p_polar_a_qk__ = reshape(k_p_polar_a_qk_,[n_q,n_k_p_r]);
if (max(std(k_p_polar_a_qk__,1,2))>1e-6); disp(sprintf(' %% Warning, flag_uniform_over_n_k_p_r==1 expected in %s',str_thisfunction)); end;

k_p_azimu_b_q_ = k_p_azimu_b_qk__(:,1+0);
k_p_polar_a_q_ = k_p_polar_a_qk__(:,1+0);
[k_p_polar_a_unique_,ij_unique_,ij_return_] = unique(k_p_polar_a_q_); index_return_ = ij_return_ - 1;
n_polar_a_unique = numel(k_p_polar_a_unique_);
l_max_max = max(l_max_);
if ~exist('sqrt_2lp1_','var'); sqrt_2lp1_=[]; end;
if ~exist('sqrt_2mp1_','var'); sqrt_2mp1_=[]; end;
if ~exist('sqrt_rat0_m_','var'); sqrt_rat0_m_=[]; end;
if ~exist('sqrt_rat3_lm__','var'); sqrt_rat3_lm__=[]; end;
if ~exist('sqrt_rat4_lm__','var'); sqrt_rat4_lm__=[]; end;
tmp_t = tic();
parameter_ylgndr = struct('type','ylgndr');
parameter_ylgndr.flag_verbose = 0;
parameter_ylgndr.flag_d = 0;
parameter_ylgndr.flag_dd = 0;
[ ...
 parameter_ylgndr ...
,d0y_jlm___ ...
,~ ...
,~ ...
,sqrt_2lp1_ ...
,sqrt_2mp1_ ...
,sqrt_rat0_m_ ...
,sqrt_rat3_lm__ ...
,sqrt_rat4_lm__ ...
] = ...
ylgndr_2( ...
 parameter_ylgndr ...
,l_max_max ...
,cos(k_p_polar_a_unique_) ...
,sqrt_2lp1_ ...
,sqrt_2mp1_ ...
,sqrt_rat0_m_ ...
,sqrt_rat3_lm__ ...
,sqrt_rat4_lm__ ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% ylgndr_2: %0.6fs',tmp_t)); end;

%%%%%%%%;
% This version involves preliminary inflation of m_val_, ;
% as well as extraneous exp evaluations, ;
% all to avoid memory movement. ;
%%%%%%%%;

tmp_t = tic();
d0y_lmj___ = permute(d0y_jlm___,1+[1,2,0]); %<-- permutation before inflation is faster. ;
d0y_jml___ = permute(d0y_lmj___(:,1+abs(-l_max_max:+l_max_max),:),1+[2,1,0]); %<-- retain unique cos(polar_a_{j}) for now.; 
d0y_jy__ = zeros(n_polar_a_unique,n_y_max);
for nl=0:numel(l_max_)-1;
l_max = l_max_(1+nl);
for l_val=0:l_max;
m_val_ = -l_val:+l_val;
d0y_jy__(:,1+l_val.^2+l_val+m_val_) = d0y_jml___(:,1+l_max_max+m_val_,1+l_val);
end;%for l_val=0:l_max;
end;%for nl=0:numel(l_max_)-1;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% d0y_jy__: %0.6fs',tmp_t)); end;

tmp_t = tic();
d0y_qy__ = d0y_jy__(1+index_return_,:);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% d0y_jy__: %0.6fs',tmp_t)); end;

tmp_t = tic();
m_val_y_ = zeros(n_y_max,1);
for nl=0:numel(l_max_)-1;
l_max = l_max_(1+nl);
for l_val=0:l_max;
m_val_ = -l_val:+l_val;
m_val_y_(1+l_val.^2+l_val+m_val_) = m_val_;
end;%for l_val=0:l_max;
end;%for nl=0:numel(l_max_)-1;
expimb_qy__ = permute(reshape(exp(+i*bsxfun(@times,m_val_y_,reshape(k_p_azimu_b_q_,[1,n_q]))),[n_y_max,n_q]),1+[1,0]);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% expimb_qy__: %0.6fs',tmp_t)); end;

tmp_t = tic();
d0Y_qy__ = d0y_qy__.*expimb_qy__ / sqrt(4*pi);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% d0Y_qy__: %0.6fs',tmp_t)); end;

tmp_t = tic();
a_k_Y_yk__ = local_yk__from_yk_(n_k_p_r,l_max_,a_k_Y_yk_);
a_k_p_qk_ = reshape(d0Y_qy__*a_k_Y_yk__,[n_qk,1]);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% a_k_p_qk_: %0.6fs',tmp_t)); end;

%{
%%%%%%%%;
% Here is a slower version. ;
%%%%%%%%;
tmp_t = tic();
d0y_lmj___ = permute(d0y_jlm___,1+[1,2,0]); %<-- permutation before inflation is faster. ;
d0y_lmq___ = d0y_lmj___(:,1+abs(-l_max_max:+l_max_max),1+index_return_);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% d0y_lmq___: %0.6fs',tmp_t)); end;

tmp_t = tic();
expimb_mq__ = exp(+i*bsxfun(@times,reshape(-l_max_max:+l_max_max,[1+2*l_max_max,1]),reshape(k_p_azimu_b_q_,[1,n_q])));
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% expimb_mq__: %0.6fs',tmp_t)); end;

tmp_t = tic();
d0Y_qml___ = permute(bsxfun(@times,d0y_lmq___,reshape(expimb_mq__,[1,1+2*l_max_max,n_q]))/sqrt(4*pi),1+[2,1,0]);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% d0Y_qml___: %0.6fs',tmp_t)); end;

tmp_t = tic();
n_y_ = (1+l_max_(:)).^2; n_y_sum = sum(n_y_); n_y_max = max(n_y_); n_y_csum_ = cumsum([0;n_y_]);
d0Y_qy__ = zeros(n_q,n_y_max);
for nl=0:numel(l_max_)-1;
l_max = l_max_(1+nl);
for l_val=0:l_max;
m_val_ = -l_val:+l_val;
d0Y_qy__(:,1+l_val.^2+l_val+m_val_) = d0Y_qml___(:,1+l_max_max+m_val_,1+l_val);
end;%for l_val=0:l_max;
end;%for nl=0:numel(l_max_)-1;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% d0Y_qy__: %0.6fs',tmp_t)); end;

tmp_t = tic();
a_k_Y_yk__ = local_yk__from_yk_(n_k_p_r,l_max_,a_k_Y_yk_);
a_k_p_qk_ = reshape(d0Y_qy__*a_k_Y_yk__,[n_qk,1]);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% a_k_p_qk_: %0.6fs',tmp_t)); end;
%}

%{
%%%%%%%%;
% Here is an even slower version. ;
%%%%%%%%;
tmp_t = tic();
d0y_lmj___ = permute(d0y_jlm___,1+[1,2,0]); %<-- permutation before inflation is faster. ;
d0y_lmq___ = d0y_lmj___(:,1+abs(-l_max_max:+l_max_max),1+index_return_);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% d0y_lmq___: %0.6fs',tmp_t)); end;

tmp_t = tic();
expimb_mq__ = exp(+i*bsxfun(@times,reshape(-l_max_max:+l_max_max,[1+2*l_max_max,1]),reshape(k_p_azimu_b_q_,[1,n_q])));
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% expimb_mq__: %0.6fs',tmp_t)); end;

tmp_t = tic();
d0Y_lmq___ = bsxfun(@times,d0y_lmq___,reshape(expimb_mq__,[1,1+2*l_max_max,n_q]))/sqrt(4*pi);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% d0Y_lmq___: %0.6fs',tmp_t)); end;

tmp_t = tic();
n_y_ = (1+l_max_(:)).^2; n_y_sum = sum(n_y_); n_y_max = max(n_y_); n_y_csum_ = cumsum([0;n_y_]);
d0Y_yq__ = zeros(n_y_max,n_q);
for nl=0:numel(l_max_)-1;
l_max = l_max_(1+nl);
for l_val=0:l_max;
m_val_ = -l_val:+l_val;
d0Y_yq__(1+l_val.^2+l_val+m_val_,:) = d0Y_lmq___(1+l_val,1+l_max_max+m_val_,:);
end;%for l_val=0:l_max;
end;%for nl=0:numel(l_max_)-1;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% d0Y_yq__: %0.6fs',tmp_t)); end;

tmp_t = tic();
a_k_Y_yk__ = local_yk__from_yk_(n_k_p_r,l_max_,a_k_Y_yk_);
a_k_p_qk_ = reshape(transpose(d0Y_yq__)*a_k_Y_yk__,[n_qk,1]);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% a_k_p_qk_: %0.6fs',tmp_t)); end;
%}

if (flag_verbose>0); disp(sprintf(' %% [finished %s] n_qk %d, n_y_sum %d',str_thisfunction,n_qk,sum(n_y_))); end;
