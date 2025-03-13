function ...
[ ...
 scaling_volumetric ...
,weight_3d_riesz_k_p_r_ ...
,weight_3d_riesz_qk_ ...
,weight_3d_riesz_qk__ ...
] = ...
get_weight_3d_riesz_2( ...
 flag_verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_2d_k_p_r_ ...
,n_w_max ...
,weight_2d_wk_ ...
,weight_3d_k_p_r_ ...
,n_q ...
,weight_3d_k_p_qk_ ...
);

na=0;
if (nargin<1+na); flag_verbose=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); k_p_r_=[]; end; na=na+1;
if (nargin<1+na); k_p_r_max=[]; end; na=na+1;
if (nargin<1+na); weight_2d_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); n_w_max=[]; end; na=na+1;
if (nargin<1+na); weight_2d_wk_=[]; end; na=na+1;
if (nargin<1+na); weight_3d_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); n_q=[]; end; na=na+1;
if (nargin<1+na); weight_3d_k_p_qk_=[]; end; na=na+1;

%%%%%%%%;
% Note that here we expect: ;
% sum(weight_2d_k_p_r_) == (pi*k_p_r_max^2), ;
% sum(weight_2d_wk_) == (pi*k_p_r_max^2)/(4*pi^2), ;
% sum(weight_3d_k_p_r_) == (4/3)*pi*k_p_r_max^3/(4*pi), ;
% sum(weight_3d_k_p_qk_) == (4/3)*pi*k_p_r_max^3, ;
% sum(weight_3d_riesz_k_p_r_) == 4*pi^2*k_p_r_max^2/(4*pi), ;
% sum(weight_3d_riesz_qk_) == 4*pi^2*k_p_r_max^2, ;
% So if principal-modes are used, then: ;
% k_p_r_max = 1.0;
% k_p_r_ = ones(pm_n_k_p_r_,1);
% sum(pm_weight_2d_k_p_r_) == (pi) --> pm_weight_2d_k_p_r_ = ones(pm_n_k_p_r,1)/max(1,pm_n_k_p_r)*pi ;
% sum(pm_weight_2d_wk_) == (pi)/(4*pi^2) = 1/(4*pi) --> pm_weight_2d_wk_ = ones(pm_n_w_sum,1)/max(1,pm_n_w_sum)/(4*pi) ;
% sum(pm_weight_3d_k_p_r_) == (4/3)*pi/(4*pi) = (1/3) --> pm_weight_3d_k_p_r_ = ones(pm_n_k_p_r,1)/max(1,pm_n_k_p_r)/3 ;
% sum(pm_weight_3d_k_p_qk_) == (4/3)*pi --> pm_weight_3d_k_p_qk_ = ones(n_q*pm_n_k_p_r,1)/max(1,n_q*pm_n_k_p_r)*((4/3)*pi) ;
% sum(pm_weight_3d_riesz_k_p_r_) == 4*pi^2/(4*pi) = pi --> pm_weight_3d_riesz_k_p_r_ = ones(pm_n_k_p_r,1)/max(1,pm_n_k_p_r)*pi ;
% sum(pm_weight_3d_riesz_qk_) == 4*pi^2 --> pm_weight_3d_riesz_qk_ = ones(n_q*pm_n_k_p_r,1)/max(1,n_q*pm_n_k_p_r)*(4*pi^2) ;
% With these conventions we see: ;
% pm_weight_3d_riesz_k_p_r_ = pm_weight_2d_k_p_r_ = ones(pm_n_k_p_r,1)/max(1,pm_n_k_p_r)*pi ;
% pm_weight_3d_riesz_qk_ = pm_weight_3d_k_p_qk_ * mean(pm_weight_2d_k_p_r_/pm_weight_3d_k_p_r_) = ones(n_q*pm_n_k_p_r,1)/max(1,n_q*pm_n_k_p_r)*((4/3)*pi)*3*pi ;
% i.e.: pm_weight_3d_riesz_qk_ = ones(n_q*pm_n_k_p_r,1)/max(1,n_q*pm_n_k_p_r)*(4*pi^2) ;
%%%%%%%%;

n_qk_csum_ = cumsum([0;n_q*ones(n_k_p_r,1)]);

%%%%%%%%;
% We construct the riesz integration-weights on the sphere. ;
% These are associated with the riesz-potential 1/k^2.5, ;
% or a weighting-function (for the squared-L2-norm) of 1/k. ;
%%%%%%%%;
weight_3d_riesz_k_p_r_ = weight_3d_k_p_r_;
weight_3d_riesz_qk_ = weight_3d_k_p_qk_;
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
weight_3d_k_p_r = weight_3d_k_p_r_(1+nk_p_r);
weight_2d_k_p_r = weight_2d_k_p_r_(1+nk_p_r);
weight_3d_riesz_k_p_r_(1+nk_p_r) = weight_3d_k_p_r_(1+nk_p_r) * weight_2d_k_p_r / max(1e-16,weight_3d_k_p_r);
tmp_index_ = n_qk_csum_(1+nk_p_r):n_qk_csum_(1+nk_p_r+1)-1;
if (flag_verbose> 1); disp(sprintf(' %% k_p_r %0.6f: sum(weight_3d_k_p_qk_(1+tmp_index_))/(weight_3d_k_p_r * 4*pi): %0.16f',k_p_r,sum(weight_3d_k_p_qk_(1+tmp_index_))/(4*pi*weight_3d_k_p_r))); end;
weight_3d_riesz_qk_(1+tmp_index_) = weight_3d_k_p_qk_(1+tmp_index_) * weight_2d_k_p_r / max(1e-16,weight_3d_k_p_r);
if (flag_verbose> 1); disp(sprintf(' %% k_p_r %0.6f: sum(weight_3d_riesz_qk_(1+tmp_index_))/(weight_2d_k_p_r * 4*pi): %0.16f',k_p_r,sum(weight_3d_riesz_qk_(1+tmp_index_))/(4*pi*weight_2d_k_p_r))); end;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
weight_3d_riesz_qk__ = reshape(weight_3d_riesz_qk_,[n_q,n_k_p_r]);

%%%%%%%%;
% Calibrate scaling factor. ;
%%%%%%%%;
term_deltafunc = sqrt(2*pi);
term_2 = (pi*k_p_r_max^2)/(4*pi^2);
if (flag_verbose>0); disp(sprintf(' %% sum(weight_2d_wk_) vs (pi*k_p_r_max^2)/(4*pi^2): %0.16f',fnorm(sum(weight_2d_wk_) - term_2))); end;
term_3 = (4/3)*pi*k_p_r_max^3;
if (flag_verbose>0); disp(sprintf(' %% sum(weight_3d_k_p_qk_) vs (4/3)*pi*k_p_r_max^3: %0.16f',fnorm(sum(weight_3d_k_p_qk_) - term_3))); end;
term_3r = (4*pi^2*k_p_r_max^2);
if (flag_verbose>0); disp(sprintf(' %% sum(weight_3d_riesz_qk_) vs 4*pi^2*k_p_r_max^2: %0.16f',fnorm(sum(weight_3d_riesz_qk_) - term_3r))); end;
scaling_volumetric = term_3r / term_2 / term_deltafunc ;
if (flag_verbose>0); disp(sprintf(' %% scaling_volumetric: %+0.6f',scaling_volumetric)); end;
if (flag_verbose>0); disp(sprintf(' %% (4*pi)^2 * sqrt(pi/2): %+0.6f',(4*pi)^2 * sqrt(pi/2))); end;
%%%%%%%%;

