function ...
[ ...
 X0 ...
,C0 ...
,weight_Y_ ...
,e_k_Y_ ...
,u_k_Y_ ...
] = ...
register_spharm_to_spharm_3( ...
 verbose ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_ ...
,b_k_Y_ ...
,weight_Y_ ...
,e_k_Y_ ...
,u_k_Y_ ...
);
% registration between molecule_A and molecule_B ;
% ;
% Input: ;
% verbose = integer verbosity_level ;
% n_k_p_r = integer maximum k ;
% k_p_r_ = real array of length n_k_p_r; k_p_r_(nk) = k_p_r_value for shell nk ;
% weight_3d_k_p_r_ = real array of length n_k_p_r: weight_3d_k_p_r_(nk) = quadrature weight for shell nk. This should be designed so that 4*pi*sum(weight_3d_k_p_r_) = volume of sphere. ;
% l_max_ = integer array of length n_k_p_r; l_max_(nk) = spherical harmonic order on shell nk; l_max_(nk) corresponds to n_lm_(nk) = (l_max_(nk)+1)^2 coefficients ;
% a_k_Y_ = complex array of length \sum_{nk} (l_max_(nk)+1)^2 ; coefficients are ordered in a row, with m varying quickly and l varying slowly ;
% b_k_Y_ = complex array of length \sum_{nk} (l_max_(nk)+1)^2 ; coefficients are ordered in a row, with m varying quickly and l varying slowly ;
%
% Output: ;
% X0 = inner product, ;
% C0 = correlation, ;
% weight_Y_ = integration-weights. ;
% e_k_Y_ = indicator vector. ;
% u_k_Y_ = vector for projection onto average. ;

na=0;
if (nargin<1+na); verbose=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); k_p_r_=[]; end; na=na+1;
if (nargin<1+na); weight_3d_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); l_max_=[]; end; na=na+1;
if (nargin<1+na); a_k_Y_=[]; end; na=na+1;
if (nargin<1+na); b_k_Y_=[]; end; na=na+1;
if (nargin<1+na); weight_Y_=[]; end; na=na+1;
if (nargin<1+na); e_k_Y_=[]; end; na=na+1;
if (nargin<1+na); u_k_Y_=[]; end; na=na+1;

if (verbose); disp(sprintf(' %% [entering register_spharm_to_spharm_3]')); end;

if (verbose>1); disp(sprintf(' %% indices for counting arrays')); end;
n_lm_ = (l_max_+1).^2;
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
l_max_max = max(l_max_);
m_max_ = -l_max_max : +l_max_max;
n_m_max = length(m_max_);
if (verbose>1); disp(sprintf(' %% l_max_max %d n_m_max %d',l_max_max,n_m_max)); end;

if ( isempty(weight_Y_) | isempty(e_k_Y_) | isempty(u_k_Y_) );
%%%%%%%%;
if (verbose>1); disp(sprintf(' %% setting up weight_Y_ for inner-product and e_k_Y_ for integration')); end;
weight_Y_row_ = zeros(n_lm_sum,1);
weight_Y_col_ = zeros(n_lm_sum,1);
weight_Y_val_ = zeros(n_lm_sum,1);
e_k_Y_ = zeros(n_lm_sum,1);
na=0;
for nk_p_r=1:n_k_p_r;
tmp_ij_ = na + (1:n_lm_(nk_p_r));
weight_Y_row_(tmp_ij_) = tmp_ij_;
weight_Y_col_(tmp_ij_) = tmp_ij_;
weight_Y_val_(tmp_ij_) = weight_3d_k_p_r_(nk_p_r);
e_k_Y_(1+na) = 1;
na=na+n_lm_(nk_p_r);
end;%for nk_p_r=1:n_k_p_r;
weight_Y_ = sparse(weight_Y_row_,weight_Y_col_,weight_Y_val_,n_lm_sum,n_lm_sum);
e_avg = ctranspose(e_k_Y_)*(weight_Y_*e_k_Y_);
u_k_Y_ = e_k_Y_./max(sqrt(e_avg),1e-12);
%%%%%%%%;
end;%if ( isempty(weight_Y_) | isempty(e_k_Y_) | isempty(u_k_Y_) );

X0 = 0;
for nk_p_r = 0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
l_max = l_max_(1+nk_p_r); n_lm = n_lm_(1+nk_p_r);
tmp_ij_ = n_lm_csum_(1+nk_p_r) + (1:n_lm);
X0 = X0 + weight_3d_k_p_r_(1+nk_p_r) * sum(conj(a_k_Y_(tmp_ij_)).*(b_k_Y_(tmp_ij_)));
end;%for nk_p_r = 0:n_k_p_r-1;
if (verbose); disp(sprintf(' %% X0 %0.16f + i%0.16f',real(X0),imag(X0))); end;

if (verbose>1); disp(sprintf(' %% centering and normalizing both a_k_Y_ and b_k_Y_')); end;
a_avg = ctranspose(u_k_Y_)*(weight_Y_*a_k_Y_);
a_k_Y_norm_ = (a_k_Y_ - a_avg*u_k_Y_);
a_std = sqrt(ctranspose(a_k_Y_norm_)*(weight_Y_*a_k_Y_norm_));
a_k_Y_norm_ = a_k_Y_norm_./max(a_std,1e-12);
b_avg = ctranspose(u_k_Y_)*(weight_Y_*b_k_Y_);
b_k_Y_norm_ = (b_k_Y_ - b_avg*u_k_Y_);
b_std = sqrt(ctranspose(b_k_Y_norm_)*(weight_Y_*b_k_Y_norm_));
b_k_Y_norm_ = b_k_Y_norm_./max(b_std,1e-12);
a_avg = ctranspose(u_k_Y_)*(weight_Y_*a_k_Y_norm_);
a_std = sqrt(ctranspose(a_k_Y_norm_)*(weight_Y_*a_k_Y_norm_));
b_avg = ctranspose(u_k_Y_)*(weight_Y_*b_k_Y_norm_);
b_std = sqrt(ctranspose(b_k_Y_norm_)*(weight_Y_*b_k_Y_norm_));
if (verbose>1); disp(sprintf(' %% a_avg %+0.6f , a_std %+0.6f , b_avg %+0.6f , b_std %+0.6f',real(a_avg),a_std,real(b_avg),b_std)); end;

C0 = 0;
C0 = ctranspose(a_k_Y_norm_)*(weight_Y_*b_k_Y_norm_);
if (verbose); disp(sprintf(' %% C0 %0.16f + i%0.16f',real(C0),imag(C0))); end;

if (verbose); disp(sprintf(' %% [finished register_spharm_to_spharm_3]')); end;
