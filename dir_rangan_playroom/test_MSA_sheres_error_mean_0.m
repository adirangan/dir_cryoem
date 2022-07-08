function C__ = test_MSA_sheres_error_mean_0(A_,dij,prctile_);
% calculates mean over percentile-range in prctile_. ;

verbose=0;

na=0;
if (nargin<1+na); A_ = []; end; na=na+1;
if (nargin<1+na); dij = []; end; na=na+1;
if (nargin<1+na); prctile_ = []; end; na=na+1;

if isempty(dij); dij=1; end;
if isempty(prctile_); prctile_ = [ 1;99]; end;

d_A_ = size(A_);
n_dim = numel(d_A_);
tmp_p_ = zeros(n_dim,1);
tmp_p_(n_dim) = dij;
tmp_p_(1:n_dim-1) = setdiff(1:n_dim,dij);

B_ = permute(A_,tmp_p_);
d_B_ = size(B_);
assert(fnorm(d_B_-d_A_(tmp_p_))==0);
n_pre = prod(d_B_(1:n_dim-1));
n_pos = d_B_(end);
B__ = reshape(B_,[n_pre,n_pos]);
C_ = zeros(n_pre,1);
for npre=0:n_pre-1;
tmp_B_ = B__(1+npre,:);
tmp_B_ = tmp_B_(find(isfinite(tmp_B_)));
tmp_B_upb = prctile(tmp_B_,max(prctile_));
tmp_B_lob = prctile(tmp_B_,min(prctile_));
tmp_ij_ = find( (tmp_B_>=tmp_B_lob) & (tmp_B_<=tmp_B_upb) );
if (verbose); disp(sprintf(' %% upb %f lob %f <-- %d',tmp_B_upb,tmp_B_lob,numel(tmp_ij_))); end;
tmp_B_mid_ = tmp_B_(tmp_ij_);
C_(1+npre) = mean(tmp_B_mid_);
end;%for npre=0:n_pre-1;
C__ = reshape(C_,d_B_(1:n_dim-1));



