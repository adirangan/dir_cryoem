function ...
[ ...
 parameter ...
,C_zero_qk__ ...
,C_zero_wk__ ...
,gamma_new_M_ ...
,Rp0_M_ ...
,tmp_F_inv_qM__ ...
,tmp_F_Mq__ ...
] = ...
test_msa_update_3( ...
 parameter ...	   
,A00_qk__ ...
,M_kM__ ...
);

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); A00_qk__=[]; end; na=na+1;
if (nargin<1+na); M_kM__=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;
if ~isfield(parameter,'n_newt'); parameter.n_newt = 2; end;
n_newt = parameter.n_newt;

n_q = size(A00_qk__,1);
q_max = (n_q-1)/2;
q_ = transpose(-q_max:+q_max);
assert(n_q==numel(q_));
n_k_p_r = size(A00_qk__,2);
assert(size(M_kM__,1)==n_k_p_r);
n_M = size(M_kM__,2);
n_gamma = max(n_q*32,1024);
gamma_ = linspace(0,2*pi,1+n_gamma);
gamma_ = reshape(gamma_(1:n_gamma),[n_gamma,1]);
dgamma = mean(diff(gamma_));
n_w = n_gamma;
F_wq__ = exp(+i*gamma_*transpose(q_));
F_inv_qw__ = ctranspose(F_wq__)/n_gamma;

L2_X00_qk__ = @(X00_qk__) sqrt(2*pi*sum(abs(X00_qk__).^2,'all'));
L2_X00_wk__ = @(X00_wk__) sqrt(dgamma*sum(abs(X00_wk__).^2,'all'));

Ap0_qk__ = bsxfun(@times,A00_qk__,i*reshape(q_,n_q,1));
App_qk__ = bsxfun(@times,Ap0_qk__,i*reshape(q_,n_q,1));

A00_wk__ = F_wq__*A00_qk__;
Ap0_wk__ = F_wq__*Ap0_qk__;
App_wk__ = F_wq__*App_qk__;

f_A00_w = @(gamma) sum(bsxfun(@times,A00_qk__,exp(+i*gamma*q_)),1);
f_Ap0_w = @(gamma) sum(bsxfun(@times,Ap0_qk__,exp(+i*gamma*q_)),1);
f_App_w = @(gamma) sum(bsxfun(@times,App_qk__,exp(+i*gamma*q_)),1);

f_A00_w_ = @(gamma_) exp(+i*gamma_(:)*transpose(q_))*A00_qk__;
f_Ap0_w_ = @(gamma_) exp(+i*gamma_(:)*transpose(q_))*Ap0_qk__;
f_App_w_ = @(gamma_) exp(+i*gamma_(:)*transpose(q_))*App_qk__;

f_R00_w_ = @(gamma_,M_k_) sum(abs(bsxfun(@minus,f_A00_w_(gamma_),M_k_)).^2,2);
f_Rp0_w_ = @(gamma_,M_k_) 2*real(sum(conj(bsxfun(@minus,f_A00_w_(gamma_),M_k_)).*f_Ap0_w_(gamma_),2));
f_Rpp_w_ = @(gamma_,M_k_) 2*real(sum(conj(f_Ap0_w_(gamma_)).*f_Ap0_w_(gamma_),2)) + 2*real(sum(conj(bsxfun(@minus,f_A00_w_(gamma_),M_k_)).*f_App_w_(gamma_),2));

Rp0_M_ = zeros(n_M,1);
gamma_new_M_ = zeros(n_M,1);
for nM=0:n_M-1;
M_k_ = reshape(M_kM__(:,1+nM),[1,n_k_p_r]);
tmp_R_w_ = sum(abs(bsxfun(@minus,A00_wk__,M_k_)).^2,2);
[~,ij_gamma_upd] = min(tmp_R_w_); index_gamma_upd = ij_gamma_upd - 1;
gamma_pos = gamma_(1+index_gamma_upd);
gamma_upd = gamma_pos + dgamma*quadratic_1d_interpolation_0(-tmp_R_w_(1+periodize(index_gamma_upd+[-1:+1],0,n_gamma)));
gamma_new = gamma_upd;
for nnewt=0:n_newt-1;
tmp_Rpp = f_Rpp_w_(gamma_new,M_k_); if abs(tmp_Rpp)<1e-16; tmp_Rpp = 1e-16*sign(tmp_Rpp); end;
gamma_new = gamma_new - f_Rp0_w_(gamma_new,M_k_)/tmp_Rpp;
end;%for nnewt=0:n_newt-1;
Rp0_M_(1+nM) = f_Rp0_w_(gamma_new,M_k_);
gamma_new_M_(1+nM) = gamma_new;
end;%for nM=0:n_M-1;
%%%%;
tmp_F_Mq__ = exp(+i*gamma_new_M_*transpose(q_));
tmp_F_inv_qM__ = pinv(tmp_F_Mq__);
C_zero_qk__ = tmp_F_inv_qM__*transpose(M_kM__);
C_zero_wk__ = F_wq__*C_zero_qk__;
