function S_p_ = interp_q_to_p_block_0(n_r,n_w_,n_A,S_q_);

if (nargin<4);
test_transforms_dr.m();
disp('returning'); return;
end;%if (nargin<4);

if (numel(unique(n_w_))> 1);
S_p_ = interp_p_to_p(n_r,n_w_,n_A,S_q_);
end;%if (numel(unique(n_w_))> 1);
if (numel(unique(n_w_))==1);
n_w = n_w_(1+0);
S_q_ = reshape(S_q_,n_w,[]);
S_p_ = ifft(S_q_,[],1)*sqrt(n_w);
S_p_ = reshape(S_p_,[numel(S_p_),1]);
end;%if (numel(unique(n_w_))==1);
