function S_q_ = interp_p_to_q_block_0(n_r,n_w_,n_A,S_p_);

if (nargin<4);
test_transforms_dr.m();
disp('returning'); return;
end;%if (nargin<4);

if (numel(unique(n_w_))> 1);
S_q_ = interp_p_to_q(n_r,n_w_,n_A,S_p_);
end;%if (numel(unique(n_w_))> 1);
if (numel(unique(n_w_))==1);
n_w = n_w_(1+0);
S_p_ = reshape(S_p_,n_w,[]);
S_q_ = fft(S_p_,[],1)/sqrt(n_w);
S_q_ = reshape(S_q_,[numel(S_q_),1]);
end;%if (numel(unique(n_w_))==1);
