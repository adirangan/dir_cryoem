function S_q_ = interp_p_to_q_block_0(n_r,n_w_,n_A,S_p_);

if (nargin<4);
test_transforms_dr.m();
disp('returning'); return;
end;%if (nargin<4);

n_S = idivide(cast(numel(S_p_),'int64'),cast(n_A,'int64'));
if numel(S_p_)~=n_A*n_S;
disp(sprintf(' %% Warning, n_A %d n_S %d in interp_p_to_q_block_0',n_A,n_S));
end;%if numel(S_p_)~=n_A*n_S;

if (numel(unique(n_w_))> 1);
S_q_ = zeros(n_A*n_S,1);
for nS=0:n_S-1;
tmp_index_ = nS*n_A + [0:n_A-1];
S_q_(1+tmp_index_) = interp_p_to_q(n_r,n_w_,n_A,S_p_(:,1+tmp_index_));
end;%for nS=0:n_S-1;
end;%if (numel(unique(n_w_))> 1);
if (numel(unique(n_w_))==1);
n_w = n_w_(1+0);
S_q_ = reshape(fft(reshape(S_p_,[n_w,n_r,n_S]),[],1+0)/max(1,sqrt(n_w)),[numel(S_p_),1]);
end;%if (numel(unique(n_w_))==1);
