function S_q_ = interp_p_to_q(n_r,n_w_,n_A,S_p_);

str_thisfunction = 'interp_p_to_q';

if (nargin<4);
test_transforms_dr.m();
disp('returning'); return;
end;%if (nargin<4);

n_S = idivide(cast(numel(S_p_),'int64'),cast(n_A,'int64'));
if numel(S_p_)~=n_A*n_S;
disp(sprintf(' %% Warning, n_A %d n_S %d in %s',n_A,n_S,str_thisfunction));
end;%if numel(S_p_)~=n_A*n_S;

%%%%%%%%%%%%%%%%;
if (numel(unique(n_w_))> 1);
S_q_ = zeros(n_A*n_S,1);
%%%%%%%%;
if n_S==1;
%%%%;
% single transformation on adaptive grid. ;
%%%%;
n_w_max = n_w_(1+n_r-1);
ic=0;
for nr=0:n_r-1;
n_w = n_w_(1+nr);
if (n_w>0);
tmp_index_ = ic + [0:n_w-1];
S_q_(1+tmp_index_) = fft(S_p_(1+tmp_index_))/max(1,sqrt(n_w));
ic = ic + n_w;
end;%if (n_w>0);
end;%for nr=0:n_r-1;
assert(ic==n_A);
%%%%;
end;%if n_S==1;
%%%%%%%%;
if n_S> 1;
for nS=0:n_S-1;
tmp_index_ = nS*n_A + [0:n_A-1];
S_q_(1+tmp_index_) = interp_p_to_q(n_r,n_w_,n_A,S_p_(1+tmp_index_));
end;%for nS=0:n_S-1;
S_q_ = reshape(S_q_,size(S_p_));
end;%if n_S> 1;
%%%%%%%%;
end;%if (numel(unique(n_w_))> 1);
%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%;
if (numel(unique(n_w_))==1);
n_w = n_w_(1+0);
S_q_ = reshape(fft(reshape(S_p_,[n_w,n_r,n_S]),[],1+0)/max(1,sqrt(n_w)),[numel(S_p_),1]);
S_q_ = reshape(S_q_,size(S_p_));
end;%if (numel(unique(n_w_))==1);
%%%%%%%%%%%%%%%%;

