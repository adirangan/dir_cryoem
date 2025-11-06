function S_p_ = interp_q_to_p(n_r,n_w_,n_A,S_q_);

str_thisfunction = 'interp_q_to_p';

if (nargin<4);
test_transforms_dr.m();
disp('returning'); return;
end;%if (nargin<4);

n_S = idivide(cast(numel(S_q_),'int64'),cast(n_A,'int64'));
if numel(S_q_)~=n_A*n_S;
disp(sprintf(' %% Warning, n_A %d n_S %d in %s',n_A,n_S,str_thisfunction));
end;%if numel(S_q_)~=n_A*n_S;

%%%%%%%%%%%%%%%%;
if (numel(unique(n_w_))> 1);
S_p_ = zeros(n_A*n_S,1);
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
S_p_(1+tmp_index_) = ifft(S_q_(1+tmp_index_))*max(0,sqrt(n_w));
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
S_p_(1+tmp_index_) = interp_q_to_p(n_r,n_w_,n_A,S_q_(1+tmp_index_));
end;%for nS=0:n_S-1;
end;%if n_S> 1;
%%%%%%%%;
end;%if (numel(unique(n_w_))> 1);
%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%;
if (numel(unique(n_w_))==1);
n_w = n_w_(1+0);
S_p_ = reshape(ifft(reshape(S_q_,[n_w,n_r,n_S]),[],1+0)*max(0,sqrt(n_w)),[numel(S_q_),1]);
end;%if (numel(unique(n_w_))==1);
%%%%%%%%%%%%%%%%;

