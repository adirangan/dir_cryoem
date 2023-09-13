function ...
R_x_p_ = ...
radon_k_p_to_x_p_xxnufft( ...
 n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,S_k_p_ ...
);
%%%%%%%%;
verbose=0;
n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
if (n_w_sum~=n_w_max*n_k_p_r); disp(sprintf(' %% Warning, incompatible n_w_ in radon_k_p_to_x_p_xxnufft')); end;
if mod(n_w_max,2)~=0; disp(sprintf(' %% Warning, incompatible n_w_max in radon_k_p_to_x_p_xxnufft')); end;
n_w_2 = n_w_max/2;
[ ...
 weight_1d_k_p_r_ ...
] = ...
get_weight_1d_0( ...
 verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
);
weight_1d_k_p_twin_j_ = [flipud(weight_1d_k_p_r_(:));weight_1d_k_p_r_];
S_k_p_wk__ = reshape(S_k_p_,[n_w_max,n_k_p_r]);
k_p_r_twin_j_ = [-flipud(k_p_r_(:));+k_p_r_(:)]; n_j = numel(k_p_r_twin_j_);
x_p_r_twin_j_ = k_p_r_twin_j_/k_p_r_max;
S_k_p_twin_wj__ = zeros(n_w_max,n_j);
for nw=0:n_w_max-1;
nw0 = nw;
nw1 = periodize(nw0+n_w_2,0,n_w_max);
S_k_p_twin_wj__(1+nw,1 + 0*n_k_p_r + [0:n_k_p_r-1]) = fliplr(S_k_p_wk__(1+nw1,:));
S_k_p_twin_wj__(1+nw,1 + 1*n_k_p_r + [0:n_k_p_r-1]) = S_k_p_wk__(1+nw0,:);
end;%for nw=0:n_w_max-1;
R_x_p_twin_wj__ = zeros(n_w_max,n_j);
for nw=0:n_w_max-1;
R_x_p_twin_wj__(1+nw,:) = xxnufft1d3(n_j,2*pi*k_p_r_twin_j_,S_k_p_twin_wj__(1+nw,:).*transpose(weight_1d_k_p_twin_j_),+1,1e-12,n_j,x_p_r_twin_j_);
end;%for nw=0:n_w_max-1;
R_x_p_wk__ = zeros(n_w_max,n_k_p_r);
R_x_p_wk__ = R_x_p_twin_wj__(:,1+n_k_p_r+[0:n_k_p_r-1]);
R_x_p_ = reshape(R_x_p_wk__,[n_w_sum,1]);
