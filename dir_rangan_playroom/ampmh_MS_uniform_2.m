function  ...
[ ...
 parameter ...
,euler_polar_a_ ...
,euler_azimu_b_ ...
,euler_gamma_z_ ...
,image_delta_x_ ...
,image_delta_y_ ...
,image_I_value_ ...
,image_X_value_ ...
,image_S_index_ ...
] = ...
ampmh_MS_uniform_2( ...
 parameter ...
,n_w_max ...
,n_S ...
,viewing_polar_a_all_ ...
,viewing_azimu_b_all_ ...
,n_M ...
,X_wSM___ ...
,delta_x_wSM___ ...
,delta_y_wSM___ ...
,gamma_z_wSM___ ...
,I_value_wSM___ ...
);
%%%%%%%%;
if isempty(parameter);
parameter = struct('type','parameter');
end;%if isempty(parameter);
%%%%%%%%;
euler_polar_a_ = zeros(n_M,1);
euler_azimu_b_ = zeros(n_M,1);
euler_gamma_z_ = zeros(n_M,1);
image_delta_x_ = zeros(n_M,1);
image_delta_y_ = zeros(n_M,1);
image_I_value_ = zeros(n_M,1);
image_X_value_ = zeros(n_M,1);
image_S_index_ = zeros(n_M,1);
%%%%%%%%;
index_nS_permutation_ = randperm(n_S)-1;
%%%%%%%%;
flag_M_used_ = zeros(n_M,1);
tmp_nS=0; nM_sum=0;
while (sum(flag_M_used_)<n_M);
index_M_unused_ = find(flag_M_used_==0)-1;
index_nS = index_nS_permutation_(1+tmp_nS);
if (ndims(X_wSM___)==3);
[image_X_value,index_wM_best] = max(real(X_wSM___(:,1+index_nS,1+index_M_unused_)),[],'all','linear'); index_wM_best = index_wM_best-1;
[nw_best,index_M_best] = ind2sub([n_w_max,numel(index_M_unused_)],1+index_wM_best); 
nw_best = nw_best-1; index_M_best = index_M_best-1;
nM_best = index_M_unused_(1+index_M_best);
gamma_z_best = 2*pi*nw_best/n_w_max; %<-- gamma_z_wSM__ unnecessary. ;
delta_x_best = delta_x_wSM___(1+nw_best,1+index_nS,1+nM_best);
delta_y_best = delta_y_wSM___(1+nw_best,1+index_nS,1+nM_best);
I_value_best = I_value_wSM___(1+nw_best,1+index_nS,1+nM_best);
end;%if (ndims(X_wSM___)==3);
if (ndims(X_wSM___)==2);
[image_X_value,index_M_best] = max(real(X_wSM___(1+index_nS,1+index_M_unused_))); index_M_best = index_M_best-1;
nM_best = index_M_unused_(1+index_M_best);
gamma_z_best = gamma_z_wSM___(1+index_nS,1+nM_best);
delta_x_best = delta_x_wSM___(1+index_nS,1+nM_best);
delta_y_best = delta_y_wSM___(1+index_nS,1+nM_best);
I_value_best = I_value_wSM___(1+index_nS,1+nM_best);
end;%if (ndims(X_wSM___)==2);
flag_M_used_(1+nM_best)=1;
euler_polar_a_(1+nM_best) = viewing_polar_a_all_(1+index_nS);
euler_azimu_b_(1+nM_best) = viewing_azimu_b_all_(1+index_nS);
euler_gamma_z_(1+nM_best) = gamma_z_best;
image_delta_x_(1+nM_best) = delta_x_best;
image_delta_y_(1+nM_best) = delta_y_best;
image_I_value_(1+nM_best) = I_value_best;
image_X_value_(1+nM_best) = image_X_value;
image_S_index_(1+nM_best) = index_nS;
tmp_nS=tmp_nS+1; if (tmp_nS>=n_S); tmp_nS=0; end;
nM_sum=nM_sum+1;
end;%while (sum(flag_M_used_)<n_M);
assert(nM_sum==n_M);
%%%%%%%%;


