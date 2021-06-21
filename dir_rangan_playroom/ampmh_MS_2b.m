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
ampmh_MS_2( ...
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
if (~isfield(parameter,'flag_euler_polar_a_restrict')); parameter.flag_euler_polar_a_restrict = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'euler_polar_a_restrict_band')); parameter.euler_polar_a_restrict_band = pi/16; end; %<-- parameter_bookmark. ;
flag_euler_polar_a_restrict = parameter.flag_euler_polar_a_restrict;
euler_polar_a_restrict_band = parameter.euler_polar_a_restrict_band;
%%%%%%%%;
if flag_euler_polar_a_restrict==1; %<-- only equatorial polar_a allowed. ;
index_polar_a_retain_ = efind(abs(pi/2-viewing_polar_a_all_)<=euler_polar_a_restrict_band);
n_polar_a_retain = numel(index_polar_a_retain_);
while (n_polar_a_retain<=0);
euler_polar_a_restrict_band = euler_polar_a_restrict_band + pi/16;
index_polar_a_retain_ = efind(abs(pi/2-viewing_polar_a_all_)<=euler_polar_a_restrict_band);
n_polar_a_retain = numel(index_polar_a_retain_);
end;%while (n_polar_a_retain<=0);
end;%if flag_euler_polar_a_restrict==1; %<-- only equatorial polar_a allowed. ;
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
if flag_euler_polar_a_restrict==0; %<-- any polar_a allowed. ;
flag_M_used_ = zeros(n_M,1);
tmp_permutation_ = randperm(n_S)-1;
nS=0; nM_sum=0;
while (sum(flag_M_used_)<n_M);
index_M_unused_ = find(flag_M_used_==0)-1;
index_nS = tmp_permutation_(1+nS);
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
nS=nS+1; if (nS>=n_S); nS=0; end;
nM_sum=nM_sum+1;
end;%while (sum(flag_M_used_)<n_M);
assert(nM_sum==n_M);
end;%if flag_euler_polar_a_restrict==0; %<-- any polar_a allowed. ;
%%%%%%%%;
if flag_euler_polar_a_restrict==1; %<-- only equatorial polar_a allowed. ;

end;%if flag_euler_polar_a_restrict==1; %<-- only equatorial polar_a allowed. ;
