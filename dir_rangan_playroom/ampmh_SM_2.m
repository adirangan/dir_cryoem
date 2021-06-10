function  ...
[ ...
 euler_polar_a_ ...
,euler_azimu_b_ ...
,euler_gamma_z_ ...
,image_delta_x_ ...
,image_delta_y_ ...
,image_I_value_ ...
,image_X_value_ ...
,image_S_index_ ...
] = ...
ampmh_SM_2( ...
 f_rand ...
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
euler_polar_a_ = zeros(n_M,1);
euler_azimu_b_ = zeros(n_M,1);
euler_gamma_z_ = zeros(n_M,1);
image_delta_x_ = zeros(n_M,1);
image_delta_y_ = zeros(n_M,1);
image_I_value_ = zeros(n_M,1);
image_X_value_ = zeros(n_M,1);
image_S_index_ = zeros(n_M,1);
%%%%%%%%;
for nM=0:n_M-1;
tmp_image_X_value_ = zeros(n_S,1);
tmp_gamma_z_S_ = zeros(n_S,1);
tmp_delta_x_S_ = zeros(n_S,1);
tmp_delta_y_S_ = zeros(n_S,1);
tmp_I_value_S_ = zeros(n_S,1);
for nS=0:n_S-1;
if (ndims(X_wSM___)==3);
[tmp_X,nw] = max(squeeze(real(X_wSM___(:,1+nS,1+nM)))); nw = nw - 1;
tmp_gamma_z = gamma_z_wSM___(1+nw,1+nS,1+nM);
tmp_delta_x = delta_x_wSM___(1+nw,1+nS,1+nM);
tmp_delta_y = delta_y_wSM___(1+nw,1+nS,1+nM);
tmp_I_value = I_value_wSM___(1+nw,1+nS,1+nM);
end;%if (ndims(X_wSM___)==3);
if (ndims(X_wSM___)==2);
tmp_X = X_wSM___(1+nS,1+nM);
tmp_gamma_z = gamma_z_wSM___(1+nS,1+nM);
tmp_delta_x = delta_x_wSM___(1+nS,1+nM);
tmp_delta_y = delta_y_wSM___(1+nS,1+nM);
tmp_I_value = I_value_wSM___(1+nS,1+nM);
end;%if (ndims(X_wSM___)==2);
tmp_image_X_value_(1+nS) = tmp_X;
tmp_gamma_z_S_(1+nS) = tmp_gamma_z;
tmp_delta_x_S_(1+nS) = tmp_delta_x;
tmp_delta_y_S_(1+nS) = tmp_delta_y;
tmp_I_value_S_(1+nS) = tmp_I_value;
end;%for nS=0:n_S-1;
if (f_rand> 0); 
tmp_X_ij_ = find(tmp_image_X_value_>=prctile(tmp_image_X_value_,100-100*f_rand)) - 1; 
tmp_X_ij_index = max(0,min(numel(tmp_X_ij_)-1,floor(numel(tmp_X_ij_)*rand())));
nS_best = tmp_X_ij_(1+tmp_X_ij_index);
end;%if (f_rand> 0); 
if (f_rand<=0); 
[~,nS_best] = max(tmp_image_X_value_); nS_best = nS_best - 1; 
end;%if (f_rand<=0); 
euler_polar_a_(1+nM) = viewing_polar_a_all_(1+nS_best);
euler_azimu_b_(1+nM) = viewing_azimu_b_all_(1+nS_best);
euler_gamma_z_(1+nM) = tmp_gamma_z_S_(1+nS_best);
image_delta_x_(1+nM) = tmp_delta_x_S_(1+nS_best);
image_delta_y_(1+nM) = tmp_delta_y_S_(1+nS_best);
image_I_value_(1+nM) = tmp_I_value_S_(1+nS_best);
image_X_value_(1+nM) = tmp_image_X_value_(1+nS_best);
image_S_index_(1+nM) = nS_best;
end;%for nM=0:n_M-1;
