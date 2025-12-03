function ...
[ ...
 V_UX_M_lwnM____ ...
] = ...
tpmh_VUXM_lwnM____4( ...
 FTK ...
,n_k_p_r ...
,n_w_ ...
,n_M ...
,M_k_q__ ...
,n_UX_rank ...
,UX_kn__ ...
,X_weight_r_ ...
);

str_thisfunction = 'tpmh_VUXM_lwnM____4';

if nargin<1;
flag_verbose=1;
if (flag_verbose>0); disp(sprintf(' %% testing %s',str_thisfunction)); end;
%%%%;
rng(0);
k_p_r_max = 48/(2*pi); k_eq_d = 1.0/(2.0*pi); str_T_vs_L = 'L';
flag_unif_vs_adap = 0; flag_tensor_vs_adap = 0;
[ ...
 n_k_all ...
,n_k_all_csum_ ...
,k_p_r_all_ ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,weight_3d_k_all_ ...
,weight_shell_k_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,k_c_0_all_ ...
,k_c_1_all_ ...
,k_c_2_all_ ...
,J_node_ ...
,J_weight_ ...
,J_chebfun_ ...
,J_polyval_ ...
,n_polar_a_k_ ...
,polar_a_ka__ ...
,n_azimu_b_ka__ ...
] = ...
sample_sphere_7( ...
 flag_verbose ...
,k_p_r_max ...
,k_eq_d ...
,str_T_vs_L ...
,flag_unif_vs_adap ...
,flag_tensor_vs_adap ...
) ;
%%%%;
delta_r_max = 0.15;
svd_eps = 1e-2;
n_delta_v_requested = 16;
FTK = ...
ampmh_FTK_1( ...
 n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,delta_r_max ...
,svd_eps ...
,n_delta_v_requested ...
);
%%%%;
n_w_max = 2*round(1+2*pi*k_p_r_max);
n_w_ = n_w_max*ones(n_k_p_r,1);
n_w_sum = sum(n_w_);
n_M = 128;
M_k_q__ = randn(n_w_sum,n_M) + i*randn(n_w_sum,n_M);
n_UX_rank = n_k_p_r-1; %<-- just to check dimension;
UX_kn__ = orth(randn(n_k_p_r,n_UX_rank));
X_weight_r_ = rand(n_k_p_r,1);
%%%%;
flag_check=0;
if flag_check;
tmp_t = tic();
V_UX_M_2_lwnM____ = tpmh_VUXM_lwnM____2(FTK,n_k_p_r,n_w_,n_M,M_k_q__,n_UX_rank,UX_kn__,X_weight_r_);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% tpmh_VUXM_lwnM____2: %fs',tmp_t)); end;
end;%if flag_check;
%%%%;
tmp_t = tic();
[ ...
 V_UX_M_4_lwnM____ ...
] = ...
tpmh_VUXM_lwnM____4( ...
 FTK ...
,n_k_p_r ...
,n_w_ ...
,n_M ...
,M_k_q__ ...
,n_UX_rank ...
,UX_kn__ ...
,X_weight_r_ ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% tpmh_VUXM_lwnM____4: %fs',tmp_t)); end;
%%%%;
tmp_t = tic();
[ ...
 V_UX_M_3_lwnM____ ...
] = ...
tpmh_VUXM_lwnM____3( ...
 FTK ...
,n_k_p_r ...
,n_w_ ...
,n_M ...
,M_k_q__ ...
,n_UX_rank ...
,UX_kn__ ...
,X_weight_r_ ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% tpmh_VUXM_lwnM____3: %fs',tmp_t)); end;
fnorm_disp(flag_verbose,'V_UX_M_3_lwnM____',V_UX_M_3_lwnM____,'V_UX_M_4_lwnM____',V_UX_M_4_lwnM____);
disp('returning');return;
end;%if nargin<1;

flag_verbose=0;
n_w_max = max(n_w_); n_w_2 = round(n_w_max/2);
l_max = max(abs(FTK.svd_l_));
%%%%%%%%;

%%%%%%%%;
tmp_t = tic();
V_r__ = reshape(FTK.svd_polyval_V_r_,[FTK.n_svd_l,n_k_p_r]);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% V_r__: %0.6fs',tmp_t)); end;
tmp_t = tic();
V_UX_lrn___ = zeros(FTK.n_svd_l,n_k_p_r,n_UX_rank);
for nUX_rank=0:n_UX_rank-1;
V_UX_lrn___(:,:,1+nUX_rank) = V_r__*diag(UX_kn__(:,1+nUX_rank).*X_weight_r_(:));
end;%for nUX_rank=0:n_UX_rank-1;
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% V_UX_lrn___: %0.6fs',tmp_t)); end;
%%%%%%%%;
%tmp_t = tic();
%index_nw_out__ = cell(1+2*l_max,1);
%index_nw_0in__ = cell(1+2*l_max,1);
%n_nw_ = zeros(1+2*l_max,1);
%for l_shift=-l_max:+l_max;
%index_nw_out_head_ = (0:n_w_2-1-max(0,+l_shift));
%index_nw_0in_head_ = [ n_w_max+l_shift:n_w_max-1 , max(0,+l_shift):n_w_2-1-max(0,-l_shift) ];
%index_nw_out_tail_ = (n_w_2+1+max(0,-l_shift):n_w_max-1);
%index_nw_0in_tail_ = [ n_w_2+1+max(0,+l_shift):n_w_max-1-max(0,-l_shift), 0:+l_shift-1 ];
%index_nw_out_ = [ index_nw_out_head_ , index_nw_out_tail_ ];
%index_nw_0in_ = [ index_nw_0in_head_ , index_nw_0in_tail_ ];
%n_nw = numel(index_nw_out_);
%n_nw_(1+l_max+l_shift) = n_nw;
%index_nw_out__{1+l_max+l_shift} = index_nw_out_;
%index_nw_0in__{1+l_max+l_shift} = index_nw_0in_;
%end;%for l_shift=-l_max:+l_max;
%tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% index_nw_xxx__: %0.6fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
assert(mod(n_w_max,2)==0); %<-- assert that n_w_max is even. ;
index_nw_zerobased_from_centered_ = [ n_w_2:n_w_max-1 , 0:n_w_2-1 ]; %<-- note that we place n_w_2 mode first, ;
index_nw_centered_from_zerobased_ = [ n_w_2:n_w_max-1 , 0:n_w_2-1 ]; %<-- note that we place first mode at n_w_2. ;
%index_nw_centered_out_start_ = zeros(1+2*l_max,1);
%index_nw_centered_out_final_ = zeros(1+2*l_max,1);
%index_nw_centered_0in_start_ = zeros(1+2*l_max,1);
%index_nw_centered_0in_final_ = zeros(1+2*l_max,1);
%for l_shift=-l_max:+l_max;
%index_nw_centered_out_start_(1+l_max+l_shift) = 1+max(0,-l_shift);
%index_nw_centered_out_final_(1+l_max+l_shift) = n_w_max-1+min(0,-l_shift);
%index_nw_centered_0in_start_(1+l_max+l_shift) = 1+max(0,+l_shift);
%index_nw_centered_0in_final_(1+l_max+l_shift) = n_w_max-1+min(0,+l_shift);
%end;%for l_shift=-l_max:+l_max;
%tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% index_nw_centered_xxx_xxxxx__: %0.6fs',tmp_t)); end;
%%%%%%%%;
%tmp_t = tic();
%M_k_q_centered_rwM___ = permute(reshape(M_k_q__,[n_w_max,n_k_p_r,n_M]),[2,1,3]);
%M_k_q_centered_rwM___ = M_k_q_centered_rwM___(:,1+index_nw_centered_from_zerobased_,:);
%tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% M_k_q_centered_rwM___: %0.6fs',tmp_t)); end;
%%%%%%%%;
index_nw_centered_padded_out_start_ = zeros(1+2*l_max,1);
index_nw_centered_padded_out_final_ = zeros(1+2*l_max,1);
index_nw_centered_padded_0in_start_ = zeros(1+2*l_max,1);
index_nw_centered_padded_0in_final_ = zeros(1+2*l_max,1);
for l_shift=-l_max:+l_max;
index_nw_centered_padded_out_start_(1+l_max+l_shift) = 1+l_max;
index_nw_centered_padded_out_final_(1+l_max+l_shift) = n_w_max-1+l_max;
index_nw_centered_padded_0in_start_(1+l_max+l_shift) = 1+l_max+l_shift;
index_nw_centered_padded_0in_final_(1+l_max+l_shift) = n_w_max-1+l_max+l_shift;
end;%for l_shift=-l_max:+l_max;
%M_k_q_centered_padded_wrM___ = cat(1+0,zeros(l_max,n_k_p_r,n_M),circshift(reshape(M_k_q__,[n_w_max,n_k_p_r,n_M]),-n_w_2,1),zeros(l_max,n_k_p_r,n_M));
%M_k_q_centered_padded_wrM___(1+l_max,:,:) = 0;
%M_k_q_centered_padded_rwM___ = permute(M_k_q_centered_padded_wrM___,[2,1,3]);
%%%%%%%%;
%{
tmp_t = tic();
V_UX_M_centered_padded_lwMn____ = zeros(FTK.n_svd_l,n_w_max+2*l_max,n_M,n_UX_rank);
%V_UX_M_centered_lwMn____ = zeros(FTK.n_svd_l,n_w_max,n_M,n_UX_rank);
for nl=0:FTK.n_svd_l-1;
l_shift = FTK.svd_l_(1+nl);
%index_nw_centered_out_start = index_nw_centered_out_start_(1+l_max+l_shift);
%index_nw_centered_out_final = index_nw_centered_out_final_(1+l_max+l_shift);
%index_nw_centered_0in_start = index_nw_centered_0in_start_(1+l_max+l_shift);
%index_nw_centered_0in_final = index_nw_centered_0in_final_(1+l_max+l_shift);
%index_nw_centered_out_ = [index_nw_centered_out_start:index_nw_centered_out_final];
%index_nw_centered_0in_ = [index_nw_centered_0in_start:index_nw_centered_0in_final];
index_nw_centered_padded_out_start = index_nw_centered_padded_out_start_(1+l_max+l_shift);
index_nw_centered_padded_out_final = index_nw_centered_padded_out_final_(1+l_max+l_shift);
index_nw_centered_padded_0in_start = index_nw_centered_padded_0in_start_(1+l_max+l_shift);
index_nw_centered_padded_0in_final = index_nw_centered_padded_0in_final_(1+l_max+l_shift);
index_nw_centered_padded_out_ = [index_nw_centered_padded_out_start:index_nw_centered_padded_out_final];
index_nw_centered_padded_0in_ = [index_nw_centered_padded_0in_start:index_nw_centered_padded_0in_final];
n_nw = numel(index_nw_centered_padded_out_);
for nUX_rank=0:n_UX_rank-1;
for nM=0:n_M-1;
V_UX_M_centered_padded_lwMn____(1+nl,1+index_nw_centered_padded_out_,1+nM,1+nUX_rank) = V_UX_lrn___(1+nl,:,1+nUX_rank)*M_k_q_centered_padded_rwM___(:,1+index_nw_centered_padded_0in_,1+nM) / n_w_max ;
%V_UX_M_centered_lwMn____(1+nl,1+index_nw_centered_out_,1+nM,1+nUX_rank) = V_UX_lrn___(1+nl,:,1+nUX_rank)*M_k_q_centered_rwM___(:,1+index_nw_centered_0in_,1+nM) / n_w_max ;
end;%for nM=0:n_M-1;
end;%for nUX_rank=0:n_UX_rank-1;
end;%for nl=0:FTK.n_svd_l-1;
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% V_UX_M_centered_padded_lwMn____: %0.6fs',tmp_t)); end;
tmp_t = tic();
V_UX_M_lwnM____ = permute(V_UX_M_centered_padded_lwMn____(:,1+l_max+index_nw_zerobased_from_centered_,:,:),[1,2,4,3]);
%V_UX_M_lwnM____ = permute(V_UX_M_centered_lwMn____(:,1+index_nw_zerobased_from_centered_,:,:),[1,2,4,3]);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% V_UX_M_lwMn____: %0.6fs',tmp_t)); end;
%}

tmp_t = tic();
V_UX_nrl___ = permute(V_UX_lrn___,[3,2,1]);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% V_UX_rnl___: %0.6fs',tmp_t)); end;
tmp_t = tic();
M_k_q_centered_padded_rwM___ = cat(1+0,zeros(l_max,n_k_p_r,n_M),circshift(reshape(M_k_q__,[n_w_max,n_k_p_r,n_M]),-n_w_2,1),zeros(l_max,n_k_p_r,n_M));
M_k_q_centered_padded_rwM___(1+l_max,:,:) = 0;
M_k_q_centered_padded_rwM___ = permute(M_k_q_centered_padded_rwM___,[2,1,3]);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% M_k_q_centered_padded_rwM___: %0.6fs',tmp_t)); end;
tmp_t = tic();
V_UX_M_centered_nwMl____ = zeros(n_UX_rank,n_w_max,n_M,FTK.n_svd_l);
for nl=0:FTK.n_svd_l-1;
l_shift = FTK.svd_l_(1+nl);
index_nw_centered_padded_out_start = index_nw_centered_padded_out_start_(1+l_max+l_shift);
index_nw_centered_padded_out_final = index_nw_centered_padded_out_final_(1+l_max+l_shift);
index_nw_centered_padded_0in_start = index_nw_centered_padded_0in_start_(1+l_max+l_shift);
index_nw_centered_padded_0in_final = index_nw_centered_padded_0in_final_(1+l_max+l_shift);
index_nw_centered_padded_out_ = [index_nw_centered_padded_out_start:index_nw_centered_padded_out_final];
index_nw_centered_padded_0in_ = [index_nw_centered_padded_0in_start:index_nw_centered_padded_0in_final];
V_UX_M_centered_nwMl____(:,2:end,:,1+nl) = pagemtimes(reshape(V_UX_nrl___(:,:,1+nl),[n_UX_rank,n_k_p_r,1])/n_w_max,reshape(M_k_q_centered_padded_rwM___(:,1+index_nw_centered_padded_0in_,:),[n_k_p_r,n_w_max-1,n_M]))  ;
end;%for nl=0:FTK.n_svd_l-1;
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% V_UX_M_centered_padded_lwMn____: %0.6fs',tmp_t)); end;
tmp_t = tic();
V_UX_M_lwnM____ = permute(V_UX_M_centered_nwMl____(:,1+index_nw_zerobased_from_centered_,:,:),[4,2,1,3]);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% V_UX_M_lwMn____: %0.6fs',tmp_t)); end;

%fnorm_disp(1,'A',V_UX_M_centered_lwMn____(:,1+[0:n_w_max-1],:,:),'B',V_UX_M_centered_padded_lwMn____(:,1+l_max+[0:n_w_max-1],:,:));


