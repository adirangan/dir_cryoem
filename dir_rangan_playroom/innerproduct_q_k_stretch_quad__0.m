function C_q_wSM___ = innerproduct_q_k_stretch_quad__0(n_r,weight_p_r_,n_w_,n_S,S_q__,S_q_rSw___,n_M,M_q__,M_q_rMw___) ;
% Assumes quasi-uniform polar-grid. ;
% Assumes that n_w_ is a list of even integers. ;
if (nargin<7);
n_r = 6;
weight_p_r_ = rand(n_r,1);
n_w_ = 2*ceil(transpose(linspace(15,50,n_r)));
n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_S = 10; S_q__ = randn(n_w_sum,n_S) + i*randn(n_w_sum,n_S);
n_M = 12; M_q__ = randn(n_w_sum,n_M) + i*randn(n_w_sum,n_M);
C_q_wSM_0___ = innerproduct_q_k_stretch_quad__0(n_r,weight_p_r_,n_w_,n_S,S_q__,[],n_M,M_q__,[]);
C_q_wSM_1___ = zeros(n_w_max,n_S,n_M);
for nS=0:n_S-1;
for nM=0:n_M-1;
C_q_wSM_1___(:,1+nS,1+nM) = innerproduct_q_k_stretch_quad_0(n_r,[],weight_p_r_,n_w_,[],S_q__(:,1+nS),M_q__(:,1+nM));
end;%for nM=0:n_M-1;
end;%for nS=0:n_S-1;
disp(sprintf(' %% C_q_wSM_0___ vs C_q_wSM_1___: %0.16f',fnorm(C_q_wSM_0___-C_q_wSM_1___)/fnorm(C_q_wSM_0___)));
disp('returning'); return;
end;%if (nargin<7);

n_w_max = max(n_w_);
if isempty(S_q_rSw___); S_q_rSw___ = permute(innerproduct_q_k_stretch_quad_stack__0(n_r,weight_p_r_,n_w_,n_S,S_q__),[1,3,2]); end;
if isempty(M_q_rMw___); M_q_rMw___ = permute(innerproduct_q_k_stretch_quad_stack__0(n_r,weight_p_r_,n_w_,n_M,M_q__),[1,3,2]); end;
C_q_SMw___ = zeros(n_S,n_M,n_w_max);
for nw=0:n_w_max-1;
C_q_SMw___(:,:,1+nw) = ctranspose(S_q_rSw___(:,:,1+nw))*M_q_rMw___(:,:,1+nw);
end;%for nw=0:n_w_max-1;
C_q_wSM___ = permute(C_q_SMw___,[3,1,2]);
