function ...
[ ...
 parameter ...
,a_UCTF_UX_Y_ync__ ...
,n_quad_from_data_q_ ...
] = ...
qbp_pm_1( ...
 parameter ...
,pm_n_k_p_r ...
,pm_l_max_ ...
,pm_n_w_ ...
,n_M ...
,UX_M_k_p_wnM__ ...
,n_CTF_rank ...
,VSCTF_Mc__ ...
,euler_polar_a_ ...
,euler_azimu_b_ ...
,euler_gamma_z_ ...
,image_I_value_...
);
%%%%%%%%;
% Applies quadrature-based back-propagation to solve for a_UX_Y_. ;
% Associates CTF_k_p_r__(:,1+nM) with image M_k_p__(:,1+nM);
% ;
% Input: ;
% quad_k_eq_d: real equatorial distance used for determining quadrature nodes on sphere (radius assumed to be 1). ;
% pm_n_UX_rank = pm_n_k_p_r: integer number of principal-volume shells retained. ;
% pm_l_max_: integer array of size pm_n_k_p_r. pm_l_max_(1+pm_nk_p_r) is the order used for a_UCTF_UX_Y_ync__ on principal-volume-shell pm_nk_p_r. ;
% pm_n_w_: integer array of size pm_n_k_p_r. pm_n_w_(1+pm_nk_p_r) is the number of inplane_gamma_z values recorded at that principal-image-ring. ;
% n_M: integer number of images. ;
% UX_M_k_p_wnM__: complex array of size (pm_n_w_sum,n_M). stack of principal-images in k_p_ format. ;
% We assume that each column of UX_M_k_p_wnM__ corresponds to a single principal-image, which is itself a stack of principal-image-rings. ;
% n_CTF_rank: integer number of CTF ranks to consider. ;
% (implicit) UCTF_kc__: real array of size (n_k_p_r,n_CTF_rank). ;
% VSCTF_Mc__: real array of size (n_M,n_CTF_rank). ;
% We assume that the full (isotropic) CTF-function CTF_k_p_r__, given by: ;
% (implicit) CTF_k_p_r__: real array of size (n_k_p_r,n_M). ;
% can be approximated via: ;
% CTF_k_p_r__ = UCTF_kc__*transpose(VSCTF_Mc__);
% which is a low-rank approximation with n_CTF_rank terms. ;
% If n_CTF_rank<=0 or VSCTF_Mc__ is empty, we assume that: ;
% n_CTF_rank=1; UCTF_kc__ = ones(pm_n_k_p_r,1); VSCTF_Mc__ = ones(n_M,1);
% euler_polar_a_: real array of size n_M. polar_a used for each image ;
% euler_azimu_b_: real array of size n_M. azimu_b used for each image ;
% euler_gamma_z_: real array of size n_M. gamma_z used for each image ;
% image_I_value_: real array of size n_M. I_value used for each image ;
% ;
% Output: ;
% a_UCTF_UX_Y_ync__: complex array of size (pm_n_lm_sum,n_UCTF_rank). output functions in k_Y_ format. ;
% This output function should approximately satisfy the least-square problem: ;
% \sum_{nCTF_rank=0}^{n_CTF_rank-1} S * [ \tau_{1+nM} * VSCTF_Mc__(1+nM,1+nCTF_rank) ] * a_UCTF_UX_Y_ync___(:,1+pm_nUX_rank,1+nCTF_rank) = [ UX_M_k_p_wnM___(:,1+pm_nUX_rank,1+nM) ] \forall nM \in [0,\ldots,n_M-1] and \forall pm_nUX_rank \in [0,pm_n_UX_rank-1]. ;
% where : ;
% \tau_{1+nM} corresponds to rotation by the viewing-angle associated with image nM, and ;
% S is the template-operator (i.e., equatorial-evaluation), and ;
% UX_M_k_p_wnM___(:,1+pm_nUX_rank,1+nM) = UX_M_k_p_wnM__(1+pm_n_w_csum_(1+pm_nUX_rank) + (0:pm_n_w_(1+pm_nUX_rank)-1),1+nM), and ;
% a_UCTF_UX_Y_ync___(:,1+pm_nUX_rank,1+nCTF_rank) = \sum_{nk_p_r=0}^{n_k_p_r} UCTF_(1+nk_p_r,1+nCTF_rank) * UX_(1+nk_p_r,1+pm_nUX_rank) * a_UX_Y__(:,1+nk_p_r). ;
% ;
% Rather than using least-squares to solve this problem, we instead generate a quadrature-grid on the sphere, ;
% map each data-point to its closest quadrature-gridpoint, and then numerically integrate to recover a_UCTF_UX_Y_ync__. ;
% ;
% To account for VSCTF_Mc__ when n_CTF_rank==1 we solve, for each local-neighborhood associated with a particular quadrature point, ;
% the following local-least-squares-problem: ;
% \argmin_{a} \sum_{nlocal} \| UX_M_I_wM__(nlocal) - VSCTF_wMc__(nlocal,1+0)*a \|^2, ;
% which has a solution: ;
% a = \{ \sum_{nlocal} VSCTF_wMc__(nlocal,1+0)*UX_M_I_wM__(nlocal) \} / \{ \sum_{nlocal} VSCTF_wMc__(nlocal,1+0)*VSCTF_wMc__(nlocal,1+0) \}
% where UX_M_I_wM__(nlocal) refers to one of the image-points UX_M_k_p_wnM__ in the local-neighborhood of the quadrature-point, ;
% and VSCTF_wMc__(nlocal,1+0) refers to the VSCTF_Mc__(?,1+0) associated with that image-point. ;
% ;
% When (n_CTF_rank> 1) & (flag_qbp_pm_solve_in_CTF_sequence==1), ;
% We solve for the dominant CTF-modes of a first, then use an approximate residual to solve for the second (and subsequent) CTF-modes. ;
% ;
% When (n_CTF_rank> 1) & (flag_qbp_pm_solve_in_CTF_sequence==0), ;
% We solve for all the CTF-modes simultaneously by using the local-least-squares-problem: ;
% \argmin_{a_} \sum_{nlocal} \| UX_M_I_wM__(nlocal) - \sum_{nCTF_rank} VSCTF_wMc__(nlocal,1+nCTF_rank)*a_(1+nCTF_rank) \|^2, ;
% which has a solution: ;
% \sum_{nlocal} VSCTF_wMc__(nlocal,1+nCTF_rank_1)*UX_M_I_wM__(nlocal) = \sum_{nCTF_rank_0} VSCTF2_cc__(1+nCTF_rank_1,1+nCTF_rank_0)*a_(1+nCTF_rank_0), ;
% where: ;
% VSCTF2_cc__(1+nCTF_rank_1,1+nCTF_rank_0) = \sum_{nlocal} VSCTF_wMc__(nlocal,1+nCTF_rank_1)*VSCTF_wMc__(nlocal,1+nCTF_rank_0). ;
% ;
%%%%%%%%;

if (nargin<1);
disp(sprintf(' %% [testing qbp_pm_1]'));
n_CTF_rank = 2;
pm_n_UX_rank = 3;
pm_l_max = 32; pm_l_max_ = pm_l_max*ones(pm_n_UX_rank,1);
pm_n_lm = (1+pm_l_max)^2; pm_n_lm_ = pm_n_lm*ones(pm_n_UX_rank,1); pm_n_lm_sum = sum(pm_n_lm_);
pm_n_w = 160*1; pm_n_w_ = pm_n_w*ones(pm_n_UX_rank,1); pm_n_w_sum = sum(pm_n_w_);
pm_k_p_r_ = ones(pm_n_UX_rank,1); pm_k_p_r_max = 1; pm_weight_k_p_r_ = ones(pm_n_UX_rank,1);
template_viewing_k_eq_d = 1/(2*pi)/2;
template_k_eq_d = -1;
quad_k_eq_d = template_viewing_k_eq_d/4;
%%%%%%%%;
rng(1);
a_true_k_Y_ync__ = crandn(pm_n_lm_sum,n_CTF_rank);
for nCTF_rank=0:n_CTF_rank-1;
[ ...
 S_k_p_wnSc___(:,:,1+nCTF_rank) ...
,pm_n_w_ ...
,~ ...
,~ ...
,n_S ...
,viewing_azimu_b_all_ ...
,viewing_polar_a_all_ ...
,viewing_weight_all_ ...
] = ...
get_template_1( ...
 0 ...
,pm_n_UX_rank ...
,pm_k_p_r_ ...
,pm_k_p_r_max ...
,pm_weight_k_p_r_ ...
,pm_l_max_ ...
,a_true_k_Y_ync__(:,1+nCTF_rank) ...
,template_viewing_k_eq_d ...
,-1 ...
,pm_n_w_ ...
);
end;%for nCTF_rank=0:n_CTF_rank-1;
%%%%%%%%;
S_CTF_k_p_wnSc___ = zeros(pm_n_w_sum,n_S,n_CTF_rank);
CTF_Sc__ = 0.5 + 0.5*rand(n_S,n_CTF_rank);
flag_decay=1;
if flag_decay;
decay_step = 0.125;
for nCTF_rank=0:n_CTF_rank-1;
CTF_Sc__(:,1+nCTF_rank) = (decay_step.^nCTF_rank)*CTF_Sc__(:,1+nCTF_rank);
end;%for nCTF_rank=0:n_CTF_rank-1;
end;%if flag_decay;
for nCTF_rank=0:n_CTF_rank-1;
S_CTF_k_p_wnSc___(:,:,1+nCTF_rank) = S_k_p_wnSc___(:,:,1+nCTF_rank)*diag(CTF_Sc__(1+nCTF_rank));
end;%for nCTF_rank=0:n_CTF_rank-1;
n_M = n_S;
N_k_p_brut_wnM__ = sum(S_CTF_k_p_wnSc___,3);
%%%%%%%%;
tmp_t = tic();
parameter = struct('type','parameter');
parameter.template_viewing_k_eq_d = template_viewing_k_eq_d;
parameter.quad_k_eq_d = quad_k_eq_d;
parameter.flag_qbp_pm_solve_in_CTF_sequence = 0;
[ ...
 parameter ...
 ,a_0qbp_k_Y_ync__ ...
 ,n_quad_from_data_q_ ...
] = ...
qbp_pm_1( ...
 parameter ...
,pm_n_UX_rank ...
,pm_l_max_ ...
,pm_n_w_ ...
,n_S ...
,N_k_p_brut_wnM__ ...
,n_CTF_rank ...
,CTF_Sc__ ...
,viewing_polar_a_all_ ...
,viewing_azimu_b_all_ ...
,zeros(n_S,1) ...
,ones(n_S,1)...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% pm_l_max %d, n_S: %d, quad_k_eq_d %0.3f, flag_qbp_pm_solve_in_CTF_sequence %d, qbp_pm_1: %0.3fs',pm_l_max,n_S,parameter.quad_k_eq_d,parameter.flag_qbp_pm_solve_in_CTF_sequence,tmp_t));
%%%%;
tmp_t = tic();
N_k_p_0qbp_wnMc___ = zeros(pm_n_w_sum,n_M,n_CTF_rank);
for nCTF_rank=0:n_CTF_rank-1;
[ ...
 N_k_p_0qbp_wnMc___(:,:,1+nCTF_rank) ...
] = ...
get_template_1( ...
 0 ...
,pm_n_UX_rank ...
,ones(pm_n_UX_rank,1) ...
,1 ...
,ones(pm_n_UX_rank,1) ...
,pm_l_max*ones(pm_n_UX_rank,1) ...
,a_0qbp_k_Y_ync__(:,1+nCTF_rank) ...
,template_viewing_k_eq_d/pm_k_p_r_max ...
,template_k_eq_d ...
,pm_n_w*ones(pm_n_UX_rank,1) ...
);
end;%for nCTF_rank=0:n_CTF_rank-1;
tmp_t = toc(tmp_t); disp(sprintf(' %% N_k_p_0qbp_wnMc___: %0.3fs',tmp_t));
N_k_p_0qbp_wnM__ = zeros(pm_n_w_sum,n_M);
for nCTF_rank=0:n_CTF_rank-1;
for nM=0:n_M-1;
N_k_p_0qbp_wnM__(:,1+nM) = N_k_p_0qbp_wnM__(:,1+nM) + N_k_p_0qbp_wnMc___(:,1+nM,1+nCTF_rank)*CTF_Sc__(1+nM,1+nCTF_rank);
end;%for nM=0:n_M-1;
end;%for nCTF_rank=0:n_CTF_rank-1;
disp(sprintf(' %% real(corr(N_k_p_0qbp_wnM__(:),N_k_p_brut_wnM__(:))) %0.16f',real(corr(N_k_p_0qbp_wnM__(:),N_k_p_brut_wnM__(:)))));
disp(sprintf(' %% N_k_p_brut_wnM__ vs N_k_p_0qbp_wnM__: %0.16f',fnorm(N_k_p_brut_wnM__ - N_k_p_0qbp_wnM__)/fnorm(N_k_p_brut_wnM__)));
%%%%%%%%;
tmp_t = tic();
parameter = struct('type','parameter');
parameter.template_viewing_k_eq_d = template_viewing_k_eq_d;
parameter.quad_k_eq_d = quad_k_eq_d;
parameter.flag_qbp_pm_solve_in_CTF_sequence = 1;
[ ...
 parameter ...
 ,a_1qbp_k_Y_ync__ ...
 ,n_quad_from_data_q_ ...
] = ...
qbp_pm_1( ...
 parameter ...
,pm_n_UX_rank ...
,pm_l_max_ ...
,pm_n_w_ ...
,n_S ...
,N_k_p_brut_wnM__...
,n_CTF_rank ...
,CTF_Sc__ ...
,viewing_polar_a_all_ ...
,viewing_azimu_b_all_ ...
,zeros(n_S,1) ...
,ones(n_S,1)...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% pm_l_max %d, n_S: %d, quad_k_eq_d %0.3f, flag_qbp_pm_solve_in_CTF_sequence %d, qbp_pm_1: %0.3fs',pm_l_max,n_S,parameter.quad_k_eq_d,parameter.flag_qbp_pm_solve_in_CTF_sequence,tmp_t));
%%%%;
tmp_t = tic();
N_k_p_1qbp_wnMc___ = zeros(pm_n_w_sum,n_M,n_CTF_rank);
for nCTF_rank=0:n_CTF_rank-1;
[ ...
 N_k_p_1qbp_wnMc___(:,:,1+nCTF_rank) ...
] = ...
get_template_1( ...
 0 ...
,pm_n_UX_rank ...
,ones(pm_n_UX_rank,1) ...
,1 ...
,ones(pm_n_UX_rank,1) ...
,pm_l_max*ones(pm_n_UX_rank,1) ...
,a_1qbp_k_Y_ync__(:,1+nCTF_rank) ...
,template_viewing_k_eq_d/pm_k_p_r_max ...
,template_k_eq_d ...
,pm_n_w*ones(pm_n_UX_rank,1) ...
);
end;%for nCTF_rank=0:n_CTF_rank-1;
tmp_t = toc(tmp_t); disp(sprintf(' %% N_k_p_1qbp_wnMc___: %0.3fs',tmp_t));
N_k_p_1qbp_wnM__ = zeros(pm_n_w_sum,n_M);
for nCTF_rank=0:n_CTF_rank-1;
for nM=0:n_M-1;
N_k_p_1qbp_wnM__(:,1+nM) = N_k_p_1qbp_wnM__(:,1+nM) + N_k_p_1qbp_wnMc___(:,1+nM,1+nCTF_rank)*CTF_Sc__(1+nM,1+nCTF_rank);
end;%for nM=0:n_M-1;
end;%for nCTF_rank=0:n_CTF_rank-1;
disp(sprintf(' %% real(corr(N_k_p_1qbp_wnM__(:),N_k_p_brut_wnM__(:))) %0.16f',real(corr(N_k_p_1qbp_wnM__(:),N_k_p_brut_wnM__(:)))));
disp(sprintf(' %% N_k_p_brut_wnM__ vs N_k_p_1qbp_wnM__: %0.16f',fnorm(N_k_p_brut_wnM__ - N_k_p_1qbp_wnM__)/fnorm(N_k_p_brut_wnM__)));
%%%%%%%%;
tmp_t = tic();
parameter = struct('type','parameter');
parameter.template_viewing_k_eq_d = template_viewing_k_eq_d;
parameter.cg_lsq_pcg_tol = 1e-4;
parameter.cg_lsq_pcg_maxit = 512;
[ ...
 parameter ...
,a_0lsq_k_Y_ync__ ...
,Residual_wnM__ ...
,pcg_flag_n_ ...
,pcg_relres_n_ ...
,pcg_iter_n_ ...
,pcg_resvec_ni__ ...
] = ...
cg_lsq_pm_2( ...
 parameter ...
,pm_n_UX_rank ...
,pm_l_max_ ...
,pm_n_w_ ...
,n_S ...
,N_k_p_brut_wnM__...
,n_CTF_rank ...
,CTF_Sc__ ...
,viewing_polar_a_all_ ...
,viewing_azimu_b_all_ ...
,zeros(n_S,1) ...
,ones(n_S,1) ...
);
%disp(sprintf(' %% pcg_flag_n_: ')); disp(pcg_flag_n_);
%disp(sprintf(' %% pcg_relres_n_: ')); disp(pcg_relres_n_);
%disp(sprintf(' %% pcg_iter_n_: ')); disp(pcg_iter_n_);
tmp_t = toc(tmp_t); disp(sprintf(' %% cg_lsq_pm_2: %0.3fs',tmp_t));
%%%%;
tmp_t = tic();
N_k_p_0lsq_wnMc___ = zeros(pm_n_w_sum,n_M,n_CTF_rank);
for nCTF_rank=0:n_CTF_rank-1;
[ ...
 N_k_p_0lsq_wnMc___(:,:,1+nCTF_rank) ...
] = ...
get_template_1( ...
 0 ...
,pm_n_UX_rank ...
,ones(pm_n_UX_rank,1) ...
,1 ...
,ones(pm_n_UX_rank,1) ...
,pm_l_max*ones(pm_n_UX_rank,1) ...
,a_0lsq_k_Y_ync__(:,1+nCTF_rank) ...
,template_viewing_k_eq_d/pm_k_p_r_max ...
,template_k_eq_d ...
,pm_n_w*ones(pm_n_UX_rank,1) ...
);
end;%for nCTF_rank=0:n_CTF_rank-1;
tmp_t = toc(tmp_t); disp(sprintf(' %% N_k_p_0lsq_wnMc___: %0.3fs',tmp_t));
N_k_p_0lsq_wnM__ = zeros(pm_n_w_sum,n_M);
for nCTF_rank=0:n_CTF_rank-1;
for nM=0:n_M-1;
N_k_p_0lsq_wnM__(:,1+nM) = N_k_p_0lsq_wnM__(:,1+nM) + N_k_p_0lsq_wnMc___(:,1+nM,1+nCTF_rank)*CTF_Sc__(1+nM,1+nCTF_rank);
end;%for nM=0:n_M-1;
end;%for nCTF_rank=0:n_CTF_rank-1;
disp(sprintf(' %% real(corr(N_k_p_0lsq_wnM__(:),N_k_p_brut_wnM__(:))) %0.16f',real(corr(N_k_p_0lsq_wnM__(:),N_k_p_brut_wnM__(:)))));
disp(sprintf(' %% N_k_p_brut_wnM__ vs N_k_p_0lsq_wnM__: %0.16f',fnorm(N_k_p_brut_wnM__ - N_k_p_0lsq_wnM__)/fnorm(N_k_p_brut_wnM__)));
%%%%%%%%;
figure(1);clf;figsml; plot(n_quad_from_data_q_,'.'); title('n_quad_from_data_q_','Interpreter','none');
figure(2);clf;figbig;
for nCTF_rank=0:n_CTF_rank-1;
subplot(n_CTF_rank,6,1+0+nCTF_rank*6);
plot(real(a_true_k_Y_ync__(:,1+nCTF_rank)),real(a_0qbp_k_Y_ync__(:,1+nCTF_rank)),'.',[-2,+2],[-2,+2],'k-'); xlabel('true');ylabel('0qbp');title(sprintf('real flag %d nCTF_rank %d',0,nCTF_rank),'Interpreter','none');axis(2*[-1,1,-1,1]);axis square;
subplot(n_CTF_rank,6,1+1+nCTF_rank*6);
plot(imag(a_true_k_Y_ync__(:,1+nCTF_rank)),imag(a_0qbp_k_Y_ync__(:,1+nCTF_rank)),'.',[-2,+2],[-2,+2],'k-'); xlabel('true');ylabel('0qbp');title(sprintf('imag flag %d nCTF_rank %d',0,nCTF_rank),'Interpreter','none');axis(2*[-1,1,-1,1]);axis square;
subplot(n_CTF_rank,6,1+2+nCTF_rank*6);
plot(real(a_true_k_Y_ync__(:,1+nCTF_rank)),real(a_1qbp_k_Y_ync__(:,1+nCTF_rank)),'.',[-2,+2],[-2,+2],'k-'); xlabel('true');ylabel('1qbp');title(sprintf('real flag %d nCTF_rank %d',1,nCTF_rank),'Interpreter','none');axis(2*[-1,1,-1,1]);axis square;
subplot(n_CTF_rank,6,1+3+nCTF_rank*6);
plot(imag(a_true_k_Y_ync__(:,1+nCTF_rank)),imag(a_1qbp_k_Y_ync__(:,1+nCTF_rank)),'.',[-2,+2],[-2,+2],'k-'); xlabel('true');ylabel('1qbp');title(sprintf('imag flag %d nCTF_rank %d',1,nCTF_rank),'Interpreter','none');axis(2*[-1,1,-1,1]);axis square;
subplot(n_CTF_rank,6,1+4+nCTF_rank*6);
plot(real(a_true_k_Y_ync__(:,1+nCTF_rank)),real(a_0lsq_k_Y_ync__(:,1+nCTF_rank)),'.',[-2,+2],[-2,+2],'k-'); xlabel('true');ylabel('0lsq');title(sprintf('real flag %d nCTF_rank %d',1,nCTF_rank),'Interpreter','none');axis(2*[-1,1,-1,1]);axis square;
subplot(n_CTF_rank,6,1+5+nCTF_rank*6);
plot(imag(a_true_k_Y_ync__(:,1+nCTF_rank)),imag(a_0lsq_k_Y_ync__(:,1+nCTF_rank)),'.',[-2,+2],[-2,+2],'k-'); xlabel('true');ylabel('0lsq');title(sprintf('imag flag %d nCTF_rank %d',1,nCTF_rank),'Interpreter','none');axis(2*[-1,1,-1,1]);axis square;
end;%for nCTF_rank=0:n_CTF_rank-1;
%%%%%%%%;
disp(sprintf(' %% corr(a_true_k_Y_ync__,a_0qbp_k_Y_ync__): %0.16f',corr(a_true_k_Y_ync__(:),a_0qbp_k_Y_ync__(:))));
disp(sprintf(' %% a_true_k_Y_ync__ vs a_0qbp_k_Y_ync__: %0.16f',fnorm(a_true_k_Y_ync__(:) - a_0qbp_k_Y_ync__(:))/fnorm(a_true_k_Y_ync__(:))));
disp(sprintf(' %% corr(a_true_k_Y_ync__,a_1qbp_k_Y_ync__): %0.16f',corr(a_true_k_Y_ync__(:),a_1qbp_k_Y_ync__(:))));
disp(sprintf(' %% a_true_k_Y_ync__ vs a_1qbp_k_Y_ync__: %0.16f',fnorm(a_true_k_Y_ync__(:) - a_1qbp_k_Y_ync__(:))/fnorm(a_true_k_Y_ync__(:))));
disp(sprintf(' %% corr(a_true_k_Y_ync__,a_0lsq_k_Y_ync__): %0.16f',corr(a_true_k_Y_ync__(:),a_0lsq_k_Y_ync__(:))));
disp(sprintf(' %% a_true_k_Y_ync__ vs a_0lsq_k_Y_ync__: %0.16f',fnorm(a_true_k_Y_ync__(:) - a_0lsq_k_Y_ync__(:))/fnorm(a_true_k_Y_ync__(:))));
%%%%%%%%;
disp('returning');return;
end;%if (nargin<1);

verbose = 0;
if (verbose>0); disp(sprintf(' %% [entering qbp_pm_1]')); end;
pm_n_UX_rank = pm_n_k_p_r;

if isempty(parameter);
parameter = struct('type','parameter');
end;%if isempty(parameter);
%%%%%%%%;
if (~isfield(parameter,'tolerance_master')); parameter.tolerance_master = 1e-2; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_qbp_pm_solve_in_CTF_sequence')); parameter.flag_qbp_pm_solve_in_CTF_sequence = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'qbp_eps')); parameter.qbp_eps = parameter.tolerance_master*1e-1; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'quad_k_eq_d')); 
pm_n_lm_max = (1+max(pm_l_max_))^2; quad_k_eq_d_A = sqrt(4*pi/max(pm_n_lm_max));
if (~isfield(parameter,'template_viewing_k_eq_d'));
quad_k_eq_d_B = 1/(2*pi)/2;
end;%if (~isfield(parameter,'template_viewing_k_eq_d'));
if ( isfield(parameter,'template_viewing_k_eq_d'));
quad_k_eq_d_B = parameter.template_viewing_k_eq_d/2;
end;%if ( isfield(parameter,'template_viewing_k_eq_d'));
quad_k_eq_d = min(quad_k_eq_d_A,quad_k_eq_d_B);
parameter.quad_k_eq_d = quad_k_eq_d; %<-- parameter_bookmark. ;
end;%if (~isfield(parameter,'quad_k_eq_d')); 
%%%%%%%%;
flag_qbp_pm_solve_in_CTF_sequence = parameter.flag_qbp_pm_solve_in_CTF_sequence;
quad_k_eq_d = parameter.quad_k_eq_d;
qbp_eps = parameter.qbp_eps;

[ ...
 parameter ...
,pm_n_UX_rank ...
,pm_n_w_ ...
,n_M ...
,UX_M_k_p_wnM__ ...
,n_CTF_rank ...
,VSCTF_Mc__ ...
] = ...
cg_lsq_pm_reduce_1( ...
 parameter ...
,pm_n_UX_rank ...
,pm_n_w_ ...
,n_M ...
,UX_M_k_p_wnM__ ...
,n_CTF_rank ...
,VSCTF_Mc__ ...
);

if (isempty(image_I_value_)); image_I_value_ = ones(n_M,1); end;

pm_n_UX_rank = pm_n_k_p_r;
pm_l_max_ = pm_l_max_(1:pm_n_UX_rank);
pm_n_lm_ = (1+pm_l_max_).^2;
pm_n_lm_max = max(pm_n_lm_);
pm_n_lm_sum = sum(pm_n_lm_);
pm_n_lm_csum_ = cumsum([0;pm_n_lm_]);

pm_n_w_max = max(pm_n_w_(1:pm_n_k_p_r));
pm_n_w_sum = sum(pm_n_w_(1:pm_n_k_p_r));
pm_n_w_csum_ = cumsum([0;pm_n_w_(1:pm_n_k_p_r)]);

if (n_CTF_rank<=0); n_CTF_rank = 1; VSCTF_Mc__ = ones(n_M,1); end;
if (isempty(VSCTF_Mc__)); n_CTF_rank = 1; VSCTF_Mc__ = ones(n_M,1); end;

[ ...
 quad_n_all ...
,quad_azimu_b_all_ ...
,quad_polar_a_all_ ...
,quad_weight_all_ ...
,quad_k_c_0_all_ ...
,quad_k_c_1_all_ ...
,quad_k_c_2_all_ ...
,~ ...
,~ ...
,~ ...
] = ...
sample_shell_5( ...
 1.0 ...
,quad_k_eq_d ...
,'L' ...
) ;
quad_k_c_qd__ = [ quad_k_c_0_all_ , quad_k_c_1_all_ , quad_k_c_2_all_ ];

flag_unique_pm_n = 0;
if (numel(unique(pm_l_max_))==1 & numel(unique(pm_n_lm_))==1 & numel(unique(pm_n_w_))==1);
flag_unique_pm_n = 1;
pm_l_max = pm_l_max_(1+0);
Ylm__ = get_Ylm__(1+pm_l_max,0:pm_l_max,quad_n_all,quad_azimu_b_all_,quad_polar_a_all_);
pm_n_lm = pm_n_lm_(1+0);
Ylm_yq__ = zeros(pm_n_lm,quad_n_all);
nml=0;
for l_val=0:pm_l_max;
for m_val=-l_val:+l_val;
Ylm_yq__(1+nml,:) = Ylm__{1+l_val}(1+l_val+m_val,:);
nml=nml+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:pm_l_max;
Ylm_w_yq__ = Ylm_yq__ * sparse(1:quad_n_all,1:quad_n_all,quad_weight_all_,quad_n_all,quad_n_all);
pm_n_w = pm_n_w_(1+0);
[ ...
 data_k_p_polar_a__ ...
,data_k_p_azimu_b__ ...
,data_k_c_0__ ...
,data_k_c_1__ ...
,data_k_c_2__ ...
] = ...
cg_rhs_1( ...
 n_M ...
,pm_n_w ...
,euler_polar_a_ ...
,euler_azimu_b_ ...
,+euler_gamma_z_ ...
);
data_k_c_wMd__ = [ data_k_c_0__(:) , data_k_c_1__(:) , data_k_c_2__(:) ];
index_quad_from_data_ = knnsearch(quad_k_c_qd__,data_k_c_wMd__,'K',1); index_quad_from_data_ = index_quad_from_data_ - 1;
quad_from_data_qwM__ = sparse(1+index_quad_from_data_,1:pm_n_w*n_M,1,quad_n_all,pm_n_w*n_M);
n_quad_from_data_q_ = quad_from_data_qwM__*ones(pm_n_w*n_M,1); %<-- number of data-points per quadrature-point. ;
data_from_quad_wMq__ = bsxfun(@rdivide,transpose(quad_from_data_qwM__),max(1,transpose(n_quad_from_data_q_)));
VSCTF_wMc__ = reshape(repmat(reshape(VSCTF_Mc__,[1,n_M,n_CTF_rank]),[pm_n_w,1,1]),[pm_n_w*n_M,n_CTF_rank]);
%%%%;
if (flag_qbp_pm_solve_in_CTF_sequence==1);
VSCTF2_qc__ = quad_from_data_qwM__*abs(VSCTF_wMc__).^2;
end;%if (flag_qbp_pm_solve_in_CTF_sequence==1);
%%%%;
if (flag_qbp_pm_solve_in_CTF_sequence==0);
tmp_t=tic();
VSCTF2_qcc___ = zeros(quad_n_all,n_CTF_rank,n_CTF_rank);
for nCTF_rank_0=0:n_CTF_rank-1;
for nCTF_rank_1=nCTF_rank_0:n_CTF_rank-1;
VSCTF2_qcc___(:,1+nCTF_rank_0,1+nCTF_rank_1) = quad_from_data_qwM__*(VSCTF_wMc__(:,1+nCTF_rank_0).*VSCTF_wMc__(:,1+nCTF_rank_1));
VSCTF2_qcc___(:,1+nCTF_rank_1,1+nCTF_rank_0) = VSCTF2_qcc___(:,1+nCTF_rank_0,1+nCTF_rank_1);
end;%for nCTF_rank_1=nCTF_rank_0:n_CTF_rank-1;
end;%for nCTF_rank_0=0:n_CTF_rank-1;
VSCTF2_ccq___ = permute(VSCTF2_qcc___,[2,3,1]);
VSCTF2_pinv_ccq___ = zeros(n_CTF_rank,n_CTF_rank,quad_n_all);
for quad_nall=0:quad_n_all-1;
VSCTF2_cc__ = VSCTF2_ccq___(:,:,1+quad_nall);
VSCTF2_fnorm = fnorm(VSCTF2_cc__);
VSCTF2_pinv__ = pinv(VSCTF2_cc__,VSCTF2_fnorm*qbp_eps);
VSCTF2_pinv_ccq___(:,:,1+quad_nall) = VSCTF2_pinv__;
end;%for quad_nall=0:quad_n_all-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% VSCTF2: %0.2fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'VSCTF2',tmp_t);
end;%if (flag_qbp_pm_solve_in_CTF_sequence==0);
%%%%;
end;%if (numel(unique(pm_l_max_))==1 & numel(unique(pm_n_lm_))==1 & numel(unique(pm_n_w_))==1);

if (~flag_unique_pm_n);
error(sprintf(' %% Error, set all values of pm_l_max, pm_n_lm and pm_n_w to be the same in cg_lsq_pm_1.'));
end;%if (~flag_unique_pm_n);

a_UCTF_UX_Y_ync__ = zeros(pm_n_lm_sum,n_CTF_rank);
for pm_nk_p_r=0:pm_n_k_p_r-1;
%%%%%%%%;
if (flag_qbp_pm_solve_in_CTF_sequence==1);
index_Y_ = pm_n_lm_csum_(1+pm_nk_p_r) + (0:pm_n_lm-1);
index_pm_nw_ = pm_n_w_csum_(1+pm_nk_p_r) + (0:pm_n_w-1);
UX_M_I_wM__ = bsxfun(@times,UX_M_k_p_wnM__(1+index_pm_nw_,:),transpose(image_I_value_));
for nCTF_rank=0:n_CTF_rank-1;
UX_M_I_VSCTF_wM__ = bsxfun(@times,UX_M_I_wM__,transpose(VSCTF_Mc__(:,1+nCTF_rank)));
quad_from_data_UX_M_I_VSCTF_normalized_q_ = (quad_from_data_qwM__ * UX_M_I_VSCTF_wM__(:))./max(1e-12,VSCTF2_qc__(:,1+nCTF_rank));
a_UCTF_UX_Y_ync__(1+index_Y_,1+nCTF_rank) = conj(Ylm_w_yq__)*quad_from_data_UX_M_I_VSCTF_normalized_q_;
UX_M_I_wM__ = UX_M_I_wM__ - bsxfun(@times,reshape(data_from_quad_wMq__*(transpose(Ylm_yq__)*a_UCTF_UX_Y_ync__(1+index_Y_,1+nCTF_rank)),[pm_n_w,n_M]),transpose(VSCTF_Mc__(:,1+nCTF_rank))); %<-- replace with approximate residual. ;
end;%for nCTF_rank=0:n_CTF_rank-1;
end;%if (flag_qbp_pm_solve_in_CTF_sequence==1);
%%%%%%%%;
if (flag_qbp_pm_solve_in_CTF_sequence==0);
index_Y_ = pm_n_lm_csum_(1+pm_nk_p_r) + (0:pm_n_lm-1);
index_pm_nw_ = pm_n_w_csum_(1+pm_nk_p_r) + (0:pm_n_w-1);
UX_M_I_wM__ = bsxfun(@times,UX_M_k_p_wnM__(1+index_pm_nw_,:),transpose(image_I_value_));
quad_from_data_UX_M_I_VSCTF_qc__ = zeros(quad_n_all,n_CTF_rank);
for nCTF_rank=0:n_CTF_rank-1;
tmp_UX_M_I_VSCTF_wM__ = bsxfun(@times,UX_M_I_wM__,transpose(VSCTF_Mc__(:,1+nCTF_rank)));
quad_from_data_UX_M_I_VSCTF_qc__(:,1+nCTF_rank) = quad_from_data_qwM__ * tmp_UX_M_I_VSCTF_wM__(:);
end;%for nCTF_rank=0:n_CTF_rank-1;
quad_from_data_UX_M_I_VSCTF_cq__ = permute(quad_from_data_UX_M_I_VSCTF_qc__,[2,1]);
pinv_quad_from_data_UX_M_I_VSCTF_cq__ = zeros(n_CTF_rank,quad_n_all);
for quad_nall=0:quad_n_all-1;
pinv_quad_from_data_UX_M_I_VSCTF_cq__(:,1+quad_nall) = VSCTF2_pinv_ccq___(:,:,1+quad_nall)*quad_from_data_UX_M_I_VSCTF_cq__(:,1+quad_nall);
end;%for quad_nall=0:quad_n_all-1;
pinv_quad_from_data_UX_M_I_VSCTF_qc__ = permute(pinv_quad_from_data_UX_M_I_VSCTF_cq__,[2,1]);
a_UCTF_UX_Y_ync__(1+index_Y_,:) = conj(Ylm_w_yq__)*pinv_quad_from_data_UX_M_I_VSCTF_qc__;
end;%if (flag_qbp_pm_solve_in_CTF_sequence==0);
%%%%%%%%;
end;%for pm_nk_p_r=0:pm_n_k_p_r-1;
  
if (verbose>0); disp(sprintf(' %% [finished qbp_pm_1]')); end;  
