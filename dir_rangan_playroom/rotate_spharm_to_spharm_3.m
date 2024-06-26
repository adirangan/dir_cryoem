function ...
[ ...
 b_ ...
,d0W__ ...
,V_lmm___ ...
,L_lm__ ...
] = ...
rotate_spharm_to_spharm_3( ...
 n_k ...
,k_ ...
,l_max_ ...
,a_ ...
,euler_ ...
,d0W__ ...
,V_lmm___ ...
,L_lm__ ...
);
% Rotates spherical harmonic expansion a_ by euler_, producing b_ ;
% Note that there needs to be a sign-flip for the interior angle. ;
% ;
% verbose = integer verbosity_level ;
% W_beta_0in__ = wignerd values. If empty these will be computed. ;
% n_k = integer maximum k ;
% k_ = real array of length n_k; k_(nk) = k_value for shell nk ;
% l_max_ = integer array of length n_k; l_max_(nk) = spherical harmonic order on shell nk; l_max_(nk) corresponds to n_lm_(nk) = (l_max_(nk)+1)^2 coefficients ;
% a_ = complex array of length \sum_{nk} (n_lm_(nk)+1)^2 ; coefficients are ordered in a row, with m varying quickly and l varying slowly ;
% euler_ equals array of euler angles; [alpha, beta, gamma]. Note that this corresponds to rotating k_ by the inverse [-gamma, -beta, -alpha]. ;
% ;
% b_ = complex array of length \sum_{nk} (n_lm_(nk)+1)^2 ; coefficients are ordered in a row, with m varying quickly and l varying slowly ;
% b_ corresponds to rotated molecule ;
% d0W__ = wignerd values (computed using wignerd_c). ;
% V_lmm___ = precomputation for wignerd_c. ;
% L_lm__ = precomputation for wignerd_c. ;

str_thisfunction = 'rotate_spharm_to_spharm_3';

if (nargin<1);
%%%%%%%%;
flag_verbose = 1; rng(0);
if (flag_verbose); disp(sprintf(' %% testing rotate_spharm_to_spharm_3')); end;
%%%%%%%%;
k_p_r_max = 48/(2*pi); k_eq_d = 1.0/(2*pi);
[ ...
 n_k_all ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,weight_shell_k_ ...
,k_c_0_all_ ...
,k_c_1_all_ ...
,k_c_2_all_ ...
] = ...
sample_shell_5( ...
 k_p_r_max ...
,k_eq_d ...
,'L' ...
) ;
%%%%%%%%;
n_k_p_r = 1;
n_k_all_csum_ = [0;n_k_all];
k_p_r_all_ = k_p_r_max*ones(n_k_all,1);
k_p_r_ = k_p_r_max;
weight_3d_k_p_r_ = 1;
weight_3d_k_all_ = weight_shell_k_;
%%%%%%%%;
l_max_upb = 36;
l_max_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
l_max_(1+nk_p_r) = max(0,min(l_max_upb,1+ceil(2*pi*k_p_r_(1+nk_p_r))));
end;%for nk_p_r=0:n_k_p_r-1;
n_lm_ = (l_max_+1).^2;
n_lm_max = max(n_lm_);
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
l_max_max = max(l_max_);
m_max_ = -l_max_max : +l_max_max;
n_m_max = length(m_max_);
Y_l_val_ = zeros(n_lm_sum,1);
Y_m_val_ = zeros(n_lm_sum,1);
Y_k_val_ = zeros(n_lm_sum,1);
for nk_p_r=0:n_k_p_r-1;
l_max = l_max_(1+nk_p_r);
tmp_l_val_ = zeros(n_lm_(1+nk_p_r),1);
tmp_m_val_ = zeros(n_lm_(1+nk_p_r),1);
na=0; 
for l_val=0:l_max;
for m_val=-l_val:+l_val;
tmp_l_val_(1+na) = l_val;
tmp_m_val_(1+na) = m_val;
na=na+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
Y_l_val_(1+tmp_index_) = tmp_l_val_;
Y_m_val_(1+tmp_index_) = tmp_m_val_;
Y_k_val_(1+tmp_index_) = k_p_r_(1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
weight_Y_ = zeros(n_lm_sum,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
weight_Y_(1+tmp_index_) = weight_3d_k_p_r_(1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
if (flag_verbose); disp(sprintf(' %% n_k_all %d, l_max_max %d',n_k_all,l_max_max)); end;
%%%%%%%%;
a_k_Y_true_ = ( randn(n_lm_sum,1) + i*randn(n_lm_sum,1) ).*exp(-(Y_l_val_+abs(Y_m_val_)).^2/(2*(l_max_max/4)^2));
tmp_t = tic;
[a_k_p_quad_] = real(convert_spharm_to_k_p_1(flag_verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_Y_true_));
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_true_ --> a_k_p_quad_ time %0.2fs',tmp_t));
tmp_t = tic;
[a_k_Y_reco_] = convert_k_p_to_spharm_1(flag_verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_p_quad_);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_reco_ time %0.2fs',tmp_t));
a_k_Y_true_ = a_k_Y_reco_;
%disp(sprintf(' %% a_k_Y_reco error: %0.16f',fnorm(a_k_Y_true_-a_k_Y_reco_)/fnorm(a_k_Y_true_)));
%%%%%%%%;
euler_d_ = [+0.4;-0.25;-0.1];
%%%%%%%%;
euler_b_ = [euler_d_(1+0);0;0];
b_k_Y_true_ = rotate_spharm_to_spharm_3(n_k_p_r,k_p_r_,l_max_,a_k_Y_true_,euler_b_*pi);
tmp_t = tic;
[b_k_p_quad_] = real(convert_spharm_to_k_p_1(flag_verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,b_k_Y_true_));
tmp_t = toc(tmp_t); disp(sprintf(' %% b_k_Y_true_ --> b_k_p_quad_ time %0.2fs',tmp_t));
%%%%%%%%;
euler_c_ = [euler_d_(1+0);euler_d_(1+1);0];
c_k_Y_true_ = rotate_spharm_to_spharm_3(n_k_p_r,k_p_r_,l_max_,a_k_Y_true_,euler_c_*pi);
tmp_t = tic;
[c_k_p_quad_] = real(convert_spharm_to_k_p_1(flag_verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,c_k_Y_true_));
tmp_t = toc(tmp_t); disp(sprintf(' %% c_k_Y_true_ --> c_k_p_quad_ time %0.2fs',tmp_t));
%%%%%%%%;
d_k_Y_true_ = rotate_spharm_to_spharm_3(n_k_p_r,k_p_r_,l_max_,a_k_Y_true_,euler_d_*pi);
tmp_t = tic;
[d_k_p_quad_] = real(convert_spharm_to_k_p_1(flag_verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,d_k_Y_true_));
tmp_t = toc(tmp_t); disp(sprintf(' %% d_k_Y_true_ --> d_k_p_quad_ time %0.2fs',tmp_t));
%%%%%%%%;
e_k_Y_true_ = a_k_Y_true_;
e_k_Y_true_ = rotate_spharm_to_spharm_3(n_k_p_r,k_p_r_,l_max_,e_k_Y_true_,[euler_d_(1+0);0;0]*pi);
e_k_Y_true_ = rotate_spharm_to_spharm_3(n_k_p_r,k_p_r_,l_max_,e_k_Y_true_,[0;euler_d_(1+1);0]*pi);
e_k_Y_true_ = rotate_spharm_to_spharm_3(n_k_p_r,k_p_r_,l_max_,e_k_Y_true_,[0;0;euler_d_(1+2)]*pi);
disp(sprintf(' %% e_k_Y_true vs d_k_Y_true error: %0.16f',fnorm(e_k_Y_true_-d_k_Y_true_)/fnorm(e_k_Y_true_)));
%%%%%%%%;
figure(1);clf;
flag_2d_vs_3d = 1;
subplot(2,2,1);
imagesc_polar_a_azimu_b_0(k_p_polar_a_all_,k_p_azimu_b_all_,real(a_k_p_quad_),4*[-1,+1],colormap(colormap_80s));
xlim([0,2*pi]); ylim([0,1*pi]);
set(gca,'XTick',linspace(0,2*pi,21),'XTickLabel',linspace(0,2,21)); xlabel('azimu/\pi'); xtickangle(90);
set(gca,'YTick',linspace(0,1*pi,11),'YTickLabel',linspace(0,1,11)); ylabel('polar/\pi');
title(sprintf('\\tau [0;0;0]\\pi'));
subplot(2,2,2);
imagesc_polar_a_azimu_b_0(k_p_polar_a_all_,k_p_azimu_b_all_,real(b_k_p_quad_),4*[-1,+1],colormap(colormap_80s));
xlim([0,2*pi]); ylim([0,1*pi]);
set(gca,'XTick',linspace(0,2*pi,21),'XTickLabel',linspace(0,2,21)); xlabel('azimu/\pi'); xtickangle(90);
set(gca,'YTick',linspace(0,1*pi,11),'YTickLabel',linspace(0,1,11)); ylabel('polar/\pi');
title(sprintf('\\tau [%+0.2f;%+0.2f;%+0.2f]\\pi',euler_b_));
subplot(2,2,3);
imagesc_polar_a_azimu_b_0(k_p_polar_a_all_,k_p_azimu_b_all_,real(c_k_p_quad_),4*[-1,+1],colormap(colormap_80s));
xlim([0,2*pi]); ylim([0,1*pi]);
set(gca,'XTick',linspace(0,2*pi,21),'XTickLabel',linspace(0,2,21)); xlabel('azimu/\pi'); xtickangle(90);
set(gca,'YTick',linspace(0,1*pi,11),'YTickLabel',linspace(0,1,11)); ylabel('polar/\pi');
title(sprintf('\\tau [%+0.2f;%+0.2f;%+0.2f]\\pi',euler_c_));
subplot(2,2,4);
imagesc_polar_a_azimu_b_0(k_p_polar_a_all_,k_p_azimu_b_all_,real(d_k_p_quad_),4*[-1,+1],colormap(colormap_80s));
xlim([0,2*pi]); ylim([0,1*pi]);
set(gca,'XTick',linspace(0,2*pi,21),'XTickLabel',linspace(0,2,21)); xlabel('azimu/\pi'); xtickangle(90);
set(gca,'YTick',linspace(0,1*pi,11),'YTickLabel',linspace(0,1,11)); ylabel('polar/\pi');
title(sprintf('\\tau [%+0.2f;%+0.2f;%+0.2f]\\pi',euler_d_));
figbig;
%%%%%%%%;
figure(2);clf;
flag_2d_vs_3d = 0;
subplot(2,2,1);
imagesc_polar_a_azimu_b_0(k_p_polar_a_all_,k_p_azimu_b_all_,real(a_k_p_quad_),4*[-1,+1],colormap(colormap_80s),flag_2d_vs_3d);
xlabel('x');ylabel('y');zlabel('z'); axis vis3d;
title(sprintf('\\tau [0;0;0]\\pi'));
subplot(2,2,2);
imagesc_polar_a_azimu_b_0(k_p_polar_a_all_,k_p_azimu_b_all_,real(b_k_p_quad_),4*[-1,+1],colormap(colormap_80s),flag_2d_vs_3d);
xlabel('x');ylabel('y');zlabel('z'); axis vis3d;
title(sprintf('\\tau [%+0.2f;%+0.2f;%+0.2f]\\pi',euler_b_));
subplot(2,2,3);
imagesc_polar_a_azimu_b_0(k_p_polar_a_all_,k_p_azimu_b_all_,real(c_k_p_quad_),4*[-1,+1],colormap(colormap_80s),flag_2d_vs_3d);
xlabel('x');ylabel('y');zlabel('z'); axis vis3d;
title(sprintf('\\tau [%+0.2f;%+0.2f;%+0.2f]\\pi',euler_c_));
subplot(2,2,4);
imagesc_polar_a_azimu_b_0(k_p_polar_a_all_,k_p_azimu_b_all_,real(d_k_p_quad_),4*[-1,+1],colormap(colormap_80s),flag_2d_vs_3d);
xlabel('x');ylabel('y');zlabel('z'); axis vis3d;
title(sprintf('\\tau [%+0.2f;%+0.2f;%+0.2f]\\pi',euler_d_));
figbig;
%%%%%%%%;
figure(3);clf;
flag_2d_vs_3d = 0;
subplot(2,3,1);
imagesc_polar_a_azimu_b_0(k_p_polar_a_all_,k_p_azimu_b_all_,k_c_0_all_,[-1,+1],colormap(colormap_80s),flag_2d_vs_3d);
xlabel('x');ylabel('y');zlabel('z'); axis vis3d;
title(sprintf('x'));
subplot(2,3,2);
imagesc_polar_a_azimu_b_0(k_p_polar_a_all_,k_p_azimu_b_all_,k_c_1_all_,[-1,+1],colormap(colormap_80s),flag_2d_vs_3d);
xlabel('x');ylabel('y');zlabel('z'); axis vis3d;
title(sprintf('y',euler_b_));
subplot(2,3,3);
imagesc_polar_a_azimu_b_0(k_p_polar_a_all_,k_p_azimu_b_all_,k_c_2_all_,[-1,+1],colormap(colormap_80s),flag_2d_vs_3d);
xlabel('x');ylabel('y');zlabel('z'); axis vis3d;
title(sprintf('z',euler_c_));
subplot(2,3,4);
imagesc_polar_a_azimu_b_0(k_p_polar_a_all_,k_p_azimu_b_all_,k_p_polar_a_all_,[0,1*pi],colormap(colormap_80s),flag_2d_vs_3d);
xlabel('x');ylabel('y');zlabel('z'); axis vis3d;
title(sprintf('polar',euler_c_));
subplot(2,3,5);
imagesc_polar_a_azimu_b_0(k_p_polar_a_all_,k_p_azimu_b_all_,k_p_azimu_b_all_,[0,2*pi],colormap(colormap_80s),flag_2d_vs_3d);
xlabel('x');ylabel('y');zlabel('z'); axis vis3d;
title(sprintf('azimu',euler_c_));
figbig;
%%%%%%%%;
disp('returning'); return;
end;%if (nargin<1);

flag_verbose=0;

na=0;
if (nargin<1+na); n_k=[]; end; na=na+1;
if (nargin<1+na); k_=[]; end; na=na+1;
if (nargin<1+na); l_max_=[]; end; na=na+1;
if (nargin<1+na); a_=[]; end; na=na+1;
if (nargin<1+na); euler_=[]; end; na=na+1;
if (nargin<1+na); d0W__=[]; end; na=na+1;
if (nargin<1+na); V_lmm___=[]; end; na=na+1;
if (nargin<1+na); L_lm__=[]; end; na=na+1;

if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

n_lm_ = (l_max_+1).^2;
k_max = k_(end);
l_max = l_max_(end);
m_max_ = -l_max : +l_max;
n_m_max = length(m_max_);

if (flag_verbose>1); disp(sprintf(' %% rotating molecule by: [%0.2f %0.2f %0.2f]',+euler_)); end;
if (flag_verbose>1); disp(sprintf(' %% rotating coordinate_frame by: [%0.2f %0.2f %0.2f]',-euler_(3),-euler_(2),-euler_(1))); end;
b_ = zeros(size(a_));
l_max_max = max(l_max_);
W_alpha_ = exp(+i*[-l_max_max:+l_max_max]*-euler_(3));
W_gamma_ = exp(+i*[-l_max_max:+l_max_max]*-euler_(1));
if ( isempty(d0W__)); d0W__ = wignerd_c(l_max_max,+euler_(2),V_lmm___,L_lm__); end;
for nk=1:n_k;
l_max = l_max_(nk); n_lm = n_lm_(nk); ix_base = sum(n_lm_(1:nk-1));
a_k_ = a_(ix_base + (1:n_lm)); b_k_ = zeros(size(a_k_));
for l_val=0:l_max;
tmp_ij_ = 1 + [0:1+2*l_val-1] + l_max_max-l_val;
W_alpha = diag(W_alpha_(tmp_ij_));
W_gamma = diag(W_gamma_(tmp_ij_));
a_k_tmp = a_k_(1+l_val*(l_val+1) + (-l_val:+l_val));
a_k_tmp = reshape(a_k_tmp,2*l_val+1,1);
b_k_tmp = W_alpha*d0W__{1+l_val}*W_gamma*a_k_tmp;
b_k_(1+l_val*(l_val+1) + (-l_val:+l_val)) = b_k_tmp;
end;%for l_val=0:l_max;
b_(ix_base + (1:n_lm)) = b_k_;
end;%for nk=1:n_k;

if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
