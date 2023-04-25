function ...
[ ...
 b_k_p_ba_ ...
,shell_scatter_from_tensor_sba__ ...
] = ...
shell_rotate_k_p_to_k_p_0( ...
 verbose ...
,n_azimu_b ...
,n_polar_a ...
,polar_a_ ...
,a_k_p_ba_ ...
,euler_ ...
,flag_flip ...
,n_order ...
);

str_thisfunction = 'shell_rotate_k_p_to_k_p_0';

if (nargin<1);
%%%%%%%%;
verbose = 1; rng(0);
if (verbose); disp(sprintf(' %% testing %s',str_thisfunction)); end;
%%%%%%%%;
k_p_r_max = 48/(2*pi); k_eq_d = 1.0/(2*pi);
flag_tensor_vs_adap = 1;
[ ...
 n_k_all ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,weight_shell_k_ ...
,k_c_0_all_ ...
,k_c_1_all_ ...
,k_c_2_all_ ...
,n_polar_a ...
,polar_a_ ...
,n_azimu_b_ ...
] = ...
sample_shell_6( ...
 k_p_r_max ...
,k_eq_d ...
,'L' ...
,flag_tensor_vs_adap ...
) ;
n_azimu_b = n_azimu_b_(1+0);
azimu_b_ = linspace(0,2*pi,1+n_azimu_b); azimu_b_ = transpose(azimu_b_(1:end-1));
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
if (verbose); disp(sprintf(' %% n_k_all %d, l_max_max %d',n_k_all,l_max_max)); end;
%%%%%%%%;
a_k_Y_true_ = ( randn(n_lm_sum,1) + i*randn(n_lm_sum,1) ).*exp(-(Y_l_val_+abs(Y_m_val_)).^2/(2*(l_max_max/4)^2));
tmp_t = tic;
[a_k_p_quad_] = real(convert_spharm_to_k_p_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_Y_true_));
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_true_ --> a_k_p_quad_ time %0.2fs',tmp_t));
tmp_t = tic;
[a_k_Y_reco_] = convert_k_p_to_spharm_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_p_quad_);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_reco_ time %0.2fs',tmp_t));
a_k_Y_true_ = a_k_Y_reco_;
disp(sprintf(' %% a_k_Y_reco error: %0.16f',fnorm(a_k_Y_true_-a_k_Y_reco_)/fnorm(a_k_Y_true_)));
%%%%%%%%;
euler_d_ = [+0.4;-0.25;-0.1];
%%%%%%%%;
euler_b_ = [euler_d_(1+0);0;0];
b_k_p_true_ = shell_rotate_k_p_to_k_p_0(1+verbose,n_azimu_b,n_polar_a,polar_a_,a_k_p_quad_,euler_b_*pi);
b_k_Y_true_ = rotate_spharm_to_spharm_2(1+verbose,[],n_k_p_r,k_p_r_,l_max_,a_k_Y_true_,euler_b_*pi);
tmp_t = tic;
[b_k_p_quad_] = real(convert_spharm_to_k_p_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,b_k_Y_true_));
tmp_t = toc(tmp_t); disp(sprintf(' %% b_k_Y_true_ --> b_k_p_quad_ time %0.2fs',tmp_t));
disp(sprintf(' %% b_k_p_true_ vs b_k_p_quad_: %0.16f',fnorm(b_k_p_true_ - b_k_p_quad_)/fnorm(b_k_p_true_))); 
%{
figure(1);clf;figmed;
subplot(1,3,1);
imagesc_polar_a_azimu_b_0( ...
 k_p_polar_a_all_ ... 
,k_p_azimu_b_all_ ... 
,real(a_k_p_quad_) ... 
);
title('a_k_p_quad_','Interpreter','none');
subplot(1,3,2);
imagesc_polar_a_azimu_b_0( ...
 k_p_polar_a_all_ ... 
,k_p_azimu_b_all_ ... 
,real(b_k_p_true_) ... 
);
title('b_k_p_true_','Interpreter','none');
subplot(1,3,3);
imagesc_polar_a_azimu_b_0( ...
 k_p_polar_a_all_ ... 
,k_p_azimu_b_all_ ... 
,real(b_k_p_quad_) ... 
);
title('b_k_p_quad_','Interpreter','none');
error('stop');
%}
%%%%%%%%;
euler_c_ = [euler_d_(1+0);euler_d_(1+1);0];
c_k_p_true_ = shell_rotate_k_p_to_k_p_0(1+verbose,n_azimu_b,n_polar_a,polar_a_,a_k_p_quad_,euler_c_*pi);
c_k_Y_true_ = rotate_spharm_to_spharm_2(1+verbose,[],n_k_p_r,k_p_r_,l_max_,a_k_Y_true_,euler_c_*pi);
tmp_t = tic;
[c_k_p_quad_] = real(convert_spharm_to_k_p_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,c_k_Y_true_));
tmp_t = toc(tmp_t); disp(sprintf(' %% c_k_Y_true_ --> c_k_p_quad_ time %0.2fs',tmp_t));
disp(sprintf(' %% c_k_p_true_ vs c_k_p_quad_: %0.16f',fnorm(c_k_p_true_ - c_k_p_quad_)/fnorm(c_k_p_true_))); 
%%%%%%%%;
d_k_p_true_ = shell_rotate_k_p_to_k_p_0(1+verbose,n_azimu_b,n_polar_a,polar_a_,a_k_p_quad_,euler_d_*pi);
d_k_Y_true_ = rotate_spharm_to_spharm_2(1+verbose,[],n_k_p_r,k_p_r_,l_max_,a_k_Y_true_,euler_d_*pi);
tmp_t = tic;
[d_k_p_quad_] = real(convert_spharm_to_k_p_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,d_k_Y_true_));
tmp_t = toc(tmp_t); disp(sprintf(' %% d_k_Y_true_ --> d_k_p_quad_ time %0.2fs',tmp_t));
disp(sprintf(' %% d_k_p_true_ vs d_k_p_quad_: %0.16f',fnorm(d_k_p_true_ - d_k_p_quad_)/fnorm(d_k_p_true_))); 
%%%%%%%%;
e_k_Y_true_ = a_k_Y_true_;
e_k_Y_true_ = rotate_spharm_to_spharm_2(1+verbose,[],n_k_p_r,k_p_r_,l_max_,e_k_Y_true_,[euler_d_(1+0);0;0]*pi);
e_k_Y_true_ = rotate_spharm_to_spharm_2(1+verbose,[],n_k_p_r,k_p_r_,l_max_,e_k_Y_true_,[0;euler_d_(1+1);0]*pi);
e_k_Y_true_ = rotate_spharm_to_spharm_2(1+verbose,[],n_k_p_r,k_p_r_,l_max_,e_k_Y_true_,[0;0;euler_d_(1+2)]*pi);
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

if (verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

na=0;
if (nargin<1+na); verbose=[]; end; na=na+1;
if (nargin<1+na); n_azimu_b=[]; end; na=na+1;
if (nargin<1+na); n_polar_a=[]; end; na=na+1;
if (nargin<1+na); polar_a_=[]; end; na=na+1;
if (nargin<1+na); a_k_p_ba_=[]; end; na=na+1;
if (nargin<1+na); euler_=[]; end; na=na+1;
if (nargin<1+na); flag_flip=[]; end; na=na+1;
if (nargin<1+na); n_order=[]; end; na=na+1;
if isempty(verbose); verbose=0; end;
if isempty(polar_a_); polar_a_ = transpose(linspace(pi,0,n_polar_a)); end;
if isempty(euler_); euler_ = [0;0;0]; end;
if isempty(flag_flip); flag_flip = 0; end;
if isempty(n_order); n_order = 5; end;

n_ba = n_azimu_b*n_polar_a;

%%%%%%%%;
% define rotation matrices. ;
%%%%%%%%;
Rz = @(azimu_b) ...
[ +cos(azimu_b) -sin(azimu_b) 0 ; ...
  +sin(azimu_b) +cos(azimu_b) 0 ; ...
   0             0            1 ; ...
] ;
%%%%%%%%;
Ry = @(polar_a) ...
[ +cos(polar_a) 0 +sin(polar_a) ; ...
   0            1  0            ; ...
  -sin(polar_a) 0 +cos(polar_a) ; ...
];
%%%%%%%%;

%%%%%%%%;
% define scatter. ;
%%%%%%%%;
azimu_b_ = linspace(0,2*pi,n_azimu_b+1); azimu_b_ = transpose(azimu_b_(1:end-1));
[azimu_b_ba__,polar_a_ba__] = ndgrid(azimu_b_,polar_a_);
c_0_ba__ = cos(azimu_b_ba__).*sin(polar_a_ba__);
c_1_ba__ = sin(azimu_b_ba__).*sin(polar_a_ba__);
c_2_ba__ = cos(polar_a_ba__);
if flag_flip; c_0_ba__ = -c_0_ba__; c_1_ba__ = -c_1_ba__; c_2_ba__ = -c_2_ba__; end;%if flag_flip;
c_3ba__ = Rz(-euler_(1+0))*Ry(-euler_(1+1))*Rz(-euler_(1+2))*[reshape(c_0_ba__,[1,n_ba]);reshape(c_1_ba__,[1,n_ba]);reshape(c_2_ba__,[1,n_ba])];
c_0_ba__ = reshape(c_3ba__(1+0,:),[n_azimu_b,n_polar_a]);
c_1_ba__ = reshape(c_3ba__(1+1,:),[n_azimu_b,n_polar_a]);
c_2_ba__ = reshape(c_3ba__(1+2,:),[n_azimu_b,n_polar_a]);
c_r01_ba__ = sqrt(c_0_ba__.^2 + c_1_ba__.^2);
c_b_ba__ = atan2(c_1_ba__,c_0_ba__);
c_a_ba__ = atan2(c_r01_ba__,c_2_ba__);
%%%%%%%%;
% construct interpolator. ;
%%%%%%%%;
flag_polar_a_ascend_vs_descend = (polar_a_(end)>=polar_a_(1));
[ ...
 shell_scatter_from_tensor_sba__ ...
] = ...
shell_k_p_scatter_from_tensor_interpolate_n_5( ...
 n_order ...
,n_azimu_b ...
,n_polar_a ...
,polar_a_ ...
,n_ba ...
,c_b_ba__(:) ...
,c_a_ba__(:) ...
,flag_polar_a_ascend_vs_descend ...
);
%%%%%%%%;
% apply interpolator. ;
%%%%%%%%;
b_k_p_ba_ = [];
if ~isempty(a_k_p_ba_);
b_k_p_ba_ = reshape(shell_scatter_from_tensor_sba__*reshape(a_k_p_ba_,[n_ba,1]),size(a_k_p_ba_));
end;%if ~isempty(a_k_p_ba_);

if (verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
