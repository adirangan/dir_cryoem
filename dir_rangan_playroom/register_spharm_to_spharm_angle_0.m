function [X0_,n_op,n_mult] = register_spharm_to_spharm_angle_0(verbose,n_k,k_,n_l_,a_,b_,n_beta,beta_,n_alpha,alpha_,n_gamma,gamma_);
% tests registration between molecule_A and molecule_B using an array of beta (fast only);
% ;
% verbose = integer verbosity_level ;
% n_k = integer maximum k ;
% k_ = real array of length n_k; k_(nk) = k_value for shell nk ;
% n_l_ = integer array of length n_k; n_l_(nk) = spherical harmonic order on shell nk; n_l_(nk) corresponds to n_lm_(nk) = (n_l_(nk)+1)^2 coefficients ;
% a_ = complex array of length \sum_{nk} (n_l_(nk)+1)^2 ; coefficients are ordered in a row, with m varying quickly and l varying slowly ;
% b_ = complex array of length \sum_{nk} (n_l_(nk)+1)^2 ; coefficients are ordered in a row, with m varying quickly and l varying slowly ;
% a_ corresponds to molecule_a, b_ to molecule_b ;
% n_beta = integer number of beta angles ;
% beta_ = real array of beta angles ;
% n_alpha = integer number of alpha angles (optional);
% alpha_ = real array of alpha angles (optional);
% n_gamma = integer number of gamma angles (optional);
% gamma_ = real array of gamma angles (optional);
% ;
% If arrays alpha_ and gamma_ are not provided we use the standard arrays: ;
% alpha_ = linspace(0,2*pi,n_m_max+1); alpha_ = alpha_(1:end-1);
% gamma_ = linspace(0,2*pi,n_m_max+1); gamma_ = gamma_(1:end-1);
% and use the standard fft to calculate X0_.; 
% ;
% However, if n_alpha, alpha_ and n_gamma, gamma_ are provided, we use these arrays instead, ;
% and use the nufft to calculate X0_. ;
% ;
% X_ = complex array of size (n_alpha,n_gamma,n_beta);
% The default values of n_alpha and n_gamma are n_m_max. ;
% X_(nalpha,ngamma,nbeta) corresponds to the innerproduct between molecule_A and molecule_B, where ;
% the latter has been rotated by euler-angles alpha,beta,gamma. ;
% Note that alpha_ and gamma_ are arrays from 0 to 2*pi, ;
% whereas beta_ is an array from -pi to pi. ;

if nargin<1;
verbose=0;
n_k=1;k_=1;n_l_=128;n_lm_=(n_l_+1).^2;
a_=randn(sum(n_lm_),1);b_=randn(sum(n_lm_),1);
n_beta=128;beta_=linspace(-pi,+pi,n_beta+1);beta_=beta_(1:end-1);
tic;
[X0_,n_op,n_mult] = register_spharm_to_spharm_angle_0(verbose,n_k,k_,n_l_,a_,b_,n_beta,beta_);
toc;
disp('returning');return;
end;%if nargin<1;

n_lm_ = (n_l_+1).^2; n_lm_csum_ = cumsum([0;n_lm_(:)]);
k_max = k_(end);
n_l_max = n_l_(end);
m_max_ = -n_l_max : +n_l_max;
n_m_max = length(m_max_);

if (nargin<=8);
fft_flag = 0; 
n_alpha = n_m_max;
n_gamma = n_m_max;
end;%if (nargin<=8);
if (nargin>=9 & nargin<=12); 
fft_flag = 1;
end;%if (nargin>=9 & nargin<=12); 

X0_ = zeros(n_alpha,n_gamma,n_beta);
verbose_tab = 0;
for nbeta = 0:n_beta-1;
beta = beta_(1+nbeta);
% Note that rotating a molecule by [+alpha,+beta,+gamma] ;
% corresponds to rotating the coordinate-frame by [-gamma,-beta,-alpha] ;
W_ = wignerd_b(n_l_max,-beta);
C_ = zeros(n_m_max,n_m_max);
n_W_ = zeros(1,1+n_l_max); for (nl=0:n_l_max); n_W_(1+nl) = numel(W_{1+nl}); end;
n_wignerd_ = cumsum([0,n_W_]);
wignerd_ = zeros(1,n_wignerd_(end));
for (nl=0:n_l_max); wignerd_(n_wignerd_(1+nl) + (1:numel(W_{1+nl}))) = W_{1+nl}(:); end;
tic; n_mult=0; n_op = 0;
for nmn = 0:n_m_max-1;
mn = -n_l_max + nmn;
for nmp = 0:n_m_max-1;
mp = -n_l_max + nmp;
C = 0;
a_base = 0;
b_base = 0;
for nk = 0:n_k-1;
k_val = k_(1+nk); k2_val = k_val.^2;
n_l = n_l_(1+nk); n_lm = (1+n_l)*(1+n_l);
% X0: ;
%ix_base = sum(n_lm_(1:1+nk-1)); a_k_ = a_(ix_base + (1:n_lm)); b_k_ = b_(ix_base + (1:n_lm));
ix = n_lm_csum_(1+nk); 
if (abs(mn)<=n_l & abs(mp)<=n_l);
%for nl = 0:n_l;
for nl = max(abs(mn),abs(mp)):n_l;
nl2 = nl*(nl+1); ix_base = ix + nl2; ix_mn = ix_base + mn; ix_mp = ix_base + mp;
n_wignerd = n_wignerd_(1+nl);
mn_flag=0; if (abs(mn)<=nl); mn_flag=1; a_mn = nl2 + mn; end;
mp_flag=0; if (abs(mp)<=nl); mp_flag=1; b_mp = nl2 + mp; end;
if (mn_flag & mp_flag);
nwignerd = n_wignerd + (nl+mn) + (nl+mp)*(1+2*nl);
wignerd = wignerd_(1+nwignerd);
a = a_(1+a_base+a_mn);
b = b_(1+b_base+b_mp);
%C = C + k2_val * conj(a)*wignerd*b;
% X0: ;
%C = C + k2_val * conj(a_k_(1+a_mn))*W_{1+nl}(1+nl+mn,1+nl+mp)*b_k_(1+b_mp);
%C = C + k2_val * conj(a)*W_{1+nl}(1+nl+mn,1+nl+mp)*b;
C = C + k2_val * conj(a)*wignerd*b;
n_mult = n_mult+3; n_op = n_op + 1;
if (verbose>2); disp(sprintf(' ')); end; %verbose
if (verbose>2); disp(sprintf(' nmn %d mn %d ',nmn-1,mn)); end; %verbose
if (verbose>2); disp(sprintf(' nmp %d mp %d ',nmp-1,mp)); end; %verbose
if (verbose>2); disp(sprintf(' nk %d k %f ',nk-1,k_val)); end; %verbose
if (verbose>2); disp(sprintf(' n_l %d n_lm %d ',n_l,n_lm)); end; %verbose
if (verbose>2); disp(sprintf(' a_base %d b_base %d ',a_base,b_base)); end; %verbose
if (verbose>2); disp(sprintf(' n_wignerd %d a_mn %d b_mp %d',n_wignerd,ix_mn-1,ix_mp-1)); end; %verbose
if (verbose>2); disp(sprintf(' nwignerd %d wignerd %f ',nwignerd,wignerd)); end; %verbose
if (verbose>2); disp(sprintf(' a %f,%f b %f,%f ',real(a_(ix_mn)),imag(a_(ix_mn)),real(b_(ix_mp)),imag(b_(ix_mp)))); end; %verbose
if (verbose>2); disp(sprintf(' C %f,%f ',real(C),imag(C))); end; %verbose
if (verbose>2); verbose_tab = verbose_tab+1; end; %verbose
if (verbose_tab>5); disp('returning'); return; end;
end;%if (mn_flag & mp_flag);
end;%for nl = 0:n_l;
end;%if (abs(mn)<=n_l & abs(mp)<=n_l);
a_base = a_base + n_lm; b_base = b_base + n_lm;
end;%for nk = 0:n_k-1;
C_(1+nmn,1+nmp) = C;
end;%for nmp = 0:n_m_max-1;
end;%for nmn = 0:n_m_max-1;
t_0 = toc; t_0_(1+nbeta) = t_0;
if (verbose); disp(sprintf(' %% t_0_(%d) %0.4f, n_mult/s %f n_op/s %f',1+nbeta,t_0,n_mult/t_0,n_op/t_0)); end;
if (fft_flag==0);
tmp_C = recenter2(squeeze(C_(:,:)));
X0_(:,:,1+nbeta) = fft2(tmp_C);
end;%if (fft_flag==0);
if (fft_flag==1);
[gamma__,alpha__] = meshgrid(gamma_,alpha_);
tmp_C = decenter2(recenter2(squeeze(C_(:,:))));
tmp = nufft2d2(n_alpha*n_gamma,alpha__(:),gamma__(:),-1,1e-12,n_m_max,n_m_max,tmp_C);
X0_(:,:,1+nbeta) = reshape(tmp,n_alpha,n_gamma);
end;%if (fft_flag==1);
end;%for nbeta = 0:n_beta-1;
