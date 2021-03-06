function X_ = test_spharm_6(n_k,k_,n_l_,a_,b_);
% tests registration between molecule_A and molecule_B using an array of beta (fast only);
% when no inputs are passed we import two spherical harmonic representations (generated by kspacegrid_to_model): ;
% molecule_A: modsph_A_ori = spiral ;
% molecule_B: modsph_B_ori = spiral with twisted tail ;
% ;
% n_k = integer maximum k ;
% k_ = integer array of length n_k; k_(nk) = k_value for shell nk ;
% n_l_ = integer array of length n_k; n_l_(nk) = spherical harmonic order on shell nk; n_l_(nk) corresponds to n_lm_(nk) = (n_l_(nk)+1)^2 coefficients ;
% a_ = complex array of length \sum_{nk} (n_lm_(nk)+1)^2 ; coefficients are ordered in a row, with m varying quickly and l varying slowly ;
% b_ = complex array of length \sum_{nk} (n_lm_(nk)+1)^2 ; coefficients are ordered in a row, with m varying quickly and l varying slowly ;
% a_ corresponds to molecule_a, b_ to molecule_b ;
% ;
% X_ = complex array of size (n_m_max,n_m_max,n_m_max) ;
% X_(nalpha,ngamma,nbeta) corresponds to the innerproduct between molecule_A and molecule_B, where ;
% the latter has been rotated by euler-angles alpha,beta,gamma. ;
% Note that alpha_ and gamma_ are arrays from 0 to 2*pi, ;
% whereas beta_ is an array from -pi to pi. ;
% ;
% test with: ;
%{
  X_ = test_spharm_6();
  %}

verbose=1;

if nargin<5;
isph_start_ = MDA_read_i4('./dir_mda6/isph_start_.mda');
nterms_sph_ = MDA_read_i4('./dir_mda6/nterms_sph_.mda');
modsph_A_ori_ = MDA_read_c16('./dir_mda6/modsph_A_ori_.mda');
modsph_B_ori_ = MDA_read_c16('./dir_mda6/modsph_B_ori_.mda');
n_k = length(isph_start_);
k_ = 1:n_k;
n_l_ = nterms_sph_;
n_lm_ = (n_l_+1).^2;
a_ = modsph_A_ori_;
b_ = modsph_B_ori_;
end;%if nargin<4;

% generating innerproduct array over beta_;
k_max = k_(end);
n_l_max = n_l_(end);
m_max_ = -n_l_max : +n_l_max;
n_m_max = length(m_max_);

beta_ = linspace(-pi,pi,n_m_max); n_beta = length(beta_);
W_ = cell(n_m_max,1);
for nbeta = 1:n_beta;
beta = beta_(nbeta);
W_{nbeta} = wignerd_b(n_l_max,beta);
end;%for nbeta = 1:n_beta;

C_ = zeros(n_m_max,n_m_max,n_beta);
for nmn = 1:n_m_max;
mn = m_max_(nmn);
for nmp = 1:n_m_max;
mp = m_max_(nmp);
C_tmp = zeros(1,n_beta);
for nk = 1:n_k;
k_val = nk;
n_l = n_l_(nk); n_lm = n_lm_(nk); ix_base = sum(n_lm_(1:nk-1));
a_k_ = a_(ix_base + (1:n_lm)); b_k_ = b_(ix_base + (1:n_lm));
if (abs(mn)<=n_l & abs(mp)<=n_l);
for nl = 0:n_l;
mn_flag=0; if (abs(mn)<=nl); mn_flag=1; ix_mn = 1+nl*(nl+1)+mn; end;
mp_flag=0; if (abs(mp)<=nl); mp_flag=1; ix_mp = 1+nl*(nl+1)+mp; end;
if (mn_flag & mp_flag);
for nbeta = 1:n_beta;
C_tmp(nbeta) = C_tmp(nbeta) + k_val^2 * conj(a_k_(ix_mn))*W_{nbeta}{1+nl}(1+nl+mn,1+nl+mp)*b_k_(ix_mp);
end;%for nbeta = 1:n_beta;
end;%if (mn_flag & mp_flag);
end;%for nl = 0:n_l;
end;%if (abs(mn)<=n_l & abs(mp)<=n_l);
end;%for nk = 1:n_k;
C_(nmn,nmp,:) = C_tmp(:);
end;%for nmp = 1:n_m_max;
end;%for nmn = 1:n_m_max;

X_ = zeros(n_m_max,n_m_max,n_beta);
for nbeta = 1:n_beta;
X_(:,:,nbeta) = recenter2(fft2(recenter2(squeeze(C_(:,:,nbeta)))));
end;%for nbeta = 1:n_beta;

prows = 5;
pcols = ceil(n_beta/prows);
for nbeta = 1:n_beta;
subplot(prows,pcols,nbeta); imagesc(real(squeeze(X_(:,:,nbeta)))); set(gca,'XTick',[],'YTick',[]);
end;%for nbeta = 1:n_beta;


