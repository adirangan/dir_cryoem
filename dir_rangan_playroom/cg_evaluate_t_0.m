function a_k_Y_ = cg_evaluate_t_0(n_polar_a,n_azimu_b,a_k_p__,l_max);
% adjoint-evaluation of a_k_p__ on tensor spherical grid. ;
% Note that a_k_Y_ is ordered linearly, with m_val varying quickly and l_val varying slowly. ;
% Note that a_k_p__ is ordered with polar_a varying quickly, and azimu_b varying slowly. ;
%%%%%%%%;
if nargin<4;
cg_evaluate_n_0();
disp('returning'); return;
end;%if nargin<4;

n_l_val = 1 + 1*l_max;
n_lm_max = (1+l_max)^2;
n_m_abs = 1 + 1*l_max;
n_m_val = 1 + 2*l_max;
assert(n_azimu_b>=n_m_val); %<-- important previously for ifft. ;
polar_a_ = transpose(linspace(0,pi,n_polar_a));
cos_polar_a_ = cos(polar_a_);
azimu_b_ = transpose(linspace(0,2*pi,n_azimu_b+1));
azimu_b_ = azimu_b_(1:end-1);

legendre_normalization__ = zeros(n_l_val,n_m_abs);
for l_val=0:l_max;
tmp_a1 = ((1+2*l_val)/(4*pi));
for m_abs=0:l_val;
tmp_a2 = exp(lfactorial(l_val-m_abs) - lfactorial(l_val+m_abs));
tmp_a3 = sqrt(tmp_a1*tmp_a2);
legendre_normalization__(1+l_val,1+m_abs) = tmp_a3;
end;%for m_abs=0:l_val;
end;%for l_val=0:l_max;

legendre_evaluate_lmj___ = zeros(n_l_val,n_m_val,n_polar_a); %<-- will include normalization. ;
for l_val=0:l_max;
tmp_n_m_abs = 1 + 1*l_val;
tmp_P__ = legendre(l_val,cos_polar_a_,'unnorm');
tmp_ZP__ = diag(legendre_normalization__(1+l_val,1:tmp_n_m_abs))*reshape(tmp_P__,[tmp_n_m_abs,n_polar_a]);
legendre_evaluate_lmj___(1+l_val,1 + l_max + [0:+1:+l_val],:) = tmp_ZP__;
legendre_evaluate_lmj___(1+l_val,1 + l_max + [0:-1:-l_val],:) = tmp_ZP__;
end;%for l_val=0:l_max;

legendre_evaluate_mlj___ = permute(legendre_evaluate_lmj___,[2,1,3]);

%%%%%%%%;
% Our main goal is to compute: ;
% a_k_Y(l_val;m_val) = a_k_Y_(1+l_val^2+l_val+m_val) = ;
% \sum_{nazimu_b=0}^{n_azimu_b-1} \sum_{npolar_a=0}}^{n_polar_a-1} .... ;
%  conj(Ylm__(l_val;m_val)(polar_a,azimu_b)) * a_k_p__(1+npolar_a,1+nazimu_b), ;
% Where polar_a = polar_a_(npolar_a), and azimu_b = azimu_b_(nazimu_b). ;
% ;
% Note that conj(Ylm__(l_val;m_val)(polar_a,azimu_b)) can be decomposed as: ;
% conj(Ylm__(l_val;m_val)) = ... ;
% legendre_normalization(l_val;m_val) ... ;
% * legendre_evaluate(l_val;m_val)(cos(polar_a_(npolar_a))) ... ;
% * exp(-i*m_val*azimu_b_(nazimu_b)). ;
% ;
% We will first calculate: ;
% spherical_harmonic_unphased___(1+l_max+m_val,1+l_val,1+nazimu_b)
% \sum_{npolar_a=0}^{n_polar_a-1} ... ;
%  legendre_normalization(l_val;m_val) ... ;
%  * legendre_evaluate(l_val;m_val;cos_polar_a_(npolar_a))) ... ;
%  * a_k_p__(1+npolar_a,1+nazimu_b) ;
% ;
% And then calculate: ;
% a_k_Y_(1+l_val^2+l_val+m_val) = ;
% a_k_Y_+(1+l_val,1+l_max+m_val) = ;
% \sum_{nazimu_b=0}^{n_azimu_b-1} ... ;
% spherical_harmonic_unphased___(1+l_max+m_val,1+l_val,1+nazimu_b) ... ;
% * exp(-i*m_val*azimu_b_(nazimu_b)). ;
%%%%%%%%;

%%%%%%%%;
% First sum over npolar_a. ;
% Note that legendre_evaluate_mlj___ already accounts for the normalization. ;
%%%%%%%%;
spherical_harmonic_unphased___ = zeros(n_m_val,n_l_val,n_azimu_b);
spherical_harmonic_unphased___ = reshape(reshape(legendre_evaluate_mlj___,n_m_val*n_l_val,n_polar_a)*a_k_p__,[n_m_val,n_l_val,n_azimu_b]);

%%%%%%%%;
% Now sum over nazimu_b. ;
% Note that FFT computes: ;
%              N                                                  ;
%      X(k) = sum  x(n)*exp(-j*2*pi*(k-1)*(n-1)/N), 1 <= k <= N.  ;
%             n=1                                                 ;
% Whereas we would like to compute: ;
% a_k_Y_(1+l_val^2+l_val+m_val) = ;
% a_k_Y_lm__(1+l_val,1+l_max+m_val) = ;
% \sum_{nazimu_b=0}^{n_azimu_b-1} ... ;
% spherical_harmonic_unphased___(1+l_max+m_val,1+l_val,1+nazimu_b) ... ;
% * exp(-i*m_val*azimu_b_(nazimu_b)). ;
% ;
% Noting that: m_val = -l_max + nm_val, ;
% where nm_val ranges across [0:(2*l_max)] ;
% and that: azimu_b_(nazimu_b) = 2*pi*nazimu_b/n_azimu_b, ;
% We can rewrite this sum as: ;
% 
% a_k_Y_lm__(1+l_val,1+nm_val) = ;
% \sum_{nazimu_b=0}^{n_azimu_b-1} ... ;
% spherical_harmonic_unphased___(1+nm_val,1+l_val,1+nazimu_b) ... ;
% * exp(-i*2*pi*(nm_val-l_max)*nazimu_b/n_azimu_b). ;
% which equals: ;
% a_k_Y_lm__(1+l_val,1+nm_val) = ;
% \sum_{nazimu_b=0}^{n_azimu_b-1} ... ;
% spherical_harmonic_unphased___(1+nm_val,1+l_val,1+nazimu_b) ... ;
% * exp(+i*2*pi*l_max*nazimu_b/n_azimu_b) ... ;
% * exp(-i*2*pi*nm_val*nazimu_b/n_azimu_b). ;
%%%%%%%%;
assert(n_azimu_b>=n_m_val);
a_k_Y_lm__ = zeros(n_l_val,n_m_val);
expil = exp(+i*2*pi*l_max/n_azimu_b);
expil__ = sparse(1:n_azimu_b,1:n_azimu_b,expil.^transpose([0:n_azimu_b-1]),n_azimu_b,n_azimu_b);
expi = exp(-i*2*pi/n_azimu_b);
expi__ = expi.^(transpose([0:n_m_val-1])*[0:n_azimu_b-1]);
for l_val=0:l_max;
a_k_Y_lm__(1+l_val,:) = sum((squeeze(spherical_harmonic_unphased___(:,1+l_val,:)).*expi__)*expil__,2);
end;%for l_val=0:l_max;
a_k_Y_ = convert_spharm__to_spharm_0(l_max,a_k_Y_lm__);
%{
for l_val=0:l_max;
a_k_Y_(1+l_val^2+l_val+[-l_val:+l_val]) = a_k_Y_lm__(1+l_val,1+l_max+[-l_val:+l_val]);
end;%for l_val=0:l_max;
 %}





