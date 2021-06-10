function [legendre_evaluate_ljm___,legendre_evaluate_mlj___,expil__,expi__] = legendre_evaluate_ljm___0(l_max,cos_polar_a_,n_azimu_b);
% evaluates associated legendre polynomials at cos_polar_a_. ;
% stores output in stacked matrix format: ;
% legendre_evaluate_ljm___(1+l_val,1+npolar_a,1+l_max+m_val) contains: ;
% Z(l_val;m_val)*P(;_val;m_val)(cos_polar_a_(1+npolar_a)). ;

n_polar_a = numel(cos_polar_a_);
n_l_val = 1 + 1*l_max;
n_lm_max = (1+l_max)^2;
n_m_abs = 1 + 1*l_max;
n_m_val = 1 + 2*l_max;

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

legendre_evaluate_ljm___ = permute(legendre_evaluate_lmj___,[1,3,2]);

if nargout>1;
legendre_evaluate_mlj___ = permute(legendre_evaluate_ljm___,[3,1,2]);
expil = exp(+i*2*pi*l_max/n_azimu_b); 
expil__ = sparse(1:n_azimu_b,1:n_azimu_b,expil.^transpose([0:n_azimu_b-1]),n_azimu_b,n_azimu_b);
expi = exp(-i*2*pi/n_azimu_b); 
expi__ = expi.^(transpose([0:n_m_val-1])*[0:n_azimu_b-1]);
end;%if nargout>1;
