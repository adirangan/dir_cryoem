function f_val_ = shevalspherep(Y_,l_max,n_polar_a,n_azimu_b,polar_a_);
%      This subroutine evaluates a spherical harmonic expansion on an ;
%      nquad x nquadm grid on the unit sphere. ;
% ;
% --------------------------------------------------------------------- ;
%      INPUT: ;
% ;
%      Y_         : coefficients of spherical harmonic exp. (packed);
%      l_max      : number of terms in the orig. expansion;
%      n_polar_a  : number of quadrature nodes in polar_a;
%      n_azimu_b  : number of quadrature nodes in azimu_b;
%      polar_a_   : Legendre nodes in polar_a (x_j = cos polar_a_j).;
% ---------------------------------------------------------------------;
%      OUTPUT:;
% ;
%      f_val_     : function value on tensor product;
%                 mesh on target sphere. azimu_b is the fast variable, polar_a slow;
% ***********************************************************************;
%      converted to OPENMP, Barnett 6/22/16;
% ;
assert(numel(Y_)>=(1+l_max)^2);
f_val_ = zeros(n_azimu_b,n_polar_a);
f_tmp_ = zeros(1+2*l_max,n_polar_a);
expi_azimu_b = 0; expi_azimu_b_k = 0;
ynm__ = zeros(1+l_max,1+l_max);

for npolar_a=1:n_polar_a;
cos_polar_a = polar_a_(npolar_a);
sin_polar_a = sqrt(1 - cos_polar_a^2);
ynm__ = ylgndr(l_max,cos_polar_a);
for m_val=-l_max:+l_max;
m_abs = abs(m_val);
for l_val=m_abs:l_max;
tab = l_val*(l_val+1) + m_val + 1;
f_tmp_(1+l_max+m_val,npolar_a) = f_tmp_(1+l_max+m_val,npolar_a) +Y_(tab)*ynm__(1+l_val,1+m_abs);
end;%for l_val=m_abs:l_max;
end;%for m_val=-l_max:+l_max;
end;%for npolar_a=1:n_polar_a;

for npolar_a = 1:n_polar_a;
for nazimu_b = 1:n_azimu_b;
f_val_(nazimu_b,npolar_a) = 0;
expi_azimu_b_k = exp(2*pi*i*(nazimu_b-1)/n_azimu_b);
expi_azimu_b = expi_azimu_b_k^(-l_max);
for m_val = -l_max:+l_max;
f_val_(nazimu_b,npolar_a) = f_val_(nazimu_b,npolar_a) + f_tmp_(1+l_max+m_val,npolar_a)*expi_azimu_b;
expi_azimu_b = expi_azimu_b*expi_azimu_b_k;
end;%for m_val = -l_max:+l_max;
end;%for nazimu_b = 1:n_azimu_b;
end;%for npolar_a = 1:n_polar_a;
