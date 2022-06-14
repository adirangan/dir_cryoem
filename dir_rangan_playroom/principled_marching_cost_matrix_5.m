function [ ...
 X__ ...
,X_weight_r_ ...
,X_ori__ ...
,X_tau__ ...
,weight_so3 ...
,n_m_max ...
,polar_a_ ...
,azimu_b_ ...
,gamma_z_ ...
] = ...
principled_marching_cost_matrix_5( ...
 n_k_p_r ...
,k_p_r_ ...
,weight_k_p_r_ ...
,l_max_ ...
,n_molecule ...
,molecule_density_ ...
,a_k_Y__ ...
,CTF_k_p_r_xcor_kk__ ...
,delta_sigma ...
,pm_delta_integral_tolerance ...
);
% more efficient calculation of volumetric cost. ;
% scales by inverse-standard-deviations. ;
% scales by cross-correlation of CTF_k_p_r_. ;
% integrates over multiple molecules. ;
% accounts for (isotropic gaussian) distribution of translations with standard-deviation delta_sigma. ;

if (nargin<10); pm_delta_integral_tolerance = 1e-2; end;
if isempty(pm_delta_integral_tolerance); pm_delta_integral_tolerance = 1e-2; end;
if (pm_delta_integral_tolerance<=0); pm_delta_integral_tolerance = 1e-2; end;

verbose=1;
if (isempty(n_molecule)); n_molecule = 1; molecule_density_ = ones(n_molecule,1); end;
if (n_molecule==0); n_molecule = 1; molecule_density_ = ones(n_molecule,1); end;
%%%%%%%%;
n_lm_ = (1+l_max_).^2;
n_lm_max = max(n_lm_);
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
l_max_max = max(l_max_); dWtdkd__l_max_max = 2*l_max_max;
m_max_ = -l_max_max : +l_max_max;
n_m_max = 1+2*l_max_max;
n_polar_a = n_m_max; polar_a_ = linspace(-pi,pi,n_polar_a+1); polar_a_ = polar_a_(1:end-1);
weight_so3 = (2*pi)*(2*pi)*4; %<-- total volume of so3. ;
a_k_Y___ = zeros(n_lm_max,n_k_p_r);
for nmolecule=0:n_molecule-1;
for nk_p_r=0:n_k_p_r-1;
tmp_ij_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
a_k_Y___(1:n_lm_(1+nk_p_r),1+nk_p_r,1+nmolecule) = a_k_Y__(1+tmp_ij_,1+nmolecule);
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nmolecule=0:n_molecule-1;
if (isempty(molecule_density_)); molecule_density_ = ones(n_molecule,1); end;
molecule_density_ = molecule_density_/sum(molecule_density_);

%%%%%%%%;
[I_pos__,I_neg_] = pm_delta_integral_1(n_k_p_r,k_p_r_,delta_sigma,l_max_max,pm_delta_integral_tolerance);
I_pos__ = I_pos__/(2*pi);
I_neg__ = (I_neg_*transpose(I_neg_))/(2*pi)^2;
%%%%%%%%;
X_weight_r_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
X_weight_r_(1+nk_p_r) = sqrt(weight_k_p_r_(1+nk_p_r));
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
X_ori__ = zeros(n_k_p_r,n_k_p_r);
X_tau__ = zeros(n_k_p_r,n_k_p_r);
for nmolecule_0=0:n_molecule-1; for nmolecule_1=0:n_molecule-1;
molecule_density_0 = molecule_density_(1+nmolecule_0); 
molecule_density_1 = molecule_density_(1+nmolecule_1); 
molecule_density = molecule_density_0*molecule_density_1;
for nk_p_r_0=0:n_k_p_r-1;for nk_p_r_1=0:n_k_p_r-1;
X_weight_r_0 = X_weight_r_(1+nk_p_r_0); %<-- volume element = inverse-variance associated with nk_p_r_0. ;
X_weight_r_1 = X_weight_r_(1+nk_p_r_1); %<-- volume element = inverse-variance associated with nk_p_r_1. ;
X_weight_r_01 = X_weight_r_0*X_weight_r_1; %<-- product of inverse-standard-deviations associated with nk_p_r_0 and nk_p_r_1. ;
X_ori__(1+nk_p_r_0,1+nk_p_r_1) = X_ori__(1+nk_p_r_0,1+nk_p_r_1) + molecule_density*register_spharm_to_spharm_2(verbose,1,1,1,l_max_max,a_k_Y___(:,1+nk_p_r_0,1+nmolecule_0),a_k_Y___(:,1+nk_p_r_1,1+nmolecule_1)) * X_weight_r_01 * CTF_k_p_r_xcor_kk__(1+nk_p_r_0,1+nk_p_r_1);
X_tau__(1+nk_p_r_0,1+nk_p_r_1) = X_tau__(1+nk_p_r_0,1+nk_p_r_1) + molecule_density*(4*pi)^2 * conj(a_k_Y___(1+0,1+nk_p_r_0,1+nmolecule_0))*a_k_Y___(1+0,1+nk_p_r_1,1+nmolecule_1) * X_weight_r_01 * CTF_k_p_r_xcor_kk__(1+nk_p_r_0,1+nk_p_r_1);
end;end;%for nk_p_r_0=0:n_k_p_r-1;for nk_p_r_1=nk_p_r_0:n_k_p_r-1;
end;end;%for nmolecule_0=0:n_molecule-1; for nmolecule_1=0:n_molecule-1;
%%%%%%%%;
X__ = real(X_ori__)*weight_so3.*I_pos__ - real(X_tau__).*I_neg__;
%%%%%%%%;

