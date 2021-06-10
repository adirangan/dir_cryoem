function [n_all,n_sub_,k_p_r_all_,azimu_b_all_,polar_a_all_,weight_all_,n_k_p_r,k_p_r_,w_k_p_r_,k_c_0_,k_c_1_,k_c_2_] = sample_sphere_5(verbose,k_p_r_max,eq_d,TorL) ;
% generate azimu_b_ and polar_a_ arrays sampled on sphere of radius k_p_r_max at equatorial_distance eq_d ;
% Note: azimu_b_ in [0,2*pi], polar_a_ in [0,1*pi];
% Note that: ;
% \int (sphere of radius K) exp(+i*k*delta) = I1 = I2 = I3 = I4, where: ;
% I1 = dblquad(@(k,polar_a) 2*pi*cos(delta.*k.*cos(polar_a)).*k.^2.*sin(polar_a),0,Kmax,0,pi);
% I2 = 4*pi*(1/delta^1)*quad(@(k) k.^2.*sin(delta*k)./k,0,Kmax);
% I3 = 4*pi*(1/delta^3)*quad(@(k) k.*sin(k),0,Kmax*delta);
% I4 = 4*pi*(1/delta^3)*(sin(Kd) - Kd*cos(Kd)); 
if (nargin<1);
verbose=1;
k_p_r_max = 7.5;
eq_d_ = 0.5.^[-2:0.25:3]; n_eq_d = length(eq_d_);
E_ = zeros(n_eq_d,1);
for neq_d = 1:n_eq_d;
eq_d = eq_d_(neq_d);
sample_sphere_5(verbose,k_p_r_max,eq_d,'L');
end;%for neq_d = 1:n_eq_d;
disp('returning'); return;
end;%if (nargin<1);

if (verbose>1); disp(sprintf(' %% [entering sample_sphere_5]')); end;
n_k_p_r = 1+ceil(k_p_r_max/eq_d); [b_jx_,b_jw_] = jacpts(n_k_p_r,0,2);
k_p_r_ = (b_jx_+1.0)*k_p_r_max/2; w_k_p_r_ = b_jw_*(k_p_r_max/2)^3;
n_all = 0; n_sub_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
[n_sub] = sample_shell_5(k_p_r,eq_d,TorL) ;
n_sub_(1+nk_p_r) = n_all;
n_all = n_all + n_sub;
end;%for nk_p_r=0:n_k_p_r-1;
k_p_r_all_ = zeros(n_all,1); azimu_b_all_ = zeros(n_all,1); polar_a_all_ = zeros(n_all,1); weight_all_ = zeros(n_all,1);
ix=0;
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
[n_sub,azimu_b_sub_,polar_a_sub_,weight_sub_] = sample_shell_5(k_p_r,eq_d,TorL) ;
ij_ = ix + (0:n_sub-1);
k_p_r_all_(1+ij_) = k_p_r*ones(n_sub,1);
azimu_b_all_(1+ij_) = azimu_b_sub_;
polar_a_all_(1+ij_) = polar_a_sub_;
weight_all_(1+ij_) = weight_sub_ * w_k_p_r_(1+nk_p_r) / k_p_r.^2;
ix = ix + n_sub;
end;%for nk_p_r=0:n_k_p_r-1;

delta_ = [0.51;0.75;1.9]; delta_012 = fnorm(delta_); delta_01 = fnorm(delta_(1:2)); 
delta_azimu_b = atan2(delta_(1+1),delta_(1+0)); delta_polar_a = atan2(delta_01,delta_(1+2));
Ix = 4*pi*(1/delta_012^3)*( sin(k_p_r_max*delta_012) - (k_p_r_max*delta_012)*cos(k_p_r_max*delta_012) );
k_c_0_all_ = k_p_r_all_.*cos(azimu_b_all_).*sin(polar_a_all_);
k_c_1_all_ = k_p_r_all_.*sin(azimu_b_all_).*sin(polar_a_all_);
k_c_2_all_ = k_p_r_all_.*cos(polar_a_all_);
f_all_ = exp(+i* ( k_c_0_all_*delta_(1+0) + k_c_1_all_*delta_(1+1) + k_c_2_all_*delta_(1+2) ));
Iq = sum(f_all_.*weight_all_);
if (verbose>0); disp(sprintf(' %% k_p_r_max: %0.2f; eq_d %0.2f; error %0.16f',k_p_r_max,eq_d,fnorm(Ix-Iq)/fnorm(Ix))); end;

if (nargout>8);
k_c_0_ = k_p_r_all_.*cos(azimu_b_all_).*sin(polar_a_all_);
k_c_1_ = k_p_r_all_.*sin(azimu_b_all_).*sin(polar_a_all_);
k_c_2_ = k_p_r_all_.*cos(polar_a_all_);
end;%if (nargout>8);

if (verbose>1); disp(sprintf(' %% [finished sample_sphere_5]')); end;


