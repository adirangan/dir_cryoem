function [I_pos__,I_neg_] = pm_delta_integral_0(n_k_p_r,k_p_r_,delta_sigma,l_max,pm_delta_integral_tolerance);
% calculates integrals associated with delta: ;
%
% Input: ;
% n_k_p_r: integer number of k-values. ;
% k_p_r_: real array of size n_k_p_r. k_values. ;
% delta_sigma: real standard-deviation for isotropic gaussian distribution fo delta-values. ;
% l_max: maximum l_val to use in expansion. ;
% pm_delta_integral_tolerance: real tolerance for warning (default 1e-2). ;
%
% Output: ;
% I_pos__ = real array of size (n_k_p_r,n_k_p_r). ;
% I_pos__(nk_p_r_0,nk_p_r_1) = \int dphi * 1/twopi/delta_sigma^2 * exp(-delta^2/2/delta_sigma^2) * delta * d_delta * d_omega * exp(+i*twopi*k_0*delta*cos(phi - omega)) * exp(-i*twopi*k_1*delta*cos(phi - omega)). ;
% where k_0 = k_p_r_(nk_p_r_0) and k_1 = k_p_r_(nk_p_r_1). ;
% I_neg_ = real array of size n_k_p_r. ;
% I_neg_(nk_p_r) = \int dphi * 1/twopi/delta_sigma^2 * exp(-delta^2/2/delta_sigma^2) * delta * d_delta * d_omega * exp(\pm i*twopi*k*delta*cos(phi - omega)) ;

if (nargin<5); pm_delta_integral_tolerance = 1e-2; end;
if isempty(pm_delta_integral_tolerance); pm_delta_integral_tolerance = 1e-2; end;
if (pm_delta_integral_tolerance<=0); pm_delta_integral_tolerance = 1e-2; end;

twopi = 2*pi;

rmu = @(d,s) exp( - d.^2 ./ (2*s.^2) ) / twopi ./ s.^2 .* d ;

I_neg_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
I_neg_(1+nk_p_r) = integral(@(d) besselj(0,twopi*k_p_r*d).*rmu(d,delta_sigma),0,+Inf)*twopi^2;
end;%for nk_p_r=0:n_k_p_r-1;

I_pos__ = zeros(n_k_p_r,n_k_p_r);
for nk_p_r_0=0:n_k_p_r-1;
for nk_p_r_1=0:n_k_p_r-1;
k_p_r_0 = k_p_r_(1+nk_p_r_0);
k_p_r_1 = k_p_r_(1+nk_p_r_1);
l_max = 32;
I_form_ = zeros(1+2*l_max,1);
for l_val=-l_max:l_max;
I_form_(1+l_max+l_val) = integral(@(d) besselj(l_val,twopi*k_p_r_0*d).*besselj(l_val,twopi*k_p_r_1*d).*rmu(d,delta_sigma),0,+Inf)*twopi^2;
end;%for l_val=0:l_max;
if (fnorm(I_form_(end))/fnorm(I_form_(1+l_max))>pm_delta_integral_tolerance); disp(sprintf(' %% Warning, increase l_max in pm_delta_integral_0: fnorm(I_form_(end))/fnorm(I_form_(1+l_max)) = %0.6f/%0.6f = %0.6f',fnorm(I_form_(end)),fnorm(I_form_(1+l_max)),fnorm(I_form_(end))/fnorm(I_form_(1+l_max)))); end;
I_form = sum(I_form_);
I_pos__(1+nk_p_r_0,1+nk_p_r_1) = I_form;
end;%for nk_p_r_1=0:n_k_p_r-1;
end;%for nk_p_r_0=0:n_k_p_r-1;


