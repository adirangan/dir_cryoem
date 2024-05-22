function ...
a_k_Y_ = ...
plane_wave_expansion_1( ...
 n_k_p_r ...
,k_p_r_ ...
,x_ ...
,l_max_ ...
);
% plane-wave expansion in terms of spherical-bessel-functions and spherical-harmonics. ;
% we assume plane-wave is of the form exp( +i * < 2*pi*k_ , x_ > ). ;

na=0;
if (nargin<1+na); n_k_p_r = []; end; na=na+1;
if (nargin<1+na); k_p_r_ = []; end; na=na+1;
if (nargin<1+na); x_ = []; end; na=na+1;
if (nargin<1+na); l_max_ = []; end; na=na+1;

verbose=0;

l_max_max = max(l_max_);
n_lm_ = (1+l_max_).^2;
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);

x_p_r = fnorm(x_); x_p_r01 = fnorm(x_(1:2)); x_p_azimu_b = atan2(x_(2),x_(1)); x_p_polar_a = atan2(x_p_r01,x_(3));
t_p_r_ = 2*pi*k_p_r_*x_p_r;
Ylm_x__ = get_Ylm__(1+l_max_max,0:l_max_max,1,x_p_azimu_b,x_p_polar_a,0);
a_k_Y_ = zeros(n_lm_sum,1);
na=0;
for nk_p_r=0:n_k_p_r-1;
t_p_r = t_p_r_(1+nk_p_r);
n_lm = n_lm_(1+nk_p_r);
l_max = l_max_(1+nk_p_r);
for l_val=0:l_max;
Ylm_d_x_ = Ylm_x__{1+l_val}(:);
if abs(t_p_r)>=1e-12; jl = besselj(l_val+0.5,t_p_r).*sqrt(pi/(2*t_p_r)); end;
if abs(t_p_r)< 1e-12; jl = 1*(l_val==0) + 0*(l_val> 0); end;
a_k_Y_(1+na+[0:1+2*l_val - 1]) = 4*pi*(i^l_val)*jl*conj(Ylm_d_x_);
na=na+1+2*l_val;
end;%for l_val=0:l_max;
end;%for nk_p_r=0:n_k_p_r-1;
assert(na==n_lm_sum);




