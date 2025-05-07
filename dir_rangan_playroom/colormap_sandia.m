function output = colormap_sandia(n_c,gamma);
if (nargin<1); n_c = 64; gamma = 1.0; end;

c_ = zeros(n_c,3);
%%%%;
t_ = transpose(linspace(0,1,n_c));
t2_ = min(2*t_,2*(1-t_));
s_ = transpose(linspace(-1,1,n_c));
s2_ = abs(s_);
sp_ = +max(0,+s_);
sn_ = +max(0,-s_);
c_ = [sn_,sp_,sn_].^(gamma);
%%%%;
output = c_;

