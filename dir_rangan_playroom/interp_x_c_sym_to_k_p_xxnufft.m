function ...
S_k_p_wk_ = ...
interp_x_c_sym_to_k_p_xxnufft( ...
 n_x0 ...
,diameter_x0_c ...
,n_x1 ...
,diameter_x1_c ...
,flag_sym ...
,S_x_c_xx_ ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,weight_2d_wk_ ...
) ;
%%%%%%%%;
% if flag_sym==1, assumes the x-values are symmetric about zero. ;
%%%%%%%%;
verbose=0;
half_diameter_x0_c = diameter_x0_c/2.0;
half_diameter_x1_c = diameter_x1_c/2.0;
x_p_r_max = max(half_diameter_x0_c,half_diameter_x1_c);
if flag_sym==0;
x_c_0_ = linspace(-half_diameter_x0_c,+half_diameter_x0_c,1+n_x1);
x_c_0_ = transpose(x_c_0_(1:n_x1));
x_c_1_ = linspace(-half_diameter_x1_c,+half_diameter_x1_c,1+n_x1);
x_c_1_ = transpose(x_c_1_(1:n_x1));
end;%if flag_sym==0;
if flag_sym==1;
x_c_0_ = linspace(-half_diameter_x0_c,+half_diameter_x0_c,n_x1);
x_c_0_ = transpose(x_c_0_(1:n_x1));
x_c_1_ = linspace(-half_diameter_x1_c,+half_diameter_x1_c,n_x1);
x_c_1_ = transpose(x_c_1_(1:n_x1));
end;%if flag_sym==1;
dx0 = mean(diff(x_c_0_));
dx1 = mean(diff(x_c_1_));
[x_c_0__,x_c_1__] = ndgrid(x_c_0_,x_c_1_);
n_xx_c = n_x0*n_x1;
weight_xx_c_ = dx0*dx1;
%%%%%%%%;
n_w_sum = sum(n_w_);
k_c_0_wk_ = zeros(n_w_sum,1);
k_c_1_wk_ = zeros(n_w_sum,1);
na=0;
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
n_w = n_w_(1+nk_p_r);
for nw=0:n_w-1;
omega = (2.0d0*pi*nw)/max(1,n_w);
k_c_0_wk_(1+na) = k_p_r*cos(omega);
k_c_1_wk_(1+na) = k_p_r*sin(omega);
na=na+1;
end;%for nw=0:n_w-1;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
eta = pi/x_p_r_max; tmp_t = tic;
S_k_p_wk_ = xxnufft2d3(n_xx_c,x_c_0__(:)*eta,x_c_1__(:)*eta,S_x_c_xx_(:).*weight_xx_c_(:),-1,1e-12,n_w_sum,2*pi*k_c_0_wk_/eta,2*pi*k_c_1_wk_/eta)/sqrt(2*pi)/sqrt(2*pi) * (sqrt(2*pi)^2) ;
%%%%%%%%;

