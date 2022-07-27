function ...
[ ...
 n_k_all ...
,n_k_all_csum_ ...
,k_p_r_all_ ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,weight_3d_k_all_ ...
,weight_shell_k_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,k_c_0_all_ ...
,k_c_1_all_ ...
,k_c_2_all_ ...
,J_node_ ...
,J_weight_ ...
,J_chebfun_ ...
,J_polyval_ ...
] = ...
sample_sphere_7( ...
 verbose ...
,k_p_r_max ...
,k_eq_d ...
,str_T_vs_L ...
,flag_unif_vs_adap ...
) ;
% generate k_p_azimu_b_ and k_p_polar_a_ arrays sampled on sphere of radius k_p_r_max at equatorial_distance k_eq_d ;
% Note: k_p_azimu_b_ in [0,2*pi], k_p_polar_a_ in [0,1*pi];
% Note that: ;
% \int (sphere of radius K) exp(+i*k*delta) = I1 = I2 = I3 = I4, where: ;
% I1 = dblquad(@(k,k_p_polar_a) 2*pi*cos(delta.*k.*cos(k_p_polar_a)).*k.^2.*sin(k_p_polar_a),0,Kmax,0,pi);
% I2 = 4*pi*(1/delta^1)*quad(@(k) k.^2.*sin(delta*k)./k,0,Kmax);
% I3 = 4*pi*(1/delta^3)*quad(@(k) k.*sin(k),0,Kmax*delta);
% I4 = 4*pi*(1/delta^3)*(sin(Kd) - Kd*cos(Kd)); 
% ;
% Inputs: ;
% verbose: integer verbosity level. ;
% k_p_r_max: real largest value of k (shell radius). ;
% k_eq_d: real equatorial distance required on each shell (for sampling). ;
% str_T_vs_L: string either 'T' for tschebycheff or 'L' [default] for legendre quadrature in the polar_a-direction (errors are comparable). ;
% flag_unif_vs_adap: integer either '1' for uniform or '0' [default] for adaptive spherical grid (in k_p_r) over each k_p_r-shell. ;
% ;
% Outputs: ;
% n_k_all: integer number of total points. ;
% n_k_all_csum_: integer array of size n_k_p_r. n_k_all_csum_(nk_p_r) is the starting index (0-based) for shell nk_p_r. ;
% k_p_r_all_: real array of size n_k_all. k-values for each point. ;
% k_p_azimu_b_all_: real array of size n_k_all. azimu_b values for each point. ;
% k_p_polar_a_all_: real array of size n_k_all. polar_a values for each point. ;
% weight_3d_k_all_: real array of size n_k_all. quadrature weights for each point: designed so that sum(weight_3d_k_all_) = volume of sphere = 4/3*pi*k_p_r_max^3. ;
% weight_shell_k_: real array of size n_k_all. quadrature weights for each shell (in sequence). ;
% n_k_p_r: integer number of shells. ;
% k_p_r_: real array of size n_k_p_r. k-values for each shell. ;
% weight_3d_k_p_r_: real array of size n_k_p_r. quadrature weights for each shell. These are designed so that 4*pi*sum(weight_3d_k_p_r_) = volume of sphere. ;
% k_c_0_all_: real array of size n_k_all: k_c_0 values for each point. ;
% k_c_1_all_: real array of size n_k_all: k_c_1 values for each point. ;
% k_c_2_all_: real array of size n_k_all: k_c_2 values for each point. ;
% J_node_: real array of size n_k_p_r: list of quadrature nodes on interval [-1,+1]. ;
% J_weight_: real array of size n_k_p_r: list of quadrature weights on interval [-1,+1]. ;
% J_chebfun_: chebfun cell array of size 1+n_k_p_r: Jv_{1+nk_p_r} is the chebfun associated with the jacobi-polynomial (0,2) of order nk_p_r. ;
% J_polyval_: matrix of size (1+n_k_p_r)*(n_k_p_r): J_polyval_(1+nk_p_r,:) is equal to J_chebfun_{1+nk_p_r} evaluted at the quadrature nodes. ;
% ;

na=0;
if (nargin<1+na); verbose=[]; end; na=na+1;
if (nargin<1+na); k_p_r_max=[]; end; na=na+1;
if (nargin<1+na); k_eq_d=[]; end; na=na+1;
if (nargin<1+na); str_T_vs_L=[]; end; na=na+1;
if (nargin<1+na); flag_unif_vs_adap=[]; end; na=na+1;

if (nargin<1);
verbose=2;
k_p_r_max = 48/(2*pi);
k_eq_d_ = 0.5.^[-2:0.25:3]; n_k_eq_d = length(k_eq_d_);
E_ = zeros(n_k_eq_d,1);
for nk_eq_d = 1:n_k_eq_d;
k_eq_d = k_eq_d_(nk_eq_d);
sample_sphere_7(verbose,k_p_r_max,k_eq_d,'T',0);
sample_sphere_7(verbose,k_p_r_max,k_eq_d,'L',0);
sample_sphere_7(verbose,k_p_r_max,k_eq_d,'T',1);
sample_sphere_7(verbose,k_p_r_max,k_eq_d,'L',1);
end;%for nk_eq_d = 1:n_k_eq_d;
disp('returning'); return;
end;%if (nargin<1);

if isempty(str_T_vs_L); str_T_vs_L = 'L'; end;
if isempty(flag_unif_vs_adap); flag_unif_vs_adap=0; end;

if (verbose>2); disp(sprintf(' %% [entering sample_sphere_7]')); end;
n_k_p_r = 1+ceil(k_p_r_max/k_eq_d); [a_jx_,a_jw_] = jacpts(n_k_p_r,0,2);
k_p_r_ = (a_jx_+1.0)*k_p_r_max/2; weight_3d_k_p_r_ = a_jw_*(k_p_r_max/2)^3;
J_node_ = a_jx_ ; J_weight_ = a_jw_ ;
J_chebfun_ = cell(1+n_k_p_r,1+n_k_p_r); J_polyval_ = zeros(1+n_k_p_r,n_k_p_r);
for nk_p_r=0:n_k_p_r;
J_chebfun = jacpoly(nk_p_r,0,2)*sqrt(2*nk_p_r+3)/sqrt(8);
J_polyval_(1+nk_p_r,:) = J_chebfun(a_jx_);
J_chebfun_{1+nk_p_r} = J_chebfun;
end;%for nk_p_r=0:n_k_p_r;
n_k_all = 0; n_k_all_csum_ = zeros(1+n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
if flag_unif_vs_adap==0;
[n_k_all_csum] = sample_shell_5(k_p_r,k_eq_d,str_T_vs_L) ;
end;%if flag_unif_vs_adap==0;
if flag_unif_vs_adap==1;
[n_k_all_csum] = sample_shell_5(k_p_r,k_eq_d*k_p_r/k_p_r_max,str_T_vs_L) ;
end;%if flag_unif_vs_adap==1;
n_k_all_csum_(1+nk_p_r) = n_k_all;
n_k_all = n_k_all + n_k_all_csum;
end;%for nk_p_r=0:n_k_p_r-1;
n_k_all_csum_(1+n_k_p_r) = n_k_all;
k_p_r_all_ = zeros(n_k_all,1); k_p_azimu_b_all_ = zeros(n_k_all,1); k_p_polar_a_all_ = zeros(n_k_all,1); weight_3d_k_all_ = zeros(n_k_all,1); weight_shell_k_ = zeros(n_k_all,1);
ix=0;
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
if flag_unif_vs_adap==0;
[n_k_all_csum,k_p_azimu_b_sub_,k_p_polar_a_sub_,weight_k_sub_] = sample_shell_5(k_p_r,k_eq_d,str_T_vs_L) ;
end;%if flag_unif_vs_adap==0;
if flag_unif_vs_adap==1;
[n_k_all_csum,k_p_azimu_b_sub_,k_p_polar_a_sub_,weight_k_sub_] = sample_shell_5(k_p_r,k_eq_d*k_p_r/k_p_r_max,str_T_vs_L) ;
end;%if flag_unif_vs_adap==1;
ij_ = ix + (0:n_k_all_csum-1);
k_p_r_all_(1+ij_) = k_p_r*ones(n_k_all_csum,1);
k_p_azimu_b_all_(1+ij_) = k_p_azimu_b_sub_;
k_p_polar_a_all_(1+ij_) = k_p_polar_a_sub_;
weight_shell_k_(1+ij_) = weight_k_sub_;
weight_3d_k_all_(1+ij_) = weight_k_sub_ * weight_3d_k_p_r_(1+nk_p_r) / k_p_r.^2;
ix = ix + n_k_all_csum;
end;%for nk_p_r=0:n_k_p_r-1;

if (verbose>1);
delta_ = [0.51;0.75;1.9]; delta_012 = fnorm(delta_); delta_01 = fnorm(delta_(1:2)); 
delta_k_p_azimu_b = atan2(delta_(1+1),delta_(1+0)); delta_k_p_polar_a = atan2(delta_01,delta_(1+2));
Ix = 4*pi*(1/delta_012^3)*( sin(k_p_r_max*delta_012) - (k_p_r_max*delta_012)*cos(k_p_r_max*delta_012) );
k_c_0_all_ = k_p_r_all_.*cos(k_p_azimu_b_all_).*sin(k_p_polar_a_all_);
k_c_1_all_ = k_p_r_all_.*sin(k_p_azimu_b_all_).*sin(k_p_polar_a_all_);
k_c_2_all_ = k_p_r_all_.*cos(k_p_polar_a_all_);
f_all_ = exp(+i* ( k_c_0_all_*delta_(1+0) + k_c_1_all_*delta_(1+1) + k_c_2_all_*delta_(1+2) ));
Iq = sum(f_all_.*weight_3d_k_all_);
disp(sprintf(' %% sample_sphere_7: k_p_r_max: %0.2f; k_eq_d %0.2f; n_k_all %.8d n_k_shell %.6d flag_unif_vs_adap %d str_T_vs_L %s; error %0.16f',k_p_r_max,k_eq_d,n_k_all,n_k_all_csum_(end)-n_k_all_csum_(end-1),flag_unif_vs_adap,str_T_vs_L,fnorm(Ix-Iq)/fnorm(Ix))); 
end;%if (verbose>0); 

if (nargout>10);
k_c_0_all_ = k_p_r_all_.*cos(k_p_azimu_b_all_).*sin(k_p_polar_a_all_);
k_c_1_all_ = k_p_r_all_.*sin(k_p_azimu_b_all_).*sin(k_p_polar_a_all_);
k_c_2_all_ = k_p_r_all_.*cos(k_p_polar_a_all_);
end;%if (nargout>10);

if (verbose>2); disp(sprintf(' %% [finished sample_sphere_7]')); end;


