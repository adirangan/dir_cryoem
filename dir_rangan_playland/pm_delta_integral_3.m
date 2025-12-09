function ...
[ ...
 I_pos__ ...
,I_neg_ ...
] = ...
pm_delta_integral_3( ...
 n_k_p_r ...
,k_p_r_ ...
,delta_sigma ...
,l_max ...
,pm_delta_integral_tolerance ...
);
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
% where k = k_p_r_(nk_p_r). ;

flag_verbose=0;
str_thisfunction = 'pm_delta_integral_3';
if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

na=0;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); k_p_r_=[]; end; na=na+1;
if (nargin<1+na); delta_sigma=[]; end; na=na+1;
if (nargin<1+na); l_max=[]; end; na=na+1;
if (nargin<1+na); pm_delta_integral_tolerance=[]; end; na=na+1;

if isempty(pm_delta_integral_tolerance); pm_delta_integral_tolerance = 1e-3; end; %<-- single-digit precision for I_neg_. ;
if (pm_delta_integral_tolerance<=0); pm_delta_integral_tolerance = 1e-3; end; %<-- single-digit precision for I_neg_. ;

%%%%%%%%;
if nargin<4;
%%%%;
flag_verbose=1;
n_k_p_r = 9; k_p_r_ = transpose(linspace(1,9,n_k_p_r)); delta_sigma = 0.05; l_max = 32;
tmp_t = tic();
[I_pos_0__,I_neg_0_] = pm_delta_integral_0(n_k_p_r,k_p_r_,delta_sigma,l_max,pm_delta_integral_tolerance);
tmp_t = toc(tmp_t); disp(sprintf(' %% pm_delta_integral_0: %0.6fs',tmp_t));
tmp_t = tic();
[I_pos_1__,I_neg_1_] = pm_delta_integral_1(n_k_p_r,k_p_r_,delta_sigma,l_max,pm_delta_integral_tolerance);
tmp_t = toc(tmp_t); disp(sprintf(' %% pm_delta_integral_1: %0.6fs',tmp_t));
tmp_t = tic();
[I_pos_2__,I_neg_2_] = pm_delta_integral_2(n_k_p_r,k_p_r_,delta_sigma,l_max,pm_delta_integral_tolerance);
tmp_t = toc(tmp_t); disp(sprintf(' %% pm_delta_integral_2: %0.6fs',tmp_t));
tmp_t = tic();
[I_pos_3__,I_neg_3_] = pm_delta_integral_3(n_k_p_r,k_p_r_,delta_sigma,l_max,pm_delta_integral_tolerance);
tmp_t = toc(tmp_t); disp(sprintf(' %% pm_delta_integral_3: %0.6fs',tmp_t));
fnorm_disp(flag_verbose,'I_pos_0__',I_pos_0__,'I_pos_1__',I_pos_1__,' %%<-- should be zero');
fnorm_disp(flag_verbose,'I_pos_0__',I_pos_0__,'I_pos_2__',I_pos_2__,' %%<-- should be zero');
fnorm_disp(flag_verbose,'I_pos_0__',I_pos_0__,'I_pos_3__',I_pos_3__,' %%<-- should be zero');
fnorm_disp(flag_verbose,'I_neg_0_',I_neg_0_,'I_neg_1_',I_neg_1_,' %%<-- should be zero');
fnorm_disp(flag_verbose,'I_neg_0_',I_neg_0_,'I_neg_2_',I_neg_1_,' %%<-- should be zero');
fnorm_disp(flag_verbose,'I_neg_0_',I_neg_0_,'I_neg_3_',I_neg_1_,' %%<-- should be zero');
%%%%;
% compare to ../dir_rangan_python/dir_pymat/test_ylgndr_2.mat. ;
%%%%;
dir_pymat = '../dir_rangan_python/dir_pymat';
fname_pymat = sprintf('%s/test_pm_delta_integral_3.mat',dir_pymat);
if ~exist(fname_pymat,'file'); disp(sprintf(' %% Warning, %s not found',fname_pymat)); end;
if  exist(fname_pymat,'file');
disp(sprintf(' %% %s found, loading',fname_pymat));
tmp_ = load(fname_pymat);
fnorm_disp(flag_verbose,'I_pos_3__',I_pos_3__,'tmp_.I_pos_3__',tmp_.I_pos_3__,' %%<-- should be <1e-6');
fnorm_disp(flag_verbose,'I_neg_3_',I_neg_3_,'tmp_.I_neg_3_',tmp_.I_neg_3_,' %%<-- should be <1e-6');
end;%if  exist(fname_pymat,'file');
%%%%;
end;%if nargin<4;
%%%%%%%%;

twopi = 2*pi;

%%%%%%%%;
% Note that, defining: ;
% I_neg = @(ks) 2*pi*(0.5*ks*sqrt(pi/2)*exp(-ks*ks/4)*(besseli(-0.5,ks*ks/4) - besseli(+0.5,ks*ks/4)));
% we see that I_neg is within ~1e-6 of 2*pi for ks<=1e-3. ;
%%%%%%%%;

ks_ = twopi*k_p_r_*delta_sigma;
ksks4_ = ks_.^2/4;
I_neg_ = zeros(n_k_p_r,1);
tmp_index_ = efind(ks_> pm_delta_integral_tolerance);
I_neg_(1+tmp_index_) = twopi.*(0.5.*ks_(1+tmp_index_).*sqrt(pi/2).*exp(-ksks4_(1+tmp_index_)).*(besseli(-0.5,ksks4_(1+tmp_index_)) - besseli(+0.5,ksks4_(1+tmp_index_))));
tmp_index_ = efind(ks_<=pm_delta_integral_tolerance);
I_neg_(1+tmp_index_) = twopi;

I_pos__ = zeros(n_k_p_r,n_k_p_r);
[k_p_r_0__,k_p_r_1__] = ndgrid(k_p_r_,k_p_r_);
tmp_z__ = twopi^2 .* k_p_r_0__.*k_p_r_1__ .* delta_sigma^2;
I_pos__ = twopi .* exp(-twopi^2.*0.5.*(k_p_r_0__.^2 + k_p_r_1__.^2).*delta_sigma^2) .* exp(tmp_z__);

if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;


