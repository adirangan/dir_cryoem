function ...
[ ...
 output ...
] = ...
I_PP_0( ...
 K ...
,delta_S_ ...
,delta_M_ ...
);
%%%%%%%%;
% calculates the integral: ;
% \int_{k=0}^{K} \int_{psi=0}^{2*pi} conj[S] * M dpsi * kdk ;
% where: ;
% S   = exp( +i*2*pi*k*delta_S * cos(psi - omega_S) ) ;
% M   = exp( +i*2*pi*k*delta_M * cos(psi - omega_M) ) ;
% where: ;
% delta_S_ = delta_S * [cos(omega_S);sin(omega_S)] ;
% delta_M_ = delta_M * [cos(omega_M);sin(omega_M)] ;
%%%%%%%%;
str_thisfunction = 'I_PP_0';

if nargin<1;
flag_verbose=1; flag_disp=0; nf=0;
if flag_verbose>0; disp(sprintf(' %% testing %s',str_thisfunction)); end;
%%%%%%%%;
% Now set up polar-quadrature-weights on disk. ;
%%%%%%%%;
k_int = 48;
k_eq_d_double = 0.25;
n_w_int = 4;
k_p_r_max = k_int/(2*pi); k_eq_d = k_eq_d_double/(2*pi); str_L = 'L';
[ ...
 n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
] = ...
get_weight_3d_1( ...
 flag_verbose ...
,k_p_r_max ...
,k_eq_d ...
,str_L ...
);
%%%%;
l_max_upb = round(2*pi*k_p_r_max);
l_max_max = min(l_max_upb,1+ceil(2*pi*k_p_r_(end)));
n_w_max = n_w_int*2*(l_max_max+1); n_w_0in = n_w_max; n_w_0in_ = n_w_max*ones(n_k_p_r,1);
[ ...
 n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_wk_ ...
,k_p_r_wk_ ...
,k_p_w_wk_ ...
,k_c_0_wk_ ...
,k_c_1_wk_ ...
] = ...
get_weight_2d_2( ...
 0*flag_verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,-1 ...
,n_w_0in_ ...
);
n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_]);
%%%%%%%%;
% Generate images, template and integral. ;
%%%%%%%%;
rng(3);
K = k_p_r_max;
delta_S_ = 1e-1*2*pi*rand(2,1); delta_S = fnorm(delta_S_); omega_S = atan2(delta_S_(1+1),delta_S_(1+0));
delta_M_ = 1e-1*2*pi*rand(2,1); delta_M = fnorm(delta_M_); omega_M = atan2(delta_M_(1+1),delta_M_(1+0));
S_k_p_wk_ = exp(2*pi*i*k_p_r_wk_.*delta_S.*cos(k_p_w_wk_ - omega_S));
M_k_p_wk_ = exp(2*pi*i*k_p_r_wk_.*delta_M.*cos(k_p_w_wk_ - omega_M));
%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figbig; p_row=2;p_col=2;np=0;
Slim_ = max(abs(S_k_p_wk_),[],'all')*[-1,+1];
Mlim_ = max(abs(M_k_p_wk_),[],'all')*[-1,+1];
subplot(p_row,p_col,1+np);np=np+1; imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(S_k_p_wk_),Slim_); axis image; axisnotick; title('real(S_k_p_wk_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1; imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(S_k_p_wk_),Slim_); axis image; axisnotick; title('imag(S_k_p_wk_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1; imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(M_k_p_wk_),Mlim_); axis image; axisnotick; title('real(M_k_p_wk_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1; imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(M_k_p_wk_),Mlim_); axis image; axisnotick; title('imag(M_k_p_wk_)','Interpreter','none');
end;%if flag_disp;
%%%%;
output_qua = sum(conj(S_k_p_wk_).*M_k_p_wk_.*weight_2d_wk_,'all')*(2*pi)^2;
[ ...
 output_mid ...
] = ...
I_PP_0( ...
 K ...
,delta_S_ ...
,delta_M_ ...
);
%%%%;
fnorm_disp(flag_verbose,'output_qua',output_qua,'output_mid',output_mid,' %<-- should be small (e.g., <1e-6 if quadrature resolved)');
%%%%;
if flag_verbose>0; disp(sprintf(' %% returning')); end; return;
end;%if nargin<1;

na=0;
if (nargin<1+na); K=[]; end; na=na+1;
if (nargin<1+na); delta_S_=[]; end; na=na+1;
if (nargin<1+na); delta_M_=[]; end; na=na+1;

if isempty(K); K = 48.0/(2*pi); end;
if isempty(delta_S_); delta_S_ = [1.0;0.0]; end;
if isempty(delta_M_); delta_M_ = [0.0;1.0]; end;

flag_1 = 1;

delta_T_ = delta_M_ - delta_S_ ;
delta_T = sqrt(delta_T_(1+0).^2 + delta_T_(1+1).^2) ;
omega_T = atan2(delta_T_(1+1),delta_T_(1+0)) ;
t = 2*pi*delta_T ;
tK = t*K ;
b0 = besselj(0,tK);
b2 = besselj(2,tK);
output = pi*K.^2*(b0 + b2);
