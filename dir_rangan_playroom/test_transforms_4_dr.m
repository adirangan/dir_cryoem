flag_verbose = 1;

%%%%%%%%;
% test bessel-function recursion. ;
%%%%%%%%;
kd = rand(); tmp_eps = 1e-4;
tmp_J0 = besselj(0,kd);
tmp_J2 = besselj(2,kd);
tmp_dJ1 = (besselj(1,kd+tmp_eps)-besselj(1,kd-tmp_eps))/(2*tmp_eps);
disp(sprintf(' %% tmp_J0 - tmp_J2 + tmp_J2 vs 2.0*tmp_dJ1 + tmp_J2: %0.16f %%<-- should be <1e-6 ',(tmp_J0 - tmp_J2 + tmp_J2) - (2.0*tmp_dJ1 + tmp_J2)));

%%%%%%%%;
% test bessel-function product. ;
% \int_{0}^{1} x J_{\nu}(\alpha x) J_{\nu}(\beta x) = ;
% \frac{\beta J_{\nu-1}(\beta) J_{\nu}(\alpha) - \alpha J_{\nu-1}(\alpha) J_{\nu}(\beta)}{\alpha^{2}-\beta^{2}} ;
%%%%%%%%;
tmp_alpha = rand();
tmp_beta = rand();
tmp_nu = 1;
tmp_I = integral(@(x) x.*besselj(tmp_nu,tmp_alpha*x).*besselj(tmp_nu,tmp_beta*x),0,1);
tmp_J = (tmp_beta*besselj(tmp_nu-1,tmp_beta)*besselj(tmp_nu,tmp_alpha) - tmp_alpha*besselj(tmp_nu-1,tmp_alpha)*besselj(tmp_nu,tmp_beta))/(tmp_alpha^2-tmp_beta^2);
if (flag_verbose>0); disp(sprintf(' %% tmp_I vs tmp_J: %0.16f %%<-- should be <1e-6 ',fnorm(tmp_I-tmp_J)/fnorm(tmp_I))); end;

%%%%%%%%;
% \int x^n J_{n-1}(x) dx = x^n J_{n}(x) + C. ;
%%%%%%%%;
tmp_n = 3; tmp_a = rand(); tmp_b = rand();
tmp_J1 = @(x) x.^tmp_n .* besselj(tmp_n,x);
tmp_J = tmp_J1(tmp_b) - tmp_J1(tmp_a);
tmp_I = integral(@(x) x.^tmp_n .* besselj(tmp_n-1,x) , tmp_a,tmp_b);
if (flag_verbose>0); disp(sprintf(' %% tmp_I vs tmp_J: %0.16f %%<-- should be <1e-6 ',fnorm(tmp_I-tmp_J)/fnorm(tmp_I))); end;

%%%%%%%%;
% test bessel-function fourier-transform. ;
% \frac{1}{2\pi} \int_{0}^{2\pi} \exp(-i kx cos(\psi-omega)) \exp(-iq(\psi-\phi)) d\psi = exp(-iq(\omega - \phi + \pi/2)) J_{q}(kx) ;
% \frac{1}{2\pi} \int_{0}^{2\pi} \exp(+i kx cos(\psi-omega)) \exp(+iq(\psi-\phi)) d\psi = exp(+iq(\omega - \phi + \pi/2)) J_{q}(kx) ;
% \frac{1}{2\pi} \int_{0}^{2\pi} \exp(-i kx cos(\psi-omega)) \exp(+iq(\psi-\phi)) d\psi = exp(+iq(\omega - \phi + \pi/2 + \pi)) J_{q}(kx) ;
% \frac{1}{2\pi} \int_{0}^{2\pi} \exp(+i kx cos(\psi-omega)) \exp(-iq(\psi-\phi)) d\psi = exp(-iq(\omega - \phi + \pi/2 + \pi)) J_{q}(kx) ;
%%%%%%%%;
tmp_omega = 2*pi*rand(); tmp_phi = 2*pi*rand();
n_kx = 13;
kx_ = transpose(linspace(0,3,n_kx));
q_ = transpose(-4:+4); n_q = numel(q_);
[kx__,q__] = ndgrid(kx_,q_);
tmp_Inn__ = zeros(n_kx,n_q); tmp_Jnn__ = zeros(n_kx,n_q);
tmp_Ipp__ = zeros(n_kx,n_q); tmp_Jpp__ = zeros(n_kx,n_q);
tmp_Inp__ = zeros(n_kx,n_q); tmp_Jnp__ = zeros(n_kx,n_q);
for nkx=0:n_kx-1;
for nq=0:n_q-1;
kx = kx_(1+nkx);
q = q_(1+nq);
tmp_Inn = integral(@(psi) exp(-i*kx*cos(psi-tmp_omega)).*exp(-i*q*(psi-tmp_phi)),0,2*pi)/(2*pi);
tmp_Jnn = exp(-i*q*(tmp_omega - tmp_phi + pi/2))*besselj(q,kx);
tmp_Inn__(1+nkx,1+nq) = tmp_Inn; tmp_Jnn__(1+nkx,1+nq) = tmp_Jnn;
tmp_Ipp = integral(@(psi) exp(+i*kx*cos(psi-tmp_omega)).*exp(+i*q*(psi-tmp_phi)),0,2*pi)/(2*pi);
tmp_Jpp = exp(+i*q*(tmp_omega - tmp_phi + pi/2))*besselj(q,kx);
tmp_Ipp__(1+nkx,1+nq) = tmp_Ipp; tmp_Jpp__(1+nkx,1+nq) = tmp_Jpp;
tmp_Inp = integral(@(psi) exp(-i*kx*cos(psi-tmp_omega)).*exp(+i*q*(psi-tmp_phi)),0,2*pi)/(2*pi);
tmp_Jnp = exp(+i*q*(tmp_omega - tmp_phi + pi/2 + pi))*besselj(q,kx);
tmp_Inp__(1+nkx,1+nq) = tmp_Inp; tmp_Jnp__(1+nkx,1+nq) = tmp_Jnp;
tmp_Ipn = integral(@(psi) exp(+i*kx*cos(psi-tmp_omega)).*exp(-i*q*(psi-tmp_phi)),0,2*pi)/(2*pi);
tmp_Jpn = exp(-i*q*(tmp_omega - tmp_phi + pi/2 + pi))*besselj(q,kx);
tmp_Ipn__(1+nkx,1+nq) = tmp_Ipn; tmp_Jpn__(1+nkx,1+nq) = tmp_Jpn;
end;%for nq=0:n_q-1;
end;%for nkx=0:n_kx-1;
if (flag_verbose>0); disp(sprintf(' %% tmp_Inn__ vs tmp_Jnn__: %0.16f %%<-- should be <1e-6 ',fnorm(tmp_Inn__-tmp_Jnn__)/fnorm(tmp_Inn__))); end;
if (flag_verbose>0); disp(sprintf(' %% tmp_Ipp__ vs tmp_Jpp__: %0.16f %%<-- should be <1e-6 ',fnorm(tmp_Ipp__-tmp_Jpp__)/fnorm(tmp_Ipp__))); end;
if (flag_verbose>0); disp(sprintf(' %% tmp_Inp__ vs tmp_Jnp__: %0.16f %%<-- should be <1e-6 ',fnorm(tmp_Inp__-tmp_Jnp__)/fnorm(tmp_Inp__))); end;
if (flag_verbose>0); disp(sprintf(' %% tmp_Ipn__ vs tmp_Jpn__: %0.16f %%<-- should be <1e-6 ',fnorm(tmp_Ipn__-tmp_Jpn__)/fnorm(tmp_Ipn__))); end;

%%%%%%%%;
% test plane-wave integration on disk. ;
% \int_{0}^{K}\int_{0}^{2\pi} exp(i 2\pi k \delta cos(\psi - \omega)) d\psi kdk = \pi K^2 (J_{0}(2\pi K\delta)+J_{2}(2\pi K\delta)) ;
%%%%%%%%;
k_p_r_max = 48/(2*pi);
delta_max = 0.1;
delta_ = transpose(linspace(0,delta_max,64)); n_delta = numel(delta_);
tmp_I_ = zeros(n_delta,1);
tmp_J_ = zeros(n_delta,1);
for ndelta=0:n_delta-1;
delta = delta_(1+ndelta);
tmp_omega = 2*pi*rand();
tmp_I = integral2(@(k_p_r,psi) k_p_r.*exp(2*pi*i*k_p_r.*delta.*cos(psi-tmp_omega)),0,k_p_r_max,0,2*pi);
tmp_J = pi*k_p_r_max^2*(besselj(0,2*pi*k_p_r_max*delta) + besselj(2,2*pi*k_p_r_max*delta));
tmp_I_(1+ndelta) = tmp_I;
tmp_J_(1+ndelta) = tmp_J;
end;%for ndelta=0:n_delta-1;
if (flag_verbose>0); disp(sprintf(' %% tmp_I_ vs tmp_J_: %0.16f %%<-- should be <1e-6 ',fnorm(tmp_I_-tmp_J_)/fnorm(tmp_I_))); end;

%%%%%%%%;
% integral of plane-wave times planar-function. ;
% (i.e., strongly anisotropic CTF). ;
%%%%%%%%;
% \int_{0}^{K}\int_{0}^{2\pi} (-2ik cos(\psi-\phi)) \cdot ... ;
%                              exp(i 2\pi k \delta cos(\psi - \omega)) d\psi kdk ... ;
%                            = \pi K^2 Kcos(\phi-\omega) (J_{1}(2\pi K\delta)+J_{3}(2\pi K\delta)) ;
%%%%%%%%;
k_p_r_max = 48/(2*pi);
delta_max = 0.1;
delta_ = transpose(linspace(0,delta_max,64)); n_delta = numel(delta_);
tmp_I_ = zeros(n_delta,1);
tmp_J_ = zeros(n_delta,1);
for ndelta=0:n_delta-1;
delta = delta_(1+ndelta);
tmp_phi = 2*pi*rand();
tmp_omega = 2*pi*rand();
tmp_I = integral2(@(k_p_r,psi) k_p_r.*(-2*i*k_p_r.*cos(psi-tmp_phi)).*exp(2*pi*i*k_p_r.*delta.*cos(psi-tmp_omega)),0,k_p_r_max,0,2*pi);
tmp_J = pi*k_p_r_max^2*k_p_r_max*cos(tmp_phi-tmp_omega)*(besselj(1,2*pi*k_p_r_max*delta) + besselj(3,2*pi*k_p_r_max*delta));
tmp_I_(1+ndelta) = tmp_I;
tmp_J_(1+ndelta) = tmp_J;
end;%for ndelta=0:n_delta-1;
if (flag_verbose>0); disp(sprintf(' %% tmp_I_ vs tmp_J_: %0.16f %%<-- should be <1e-6 ',fnorm(tmp_I_-tmp_J_)/fnorm(tmp_I_))); end;

%%%%%%%%;
% integral of plane-wave times planar-function times plane-wave. ;
%%%%%%%%;
% C = -2i\vk \dot \vn = -2i k cos(\psi-\phi) ;
% S = exp(i 2\pi \vk \dot \vd_{S}) = exp(i 2\pi k \delta_{S} cos(\psi - \omega_{S})) ;
% M = exp(i 2\pi \vk \dot \vd_{M}) = exp(i 2\pi k \delta_{M} cos(\psi - \omega_{M})) ;
% \vd_{T} = \vd_{M} - \vd_{S} = \delta_{T} [ \cos(\omega_{T}) ; \sin(\omega_{T}) ] ;
% \int_{0}^{K}\int_{0}^{2\pi} [C\cdot S]^{\dagger} M d\psi kdk ... ;
%                            = - \pi K^2 Kcos(\phi-\omega_{T}) (J_{1}(2\pi K\delta_{T})+J_{3}(2\pi K\delta_{T})) ;
%%%%%%%%;
k_p_r_max = 48/(2*pi);
delta_max = 0.1;
tmp_delta_S = rand()*delta_max;
tmp_omega_S = 2*pi*rand();
tmp_delta_S_ = tmp_delta_S * [cos(tmp_omega_S) ; sin(tmp_omega_S)];
tmp_delta_M = rand()*delta_max;
tmp_omega_M = 2*pi*rand();
tmp_delta_M_ = tmp_delta_M * [cos(tmp_omega_M) ; sin(tmp_omega_M)];
tmp_delta_T_ = tmp_delta_M_ - tmp_delta_S_;
tmp_delta_T = fnorm(tmp_delta_T_);
tmp_omega_T = atan2(tmp_delta_T_(1+1),tmp_delta_T_(1+0));
tmp_phi = 2*pi*rand();
tmp_I = integral2( ...
@(k_p_r,psi) k_p_r ...
 .* conj((-2*i*k_p_r.*cos(psi-tmp_phi)).*exp(2*pi*i*k_p_r.*tmp_delta_S.*cos(psi-tmp_omega_S))) ...
 .* exp(2*pi*i*k_p_r.*tmp_delta_M.*cos(psi-tmp_omega_M)) ...
,0,k_p_r_max,0,2*pi);
tmp_J = -pi*k_p_r_max^2*k_p_r_max*cos(tmp_phi-tmp_omega_T)*(besselj(1,2*pi*k_p_r_max*tmp_delta_T) + besselj(3,2*pi*k_p_r_max*tmp_delta_T));
if (flag_verbose>0); disp(sprintf(' %% tmp_I vs tmp_J: %0.16f %%<-- should be <1e-6 ',fnorm(tmp_I-tmp_J)/fnorm(tmp_I))); end;

%%%%%%%%;
% integral of plane-wave times planar-function times plane-wave times planar-function. ;
%%%%%%%%;
% C_{S} = 2\vk \dot \vn_{S} = 2 k cos(\psi-\phi_{S}) ;
% S = C_{S} \cdot exp(i 2\pi \vk \dot \vd_{S}) = 2 k cos(\psi-\phi_{S}) \cdot exp(i 2\pi k \delta_{S} cos(\psi - \omega_{S})) ;
% C_{M} = -2\vk \dot \vn_{M} = -2i k cos(\psi-\phi_{M}) ;
% M = exp(i 2\pi \vk \dot \vd_{M}) = exp(i 2\pi k \delta_{M} cos(\psi - \omega_{M})) ;
% \vd_{T} = \vd_{M} - \vd_{S} = \delta_{T} [ \cos(\omega_{T}) ; \sin(\omega_{T}) ] ;
% phi_{+} = \phi_{S} + \phi_{M} ;
% phi_{-} = \phi_{S} - \phi_{M} ;
% \int_{0}^{K}\int_{0}^{2\pi} [C\cdot S]^{\dagger} M d\psi kdk ... ;
%                            = 4\pi cos(\phi_{-}) K^4 (J0/4 + J2/6 - J4/12) - 4\pi cos(\phi_{+} - 2\omega_{T}) K^4 (J2/6 + J4/6) ;
% where each J is evaluated at (2\pi K\delta_{T}) ;
%%%%%%%%;
k_p_r_max = 48/(2*pi);
delta_max = 0.1;
tmp_delta_S = rand()*delta_max;
tmp_omega_S = 2*pi*rand();
tmp_delta_S_ = tmp_delta_S * [cos(tmp_omega_S) ; sin(tmp_omega_S)];
tmp_delta_M = rand()*delta_max;
tmp_omega_M = 2*pi*rand();
tmp_delta_M_ = tmp_delta_M * [cos(tmp_omega_M) ; sin(tmp_omega_M)];
tmp_delta_T_ = tmp_delta_M_ - tmp_delta_S_;
tmp_delta_T = fnorm(tmp_delta_T_);
tmp_omega_T = atan2(tmp_delta_T_(1+1),tmp_delta_T_(1+0));
tmp_phi_S = 2*pi*rand();
tmp_phi_M = 2*pi*rand();
tmp_phi_pos = tmp_phi_S + tmp_phi_M;
tmp_phi_neg = tmp_phi_S - tmp_phi_M;
tmp_I = integral2( ...
@(k_p_r,psi) k_p_r ...
 .* conj((2*k_p_r.*cos(psi-tmp_phi_M)).*(2*k_p_r.*cos(psi-tmp_phi_S)).*exp(2*pi*i*k_p_r.*tmp_delta_S.*cos(psi-tmp_omega_S))) ...
 .* exp(2*pi*i*k_p_r.*tmp_delta_M.*cos(psi-tmp_omega_M)) ...
,0,k_p_r_max,0,2*pi);
tmp_kd = 2*pi*k_p_r_max*tmp_delta_T;
tmp_J = ...
+4*pi*cos(tmp_phi_neg)*k_p_r_max^4*(besselj(0,tmp_kd)/4 + besselj(2,tmp_kd)/6 - besselj(4,tmp_kd)/12) ...
-4*pi*cos(tmp_phi_pos-2*tmp_omega_T)*k_p_r_max^4*(besselj(2,tmp_kd)/6 + besselj(4,tmp_kd)/6) ...
;
if (flag_verbose>0); disp(sprintf(' %% tmp_I vs tmp_J: %0.16f %%<-- should be <1e-6 ',fnorm(tmp_I-tmp_J)/fnorm(tmp_I))); end;

%%%%%%%%;
% Note that, with this notation: ;
% The norm of (C_{M}\otimes R_{\gamma}\circ S) depends on the rotation-angle \gamma. ;
% Note also that R_{-\gamma}\circ C_{M} = 2*k_p_r*cos(\psi+\gamma - tmp_phi_M). ;
% I.e.: ;
% \int_{0}^{K}\int_{0}^{2\pi} [C\cdot S]^{\dagger} [C\cdot S] d\psi kdk ... ;
%                            = 2\pi K^6 / 3 (1 + 2\cos(\phi_{-}+\gamma)^2) ;
%%%%%%%%;
k_p_r_max = 48/(2*pi);
delta_max = 0.1;
tmp_delta_S = rand()*delta_max;
tmp_omega_S = 2*pi*rand();
tmp_delta_S_ = tmp_delta_S * [cos(tmp_omega_S) ; sin(tmp_omega_S)];
tmp_delta_M = rand()*delta_max;
tmp_omega_M = 2*pi*rand();
tmp_delta_M_ = tmp_delta_M * [cos(tmp_omega_M) ; sin(tmp_omega_M)];
tmp_gamma = 2*pi*rand();
tmp_phi_S = 2*pi*rand();
tmp_phi_M = 2*pi*rand();
tmp_phi_neg = tmp_phi_S - tmp_phi_M;
tmp_I = integral2( ...
@(k_p_r,psi) k_p_r ...
 .* conj((2*k_p_r.*cos(psi+tmp_gamma-tmp_phi_M)).*(2*k_p_r.*cos(psi-tmp_phi_S)).*exp(2*pi*i*k_p_r.*tmp_delta_S.*cos(psi-tmp_omega_S))) ...
 .*     ((2*k_p_r.*cos(psi+tmp_gamma-tmp_phi_M)).*(2*k_p_r.*cos(psi-tmp_phi_S)).*exp(2*pi*i*k_p_r.*tmp_delta_S.*cos(psi-tmp_omega_S))) ...
,0,k_p_r_max,0,2*pi);
tmp_J = 2*pi*k_p_r_max^6 * (1/3)*(1 + 2*cos(tmp_phi_neg + tmp_gamma)*cos(tmp_phi_neg + tmp_gamma)) ;
if (flag_verbose>0); disp(sprintf(' %% tmp_I vs tmp_J: %0.16f %%<-- should be <1e-6 ',fnorm(tmp_I-tmp_J)/fnorm(tmp_I))); end;

%%%%%%%%;
% integral of plane-wave times bessel times plane-wave. ;
% (i.e., isotropic CTF). ;
%%%%%%%%;
% C = J_{0}(\alpha k) ;
% S = exp(i 2\pi \vk \dot \vd_{S}) = exp(i 2\pi k \delta_{S} cos(\psi - \omega_{S})) ;
% M = exp(i 2\pi \vk \dot \vd_{M}) = exp(i 2\pi k \delta_{M} cos(\psi - \omega_{M})) ;
% \vd_{T} = \vd_{M} - \vd_{S} = \delta_{T} [ \cos(\omega_{T}) ; \sin(\omega_{T}) ] ;
% a = \alpha K ;
% b = 2\pi K \delta_{T} ;
% \int_{0}^{K}\int_{0}^{2\pi} [C\cdot S]^{\dagger} M d\psi kdk ... ;
%                            = 2\pi K^2 \frac{ bJ_{-1}(b)J_{0}(a) - aJ_{-1}(a)J_{0}(b) }{ a^2 - b^2 } ;
%%%%%%%%;
k_p_r_max = 48/(2*pi);
delta_max = 0.1;
tmp_delta_S = rand()*delta_max;
tmp_omega_S = 2*pi*rand();
tmp_delta_S_ = tmp_delta_S * [cos(tmp_omega_S) ; sin(tmp_omega_S)];
tmp_delta_M = rand()*delta_max;
tmp_omega_M = 2*pi*rand();
tmp_delta_M_ = tmp_delta_M * [cos(tmp_omega_M) ; sin(tmp_omega_M)];
tmp_delta_T_ = tmp_delta_M_ - tmp_delta_S_;
tmp_delta_T = fnorm(tmp_delta_T_);
tmp_omega_T = atan2(tmp_delta_T_(1+1),tmp_delta_T_(1+0));
tmp_alpha = rand();
tmp_a = tmp_alpha*k_p_r_max;
tmp_b = 2*pi*k_p_r_max*tmp_delta_T;
tmp_I = integral2( ...
@(k_p_r,psi) k_p_r ...
 .* conj(besselj(0,tmp_alpha*k_p_r).*exp(2*pi*i*k_p_r.*tmp_delta_S.*cos(psi-tmp_omega_S))) ...
 .* exp(2*pi*i*k_p_r.*tmp_delta_M.*cos(psi-tmp_omega_M)) ...
,0,k_p_r_max,0,2*pi);
tmp_c = max(tmp_a,tmp_b); tmp_d = min(tmp_a,tmp_b);
tmp_J = 2*pi*k_p_r_max^2 * (tmp_d*besselj(-1,tmp_d)*besselj(0,tmp_c) - tmp_c*besselj(-1,tmp_c)*besselj(0,tmp_d))/max(1e-12,tmp_c^2-tmp_d^2);
if (flag_verbose>0); disp(sprintf(' %% tmp_I vs tmp_J: %0.16f %%<-- should be <1e-6 ',fnorm(tmp_I-tmp_J)/fnorm(tmp_I))); end;

%%%%%%%%;
% Now set up polar-quadrature-weights. ;
%%%%%%%%;
k_p_r_max = 48/(2*pi); k_eq_d = 1.0/(2*pi); TorL = 'L';
[ ...
 n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
] = ...
get_weight_3d_1( ...
 flag_verbose ...
,k_p_r_max ...
,k_eq_d ...
,TorL ...
);
%%%%;
l_max_upb = round(2*pi*k_p_r_max);
l_max_max = min(l_max_upb,1+ceil(2*pi*k_p_r_(end)));
n_w_max = 2*(l_max_max+1); n_w_0in_ = n_w_max*ones(n_k_p_r,1);
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
 flag_verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,-1 ...
,n_w_0in_ ...
);
%%%%%%%%;
n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_]);

%%%%%%%%;
% Now use polar-quadrature to estimate integral of plane-wave times bessel times plane-wave. ;
%%%%%%%%;
delta_max = 0.1;
tmp_delta_S = rand()*delta_max;
tmp_omega_S = 2*pi*rand();
tmp_delta_S_ = tmp_delta_S * [cos(tmp_omega_S) ; sin(tmp_omega_S)];
tmp_delta_M = rand()*delta_max;
tmp_omega_M = 2*pi*rand();
tmp_delta_M_ = tmp_delta_M * [cos(tmp_omega_M) ; sin(tmp_omega_M)];
tmp_delta_T_ = tmp_delta_M_ - tmp_delta_S_;
tmp_delta_T = fnorm(tmp_delta_T_);
tmp_omega_T = atan2(tmp_delta_T_(1+1),tmp_delta_T_(1+0));
tmp_alpha = rand();
tmp_a = tmp_alpha*k_p_r_max;
tmp_b = 2*pi*k_p_r_max*tmp_delta_T;
tmp_I = integral2( ...
@(k_p_r,psi) k_p_r ...
 .* conj(besselj(0,tmp_alpha*k_p_r).*exp(2*pi*i*k_p_r.*tmp_delta_S.*cos(psi-tmp_omega_S))) ...
 .* exp(2*pi*i*k_p_r.*tmp_delta_M.*cos(psi-tmp_omega_M)) ...
,0,k_p_r_max,0,2*pi);
tmp_c = max(tmp_a,tmp_b); tmp_d = min(tmp_a,tmp_b);
tmp_J = 2*pi*k_p_r_max^2 * (tmp_d*besselj(-1,tmp_d)*besselj(0,tmp_c) - tmp_c*besselj(-1,tmp_c)*besselj(0,tmp_d))/max(1e-12,tmp_c^2-tmp_d^2);
CTF_wk_ = besselj(0,tmp_alpha*k_p_r_wk_);
S_wk_ = exp(2*pi*i*k_p_r_wk_.*tmp_delta_S.*cos(k_p_w_wk_-tmp_omega_S));
M_wk_ = exp(2*pi*i*k_p_r_wk_.*tmp_delta_M.*cos(k_p_w_wk_-tmp_omega_M));
tmp_Q = sum(conj(CTF_wk_.*S_wk_).*M_wk_.*weight_2d_wk_,'all')*(2*pi)^2;
if (flag_verbose>0); disp(sprintf(' %% tmp_I vs tmp_J: %0.16f %%<-- should be <1e-6 ',fnorm(tmp_I-tmp_J)/fnorm(tmp_I))); end;
if (flag_verbose>0); disp(sprintf(' %% tmp_Q vs tmp_J: %0.16f %%<-- should be <1e-2 ',fnorm(tmp_Q-tmp_J)/fnorm(tmp_Q))); end;

flag_CTF_use = 1;
flag_delta_use = 1;
n_S = 7; n_M = 5;
if (flag_verbose); disp(sprintf(' %% testing ampmh_X_wSM___8.m using n_S %d templates and n_M %d images. ',n_S,n_M)); end;
delta_max = 0.1;
delta_S_S_ = rand(n_S,1)*delta_max;
omega_S_S_ = 2*pi*rand(n_S,1);
delta_S_2S__ = zeros(2,n_S);
S_k_p_wkS__ = zeros(n_w_sum,n_S);
S_k_q_wkS__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
tmp_delta_S = delta_S_S_(1+nS);
tmp_omega_S = omega_S_S_(1+nS);
tmp_delta_S_ = tmp_delta_S * [cos(tmp_omega_S) ; sin(tmp_omega_S)];
delta_S_2S__(:,1+nS) = tmp_delta_S_;
S_k_p_wkS__(:,1+nS) = exp(2*pi*i*k_p_r_wk_.*tmp_delta_S.*cos(k_p_w_wk_-tmp_omega_S));;
S_k_q_wkS__(:,1+nS) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,S_k_p_wkS__(:,1+nS));
end;%for nS=0:n_S-1;
delta_M_M_ = rand(n_M,1)*delta_max;
omega_M_M_ = 2*pi*rand(n_M,1);
delta_M_M_(1+0) = delta_S_S_(1+0); omega_M_M_(1+0) = omega_S_S_(1+0); %<-- ensure that the first image matches the first template. ;
delta_M_2M__ = zeros(2,n_M);
M_k_p_wkM__ = zeros(n_w_sum,n_M);
M_k_q_wkM__ = zeros(n_w_sum,n_M);
for nM=0:n_M-1;
tmp_delta_M = delta_M_M_(1+nM);
tmp_omega_M = omega_M_M_(1+nM);
tmp_delta_M_ = tmp_delta_M * [cos(tmp_omega_M) ; sin(tmp_omega_M)];
delta_M_2M__(:,1+nM) = tmp_delta_M_;
M_k_p_wkM__(:,1+nM) = exp(2*pi*i*k_p_r_wk_.*tmp_delta_M.*cos(k_p_w_wk_-tmp_omega_M));;
M_k_q_wkM__(:,1+nM) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_k_p_wkM__(:,1+nM));
end;%for nM=0:n_M-1;
CTF_alpha = flag_CTF_use*1.5;
CTF_k_p_r_k_ = besselj(0,CTF_alpha*k_p_r_);
CTF_k_p_r_wk_ = besselj(0,CTF_alpha*k_p_r_wk_); CTF_wk_ = CTF_k_p_r_wk_;
%%%%%%%%;
% Prepare FTK. ;
%%%%%%%%;
delta_r_max = flag_delta_use*0.05; n_delta_v_requested = flag_delta_use*16;
delta_r_p = 0.05;
delta_r_s = delta_r_max/sqrt(2*log(1/delta_r_p));
delta_r_N = delta_r_max * (2*pi*k_p_r_max) / (pi*sqrt(2));
svd_eps = 1e-3;
tmp_t = tic();
FTK = ampmh_FTK_1(n_k_p_r,k_p_r_,k_p_r_max,delta_r_max,svd_eps,n_delta_v_requested);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% FTK: %0.3fs',tmp_t)); end;
disp(sprintf(' %% p-val %0.4f delta_r_max %0.6f sigma %0.4f N_pixel %0.4f --> FTK.n_svd_l %d, n_delta_v_requested %d',delta_r_p,delta_r_max,delta_r_s,delta_r_N,FTK.n_svd_l,n_delta_v_requested));
%%%%%%%%;
% Prepare pricipal-modes. ;
%%%%%%%%;
n_UX_rank = n_k_p_r-1; %<-- just to check dimensions.; 
UX_kn__ = zeros(n_k_p_r,n_UX_rank);
X_weight_r_ = zeros(n_k_p_r,1);
[ ...
 X_kk__ ...
,X_weight_r_ ...
] = ...
principled_marching_empirical_cost_matrix_0( ...
 n_k_p_r ...
,k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_ ...
,n_M ...
,M_k_p_wkM__ ...
);
[tmp_UX_kn__,tmp_SX_k__,tmp_VX_kn__] = svds(X_kk__,n_UX_rank);
UX_kn__(:,1:n_UX_rank) = tmp_UX_kn__(:,1:n_UX_rank);
pm_n_UX_rank = n_UX_rank; %<-- just to check dimension. ;
%%%%%%%%;
% Prepare quasi-images. ;
%%%%%%%%;
n_UX_rank = n_k_p_r;
tmp_t = tic();
svd_VUXM_lwnM____ = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,n_M,M_k_q_wkM__,pm_n_UX_rank,UX_kn__,X_weight_r_);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% svd_VUXM_lwnM____: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Now calculate norms of the non-translated images. ;
%%%%%%%%;
tmp_t = tic();
UX_M_l2_dM__ = ampmh_UX_M_l2_dM__1(FTK,n_w_,n_M,pm_n_UX_rank,svd_VUXM_lwnM____);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% UX_M_l2_dM__: %0.3fs',tmp_t)); end;
disp(sprintf(' %% average l2-norm of templates: %0.16f %%<-- should be 1 ',mean(UX_M_l2_dM__(:))/(pi*k_p_r_max^2)));
%%%%%%%%;
% Prepare CTF_UX_S_k_q_wnS__. ;
%%%%%%%%;
S_k_q_wSk___ = permute(reshape(S_k_q_wkS__,[n_w_max,n_k_p_r,n_S]),[1,3,2]);
CTF_UX_S_k_q_wnS__ = reshape(permute(reshape(reshape(S_k_q_wSk___,[n_w_max*n_S,n_k_p_r])*diag(X_weight_r_.*CTF_k_p_r_k_)*UX_kn__,[n_w_max,n_S,pm_n_UX_rank]),[1,3,2]),[n_w_max*pm_n_UX_rank,n_S]);
CTF_UX_S_l2_S_ = sum(abs(CTF_UX_S_k_q_wnS__).^2,1)/n_w_max;
%%%%%%%%;
% Calculate ampmh_X_wSM___8. ;
%%%%%%%%;
parameter = struct('type','parameter');
parameter.flag_compute_I_value = 0;
tmp_t = tic();
[ ...
 parameter ...
,X_wSM_ampm___ ...
,delta_x_wSM___ ...
,delta_y_wSM___ ...
,gamma_z_wSM___ ...
,I_value_wSM___ ...
] = ...
ampmh_X_wSM___8( ...
 parameter ...
,FTK ...
,n_w_max ...
,pm_n_UX_rank ...
,n_S ...
,CTF_UX_S_k_q_wnS__ ...
,CTF_UX_S_l2_S_ ...
,n_M ...
,svd_VUXM_lwnM____ ...
,UX_M_l2_dM__ ...
);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% X_wSM___: %0.3fs',tmp_t)); end;
X_wSM_ampm___ = real(X_wSM_ampm___);
%%%%%%%%;
% Calculate true landscape of innerproducts for the same set of translations. ;
%%%%%%%%;
X_wSM_bess___ = zeros(n_w_max,n_S,n_M);
X_wSM_quad___ = zeros(n_w_max,n_S,n_M);
tmp_a = CTF_alpha*k_p_r_max;
for nS=0:n_S-1;
if (flag_verbose>0); disp(sprintf(' %% nS %d/%d',nS,n_S)); end;
tmp_delta_S = delta_S_S_(1+nS);
tmp_omega_S = omega_S_S_(1+nS);
tmp_delta_S_ = tmp_delta_S * [cos(tmp_omega_S) ; sin(tmp_omega_S)];
for nM=0:n_M-1;
tmp_delta_M = delta_M_M_(1+nM);
tmp_omega_M = omega_M_M_(1+nM);
tmp_delta_M_ = tmp_delta_M * [cos(tmp_omega_M) ; sin(tmp_omega_M)];
for nw=0:n_w_max-1;
gamma_z = (2*pi*nw)/n_w_max; assert(gamma_z==gamma_z_wSM___(1+nw,1+nS,1+nM));
cc = cos(+gamma_z); sc = sin(+gamma_z);
Rz = [+cc , -sc ; +sc , +cc];
delta_x = delta_x_wSM___(1+nw,1+nS,1+nM);
delta_y = delta_y_wSM___(1+nw,1+nS,1+nM);
tmp_delta_d_ = [delta_x;delta_y];
tmp_delta_d = fnorm(tmp_delta_d_);
tmp_omega_d = atan2(delta_y,delta_x);
tmp_delta_TM_ = tmp_delta_M_ - tmp_delta_d_; %<-- multiply image by exp(-2*pi*i*dot(k_,delta_)). ;
tmp_delta_TM = fnorm(tmp_delta_TM_);
tmp_omega_TM = atan2(tmp_delta_TM_(1+1),tmp_delta_TM_(1+0));
tmp_delta_RS_ = Rz*tmp_delta_S_; %<-- rotate template by +gamma_z. ;
tmp_delta_RS = fnorm(tmp_delta_RS_);
tmp_omega_RS = atan2(tmp_delta_RS_(1+1),tmp_delta_RS_(1+0));
%%%%;
%{
tmp_I = integral2( ...
@(k_p_r,psi) k_p_r ...
 .* conj(besselj(0,CTF_alpha*k_p_r).*exp(2*pi*i*k_p_r.*tmp_delta_RS.*cos(psi-tmp_omega_RS))) ...
 .* exp(2*pi*i*k_p_r.*tmp_delta_TM.*cos(psi-tmp_omega_TM)) ...
,0,k_p_r_max,0,2*pi);
%}
%%%%;
RS_wk_ = exp(2*pi*i*k_p_r_wk_.*tmp_delta_RS.*cos(k_p_w_wk_-tmp_omega_RS));
TM_wk_ = exp(2*pi*i*k_p_r_wk_.*tmp_delta_TM.*cos(k_p_w_wk_-tmp_omega_TM));
tmp_Q = sum(conj(CTF_wk_.*RS_wk_).*TM_wk_.*weight_2d_wk_,'all')*(2*pi)^2;
tmp_RS_l2 = sum(conj(CTF_wk_.*RS_wk_).*(CTF_wk_.*RS_wk_).*weight_2d_wk_,'all')*(2*pi)^2;
tmp_TM_l2 = sum(conj(TM_wk_).*(TM_wk_).*weight_2d_wk_,'all')*(2*pi)^2;
%%%%;
tmp_delta_RSTM_ = tmp_delta_TM_ - tmp_delta_RS_;
tmp_delta_RSTM = fnorm(tmp_delta_RSTM_);
tmp_b = 2*pi*k_p_r_max*tmp_delta_RSTM;
tmp_c = max(tmp_a,tmp_b); tmp_d = min(tmp_a,tmp_b);
tmp_J = sqrt(tmp_RS_l2*tmp_TM_l2); %<-- default if formula fails. ;
if (abs(tmp_c-tmp_d)>1e-12);
tmp_J = 2*pi*k_p_r_max^2 * (tmp_d*besselj(-1,tmp_d)*besselj(0,tmp_c) - tmp_c*besselj(-1,tmp_c)*besselj(0,tmp_d))/max(1e-12,tmp_c^2-tmp_d^2);
end;%if (abs(tmp_c-tmp_d)>1e-12);
%%%%;
X_wSM_bess___(1+nw,1+nS,1+nM) = tmp_J/max(1e-12,sqrt(tmp_RS_l2*tmp_TM_l2));
X_wSM_quad___(1+nw,1+nS,1+nM) = tmp_Q/max(1e-12,sqrt(tmp_RS_l2*tmp_TM_l2));
end;%for nw=0:n_w_max-1;
end;%for nM=0:n_M-1;
end;%for nS=0:n_S-1;
if (flag_verbose>0); disp(sprintf(' %% X_wSM_bess___ vs X_wSM_quad___: %0.16f %%<-- should be <1e-2 ',fnorm(X_wSM_bess___-X_wSM_quad___)/fnorm(X_wSM_bess___))); end;
if (flag_verbose>0); disp(sprintf(' %% X_wSM_bess___ vs X_wSM_ampm___: %0.16f %%<-- should be <1e-2 ',fnorm(X_wSM_bess___-X_wSM_ampm___)/fnorm(X_wSM_bess___))); end;

%plot(real(X_wSM_bess___(:)),real(X_wSM_ampm___(:)),'.');
