%%%%%%%%;
% Paraxial approximation to shrodinger-equation. ;
% Beam headed in z-direction. ;
% Using polar (i.e., radial) description of x-y-plane for now. ;
%%%%%%%%;

flag_verbose=1; nf=0;
flag_disp = 1; nf=0; flag_replot = 1;

%%%%%%%%;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
%%%%%%%%;
dir_base = sprintf('/%s/rangan/dir_cryoem/dir_rangan_besselbeam/dir_jpg',string_root);

%%%%%%%%;
% Simple test of J0 integral. ;
%%%%%%%%;
kx = 1.75;
g = @(w_) exp(-i*2*pi*kx*cos(w_));
I_quad = integral(g,0,2*pi);
I_form = besselj(0,2*pi*kx)*(2*pi);
fnorm_disp(flag_verbose,'I_form',I_form,'I_quad',I_quad,' %<-- should be zero.');
%%%%%%%%;

%%%%%%%%;
% Simple test of Jq integral (frwd fft). ;
%%%%%%%%;
for nq=0:8-1;
q = nq; omega = pi/3;
kx = 1.75;
g = @(w_) exp(-i*2*pi*kx*cos(omega - w_)) .* exp(-i*q*w_) ;
I_quad = integral(g,0,2*pi);
I_form = (-i).^q .* exp(-i*q*omega) .* besselj(q,2*pi*kx)*(2*pi);
fnorm_disp(flag_verbose,'I_form',I_form,'I_quad',I_quad,sprintf(' %%<-- q %d: error should be zero.',q));
end;%for nq=0:8-1;
%%%%%%%%;

%%%%%%%%;
% Simple test of Jq integral (back fft). ;
%%%%%%%%;
for nq=0:8-1;
q = nq; omega = pi/3;
kx = 1.75;
g = @(w_) exp(-i*2*pi*kx*cos(omega - w_)) .* exp(+i*q*w_) ;
I_quad = integral(g,0,2*pi);
I_form = (-i).^q .* exp(+i*q*omega) .* besselj(q,2*pi*kx)*(2*pi);
fnorm_disp(flag_verbose,'I_form',I_form,'I_quad',I_quad,sprintf(' %%<-- q %d: error should be zero.',q));
end;%for nq=0:8-1;
%%%%%%%%;

%%%%%%%%;
% Here we test a version of the integral in Gradshteyn and Ryzhik: 2.33 (just to make sure there are no typos). ;
%%%%%%%%;
a = 3.0;
b = -0.5;
g = @(x_) exp(-a*x_.^2 + b*x_);
I_quad = integral(g,-10,+10);
I_form = sqrt(pi/max(1e-12,a)) * exp(b^2/max(1e-12,4*a));
fnorm_disp(flag_verbose,'I_form',I_form,'I_quad',I_quad,' %<-- should be zero if integral is resolved.');
%%%%%%%%;

%%%%%%%%;
% Here we test another integral. ;
%%%%%%%%;
eta = 2.5; k = 0.1; x_min = -10*4; x_max = +10*4;
g0 = @(x_) exp(-i*(eta.^2*x_.^2 + 2*pi*k*x_));
%g0 = @(x_) exp(-i*(-eta.^2*x_.^2 + 2*pi*k*x_));
I_quad = integral(g0,x_min,x_max);
g1 = @(x_) exp(+i*(pi*k/eta)^2) * exp(-i*(eta*x_ + pi*k/eta).^2);
I_qua1 = integral(g1,x_min,x_max);
fnorm_disp(flag_verbose,'I_qua1',I_qua1,'I_quad',I_quad,' %<-- should be zero (change of variables).');
g2 = @(y_) exp(+i*(pi*k/eta)^2) * exp(-i*y_.^2);
y_min = x_min*eta - (pi*k/eta); y_max = x_max*eta - (pi*k/eta);
I_qua2 = integral(g2,y_min,y_max) / eta;
fnorm_disp(flag_verbose,'I_qua2',I_qua2,'I_quad',I_quad,' %<-- should be zero (change of variables).');
I_form = exp(-i*pi/4)*sqrt(pi)/eta * exp(+i*(pi*k/eta)^2); %<-- principal-value of oscillatory integral. ;
%I_form = i*exp(-i*pi/4)*sqrt(pi)/eta * exp(-i*(pi*k/eta)^2); %<-- principal-value of oscillatory integral. ;
fnorm_disp(flag_verbose,'I_form',I_form,'I_quad',I_quad,' %<-- should be zero if integral is resolved.');
%x_ = linspace(x_min,x_max,1024*8); plot(x_,g0(x_),'o-');
%%%%%%%%;

%%%%%%%%;
% Consequently, the (full-space) integral of: ;
% exp(-i*(eta.^2*x_0_.^2 + 2*pi*k_0*x_0_)) .* exp(-i*(eta.^2*x_1_.^2 + 2*pi*k_1*x_1_)) dx_0 dx_1 ;
% --> exp(-i*pi/4)*sqrt(pi)/eta * exp(+i*(pi*k_0/eta)^2) .* exp(-i*pi/4)*sqrt(pi)/eta * exp(+i*(pi*k_1/eta)^2) ;
% --> exp(-i*pi/2) * pi/eta^2 * exp( +i * (pi/eta)^2 * (k_0^2 + k_1^2) ) ;
%%%%%%%%;
eta = 2.50; k = 0.1; x_min = 0; x_max = 19;
h0 = @(x_) x_ .* exp(-i*eta^2*x_.^2) .* besselj(0,2*pi*k*x_) * (2*pi);
I_quad = integral(h0,x_min,x_max);
%x_ = linspace(x_min,x_max,1024*8); plot(x_,h0(x_),'o-');
h1 = @(x_) +1/(2*-i*eta^2) * exp(-i*eta^2*x_.^2) .* besselj(0,2*pi*k*x_) * (2*pi);
h2 = @(x_) -1/(2*-i*eta^2) * exp(-i*eta^2*x_.^2) .* besselj(1,2*pi*k*x_) * (2*pi*k)*(2*pi);
I_qua2 = h1(x_max) - h1(x_min) - integral(h2,x_min,x_max);
%x_ = linspace(x_min,x_max,1024*8); plot(x_,h1(x_),'o-');
fnorm_disp(flag_verbose,'I_qua2',I_qua2,'I_quad',I_quad,' %<-- should be zero (integration by parts).');
I_form = exp(-i*pi/2) * pi/eta^2 * exp( +i * (pi/eta)^2 * (k^2) ) ;
fnorm_disp(flag_verbose,'I_form',I_form,'I_quad',I_quad,' %<-- should be zero if integral is resolved.');
%%%%%%%%;

%%%%%%%%;
% Set up spatial-grid for testing. ;
%%%%%%%%;
half_diameter_x_c = 1.0d0;
diameter_x_c = 2.0d0*half_diameter_x_c;
%%%%%%%%;
n_x_M_u = 128;
x_c_0_ = -half_diameter_x_c + transpose([0:n_x_M_u-1]/n_x_M_u)*diameter_x_c;
x_c_1_ = -half_diameter_x_c + transpose([0:n_x_M_u-1]/n_x_M_u)*diameter_x_c;
x_c_0_ = transpose(linspace(-half_diameter_x_c,+half_diameter_x_c,n_x_M_u+1)); x_c_0_ = x_c_0_(1:n_x_M_u);
x_c_1_ = transpose(linspace(-half_diameter_x_c,+half_diameter_x_c,n_x_M_u+1)); x_c_1_ = x_c_1_(1:n_x_M_u);
%x_c_0_ = transpose(linspace(-half_diameter_x_c,+half_diameter_x_c,n_x_M_u));
%x_c_1_ = transpose(linspace(-half_diameter_x_c,+half_diameter_x_c,n_x_M_u));
[x_c_0__,x_c_1__] = ndgrid(x_c_0_,x_c_1_);
dx = mean(diff(x_c_0_));
n_xx_M_u = n_x_M_u^2;

%%%%%%%%;
% Now set up polar grids in real-space and fourier-space. ;
%%%%%%%%;
x_res = 2;
x_int = 48*x_res;
x_eq_d_double = 0.5;
n_t_int = 1;
%%%%%%%%;
k_res = 2; %<-- be careful with the aliasing associated with n_x_M_u. ;
k_int = 48*k_res;
k_eq_d_double = 0.5;
n_w_int = n_t_int;
%%%%%%%%;
x_p_r_max = half_diameter_x_c;
k_p_r_max = k_int/(2*pi);
x_eq_d = x_eq_d_double/(2*pi)/k_p_r_max * 1.05; %<-- to check dimensions. ;
k_eq_d = k_eq_d_double/(2*pi)/x_p_r_max;
%%%%%%%%;

%%%%%%%%;
[ ...
 n_x_p_r ...
,x_p_r_ ...
,weight_3d_x_p_r_ ...
] = ...
get_weight_3d_1( ...
 max(0,flag_verbose-1) ...
,x_p_r_max ...
,x_eq_d ...
);
fnorm_disp(flag_verbose,'sum(weight_3d_x_p_r_)*4*pi',sum(weight_3d_x_p_r_)*4*pi,'volume',4/3*pi*x_p_r_max^3,' %<-- should be small');
%%%%%%%%;
n_t_max = n_t_int*2*(x_int+1) + 2;  %<-- to check dimensions. ;
n_t_0in_ = n_t_max*ones(n_x_p_r,1);
%%%%%%%%;
[ ...
 n_t_ ...
,weight_2d_x_p_r_ ...
,weight_2d_x_p_tx_ ...
,x_p_r_tx_ ...
,x_p_t_tx_ ...
,x_c_0_tx_ ...
,x_c_1_tx_ ...
] = ...
get_weight_2d_2( ...
 max(0,flag_verbose-1) ...
,n_x_p_r ...
,x_p_r_ ...
,x_p_r_max ...
,-1 ...
,n_t_0in_ ...
);
fnorm_disp(flag_verbose,'sum(weight_2d_x_p_r_)',sum(weight_2d_x_p_r_),'area',pi*x_p_r_max^2,' %<-- should be small');
n_t_sum = sum(n_t_); n_t_csum_ = cumsum([0;n_t_]);
%%%%%%%%;

%%%%%%%%;
[ ...
 n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
] = ...
get_weight_3d_1( ...
 max(0,flag_verbose-1) ...
,k_p_r_max ...
,k_eq_d ...
);
fnorm_disp(flag_verbose,'sum(weight_3d_k_p_r_)*4*pi',sum(weight_3d_k_p_r_)*4*pi,'volume',4/3*pi*k_p_r_max^3,' %<-- should be small');
%%%%%%%%;
%n_w_max = 12; 
n_w_max = n_w_int*2*(k_int+1);
n_w_0in_ = n_w_max*ones(n_k_p_r,1);
%%%%%%%%;
[ ...
 n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_k_p_wk_ ...
,k_p_r_wk_ ...
,k_p_w_wk_ ...
,k_c_0_wk_ ...
,k_c_1_wk_ ...
] = ...
get_weight_2d_2( ...
 max(0,flag_verbose-1) ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,-1 ...
,n_w_0in_ ...
);
fnorm_disp(flag_verbose,'sum(weight_2d_k_p_r_)',sum(weight_2d_k_p_r_),'area',pi*k_p_r_max^2,' %<-- should be small');
n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_]);
%%%%%%%%;

%%%%%%%%;
% Now test gaussian. ;
%%%%%%%%;
sigma_x_c = 0.0625;
sigma_k_p = 1/sigma_x_c;
delta_ = 3*0.75*[+0.1,-0.2];
M_x_c_form_ = 1/(sqrt(2*pi)*sigma_x_c)^2 * exp( -( (x_c_0__-delta_(1+0)).^2 + (x_c_1__-delta_(1+1)).^2 ) / (2*sigma_x_c^2) );
M_x_c_form_l2 = sum(M_x_c_form_.^2,'all')*dx^2;
disp(sprintf(' %% sum(M_x_c_form_*dx^2,''all'') = %0.16f',sum(M_x_c_form_*dx^2,'all')));
disp(sprintf(' %% M_x_c_form_l2 = %0.16f',M_x_c_form_l2));
M_x_p_form_ = 1/(sqrt(2*pi)*sigma_x_c)^2 * exp( -( (x_c_0_tx_-delta_(1+0)).^2 + (x_c_1_tx_-delta_(1+1)).^2 ) / (2*sigma_x_c^2) );
M_x_p_form_l2 = sum(abs(M_x_p_form_).^2 .* weight_2d_x_p_tx_) * (2*pi)^2;
disp(sprintf(' %% M_x_p_form_l2 = %0.16f',M_x_p_form_l2));
M_k_p_quad_ = interp_x_c_to_k_p_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,M_x_c_form_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_M_u^2)*dx^2;
M_k_p_quad_l2 = sum(abs(M_k_p_quad_).^2 .* weight_2d_k_p_wk_) * (2*pi)^2;
disp(sprintf(' %% M_k_p_quad_l2 = %0.16f',M_k_p_quad_l2));
M_k_p_form_ = exp( -( (2*pi*k_c_0_wk_).^2 + (2*pi*k_c_1_wk_).^2 ) / (2/sigma_x_c^2) ) .* exp( - 2*pi*i*( k_c_0_wk_*delta_(1+0) + k_c_1_wk_*delta_(1+1) ) );
M_k_p_form_l2 = sum(abs(M_k_p_form_).^2 .* weight_2d_k_p_wk_) * (2*pi)^2;
disp(sprintf(' %% M_k_p_form_l2 = %0.16f',M_k_p_form_l2));
fnorm_disp(flag_verbose,'M_k_p_form_',M_k_p_form_,'M_k_p_quad_',M_k_p_quad_,' %<-- should be small if no aliasing');
M_k_p_reco_ = interp_x_p_to_k_p_xxnufft(n_k_p_r,k_p_r_,k_p_r_max,n_w_,weight_2d_k_p_wk_,n_x_p_r,x_p_r_,x_p_r_max,n_t_,weight_2d_x_p_tx_,M_x_p_form_)*(2*pi)^2;
M_k_p_reco_l2 = sum(abs(M_k_p_reco_).^2 .* weight_2d_k_p_wk_) * (2*pi)^2;
disp(sprintf(' %% M_k_p_reco_l2 = %0.16f',M_k_p_reco_l2));
M_x_c_reco_ = interp_k_p_to_x_c_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,M_k_p_quad_.*weight_2d_k_p_wk_*(2*pi)^2)*sqrt(n_x_M_u^2) * n_w_sum;
M_x_c_reco_l2 = sum(abs(M_x_c_reco_).^2,'all')*dx^2;
disp(sprintf(' %% M_x_c_reco_l2 = %0.16f',M_x_c_reco_l2));
fnorm_disp(flag_verbose,'M_x_c_form_',M_x_c_form_,'M_x_c_reco_',M_x_c_reco_,' %<-- should be small if no aliasing');
M_x_p_reco_ = interp_k_p_to_x_p_xxnufft(n_x_p_r,x_p_r_,x_p_r_max,n_t_,weight_2d_x_p_tx_,n_k_p_r,k_p_r_,k_p_r_max,n_w_,weight_2d_k_p_wk_,M_k_p_form_)*(2*pi)^2;
M_x_p_reco_l2 = sum(abs(M_x_p_reco_).^2 .* weight_2d_x_p_tx_) * (2*pi)^2;
disp(sprintf(' %% M_x_p_reco_l2 = %0.16f',M_x_p_reco_l2));

%%%%%%%%;
% Now test fourier-bessel representation. ;
%%%%%%%%;
M_x_q_quad_ = interp_p_to_q(n_x_p_r,n_t_,n_t_sum,M_x_p_form_);
M_x_q_quad_l2 = sum(abs(M_x_q_quad_).^2 .* weight_2d_x_p_tx_) * (2*pi)^2;
disp(sprintf(' %% M_x_q_quad_l2 = %0.16f',M_x_q_quad_l2));
M_k_q_quad_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_k_p_form_);
M_k_q_quad_l2 = sum(abs(M_k_q_quad_).^2 .* weight_2d_k_p_wk_) * (2*pi)^2;
disp(sprintf(' %% M_k_q_quad_l2 = %0.16f',M_k_q_quad_l2));
%%%%%%%%;
%imagesc_q(n_x_p_r,x_p_r_,n_t_,n_t_sum,real(M_x_q_quad_));
%imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(M_k_q_quad_));
%%%%%%%%;
M_x_q_quad__ = reshape(M_x_q_quad_,[n_t_max,n_x_p_r]);
M_k_q_quad__ = reshape(M_k_q_quad_,[n_w_max,n_k_p_r]);
%%%%%%%%;
for q_val=-6:+6;
nt = q_val; if (nt<0); nt=nt+n_t_max; end;
nw = q_val; if (nw<0); nw=nw+n_w_max; end;
M_x_q_x_ = reshape(M_x_q_quad__(1+nt,:),[n_x_p_r,1]);
M_k_q_k_ = reshape(M_k_q_quad__(1+nw,:),[n_k_p_r,1]);
N_x_q_x_ = zeros(n_x_p_r,1);
for nx_p_r=0:n_x_p_r-1;
x_p_r = x_p_r_(1+nx_p_r);
N_x_q_x_(1+nx_p_r) = sum( (+i).^q_val .* besselj(q_val,2*pi*x_p_r*k_p_r_) .* M_k_q_k_ .* weight_2d_k_p_r_ ); %<-- note that weight_2d_k_p_r_ accounts for radial weighting. ;
end;%for nx_p_r=0:n_x_p_r-1;
fnorm_disp(flag_verbose,'M_x_q_x_',M_x_q_x_,'N_x_q_x_',N_x_q_x_,sprintf(' %%<-- q_val %+d: should be small',q_val));
N_k_q_k_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
N_k_q_k_(1+nk_p_r) = sum( (-i).^q_val .* besselj(q_val,2*pi*k_p_r*x_p_r_) .* M_x_q_x_ .* weight_2d_x_p_r_ ); %<-- note that weight_2d_x_p_r_ accounts for radial weighting. ;
end;%for nk_p_r=0:n_k_p_r-1;
fnorm_disp(flag_verbose,'M_k_q_k_',M_k_q_k_,'N_k_q_k_',N_k_q_k_,sprintf(' %%<-- q_val %+d: should be small',q_val));
end;%for q_val=-6:+6;
%%%%%%%%;

%%%%%%%%;
if flag_disp>0;
figure(1+nf);nf=nf+1;figbig;fig80s;
p_row = 2; p_col = 4; np=0;
M_x_p_lim_ = reshape(prctile(real(M_x_p_form_),[  0,100],'all'),[1,2]);
M_k_p_lim_ = reshape(prctile(real(M_k_p_form_),[  0,100],'all'),[1,2]);
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,real(M_x_c_form_),M_x_p_lim_);axis image;axisnotick; title('M_x_c_form_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_p(n_x_p_r,x_p_r_,n_t_,n_t_sum,real(M_x_p_form_),M_x_p_lim_);axis image;axisnotick; title('M_x_p_form_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_p(n_x_p_r,x_p_r_,n_t_,n_t_sum,real(M_x_p_reco_),M_x_p_lim_);axis image;axisnotick; title('M_x_p_reco_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(M_k_p_quad_),M_k_p_lim_);axis image;axisnotick; title('M_k_p_quad_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(M_k_p_form_),M_k_p_lim_);axis image;axisnotick; title('M_k_p_form_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(M_k_p_reco_),M_k_p_lim_);axis image;axisnotick; title('M_k_p_reco_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,real(M_x_c_reco_),M_x_p_lim_);axis image;axisnotick; title('M_x_c_reco_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,real(M_x_c_reco_-M_x_c_form_),M_x_p_lim_);axis image;axisnotick; title('difference','Interpreter','none');
end;%if flag_disp>0;
%%%%%%%%;

disp('returning'); return;

%%%%%%%%;
% simple plane-wave incident on lens. ;
%%%%%%%%;
eta = 2.50; %<-- kz/(2*f) = z-frequency / (2*focal_length). ;
%%%%%%%%;
str_type = 'ringlens';
%%%%;
if strcmp(str_type,'quadlens');
M_k_p_r_ = exp(-i*pi/2) * pi/max(1e-12,eta^2) * exp( + i * (pi/max(1e-12,eta))^2 * k_p_r_.^2 );
M_k_p_wk_ = reshape(transpose(repmat(M_k_p_r_,[1,n_w_max])),[n_w_sum,1]);
end;%if strcmp(str_type,'quadlens');
%%%%;
if strcmp(str_type,'ringlens');
ring_inn = 0.35; ring_ext = 0.36; %<-- interior and exterior radii of annulus. ;
%ring_inn = 0.00; ring_ext = 0.55; %<-- interior and exterior radii of annulus. ;
M_k_p_r_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
x_min = ring_inn; x_max = ring_ext;
h0 = @(x_) x_ .* exp(+i*eta^2*x_.^2) .* besselj(0,2*pi*k_p_r*x_) * (2*pi);
I_quad = integral(h0,x_min,x_max);
M_k_p_r_(1+nk_p_r) = I_quad;
end;%for nk_p_r=0:n_k_p_r-1;
M_k_p_wk_ = reshape(transpose(repmat(M_k_p_r_,[1,n_w_max])),[n_w_sum,1]);
end;%if strcmp(str_type,'quadlens');
%%%%;

%%%%%%%%;
% initial radial waveform. ;
%%%%%%%%;
flag_check=0;
if flag_check;
n_x_p_r = n_k_p_r; x_p_r_ = linspace(0,x_p_r_max,n_x_p_r);
M_x_p_r_ = zeros(n_x_p_r,1);
for nx_p_r=0:n_x_p_r-1;
x_p_r = x_p_r_(1+nx_p_r);
tmp_I = sum(2*pi*besselj(0,2*pi*k_p_r_*x_p_r).*M_k_p_r_.*k_p_r_.*weight_2d_k_p_r_);
M_x_p_r_(1+nx_p_r) = tmp_I;
end;%for nx_p_r=0:n_x_p_r-1;
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,2,1); cla; hold on;
plot(x_p_r_,abs(M_x_p_r_),'ko-');
plot(x_p_r_,real(M_x_p_r_),'r-');
plot(x_p_r_,imag(M_x_p_r_),'g-');
title('abs(M_x_p_r_)','Interpreter','none');
subplot(1,2,2); cla; hold on;
plot(k_p_r_,abs(M_k_p_r_),'ko-');
plot(k_p_r_,real(M_k_p_r_),'r-');
plot(k_p_r_,imag(M_k_p_r_),'g-');
title('abs(M_k_p_r_)','Interpreter','none');
disp('returning'); return;
end;%if flag_check;
%%%%%%%%;

%%%%%%%%;
% evolve radial waveform. ;
%%%%%%%%;
flag_check=1;
if flag_check;
%%%%;
n_x_p_r = 1024; x_p_r_ = linspace(0,x_p_r_max,n_x_p_r);
x_p_r_from_k_p_r_xk__ = zeros(n_x_p_r,n_k_p_r);
for nx_p_r=0:n_x_p_r-1;
x_p_r = x_p_r_(1+nx_p_r);
x_p_r_from_k_p_r_xk__(1+nx_p_r,:) = 2*pi*besselj(0,2*pi*k_p_r_*x_p_r).*k_p_r_.*weight_2d_k_p_r_;
end;%for nx_p_r=0:n_x_p_r-1;
%%%%;
lambda = 1.0; %<-- here measure z in terms of lambda. ;
dzl = 1.0/8192.0;
zl_max = 512*dzl;
n_zl = zl_max/dzl;
P_k_p_r_ = exp(-i*pi*lambda*dzl*k_p_r_.^2); %<-- radial fourier propagator ;
M_k_p_r_kt__ = zeros(n_k_p_r,n_zl);
M_k_p_r_kt__(:,1+0) = M_k_p_r_;
for nzl=1:n_zl-1;
M_k_p_r_kt__(:,1+nzl+0) = M_k_p_r_kt__(:,1+nzl-1).*P_k_p_r_;
end;%for nzl=1:n_zl-1;
%%%%;
M_x_p_r_xt__ = x_p_r_from_k_p_r_xk__*M_k_p_r_kt__;
%%%%;
figure(1+nf);nf=nf+1;clf;figmed;
subplot_{1} = subplot(2,1,1);
tmp_xt__ = abs(M_x_p_r_xt__);
tmp_xt__ = [flipud(tmp_xt__) ; tmp_xt__];
imagesc(tmp_xt__); set(gca,'ydir','normal'); axisnotick; title('abs(M_x_p_r_xt__)','Interpreter','none');
subplot_{2} = subplot(2,1,2);
tmp_xt__ = angle(M_x_p_r_xt__);
tmp_xt__ = [flipud(tmp_xt__) ; tmp_xt__];
imagesc(tmp_xt__,[-pi,+pi]); set(gca,'ydir','normal'); axisnotick; title('angle(M_x_p_r_xt__)','Interpreter','none');
colormap(subplot_{2},colormap('hsv'));
colormap(subplot_{1},colormap_81s());
%%%%;
disp('returning'); return;
end;%if flag_check;
%%%%%%%%;

%%%%%%%%;
lambda = 1.0; %<-- here measure z in terms of lambda. ;
dzl = 1.0/64.0;
zl_max = 2.0;
n_zl = zl_max/dzl;
tmp_M_k_p_wk_ = M_k_p_wk_;
%%%%;
figure(1+nf);nf=nf+1;clf;figsml;
fontsize_use = 18;
gamma_use = 2.0;
colormap(colormap_pm().^gamma_use);
clim_ = [];
%%%%;
for nzl=0:n_zl-1;
tmp_M_k_p_wk_ = tmp_M_k_p_wk_.*exp(-i*pi*lambda*dzl*k_p_r_wk_.^2); %<-- paraxial fourier propagator ;
tmp_N_x_c_xx__ = ...
interp_k_p_to_x_c_xxnufft( ...
 n_x_M_u ...
,diameter_x_c ...
,n_x_M_u ...
,diameter_x_c ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,tmp_M_k_p_wk_.*weight_2d_k_p_wk_*(4*pi^2) ...
)/(sqrt(n_x_M_u^2)*dx^2)*n_w_max^2*4 ;
if isempty(clim_); clim_ = [-1,+1]*prctile(abs(real(tmp_N_x_c_xx__)),100,'all'); end;
subplot(1,1,1); cla;
imagesc_c(n_x_M_u,x_u_0_,n_x_M_u,x_u_1_,real(tmp_N_x_c_xx__),clim_,colormap_pm.^gamma_use); axis image;
set(gca,'FontSize',fontsize_use);
set(gca,'XTick',[-half_diameter_x_c,0,+half_diameter_x_c]);
set(gca,'YTick',[-half_diameter_x_c,0,+half_diameter_x_c]);
set(gca,'TickLength',[0,0]);
xlabel('$x_{1}$','Interpreter','latex');ylabel('$x_{2}$','Interpreter','latex');title('$\Re\{A(\vec{x})\}$','Interpreter','latex');
drawnow();
end;%for nzl=0:n_zl-1;
%%%%;



