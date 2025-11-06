function ...
S_k_p_ = ...
interp_x_c_to_k_p_xxnufft( ...
 n_x1 ...
,diameter_x1_c ...
,n_x2 ...
,diameter_x2_c ...
,S_x_c_ ...
,n_r ...
,grid_k_p_ ...
,n_w_ ...
) ;
%%%%%%%%;
% Assuming that S_x_c_ is the original image ordered as: ;
% S_x_c_(nx1 + nx2*n_x1) = pixel nx1,nx2 ;
% we use xxnufft2d2 to calculate ;
% S_k_p_(na) = S_k_p_(nw + n_w_csum_(nr)) ;
% given by: ;
% S_k_p_(na) = 1/Z * \sum_{nx1=0,nx2=0}^{nx1=n_x1-1,nx2=n_x2-1}  ;
%              S_x_c_(nx1 + nx2*n_x1)  ;
%              * exp(-i * (x1_(nx1)*k1+x2_(nx2)*k2)) ;
% where:  ;
% Z = sqrt(n_x1*n_x2) ;
% and: ;
% k1 = 2*pi*grid_k_p_(nr)*cos(2*pi*nw/n_w_(nr)) ;
% k2 = 2*pi*grid_k_p_(nr)*sin(2*pi*nw/n_w_(nr)) ;
% and ;
% x1_(nx1) = half_diameter_x1_c*(-1 + 2*nx1/n_x1) ;
% x2_(nx2) = half_diameter_x2_c*(-1 + 2*nx2/n_x2) ;
% or equivalently: ;
% x1_(nx1) = diameter_x1_c/n_x1*(nx1-n_x1/2) ;
% x2_(nx2) = diameter_x2_c/n_x2*(nx2-n_x2/2) ;
% where: ;
% nx1-n_x1/2 runs from -n_x1/2 up to +n_x1/2-1 ;
% nx2-n_x2/2 runs from -n_x2/2 up to +n_x2/2-1 ;
% %%%%%%%% ;
% Note that grid_k_p_ is assumed to provide the k-magnitude ;
% in terms of ordinary frequency (not angular frequency). ;
% Hence the 2*pi in the definitions of k1 and k2. ; ;
% %%%%%%%% ;
% We use use the type-2 xxnufft: ;
% int xxnufft2d2(BIGINT nj,FLT* xj,FLT *yj,CPX* cj,int iflag,FLT eps, ;
%                BIGINT ms, BIGINT mt, CPX* fk) ;
% cj[j] =  SUM   fk[k1,k2] exp(+/-i (k1 xj[j] + k2 yj[j]))      for j = 0,...,nj-1 ;
%         k1,k2 ;
% where sum is over -ms/2 <= k1 <= (ms-1)/2, -mt/2 <= k2 <= (mt-1)/2, ;
% Inputs: ;
% nj     number of targets (int64, aka BIGINT) ;
% xj,yj     x,y locations of targets (each a size-nj FLT array) in [-3pi,3pi] ;
% fk     FLT complex array of Fourier transform values (size ms*mt, ;
%        increasing fast in ms then slow in mt, ie Fortran ordering). ;
% iflag  if >=0, uses + sign in exponential, otherwise - sign (int) ;
% eps    precision requested (>1e-16) ;
% ms,mt  numbers of Fourier modes given in x and y (int64) ;
%        each may be even or odd; ;
%        in either case the mode range is integers lying in [-m/2, (m-1)/2]. ;
% Outputs: ;
% cj     size-nj complex FLT array of target values ;
%        (ie, stored as 2*nj FLTs interleaving Re, Im). ;
% returned value - 0 if success, else see ../docs/usage.rst ;
% %%%%%%%% ;
% Because of the conventions we make regarding the image and frequencies, ;
% we call the type-2 xxnufft using: ;
% nj = n_w_sum = n_A ;
% xj = k1_*diameter_x1_c/n_x1 = 2*pi*grid_k_p_(nr)*cos(2*pi*nw/n_w_(nr))*diameter_x1_c/n_x1 ;
% yj = k2_*diameter_x2_c/n_x2 = 2*pi*grid_k_p_(nr)*sin(2*pi*nw/n_w_(nr))*diameter_x2_c/n_x2 ;
% cj = S_k_p_ ;
% iflag = -1 ;
% eps = 1e-6 ;
% ms = n_x1 ;
% mt = n_x2 ;
% fk = S_x_c_ ;
%%%%%%%%;

if (nargin<1);
verbose=1;
%%%%%%%%;
% First define integral of <f,f>. ;
%%%%%%%%;
h2d_ = @(kd) 4*pi^2*(besselj(0,kd) + besselj(2,kd)); % calculates <f_j,f_k>, normalized so that <f,f> = (4*pi^2);
dh2d_ = @(kd) 4*pi^3*(besselj(-1,kd) - besselj(+3,kd));
h3d_ = @(kd) 4*pi*( sin(kd) - (kd).*cos(kd) ) ./ kd.^3 ; % calculates <f_j,f_k>, normalized so that <f,f> = 4*pi/3;
dh3d_ = @(kd) 12*pi*( (kd.^2/3 - 1) .* sin(kd) + (kd).*cos(kd) ) ./ kd.^4 ;
%%%%%%%%;
% Define spatial grid. ;
%%%%%%%%;
half_diameter_x_c = 1.0d0;
x_p_r_max = half_diameter_x_c;
diameter_x_c = 2*half_diameter_x_c;
n_x_M_u = 128;
x_c_0_ = -half_diameter_x_c + transpose([0:n_x_M_u-1]/n_x_M_u)*diameter_x_c;
x_c_1_ = -half_diameter_x_c + transpose([0:n_x_M_u-1]/n_x_M_u)*diameter_x_c;
x_c_0_ = transpose(linspace(-x_p_r_max,+x_p_r_max,n_x_M_u+1)); x_c_0_ = x_c_0_(1:n_x_M_u);
x_c_1_ = transpose(linspace(-x_p_r_max,+x_p_r_max,n_x_M_u+1)); x_c_1_ = x_c_1_(1:n_x_M_u);
%x_c_0_ = transpose(linspace(-x_p_r_max,+x_p_r_max,n_x_M_u));
%x_c_1_ = transpose(linspace(-x_p_r_max,+x_p_r_max,n_x_M_u));
[x_c_0__,x_c_1__] = ndgrid(x_c_0_,x_c_1_);
dx = mean(diff(x_c_0_));
n_xx_M_u = n_x_M_u^2;
%%%%%%%%;
% Now set up and test k-quadrature on sphere. ;
%%%%%%%%;
k_int = 2*48; k_eq_d_double = 1.0; n_w_int = 1;
flag_verbose=0; k_p_r_max = k_int/(2*pi); k_eq_d = k_eq_d_double/(2*pi); str_L = 'L';
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
n_w_max = n_w_int*2*(k_int+1);
n_w_0in_ = n_w_max*ones(n_k_p_r,1);
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
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);
%%%%%%%%;
n_source = 4;
rng(0);
delta_a_c__ = zeros(2,n_source);
delta_b_c__ = zeros(2,n_source);
for nsource=0:n_source-1;
delta_a_c_ = 0.125*(2*rand(2,1)-1);
delta_a_c__(:,1+nsource) = delta_a_c_;
delta_b_c_ = 0.125*(2*rand(2,1)-1);
delta_b_c__(:,1+nsource) = delta_b_c_;
end;%for nsource=0:n_source-1;
a_k_p_form_ = zeros(n_w_sum,1);
b_k_p_form_ = zeros(n_w_sum,1);
for nsource=0:n_source-1;
delta_a_c_ = delta_a_c__(:,1+nsource);
a_k_p_form_ = a_k_p_form_ + exp(+i*2*pi*(k_c_0_wk_*delta_a_c_(1+0) + k_c_1_wk_*delta_a_c_(1+1)));
delta_b_c_ = delta_b_c__(:,1+nsource);
b_k_p_form_ = b_k_p_form_ + exp(+i*2*pi*(k_c_0_wk_*delta_b_c_(1+0) + k_c_1_wk_*delta_b_c_(1+1)));
end;%for nsource=0:n_source-1;
I_a_quad = sum(a_k_p_form_.*weight_2d_wk_)*(2*pi)^2;
I_b_quad = sum(b_k_p_form_.*weight_2d_wk_)*(2*pi)^2;
I_a_form = 0;
I_b_form = 0;
for nsource=0:n_source-1;
delta_a_c_ = delta_a_c__(:,1+nsource);
delta_b_c_ = delta_b_c__(:,1+nsource);
tmp_kd = 2*pi*k_p_r_max*fnorm(delta_a_c_);
I_a_form = I_a_form + h2d_(tmp_kd)*k_p_r_max^2/(4*pi);
tmp_kd = 2*pi*k_p_r_max*fnorm(delta_b_c_);
I_b_form = I_b_form + h2d_(tmp_kd)*k_p_r_max^2/(4*pi);
end;%for nsource=0:n_source-1;
disp(sprintf(' %% I_a_form vs I_a_quad %0.16f',fnorm(I_a_form-I_a_quad)/fnorm(I_a_form)));
disp(sprintf(' %% I_b_form vs I_b_quad %0.16f',fnorm(I_b_form-I_b_quad)/fnorm(I_b_form)));
%%%%%%%%;
a_k_p_form_l2 = sum(abs(a_k_p_form_).^2.*weight_2d_wk_)*(2*pi)^2;
a_x_c_quad_ = interp_k_p_to_x_c_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,a_k_p_form_.*weight_2d_wk_*(2*pi)^2)*sqrt(n_x_M_u^2) * n_w_sum;
a_x_c_quad_l2 = sum(abs(a_x_c_quad_).^2,'all')*dx^2;
b_k_p_form_l2 = sum(abs(b_k_p_form_).^2.*weight_2d_wk_)*(2*pi)^2;
b_x_c_quad_ = interp_k_p_to_x_c_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,b_k_p_form_.*weight_2d_wk_*(2*pi)^2)*sqrt(n_x_M_u^2) * n_w_sum;
b_x_c_quad_l2 = sum(abs(b_x_c_quad_).^2,'all')*dx^2;
disp(sprintf(' %% a_k_p_form_l2 vs a_x_c_quad_l2: %0.16f',fnorm(a_k_p_form_l2 - a_x_c_quad_l2)/fnorm(a_k_p_form_l2)));
disp(sprintf(' %% b_k_p_form_l2 vs b_x_c_quad_l2: %0.16f',fnorm(b_k_p_form_l2 - b_x_c_quad_l2)/fnorm(b_k_p_form_l2)));
%%%%%%%%;
a_x_c_form_ = zeros(n_x_M_u,n_x_M_u);
b_x_c_form_ = zeros(n_x_M_u,n_x_M_u);
for nsource=0:n_source-1;
delta_a_c_ = delta_a_c__(:,1+nsource);
delta_b_c_ = delta_b_c__(:,1+nsource);
tmp_kd__ = 2*pi * k_p_r_max * sqrt( (x_c_0__ + delta_a_c_(1+0)).^2 + (x_c_1__ + delta_a_c_(1+1)).^2 ) ;
a_x_c_form_ = a_x_c_form_ + h2d_(tmp_kd__)*k_p_r_max^2/(4*pi);
tmp_kd__ = 2*pi * k_p_r_max * sqrt( (x_c_0__ + delta_b_c_(1+0)).^2 + (x_c_1__ + delta_b_c_(1+1)).^2 ) ;
b_x_c_form_ = b_x_c_form_ + h2d_(tmp_kd__)*k_p_r_max^2/(4*pi);
end;%for nsource=0:n_source-1;
disp(sprintf(' %% a_x_c_form_(:) vs a_x_c_quad_(:): %0.16f',fnorm(a_x_c_form_(:) - a_x_c_quad_(:))/fnorm(a_x_c_form_(:))));
disp(sprintf(' %% b_x_c_form_(:) vs b_x_c_quad_(:): %0.16f',fnorm(b_x_c_form_(:) - b_x_c_quad_(:))/fnorm(b_x_c_form_(:))));
a_x_c_form_l2 = sum(abs(a_x_c_form_(:)).^2)*dx^2;
b_x_c_form_l2 = sum(abs(b_x_c_form_(:)).^2)*dx^2;
a_x_c_quad_l2 = sum(abs(a_x_c_quad_(:)).^2)*dx^2;
b_x_c_quad_l2 = sum(abs(b_x_c_quad_(:)).^2)*dx^2;
disp(sprintf(' %% a_x_c_form_l2 vs a_x_c_quad_l2: %0.16f',fnorm(a_x_c_form_l2 - a_x_c_quad_l2)/fnorm(a_x_c_form_l2)));
disp(sprintf(' %% b_x_c_form_l2 vs b_x_c_quad_l2: %0.16f',fnorm(b_x_c_form_l2 - b_x_c_quad_l2)/fnorm(b_x_c_form_l2)));
%%%%%%%%;
flag_plot=0;
if flag_plot;
subplot(1,3,1);imagesc(reshape(real(a_x_c_form_),[n_x_M_u,n_x_M_u])); axis image; axisnotick; colorbar;
subplot(1,3,2);imagesc(reshape(real(a_x_c_quad_),[n_x_M_u,n_x_M_u])); axis image; axisnotick; colorbar;
subplot(1,3,3);imagesc(reshape(real(a_x_c_form_ - a_x_c_quad_),[n_x_M_u,n_x_M_u])); axis image; axisnotick; colorbar;
end;%if flag_plot;
%%%%%%%%;
a_k_p_quad_ = interp_x_c_to_k_p_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,a_x_c_form_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_M_u^2)*dx^2;
a_k_p_quad_l2 = sum(abs(a_k_p_quad_).^2.*weight_2d_wk_)*(2*pi)^2;
b_k_p_quad_ = interp_x_c_to_k_p_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,b_x_c_form_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_M_u^2)*dx^2;
b_k_p_quad_l2 = sum(abs(b_k_p_quad_).^2.*weight_2d_wk_)*(2*pi)^2;
disp(sprintf(' %% a_x_c_form_l2 vs a_k_p_quad_l2: %0.16f',fnorm(a_x_c_form_l2 - a_k_p_quad_l2)/fnorm(a_x_c_form_l2)));
disp(sprintf(' %% b_x_c_form_l2 vs b_k_p_quad_l2: %0.16f',fnorm(b_x_c_form_l2 - b_k_p_quad_l2)/fnorm(b_x_c_form_l2)));
disp(sprintf(' %% a_k_p_form_l2 vs a_k_p_quad_l2: %0.16f',fnorm(a_k_p_form_l2 - a_k_p_quad_l2)/fnorm(a_k_p_form_l2)));
disp(sprintf(' %% b_k_p_form_l2 vs b_k_p_quad_l2: %0.16f',fnorm(b_k_p_form_l2 - b_k_p_quad_l2)/fnorm(b_k_p_form_l2)));
%%%%%%%%;
% Now test gaussian. ;
%%%%%%%%;
tmp_sigma_x_c = 0.0625;
tmp_sigma_k_p = 1/tmp_sigma_x_c;
tmp_delta_ = 0.75*[+0.1,-0.2];
tmp_M_x_c_ = 1/(sqrt(2*pi)*tmp_sigma_x_c)^2 * exp( -( (x_c_0__-tmp_delta_(1+0)).^2 + (x_c_1__-tmp_delta_(1+1)).^2 ) / (2*tmp_sigma_x_c^2) );
tmp_M_x_c_l2 = sum(tmp_M_x_c_.^2,'all')*dx^2;
disp(sprintf(' %% sum(tmp_M_x_c_*dx^2,''all'') = %0.16f',sum(tmp_M_x_c_*dx^2,'all')));
disp(sprintf(' %% tmp_M_x_c_l2 = %0.16f',tmp_M_x_c_l2));
tmp_M_k_p_ = interp_x_c_to_k_p_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,tmp_M_x_c_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_M_u^2)*dx^2;
tmp_M_k_p_l2 = sum(abs(tmp_M_k_p_).^2 .* weight_2d_wk_) * (2*pi)^2;
disp(sprintf(' %% tmp_M_k_p_l2 = %0.16f',tmp_M_k_p_l2));
tmp_M_k_p_form_ = zeros(n_w_sum,1);
na=0;
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
n_w = n_w_(1+nk_p_r);
for nw=0:n_w-1;
k_x_c_0 = k_p_r*cos(2*pi*nw/n_w);
k_x_c_1 = k_p_r*sin(2*pi*nw/n_w);
tmp_M_k_p_form_(1+na) = exp( -( (2*pi*k_x_c_0).^2 + (2*pi*k_x_c_1).^2 ) / (2/tmp_sigma_x_c^2) ) .* exp( - 2*pi*i*( k_x_c_0*tmp_delta_(1+0) + k_x_c_1*tmp_delta_(1+1) ) );
na=na+1;
end;%for nw=0:n_w-1;
end;%for nk_p_r=0:n_k_p_r-1;
tmp_M_k_p_form_l2 = sum(abs(tmp_M_k_p_form_).^2 .* weight_2d_wk_) * (2*pi)^2;
disp(sprintf(' %% tmp_M_k_p_form_l2 = %0.16f',tmp_M_k_p_form_l2));
disp(sprintf(' %% tmp_M_k_p_ vs tmp_M_k_p_form: %0.16f',fnorm(tmp_M_k_p_ - tmp_M_k_p_form_)/fnorm(tmp_M_k_p_)));
tmp_M_x_c_reco_ = interp_k_p_to_x_c_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_M_k_p_.*weight_2d_wk_*(2*pi)^2)*sqrt(n_x_M_u^2) * n_w_sum;
tmp_M_x_c_reco_l2 = sum(abs(tmp_M_x_c_reco_).^2,'all')*dx^2;
disp(sprintf(' %% tmp_M_x_c_reco_l2 = %0.16f',tmp_M_x_c_reco_l2));
disp(sprintf(' %% tmp_M_x_c_ vs tmp_M_x_c_reco: %0.16f',fnorm(tmp_M_x_c_ - tmp_M_x_c_reco_)/fnorm(tmp_M_x_c_)));
figure(1);figbig;fig80s;
p_row = 2; p_col = 3; np=0;
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,tmp_M_x_c_);axis image;axisnotick; title('tmp_M_x_c_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(tmp_M_k_p_));axis image;axisnotick; title('tmp_M_k_p_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(tmp_M_k_p_form_));axis image;axisnotick; title('tmp_M_k_p_form_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,real(tmp_M_x_c_reco_));axis image;axisnotick; title('tmp_M_x_c_reco_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,real(tmp_M_x_c_reco_-tmp_M_x_c_));axis image;axisnotick; title('difference','Interpreter','none');
disp('returning'); return;
end;%if (nargin<1);

flag_verbose=0;
iflag = -1;
eps = 1.0d-6;
flag_2_versus_3 = 0; %<-- logical: compare type-2 and type-3 results. ;
dx1 = diameter_x1_c/max(1,n_x1) ;
dx2 = diameter_x2_c/max(1,n_x2) ;
n_A = sum(n_w_) ;
k1_ = zeros(n_A,1);
k2_ = zeros(n_A,1);

if (flag_2_versus_3);
tmp_S_k_p_ = zeros(n_A,1);
x1_ = zeros(n_x1*n_x2,1);
x2_ = zeros(n_x1*n_x2,1);
end; %if (flag_2_versus_3);

S_k_p_ = zeros(n_A,1);
na=0;
for nr=0:n_r-1;
r = 2*pi*grid_k_p_(1+nr);
for nw=0:n_w_(1+nr)-1;
omega = (2*pi*nw)/max(1,n_w_(1+nr));
k1_(1+na) = r*cos(omega)*dx1;
k2_(1+na) = r*sin(omega)*dx2;
na = na+1;
end;%for nw=0:n_w_(1+nr)-1;
end;%for nr=0:n_r-1;
if (na~=n_A);
disp(sprintf(' %% Warning, na~=n_A in interp_x_c_to_k_p_xxnufft.m'));
end;%if (na~=n_A);
[S_k_p_,ier] = xxnufft2d2(n_A,k1_,k2_,iflag,eps,n_x1,n_x2,S_x_c_);
if (ier~=0); disp(sprintf(' %% Warning! precision out of range in xxnufft2d2 in interp_x_c_to_k_p_xxnufft.m')); end;
if (flag_2_versus_3);
for nx1=0:n_x1-1;
for nx2=0:n_x2-1;
x1_(1+nx1+nx2*n_x1) = nx1-(1.0d0*n_x1)/2.0d0;
x2_(1+nx1+nx2*n_x1) = nx2-(1.0d0*n_x2)/2.0d0;
end%for nx2=0:n_x2-1;
end%for nx1=0:n_x1-1;
[tmp_S_k_p_,ier] = xxnufft2d3(n_x1*n_x2,x1_,x2_,S_x_c_,iflag,eps,n_A,k1_,k2_);
if (ier~=0); disp(sprintf(' %% Warning! precision out of range in xxnufft2d3 in interp_x_c_to_k_p_xxnufft.m')); end;
print_sub(n_x1*n_x2,x1_,' x1_: ');
print_sub(n_x1*n_x2,x2_,' x2_: ');
print_sub(n_A,    S_k_p_,' 2: S_k_p_: ');
print_sub(n_A,tmp_S_k_p_,' 3: S_k_p_: ');
disp(sprintf(' %% type_2_versus_3: %0.16f',norm(S_k_p_-tmp_S_k_p_,'fro')/max(norm(S_k_p_,'fro'),norm(tmp_S_k_p_,'fro'))));
end;%if (flag_2_versus_3);
S_k_p_ = S_k_p_/sqrt(n_x1*n_x2);

