nf=0;

%%%%%%%%;
n_pixel = 7;
y_ = randn(n_pixel,1) + i*randn(n_pixel,1);
x_ = randn(n_pixel,1) + i*randn(n_pixel,1);
lambda = 0.3;
u = -1.2; v = 3.4;
mu = u + i*v;
N = 0.7;
z_ = N*x_ + mu;
e_ = ones(n_pixel,1);
I_all_0 = sum(abs(y_ - z_).^2,'all');
I_yy = sum(conj(y_).*y_,'all');
I_xy = sum(conj(x_).*y_,'all');
I_1y = sum(conj(e_).*y_,'all');
I_yx = sum(conj(y_).*x_,'all');
I_xx = sum(conj(x_).*x_,'all');
I_1x = sum(conj(e_).*x_,'all');
I_y1 = sum(conj(y_).*e_,'all');
I_x1 = sum(conj(x_).*e_,'all');
I_11 = sum(conj(e_).*e_,'all');
I_all_1 = I_yy + N^2*I_xx - 2*N*real(I_xy) + u^2*I_11 + v^2*I_11 - 2*u*(real(I_1y) - N*real(I_1x)) - 2*v*(imag(I_1y) - N*imag(I_1x));
disp(sprintf(' %% I_all_0 vs I_all_1: %0.16f',fnorm(I_all_0-I_all_1)/max(1e-12,fnorm(I_all_0))));
%%%%%%%%;

%%%%%%%%;
a = 3.0;
b = -0.5;
g = @(x_) exp(-a*x_.^2 + b*x_);
I_quad = integral(g,-10,+10);
I_form = sqrt(pi/a) * exp(b^2/(4*a));
disp(sprintf(' %% I_form vs I_quad: %0.16f',fnorm(I_form-I_quad)/max(1e-12,fnorm(I_form))));
%%%%%%%%;

%%%%%%%%;
a = 3.0;
b = -0.5;
g = @(x_) exp(-a*x_.^2 + 2*b*x_);
I_quad = integral(g,-10,+10);
I_form = sqrt(pi/a) * exp(b^2/a);
disp(sprintf(' %% I_form vs I_quad: %0.16f',fnorm(I_form-I_quad)/max(1e-12,fnorm(I_form))));
%%%%%%%%;
 
%%%%%%%%;
a = 3.0;
b = -0.5;
lambda = 1.75;
g = @(x_) exp(-(a*x_.^2 + 2*b*x_)/(2*lambda^2));
I_quad = integral(g,-10,+10);
I_form = sqrt(2*pi*lambda^2/a) * exp(b^2/(a*2*lambda^2));
disp(sprintf(' %% I_form vs I_quad: %0.16f',fnorm(I_form-I_quad)/max(1e-12,fnorm(I_form))));
%%%%%%%%;

%%%%%%%%;
a = 3.0;
b = 2.25*i;
h = @(x_) exp(-a*x_.^2 + i*b*x_);
I_quad = integral(h,-10,+10);
I_form = sqrt(pi/a) * exp(-b^2/(4*a));
disp(sprintf(' %% I_form vs I_quad: %0.16f',fnorm(I_form-I_quad)/max(1e-12,fnorm(I_form))));
%%%%%%%%;

%%%%%%%%;
% perform integral calculation for complex constant offset. ;
%%%%%%%%;
n_pixel = 3;
y_ = 0.25*randn(n_pixel,1) + i*0.25*randn(n_pixel,1);
y0 = y_(1+0);
y1 = y_(1+1);
y2 = y_(1+2);
x_ = 0.25*randn(n_pixel,1) + i*0.25*randn(n_pixel,1);
x0 = x_(1+0);
x1 = x_(1+1);
x2 = x_(1+2);
lambda = 0.3;
u = -0.2; v = 0.45;
mu = u + i*v;
N = 0.37;
z_ = N*x_ + mu;
e_ = ones(n_pixel,1);
I_all_0 = @(u_,v_,N_) ...
  + conj(y0 - N_*x0 - (u_ + i*v_)).*(y0 - N_*x0 - (u_ + i*v_)) ...
  + conj(y1 - N_*x1 - (u_ + i*v_)).*(y1 - N_*x1 - (u_ + i*v_)) ...
  + conj(y2 - N_*x2 - (u_ + i*v_)).*(y2 - N_*x2 - (u_ + i*v_)) ...
  ;
I_yy = sum(conj(y_).*y_,'all');
I_xy = sum(conj(x_).*y_,'all');
I_1y = sum(conj(e_).*y_,'all');
I_yx = sum(conj(y_).*x_,'all');
I_xx = sum(conj(x_).*x_,'all');
I_1x = sum(conj(e_).*x_,'all');
I_y1 = sum(conj(y_).*e_,'all');
I_x1 = sum(conj(x_).*e_,'all');
I_11 = sum(conj(e_).*e_,'all');
I_A = -real(I_1x)^2 - imag(I_1x)^2 + I_xx*I_11;
I_B = -real(I_1x)*real(I_1y) - imag(I_1x)*imag(I_1y) + real(I_xy)*I_11;
I_C = +real(I_1y)^2 + imag(I_1y)^2 - I_yy*I_11;
I_tmp_0 = @(N) ...
  (I_A*N.^2 -2*I_B*N - I_C)./I_11 ...
  ;
exp_all_0 = @(u_,v_,N_) exp(-I_all_0(u_,v_,N_)./(2*lambda.^2)) .* (2*pi*lambda.^2).^(-n_pixel/2) ;
exp_all_1 = @(u_,v_) exp_all_0(u_,v_,N);
exp_tmp_0 = @(N_) exp(-I_tmp_0(N_)./(2*lambda.^2)) .* (2*pi*lambda.^2).^(-n_pixel/2) .* (2*pi*lambda.^2./I_11) ;
exp_tmp_1 = exp_tmp_0(N);
I_quad_1 = integral2(exp_all_1,-10,+10,-10,+10);
I_form_1 = exp_tmp_1;
disp(sprintf(' %% I_form_1 vs I_quad_1: %0.16f',fnorm(I_form_1-I_quad_1)/max(1e-12,fnorm(I_form_1))));
%%%%%%%%;
I_quad_0 = integral(exp_tmp_0,-10,+10);
I_form_0 = sqrt(2*pi*lambda^2*I_11/I_A)*exp((I_B^2/(2*lambda^2*I_11*I_A)))*exp(I_C/(2*lambda^2*I_11)) .* (2*pi*lambda.^2).^(-n_pixel/2) .* (2*pi*lambda.^2./I_11);
disp(sprintf(' %% I_form_0 vs I_quad_0: %0.16f',fnorm(I_form_0-I_quad_0)/max(1e-12,fnorm(I_form_0))));
%%%%%%%%;
exp_tmp_0 = @(N_) exp(-(I_A*N_.^2 - 2*I_B*N_ - I_C)/(2*lambda^2*I_11)) .* (2*pi*lambda.^2).^(-n_pixel/2) .* (2*pi*lambda.^2./I_11) ;
I_quad_0 = integral(exp_tmp_0,-10,+10);
I_form_0 = sqrt(2*pi*lambda^2*I_11/I_A)*exp(I_B^2/I_A/(2*lambda^2*I_11))*exp(I_C/(2*lambda^2*I_11)) .* (2*pi*lambda.^2).^(-n_pixel/2) .* (2*pi*lambda.^2./I_11) ;
disp(sprintf(' %% I_form_0 vs I_quad_0: %0.16f',fnorm(I_form_0-I_quad_0)/max(1e-12,fnorm(I_form_0))));
I_form_0 = (2*pi).^(3/2-n_pixel/2)*(I_11*I_A).^(-1/2)*lambda.^(3-n_pixel) * exp(I_B^2/I_A/(2*lambda^2*I_11))*exp(I_C/(2*lambda^2*I_11)) ;
disp(sprintf(' %% I_form_0 vs I_quad_0: %0.16f',fnorm(I_form_0-I_quad_0)/max(1e-12,fnorm(I_form_0))));
%%%%%%%%;

%%%%%%%%;
n = -5.5; a = 3.25;
egamma = @(lambda) lambda.^n .* exp(-a./lambda.^2);
x_max = 64;
I_quad = integral(egamma,1e-9,x_max,'AbsTol',1e-9,'RelTol',1e-9,'WayPoints',[1e-9:1:x_max]);
I_form = 0.5 * a.^(+n/2+1/2) * gamma(-n/2 - 1/2);
%x_ = linspace(0,x_max,1024); plot(x_,egamma(x_),'r-','LineWidth',3);
disp(sprintf(' %% I_form vs I_quad: %0.16f',fnorm(I_form-I_quad)/max(1e-12,fnorm(I_form))));
%%%%%%%%;

%%%%%%%%;
n_pixel = 7; %<-- should be more than 3+1. ;
tmp_factor = (I_B^2/I_A + I_C)/(2*I_11);
exp_2 = @(lambda_) (2*pi)^(3/2-n_pixel/2) * (I_11*I_A).^(-1/2) * (lambda_).^(3-n_pixel) .* exp( tmp_factor ./ lambda_.^2 );
if tmp_factor> 0;
disp(sprintf(' %% Warning, tmp_factor> 0 (integral not defined);',tmp_factor));
end;%if tmp_factor> 0;
x_max = -tmp_factor*128;
I_quad = integral(exp_2,1e-9,x_max,'AbsTol',1e-9,'RelTol',1e-9,'WayPoints',[1e-9:x_max/1024:x_max]);
I_form = (2*pi)^(3/2-n_pixel/2) * (I_11*I_A).^(-1/2) * 0.5 * (-tmp_factor).^(2 - n_pixel/2) * gamma(n_pixel/2 - 2);
disp(sprintf(' %% I_form vs I_quad: %0.16f',fnorm(I_form-I_quad)/max(1e-12,fnorm(I_form))));
lI_form = (3/2-n_pixel/2)*log(2*pi) -0.5*log(I_11*I_A) - log(2) - (n_pixel/2-2)*log(-tmp_factor) + gammaln(n_pixel/2 - 2);
I_form = exp(lI_form);
disp(sprintf(' %% I_form vs I_quad: %0.16f',fnorm(I_form-I_quad)/max(1e-12,fnorm(I_form))));
%%%%%%%%;

%%%%%%%%;
% Now look at the sinc. ;
%%%%%%%%;
half_diameter_x = 1.0;
full_diameter_x = 2*half_diameter_x;
n_x_0 = 128+0; %<-- just to check dimensions. ;
n_x_1 = 128+2; %<-- just to check dimensions. ;
x_0_ = linspace(-half_diameter_x,+half_diameter_x,n_x_0+1); x_0_ = transpose(x_0_(1:n_x_0));
dx_0 = mean(diff(x_0_));
x_1_ = linspace(-half_diameter_x,+half_diameter_x,n_x_1+1); x_1_ = transpose(x_1_(1:n_x_1));
dx_1 = mean(diff(x_1_));
[x_0__,x_1__] = ndgrid(x_0_,x_1_);
n_xx = n_x_0*n_x_1; dxx = dx_0*dx_1;
%%%%%%%%;
sinc_0 = max(1e-12,dx_0*(32+1.0));
sinc_1 = max(1e-12,dx_1*(32+1.0));
%%%%%%%%;
half_diameter_k = 4*2*pi/sinc_0;
full_diameter_k = 2*half_diameter_k;
n_k_0 = 128+4; %<-- just to check dimensions. ;
n_k_1 = 128+6; %<-- just to check dimensions. ;
k_0_ = linspace(-half_diameter_k,+half_diameter_k,n_k_0+1); k_0_ = transpose(k_0_(1:n_k_0));
dk_0 = mean(diff(k_0_));
k_1_ = linspace(-half_diameter_k,+half_diameter_k,n_k_1+1); k_1_ = transpose(k_1_(1:n_k_1));
dk_1 = mean(diff(k_1_));
[k_0__,k_1__] = ndgrid(k_0_,k_1_);
n_kk = n_k_0*n_k_1; dkk = dk_0*dk_1;
%%%%%%%%;
sigma_x = sinc_0/3;
f_x = @(x) 1/sqrt(2*pi)/sigma_x*exp(-x_0_.^2/(2*sigma_x^2));
f_x_ = f_x(x_0_);
f_k_ = xxnufft1d3(n_x_0,x_0_,f_x_,-1,1e-6,n_k_0,k_0_)*dx_0;
disp(sprintf(' %% fnorm(imag(f_k_)) %0.16f',fnorm(imag(f_k_))));
f_k_ = real(f_k_);
sigma_k = 1./sigma_x;
g_k_ = exp(-k_0_.^2/(2*sigma_k^2));
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,2,1);
plot(x_0_,f_x_,'kx-'); xlabel('x'); ylabel('f(x)');
subplot(1,2,2);
plot(k_0_,f_k_,'kx-',k_0_,g_k_,'ro-'); xlabel('k'); ylabel('f(k)');
sgtitle('gaussian');
disp(sprintf(' %% sum(f_x_.^2)*dx_0: %0.16f %%<-- energy in real-space',sum(f_x_.^2,'all')*dx_0));
disp(sprintf(' %% sum(f_k_.^2)*dk_0/(2*pi): %0.16f %%<-- energy in fourier-space',sum(f_k_.^2,'all')*dk_0/(2*pi)));
%%%%%%%%;
f_x = @(x) cast(abs(x-0)<sinc_0/2,'double');
f_x_ = f_x(x_0_);
f_k_ = xxnufft1d3(n_x_0,x_0_,f_x_,-1,1e-6,n_k_0,k_0_)*dx_0;
disp(sprintf(' %% fnorm(imag(f_k_)) %0.16f',fnorm(imag(f_k_))));
f_k_ = real(f_k_);
s_k_ = (sinc_0).*sinc(sinc_0/(2*pi)*k_0_);
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,2,1);
plot(x_0_,f_x_,'kx-'); xlabel('x'); ylabel('f(x)');
subplot(1,2,2);
plot(k_0_,f_k_,'kx-',k_0_,s_k_,'ro-'); xlabel('k'); ylabel('f(k)');
sgtitle('square-wave');
disp(sprintf(' %% sum(f_x_.^2)*dx_0: %0.16f %%<-- energy in real-space',sum(f_x_.^2,'all')*dx_0));
disp(sprintf(' %% sum(f_k_.^2)*dk_0/(2*pi): %0.16f %%<-- not all energy captured in fourier-space',sum(f_k_.^2,'all')*dk_0/(2*pi)));
%%%%%%%%;
f_xx = @(x_0,x_1) cast((abs(x_0-0)<sinc_0/2) .* (abs(x_1-0)<sinc_1/2),'double') ;
f_xx__ = f_xx(x_0__,x_1__);
%%%%%%%%;
f_kk__ = reshape(xxnufft2d3(n_xx,x_0__(:),x_1__(:),f_xx__(:),-1,1e-6,n_kk,k_0__(:),k_1__(:))*dxx,[n_k_0,n_k_1]);
disp(sprintf(' %% fnorm(imag(f_kk__)) %0.16f',fnorm(imag(f_kk__))));
f_kk__ = real(f_kk__);
g_kk__ = (sinc_0).*sinc(sinc_0/(2*pi)*k_0__).*(sinc_1).*sinc(sinc_1/(2*pi)*k_1__);
disp(sprintf(' %% sum(f_xx__.^2)*dxx: %0.16f %%<-- energy in real-space',sum(f_xx__.^2,'all')*dxx));
disp(sprintf(' %% sum(f_kk__.^2)*dkk/(2*pi)^2: %0.16f %%<-- not all energy captured in fourier-space',sum(f_kk__.^2,'all')*dkk/(2*pi)^2));
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,3,1);
imagesc_c(n_x_0,x_0_,n_x_1,x_1_,f_xx__,[0,1],colormap_beach());
xlim(half_diameter_x*[-1,+1]); xlabel('x_0','Interpreter','none');
ylim(half_diameter_x*[-1,+1]); ylabel('x_1','Interpreter','none');
axisnotick; title('x-space');
subplot(1,3,2);
imagesc_c(n_k_0,k_0_,n_k_1,k_1_,f_kk__,sinc_0*sinc_1*[-1,+1],colormap_beach());
xlim(half_diameter_k*[-1,+1]); xlabel('k_0','Interpreter','none');
ylim(half_diameter_k*[-1,+1]); ylabel('k_1','Interpreter','none');
axisnotick; title('k-space nufft');
subplot(1,3,3);
imagesc_c(n_k_0,k_0_,n_k_1,k_1_,g_kk__,sinc_0*sinc_1*[-1,+1],colormap_beach());
xlim(half_diameter_k*[-1,+1]); xlabel('k_0','Interpreter','none');
ylim(half_diameter_k*[-1,+1]); ylabel('k_1','Interpreter','none');
axisnotick; title('k-space formula');
%%%%%%%%;

%%%%%%%%;
% Now repeat integral calculation for generic real-valued offset w. ;
%%%%%%%%;
n_pixel = 3;
y_ = 0.25*randn(n_pixel,1) + i*0.25*randn(n_pixel,1);
y0 = y_(1+0);
y1 = y_(1+1);
y2 = y_(1+2);
x_ = 0.25*randn(n_pixel,1) + i*0.25*randn(n_pixel,1);
x0 = x_(1+0);
x1 = x_(1+1);
x2 = x_(1+2);
w_ = 0.25*randn(n_pixel,1);
w0 = w_(1+0);
w1 = w_(1+1);
w2 = w_(1+2);
lambda = 0.3;
u = -0.2;
N = 0.37;
z_ = N*x_ + u*w_;
e_ = ones(n_pixel,1);
I_all_0 = @(u_,N_) ...
  + conj(y0 - N_*x0 - u_*w0).*(y0 - N_*x0 - u_*w0) ...
  + conj(y1 - N_*x1 - u_*w1).*(y1 - N_*x1 - u_*w1) ...
  + conj(y2 - N_*x2 - u_*w2).*(y2 - N_*x2 - u_*w2) ...
  ;
I_yy = sum(conj(y_).*y_,'all');
I_xy = sum(conj(x_).*y_,'all');
I_wy = sum(conj(w_).*y_,'all');
I_yx = sum(conj(y_).*x_,'all');
I_xx = sum(conj(x_).*x_,'all');
I_wx = sum(conj(w_).*x_,'all');
I_yw = sum(conj(y_).*w_,'all');
I_xw = sum(conj(x_).*w_,'all');
I_ww = sum(conj(w_).*w_,'all');
I_A = -real(I_wx)^2 + I_xx*I_ww;
I_B = -real(I_wx)*real(I_wy) + real(I_xy)*I_ww;
I_C = +real(I_wy)^2 - I_yy*I_ww;
I_tmp_0 = @(N) ...
  (I_A*N.^2 -2*I_B*N - I_C)./I_ww ...
  ;
exp_all_0 = @(u_,N_) exp(-I_all_0(u_,N_)./(2*lambda.^2)) .* (2*pi*lambda.^2).^(-n_pixel/2) ;
exp_all_1 = @(u_) exp_all_0(u_,N);
exp_tmp_0 = @(N_) exp(-I_tmp_0(N_)./(2*lambda.^2)) .* (2*pi*lambda.^2).^(-n_pixel/2) .* sqrt(2*pi*lambda.^2./I_ww) ;
exp_tmp_1 = exp_tmp_0(N);
I_quad_1 = integral(exp_all_1,-10,+10);
I_form_1 = exp_tmp_1;
disp(sprintf(' %% I_form_1 vs I_quad_1: %0.16f',fnorm(I_form_1-I_quad_1)/max(1e-12,fnorm(I_form_1))));
%%%%%%%%;
I_quad_0 = integral(exp_tmp_0,-10,+10);
I_form_0 = sqrt(2*pi*lambda^2*I_ww/I_A)*exp((I_B^2/(2*lambda^2*I_ww*I_A)))*exp(I_C/(2*lambda^2*I_ww)) .* (2*pi*lambda.^2).^(-n_pixel/2) .* sqrt(2*pi*lambda.^2./I_ww);
disp(sprintf(' %% I_form_0 vs I_quad_0: %0.16f',fnorm(I_form_0-I_quad_0)/max(1e-12,fnorm(I_form_0))));
%%%%%%%%;
exp_tmp_0 = @(N_) exp(-(I_A*N_.^2 - 2*I_B*N_ - I_C)/(2*lambda^2*I_ww)) .* (2*pi*lambda.^2).^(-n_pixel/2) .* sqrt(2*pi*lambda.^2./I_ww) ;
I_quad_0 = integral(exp_tmp_0,-10,+10);
I_form_0 = sqrt(2*pi*lambda^2*I_ww/I_A)*exp(I_B^2/I_A/(2*lambda^2*I_ww))*exp(I_C/(2*lambda^2*I_ww)) .* (2*pi*lambda.^2).^(-n_pixel/2) .* sqrt(2*pi*lambda.^2./I_ww) ;
disp(sprintf(' %% I_form_0 vs I_quad_0: %0.16f',fnorm(I_form_0-I_quad_0)/max(1e-12,fnorm(I_form_0))));
I_form_0 = (2*pi).^(1-n_pixel/2)*(I_A).^(-1/2)*lambda.^(2-n_pixel) * exp(I_B^2/I_A/(2*lambda^2*I_ww))*exp(I_C/(2*lambda^2*I_ww)) ;
disp(sprintf(' %% I_form_0 vs I_quad_0: %0.16f',fnorm(I_form_0-I_quad_0)/max(1e-12,fnorm(I_form_0))));
%%%%%%%%;

%%%%%%%%;
n_pixel = 7; %<-- should be more than 3+1. ;
tmp_factor = (I_B^2/I_A + I_C)/(2*I_ww);
if tmp_factor> 0;
disp(sprintf(' %% Warning, tmp_factor> 0 (integral not defined);',tmp_factor));
end;%if tmp_factor> 0;
exp_2 = @(lambda_) (2*pi)^(1-n_pixel/2) * (I_A).^(-1/2) * (lambda_).^(2-n_pixel) .* exp( tmp_factor ./ lambda_.^2 );
x_max = -tmp_factor*64;
I_quad = integral(exp_2,1e-9,x_max,'AbsTol',1e-9,'RelTol',1e-9,'WayPoints',[1e-9:x_max/1024:x_max]);
I_form = (2*pi)^(1-n_pixel/2) * (I_A).^(-1/2) * 0.5 * (-tmp_factor).^(1.5 - n_pixel/2) * gamma(n_pixel/2 - 1.5);
disp(sprintf(' %% I_form vs I_quad: %0.16f',fnorm(I_form-I_quad)/max(1e-12,fnorm(I_form))));
lI_form = (1-n_pixel/2)*log(2*pi) -0.5*log(I_A) - log(2) - (n_pixel/2-1.5)*log(-tmp_factor) + gammaln(n_pixel/2 - 1.5);
I_form = exp(lI_form);
disp(sprintf(' %% I_form vs I_quad: %0.16f',fnorm(I_form-I_quad)/max(1e-12,fnorm(I_form))));
%%%%%%%%;


 








