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
n_pixel = 7;
y_ = 0.25*randn(n_pixel,1) + i*0.25*randn(n_pixel,1);
y0 = y_(1+0);
y1 = y_(1+1);
y2 = y_(1+2);
y3 = y_(1+3);
y4 = y_(1+4);
y5 = y_(1+5);
y6 = y_(1+6);
x_ = 0.25*randn(n_pixel,1) + i*0.25*randn(n_pixel,1);
x0 = x_(1+0);
x1 = x_(1+1);
x2 = x_(1+2);
x3 = x_(1+3);
x4 = x_(1+4);
x5 = x_(1+5);
x6 = x_(1+6);
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
  + conj(y3 - N_*x3 - (u_ + i*v_)).*(y3 - N_*x3 - (u_ + i*v_)) ...
  + conj(y4 - N_*x4 - (u_ + i*v_)).*(y4 - N_*x4 - (u_ + i*v_)) ...
  + conj(y5 - N_*x5 - (u_ + i*v_)).*(y5 - N_*x5 - (u_ + i*v_)) ...
  + conj(y6 - N_*x6 - (u_ + i*v_)).*(y6 - N_*x6 - (u_ + i*v_)) ...
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
x_max = 32;
I_quad = integral(exp_2,1e-9,x_max,'AbsTol',1e-9,'RelTol',1e-9,'WayPoints',[1e-9:1:x_max]);
I_form = (2*pi)^(3/2-n_pixel/2) * (I_11*I_A).^(-1/2) * 0.5 * (-tmp_factor).^(2 - n_pixel/2) * gamma(n_pixel/2 - 2);
disp(sprintf(' %% I_form vs I_quad: %0.16f',fnorm(I_form-I_quad)/max(1e-12,fnorm(I_form))));
lI_form = (3/2-n_pixel/2)*log(2*pi) -0.5*log(I_11*I_A) - log(2) - (n_pixel/2-2)*log(-tmp_factor) + gammaln(n_pixel/2 - 2);
I_form = exp(lI_form);
disp(sprintf(' %% I_form vs I_quad: %0.16f',fnorm(I_form-I_quad)/max(1e-12,fnorm(I_form))));
%%%%%%%%;

%%%%%%%%;
n_pixel = 9; %<-- should be more than 3+1. ;
tmp_factor = -(I_B^2/I_A + I_C)/(2*I_11);
exp_2 = @(lambda_) (2*pi)^(3/2-n_pixel/2) * (I_11*I_A).^(-1/2) * (lambda_).^(3-n_pixel) .* exp( -tmp_factor ./ lambda_.^2 );
x_max = 32;
I_quad = integral(exp_2,1e-9,x_max,'AbsTol',1e-9,'RelTol',1e-9,'WayPoints',[1e-9:1:x_max]);
I_form = (2*pi)^(3/2-n_pixel/2) * (I_11*I_A).^(-1/2) * 0.5 * (tmp_factor).^(2 - n_pixel/2) * gamma(n_pixel/2 - 2);
disp(sprintf(' %% I_form vs I_quad: %0.16f',fnorm(I_form-I_quad)/max(1e-12,fnorm(I_form))));
lI_form = (3/2-n_pixel/2)*log(2*pi) -0.5*log(I_11*I_A) - log(2) - (n_pixel/2-2)*log(tmp_factor) + gammaln(n_pixel/2 - 2);
I_form = exp(lI_form);
disp(sprintf(' %% I_form vs I_quad: %0.16f',fnorm(I_form-I_quad)/max(1e-12,fnorm(I_form))));
%%%%%%%%;










