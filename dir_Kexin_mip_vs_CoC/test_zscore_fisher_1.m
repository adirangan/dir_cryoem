%%%%%%%%;
% treating l as the values to be averaged, rather than a distribution. ;
%%%%%%%%;
clear;
nf=0;
n_x = 256+1;
n_y = 256+1;
x_ = transpose(linspace(-2,+2,n_x));
dx = mean(diff(x_));
y_ = transpose(linspace(-2,+2,n_y));
dy = mean(diff(y_));
[x__,y__] = ndgrid(x_,y_);
sigma_x = 0.0850;
sigma_y = 0.1025;
rbar = 7.5;
nll = @(x,y) 0.5 * (x.*x/sigma_x^2 + y.*y/sigma_y^2) ;
l = @(x,y) rbar * 1/(2*pi)/sigma_x/sigma_y .* exp(-nll(x,y)) ;
l__ = l(x__,y__);
figure(1+nf);nf=nf+1;clf;fig80s;
imagesc(l__);
axis image; axisnotick; xlabel('x'); ylabel('y');
dom0 = sum(l__.^0,'all')*dx*dy;
disp(sprintf(' %% dom0 %0.16f',dom0));
avg0 = sum(l__.^1,'all')*dx*dy / dom0;
disp(sprintf(' %% avg0 %0.16f',avg0));
avg1 = rbar / dom0;
disp(sprintf(' %% avg1 %0.16f',avg1));
var0 = sum((l__ - avg0).^2,'all')*dx*dy / dom0;
disp(sprintf(' %% var0 %0.16f',var0));
var1 = sum(l__.^2,'all')*dx*dy / dom0 - avg1^2;
disp(sprintf(' %% var1 %0.16f',var1));
var2 = (rbar^2/dom0) * (1/pi/(2*sigma_x)/(2*sigma_y) - 1/dom0);
disp(sprintf(' %% var2 %0.16f',var2));
max0 = max(l__,[],'all');
disp(sprintf(' %% max0 %0.16f',max0));
max1 = rbar/(2*pi)/sigma_x/sigma_y;
disp(sprintf(' %% max1 %0.16f',max1));

nx_mid = ceil(n_x/2);
ny_mid = ceil(n_y/2);
dxxdl = (l__(nx_mid+1,ny_mid) - 2*l__(nx_mid+0,ny_mid) + l__(nx_mid-1,ny_mid))/max(1e-12,dx*dx);
dyydl = (l__(nx_mid,ny_mid+1) - 2*l__(nx_mid,ny_mid+0) + l__(nx_mid,ny_mid-1))/max(1e-12,dy*dy);
disp(sprintf(' %% dxxdl %0.16f',dxxdl));
disp(sprintf(' %% -rbar/sigma_x.^2 / (2*pi)/sigma_x/sigma_y %0.16f',-rbar/sigma_x.^2 / (2*pi)/sigma_x/sigma_y));
disp(sprintf(' %% dyydl %0.16f',dyydl));
disp(sprintf(' %% -rbar/sigma_y.^2 / (2*pi)/sigma_x/sigma_y %0.16f',-rbar/sigma_y.^2 / (2*pi)/sigma_x/sigma_y));







