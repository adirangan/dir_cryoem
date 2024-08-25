clear;
nf=0;
n_x = 128+1;
n_y = 128+1;
x_ = transpose(linspace(-2,+2,n_x));
dx = mean(diff(x_));
y_ = transpose(linspace(-2,+2,n_y));
dy = mean(diff(y_));
[x__,y__] = ndgrid(x_,y_);
sigma_x = 0.150;
sigma_y = 0.325;
rbar = 7.5;
nll = @(x,y) 0.5 * (x.*x/sigma_x^2 + y.*y/sigma_y^2) ;
l = @(x,y) rbar * 1/(2*pi)/sigma_x/sigma_y .* exp(-nll(x,y)) ;
figure(1+nf);nf=nf+1;clf;fig80s;
imagesc(l(x__,y__));
axis image; axisnotick; xlabel('x'); ylabel('y');
avg = sum((x__.^2 + y__.^2).^0.*l(x__,y__),'all')*dx*dy;
disp(sprintf(' %% avg %0.16f',avg));
disp(sprintf(' %% rbar %0.16f',rbar));
var = sum((x__.^2 + y__.^2).^1.*l(x__,y__),'all')*dx*dy;
disp(sprintf(' %% var %0.16f',var));
disp(sprintf(' %% rbar^1*(sigma_x.^2 + sigma_y.^2) %0.16f',rbar^1*(sigma_x.^2 + sigma_y.^2)));
nx_mid = ceil(n_x/2);
ny_mid = ceil(n_y/2);
l__ = l(x__,y__);
dxxdl = (l__(nx_mid+1,ny_mid) - 2*l__(nx_mid+0,ny_mid) + l__(nx_mid-1,ny_mid))/max(1e-12,dx*dx);
dyydl = (l__(nx_mid,ny_mid+1) - 2*l__(nx_mid,ny_mid+0) + l__(nx_mid,ny_mid-1))/max(1e-12,dy*dy);
disp(sprintf(' %% dxxdl %0.16f',dxxdl));
disp(sprintf(' %% -rbar/sigma_x.^2 / (2*pi)/sigma_x/sigma_y %0.16f',-rbar/sigma_x.^2 / (2*pi)/sigma_x/sigma_y));
disp(sprintf(' %% dyydl %0.16f',dyydl));
disp(sprintf(' %% -rbar/sigma_y.^2 / (2*pi)/sigma_x/sigma_y %0.16f',-rbar/sigma_y.^2 / (2*pi)/sigma_x/sigma_y));







