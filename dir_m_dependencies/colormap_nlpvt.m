function output = colormap_nlpvt(n_c,step_threshold,gamma_r,gamma_g,gamma_b);
% threshold at step. ;
if (nargin<1); n_c = 64; end;
if (nargin<2); step_threshold = 1/9; end;
if (nargin<3); gamma_r = 0.50; end;
if (nargin<4); gamma_g = 1.25; end;
if (nargin<5); gamma_b = 0.25; end;

if n_c<0;
%%%%%%%%;
figure(1);clf;
gamma_r = 0.50;
n_g = 1+8;
gamma_g_ = linspace(0,2,n_g);
gamma_b_ = linspace(0,2,n_g);
n_c = 64;
p_ = 0:n_c-1;
x_ = [0+p_;1+p_;1+p_;0+p_;0+p_];
y_ = repmat([0;0;1;1;0],[1,n_c]);
for ng1=0:n_g-1; for ng2=0:n_g-1;
gamma_g = gamma_g_(1+ng2);
gamma_b = gamma_b_(1+ng1);
subplot(n_g,n_g,1+ng1+ng2*n_g);
c_ = reshape(colormap_nlpvt(n_c,1/9,gamma_r,gamma_g,gamma_b),[1,n_c,3]);
p_ = patch(x_,y_,c_); set(p_,'EdgeColor','none');
title(sprintf('%0.2f %0.2f %0.2f',gamma_r,gamma_g,gamma_b));
axisnotick;
end;end;%for ng1=0:n_g-1; for ng2=0:n_g-1;
sgtitle(sprintf('gamma_r %0.2f',gamma_r));
figbig;
%%%%%%%%;
figure(2);clf;
gamma_b = 0.50;
n_g = 1+8;
gamma_r_ = linspace(0,2,n_g);
gamma_g_ = linspace(0,2,n_g);
n_c = 64;
p_ = 0:n_c-1;
x_ = [0+p_;1+p_;1+p_;0+p_;0+p_];
y_ = repmat([0;0;1;1;0],[1,n_c]);
for ng1=0:n_g-1; for ng2=0:n_g-1;
gamma_r = gamma_r_(1+ng1);
gamma_g = gamma_g_(1+ng2);
subplot(n_g,n_g,1+ng1+ng2*n_g);
c_ = reshape(colormap_nlpvt(n_c,1/9,gamma_r,gamma_g,gamma_b),[1,n_c,3]);
p_ = patch(x_,y_,c_); set(p_,'EdgeColor','none');
title(sprintf('%0.2f %0.2f %0.2f',gamma_r,gamma_g,gamma_b));
axisnotick;
end;end;%for ng1=0:n_g-1; for ng2=0:n_g-1;
sgtitle(sprintf('gamma_b %0.2f',gamma_b));
figbig;
%%%%%%%%;
disp('returning'); return;
end;%if n_c<0;

tmp_r_ = transpose(linspace(1,0,n_c)).^gamma_r;
tmp_g_ = transpose(linspace(1,0,n_c)).^gamma_g;
tmp_b_ = transpose(linspace(1,0,n_c)).^gamma_b;
tmp_index_ = efind(linspace(0,1,n_c)<step_threshold);
tmp_r_(1+tmp_index_) = tmp_r_(1+tmp_index_).^4.00;
tmp_g_(1+tmp_index_) = tmp_g_(1+tmp_index_).^0.25;
tmp_b_(1+tmp_index_) = tmp_b_(1+tmp_index_).^0.25;
c_ = [ tmp_r_ , tmp_g_ , tmp_b_ ];
output = c_;
