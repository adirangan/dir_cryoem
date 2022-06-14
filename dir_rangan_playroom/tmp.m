%%%%%%%%;
% is max-likelihood with angle-distributions per image convex? ;
%%%%%%%%;
pmxgf = @(a,b,x,m,s) exp( -(a.*x + b - m).^2 ./ (2.*s.^2) ) / sqrt(2*pi) ./ s;
pmgf = @(a,b,xp,xq,p,m,s) p.*pmxgf(a,b,xp,m,s) + (1-p).*pmxgf(a,b,xq,m,s);
pm1m2gf = @(a,b,x1p,x1q,p1,m1,x2p,x2q,p2,m2,s) pmgf(a,b,x1p,x1q,p1,m1,s).*pmgf(a,b,x2p,x2q,p2,m2,s);
pm1m2m3gf = @(a,b,x1p,x1q,p1,m1,x2p,x2q,p2,m2,s) pmgf(a,b,x1p,x1q,p1,m1,s).*pmgf(a,b,x2p,x2q,p2,m2,s);
n_a = 128;
n_b = 128;
a_ = linspace(-5,+5,n_a);
b_ = linspace(-5,+5,n_a);
[a__,b__] = ndgrid(a_,b_);
n_contour = 32;

rseed_ = 48 + (0:48-1); n_rseed = numel(rseed_);
figure(1);clf;figbig;fig80s;
p_row = 6; p_col = ceil(n_rseed/p_row); np=0;
for nrseed=0:n_rseed-1;
rseed = rseed_(1+nrseed);
rng(rseed);
s=1.0;
m1 = 2*rand()-1; m2 = 2*rand()-1; m3 = 2*rand()-1;
x1p = 2*rand()-1; x1q = 2*rand()-1;
x2p = 2*rand()-1; x2q = 2*rand()-1;
x3p = 2*rand()-1; x3q = 2*rand()-1;
p1 = rand(); p2 = rand(); p3 = rand();
pm1m2gf__ = pm1m2gf(a__,b__,x1p,x1q,p1,m1,x2p,x2q,p1,m2,s);
nlpm1m2gf__ = -log(pm1m2gf__);
nlplim_ = prctile(nlpm1m2gf__,linspace(0,50,n_contour),'all');
subplot(p_row,p_col,1+nrseed);
contour(nlpm1m2gf__,nlplim_);
axis image; axisnotick;
drawnow();
end;%for nrseed=0:n_rseed-1;





