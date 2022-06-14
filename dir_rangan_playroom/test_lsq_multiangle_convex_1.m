%%%%%%%%;
% is max-likelihood with angle-distributions per image convex? ;
%%%%%%%%;
pmxgf = @(a,b,x,m,s) exp( -(a.*x + b - m).^2 ./ (2.*s.^2) ) / sqrt(2*pi) ./ s;
pmgf = @(a,b,xp,xq,p,m,s) p.*pmxgf(a,b,xp,m,s) + (1-p).*pmxgf(a,b,xq,m,s);
n_a = 128;
n_b = 128;
vlim = 5;
a_ = linspace(-vlim,+vlim,n_a);
b_ = linspace(-vlim,+vlim,n_a);
[a__,b__] = ndgrid(a_,b_);
n_contour = 16;
n_image = 1;

rseed_ = 48 + (0:48-1); n_rseed = numel(rseed_);
figure(1);clf;figbig;fig80s;
p_row = 6; p_col = ceil(n_rseed/p_row); np=0;
for nrseed=0:n_rseed-1;
rseed = rseed_(1+nrseed);
rng(rseed);
s=1.0;
m_ = 2*rand(n_image,1)-1;
xp_ = 2*rand(n_image,1)-1;
xq_ = 2*rand(n_image,1)-1;
p_ = rand(n_image,1);
%%%%;
pmmgf__ = ones(n_a,n_b);
for nimage=0:n_image-1;
pmmgf__ = pmmgf__.*pmgf(a__,b__,xp_(1+nimage),xq_(1+nimage),p_(1+nimage),m_(1+nimage),s);
end;%for nimage=0:n_image-1;
%%%%;
nlpmmgf__ = -log(pmmgf__);
nlplim_ = prctile(nlpmmgf__,linspace(0,15,n_contour),'all');
subplot(p_row,p_col,1+nrseed);
contour(nlpmmgf__,nlplim_);
axis image; axisnotick;
drawnow();
end;%for nrseed=0:n_rseed-1;





