function ...
imagesc_Y( ...
 k_p_r_max ...
,n_r ...
,k_p_r_ ...
,l_max_ ...
,S_Y_ ...
,clim_ ...
,cmap_ ...
);
% spherical harmonic imagesc;
if (nargin<6); clim_ = []; end;
if (nargin<7); cmap_ = []; end;
if isempty(clim_); clim_ = mean(S_Y_) + 2.5*std(S_Y_,1)*[-1,+1]; end;
if isempty(cmap_); cmap_ = colormap_beach(); end;
n_c = size(cmap_,1);
if isempty(clim_); clim_ = mean(S_Y_) + std(S_Y_)*1.5*[-1,1]; end;%if isempty(clim_);
l_max_ = l_max_(1:n_r);
n_l = length(l_max_);
n_lm_ = (l_max_+1).^2;
n_lm_csum_ = cumsum([0;n_lm_(:)]);
n_lm_sum = sum(n_lm_);
X_ = zeros(5,n_lm_sum);
Y_ = zeros(5,n_lm_sum);
Z_ = zeros(5,n_lm_sum);
C_ = zeros(1,n_lm_sum,3);
ix=0;
for nr=0:n_r-1;
k_p_r = k_p_r_(1+nr);
if (nr==0); k_p_r_pre = 0.0d0; else k_p_r_pre = k_p_r_(1+nr-1); end;
if (nr==n_r-1); k_p_r_pos = k_p_r_max; else k_p_r_pos = k_p_r_(1+nr+1); end;
l_max_max = l_max_(1+nr);
for l_val=0:l_max_max;
n_m = 2*l_val+1;
for nm=0:n_m-1;
m_val = -l_val+nm;
nc = max(1,min(n_c,floor(n_c*(S_Y_(1+ix) - min(clim_))/diff(clim_))));
C_(1,1+ix,1:3) = cmap_(nc,:);
l0_(1+ix) = l_val-0.5; l1_(1+ix) = l_val+0.5;
m0_(1+ix) = m_val-0.5; m1_(1+ix) = m_val+0.5;
r0_(1+ix) = k_p_r_pre; r1_(1+ix) = k_p_r_pos; rh_(1+ix) = k_p_r;
ix=ix+1;
end;%for nm=0:n_m-1;
end;%for l_val=0:l_max_max;
end;%for nr=0:n_r-1;
assert(ix==n_lm_sum);
X_ = [l0_ ; l0_ ; l1_ ; l1_ ; l0_ ];
Y_ = [m0_ ; m1_ ; m1_ ; m0_ ; m0_ ];
%Z_ = [r0_ ; r0_ ; r1_ ; r1_ ; r0_ ];
Z_ = [rh_ ; rh_ ; rh_ ; rh_ ; rh_ ];
hold on;
for nr=0:n_r-1;
ij_ = n_lm_csum_(1+nr) + (1:n_lm_(1+nr));
%p=patch(X_(:,ij_),Y_(:,ij_),C_(1,ij_,:)); set(p,'EdgeColor','none');
p=patch(X_(:,ij_),Y_(:,ij_),Z_(:,ij_),C_(1,ij_,:)); set(p,'EdgeColor','none');
end;%for nr=0:n_r-1;
hold off;
