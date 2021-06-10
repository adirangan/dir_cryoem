function ...
[ ...
 template_wkSa___ ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_k_all_ ...
,n_viewing_all ...
,viewing_azimu_b_all_ ...
,viewing_polar_a_all_ ...
,viewing_weight_all_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,template_k_c_0__ ...
,template_k_c_1__ ...
,template_k_c_2__ ...
] = ...
get_template__1( ...
 verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_k_p_r_ ...
,l_max_ ...
,n_a ...
,a_k_Y_ya__ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,n_w_0in_ ...
);
% uses spherical-harmonic-expansions a_k_Y_ya__ to evaluate templates on a collection of points on spherical shells determined by k_p_r_. ;
% ;
% inputs: ;
% ;
% verbose = integer verbosity_level. ;
% n_k_p_r = integer maximum number of shells. ;
% k_p_r_ = real array of length n_k_p_r; k_p_r_(nk_p_r) = k_p_r_value for shell nk_p_r. ;
% k_p_r_max = real maximum k-value. ;
% weight_k_p_r_ = real array of length n_k_p_r; radial quadrature weights. ;
% l_max_ = integer array of length n_k_p_r; l_max_(nk_p_r) = spherical harmonic order on shell nk_p_r; l_max_(nk_p_r) corresponds to n_lm_(nk_p_r) = (l_max_(nk_p_r)+1)^2 coefficients. ;
% n_a = integer number of a_k_Y_. ;
% a_k_Y_ya__ = complex array of size ( \sum_{nk_p_r} (n_lm_(nk_p_r)+1)^2 , n_a ). ;
%              each a_k_Y_ = a_k_Y_ya__(:,1+na) is an expansion, ;
%              coefficients are ordered in a row, with m varying quickly, l varying slowly and k varying most slowly. ;
% viewing_k_eq_d = real equatorial-distance used for sampling viewing angles and templates. ;
% template_k_eq_d = real equatorial-distance used for sampling inplane-shifts along each template. ;
% n_w_0in_ = integer array of length n_k_p_r; used if template_k_eq_d <=0; desired n_w_ for templates. ;
% ;
% outputs: ;
% ;
% template_wkSa___ = complex array of templates for each viewing angle and expansion. ;
%                    each template_ = template_wkSa___(:,1+nS,1+na) is a template, ;
%                    coefficients are ordered in a row, with nw varying quickly and nk varying slowly. ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

if nargin<10;
verbose=0;
n_k_p_r = 8;
k_p_r_ = transpose(1:8);
k_p_r_max = 9;
weight_k_p_r_ = transpose(linspace(2,3,n_k_p_r));
l_max_ = transpose(2:9);
n_lm_sum = sum((1+l_max_).^2);
n_a = 7;
a_k_Y_ya__ = crandn(n_lm_sum,n_a);
viewing_k_eq_d = 1;
template_k_eq_d = 1;
%%%%%%%%;
tmp_t = tic();
[ ...
 template_A_wkSa___ ...
] = ...
get_template__1( ...
 verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_k_p_r_ ...
,l_max_ ...
,n_a ...
,a_k_Y_ya__ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,[] ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% together: %0.2fs',tmp_t));
%%%%%%%%;
template_B_wkSa___ = zeros(size(template_A_wkSa___));
tmp_t=tic();
for na=0:n_a-1;
[ ...
 template_B_wkSa___(:,:,1+na) ...
] = ...
get_template_1( ...
 verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_k_p_r_ ...
,l_max_ ...
,a_k_Y_ya__(:,1+na) ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,[] ...
);
end;%for na=0:n_a-1;
tmp_t = toc(tmp_t); disp(sprintf(' %% individual: %0.2fs',tmp_t));
%%%%%%%%;
disp(sprintf(' %% individual vs together: %0.16f', fnorm(template_B_wkSa___ - template_A_wkSa___)/fnorm(template_B_wkSa___)));
disp(sprintf(' %% returning')); return;
end;% if nargin<11;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

n_lm_ = (l_max_+1).^2;
n_lm_max = max(n_lm_);
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
l_max_max = max(l_max_);
m_max_ = -l_max_max : +l_max_max;
n_m_max = length(m_max_);
if (verbose); disp(sprintf(' %% n_k_p_r %d l_max_max %d n_lm_sum %d',n_k_p_r,l_max_max,n_lm_sum)); end;
if (isempty(a_k_Y_ya__)); a_k_Y_ya__ = zeros(n_lm_sum,n_a); end;

%%%%%%%%;
if (verbose); disp(sprintf(' %% First determine the viewing angles on the outermost shell.')); end;
%%%%%%%%;
[ ...
 n_viewing_all ...
,viewing_azimu_b_all_ ...
,viewing_polar_a_all_ ...
,viewing_weight_all_ ...
,~ ...
,~ ...
,~ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
] = ...
sample_shell_5( ...
 k_p_r_max ...
,viewing_k_eq_d ...
,'L' ...
) ; %<-- obtain viewing angles on outer shell. ;
n_viewing_azimu_b_sum = sum(n_viewing_azimu_b_);
n_viewing_azimu_b_csum_ = cumsum([0;n_viewing_azimu_b_]);
if (verbose); disp(sprintf(' %% n_viewing_all %d n_viewing_polar_a %d n_viewing_azimu_b_sum %d',n_viewing_all,n_viewing_polar_a,n_viewing_azimu_b_sum)); end;

%%%%%%%%;
if (verbose); disp(sprintf(' %% Now determine the points along each equatorial plane (i.e., the points for each template).')); end;
% Note that the radial arrangement of these points is determined by the radii of the shells (i..e, the input k_p_r_). ;
%%%%%%%%;
n_w_ = zeros(n_k_p_r,1);
if (template_k_eq_d>0);
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
n_equator = 3+round(2*pi*k_p_r/template_k_eq_d);
n_polar_a = 3+round(n_equator/2);
n_w_(1+nk_p_r) = 2*n_polar_a;
end;%for nk_p_r=0:n_k_p_r-1;
end;%if (template_k_eq_d>0);
if (template_k_eq_d<=0);
n_w_ = n_w_0in_;
assert(numel(n_w_)==n_k_p_r); assert(min(n_w_)>0);
end;%if (template_k_eq_d<=0);
n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);
if (verbose); disp(sprintf(' %% n_w_max %d n_w_sum %d',n_w_max,n_w_sum)); end;
%%%%%%%%;
if (verbose); disp(sprintf(' %% Set up integration weights for the templates.')); end;
%%%%%%%%;
tmp_P_ = zeros(n_k_p_r,n_k_p_r); %<-- polynomials of order 0:n_k_p_r-1 evaluated on k_p_r_/k_p_r_max. ;
tmp_I_ = zeros(n_k_p_r,1); %<-- integrals of those polynomials on the 2d-disc of radius 1. ;
for nk_p_r=0:n_k_p_r-1;
tmp_x = @(x) x.^nk_p_r;
tmp_P_(1+nk_p_r,:) = tmp_x(k_p_r_/k_p_r_max);
tmp_I_(1+nk_p_r) = 2*pi*1/(nk_p_r+2);
end;%for nk_p_r=0:n_k_p_r-1;
tmp_W_ = pinv(tmp_P_,1e-6)*tmp_I_;
if (verbose>1); disp(sprintf(' %% weight error: %0.16f',fnorm(tmp_P_*tmp_W_ - tmp_I_)/fnorm(tmp_I_))); end;
weight_2d_k_p_r_ = tmp_W_*k_p_r_max^2;
weight_2d_k_all_ = zeros(n_w_sum,1);
for nk_p_r=0:n_k_p_r-1;
tmp_ij_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
weight_2d_k_all_(1+tmp_ij_) = weight_2d_k_p_r_(1+nk_p_r) / max(1,n_w_(1+nk_p_r)) / (2*pi)^2;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
% Set up inner gamma_z for the templates. ;
%%%%%%%%;
gamma_z_all_ = zeros(n_w_sum,1);
for nk_p_r=0:n_k_p_r-1;
tmp_ij_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
tmp_val_ = transpose(linspace(0,2*pi,n_w_(1+nk_p_r)+1));
gamma_z_all_(1+tmp_ij_) = tmp_val_(1:n_w_(1+nk_p_r));
end;%for nk_p_r=0:n_k_p_r-1;
cc_ = cos(gamma_z_all_);
sc_ = sin(gamma_z_all_);

%%%%%%%%;
% The general formula used here is as follows. ;
% let sa and ca be sin(polar_a) and cos(polar_a), respectively. ;
% let sb and cb be sin(azimu_b) and cos(azimu_b), respectively. ;
% let sc and cc be sin(gamma_z) and cos(gamma_z), respectively. ;
% And rotation by azimu_b about the +z-axis is represented as: ;
% Rz(azimu_b) = ;
% [ +cb -sb 0 ] ;
% [ +sb +cb 0 ] ;
% [  0   0  1 ] ;
% And rotation by polar_a about the +y-axis is represented as: ;
% Ry(polar_a) = ;
% [ +ca 0 +sa ] ;
% [  0  1  0  ] ;
% [ -sa 0 +ca ] ;
% And rotation by gamma_z about the +z-axis is represented as: ;
% Rz(gamma_z) = ;
% [ +cc -sc 0 ] ;
% [ +sc +cc 0 ] ;
% [  0   0  1 ] ;
% Which, collectively, implies that under the transform: ;
% Rz(azimu_b) * Ry(polar_a) * Rz(gamma_z), ;
% Which is the same as: ;
% [ +cb -sb 0 ] [ +ca*cc -ca*sc +sa ]   [ +cb*ca*cc - sb*sc , -cb*ca*sc -sb*cc , +cb*sa ];
% [ +sb +cb 0 ] [ +sc    +cc    0   ] = [ +sb*ca*cc + cb*sc , -sb*ca*sc +cb*cc , +sb*sa ];
% [  0   0  1 ] [ -sa*cc +sa*sc +ca ]   [ -sa*cc            , +sa*sc           , +ca    ];
% the point [1;0;0] is mapped to: ;
% [ template_k_c_0 ; template_k_c_1 ; template_k_c_2 ] = [ +cb*ca*cc - sb*sc ; +sb*ca*cc + cb*sc ; -sa*cc ];
%%%%%%%%;

if (verbose); disp(sprintf(' %% template: (%d,%d)=%d (%0.2f GB)',n_w_sum,n_viewing_all,n_w_sum*n_viewing_all,n_w_sum*n_viewing_all*16/1e9)); end;

%%%%%%%%;
if (verbose); disp(sprintf(' %% Now construct array of k_c_?_ values for the templates.')); end;
%%%%%%%%;
template_k_p_r_ = zeros(n_w_sum,1);
for nk_p_r=0:n_k_p_r-1;
tmp_ij_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
template_k_p_r_(1+tmp_ij_) = k_p_r_(1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
template_k_c_0__ = zeros(n_w_sum,n_viewing_all);
template_k_c_1__ = zeros(n_w_sum,n_viewing_all);
template_k_c_2__ = zeros(n_w_sum,n_viewing_all);
for nviewing_all=0:n_viewing_all-1;
viewing_polar_a = viewing_polar_a_all_(1+nviewing_all); ca = cos(viewing_polar_a); sa = sin(viewing_polar_a);
viewing_azimu_b = viewing_azimu_b_all_(1+nviewing_all); cb = cos(viewing_azimu_b); sb = sin(viewing_azimu_b);
template_k_c_0__(:,1+nviewing_all) = (+cb*ca*cc_ - sb*sc_).*template_k_p_r_;
template_k_c_1__(:,1+nviewing_all) = (+sb*ca*cc_ + cb*sc_).*template_k_p_r_;
template_k_c_2__(:,1+nviewing_all) = (-sa*cc_            ).*template_k_p_r_;
end;%for nviewing_all=0:n_viewing_all-1;
template_azimu_b__ = atan2(template_k_c_1__,template_k_c_0__);
expi_template_azimu_b__ = exp(i*template_azimu_b__);
clear template_azimu_b__;
if nargout<11; clear template_k_c_0__ template_k_c_1__ template_k_c_2__ ; end;
%%%%%%%%;
% We also use a condensed array, called condense_k_c_2__, which only depends on the polar_a, and not on azimu_b. ;
%%%%%%%%;
condense_k_c_2__ = zeros(n_w_sum,n_viewing_polar_a);
for nviewing_polar_a=0:n_viewing_polar_a-1;
viewing_polar_a = viewing_polar_a_(1+nviewing_polar_a); ca = cos(viewing_polar_a); sa = sin(viewing_polar_a);
condense_k_c_2__(:,1+nviewing_polar_a) = -sa*cc_ ;
end;%for nviewing_polar_a=0:n_viewing_all-1;

%%%%%%%%;
if (verbose); disp(sprintf(' %% Now evaluate associated legendre polynomials at the varous k_c_2 values.')); end;
% Here legendre_evaluate_{1+nk_p_r}{1+l_val}(1+l_val+m_val,1+nw,1+nviewing_polar_a) contains ;
% the associated legendre-function of degree l_val and order abs(m_val) (ranging from 0 to +l_val) ;
% evaluated at the k_c_2 value stored in condense_k_c_2__(1+n_w_csum_(1+nk_p_r)+nw,1+nviewing_polar_a). ;
% Note that this is associated with ring/shell nk_p_r and viewing_polar_a_(1+nviewing_polar_a). ;
% The legendre_normalization_{1+nk_p_r}{1+l_val}(1+abs(m_val)) contains ;
% The normalization coefficient for the spherical harmonics associated with l_val and m_val. ;
% Note that this is somewhat redundant (as it does not depend explicitly on the shell). ;
%%%%%%%%;
legendre_evaluate_ = cell(n_k_p_r,1);
legendre_normalization_ = cell(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
l_max = l_max_(1+nk_p_r);
template_ij_ = n_w_csum_(1+nk_p_r) + (1:n_w_(1+nk_p_r));
legendre_evaluate_{1+nk_p_r} = cell(l_max+1,1);
for l_val=0:l_max;
tmp_P__ = zeros(1+1*l_val,n_w_(1+nk_p_r),n_viewing_polar_a);
tmp_a1 = ((1+2*l_val)/(4*pi));
m_val_ = -l_val:+l_val;
tmp_a2_ = exp(lfactorial(l_val-abs(m_val_)) - lfactorial(l_val+abs(m_val_)));
tmp_a3_ = sqrt(tmp_a1*tmp_a2_);
tmp_P__ = legendre(l_val,condense_k_c_2__(template_ij_,:),'unnorm');
tmp_t = tic;
legendre_evaluate_{1+nk_p_r}{1+l_val} = reshape(tmp_P__,[1+1*l_val,n_w_(1+nk_p_r),n_viewing_polar_a]);
tmp_t = toc(tmp_t);
if (verbose>1); disp(sprintf(' %% nk_p_r %d/%d l_val %d/%d legendre_evaluate(%d,%d) %0.2fs',nk_p_r,n_k_p_r,l_val,l_max,n_w_(1+nk_p_r),n_viewing_polar_a,tmp_t)); end;
legendre_normalization_{1+nk_p_r}{1+l_val} = tmp_a3_;
end;%for l_val=0:l_max;
end;%for nk_p_r=0:n_k_p_r-1;

%%%%%%%%;
if (verbose); disp(sprintf(' %% Now accumulate the legendre_evaluates over l_val, for each m_val.')); end;
% We account for the normalization coefficients here, ;
% so that later we can apply the complex exponential to produce the spherical-harmonics. ;
% More specifically: ;
% spherical_harmonic_unphased_kmwpa_____{1+nk_p_r}(1+l_max_(1+nk_p_r)+m_val,1+nw,1+nviewing_polar_a,1+na) ;
% contains the sum: ;
% \sum_{l_val=0}^{l_max_(1+nk_p_r)} ... ;
%             legendre_normalization_{1+nk_p_r}{1+l_val}(1+l_val+m_val) ... ;
%           * legendre_evaluate_{1+nk_p_r}{1+l_val}(1+l_val+m_val,1+nw,1+nviewing_polar_a) ... ;
%           * a_k_Y_ya__(1+n_lm_csum_(1+nk_p_r)+(1+l_val-1)^2+l_val+m_val,1+na). ;
% Note that this is 'unphased', in the sense that the full evaluate requires: ;
% spherical_harmonic_evaluate_kwSa___(1+n_w_csum_(1+nk_p_r)+nw,1+nviewing_all,1+na) = ... ;
% \sum_{m_val=-l_max_(1+nk_p_r)}^{+l_max_(1+nk_p_r)} ... ;
%             spherical_harmonic_unphased_kmwpa_____{1+nk_p_r}(1+l_max_(1+nk_p_r)+m_val,1+nw,1+nviewing_polar_a,1+na) ... ;
%           * exp(+i*m_val*template_azimu_b__(1+n_w_csum_(1+nk_p_r)+nw,1+nviewing_all)). ;
% Note that this final exponential can be calculated as: ;
% expi_template_azimu_b__(1+n_w_csum_(1+nk_p_r)+nw,1+nviewing_all)^m_val. ;
%%%%%%%%;
spherical_harmonic_unphased_kmwpa_____ = cell(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
n_w = n_w_(1+nk_p_r);
spherical_harmonic_unphased_kmwpa_____{1+nk_p_r} = zeros(1+2*l_max_(1+nk_p_r),n_w,n_viewing_polar_a,n_a);
l_max = l_max_(1+nk_p_r);
for m_val=-l_max:l_max;
for l_val=abs(m_val):l_max;
spherical_harmonic_unphased_kmwpa_____{1+nk_p_r}(1+l_max_(1+nk_p_r)+m_val,:,:,:) = spherical_harmonic_unphased_kmwpa_____{1+nk_p_r}(1+l_max_(1+nk_p_r)+m_val,:,:,:) + legendre_normalization_{1+nk_p_r}{1+l_val}(1+l_val+m_val) * bsxfun(@times,reshape(legendre_evaluate_{1+nk_p_r}{1+l_val}(1+abs(m_val),:,:),[1,n_w,n_viewing_polar_a,1]),reshape(a_k_Y_ya__(1+n_lm_csum_(1+nk_p_r)+(1+l_val-1)^2+l_val+m_val,:),[1,1,1,n_a]));
end;%for l_val=abs(m_val):l_max;
end;%for m_val=-l_max:l_max;
end;%for nk_p_r=0:n_k_p_r-1;

%%%%%%%%;
if (verbose); disp(sprintf(' %% now perform the final sum over m_val.')); end;
%%%%%%%%;
spherical_harmonic_evaluate_wkSa___ = zeros(n_w_sum,n_viewing_all,n_a);
for nk_p_r=0:n_k_p_r-1;
n_w = n_w_(1+nk_p_r);
tmp_spherical_harmonic_evaluate_wSa___ = zeros(n_w,n_viewing_all,n_a);
l_max = l_max_(1+nk_p_r);
nviewing_all=0;
for nviewing_polar_a=0:n_viewing_polar_a-1;
n_viewing_azimu_b = n_viewing_azimu_b_(1+nviewing_polar_a);
%viewing_polar_a = viewing_polar_a_(1+nviewing_polar_a);
for nviewing_azimu_b=0:n_viewing_azimu_b-1;
%viewing_azimu_b = 2*pi*nviewing_azimu_b/max(1,n_viewing_azimu_b);
tmp_expi_sub = expi_template_azimu_b__(1+n_w_csum_(1+nk_p_r)+(0:n_w-1),1+nviewing_all);
tmp_expi_pre = tmp_expi_sub.^(-l_max);
tmp_expi_pos = tmp_expi_pre;
tmp_sum_wa__ = zeros(n_w,n_a);
for m_val=-l_max:+l_max;
tmp_sum_wa__ = tmp_sum_wa__ + bsxfun(@times,reshape(spherical_harmonic_unphased_kmwpa_____{1+nk_p_r}(1+l_max+m_val,:,1+nviewing_polar_a,:),[n_w,n_a]),tmp_expi_pos);
tmp_expi_pos = tmp_expi_pos.*tmp_expi_sub;
end;%for m_val=-l_max:+l_max;
tmp_spherical_harmonic_evaluate_wSa___(:,1+nviewing_all,:) = tmp_sum_wa__;
nviewing_all=nviewing_all+1;
end;%for nviewing_azimu_b=0:n_viewing_azimu_b-1;
end;%for nviewing_polar_a=0:n_viewing_polar_a-1;
spherical_harmonic_evaluate_wkSa___(1+n_w_csum_(1+nk_p_r)+(0:n_w-1),:,:) = tmp_spherical_harmonic_evaluate_wSa___;
end;%for nk_p_r=0:n_k_p_r-1;

template_wkSa___ = spherical_harmonic_evaluate_wkSa___;
