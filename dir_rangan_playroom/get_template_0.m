function ...
[ ...
 template__ ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_k_all_ ...
,n_viewing_all ...
,viewing_azimu_b_all_ ...
,viewing_polar_a_all_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,template_k_c_0__ ...
,template_k_c_1__ ...
,template_k_c_2__ ...
] = ...
get_template_0( ...
 verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_k_p_r_ ...
,l_max_ ...
,a_k_Y_ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,n_w_0in_ ...
);
% uses spherical-harmonic-expansion a_k_Y_ to evaluate templates on a collection of points on spherical shells determined by k_p_r_. ;
% ;
% inputs: ;
% ;
% verbose = integer verbosity_level. ;
% n_k_p_r = integer maximum number of shells. ;
% k_p_r_ = real array of length n_k_p_r; k_p_r_(nk_p_r) = k_p_r_value for shell nk_p_r. ;
% k_p_r_max = real maximum k-value. ;
% weight_k_p_r_ = real array of length n_k_p_r; radial quadrature weights. ;
% l_max_ = integer array of length n_k_p_r; l_max_(nk_p_r) = spherical harmonic order on shell nk_p_r; l_max_(nk_p_r) corresponds to n_lm_(nk_p_r) = (l_max_(nk_p_r)+1)^2 coefficients. ;
% a_k_Y_ = complex array of length \sum_{nk_p_r} (n_lm_(nk_p_r)+1)^2 ; coefficients are ordered in a row, with m varying quickly, l varying slowly and k varying most slowly. ;
% viewing_k_eq_d = real equatorial-distance used for sampling viewing angles and templates. ;
% template_k_eq_d = real equatorial-distance used for sampling viewing angles and templates. ;
% n_w_0in_ = integer array of length n_k_p_r; used if template_k_eq_d <=0; desired n_w_ for templates. ;
% ;
% outputs: ;
% ;
% template__ = complex array of templates for each viewing angle. ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

if nargin<8;

%%%%%%%%;
% First define integral of <f,f>. ;
%%%%%%%%;
h2d_ = @(kd) 4*pi^2*(besselj(0,kd) + besselj(2,kd)); % calculates <f_j,f_k>, normalized so that <f,f> = (4*pi^2);
dh2d_ = @(kd) 4*pi^3*(besselj(-1,kd) - besselj(+3,kd));
h3d_ = @(kd) 4*pi*( sin(kd) - (kd).*cos(kd) ) ./ kd.^3 ; % calculates <f_j,f_k>, normalized so that <f,f> = 4*pi/3;
dh3d_ = @(kd) 12*pi*( (kd.^2/3 - 1) .* sin(kd) + (kd).*cos(kd) ) ./ kd.^4 ;

%%%%%%%%;
% Now set up and test k-quadrature on sphere. ;
%%%%%%%%;
verbose=0; k_p_r_max = 9.0d0; k_eq_d = 1.0/32*k_p_r_max*0.95; TorL = 'L';
[n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_k_p_r_,k_c_0_all_,k_c_1_all_,k_c_2_all_] = sample_sphere_6(verbose,k_p_r_max,k_eq_d,'L') ;
delta_a_c_ = [+0.15;-0.25;+0.35];
delta_a_p_r = sqrt(delta_a_c_(1+0)^2 + delta_a_c_(1+1)^2 + delta_a_c_(1+2)^2);
%disp(sprintf(' %% 2*pi*k_p_r_max*delta_a_p_r = %0.2f',2*pi*k_p_r_max*delta_a_p_r));
a_k_p_form_ = exp(+i*2*pi*(k_c_0_all_*delta_a_c_(1+0) + k_c_1_all_*delta_a_c_(1+1) + k_c_2_all_*delta_a_c_(1+2)));
I_quad = sum(a_k_p_form_.*weight_k_all_);
I_form = h3d_(2*pi*k_p_r_max*fnorm(delta_a_c_))*k_p_r_max^3;
disp(sprintf(' %% I_form vs I_quad %0.16f',fnorm(I_form-I_quad)/fnorm(I_form)));

%%%%%%%%;
% Now convert to a_k_Y_ ; 
%%%%%%%%;
verbose=0;
l_max_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
l_max_(1+nk_p_r) = 1+ceil(2*pi*k_p_r_(1+nk_p_r));
%l_max_(1+nk_p_r) = 1+ceil(1*pi*k_p_r_(1+nk_p_r));
end;%for nk_p_r=0:n_k_p_r-1;
n_lm_ = (l_max_+1).^2;
n_lm_max = max(n_lm_);
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
l_max_max = max(l_max_); dWtdkd__l_max_max = 2*l_max_max;
m_max_ = -l_max_max : +l_max_max;
n_m_max = length(m_max_);
disp(sprintf(' %% n_k_p_r %d l_max_max %d n_lm_sum %d',n_k_p_r,l_max_max,n_lm_sum));
Y_l_val_ = zeros(n_lm_sum,1);
Y_m_val_ = zeros(n_lm_sum,1);
for nk_p_r=0:n_k_p_r-1;
l_max = l_max_(1+nk_p_r);
tmp_l_val_ = zeros(n_lm_(1+nk_p_r),1);
tmp_m_val_ = zeros(n_lm_(1+nk_p_r),1);
na=0; 
for l_val=0:l_max;
for m_val=-l_val:+l_val;
tmp_l_val_(1+na) = l_val;
tmp_m_val_(1+na) = m_val;
na=na+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max;
tmp_ij_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
Y_l_val_(1+tmp_ij_) = tmp_l_val_;
Y_m_val_(1+tmp_ij_) = tmp_m_val_;
end;%for nk_p_r=0:n_k_p_r-1;
weight_Y_ = zeros(n_lm_sum,1);
for nk_p_r=0:n_k_p_r-1;
tmp_ij_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
weight_Y_(1+tmp_ij_) = weight_k_p_r_(1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
tmp_t = tic;
[a_k_Y_quad_] = convert_k_p_to_spharm_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_k_p_r_,l_max_,a_k_p_form_);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_quad_ time %0.2fs',tmp_t));
disp(sprintf(' %% Here we should ensure that a_k_Y_ on the outer shells has decayed to the desired precision.'));
disp(sprintf(' %% Moreover, we should ensure that a_k_Y_(l,m) has decayed for large l,m at each shell.'));
flag_plot=0;
if flag_plot;
%imagesc_Y(k_p_r_max,n_k_p_r,k_p_r_,l_max_,log10(abs(a_k_Y_quad_)),[-10,0],colormap_beach());
subplot(1,2,1); plot(Y_l_val_,log10(abs(a_k_Y_quad_)),'.'); xlabel('l'); ylabel('log10(abs(a_k_Y_quad_))','Interpreter','none');
subplot(1,2,2); plot(Y_m_val_,log10(abs(a_k_Y_quad_)),'.'); xlabel('m'); ylabel('log10(abs(a_k_Y_quad_))','Interpreter','none');
end;%if flag_plot;
tmp_t = tic;
[a_k_p_quad_] = convert_spharm_to_k_p_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_k_p_r_,l_max_,a_k_Y_quad_);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_quad --> a_k_p_quad_ time %0.2fs',tmp_t));
disp(sprintf(' %% nufft3d3: a_k_p_quad error: %0.16f',fnorm(a_k_p_form_-a_k_p_quad_)/fnorm(a_k_p_form_)));
%%%%%%%%;
a_k_Y_quad__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_ij1_ = (0:n_lm_(1+nk_p_r)-1);
tmp_ij2_ = n_lm_csum_(1+nk_p_r) + tmp_ij1_;
a_k_Y_quad__(1+tmp_ij1_,1+nk_p_r) = a_k_Y_quad_(1+tmp_ij2_);
end;%for nk_p_r=0:n_k_p_r-1;

%%%%%%%%;
% Now generate templates. ;
%%%%%%%%;
verbose=1;
a_k_Y_ = a_k_Y_quad_;
template_k_eq_d = k_eq_d;
viewing_k_eq_d = 5*k_eq_d;
[template_quad__,n_w_,weight_2d_k_p_r_,weight_2d_k_all_,n_viewing_all,viewing_azimu_b_all_,viewing_polar_a_all_,n_viewing_polar_a,viewing_polar_a_,n_viewing_azimu_b_,template_k_c_0__,template_k_c_1__,template_k_c_2__] = get_template_0(verbose,n_k_p_r,k_p_r_,k_p_r_max,weight_k_p_r_,l_max_,a_k_Y_,viewing_k_eq_d,template_k_eq_d);
template_form__ = exp(+i*2*pi*(template_k_c_0__*delta_a_c_(1+0) + template_k_c_1__*delta_a_c_(1+1) + template_k_c_2__*delta_a_c_(1+2)));
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);
flag_plot=0;
if flag_plot;
nviewing_all = 12;
subplot(1,2,1); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(template_quad__(:,nviewing_all)),[-1,+1],colormap_beach()); title('quad');
subplot(1,2,2); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(template_form__(:,nviewing_all)),[-1,+1],colormap_beach()); title('form');
end;%if flag_plot;
disp(sprintf(' %% testing template evaluation: template_form__ vs template_quad__: %0.16f',fnorm(template_form__-template_quad__)/fnorm(template_form__)));
%%%%%%%%;
% Now test k-quadrature on disc. ;
%%%%%%%%;
template_2d_k_c_0_ = zeros(n_w_sum,1);
template_2d_k_c_1_ = zeros(n_w_sum,1);
template_2d_k_p_r_ = zeros(n_w_sum,1);
na=0;
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
for nw=0:n_w_(1+nk_p_r)-1;
gamma_z = 2*pi*nw/max(1,n_w_(1+nk_p_r));
template_2d_k_c_0_(1+na) = k_p_r*cos(gamma_z);
template_2d_k_c_1_(1+na) = k_p_r*sin(gamma_z);
template_2d_k_p_r_(1+na) = k_p_r;
na=na+1;
end;%for nw=0:n_w_(1+nk_p_r)-1;
end;%for nk_p_r=0:n_k_p_r-1;
assert(na==n_w_sum);
S_k_p_form_ = exp(+i*2*pi*(template_2d_k_c_0_*delta_a_c_(1+0) + template_2d_k_c_1_*delta_a_c_(1+1)));
%S_k_p_form_ = template_2d_k_p_r_.^3;
I_quad = sum(S_k_p_form_.*weight_2d_k_all_)*(2*pi)^2;
I_form = h2d_(2*pi*k_p_r_max*fnorm(delta_a_c_(1+(0:1))))/(2*pi)^2 * (pi*k_p_r_max^2);
%I_form = 2*pi*k_p_r_max^5/5;
disp(sprintf(' %% testing template quadrature: I_form vs I_quad: %0.16f',fnorm(I_form-I_quad)/fnorm(I_form)));

disp(sprintf(' %% returning')); return;
end;% if nargin<8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

n_lm_ = (l_max_+1).^2;
n_lm_max = max(n_lm_);
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
l_max_max = max(l_max_);
m_max_ = -l_max_max : +l_max_max;
n_m_max = length(m_max_);
if (verbose); disp(sprintf(' %% n_k_p_r %d l_max_max %d n_lm_sum %d',n_k_p_r,l_max_max,n_lm_sum)); end;
if (isempty(a_k_Y_)); a_k_Y_ = zeros(n_lm_sum,1); end;

%%%%%%%%;
if (verbose); disp(sprintf(' %% First determine the viewing angles on the outermost shell.')); end;
%%%%%%%%;
[n_viewing_all,viewing_azimu_b_all_,viewing_polar_a_all_,~,~,~,~,n_viewing_polar_a,viewing_polar_a_,n_viewing_azimu_b_] = sample_shell_5(k_p_r_max,viewing_k_eq_d,'L') ; %<-- obtain viewing angles on outer shell. ;
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
% spherical_harmonic_unphased_{1+nk_p_r}(1+l_max_(1+nk_p_r)+m_val,1+nw,1+nviewing_polar_a) ;
% contains the sum: ;
% \sum_{l_val=0}^{l_max_(1+nk_p_r)} ... ;
%             legendre_normalization_{1+nk_p_r}{1+l_val}(1+l_val+m_val) ... ;
%           * legendre_evaluate_{1+nk_p_r}{1+l_val}(1+l_val+m_val,1+nw,1+nviewing_polar_a) ... ;
%           * a_k_Y_(1+n_lm_csum_(1+nk_p_r)+(1+l_val-1)^2+l_val+m_val). ;
% Note that this is 'unphased', in the sense that the full evaluate requires: ;
% spherical_harmonic_evaluate__(1+n_w_csum_(1+nk_p_r)+nw,1+nviewing_all) = ... ;
% \sum_{m_val=-l_max_(1+nk_p_r)}^{+l_max_(1+nk_p_r)} ... ;
%             spherical_harmonic_unphased_{1+nk_p_r}(1+l_max_(1+nk_p_r)+m_val,1+nw,1+nviewing_polar_a) ... ;
%           * exp(+i*m_val*template_azimu_b__(1+n_w_csum_(1+nk_p_r)+nw,1+nviewing_all)). ;
% Note that this final exponential can be calculated as: ;
% expi_template_azimu_b__(1+n_w_csum_(1+nk_p_r)+nw,1+nviewing_all)^m_val. ;
%%%%%%%%;
spherical_harmonic_unphased_ = cell(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
spherical_harmonic_unphased_{1+nk_p_r} = zeros(1+2*l_max_(1+nk_p_r),n_w_(1+nk_p_r),n_viewing_polar_a);
l_max = l_max_(1+nk_p_r);
for m_val=-l_max:l_max;
for l_val=abs(m_val):l_max;
spherical_harmonic_unphased_{1+nk_p_r}(1+l_max_(1+nk_p_r)+m_val,:,:) = spherical_harmonic_unphased_{1+nk_p_r}(1+l_max_(1+nk_p_r)+m_val,:,:) + legendre_normalization_{1+nk_p_r}{1+l_val}(1+l_val+m_val) * legendre_evaluate_{1+nk_p_r}{1+l_val}(1+abs(m_val),:,:) * a_k_Y_(1+n_lm_csum_(1+nk_p_r)+(1+l_val-1)^2+l_val+m_val);
end;%for l_val=abs(m_val):l_max;
end;%for m_val=-l_max:l_max;
end;%for nk_p_r=0:n_k_p_r-1;

%%%%%%%%;
if (verbose); disp(sprintf(' %% now perform the final sum over m_val.')); end;
%%%%%%%%;
spherical_harmonic_evaluate__ = zeros(n_w_sum,n_viewing_all);
for nk_p_r=0:n_k_p_r-1;
tmp_spherical_harmonic_evaluate__ = zeros(n_w_(1+nk_p_r),n_viewing_all);
l_max = l_max_(1+nk_p_r);
nviewing_all=0;
for nviewing_polar_a=0:n_viewing_polar_a-1;
n_viewing_azimu_b = n_viewing_azimu_b_(1+nviewing_polar_a);
%viewing_polar_a = viewing_polar_a_(1+nviewing_polar_a);
for nviewing_azimu_b=0:n_viewing_azimu_b-1;
%viewing_azimu_b = 2*pi*nviewing_azimu_b/max(1,n_viewing_azimu_b);
tmp_expi_sub = expi_template_azimu_b__(1+n_w_csum_(1+nk_p_r)+(0:n_w_(1+nk_p_r)-1),1+nviewing_all);
tmp_expi_pre = tmp_expi_sub.^(-l_max);
tmp_expi_pos = tmp_expi_pre;
tmp_sum_ = zeros(n_w_(1+nk_p_r),1);
for m_val=-l_max:+l_max;
tmp_sum_ = tmp_sum_ + transpose(spherical_harmonic_unphased_{1+nk_p_r}(1+l_max+m_val,:,1+nviewing_polar_a)).*tmp_expi_pos;
tmp_expi_pos = tmp_expi_pos.*tmp_expi_sub;
end;%for m_val=-l_max:+l_max;
tmp_spherical_harmonic_evaluate__(:,1+nviewing_all) = tmp_sum_;
nviewing_all=nviewing_all+1;
end;%for nviewing_azimu_b=0:n_viewing_azimu_b-1;
end;%for nviewing_polar_a=0:n_viewing_polar_a-1;
spherical_harmonic_evaluate__(1+n_w_csum_(1+nk_p_r)+(0:n_w_(1+nk_p_r)-1),:) = tmp_spherical_harmonic_evaluate__;
end;%for nk_p_r=0:n_k_p_r-1;

template__ = spherical_harmonic_evaluate__;
