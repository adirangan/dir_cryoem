function ...
[ ...
 template_waS___ ...
,n_w ...
,n_viewing_all ...
,viewing_azimu_b_all_ ...
,viewing_polar_a_all_ ...
] = ...
pm_template_2( ...
 verbose ...
,l_max ...
,n_a ...
,a_k_Y_ya__ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,n_w_0in ...
);
% uses spherical-harmonic-expansions a_k_Y_ya__ to evaluate templates on a collection of points on spherical shells. ;
% each spherical-shell has the same resolution, determiend by viewing_k_eq_d and template_k_eq_d and/or n_w_max. ;
% ;
% inputs: ;
% ;
% verbose = integer verbosity_level. ;
% l_max = spherical harmonic order on each shell. ;
%         corresponds to n_lm = (l_max+1)^2 coefficients. ;
% n_a = integer number of shells. ;
% a_k_Y_ya__ = complex array of size (n_lm,n_a). ;
%              coefficients are ordered in a row, with m varying quickly, l varying slowly. ;
% viewing_k_eq_d = real equatorial-distance used for sampling viewing angles and templates. ;
% template_k_eq_d = real equatorial-distance used for sampling inplane-shifts along each template. ;
% n_w_0in = integer. used if template_k_eq_d <=0; desired n_w for templates. ;
% ;
% outputs: ;
% ;
% template_waS___ = complex array of templates for each viewing angle. ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

if nargin<7;
verbose=0;
n_k_p_r = 18;
k_p_r_max = 1;
k_p_r_ = ones(n_k_p_r,1);
weight_k_p_r_ = ones(n_k_p_r,1);
l_max = 9;
n_lm = (l_max+1)^2;
l_max_ = l_max*ones(n_k_p_r,1);
n_lm_sum = sum((1+l_max_).^2);
a_k_Y_ = crandn(n_lm_sum,1);
viewing_k_eq_d = 1/pi;
template_k_eq_d = -1;
n_w_max = 32;
%%%%%%%%;
tmp_t = tic();
[ ...
 template_waS___ ...
,~ ...
,n_S ...
,~ ...
,~ ...
] = ...
pm_template_2( ...
 verbose ...
,l_max ...
,n_k_p_r ...
,reshape(a_k_Y_,[n_lm,n_k_p_r]) ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,n_w_max ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% pm_template_2: %0.2fs',tmp_t));
template_waS__ = reshape(template_waS___,[n_w_max*n_k_p_r,n_S]);
%%%%%%%%;
tmp_t = tic();
[ ...
 template_wkS__ ...
] = ...
get_template_1( ...
 verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_k_p_r_ ...
,l_max_ ...
,a_k_Y_ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,n_w_max*ones(n_k_p_r,1) ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% get_template_1: %0.2fs',tmp_t));
%%%%%%%%;
disp(sprintf(' %% template_wkS__ vs template_waS__: %0.16f', fnorm(template_wkS__ - template_waS__)/fnorm(template_wkS__)));
disp(sprintf(' %% returning')); return;
end;% if nargin<7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

n_lm = (l_max+1)^2;
m_max_ = -l_max : +l_max;
n_m_max = length(m_max_); %<-- 2*l_max+1;
if (verbose); disp(sprintf(' %% l_max %d n_lm %d',l_max,n_lm)); end;
if (isempty(a_k_Y_ya__)); a_k_Y_ya__ = zeros(n_lm,1); end;

%%%%%%%%;
if (verbose); disp(sprintf(' %% First determine the viewing angles.')); end;
%%%%%%%%;
k_p_r = 1;
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
 k_p_r ...
,viewing_k_eq_d ...
,'L' ...
) ; %<-- obtain viewing angles. ;
n_viewing_azimu_b_sum = sum(n_viewing_azimu_b_);
n_viewing_azimu_b_csum_ = cumsum([0;n_viewing_azimu_b_]);
if (verbose); disp(sprintf(' %% n_viewing_all %d n_viewing_polar_a %d n_viewing_azimu_b_sum %d',n_viewing_all,n_viewing_polar_a,n_viewing_azimu_b_sum)); end;

%%%%%%%%;
if (verbose); disp(sprintf(' %% Now determine the points along each equatorial plane (i.e., the points for each template).')); end;
%%%%%%%%;
n_w = 0;
if (template_k_eq_d>0);
k_p_r = 1;
n_equator = 3+round(2*pi*k_p_r/template_k_eq_d);
n_polar_a = 3+round(n_equator/2);
n_w = 2*n_polar_a;
end;%if (template_k_eq_d>0);
if (template_k_eq_d<=0);
n_w = max(6,n_w_0in);
end;%if (template_k_eq_d<=0);
if (verbose); disp(sprintf(' %% n_w %d',n_w)); end;
%%%%%%%%;
% Set up inner gamma_z for the templates. ;
%%%%%%%%;
gamma_z_ = zeros(n_w,1); gamma_z_ = transpose(linspace(0,2*pi,n_w+1)); gamma_z_ = gamma_z_(1:n_w);
cc_ = cos(gamma_z_); sc_ = sin(gamma_z_);

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

if (verbose); disp(sprintf(' %% template: (%d,%d,%d)=%d (%0.2f GB)',n_w,n_viewing_all,n_a,n_w*n_viewing_all*n_a,n_w*n_viewing_all*16/1e9)); end;

%%%%%%%%;
if (verbose); disp(sprintf(' %% Now construct array of k_c_?_ values for the templates.')); end;
%%%%%%%%;
template_k_p_r = 1;
%%%%%%%%;
template_k_c_0__ = zeros(n_w,n_viewing_all);
template_k_c_1__ = zeros(n_w,n_viewing_all);
template_k_c_2__ = zeros(n_w,n_viewing_all);
if (verbose); disp(sprintf(' %% template_k_c___: (%d,%d)=%d (%0.2f GB)',n_w,n_viewing_all,n_w*n_viewing_all,n_w*n_viewing_all*8/1e9)); end;
for nviewing_all=0:n_viewing_all-1;
viewing_polar_a = viewing_polar_a_all_(1+nviewing_all); ca = cos(viewing_polar_a); sa = sin(viewing_polar_a);
viewing_azimu_b = viewing_azimu_b_all_(1+nviewing_all); cb = cos(viewing_azimu_b); sb = sin(viewing_azimu_b);
template_k_c_0__(:,1+nviewing_all) = (+cb*ca*cc_ - sb*sc_).*template_k_p_r;
template_k_c_1__(:,1+nviewing_all) = (+sb*ca*cc_ + cb*sc_).*template_k_p_r;
template_k_c_2__(:,1+nviewing_all) = (-sa*cc_            ).*template_k_p_r;
end;%for nviewing_all=0:n_viewing_all-1;
template_azimu_b__ = atan2(template_k_c_1__,template_k_c_0__);
expi_template_azimu_b__ = exp(i*template_azimu_b__);
clear template_azimu_b__;
clear template_k_c_0__ template_k_c_1__ template_k_c_2__ ;
%%%%%%%%;
% We also use a condensed array, called condense_k_c_2__, which only depends on the polar_a, and not on azimu_b. ;
%%%%%%%%;
condense_k_c_2__ = zeros(n_w,n_viewing_polar_a);
if (verbose); disp(sprintf(' %% condense_k_c_2__: (%d,%d)=%d (%0.2f GB)',n_w,n_viewing_polar_a,n_w*n_viewing_polar_a,n_w*n_viewing_polar_a*8/1e9)); end;
for nviewing_polar_a=0:n_viewing_polar_a-1;
viewing_polar_a = viewing_polar_a_(1+nviewing_polar_a); ca = cos(viewing_polar_a); sa = sin(viewing_polar_a);
condense_k_c_2__(:,1+nviewing_polar_a) = -sa*cc_ ;
end;%for nviewing_polar_a=0:n_viewing_all-1;

%%%%%%%%;
if (verbose); disp(sprintf(' %% Now evaluate associated legendre polynomials at the varous k_c_2 values.')); end;
% Here legendre_evaluate_lmwS____{1+l_val}(1+l_val+m_val,1+nw,1+nviewing_polar_a) contains ;
% the associated legendre-function of degree l_val and order abs(m_val) (ranging from 0 to +l_val) ;
% evaluated at the k_c_2 value stored in condense_k_c_2__(1+nw,1+nviewing_polar_a). ;
% Note that this is associated with viewing_polar_a_(1+nviewing_polar_a). ;
% The legendre_normalization_{1+l_val}(1+abs(m_val)) contains ;
% The normalization coefficient for the spherical harmonics associated with l_val and m_val. ;
% Note that this is somewhat redundant (as it does not depend explicitly on the shell). ;
%%%%%%%%;
legendre_evaluate_lmwS____ = cell(l_max+1,1);
legendre_normalization_ = cell(l_max+1,1);
for l_val=0:l_max;
tmp_P__ = zeros(1+1*l_val,n_w,n_viewing_polar_a);
tmp_a1 = ((1+2*l_val)/(4*pi));
m_val_ = -l_val:+l_val;
tmp_a2_ = exp(lfactorial(l_val-abs(m_val_)) - lfactorial(l_val+abs(m_val_)));
tmp_a3_ = sqrt(tmp_a1*tmp_a2_);
tmp_P__ = legendre(l_val,condense_k_c_2__,'unnorm');
tmp_t = tic;
legendre_evaluate_lmwS____{1+l_val} = reshape(tmp_P__,[1+1*l_val,n_w,n_viewing_polar_a]);
tmp_t = toc(tmp_t);
if (verbose>1); disp(sprintf(' %% l_val %d/%d legendre_evaluate(%d,%d) %0.2fs',l_val,l_max,n_w,n_viewing_polar_a,tmp_t)); end;
legendre_normalization_{1+l_val} = tmp_a3_;
end;%for l_val=0:l_max;
%%%%%%%%;
legendre_evaluate_normalized_lmwp____ = zeros(l_max,n_m_max,n_w,n_viewing_polar_a);
if (verbose); disp(sprintf(' %% legendre_evaluate_normalized_lmwp____: (%d,%d,%d,%d)=%d (%0.2f GB)',(1+l_max),n_m_max,n_w,n_viewing_polar_a,(1+l_max)*n_m_max*n_w*n_viewing_polar_a,(1+l_max)*n_m_max*n_w*n_viewing_polar_a*8/1e9)); end;
for l_val=0:l_max;
index_m_out_ = l_max + [-l_val:+l_val];
index_m_0in_ = [l_val:-1:1,0:l_val];
legendre_evaluate_normalized_lmwp____(1+l_val,1+index_m_out_,:,:) = bsxfun(@times,legendre_evaluate_lmwS____{1+l_val}(1+index_m_0in_,:,:),reshape(legendre_normalization_{1+l_val},[1+2*l_val,1,1]));
end;%for l_val=0:l_max;
%%%%%%%%;
legendre_evaluate_normalized_lwpm___ = reshape(permute(legendre_evaluate_normalized_lmwp____,[1,3,4,2]),[1+l_max,n_w*n_viewing_polar_a,n_m_max]);
if (verbose); disp(sprintf(' %% legendre_evaluate_normalized_lwpm____: (%d,%d,%d,%d)=%d (%0.2f GB)',(1+l_max),n_w,n_viewing_polar_a,n_m_max,(1+l_max)*n_w*n_viewing_polar_a*n_m_max,(1+l_max)*n_w*n_viewing_polar_a*n_m_max*8/1e9)); end;

%%%%%%%%;
% unroll a_k_Y_ya__. ;
%%%%%%%%;
a_k_Y_lma___ = zeros(1+l_max,n_m_max,n_a);
if (verbose); disp(sprintf(' %% a_k_Y_lma___: (%d,%d,%d)=%d (%0.2f GB)',(1+l_max),n_m_max,n_a,(1+l_max)*n_m_max*n_a,(1+l_max)*n_m_max*n_a*16/1e9)); end;
for l_val=0:l_max;
index_m_0in_ = (1+l_val-1)^2+l_val+[-l_val:+l_val];
index_m_out_ = l_max + [-l_val:+l_val];
a_k_Y_lma___(1+l_val,1+index_m_out_,:) = a_k_Y_ya__(1+index_m_0in_,:);
end;%for l_val=0:l_max;
a_k_Y_lam___ = permute(a_k_Y_lma___,[1,3,2]);
if (verbose); disp(sprintf(' %% a_k_Y_lam___: (%d,%d,%d)=%d (%0.2f GB)',(1+l_max),n_a,n_m_max,(1+l_max)*n_a*n_m_max,(1+l_max)*n_a*n_m_max*16/1e9)); end;

%%%%%%%%;
if (verbose); disp(sprintf(' %% Now accumulate the legendre_evaluates over l_val, for each m_val.')); end;
% We account for the normalization coefficients here, ;
% so that later we can apply the complex exponential to produce the spherical-harmonics. ;
% More specifically, for each na: ;
% spherical_harmonic_unphased_(1+l_max+m_val,1+nw,1+nviewing_polar_a) ;
% contains the sum: ;
% \sum_{l_val=0}^{l_max} ... ;
%             legendre_normalization_{1+l_val}(1+l_val+m_val) ... ;
%           * legendre_evaluate_lmwS____{1+l_val}(1+l_val+m_val,1+nw,1+nviewing_polar_a) ... ;
%           * a_k_Y_(1+(1+l_val-1)^2+l_val+m_val). ;
% Note that this is 'unphased', in the sense that the full evaluate requires: ;
% spherical_harmonic_evaluate__(1+nw,1+nviewing_all) = ... ;
% \sum_{m_val=-l_max}^{+l_max} ... ;
%             spherical_harmonic_unphased_(1+l_max+m_val,1+nw,1+nviewing_polar_a) ... ;
%           * exp(+i*m_val*template_azimu_b__(1+nw,1+nviewing_all)). ;
% Note that this final exponential can be calculated as: ;
% expi_template_azimu_b__(1+nw,1+nviewing_all)^m_val. ;
%%%%%%%%;
%{
spherical_harmonic_unphased_mwpa____ = zeros(1+2*l_max,n_w,n_viewing_polar_a,n_a);
for m_val=-l_max:l_max;
for l_val=abs(m_val):l_max;
spherical_harmonic_unphased_mwpa____(1+l_max+m_val,:,:,:) = spherical_harmonic_unphased_mwpa____(1+l_max+m_val,:,:,:) + bsxfun(@times,reshape(legendre_evaluate_lmwS____{1+l_val}(1+abs(m_val),:,:),[1,n_w,n_viewing_polar_a,1]),reshape(legendre_normalization_{1+l_val}(1+l_val+m_val)*a_k_Y_ya__(1+(1+l_val-1)^2+l_val+m_val,:),[1,1,1,n_a]));
end;%for l_val=abs(m_val):l_max;
end;%for m_val=-l_max:l_max;
%}
%%%%%%%%;

n_a_per_a_batch = max(1,ceil( 0.5e9 / (n_w*n_viewing_polar_a*(1+2*l_max)*8) )); %<-- 0.5GB limit. ;
n_a_batch = ceil(n_a/n_a_per_a_batch);
if (verbose); disp(sprintf(' %% n_a_per_a_batch %d, n_a_batch %d',n_a_per_a_batch,n_a_batch)); end;
for na_batch=0:n_a_batch-1;
index_a_ = n_a_per_a_batch*na_batch + [0:n_a_per_a_batch-1];
index_a_ = index_a_(find(index_a_<n_a));
tmp_n_a = numel(index_a_);
if (tmp_n_a>0);
tmp_a_k_Y_lam___ = a_k_Y_lam___(:,1+index_a_,:);
  
spherical_harmonic_unphased_wpam____ = zeros(n_w,n_viewing_polar_a,tmp_n_a,1+2*l_max);
if (verbose);
disp(sprintf(' %% spherical_harmonic_unphased_wpam____: (%d,%d,%d,%d)=%d (%0.2f GB)' ...
,n_w,n_viewing_polar_a,tmp_n_a,1+2*l_max ...
,n_w*n_viewing_polar_a*tmp_n_a*(1+2*l_max) ...
,n_w*n_viewing_polar_a*tmp_n_a*(1+2*l_max)*8/1e9 ...
)); end;
for nm=0:n_m_max-1;
spherical_harmonic_unphased_wpam____(:,:,:,1+nm) = reshape(transpose(legendre_evaluate_normalized_lwpm___(:,:,1+nm))*tmp_a_k_Y_lam___(:,:,1+nm),[n_w,n_viewing_polar_a,tmp_n_a]);
end;%for nm=0:n_m_max-1;
spherical_harmonic_unphased_mawp____ = permute(spherical_harmonic_unphased_wpam____,[4,3,1,2]);
if (verbose); 
disp(sprintf(' %% spherical_harmonic_unphased_mawp____: (%d,%d,%d,%d)=%d (%0.2f GB)' ...
,(1+2*l_max),tmp_n_a,n_w,n_viewing_polar_a ...
,(1+2*l_max)*tmp_n_a*n_w*n_viewing_polar_a ...
,(1+2*l_max)*tmp_n_a*n_w*n_viewing_polar_a*8/1e9 ...
)); end;

%%%%%%%%;
if (verbose); disp(sprintf(' %% now perform the final sum over m_val.')); end;
%%%%%%%%;
%{
spherical_harmonic_evaluate_waS___ = zeros(n_w,n_a,n_viewing_all);
nviewing_all=0;
for nviewing_polar_a=0:n_viewing_polar_a-1;
n_viewing_azimu_b = n_viewing_azimu_b_(1+nviewing_polar_a);
for nviewing_azimu_b=0:n_viewing_azimu_b-1;
tmp_expi_sub = expi_template_azimu_b__(:,1+nviewing_all);
tmp_expi_pre = tmp_expi_sub.^(-l_max);
tmp_expi_pos = tmp_expi_pre;
tmp_sum_wa__ = zeros(n_w,n_a);
for m_val=-l_max:+l_max;
tmp_sum_wa__ = tmp_sum_wa__ + bsxfun(@times,reshape(spherical_harmonic_unphased_mwpa____(1+l_max+m_val,:,1+nviewing_polar_a,:),[n_w,n_a]),tmp_expi_pos);
tmp_expi_pos = tmp_expi_pos.*tmp_expi_sub;
end;%for m_val=-l_max:+l_max;
spherical_harmonic_evaluate_waS___(:,:,1+nviewing_all) = tmp_sum_wa__;
nviewing_all=nviewing_all+1;
end;%for nviewing_azimu_b=0:n_viewing_azimu_b-1;
end;%for nviewing_polar_a=0:n_viewing_polar_a-1;
%}

spherical_harmonic_evaluate_waS___ = zeros(n_w,tmp_n_a,n_viewing_all);
if (verbose);
disp(sprintf(' %% spherical_harmonic_evaluate_waS___: (%d,%d,%d)=%d (%0.2f GB)' ...
,n_w,tmp_n_a,n_viewing_all ...
,n_w*tmp_n_a*n_viewing_all ...
,n_w*tmp_n_a*n_viewing_all*16/1e9 ...
)); end;
nviewing_all=0;
for nviewing_polar_a=0:n_viewing_polar_a-1;
n_viewing_azimu_b = n_viewing_azimu_b_(1+nviewing_polar_a);
for nviewing_azimu_b=0:n_viewing_azimu_b-1;
tmp_expi_sub = expi_template_azimu_b__(:,1+nviewing_all);
tmp_expi_pre = tmp_expi_sub.^(-l_max);
tmp_expi_pos = tmp_expi_pre;
%%%%;
tmp_expi_wm__ = zeros(n_w,n_m_max); tmp_expi_wm__(:,1) = tmp_expi_pos;
for nm=1:n_m_max-1;
tmp_expi_wm__(:,1+nm) = tmp_expi_wm__(:,1+nm-1).*tmp_expi_sub;
end;%for nm=1:n_m_max-1;
%%%%;
tmp_sum_wa__ = zeros(n_w,tmp_n_a);
for nw=0:n_w-1;
tmp_sum_wa__(1+nw,:) = tmp_expi_wm__(1+nw,:)*spherical_harmonic_unphased_mawp____(:,:,1+nw,1+nviewing_polar_a);
end;%for nw=0:n_w-1;
spherical_harmonic_evaluate_waS___(:,:,1+nviewing_all) = tmp_sum_wa__;
%%%%;
%tmp_expi_mw__ = transpose(tmp_expi_wm__);
%spherical_harmonic_evaluate_waS___(:,:,1+nviewing_all) = transpose(squeeze(sum(bsxfun(@times,reshape(tmp_expi_mw__,[n_m_max,1,n_w]),spherical_harmonic_unphased_mawp____(:,:,:,1+nviewing_polar_a)),1)));
%%%%;
nviewing_all=nviewing_all+1;
end;%for nviewing_azimu_b=0:n_viewing_azimu_b-1;
end;%for nviewing_polar_a=0:n_viewing_polar_a-1;

template_waS___(:,1+index_a_,:) = spherical_harmonic_evaluate_waS___;

end;%if (tmp_n_a>0);
end;%for na_batch=0:n_a_batch-1;
