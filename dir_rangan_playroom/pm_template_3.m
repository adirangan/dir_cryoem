function ...
[ ...
 template_wkS___ ...
,n_w ...
,n_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,sqrt_2lp1_ ...
,sqrt_2mp1_ ...
,sqrt_rat0_m_ ...
,sqrt_rat3_lm__ ...
,sqrt_rat4_lm__ ...
] = ...
pm_template_3( ...
 flag_verbose ...
,l_max ...
,n_k ...
,a_k_Y_yk__ ...
,viewing_euler_k_eq_d ...
,template_inplane_k_eq_d ...
,n_w_input ...
,n_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,sqrt_2lp1_ ...
,sqrt_2mp1_ ...
,sqrt_rat0_m_ ...
,sqrt_rat3_lm__ ...
,sqrt_rat4_lm__ ...
);
% uses spherical-harmonic-expansions a_k_Y_yk__ to evaluate templates on a collection of points on spherical shells. ;
% each spherical-shell has the same resolution, determined by viewing_euler_k_eq_d and template_inplane_k_eq_d and/or n_w_max. ;
% ;
% inputs: ;
% ;
% flag_verbose = integer verbosity_level. ;
% l_max = spherical harmonic order on each shell. ;
%         corresponds to n_y = (l_max+1)^2 coefficients. ;
% n_k = integer number of shells. ;
% a_k_Y_yk__ = complex array of size (n_y,n_k). ;
%              coefficients are ordered in a row, with m varying quickly, l varying slowly. ;
% viewing_euler_k_eq_d = real equatorial-distance used for sampling viewing angles and templates. ;
% template_inplane_k_eq_d = real equatorial-distance used for sampling inplane-shifts along each template. ;
% n_w_input = integer. used if template_inplane_k_eq_d <=0; desired n_w for templates (actual n_w:=max(6,n_w_input)). ;
% n_S = integer. number of viewing angles (i.e., number of templates) .;
% viewing_azimu_b_S_ = real array of size (n_S,1). ;
%                        azimu_b values for each template. ;
% viewing_polar_a_S_ = real array of size (n_S,1). ;
%                        polar_a values for each template. ;
% viewing_weight_S_ = real array of size (n_S,1). ;
%                       integration weight (on shell of radius 1) for each template. ;
% n_viewing_polar_a = integer. number of distinct polar_a across the viewing angles. ;
% viewing_polar_a_ = real array of size (n_viewing_polar_a,1). ;
%                    polar_a values for each viewing_polar_a_. ;
% n_viewing_azimu_b_ = integer array of size (n_viewing_polar_a,1). ;
%                      number of azimu_b values for each polar_a. ;
%                      These azimu_b values are assumed to be equispaced on [0,2*pi). ;
% ;
% outputs: ;
% ;
% template_wkS___ = complex array of templates for each viewing angle. ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

str_thisfunction = 'pm_template_3';

if nargin<1;
flag_verbose = 1;
n_k_p_r = 49;%n_k_p_r = 7;
k_p_r_max = 1;
k_p_r_ = ones(n_k_p_r,1);
weight_k_p_r_ = ones(n_k_p_r,1);
l_max = 80;%l_max = 13;
n_y = (l_max+1)^2;
l_max_ = l_max*ones(n_k_p_r,1);
n_y_sum = sum((1+l_max_).^2);
a_k_Y_ = (mod(transpose([0:n_y_sum-1]),89)-44)/89 + i*(mod(transpose([0:n_y_sum-1]),97)-48)/97;
viewing_euler_k_eq_d = 1/(2*pi);
template_inplane_k_eq_d = -1;
n_w_max = max(6,98);
%%%%%%%%;
tmp_t = tic();
[ ...
 template2_wkS___ ...
,n_w_ ...
,n_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
] = ...
pm_template_2( ...
 flag_verbose ...
,l_max ...
,n_k_p_r ...
,reshape(a_k_Y_,[n_y,n_k_p_r]) ...
,viewing_euler_k_eq_d ...
,template_inplane_k_eq_d ...
,n_w_max ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% pm_template_3: %0.2fs',tmp_t));
disp(sprintf(' %% n_S %d',n_S));
template2_wkS__ = reshape(template2_wkS___,[n_w_max*n_k_p_r,n_S]);
%%%%%%%%;
tmp_t = tic();
[ ...
 template3_wkS___ ...
,n_w_ ...
,n_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
] = ...
pm_template_3( ...
 flag_verbose ...
,l_max ...
,n_k_p_r ...
,reshape(a_k_Y_,[n_y,n_k_p_r]) ...
,viewing_euler_k_eq_d ...
,template_inplane_k_eq_d ...
,n_w_max ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% pm_template_3: %0.2fs',tmp_t));
disp(sprintf(' %% n_S %d',n_S));
template3_wkS__ = reshape(template3_wkS___,[n_w_max*n_k_p_r,n_S]);
%%%%%%%%;
fnorm_disp(flag_verbose,'template2_wkS__',template2_wkS__,'template3_wkS__',template3_wkS__,' %%<-- should be zero');
disp(sprintf(' %% returning')); return;
end;% if nargin<7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

na=0;
if (nargin<1+na); flag_verbose=[]; end; na=na+1;
if (nargin<1+na); l_max=[]; end; na=na+1;
if (nargin<1+na); n_k=[]; end; na=na+1;
if (nargin<1+na); a_k_Y_yk__=[]; end; na=na+1;
if (nargin<1+na); viewing_euler_k_eq_d=[]; end; na=na+1;
if (nargin<1+na); template_inplane_k_eq_d=[]; end; na=na+1;
if (nargin<1+na); n_w_input=[]; end; na=na+1;
if (nargin<1+na); n_S=[]; end; na=na+1;
if (nargin<1+na); viewing_azimu_b_S_=[]; end; na=na+1;
if (nargin<1+na); viewing_polar_a_S_=[]; end; na=na+1;
if (nargin<1+na); viewing_weight_S_=[]; end; na=na+1;
if (nargin<1+na); n_viewing_polar_a=[]; end; na=na+1;
if (nargin<1+na); viewing_polar_a_=[]; end; na=na+1;
if (nargin<1+na); n_viewing_azimu_b_=[]; end; na=na+1;
if (nargin<1+na); sqrt_2lp1_ = []; end; na=na+1;
if (nargin<1+na); sqrt_2mp1_ = []; end; na=na+1;
if (nargin<1+na); sqrt_rat0_m_ = []; end; na=na+1;
if (nargin<1+na); sqrt_rat3_lm__ = []; end; na=na+1;
if (nargin<1+na); sqrt_rat4_lm__ = []; end; na=na+1;

if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

memory_limit_GB = 0.5; %<-- 0.5GB for n_k_per_k_batch. ;

n_y = (l_max+1)^2;
m_max_ = -l_max : +l_max;
n_m_max = length(m_max_); %<-- 2*l_max+1;
if (flag_verbose); disp(sprintf(' %% l_max %d n_y %d',l_max,n_y)); end;
if (isempty(a_k_Y_yk__)); a_k_Y_yk__ = zeros(n_y,1); end;

%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% First determine the viewing angles.')); end;
%%%%%%%%;
tmp_t=tic();
if isempty(n_S);
k_p_r = 1;
[ ...
 n_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,viewing_weight_S_ ...
,~ ...
,~ ...
,~ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
] = ...
sample_shell_5( ...
 k_p_r ...
,viewing_euler_k_eq_d ...
,'L' ...
) ; %<-- obtain viewing angles. ;
end;%if isempty(n_S);
%%%%;
if isempty(viewing_weight_S_); viewing_weight_S_ = ones(n_S,1); end;
%%%%;
if isempty(n_viewing_polar_a);
viewing_polar_a_ = unique(viewing_polar_a_S_); n_viewing_polar_a = numel(viewing_polar_a_);
n_viewing_azimu_b_ = zeros(n_viewing_polar_a,1);
for nviewing_polar_a=0:n_viewing_polar_a-1;
viewing_polar_a = viewing_polar_a_(1+nviewing_polar_a);
n_viewing_azimu_b_(1+nviewing_polar_a) = numel(efind(abs(viewing_polar_a_S_-viewing_polar_a)<1e-12));
end;%for nviewing_polar_a=0:n_viewing_polar_a-1;
end;%if isempty(n_viewing_polar_a);
%%%%;
n_viewing_azimu_b_sum = sum(n_viewing_azimu_b_);
n_viewing_azimu_b_csum_ = cumsum([0;n_viewing_azimu_b_]);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% determine viewing-angles: %0.6fs',tmp_t)); end;
if (flag_verbose>0); disp(sprintf(' %% n_S %d n_viewing_polar_a %d n_viewing_azimu_b_sum %d',n_S,n_viewing_polar_a,n_viewing_azimu_b_sum)); end;
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% Now determine the points along each equatorial plane (i.e., the points for each template).')); end;
%%%%%%%%;
tmp_t=tic();
n_w = 0;
if (template_inplane_k_eq_d>0);
k_p_r = 1;
n_equator = 3+round(2*pi*k_p_r/template_inplane_k_eq_d);
n_polar_a = 3+round(n_equator/2);
n_w = 2*n_polar_a;
end;%if (template_inplane_k_eq_d>0);
if (template_inplane_k_eq_d<=0);
n_w = max(6,n_w_input);
end;%if (template_inplane_k_eq_d<=0);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% determine n_w: %0.6fs',tmp_t)); end;
if (flag_verbose>0); disp(sprintf(' %% n_w %d',n_w)); end;
%%%%%%%%;
% Set up inner gamma_z for the templates. ;
%%%%%%%%;
gamma_z_ = zeros(n_w,1); gamma_z_ = transpose(linspace(0,2*pi,n_w+1)); gamma_z_ = gamma_z_(1:n_w);
cc_ = cos(gamma_z_); sc_ = sin(gamma_z_);
%%%%%%%%;
% initialize template_wkS___. ;
%%%%%%%%;
template_wkS___ = zeros(n_w,n_k,n_S);

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

if (flag_verbose>0); disp(sprintf(' %% template: (%d,%d,%d)=%d (%0.2f GB)',n_w,n_S,n_k,n_w*n_S*n_k,n_w*n_S*n_k*16/1e9)); end;

%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% Now construct array of k_c_?_ values for the templates.')); end;
%%%%%%%%;
tmp_t=tic();
template_k_p_r = 1;
%%%%%%%%;
template_k_c_0__ = zeros(n_w,n_S);
template_k_c_1__ = zeros(n_w,n_S);
%template_k_c_2__ = zeros(n_w,n_S);
if (flag_verbose>0); disp(sprintf(' %% template_k_c___: (%d,%d)=%d (%0.2f GB)',n_w,n_S,n_w*n_S,n_w*n_S*8/1e9)); end;
ca_ = cos(viewing_polar_a_S_); sa_ = sin(viewing_polar_a_S_);
cb_ = cos(viewing_azimu_b_S_); sb_ = sin(viewing_azimu_b_S_);
ca__ = repmat(reshape(ca_,[1,n_S]),[n_w,1]); sa__ = repmat(reshape(sa_,[1,n_S]),[n_w,1]);
cb__ = repmat(reshape(cb_,[1,n_S]),[n_w,1]); sb__ = repmat(reshape(sb_,[1,n_S]),[n_w,1]);
cc__ = repmat(reshape(cc_,[n_w,1]),[1,n_S]); sc__ = repmat(reshape(sc_,[n_w,1]),[1,n_S]);
template_k_c_0__ = (+cb__.*ca__.*cc__ - sb__.*sc__).*template_k_p_r;
template_k_c_1__ = (+sb__.*ca__.*cc__ + cb__.*sc__).*template_k_p_r;
%template_k_c_2__ = (-sa__.*cc__                   ).*template_k_p_r;
template_azimu_b__ = atan2(template_k_c_1__,template_k_c_0__);
expi_template_azimu_b__ = exp(i*template_azimu_b__);
clear template_azimu_b__;
clear template_k_c_0__ template_k_c_1__ template_k_c_2__ ;
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% expi_template_azimu_b__: %0.6fs',tmp_t)); end;
%%%%%%%%;
% We also use a condensed array, called condense_k_c_2__, which only depends on the polar_a, and not on azimu_b. ;
%%%%%%%%;
tmp_t=tic();
condense_k_c_2__ = zeros(n_w,n_viewing_polar_a);
if (flag_verbose>0); disp(sprintf(' %% condense_k_c_2__: (%d,%d)=%d (%0.2f GB)',n_w,n_viewing_polar_a,n_w*n_viewing_polar_a,n_w*n_viewing_polar_a*8/1e9)); end;
condense_k_c_2__ = -bsxfun(@times,reshape(sin(viewing_polar_a_),[1,n_viewing_polar_a]),reshape(cc_,[n_w,1]));
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% condense_k_c_2__: %0.6fs',tmp_t)); end;
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% Now evaluate associated legendre polynomials at the various k_c_2 values.')); end;
%%%%%%%%;
% we use ylgndr_2 (i.e., yrecursion.f). ;
%%%%%%%%;
tmp_t=tic();
[ ...
 ~ ...
,y_jlm___ ...
,~ ...
,~ ...
,sqrt_2lp1_ ...
,sqrt_2mp1_ ...
,sqrt_rat0_m_ ...
,sqrt_rat3_lm__ ...
,sqrt_rat4_lm__ ...
] = ...
ylgndr_2( ...
 struct('flag_verbose',0,'flag_d',0,'flag_dd',0) ...
,l_max ...
,condense_k_c_2__(:) ...
,sqrt_2lp1_ ...
,sqrt_2mp1_ ...
,sqrt_rat0_m_ ...
,sqrt_rat3_lm__ ...
,sqrt_rat4_lm__ ...
);
y_lmwa____ = reshape(permute(y_jlm___,1+[1,2,0]),[1+l_max,1+l_max,n_w,n_viewing_polar_a])/sqrt(4*pi);
legendre_use_evaluate_normalized_lmwa____ = zeros(1+l_max,n_m_max,n_w,n_viewing_polar_a);
legendre_use_evaluate_normalized_lmwa____(:,1+l_max+[0:+1:+l_max],:,:) = y_lmwa____;
legendre_use_evaluate_normalized_lmwa____(:,1+l_max+[0:-1:-l_max],:,:) = y_lmwa____;
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% legendre_use_evaluate_normalize_lmwa____: %0.6fs',tmp_t)); end;
%%%%%%%%;

%%%%%%%%;
tmp_t=tic();
legendre_evaluate_normalized_lwam___ = reshape(permute(legendre_use_evaluate_normalized_lmwa____,1+[0,2,3,1]),[1+l_max,n_w*n_viewing_polar_a,n_m_max]);
if (flag_verbose>0); disp(sprintf(' %% legendre_evaluate_normalized_lwam____: (%d,%d,%d,%d)=%d (%0.2f GB)',(1+l_max),n_w,n_viewing_polar_a,n_m_max,(1+l_max)*n_w*n_viewing_polar_a*n_m_max,(1+l_max)*n_w*n_viewing_polar_a*n_m_max*8/1e9)); end;
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% legendre_use_evaluate_normalize_lwam____: %0.6fs',tmp_t)); end;
%%%%%%%%;

%%%%%%%%;
% unroll a_k_Y_yk__. ;
%%%%%%%%;
tmp_t=tic();
a_k_Y_lmk___ = zeros(1+l_max,n_m_max,n_k);
if (flag_verbose>0); disp(sprintf(' %% a_k_Y_lmk___: (%d,%d,%d)=%d (%0.2f GB)',(1+l_max),n_m_max,n_k,(1+l_max)*n_m_max*n_k,(1+l_max)*n_m_max*n_k*16/1e9)); end;
for l_val=0:l_max;
index_m_0in_ = (1+l_val-1)^2+l_val+[-l_val:+l_val];
index_m_out_ = l_max + [-l_val:+l_val];
a_k_Y_lmk___(1+l_val,1+index_m_out_,:) = a_k_Y_yk__(1+index_m_0in_,:);
end;%for l_val=0:l_max;
a_k_Y_lkm___ = permute(a_k_Y_lmk___,1+[0,2,1]);
if (flag_verbose>0); disp(sprintf(' %% a_k_Y_lkm___: (%d,%d,%d)=%d (%0.2f GB)',(1+l_max),n_k,n_m_max,(1+l_max)*n_k*n_m_max,(1+l_max)*n_k*n_m_max*16/1e9)); end;
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% a_k_Y_lkm___: %0.6fs',tmp_t)); end;

%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% Now accumulate the legendre_evaluates over l_val, for each m_val.')); end;
% We account for the normalization coefficients here, ;
% so that later we can apply the complex exponential to produce the spherical-harmonics. ;
% More specifically, for each nk: ;
% spherical_harmonic_unphased_(1+l_max+m_val,1+nw,1+nviewing_polar_a) ;
% contains the sum: ;
% \sum_{l_val=0}^{l_max} ... ;
%             legendre_normalization_{1+l_val}(1+l_val+m_val) ... ;
%           * legendre_evaluate_lmwS____{1+l_val}(1+l_val+m_val,1+nw,1+nviewing_polar_a) ... ;
%           * a_k_Y_(1+(1+l_val-1)^2+l_val+m_val). ;
% Note that this is 'unphased', in the sense that the full evaluate requires: ;
% spherical_harmonic_evaluate__(1+nw,1+nS) = ... ;
% \sum_{m_val=-l_max}^{+l_max} ... ;
%             spherical_harmonic_unphased_(1+l_max+m_val,1+nw,1+nviewing_polar_a) ... ;
%           * exp(+i*m_val*template_azimu_b__(1+nw,1+nS)). ;
% Note that this final exponential can be calculated as: ;
% expi_template_azimu_b__(1+nw,1+nS)^m_val. ;
%%%%%%%%;
%{
spherical_harmonic_unphased_mwak____ = zeros(1+2*l_max,n_w,n_viewing_polar_a,n_k);
for m_val=-l_max:l_max;
for l_val=abs(m_val):l_max;
spherical_harmonic_unphased_mwak____(1+l_max+m_val,:,:,:) = spherical_harmonic_unphased_mwak____(1+l_max+m_val,:,:,:) + bsxfun(@times,reshape(legendre_evaluate_lmwS____{1+l_val}(1+abs(m_val),:,:),[1,n_w,n_viewing_polar_a,1]),reshape(legendre_normalization_{1+l_val}(1+l_val+m_val)*a_k_Y_yk__(1+(1+l_val-1)^2+l_val+m_val,:),[1,1,1,n_k]));
end;%for l_val=abs(m_val):l_max;
end;%for m_val=-l_max:l_max;
%}
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
n_k_per_k_batch = min(n_k,max(1,ceil( memory_limit_GB * 1e9 / max(1,n_w*n_viewing_polar_a*(1+2*l_max)*8) ))); %<-- memory_limit_GB. ;
if (flag_verbose>2); n_k_per_k_batch=min(8,n_k_per_k_batch); end;
n_k_batch = ceil(n_k/max(1,n_k_per_k_batch));
if (flag_verbose>0); disp(sprintf(' %% n_k_per_k_batch %d, n_k_batch %d',n_k_per_k_batch,n_k_batch)); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for nk_batch=0:n_k_batch-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
index_k_ = n_k_per_k_batch*nk_batch + [0:n_k_per_k_batch-1];
index_k_ = index_k_(find(index_k_<n_k));
n_k_sub = numel(index_k_);
%%%%%%%%%%%%%%%%;
if (n_k_sub>0);
%%%%%%%%%%%%%%%%;

tmp_a_k_Y_lkm___ = a_k_Y_lkm___(:,1+index_k_,:);
  
tmp_t=tic();
spherical_harmonic_unphased_wakm____ = zeros(n_w,n_viewing_polar_a,n_k_sub,1+2*l_max);
if (flag_verbose>0);
disp(sprintf(' %% spherical_harmonic_unphased_wakm____: (%d,%d,%d,%d)=%d (%0.2f GB)' ...
,n_w,n_viewing_polar_a,n_k_sub,1+2*l_max ...
,n_w*n_viewing_polar_a*n_k_sub*(1+2*l_max) ...
,n_w*n_viewing_polar_a*n_k_sub*(1+2*l_max)*8/1e9 ...
)); end;
spherical_harmonic_unphased_wakm____ = reshape(pagemtimes(legendre_evaluate_normalized_lwam___,'transpose',tmp_a_k_Y_lkm___,'none'),[n_w,n_viewing_polar_a,n_k_sub,1+2*l_max]);
spherical_harmonic_unphased_mkwa____ = permute(spherical_harmonic_unphased_wakm____,1+[3,2,0,1]);
if (flag_verbose>0); 
disp(sprintf(' %% spherical_harmonic_unphased_mkwa____: (%d,%d,%d,%d)=%d (%0.2f GB)' ...
,(1+2*l_max),n_k_sub,n_w,n_viewing_polar_a ...
,(1+2*l_max)*n_k_sub*n_w*n_viewing_polar_a ...
,(1+2*l_max)*n_k_sub*n_w*n_viewing_polar_a*8/1e9 ...
)); end;
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% spherical_harmonic_unphased_mkwa____: %0.6fs',tmp_t)); end;

%%%%%%%%;
% Here we rewrite the final sum from pm_template_2 to use pagemtimes. ;
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% now perform the final sum over m_val.')); end;
%%%%%%%%;
tmp_t=tic();
spherical_harmonic_evaluate_wkS___ = zeros(n_w,n_k_sub,n_S);
if (flag_verbose>0);
disp(sprintf(' %% spherical_harmonic_evaluate_wkS___: (%d,%d,%d)=%d (%0.2f GB)' ...
,n_w,n_k_sub,n_S ...
,n_w*n_k_sub*n_S ...
,n_w*n_k_sub*n_S*16/1e9 ...
)); end;
n_1 = 1;
nS=0;
for nviewing_polar_a=0:n_viewing_polar_a-1;
tmp_shu_mkw1____ = reshape(spherical_harmonic_unphased_mkwa____(:,:,:,1+nviewing_polar_a),[n_m_max,n_k_sub,n_w,n_1]);
n_viewing_azimu_b = n_viewing_azimu_b_(1+nviewing_polar_a);
n_S_sub = n_viewing_azimu_b;
index_nS_sub_ = nS + [0:n_S_sub-1];
tmp_expi_sub_wSm___ = bsxfun(@power,expi_template_azimu_b__(:,1+index_nS_sub_),reshape([-l_max:+l_max],[1,1,n_m_max]));
tmp_expi_sub_1mwS____ = reshape(permute(tmp_expi_sub_wSm___,1+[2,0,1]),[n_1,n_m_max,n_w,n_S_sub]);
spherical_harmonic_evaluate_wkS___(:,:,1+index_nS_sub_) = permute(reshape(pagemtimes(tmp_expi_sub_1mwS____,tmp_shu_mkw1____),[n_k_sub,n_w,n_S_sub]),1+[1,0,2]);
nS=nS+n_S_sub;
end;%for nviewing_polar_a=0:n_viewing_polar_a-1;
assert(nS==n_S);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% spherical_harmonic_evaluate_wkS___: %0.6fs',tmp_t)); end;

template_wkS___(:,1+index_k_,:) = spherical_harmonic_evaluate_wkS___;

%%%%%%%%%%%%%%%%;
end;%if (n_k_sub>0);
%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for nk_batch=0:n_k_batch-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% finished batches. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
