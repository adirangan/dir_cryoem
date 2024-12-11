function ...
[ ...
 template_wkS__ ...
,n_w ...
,n_viewing_all ...
,viewing_azimu_b_all_ ...
,viewing_polar_a_all_ ...
,viewing_weight_all_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,viewing_gamma_z ...
,dtemplateda_wkS__ ...
,dtemplatedb_wkS__ ...
,dtemplatedc_wkS__ ...
,ddtemplatedaa_wkS__ ...
,ddtemplatedab_wkS__ ...
,ddtemplatedac_wkS__ ...
,ddtemplatedbb_wkS__ ...
,ddtemplatedbc_wkS__ ...
,ddtemplatedcc_wkS__ ...
] = ...
fft_template_4( ...
 verbose ...
,n_x_u ...
,x_p_r_max ...
,a_x_u_ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,n_w_0in ...
,n_viewing_all ...
,viewing_azimu_b_all_ ...
,viewing_polar_a_all_ ...
,viewing_weight_all_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,viewing_gamma_z ...
);
% applies fft to a_x_u_ to evaluate templates on a collection of points on spherical shells. ;
% each spherical-shell has the same resolution, determined by viewing_k_eq_d and template_k_eq_d and/or n_w_max. ;
% ;
% inputs: ;
% ;
% verbose = integer verbosity_level. ;
% n_x_u = integer number of voxels on a side. ;
% x_p_r_max = real half-side-length of box. ;
%             We assume that the voxels extend from x_p_r_max*[-1,+1]. ;
% a_x_u_ = complex array of size (n_x_u,n_x_u,n_x_u). ;
% viewing_k_eq_d = real equatorial-distance used for sampling viewing angles and templates. ;
% template_k_eq_d = real equatorial-distance used for sampling inplane-shifts along each template. ;
% n_w_0in = integer. used if template_k_eq_d <=0; desired n_w for templates. ;
% n_viewing_all = integer. number of viewing angles (i.e., number of templates) .;
% viewing_azimu_b_all_ = real array of size (n_viewing_all,1). ;
%                        azimu_b values for each template. ;
% viewing_polar_a_all_ = real array of size (n_viewing_all,1). ;
%                        polar_a values for each template. ;
% viewing_weight_all_ = real array of size (n_viewing_all,1). ;
%                       integration weight (on shell of radius 1) for each template. ;
% n_viewing_polar_a = integer. number of distinct polar_a across the viewing angles. ;
% viewing_polar_a_ = real array of size (n_viewing_polar_a,1). ;
%                    polar_a values for each viewing_polar_a_. ;
% n_viewing_azimu_b_ = integer array of size (n_viewing_polar_a,1). ;
%                      number of azimu_b values for each polar_a. ;
%                      These azimu_b values are assumed to be equispaced on [0,2*pi). ;
% viewing_gamma_z = real gamma_z value for each template. (typically 0.0). ;
% ;
% outputs: ;
% ;
% template_wkS__ = complex array of templates. ;
%                  template_wkS__(1+nw+nk_p_r*n_w_max,1+nS) ;
%                  stores template value for angle-index nw, radial-index nk_p_r, ;
%                  and viewing_azimu_b = viewing_azimu_b_all_(1+nS). ;
%                  and viewing_polar_a = viewing_polar_a_all_(1+nS). ;
% dtemplatedx_wkS__ = complex array analogous to template_wkS__. ;
%                     stores first-derivative of template with respect to: ;
%                     x==a: polar_a ; %<-- note that the first-derivative with respect to polar_a has a different sign than wignerd_c produces. ;
%                     x==b: azimu_b ;
%                     x==c: gamma_z ;
% ddtemplatedxy_wkS__ = complex array analogous to template_wkS__. ;
%                       stores second-derivative of template with respect to x and y (see above). ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

str_thisfunction = 'fft_template_4';

if nargin<1;
verbose = 2; nf=1;
if (verbose); disp(sprintf(' %% testing %s',str_thisfunction)); end;
if (verbose); disp(sprintf(' %% calling test_CryoLike_fft_template_20241210')); end;
test_slice_vs_volume_integral_5;
disp(sprintf(' %% returning')); return;
end;% if nargin<7;
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

na=0;
if (nargin<1+na); verbose=[]; end; na=na+1;
if (nargin<1+na); n_x_u=[]; end; na=na+1;
if (nargin<1+na); x_p_r_max=[]; end; na=na+1;
if (nargin<1+na); a_x_u_=[]; end; na=na+1;
if (nargin<1+na); viewing_k_eq_d=[]; end; na=na+1;
if (nargin<1+na); template_k_eq_d=[]; end; na=na+1;
if (nargin<1+na); n_w_0in=[]; end; na=na+1;
if (nargin<1+na); n_viewing_all=[]; end; na=na+1;
if (nargin<1+na); viewing_azimu_b_all_=[]; end; na=na+1;
if (nargin<1+na); viewing_polar_a_all_=[]; end; na=na+1;
if (nargin<1+na); viewing_weight_all_=[]; end; na=na+1;
if (nargin<1+na); n_viewing_polar_a=[]; end; na=na+1;
if (nargin<1+na); viewing_polar_a_=[]; end; na=na+1;
if (nargin<1+na); n_viewing_azimu_b_=[]; end; na=na+1;
if (nargin<1+na); viewing_gamma_z=[]; end; na=na+1;

flag_d0 = 1;
flag_d1 = (nargout>=11);
flag_d2 = (nargout>=14);

if isempty(x_p_r_max); x_p_r_max = 1.0; end;
n_xxx_u = n_x_u^3;
if (verbose); disp(sprintf(' %% n_x_u %d',n_x_u)); end;
if (isempty(a_x_u_)); a_x_u_ = zeros(n_xxx_u,1); end;
if (verbose>2); fprintf(1,' %% \t\n'); darray_printf_margin(real(a_x_u_),1,n_xxx_u,' %% a_x_u_real___: '); end;
if (verbose>2); fprintf(1,' %% \t\n'); darray_printf_margin(imag(a_x_u_),1,n_xxx_u,' %% a_x_u_imag___: '); end;

%%%%%%%%;
if (verbose); disp(sprintf(' %% First determine the viewing angles.')); end;
%%%%%%%%;
if isempty(n_viewing_all);
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
end;%if isempty(n_viewing_all);
%%%%;
if isempty(viewing_weight_all_); viewing_weight_all_ = ones(n_viewing_all,1); end;
%%%%;
if isempty(n_viewing_polar_a);
viewing_polar_a_ = unique(viewing_polar_a_all_); n_viewing_polar_a = numel(viewing_polar_a_);
n_viewing_azimu_b_ = zeros(n_viewing_polar_a,1);
for nviewing_polar_a=0:n_viewing_polar_a-1;
viewing_polar_a = viewing_polar_a_(1+nviewing_polar_a);
n_viewing_azimu_b_(1+nviewing_polar_a) = numel(efind(abs(viewing_polar_a_all_-viewing_polar_a)<1e-12));
end;%for nviewing_polar_a=0:n_viewing_polar_a-1;
end;%if isempty(n_viewing_polar_a);
%%%%;
n_viewing_azimu_b_sum = sum(n_viewing_azimu_b_);
n_viewing_azimu_b_csum_ = cumsum([0;n_viewing_azimu_b_]);
if (verbose); disp(sprintf(' %% n_viewing_all %d n_viewing_polar_a %d n_viewing_azimu_b_sum %d',n_viewing_all,n_viewing_polar_a,n_viewing_azimu_b_sum)); end;
if (verbose>2); fprintf(1,' %% \t\n'); darray_printf_margin(viewing_azimu_b_all_,1,n_viewing_all,' %% viewing_azimu_b_all_: '); end;
if (verbose>2); fprintf(1,' %% \t\n'); darray_printf_margin(viewing_polar_a_all_,1,n_viewing_all,' %% viewing_polar_a_all_: '); end;
if (verbose>2); fprintf(1,' %% \t\n'); darray_printf_margin(viewing_weight_all_,1,n_viewing_all,' %% viewing_weight_all_: '); end;
if (verbose>2); fprintf(1,' %% \t\n'); darray_printf_margin(viewing_polar_a_,1,n_viewing_polar_a,' %% viewing_polar_a_: '); end;
if (verbose>2); fprintf(1,' %% \t\n'); darray_printf_margin(n_viewing_azimu_b_,1,n_viewing_polar_a,' %% n_viewing_azimu_b_: '); end;
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
%n_w = max(6,n_w_0in);
n_w = n_w_0in; %<-- no minimum. ;
end;%if (template_k_eq_d<=0);
if (verbose); disp(sprintf(' %% n_w %d',n_w)); end;
n_w_max = n_w; n_w_sum = n_w_max*n_k_p_r; n_w_ = n_w_max*ones(n_k_p_r,1); n_w_csum_ = cumsum([0;n_w_]);
%%%%%%%%;

%%%%%%%%;
% Set up inner gamma_z for the templates. ;
%%%%%%%%;
gamma_z_ = zeros(n_w,1); gamma_z_ = transpose(linspace(0,2*pi,n_w+1)); gamma_z_ = gamma_z_(1:n_w);
cc_ = cos(gamma_z_); sc_ = sin(gamma_z_);
if (verbose>2); fprintf(1,' %% \t\n'); darray_printf_margin(gamma_z_,1,n_w,' %% gamma_z_: '); end;
if (verbose>2); fprintf(1,' %% \t\n'); darray_printf_margin(cc_,1,n_w,' %% cc_: '); end;
if (verbose>2); fprintf(1,' %% \t\n'); darray_printf_margin(sc_,1,n_w,' %% sc_: '); end;
%%%%%%%%;

if isempty(a_k_Y_lkm___);
%%%%%%%%;
% Now we unroll the a_k_Y_yk__. ;
%%%%%%%%;
tmp_t = tic();
a_k_Y_lmk___ = zeros(1+l_max_max,n_m_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_a_k_Y_lm_ = a_k_Y_yk__(:,1+nk_p_r);
tmp_a_k_Y_lm__ = zeros(1+l_max_max,n_m_max);
l_max = l_max_(1+nk_p_r);
na=0;
for l_val=0:l_max;
for m_val=-l_val:+l_val;
assert(na==l_val*(l_val+1)+m_val);
tmp_a_k_Y_lm__(1+l_val,1+l_max_max+m_val) = tmp_a_k_Y_lm_(1+na);
na=na+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max;
a_k_Y_lmk___(:,:,1+nk_p_r) = tmp_a_k_Y_lm__;
end;%for nk_p_r=0:n_k_p_r-1;
a_k_Y_lkm___ = permute(a_k_Y_lmk___,[1,3,2]);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% a_k_Y_lkm___: %0.2fs',tmp_t)); end;
%%%%%%%%;
end;%if isempty(a_k_Y_lkm___);

if flag_d0; template_wkS__ = zeros(n_w_sum,n_viewing_all); end;
if flag_d1; dtemplateda_wkS__ = zeros(n_w_sum,n_viewing_all); end;
if flag_d1; dtemplatedb_wkS__ = zeros(n_w_sum,n_viewing_all); end;
if flag_d1; dtemplatedc_wkS__ = zeros(n_w_sum,n_viewing_all); end;
if flag_d2; ddtemplatedaa_wkS__ = zeros(n_w_sum,n_viewing_all); end;
if flag_d2; ddtemplatedab_wkS__ = zeros(n_w_sum,n_viewing_all); end;
if flag_d2; ddtemplatedac_wkS__ = zeros(n_w_sum,n_viewing_all); end;
if flag_d2; ddtemplatedbb_wkS__ = zeros(n_w_sum,n_viewing_all); end;
if flag_d2; ddtemplatedbc_wkS__ = zeros(n_w_sum,n_viewing_all); end;
if flag_d2; ddtemplatedcc_wkS__ = zeros(n_w_sum,n_viewing_all); end;
flag_d0W_use = 1;
if flag_d0 & isempty(d0W_betazeta_mlma____); d0W_betazeta_mlma____ = zeros(n_m_max,1+l_max_max,n_m_max,n_viewing_polar_a); flag_d0W_use = 0; end;
flag_d1W_use = 1;
if flag_d1 & isempty(d1W_betazeta_mlma____); d1W_betazeta_mlma____ = zeros(n_m_max,1+l_max_max,n_m_max,n_viewing_polar_a); flag_d1W_use = 0; end;
flag_d2W_use = 1;
if flag_d2 & isempty(d2W_betazeta_mlma____); d2W_betazeta_mlma____ = zeros(n_m_max,1+l_max_max,n_m_max,n_viewing_polar_a); flag_d2W_use = 0; end;
%%%%%%%%;
for nviewing_polar_a=0:n_viewing_polar_a-1;
%%%%%%%%;
viewing_polar_a = viewing_polar_a_(1+nviewing_polar_a);
tmp_index_ = efind(abs(viewing_polar_a_all_-viewing_polar_a)<1e-12);
n_viewing_azimu_b = n_viewing_azimu_b_(1+nviewing_polar_a);
assert(numel(tmp_index_)==n_viewing_azimu_b);
viewing_azimu_b_ = viewing_azimu_b_all_(1+tmp_index_);
if (verbose); disp(sprintf(' %% n_viewing_azimu_b %d, size(tmp_index_) [%d,%d]',n_viewing_azimu_b,size(tmp_index_))); end;
if (verbose); fprintf(1,' %% \t\n'); darray_printf_margin(tmp_index_,1,n_viewing_azimu_b,' %% tmp_index_: '); end;
if isempty(viewing_gamma_z); viewing_gamma_z = 0.0; end;
d0W_betazeta_mlm___ = []; if flag_d0 & flag_d0W_use; d0W_betazeta_mlm___ = d0W_betazeta_mlma____(:,:,:,1+nviewing_polar_a); end;
d1W_betazeta_mlm___ = []; if flag_d1 & flag_d1W_use; d1W_betazeta_mlm___ = d1W_betazeta_mlma____(:,:,:,1+nviewing_polar_a); end;
d2W_betazeta_mlm___ = []; if flag_d2 & flag_d2W_use; d2W_betazeta_mlm___ = d2W_betazeta_mlma____(:,:,:,1+nviewing_polar_a); end;
%%%%;
if  flag_d0 & ~flag_d1 & ~flag_d2;
[ ...
 template_wkb__ ...
,V_lmm___ ...
,L_lm__ ...
,d0W_betazeta_mlm___ ...
] = ...
sph_template_single_polar_a_3( ...
 verbose ...
,l_max ...
,n_k_p_r ...
,a_k_Y_lkm___ ...
,n_w ...
,viewing_polar_a ...
,n_viewing_azimu_b ...
,viewing_azimu_b_ ...
,viewing_gamma_z ...
,V_lmm___ ...
,L_lm__ ...
,d0W_betazeta_mlm___ ...
);
if (verbose);
disp(sprintf(' %% flag_d0 %d(%d) flag_d1 %d(%d) flag_d2 %d(%d)',flag_d0,flag_d0W_use,flag_d1,flag_d1W_use,flag_d2,flag_d2W_use));
disp(sprintf(' %% size(template_wkb__): [%d,%d]',size(template_wkb__)));
disp(sprintf(' %% size(d0W_betazeta_mlm___): [%d,%d,%d]',size(d0W_betazeta_mlm___)));
end;%if (verbose);
end;%if  flag_d0 &  flag_d1 &  flag_d2;
%%%%;
if  flag_d0 &  flag_d1 & ~flag_d2;
[ ...
 template_wkb__ ...
,V_lmm___ ...
,L_lm__ ...
,d0W_betazeta_mlm___ ...
,dtemplateda_wkb__ ...
,dtemplatedb_wkb__ ...
,dtemplatedc_wkb__ ...
,d1W_betazeta_mlm___ ...
] = ...
sph_template_single_polar_a_3( ...
 verbose ...
,l_max ...
,n_k_p_r ...
,a_k_Y_lkm___ ...
,n_w ...
,viewing_polar_a ...
,n_viewing_azimu_b ...
,viewing_azimu_b_ ...
,viewing_gamma_z ...
,V_lmm___ ...
,L_lm__ ...
,d0W_betazeta_mlm___ ...
,d1W_betazeta_mlm___ ...
);
if (verbose);
disp(sprintf(' %% flag_d0 %d(%d) flag_d1 %d(%d) flag_d2 %d(%d)',flag_d0,flag_d0W_use,flag_d1,flag_d1W_use,flag_d2,flag_d2W_use));
disp(sprintf(' %% size(template_wkb__): [%d,%d]',size(template_wkb__)));
disp(sprintf(' %% size(d0W_betazeta_mlm___): [%d,%d,%d]',size(d0W_betazeta_mlm___)));
disp(sprintf(' %% size(dtemplateda_wkb__): [%d,%d]',size(dtemplateda_wkb__)));
disp(sprintf(' %% size(dtemplatedb_wkb__): [%d,%d]',size(dtemplatedb_wkb__)));
disp(sprintf(' %% size(dtemplatedc_wkb__): [%d,%d]',size(dtemplatedc_wkb__)));
disp(sprintf(' %% size(d1W_betazeta_mlm___): [%d,%d,%d]',size(d1W_betazeta_mlm___)));
end;%if (verbose);
end;%if  flag_d0 &  flag_d1 &  flag_d2;
%%%%;
if  flag_d0 &  flag_d1 &  flag_d2;
[ ...
 template_wkb__ ...
,V_lmm___ ...
,L_lm__ ...
,d0W_betazeta_mlm___ ...
,dtemplateda_wkb__ ...
,dtemplatedb_wkb__ ...
,dtemplatedc_wkb__ ...
,d1W_betazeta_mlm___ ...
,ddtemplatedaa_wkb__ ...
,ddtemplatedab_wkb__ ...
,ddtemplatedac_wkb__ ...
,ddtemplatedbb_wkb__ ...
,ddtemplatedbc_wkb__ ...
,ddtemplatedcc_wkb__ ...
,d2W_betazeta_mlm___ ...
] = ...
sph_template_single_polar_a_3( ...
 verbose ...
,l_max ...
,n_k_p_r ...
,a_k_Y_lkm___ ...
,n_w ...
,viewing_polar_a ...
,n_viewing_azimu_b ...
,viewing_azimu_b_ ...
,viewing_gamma_z ...
,V_lmm___ ...
,L_lm__ ...
,d0W_betazeta_mlm___ ...
,d1W_betazeta_mlm___ ...
,d2W_betazeta_mlm___ ...
);
if (verbose);
disp(sprintf(' %% flag_d0 %d(%d) flag_d1 %d(%d) flag_d2 %d(%d)',flag_d0,flag_d0W_use,flag_d1,flag_d1W_use,flag_d2,flag_d2W_use));
disp(sprintf(' %% size(template_wkb__): [%d,%d]',size(template_wkb__)));
disp(sprintf(' %% size(d0W_betazeta_mlm___): [%d,%d,%d]',size(d0W_betazeta_mlm___)));
disp(sprintf(' %% size(dtemplateda_wkb__): [%d,%d]',size(dtemplateda_wkb__)));
disp(sprintf(' %% size(dtemplatedb_wkb__): [%d,%d]',size(dtemplatedb_wkb__)));
disp(sprintf(' %% size(dtemplatedc_wkb__): [%d,%d]',size(dtemplatedc_wkb__)));
disp(sprintf(' %% size(d1W_betazeta_mlm___): [%d,%d,%d]',size(d1W_betazeta_mlm___)));
disp(sprintf(' %% size(ddtemplatedaa_wkb__): [%d,%d]',size(ddtemplatedaa_wkb__)));
disp(sprintf(' %% size(ddtemplatedab_wkb__): [%d,%d]',size(ddtemplatedab_wkb__)));
disp(sprintf(' %% size(ddtemplatedac_wkb__): [%d,%d]',size(ddtemplatedac_wkb__)));
disp(sprintf(' %% size(ddtemplatedbb_wkb__): [%d,%d]',size(ddtemplatedbb_wkb__)));
disp(sprintf(' %% size(ddtemplatedbc_wkb__): [%d,%d]',size(ddtemplatedbc_wkb__)));
disp(sprintf(' %% size(ddtemplatedcc_wkb__): [%d,%d]',size(ddtemplatedcc_wkb__)));
disp(sprintf(' %% size(d2W_betazeta_mlm___): [%d,%d,%d]',size(d2W_betazeta_mlm___)));
end;%if (verbose);
end;%if  flag_d0 &  flag_d1 &  flag_d2;
%%%%;
if flag_d0 & ~flag_d0W_use; d0W_betazeta_mlma____(:,:,:,1+nviewing_polar_a) = d0W_betazeta_mlm___; clear d0W_betazeta_mlm___; end;
if flag_d1 & ~flag_d1W_use; d1W_betazeta_mlma____(:,:,:,1+nviewing_polar_a) = d1W_betazeta_mlm___; clear d1W_betazeta_mlm___; end;
if flag_d2 & ~flag_d2W_use; d2W_betazeta_mlma____(:,:,:,1+nviewing_polar_a) = d2W_betazeta_mlm___; clear d2W_betazeta_mlm___; end;
if flag_d0; template_wkS__(:,1+tmp_index_) = template_wkb__; end;
if flag_d1; dtemplateda_wkS__(:,1+tmp_index_) = dtemplateda_wkb__; end;
if flag_d1; dtemplatedb_wkS__(:,1+tmp_index_) = dtemplatedb_wkb__; end;
if flag_d1; dtemplatedc_wkS__(:,1+tmp_index_) = dtemplatedc_wkb__; end;
if flag_d2; ddtemplatedaa_wkS__(:,1+tmp_index_) = ddtemplatedaa_wkb__; end;
if flag_d2; ddtemplatedab_wkS__(:,1+tmp_index_) = ddtemplatedab_wkb__; end;
if flag_d2; ddtemplatedac_wkS__(:,1+tmp_index_) = ddtemplatedac_wkb__; end;
if flag_d2; ddtemplatedbb_wkS__(:,1+tmp_index_) = ddtemplatedbb_wkb__; end;
if flag_d2; ddtemplatedbc_wkS__(:,1+tmp_index_) = ddtemplatedbc_wkb__; end;
if flag_d2; ddtemplatedcc_wkS__(:,1+tmp_index_) = ddtemplatedcc_wkb__; end;
if (verbose);
if flag_d0; disp(sprintf(' %% size(template_wkS__): [%d,%d]',size(template_wkS__))); end;
if flag_d1; disp(sprintf(' %% size(dtemplateda_wkS__): [%d,%d]',size(dtemplateda_wkS__))); end;
if flag_d1; disp(sprintf(' %% size(dtemplatedb_wkS__): [%d,%d]',size(dtemplatedb_wkS__))); end;
if flag_d1; disp(sprintf(' %% size(dtemplatedc_wkS__): [%d,%d]',size(dtemplatedc_wkS__))); end;
if flag_d2; disp(sprintf(' %% size(ddtemplatedaa_wkS__): [%d,%d]',size(ddtemplatedaa_wkS__))); end;
if flag_d2; disp(sprintf(' %% size(ddtemplatedab_wkS__): [%d,%d]',size(ddtemplatedab_wkS__))); end;
if flag_d2; disp(sprintf(' %% size(ddtemplatedac_wkS__): [%d,%d]',size(ddtemplatedac_wkS__))); end;
if flag_d2; disp(sprintf(' %% size(ddtemplatedbb_wkS__): [%d,%d]',size(ddtemplatedbb_wkS__))); end;
if flag_d2; disp(sprintf(' %% size(ddtemplatedbc_wkS__): [%d,%d]',size(ddtemplatedbc_wkS__))); end;
if flag_d2; disp(sprintf(' %% size(ddtemplatedcc_wkS__): [%d,%d]',size(ddtemplatedcc_wkS__))); end;
end;%if (verbose);
%%%%%%%%;
end;%for nviewing_polar_a=0:n_viewing_polar_a-1;
%%%%%%%%;

