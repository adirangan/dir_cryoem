function ...
[ ...
 parameter ...
,template_wkS__ ...
,a_k_u_fft3d0___ ...
,n_w ...
,n_viewing_all ...
,viewing_azimu_b_all_ ...
,viewing_polar_a_all_ ...
,viewing_weight_all_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
] = ...
fft_template_4( ...
 parameter ...
,n_x_u ...
,x_p_r_max ...
,a_x_u_form_ ...
,k_p_r_max ...
,k_eq_d ...
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
);
% Considering applying fft to a_x_u_form_ to evaluate templates on a collection of points on spherical shells. ;
% each spherical-shell has the same resolution, determined by viewing_k_eq_d and template_k_eq_d and/or n_w_max. ;
% Warning, not implemented yet. ;
% ;
% inputs: ;
% ;
% parameter = struct holding parameters. ;
% n_x_u = integer number of voxels on a side. ;
% x_p_r_max = real half-side-length of box. ;
%             We assume that the voxels extend from x_p_r_max*[-1,+1]. ;
% a_x_u_form_ = complex array of size (n_x_u,n_x_u,n_x_u). ;
% k_p_r_max = real maximum frequency. ;
% k_eq_d = real distance between points along equator on largest frequency-shell. ;
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
% ;
% outputs: ;
% ;
% template_wkS__ = complex array of templates. ;
%                  template_wkS__(1+nw+nk_p_r*n_w_max,1+nS) ;
%                  stores template value for angle-index nw, radial-index nk_p_r, ;
%                  and viewing_azimu_b = viewing_azimu_b_all_(1+nS). ;
%                  and viewing_polar_a = viewing_polar_a_all_(1+nS). ;
% a_k_u_fft3d0___ = complex array of size (n_k_u_pad,n_k_u_pad,n_k_u_pad) storing the fourier-transform of a_x_u_form_. ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

str_thisfunction = 'fft_template_4';

if nargin<1;
flag_verbose = 2; nf=1;
if (flag_verbose>0); disp(sprintf(' %% testing %s',str_thisfunction)); end;
if (flag_verbose>0); disp(sprintf(' %% calling test_CryoLike_fft_template_20241210')); end;
disp(sprintf(' %% returning')); return;
end;% if nargin<1;
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_x_u=[]; end; na=na+1;
if (nargin<1+na); x_p_r_max=[]; end; na=na+1;
if (nargin<1+na); a_x_u_form_=[]; end; na=na+1;
if (nargin<1+na); k_p_r_max=[]; end; na=na+1;
if (nargin<1+na); k_eq_d=[]; end; na=na+1;
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

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
flag_verbose=parameter.flag_verbose;
if ~isfield(parameter,'tolerance_master'); parameter.tolerance_master=1e-2; end;
tolerance_master=parameter.tolerance_master;
if (~isfield(parameter,'GB_max')); parameter.GB_max = 32; end; %<-- parameter_bookmark. ;
GB_max = parameter.GB_max;
if (~isfield(parameter,'n_order')); parameter.n_order = 3; end; %<-- parameter_bookmark. ;
n_order = parameter.n_order;
if (~isfield(parameter,'n_S_per_Sbatch')); parameter.n_S_per_Sbatch = 128; end; %<-- parameter_bookmark. ;
n_S_per_Sbatch = parameter.n_S_per_Sbatch;
if (~isfield(parameter,'n_x_u_pad_factor')); parameter.n_x_u_pad_factor = 1.0; end; %<-- parameter_bookmark. ;
n_x_u_pad_factor = parameter.n_x_u_pad_factor;

if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

half_diameter_x_c = x_p_r_max;
diameter_x_c = 2.0d0*half_diameter_x_c;
%x_u_0_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u+1); x_u_0_ = transpose(x_u_0_(1:n_x_u));
%x_u_1_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u+1); x_u_1_ = transpose(x_u_1_(1:n_x_u));
%x_u_2_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u+1); x_u_2_ = transpose(x_u_2_(1:n_x_u));
x_u_0_ = transpose(linspace(-x_p_r_max,+x_p_r_max,n_x_u)); d_x_0 = mean(diff(x_u_0_));
x_u_1_ = transpose(linspace(-x_p_r_max,+x_p_r_max,n_x_u)); d_x_1 = mean(diff(x_u_1_));
x_u_2_ = transpose(linspace(-x_p_r_max,+x_p_r_max,n_x_u)); d_x_2 = mean(diff(x_u_2_));
[x_u_0___,x_u_1___,x_u_2___] = ndgrid(x_u_0_,x_u_1_,x_u_2_); n_xxx_u = n_x_u^3;
weight_xxx_u_ = d_x_0 * d_x_1 * d_x_2 ;

%%%%%%%%;
% First determine radii of spherical grid. ;
%%%%%%%%;
str_TorL = 'L';
[ ...
 n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
] = ...
get_weight_3d_1( ...
 flag_verbose ...
,k_p_r_max ...
,k_eq_d ...
,str_TorL ...
);

%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% Then determine the viewing angles.')); end;
%%%%%%%%;
if isempty(n_viewing_all);
k_p_r_1 = 1;
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
 k_p_r_1 ...
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
if (flag_verbose>0); disp(sprintf(' %% n_viewing_all %d n_viewing_polar_a %d n_viewing_azimu_b_sum %d',n_viewing_all,n_viewing_polar_a,n_viewing_azimu_b_sum)); end;
if (flag_verbose>2); fprintf(1,' %% \t\n'); darray_printf_margin(viewing_azimu_b_all_,1,n_viewing_all,' %% viewing_azimu_b_all_: '); end;
if (flag_verbose>2); fprintf(1,' %% \t\n'); darray_printf_margin(viewing_polar_a_all_,1,n_viewing_all,' %% viewing_polar_a_all_: '); end;
if (flag_verbose>2); fprintf(1,' %% \t\n'); darray_printf_margin(viewing_weight_all_,1,n_viewing_all,' %% viewing_weight_all_: '); end;
if (flag_verbose>2); fprintf(1,' %% \t\n'); darray_printf_margin(viewing_polar_a_,1,n_viewing_polar_a,' %% viewing_polar_a_: '); end;
if (flag_verbose>2); fprintf(1,' %% \t\n'); darray_printf_margin(n_viewing_azimu_b_,1,n_viewing_polar_a,' %% n_viewing_azimu_b_: '); end;
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% Now determine the points along each equatorial plane (i.e., the points for each template).')); end;
%%%%%%%%;
n_w = 0;
if (template_k_eq_d>0);
k_p_r_1 = 1;
n_equator = 3+round(2*pi*k_p_r_1/template_k_eq_d);
n_polar_a = 3+round(n_equator/2);
n_w = 2*n_polar_a;
end;%if (template_k_eq_d>0);
if (template_k_eq_d<=0);
%n_w = max(6,n_w_0in);
n_w = n_w_0in; %<-- no minimum. ;
end;%if (template_k_eq_d<=0);
if (flag_verbose>0); disp(sprintf(' %% n_w %d',n_w)); end;
n_w_max = n_w; n_w_sum = n_w_max*n_k_p_r; n_w_ = n_w_max*ones(n_k_p_r,1); n_w_csum_ = cumsum([0;n_w_]);
%%%%%%%%;

%%%%%%%%;
% Set up inner gamma_z for the templates. ;
%%%%%%%%;
gamma_z_ = zeros(n_w,1); gamma_z_ = transpose(linspace(0,2*pi,n_w+1)); gamma_z_ = gamma_z_(1:n_w);
cc_ = cos(gamma_z_); sc_ = sin(gamma_z_);
if (flag_verbose>2); fprintf(1,' %% \t\n'); darray_printf_margin(gamma_z_,1,n_w,' %% gamma_z_: '); end;
if (flag_verbose>2); fprintf(1,' %% \t\n'); darray_printf_margin(cc_,1,n_w,' %% cc_: '); end;
if (flag_verbose>2); fprintf(1,' %% \t\n'); darray_printf_margin(sc_,1,n_w,' %% sc_: '); end;
%%%%%%%%;

%%%%;
% copy parameters. ;
%%%%;
n_S = n_viewing_all;
viewing_polar_a_S_ = viewing_polar_a_all_;
viewing_azimu_b_S_ = viewing_azimu_b_all_;
viewing_gamma_z_S_ = zeros(n_S,1);
%%%%%%%%;
% Note that the convention in cg_rhs_2 is to ;
% subtract viewing_gamma_z from inplane_gamma_z. ;
%%%%%%%%;
[ ...
 k_p_polar_a_wS__ ...
,k_p_azimu_b_wS__ ...
,k_c_0_wS__ ...
,k_c_1_wS__ ...
,k_c_2_wS__ ...
] = ...
cg_rhs_2( ...
 n_S ...
,n_w_max ...
,viewing_polar_a_S_ ...
,viewing_azimu_b_S_ ...
,viewing_gamma_z_S_ ...
);
%%%%%%%%;
k_c_0_wkS___ = bsxfun(@times,reshape(k_p_r_,[1,n_k_p_r,1]),reshape(k_c_0_wS__,[n_w_max,1,n_S]));
k_c_1_wkS___ = bsxfun(@times,reshape(k_p_r_,[1,n_k_p_r,1]),reshape(k_c_1_wS__,[n_w_max,1,n_S]));
k_c_2_wkS___ = bsxfun(@times,reshape(k_p_r_,[1,n_k_p_r,1]),reshape(k_c_2_wS__,[n_w_max,1,n_S]));

%%%%%%%%;
% With this in mind we check k_p_r_max. ;
%%%%%%%%;
flag_check_k = (k_p_r_max< (n_x_u-1)/(4*x_p_r_max)) ;
if (flag_verbose>0); disp(sprintf(' %% flag_check_k %d: k_p_r_max %0.6f (n_x_u-1)/(4*x_p_r_max) %0.6f',flag_check_k,k_p_r_max,(n_x_u-1)/(4*x_p_r_max))); end;
if ~flag_check_k; disp(sprintf(' %% Warning: flag_check_k %d: k_p_r_max %0.6f >= (n_x_u-1)/(4*x_p_r_max) %0.6f',flag_check_k,k_p_r_max,(n_x_u-1)/(4*x_p_r_max))); end;

%%%%%%%%;
% Now we create a_k_u_fft3d0___. ;
%%%%%%%%;
n_x_u_pad = n_x_u_pad_factor*n_x_u;
tmp_2pik_fft1d0_ = transpose(periodize([0:n_x_u_pad-1],-n_x_u_pad/2,+n_x_u_pad/2))*(2*pi)/(2*x_p_r_max)*(n_x_u-1)/(n_x_u_pad);
tmp_2pik_fft1d0_lim_ = [min(tmp_2pik_fft1d0_),max(tmp_2pik_fft1d0_)];
[tmp_2pik_0_fft3d0___,tmp_2pik_1_fft3d0___,tmp_2pik_2_fft3d0___] = ndgrid(tmp_2pik_fft1d0_,tmp_2pik_fft1d0_,tmp_2pik_fft1d0_);
tmp_isgn = -1;
if tmp_isgn == -1; a_k_u_fft3d0___ =       1.0  *fftshift( fftn(reshape(a_x_u_form_,[n_x_u,n_x_u,n_x_u]),[n_x_u_pad,n_x_u_pad,n_x_u_pad]).*exp(-i*tmp_isgn*(tmp_2pik_0_fft3d0___+tmp_2pik_1_fft3d0___+tmp_2pik_2_fft3d0___)*x_p_r_max)); end;
if tmp_isgn == +1; a_k_u_fft3d0___ = n_x_u_pad^3*fftshift(ifftn(reshape(a_x_u_form_,[n_x_u,n_x_u,n_x_u]),[n_x_u_pad,n_x_u_pad,n_x_u_pad]).*exp(-i*tmp_isgn*(tmp_2pik_0_fft3d0___+tmp_2pik_1_fft3d0___+tmp_2_2pik_fft3d0___)*x_p_r_max)); end;

%%%%%%%%;
% Now we create the interpolation operator. ;
%%%%%%%%;
GB_per_Sbatch = (1+n_order)^4 * n_w_max*n_k_p_r*n_S_per_Sbatch*16/1e9; %<-- n_nodex^4*n_scatter used for flag_dd. ;
if (flag_verbose>0); disp(sprintf(' %% GB_per_Sbatch %0.2f',GB_per_Sbatch)); end;
if (GB_per_Sbatch>=GB_max); disp(sprintf(' %% Warning, GB_per_Sbatch %d',GB_per_Sbatch)); end;
scatter_from_tensor_s012__ = sparse([],[],[],n_w_max*n_k_p_r*0,n_x_u_pad^3,0);
n_Sbatch = ceil(n_S/n_S_per_Sbatch);
if (flag_verbose>0); disp(sprintf(' %% n_Sbatch %d',n_Sbatch)); end;
for nSbatch=0:n_Sbatch-1;
index_S_in_Sbatch_ = nSbatch*n_S_per_Sbatch + (0:n_S_per_Sbatch-1);
index_S_in_Sbatch_ = index_S_in_Sbatch_(find(index_S_in_Sbatch_<n_S)); n_S_sub = numel(index_S_in_Sbatch_);
if (flag_verbose>0); disp(sprintf(' %% nSbatch %d/%d index_S_in_Sbatch_ %d-->%d',nSbatch,n_Sbatch,index_S_in_Sbatch_(1+0),index_S_in_Sbatch_(1+n_S_sub-1))); end;
if (flag_verbose>0 & mod(nSbatch,32)==0); disp(sprintf(' %% nSbatch %d/%d index_S_in_Sbatch_ %d-->%d',nSbatch,n_Sbatch,index_S_in_Sbatch_(1+0),index_S_in_Sbatch_(1+n_S_sub-1))); end;
if (n_S_sub>0);
tmp_t = tic();
k_c_0_wkS_sub___ = k_c_0_wkS___(:,:,1+index_S_in_Sbatch_);
k_c_1_wkS_sub___ = k_c_1_wkS___(:,:,1+index_S_in_Sbatch_);
k_c_2_wkS_sub___ = k_c_2_wkS___(:,:,1+index_S_in_Sbatch_);
[ ...
 scatter_from_tensor_sub_s012__ ...
] = ...
volume_k_c_scatter_from_tensor_interpolate_n_8( ...
 n_order ...
,n_x_u_pad ...
,n_x_u_pad ...
,n_x_u_pad ...
,tmp_2pik_fft1d0_lim_ ...
,tmp_2pik_fft1d0_lim_ ...
,tmp_2pik_fft1d0_lim_ ...
,n_w_max*n_k_p_r*n_S_sub ...
,2*pi*k_c_0_wkS_sub___(:) ...
,2*pi*k_c_1_wkS_sub___(:) ...
,2*pi*k_c_2_wkS_sub___(:) ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% nSbatch %d/%d: scatter_from_tensor_sub_s012__: %0.6fs',nSbatch,n_Sbatch,tmp_t)); end;
scatter_from_tensor_s012__ = cat(1,scatter_from_tensor_s012__,scatter_from_tensor_sub_s012__);
end;%if (n_S_sub>0);
end;%for nSbatch=0:n_Sbatch-1;

%%%%%%%%;
% Now use the interpolation operator to build the templates. ;
%%%%%%%%;
tmp_t = tic;
S_k_p_wkS__ = reshape(scatter_from_tensor_s012__*a_k_u_fft3d0___(:).*weight_xxx_u_(:),[n_w_sum,n_S]);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% %% scatter_from_tensor_s012__: S_3_k_p_wkS__ time %0.2fs',tmp_t)); end;
%%%%;
for nS=0:n_S-1;
S_k_p_wkS__(:,1+nS) = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_wkS__(:,1+nS),-viewing_gamma_z_S_(1+nS));
end;%for nS=0:n_S-1;
%%%%;
template_wkS__ = S_k_p_wkS__;

if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
