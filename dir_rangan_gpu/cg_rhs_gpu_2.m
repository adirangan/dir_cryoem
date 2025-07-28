function ...
[ ...
 k_p_polar_a_gpu_wM__ ...
,k_p_azimu_b_gpu_wM__ ...
,k_c_0_gpu_wM__ ...
,k_c_1_gpu_wM__ ...
,k_c_2_gpu_wM__ ...
,k_p_r01_gpu_wM__ ...
,dtau_k_p_polar_a_gpu_wM3___ ...
,dtau_k_p_azimu_b_gpu_wM3___ ...
,dtau_k_c_0_gpu_wM3___ ...
,dtau_k_c_1_gpu_wM3___ ...
,dtau_k_c_2_gpu_wM3___ ...
,dtau_k_p_r01_gpu_wM3___ ...
,dtau_dtau_k_p_polar_a_gpu_wM33____ ...
,dtau_dtau_k_p_azimu_b_gpu_wM33____ ...
,dtau_dtau_k_c_0_gpu_wM33____ ...
,dtau_dtau_k_c_1_gpu_wM33____ ...
,dtau_dtau_k_c_2_gpu_wM33____ ...
,dtau_dtau_k_p_r01_gpu_wM33____ ...
] = ...
cg_rhs_gpu_2( ...
 n_M ...
,n_w ...
,viewing_polar_a_M_ ...
,viewing_azimu_b_M_ ...
,viewing_gamma_z_M_ ...
,flag_cleanup ...
);
%%%%%%%%;
% returns right-hand-side for a least-squares problem, ;
% comprising the values of M_k_p_gpu_wM__, along with their locations (in k_p_polar_a and k_p_azimu_b_). ;
% In this context, M_k_p_gpu_wM__ refers to a complex array of size n_w-by-n_M. ;
% M_k_p_gpu_wM__(1+nw,1+nM) contains the (complex) value of image-ring nM at angle gamma_z_(1+nw) = 2*pi*nw/n_w. ;
%%%%;
% Inputs: ;
% n_M = integer number of image-rings (each image is a single ring). ;
% n_w = integer number of gamma_z values for each image-ring. ;
%       We assume gamma_z_ = linspace(0,2*pi,n_w+1); gamma_z_ = gamma_z(1:end-1);
% viewing_polar_a_M_ = double array of size n_M. viewing_polar_a_M_(1+nM) contains the viewing_polar_a for image-ring nM. ;
% viewing_azimu_b_M_ = double array of size n_M. viewing_azimu_b_M_(1+nM) contains the viewing_azimu_b for image-ring nM. ;
% viewing_gamma_z_M_ = double array of size n_M. viewing_gamma_z_M_(1+nM) contains the viewing_gamma_z for image-ring nM. ;
% flag_cleanup = integer (default 1). if set will threshold / cleanup large derivative. ;
%%%%;
% Outputs: ;
% k_p_polar_a_gpu_wM__ = double array of size n_w-by-n_M. k_p_polar_a_gpu_wM__(1+nw,1+nM) lists the k_p_polar_a for point nw of image-ring nM. ;
% k_p_azimu_b_gpu_wM__ = double array of size n_w-by-n_M. k_p_azimu_b_gpu_wM__(1+nw,1+nM) lists the k_p_azimu_b for point nw of image-ring nM. ;
% k_c_0_gpu_wM__ = double array of size n_w-by-n_M. k_c_0_gpu_wM__(1+nw,1+nM) lists the k_c_0 for point nw of image-ring nM. ;
% k_c_1_gpu_wM__ = double array of size n_w-by-n_M. k_c_0_gpu_wM__(1+nw,1+nM) lists the k_c_1 for point nw of image-ring nM. ;
% k_c_2_gpu_wM__ = double array of size n_w-by-n_M. k_c_0_gpu_wM__(1+nw,1+nM) lists the k_c_2 for point nw of image-ring nM. ;
% ;
% first-derivates (e.g., dtau_k_c_0_gpu_wM3___) are organized such that: ;
% dtau_k_c_0_gpu_wM3___(:,:,1+0) is the derivative w.r.t. viewing_polar_a. ;
% dtau_k_c_0_gpu_wM3___(:,:,1+1) is the derivative w.r.t. viewing_azimu_b. ;
% dtau_k_c_0_gpu_wM3___(:,:,1+2) is the derivative w.r.t. viewing_gamma_z. ;
% ;
% second-derivates (e.g., dtau_dtau_k_c_0_gpu_wM33____) are organized such that: ;
% dtau_dtau_k_c_0_gpu_wM33____(:,:,1+0,1+0) is the derivative w.r.t. viewing_polar_a and viewing_polar_a. ;
% dtau_dtau_k_c_0_gpu_wM33____(:,:,1+0,1+1) is the derivative w.r.t. viewing_polar_a and viewing_azimu_b. ;
% dtau_dtau_k_c_0_gpu_wM33____(:,:,1+0,1+2) is the derivative w.r.t. viewing_polar_a and viewing_gamma_z. ;
% dtau_dtau_k_c_0_gpu_wM33____(:,:,1+1,1+0) is the derivative w.r.t. viewing_azimu_b and viewing_polar_a. ;
% dtau_dtau_k_c_0_gpu_wM33____(:,:,1+1,1+1) is the derivative w.r.t. viewing_azimu_b and viewing_azimu_b. ;
% dtau_dtau_k_c_0_gpu_wM33____(:,:,1+1,1+2) is the derivative w.r.t. viewing_azimu_b and viewing_gamma_z. ;
% dtau_dtau_k_c_0_gpu_wM33____(:,:,1+2,1+0) is the derivative w.r.t. viewing_gamma_z and viewing_polar_a. ;
% dtau_dtau_k_c_0_gpu_wM33____(:,:,1+2,1+1) is the derivative w.r.t. viewing_gamma_z and viewing_azimu_b. ;
% dtau_dtau_k_c_0_gpu_wM33____(:,:,1+2,1+2) is the derivative w.r.t. viewing_gamma_z and viewing_gamma_z. ;
%%%%;
% For this calculation we use the same assumptions made in get_template_0.m: ;
% The general formula used here is as follows. ;
% let sa and ca be sin(viewing_polar_a) and cos(viewing_polar_a), respectively. ;
% let sb and cb be sin(viewing_azimu_b) and cos(viewing_azimu_b), respectively. ;
% let sc and cc be sin(inplane_gamma_z) and cos(inplane_gamma_z), respectively. ;
% And rotation by viewing_azimu_b about the +z-axis is represented as: ;
% Rz(viewing_azimu_b) = ;
% [ +cb -sb 0 ] ;
% [ +sb +cb 0 ] ;
% [  0   0  1 ] ;
% And rotation by viewing_polar_a about the +y-axis is represented as: ;
% Ry(viewing_polar_a) = ;
% [ +ca 0 +sa ] ;
% [  0  1  0  ] ;
% [ -sa 0 +ca ] ;
% And rotation by inplane_gamma_z about the +z-axis is represented as: ;
% Rz(inplane_gamma_z) = ;
% [ +cc -sc 0 ] ;
% [ +sc +cc 0 ] ;
% [  0   0  1 ] ;
% Which, collectively, implies that under the transform: ;
% Rz(viewing_azimu_b) * Ry(viewing_polar_a) * Rz(inplane_gamma_z), ;
% Which is the same as: ;
% [ +cb -sb 0 ] [ +ca*cc -ca*sc +sa ]   [ +cb*ca*cc - sb*sc , -cb*ca*sc -sb*cc , +cb*sa ];
% [ +sb +cb 0 ] [ +sc    +cc    0   ] = [ +sb*ca*cc + cb*sc , -sb*ca*sc +cb*cc , +sb*sa ];
% [  0   0  1 ] [ -sa*cc +sa*sc +ca ]   [ -sa*cc            , +sa*sc           , +ca    ];
% the point [1;0;0] is mapped to: ;
% [ template_k_c_0 ; template_k_c_1 ; template_k_c_2 ] = [ +cb*ca*cc - sb*sc ; +sb*ca*cc + cb*sc ; -sa*cc ];
%%%%;
% In other words: 
% the point M_k_p_w_(1+nw)==M_k_p_gpu_wM__(1+nw,1+nM), ;
% corresponding to the nw-th point on the nM-th image-ring, ;
% is generated using inplane_gamma_z = 2*pi*nw/n_w, ;
% as well as viewing_polar_a and viewing_azimu_b. ;
% Thus, the k_c_ location for that point (in 3d) is: ;
% k_c_ = [ +cb*ca*cc - sb*sc ; +sb*ca*cc + cb*sc ; -sa*cc ]. ;
% Correspondingly, the k_p_r01, k_p_polar_a and k_p_azimu_b for that point can be calculated as: ;
% k_p_r01 = sqrt(k_c_(1+0)^2 + k_c_(1+1)^2) = ;
% k_p_polar_a = atan2(k_p_r01,k_c_(1+2));
% k_p_azimu_b = atan2(k_c_(1+1),k_c_(1+0));
%%%%;
% Here we adopt the convention that the input viewing_gamma_z for each image-ring ;
% should be applied as an inplane rotation of the image-ring itself before calculating k_c_. ;
% Thus, the viewing_gamma_z is simply subtracted from the inplane_gamma_z associated with each point nw in that image-ring. ;
%%%%;
% Note that: ;
% cos(template_polar_a) = -sa*cc ;
% template_azimu_b = atan2( +sb*ca*cc + cb*sc , *cb*ca*cc - sb*sc );
% implying that: ;
% \partial_{viewing_polar_a} cos(template_polar_a) = -ca*cc ;
% \partial_{viewing_azimu_b} cos(template_polar_a) = 0 ;
% \partial_{viewing_gamma_z} cos(template_polar_a) = +sa*sc ;
% and since \partial_{y/x} atan(y/x) = 1/(1+(y/x)^2) = x^2/(x^2+y^2) ;
% and \partial_{z} y(z)/x(z) = (y'x-yx')/x^2, ;
% \partial_{viewing_polar_a} template_azimu_b = (\partial_{viewing_polar_a} y * x - y * \partial_{viewing_polar_a} x ) / (x^2 + y^2) ;
% \partial_{viewing_azimu_b} template_azimu_b = (\partial_{viewing_azimu_b} y * x - y * \partial_{viewing_azimu_b} x ) / (x^2 + y^2) ;
% \partial_{viewing_gamma_z} template_azimu_b = (\partial_{viewing_gamma_z} y * x - y * \partial_{viewing_gamma_z} x ) / (x^2 + y^2) ;
% where:
% \partial_{viewing_polar_a} x = -cb*sa*cc ;
% \partial_{viewing_polar_a} y = -sb*sa*cc ;
% \partial_{viewing_azimu_b} x = -sb*ca*cc - cb*sc ;
% \partial_{viewing_azimu_b} y = +cb*ca*cc - sb*sc ;
% \partial_{viewing_gamma_z} x = -cb*ca*sc - sb*cc ;
% \partial_{viewing_gamma_z} y = -sb*ca*sc + cb*cc ;
%%%%%%%%;

str_thisfunction = 'cg_rhs_gpu_2';
flag_verbose=0;

if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

na=0;
if (nargin<1+na); n_M=[]; end; na=na+1;
if (nargin<1+na); n_w=[]; end; na=na+1;
if (nargin<1+na); viewing_polar_a_M_=[]; end; na=na+1;
if (nargin<1+na); viewing_azimu_b_M_=[]; end; na=na+1;
if (nargin<1+na); viewing_gamma_z_M_=[]; end; na=na+1;
if (nargin<1+na); flag_cleanup=[]; end; na=na+1;

if isempty(viewing_gamma_z_M_); viewing_gamma_z_M_=zeros(n_M,1); end;


f_zero = gpuArray( single(0.0));
viewing_polar_a_gpu_M_ = gpuArray( viewing_polar_a_M_ );
viewing_azimu_b_gpu_M_ = gpuArray( viewing_azimu_b_M_ );
viewing_gamma_z_gpu_M_ = gpuArray( viewing_gamma_z_M_ );

flag_1 = 1;
flag_d = (nargout>=1+1*6) ;
flag_dd = (nargout>=1+2*6) ;
tolerance_cg_rhs = 1e-16;
tolerance_pole = sqrt(tolerance_cg_rhs);
tolerance_cg_rhs_upb = inv(sqrt(1e-16));
if isempty(flag_cleanup); flag_cleanup=1; end;

%%%%%%%%;
if flag_1;
%%%%;
sa_gpu_wM__ = repmat(reshape(sin(viewing_polar_a_gpu_M_),[1,n_M]),[n_w,1]); ca_gpu_wM__ = repmat(reshape(cos(viewing_polar_a_gpu_M_),[1,n_M]),[n_w,1]);
sb_gpu_wM__ = repmat(reshape(sin(viewing_azimu_b_gpu_M_),[1,n_M]),[n_w,1]); cb_gpu_wM__ = repmat(reshape(cos(viewing_azimu_b_gpu_M_),[1,n_M]),[n_w,1]);
viewing_gamma_z_gpu_M_ = reshape(viewing_gamma_z_gpu_M_,[1,n_M]); inplane_gamma_z_gpu_w_ = gpuArray( reshape(2*pi*[0:n_w-1]/max(1,n_w),[n_w,1]) );
combine_gamma_z_gpu_wM__ = repmat(inplane_gamma_z_gpu_w_,[1,n_M]) - repmat(viewing_gamma_z_gpu_M_,[n_w,1]);
sc_gpu_wM__ = sin(combine_gamma_z_gpu_wM__); cc_gpu_wM__ = cos(combine_gamma_z_gpu_wM__);
k_c_0_gpu_wM__ = +cb_gpu_wM__.*ca_gpu_wM__.*cc_gpu_wM__ - sb_gpu_wM__.*sc_gpu_wM__ ; 
k_c_1_gpu_wM__ = +sb_gpu_wM__.*ca_gpu_wM__.*cc_gpu_wM__ + cb_gpu_wM__.*sc_gpu_wM__ ; 
k_c_2_gpu_wM__ = -sa_gpu_wM__.*cc_gpu_wM__ ;
k_p_r01_bkp_gpu_wM__ = sqrt(k_c_0_gpu_wM__.^2 + k_c_1_gpu_wM__.^2);
k_p_r01_gpu_wM__ = sqrt((ca_gpu_wM__.*cc_gpu_wM__).^2 + (sc_gpu_wM__).^2);
fnorm_disp(flag_verbose,'k_p_r01_gpu_wM__',k_p_r01_gpu_wM__,'k_p_r01_bkp_gpu_wM__',k_p_r01_bkp_gpu_wM__,' %<-- should be <1e-12');
k_p_polar_a_gpu_wM__ = atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__);
k_p_azimu_b_gpu_wM__ = atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__);
%%%%%%%%;
% cleanup. ;
%%%%%%%%;
for nM=0:n_M-1;
if abs(periodize(viewing_polar_a_gpu_M_-pi/2,-pi/2,+pi/2)<tolerance_pole);
tmp_index_0_pos_ = efind(k_c_0_gpu_wM__(:,1+nM)> +tolerance_pole);
tmp_index_0_neg_ = efind(k_c_0_gpu_wM__(:,1+nM)< -tolerance_pole);
flag_0_use = fnorm(k_c_0_gpu_wM__(:,1+nM))>=0.25;
tmp_index_1_pos_ = efind(k_c_1_gpu_wM__(:,1+nM)> +tolerance_pole);
tmp_index_1_neg_ = efind(k_c_1_gpu_wM__(:,1+nM)< -tolerance_pole);
flag_1_use = fnorm(k_c_1_gpu_wM__(:,1+nM))>=0.25;
tmp_index_2_pos_ = efind(abs(k_c_2_gpu_wM__(:,1+nM)-1)< tolerance_pole);
tmp_index_2_neg_ = efind(abs(k_c_2_gpu_wM__(:,1+nM)+1)< tolerance_pole);
assert(flag_0_use | flag_1_use);
if flag_0_use
tmp_azimu_b_pos = mean(k_p_azimu_b_gpu_wM__(1+tmp_index_0_pos_,1+nM));
tmp_azimu_b_neg = mean(k_p_azimu_b_gpu_wM__(1+tmp_index_0_neg_,1+nM));
else
tmp_azimu_b_pos = mean(k_p_azimu_b_gpu_wM__(1+tmp_index_1_pos_,1+nM));
tmp_azimu_b_neg = mean(k_p_azimu_b_gpu_wM__(1+tmp_index_1_neg_,1+nM));
end;%if flag_0_use;
tmp_azimu_b_mid = periodize(0.5*(tmp_azimu_b_pos + tmp_azimu_b_neg),0,2*pi);
k_p_azimu_b_gpu_wM__(1+tmp_index_2_pos_,1+nM) = periodize(tmp_azimu_b_mid + 0*pi,0,2*pi);
k_p_azimu_b_gpu_wM__(1+tmp_index_2_neg_,1+nM) = periodize(tmp_azimu_b_mid + 1*pi,0,2*pi);
end;%if abs(periodize(viewing_polar_a_gpu_M_-pi/2,-pi/2,+pi/2)<tolerance_pole);
end;%for nM=0:n_M-1;
%%%%;
end;%if flag_1;
%%%%%%%%;

%%%%%%%%;
if flag_d;
%%%%;
r01 = @(r0,r1) sqrt(r0.^2 + r1.^2) ;
d0_r01 = @(r0,r1) r0./max(tolerance_cg_rhs,r01(r0,r1));
d1_r01 = @(r0,r1) r1./max(tolerance_cg_rhs,r01(r0,r1));
%z = @(y,x) y./x ;
dx_z = @(y,x) -y./max(tolerance_cg_rhs^2,x.^2) ;
%dy_z = @(y,x) +1./x ; %<-- Warning: signed divide. ;
dz_atan2 = @(y,x) x.^2 ./ max(tolerance_cg_rhs^2,x.^2 + y.^2) ;
dz_atan2_times_dy_z = @(y,x) +x.^1 ./ max(tolerance_cg_rhs^2,x.^2 + y.^2) ; %<-- avoids signed divide. ;
dz_atan2_times_dx_z = @(y,x) -y.^1 ./ max(tolerance_cg_rhs^2,x.^2 + y.^2) ; %<-- avoids signed divide. ;
%dx_atan2 = @(y,x) dz_atan2(y,x).*dx_z(y,x) ; %<-- uses signed divide. ;
dx_atan2 = @(y,x) dz_atan2_times_dx_z(y,x) ; %<-- avoids signed divide. ;
%dy_atan2 = @(y,x) dz_atan2(y,x).*dy_z(y,x) ; %<-- uses signed divide. ;
dy_atan2 = @(y,x) dz_atan2_times_dy_z(y,x) ; %<-- avoids signed divide. ;
da_sa_gpu_wM__ = repmat(reshape(+cos(viewing_polar_a_gpu_M_),[1,n_M]),[n_w,1]); da_ca_gpu_wM__ = repmat(reshape(-sin(viewing_polar_a_gpu_M_),[1,n_M]),[n_w,1]);
db_sb_gpu_wM__ = repmat(reshape(+cos(viewing_azimu_b_gpu_M_),[1,n_M]),[n_w,1]); db_cb_gpu_wM__ = repmat(reshape(-sin(viewing_azimu_b_gpu_M_),[1,n_M]),[n_w,1]);
dc_sc_gpu_wM__ = -cos(combine_gamma_z_gpu_wM__); dc_cc_gpu_wM__ = +sin(combine_gamma_z_gpu_wM__); %<-- note sign switch. ;
da_k_c_0_gpu_wM__ = +cb_gpu_wM__.*da_ca_gpu_wM__.*cc_gpu_wM__ ;
db_k_c_0_gpu_wM__ = +db_cb_gpu_wM__.*ca_gpu_wM__.*cc_gpu_wM__ - db_sb_gpu_wM__.*sc_gpu_wM__ ; 
dc_k_c_0_gpu_wM__ = +cb_gpu_wM__.*ca_gpu_wM__.*dc_cc_gpu_wM__ - sb_gpu_wM__.*dc_sc_gpu_wM__ ;
da_k_c_1_gpu_wM__ = +sb_gpu_wM__.*da_ca_gpu_wM__.*cc_gpu_wM__ ;
db_k_c_1_gpu_wM__ = +db_sb_gpu_wM__.*ca_gpu_wM__.*cc_gpu_wM__ + db_cb_gpu_wM__.*sc_gpu_wM__ ; 
dc_k_c_1_gpu_wM__ = +sb_gpu_wM__.*ca_gpu_wM__.*dc_cc_gpu_wM__ + cb_gpu_wM__.*dc_sc_gpu_wM__ ; 
da_k_c_2_gpu_wM__ = -da_sa_gpu_wM__.*cc_gpu_wM__ ;
db_k_c_2_gpu_wM__ = gpuArray( zeros(n_w,n_M) );
dc_k_c_2_gpu_wM__ = -sa_gpu_wM__.*dc_cc_gpu_wM__ ;
%da_k_p_r01_gpu_wM__ = (k_c_0_gpu_wM__.*da_k_c_0_gpu_wM__ + k_c_1_gpu_wM__.*da_k_c_1_gpu_wM__)./k_p_r01_gpu_wM__;
%db_k_p_r01_gpu_wM__ = (k_c_0_gpu_wM__.*db_k_c_0_gpu_wM__ + k_c_1_gpu_wM__.*db_k_c_1_gpu_wM__)./k_p_r01_gpu_wM__;
%dc_k_p_r01_gpu_wM__ = (k_c_0_gpu_wM__.*dc_k_c_0_gpu_wM__ + k_c_1_gpu_wM__.*dc_k_c_1_gpu_wM__)./k_p_r01_gpu_wM__;
da_k_p_r01_gpu_wM__ = d0_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*da_k_c_0_gpu_wM__ + d1_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*da_k_c_1_gpu_wM__ ;
db_k_p_r01_gpu_wM__ = d0_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*db_k_c_0_gpu_wM__ + d1_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*db_k_c_1_gpu_wM__ ;
dc_k_p_r01_gpu_wM__ = d0_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*dc_k_c_0_gpu_wM__ + d1_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*dc_k_c_1_gpu_wM__ ;
da_k_p_polar_a_gpu_wM__ = dy_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*da_k_p_r01_gpu_wM__ + dx_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*da_k_c_2_gpu_wM__;
db_k_p_polar_a_gpu_wM__ = dy_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*db_k_p_r01_gpu_wM__ + dx_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*db_k_c_2_gpu_wM__;
dc_k_p_polar_a_gpu_wM__ = dy_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*dc_k_p_r01_gpu_wM__ + dx_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*dc_k_c_2_gpu_wM__;
da_k_p_azimu_b_gpu_wM__ = dy_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*da_k_c_1_gpu_wM__ + dx_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*da_k_c_0_gpu_wM__;
db_k_p_azimu_b_gpu_wM__ = dy_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*db_k_c_1_gpu_wM__ + dx_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*db_k_c_0_gpu_wM__;
dc_k_p_azimu_b_gpu_wM__ = dy_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*dc_k_c_1_gpu_wM__ + dx_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*dc_k_c_0_gpu_wM__;
%%%%%%%%;
% cleanup. ;
%%%%%%%%;
if flag_cleanup;
%%;
for ntype=0:3-1;
na=0;
if ntype==na; tmp_gpu_wM__ = da_k_p_r01_gpu_wM__; tmp_str_ = 'da_k_p_r01_gpu_wM__'; end; na=na+1;
if ntype==na; tmp_gpu_wM__ = db_k_p_r01_gpu_wM__; tmp_str_ = 'db_k_p_r01_gpu_wM__'; end; na=na+1;
if ntype==na; tmp_gpu_wM__ = dc_k_p_r01_gpu_wM__; tmp_str_ = 'dc_k_p_r01_gpu_wM__'; end; na=na+1;
tmp_index_ = efind(abs(tmp_gpu_wM__)> tolerance_cg_rhs_upb); tmp_n_index = numel(tmp_index_);
if tmp_n_index> 0;
if (flag_verbose>1); disp(sprintf(' %% cleanup %s: tmp_n_index %d',tmp_str_,tmp_n_index)); end;
tmp_clean_gpu_wM__ = 0.5*(circshift(tmp_gpu_wM__,+1,1)+circshift(tmp_gpu_wM__,-1,1));
tmp_gpu_wM__(1+tmp_index_) = tmp_clean_gpu_wM__(1+tmp_index_);
end;%if tmp_n_index> 0;
na=0;
if ntype==na; da_k_p_r01_gpu_wM__ = tmp_gpu_wM__; end; na=na+1;
if ntype==na; db_k_p_r01_gpu_wM__ = tmp_gpu_wM__; end; na=na+1;
if ntype==na; dc_k_p_r01_gpu_wM__ = tmp_gpu_wM__; end; na=na+1;
end;%for ntype=0:3-1;
%%;
for ntype=0:3-1;
na=0;
if ntype==na; tmp_gpu_wM__ = da_k_p_polar_a_gpu_wM__; tmp_str_ = 'da_k_p_polar_a_gpu_wM__'; end; na=na+1;
if ntype==na; tmp_gpu_wM__ = db_k_p_polar_a_gpu_wM__; tmp_str_ = 'db_k_p_polar_a_gpu_wM__'; end; na=na+1;
if ntype==na; tmp_gpu_wM__ = dc_k_p_polar_a_gpu_wM__; tmp_str_ = 'dc_k_p_polar_a_gpu_wM__'; end; na=na+1;
tmp_index_ = efind(abs(tmp_gpu_wM__)> tolerance_cg_rhs_upb); tmp_n_index = numel(tmp_index_);
if tmp_n_index> 0;
if (flag_verbose>1); disp(sprintf(' %% cleanup %s: tmp_n_index %d',tmp_str_,tmp_n_index)); end;
tmp_clean_gpu_wM__ = 0.5*(circshift(tmp_gpu_wM__,+1,1)+circshift(tmp_gpu_wM__,-1,1));
tmp_gpu_wM__(1+tmp_index_) = tmp_clean_gpu_wM__(1+tmp_index_);
end;%if tmp_n_index> 0;
na=0;
if ntype==na; da_k_p_polar_a_gpu_wM__ = tmp_gpu_wM__; end; na=na+1;
if ntype==na; db_k_p_polar_a_gpu_wM__ = tmp_gpu_wM__; end; na=na+1;
if ntype==na; dc_k_p_polar_a_gpu_wM__ = tmp_gpu_wM__; end; na=na+1;
end;%for ntype=0:3-1;
%%;
for ntype=0:3-1;
na=0;
if ntype==na; tmp_gpu_wM__ = da_k_p_azimu_b_gpu_wM__; tmp_str_ = 'da_k_p_azimu_b_gpu_wM__'; end; na=na+1;
if ntype==na; tmp_gpu_wM__ = db_k_p_azimu_b_gpu_wM__; tmp_str_ = 'db_k_p_azimu_b_gpu_wM__'; end; na=na+1;
if ntype==na; tmp_gpu_wM__ = dc_k_p_azimu_b_gpu_wM__; tmp_str_ = 'dc_k_p_azimu_b_gpu_wM__'; end; na=na+1;
tmp_index_ = efind(abs(tmp_gpu_wM__)> tolerance_cg_rhs_upb); tmp_n_index = numel(tmp_index_);
if tmp_n_index> 0;
if (flag_verbose>1); disp(sprintf(' %% cleanup %s: tmp_n_index %d',tmp_str_,tmp_n_index)); end;
tmp_clean_gpu_wM__ = 0.5*(circshift(tmp_gpu_wM__,+1,1)+circshift(tmp_gpu_wM__,-1,1));
tmp_gpu_wM__(1+tmp_index_) = tmp_clean_gpu_wM__(1+tmp_index_);
end;%if tmp_n_index> 0;
na=0;
if ntype==na; da_k_p_azimu_b_gpu_wM__ = tmp_gpu_wM__; end; na=na+1;
if ntype==na; db_k_p_azimu_b_gpu_wM__ = tmp_gpu_wM__; end; na=na+1;
if ntype==na; dc_k_p_azimu_b_gpu_wM__ = tmp_gpu_wM__; end; na=na+1;
end;%for ntype=0:3-1;
%%;
end;%if flag_cleanup;
%%%%%%%%;
dtau_k_c_0_gpu_wM3___(:,:,1+0) = da_k_c_0_gpu_wM__;
dtau_k_c_0_gpu_wM3___(:,:,1+1) = db_k_c_0_gpu_wM__;
dtau_k_c_0_gpu_wM3___(:,:,1+2) = dc_k_c_0_gpu_wM__;
dtau_k_c_1_gpu_wM3___(:,:,1+0) = da_k_c_1_gpu_wM__;
dtau_k_c_1_gpu_wM3___(:,:,1+1) = db_k_c_1_gpu_wM__;
dtau_k_c_1_gpu_wM3___(:,:,1+2) = dc_k_c_1_gpu_wM__;
dtau_k_c_2_gpu_wM3___(:,:,1+0) = da_k_c_2_gpu_wM__;
dtau_k_c_2_gpu_wM3___(:,:,1+1) = db_k_c_2_gpu_wM__;
dtau_k_c_2_gpu_wM3___(:,:,1+2) = dc_k_c_2_gpu_wM__;
dtau_k_p_r01_gpu_wM3___(:,:,1+0) = da_k_p_r01_gpu_wM__;
dtau_k_p_r01_gpu_wM3___(:,:,1+1) = db_k_p_r01_gpu_wM__;
dtau_k_p_r01_gpu_wM3___(:,:,1+2) = dc_k_p_r01_gpu_wM__;
dtau_k_p_polar_a_gpu_wM3___(:,:,1+0) = da_k_p_polar_a_gpu_wM__;
dtau_k_p_polar_a_gpu_wM3___(:,:,1+1) = db_k_p_polar_a_gpu_wM__;
dtau_k_p_polar_a_gpu_wM3___(:,:,1+2) = dc_k_p_polar_a_gpu_wM__;
dtau_k_p_azimu_b_gpu_wM3___(:,:,1+0) = da_k_p_azimu_b_gpu_wM__;
dtau_k_p_azimu_b_gpu_wM3___(:,:,1+1) = db_k_p_azimu_b_gpu_wM__;
dtau_k_p_azimu_b_gpu_wM3___(:,:,1+2) = dc_k_p_azimu_b_gpu_wM__;
%%%%;
end;%if flag_d;
%%%%%%%%;

%%%%%%%%;
if flag_dd;
%%%%;
r01 = @(r0,r1) sqrt(r0.^2 + r1.^2) ;
d0_d0_r01 = @(r0,r1) 1./max(tolerance_cg_rhs,r01(r0,r1)) - r0.^2./max(tolerance_cg_rhs,r01(r0,r1)).^3 ;
d1_d0_r01 = @(r0,r1) - (r0.*r1)./max(tolerance_cg_rhs,r01(r0,r1)).^3 ;
d0_d1_r01 = @(r0,r1) - (r1.*r0)./max(tolerance_cg_rhs,r01(r0,r1)).^3 ;
d1_d1_r01 = @(r0,r1) 1./max(tolerance_cg_rhs,r01(r0,r1)) - r1.^2./max(tolerance_cg_rhs,r01(r0,r1)).^3 ;
%dx_dx_z = @(y,x) +2*y./x.^3 ; %<-- Warning: signed divide. ;
%dz_atan2_times_dx_dx_z = @(y,x) +2*y.*x.^1 ./ max(tolerance_cg_rhs^2,x.^2 + y.^2) ./ max(tolerance_cg_rhs^2,x.^2) ; %<-- avoids signed divide. ;
dz_atan2_times_dx_dx_z = @(y,x) +2*y.*x.^1 ./ max(tolerance_cg_rhs^3,(x.^2 + y.^2).*x.^2) ; %<-- avoids signed divide. ;
dy_dx_z = @(y,x) -1./max(tolerance_cg_rhs^2,x.^2) ;
dz_atan2_times_dy_dx_z = @(y,x) -1 ./ max(tolerance_cg_rhs^2,x.^2 + y.^2) ; %<-- avoids signed divide. ;
dx_dy_z = @(y,x) -1./max(tolerance_cg_rhs^2,x.^2) ;
dz_atan2_times_dx_dy_z = @(y,x) -1 ./ max(tolerance_cg_rhs^2,x.^2 + y.^2) ; %<-- avoids signed divide. ;
dy_dy_z = @(y,x) gpuArray( zeros(size(y)) );
dz_atan2_times_dy_dy_z = @(y,x) gpuArray( zeros(size(y)) ) ; %<-- avoids signed divide. ;
dz_dz_atan2 = @(y,x) -2*(x.^3.*y) ./ max(tolerance_cg_rhs^2,x.^2 + y.^2).^2 ;
dz_dz_atan2_times_dy_z = @(y,x) -2*(x.^2.*y.^1) ./ max(tolerance_cg_rhs^2,x.^2 + y.^2).^2 ; %<-- avoids signed divide. ;
dz_dz_atan2_times_dx_z = @(y,x) +2*(x.^1.*y.^2) ./ max(tolerance_cg_rhs^2,x.^2 + y.^2).^2 ; %<-- avoids signed divide. ;
dz_dz_atan2_times_dy_z_times_dy_z = @(y,x) -2*(x.^1.*y.^1) ./ max(tolerance_cg_rhs^2,x.^2 + y.^2).^2 ; %<-- avoids signed divide. ;
dz_dz_atan2_times_dx_z_times_dy_z = @(y,x) +2*(x.^0.*y.^2) ./ max(tolerance_cg_rhs^2,x.^2 + y.^2).^2 ; %<-- avoids signed divide. ;
dz_dz_atan2_times_dy_z_times_dx_z = @(y,x) +2*(x.^0.*y.^2) ./ max(tolerance_cg_rhs^2,x.^2 + y.^2).^2 ; %<-- avoids signed divide. ;
%dz_dz_atan2_times_dx_z_times_dx_z = @(y,x) -2*(x.^1.*y.^3) ./ max(tolerance_cg_rhs^2,x.^2 + y.^2).^2 ./ max(tolerance_cg_rhs^2,x.^2); %<-- avoids signed divide. ;
dz_dz_atan2_times_dx_z_times_dx_z = @(y,x) -2*(x.^1.*y.^3) ./ max(tolerance_cg_rhs^5,(x.^2 + y.^2).^2.*x.^2); %<-- avoids signed divide. ;
%dx_dx_atan2 = @(y,x) dz_dz_atan2(y,x).*dx_z(y,x).*dx_z(y,x) + dz_atan2(y,x).*dx_dx_z(y,x) ; %<-- uses signed divide. ;
dx_dx_atan2 = @(y,x) dz_dz_atan2_times_dx_z_times_dx_z(y,x) + dz_atan2_times_dx_dx_z(y,x) ; %<-- avoids signed divide. ;
%dy_dx_atan2 = @(y,x) dz_dz_atan2(y,x).*dy_z(y,x).*dx_z(y,x) + dz_atan2(y,x).*dy_dx_z(y,x) ; %<-- uses signed divide. ;
dy_dx_atan2 = @(y,x) dz_dz_atan2_times_dy_z_times_dx_z(y,x) + dz_atan2_times_dy_dx_z(y,x) ; %<-- avoids signed divide. ;
%dx_dy_atan2 = @(y,x) dz_dz_atan2(y,x).*dx_z(y,x).*dy_z(y,x) + dz_atan2(y,x).*dx_dy_z(y,x) ; %<-- uses signed divide. ;
dx_dy_atan2 = @(y,x) dz_dz_atan2_times_dx_z_times_dy_z(y,x) + dz_atan2_times_dx_dy_z(y,x) ; %<-- avoids signed divide. ;
%dy_dy_atan2 = @(y,x) dz_dz_atan2(y,x).*dy_z(y,x).*dy_z(y,x) + dz_atan2(y,x).*dy_dy_z(y,x) ; %<-- uses signed divide. ;
dy_dy_atan2 = @(y,x) dz_dz_atan2_times_dy_z_times_dy_z(y,x) + dz_atan2_times_dy_dy_z(y,x) ; %<-- avoids signed divide. ;
%%;
da_da_sa_gpu_wM__ = repmat(reshape(-sin(viewing_polar_a_gpu_M_),[1,n_M]),[n_w,1]);
da_da_ca_gpu_wM__ = repmat(reshape(-cos(viewing_polar_a_gpu_M_),[1,n_M]),[n_w,1]);
db_db_sb_gpu_wM__ = repmat(reshape(-sin(viewing_azimu_b_gpu_M_),[1,n_M]),[n_w,1]);
db_db_cb_gpu_wM__ = repmat(reshape(-cos(viewing_azimu_b_gpu_M_),[1,n_M]),[n_w,1]);
dc_dc_sc_gpu_wM__ = -sin(combine_gamma_z_gpu_wM__); %<-- note sign switch. ;
dc_dc_cc_gpu_wM__ = -cos(combine_gamma_z_gpu_wM__); %<-- note sign switch. ;
da_da_k_c_0_gpu_wM__ = +cb_gpu_wM__.*da_da_ca_gpu_wM__.*cc_gpu_wM__ ;
db_da_k_c_0_gpu_wM__ = +db_cb_gpu_wM__.*da_ca_gpu_wM__.*cc_gpu_wM__ ;
dc_da_k_c_0_gpu_wM__ = +cb_gpu_wM__.*da_ca_gpu_wM__.*dc_cc_gpu_wM__ ;
da_db_k_c_0_gpu_wM__ = +db_cb_gpu_wM__.*da_ca_gpu_wM__.*cc_gpu_wM__ ;
db_db_k_c_0_gpu_wM__ = +db_db_cb_gpu_wM__.*ca_gpu_wM__.*cc_gpu_wM__ - db_db_sb_gpu_wM__.*sc_gpu_wM__ ; 
dc_db_k_c_0_gpu_wM__ = +db_cb_gpu_wM__.*ca_gpu_wM__.*dc_cc_gpu_wM__ - db_sb_gpu_wM__.*dc_sc_gpu_wM__ ; 
da_dc_k_c_0_gpu_wM__ = +cb_gpu_wM__.*da_ca_gpu_wM__.*dc_cc_gpu_wM__ ;
db_dc_k_c_0_gpu_wM__ = +db_cb_gpu_wM__.*ca_gpu_wM__.*dc_cc_gpu_wM__ - db_sb_gpu_wM__.*dc_sc_gpu_wM__ ;
dc_dc_k_c_0_gpu_wM__ = +cb_gpu_wM__.*ca_gpu_wM__.*dc_dc_cc_gpu_wM__ - sb_gpu_wM__.*dc_dc_sc_gpu_wM__ ;
da_da_k_c_1_gpu_wM__ = +sb_gpu_wM__.*da_da_ca_gpu_wM__.*cc_gpu_wM__ ;
db_da_k_c_1_gpu_wM__ = +db_sb_gpu_wM__.*da_ca_gpu_wM__.*cc_gpu_wM__ ;
dc_da_k_c_1_gpu_wM__ = +sb_gpu_wM__.*da_ca_gpu_wM__.*dc_cc_gpu_wM__ ;
da_db_k_c_1_gpu_wM__ = +db_sb_gpu_wM__.*da_ca_gpu_wM__.*cc_gpu_wM__ ;
db_db_k_c_1_gpu_wM__ = +db_db_sb_gpu_wM__.*ca_gpu_wM__.*cc_gpu_wM__ + db_db_cb_gpu_wM__.*sc_gpu_wM__ ; 
dc_db_k_c_1_gpu_wM__ = +db_sb_gpu_wM__.*ca_gpu_wM__.*dc_cc_gpu_wM__ + db_cb_gpu_wM__.*dc_sc_gpu_wM__ ; 
da_dc_k_c_1_gpu_wM__ = +sb_gpu_wM__.*da_ca_gpu_wM__.*dc_cc_gpu_wM__ ;
db_dc_k_c_1_gpu_wM__ = +db_sb_gpu_wM__.*ca_gpu_wM__.*dc_cc_gpu_wM__ + db_cb_gpu_wM__.*dc_sc_gpu_wM__ ; 
dc_dc_k_c_1_gpu_wM__ = +sb_gpu_wM__.*ca_gpu_wM__.*dc_dc_cc_gpu_wM__ + cb_gpu_wM__.*dc_dc_sc_gpu_wM__ ; 
da_da_k_c_2_gpu_wM__ = -da_da_sa_gpu_wM__.*cc_gpu_wM__ ;
db_da_k_c_2_gpu_wM__ = gpuArray( zeros(n_w,n_M) );
dc_da_k_c_2_gpu_wM__ = -da_sa_gpu_wM__.*dc_cc_gpu_wM__ ;
da_db_k_c_2_gpu_wM__ = gpuArray( zeros(n_w,n_M) );
db_db_k_c_2_gpu_wM__ = gpuArray( zeros(n_w,n_M) );
dc_db_k_c_2_gpu_wM__ = gpuArray( zeros(n_w,n_M) );
da_dc_k_c_2_gpu_wM__ = -da_sa_gpu_wM__.*dc_cc_gpu_wM__ ;
db_dc_k_c_2_gpu_wM__ = gpuArray( zeros(n_w,n_M) );
dc_dc_k_c_2_gpu_wM__ = -sa_gpu_wM__.*dc_dc_cc_gpu_wM__ ;
%%;
da_da_k_p_r01_gpu_wM__ = ...
+ d0_d0_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*da_k_c_0_gpu_wM__.*da_k_c_0_gpu_wM__ ...
+ d1_d0_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*da_k_c_1_gpu_wM__.*da_k_c_0_gpu_wM__ ...
+ d0_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*da_da_k_c_0_gpu_wM__ ...
+ d0_d1_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*da_k_c_0_gpu_wM__.*da_k_c_1_gpu_wM__ ...
+ d1_d1_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*da_k_c_1_gpu_wM__.*da_k_c_1_gpu_wM__ ...
+ d1_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*da_da_k_c_1_gpu_wM__ ...
;
db_da_k_p_r01_gpu_wM__ = ...
+ d0_d0_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*db_k_c_0_gpu_wM__.*da_k_c_0_gpu_wM__ ...
+ d1_d0_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*db_k_c_1_gpu_wM__.*da_k_c_0_gpu_wM__ ...
+ d0_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*db_da_k_c_0_gpu_wM__ ...
+ d0_d1_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*db_k_c_0_gpu_wM__.*da_k_c_1_gpu_wM__ ...
+ d1_d1_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*db_k_c_1_gpu_wM__.*da_k_c_1_gpu_wM__ ...
+ d1_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*db_da_k_c_1_gpu_wM__ ...
;
dc_da_k_p_r01_gpu_wM__ = ...
+ d0_d0_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*dc_k_c_0_gpu_wM__.*da_k_c_0_gpu_wM__ ...
+ d1_d0_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*dc_k_c_1_gpu_wM__.*da_k_c_0_gpu_wM__ ...
+ d0_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*dc_da_k_c_0_gpu_wM__ ...
+ d0_d1_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*dc_k_c_0_gpu_wM__.*da_k_c_1_gpu_wM__ ...
+ d1_d1_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*dc_k_c_1_gpu_wM__.*da_k_c_1_gpu_wM__ ...
+ d1_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*dc_da_k_c_1_gpu_wM__ ...
;
%%;
da_db_k_p_r01_gpu_wM__ =  ...
+ d0_d0_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*da_k_c_0_gpu_wM__.*db_k_c_0_gpu_wM__ ...
+ d1_d0_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*da_k_c_1_gpu_wM__.*db_k_c_0_gpu_wM__ ...
+ d0_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*da_db_k_c_0_gpu_wM__ ...
+ d0_d1_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*da_k_c_0_gpu_wM__.*db_k_c_1_gpu_wM__ ...
+ d1_d1_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*da_k_c_1_gpu_wM__.*db_k_c_1_gpu_wM__ ...
+ d1_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*da_db_k_c_1_gpu_wM__ ...
;
db_db_k_p_r01_gpu_wM__ =  ...
+ d0_d0_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*db_k_c_0_gpu_wM__.*db_k_c_0_gpu_wM__ ...
+ d1_d0_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*db_k_c_1_gpu_wM__.*db_k_c_0_gpu_wM__ ...
+ d0_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*db_db_k_c_0_gpu_wM__ ...
+ d0_d1_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*db_k_c_0_gpu_wM__.*db_k_c_1_gpu_wM__ ...
+ d1_d1_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*db_k_c_1_gpu_wM__.*db_k_c_1_gpu_wM__ ...
+ d1_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*db_db_k_c_1_gpu_wM__ ...
;
dc_db_k_p_r01_gpu_wM__ =  ...
+ d0_d0_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*dc_k_c_0_gpu_wM__.*db_k_c_0_gpu_wM__ ...
+ d1_d0_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*dc_k_c_1_gpu_wM__.*db_k_c_0_gpu_wM__ ...
+ d0_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*dc_db_k_c_0_gpu_wM__ ...
+ d0_d1_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*dc_k_c_0_gpu_wM__.*db_k_c_1_gpu_wM__ ...
+ d1_d1_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*dc_k_c_1_gpu_wM__.*db_k_c_1_gpu_wM__ ...
+ d1_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*dc_db_k_c_1_gpu_wM__ ...
;
%%;
da_dc_k_p_r01_gpu_wM__ = ...
+ d0_d0_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*da_k_c_0_gpu_wM__.*dc_k_c_0_gpu_wM__ ...
+ d1_d0_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*da_k_c_1_gpu_wM__.*dc_k_c_0_gpu_wM__ ...
+ d0_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*da_dc_k_c_0_gpu_wM__ ...
+ d0_d1_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*da_k_c_0_gpu_wM__.*dc_k_c_1_gpu_wM__ ...
+ d1_d1_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*da_k_c_1_gpu_wM__.*dc_k_c_1_gpu_wM__ ...
+ d1_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*da_dc_k_c_1_gpu_wM__ ...
;
db_dc_k_p_r01_gpu_wM__ = ...
+ d0_d0_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*db_k_c_0_gpu_wM__.*dc_k_c_0_gpu_wM__ ...
+ d1_d0_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*db_k_c_1_gpu_wM__.*dc_k_c_0_gpu_wM__ ...
+ d0_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*db_dc_k_c_0_gpu_wM__ ...
+ d0_d1_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*db_k_c_0_gpu_wM__.*dc_k_c_1_gpu_wM__ ...
+ d1_d1_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*db_k_c_1_gpu_wM__.*dc_k_c_1_gpu_wM__ ...
+ d1_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*db_dc_k_c_1_gpu_wM__ ...
;
dc_dc_k_p_r01_gpu_wM__ = ...
+ d0_d0_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*dc_k_c_0_gpu_wM__.*dc_k_c_0_gpu_wM__ ...
+ d1_d0_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*dc_k_c_1_gpu_wM__.*dc_k_c_0_gpu_wM__ ...
+ d0_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*dc_dc_k_c_0_gpu_wM__ ...
+ d0_d1_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*dc_k_c_0_gpu_wM__.*dc_k_c_1_gpu_wM__ ...
+ d1_d1_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*dc_k_c_1_gpu_wM__.*dc_k_c_1_gpu_wM__ ...
+ d1_r01(k_c_0_gpu_wM__,k_c_1_gpu_wM__).*dc_dc_k_c_1_gpu_wM__ ...
;
%%;
da_da_k_p_polar_a_gpu_wM__ = ...
  + dy_dy_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*da_k_p_r01_gpu_wM__.*da_k_p_r01_gpu_wM__ ...
  + dx_dy_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*da_k_c_2_gpu_wM__.*da_k_p_r01_gpu_wM__ ...
  + dy_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*da_da_k_p_r01_gpu_wM__ ...
  + dy_dx_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*da_k_p_r01_gpu_wM__.*da_k_c_2_gpu_wM__ ...
  + dx_dx_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*da_k_c_2_gpu_wM__.*da_k_c_2_gpu_wM__ ...
  + dx_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*da_da_k_c_2_gpu_wM__ ...
  ;
db_da_k_p_polar_a_gpu_wM__ = ...
  + dy_dy_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*db_k_p_r01_gpu_wM__.*da_k_p_r01_gpu_wM__ ...
  + dx_dy_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*db_k_c_2_gpu_wM__.*da_k_p_r01_gpu_wM__ ...
  + dy_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*db_da_k_p_r01_gpu_wM__ ...
  + dy_dx_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*db_k_p_r01_gpu_wM__.*da_k_c_2_gpu_wM__ ...
  + dx_dx_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*db_k_c_2_gpu_wM__.*da_k_c_2_gpu_wM__ ...
  + dx_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*db_da_k_c_2_gpu_wM__ ...
  ;
dc_da_k_p_polar_a_gpu_wM__ = ...
  + dy_dy_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*dc_k_p_r01_gpu_wM__.*da_k_p_r01_gpu_wM__ ...
  + dx_dy_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*dc_k_c_2_gpu_wM__.*da_k_p_r01_gpu_wM__ ...
  + dy_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*dc_da_k_p_r01_gpu_wM__ ...
  + dy_dx_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*dc_k_p_r01_gpu_wM__.*da_k_c_2_gpu_wM__ ...
  + dx_dx_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*dc_k_c_2_gpu_wM__.*da_k_c_2_gpu_wM__ ...
  + dx_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*dc_da_k_c_2_gpu_wM__ ...
  ;
da_db_k_p_polar_a_gpu_wM__ = ...
  + dy_dy_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*da_k_p_r01_gpu_wM__.*db_k_p_r01_gpu_wM__ ...
  + dx_dy_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*da_k_c_2_gpu_wM__.*db_k_p_r01_gpu_wM__ ...
  + dy_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*da_db_k_p_r01_gpu_wM__ ...
  + dy_dx_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*da_k_p_r01_gpu_wM__.*db_k_c_2_gpu_wM__ ...
  + dx_dx_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*da_k_c_2_gpu_wM__.*db_k_c_2_gpu_wM__ ...
  + dx_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*da_db_k_c_2_gpu_wM__ ...
  ;
db_db_k_p_polar_a_gpu_wM__ = ...
  + dy_dy_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*db_k_p_r01_gpu_wM__.*db_k_p_r01_gpu_wM__ ...
  + dx_dy_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*db_k_c_2_gpu_wM__.*db_k_p_r01_gpu_wM__ ...
  + dy_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*db_db_k_p_r01_gpu_wM__ ...
  + dy_dx_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*db_k_p_r01_gpu_wM__.*db_k_c_2_gpu_wM__ ...
  + dx_dx_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*db_k_c_2_gpu_wM__.*db_k_c_2_gpu_wM__ ...
  + dx_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*db_db_k_c_2_gpu_wM__ ...
  ;
dc_db_k_p_polar_a_gpu_wM__ = ...
  + dy_dy_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*dc_k_p_r01_gpu_wM__.*db_k_p_r01_gpu_wM__ ...
  + dx_dy_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*dc_k_c_2_gpu_wM__.*db_k_p_r01_gpu_wM__ ...
  + dy_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*dc_db_k_p_r01_gpu_wM__ ...
  + dy_dx_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*dc_k_p_r01_gpu_wM__.*db_k_c_2_gpu_wM__ ...
  + dx_dx_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*dc_k_c_2_gpu_wM__.*db_k_c_2_gpu_wM__ ...
  + dx_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*dc_db_k_c_2_gpu_wM__ ...
  ;
da_dc_k_p_polar_a_gpu_wM__ = ...
  + dy_dy_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*da_k_p_r01_gpu_wM__.*dc_k_p_r01_gpu_wM__ ...
  + dx_dy_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*da_k_c_2_gpu_wM__.*dc_k_p_r01_gpu_wM__ ...
  + dy_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*da_dc_k_p_r01_gpu_wM__ ...
  + dy_dx_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*da_k_p_r01_gpu_wM__.*dc_k_c_2_gpu_wM__ ...
  + dx_dx_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*da_k_c_2_gpu_wM__.*dc_k_c_2_gpu_wM__ ...
  + dx_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*da_dc_k_c_2_gpu_wM__ ...
  ;
db_dc_k_p_polar_a_gpu_wM__ = ...
  + dy_dy_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*db_k_p_r01_gpu_wM__.*dc_k_p_r01_gpu_wM__ ...
  + dx_dy_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*db_k_c_2_gpu_wM__.*dc_k_p_r01_gpu_wM__ ...
  + dy_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*db_dc_k_p_r01_gpu_wM__ ...
  + dy_dx_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*db_k_p_r01_gpu_wM__.*dc_k_c_2_gpu_wM__ ...
  + dx_dx_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*db_k_c_2_gpu_wM__.*dc_k_c_2_gpu_wM__ ...
  + dx_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*db_dc_k_c_2_gpu_wM__ ...
  ;
dc_dc_k_p_polar_a_gpu_wM__ = ...
  + dy_dy_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*dc_k_p_r01_gpu_wM__.*dc_k_p_r01_gpu_wM__ ...
  + dx_dy_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*dc_k_c_2_gpu_wM__.*dc_k_p_r01_gpu_wM__ ...
  + dy_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*dc_dc_k_p_r01_gpu_wM__ ...
  + dy_dx_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*dc_k_p_r01_gpu_wM__.*dc_k_c_2_gpu_wM__ ...
  + dx_dx_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*dc_k_c_2_gpu_wM__.*dc_k_c_2_gpu_wM__ ...
  + dx_atan2(k_p_r01_gpu_wM__,k_c_2_gpu_wM__).*dc_dc_k_c_2_gpu_wM__ ...
  ;
%%;
da_da_k_p_azimu_b_gpu_wM__ = ...
  + dy_dy_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*da_k_c_1_gpu_wM__.*da_k_c_1_gpu_wM__ ...
  + dx_dy_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*da_k_c_0_gpu_wM__.*da_k_c_1_gpu_wM__ ...
  + dy_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*da_da_k_c_1_gpu_wM__ ...
  + dy_dx_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*da_k_c_1_gpu_wM__.*da_k_c_0_gpu_wM__ ...
  + dx_dx_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*da_k_c_0_gpu_wM__.*da_k_c_0_gpu_wM__ ...
  + dx_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*da_da_k_c_0_gpu_wM__ ...
  ;
db_da_k_p_azimu_b_gpu_wM__ = ...
  + dy_dy_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*db_k_c_1_gpu_wM__.*da_k_c_1_gpu_wM__ ...
  + dx_dy_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*db_k_c_0_gpu_wM__.*da_k_c_1_gpu_wM__ ...
  + dy_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*db_da_k_c_1_gpu_wM__ ...
  + dy_dx_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*db_k_c_1_gpu_wM__.*da_k_c_0_gpu_wM__ ...
  + dx_dx_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*db_k_c_0_gpu_wM__.*da_k_c_0_gpu_wM__ ...
  + dx_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*db_da_k_c_0_gpu_wM__ ...
  ;
dc_da_k_p_azimu_b_gpu_wM__ = ...
  + dy_dy_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*dc_k_c_1_gpu_wM__.*da_k_c_1_gpu_wM__ ...
  + dx_dy_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*dc_k_c_0_gpu_wM__.*da_k_c_1_gpu_wM__ ...
  + dy_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*dc_da_k_c_1_gpu_wM__ ...
  + dy_dx_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*dc_k_c_1_gpu_wM__.*da_k_c_0_gpu_wM__ ...
  + dx_dx_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*dc_k_c_0_gpu_wM__.*da_k_c_0_gpu_wM__ ...
  + dx_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*dc_da_k_c_0_gpu_wM__ ...
  ;
da_db_k_p_azimu_b_gpu_wM__ = ...
  + dy_dy_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*da_k_c_1_gpu_wM__.*db_k_c_1_gpu_wM__ ...
  + dx_dy_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*da_k_c_0_gpu_wM__.*db_k_c_1_gpu_wM__ ...
  + dy_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*da_db_k_c_1_gpu_wM__ ...
  + dy_dx_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*da_k_c_1_gpu_wM__.*db_k_c_0_gpu_wM__ ...
  + dx_dx_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*da_k_c_0_gpu_wM__.*db_k_c_0_gpu_wM__ ...
  + dx_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*da_db_k_c_0_gpu_wM__ ...
  ;
db_db_k_p_azimu_b_gpu_wM__ = ...
  + dy_dy_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*db_k_c_1_gpu_wM__.*db_k_c_1_gpu_wM__ ...
  + dx_dy_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*db_k_c_0_gpu_wM__.*db_k_c_1_gpu_wM__ ...
  + dy_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*db_db_k_c_1_gpu_wM__ ...
  + dy_dx_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*db_k_c_1_gpu_wM__.*db_k_c_0_gpu_wM__ ...
  + dx_dx_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*db_k_c_0_gpu_wM__.*db_k_c_0_gpu_wM__ ...
  + dx_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*db_db_k_c_0_gpu_wM__ ...
  ;
dc_db_k_p_azimu_b_gpu_wM__ = ...
  + dy_dy_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*dc_k_c_1_gpu_wM__.*db_k_c_1_gpu_wM__ ...
  + dx_dy_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*dc_k_c_0_gpu_wM__.*db_k_c_1_gpu_wM__ ...
  + dy_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*dc_db_k_c_1_gpu_wM__ ...
  + dy_dx_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*dc_k_c_1_gpu_wM__.*db_k_c_0_gpu_wM__ ...
  + dx_dx_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*dc_k_c_0_gpu_wM__.*db_k_c_0_gpu_wM__ ...
  + dx_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*dc_db_k_c_0_gpu_wM__ ...
  ;
da_dc_k_p_azimu_b_gpu_wM__ = ...
  + dy_dy_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*da_k_c_1_gpu_wM__.*dc_k_c_1_gpu_wM__ ...
  + dx_dy_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*da_k_c_0_gpu_wM__.*dc_k_c_1_gpu_wM__ ...
  + dy_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*da_dc_k_c_1_gpu_wM__ ...
  + dy_dx_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*da_k_c_1_gpu_wM__.*dc_k_c_0_gpu_wM__ ...
  + dx_dx_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*da_k_c_0_gpu_wM__.*dc_k_c_0_gpu_wM__ ...
  + dx_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*da_dc_k_c_0_gpu_wM__ ...
  ;
db_dc_k_p_azimu_b_gpu_wM__ = ...
  + dy_dy_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*db_k_c_1_gpu_wM__.*dc_k_c_1_gpu_wM__ ...
  + dx_dy_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*db_k_c_0_gpu_wM__.*dc_k_c_1_gpu_wM__ ...
  + dy_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*db_dc_k_c_1_gpu_wM__ ...
  + dy_dx_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*db_k_c_1_gpu_wM__.*dc_k_c_0_gpu_wM__ ...
  + dx_dx_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*db_k_c_0_gpu_wM__.*dc_k_c_0_gpu_wM__ ...
  + dx_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*db_dc_k_c_0_gpu_wM__ ...
  ;
dc_dc_k_p_azimu_b_gpu_wM__ = ...
  + dy_dy_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*dc_k_c_1_gpu_wM__.*dc_k_c_1_gpu_wM__ ...
  + dx_dy_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*dc_k_c_0_gpu_wM__.*dc_k_c_1_gpu_wM__ ...
  + dy_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*dc_dc_k_c_1_gpu_wM__ ...
  + dy_dx_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*dc_k_c_1_gpu_wM__.*dc_k_c_0_gpu_wM__ ...
  + dx_dx_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*dc_k_c_0_gpu_wM__.*dc_k_c_0_gpu_wM__ ...
  + dx_atan2(k_c_1_gpu_wM__,k_c_0_gpu_wM__).*dc_dc_k_c_0_gpu_wM__ ...
  ;
%%;
% cleanup. ;
%%;
if flag_cleanup;
%%;
for ntype=0:9-1;
na=0;
if ntype==na; tmp_gpu_wM__ = da_da_k_p_r01_gpu_wM__; tmp_str_ = 'da_da_k_p_r01_gpu_wM__'; end; na=na+1;
if ntype==na; tmp_gpu_wM__ = da_db_k_p_r01_gpu_wM__; tmp_str_ = 'da_db_k_p_r01_gpu_wM__'; end; na=na+1;
if ntype==na; tmp_gpu_wM__ = da_dc_k_p_r01_gpu_wM__; tmp_str_ = 'da_dc_k_p_r01_gpu_wM__'; end; na=na+1;
if ntype==na; tmp_gpu_wM__ = db_da_k_p_r01_gpu_wM__; tmp_str_ = 'db_da_k_p_r01_gpu_wM__'; end; na=na+1;
if ntype==na; tmp_gpu_wM__ = db_db_k_p_r01_gpu_wM__; tmp_str_ = 'db_db_k_p_r01_gpu_wM__'; end; na=na+1;
if ntype==na; tmp_gpu_wM__ = db_dc_k_p_r01_gpu_wM__; tmp_str_ = 'db_dc_k_p_r01_gpu_wM__'; end; na=na+1;
if ntype==na; tmp_gpu_wM__ = dc_da_k_p_r01_gpu_wM__; tmp_str_ = 'dc_da_k_p_r01_gpu_wM__'; end; na=na+1;
if ntype==na; tmp_gpu_wM__ = dc_db_k_p_r01_gpu_wM__; tmp_str_ = 'dc_db_k_p_r01_gpu_wM__'; end; na=na+1;
if ntype==na; tmp_gpu_wM__ = dc_dc_k_p_r01_gpu_wM__; tmp_str_ = 'dc_dc_k_p_r01_gpu_wM__'; end; na=na+1;
tmp_index_ = efind(abs(tmp_gpu_wM__)> tolerance_cg_rhs_upb); tmp_n_index = numel(tmp_index_);
if tmp_n_index> 0;
if (flag_verbose>1); disp(sprintf(' %% cleanup %s: tmp_n_index %d',tmp_str_,tmp_n_index)); end;
tmp_clean_gpu_wM__ = 0.5*(circshift(tmp_gpu_wM__,+1,1)+circshift(tmp_gpu_wM__,-1,1));
tmp_gpu_wM__(1+tmp_index_) = tmp_clean_gpu_wM__(1+tmp_index_);
end;%if tmp_n_index> 0;
na=0;
if ntype==na; da_da_k_p_r01_gpu_wM__ = tmp_gpu_wM__; end; na=na+1;
if ntype==na; da_db_k_p_r01_gpu_wM__ = tmp_gpu_wM__; end; na=na+1;
if ntype==na; da_dc_k_p_r01_gpu_wM__ = tmp_gpu_wM__; end; na=na+1;
if ntype==na; db_da_k_p_r01_gpu_wM__ = tmp_gpu_wM__; end; na=na+1;
if ntype==na; db_db_k_p_r01_gpu_wM__ = tmp_gpu_wM__; end; na=na+1;
if ntype==na; db_dc_k_p_r01_gpu_wM__ = tmp_gpu_wM__; end; na=na+1;
if ntype==na; dc_da_k_p_r01_gpu_wM__ = tmp_gpu_wM__; end; na=na+1;
if ntype==na; dc_db_k_p_r01_gpu_wM__ = tmp_gpu_wM__; end; na=na+1;
if ntype==na; dc_dc_k_p_r01_gpu_wM__ = tmp_gpu_wM__; end; na=na+1;
end;%for ntype=0:9-1;
%%;
for ntype=0:9-1;
na=0;
if ntype==na; tmp_gpu_wM__ = da_da_k_p_polar_a_gpu_wM__; tmp_str_ = 'da_da_k_p_polar_a_gpu_wM__'; end; na=na+1;
if ntype==na; tmp_gpu_wM__ = da_db_k_p_polar_a_gpu_wM__; tmp_str_ = 'da_db_k_p_polar_a_gpu_wM__'; end; na=na+1;
if ntype==na; tmp_gpu_wM__ = da_dc_k_p_polar_a_gpu_wM__; tmp_str_ = 'da_dc_k_p_polar_a_gpu_wM__'; end; na=na+1;
if ntype==na; tmp_gpu_wM__ = db_da_k_p_polar_a_gpu_wM__; tmp_str_ = 'db_da_k_p_polar_a_gpu_wM__'; end; na=na+1;
if ntype==na; tmp_gpu_wM__ = db_db_k_p_polar_a_gpu_wM__; tmp_str_ = 'db_db_k_p_polar_a_gpu_wM__'; end; na=na+1;
if ntype==na; tmp_gpu_wM__ = db_dc_k_p_polar_a_gpu_wM__; tmp_str_ = 'db_dc_k_p_polar_a_gpu_wM__'; end; na=na+1;
if ntype==na; tmp_gpu_wM__ = dc_da_k_p_polar_a_gpu_wM__; tmp_str_ = 'dc_da_k_p_polar_a_gpu_wM__'; end; na=na+1;
if ntype==na; tmp_gpu_wM__ = dc_db_k_p_polar_a_gpu_wM__; tmp_str_ = 'dc_db_k_p_polar_a_gpu_wM__'; end; na=na+1;
if ntype==na; tmp_gpu_wM__ = dc_dc_k_p_polar_a_gpu_wM__; tmp_str_ = 'dc_dc_k_p_polar_a_gpu_wM__'; end; na=na+1;
tmp_index_ = efind(abs(tmp_gpu_wM__)> tolerance_cg_rhs_upb); tmp_n_index = numel(tmp_index_);
if tmp_n_index> 0;
if (flag_verbose>1); disp(sprintf(' %% cleanup %s: tmp_n_index %d',tmp_str_,tmp_n_index)); end;
tmp_clean_gpu_wM__ = 0.5*(circshift(tmp_gpu_wM__,+1,1)+circshift(tmp_gpu_wM__,-1,1));
tmp_gpu_wM__(1+tmp_index_) = tmp_clean_gpu_wM__(1+tmp_index_);
end;%if tmp_n_index> 0;
na=0;
if ntype==na; da_da_k_p_polar_a_gpu_wM__ = tmp_gpu_wM__; end; na=na+1;
if ntype==na; da_db_k_p_polar_a_gpu_wM__ = tmp_gpu_wM__; end; na=na+1;
if ntype==na; da_dc_k_p_polar_a_gpu_wM__ = tmp_gpu_wM__; end; na=na+1;
if ntype==na; db_da_k_p_polar_a_gpu_wM__ = tmp_gpu_wM__; end; na=na+1;
if ntype==na; db_db_k_p_polar_a_gpu_wM__ = tmp_gpu_wM__; end; na=na+1;
if ntype==na; db_dc_k_p_polar_a_gpu_wM__ = tmp_gpu_wM__; end; na=na+1;
if ntype==na; dc_da_k_p_polar_a_gpu_wM__ = tmp_gpu_wM__; end; na=na+1;
if ntype==na; dc_db_k_p_polar_a_gpu_wM__ = tmp_gpu_wM__; end; na=na+1;
if ntype==na; dc_dc_k_p_polar_a_gpu_wM__ = tmp_gpu_wM__; end; na=na+1;
end;%for ntype=0:9-1;
%%;
for ntype=0:9-1;
na=0;
if ntype==na; tmp_gpu_wM__ = da_da_k_p_azimu_b_gpu_wM__; tmp_str_ = 'da_da_k_p_azimu_b_gpu_wM__'; end; na=na+1;
if ntype==na; tmp_gpu_wM__ = da_db_k_p_azimu_b_gpu_wM__; tmp_str_ = 'da_db_k_p_azimu_b_gpu_wM__'; end; na=na+1;
if ntype==na; tmp_gpu_wM__ = da_dc_k_p_azimu_b_gpu_wM__; tmp_str_ = 'da_dc_k_p_azimu_b_gpu_wM__'; end; na=na+1;
if ntype==na; tmp_gpu_wM__ = db_da_k_p_azimu_b_gpu_wM__; tmp_str_ = 'db_da_k_p_azimu_b_gpu_wM__'; end; na=na+1;
if ntype==na; tmp_gpu_wM__ = db_db_k_p_azimu_b_gpu_wM__; tmp_str_ = 'db_db_k_p_azimu_b_gpu_wM__'; end; na=na+1;
if ntype==na; tmp_gpu_wM__ = db_dc_k_p_azimu_b_gpu_wM__; tmp_str_ = 'db_dc_k_p_azimu_b_gpu_wM__'; end; na=na+1;
if ntype==na; tmp_gpu_wM__ = dc_da_k_p_azimu_b_gpu_wM__; tmp_str_ = 'dc_da_k_p_azimu_b_gpu_wM__'; end; na=na+1;
if ntype==na; tmp_gpu_wM__ = dc_db_k_p_azimu_b_gpu_wM__; tmp_str_ = 'dc_db_k_p_azimu_b_gpu_wM__'; end; na=na+1;
if ntype==na; tmp_gpu_wM__ = dc_dc_k_p_azimu_b_gpu_wM__; tmp_str_ = 'dc_dc_k_p_azimu_b_gpu_wM__'; end; na=na+1;
tmp_index_ = efind(abs(tmp_gpu_wM__)> tolerance_cg_rhs_upb); tmp_n_index = numel(tmp_index_);
if tmp_n_index> 0;
if (flag_verbose>1); disp(sprintf(' %% cleanup %s: tmp_n_index %d',tmp_str_,tmp_n_index)); end;
tmp_clean_gpu_wM__ = 0.5*(circshift(tmp_gpu_wM__,+1,1)+circshift(tmp_gpu_wM__,-1,1));
tmp_gpu_wM__(1+tmp_index_) = tmp_clean_gpu_wM__(1+tmp_index_);
end;%if tmp_n_index> 0;
na=0;
if ntype==na; da_da_k_p_azimu_b_gpu_wM__ = tmp_gpu_wM__; end; na=na+1;
if ntype==na; da_db_k_p_azimu_b_gpu_wM__ = tmp_gpu_wM__; end; na=na+1;
if ntype==na; da_dc_k_p_azimu_b_gpu_wM__ = tmp_gpu_wM__; end; na=na+1;
if ntype==na; db_da_k_p_azimu_b_gpu_wM__ = tmp_gpu_wM__; end; na=na+1;
if ntype==na; db_db_k_p_azimu_b_gpu_wM__ = tmp_gpu_wM__; end; na=na+1;
if ntype==na; db_dc_k_p_azimu_b_gpu_wM__ = tmp_gpu_wM__; end; na=na+1;
if ntype==na; dc_da_k_p_azimu_b_gpu_wM__ = tmp_gpu_wM__; end; na=na+1;
if ntype==na; dc_db_k_p_azimu_b_gpu_wM__ = tmp_gpu_wM__; end; na=na+1;
if ntype==na; dc_dc_k_p_azimu_b_gpu_wM__ = tmp_gpu_wM__; end; na=na+1;
end;%for ntype=0:9-1;
%%;
end;%if flag_cleanup;
%%;
dtau_dtau_k_c_0_gpu_wM33____(:,:,1+0,1+0) = da_da_k_c_0_gpu_wM__;
dtau_dtau_k_c_0_gpu_wM33____(:,:,1+0,1+1) = db_da_k_c_0_gpu_wM__;
dtau_dtau_k_c_0_gpu_wM33____(:,:,1+0,1+2) = dc_da_k_c_0_gpu_wM__;
dtau_dtau_k_c_0_gpu_wM33____(:,:,1+1,1+0) = da_db_k_c_0_gpu_wM__;
dtau_dtau_k_c_0_gpu_wM33____(:,:,1+1,1+1) = db_db_k_c_0_gpu_wM__;
dtau_dtau_k_c_0_gpu_wM33____(:,:,1+1,1+2) = dc_db_k_c_0_gpu_wM__;
dtau_dtau_k_c_0_gpu_wM33____(:,:,1+2,1+0) = da_dc_k_c_0_gpu_wM__;
dtau_dtau_k_c_0_gpu_wM33____(:,:,1+2,1+1) = db_dc_k_c_0_gpu_wM__;
dtau_dtau_k_c_0_gpu_wM33____(:,:,1+2,1+2) = dc_dc_k_c_0_gpu_wM__;
dtau_dtau_k_c_1_gpu_wM33____(:,:,1+0,1+0) = da_da_k_c_1_gpu_wM__;
dtau_dtau_k_c_1_gpu_wM33____(:,:,1+0,1+1) = db_da_k_c_1_gpu_wM__;
dtau_dtau_k_c_1_gpu_wM33____(:,:,1+0,1+2) = dc_da_k_c_1_gpu_wM__;
dtau_dtau_k_c_1_gpu_wM33____(:,:,1+1,1+0) = da_db_k_c_1_gpu_wM__;
dtau_dtau_k_c_1_gpu_wM33____(:,:,1+1,1+1) = db_db_k_c_1_gpu_wM__;
dtau_dtau_k_c_1_gpu_wM33____(:,:,1+1,1+2) = dc_db_k_c_1_gpu_wM__;
dtau_dtau_k_c_1_gpu_wM33____(:,:,1+2,1+0) = da_dc_k_c_1_gpu_wM__;
dtau_dtau_k_c_1_gpu_wM33____(:,:,1+2,1+1) = db_dc_k_c_1_gpu_wM__;
dtau_dtau_k_c_1_gpu_wM33____(:,:,1+2,1+2) = dc_dc_k_c_1_gpu_wM__;
dtau_dtau_k_c_2_gpu_wM33____(:,:,1+0,1+0) = da_da_k_c_2_gpu_wM__;
dtau_dtau_k_c_2_gpu_wM33____(:,:,1+0,1+1) = db_da_k_c_2_gpu_wM__;
dtau_dtau_k_c_2_gpu_wM33____(:,:,1+0,1+2) = dc_da_k_c_2_gpu_wM__;
dtau_dtau_k_c_2_gpu_wM33____(:,:,1+1,1+0) = da_db_k_c_2_gpu_wM__;
dtau_dtau_k_c_2_gpu_wM33____(:,:,1+1,1+1) = db_db_k_c_2_gpu_wM__;
dtau_dtau_k_c_2_gpu_wM33____(:,:,1+1,1+2) = dc_db_k_c_2_gpu_wM__;
dtau_dtau_k_c_2_gpu_wM33____(:,:,1+2,1+0) = da_dc_k_c_2_gpu_wM__;
dtau_dtau_k_c_2_gpu_wM33____(:,:,1+2,1+1) = db_dc_k_c_2_gpu_wM__;
dtau_dtau_k_c_2_gpu_wM33____(:,:,1+2,1+2) = dc_dc_k_c_2_gpu_wM__;
%%;
dtau_dtau_k_p_r01_gpu_wM33____(:,:,1+0,1+0) = da_da_k_p_r01_gpu_wM__;
dtau_dtau_k_p_r01_gpu_wM33____(:,:,1+0,1+1) = db_da_k_p_r01_gpu_wM__;
dtau_dtau_k_p_r01_gpu_wM33____(:,:,1+0,1+2) = dc_da_k_p_r01_gpu_wM__;
dtau_dtau_k_p_r01_gpu_wM33____(:,:,1+1,1+0) = da_db_k_p_r01_gpu_wM__;
dtau_dtau_k_p_r01_gpu_wM33____(:,:,1+1,1+1) = db_db_k_p_r01_gpu_wM__;
dtau_dtau_k_p_r01_gpu_wM33____(:,:,1+1,1+2) = dc_db_k_p_r01_gpu_wM__;
dtau_dtau_k_p_r01_gpu_wM33____(:,:,1+2,1+0) = da_dc_k_p_r01_gpu_wM__;
dtau_dtau_k_p_r01_gpu_wM33____(:,:,1+2,1+1) = db_dc_k_p_r01_gpu_wM__;
dtau_dtau_k_p_r01_gpu_wM33____(:,:,1+2,1+2) = dc_dc_k_p_r01_gpu_wM__;
%%;
dtau_dtau_k_p_polar_a_gpu_wM33____(:,:,1+0,1+0) = da_da_k_p_polar_a_gpu_wM__;
dtau_dtau_k_p_polar_a_gpu_wM33____(:,:,1+0,1+1) = db_da_k_p_polar_a_gpu_wM__;
dtau_dtau_k_p_polar_a_gpu_wM33____(:,:,1+0,1+2) = dc_da_k_p_polar_a_gpu_wM__;
dtau_dtau_k_p_polar_a_gpu_wM33____(:,:,1+1,1+0) = da_db_k_p_polar_a_gpu_wM__;
dtau_dtau_k_p_polar_a_gpu_wM33____(:,:,1+1,1+1) = db_db_k_p_polar_a_gpu_wM__;
dtau_dtau_k_p_polar_a_gpu_wM33____(:,:,1+1,1+2) = dc_db_k_p_polar_a_gpu_wM__;
dtau_dtau_k_p_polar_a_gpu_wM33____(:,:,1+2,1+0) = da_dc_k_p_polar_a_gpu_wM__;
dtau_dtau_k_p_polar_a_gpu_wM33____(:,:,1+2,1+1) = db_dc_k_p_polar_a_gpu_wM__;
dtau_dtau_k_p_polar_a_gpu_wM33____(:,:,1+2,1+2) = dc_dc_k_p_polar_a_gpu_wM__;
dtau_dtau_k_p_azimu_b_gpu_wM33____(:,:,1+0,1+0) = da_da_k_p_azimu_b_gpu_wM__;
dtau_dtau_k_p_azimu_b_gpu_wM33____(:,:,1+0,1+1) = db_da_k_p_azimu_b_gpu_wM__;
dtau_dtau_k_p_azimu_b_gpu_wM33____(:,:,1+0,1+2) = dc_da_k_p_azimu_b_gpu_wM__;
dtau_dtau_k_p_azimu_b_gpu_wM33____(:,:,1+1,1+0) = da_db_k_p_azimu_b_gpu_wM__;
dtau_dtau_k_p_azimu_b_gpu_wM33____(:,:,1+1,1+1) = db_db_k_p_azimu_b_gpu_wM__;
dtau_dtau_k_p_azimu_b_gpu_wM33____(:,:,1+1,1+2) = dc_db_k_p_azimu_b_gpu_wM__;
dtau_dtau_k_p_azimu_b_gpu_wM33____(:,:,1+2,1+0) = da_dc_k_p_azimu_b_gpu_wM__;
dtau_dtau_k_p_azimu_b_gpu_wM33____(:,:,1+2,1+1) = db_dc_k_p_azimu_b_gpu_wM__;
dtau_dtau_k_p_azimu_b_gpu_wM33____(:,:,1+2,1+2) = dc_dc_k_p_azimu_b_gpu_wM__;
%%%%;
end;%if flag_dd;
%%%%%%%%;

%%%%%%%%;
if flag_1;
if sum(~isfinite(k_p_polar_a_gpu_wM__),'all')> 0; disp(sprintf(' %% Warning, k_p_polar_a_gpu_wM__ not finite in %s',str_thisfunction)); end;
if sum(~isfinite(k_p_azimu_b_gpu_wM__),'all')> 0; disp(sprintf(' %% Warning, k_p_azimu_b_gpu_wM__ not finite in %s',str_thisfunction)); end;
if sum(~isfinite(k_c_0_gpu_wM__),'all')> 0; disp(sprintf(' %% Warning, k_c_0_gpu_wM__ not finite in %s',str_thisfunction)); end;
if sum(~isfinite(k_c_1_gpu_wM__),'all')> 0; disp(sprintf(' %% Warning, k_c_1_gpu_wM__ not finite in %s',str_thisfunction)); end;
if sum(~isfinite(k_c_2_gpu_wM__),'all')> 0; disp(sprintf(' %% Warning, k_c_2_gpu_wM__ not finite in %s',str_thisfunction)); end;
if sum(~isfinite(k_p_r01_gpu_wM__),'all')> 0; disp(sprintf(' %% Warning, k_p_r01_gpu_wM__ not finite in %s',str_thisfunction)); end;
end;%if flag_1;
if flag_d;
if sum(~isfinite(dtau_k_p_polar_a_gpu_wM3___),'all')> 0; disp(sprintf(' %% Warning, dtau_k_p_polar_a_gpu_wM3___ not finite in %s',str_thisfunction)); end;
if sum(~isfinite(dtau_k_p_azimu_b_gpu_wM3___),'all')> 0; disp(sprintf(' %% Warning, dtau_k_p_azimu_b_gpu_wM3___ not finite in %s',str_thisfunction)); end;
if sum(~isfinite(dtau_k_c_0_gpu_wM3___),'all')> 0; disp(sprintf(' %% Warning, dtau_k_c_0_gpu_wM3___ not finite in %s',str_thisfunction)); end;
if sum(~isfinite(dtau_k_c_1_gpu_wM3___),'all')> 0; disp(sprintf(' %% Warning, dtau_k_c_1_gpu_wM3___ not finite in %s',str_thisfunction)); end;
if sum(~isfinite(dtau_k_c_2_gpu_wM3___),'all')> 0; disp(sprintf(' %% Warning, dtau_k_c_2_gpu_wM3___ not finite in %s',str_thisfunction)); end;
if sum(~isfinite(dtau_k_p_r01_gpu_wM3___),'all')> 0; disp(sprintf(' %% Warning, dtau_k_p_r01_gpu_wM3___ not finite in %s',str_thisfunction)); end;
end;%if flag_d;
if flag_dd;
if sum(~isfinite(dtau_dtau_k_p_polar_a_gpu_wM33____),'all')> 0; disp(sprintf(' %% Warning, dtau_dtau_k_p_polar_a_gpu_wM33____ not finite in %s',str_thisfunction)); end;
if sum(~isfinite(dtau_dtau_k_p_azimu_b_gpu_wM33____),'all')> 0; disp(sprintf(' %% Warning, dtau_dtau_k_p_azimu_b_gpu_wM33____ not finite in %s',str_thisfunction)); end;
if sum(~isfinite(dtau_dtau_k_c_0_gpu_wM33____),'all')> 0; disp(sprintf(' %% Warning, dtau_dtau_k_c_0_gpu_wM33____ not finite in %s',str_thisfunction)); end;
if sum(~isfinite(dtau_dtau_k_c_1_gpu_wM33____),'all')> 0; disp(sprintf(' %% Warning, dtau_dtau_k_c_1_gpu_wM33____ not finite in %s',str_thisfunction)); end;
if sum(~isfinite(dtau_dtau_k_c_2_gpu_wM33____),'all')> 0; disp(sprintf(' %% Warning, dtau_dtau_k_c_2_gpu_wM33____ not finite in %s',str_thisfunction)); end;
if sum(~isfinite(dtau_dtau_k_p_r01_gpu_wM33____),'all')> 0; disp(sprintf(' %% Warning, dtau_dtau_k_p_r01_gpu_wM33____ not finite in %s',str_thisfunction)); end;
end;%if flag_dd;
%%%%%%%%;

if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

function ...
disp_error(flag_verbose,str_0_,d_0_,str_1_,d_1_);
if (flag_verbose>0);
f_numerator = fnorm(d_0_-d_1_);
f_denomator = min(fnorm(d_0_),fnorm(d_1_));
r_error = f_numerator/max(1e-12,f_denomator);
n_header = 64;
str_prefix = sprintf(' %% %s vs %s: ',str_0_,str_1_);
str_header = strcat(str_prefix,' '*ones(1,n_header-numel(str_prefix)));
disp(sprintf('%snumerator %16.8f denomator %16.8f (r_error %16.8f)',str_header,f_numerator,f_denomator,r_error));
end;%if (flag_verbose>0);
