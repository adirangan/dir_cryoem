function imagesc_S_k_p_3d_2( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,weight_2d_k_all_ ...
,n_S ...
,S_k_p_wkS__ ...
,template_viewing_azimu_b_S_ ...
,template_viewing_polar_a_S_ ...
);

str_thisfunction = 'imagesc_S_k_p_3d_2';

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); k_p_r_=[]; end; na=na+1;
if (nargin<1+na); k_p_r_max=[]; end; na=na+1;
if (nargin<1+na); n_w_=[]; end; na=na+1;
if (nargin<1+na); weight_2d_k_all_=[]; end; na=na+1;
if (nargin<1+na); n_S=[]; end; na=na+1;
if (nargin<1+na); S_k_p_wkS__=[]; end; na=na+1;
if (nargin<1+na); template_viewing_azimu_b_S_=[]; end; na=na+1;
if (nargin<1+na); template_viewing_polar_a_S_=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;
if ~isfield(parameter,'tolerance_master'); parameter.tolerance_master = 1e-2; end;
tolerance_master = parameter.tolerance_master;
if ~isfield(parameter,'tolerance_stack'); parameter.tolerance_stack = 1e-3; end;
tolerance_stack = parameter.tolerance_stack;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
flag_verbose = parameter.flag_verbose;
if ~isfield(parameter,'flag_percent_use'); parameter.flag_percent_use = 1; end;
flag_percent_use = parameter.flag_percent_use;
if ~isfield(parameter,'percent_threshold'); parameter.percent_threshold = 98.5; end;
percent_threshold = parameter.percent_threshold;
if ~isfield(parameter,'n_contour'); parameter.n_contour = 16; end;
n_contour = parameter.n_contour;
if ~isfield(parameter,'vlim_'); parameter.vlim_ = []; end;
vlim_ = parameter.vlim_;
if ~isfield(parameter,'vval_'); parameter.vval_ = []; end;
vval_ = parameter.vval_;
if ~isfield(parameter,'flag_k_c_interp'); parameter.flag_k_c_interp = 1; end;
flag_k_c_interp = parameter.flag_k_c_interp;
if ~isfield(parameter,'n_order'); parameter.n_order = 5; end;
n_order = parameter.n_order;
if ~isfield(parameter,'flag_pad_xxnufft'); parameter.flag_pad_xxnufft = 0; end;
flag_pad_xxnufft = parameter.flag_pad_xxnufft;
if ~isfield(parameter,'n_x_c'); parameter.n_x_c = 128; end;
n_x_c = parameter.n_x_c;
n_k_c = n_x_c;
if ~isfield(parameter,'diameter_x_c'); parameter.diameter_x_c = 2.0; end;
diameter_x_c = parameter.diameter_x_c;
half_diameter_x_c = diameter_x_c/2.0;
if ~isfield(parameter,'c_use__'); parameter.c_use__ = colormap_80s; end;
c_use__ = parameter.c_use__;
n_c_use = size(c_use__,1);
if ~isfield(parameter,'k_p_r_max_use'); parameter.k_p_r_max_use = k_p_r_max; end;
k_p_r_max_use = parameter.k_p_r_max_use;
if ~isfield(parameter,'flag_patch_vs_line'); parameter.flag_patch_vs_line = 1; end;
flag_patch_vs_line = parameter.flag_patch_vs_line;
if ~isfield(parameter,'linewidth_3d_use'); parameter.linewidth_3d_use = 3; end;
linewidth_3d_use = parameter.linewidth_3d_use;
if ~isfield(parameter,'linewidth_contour_use'); parameter.linewidth_contour_use = 1; end;
linewidth_contour_use = parameter.linewidth_contour_use;
if ~isfield(parameter,'linewidth_grid_use'); parameter.linewidth_grid_use = 1; end;
linewidth_grid_use = parameter.linewidth_grid_use;
if ~isfield(parameter,'shadow_offset'); parameter.shadow_offset = 0.5; end;
shadow_offset = parameter.shadow_offset;
if ~isfield(parameter,'shadow_gamma'); parameter.shadow_gamma = 0.5; end;
shadow_gamma = parameter.shadow_gamma;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

if isempty(n_w_); n_w_ = min(64,n_k_c)*ones(n_k_p_r,1); end;
n_w_ = n_w_(:); n_w_max = max(n_w_); n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_(:)]);
x_c_0_lim_ = half_diameter_x_c*[-1,+1];
x_c_1_lim_ = half_diameter_x_c*[-1,+1];
x_c_2_lim_ = half_diameter_x_c*[-1,+1];
dx = diameter_x_c/n_x_c;
k_c_0_lim_ = k_p_r_max*[-1,+1];
k_c_1_lim_ = k_p_r_max*[-1,+1];
k_c_2_lim_ = k_p_r_max*[-1,+1];
dk = 2*k_p_r_max/n_k_c;
x_c_0_ = transpose(linspace(-half_diameter_x_c,+half_diameter_x_c,n_x_c));
x_c_1_ = transpose(linspace(-half_diameter_x_c,+half_diameter_x_c,n_x_c));
x_c_2_ = transpose(linspace(-half_diameter_x_c,+half_diameter_x_c,n_x_c));
k_c_0_ = transpose(linspace(-k_p_r_max,+k_p_r_max,n_k_c));
k_c_1_ = transpose(linspace(-k_p_r_max,+k_p_r_max,n_k_c));
k_c_2_ = transpose(linspace(-k_p_r_max,+k_p_r_max,n_k_c));
[x_c_0__,x_c_1__] = ndgrid(x_c_0_,x_c_1_);
x_c_r__ = sqrt(x_c_0__.^2 + x_c_1__.^2);
[k_c_0__,k_c_1__] = ndgrid(k_c_0_,k_c_1_);
k_c_r__ = sqrt(k_c_0__.^2 + k_c_1__.^2);
k_c_w__ = atan2(k_c_1__,k_c_0__);
%%%%%%%%%%%%%%%%;
if flag_pad_xxnufft==1;
n_pad = 3; n_rep = (n_pad-1)/2;
n_x_c_pad = n_pad*n_x_c;
diameter_x_c_pad = n_pad*diameter_x_c;
half_diameter_x_c_pad = 0.5*diameter_x_c_pad;
x_c_pad_0_lim_ = half_diameter_x_c_pad*[-1,+1];
x_c_pad_1_lim_ = half_diameter_x_c_pad*[-1,+1];
x_c_pad_2_lim_ = half_diameter_x_c_pad*[-1,+1];
dx_c_pad = diameter_x_c_pad/n_x_c_pad;
x_c_pad_0_ = transpose(linspace(-half_diameter_x_c_pad,+half_diameter_x_c_pad,n_x_c_pad));
x_c_pad_1_ = transpose(linspace(-half_diameter_x_c_pad,+half_diameter_x_c_pad,n_x_c_pad));
x_c_pad_2_ = transpose(linspace(-half_diameter_x_c_pad,+half_diameter_x_c_pad,n_x_c_pad));
[x_c_pad_0__,x_c_pad_1__] = ndgrid(x_c_pad_0_,x_c_pad_1_);
x_c_pad_r__ = sqrt(x_c_pad_0__.^2 + x_c_pad_1__.^2);
sigma_pad = 0.25*half_diameter_x_c;
x_c_app__ = exp(-max(0,x_c_r__-half_diameter_x_c).^2/(2*sigma_pad^2));
x_c_pad_app__ = exp(-max(0,x_c_pad_r__-half_diameter_x_c).^2/(2*sigma_pad^2));
end;%if flag_pad_xxnufft==1;
%%%%%%%%%%%%%%%%;
if flag_k_c_interp;
[ ...
 scatter_from_tensor_swk__ ...
] = ...
disk_k_p_scatter_from_tensor_interpolate_n_4( ...
 n_order ...
,n_w_max ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_k_c^2 ...
,k_c_w__(:) ...
,k_c_r__(:) ...
);
end;%if flag_k_c_interp;
%%%%%%%%%%%%%%%%;

if isempty(template_viewing_azimu_b_S_); template_viewing_azimu_b_S_ = zeros(n_S,1); end;
if isempty(template_viewing_polar_a_S_); template_viewing_polar_a_S_ = zeros(n_S,1); end;

Rz = @(azimu_b) [ +cos(azimu_b) -sin(azimu_b)  0 ; +sin(azimu_b) +cos(azimu_b)  0 ;  0  0  1 ] ;
Ry = @(polar_a) [ +cos(polar_a)  0 +sin(polar_a) ;  0  1  0 ; -sin(polar_a)  0 +cos(polar_a) ] ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if n_contour<=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

k_p_r_pre_wk_ = zeros(n_w_sum,1);
k_p_r_mid_wk_ = zeros(n_w_sum,1);
k_p_r_pos_wk_ = zeros(n_w_sum,1);
k_p_w_pre_wk_ = zeros(n_w_sum,1);
k_p_w_mid_wk_ = zeros(n_w_sum,1);
k_p_w_pos_wk_ = zeros(n_w_sum,1);
for nk_p_r=0:n_k_p_r-1;
k_p_r_pre =         0; if nk_p_r>        0; k_p_r_pre = 0.5*(k_p_r_(1+nk_p_r-1) + k_p_r_(1+nk_p_r+0)); end;
k_p_r_mid = k_p_r_(1+nk_p_r);
k_p_r_pos = k_p_r_max; if nk_p_r<n_k_p_r-1; k_p_r_pos = 0.5*(k_p_r_(1+nk_p_r+1) + k_p_r_(1+nk_p_r+0)); end;
n_w = n_w_(1+nk_p_r);
tmp_psi_ = 2*pi*transpose(0:n_w-1)/max(1,n_w); dpsi = mean(diff(tmp_psi_));
tmp_index_ = n_w_csum_(1+nk_p_r):n_w_csum_(1+nk_p_r+1)-1;
k_p_r_pre_wk_(1+tmp_index_) = k_p_r_pre;
k_p_r_mid_wk_(1+tmp_index_) = k_p_r_mid;
k_p_r_pos_wk_(1+tmp_index_) = k_p_r_pos;
k_p_w_pre_wk_(1+tmp_index_) = tmp_psi_ - 0.5*dpsi;
k_p_w_mid_wk_(1+tmp_index_) = tmp_psi_ + 0.0*dpsi;
k_p_w_pos_wk_(1+tmp_index_) = tmp_psi_ + 0.5*dpsi;
end;%for nk_p_r=0:n_k_p_r-1;
k_c_0_ne_wk_ = k_p_r_pos_wk_.*cos(k_p_w_pos_wk_);
k_c_0_nw_wk_ = k_p_r_pre_wk_.*cos(k_p_w_pos_wk_);
k_c_0_sw_wk_ = k_p_r_pre_wk_.*cos(k_p_w_pre_wk_);
k_c_0_se_wk_ = k_p_r_pos_wk_.*cos(k_p_w_pre_wk_);
k_c_1_ne_wk_ = k_p_r_pos_wk_.*sin(k_p_w_pos_wk_);
k_c_1_nw_wk_ = k_p_r_pre_wk_.*sin(k_p_w_pos_wk_);
k_c_1_sw_wk_ = k_p_r_pre_wk_.*sin(k_p_w_pre_wk_);
k_c_1_se_wk_ = k_p_r_pos_wk_.*sin(k_p_w_pre_wk_);
k_c_ne_wk3__ = cat(2,k_c_0_ne_wk_,k_c_1_ne_wk_,zeros(n_w_sum,1));
k_c_nw_wk3__ = cat(2,k_c_0_nw_wk_,k_c_1_nw_wk_,zeros(n_w_sum,1));
k_c_sw_wk3__ = cat(2,k_c_0_sw_wk_,k_c_1_sw_wk_,zeros(n_w_sum,1));
k_c_se_wk3__ = cat(2,k_c_0_se_wk_,k_c_1_se_wk_,zeros(n_w_sum,1));

n_1 = 1; n_3 = 3;
for nS=0:n_S-1;
%%%%%%%%;
template_viewing_azimu_b = template_viewing_azimu_b_S_(1+nS);
template_viewing_polar_a = template_viewing_polar_a_S_(1+nS);
S_k_p_wk_ = S_k_p_wkS__(:,1+nS);
nc_use_ = max(0,min(n_c_use-1,floor(n_c_use*(real(S_k_p_wk_)-min(vlim_))/diff(vlim_))));
tmp_c1_1wk3___ = reshape(c_use__(1+nc_use_,:),[1,n_w_sum,n_3]);
tmp_cg_1wk3___ = (shadow_offset + (1-shadow_offset)*tmp_c1_1wk3___).^shadow_gamma;
tmp_k_c_ne_wk3__ = k_c_ne_wk3__*transpose(Ry(template_viewing_polar_a))*transpose(Rz(template_viewing_azimu_b));
tmp_k_c_nw_wk3__ = k_c_nw_wk3__*transpose(Ry(template_viewing_polar_a))*transpose(Rz(template_viewing_azimu_b));
tmp_k_c_sw_wk3__ = k_c_sw_wk3__*transpose(Ry(template_viewing_polar_a))*transpose(Rz(template_viewing_azimu_b));
tmp_k_c_se_wk3__ = k_c_se_wk3__*transpose(Ry(template_viewing_polar_a))*transpose(Rz(template_viewing_azimu_b));
tmp_k_c_0_ne_wk_ = tmp_k_c_ne_wk3__(:,1+0); tmp_k_c_1_ne_wk_ = tmp_k_c_ne_wk3__(:,1+1); tmp_k_c_2_ne_wk_ = tmp_k_c_ne_wk3__(:,1+2);
tmp_k_c_0_nw_wk_ = tmp_k_c_nw_wk3__(:,1+0); tmp_k_c_1_nw_wk_ = tmp_k_c_nw_wk3__(:,1+1); tmp_k_c_2_nw_wk_ = tmp_k_c_nw_wk3__(:,1+2);
tmp_k_c_0_sw_wk_ = tmp_k_c_sw_wk3__(:,1+0); tmp_k_c_1_sw_wk_ = tmp_k_c_sw_wk3__(:,1+1); tmp_k_c_2_sw_wk_ = tmp_k_c_sw_wk3__(:,1+2);
tmp_k_c_0_se_wk_ = tmp_k_c_se_wk3__(:,1+0); tmp_k_c_1_se_wk_ = tmp_k_c_se_wk3__(:,1+1); tmp_k_c_2_se_wk_ = tmp_k_c_se_wk3__(:,1+2);
tmp_k_c_0_5wk__ = permute(cat(2,tmp_k_c_0_ne_wk_,tmp_k_c_0_nw_wk_,tmp_k_c_0_sw_wk_,tmp_k_c_0_se_wk_,tmp_k_c_0_ne_wk_),[2,1]);
tmp_k_c_1_5wk__ = permute(cat(2,tmp_k_c_1_ne_wk_,tmp_k_c_1_nw_wk_,tmp_k_c_1_sw_wk_,tmp_k_c_1_se_wk_,tmp_k_c_1_ne_wk_),[2,1]);
tmp_k_c_2_5wk__ = permute(cat(2,tmp_k_c_2_ne_wk_,tmp_k_c_2_nw_wk_,tmp_k_c_2_sw_wk_,tmp_k_c_2_se_wk_,tmp_k_c_2_ne_wk_),[2,1]);
%%%%;
hold on;
if flag_patch_vs_line==1;
p = patch( +0*k_p_r_max_use + 1*tmp_k_c_0_5wk__ , +0*k_p_r_max_use + 1*tmp_k_c_1_5wk__ , +0*k_p_r_max_use + 1*tmp_k_c_2_5wk__ , tmp_c1_1wk3___ ); set(p,'LineStyle','none');
p = patch( +1*k_p_r_max_use + 0*tmp_k_c_0_5wk__ , +0*k_p_r_max_use + 1*tmp_k_c_1_5wk__ , +0*k_p_r_max_use + 1*tmp_k_c_2_5wk__ , tmp_cg_1wk3___ ); set(p,'LineStyle','none');
p = patch( +0*k_p_r_max_use + 1*tmp_k_c_0_5wk__ , +1*k_p_r_max_use + 0*tmp_k_c_1_5wk__ , +0*k_p_r_max_use + 1*tmp_k_c_2_5wk__ , tmp_cg_1wk3___ ); set(p,'LineStyle','none');
p = patch( +0*k_p_r_max_use + 1*tmp_k_c_0_5wk__ , +0*k_p_r_max_use + 1*tmp_k_c_1_5wk__ , -1*k_p_r_max_use + 0*tmp_k_c_2_5wk__ , tmp_cg_1wk3___ ); set(p,'LineStyle','none');
end;%if flag_patch_vs_line==1;
if flag_patch_vs_line==0;
%<-- do nothing. ;
end;%if flag_patch_vs_line==0;
hold off;
%%%%;
%%%%%%%%;
end;%for nS=0:n_S-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if n_contour<=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if n_contour>=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

S_k_c_01S___ = zeros(n_k_c,n_k_c,n_S);
for nS=0:n_S-1;
S_k_p_ = S_k_p_wkS__(:,1+nS);
S_k_p_l2 = sum(abs(S_k_p_).^2.*weight_2d_k_all_)*(2*pi)^2;
%%%%%%%%%%%%%%%%;
if flag_k_c_interp==0;
%%%%%%%%%%%%%%%%;
S_x_c_ = interp_k_p_to_x_c_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,n_k_p_r,k_p_r_,n_w_,S_k_p_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_c^2) * n_w_sum;
S_x_c__ = reshape(S_x_c_,[n_x_c,n_x_c]);
S_x_c_l2 = sum(abs(S_x_c__).^2,'all')*dx^2;
if flag_pad_xxnufft==1;
eta_pad = pi/diameter_x_c_pad;
S_x_c_pad__ = ...
[ ...
 repmat(S_x_c__(end,end),[n_rep*n_x_c,n_rep*n_x_c]) , repmat(S_x_c__(end,  :),[n_rep*n_x_c,          1]) , repmat(S_x_c__(  1,end),[n_rep*n_x_c,n_rep*n_x_c]) , 
 repmat(S_x_c__(  :,end),[          1,n_rep*n_x_c]) , repmat(S_x_c__(  :,  :),[          1,          1]) , repmat(S_x_c__(  :,  1),[          1,n_rep*n_x_c]) ,
% repmat(S_x_c__(  :,end),[          1,n_rep*n_x_c]) , repmat(S_x_c__(  1,  1),[    1*n_x_c,    1*n_x_c]) , repmat(S_x_c__(  :,  1),[          1,n_rep*n_x_c]) , 
 repmat(S_x_c__(end,  1),[n_rep*n_x_c,n_rep*n_x_c]) , repmat(S_x_c__(  1,  :),[n_rep*n_x_c,          1]) , repmat(S_x_c__(  1,  1),[n_rep*n_x_c,n_rep*n_x_c]) , 
].*x_c_pad_app__;
S_k_c__ = reshape( ...
		   xxnufft2d3( ...
			       n_x_c_pad^2 ...
			       ,x_c_pad_0__(:)*eta_pad ...
			       ,x_c_pad_1__(:)*eta_pad ...
			       ,S_x_c_pad__(:)*dx_c_pad^2 ...
			       ,-1 ...
			       ,1e-7 ...
			       ,n_k_c^2 ...
			       ,2*pi*k_c_0__(:)/eta_pad ...
			       ,2*pi*k_c_1__(:)/eta_pad ...
			       ) ...
		   ,[n_k_c,n_k_c] ...
		   );
end;%if flag_pad_xxnufft==1;
if flag_pad_xxnufft==0;
eta = pi/diameter_x_c;
S_k_c__ = reshape(xxnufft2d3(n_x_c^2,x_c_0__(:)*eta,x_c_1__(:)*eta,S_x_c__(:)*dx^2,-1,1e-7,n_k_c^2,2*pi*k_c_0__(:)/eta,2*pi*k_c_1__(:)/eta),[n_k_c,n_k_c]);
end;%if flag_pad_xxnufft==1;
S_k_c_l2 = sum(abs(S_k_c__).^2,'all')*dk^2;
%%%%%%%%;
flag_check=0;
if flag_check | flag_verbose>1;
disp(sprintf(' %% S_k_p_l2: %0.6f',S_k_p_l2));
disp(sprintf(' %% S_x_c_l2: %0.6f',S_x_c_l2));
disp(sprintf(' %% S_x_c_l2/S_k_p_l2: %0.6f',S_x_c_l2/S_k_p_l2));
disp(sprintf(' %% S_k_c_l2: %0.6f',S_k_c_l2));
end;%if flag_check | flag_verbose>1;
if flag_check;
figure(1);clf;figmed;
subplot(1,3,1); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(S_k_p_),vlim_,colormap_80s); axis image; axisnotick;
title('S_k_p_','Interpreter','none');
if flag_pad_xxnufft==0;
subplot(1,3,2); imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,(real(S_x_c__)),[],colormap_beach); axis(half_diameter_x_c*[-1,1,-1,1]); axis equal; axisnotick;
end;%if flag_pad_xxnufft==0;
title('(S_x_c__)','Interpreter','none');
if flag_pad_xxnufft==1;
subplot(1,3,2); imagesc_c(n_x_c_pad,x_c_pad_0_,n_x_c_pad,x_c_pad_1_,(real(S_x_c_pad__)),[],colormap_beach); axis(half_diameter_x_c_pad*[-1,1,-1,1]); axis equal; axisnotick;
title('(S_x_c_pad__)','Interpreter','none');
end;%if flag_pad_xxnufft==1;
subplot(1,3,3); imagesc_c(n_k_c,k_c_0_,n_k_c,k_c_1_,real(S_k_c__),vlim_,colormap_80s); axis(k_p_r_max*[-1,1,-1,1]); axis equal; axisnotick;
title('S_k_c__','Interpreter','none');
error('stopping');
end;%if flag_check;
%%%%%%%%%%%%%%%%;
end;%if flag_k_c_interp==0;
%%%%%%%%%%%%%%%%;
if flag_k_c_interp==1;
%%%%%%%%%%%%%%%%;
S_k_c__ = reshape(scatter_from_tensor_swk__*S_k_p_,[n_k_c,n_k_c]);
S_k_c_l2 = sum(abs(S_k_c__).^2,'all')*dk^2;
flag_check=0;
if flag_check | flag_verbose>1;
disp(sprintf(' %% S_k_p_l2: %0.6f',S_k_p_l2));
disp(sprintf(' %% S_k_c_l2: %0.6f',S_k_c_l2));
end;%if flag_check | flag_verbose>1;
if flag_check;
figure(1);clf;figmed;
subplot(1,2,1); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(S_k_p_),vlim_,colormap_80s); axis image; axisnotick;
title('S_k_p_','Interpreter','none');
subplot(1,2,2); imagesc_c(n_k_c,k_c_0_,n_k_c,k_c_1_,real(S_k_c__),vlim_,colormap_80s); axis(k_p_r_max*[-1,1,-1,1]); axis equal; axisnotick;
title('S_k_c__','Interpreter','none');
error('stopping');
end;%if flag_check;
%%%%%%%%%%%%%%%%;
end;%if flag_k_c_interp==1;
%%%%%%%%%%%%%%%%;
S_k_c_01S___(:,:,1+nS) = S_k_c__;
end;%for nS=0:n_S-1;
if flag_percent_use==1;
vval = prctile(abs(S_k_c_01S___),percent_threshold,'all'); vlim_ = vval*[-1,+1];
vval_ = transpose(linspace(min(vlim_),max(vlim_),n_contour));
end;%if flag_percent_use==1;
if flag_percent_use==0;
n_contour = numel(vval_);
end;%if flag_percent_use==0;

flag_check=0;
if flag_check;
figure(1);clf;figsml;
hold on;
for nS=0:n_S-1;
template_viewing_azimu_b = 2*pi*nS/n_S;
tmp_ = Rz(template_viewing_azimu_b)*[0;1;0];
plot3(tmp_(1),tmp_(2),tmp_(3),'o');
end;%for nS=0:n_S-1;
axis vis3d;
hold off;
error('stopping');
end;%if flag_check;

%plot(1:n_S,template_viewing_azimu_b_S_,'o-',1:n_S,2*pi*(0:n_S-1)/n_S,'x-'); error('stopping');

if flag_check;
figure(1);clf;figsml;
hold on;
for nS=0:n_S-1;
template_viewing_azimu_b = template_viewing_azimu_b_S_(1+nS);
template_viewing_polar_a = template_viewing_polar_a_S_(1+nS);
%tmp_ = Rz(template_viewing_azimu_b)*Ry(template_viewing_polar_a)*[0;1;0];
tmp_ = Rz(template_viewing_azimu_b)*[0;1;0];
plot3(tmp_(1),tmp_(2),tmp_(3),'o');
end;%for nS=0:n_S-1;
axis vis3d;
hold off;
error('stopping');
end;%if flag_check;

for nS=0:n_S-1;
%%%%%%%%;
template_viewing_azimu_b = template_viewing_azimu_b_S_(1+nS);
template_viewing_polar_a = template_viewing_polar_a_S_(1+nS);
S_k_c__ = S_k_c_01S___(:,:,1+nS);
[contour_cut__] = contourc(transpose(real(S_k_c__)),vval_);
[n_item,val_i_,length_i_,k_c_0_il__,k_c_1_il__] = cell_from_contour_0(contour_cut__);
nitem_ = 0:n_item-1; if (mod(n_item,2)==1); nitem_ = [nitem_,floor(n_item/2)]; end; tmp_n2 = numel(nitem_)/2; tmp_n = 2*tmp_n2;
nitem_ = reshape(nitem_,[tmp_n2,2]);
nitem_(:,1) = flipud(nitem_(:,1));
nitem_ = reshape(transpose(nitem_),1,2*tmp_n2);
nitem_ = unique(nitem_,'stable');
%%%%;
hold on;
%%;
tmp_k_upd = @(k) k_p_r_max * (k-0.5 - n_k_c/2)/(n_k_c/2);
%%;
for sign_val=[-1,+1];
tmpn=0;
for nitem=nitem_;
val = val_i_(1+nitem);
nc_use = max(0,min(n_c_use-1,floor(n_c_use*(val-min(vlim_))/diff(vlim_))));
tmp_length = length_i_(1+nitem);
tmp_k_c_0_ = tmp_k_upd(k_c_0_il__{1+nitem});
tmp_k_c_1_ = tmp_k_upd(k_c_1_il__{1+nitem});
tmp_k_c_2_ = zeros(tmp_length,1) + sign_val*tolerance_stack*k_p_r_max_use*tmpn/tmp_n;
tmp_k_c_ld__ = transpose(Rz(template_viewing_azimu_b)*Ry(template_viewing_polar_a)*transpose([ tmp_k_c_0_ , tmp_k_c_1_ , tmp_k_c_2_ ]));
tmp_k_c_0_ = tmp_k_c_ld__(:,1+0);
tmp_k_c_1_ = tmp_k_c_ld__(:,1+1);
tmp_k_c_2_ = tmp_k_c_ld__(:,1+2);
if flag_patch_vs_line==1;
p = patch( tmp_k_c_0_ , tmp_k_c_1_ , tmp_k_c_2_ , c_use__(1+nc_use,:) ); set(p,'LineStyle','none');
end;%if flag_patch_vs_line==1;
if flag_patch_vs_line==0;
l = line( tmp_k_c_0_ , tmp_k_c_1_ , tmp_k_c_2_ ); set(l,'Color',c_use__(1+nc_use,:),'LineWidth',linewidth_3d_use);
end;%if flag_patch_vs_line==0;
tmpn=tmpn+1;
end;%for nitem=nitem_;
end;%for sign_val=[-1,+1];
%%;
tmpn=0;
for nitem=nitem_;
val = val_i_(1+nitem);
nc_use = max(0,min(n_c_use-1,floor(n_c_use*(val-min(vlim_))/diff(vlim_))));
tmp_length = length_i_(1+nitem);
tmp_k_c_0_ = tmp_k_upd(k_c_0_il__{1+nitem});
tmp_k_c_1_ = tmp_k_upd(k_c_1_il__{1+nitem});
tmp_k_c_2_ = zeros(tmp_length,1);
tmp_k_c_ld__ = transpose(Rz(template_viewing_azimu_b)*Ry(template_viewing_polar_a)*transpose([ tmp_k_c_0_ , tmp_k_c_1_ , tmp_k_c_2_ ]));
tmp_k_c_0_ = tmp_k_c_ld__(:,1+0);
tmp_k_c_1_ = tmp_k_c_ld__(:,1+1);
tmp_k_c_2_ = tmp_k_c_ld__(:,1+2);
tmp_c_ = c_use__(1+nc_use,:);
tmp_c_ = (shadow_offset + (1-shadow_offset)*tmp_c_).^shadow_gamma;
l = line( +k_p_r_max_use*ones(tmp_length,1) , tmp_k_c_1_ , tmp_k_c_2_ ); set(l,'LineWidth',linewidth_contour_use,'Color',tmp_c_);
l = line( tmp_k_c_0_ , +k_p_r_max_use*ones(tmp_length,1) , tmp_k_c_2_ ); set(l,'LineWidth',linewidth_contour_use,'Color',tmp_c_);
l = line( tmp_k_c_0_ , tmp_k_c_1_ , -k_p_r_max_use*ones(tmp_length,1) ); set(l,'LineWidth',linewidth_contour_use,'Color',tmp_c_);
tmpn=tmpn+1;
end;%for nitem=nitem_;
%%;
hold off;
%%%%;
%%%%%%%%;
end;%for nS=0:n_S-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if n_contour>=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% plot grid backdrop. ;
%%%%%%%%;
n_grid = 1+8; color_use = 0.65*[1,1,1];
hold on;
for ngrid=0:n_grid-1;
tmp_k = -1+ngrid/max(1,n_grid-1)*2;
line(k_p_r_max_use*[+1,+1],k_p_r_max_use*[-1,+1],k_p_r_max_use*tmp_k*[+1,+1],'LineWidth',linewidth_grid_use,'Color',color_use);
line(k_p_r_max_use*[+1,+1],k_p_r_max_use*tmp_k*[+1,+1],k_p_r_max_use*[-1,+1],'LineWidth',linewidth_grid_use,'Color',color_use);
line(k_p_r_max_use*[-1,+1],k_p_r_max_use*[+1,+1],k_p_r_max_use*tmp_k*[+1,+1],'LineWidth',linewidth_grid_use,'Color',color_use);
line(k_p_r_max_use*tmp_k*[+1,+1],k_p_r_max_use*[+1,+1],k_p_r_max_use*[-1,+1],'LineWidth',linewidth_grid_use,'Color',color_use);
line(k_p_r_max_use*[-1,+1],k_p_r_max_use*tmp_k*[+1,+1],k_p_r_max_use*[-1,-1],'LineWidth',linewidth_grid_use,'Color',color_use);
line(k_p_r_max_use*tmp_k*[+1,+1],k_p_r_max_use*[-1,+1],k_p_r_max_use*[-1,-1],'LineWidth',linewidth_grid_use,'Color',color_use);
end;%for ngrid=0:n_grid-1;
hold off;

xlim((1+tolerance_stack)*k_p_r_max_use*[-1,+1]); ylim((1+tolerance_stack)*k_p_r_max_use*[-1,+1]); zlim((1+tolerance_stack)*k_p_r_max_use*[-1,+1]); axis vis3d;
view([-65,20]);

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
