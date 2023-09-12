function ...
[ ...
 parameter ...
] = ...
imagesc_p_strip_0( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,n_w_sum ...
,S_k_p_ ...
,n_gamma_z ...
,gamma_z_ ...
);

str_thisfunction = 'imagesc_p_strip_0';

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); k_p_r_=[]; end; na=na+1;
if (nargin<1+na); k_p_r_max=[]; end; na=na+1;
if (nargin<1+na); n_w_=[]; end; na=na+1;
if (nargin<1+na); n_w_sum=[]; end; na=na+1;
if (nargin<1+na); S_k_p_=[]; end; na=na+1;
if (nargin<1+na); n_gamma_z=[]; end; na=na+1;
if (nargin<1+na); gamma_z_=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
flag_verbose = parameter.flag_verbose;
if ~isfield(parameter,'strip_width'); parameter.strip_width = 0.125; end;
strip_width = parameter.strip_width;
if ~isfield(parameter,'clim_'); parameter.clim_ = prctile(abs(S_k_p_),100,'all')*[-1,+1]; end;
clim_ = parameter.clim_;
if ~isfield(parameter,'c_use__'); parameter.c_use__ = colormap_beach; end;
c_use__ = parameter.c_use__;
if ~isfield(parameter,'k_scale'); parameter.k_scale = 1.0; end;
k_scale = parameter.k_scale;
if ~isfield(parameter,'k_0_offset'); parameter.k_0_offset = 0.0; end;
k_0_offset = parameter.k_0_offset;
if ~isfield(parameter,'k_1_offset'); parameter.k_1_offset = 0.0; end;
k_1_offset = parameter.k_1_offset;
if ~isfield(parameter,'gamma_offset'); parameter.gamma_offset = 0.0; end;
gamma_offset = parameter.gamma_offset;
if ~isfield(parameter,'gamma_offset_'); if ~isempty(n_gamma_z); parameter.gamma_offset_ = 0.0*gamma_z_; else parameter.gamma_offset_ = 0.0*zeros(1,n_w_sum); end; end;
gamma_offset_ = parameter.gamma_offset_;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

n_c_use = size(c_use__,1);
n_w_csum_ = cumsum([0;n_w_(:)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if ~isempty(n_gamma_z);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%;
for ngamma_z=0:n_gamma_z-1;
%%%%%%%%;
gamma_z = gamma_z_(1+ngamma_z); gamma_z = periodize(gamma_z,0,2*pi);
%%%%;
S_pos_k_ = zeros(n_k_p_r,1);
S_neg_k_ = zeros(n_k_p_r,1);
r_pos_k_ = zeros(n_k_p_r,1);
r_neg_k_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
n_w = n_w_(1+nk_p_r);
S_w_ = S_k_p_(1+n_w_csum_(1+nk_p_r)+[0:n_w-1]);
nw_pos = periodize(round(n_w*gamma_z/(2*pi)),0,n_w);
nw_neg = periodize(nw_pos+n_w/2,0,n_w); %<-- assume n_w is divisible by 2. ;
S_pos_k_(1+nk_p_r) = S_w_(1+nw_pos);
r_pos_k_(1+nk_p_r) = k_p_r;
S_neg_k_(1+nk_p_r) = S_w_(1+nw_neg);
r_neg_k_(1+nk_p_r) = k_p_r;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%;
S_k_ = [+flipud(S_neg_k_);+S_pos_k_];
r_neg_k_ = -flipud(r_neg_k_);
r_neg_lob_k_ = [-k_p_r_max;+0.5*r_neg_k_(1:end-1)+0.5*r_neg_k_(2:end)];
r_neg_upb_k_ = [+0.5*r_neg_k_(1:end-1)+0.5*r_neg_k_(2:end);-0.0];
r_pos_lob_k_ = [+0.0;+0.5*r_pos_k_(1:end-1)+0.5*r_pos_k_(2:end)];
r_pos_upb_k_ = [+0.5*r_pos_k_(1:end-1)+0.5*r_pos_k_(2:end);+k_p_r_max];
r_lob_k_ = [r_neg_lob_k_;r_pos_lob_k_];
r_upb_k_ = [r_neg_upb_k_;r_pos_upb_k_];
k_0_k__ = [ r_lob_k_ , r_upb_k_ , r_upb_k_ , r_lob_k_ ];
k_0_k__ = transpose(k_0_k__);
k_1_k__ = k_p_r_max*strip_width*[ -ones(2*n_k_p_r,1) , -ones(2*n_k_p_r,1) , +ones(2*n_k_p_r,1) , +ones(2*n_k_p_r,1) ];
k_1_k__ = transpose(k_1_k__);
R__ = [ ...
 +cos(gamma_z+gamma_offset+gamma_offset_(1+ngamma_z)) ...
,-sin(gamma_z+gamma_offset+gamma_offset_(1+ngamma_z)) ...
;+sin(gamma_z+gamma_offset+gamma_offset_(1+ngamma_z)) ...
,+cos(gamma_z+gamma_offset+gamma_offset_(1+ngamma_z)) ...
];
k_k__ = R__ * [ reshape(k_0_k__,[1,4*2*n_k_p_r]) ; reshape(k_1_k__,[1,4*2*n_k_p_r]) ];
k_0_k__ = reshape(k_k__(1+0,:),[4,2*n_k_p_r]);
k_1_k__ = reshape(k_k__(1+1,:),[4,2*n_k_p_r]);
nc_ = max(0,min(n_c_use-1,floor(n_c_use*(real(S_k_)-min(clim_))/diff(clim_))));
c_ = reshape(c_use__(1+nc_,:),[1,2*n_k_p_r,3]);
patch(k_0_k__+k_0_offset,k_1_k__+k_1_offset,c_,'LineStyle','none');
%%%%%%%%;
end;%for ngamma_z=0:n_gamma_z-1;
%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if ~isempty(n_gamma_z);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if  isempty(n_gamma_z);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
tmp_k0_ = zeros(1,n_w_sum); tmp_k1_ = zeros(1,n_w_sum);
tmp_w0_ = zeros(1,n_w_sum); tmp_w1_ = zeros(1,n_w_sum);
c_ = zeros(1,n_w_sum,3);
ic=0;
for nk_p_r=0:n_k_p_r-1;
tmp_k = k_p_r_(1+nk_p_r);
if nk_p_r==        0; tmp_k_pre = k_p_r_(1+0); else tmp_k_pre = 0.5*(k_p_r_(1+nk_p_r-1) + k_p_r_(1+nk_p_r)); end;
if nk_p_r==n_k_p_r-1; tmp_k_pos = k_p_r_(end); else tmp_k_pos = 0.5*(k_p_r_(1+nk_p_r+1) + k_p_r_(1+nk_p_r)); end;
n_w = n_w_(1+nk_p_r);
dw = 2*pi/n_w;
for nw=0:n_w_(1+nk_p_r)-1;
tmp_w_pre = nw*dw - 0.5*dw;
tmp_w_pos = nw*dw + 0.5*dw;
tmp_k0_(1+ic) = tmp_k_pre; tmp_k1_(1+ic) = tmp_k_pos;
tmp_w0_(1+ic) = tmp_w_pre; tmp_w1_(1+ic) = tmp_w_pos;
nc = max(0,min(n_c_use-1,floor(n_c_use*(real(S_k_p_(1+ic))-min(clim_))/diff(clim_))));
c_(1,1+ic,:) = c_use__(1+nc,:);
ic=ic+1;
end;%for nw=0:n_w_(1+nk_p_r)-1;
end;%for nk_p_r=0:n_k_p_r-1;
tmp_w0_ = tmp_w0_ + gamma_offset + gamma_offset_;
tmp_w1_ = tmp_w1_ + gamma_offset + gamma_offset_;
k_0_k__ = [tmp_k0_.*cos(tmp_w0_) ; tmp_k1_.*cos(tmp_w0_) ; tmp_k1_.*cos(tmp_w1_) ; tmp_k0_.*cos(tmp_w1_) ; tmp_k0_.*cos(tmp_w0_) ];
k_1_k__ = [tmp_k0_.*sin(tmp_w0_) ; tmp_k1_.*sin(tmp_w0_) ; tmp_k1_.*sin(tmp_w1_) ; tmp_k0_.*sin(tmp_w1_) ; tmp_k0_.*sin(tmp_w0_) ];
p=patch(k_0_k__+k_0_offset,k_1_k__+k_1_offset,c_,'LineStyle','none');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if  isempty(n_gamma_z);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
