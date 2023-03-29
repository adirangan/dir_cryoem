function imagesc_S_k_p_3d_0( ...
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

str_thisfunction = 'imagesc_S_k_p_3d_0';

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
if ~isfield(parameter,'tolerance_stack'); parameter.tolerance_stack = 1e-2; end;
tolerance_stack = parameter.tolerance_stack;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
flag_verbose = parameter.flag_verbose;
if ~isfield(parameter,'percent_threshold'); parameter.percent_threshold = 98.5; end;
percent_threshold = parameter.percent_threshold;
if ~isfield(parameter,'n_contour'); parameter.n_contour = 4; end;
n_contour = parameter.n_contour;
if ~isfield(parameter,'n_x_c'); parameter.n_x_c = 128; end;
n_x_c = parameter.n_x_c;
n_k_c = n_x_c;
if ~isfield(parameter,'diameter_x_c'); parameter.diameter_x_c = 2.0; end;
diameter_x_c = parameter.diameter_x_c;
half_diameter_x_c = diameter_x_c/2.0;
if ~isfield(parameter,'c_use__'); parameter.c_use__ = colormap_80s; end;
c_use__ = parameter.c_use__;
n_c_use = size(c_use__,1);

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

n_w_ = n_w_(:); n_w_max = max(n_w_); n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_(:)]);
x_0_lim_ = half_diameter_x_c*[-1,+1];
x_1_lim_ = half_diameter_x_c*[-1,+1];
x_2_lim_ = half_diameter_x_c*[-1,+1];
k_0_lim_ = k_p_r_max*[-1,+1];
k_1_lim_ = k_p_r_max*[-1,+1];
k_2_lim_ = k_p_r_max*[-1,+1];

if isempty(template_viewing_azimu_b_S_); template_viewing_azimu_b_S_ = zeros(n_S,1); end;
if isempty(template_viewing_polar_a_S_); template_viewing_polar_a_S_ = zeros(n_S,1); end;

Rz = @(azimu_b) [ +cos(azimu_b) -sin(azimu_b)  0 ; +sin(azimu_b) +cos(azimu_b)  0 ;  0  0  1 ] ;
Ry = @(polar_a) [ +cos(polar_a)  0 +sin(polar_a) ;  0  1  0 ; -sin(polar_a)  0 +cos(polar_a) ] ;

S_k_c_01S___ = zeros(n_k_c,n_k_c,n_S);
for nS=0:n_S-1;
S_k_p_ = S_k_p_wkS__(:,1+nS);
S_x_c_ = interp_k_p_to_x_c_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,n_k_p_r,k_p_r_,n_w_,S_k_p_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_c^2) * n_w_sum;
S_x_c__ = reshape(S_x_c_,[n_x_c,n_x_c]);
S_k_c__ = fftshift(fft2(fftshift(S_x_c__)))/n_k_c;
S_k_c_01S___(:,:,1+nS) = S_k_c__;
end;%for nS=0:n_S-1;
vval = prctile(abs(S_k_c_01S___),percent_threshold,'all'); vlim_ = vval*[-1,+1];

%{
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
%}

%plot(1:n_S,template_viewing_azimu_b_S_,'o-',1:n_S,2*pi*(0:n_S-1)/n_S,'x-'); error('stopping');

%{
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
%}

for nS=0:n_S-1;
%%%%%%%%;
template_viewing_azimu_b = template_viewing_azimu_b_S_(1+nS);
template_viewing_polar_a = template_viewing_polar_a_S_(1+nS);
S_k_c__ = S_k_c_01S___(:,:,1+nS);
[contour_cut__] = contourc(transpose(real(S_k_c__)),linspace(min(vlim_),max(vlim_),n_contour));
[n_item,val_i_,length_i_,k_c_0_il__,k_c_1_il__] = cell_from_contour_0(contour_cut__);
nitem_ = 0:n_item-1; if (mod(n_item,2)==1); nitem_ = [nitem_,floor(n_item/2)]; end; tmp_n2 = numel(nitem_)/2; tmp_n = 2*tmp_n2;
nitem_ = reshape(nitem_,[tmp_n2,2]);
nitem_(:,1) = flipud(nitem_(:,1));
nitem_ = reshape(transpose(nitem_),1,2*tmp_n2);
nitem_ = unique(nitem_,'stable');
%%%%;
hold on;
%%;
tmpn=0;
for nitem=nitem_;
val = val_i_(1+nitem);
nc_use = max(0,min(n_c_use-1,floor(n_c_use*(val-min(vlim_))/diff(vlim_))));
tmp_length = length_i_(1+nitem);
tmp_k_c_0_ = k_c_0_il__{1+nitem}; tmp_k_c_0_ = k_p_r_max*(tmp_k_c_0_ - n_k_c/2 - 0.5)*2/n_k_c;
tmp_k_c_1_ = k_c_1_il__{1+nitem}; tmp_k_c_1_ = k_p_r_max*(tmp_k_c_1_ - n_k_c/2 - 0.5)*2/n_k_c;
tmp_k_c_2_ = zeros(tmp_length,1) + tolerance_stack*k_p_r_max*tmpn/tmp_n;
tmp_k_c_ld__ = transpose(Rz(template_viewing_azimu_b)*Ry(template_viewing_polar_a)*transpose([ tmp_k_c_0_ , tmp_k_c_1_ , tmp_k_c_2_ ]));
tmp_k_c_0_ = tmp_k_c_ld__(:,1+0);
tmp_k_c_1_ = tmp_k_c_ld__(:,1+1);
tmp_k_c_2_ = tmp_k_c_ld__(:,1+2);
p = patch( tmp_k_c_0_ , tmp_k_c_1_ , tmp_k_c_2_ , c_use__(1+nc_use,:) );
set(p,'LineStyle','none');
tmpn=tmpn+1;
end;%for nitem=nitem_;
%%;
tmpn=0;
for nitem=nitem_;
val = val_i_(1+nitem);
nc_use = max(0,min(n_c_use-1,floor(n_c_use*(val-min(vlim_))/diff(vlim_))));
tmp_length = length_i_(1+nitem);
tmp_k_c_0_ = k_c_0_il__{1+nitem}; tmp_k_c_0_ = k_p_r_max*(tmp_k_c_0_ - n_k_c/2 - 0.5)*2/n_k_c;
tmp_k_c_1_ = k_c_1_il__{1+nitem}; tmp_k_c_1_ = k_p_r_max*(tmp_k_c_1_ - n_k_c/2 - 0.5)*2/n_k_c;
tmp_k_c_2_ = zeros(tmp_length,1) - tolerance_stack*k_p_r_max*tmpn/tmp_n;
tmp_k_c_ld__ = transpose(Rz(template_viewing_azimu_b)*Ry(template_viewing_polar_a)*transpose([ tmp_k_c_0_ , tmp_k_c_1_ , tmp_k_c_2_ ]));
tmp_k_c_0_ = tmp_k_c_ld__(:,1+0);
tmp_k_c_1_ = tmp_k_c_ld__(:,1+1);
tmp_k_c_2_ = tmp_k_c_ld__(:,1+2);
p = patch( tmp_k_c_0_ , tmp_k_c_1_ , tmp_k_c_2_ , c_use__(1+nc_use,:) );
set(p,'LineStyle','none');
tmpn=tmpn+1;
end;%for nitem=nitem_;
%%;
hold off;
%%%%;
%%%%%%%%;
end;%for nS=0:n_S-1;

xlim(k_p_r_max*[-1,+1]); ylim(k_p_r_max*[-1,+1]); zlim(k_p_r_max*[-1,+1]); axis vis3d;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
