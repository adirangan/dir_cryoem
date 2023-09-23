function imagesc_S_k_p_3d_belt_3( ...
 parameter ...
,n_S ...
,template_viewing_azimu_b_S_ ...
,template_viewing_polar_a_S_ ...
);

str_thisfunction = 'imagesc_S_k_p_3d_2';

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_S=[]; end; na=na+1;
if (nargin<1+na); template_viewing_azimu_b_S_=[]; end; na=na+1;
if (nargin<1+na); template_viewing_polar_a_S_=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;
if ~isfield(parameter,'tolerance_master'); parameter.tolerance_master = 1e-2; end;
tolerance_master = parameter.tolerance_master;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
flag_verbose = parameter.flag_verbose;
if ~isfield(parameter,'k_p_r_max'); parameter.k_p_r_max = 1.0; end;
k_p_r_max = parameter.k_p_r_max;
if ~isfield(parameter,'n_gamma_z'); parameter.n_gamma_z = 1024; end;
n_gamma_z = parameter.n_gamma_z;
if ~isfield(parameter,'flag_fill'); parameter.flag_fill = 1; end;
flag_fill = parameter.flag_fill;
if ~isfield(parameter,'c_fill_use_'); parameter.c_fill_use_ = [1,0,1]; end;
c_fill_use_ = parameter.c_fill_use_;
if ~isfield(parameter,'linewidth_use'); parameter.linewidth_use = 1.0; end;
linewidth_use = parameter.linewidth_use;
if ~isfield(parameter,'c_line_use_'); parameter.c_line_use_ = [0,1,0]; end;
c_line_use_ = parameter.c_line_use_;
if ~isfield(parameter,'markersize_use'); parameter.markersize_use = 8.0; end;
markersize_use = parameter.markersize_use;
if ~isfield(parameter,'c_mark_use_'); parameter.c_mark_use_ = [0,1,1]; end;
c_mark_use_ = parameter.c_mark_use_;


if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

if isempty(n_S); n_S = 1; end;
if isempty(template_viewing_azimu_b_S_); template_viewing_azimu_b_S_ = zeros(n_S,1); end;
if isempty(template_viewing_polar_a_S_); template_viewing_polar_a_S_ = zeros(n_S,1); end;

Rz = @(azimu_b) [ +cos(azimu_b) -sin(azimu_b)  0 ; +sin(azimu_b) +cos(azimu_b)  0 ;  0  0  1 ] ;
Ry = @(polar_a) [ +cos(polar_a)  0 +sin(polar_a) ;  0  1  0 ; -sin(polar_a)  0 +cos(polar_a) ] ;
gamma_z_ = linspace(0,2*pi,1+n_gamma_z); gamma_z_ = transpose(gamma_z_(1:n_gamma_z));

k_pre_0_z_ = k_p_r_max*cos(gamma_z_);
k_pre_1_z_ = k_p_r_max*sin(gamma_z_);
k_pre_2_z_ = zeros(n_gamma_z,1);
hold on;
for nS=0:n_S-1;
template_viewing_azimu_b = template_viewing_azimu_b_S_(1+nS);
template_viewing_polar_a = template_viewing_polar_a_S_(1+nS);
tmp__ = transpose([k_pre_0_z_,k_pre_1_z_,k_pre_2_z_]);
tmp__ = Ry(template_viewing_polar_a)*tmp__;
tmp__ = Rz(template_viewing_azimu_b)*tmp__;
tmp__ = transpose(tmp__);
k_pos_0_z_ = tmp__(:,1+0); k_pos_1_z_ = tmp__(:,1+1); k_pos_2_z_ = tmp__(:,1+2);
tmp_p=patch(k_pos_0_z_,k_pos_1_z_,k_pos_2_z_,c_fill_use_);
if flag_fill==0; set(tmp_p,'FaceColor','none'); end;
if linewidth_use==0; set(tmp_p,'LineStyle','none'); end;
if linewidth_use> 0; set(tmp_p,'LineWidth',linewidth_use); set(tmp_p,'EdgeColor',c_line_use_); end;
tmp__ = transpose([0,0,k_p_r_max*1.25]);
tmp__ = Ry(template_viewing_polar_a)*tmp__;
tmp__ = Rz(template_viewing_azimu_b)*tmp__;
tmp__ = transpose(tmp__);
k_pos_0_z_ = tmp__(:,1+0); k_pos_1_z_ = tmp__(:,1+1); k_pos_2_z_ = tmp__(:,1+2);
plot3(k_pos_0_z_,k_pos_1_z_,k_pos_2_z_,'.','MarkerSize',markersize_use,'MarkerEdgeColor',c_mark_use_);
end;%for nS=0:n_S-1;
hold off;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
