function ...
[ ...
 parameter ...
,template_tree ...
] = ...
get_template_tree_1( ...
 parameter ...
,viewing_k_eq_d_min ...
,n_level ...
,k_p_r_max ...
);

str_thisfunction = 'get_template_tree_1';

if nargin<1;
disp(sprintf(' %% testing %s',str_thisfunction));
parameter = struct('type','parameter');
parameter.flag_verbose = 1;
viewing_k_eq_d_min = 1.0; %<-- typical viewing_k_eq_d_min for coarse reconstruction. ;
n_level = 5;
k_p_r_max = 48/(2*pi); %<-- typical k_p_r_max for coarse reconstruction. ;
[ ...
 ~ ...
,template_tree ...
] = ...
get_template_tree_1( ...
 parameter ...
,viewing_k_eq_d_min ...
,n_level ...
,k_p_r_max ...
);
disp('returning'); return;
end;%if nargin<1;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); viewing_k_eq_d_min=[]; end; na=na+1;
if (nargin<1+na); n_level=[]; end; na=na+1;
if (nargin<1+na); k_p_r_max=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
flag_verbose=parameter.flag_verbose;

if isempty(viewing_k_eq_d_min); viewing_k_eq_d_min = 1.0/(2*pi); end;
if isempty(n_level); n_level = 5; end;
if isempty(k_p_r_max); k_p_r_max = 1.0; end;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

template_tree = struct('type','template_tree');
template_tree.k_p_r_max = k_p_r_max;
%%%%%%%%;
% build a hierarchy of templates, with templates defined on the finest level. ;
%%%%%%%%;
template_tree.n_level = n_level;
template_tree.viewing_k_eq_d_max = max(2*pi,viewing_k_eq_d_min);
template_tree.viewing_k_eq_d_min = viewing_k_eq_d_min;
template_tree.viewing_k_eq_d_ = zeros(n_level,1);
template_tree.n_S_ = zeros(n_level,1);
template_tree.viewing_azimu_b_all__ = cell(n_level,1);
template_tree.viewing_polar_a_all__ = cell(n_level,1);
template_tree.viewing_pole_k_c_0_all__ = cell(n_level,1);
template_tree.viewing_pole_k_c_1_all__ = cell(n_level,1);
template_tree.viewing_pole_k_c_2_all__ = cell(n_level,1);
template_tree.viewing_pole_k_c_Sd___ = cell(n_level,1);
%%%%%%%%;
% Initialize template_tree. ;
%%%%%%%%;
template_tree.viewing_k_eq_d_ = transpose(exp(linspace(log(template_tree.viewing_k_eq_d_max),log(template_tree.viewing_k_eq_d_min),n_level)));
for nlevel=0:n_level-1;
[ ...
 template_tree.n_S_(1+nlevel) ...
,template_tree.viewing_azimu_b_all__{1+nlevel} ...
,template_tree.viewing_polar_a_all__{1+nlevel} ...
,~ ...
,template_tree.viewing_pole_k_c_0_all__{1+nlevel} ...
,template_tree.viewing_pole_k_c_1_all__{1+nlevel} ...
,template_tree.viewing_pole_k_c_2_all__{1+nlevel} ...
, ...
] = ...
sample_shell_5( ...
 template_tree.k_p_r_max ...
,template_tree.viewing_k_eq_d_(1+nlevel) ...
,'L' ...
) ; %<-- obtain viewing angles on outer shell. ;
template_tree.viewing_pole_k_c_Sd___{1+nlevel} = ...
[ ...
 template_tree.viewing_pole_k_c_0_all__{1+nlevel} ...
,template_tree.viewing_pole_k_c_1_all__{1+nlevel} ...
,template_tree.viewing_pole_k_c_2_all__{1+nlevel} ...
];
end;%for nlevel=0:n_level-1;
%%%%%%%%;
template_tree.n_S_max = template_tree.n_S_(n_level);
n_S_max = template_tree.n_S_max;
template_tree.shift_z__ = get_template_inplane_shift_0(n_S_max,template_tree.viewing_polar_a_all__{n_level},template_tree.viewing_azimu_b_all__{n_level});
%%%%%%%%;
% Find representatives at finest level. ;
%%%%%%%%;
template_tree.representative_of_level_ = -ones(n_S_max,1);
template_tree.index_representative_fine_from_coarse__ = cell(n_level,1);
template_tree.index_representative_coarse_from_fine__ = cell(n_level,1);
tmp_index_representative_fine_remaining_ = 0:template_tree.n_S_(n_level)-1;
for nlevel=0:n_level-1;
n_S = template_tree.n_S_(1+nlevel);
[tmp_index_] = ...
knnsearch( ...
 template_tree.viewing_pole_k_c_Sd___{n_level}(1+tmp_index_representative_fine_remaining_,:) ...
,template_tree.viewing_pole_k_c_Sd___{1+nlevel} ...
);
tmp_index_ = tmp_index_ - 1;
if (numel(tmp_index_)<=0); error(sprintf(' %% Error: no representatives remaining at nlevel %d/%d, reduce n_level in get_template_tree_0',nlevel,n_level)); end;
tmp_index_representative_fine_ = tmp_index_representative_fine_remaining_(1+tmp_index_);
template_tree.index_representative_fine_from_coarse__{1+nlevel} = tmp_index_representative_fine_;
assert(numel(find(template_tree.representative_of_level_(1+tmp_index_representative_fine_)==-1))==numel(tmp_index_representative_fine_));
template_tree.representative_of_level_(1+tmp_index_representative_fine_) = nlevel;
tmp_index_representative_coarse_from_fine_ = -ones(n_S_max,1);
tmp_index_representative_coarse_from_fine_(1+tmp_index_representative_fine_) = 0:n_S-1;
template_tree.index_representative_coarse_from_fine__{1+nlevel} = tmp_index_representative_coarse_from_fine_;
tmp_index_representative_fine_remaining_ = setdiff(tmp_index_representative_fine_remaining_,tmp_index_representative_fine_);
end;%for nlevel=0:n_level-1;
%%%%%%%%;
% Find neighbors from one level to the next. ;
%%%%%%%%;
template_tree.index_neighbor_fine_from_coarse___ = cell(n_level,1);
template_tree.index_neighbor_coarse_from_fine__ = cell(n_level,1);
for nlevel=0:n_level-2;
n_S = template_tree.n_S_(1+nlevel);
[tmp_index_neighbor_coarse_from_fine_] = ...
knnsearch( ...
 template_tree.viewing_pole_k_c_Sd___{1+nlevel+0} ...
,template_tree.viewing_pole_k_c_Sd___{1+nlevel+1} ...
);
tmp_index_neighbor_coarse_from_fine_ = tmp_index_neighbor_coarse_from_fine_ - 1;
tmp_index_neighbor_fine_from_coarse__ = cell(n_S,1);
for nS=0:n_S-1;
tmp_index_neighbor_fine_from_coarse__{1+nS} = efind(tmp_index_neighbor_coarse_from_fine_==nS);
end;%for nS=0:n_S-1;
template_tree.index_neighbor_fine_from_coarse___{1+nlevel} = tmp_index_neighbor_fine_from_coarse__;
template_tree.index_neighbor_coarse_from_fine__{1+nlevel} = tmp_index_neighbor_coarse_from_fine_;
end;%for nlevel=0:n_level-2;
%%%%%%%%;
% Build neighbor representative list. ;
%%%%%%%%;
template_tree.index_neighbor_representative_fine_from_coarse___ = cell(n_level,1);
template_tree.index_neighbor_representative_coarse_from_fine___ = cell(n_level,1);
for nlevel=0:n_level-2;
n_S = template_tree.n_S_(1+nlevel);
tmp_index_neighbor_representative_fine_from_coarse__ = cell(n_S,1);
tmp_index_neighbor_representative_coarse_from_fine__ = cell(n_S,1);
tmp_index_neighbor_fine_from_coarse__ = template_tree.index_neighbor_fine_from_coarse___{1+nlevel};
tmp_index_representative_fine_from_coarse_ = template_tree.index_representative_fine_from_coarse__{1+nlevel+1};
for nS=0:n_S-1;
tmp_index_neighbor_representative_fine_from_coarse__{1+nS} = tmp_index_representative_fine_from_coarse_(1+tmp_index_neighbor_fine_from_coarse__{1+nS});
tmp_index_neighbor_representative_coarse_from_fine__{1+nS} = -ones(n_S_max,1);
tmp_index_neighbor_representative_coarse_from_fine__{1+nS}(1+tmp_index_neighbor_representative_fine_from_coarse__{1+nS}) = nS;
end;%for nS=0:n_S-1;
template_tree.index_neighbor_representative_fine_from_coarse___{1+nlevel} = tmp_index_neighbor_representative_fine_from_coarse__;
template_tree.index_neighbor_representative_coarse_from_fine___{1+nlevel} = tmp_index_neighbor_representative_coarse_from_fine__;
end;%for nlevel=0:n_level-2;
%%%%%%%%;

%%%%%%%%;
% check index values. ;
%%%%%%%%;
for nlevel=0:n_level-2;
n_S = template_tree.n_S_(1+nlevel);
tmp_index_representative_fine_from_coarse_ = template_tree.index_representative_fine_from_coarse__{1+nlevel};
tmp_index_representative_coarse_from_fine_ = template_tree.index_representative_coarse_from_fine__{1+nlevel};
assert(numel(union(tmp_index_representative_coarse_from_fine_(1+tmp_index_representative_fine_from_coarse_(1+(0:n_S-1))),transpose(0:n_S-1)))<=n_S);
assert(numel(find(template_tree.representative_of_level_(1+tmp_index_representative_fine_from_coarse_(1+(0:n_S-1)))==nlevel))==n_S);
end;%for nlevel=0:n_level-2;
%%%%%%%%;
for nlevel=0:n_level-2;
n_S = template_tree.n_S_(1+nlevel);
tmp_index_neighbor_coarse_from_fine_ = template_tree.index_neighbor_coarse_from_fine__{1+nlevel};
tmp_index_neighbor_fine_from_coarse__ = template_tree.index_neighbor_fine_from_coarse___{1+nlevel};
for nS=0:n_S-1;
assert(numel(find(tmp_index_neighbor_coarse_from_fine_(1+tmp_index_neighbor_fine_from_coarse__{1+nS})==nS))==numel(tmp_index_neighbor_fine_from_coarse__{1+nS}));
end;%for nS=0:n_S-1;
end;%for nlevel=0:n_level-2;
%%%%%%%%;
for nlevel=0:n_level-2;
n_S = template_tree.n_S_(1+nlevel);
tmp_index_neighbor_representative_fine_from_coarse__ = template_tree.index_neighbor_representative_fine_from_coarse___{1+nlevel};
tmp_index_neighbor_representative_coarse_from_fine__ = template_tree.index_neighbor_representative_coarse_from_fine___{1+nlevel};
for nS=0:n_S-1;
assert(numel(find(tmp_index_neighbor_representative_coarse_from_fine__{1+nS}(1+tmp_index_neighbor_representative_fine_from_coarse__{1+nS})==nS))==numel(tmp_index_neighbor_representative_fine_from_coarse__{1+nS}));
end;%for nS=0:n_S-1;
end;%for nlevel=0:n_level-2;
%%%%%%%%;

%%%%%%%%;
% For each node on each level, find the 8 nearest neighbors. ;
%%%%%%%%;
%{
template_tree.index_neighbor8_coarse_from_coarse___ = cell(n_level,1);
for nlevel=0:n_level-1;
n_S = template_tree.n_S_(1+nlevel);
[tmp_index_neighbor8_coarse_from_coarse__] = ...
knnsearch( ...
 template_tree.viewing_pole_k_c_Sd___{1+nlevel+0} ...
,template_tree.viewing_pole_k_c_Sd___{1+nlevel+0} ...
,'K' ...
,8 ...
);
tmp_index_neighbor8_coarse_from_coarse__ = tmp_index_neighbor8_coarse_from_coarse__ - 1;
template_tree.index_neighbor8_coarse_from_coarse___{1+nlevel} = tmp_index_neighbor8_coarse_from_coarse__;
end;%for nlevel=0:n_level-1;
%}
%%%%%%%%;
% For each node on each level, find the maximum radius of the neighbor-list (at the next-finer level). ;
%%%%%%%%;
template_tree.index_neighbor_fine_from_coarse_avg___ = cell(n_level,1);
template_tree.index_neighbor_fine_from_coarse_rad__ = cell(n_level,1);
for nlevel=0:n_level-2;
n_S = template_tree.n_S_(1+nlevel);
tmp_index_neighbor_fine_from_coarse__ = template_tree.index_neighbor_fine_from_coarse___{1+nlevel};
tmp_viewing_pole_k_c_avg__ = zeros(n_S,3);
tmp_viewing_pole_k_c_rad_ = zeros(n_S,1);
for nS=0:n_S-1;
tmp_index_neighbor_fine_from_coarse_ = tmp_index_neighbor_fine_from_coarse__{1+nS};
tmp_viewing_pole_k_c_Sd__ = template_tree.viewing_pole_k_c_Sd___{1+nlevel+1}(1+tmp_index_neighbor_fine_from_coarse_,:);
tmp_viewing_pole_k_c_avg_d_ = mean(tmp_viewing_pole_k_c_Sd__,1);
tmp_viewing_pole_k_c_avg__(1+nS,:) = tmp_viewing_pole_k_c_avg_d_;
tmp_viewing_pole_k_c_cnt_Sd__ = bsxfun(@minus,tmp_viewing_pole_k_c_Sd__,tmp_viewing_pole_k_c_avg_d_);
tmp_viewing_pole_k_c_rad = max(sqrt(sum(tmp_viewing_pole_k_c_cnt_Sd__.^2,2)));
tmp_viewing_pole_k_c_rad_(1+nS) = tmp_viewing_pole_k_c_rad;
end;%for nS=0:n_S-1;
template_tree.index_neighbor_fine_from_coarse_avg___{1+nlevel} = tmp_viewing_pole_k_c_avg__;
template_tree.index_neighbor_fine_from_coarse_rad__{1+nlevel} = tmp_viewing_pole_k_c_rad_;
end;%for nlevel=0:n_level-2;
%%%%%%%%;
% For each node on each level, find a list of elements (at the next-finer level) within 2 radii. ;
%%%%%%%%;
template_tree.index_neighbor_fine_from_coarse_exp___ = cell(n_level,1);
for nlevel=0:n_level-2;
n_S = template_tree.n_S_(1+nlevel);
tmp_index_neighbor_fine_from_coarse_exp__ = cell(n_S,1);
for nS=0:n_S-1;
tmp_index_neighbor_fine_from_coarse_avg_d_ = template_tree.index_neighbor_fine_from_coarse_avg___{1+nlevel}(1+nS,:);
tmp_index_neighbor_fine_from_coarse_rad = template_tree.index_neighbor_fine_from_coarse_rad__{1+nlevel}(1+nS);
tmp_viewing_pole_k_c_Sd__ = template_tree.viewing_pole_k_c_Sd___{1+nlevel+1};
tmp_viewing_pole_k_c_cnt_Sd__ = bsxfun(@minus,tmp_viewing_pole_k_c_Sd__,tmp_index_neighbor_fine_from_coarse_avg_d_);
tmp_viewing_pole_k_c_rad_S_ = sqrt(sum(tmp_viewing_pole_k_c_cnt_Sd__.^2,2));
tmp_index_neighbor_fine_from_coarse_exp_ = efind(tmp_viewing_pole_k_c_rad_S_<=2*tmp_index_neighbor_fine_from_coarse_rad);
tmp_index_neighbor_fine_from_coarse_exp__{1+nS} = tmp_index_neighbor_fine_from_coarse_exp_;
end;%for nS=0:n_S-1;
template_tree.index_neighbor_fine_from_coarse_exp___{1+nlevel} = tmp_index_neighbor_fine_from_coarse_exp__;
end;%for nlevel=0:n_level-2;
%%%%%%%%;

%%%%%%%%;
% Visualize: ;
%%%%%%%%;
nf=0;
flag_plot=flag_verbose;
if flag_plot;
figure(1+nf);nf=nf+1;clf;figsml;
imagesc_polar_a_azimu_b_0(template_tree.viewing_polar_a_all__{n_level},template_tree.viewing_azimu_b_all__{n_level},template_tree.representative_of_level_,[0,n_level-1],colormap('lines'),1);
xlim([0,2*pi]);ylim([0,1*pi]); axisnotick;
end;%if flag_plot;
%%%%%%%%;
% Visualize second-to-finest-level: ;
%%%%%%%%;
flag_plot=flag_verbose;
if flag_plot;
figure(1+nf);nf=nf+1;clf;hold on; figsml;
tmp_nlevel = n_level-2;
tmp_n_S = template_tree.n_S_(1+tmp_nlevel);
c_beach__ = colormap_beach(); n_c_beach = size(c_beach__,1);
str_symbol_ = {'o','^','s','p','h'}; markersize_use = 25;
hold on;
plot3( ...
 template_tree.viewing_pole_k_c_0_all__{n_level}(:) ...
,template_tree.viewing_pole_k_c_1_all__{n_level}(:) ...
,template_tree.viewing_pole_k_c_2_all__{n_level}(:) ...
,'.' ...
,'Color' , 'k' ...
);
for tmp_nS=0:tmp_n_S-1;
str_symbol = str_symbol_{1+mod(tmp_nS,5)};
nc = max(0,min(n_c_beach-1,floor(n_c_beach*tmp_nS/tmp_n_S)));
tmp_index_ = template_tree.index_neighbor_fine_from_coarse___{1+tmp_nlevel}{1+tmp_nS};
plot3( ...
 template_tree.viewing_pole_k_c_0_all__{n_level}(1+tmp_index_) ...
,template_tree.viewing_pole_k_c_1_all__{n_level}(1+tmp_index_) ...
,template_tree.viewing_pole_k_c_2_all__{n_level}(1+tmp_index_) ...
,str_symbol ...
,'Color' , c_beach__(1+nc,:) ...
,'MarkerSize', markersize_use ...
,'MarkerFaceColor', c_beach__(1+nc,:) ...
,'MarkerEdgeColor', 'k' ...
);
end;%for tmp_nS=0:tmp_n_S-1;
hold off;
axis equal; axis vis3d;
end;%if flag_plot;
%%%%%%%%;
% Visualize expanded neighborhoods at second-to-finest-level: ;
%%%%%%%%;
flag_plot=flag_verbose;
if flag_plot;
figure(1+nf);nf=nf+1;clf;hold on; figsml;
tmp_nlevel = n_level-2;
tmp_n_S = template_tree.n_S_(1+tmp_nlevel);
c_beach__ = colormap_beach(); n_c_beach = size(c_beach__,1);
str_symbol_ = {'o','^','s','p','h'}; markersize_big = 25; markersize_sml = 15;
tmp_nS0 = floor(0.3*tmp_n_S/2);
str_symbol = str_symbol_{1+mod(tmp_nS0,5)};
nc = max(0,min(n_c_beach-1,floor(n_c_beach*tmp_nS0/tmp_n_S)));
hold on;
plot3( ...
 template_tree.viewing_pole_k_c_0_all__{n_level}(:) ...
,template_tree.viewing_pole_k_c_1_all__{n_level}(:) ...
,template_tree.viewing_pole_k_c_2_all__{n_level}(:) ...
,'.' ...
,'Color' , 'k' ...
);
tmp_index_ = ...
setdiff( ...
 template_tree.index_neighbor_fine_from_coarse_exp___{1+tmp_nlevel}{1+tmp_nS0} ...
,template_tree.index_neighbor_fine_from_coarse___{1+tmp_nlevel}{1+tmp_nS0} ...
);
plot3( ...
 template_tree.viewing_pole_k_c_0_all__{n_level}(1+tmp_index_) ...
,template_tree.viewing_pole_k_c_1_all__{n_level}(1+tmp_index_) ...
,template_tree.viewing_pole_k_c_2_all__{n_level}(1+tmp_index_) ...
,str_symbol ...
,'Color' , c_beach__(1+nc,:) ...
,'MarkerSize', markersize_sml ...
,'MarkerFaceColor', 0.5*c_beach__(1+nc,:) ...
,'MarkerEdgeColor', 'k' ...
);
tmp_index_ = template_tree.index_neighbor_fine_from_coarse___{1+tmp_nlevel}{1+tmp_nS0};
plot3( ...
 template_tree.viewing_pole_k_c_0_all__{n_level}(1+tmp_index_) ...
,template_tree.viewing_pole_k_c_1_all__{n_level}(1+tmp_index_) ...
,template_tree.viewing_pole_k_c_2_all__{n_level}(1+tmp_index_) ...
,str_symbol ...
,'Color' , c_beach__(1+nc,:) ...
,'MarkerSize', markersize_big ...
,'MarkerFaceColor', c_beach__(1+nc,:) ...
,'MarkerEdgeColor', 'k' ...
);
hold off;
axis equal; axis vis3d;
end;%if flag_plot;
%%%%%%%%;

%%%%%%%%;
% Now, at each level, compose a list of nodes which have an azimuthal-angle within pi/sqrt(3) + [0,pi). ;
% These irrational endpoints are chosen so that it is unlikely that any points will lie on the boundaries. ;
% Each such 'pnode' will represent a pair of antipodal points within the projective sphere. ;
%%%%%%%%;
template_tree.n_P_ = zeros(n_level,1);
template_tree.index_nS_from_nP__ = cell(n_level,1);
for nlevel=0:n_level-1;
tmp_viewing_azimu_b_all_ = periodize(template_tree.viewing_azimu_b_all__{1+nlevel},0,2*pi);
index_nS_from_nP_ = efind( (tmp_viewing_azimu_b_all_>= pi/sqrt(3)) & (tmp_viewing_azimu_b_all_< pi/sqrt(3) + pi ));
n_P = numel(index_nS_from_nP_);
template_tree.index_nS_from_nP__{1+nlevel} = index_nS_from_nP_;
template_tree.n_P_(1+nlevel) = n_P;
end;%for nlevel=0:n_level-1;
%%%%%%%%;
% Now, for each pnode find the neighbor nodes (not pnodes) at the next level. ;
%%%%%%%%;
template_tree.index_neighbor_fine_from_pcoarse___ = cell(n_level,1);
template_tree.index_neighbor_pcoarse_from_fine__ = cell(n_level,1);
for nlevel=0:n_level-2;
n_P = template_tree.n_P_(1+nlevel);
index_nS_from_nP_ = template_tree.index_nS_from_nP__{1+nlevel};
tmp_viewing_pole_k_c_Pd__ = template_tree.viewing_pole_k_c_Sd___{1+nlevel+0}(1+index_nS_from_nP_,:);
tmp_viewing_pole_k_c_Sd__ = template_tree.viewing_pole_k_c_Sd___{1+nlevel+1};
[tmp_index_neighbor_pcoarse_from_fine_] = ...
knnsearch( ...
[+tmp_viewing_pole_k_c_Pd__ ; -tmp_viewing_pole_k_c_Pd__] ...
,tmp_viewing_pole_k_c_Sd__ ...
);
tmp_index_neighbor_pcoarse_from_fine_ = tmp_index_neighbor_pcoarse_from_fine_ - 1;
tmp_index_neighbor_pcoarse_from_fine_ = periodize(tmp_index_neighbor_pcoarse_from_fine_,0,n_P);
tmp_index_neighbor_fine_from_pcoarse__ = cell(n_P,1);
for nP=0:n_P-1;
tmp_index_neighbor_fine_from_pcoarse__{1+nP} = efind(tmp_index_neighbor_pcoarse_from_fine_==nP);
end;%for nP=0:n_P-1;
template_tree.index_neighbor_fine_from_pcoarse___{1+nlevel} = tmp_index_neighbor_fine_from_pcoarse__;
template_tree.index_neighbor_pcoarse_from_fine__{1+nlevel} = tmp_index_neighbor_pcoarse_from_fine_;
end;%for nlevel=0:n_level-2;
%%%%%%%%;
% For each pnode on each level, find a list of nodes (not pnodes) at the next-finer level within 2 radii. ;
%%%%%%%%;
template_tree.index_neighbor_fine_from_pcoarse_exp___ = cell(n_level,1);
for nlevel=0:n_level-2;
n_P = template_tree.n_P_(1+nlevel);
index_nS_from_nP_ = template_tree.index_nS_from_nP__{1+nlevel};
tmp_index_neighbor_fine_from_pcoarse_exp__ = cell(n_P,1);
for nP=0:n_P-1;
nS = index_nS_from_nP_(1+nP);
tmp_index_neighbor_fine_from_coarse_avg_d_ = template_tree.index_neighbor_fine_from_coarse_avg___{1+nlevel}(1+nS,:);
tmp_index_neighbor_fine_from_coarse_rad = template_tree.index_neighbor_fine_from_coarse_rad__{1+nlevel}(1+nS);
tmp_viewing_pole_k_c_Sd__ = template_tree.viewing_pole_k_c_Sd___{1+nlevel+1};
tmp_index_neighbor_fine_from_pcoarse_exp_ = [];
tmp_viewing_pole_k_c_cnt_Sd__ = bsxfun(@minus,tmp_viewing_pole_k_c_Sd__,+tmp_index_neighbor_fine_from_coarse_avg_d_);
tmp_viewing_pole_k_c_rad_S_ = sqrt(sum(tmp_viewing_pole_k_c_cnt_Sd__.^2,2));
tmp_index_neighbor_fine_from_pcoarse_exp_ = union(tmp_index_neighbor_fine_from_pcoarse_exp_,efind(tmp_viewing_pole_k_c_rad_S_<=2*tmp_index_neighbor_fine_from_coarse_rad));
tmp_viewing_pole_k_c_cnt_Sd__ = bsxfun(@minus,tmp_viewing_pole_k_c_Sd__,-tmp_index_neighbor_fine_from_coarse_avg_d_);
tmp_viewing_pole_k_c_rad_S_ = sqrt(sum(tmp_viewing_pole_k_c_cnt_Sd__.^2,2));
tmp_index_neighbor_fine_from_pcoarse_exp_ = union(tmp_index_neighbor_fine_from_pcoarse_exp_,efind(tmp_viewing_pole_k_c_rad_S_<=2*tmp_index_neighbor_fine_from_coarse_rad));
tmp_index_neighbor_fine_from_pcoarse_exp__{1+nP} = tmp_index_neighbor_fine_from_pcoarse_exp_;
end;%for nP=0:n_P-1;
template_tree.index_neighbor_fine_from_pcoarse_exp___{1+nlevel} = tmp_index_neighbor_fine_from_pcoarse_exp__;
end;%for nlevel=0:n_level-2;
%%%%%%%%;

%%%%%%%%;
% Visualize second-to-finest-level: ;
%%%%%%%%;
flag_plot=flag_verbose;
if flag_plot;
figure(1+nf);nf=nf+1;clf;hold on; figsml;
tmp_nlevel = n_level-2;
tmp_n_P = template_tree.n_P_(1+tmp_nlevel);
c_beach__ = colormap_beach(); n_c_beach = size(c_beach__,1);
str_symbol_ = {'o','^','s','p','h'}; markersize_use = 25;
hold on;
plot3( ...
 template_tree.viewing_pole_k_c_0_all__{n_level}(:) ...
,template_tree.viewing_pole_k_c_1_all__{n_level}(:) ...
,template_tree.viewing_pole_k_c_2_all__{n_level}(:) ...
,'.' ...
,'Color' , 'k' ...
);
for tmp_nP=0:tmp_n_P-1;
str_symbol = str_symbol_{1+mod(tmp_nP,5)};
nc = max(0,min(n_c_beach-1,floor(n_c_beach*tmp_nP/tmp_n_P)));
tmp_index_ = template_tree.index_neighbor_fine_from_pcoarse___{1+tmp_nlevel}{1+tmp_nP};
plot3( ...
 template_tree.viewing_pole_k_c_0_all__{n_level}(1+tmp_index_) ...
,template_tree.viewing_pole_k_c_1_all__{n_level}(1+tmp_index_) ...
,template_tree.viewing_pole_k_c_2_all__{n_level}(1+tmp_index_) ...
,str_symbol ...
,'Color' , c_beach__(1+nc,:) ...
,'MarkerSize', markersize_use ...
,'MarkerFaceColor', c_beach__(1+nc,:) ...
,'MarkerEdgeColor', 'k' ...
);
end;%for tmp_nP=0:tmp_n_P-1;
hold off;
axis equal; axis vis3d;
end;%if flag_plot;
%%%%%%%%;
% Visualize expanded neighborhoods at second-to-finest-level: ;
%%%%%%%%;
flag_plot=flag_verbose;
if flag_plot;
figure(1+nf);nf=nf+1;clf;hold on; figbig;
tmp_nlevel = n_level-2;
tmp_n_P = template_tree.n_P_(1+tmp_nlevel);
c_beach__ = colormap_beach(); n_c_beach = size(c_beach__,1);
str_symbol_ = {'o','^','s','p','h'}; markersize_big = 25; markersize_sml = 15;
for np=0:6-1;
subplot(2,3,1+np);
tmp_nP0 = max(0,min(tmp_n_P-1,floor(tmp_n_P*np/(6-1))));
str_symbol = str_symbol_{1+mod(tmp_nP0,5)};
nc = max(0,min(n_c_beach-1,floor(n_c_beach*tmp_nP0/tmp_n_P)));
hold on;
plot3( ...
 template_tree.viewing_pole_k_c_0_all__{n_level}(:) ...
,template_tree.viewing_pole_k_c_1_all__{n_level}(:) ...
,template_tree.viewing_pole_k_c_2_all__{n_level}(:) ...
,'.' ...
,'Color' , 'k' ...
);
tmp_index_ = ...
setdiff( ...
 template_tree.index_neighbor_fine_from_pcoarse_exp___{1+tmp_nlevel}{1+tmp_nP0} ...
,template_tree.index_neighbor_fine_from_pcoarse___{1+tmp_nlevel}{1+tmp_nP0} ...
);
plot3( ...
 template_tree.viewing_pole_k_c_0_all__{n_level}(1+tmp_index_) ...
,template_tree.viewing_pole_k_c_1_all__{n_level}(1+tmp_index_) ...
,template_tree.viewing_pole_k_c_2_all__{n_level}(1+tmp_index_) ...
,str_symbol ...
,'Color' , c_beach__(1+nc,:) ...
,'MarkerSize', markersize_sml ...
,'MarkerFaceColor', 0.5*c_beach__(1+nc,:) ...
,'MarkerEdgeColor', 'k' ...
);
tmp_index_ = template_tree.index_neighbor_fine_from_pcoarse___{1+tmp_nlevel}{1+tmp_nP0};
plot3( ...
 template_tree.viewing_pole_k_c_0_all__{n_level}(1+tmp_index_) ...
,template_tree.viewing_pole_k_c_1_all__{n_level}(1+tmp_index_) ...
,template_tree.viewing_pole_k_c_2_all__{n_level}(1+tmp_index_) ...
,str_symbol ...
,'Color' , c_beach__(1+nc,:) ...
,'MarkerSize', markersize_big ...
,'MarkerFaceColor', c_beach__(1+nc,:) ...
,'MarkerEdgeColor', 'k' ...
);
hold off;
axis equal; axis vis3d;
end;%for np=0:6-1;
end;%if flag_plot;
%%%%%%%%;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

