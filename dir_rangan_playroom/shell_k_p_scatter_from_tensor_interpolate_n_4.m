function ...
[ ...
 scatter_from_tensor_sba__ ...
,scatter_from_tensor_db1da0_sba__ ...
,scatter_from_tensor_db0da1_sba__ ...
,scatter_from_tensor_db1da1_sba__ ...
,scatter_from_tensor_db2da0_sba__ ...
,scatter_from_tensor_db0da2_sba__ ...
] = ...
shell_k_p_scatter_from_tensor_interpolate_n_4( ...
 n_order ...
,n_azimu_b ...
,n_polar_a ...
,n_scatter ...
,azimu_b_scatter_ ...
,polar_a_scatter_ ...
,flag_polar_a_ascend_vs_descend ...
);
% returns sparse matrix encoding the n_order interpolation operator between a tensor_grid ;
% using n_azimu_b azimu_b (unformly from 0 to 2*pi with periodic boundary). ;
% and n_polar_a polar_a (uniformly from pi to 0 inclusive of endpoints) ;
% Note that (assuming flag_polar_a_ascend_vs_descend==0) polar_a_ descends. ;
% Note also that a shell discretization (e.g., using sample_shell_6) ;
% with str_T_vs_L=='T' will generate uniformly spaced polar_a_, ;
% however, this discretization will not include the endpoints (i.e., 0 and pi). ;
% Thus, the interpolation here will not be accurate. ;
% To correct this, see shell_k_p_scatter_from_tensor_interpolate_n_5.m ;
%%%%;
% We presume these points are associated with a (n_azimu_b)-by-(n_polar_a) matrix ;
% a_k_p_ba__ with rows corresponding to azimu_b_ and columns corresponding to polar_a_. ;
%%%%;
% Note that this particular ordering is intended to align with sample_shell_6. ;
% i.e.: (assuming flag_polar_a_ascend_vs_descend==0) ;
% [lgnd_node_,lgnd_weight_] = legpts(n_polar_a);
% polar_a_ = acos(lgnd_node_);
%%%%;
% The n_scatter points to be interpolated have coordinates stored in ;
% azimu_b_scatter_ and polar_a_scatter_. ;
% We presume these points are associated with a n_scatter-element list a_k_all_. ;
% ;
% The (n_scatter)-by-(n_azimu_b*n_polar_a) scatter_from_tensor_sba__ matrix stores the ;
% interpolation weights linking a_k_all_ to the unrolled a_k_p_ba__(:). ;
%%%%;
% If requested, we also calculate the scatter_from_tensor_sba__ derivatives: ;
% scatter_from_tensor_db1da0_sba__; %<-- azimu_b 1-derivate. ;
% scatter_from_tensor_db0da1_sba__; %<-- polar_a 1-derivate. ;
% scatter_from_tensor_db1da1_sba__; %<-- azimu_b 1-derivate, polar_a 1-derivate. ;
% scatter_from_tensor_db2da0_sba__; %<-- azimu_b 2-derivates. ;
% scatter_from_tensor_db0da2_sba__; %<-- polar_a 2-derivates. ;
%%%%%%%%;

str_thisfunction = 'shell_k_p_scatter_from_tensor_interpolate_n_4';

if (nargin<6);
rng(1);
%%%%%%%%;
flag_check=1;
if flag_check;
disp(sprintf(' %% testing %s',str_thisfunction));
n_azimu_b = 128; n_polar_a = 65;
azimu_b_ = transpose(linspace(0,2*pi,n_azimu_b+1)); azimu_b_ = azimu_b_(1:end-1);
for flag_polar_a_ascend_vs_descend=[0,1];
disp(sprintf(' %% flag_polar_a_ascend_vs_descend %d:',flag_polar_a_ascend_vs_descend));
polar_a_ = transpose(linspace(0,1*pi,n_polar_a));
if flag_polar_a_ascend_vs_descend==1; end; %<-- traditional order. ;
if flag_polar_a_ascend_vs_descend==0; polar_a_ = flipud(polar_a_); end; %<-- reverse order. ;
n_grid = n_azimu_b*n_polar_a;
[azimu_b_ba__,polar_a_ba__] = ndgrid(azimu_b_,polar_a_);
n_scat = 1024;
polar_a_scat_ = 1*pi*rand(n_scat,1);
azimu_b_scat_ = 2*pi*rand(n_scat,1);
l_max = 4;
[Ylm_grid__] = get_Ylm__(1+l_max,0:l_max,n_grid,azimu_b_ba__,polar_a_ba__);
[Ylm_scat__] = get_Ylm__(1+l_max,0:l_max,n_scat,azimu_b_scat_,polar_a_scat_);
for n_order=[3,5,7];
scatter_from_tensor_sba__ = shell_k_p_scatter_from_tensor_interpolate_n_4(n_order,n_azimu_b,n_polar_a,n_scat,azimu_b_scat_,polar_a_scat_,flag_polar_a_ascend_vs_descend);
Ylm_pint__ = cell(1+l_max,1);
for l_val=0:l_max;
for m_val=-l_val:+l_val;
Ylm_pint__{1+l_val}(1+l_val+m_val,:) = scatter_from_tensor_sba__*transpose(Ylm_grid__{1+l_val}(1+l_val+m_val,:));
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max;
disp(sprintf(' %% interpolation relative error: n_order %d ',n_order));
for l_val=0:l_max;
for m_val=-l_val:+l_val;
disp(sprintf(' %% l_val %d m_val %+d : %0.16f',l_val,m_val,fnorm(Ylm_scat__{1+l_val}(1+l_val+m_val,:) - Ylm_pint__{1+l_val}(1+l_val+m_val,:))/fnorm(Ylm_scat__{1+l_val}(1+l_val+m_val,:))));
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max;
end;%for n_order=[3,5,7];
end;%for flag_polar_a_ascend_vs_descend=[0,1];
end;%if flag_check;
%%%%%%%%;
flag_check=1;
if flag_check;
disp(sprintf(' %% testing accuracy of scatter_from_tensor_dbndam_sba__:'));
n_azimu_b = 128; n_polar_a = 65;
azimu_b_ = transpose(linspace(0,2*pi,n_azimu_b+1)); azimu_b_ = azimu_b_(1:end-1);
for flag_polar_a_ascend_vs_descend=[0,1];
disp(sprintf(' %% flag_polar_a_ascend_vs_descend %d:',flag_polar_a_ascend_vs_descend));
polar_a_ = transpose(linspace(0,1*pi,n_polar_a));
if flag_polar_a_ascend_vs_descend==1; end; %<-- traditional order. ;
if flag_polar_a_ascend_vs_descend==0; polar_a_ = flipud(polar_a_); end; %<-- reverse order. ;
n_grid = n_azimu_b*n_polar_a;
[azimu_b_ba__,polar_a_ba__] = ndgrid(azimu_b_,polar_a_);
n_scat = 1024;
polar_a_scat_ = 1*pi*rand(n_scat,1);
azimu_b_scat_ = 2*pi*rand(n_scat,1);
f_db0da0 = @(azimu_b,polar_a) (+ 1*cos(3*azimu_b)) .* (+ 1*sin(4*polar_a)) ;
f_db1da0 = @(azimu_b,polar_a) (- 3*sin(3*azimu_b)) .* (+ 1*sin(4*polar_a)) ;
f_db0da1 = @(azimu_b,polar_a) (+ 1*cos(3*azimu_b)) .* (+ 4*cos(4*polar_a)) ;
f_db1da1 = @(azimu_b,polar_a) (- 3*sin(3*azimu_b)) .* (+ 4*cos(4*polar_a)) ;
f_db2da0 = @(azimu_b,polar_a) (- 9*cos(3*azimu_b)) .* (+ 1*sin(4*polar_a)) ;
f_db0da2 = @(azimu_b,polar_a) (+ 1*cos(3*azimu_b)) .* (-16*sin(4*polar_a)) ;
f_db0da0_grid__ = f_db0da0(azimu_b_ba__,polar_a_ba__);
f_db1da0_grid__ = f_db1da0(azimu_b_ba__,polar_a_ba__);
f_db0da1_grid__ = f_db0da1(azimu_b_ba__,polar_a_ba__);
f_db1da1_grid__ = f_db1da1(azimu_b_ba__,polar_a_ba__);
f_db2da0_grid__ = f_db2da0(azimu_b_ba__,polar_a_ba__);
f_db0da2_grid__ = f_db0da2(azimu_b_ba__,polar_a_ba__);
f_db0da0_scat_ = f_db0da0(azimu_b_scat_,polar_a_scat_);
f_db1da0_scat_ = f_db1da0(azimu_b_scat_,polar_a_scat_);
f_db0da1_scat_ = f_db0da1(azimu_b_scat_,polar_a_scat_);
f_db1da1_scat_ = f_db1da1(azimu_b_scat_,polar_a_scat_);
f_db2da0_scat_ = f_db2da0(azimu_b_scat_,polar_a_scat_);
f_db0da2_scat_ = f_db0da2(azimu_b_scat_,polar_a_scat_);
for n_order=[3,5,7];
[ ...
 scatter_from_tensor_sba__ ...
,scatter_from_tensor_db1da0_sba__ ...
,scatter_from_tensor_db0da1_sba__ ...
,scatter_from_tensor_db1da1_sba__ ...
,scatter_from_tensor_db2da0_sba__ ...
,scatter_from_tensor_db0da2_sba__ ...
] = ...
shell_k_p_scatter_from_tensor_interpolate_n_4( ...
 n_order ...
,n_azimu_b ...
,n_polar_a ...
,n_scat ...
,azimu_b_scat_ ...
,polar_a_scat_ ...
,flag_polar_a_ascend_vs_descend ...
);
f_db0da0_pint_ = scatter_from_tensor_sba__*f_db0da0_grid__(:);
f_db1da0_pint_ = scatter_from_tensor_db1da0_sba__*f_db0da0_grid__(:);
f_db0da1_pint_ = scatter_from_tensor_db0da1_sba__*f_db0da0_grid__(:);
f_db1da1_pint_ = scatter_from_tensor_db1da1_sba__*f_db0da0_grid__(:);
f_db2da0_pint_ = scatter_from_tensor_db2da0_sba__*f_db0da0_grid__(:);
f_db0da2_pint_ = scatter_from_tensor_db0da2_sba__*f_db0da0_grid__(:);
disp(sprintf(' %% n_order %d: f_db0da0_scat_ vs f_db0da0_pint_: %0.16f',n_order,fnorm(f_db0da0_scat_-f_db0da0_pint_)/fnorm(f_db0da0_scat_)));
disp(sprintf(' %% n_order %d: f_db1da0_scat_ vs f_db1da0_pint_: %0.16f',n_order,fnorm(f_db1da0_scat_-f_db1da0_pint_)/fnorm(f_db1da0_scat_)));
disp(sprintf(' %% n_order %d: f_db0da1_scat_ vs f_db0da1_pint_: %0.16f',n_order,fnorm(f_db0da1_scat_-f_db0da1_pint_)/fnorm(f_db0da1_scat_)));
disp(sprintf(' %% n_order %d: f_db1da1_scat_ vs f_db1da1_pint_: %0.16f',n_order,fnorm(f_db1da1_scat_-f_db1da1_pint_)/fnorm(f_db1da1_scat_)));
disp(sprintf(' %% n_order %d: f_db2da0_scat_ vs f_db2da0_pint_: %0.16f',n_order,fnorm(f_db2da0_scat_-f_db2da0_pint_)/fnorm(f_db2da0_scat_)));
disp(sprintf(' %% n_order %d: f_db0da2_scat_ vs f_db0da2_pint_: %0.16f',n_order,fnorm(f_db0da2_scat_-f_db0da2_pint_)/fnorm(f_db0da2_scat_)));
end;%for n_order=[3,5,7];
end;%for flag_polar_a_ascend_vs_descend=[0,1];
end;%if flag_check;
%%%%%%%%;
flag_check=0;
if flag_check;
disp(sprintf(' %% testing timing:'));
n_order = 5; 
n_azimu_b = 128; n_polar_a = 65;
n_scat = 1024*16; polar_a_scat_ = 1*pi*rand(n_scat,1); azimu_b_scat_ = 2*pi*rand(n_scat,1);
scatter_from_tensor_sba__ = shell_k_p_scatter_from_tensor_interpolate_n_4(n_order,n_azimu_b,n_polar_a,n_scat,azimu_b_scat_,polar_a_scat_);
tensor_from_scatter_bas__ = transpose(scatter_from_tensor_sba__);
tmp_a_k_p_grid_ = randn(n_azimu_b*n_polar_a,1) + i*randn(n_azimu_b*n_polar_a,1);
tmp_a_k_p_scat_ = randn(n_scat,1) + i*randn(n_scat,1);
n_iteration = 256;
tmp_t=tic;
for niteration=1:n_iteration; 
tmp_a_k_p_scat_out_ = scatter_from_tensor_sba__*tmp_a_k_p_grid_;
end;%for niteration=1:n_iteration; 
tmp_t=toc(tmp_t); 
disp(sprintf(' %% n_order %d n_azimu_b %d n_polar_a %d n_scatter %d n_iteration %d scatter_from_tensor_sba__ %0.2fs',n_order,n_azimu_b,n_polar_a,n_scat,n_iteration,tmp_t));
tmp_t=tic;
for niteration=1:n_iteration; 
tmp_a_k_p_grid_out_ = tensor_from_scatter_bas__*tmp_a_k_p_scat_;
end;%for niteration=1:n_iteration; 
tmp_t=toc(tmp_t); 
disp(sprintf(' %% n_order %d n_azimu_b %d n_polar_a %d n_scatter %d n_iteration %d: tensor_from_scatter_bas__ %0.2fs',n_order,n_azimu_b,n_polar_a,n_scat,n_iteration,tmp_t));
end;%if flag_check;
%%%%%%%%;
disp(sprintf(' %% returning')); return;
end;%if (nargin<6);

na=0;
if (nargin<1+na); n_order=[]; end; na=na+1;
if (nargin<1+na); n_azimu_b=[]; end; na=na+1;
if (nargin<1+na); n_polar_a=[]; end; na=na+1;
if (nargin<1+na); n_scatter=[]; end; na=na+1;
if (nargin<1+na); azimu_b_scatter_=[]; end; na=na+1;
if (nargin<1+na); polar_a_scatter_=[]; end; na=na+1;
if (nargin<1+na); flag_polar_a_ascend_vs_descend=[]; end; na=na+1;
if isempty(flag_polar_a_ascend_vs_descend); flag_polar_a_ascend_vs_descend = 0; end;

if (n_order>min(n_azimu_b,n_polar_a)); disp(sprintf(' %% Warning, n_order %d > n_azimu_b %d n_polar_a %d',n_order,n_azimu_b,n_polar_a)); end;
n_order = min(n_order,min(n_azimu_b,n_polar_a));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% calculate scatter_from_tensor_sba__. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
sign_polar_a = +1;
if flag_polar_a_ascend_vs_descend==0;
sign_polar_a = -1; polar_a_scatter_ = pi-polar_a_scatter_; %<-- flip polar_a_scatter_ to account for reverse ordering of polar_a_. ;
end;%if flag_polar_a_ascend_vs_descend==0;
%%%%%%%%;
spacing_azimu_b = 2*pi/n_azimu_b;
spacing_polar_a = pi/(n_polar_a-1);
azimu_b_scatter_rescale_ = periodize(azimu_b_scatter_,0,2*pi)/spacing_azimu_b;
polar_a_scatter_rescale_ = max(0,min(n_polar_a-1,periodize(polar_a_scatter_,0,1*pi)/spacing_polar_a));
n_half_order = floor(n_order/2);
azimu_b_scatter_floor_ = floor(azimu_b_scatter_rescale_)-n_half_order; %<-- do not loop around boundary just yet. ;
polar_a_scatter_floor_ = min(n_polar_a-n_order,max(0,floor(polar_a_scatter_rescale_)-n_half_order));
azimu_b_scatter_shift_ = azimu_b_scatter_rescale_ - azimu_b_scatter_floor_;
polar_a_scatter_shift_ = polar_a_scatter_rescale_ - polar_a_scatter_floor_;
%%%%%%%%;
node_x_ = transpose(0:n_order-1);
weight_denominator_ = prod( repmat(node_x_,1,n_order) - repmat(transpose(node_x_),n_order,1) + eye(n_order,n_order) , 2 ) ;
weight_numerator_polar_a__ = zeros(n_order,n_scatter);
weight_numerator_azimu_b__ = zeros(n_order,n_scatter);
index_col_azimu_b__ = zeros(n_order,n_scatter);
index_col_polar_a__ = zeros(n_order,n_scatter);
index_row__ = repmat(0:n_scatter-1,n_order^2,1);
for nscatter=0:n_scatter-1;
weight_numerator_azimu_b__(:,1+nscatter) = prod( azimu_b_scatter_shift_(1+nscatter) - repmat(node_x_,1,n_order) + diag(1-(azimu_b_scatter_shift_(1+nscatter)-node_x_)) , 1 ) ;
weight_numerator_polar_a__(:,1+nscatter) = prod( polar_a_scatter_shift_(1+nscatter) - repmat(node_x_,1,n_order) + diag(1-(polar_a_scatter_shift_(1+nscatter)-node_x_)) , 1);
end;%for nscatter=0:n_scatter-1;
%%%%%%%%;
index_col_azimu_b__ = bsxfun(@plus,reshape(azimu_b_scatter_floor_,[1,n_scatter]),transpose([0:n_order-1]));
index_col_polar_a__ = bsxfun(@plus,reshape(polar_a_scatter_floor_,[1,n_scatter]),transpose([0:n_order-1]));
%%%%%%%%;
% periodize azimu_b. ;
%%%%%%%%;
tmp_index_ = efind(index_col_azimu_b__>=n_azimu_b); 
index_col_azimu_b__(1+tmp_index_) = index_col_azimu_b__(1+tmp_index_) - n_azimu_b;
tmp_index_ = efind(index_col_azimu_b__< 0        ); 
index_col_azimu_b__(1+tmp_index_) = index_col_azimu_b__(1+tmp_index_) + n_azimu_b;
%%%%%%%%;
weight_azimu_b__ = diag(1./weight_denominator_)*weight_numerator_azimu_b__;
weight_polar_a__ = diag(1./weight_denominator_)*weight_numerator_polar_a__;
weight__ = zeros(n_order^2,n_scatter);
index_col__ = zeros(n_order^2,n_scatter);
for nscatter=0:n_scatter-1;
weight__(:,1+nscatter) = reshape(weight_azimu_b__(:,1+nscatter)*transpose(weight_polar_a__(:,1+nscatter)),[n_order^2,1]);
index_col__(:,1+nscatter) = reshape(bsxfun(@plus,index_col_azimu_b__(:,1+nscatter),transpose(index_col_polar_a__(:,1+nscatter))*n_azimu_b),[n_order^2,1]);
end;%for nscatter=0:n_scatter-1;
%%%%%%%%;
scatter_from_tensor_sba__ = sparse(1+index_row__(:),1+index_col__(:),weight__(:),n_scatter,n_azimu_b*n_polar_a);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if nargout> 1;
%%%%%%%%;
% calculate scatter_from_tensor_dbndam_sba__. ;
%%%%%%%%;
weight_numerator_db1da0_azimu_b__ = zeros(n_order,n_scatter);
weight_numerator_db0da1_polar_a__ = zeros(n_order,n_scatter);
%%%%;
for nscatter=0:n_scatter-1;
%%%%;
tmp_b = azimu_b_scatter_shift_(1+nscatter);
tmp_mat0__ = tmp_b - repmat(node_x_,1,n_order) + diag(1-(tmp_b-node_x_));
tmp_weight_numerator_db1da0_azimu_b_ = zeros(n_order,1);
for norder0=1:n_order-1;
tmp_mat1__ = circshift(diag(1-(tmp_b-circshift(node_x_,-norder0))),+norder0);
tmp_weight_numerator_db1da0_azimu_b_ = tmp_weight_numerator_db1da0_azimu_b_ + transpose(prod( tmp_mat0__ + tmp_mat1__ , 1 )) ;
end;%for norder0=1:n_order-1;
weight_numerator_db1da0_azimu_b__(:,1+nscatter) = tmp_weight_numerator_db1da0_azimu_b_;
clear tmp_b tmp_mat0__ tmp_mat1__ tmp_weight_numerator_db1da0_azimu_b_;
%%%%;
tmp_a = polar_a_scatter_shift_(1+nscatter);
tmp_mat0__ = tmp_a - repmat(node_x_,1,n_order) + diag(1-(tmp_a-node_x_));
tmp_weight_numerator_db0da1_polar_a_ = zeros(n_order,1);
for norder0=1:n_order-1;
tmp_mat1__ = circshift(diag(1-(tmp_a-circshift(node_x_,-norder0))),+norder0);
tmp_weight_numerator_db0da1_polar_a_ = tmp_weight_numerator_db0da1_polar_a_ + transpose(prod( tmp_mat0__ + tmp_mat1__ , 1 )) ;
end;%for norder0=1:n_order-1;
weight_numerator_db0da1_polar_a__(:,1+nscatter) = tmp_weight_numerator_db0da1_polar_a_;
clear tmp_a tmp_mat0__ tmp_mat1__ tmp_weight_numerator_db0da1_polar_a_;
%%%%;
end;%for nscatter=0:n_scatter-1;
weight_db1da0_azimu_b__ = diag(1./weight_denominator_)*weight_numerator_db1da0_azimu_b__/max(1e-12,spacing_azimu_b);
weight_db0da1_polar_a__ = diag(1./weight_denominator_)*weight_numerator_db0da1_polar_a__/max(1e-12,spacing_polar_a)*sign_polar_a;
%%%%;
weight_db1da0__ = zeros(n_order^2,n_scatter);
for nscatter=0:n_scatter-1;
weight_db1da0__(:,1+nscatter) = reshape(weight_db1da0_azimu_b__(:,1+nscatter)*transpose(weight_polar_a__(:,1+nscatter)),[n_order^2,1]);
end;%for nscatter=0:n_scatter-1;
scatter_from_tensor_db1da0_sba__ = sparse(1+index_row__(:),1+index_col__(:),weight_db1da0__(:),n_scatter,n_azimu_b*n_polar_a);
%%%%;
weight_db0da1__ = zeros(n_order^2,n_scatter);
for nscatter=0:n_scatter-1;
weight_db0da1__(:,1+nscatter) = reshape(weight_azimu_b__(:,1+nscatter)*transpose(weight_db0da1_polar_a__(:,1+nscatter)),[n_order^2,1]);
end;%for nscatter=0:n_scatter-1;
scatter_from_tensor_db0da1_sba__ = sparse(1+index_row__(:),1+index_col__(:),weight_db0da1__(:),n_scatter,n_azimu_b*n_polar_a);
%%%%;
weight_db1da1__ = zeros(n_order^2,n_scatter);
for nscatter=0:n_scatter-1;
weight_db1da1__(:,1+nscatter) = reshape(weight_db1da0_azimu_b__(:,1+nscatter)*transpose(weight_db0da1_polar_a__(:,1+nscatter)),[n_order^2,1]);
end;%for nscatter=0:n_scatter-1;
scatter_from_tensor_db1da1_sba__ = sparse(1+index_row__(:),1+index_col__(:),weight_db1da1__(:),n_scatter,n_azimu_b*n_polar_a);
%%%%%%%%;
weight_numerator_db2da0_azimu_b__ = zeros(n_order,n_scatter);
weight_numerator_db0da2_polar_a__ = zeros(n_order,n_scatter);
%%%%;
for nscatter=0:n_scatter-1;
%%%%;
tmp_b = azimu_b_scatter_shift_(1+nscatter);
tmp_mat0__ = tmp_b - repmat(node_x_,1,n_order) + diag(1-(tmp_b-node_x_));
tmp_weight_numerator_db2da0_azimu_b_ = zeros(n_order,1);
for norder0=1:n_order-1;
tmp_mat1__ = circshift(diag(1-(tmp_b-circshift(node_x_,-norder0))),+norder0);
for norder1=1+norder0:n_order-1;
tmp_mat2__ = circshift(diag(1-(tmp_b-circshift(node_x_,-norder1))),+norder1);
tmp_weight_numerator_db2da0_azimu_b_ = tmp_weight_numerator_db2da0_azimu_b_ + transpose(prod( tmp_mat0__ + tmp_mat1__ + tmp_mat2__, 1 )) ;
end;%for norder1=1+norder0:n_order-1;
end;%for norder0=1:n_order-1;
weight_numerator_db2da0_azimu_b__(:,1+nscatter) = tmp_weight_numerator_db2da0_azimu_b_;
clear tmp_b tmp_mat0__ tmp_mat1__ tmp_weight_numerator_db2da0_azimu_b_;
%%%%;
tmp_a = polar_a_scatter_shift_(1+nscatter);
tmp_mat0__ = tmp_a - repmat(node_x_,1,n_order) + diag(1-(tmp_a-node_x_));
tmp_weight_numerator_db0da2_polar_a_ = zeros(n_order,1);
for norder0=1:n_order-1;
tmp_mat1__ = circshift(diag(1-(tmp_a-circshift(node_x_,-norder0))),+norder0);
for norder1=1+norder0:n_order-1;
tmp_mat2__ = circshift(diag(1-(tmp_a-circshift(node_x_,-norder1))),+norder1);
tmp_weight_numerator_db0da2_polar_a_ = tmp_weight_numerator_db0da2_polar_a_ + transpose(prod( tmp_mat0__ + tmp_mat1__ + tmp_mat2__, 1 )) ;
end;%for norder1=1+norder0:n_order-1;
end;%for norder0=1:n_order-1;
weight_numerator_db0da2_polar_a__(:,1+nscatter) = tmp_weight_numerator_db0da2_polar_a_;
clear tmp_a tmp_mat0__ tmp_mat1__ tmp_weight_numerator_db0da2_polar_a_;
%%%%;
end;%for nscatter=0:n_scatter-1;
weight_db2da0_azimu_b__ = diag(1./weight_denominator_)*weight_numerator_db2da0_azimu_b__/max(1e-12,spacing_azimu_b)/max(1e-12,spacing_azimu_b)*2;
weight_db0da2_polar_a__ = diag(1./weight_denominator_)*weight_numerator_db0da2_polar_a__/max(1e-12,spacing_polar_a)/max(1e-12,spacing_polar_a)*2;
%%%%;
weight_db2da0__ = zeros(n_order^2,n_scatter);
for nscatter=0:n_scatter-1;
weight_db2da0__(:,1+nscatter) = reshape(weight_db2da0_azimu_b__(:,1+nscatter)*transpose(weight_polar_a__(:,1+nscatter)),[n_order^2,1]);
end;%for nscatter=0:n_scatter-1;
scatter_from_tensor_db2da0_sba__ = sparse(1+index_row__(:),1+index_col__(:),weight_db2da0__(:),n_scatter,n_azimu_b*n_polar_a);
%%%%;
weight_db0da2__ = zeros(n_order^2,n_scatter);
for nscatter=0:n_scatter-1;
weight_db0da2__(:,1+nscatter) = reshape(weight_azimu_b__(:,1+nscatter)*transpose(weight_db0da2_polar_a__(:,1+nscatter)),[n_order^2,1]);
end;%for nscatter=0:n_scatter-1;
scatter_from_tensor_db0da2_sba__ = sparse(1+index_row__(:),1+index_col__(:),weight_db0da2__(:),n_scatter,n_azimu_b*n_polar_a);
%%%%%%%%;
end;%if nargout> 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

