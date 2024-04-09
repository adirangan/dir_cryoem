function ...
[ ...
 scatter_from_tensor_sba__ ...
] = ...
shell_k_p_scatter_from_tensor_interpolate_n_5( ...
 n_order ...
,n_azimu_b ...
,n_polar_a ...
,polar_a_ ...
,n_scatter ...
,azimu_b_scatter_ ...
,polar_a_scatter_ ...
,flag_polar_a_ascend_vs_descend ...
);
% returns sparse matrix encoding the n_order interpolation operator between a tensor_grid ;
% using n_azimu_b azimu_b (unformly from 0 to 2*pi with periodic boundary). ;
% and n_polar_a polar_a values (not necessarily uniform from pi to 0) ;
% Note that (assuming flag_polar_a_ascend_vs_descend==0) polar_a_ descends. ;
%%%%%%%%;
% Note that a shell discretization (e.g., using sample_shell_6) ;
% with str_T_vs_L=='T' will generate uniformly spaced polar_a_, ;
% however, this discretization will not include the endpoints (i.e., 0 and pi). ;
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

str_thisfunction = 'shell_k_p_scatter_from_tensor_interpolate_n_5';

if (nargin<6);
rng(1);
%%%%%%%%;
flag_check=1;
if flag_check;
flag_disp = 1; nf=0;
disp(sprintf(' %% testing %s',str_thisfunction));
n_azimu_b = 128; n_polar_a = 65;
azimu_b_ = transpose(linspace(0,2*pi,n_azimu_b+1)); azimu_b_ = azimu_b_(1:end-1);
for flag_polar_a_ascend_vs_descend=[0,1];
disp(sprintf(' %% flag_polar_a_ascend_vs_descend %d:',flag_polar_a_ascend_vs_descend));
[ ...
 n_all ...
,azimu_b_all_ ...
,polar_a_all_ ...
,weight_all_ ...
,k_c_0_all_ ...
,k_c_1_all_ ...
,k_c_2_all_ ...
,n_polar_a ...
,polar_a_ ...
,n_azimu_b_ ...
] = ...
sample_shell_6( ...
 8.0 ...
,1.0 ...
,'L' ...
,1 ...
) ;
n_azimu_b = n_azimu_b_(1+0); azimu_b_ = transpose(linspace(0,2*pi,n_azimu_b+1)); azimu_b_ = azimu_b_(1:end-1);
if flag_polar_a_ascend_vs_descend==1; polar_a_ = flipud(polar_a_); end; %<-- traditional order. ;
if flag_polar_a_ascend_vs_descend==0; end; %<-- reverse order. ;
disp(sprintf(' %% n_azimu_b %d n_polar_a %d polar_a_: [%+0.2f,%+0.2f]',n_azimu_b,n_polar_a,polar_a_(1),polar_a_(end)));
%%%%%%%%;
n_grid = n_azimu_b*n_polar_a;
[azimu_b_ba__,polar_a_ba__] = ndgrid(azimu_b_,polar_a_);
n_scat = 1024;
polar_a_scat_ = 1*pi*rand(n_scat,1);
azimu_b_scat_ = 2*pi*rand(n_scat,1);
l_max = 4;
[Ylm_grid__] = get_Ylm__(1+l_max,0:l_max,n_grid,azimu_b_ba__,polar_a_ba__);
[Ylm_scat__] = get_Ylm__(1+l_max,0:l_max,n_scat,azimu_b_scat_,polar_a_scat_);
for n_order=[3,5,7];
scatter_from_tensor_sba__ = ...
shell_k_p_scatter_from_tensor_interpolate_n_5( ...
 n_order ...
,n_azimu_b ...
,n_polar_a ...
,polar_a_ ...
,n_scat ...
,azimu_b_scat_ ...
,polar_a_scat_ ...
,flag_polar_a_ascend_vs_descend ...
);
if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;
slim_ = max(abs(scatter_from_tensor_sba__),[],'all');
tmp_index_ = efind(scatter_from_tensor_sba__);
plot(sort(scatter_from_tensor_sba__(1+tmp_index_),'descend'),'.');
xlabel('nodes sorted by weight'); ylabel('weight');
title(sprintf('n_order %d',n_order));
end;%if flag_disp;
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

%%%%%%%%;
% flip polar_a_ if necessary. ;
%%%%%%%%;
sign_polar_a = +1;
if flag_polar_a_ascend_vs_descend==0;
sign_polar_a = -1; polar_a_scatter_ = pi-polar_a_scatter_; %<-- flip polar_a_scatter_ to account for reverse ordering of polar_a_. ;
end;%if flag_polar_a_ascend_vs_descend==0;

%%%%%%%%;
% transform polar_a_ to be uniform. ;
%%%%%%%%;
if flag_polar_a_ascend_vs_descend==1; polar_a_uni_ = transpose(0:+1:n_polar_a-1); end; %<-- uniform transformation. ;
if flag_polar_a_ascend_vs_descend==0; polar_a_uni_ = transpose(n_polar_a-1:-1:0); end; %<-- uniform transformation. ;
polar_a_scatter_ = interp1(polar_a_(:),polar_a_uni_(:),polar_a_scatter_,'spline','extrap'); %<-- use interpolation to uniformize polar_a_. ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% calculate scatter_from_tensor_sba__. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
spacing_azimu_b = 2*pi/n_azimu_b;
spacing_polar_a = 1.0;%spacing_polar_a = pi/(n_polar_a-1);
azimu_b_scatter_rescale_ = periodize(azimu_b_scatter_,0,2*pi)/spacing_azimu_b;
%polar_a_scatter_rescale_ = max(0,min(n_polar_a-1,periodize(polar_a_scatter_,0,1*pi)/spacing_polar_a));
polar_a_scatter_rescale_ = max(0-1,min(n_polar_a-1+1,polar_a_scatter_/spacing_polar_a));
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


