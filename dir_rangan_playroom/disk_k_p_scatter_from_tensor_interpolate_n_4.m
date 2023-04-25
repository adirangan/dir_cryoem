function ...
[ ...
 scatter_from_tensor_swk__ ...
,scatter_from_tensor_dw1dk0_swk__ ...
,scatter_from_tensor_dw0dk1_swk__ ...
,scatter_from_tensor_dw1dk1_swk__ ...
,scatter_from_tensor_dw2dk0_swk__ ...
,scatter_from_tensor_dw0dk2_swk__ ...
] = ...
disk_k_p_scatter_from_tensor_interpolate_n_4( ...
 n_order ...
,n_gamma_z ...
,n_k_p_rad ...
,k_p_rad_ ...
,k_p_rad_max ...
,n_scatter ...
,gamma_z_scatter_ ...
,k_p_rad_scatter_ ...
);
% returns sparse matrix encoding the n_order interpolation operator between a tensor_grid ;
% using n_gamma_z gamma_z (unformly from 0 to 2*pi with periodic boundary). ;
% and n_k_p_rad k_p_rad values stored in k_p_rad_. ;
%%%%;
% We presume these points are associated with a (n_gamma_z)-by-(n_k_p_rad) matrix ;
% S_k_p_wk__ with rows corresponding to gamma_z_ and columns corresponding to k_p_rad_. ;
%%%%;
% The n_scatter points to be interpolated have coordinates stored in ;
% gamma_z_scatter_ and k_p_rad_scatter_. ;
% We presume these points are associated with a n_scatter-element list S_k_all_. ;
% ;
% The (n_scatter)-by-(n_gamma_z*n_k_p_rad) scatter_from_tensor_swk__ matrix stores the ;
% interpolation weights linking S_k_all_ to the unrolled S_k_p_wk__(:). ;
%%%%;
% If requested, we also calculate the scatter_from_tensor_swk__ derivatives: ;
% scatter_from_tensor_dw1dk0_swk__; %<-- gamma_z 1-derivate. ;
% scatter_from_tensor_dw0dk1_swk__; %<-- k_p_rad 1-derivate. ;
% scatter_from_tensor_dw1dk1_swk__; %<-- gamma_z 1-derivate, k_p_rad 1-derivate. ;
% scatter_from_tensor_dw2dk0_swk__; %<-- gamma_z 2-derivates. ;
% scatter_from_tensor_dw0dk2_swk__; %<-- k_p_rad 2-derivates. ;
% Because we uniformize the k_p_rad_ below (and do not later account for this transformation), ;
% the k_p_rad derivatives will only be accurate if the original k_p_rad_ are uniform. ;
%%%%%%%%;

str_thisfunction = 'disk_k_p_scatter_from_tensor_interpolate_n_4';

if (nargin<1);
rng(1);
%%%%%%%%;
flag_check=1;
if flag_check;
disp(sprintf(' %% testing %s',str_thisfunction));
n_gamma_z = 128; n_k_p_rad = 65; n_order = 7;
for norder=3:2:n_order;
gamma_z_ = transpose(linspace(0,2*pi,n_gamma_z+1)); gamma_z_ = gamma_z_(1:end-1);
k_p_rad_max = 1.0;
k_p_rad_ = k_p_rad_max * transpose(0.5 + 0.5*cos(linspace(pi,0,n_k_p_rad))); %<-- nonuniform. ;
n_grid = n_gamma_z*n_k_p_rad;
[gamma_z_wk__,k_p_rad_wk__] = ndgrid(gamma_z_,k_p_rad_);
n_scat = 1024;
k_p_rad_scat_ = k_p_rad_max * rand(n_scat,1);
gamma_z_scat_ = 2*pi*rand(n_scat,1);
f_grid__ = zeros(n_gamma_z,n_k_p_rad);
f_scat_ = zeros(n_scat,1);
for nsource=0:8-1;
rng(nsource);
delta_ = 10*randn(2,1); delta_0 = delta_(1+0); delta_1 = delta_(1+1); %delta_0 = +15.0; delta_1 = -12.5;
phi = 2*pi*rand();
f = @(w,k) exp( i * k.*cos(w)*delta_0 + i * k.*sin(w)*delta_1 + i*phi);
f_grid__ = f_grid__ + f(gamma_z_wk__,k_p_rad_wk__);
f_scat_ = f_scat_ + f(gamma_z_scat_,k_p_rad_scat_);
end;%for nsource=0:8-1;
figure(1);clf;figsml;imagesc_p(n_k_p_rad,k_p_rad_,n_gamma_z*ones(n_k_p_rad,1),n_grid,real(f_grid__(:)));
scatter_from_tensor_swk__ = ...
disk_k_p_scatter_from_tensor_interpolate_n_4( ...
 norder ...
,n_gamma_z ...
,n_k_p_rad ...
,k_p_rad_ ...
,k_p_rad_max ...
,n_scat ...
,gamma_z_scat_ ...
,k_p_rad_scat_ ...
);
f_pint_ = scatter_from_tensor_swk__*f_grid__(:);
disp(sprintf(' %% norder %d: f_scat_ vs f_pint_: %0.16f',norder,fnorm(f_scat_-f_pint_)/fnorm(f_scat_)));
end;%for norder=3:2:n_order;
clear n_gamma_z n_k_p_rad n_order gamma_z_ ;
clear k_p_rad_max k_p_rad_ n_grid ;
clear gamma_z_wk__ k_p_rad_wk__ ;
clear n_scat k_p_rad_scat_ gamma_z_scat_ ;
clear delta_0 delta_1 f f_grid__ f_scat_ ;
clear scatter_from_tensor_swk__ f_pint_ ;
end;%if flag_check;
%%%%%%%%;
flag_check=1;
if flag_check;
disp(sprintf(' %% testing accuracy of scatter_from_tensor_dwndkm_swk__:'));
n_gamma_z = 128; n_k_p_rad = 65; n_order = 5;
gamma_z_ = transpose(linspace(0,2*pi,n_gamma_z+1)); gamma_z_ = gamma_z_(1:end-1);
k_p_rad_max = 1.0;
%k_p_rad_ = k_p_rad_max * transpose(0.5 + 0.5*cos(linspace(pi,0,n_k_p_rad))); %<-- nonuniform. ;
k_p_rad_ = k_p_rad_max * transpose(linspace(0,1,n_k_p_rad));
n_grid = n_gamma_z*n_k_p_rad;
[gamma_z_wk__,k_p_rad_wk__] = ndgrid(gamma_z_,k_p_rad_);
n_scat = 1024;
k_p_rad_scat_ = k_p_rad_max * rand(n_scat,1);
gamma_z_scat_ = 2*pi*rand(n_scat,1);
f_dw0dk0 = @(gamma_z,k_p_rad) (+ 1*cos(3*gamma_z)) .* (+ 1*sin(4*k_p_rad)) ;
f_dw1dk0 = @(gamma_z,k_p_rad) (- 3*sin(3*gamma_z)) .* (+ 1*sin(4*k_p_rad)) ;
f_dw0dk1 = @(gamma_z,k_p_rad) (+ 1*cos(3*gamma_z)) .* (+ 4*cos(4*k_p_rad)) ;
f_dw1dk1 = @(gamma_z,k_p_rad) (- 3*sin(3*gamma_z)) .* (+ 4*cos(4*k_p_rad)) ;
f_dw2dk0 = @(gamma_z,k_p_rad) (- 9*cos(3*gamma_z)) .* (+ 1*sin(4*k_p_rad)) ;
f_dw0dk2 = @(gamma_z,k_p_rad) (+ 1*cos(3*gamma_z)) .* (-16*sin(4*k_p_rad)) ;
f_dw0dk0_grid__ = f_dw0dk0(gamma_z_wk__,k_p_rad_wk__);
f_dw1dk0_grid__ = f_dw1dk0(gamma_z_wk__,k_p_rad_wk__);
f_dw0dk1_grid__ = f_dw0dk1(gamma_z_wk__,k_p_rad_wk__);
f_dw1dk1_grid__ = f_dw1dk1(gamma_z_wk__,k_p_rad_wk__);
f_dw2dk0_grid__ = f_dw2dk0(gamma_z_wk__,k_p_rad_wk__);
f_dw0dk2_grid__ = f_dw0dk2(gamma_z_wk__,k_p_rad_wk__);
f_dw0dk0_scat_ = f_dw0dk0(gamma_z_scat_,k_p_rad_scat_);
f_dw1dk0_scat_ = f_dw1dk0(gamma_z_scat_,k_p_rad_scat_);
f_dw0dk1_scat_ = f_dw0dk1(gamma_z_scat_,k_p_rad_scat_);
f_dw1dk1_scat_ = f_dw1dk1(gamma_z_scat_,k_p_rad_scat_);
f_dw2dk0_scat_ = f_dw2dk0(gamma_z_scat_,k_p_rad_scat_);
f_dw0dk2_scat_ = f_dw0dk2(gamma_z_scat_,k_p_rad_scat_);
for n_order=[3,5,7];
[ ...
 scatter_from_tensor_swk__ ...
,scatter_from_tensor_dw1dk0_swk__ ...
,scatter_from_tensor_dw0dk1_swk__ ...
,scatter_from_tensor_dw1dk1_swk__ ...
,scatter_from_tensor_dw2dk0_swk__ ...
,scatter_from_tensor_dw0dk2_swk__ ...
] = ...
disk_k_p_scatter_from_tensor_interpolate_n_4( ...
 n_order ...
,n_gamma_z ...
,n_k_p_rad ...
,k_p_rad_ ...
,k_p_rad_max ...
,n_scat ...
,gamma_z_scat_ ...
,k_p_rad_scat_ ...
);
f_dw0dk0_pint_ = scatter_from_tensor_swk__*f_dw0dk0_grid__(:);
f_dw1dk0_pint_ = scatter_from_tensor_dw1dk0_swk__*f_dw0dk0_grid__(:);
f_dw0dk1_pint_ = scatter_from_tensor_dw0dk1_swk__*f_dw0dk0_grid__(:);
f_dw1dk1_pint_ = scatter_from_tensor_dw1dk1_swk__*f_dw0dk0_grid__(:);
f_dw2dk0_pint_ = scatter_from_tensor_dw2dk0_swk__*f_dw0dk0_grid__(:);
f_dw0dk2_pint_ = scatter_from_tensor_dw0dk2_swk__*f_dw0dk0_grid__(:);
disp(sprintf(' %% n_order %d: f_dw0dk0_scat_ vs f_dw0dk0_pint_: %0.16f',n_order,fnorm(f_dw0dk0_scat_-f_dw0dk0_pint_)/fnorm(f_dw0dk0_scat_)));
disp(sprintf(' %% n_order %d: f_dw1dk0_scat_ vs f_dw1dk0_pint_: %0.16f',n_order,fnorm(f_dw1dk0_scat_-f_dw1dk0_pint_)/fnorm(f_dw1dk0_scat_)));
disp(sprintf(' %% n_order %d: f_dw0dk1_scat_ vs f_dw0dk1_pint_: %0.16f',n_order,fnorm(f_dw0dk1_scat_-f_dw0dk1_pint_)/fnorm(f_dw0dk1_scat_)));
disp(sprintf(' %% n_order %d: f_dw1dk1_scat_ vs f_dw1dk1_pint_: %0.16f',n_order,fnorm(f_dw1dk1_scat_-f_dw1dk1_pint_)/fnorm(f_dw1dk1_scat_)));
disp(sprintf(' %% n_order %d: f_dw2dk0_scat_ vs f_dw2dk0_pint_: %0.16f',n_order,fnorm(f_dw2dk0_scat_-f_dw2dk0_pint_)/fnorm(f_dw2dk0_scat_)));
disp(sprintf(' %% n_order %d: f_dw0dk2_scat_ vs f_dw0dk2_pint_: %0.16f',n_order,fnorm(f_dw0dk2_scat_-f_dw0dk2_pint_)/fnorm(f_dw0dk2_scat_)));
end;%for n_order=[3,5,7];
end;%if flag_check;
%%%%%%%%;
disp(sprintf(' %% returning')); return;
end;%if (nargin<1);

na=0;
if (nargin<1+na); n_order=[]; end; na=na+1;
if (nargin<1+na); n_gamma_z=[]; end; na=na+1;
if (nargin<1+na); n_k_p_rad=[]; end; na=na+1;
if (nargin<1+na); k_p_rad_=[]; end; na=na+1;
if (nargin<1+na); k_p_rad_max=[]; end; na=na+1;
if (nargin<1+na); n_scatter=[]; end; na=na+1;
if (nargin<1+na); gamma_z_scatter_=[]; end; na=na+1;
if (nargin<1+na); k_p_rad_scatter_=[]; end; na=na+1;

if (n_order>min(n_gamma_z,n_k_p_rad)); disp(sprintf(' %% Warning, n_order %d > n_gamma_z %d n_k_p_rad %d',n_order,n_gamma_z,n_k_p_rad)); end;
n_order = min(n_order,min(n_gamma_z,n_k_p_rad));

%%%%%%%%;
% transform k_p_rad_ to be uniform. ;
%%%%%%%%;
k_p_rad_uni_ = transpose(0:n_k_p_rad-1); %<-- uniform transformation. ;
k_p_rad_scatter_ = interp1(k_p_rad_(:),k_p_rad_uni_(:),k_p_rad_scatter_,'spline','extrap'); %<-- use interpolation to uniformize k_p_rad_. ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% calculate scatter_from_tensor_swk__. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
spacing_gamma_z = 2*pi/n_gamma_z;
spacing_k_p_rad = 1.0;
gamma_z_scatter_rescale_ = periodize(gamma_z_scatter_,0,2*pi)/spacing_gamma_z;
k_p_rad_scatter_rescale_ = max(0-1,min(n_k_p_rad-1+1,k_p_rad_scatter_/spacing_k_p_rad));
n_half_order = floor(n_order/2);
gamma_z_scatter_floor_ = floor(gamma_z_scatter_rescale_)-n_half_order; %<-- do not loop around boundary just yet. ;
k_p_rad_scatter_floor_ = min(n_k_p_rad-n_order,max(0,floor(k_p_rad_scatter_rescale_)-n_half_order));
gamma_z_scatter_shift_ = gamma_z_scatter_rescale_ - gamma_z_scatter_floor_;
k_p_rad_scatter_shift_ = k_p_rad_scatter_rescale_ - k_p_rad_scatter_floor_;
%%%%%%%%;
node_x_ = transpose(0:n_order-1);
weight_denominator_ = prod( repmat(node_x_,1,n_order) - repmat(transpose(node_x_),n_order,1) + eye(n_order,n_order) , 2 ) ;
weight_numerator_k_p_rad__ = zeros(n_order,n_scatter);
weight_numerator_gamma_z__ = zeros(n_order,n_scatter);
index_col_gamma_z__ = zeros(n_order,n_scatter);
index_col_k_p_rad__ = zeros(n_order,n_scatter);
index_row__ = repmat(0:n_scatter-1,n_order^2,1);
for nscatter=0:n_scatter-1;
weight_numerator_gamma_z__(:,1+nscatter) = prod( gamma_z_scatter_shift_(1+nscatter) - repmat(node_x_,1,n_order) + diag(1-(gamma_z_scatter_shift_(1+nscatter)-node_x_)) , 1 ) ;
weight_numerator_k_p_rad__(:,1+nscatter) = prod( k_p_rad_scatter_shift_(1+nscatter) - repmat(node_x_,1,n_order) + diag(1-(k_p_rad_scatter_shift_(1+nscatter)-node_x_)) , 1);
end;%for nscatter=0:n_scatter-1;
%%%%%%%%;
index_col_gamma_z__ = bsxfun(@plus,reshape(gamma_z_scatter_floor_,[1,n_scatter]),transpose([0:n_order-1]));
index_col_k_p_rad__ = bsxfun(@plus,reshape(k_p_rad_scatter_floor_,[1,n_scatter]),transpose([0:n_order-1]));
%%%%%%%%;
% periodize gamma_z. ;
%%%%%%%%;
tmp_index_ = efind(index_col_gamma_z__>=n_gamma_z); 
index_col_gamma_z__(1+tmp_index_) = index_col_gamma_z__(1+tmp_index_) - n_gamma_z;
tmp_index_ = efind(index_col_gamma_z__< 0        ); 
index_col_gamma_z__(1+tmp_index_) = index_col_gamma_z__(1+tmp_index_) + n_gamma_z;
%%%%%%%%;
weight_gamma_z__ = diag(1./weight_denominator_)*weight_numerator_gamma_z__;
weight_k_p_rad__ = diag(1./weight_denominator_)*weight_numerator_k_p_rad__;
weight__ = zeros(n_order^2,n_scatter);
index_col__ = zeros(n_order^2,n_scatter);
for nscatter=0:n_scatter-1;
weight__(:,1+nscatter) = reshape(weight_gamma_z__(:,1+nscatter)*transpose(weight_k_p_rad__(:,1+nscatter)),[n_order^2,1]);
index_col__(:,1+nscatter) = reshape(bsxfun(@plus,index_col_gamma_z__(:,1+nscatter),transpose(index_col_k_p_rad__(:,1+nscatter))*n_gamma_z),[n_order^2,1]);
end;%for nscatter=0:n_scatter-1;
%%%%%%%%;
scatter_from_tensor_swk__ = sparse(1+index_row__(:),1+index_col__(:),weight__(:),n_scatter,n_gamma_z*n_k_p_rad);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if nargout> 1;
%%%%%%%%;
% calculate scatter_from_tensor_dwndkm_swk__. ;
%%%%%%%%;
weight_numerator_dw1dk0_gamma_z__ = zeros(n_order,n_scatter);
weight_numerator_dw0dk1_k_p_rad__ = zeros(n_order,n_scatter);
%%%%;
for nscatter=0:n_scatter-1;
%%%%;
tmp_w = gamma_z_scatter_shift_(1+nscatter);
tmp_mat0__ = tmp_w - repmat(node_x_,1,n_order) + diag(1-(tmp_w-node_x_));
tmp_weight_numerator_dw1dk0_gamma_z_ = zeros(n_order,1);
for norder0=1:n_order-1;
tmp_mat1__ = circshift(diag(1-(tmp_w-circshift(node_x_,-norder0))),+norder0);
tmp_weight_numerator_dw1dk0_gamma_z_ = tmp_weight_numerator_dw1dk0_gamma_z_ + transpose(prod( tmp_mat0__ + tmp_mat1__ , 1 )) ;
end;%for norder0=1:n_order-1;
weight_numerator_dw1dk0_gamma_z__(:,1+nscatter) = tmp_weight_numerator_dw1dk0_gamma_z_;
clear tmp_w tmp_mat0__ tmp_mat1__ tmp_weight_numerator_dw1dk0_gamma_z_;
%%%%;
tmp_k = k_p_rad_scatter_shift_(1+nscatter);
tmp_mat0__ = tmp_k - repmat(node_x_,1,n_order) + diag(1-(tmp_k-node_x_));
tmp_weight_numerator_dw0dk1_k_p_rad_ = zeros(n_order,1);
for norder0=1:n_order-1;
tmp_mat1__ = circshift(diag(1-(tmp_k-circshift(node_x_,-norder0))),+norder0);
tmp_weight_numerator_dw0dk1_k_p_rad_ = tmp_weight_numerator_dw0dk1_k_p_rad_ + transpose(prod( tmp_mat0__ + tmp_mat1__ , 1 )) ;
end;%for norder0=1:n_order-1;
weight_numerator_dw0dk1_k_p_rad__(:,1+nscatter) = tmp_weight_numerator_dw0dk1_k_p_rad_;
clear tmp_k tmp_mat0__ tmp_mat1__ tmp_weight_numerator_dw0dk1_k_p_rad_;
%%%%;
end;%for nscatter=0:n_scatter-1;
weight_dw1dk0_gamma_z__ = diag(1./weight_denominator_)*weight_numerator_dw1dk0_gamma_z__/max(1e-12,spacing_gamma_z);
weight_dw0dk1_k_p_rad__ = diag(1./weight_denominator_)*weight_numerator_dw0dk1_k_p_rad__/max(1e-12,spacing_k_p_rad)/(k_p_rad_max/n_k_p_rad); %<-- correction. ;
%%%%;
weight_dw1dk0__ = zeros(n_order^2,n_scatter);
for nscatter=0:n_scatter-1;
weight_dw1dk0__(:,1+nscatter) = reshape(weight_dw1dk0_gamma_z__(:,1+nscatter)*transpose(weight_k_p_rad__(:,1+nscatter)),[n_order^2,1]);
end;%for nscatter=0:n_scatter-1;
scatter_from_tensor_dw1dk0_swk__ = sparse(1+index_row__(:),1+index_col__(:),weight_dw1dk0__(:),n_scatter,n_gamma_z*n_k_p_rad);
%%%%;
weight_dw0dk1__ = zeros(n_order^2,n_scatter);
for nscatter=0:n_scatter-1;
weight_dw0dk1__(:,1+nscatter) = reshape(weight_gamma_z__(:,1+nscatter)*transpose(weight_dw0dk1_k_p_rad__(:,1+nscatter)),[n_order^2,1]);
end;%for nscatter=0:n_scatter-1;
scatter_from_tensor_dw0dk1_swk__ = sparse(1+index_row__(:),1+index_col__(:),weight_dw0dk1__(:),n_scatter,n_gamma_z*n_k_p_rad);
%%%%;
weight_dw1dk1__ = zeros(n_order^2,n_scatter);
for nscatter=0:n_scatter-1;
weight_dw1dk1__(:,1+nscatter) = reshape(weight_dw1dk0_gamma_z__(:,1+nscatter)*transpose(weight_dw0dk1_k_p_rad__(:,1+nscatter)),[n_order^2,1]);
end;%for nscatter=0:n_scatter-1;
scatter_from_tensor_dw1dk1_swk__ = sparse(1+index_row__(:),1+index_col__(:),weight_dw1dk1__(:),n_scatter,n_gamma_z*n_k_p_rad);
%%%%%%%%;
weight_numerator_dw2dk0_gamma_z__ = zeros(n_order,n_scatter);
weight_numerator_dw0dk2_k_p_rad__ = zeros(n_order,n_scatter);
%%%%;
for nscatter=0:n_scatter-1;
%%%%;
tmp_w = gamma_z_scatter_shift_(1+nscatter);
tmp_mat0__ = tmp_w - repmat(node_x_,1,n_order) + diag(1-(tmp_w-node_x_));
tmp_weight_numerator_dw2dk0_gamma_z_ = zeros(n_order,1);
for norder0=1:n_order-1;
tmp_mat1__ = circshift(diag(1-(tmp_w-circshift(node_x_,-norder0))),+norder0);
for norder1=1+norder0:n_order-1;
tmp_mat2__ = circshift(diag(1-(tmp_w-circshift(node_x_,-norder1))),+norder1);
tmp_weight_numerator_dw2dk0_gamma_z_ = tmp_weight_numerator_dw2dk0_gamma_z_ + transpose(prod( tmp_mat0__ + tmp_mat1__ + tmp_mat2__, 1 )) ;
end;%for norder1=1+norder0:n_order-1;
end;%for norder0=1:n_order-1;
weight_numerator_dw2dk0_gamma_z__(:,1+nscatter) = tmp_weight_numerator_dw2dk0_gamma_z_;
clear tmp_w tmp_mat0__ tmp_mat1__ tmp_weight_numerator_dw2dk0_gamma_z_;
%%%%;
tmp_k = k_p_rad_scatter_shift_(1+nscatter);
tmp_mat0__ = tmp_k - repmat(node_x_,1,n_order) + diag(1-(tmp_k-node_x_));
tmp_weight_numerator_dw0dk2_k_p_rad_ = zeros(n_order,1);
for norder0=1:n_order-1;
tmp_mat1__ = circshift(diag(1-(tmp_k-circshift(node_x_,-norder0))),+norder0);
for norder1=1+norder0:n_order-1;
tmp_mat2__ = circshift(diag(1-(tmp_k-circshift(node_x_,-norder1))),+norder1);
tmp_weight_numerator_dw0dk2_k_p_rad_ = tmp_weight_numerator_dw0dk2_k_p_rad_ + transpose(prod( tmp_mat0__ + tmp_mat1__ + tmp_mat2__, 1 )) ;
end;%for norder1=1+norder0:n_order-1;
end;%for norder0=1:n_order-1;
weight_numerator_dw0dk2_k_p_rad__(:,1+nscatter) = tmp_weight_numerator_dw0dk2_k_p_rad_;
clear tmp_k tmp_mat0__ tmp_mat1__ tmp_weight_numerator_dw0dk2_k_p_rad_;
%%%%;
end;%for nscatter=0:n_scatter-1;
weight_dw2dk0_gamma_z__ = diag(1./weight_denominator_)*weight_numerator_dw2dk0_gamma_z__/max(1e-12,spacing_gamma_z)/max(1e-12,spacing_gamma_z)*2;
weight_dw0dk2_k_p_rad__ = diag(1./weight_denominator_)*weight_numerator_dw0dk2_k_p_rad__/max(1e-12,spacing_k_p_rad)/max(1e-12,spacing_k_p_rad)*2/(k_p_rad_max/n_k_p_rad)^2; %<-- correction. ;
%%%%;
weight_dw2dk0__ = zeros(n_order^2,n_scatter);
for nscatter=0:n_scatter-1;
weight_dw2dk0__(:,1+nscatter) = reshape(weight_dw2dk0_gamma_z__(:,1+nscatter)*transpose(weight_k_p_rad__(:,1+nscatter)),[n_order^2,1]);
end;%for nscatter=0:n_scatter-1;
scatter_from_tensor_dw2dk0_swk__ = sparse(1+index_row__(:),1+index_col__(:),weight_dw2dk0__(:),n_scatter,n_gamma_z*n_k_p_rad);
%%%%;
weight_dw0dk2__ = zeros(n_order^2,n_scatter);
for nscatter=0:n_scatter-1;
weight_dw0dk2__(:,1+nscatter) = reshape(weight_gamma_z__(:,1+nscatter)*transpose(weight_dw0dk2_k_p_rad__(:,1+nscatter)),[n_order^2,1]);
end;%for nscatter=0:n_scatter-1;
scatter_from_tensor_dw0dk2_swk__ = sparse(1+index_row__(:),1+index_col__(:),weight_dw0dk2__(:),n_scatter,n_gamma_z*n_k_p_rad);
%%%%%%%%;
end;%if nargout> 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

