function ...
[ ...
 scatter_from_tensor_s0__ ...
 dscatter_from_tensor_s0__ ...
 ddscatter_from_tensor_s0__ ...
] = ...
interval_periodic_k_c_scatter_from_tensor_interpolate_n_9( ...
 n_order ...
,n_k_c_0 ...
,k_c_0_lim_ ...
,n_scatter ...
,k_s_0_ ...
);
%%%%%%%%;
% Here we interpret the tensor_grid on a periodic array. ;
%%%%%%%%;
% returns sparse matrix encoding the n_order interpolation operator between a tensor_grid ;
% using k_c_?_. ;
%%%%%%%%;
% We presume these points are associated with an (n_k_c_0) array: ;
% a_k_c_0_. ;
%%%%;
% The n_scatter points to be interpolated have coordinates stored in ;
% k_s_0_, ;
% ;
% The (n_scatter)-by-(n_k_0) scatter_from_tensor_s0__ matrix stores the ;
% interpolation weights linking a_k_c_0_ to the unrolled a_k_s_. ;
%%%%;

str_thisfunction = 'interval_periodic_k_c_scatter_from_tensor_interpolate_n_9';

if (nargin<5);
rng(1);
%%%%%%%%;
flag_verbose = 1;
flag_disp = 1; nf=0;
disp(sprintf(' %% testing %s',str_thisfunction));
n_k_c_0 = 64+0;
k_p_r_max = 48/(2*pi);
k_c_0_lim_ = [0,2*pi];
k_c_0_ = linspace(k_c_0_lim_(1),k_c_0_lim_(2),n_k_c_0+1); k_c_0_ = transpose(k_c_0_(1:n_k_c_0));
k_c_0_0_ = k_c_0_;
delta_ = [5.0];
f_c_0_ = exp(+i*(k_c_0_0_*delta_(1+0)));
df_c_0_ = (+i*delta_(1+0))*exp(+i*(k_c_0_0_*delta_(1+0)));
ddf_c_0_ = (+i*delta_(1+0))*(+i*delta_(1+0))*exp(+i*(k_c_0_0_*delta_(1+0)));
n_scatter = 32*8;
k_s_0_ = 2*pi*rand(n_scatter,1);
n_order = 5;
[ ...
 scatter_from_tensor_s0__ ...
 dscatter_from_tensor_s0__ ...
 ddscatter_from_tensor_s0__ ...
] = ...
interval_periodic_k_c_scatter_from_tensor_interpolate_n_9( ...
 n_order ...
,n_k_c_0 ...
,k_c_0_lim_ ...
,n_scatter ...
,k_s_0_ ...
);
f_s_form_ = exp(+i*(k_s_0_*delta_(1+0)));
f_s_inte_ = scatter_from_tensor_s0__*f_c_0_(:);
df_s_form_ = (+i*delta_(1+0))*exp(+i*(k_s_0_*delta_(1+0)));
df_s_inte_ = dscatter_from_tensor_s0__*f_c_0_(:);
ddf_s_form_ = (+i*delta_(1+0))*(+i*delta_(1+0))*exp(+i*(k_s_0_*delta_(1+0)));
ddf_s_inte_ = ddscatter_from_tensor_s0__*f_c_0_(:);
%%%%;
figure(1+nf);clf;figbig;
subplot(3,1,1)
hold on;
plot(k_c_0_,real(f_c_0_),'co');
plot(k_s_0_,real(f_s_form_),'go');
plot(k_s_0_,real(f_s_inte_),'rx');
hold off; xlim(k_c_0_lim_);
xlabel('k_c_0_','Interpreter','none');
ylabel('f');
subplot(3,1,2)
hold on;
plot(k_c_0_,real(df_c_0_),'co');
plot(k_s_0_,real(df_s_form_),'go');
plot(k_s_0_,real(df_s_inte_),'rx');
hold off; xlim(k_c_0_lim_);
xlabel('k_c_0_','Interpreter','none');
ylabel('df');
subplot(3,1,3)
hold on;
plot(k_c_0_,real(ddf_c_0_),'co');
plot(k_s_0_,real(ddf_s_form_),'go');
plot(k_s_0_,real(ddf_s_inte_),'rx');
hold off; xlim(k_c_0_lim_);
xlabel('k_c_0_','Interpreter','none');
ylabel('ddf');
%%%%;
fnorm_disp(flag_verbose,'f_s_form_',f_s_form_,'f_s_inte_',f_s_inte_);
fnorm_disp(flag_verbose,'df_s_form_',df_s_form_,'df_s_inte_',df_s_inte_);
fnorm_disp(flag_verbose,'ddf_s_form_',ddf_s_form_,'ddf_s_inte_',ddf_s_inte_);
%%%%%%%%;
disp(sprintf(' %% returning')); return;
end;%if (nargin<6);

flag_verbose = 0;
flag_disp = 0; nf=0;
flag_check = 0;
flag_d = (nargout>=2);
flag_dd = (nargout>=3);

if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

na=0;
if (nargin<1+na); n_order=[]; end; na=na+1;
if (nargin<1+na); n_k_c_0=[]; end; na=na+1;
if (nargin<1+na); k_c_0_lim_=[]; end; na=na+1;
if (nargin<1+na); n_scatter=[]; end; na=na+1;
if (nargin<1+na); k_s_0_=[]; end; na=na+1;

%%%%%%%%;
n_part_order = round((n_order-1)/2);
if (flag_verbose>0); disp(sprintf(' %% n_order %d n_part_order %d',n_order,n_part_order)); end;
%%%%;
dk_c_0 = diff(k_c_0_lim_)/max(1,n_k_c_0-0); %<-- interpreted as periodic. ;
%%%%;
k_c_0_start = k_c_0_lim_(1); k_c_0_final = k_c_0_lim_(2);
%%%%;
if (flag_verbose>0);
disp(sprintf(' %% k_c_0_start %+0.6f k_c_0_final %+0.6f: n_k_c_0 %d --> dk_c_0 %+0.6f',k_c_0_start,k_c_0_final,n_k_c_0,dk_c_0));
end;%if (flag_verbose>0);
%%%%%%%%;

%%%%%%%%;
parameter = struct('type','parameter');
parameter.flag_verbose = flag_verbose;
parameter.n_order = n_order;
parameter.flag_d = flag_d;
parameter.flag_dd = flag_dd;
[ ...
 parameter ...
,n_part_order ...
,k_s_0_rescale_ ...
,k_s_0_rescale_start_ ...
,k_s_0_rescale_shift_ ...
,node_x_ ...
,n_nodex ...
,index_col_0_xs__ ...
,index_row_xs__ ...
,weight_denominator_x_ ...
,weight_numerator_0_xs__ ...
,weight_dnumerator_0_xs__ ...
,weight_ddnumerator_0_xs__ ...
] = ...
scatter_from_tensor_helper_n_9( ...
 parameter ...
,n_scatter ...
,k_s_0_ ...
,k_c_0_start ...
,dk_c_0 ...
);
%%%%%%%%;
index_col_0_xs__ = periodize(index_col_0_xs__,0,n_k_c_0);
%%%%%%%%;
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figsml;
plot(index_col_0_xs__(:),'.');
ylabel('index_col_0_xs__','Interpreter','none');
end;%if (flag_disp>1);
%%%%%%%%;

%%%%%%%%;
weight_0_xs__ = bsxfun(@rdivide,weight_numerator_0_xs__,reshape(weight_denominator_x_,[n_nodex,1]));
weight_xs__ = weight_0_xs__ ;
index_col_xs__ = index_col_0_xs__;
scatter_from_tensor_s0__ = sparse(1+index_row_xs__(:),1+index_col_xs__(:),weight_xs__(:),n_scatter,n_k_c_0);
%%%%%%%%;
if flag_d;
dweight_0_xs__ = bsxfun(@rdivide,weight_dnumerator_0_xs__,reshape(weight_denominator_x_,[n_nodex,1])) / dk_c_0 ;
dweight_xs__ = dweight_0_xs__ ;
index_col_xs__ = index_col_0_xs__;
dscatter_from_tensor_s0__ = sparse(1+index_row_xs__(:),1+index_col_xs__(:),dweight_xs__(:),n_scatter,n_k_c_0);
end;%if flag_d;
%%%%%%%%;
if flag_dd;
ddweight_0_xs__ = bsxfun(@rdivide,weight_ddnumerator_0_xs__,reshape(weight_denominator_x_,[n_nodex,1])) / dk_c_0^2 ;
ddweight_xs__ = ddweight_0_xs__ ;
index_col_xs__ = index_col_0_xs__;
ddscatter_from_tensor_s0__ = sparse(1+index_row_xs__(:),1+index_col_xs__(:),ddweight_xs__(:),n_scatter,n_k_c_0);
end;%if flag_dd;
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;


