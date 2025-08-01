function ...
[ ...
 parameter ...
,n_part_order ...
,k_s_0_rescale_gpu_ ...
,k_s_0_rescale_gpu_start_ ...
,k_s_0_rescale_gpu_shift_ ...
,node_x_ ...
,n_nodex ...
,index_col_0_gpu_xs__ ...
,index_row_gpu_xs__ ...
,weight_denominator_0_gpu_x_ ...
,weight_numerator_0_gpu_xs__ ...
,weight_dnumerator_0_gpu_xs__ ...
,weight_ddnumerator_0_gpu_xs__ ...
] = ...
scatter_from_tensor_helper_n_gpu_9( ...
 parameter ...
,n_scatter ...
,k_s_0_ ...
,k_c_0_start ...
,dk_c_0 ...
);

str_thisfunction = 'scatter_from_tensor_helper_n_gpu_9';

if nargin<1;
disp(sprintf(' %% testing %s',str_thisfunction));
disp(sprintf(' %% see interval_k_c_scatter_from_tensor_interpolate_n_9'));
disp('returning'); return;
end;%if nargin<1;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_scatter=[]; end; na=na+1;
if (nargin<1+na); k_s_0_=[]; end; na=na+1;
if (nargin<1+na); k_c_0_start=[]; end; na=na+1;
if (nargin<1+na); dk_c_0=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
flag_verbose=parameter.flag_verbose;
if ~isfield(parameter,'n_order'); parameter.n_order=5; end;
n_order=parameter.n_order;
if ~isfield(parameter,'gpu_batch_numel'); parameter.gpu_batch_numel=1e9; end;
gpu_batch_numel=parameter.gpu_batch_numel;
if ~isfield(parameter,'flag_d'); parameter.flag_d=1; end;
flag_d=parameter.flag_d;
if ~isfield(parameter,'flag_dd'); parameter.flag_dd=1; end;
flag_dd=parameter.flag_dd;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

n_part_order=[];
k_s_0_rescale_gpu_ = gpuArray( [] );
k_s_0_rescale_gpu_start_ = gpuArray( [] );
k_s_0_rescale_gpu_shift_ = gpuArray( [] );
node_x_=[];
n_nodex=[];
index_col_0_gpu_xs__ = gpuArray( [] );
index_row_gpu_xs__ = gpuArray( [] );
weight_denominator_0_gpu_x_ = gpuArray( [] );
weight_numerator_0_gpu_xs__ = gpuArray( [] );
weight_dnumerator_0_gpu_xs__ = gpuArray( [] );
weight_ddnumerator_0_gpu_xs__ = gpuArray( [] );

n_part_order = round((n_order-1)/2);
%%%%;
k_s_0_rescale_gpu_ = gpuArray( (k_s_0_ - k_c_0_start)/dk_c_0 ); %<-- unprotected (signed) divide, assume dk_c_0 is not close to 0. ;
k_s_0_rescale_gpu_start_ = floor(k_s_0_rescale_gpu_) - n_part_order;
k_s_0_rescale_gpu_shift_ = k_s_0_rescale_gpu_ - floor(k_s_0_rescale_gpu_) - 0.5 ;
%%%%;
if (flag_verbose>0);
disp(sprintf(' %% k_s_0_rescale_gpu_ in [%+0.6f , %+0.6f]',prctile(k_s_0_rescale_gpu_,[  0,100])));
end;%if (flag_verbose>0);
%%%%;
node_x_ = transpose(-n_part_order-0.5:1:+n_part_order+0.5);
node_gpu_x_ = gpuArray( node_x_ ) ;
if (flag_verbose>0); disp(sprintf(' %% node_x_: %s',num2str(transpose(node_x_),'%+0.2f '))); end;
n_nodex = numel(node_x_);
index_col_0_gpu_xs__ = gpuArray( zeros(n_nodex,n_scatter) );
index_col_0_gpu_xs__ = bsxfun(@plus,reshape(k_s_0_rescale_gpu_start_,[1,n_scatter]),transpose([0:n_nodex-1]));
index_row_gpu_xs__ = gpuArray( repmat(0:n_scatter-1,[n_nodex^1,1]) );
%%%%;
if (flag_verbose>0); disp(sprintf(' %% index_col_0_gpu_xs__ in [%.4d,%.4d]',min(index_col_0_gpu_xs__,[],'all'),max(index_col_0_gpu_xs__,[],'all'))); end;
%%%%%%%%;
weight_denominator_0_gpu_x_ = prod( repmat(node_gpu_x_,1,n_nodex) - repmat(transpose(node_gpu_x_),n_nodex,1) + eye(n_nodex,n_nodex) , 2 ) ;
if (flag_verbose>0); disp(sprintf(' %% weight_denominator_0_gpu_x_: %s',num2str(transpose(weight_denominator_0_gpu_x_),'%+0.2f '))); end;
%%%%;
%%%%;
weight_numerator_0_gpu_xs__ = gpuArray( zeros(n_nodex,n_scatter) );
tmp_t = tic();
weight_numerator_0_gpu_xxs___ = bsxfun(@minus,reshape(k_s_0_rescale_gpu_shift_,[1,1,n_scatter]),repmat(node_gpu_x_,[1,n_nodex,1]));
for nnodex=0:n_nodex-1; weight_numerator_0_gpu_xxs___(1+nnodex,1+nnodex,:) = 1; end;
weight_numerator_0_gpu_xs__ = reshape(prod(weight_numerator_0_gpu_xxs___,1),[n_nodex,n_scatter]);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% weight_numerator_0_gpu_xs__: time %0.6fs',tmp_t)); end;
%%%%;
if flag_d;
tmp_t = tic();
n_scatter_batch_size = max(1,min(n_scatter,floor(gpu_batch_numel/max(1,n_nodex^3))));
n_scatter_batch = ceil(n_scatter/n_scatter_batch_size);
if (flag_verbose> 1); disp(sprintf(' %% n_scatter %d n_nodex %d n_scatter_batch_size %d n_scatter_batch %d',n_scatter,n_nodex,n_scatter_batch_size,n_scatter_batch)); end;
for nscatter_batch=0:n_scatter_batch-1;
tmp_index_start = max(0,min(n_scatter-1,(nscatter_batch+0)*n_scatter_batch_size + 0));
tmp_index_final = max(0,min(n_scatter-1,(nscatter_batch+1)*n_scatter_batch_size - 1));
tmp_index_ = [tmp_index_start:tmp_index_final]; n_scatter_sub = numel(tmp_index_);
k_s_0_rescale_gpu_shift_sub_ = k_s_0_rescale_gpu_shift_(1+tmp_index_);
weight_dnumerator_0_gpu_sub_xs__ = local_d(n_nodex,node_gpu_x_,n_scatter_sub,k_s_0_rescale_gpu_shift_sub_);
weight_dnumerator_0_gpu_xs__(:,1+tmp_index_) = weight_dnumerator_0_gpu_sub_xs__;
end;%for nscatter_batch=0:n_scatter_batch-1;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% weight_dnumerator_0_gpu_xs__: time %0.6fs',tmp_t)); end;
end;%if flag_d;
%%%%;
if flag_dd;
tmp_t = tic();
n_scatter_batch_size = max(1,min(n_scatter,floor(gpu_batch_numel/max(1,n_nodex^4))));
n_scatter_batch = ceil(n_scatter/n_scatter_batch_size);
if (flag_verbose> 1); disp(sprintf(' %% n_scatter %d n_nodex %d n_scatter_batch_size %d n_scatter_batch %d',n_scatter,n_nodex,n_scatter_batch_size,n_scatter_batch)); end;
for nscatter_batch=0:n_scatter_batch-1;
tmp_index_start = max(0,min(n_scatter-1,(nscatter_batch+0)*n_scatter_batch_size + 0));
tmp_index_final = max(0,min(n_scatter-1,(nscatter_batch+1)*n_scatter_batch_size - 1));
tmp_index_ = [tmp_index_start:tmp_index_final]; n_scatter_sub = numel(tmp_index_);
k_s_0_rescale_gpu_shift_sub_ = k_s_0_rescale_gpu_shift_(1+tmp_index_);
weight_ddnumerator_0_gpu_sub_xs__ = local_dd(n_nodex,node_gpu_x_,n_scatter_sub,k_s_0_rescale_gpu_shift_sub_);
weight_ddnumerator_0_gpu_xs__(:,1+tmp_index_) = weight_ddnumerator_0_gpu_sub_xs__;
end;%for nscatter_batch=0:n_scatter_batch-1;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% weight_ddnumerator_0_gpu_xs__: time %0.6fs',tmp_t)); end;
end;%if flag_dd;
%%%%%%%%;

if (flag_verbose> 1);
disp(sprintf(' %% n_scatter %d',n_scatter));
disp(sprintf(' %% k_s_0_ %s',num2str(transpose(k_s_0_),'%0.2f ')));
disp(sprintf(' %% k_c_0_start %+0.6f',k_c_0_start));
disp(sprintf(' %% dk_c_0 %+0.6f',dk_c_0));
disp(sprintf(' %% n_part_order %d',n_part_order));
disp(sprintf(' %% k_s_0_rescale_gpu_ %s',num2str(transpose(k_s_0_rescale_gpu_),'%0.2f ')));
disp(sprintf(' %% k_s_0_rescale_gpu_start_ %s',num2str(transpose(k_s_0_rescale_gpu_start_),'%0.2f ')));
disp(sprintf(' %% k_s_0_rescale_gpu_shift_ %s',num2str(transpose(k_s_0_rescale_gpu_shift_),'%0.2f ')));
disp(sprintf(' %% node_gpu_x_ %s',num2str(transpose(node_gpu_x_),'%0.2f ')));
disp(sprintf(' %% n_nodex %d',n_nodex));
disp(sprintf(' %% index_col_0_gpu_xs__:')); disp(index_col_0_gpu_xs__);
disp(sprintf(' %% index_row_gpu_xs__:')); disp(index_row_gpu_xs__);
disp(sprintf(' %% weight_denominator_0_gpu_x_ %s',num2str(transpose(weight_denominator_0_gpu_x_),'%0.2f ')));
disp(sprintf(' %% weight_numerator_0_gpu_xs__:')); disp(weight_numerator_0_gpu_xs__);
if flag_d; disp(sprintf(' %% weight_dnumerator_0_gpu_xs__:')); disp(weight_dnumerator_0_gpu_xs__); end;
if flag_dd; disp(sprintf(' %% weight_ddnumerator_0_gpu_xs__:')); disp(weight_ddnumerator_0_gpu_xs__); end;
end;%if (flag_verbose> 1);

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

function ...
weight_dnumerator_0_gpu_xs__ = local_d(n_nodex,node_gpu_x_,n_scatter,k_s_0_rescale_gpu_shift_);
weight_dnumerator_0_gpu_xxxs____ = bsxfun(@minus,reshape(k_s_0_rescale_gpu_shift_,[1,1,1,n_scatter]),repmat(node_gpu_x_,[1,n_nodex,n_nodex,1]));
for nnodex_A=0:n_nodex-1;
weight_dnumerator_0_gpu_xxxs____(1+nnodex_A,1+nnodex_A,:,:) = 1;
end;%for nnodex_A=0:n_nodex-1;
for nnodex_B=0:n_nodex-1;
weight_dnumerator_0_gpu_xxxs____(1+nnodex_B,:,1+nnodex_B,:) = 1;
end;%for nnodex_B=0:n_nodex-1;
for nnodex_AB=0:n_nodex-1;
weight_dnumerator_0_gpu_xxxs____(1+nnodex_AB,1+nnodex_AB,1+nnodex_AB,:) = 0;
end;%for nnodex_AB=0:n_nodex-1;
weight_dnumerator_0_gpu_xs__ = reshape(sum(prod(weight_dnumerator_0_gpu_xxxs____,[1]),[2]),[n_nodex,n_scatter]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

function ...
weight_ddnumerator_0_gpu_xs__ = local_dd(n_nodex,node_gpu_x_,n_scatter,k_s_0_rescale_gpu_shift_);
weight_ddnumerator_0_gpu_xxxxs_____ = bsxfun( ...
					      @minus ...
					      ,reshape(k_s_0_rescale_gpu_shift_,[1,1,1,1,n_scatter]) ...
					      ,repmat(node_gpu_x_,[1,n_nodex,n_nodex,n_nodex,1]) ...
					    );
for nnodex_A=0:n_nodex-1;
weight_ddnumerator_0_gpu_xxxxs_____(1+nnodex_A,1+nnodex_A,:,:,:) = 1;
end;%for nnodex_A=0:n_nodex-1;
for nnodex_B=0:n_nodex-1;
weight_ddnumerator_0_gpu_xxxxs_____(1+nnodex_B,:,1+nnodex_B,:,:) = 1;
end;%for nnodex_B=0:n_nodex-1;
for nnodex_C=0:n_nodex-1;
weight_ddnumerator_0_gpu_xxxxs_____(1+nnodex_C,:,:,1+nnodex_C,:) = 1;
end;%for nnodex_C=0:n_nodex-1;
for nnodex_AB=0:n_nodex-1;
weight_ddnumerator_0_gpu_xxxxs_____(1+nnodex_AB,1+nnodex_AB,1+nnodex_AB,:,:) = 0;
end;%for nnodex_AB=0:n_nodex-1;
for nnodex_AC=0:n_nodex-1;
weight_ddnumerator_0_gpu_xxxxs_____(1+nnodex_AC,1+nnodex_AC,:,1+nnodex_AC,:) = 0;
end;%for nnodex_AC=0:n_nodex-1;
for nnodex_BC=0:n_nodex-1;
weight_ddnumerator_0_gpu_xxxxs_____(1+nnodex_BC,:,1+nnodex_BC,1+nnodex_BC,:) = 0;
end;%for nnodex_BC=0:n_nodex-1;
for nnodex_ABC=0:n_nodex-1;
weight_ddnumerator_0_gpu_xxxxs_____(1+nnodex_ABC,1+nnodex_ABC,1+nnodex_ABC,1+nnodex_ABC,:) = 0;
end;%for nnodex_ABC=0:n_nodex-1;
weight_ddnumerator_0_gpu_xs__ = reshape(sum(prod(weight_ddnumerator_0_gpu_xxxxs_____,[1]),[2,3]),[n_nodex,n_scatter]);

