function ...
[ ...
 scatter_from_tensor_s0__ ...
 dscatter_from_tensor_s0__ ...
 ddscatter_from_tensor_s0__ ...
] = ...
interval_k_c_scatter_from_tensor_interpolate_n_8( ...
 n_order ...
,n_k_c_0 ...
,k_c_0_lim_ ...
,n_scatter ...
,k_s_0_ ...
);
%%%%%%%%;
% Here we limit the scattered points to lie in the box defined by the tensor_grid. ;
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

str_thisfunction = 'interval_k_c_scatter_from_tensor_interpolate_n_8';

if (nargin<5);
rng(1);
%%%%%%%%;
flag_verbose = 1;
flag_disp = 1; nf=0;
disp(sprintf(' %% testing %s',str_thisfunction));
n_k_c_0 = 64+0;
k_p_r_max = 48/(2*pi);
k_c_0_lim_ = 1.0*k_p_r_max*[-1,+1];
k_c_0_ = transpose(linspace(k_c_0_lim_(1),k_c_0_lim_(2),n_k_c_0));
k_c_0_0_ = k_c_0_;
delta_ = [2.73];
f_c_0_ = exp(+i*(k_c_0_0_*delta_(1+0)));
df_c_0_ = (+i*delta_(1+0))*exp(+i*(k_c_0_0_*delta_(1+0)));
ddf_c_0_ = (+i*delta_(1+0))*(+i*delta_(1+0))*exp(+i*(k_c_0_0_*delta_(1+0)));
n_scatter = 32;
k_s_0_ = 0.5*k_p_r_max*(2*rand(n_scatter,1)-1);
n_order = 5;
[ ...
 scatter_from_tensor_s0__ ...
 dscatter_from_tensor_s0__ ...
 ddscatter_from_tensor_s0__ ...
] = ...
interval_k_c_scatter_from_tensor_interpolate_n_8( ...
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
if (n_k_c_0<=2*n_part_order+1); disp(sprintf(' %% Warning, n_k_c_0 %d, 2*n_part_order+1 %d in %s',n_k_c_0,2*n_part_order+1,str_thisfunction)); end;
%%%%;
dk_c_0 = diff(k_c_0_lim_)/max(1,n_k_c_0-1);
%%%%;
k_c_0_lim_use_ = [k_c_0_lim_(1)+n_part_order*dk_c_0,k_c_0_lim_(2)-n_part_order*dk_c_0];
%%%%;
if numel(efind(k_s_0_<=min(k_c_0_lim_use_)))>0; disp(sprintf(' %% Warning, k_s_0_ out of bounds in %s',str_thisfunction)); end;
if numel(efind(k_s_0_>=max(k_c_0_lim_use_)))>0; disp(sprintf(' %% Warning, k_s_0_ out of bounds in %s',str_thisfunction)); end;
k_s_0_ = max(min(k_c_0_lim_use_),min(max(k_c_0_lim_use_),k_s_0_));
%%%%;
k_c_0_start = k_c_0_lim_(1); k_c_0_final = k_c_0_lim_(2);
%%%%;
if (flag_verbose>0);
disp(sprintf(' %% k_c_0_start %+0.6f k_c_0_final %+0.6f: n_k_c_0 %d --> dk_c_0 %+0.6f',k_c_0_start,k_c_0_final,n_k_c_0,dk_c_0));
end;%if (flag_verbose>0);
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% calculate scatter_from_tensor_sba__. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
k_s_0_rescale_ = (k_s_0_ - k_c_0_start)/max(1e-12,dk_c_0);
%%%%;
if (flag_verbose>0);
disp(sprintf(' %% k_s_0_rescale_ in [%+0.6f , %+0.6f]',prctile(k_s_0_rescale_,[  0,100])));
end;%if (flag_verbose>0);
%%%%;
k_s_0_rescale_start_ = floor(k_s_0_rescale_) - n_part_order;
node_x_ = transpose(-n_part_order-0.5:1:+n_part_order+0.5);
if (flag_verbose>0); disp(sprintf(' %% node_x_: %s',num2str(transpose(node_x_),'%+0.2f '))); end;
n_nodex = numel(node_x_);
weight_denominator_x_ = prod( repmat(node_x_,1,n_nodex) - repmat(transpose(node_x_),n_nodex,1) + eye(n_nodex,n_nodex) , 2 ) ;
if (flag_verbose>0); disp(sprintf(' %% weight_denominator_x_: %s',num2str(transpose(weight_denominator_x_),'%+0.2f '))); end;
%%%%;
if (flag_verbose>1);
for tmp_x = [-1.5,-0.5,0,+0.5,+1.75];
disp(sprintf(' %% tmp_x: %0.6f',tmp_x));
tmp_d_ = prod( tmp_x - repmat(node_x_,1,n_nodex) + diag(1-(tmp_x-node_x_)) , 1 );
for nnodex=0:n_nodex-1;
tmp_d = prod(tmp_x - [node_x_(1 + [0:nnodex-1]);node_x_(1+[nnodex+1:n_nodex-1])]);
disp(sprintf(' %% nnodex %d/%d: %+0.4f vs %+0.4f',nnodex,n_nodex,tmp_d,tmp_d_(1+nnodex)));
end;%for nnodex=0:n_nodex-1;
end;%for tmp_x = [-1.5,-0.5,0,+0.5,+1.75];
end;%if (flag_verbose>1);
%%%%;
k_s_0_rescale_shift_ = k_s_0_rescale_ - k_s_0_rescale_start_ - n_part_order - 0.5 ;
%%%%;
weight_numerator_0_xs__ = zeros(n_nodex,n_scatter);
index_col_0_xs__ = zeros(n_nodex,n_scatter);
index_row_xs__ = repmat(0:n_scatter-1,[n_nodex^1,1]);
tmp_t = tic();
weight_numerator_0_xxs___ = bsxfun(@minus,reshape(k_s_0_rescale_shift_,[1,1,n_scatter]),repmat(node_x_,[1,n_nodex,1]));
for nnodex=0:n_nodex-1; weight_numerator_0_xxs___(1+nnodex,1+nnodex,:) = 1; end;
weight_numerator_0_xs__ = reshape(prod(weight_numerator_0_xxs___,1),[n_nodex,n_scatter]);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% weight_numerator_0_xs__: %0.6fs',tmp_t)); end;
%%%%;
% slower version. ;
%%%%;
flag_check=0;
if flag_check;
tmp_t = tic();
for nscatter=0:n_scatter-1;
weight_numerator_0_xs__(:,1+nscatter) = prod( k_s_0_rescale_shift_(1+nscatter) - repmat(node_x_,1,n_nodex) + diag(1-(k_s_0_rescale_shift_(1+nscatter)-node_x_)) , 1 ) ;
end;%for nscatter=0:n_scatter-1;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% weight_numerator_0_xs__: %0.6fs',tmp_t)); end;
end;%if flag_check;
%%%%;
if flag_d;
tmp_t = tic();
weight_dnumerator_0_xxxs____ = bsxfun(@minus,reshape(k_s_0_rescale_shift_,[1,1,1,n_scatter]),repmat(node_x_,[1,n_nodex,n_nodex,1]));
for nnodex_A=0:n_nodex-1;
weight_dnumerator_0_xxxs____(1+nnodex_A,1+nnodex_A,:,:) = 1;
end;%for nnodex_A=0:n_nodex-1;
for nnodex_B=0:n_nodex-1;
weight_dnumerator_0_xxxs____(1+nnodex_B,:,1+nnodex_B,:) = 1;
end;%for nnodex_B=0:n_nodex-1;
for nnodex_AB=0:n_nodex-1;
weight_dnumerator_0_xxxs____(1+nnodex_AB,1+nnodex_AB,1+nnodex_AB,:) = 0;
end;%for nnodex_AB=0:n_nodex-1;
weight_dnumerator_0_xs__ = reshape(sum(prod(weight_dnumerator_0_xxxs____,[1]),[2]),[n_nodex,n_scatter]);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% weight_dnumerator_0_xs__: %0.6fs',tmp_t)); end;
end;%if flag_d;
%%%%;
if flag_dd;
tmp_t = tic();
weight_ddnumerator_0_xxxxs_____ = bsxfun(@minus,reshape(k_s_0_rescale_shift_,[1,1,1,1,n_scatter]),repmat(node_x_,[1,n_nodex,n_nodex,n_nodex,1]));
for nnodex_A=0:n_nodex-1;
weight_ddnumerator_0_xxxxs_____(1+nnodex_A,1+nnodex_A,:,:,:) = 1;
end;%for nnodex_A=0:n_nodex-1;
for nnodex_B=0:n_nodex-1;
weight_ddnumerator_0_xxxxs_____(1+nnodex_B,:,1+nnodex_B,:,:) = 1;
end;%for nnodex_B=0:n_nodex-1;
for nnodex_C=0:n_nodex-1;
weight_ddnumerator_0_xxxxs_____(1+nnodex_C,:,:,1+nnodex_C,:) = 1;
end;%for nnodex_C=0:n_nodex-1;
for nnodex_AB=0:n_nodex-1;
weight_ddnumerator_0_xxxxs_____(1+nnodex_AB,1+nnodex_AB,1+nnodex_AB,:,:) = 0;
end;%for nnodex_AB=0:n_nodex-1;
for nnodex_AC=0:n_nodex-1;
weight_ddnumerator_0_xxxxs_____(1+nnodex_AC,1+nnodex_AC,:,1+nnodex_AC,:) = 0;
end;%for nnodex_AC=0:n_nodex-1;
for nnodex_BC=0:n_nodex-1;
weight_ddnumerator_0_xxxxs_____(1+nnodex_BC,:,1+nnodex_BC,1+nnodex_BC,:) = 0;
end;%for nnodex_BC=0:n_nodex-1;
for nnodex_ABC=0:n_nodex-1;
weight_ddnumerator_0_xxxxs_____(1+nnodex_ABC,1+nnodex_ABC,1+nnodex_ABC,1+nnodex_ABC,:) = 0;
end;%for nnodex_ABC=0:n_nodex-1;
weight_ddnumerator_0_xs__ = reshape(sum(prod(weight_ddnumerator_0_xxxxs_____,[1]),[2,3]),[n_nodex,n_scatter]);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% weight_ddnumerator_0_xs__: %0.6fs',tmp_t)); end;
end;%if flag_dd;
%%%%%%%%;
index_col_0_xs__ = bsxfun(@plus,reshape(k_s_0_rescale_start_,[1,n_scatter]),transpose([0:n_nodex-1]));
%%%%;
if (flag_verbose>0); disp(sprintf(' %% n_k_c_0 %d: index_col_0_xs__ in [%.4d,%.4d]',n_k_c_0,min(index_col_0_xs__,[],'all'),max(index_col_0_xs__,[],'all'))); end;
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figsml;
plot(index_col_0_xs__(:),'.');
ylabel('index_col_0_xs__','Interpreter','none');
end;%if (flag_disp>1);
%%%%%%%%;

%%%%;
if (flag_verbose>1);
%%%%;
[nc_0_min,tmp_ij] = min(index_col_0_xs__(:)); tmp_index = tmp_ij-1;
tmp_index_x = mod(tmp_index,n_nodex);
tmp_index_s = (tmp_index-tmp_index_x)/max(1,n_nodex);
disp(sprintf(' %% tmp_index %d (%d,%d) <-- nc_0_min %d',tmp_index,tmp_index_x,tmp_index_s,nc_0_min));
disp(sprintf(' %% index_col_0_xs__(:,1+tmp_index_s): %s',num2str(transpose(index_col_0_xs__(:,1+tmp_index_s)),' %+d')));
%%%%;
disp('returning'); return;
end;%if (flag_verbose>1);
%%%%;

%%%%%%%%;
weight_0_xs__ = bsxfun(@rdivide,weight_numerator_0_xs__,reshape(weight_denominator_x_,[n_nodex,1]));
weight_xs__ = weight_0_xs__ ;
index_col_xs__ = index_col_0_xs__;
scatter_from_tensor_s0__ = sparse(1+index_row_xs__(:),1+index_col_xs__(:),weight_xs__(:),n_scatter,n_k_c_0);
%%%%%%%%;
if flag_d;
dweight_0_xs__ = bsxfun(@rdivide,weight_dnumerator_0_xs__,reshape(weight_denominator_x_,[n_nodex,1])) / max(1e-12,dk_c_0) ;
dweight_xs__ = dweight_0_xs__ ;
index_col_xs__ = index_col_0_xs__;
dscatter_from_tensor_s0__ = sparse(1+index_row_xs__(:),1+index_col_xs__(:),dweight_xs__(:),n_scatter,n_k_c_0);
end;%if flag_d;
%%%%%%%%;
if flag_dd;
ddweight_0_xs__ = bsxfun(@rdivide,weight_ddnumerator_0_xs__,reshape(weight_denominator_x_,[n_nodex,1])) / max(1e-12,dk_c_0)^2 ;
ddweight_xs__ = ddweight_0_xs__ ;
index_col_xs__ = index_col_0_xs__;
ddscatter_from_tensor_s0__ = sparse(1+index_row_xs__(:),1+index_col_xs__(:),ddweight_xs__(:),n_scatter,n_k_c_0);
end;%if flag_dd;
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;


