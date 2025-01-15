function ...
[ ...
 scatter_from_tensor_s01__ ...
 dscatterd0_from_tensor_s01__ ...
 dscatterd1_from_tensor_s01__ ...
 ddscatterd00_from_tensor_s01__ ...
 ddscatterd01_from_tensor_s01__ ...
 ddscatterd11_from_tensor_s01__ ...
] = ...
cylinder_k_c_scatter_from_tensor_interpolate_n_9( ...
 n_order ...
,n_k_c_0 ...
,n_k_c_1 ...
,k_c_0_lim_ ...
,k_c_1_lim_ ...
,n_scatter ...
,k_s_0_ ...
,k_s_1_ ...
);
%%%%%%%%;
% We interpret the k_c_0_ array as a standard tensor array (e.g., height). ;
% We interpret the k_c_1_ array as periodic (e.g., angle). ; 
%%%%%%%%;
% returns sparse matrix encoding the n_order interpolation operator between a tensor_grid ;
% using k_c_?_. ;
%%%%%%%%;
% We presume these points are associated with an (n_k_c_0)-by-(n_k_c_1) array: ;
% a_k_c_01__. ;
%%%%;
% The n_scatter points to be interpolated have coordinates stored in ;
% k_s_0_, ;
% k_s_1_, ;
% ;
% The (n_scatter)-by-(n_k_0*n_k_1) scatter_from_tensor_s01__ matrix stores the ;
% interpolation weights linking a_k_c_01_ to the unrolled a_k_s_. ;
%%%%;

str_thisfunction = 'cylinder_k_c_scatter_from_tensor_interpolate_n_9';

if (nargin<7);
rng(1);
%%%%%%%%;
flag_verbose = 1;
flag_disp = 1; nf=0;
disp(sprintf(' %% testing %s',str_thisfunction));
n_k_c_0 = 64+0;
n_k_c_1 = 64+1;
k_p_r_max = 48/(2*pi);
k_c_0_lim_ = 1.0*k_p_r_max*[-1,+1];
k_c_1_lim_ = [0,2*pi];
k_c_0_ = transpose(linspace(k_c_0_lim_(1),k_c_0_lim_(2),n_k_c_0));
k_c_1_ = linspace(k_c_1_lim_(1),k_c_1_lim_(2),n_k_c_1+1); k_c_1_ = transpose(k_c_1_(1:n_k_c_1));
[k_c_0_01__,k_c_1_01__] = ndgrid(k_c_0_,k_c_1_);
delta_ = [0.6;3.0];
f_c_01__ = exp(+i*(k_c_0_01__*delta_(1+0) + k_c_1_01__*delta_(1+1)));
dfd0_c_01__ = (+i*delta_(1+0))*exp(+i*(k_c_0_01__*delta_(1+0) + k_c_1_01__*delta_(1+1)));
dfd1_c_01__ = (+i*delta_(1+1))*exp(+i*(k_c_0_01__*delta_(1+0) + k_c_1_01__*delta_(1+1)));
ddfd00_c_01__ = (+i*delta_(1+0))*(+i*delta_(1+0))*exp(+i*(k_c_0_01__*delta_(1+0) + k_c_1_01__*delta_(1+1)));
ddfd01_c_01__ = (+i*delta_(1+0))*(+i*delta_(1+1))*exp(+i*(k_c_0_01__*delta_(1+0) + k_c_1_01__*delta_(1+1)));
ddfd11_c_01__ = (+i*delta_(1+1))*(+i*delta_(1+1))*exp(+i*(k_c_0_01__*delta_(1+0) + k_c_1_01__*delta_(1+1)));
n_scatter = 1024;
k_s_0_ = 0.5*k_p_r_max*(2*rand(n_scatter,1)-1);
k_s_1_ = 2*pi*rand(n_scatter,1);
n_order = 5;
[ ...
 scatter_from_tensor_s01__ ...
 dscatterd0_from_tensor_s01__ ...
 dscatterd1_from_tensor_s01__ ...
 ddscatterd00_from_tensor_s01__ ...
 ddscatterd01_from_tensor_s01__ ...
 ddscatterd11_from_tensor_s01__ ...
] = ...
cylinder_k_c_scatter_from_tensor_interpolate_n_9( ...
 n_order ...
,n_k_c_0 ...
,n_k_c_1 ...
,k_c_0_lim_ ...
,k_c_1_lim_ ...
,n_scatter ...
,k_s_0_ ...
,k_s_1_ ...
);
f_s_form_ = exp(+i*(k_s_0_*delta_(1+0) + k_s_1_*delta_(1+1)));
dfd0_s_form_ = (+i*delta_(1+0))*exp(+i*(k_s_0_*delta_(1+0) + k_s_1_*delta_(1+1)));
dfd1_s_form_ = (+i*delta_(1+1))*exp(+i*(k_s_0_*delta_(1+0) + k_s_1_*delta_(1+1)));
ddfd00_s_form_ = (+i*delta_(1+0))*(+i*delta_(1+0))*exp(+i*(k_s_0_*delta_(1+0) + k_s_1_*delta_(1+1)));
ddfd01_s_form_ = (+i*delta_(1+0))*(+i*delta_(1+1))*exp(+i*(k_s_0_*delta_(1+0) + k_s_1_*delta_(1+1)));
ddfd11_s_form_ = (+i*delta_(1+1))*(+i*delta_(1+1))*exp(+i*(k_s_0_*delta_(1+0) + k_s_1_*delta_(1+1)));
f_s_inte_ = scatter_from_tensor_s01__*f_c_01__(:);
dfd0_s_inte_ = dscatterd0_from_tensor_s01__*f_c_01__(:);
dfd1_s_inte_ = dscatterd1_from_tensor_s01__*f_c_01__(:);
ddfd00_s_inte_ = ddscatterd00_from_tensor_s01__*f_c_01__(:);
ddfd01_s_inte_ = ddscatterd01_from_tensor_s01__*f_c_01__(:);
ddfd11_s_inte_ = ddscatterd11_from_tensor_s01__*f_c_01__(:);
%%%%;
figure(1+nf);clf;figbig;
p_row = 2; p_col = 3; np=0;
subplot(p_row,p_col,1+np);np=np+1;
hold on;
plot3(k_c_0_01__(:),k_c_1_01__(:),real(f_c_01__(:)),'co');
plot3(k_s_0_(:),k_s_1_(:),real(f_s_form_),'go');
plot3(k_s_0_(:),k_s_1_(:),real(f_s_inte_),'rx');
hold off;
xlim(k_c_0_lim_); xlabel('k_c_0_','Interpreter','none');
ylim(k_c_1_lim_); ylabel('k_c_1_','Interpreter','none');
zlabel('f');
axis vis3d;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
hold on;
plot3(k_c_0_01__(:),k_c_1_01__(:),real(dfd0_c_01__(:)),'co');
plot3(k_s_0_(:),k_s_1_(:),real(dfd0_s_form_),'go');
plot3(k_s_0_(:),k_s_1_(:),real(dfd0_s_inte_),'rx');
hold off;
xlim(k_c_0_lim_); xlabel('k_c_0_','Interpreter','none');
ylim(k_c_1_lim_); ylabel('k_c_1_','Interpreter','none');
zlabel('dfd0');
axis vis3d;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
hold on;
plot3(k_c_0_01__(:),k_c_1_01__(:),real(dfd1_c_01__(:)),'co');
plot3(k_s_0_(:),k_s_1_(:),real(dfd1_s_form_),'go');
plot3(k_s_0_(:),k_s_1_(:),real(dfd1_s_inte_),'rx');
hold off;
xlim(k_c_0_lim_); xlabel('k_c_0_','Interpreter','none');
ylim(k_c_1_lim_); ylabel('k_c_1_','Interpreter','none');
zlabel('dfd1');
axis vis3d;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
hold on;
plot3(k_c_0_01__(:),k_c_1_01__(:),real(ddfd00_c_01__(:)),'co');
plot3(k_s_0_(:),k_s_1_(:),real(ddfd00_s_form_),'go');
plot3(k_s_0_(:),k_s_1_(:),real(ddfd00_s_inte_),'rx');
hold off;
xlim(k_c_0_lim_); xlabel('k_c_0_','Interpreter','none');
ylim(k_c_1_lim_); ylabel('k_c_1_','Interpreter','none');
zlabel('ddfd00');
axis vis3d;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
hold on;
plot3(k_c_0_01__(:),k_c_1_01__(:),real(ddfd01_c_01__(:)),'co');
plot3(k_s_0_(:),k_s_1_(:),real(ddfd01_s_form_),'go');
plot3(k_s_0_(:),k_s_1_(:),real(ddfd01_s_inte_),'rx');
hold off;
xlim(k_c_0_lim_); xlabel('k_c_0_','Interpreter','none');
ylim(k_c_1_lim_); ylabel('k_c_1_','Interpreter','none');
zlabel('ddfd01');
axis vis3d;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
hold on;
plot3(k_c_0_01__(:),k_c_1_01__(:),real(ddfd11_c_01__(:)),'co');
plot3(k_s_0_(:),k_s_1_(:),real(ddfd11_s_form_),'go');
plot3(k_s_0_(:),k_s_1_(:),real(ddfd11_s_inte_),'rx');
hold off;
xlim(k_c_0_lim_); xlabel('k_c_0_','Interpreter','none');
ylim(k_c_1_lim_); ylabel('k_c_1_','Interpreter','none');
zlabel('ddfd11');
axis vis3d;
%%%%%%%%;
fnorm_disp(flag_verbose,'f_s_form_',f_s_form_,'f_s_inte_',f_s_inte_);
fnorm_disp(flag_verbose,'dfd0_s_form_',dfd0_s_form_,'dfd0_s_inte_',dfd0_s_inte_);
fnorm_disp(flag_verbose,'dfd1_s_form_',dfd1_s_form_,'dfd1_s_inte_',dfd1_s_inte_);
fnorm_disp(flag_verbose,'ddfd00_s_form_',ddfd00_s_form_,'ddfd00_s_inte_',ddfd00_s_inte_);
fnorm_disp(flag_verbose,'ddfd01_s_form_',ddfd01_s_form_,'ddfd01_s_inte_',ddfd01_s_inte_);
fnorm_disp(flag_verbose,'ddfd11_s_form_',ddfd11_s_form_,'ddfd11_s_inte_',ddfd11_s_inte_);
%%%%%%%%;
disp(sprintf(' %% returning')); return;
end;%if (nargin<6);

flag_verbose = 0;
flag_disp = 0; nf=0;
flag_check = 0;
flag_d = (nargout>=2);
flag_dd = (nargout>=4);

if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

na=0;
if (nargin<1+na); n_order=[]; end; na=na+1;
if (nargin<1+na); n_k_c_0=[]; end; na=na+1;
if (nargin<1+na); n_k_c_1=[]; end; na=na+1;
if (nargin<1+na); k_c_0_lim_=[]; end; na=na+1;
if (nargin<1+na); k_c_1_lim_=[]; end; na=na+1;
if (nargin<1+na); n_scatter=[]; end; na=na+1;
if (nargin<1+na); k_s_0_=[]; end; na=na+1;
if (nargin<1+na); k_s_1_=[]; end; na=na+1;

%%%%%%%%;
n_part_order = round((n_order-1)/2);
if (flag_verbose>0); disp(sprintf(' %% n_order %d n_part_order %d',n_order,n_part_order)); end;
%%%%;
if (n_k_c_0<=2*n_part_order+1); disp(sprintf(' %% Warning, n_k_c_0 %d, 2*n_part_order+1 %d in %s',n_k_c_0,2*n_part_order+1,str_thisfunction)); end;
%%%%;
dk_c_0 = diff(k_c_0_lim_)/max(1,n_k_c_0-1); %<-- interpreted as interval. ;
dk_c_1 = diff(k_c_1_lim_)/max(1,n_k_c_1-0); %<-- interpreted as periodic array. ;
%%%%;
k_c_0_lim_use_ = [k_c_0_lim_(1)+n_part_order*dk_c_0,k_c_0_lim_(2)-n_part_order*dk_c_0];
%%%%;
if numel(efind(k_s_0_<=min(k_c_0_lim_use_)))>0; disp(sprintf(' %% Warning, k_s_0_ out of bounds in %s',str_thisfunction)); end;
if numel(efind(k_s_0_>=max(k_c_0_lim_use_)))>0; disp(sprintf(' %% Warning, k_s_0_ out of bounds in %s',str_thisfunction)); end;
k_s_0_ = max(min(k_c_0_lim_use_),min(max(k_c_0_lim_use_),k_s_0_));
%%%%;
k_c_0_start = k_c_0_lim_(1); k_c_0_final = k_c_0_lim_(2);
k_c_1_start = k_c_1_lim_(1); k_c_1_final = k_c_1_lim_(2);
%%%%;
if (flag_verbose>0);
disp(sprintf(' %% k_c_0_start %+0.6f k_c_0_final %+0.6f: n_k_c_0 %d --> dk_c_0 %+0.6f',k_c_0_start,k_c_0_final,n_k_c_0,dk_c_0));
disp(sprintf(' %% k_c_1_start %+0.6f k_c_1_final %+0.6f: n_k_c_1 %d --> dk_c_1 %+0.6f',k_c_1_start,k_c_1_final,n_k_c_1,dk_c_1));
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
parameter = struct('type','parameter');
parameter.flag_verbose = flag_verbose;
parameter.n_order = n_order;
parameter.flag_d = flag_d;
parameter.flag_dd = flag_dd;
[ ...
 parameter ...
,n_part_order ...
,k_s_1_rescale_ ...
,k_s_1_rescale_start_ ...
,k_s_1_rescale_shift_ ...
,node_x_ ...
,n_nodex ...
,index_col_1_xs__ ...
,index_row_xs__ ...
,weight_denominator_x_ ...
,weight_numerator_1_xs__ ...
,weight_dnumerator_1_xs__ ...
,weight_ddnumerator_1_xs__ ...
] = ...
scatter_from_tensor_helper_n_9( ...
 parameter ...
,n_scatter ...
,k_s_1_ ...
,k_c_1_start ...
,dk_c_1 ...
);
%%%%%%%%;
index_col_1_xs__ = periodize(index_col_1_xs__,0,n_k_c_1);
%%%%%%%%;
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figsml;
plot3(index_col_0_xs__(:),index_col_1_xs__(:),index_col_2_xs__(:),'.');
xlabel('index_col_0_xs__','Interpreter','none');
ylabel('index_col_1_xs__','Interpreter','none');
axis equal; axis vis3d;
end;%if (flag_disp>1);
%%%%%%%%;
index_row_xxs__ = repmat(0:n_scatter-1,[n_nodex^2,1]);

%%%%;
if (flag_verbose>1);
%%%%;
[nc_0_min,tmp_ij] = min(index_col_0_xs__(:)); tmp_index = tmp_ij-1;
tmp_index_x = mod(tmp_index,n_nodex);
tmp_index_s = (tmp_index-tmp_index_x)/n_nodex;
disp(sprintf(' %% tmp_index %d (%d,%d) <-- nc_0_min %d',tmp_index,tmp_index_x,tmp_index_s,nc_0_min));
disp(sprintf(' %% index_col_0_xs__(:,1+tmp_index_s): %s',num2str(transpose(index_col_0_xs__(:,1+tmp_index_s)),' %+d')));
disp(sprintf(' %% index_col_1_xs__(:,1+tmp_index_s): %s',num2str(transpose(index_col_1_xs__(:,1+tmp_index_s)),' %+d')));
%%%%;
disp('returning'); return;
end;%if (flag_verbose>1);
%%%%;

%%%%%%%%;
weight_0_xs__ = bsxfun(@rdivide,weight_numerator_0_xs__,reshape(weight_denominator_x_,[n_nodex,1]));
weight_1_xs__ = bsxfun(@rdivide,weight_numerator_1_xs__,reshape(weight_denominator_x_,[n_nodex,1]));
weight_xxs__ = ...
reshape( ...
bsxfun(@times ...
,reshape(weight_0_xs__,[n_nodex,1,1,n_scatter]) ...
,reshape(weight_1_xs__,[1,n_nodex,1,n_scatter]) ...
) ...
,[n_nodex^2,n_scatter] ...
);
index_col_xxs__ =  ...
reshape( ...
bsxfun(@plus ...
,reshape(index_col_0_xs__,[n_nodex,1,1,n_scatter]) ...
,reshape(index_col_1_xs__,[1,n_nodex,1,n_scatter])*n_k_c_0 ...
) ...
,[n_nodex^2,n_scatter] ...
);
scatter_from_tensor_s01__ = sparse(1+index_row_xxs__(:),1+index_col_xxs__(:),weight_xxs__(:),n_scatter,n_k_c_0*n_k_c_1);
%%%%%%%%;
if flag_d;
dweight_0_xs__ = bsxfun(@rdivide,weight_dnumerator_0_xs__,reshape(weight_denominator_x_,[n_nodex,1])) / dk_c_0 ;
dweight_1_xs__ = bsxfun(@rdivide,weight_dnumerator_1_xs__,reshape(weight_denominator_x_,[n_nodex,1])) / dk_c_1 ;
dweightd0_xxs__ = ...
reshape( ...
bsxfun(@times ...
,reshape(dweight_0_xs__,[n_nodex,1,1,n_scatter]) ...
,reshape( weight_1_xs__,[1,n_nodex,1,n_scatter]) ...
) ...
,[n_nodex^2,n_scatter] ...
);
dweightd1_xxs__ = ...
reshape( ...
bsxfun(@times ...
,reshape( weight_0_xs__,[n_nodex,1,1,n_scatter]) ...
,reshape(dweight_1_xs__,[1,n_nodex,1,n_scatter]) ...
) ...
,[n_nodex^2,n_scatter] ...
);
index_col_xxs__ =  ...
reshape( ...
bsxfun(@plus ...
,reshape(index_col_0_xs__,[n_nodex,1,1,n_scatter]) ...
,reshape(index_col_1_xs__,[1,n_nodex,1,n_scatter])*n_k_c_0 ...
) ...
,[n_nodex^2,n_scatter] ...
);
dscatterd0_from_tensor_s01__ = sparse(1+index_row_xxs__(:),1+index_col_xxs__(:),dweightd0_xxs__(:),n_scatter,n_k_c_0*n_k_c_1);
dscatterd1_from_tensor_s01__ = sparse(1+index_row_xxs__(:),1+index_col_xxs__(:),dweightd1_xxs__(:),n_scatter,n_k_c_0*n_k_c_1);
end;%if flag_d;
%%%%%%%%;
if flag_dd;
ddweight_0_xs__ = bsxfun(@rdivide,weight_ddnumerator_0_xs__,reshape(weight_denominator_x_,[n_nodex,1])) / dk_c_0^2 ;
ddweight_1_xs__ = bsxfun(@rdivide,weight_ddnumerator_1_xs__,reshape(weight_denominator_x_,[n_nodex,1])) / dk_c_1^2 ;
ddweightd00_xxs__ = ...
reshape( ...
bsxfun(@times ...
,reshape(ddweight_0_xs__,[n_nodex,1,1,n_scatter]) ...
,reshape(  weight_1_xs__,[1,n_nodex,1,n_scatter]) ...
) ...
,[n_nodex^2,n_scatter] ...
);
ddweightd01_xxs__ = ...
reshape( ...
bsxfun(@times ...
,reshape( dweight_0_xs__,[n_nodex,1,1,n_scatter]) ...
,reshape( dweight_1_xs__,[1,n_nodex,1,n_scatter]) ...
) ...
,[n_nodex^2,n_scatter] ...
);
ddweightd11_xxs__ = ...
reshape( ...
bsxfun(@times ...
,reshape(  weight_0_xs__,[n_nodex,1,1,n_scatter]) ...
,reshape(ddweight_1_xs__,[1,n_nodex,1,n_scatter]) ...
) ...
,[n_nodex^2,n_scatter] ...
);
index_col_xxs__ =  ...
reshape( ...
bsxfun(@plus ...
,reshape(index_col_0_xs__,[n_nodex,1,1,n_scatter]) ...
,reshape(index_col_1_xs__,[1,n_nodex,1,n_scatter])*n_k_c_0 ...
) ...
,[n_nodex^2,n_scatter] ...
);
ddscatterd00_from_tensor_s01__ = sparse(1+index_row_xxs__(:),1+index_col_xxs__(:),ddweightd00_xxs__(:),n_scatter,n_k_c_0*n_k_c_1);
ddscatterd01_from_tensor_s01__ = sparse(1+index_row_xxs__(:),1+index_col_xxs__(:),ddweightd01_xxs__(:),n_scatter,n_k_c_0*n_k_c_1);
ddscatterd11_from_tensor_s01__ = sparse(1+index_row_xxs__(:),1+index_col_xxs__(:),ddweightd11_xxs__(:),n_scatter,n_k_c_0*n_k_c_1);
end;%if flag_dd;
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;


