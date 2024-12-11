function ...
[ ...
 scatter_from_tensor_s012__ ...
 dscatterd0_from_tensor_s012__ ...
 dscatterd1_from_tensor_s012__ ...
 dscatterd2_from_tensor_s012__ ...
 ddscatterd00_from_tensor_s012__ ...
 ddscatterd01_from_tensor_s012__ ...
 ddscatterd02_from_tensor_s012__ ...
 ddscatterd11_from_tensor_s012__ ...
 ddscatterd12_from_tensor_s012__ ...
 ddscatterd22_from_tensor_s012__ ...
] = ...
volume_k_c_scatter_from_tensor_interpolate_n_8( ...
 n_order ...
,n_k_c_0 ...
,n_k_c_1 ...
,n_k_c_2 ...
,k_c_0_lim_ ...
,k_c_1_lim_ ...
,k_c_2_lim_ ...
,n_scatter ...
,k_s_0_ ...
,k_s_1_ ...
,k_s_2_ ...
);
%%%%%%%%;
% Here we limit the scattered points to lie in the box defined by the tensor_grid. ;
%%%%%%%%;
% returns sparse matrix encoding the n_order interpolation operator between a tensor_grid ;
% using k_c_?_. ;
%%%%%%%%;
% We presume these points are associated with an (n_k_c_0)-by-(n_k_c_1)-by-(n_k_c_2) array: ;
% a_k_c_012___. ;
%%%%;
% The n_scatter points to be interpolated have coordinates stored in ;
% k_s_0_, ;
% k_s_1_, ;
% k_s_2_, ;
% ;
% The (n_scatter)-by-(n_k_0*n_k_1*n_k_2) scatter_from_tensor_s012__ matrix stores the ;
% interpolation weights linking a_k_c_012_ to the unrolled a_k_s_. ;
%%%%;

str_thisfunction = 'volume_k_c_scatter_from_tensor_interpolate_n_8';

if (nargin<6);
rng(1);
%%%%%%%%;
flag_verbose = 1;
flag_disp = 1; nf=0;
disp(sprintf(' %% testing %s',str_thisfunction));
n_k_c_0 = 64+0;
n_k_c_1 = 64+1;
n_k_c_2 = 64+2;
k_p_r_max = 48/(2*pi);
k_c_0_lim_ = 1.0*k_p_r_max*[-1,+1];
k_c_1_lim_ = 1.5*k_p_r_max*[-1,+1];
k_c_2_lim_ = 2.0*k_p_r_max*[-1,+1];
k_c_0_ = transpose(linspace(k_c_0_lim_(1),k_c_0_lim_(2),n_k_c_0));
k_c_1_ = transpose(linspace(k_c_1_lim_(1),k_c_1_lim_(2),n_k_c_1));
k_c_2_ = transpose(linspace(k_c_2_lim_(1),k_c_2_lim_(2),n_k_c_2));
[k_c_0_012___,k_c_1_012___,k_c_2_012___] = ndgrid(k_c_0_,k_c_1_,k_c_2_);
delta_ = [0.3;1.5;2.0];
f_c_012___ = exp(+i*(k_c_0_012___*delta_(1+0) + k_c_1_012___*delta_(1+1) + k_c_2_012___*delta_(1+2)));
dfd0_c_012___ = (+i*delta_(1+0))*exp(+i*(k_c_0_012___*delta_(1+0) + k_c_1_012___*delta_(1+1) + k_c_2_012___*delta_(1+2)));
dfd1_c_012___ = (+i*delta_(1+1))*exp(+i*(k_c_0_012___*delta_(1+0) + k_c_1_012___*delta_(1+1) + k_c_2_012___*delta_(1+2)));
dfd2_c_012___ = (+i*delta_(1+2))*exp(+i*(k_c_0_012___*delta_(1+0) + k_c_1_012___*delta_(1+1) + k_c_2_012___*delta_(1+2)));
ddfd00_c_012___ = (+i*delta_(1+0))*(+i*delta_(1+0))*exp(+i*(k_c_0_012___*delta_(1+0) + k_c_1_012___*delta_(1+1) + k_c_2_012___*delta_(1+2)));
ddfd01_c_012___ = (+i*delta_(1+0))*(+i*delta_(1+1))*exp(+i*(k_c_0_012___*delta_(1+0) + k_c_1_012___*delta_(1+1) + k_c_2_012___*delta_(1+2)));
ddfd02_c_012___ = (+i*delta_(1+0))*(+i*delta_(1+2))*exp(+i*(k_c_0_012___*delta_(1+0) + k_c_1_012___*delta_(1+1) + k_c_2_012___*delta_(1+2)));
ddfd11_c_012___ = (+i*delta_(1+1))*(+i*delta_(1+1))*exp(+i*(k_c_0_012___*delta_(1+0) + k_c_1_012___*delta_(1+1) + k_c_2_012___*delta_(1+2)));
ddfd12_c_012___ = (+i*delta_(1+1))*(+i*delta_(1+2))*exp(+i*(k_c_0_012___*delta_(1+0) + k_c_1_012___*delta_(1+1) + k_c_2_012___*delta_(1+2)));
ddfd22_c_012___ = (+i*delta_(1+2))*(+i*delta_(1+2))*exp(+i*(k_c_0_012___*delta_(1+0) + k_c_1_012___*delta_(1+1) + k_c_2_012___*delta_(1+2)));
n_scatter = 1024;
k_s_0_ = 0.5*k_p_r_max*(2*rand(n_scatter,1)-1);
k_s_1_ = 0.5*k_p_r_max*(2*rand(n_scatter,1)-1);
k_s_2_ = 0.5*k_p_r_max*(2*rand(n_scatter,1)-1);
n_order = 5;
[ ...
 scatter_from_tensor_s012__ ...
 dscatterd0_from_tensor_s012__ ...
 dscatterd1_from_tensor_s012__ ...
 dscatterd2_from_tensor_s012__ ...
 ddscatterd00_from_tensor_s012__ ...
 ddscatterd01_from_tensor_s012__ ...
 ddscatterd02_from_tensor_s012__ ...
 ddscatterd11_from_tensor_s012__ ...
 ddscatterd12_from_tensor_s012__ ...
 ddscatterd22_from_tensor_s012__ ...
] = ...
volume_k_c_scatter_from_tensor_interpolate_n_8( ...
 n_order ...
,n_k_c_0 ...
,n_k_c_1 ...
,n_k_c_2 ...
,k_c_0_lim_ ...
,k_c_1_lim_ ...
,k_c_2_lim_ ...
,n_scatter ...
,k_s_0_ ...
,k_s_1_ ...
,k_s_2_ ...
);
f_s_form_ = exp(+i*(k_s_0_*delta_(1+0) + k_s_1_*delta_(1+1) + k_s_2_*delta_(1+2)));
dfd0_s_form_ = (+i*delta_(1+0))*exp(+i*(k_s_0_*delta_(1+0) + k_s_1_*delta_(1+1) + k_s_2_*delta_(1+2)));
dfd1_s_form_ = (+i*delta_(1+1))*exp(+i*(k_s_0_*delta_(1+0) + k_s_1_*delta_(1+1) + k_s_2_*delta_(1+2)));
dfd2_s_form_ = (+i*delta_(1+2))*exp(+i*(k_s_0_*delta_(1+0) + k_s_1_*delta_(1+1) + k_s_2_*delta_(1+2)));
ddfd00_s_form_ = (+i*delta_(1+0))*(+i*delta_(1+0))*exp(+i*(k_s_0_*delta_(1+0) + k_s_1_*delta_(1+1) + k_s_2_*delta_(1+2)));
ddfd01_s_form_ = (+i*delta_(1+0))*(+i*delta_(1+1))*exp(+i*(k_s_0_*delta_(1+0) + k_s_1_*delta_(1+1) + k_s_2_*delta_(1+2)));
ddfd02_s_form_ = (+i*delta_(1+0))*(+i*delta_(1+2))*exp(+i*(k_s_0_*delta_(1+0) + k_s_1_*delta_(1+1) + k_s_2_*delta_(1+2)));
ddfd11_s_form_ = (+i*delta_(1+1))*(+i*delta_(1+1))*exp(+i*(k_s_0_*delta_(1+0) + k_s_1_*delta_(1+1) + k_s_2_*delta_(1+2)));
ddfd12_s_form_ = (+i*delta_(1+1))*(+i*delta_(1+2))*exp(+i*(k_s_0_*delta_(1+0) + k_s_1_*delta_(1+1) + k_s_2_*delta_(1+2)));
ddfd22_s_form_ = (+i*delta_(1+2))*(+i*delta_(1+2))*exp(+i*(k_s_0_*delta_(1+0) + k_s_1_*delta_(1+1) + k_s_2_*delta_(1+2)));
f_s_inte_ = scatter_from_tensor_s012__*f_c_012___(:);
dfd0_s_inte_ = dscatterd0_from_tensor_s012__*f_c_012___(:);
dfd1_s_inte_ = dscatterd1_from_tensor_s012__*f_c_012___(:);
dfd2_s_inte_ = dscatterd2_from_tensor_s012__*f_c_012___(:);
ddfd00_s_inte_ = ddscatterd00_from_tensor_s012__*f_c_012___(:);
ddfd01_s_inte_ = ddscatterd01_from_tensor_s012__*f_c_012___(:);
ddfd02_s_inte_ = ddscatterd02_from_tensor_s012__*f_c_012___(:);
ddfd11_s_inte_ = ddscatterd11_from_tensor_s012__*f_c_012___(:);
ddfd12_s_inte_ = ddscatterd12_from_tensor_s012__*f_c_012___(:);
ddfd22_s_inte_ = ddscatterd22_from_tensor_s012__*f_c_012___(:);
%%%%;
figure(1+nf);clf;figmed;
p_row = 2; p_col = 5; np=0;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
hold on;
plot([-1,+1],[-1,+1],'k-');
plot(real(f_s_form_),real(f_s_inte_),'ro');
plot(imag(f_s_form_),imag(f_s_inte_),'ro');
hold off;
xlabel('f_s_form_','Interpreter','none');
ylabel('f_s_inte_','Interpreter','none');
grid on;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
hold on;
plot([-1,+1],[-1,+1],'k-');
plot(real(dfd0_s_form_),real(dfd0_s_inte_),'ro');
plot(imag(dfd0_s_form_),imag(dfd0_s_inte_),'ro');
hold off;
xlabel('dfd0_s_form_','Interpreter','none');
ylabel('dfd0_s_inte_','Interpreter','none');
grid on;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
hold on;
plot([-1,+1],[-1,+1],'k-');
plot(real(dfd1_s_form_),real(dfd1_s_inte_),'ro');
plot(imag(dfd1_s_form_),imag(dfd1_s_inte_),'ro');
hold off;
xlabel('dfd1_s_form_','Interpreter','none');
ylabel('dfd1_s_inte_','Interpreter','none');
grid on;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
hold on;
plot([-1,+1],[-1,+1],'k-');
plot(real(dfd2_s_form_),real(dfd2_s_inte_),'ro');
plot(imag(dfd2_s_form_),imag(dfd2_s_inte_),'ro');
hold off;
xlabel('dfd2_s_form_','interpreter','none');
ylabel('dfd2_s_inte_','interpreter','none');
grid on;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
hold on;
plot([-1,+1],[-1,+1],'k-');
plot(real(ddfd00_s_form_),real(ddfd00_s_inte_),'ro');
plot(imag(ddfd00_s_form_),imag(ddfd00_s_inte_),'ro');
hold off;
xlabel('ddfd00_s_form_','interpreter','none');
ylabel('ddfd00_s_inte_','interpreter','none');
grid on;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
hold on;
plot([-1,+1],[-1,+1],'k-');
plot(real(ddfd01_s_form_),real(ddfd01_s_inte_),'ro');
plot(imag(ddfd01_s_form_),imag(ddfd01_s_inte_),'ro');
hold off;
xlabel('ddfd01_s_form_','interpreter','none');
ylabel('ddfd01_s_inte_','interpreter','none');
grid on;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
hold on;
plot([-1,+1],[-1,+1],'k-');
plot(real(ddfd11_s_form_),real(ddfd11_s_inte_),'ro');
plot(imag(ddfd11_s_form_),imag(ddfd11_s_inte_),'ro');
hold off;
xlabel('ddfd11_s_form_','interpreter','none');
ylabel('ddfd11_s_inte_','interpreter','none');
grid on;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
hold on;
plot([-1,+1],[-1,+1],'k-');
plot(real(ddfd11_s_form_),real(ddfd11_s_inte_),'ro');
plot(imag(ddfd11_s_form_),imag(ddfd11_s_inte_),'ro');
hold off;
xlabel('ddfd11_s_form_','interpreter','none');
ylabel('ddfd11_s_inte_','interpreter','none');
grid on;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
hold on;
plot([-1,+1],[-1,+1],'k-');
plot(real(ddfd12_s_form_),real(ddfd12_s_inte_),'ro');
plot(imag(ddfd12_s_form_),imag(ddfd12_s_inte_),'ro');
hold off;
xlabel('ddfd12_s_form_','interpreter','none');
ylabel('ddfd12_s_inte_','interpreter','none');
grid on;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
hold on;
plot([-1,+1],[-1,+1],'k-');
plot(real(ddfd22_s_form_),real(ddfd22_s_inte_),'ro');
plot(imag(ddfd22_s_form_),imag(ddfd22_s_inte_),'ro');
hold off;
xlabel('ddfd22_s_form_','interpreter','none');
ylabel('ddfd22_s_inte_','interpreter','none');
grid on;
%%%%;
fnorm_disp(flag_verbose,'f_s_form_',f_s_form_,'f_s_inte_',f_s_inte_);
fnorm_disp(flag_verbose,'dfd0_s_form_',dfd0_s_form_,'dfd0_s_inte_',dfd0_s_inte_);
fnorm_disp(flag_verbose,'dfd1_s_form_',dfd1_s_form_,'dfd1_s_inte_',dfd1_s_inte_);
fnorm_disp(flag_verbose,'dfd2_s_form_',dfd2_s_form_,'dfd2_s_inte_',dfd2_s_inte_);
fnorm_disp(flag_verbose,'ddfd00_s_form_',ddfd00_s_form_,'ddfd00_s_inte_',ddfd00_s_inte_);
fnorm_disp(flag_verbose,'ddfd01_s_form_',ddfd01_s_form_,'ddfd01_s_inte_',ddfd01_s_inte_);
fnorm_disp(flag_verbose,'ddfd02_s_form_',ddfd02_s_form_,'ddfd02_s_inte_',ddfd02_s_inte_);
fnorm_disp(flag_verbose,'ddfd11_s_form_',ddfd11_s_form_,'ddfd11_s_inte_',ddfd11_s_inte_);
fnorm_disp(flag_verbose,'ddfd12_s_form_',ddfd12_s_form_,'ddfd12_s_inte_',ddfd12_s_inte_);
fnorm_disp(flag_verbose,'ddfd22_s_form_',ddfd22_s_form_,'ddfd22_s_inte_',ddfd22_s_inte_);
%%%%%%%%;
disp(sprintf(' %% returning')); return;
end;%if (nargin<6);

flag_verbose = 0;
flag_disp = 0; nf=0;
flag_check = 0;
flag_d = (nargout>=2);
flag_dd = (nargout>=5);

if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

na=0;
if (nargin<1+na); n_order=[]; end; na=na+1;
if (nargin<1+na); n_k_c_0=[]; end; na=na+1;
if (nargin<1+na); n_k_c_1=[]; end; na=na+1;
if (nargin<1+na); n_k_c_2=[]; end; na=na+1;
if (nargin<1+na); k_c_0_lim_=[]; end; na=na+1;
if (nargin<1+na); k_c_1_lim_=[]; end; na=na+1;
if (nargin<1+na); k_c_2_lim_=[]; end; na=na+1;
if (nargin<1+na); n_scatter=[]; end; na=na+1;
if (nargin<1+na); k_s_0_=[]; end; na=na+1;
if (nargin<1+na); k_s_1_=[]; end; na=na+1;
if (nargin<1+na); k_s_2_=[]; end; na=na+1;

%%%%%%%%;
n_part_order = round((n_order-1)/2);
if (flag_verbose>0); disp(sprintf(' %% n_order %d n_part_order %d',n_order,n_part_order)); end;
%%%%;
if (n_k_c_0<=2*n_part_order+1); disp(sprintf(' %% warning, n_k_c_0 %d, 2*n_part_order+1 %d in %s',n_k_c_0,2*n_part_order+1,str_thisfunction)); end;
if (n_k_c_1<=2*n_part_order+1); disp(sprintf(' %% warning, n_k_c_1 %d, 2*n_part_order+1 %d in %s',n_k_c_1,2*n_part_order+1,str_thisfunction)); end;
if (n_k_c_2<=2*n_part_order+1); disp(sprintf(' %% warning, n_k_c_2 %d, 2*n_part_order+1 %d in %s',n_k_c_2,2*n_part_order+1,str_thisfunction)); end;
%%%%;
dk_c_0 = diff(k_c_0_lim_)/max(1,n_k_c_0-1);
dk_c_1 = diff(k_c_1_lim_)/max(1,n_k_c_1-1);
dk_c_2 = diff(k_c_2_lim_)/max(1,n_k_c_2-1);
%%%%;
k_c_0_lim_use_ = [k_c_0_lim_(1)+n_part_order*dk_c_0,k_c_0_lim_(2)-n_part_order*dk_c_0];
k_c_1_lim_use_ = [k_c_1_lim_(1)+n_part_order*dk_c_1,k_c_1_lim_(2)-n_part_order*dk_c_1];
k_c_2_lim_use_ = [k_c_2_lim_(1)+n_part_order*dk_c_2,k_c_2_lim_(2)-n_part_order*dk_c_2];
%%%%;
if numel(efind(k_s_0_<=min(k_c_0_lim_use_)))>0; disp(sprintf(' %% Warning, k_s_0_ out of bounds in %s',str_thisfunction)); end;
if numel(efind(k_s_0_>=max(k_c_0_lim_use_)))>0; disp(sprintf(' %% Warning, k_s_0_ out of bounds in %s',str_thisfunction)); end;
if numel(efind(k_s_1_<=min(k_c_1_lim_use_)))>0; disp(sprintf(' %% Warning, k_s_1_ out of bounds in %s',str_thisfunction)); end;
if numel(efind(k_s_1_>=max(k_c_1_lim_use_)))>0; disp(sprintf(' %% Warning, k_s_1_ out of bounds in %s',str_thisfunction)); end;
if numel(efind(k_s_2_<=min(k_c_2_lim_use_)))>0; disp(sprintf(' %% Warning, k_s_2_ out of bounds in %s',str_thisfunction)); end;
if numel(efind(k_s_2_>=max(k_c_2_lim_use_)))>0; disp(sprintf(' %% Warning, k_s_2_ out of bounds in %s',str_thisfunction)); end;
k_s_0_ = max(min(k_c_0_lim_use_),min(max(k_c_0_lim_use_),k_s_0_));
k_s_1_ = max(min(k_c_1_lim_use_),min(max(k_c_1_lim_use_),k_s_1_));
k_s_2_ = max(min(k_c_2_lim_use_),min(max(k_c_2_lim_use_),k_s_2_));
%%%%;
k_c_0_start = k_c_0_lim_(1); k_c_0_final = k_c_0_lim_(2);
k_c_1_start = k_c_1_lim_(1); k_c_1_final = k_c_1_lim_(2);
k_c_2_start = k_c_2_lim_(1); k_c_2_final = k_c_2_lim_(2);
%%%%;
if (flag_verbose>0);
disp(sprintf(' %% k_c_0_start %+0.6f k_c_0_final %+0.6f: n_k_c_0 %d --> dk_c_0 %+0.6f',k_c_0_start,k_c_0_final,n_k_c_0,dk_c_0));
disp(sprintf(' %% k_c_1_start %+0.6f k_c_1_final %+0.6f: n_k_c_1 %d --> dk_c_1 %+0.6f',k_c_1_start,k_c_1_final,n_k_c_1,dk_c_1));
disp(sprintf(' %% k_c_2_start %+0.6f k_c_2_final %+0.6f: n_k_c_2 %d --> dk_c_2 %+0.6f',k_c_2_start,k_c_2_final,n_k_c_2,dk_c_2));
end;%if (flag_verbose>0);
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% calculate scatter_from_tensor_sba__. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
k_s_0_rescale_ = (k_s_0_ - k_c_0_start)/max(1e-12,dk_c_0);
k_s_1_rescale_ = (k_s_1_ - k_c_1_start)/max(1e-12,dk_c_1);
k_s_2_rescale_ = (k_s_2_ - k_c_2_start)/max(1e-12,dk_c_2);
%%%%;
if (flag_verbose>0);
disp(sprintf(' %% k_s_0_rescale_ in [%+0.6f , %+0.6f]',prctile(k_s_0_rescale_,[  0,100])));
disp(sprintf(' %% k_s_1_rescale_ in [%+0.6f , %+0.6f]',prctile(k_s_1_rescale_,[  0,100])));
disp(sprintf(' %% k_s_2_rescale_ in [%+0.6f , %+0.6f]',prctile(k_s_2_rescale_,[  0,100])));
end;%if (flag_verbose>0);
%%%%;
k_s_0_rescale_start_ = floor(k_s_0_rescale_) - n_part_order;
k_s_1_rescale_start_ = floor(k_s_1_rescale_) - n_part_order;
k_s_2_rescale_start_ = floor(k_s_2_rescale_) - n_part_order;
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
k_s_1_rescale_shift_ = k_s_1_rescale_ - k_s_1_rescale_start_ - n_part_order - 0.5 ;
k_s_2_rescale_shift_ = k_s_2_rescale_ - k_s_2_rescale_start_ - n_part_order - 0.5 ;
%%%%;
weight_numerator_0_xs__ = zeros(n_nodex,n_scatter);
weight_numerator_1_xs__ = zeros(n_nodex,n_scatter);
weight_numerator_2_xs__ = zeros(n_nodex,n_scatter);
index_col_0_xs__ = zeros(n_nodex,n_scatter);
index_col_1_xs__ = zeros(n_nodex,n_scatter);
index_col_2_xs__ = zeros(n_nodex,n_scatter);
index_row_xxxs__ = repmat(0:n_scatter-1,[n_nodex^3,1]);
%%%%;
tmp_t = tic();
weight_numerator_0_xxs___ = bsxfun(@minus,reshape(k_s_0_rescale_shift_,[1,1,n_scatter]),repmat(node_x_,[1,n_nodex,1]));
for nnodex=0:n_nodex-1; weight_numerator_0_xxs___(1+nnodex,1+nnodex,:) = 1; end;
weight_numerator_0_xs__ = reshape(prod(weight_numerator_0_xxs___,1),[n_nodex,n_scatter]);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% weight_numerator_0_xs__: %0.6fs',tmp_t)); end;
tmp_t = tic();
weight_numerator_1_xxs___ = bsxfun(@minus,reshape(k_s_1_rescale_shift_,[1,1,n_scatter]),repmat(node_x_,[1,n_nodex,1]));
for nnodex=0:n_nodex-1; weight_numerator_1_xxs___(1+nnodex,1+nnodex,:) = 1; end;
weight_numerator_1_xs__ = reshape(prod(weight_numerator_1_xxs___,1),[n_nodex,n_scatter]);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% weight_numerator_1_xs__: %0.6fs',tmp_t)); end;
tmp_t = tic();
weight_numerator_2_xxs___ = bsxfun(@minus,reshape(k_s_2_rescale_shift_,[1,1,n_scatter]),repmat(node_x_,[1,n_nodex,1]));
for nnodex=0:n_nodex-1; weight_numerator_2_xxs___(1+nnodex,1+nnodex,:) = 1; end;
weight_numerator_2_xs__ = reshape(prod(weight_numerator_2_xxs___,1),[n_nodex,n_scatter]);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% weight_numerator_2_xs__: %0.6fs',tmp_t)); end;
%%%%;
%%%%;
if flag_d;
tmp_t = tic();
weight_dnumerator_0_xxxs____ = bsxfun(@minus,reshape(k_s_0_rescale_shift_,[1,1,1,n_scatter]),repmat(node_x_,[1,n_nodex,n_nodex,1]));
for nnodex_a=0:n_nodex-1;
weight_dnumerator_0_xxxs____(1+nnodex_a,1+nnodex_a,:,:) = 1;
end;%for nnodex_a=0:n_nodex-1;
for nnodex_b=0:n_nodex-1;
weight_dnumerator_0_xxxs____(1+nnodex_b,:,1+nnodex_b,:) = 1;
end;%for nnodex_b=0:n_nodex-1;
for nnodex_ab=0:n_nodex-1;
weight_dnumerator_0_xxxs____(1+nnodex_ab,1+nnodex_ab,1+nnodex_ab,:) = 0;
end;%for nnodex_ab=0:n_nodex-1;
weight_dnumerator_0_xs__ = reshape(sum(prod(weight_dnumerator_0_xxxs____,[1]),[2]),[n_nodex,n_scatter]);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% weight_dnumerator_0_xs__: %0.6fs',tmp_t)); end;
%%%%;
tmp_t = tic();
weight_dnumerator_1_xxxs____ = bsxfun(@minus,reshape(k_s_1_rescale_shift_,[1,1,1,n_scatter]),repmat(node_x_,[1,n_nodex,n_nodex,1]));
for nnodex_a=0:n_nodex-1;
weight_dnumerator_1_xxxs____(1+nnodex_a,1+nnodex_a,:,:) = 1;
end;%for nnodex_a=0:n_nodex-1;
for nnodex_b=0:n_nodex-1;
weight_dnumerator_1_xxxs____(1+nnodex_b,:,1+nnodex_b,:) = 1;
end;%for nnodex_b=0:n_nodex-1;
for nnodex_ab=0:n_nodex-1;
weight_dnumerator_1_xxxs____(1+nnodex_ab,1+nnodex_ab,1+nnodex_ab,:) = 0;
end;%for nnodex_ab=0:n_nodex-1;
weight_dnumerator_1_xs__ = reshape(sum(prod(weight_dnumerator_1_xxxs____,[1]),[2]),[n_nodex,n_scatter]);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% weight_dnumerator_1_xs__: %0.6fs',tmp_t)); end;
%%%%;
tmp_t = tic();
weight_dnumerator_2_xxxs____ = bsxfun(@minus,reshape(k_s_2_rescale_shift_,[1,1,1,n_scatter]),repmat(node_x_,[1,n_nodex,n_nodex,1]));
for nnodex_a=0:n_nodex-1;
weight_dnumerator_2_xxxs____(1+nnodex_a,1+nnodex_a,:,:) = 1;
end;%for nnodex_a=0:n_nodex-1;
for nnodex_b=0:n_nodex-1;
weight_dnumerator_2_xxxs____(1+nnodex_b,:,1+nnodex_b,:) = 1;
end;%for nnodex_b=0:n_nodex-1;
for nnodex_ab=0:n_nodex-1;
weight_dnumerator_2_xxxs____(1+nnodex_ab,1+nnodex_ab,1+nnodex_ab,:) = 0;
end;%for nnodex_ab=0:n_nodex-1;
weight_dnumerator_2_xs__ = reshape(sum(prod(weight_dnumerator_2_xxxs____,[1]),[2]),[n_nodex,n_scatter]);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% weight_dnumerator_2_xs__: %0.6fs',tmp_t)); end;
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
%%%%;
tmp_t = tic();
weight_ddnumerator_1_xxxxs_____ = bsxfun(@minus,reshape(k_s_1_rescale_shift_,[1,1,1,1,n_scatter]),repmat(node_x_,[1,n_nodex,n_nodex,n_nodex,1]));
for nnodex_A=0:n_nodex-1;
weight_ddnumerator_1_xxxxs_____(1+nnodex_A,1+nnodex_A,:,:,:) = 1;
end;%for nnodex_A=0:n_nodex-1;
for nnodex_B=0:n_nodex-1;
weight_ddnumerator_1_xxxxs_____(1+nnodex_B,:,1+nnodex_B,:,:) = 1;
end;%for nnodex_B=0:n_nodex-1;
for nnodex_C=0:n_nodex-1;
weight_ddnumerator_1_xxxxs_____(1+nnodex_C,:,:,1+nnodex_C,:) = 1;
end;%for nnodex_C=0:n_nodex-1;
for nnodex_AB=0:n_nodex-1;
weight_ddnumerator_1_xxxxs_____(1+nnodex_AB,1+nnodex_AB,1+nnodex_AB,:,:) = 0;
end;%for nnodex_AB=0:n_nodex-1;
for nnodex_AC=0:n_nodex-1;
weight_ddnumerator_1_xxxxs_____(1+nnodex_AC,1+nnodex_AC,:,1+nnodex_AC,:) = 0;
end;%for nnodex_AC=0:n_nodex-1;
for nnodex_BC=0:n_nodex-1;
weight_ddnumerator_1_xxxxs_____(1+nnodex_BC,:,1+nnodex_BC,1+nnodex_BC,:) = 0;
end;%for nnodex_BC=0:n_nodex-1;
for nnodex_ABC=0:n_nodex-1;
weight_ddnumerator_1_xxxxs_____(1+nnodex_ABC,1+nnodex_ABC,1+nnodex_ABC,1+nnodex_ABC,:) = 0;
end;%for nnodex_ABC=0:n_nodex-1;
weight_ddnumerator_1_xs__ = reshape(sum(prod(weight_ddnumerator_1_xxxxs_____,[1]),[2,3]),[n_nodex,n_scatter]);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% weight_ddnumerator_1_xs__: %0.6fs',tmp_t)); end;
%%%%;
tmp_t = tic();
weight_ddnumerator_2_xxxxs_____ = bsxfun(@minus,reshape(k_s_2_rescale_shift_,[1,1,1,1,n_scatter]),repmat(node_x_,[1,n_nodex,n_nodex,n_nodex,1]));
for nnodex_A=0:n_nodex-1;
weight_ddnumerator_2_xxxxs_____(1+nnodex_A,1+nnodex_A,:,:,:) = 1;
end;%for nnodex_A=0:n_nodex-1;
for nnodex_B=0:n_nodex-1;
weight_ddnumerator_2_xxxxs_____(1+nnodex_B,:,1+nnodex_B,:,:) = 1;
end;%for nnodex_B=0:n_nodex-1;
for nnodex_C=0:n_nodex-1;
weight_ddnumerator_2_xxxxs_____(1+nnodex_C,:,:,1+nnodex_C,:) = 1;
end;%for nnodex_C=0:n_nodex-1;
for nnodex_AB=0:n_nodex-1;
weight_ddnumerator_2_xxxxs_____(1+nnodex_AB,1+nnodex_AB,1+nnodex_AB,:,:) = 0;
end;%for nnodex_AB=0:n_nodex-1;
for nnodex_AC=0:n_nodex-1;
weight_ddnumerator_2_xxxxs_____(1+nnodex_AC,1+nnodex_AC,:,1+nnodex_AC,:) = 0;
end;%for nnodex_AC=0:n_nodex-1;
for nnodex_BC=0:n_nodex-1;
weight_ddnumerator_2_xxxxs_____(1+nnodex_BC,:,1+nnodex_BC,1+nnodex_BC,:) = 0;
end;%for nnodex_BC=0:n_nodex-1;
for nnodex_ABC=0:n_nodex-1;
weight_ddnumerator_2_xxxxs_____(1+nnodex_ABC,1+nnodex_ABC,1+nnodex_ABC,1+nnodex_ABC,:) = 0;
end;%for nnodex_ABC=0:n_nodex-1;
weight_ddnumerator_2_xs__ = reshape(sum(prod(weight_ddnumerator_2_xxxxs_____,[1]),[2,3]),[n_nodex,n_scatter]);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% weight_ddnumerator_2_xs__: %0.6fs',tmp_t)); end;
end;%if flag_dd;
%%%%%%%%;
index_col_0_xs__ = bsxfun(@plus,reshape(k_s_0_rescale_start_,[1,n_scatter]),transpose([0:n_nodex-1]));
index_col_1_xs__ = bsxfun(@plus,reshape(k_s_1_rescale_start_,[1,n_scatter]),transpose([0:n_nodex-1]));
index_col_2_xs__ = bsxfun(@plus,reshape(k_s_2_rescale_start_,[1,n_scatter]),transpose([0:n_nodex-1]));
%%%%;
if (flag_verbose>0); disp(sprintf(' %% n_k_c_0 %d: index_col_0_xs__ in [%.4d,%.4d]',n_k_c_0,min(index_col_0_xs__,[],'all'),max(index_col_0_xs__,[],'all'))); end;
if (flag_verbose>0); disp(sprintf(' %% n_k_c_1 %d: index_col_1_xs__ in [%.4d,%.4d]',n_k_c_1,min(index_col_1_xs__,[],'all'),max(index_col_1_xs__,[],'all'))); end;
if (flag_verbose>0); disp(sprintf(' %% n_k_c_2 %d: index_col_2_xs__ in [%.4d,%.4d]',n_k_c_2,min(index_col_2_xs__,[],'all'),max(index_col_2_xs__,[],'all'))); end;
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figsml;
plot3(index_col_0_xs__(:),index_col_1_xs__(:),index_col_2_xs__(:),'.');
xlabel('index_col_0_xs__','Interpreter','none');
ylabel('index_col_1_xs__','Interpreter','none');
zlabel('index_col_2_xs__','Interpreter','none');
axis equal; axis vis3d;
end;%if (flag_disp>1);
%%%%%%%%;

%%%%;
if (flag_verbose>1);
%%%%;
[nc_0_min,tmp_ij] = min(index_col_0_xs__(:)); tmp_index = tmp_ij-1;
tmp_index_x = mod(tmp_index,n_nodex);
tmp_index_s = (tmp_index-tmp_index_x)/n_nodex;
disp(sprintf(' %% tmp_index %d (%d,%d) <-- nc_0_min %d',tmp_index,tmp_index_x,tmp_index_s,nc_0_min));
disp(sprintf(' %% index_col_0_xs__(:,1+tmp_index_s): %s',num2str(transpose(index_col_0_xs__(:,1+tmp_index_s)),' %+d')));
disp(sprintf(' %% index_col_1_xs__(:,1+tmp_index_s): %s',num2str(transpose(index_col_1_xs__(:,1+tmp_index_s)),' %+d')));
disp(sprintf(' %% index_col_2_xs__(:,1+tmp_index_s): %s',num2str(transpose(index_col_2_xs__(:,1+tmp_index_s)),' %+d')));
%%%%;
disp('returning'); return;
end;%if (flag_verbose>1);
%%%%;

%%%%%%%%;
weight_0_xs__ = bsxfun(@rdivide,weight_numerator_0_xs__,reshape(weight_denominator_x_,[n_nodex,1]));
weight_1_xs__ = bsxfun(@rdivide,weight_numerator_1_xs__,reshape(weight_denominator_x_,[n_nodex,1]));
weight_2_xs__ = bsxfun(@rdivide,weight_numerator_2_xs__,reshape(weight_denominator_x_,[n_nodex,1]));
weight_xxxs__ = ...
reshape( ...
bsxfun(@times ...
,bsxfun(@times ...
,reshape(weight_0_xs__,[n_nodex,1,1,n_scatter]) ...
,reshape(weight_1_xs__,[1,n_nodex,1,n_scatter]) ...
) ...
,reshape(weight_2_xs__,[1,1,n_nodex,n_scatter]) ...
) ...
,[n_nodex^3,n_scatter] ...
);
index_col_xxxs__ =  ...
reshape( ...
bsxfun(@plus ...
,bsxfun(@plus ...
,reshape(index_col_0_xs__,[n_nodex,1,1,n_scatter]) ...
,reshape(index_col_1_xs__,[1,n_nodex,1,n_scatter])*n_k_c_0 ...
) ...
,reshape(index_col_2_xs__,[1,1,n_nodex,n_scatter])*n_k_c_0*n_k_c_1 ...
) ...
,[n_nodex^3,n_scatter] ...
);
scatter_from_tensor_s012__ = sparse(1+index_row_xxxs__(:),1+index_col_xxxs__(:),weight_xxxs__(:),n_scatter,n_k_c_0*n_k_c_1*n_k_c_2);
%%%%%%%%;
if flag_d;
dweight_0_xs__ = bsxfun(@rdivide,weight_dnumerator_0_xs__,reshape(weight_denominator_x_,[n_nodex,1])) / dk_c_0 ;
dweight_1_xs__ = bsxfun(@rdivide,weight_dnumerator_1_xs__,reshape(weight_denominator_x_,[n_nodex,1])) / dk_c_1 ;
dweight_2_xs__ = bsxfun(@rdivide,weight_dnumerator_2_xs__,reshape(weight_denominator_x_,[n_nodex,1])) / dk_c_2 ;
dweightd0_xxxs__ = ...
reshape( ...
bsxfun(@times ...
,bsxfun(@times ...
,reshape(dweight_0_xs__,[n_nodex,1,1,n_scatter]) ...
,reshape( weight_1_xs__,[1,n_nodex,1,n_scatter]) ...
) ...
,reshape( weight_2_xs__,[1,1,n_nodex,n_scatter]) ...
) ...
,[n_nodex^3,n_scatter] ...
);
dweightd1_xxxs__ = ...
reshape( ...
bsxfun(@times ...
,bsxfun(@times ...
,reshape( weight_0_xs__,[n_nodex,1,1,n_scatter]) ...
,reshape(dweight_1_xs__,[1,n_nodex,1,n_scatter]) ...
) ...
,reshape( weight_2_xs__,[1,1,n_nodex,n_scatter]) ...
) ...
,[n_nodex^3,n_scatter] ...
);
dweightd2_xxxs__ = ...
reshape( ...
bsxfun(@times ...
,bsxfun(@times ...
,reshape( weight_0_xs__,[n_nodex,1,1,n_scatter]) ...
,reshape( weight_1_xs__,[1,n_nodex,1,n_scatter]) ...
) ...
,reshape(dweight_2_xs__,[1,1,n_nodex,n_scatter]) ...
) ...
,[n_nodex^3,n_scatter] ...
);
index_col_xxxs__ =  ...
reshape( ...
bsxfun(@plus ...
,bsxfun(@plus ...
,reshape(index_col_0_xs__,[n_nodex,1,1,n_scatter]) ...
,reshape(index_col_1_xs__,[1,n_nodex,1,n_scatter])*n_k_c_0 ...
) ...
,reshape(index_col_2_xs__,[1,1,n_nodex,n_scatter])*n_k_c_0*n_k_c_1 ...
) ...
,[n_nodex^3,n_scatter] ...
);
dscatterd0_from_tensor_s012__ = sparse(1+index_row_xxxs__(:),1+index_col_xxxs__(:),dweightd0_xxxs__(:),n_scatter,n_k_c_0*n_k_c_1*n_k_c_2);
dscatterd1_from_tensor_s012__ = sparse(1+index_row_xxxs__(:),1+index_col_xxxs__(:),dweightd1_xxxs__(:),n_scatter,n_k_c_0*n_k_c_1*n_k_c_2);
dscatterd2_from_tensor_s012__ = sparse(1+index_row_xxxs__(:),1+index_col_xxxs__(:),dweightd2_xxxs__(:),n_scatter,n_k_c_0*n_k_c_1*n_k_c_2);
end;%if flag_d;
%%%%%%%%;
if flag_dd;
ddweight_0_xs__ = bsxfun(@rdivide,weight_ddnumerator_0_xs__,reshape(weight_denominator_x_,[n_nodex,1])) / dk_c_0^2 ;
ddweight_1_xs__ = bsxfun(@rdivide,weight_ddnumerator_1_xs__,reshape(weight_denominator_x_,[n_nodex,1])) / dk_c_1^2 ;
ddweight_2_xs__ = bsxfun(@rdivide,weight_ddnumerator_2_xs__,reshape(weight_denominator_x_,[n_nodex,1])) / dk_c_2^2 ;
ddweightd00_xxxs__ = ...
reshape( ...
bsxfun(@times ...
,bsxfun(@times ...
,reshape(ddweight_0_xs__,[n_nodex,1,1,n_scatter]) ...
,reshape(  weight_1_xs__,[1,n_nodex,1,n_scatter]) ...
) ...
,reshape(  weight_2_xs__,[1,1,n_nodex,n_scatter]) ...
) ...
,[n_nodex^3,n_scatter] ...
);
ddweightd01_xxxs__ = ...
reshape( ...
bsxfun(@times ...
,bsxfun(@times ...
,reshape( dweight_0_xs__,[n_nodex,1,1,n_scatter]) ...
,reshape( dweight_1_xs__,[1,n_nodex,1,n_scatter]) ...
) ...
,reshape(  weight_2_xs__,[1,1,n_nodex,n_scatter]) ...
) ...
,[n_nodex^3,n_scatter] ...
);
ddweightd02_xxxs__ = ...
reshape( ...
bsxfun(@times ...
,bsxfun(@times ...
,reshape( dweight_0_xs__,[n_nodex,1,1,n_scatter]) ...
,reshape(  weight_1_xs__,[1,n_nodex,1,n_scatter]) ...
) ...
,reshape( dweight_2_xs__,[1,1,n_nodex,n_scatter]) ...
) ...
,[n_nodex^3,n_scatter] ...
);
ddweightd11_xxxs__ = ...
reshape( ...
bsxfun(@times ...
,bsxfun(@times ...
,reshape(  weight_0_xs__,[n_nodex,1,1,n_scatter]) ...
,reshape(ddweight_1_xs__,[1,n_nodex,1,n_scatter]) ...
) ...
,reshape(  weight_2_xs__,[1,1,n_nodex,n_scatter]) ...
) ...
,[n_nodex^3,n_scatter] ...
);
ddweightd12_xxxs__ = ...
reshape( ...
bsxfun(@times ...
,bsxfun(@times ...
,reshape(  weight_0_xs__,[n_nodex,1,1,n_scatter]) ...
,reshape( dweight_1_xs__,[1,n_nodex,1,n_scatter]) ...
) ...
,reshape( dweight_2_xs__,[1,1,n_nodex,n_scatter]) ...
) ...
,[n_nodex^3,n_scatter] ...
);
ddweightd22_xxxs__ = ...
reshape( ...
bsxfun(@times ...
,bsxfun(@times ...
,reshape(  weight_0_xs__,[n_nodex,1,1,n_scatter]) ...
,reshape(  weight_1_xs__,[1,n_nodex,1,n_scatter]) ...
) ...
,reshape(ddweight_2_xs__,[1,1,n_nodex,n_scatter]) ...
) ...
,[n_nodex^3,n_scatter] ...
);
index_col_xxxs__ =  ...
reshape( ...
bsxfun(@plus ...
,bsxfun(@plus ...
,reshape(index_col_0_xs__,[n_nodex,1,1,n_scatter]) ...
,reshape(index_col_1_xs__,[1,n_nodex,1,n_scatter])*n_k_c_0 ...
) ...
,reshape(index_col_2_xs__,[1,1,n_nodex,n_scatter])*n_k_c_0*n_k_c_1 ...
) ...
,[n_nodex^3,n_scatter] ...
);
ddscatterd00_from_tensor_s012__ = sparse(1+index_row_xxxs__(:),1+index_col_xxxs__(:),ddweightd00_xxxs__(:),n_scatter,n_k_c_0*n_k_c_1*n_k_c_2);
ddscatterd01_from_tensor_s012__ = sparse(1+index_row_xxxs__(:),1+index_col_xxxs__(:),ddweightd01_xxxs__(:),n_scatter,n_k_c_0*n_k_c_1*n_k_c_2);
ddscatterd02_from_tensor_s012__ = sparse(1+index_row_xxxs__(:),1+index_col_xxxs__(:),ddweightd02_xxxs__(:),n_scatter,n_k_c_0*n_k_c_1*n_k_c_2);
ddscatterd11_from_tensor_s012__ = sparse(1+index_row_xxxs__(:),1+index_col_xxxs__(:),ddweightd11_xxxs__(:),n_scatter,n_k_c_0*n_k_c_1*n_k_c_2);
ddscatterd12_from_tensor_s012__ = sparse(1+index_row_xxxs__(:),1+index_col_xxxs__(:),ddweightd12_xxxs__(:),n_scatter,n_k_c_0*n_k_c_1*n_k_c_2);
ddscatterd22_from_tensor_s012__ = sparse(1+index_row_xxxs__(:),1+index_col_xxxs__(:),ddweightd22_xxxs__(:),n_scatter,n_k_c_0*n_k_c_1*n_k_c_2);
end;%if flag_dd;
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;


