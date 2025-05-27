function ...
[ ...
 scatter_from_tensor_s10__ ...
 dscatterd0_from_tensor_s10__ ...
 dscatterd1_from_tensor_s10__ ...
 ddscatterd00_from_tensor_s10__ ...
 ddscatterd01_from_tensor_s10__ ...
 ddscatterd11_from_tensor_s10__ ...
] = ...
cylinder_k_c_scatter_from_tensor_adaptive_interpolate_n_9( ...
 n_order ...
,n_k_c_0 ...
,n_k_c_1_ ...
,k_c_0_lim_ ...
,k_c_1_lim__ ...
,n_scatter ...
,k_s_0_ ...
,k_s_1_ ...
);
%%%%%%%%;
% We interpret the k_c_0_ array as a standard tensor array (e.g., height). ;
% We interpret the k_c_1_ array as periodic (e.g., angle). ; 
% However, the number of k_c_1_ points depends on nk_c_0. ;
% For this we store: ;
% n_k_c_1 :=n_k_c_1_(1+nk_c_0) ;
% k_c_1_lim_ := k_c_1_lim__(1+nk_c_0,:) ;
%%%%%%%%;
% returns sparse matrix encoding the n_order interpolation operator between a tensor_grid ;
% using k_c_?_. ;
%%%%%%%%;
% We presume these points are associated with an unrolled array: ;
% a_k_c_10_. ;
%%%%;
% The n_scatter points to be interpolated have coordinates stored in ;
% k_s_0_, ;
% k_s_1_, ;
% ;
% The (n_scatter)-by-(n_unrolled) scatter_from_tensor_s10__ matrix stores the ;
% interpolation weights linking a_k_c_10_ to the unrolled a_k_s_. ;
%%%%;

str_thisfunction = 'cylinder_k_c_scatter_from_tensor_adaptive_interpolate_n_9';

if (nargin<7);
rng(1);
%%%%%%%%;
flag_verbose = 1;
flag_disp = 1; nf=0;
disp(sprintf(' %% testing %s',str_thisfunction));
%%%%%%%%;
disp(sprintf(' %% first testing uniform n_k_c_1_'));
n_k_c_0 = 64+0;
n_k_c_1 = 64+1;
k_p_r_max = 48/(2*pi);
k_c_0_lim_ = 1.0*k_p_r_max*[-1,+1];
k_c_1_lim_ = [0,2*pi];
k_c_0_ = transpose(linspace(k_c_0_lim_(1),k_c_0_lim_(2),n_k_c_0));
k_c_1_ = linspace(k_c_1_lim_(1),k_c_1_lim_(2),n_k_c_1+1); k_c_1_ = transpose(k_c_1_(1:n_k_c_1));
[k_c_1_10__,k_c_0_10__] = ndgrid(k_c_1_,k_c_0_);
delta_ = [0.6;3.0];
f_c_10__ = exp(+i*(k_c_0_10__*delta_(1+0) + k_c_1_10__*delta_(1+1)));
dfd0_c_10__ = (+i*delta_(1+0))*exp(+i*(k_c_0_10__*delta_(1+0) + k_c_1_10__*delta_(1+1)));
dfd1_c_10__ = (+i*delta_(1+1))*exp(+i*(k_c_0_10__*delta_(1+0) + k_c_1_10__*delta_(1+1)));
ddfd00_c_10__ = (+i*delta_(1+0))*(+i*delta_(1+0))*exp(+i*(k_c_0_10__*delta_(1+0) + k_c_1_10__*delta_(1+1)));
ddfd01_c_10__ = (+i*delta_(1+0))*(+i*delta_(1+1))*exp(+i*(k_c_0_10__*delta_(1+0) + k_c_1_10__*delta_(1+1)));
ddfd11_c_10__ = (+i*delta_(1+1))*(+i*delta_(1+1))*exp(+i*(k_c_0_10__*delta_(1+0) + k_c_1_10__*delta_(1+1)));
n_scatter = 1024;
k_s_0_ = 0.5*k_p_r_max*(2*rand(n_scatter,1)-1);
k_s_1_ = 2*pi*rand(n_scatter,1);
tmp_ij_ = randperm(n_scatter,floor(n_scatter/4));
k_s_1_(tmp_ij_) = 2*pi*round(4*rand(numel(tmp_ij_),1))/4; %<-- ensure that some values of k_s_1_ are at multiples of pi/2. ;
n_order = 5;
[ ...
 scatter_from_tensor_s10__ ...
 dscatterd0_from_tensor_s10__ ...
 dscatterd1_from_tensor_s10__ ...
 ddscatterd00_from_tensor_s10__ ...
 ddscatterd01_from_tensor_s10__ ...
 ddscatterd11_from_tensor_s10__ ...
] = ...
cylinder_k_c_scatter_from_tensor_adaptive_interpolate_n_9( ...
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
f_s_inte_ = scatter_from_tensor_s10__*f_c_10__(:);
dfd0_s_inte_ = dscatterd0_from_tensor_s10__*f_c_10__(:);
dfd1_s_inte_ = dscatterd1_from_tensor_s10__*f_c_10__(:);
ddfd00_s_inte_ = ddscatterd00_from_tensor_s10__*f_c_10__(:);
ddfd01_s_inte_ = ddscatterd01_from_tensor_s10__*f_c_10__(:);
ddfd11_s_inte_ = ddscatterd11_from_tensor_s10__*f_c_10__(:);
%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 2; p_col = 3; np=0;
subplot(p_row,p_col,1+np);np=np+1;
hold on;
plot3(k_c_0_10__(:),k_c_1_10__(:),real(f_c_10__(:)),'co');
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
plot3(k_c_0_10__(:),k_c_1_10__(:),real(dfd0_c_10__(:)),'co');
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
plot3(k_c_0_10__(:),k_c_1_10__(:),real(dfd1_c_10__(:)),'co');
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
plot3(k_c_0_10__(:),k_c_1_10__(:),real(ddfd00_c_10__(:)),'co');
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
plot3(k_c_0_10__(:),k_c_1_10__(:),real(ddfd01_c_10__(:)),'co');
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
plot3(k_c_0_10__(:),k_c_1_10__(:),real(ddfd11_c_10__(:)),'co');
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
%%%%%%%%%%%%%%%%;
%%%%%%%%;
disp(sprintf(' %% now testing nonuniform n_k_c_1_'));
n_k_c_0 = 64+0;
n_k_c_1_ = 32+ceil(32*rand(n_k_c_0,1));
n_k_c_1_max = max(n_k_c_1_);
n_k_c_1_sum = sum(n_k_c_1_);
n_k_c_1_csum_ = cumsum([0;n_k_c_1_]);
k_p_r_max = 48/(2*pi);
k_c_0_lim_ = 1.0*k_p_r_max*[-1,+1];
k_c_1_lim__ = zeros(n_k_c_0,2);
for nk_c_0=0:n_k_c_0-1;
k_c_1_lim__(1+nk_c_0,:) = [0,2*pi] + pi/16*randn;
end;%for nk_c_0=0:n_k_c_0-1;
k_c_0_ = transpose(linspace(k_c_0_lim_(1),k_c_0_lim_(2),n_k_c_0));
k_c_1_10_ = zeros(n_k_c_1_sum,1);
k_c_0_10_ = zeros(n_k_c_1_sum,1);
for nk_c_0=0:n_k_c_0-1;
k_c_0 = k_c_0_(1+nk_c_0);
n_k_c_1 = n_k_c_1_(1+nk_c_0);
n_k_c_1_csum = n_k_c_1_csum_(1+nk_c_0);
k_c_1_lim_ = k_c_1_lim__(1+nk_c_0,:);
k_c_1_start = k_c_1_lim_(1+0); k_c_1_final = k_c_1_lim_(1+1);
k_c_1_ = linspace(k_c_1_start,k_c_1_final,n_k_c_1+1); k_c_1_ = transpose(k_c_1_(1:n_k_c_1));
k_c_1_10_(1+n_k_c_1_csum+[0:n_k_c_1-1]) = k_c_1_;
k_c_0_10_(1+n_k_c_1_csum+[0:n_k_c_1-1]) = k_c_0*ones(n_k_c_1,1);
end;%for nk_c_0=0:n_k_c_0-1;
delta_ = [0.6;3.0];
f_c_10_ = exp(+i*(k_c_0_10_*delta_(1+0) + k_c_1_10_*delta_(1+1)));
dfd0_c_10_ = (+i*delta_(1+0))*exp(+i*(k_c_0_10_*delta_(1+0) + k_c_1_10_*delta_(1+1)));
dfd1_c_10_ = (+i*delta_(1+1))*exp(+i*(k_c_0_10_*delta_(1+0) + k_c_1_10_*delta_(1+1)));
ddfd00_c_10_ = (+i*delta_(1+0))*(+i*delta_(1+0))*exp(+i*(k_c_0_10_*delta_(1+0) + k_c_1_10_*delta_(1+1)));
ddfd01_c_10_ = (+i*delta_(1+0))*(+i*delta_(1+1))*exp(+i*(k_c_0_10_*delta_(1+0) + k_c_1_10_*delta_(1+1)));
ddfd11_c_10_ = (+i*delta_(1+1))*(+i*delta_(1+1))*exp(+i*(k_c_0_10_*delta_(1+0) + k_c_1_10_*delta_(1+1)));
n_scatter = 1024;
k_s_0_ = 0.5*k_p_r_max*(2*rand(n_scatter,1)-1);
k_s_1_ = 2*pi*rand(n_scatter,1);
tmp_ij_ = randperm(n_scatter,floor(n_scatter/4));
k_s_1_(tmp_ij_) = 2*pi*round(4*rand(numel(tmp_ij_),1))/4; %<-- ensure that some values of k_s_1_ are at multiples of pi/2. ;
n_order = 5;
[ ...
 scatter_from_tensor_s10__ ...
 dscatterd0_from_tensor_s10__ ...
 dscatterd1_from_tensor_s10__ ...
 ddscatterd00_from_tensor_s10__ ...
 ddscatterd01_from_tensor_s10__ ...
 ddscatterd11_from_tensor_s10__ ...
] = ...
cylinder_k_c_scatter_from_tensor_adaptive_interpolate_n_9( ...
 n_order ...
,n_k_c_0 ...
,n_k_c_1_ ...
,k_c_0_lim_ ...
,k_c_1_lim__ ...
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
f_s_inte_ = scatter_from_tensor_s10__*f_c_10_(:);
dfd0_s_inte_ = dscatterd0_from_tensor_s10__*f_c_10_(:);
dfd1_s_inte_ = dscatterd1_from_tensor_s10__*f_c_10_(:);
ddfd00_s_inte_ = ddscatterd00_from_tensor_s10__*f_c_10_(:);
ddfd01_s_inte_ = ddscatterd01_from_tensor_s10__*f_c_10_(:);
ddfd11_s_inte_ = ddscatterd11_from_tensor_s10__*f_c_10_(:);
%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 2; p_col = 3; np=0;
subplot(p_row,p_col,1+np);np=np+1;
hold on;
plot3(k_c_0_10_(:),k_c_1_10_(:),real(f_c_10_(:)),'co');
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
plot3(k_c_0_10_(:),k_c_1_10_(:),real(dfd0_c_10_(:)),'co');
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
plot3(k_c_0_10_(:),k_c_1_10_(:),real(dfd1_c_10_(:)),'co');
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
plot3(k_c_0_10_(:),k_c_1_10_(:),real(ddfd00_c_10_(:)),'co');
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
plot3(k_c_0_10_(:),k_c_1_10_(:),real(ddfd01_c_10_(:)),'co');
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
plot3(k_c_0_10_(:),k_c_1_10_(:),real(ddfd11_c_10_(:)),'co');
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
if (nargin<1+na); n_k_c_1_=[]; end; na=na+1;
if (nargin<1+na); k_c_0_lim_=[]; end; na=na+1;
if (nargin<1+na); k_c_1_lim__=[]; end; na=na+1;
if (nargin<1+na); n_scatter=[]; end; na=na+1;
if (nargin<1+na); k_s_0_=[]; end; na=na+1;
if (nargin<1+na); k_s_1_=[]; end; na=na+1;

if numel(n_k_c_1_)==1;
n_k_c_1_ = repmat(n_k_c_1_(1+0),[n_k_c_0,1]);
end;%if numel(n_k_c_1_)==1;

n_k_c_1_max = max(n_k_c_1_);
n_k_c_1_sum = sum(n_k_c_1_);
n_k_c_1_csum_ = cumsum([0;n_k_c_1_]);
if (flag_verbose>0); disp(sprintf(' %% n_k_c_1_max %d n_k_c_1_sum %d',n_k_c_1_max,n_k_c_1_sum)); end;
if (flag_verbose>0); disp(sprintf(' %% n_k_c_1_: %s',num2str(transpose(n_k_c_1_),'%d '))); end;
if (flag_verbose>0); disp(sprintf(' %% n_k_c_1_csum_: %s',num2str(transpose(n_k_c_1_csum_),'%d '))); end;

if numel(k_c_1_lim__)==2;
k_c_1_lim__ = repmat([k_c_1_lim__(1+0),k_c_1_lim__(1+1)],[n_k_c_0,1]);
end;%if numel(k_c_1_lim__)==2;
if (flag_verbose>1); disp(sprintf(' %% k_c_1_lim__:')); end;
if (flag_verbose>1); disp(transpose(k_c_1_lim__)); end;

%%%%%%%%;
n_part_order = round((n_order-1)/2);
if (flag_verbose>0); disp(sprintf(' %% n_order %d n_part_order %d',n_order,n_part_order)); end;
%%%%;
if (n_k_c_0<=2*n_part_order+1); disp(sprintf(' %% Warning, n_k_c_0 %d, 2*n_part_order+1 %d in %s',n_k_c_0,2*n_part_order+1,str_thisfunction)); end;
%%%%;
dk_c_0 = diff(k_c_0_lim_)/max(1,n_k_c_0-1); %<-- interpreted as interval. ;
dk_c_1_ = zeros(n_k_c_0,1);
for nk_c_0=0:n_k_c_0-1;
n_k_c_1 = n_k_c_1_(1+nk_c_0);
dk_c_1_(1+nk_c_0) = diff(k_c_1_lim__(1+nk_c_0,:))/max(1,n_k_c_1-0); %<-- interpreted as periodic array. ;
end;%for nk_c_0=0:n_k_c_0-1;
if (flag_verbose>0); disp(sprintf(' %% dk_c_1_: %s',num2str(transpose(dk_c_1_),'%0.4f '))); end;
%%%%;
k_c_0_lim_use_ = [k_c_0_lim_(1)+n_part_order*dk_c_0,k_c_0_lim_(2)-n_part_order*dk_c_0];
%%%%;
if numel(efind(k_s_0_<=min(k_c_0_lim_use_)))>0; disp(sprintf(' %% Warning, k_s_0_ out of bounds in %s',str_thisfunction)); end;
if numel(efind(k_s_0_>=max(k_c_0_lim_use_)))>0; disp(sprintf(' %% Warning, k_s_0_ out of bounds in %s',str_thisfunction)); end;
k_s_0_ = max(min(k_c_0_lim_use_),min(max(k_c_0_lim_use_),k_s_0_));
%%%%;
k_c_0_start = k_c_0_lim_(1); k_c_0_final = k_c_0_lim_(2);
k_c_1_start_ = k_c_1_lim__(:,1+0); k_c_1_final_ = k_c_1_lim__(:,1+1);
%%%%;
if (flag_verbose>1);
disp(sprintf(' %% k_c_0_start %+0.6f k_c_0_final %+0.6f: n_k_c_0 %d --> dk_c_0 %+0.6f',k_c_0_start,k_c_0_final,n_k_c_0,dk_c_0));
for nk_c_0=0:n_k_c_0-1;
n_k_c_1 = k_c_1_start_(1+nk_c_0);
dk_c_1 = dk_c_1_(1+nk_c_0);
k_c_1_start = k_c_1_start_(1+nk_c_0);
k_c_1_final = k_c_1_final_(1+nk_c_0);
disp(sprintf(' %% k_c_1_start %+0.6f k_c_1_final %+0.6f: n_k_c_1 %d --> dk_c_1 %+0.6f',k_c_1_start,k_c_1_final,n_k_c_1,dk_c_1));
end;%for nk_c_0=0:n_k_c_0-1;
end;%if (flag_verbose>1);
%%%%%%%%;

%%%%%%%%;
parameter = struct('type','parameter');
parameter.flag_verbose = 0*flag_verbose;
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
,weight_denominator_0_x_ ...
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

nx_col_0_xs__ = repmat(transpose(0:n_nodex-1),[1,n_scatter]);
ns_col_0_xs__ = repmat(0:n_scatter-1,[n_nodex,1]);
[index_col_0_sort_xs_,tmp_ij_] = sort(index_col_0_xs__(:),'ascend');
n_xs_per_nk_c_0_ = full(sparse(1+index_col_0_sort_xs_,1+0,1,n_k_c_0,1));
n_xs_max = max(n_xs_per_nk_c_0_);
n_xs_sum = sum(n_xs_per_nk_c_0_);
assert(n_xs_sum==n_nodex*n_scatter);
n_xs_csum_ = cumsum([0;n_xs_per_nk_c_0_]);
if (flag_verbose>0); disp(sprintf(' %% n_xs_per_nk_c_0_: %s',num2str(transpose(n_xs_per_nk_c_0_),'%d '))); end;
if (flag_verbose>0); disp(sprintf(' %% n_xs_csum_: %s',num2str(transpose(n_xs_csum_),'%d '))); end;
nxs_from_nk_c_0__ = cell(n_k_c_0,1);
nx_from_nk_c_0__ = cell(n_k_c_0,1);
ns_from_nk_c_0__ = cell(n_k_c_0,1);
weight_denominator_0_from_nk_c_0__ = cell(n_k_c_0,1);
weight_numerator_0_from_nk_c_0__ = cell(n_k_c_0,1);
if flag_d; weight_dnumerator_0_from_nk_c_0__ = cell(n_k_c_0,1); end;
if flag_dd; weight_ddnumerator_0_from_nk_c_0__ = cell(n_k_c_0,1); end;
tab=0;
for nk_c_0=0:n_k_c_0-1;
n_xs_sub = n_xs_per_nk_c_0_(1+nk_c_0);
nxs_from_nk_c_0_ = tmp_ij_(1+tab+[0:n_xs_sub-1])-1;
nxs_from_nk_c_0__{1+nk_c_0} = nxs_from_nk_c_0_;
nx_from_nk_c_0__{1+nk_c_0} = nx_col_0_xs__(1+nxs_from_nk_c_0_);
ns_from_nk_c_0__{1+nk_c_0} = ns_col_0_xs__(1+nxs_from_nk_c_0_);
weight_denominator_0_from_nk_c_0__{1+nk_c_0} = weight_denominator_0_x_(1+nx_from_nk_c_0__{1+nk_c_0});
weight_numerator_0_from_nk_c_0__{1+nk_c_0} = weight_numerator_0_xs__(1+nxs_from_nk_c_0_);
if flag_d; weight_dnumerator_0_from_nk_c_0__{1+nk_c_0} = weight_dnumerator_0_xs__(1+nxs_from_nk_c_0_); end;
if flag_dd; weight_ddnumerator_0_from_nk_c_0__{1+nk_c_0} = weight_ddnumerator_0_xs__(1+nxs_from_nk_c_0_); end;
tab = tab+n_xs_sub;
end;%for nk_c_0=0:n_k_c_0-1;
assert(tab==n_nodex*n_scatter);
%%%%%%%%;
index_col_xxs_ = zeros(n_nodex^2*n_scatter,1);
index_row_xxs_ = zeros(n_nodex^2*n_scatter,1);
weight_0_xxs_ = zeros(n_nodex^2*n_scatter,1);
weight_1_xxs_ = zeros(n_nodex^2*n_scatter,1);
if flag_d;
dweight_0_xxs_ = zeros(n_nodex^2*n_scatter,1);
dweight_1_xxs_ = zeros(n_nodex^2*n_scatter,1);
end;%if flag_d;
if flag_dd;
ddweight_0_xxs_ = zeros(n_nodex^2*n_scatter,1);
ddweight_1_xxs_ = zeros(n_nodex^2*n_scatter,1);
end;%if flag_dd;
for nk_c_0=0:n_k_c_0-1;
n_xs_sub = n_xs_per_nk_c_0_(1+nk_c_0);
n_xs_csum = n_xs_csum_(1+nk_c_0);
n_k_c_1 = n_k_c_1_(1+nk_c_0);
n_k_c_1_csum = n_k_c_1_csum_(1+nk_c_0);
dk_c_1 = dk_c_1_(1+nk_c_0);
k_c_1_start = k_c_1_start_(1+nk_c_0);
k_c_1_final = k_c_1_final_(1+nk_c_0);
nxs_from_nk_c_0_ = nxs_from_nk_c_0__{1+nk_c_0};
nx_from_nk_c_0_ = nx_from_nk_c_0__{1+nk_c_0};
ns_from_nk_c_0_ = ns_from_nk_c_0__{1+nk_c_0};
k_xs_sub_1_ = k_s_1_(1+ns_from_nk_c_0_);
if (flag_verbose>1);
disp(sprintf(' %% nk_c_0 %d/%d n_xs_sub %d n_xs_csum %d n_k_c_1 %d n_k_c_1_csum %d',nk_c_0,n_k_c_0,n_xs_sub,n_xs_csum,n_k_c_1,n_k_c_1_csum));
disp(sprintf(' %% dk_c_1 %0.6f k_c_1_start %0.6f k_c_1_final %0.6f ',dk_c_1,k_c_1_start,k_c_1_final));
end;%if (flag_verbose>1);
if n_xs_sub> 0;
parameter = struct('type','parameter');
parameter.flag_verbose = 0*flag_verbose;
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
,index_col_1_xxs_sub__ ...
,index_row_xxs_sub__ ...
,weight_denominator_1_x_ ...
,weight_numerator_1_xxs_sub__ ...
,weight_dnumerator_1_xxs_sub__ ...
,weight_ddnumerator_1_xxs_sub__ ...
] = ...
scatter_from_tensor_helper_n_9( ...
 parameter ...
,n_xs_sub ...
,k_xs_sub_1_ ...
,k_c_1_start ...
,dk_c_1 ...
);
%%%%%%%%;
tmp_index_out_ = n_nodex*n_xs_csum+[0:n_nodex*n_xs_sub-1];
index_col_1_xxs_sub__ = periodize(index_col_1_xxs_sub__,0,n_k_c_1);
index_col_xxs_sub_ = reshape(index_col_1_xxs_sub__,[n_nodex*n_xs_sub,1]) + n_k_c_1_csum;
index_row_xxs_sub_ = reshape(repmat(reshape(ns_from_nk_c_0_,[1,n_xs_sub]),[n_nodex,1]),[n_nodex*n_xs_sub,1]);
index_col_xxs_(1+tmp_index_out_) = index_col_xxs_sub_;
index_row_xxs_(1+tmp_index_out_) = index_row_xxs_sub_;
weight_denominator_0_from_nk_c_0_ = weight_denominator_0_from_nk_c_0__{1+nk_c_0};
weight_numerator_0_from_nk_c_0_ = weight_numerator_0_from_nk_c_0__{1+nk_c_0};
if flag_d; weight_dnumerator_0_from_nk_c_0_ = weight_dnumerator_0_from_nk_c_0__{1+nk_c_0}; end;
if flag_dd; weight_ddnumerator_0_from_nk_c_0_ = weight_ddnumerator_0_from_nk_c_0__{1+nk_c_0}; end;
weight_0_xxs_sub_ = reshape(repmat(reshape(weight_numerator_0_from_nk_c_0_./weight_denominator_0_from_nk_c_0_,[1,n_xs_sub]),[n_nodex,1]),[n_nodex*n_xs_sub,1]);
if flag_d; dweight_0_xxs_sub_ = reshape(repmat(reshape(weight_dnumerator_0_from_nk_c_0_./weight_denominator_0_from_nk_c_0_,[1,n_xs_sub]),[n_nodex,1]),[n_nodex*n_xs_sub,1]) / dk_c_0 ; end; %<-- unprotected (signed) divide, assume dk_c_0 is not close to 0. ;
if flag_dd; ddweight_0_xxs_sub_ = reshape(repmat(reshape(weight_ddnumerator_0_from_nk_c_0_./weight_denominator_0_from_nk_c_0_,[1,n_xs_sub]),[n_nodex,1]),[n_nodex*n_xs_sub,1]) / dk_c_0^2 ; end; %<-- unprotected (signed) divide, assume dk_c_0 is not close to 0. ;
weight_0_xxs_(1+tmp_index_out_) = weight_0_xxs_sub_;
if flag_d; dweight_0_xxs_(1+tmp_index_out_) = dweight_0_xxs_sub_; end;
if flag_dd; ddweight_0_xxs_(1+tmp_index_out_) = ddweight_0_xxs_sub_; end;
weight_1_xxs_sub__ = bsxfun(@rdivide,weight_numerator_1_xxs_sub__,reshape(weight_denominator_1_x_,[n_nodex,1]));
if flag_d; dweight_1_xxs_sub__ = bsxfun(@rdivide,weight_dnumerator_1_xxs_sub__,reshape(weight_denominator_1_x_,[n_nodex,1])) / dk_c_1 ; end; %<-- unprotected (signed) divide, assume dk_c_1 is not close to 0. ;
if flag_dd; ddweight_1_xxs_sub__ = bsxfun(@rdivide,weight_ddnumerator_1_xxs_sub__,reshape(weight_denominator_1_x_,[n_nodex,1])) / dk_c_1^2 ; end; %<-- unprotected (signed) divide, assume dk_c_1 is not close to 0. ;
weight_1_xxs_(1+tmp_index_out_) = weight_1_xxs_sub__(:);
if flag_d; dweight_1_xxs_(1+tmp_index_out_) = dweight_1_xxs_sub__(:); end;
if flag_dd; ddweight_1_xxs_(1+tmp_index_out_) = ddweight_1_xxs_sub__(:); end;
if flag_verbose>2;
disp(sprintf('tmp_index_out_:')); disp(transpose(tmp_index_out_));
disp(sprintf('index_col_1_xxs_sub__:')); disp(index_col_1_xxs_sub__);
disp(sprintf('index_col_xxs_sub_:')); disp(index_col_xxs_sub_);
disp(sprintf('index_row_xxs_sub_:')); disp(index_row_xxs_sub_);
disp(sprintf('weight_denominator_0_from_nk_c_0_:')); disp(transpose(weight_denominator_0_from_nk_c_0_));
disp(sprintf('weight_numerator_0_from_nk_c_0_:')); disp(transpose(weight_numerator_0_from_nk_c_0_));
disp(sprintf('weight_0_xxs_sub_:')); disp(transpose(weight_0_xxs_sub_));
disp(sprintf('weight_1_xxs_sub__:')); disp(weight_1_xxs_sub__);
end;%if flag_verbose>2;
end;%if n_xs_sub> 0;
end;%for nk_c_0=0:n_k_c_0-1;

if flag_verbose>2;
disp(sprintf(' %% n_scatter %d n_k_c_1_sum %d',n_scatter,n_k_c_1_sum));
disp(sprintf('index_row_xxs_:')); disp(index_row_xxs_);
disp(sprintf('index_col_xxs_:')); disp(index_col_xxs_);
disp(sprintf('weight_0_xxs_:')); disp(weight_0_xxs_);
disp(sprintf('weight_1_xxs_:')); disp(weight_1_xxs_);
end;%if flag_verbose>2;

if sum(~isfinite(weight_0_xxs_))> 0;
disp(sprintf(' %% Warning: sum(~isfinite(weight_0_xxs_)) %d in %s',sum(~isfinite(weight_0_xxs_)),str_thisfunction));
end;%if sum(~isfinite(weight_0_xxs_))> 0;
if sum(~isfinite(weight_1_xxs_))> 0;
disp(sprintf(' %% Warning: sum(~isfinite(weight_1_xxs_)) %d in %s',sum(~isfinite(weight_1_xxs_)),str_thisfunction));
end;%if sum(~isfinite(weight_1_xxs_))> 0;
if flag_d;
if sum(~isfinite(dweight_0_xxs_))> 0;
disp(sprintf(' %% Warning: sum(~isfinite(dweight_0_xxs_)) %d in %s',sum(~isfinite(dweight_0_xxs_)),str_thisfunction));
end;%if sum(~isfinite(dweight_0_xxs_))> 0;
if sum(~isfinite(dweight_1_xxs_))> 0;
disp(sprintf(' %% Warning: sum(~isfinite(dweight_1_xxs_)) %d in %s',sum(~isfinite(dweight_1_xxs_)),str_thisfunction));
end;%if sum(~isfinite(dweight_1_xxs_))> 0;
end;%if flag_d;
if flag_dd;
if sum(~isfinite(ddweight_0_xxs_))> 0;
disp(sprintf(' %% Warning: sum(~isfinite(ddweight_0_xxs_)) %d in %s',sum(~isfinite(ddweight_0_xxs_)),str_thisfunction));
end;%if sum(~isfinite(ddweight_0_xxs_))> 0;
if sum(~isfinite(ddweight_1_xxs_))> 0;
disp(sprintf(' %% Warning: sum(~isfinite(ddweight_1_xxs_)) %d in %s',sum(~isfinite(ddweight_1_xxs_)),str_thisfunction));
end;%if sum(~isfinite(ddweight_1_xxs_))> 0;
end;%if flag_dd;

scatter_from_tensor_s10__ = sparse(1+index_row_xxs_,1+index_col_xxs_,weight_0_xxs_.*weight_1_xxs_,n_scatter,n_k_c_1_sum);
if flag_d;
dscatterd0_from_tensor_s10__ = sparse(1+index_row_xxs_,1+index_col_xxs_,dweight_0_xxs_.*weight_1_xxs_,n_scatter,n_k_c_1_sum);
dscatterd1_from_tensor_s10__ = sparse(1+index_row_xxs_,1+index_col_xxs_,weight_0_xxs_.*dweight_1_xxs_,n_scatter,n_k_c_1_sum);
end;%if flag_d;
if flag_dd;
ddscatterd00_from_tensor_s10__ = sparse(1+index_row_xxs_,1+index_col_xxs_,ddweight_0_xxs_.*weight_1_xxs_,n_scatter,n_k_c_1_sum);
ddscatterd01_from_tensor_s10__ = sparse(1+index_row_xxs_,1+index_col_xxs_,dweight_0_xxs_.*dweight_1_xxs_,n_scatter,n_k_c_1_sum);
ddscatterd11_from_tensor_s10__ = sparse(1+index_row_xxs_,1+index_col_xxs_,weight_0_xxs_.*ddweight_1_xxs_,n_scatter,n_k_c_1_sum);
end;%if flag_dd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
if (flag_verbose>0); disp(sprintf(' %% ')); end;
if (flag_verbose>0); disp(sprintf(' %% returning')); end; return;


