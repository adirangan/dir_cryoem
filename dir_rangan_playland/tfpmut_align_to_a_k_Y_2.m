function ...
[ ...
 parameter ...
] = ...
tfpmut_align_to_a_k_Y_2( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_ ...
,weight_2d_k_p_wk_ ...
,l_max_ ...
,n_CTF ...
,CTF_k_p_wkC__ ...
,index_nCTF_from_nM_ ...
,n_M ...
,M_k_p_wkM__ ...
,a_k_Y_yk_ ...
,n_iteration ...
,euler_polar_a_calc_Mi__ ...
,euler_azimu_b_calc_Mi__ ...
,euler_gamma_z_calc_Mi__ ...
,image_delta_x_calc_Mi__ ...
,image_delta_y_calc_Mi__ ...
,euler_polar_a_true_M_ ...
,euler_azimu_b_true_M_ ...
,euler_gamma_z_true_M_ ...
,image_delta_x_true_M_ ...
,image_delta_y_true_M_ ...
);

str_thisfunction = 'tfpmut_align_to_a_k_Y_2';

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); k_p_r_=[]; end; na=na+1;
if (nargin<1+na); k_p_r_max=[]; end; na=na+1;
if (nargin<1+na); weight_3d_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); weight_2d_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); n_w_=[]; end; na=na+1;
if (nargin<1+na); weight_2d_k_p_wk_=[]; end; na=na+1;
if (nargin<1+na); l_max_=[]; end; na=na+1;
if (nargin<1+na); n_CTF=[]; end; na=na+1;
if (nargin<1+na); CTF_k_p_wkC__=[]; end; na=na+1;
if (nargin<1+na); index_nCTF_from_nM_=[]; end; na=na+1;
if (nargin<1+na); n_M=[]; end; na=na+1;
if (nargin<1+na); M_k_p_wkM__=[]; end; na=na+1;
if (nargin<1+na); a_k_Y_yk_=[]; end; na=na+1;
if (nargin<1+na); n_iteration=[]; end; na=na+1;
if (nargin<1+na); euler_polar_a_calc_Mi__=[]; end; na=na+1;
if (nargin<1+na); euler_azimu_b_calc_Mi__=[]; end; na=na+1;
if (nargin<1+na); euler_gamma_z_calc_Mi__=[]; end; na=na+1;
if (nargin<1+na); image_delta_x_calc_Mi__=[]; end; na=na+1;
if (nargin<1+na); image_delta_y_calc_Mi__=[]; end; na=na+1;
if (nargin<1+na); euler_polar_a_true_M_=[]; end; na=na+1;
if (nargin<1+na); euler_azimu_b_true_M_=[]; end; na=na+1;
if (nargin<1+na); euler_gamma_z_true_M_=[]; end; na=na+1;
if (nargin<1+na); image_delta_x_true_M_=[]; end; na=na+1;
if (nargin<1+na); image_delta_y_true_M_=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;%if isempty(parameter);
%%%%%%%%;
if (~isfield(parameter,'flag_verbose')); parameter.flag_verbose = 0; end; %<-- parameter_bookmark. ;
flag_verbose = parameter.flag_verbose;
if (~isfield(parameter,'tolerance_master')); parameter.tolerance_master = 1e-2; end; %<-- parameter_bookmark. ;
tolerance_master = parameter.tolerance_master;
if (~isfield(parameter,'flag_qbp_vs_lsq')); parameter.flag_qbp_vs_lsq = 1; end; %<-- parameter_bookmark. ;
flag_qbp_vs_lsq = parameter.flag_qbp_vs_lsq;
if (~isfield(parameter,'cg_lsq_n_order')); parameter.cg_lsq_n_order = 5; end; %<-- parameter_bookmark. ;
cg_lsq_n_order = parameter.cg_lsq_n_order;
if (~isfield(parameter,'qbp_eps')); parameter.qbp_eps = parameter.tolerance_master; end; %<-- parameter_bookmark. ;
qbp_eps = parameter.qbp_eps;

if (flag_verbose>0); disp(sprintf(' %% [entering tfpmut_align_to_a_k_Y_2]')); end;

if ( isempty(n_iteration)); n_iteration=0; end;
if (n_iteration<=0); n_iteration = size(euler_polar_a_calc_Mi__,2); end;

n_w_ = n_w_(1:n_k_p_r); n_w_sum = sum(n_w_); n_w_max = max(n_w_); n_w_csum_ = cumsum([0;n_w_]);
if (std(diff(n_w_csum_),1)>1e-6); disp(sprintf(' %% Warning, flag_uniform_over_n_k_p_r==1 expected in %s',str_thisfunction)); end;
n_w_max = n_w_csum_(1+1);
assert(n_w_sum==n_w_max*n_k_p_r);
if mod(n_w_max,2)~=0; disp(sprintf(' %% Warning, n_w_max %d in %s',n_w_max,str_thisfunction)); end;

%%%%%%%%;
% construct CTF of same size as images. ;
%%%%%%%%
if (size(CTF_k_p_wkC__,1)==n_k_p_r);
n_CTF = size(CTF_k_p_wkC__,2);
CTF_k_p_r_kC__ = CTF_k_p_wkC__;
CTF_k_p_wkC__ = reshape(repmat(reshape(CTF_k_p_r_kC__,[1,n_k_p_r,n_CTF]),[n_w_max,1,1]),[n_w_sum,n_CTF]);
end;%if (size(CTF_k_p_wkC__,1)==n_k_p_r);
%%%%%%%%;
CTF_k_p_r_kC__ = reshape(mean(reshape(CTF_k_p_wkC__,[n_w_max,n_k_p_r,n_CTF]),1+0),[n_k_p_r,n_CTF]);
CTF_k_p_wkM__ = CTF_k_p_wkC__(:,1+index_nCTF_from_nM_);
%%%%%%%%;

flag_found=0;
if (~isempty(parameter.fname_align_a_k_Y_pre));
tmp_fname_mat = sprintf('%s.mat',parameter.fname_align_a_k_Y_pre);
if  exist(tmp_fname_mat,'file');
disp(sprintf(' %% %s found, loading',tmp_fname_mat));
load(tmp_fname_mat ...
     ,'X_best_i_','flag_flip_i_' ...
     ,'polar_a_best_i_' ...
     ,'azimu_b_best_i_' ...
     ,'gamma_z_best_i_' ...
     ,'delta_best_3i__' ...
    );
flag_found=1;
end;%if  exist(tmp_fname_mat,'file');
end;%if (~isempty(parameter.fname_align_a_k_Y_pre));

if ~flag_found;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
X_best_i_ = zeros(n_iteration,1);
flag_flip_i_ = zeros(n_iteration,1);
polar_a_best_i_ = zeros(n_iteration,1);
azimu_b_best_i_ = zeros(n_iteration,1);
gamma_z_best_i_ = zeros(n_iteration,1);
delta_best_3i__ = zeros(3,n_iteration);
%%%%%%%%;
for niteration=0:n_iteration-1;
%%%%%%%%;
tmp_euler_polar_a_M_ = euler_polar_a_calc_Mi__(:,1+niteration);
tmp_euler_azimu_b_M_ = euler_azimu_b_calc_Mi__(:,1+niteration);
tmp_euler_gamma_z_M_ = euler_gamma_z_calc_Mi__(:,1+niteration);
tmp_image_delta_x_M_ = image_delta_x_calc_Mi__(:,1+niteration);
tmp_image_delta_y_M_ = image_delta_y_calc_Mi__(:,1+niteration);
%%%%%%%%;
if flag_qbp_vs_lsq==0;
disp(sprintf(' %% Warning, flag_qbp_vs_lsq==%d not implemented in %s',flag_qbp_vs_lsq,str_thisfunction));
end;%if flag_qbp_vs_lsq==0;
%%%%%%%%;
if flag_qbp_vs_lsq==1;
tmp_t = tic;
c_k_Y_yk_ = ...
qbp_uniform_over_n_k_p_r_10( ...
 qbp_eps ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,n_M ...
,M_k_p_wkM__ ...
,index_nCTF_from_nM_ ...
,CTF_k_p_wkC__ ...
,tmp_euler_polar_a_M_ ...
,tmp_euler_azimu_b_M_ ...
,tmp_euler_gamma_z_M_ ...
,tmp_image_delta_x_M_ ...
,tmp_image_delta_y_M_ ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% niteration %.2d/%.2d: M_k_p_wkM__ --> c_k_Y_yk_ time %0.2fs',niteration,n_iteration,tmp_t)); end;
end;%if flag_qbp_vs_lsq==1;
%%%%%%%%;
tmp_t = tic;
if ~exist('d_mmlb____','var'); d_mmlb____=[]; end;
if ~exist('d_mmzb____','var'); d_mmzb____=[]; end;
if ~exist('n_UX_rank','var'); n_UX_rank=[]; end;
if ~exist('U_lz__','var'); U_lz__=[]; end;
N_wavelength = 0.0;
parameter.flag_verbose = max(0,flag_verbose-1);
[ ...
 parameter ...
,d_mmlb____ ....
,d_mmzb____ ....
,X_best ...
,flag_flip ...
,polar_a_best ...
,azimu_b_best ...
,gamma_z_best ...
,delta_best_ ...
,c_k_Y_best_yk_ ...
,~ ...
,~ ...
,~ ...
] = ...
register_spharm_to_spharm_3_stripped_wrap_0( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,N_wavelength ...
,l_max_ ...
,n_UX_rank ...
,U_lz__ ...
,d_mmlb____ ...
,d_mmzb____ ...
,a_k_Y_yk_ ...
,c_k_Y_yk_ ...
);
parameter.flag_verbose = flag_verbose;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% niteration %.2d/%.2d: X_best time %0.2fs',niteration,n_iteration,tmp_t)); end;
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% niteration %.2d/%.2d: X_best %0.4f flag_flip %d',niteration,n_iteration,X_best,flag_flip)); end;
X_best_i_(1+niteration) = X_best;
flag_flip_i_(1+niteration) = flag_flip;
polar_a_best_i_(1+niteration) = polar_a_best;
azimu_b_best_i_(1+niteration) = azimu_b_best;
gamma_z_best_i_(1+niteration) = gamma_z_best;
delta_best_3i__(:,1+niteration) = delta_best_;
%%%%%%%%;
end;%for niteration=0:n_iteration-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if ~flag_found;

flag_found=0;
if (~isempty(parameter.fname_align_a_k_Y_pre));
tmp_fname_mat = sprintf('%s.mat',parameter.fname_align_a_k_Y_pre);
if (~exist(tmp_fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',tmp_fname_mat));
save(tmp_fname_mat ...
     ,'X_best_i_','flag_flip_i_' ...
     ,'polar_a_best_i_' ...
     ,'azimu_b_best_i_' ...
     ,'gamma_z_best_i_' ...
     ,'delta_best_3i__' ...
    );
end;%if (~exist(tmp_fname_mat,'file'));
flag_found=1;
end;%if (~isempty(parameter.fname_align_a_k_Y_pre));

if flag_found;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (~isempty(parameter.fname_align_a_k_Y_pre));
tmp_fname_fig_jpg = sprintf('%s.jpg',parameter.fname_align_a_k_Y_pre);
if (~exist(tmp_fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',tmp_fname_fig_jpg));
delta_sigma = 1.0 * std([image_delta_x_true_M_(1:n_M);image_delta_y_true_M_(1:n_M)],1); %<-- no reduction. ;
dscale = 2.5;
c2d_euler__ = colormap_polar_a_azimu_b_2d(+euler_polar_a_true_M_(1:n_M),+euler_azimu_b_true_M_(1:n_M),0.35);
markersize_euler = 25;
c2d_delta__ = colormap_gaussian_2d(+image_delta_x_true_M_(1:n_M),+image_delta_y_true_M_(1:n_M),dscale*delta_sigma,0.35);
markersize_delta = 25;
figure(3); clf; figbig; 
p_row = 4; p_col = 3*ceil((1+n_iteration)/p_row); np=0;
%%%%%%%%;
subplot(p_row,p_col,1+np+[0,1]); np=np+2;
hold on;
patch(2*pi*[+0;+1;+1;+0],1*pi*[+0;+0;+1;+1],0.65*[1,1,1],'EdgeColor','none');
scatter(0*pi+euler_azimu_b_true_M_(1:n_M),1*pi-euler_polar_a_true_M_(1:n_M),markersize_euler,c2d_euler__(1:n_M,:),'filled');
hold off;
axisnotick;
axis([0,2*pi,0,1*pi]); axis equal; xlabel('azimu_b'); ylabel('polar_a'); title('true view');
%%%%%%%%;
subplot(p_row,p_col,1+np); np=np+1;
hold on;
patch(dscale*delta_sigma*[-1;+1;+1;-1],dscale*delta_sigma*[-1;-1;+1;+1],'k');
scatter(+image_delta_x_true_M_(1:n_M),+image_delta_y_true_M_(1:n_M),markersize_delta,c2d_delta__(1:n_M,:),'filled');
hold off;
axisnotick;
axis(dscale*delta_sigma*[-1,+1,-1,+1]); axis square; xlabel('dx'); ylabel('dy'); title('true delta');
%%%%%%%%;
for niteration=0:n_iteration-1;
%%%%%%%%;
X_best = X_best_i_(1+niteration);
flag_flip = flag_flip_i_(1+niteration);
polar_a_best = polar_a_best_i_(1+niteration);
azimu_b_best = azimu_b_best_i_(1+niteration);
gamma_z_best = gamma_z_best_i_(1+niteration);
delta_best_ = delta_best_3i__(:,1+niteration);
euler_best_ = [+azimu_b_best,+polar_a_best,+gamma_z_best];
R_best_ = euler_to_R(euler_best_);
%%%%%%%%;
tmp_euler_azimu_b_ = euler_azimu_b_calc_Mi__(:,1+niteration);
tmp_euler_polar_a_ = euler_polar_a_calc_Mi__(:,1+niteration);
tmp_k_c_0_ = sin(tmp_euler_polar_a_).*cos(tmp_euler_azimu_b_);
tmp_k_c_1_ = sin(tmp_euler_polar_a_).*sin(tmp_euler_azimu_b_);
tmp_k_c_2_ = cos(tmp_euler_polar_a_);
tmp__ = inv(R_best_) * transpose([tmp_k_c_0_(:) , tmp_k_c_1_(:) , tmp_k_c_2_(:)]);
tmp_k_c_0_(:) = transpose(tmp__(1+0,:));
tmp_k_c_1_(:) = transpose(tmp__(1+1,:));
tmp_k_c_2_(:) = transpose(tmp__(1+2,:));
tmp_k_c_01_ = sqrt(tmp_k_c_0_.^2 + tmp_k_c_1_.^2);
tmp_euler_azimu_b_ = periodize(atan2(tmp_k_c_1_,tmp_k_c_0_),0,2*pi);
tmp_euler_polar_a_ = atan2(tmp_k_c_01_,tmp_k_c_2_);
%%%%%%%%;
subplot(p_row,p_col,1+np+[0,1]); np=np+2; cla;
hold on;
patch(2*pi*[+0;+1;+1;+0],1*pi*[+0;+0;+1;+1],0.65*[1,1,1],'EdgeColor','none');
scatter(0*pi+tmp_euler_azimu_b_,1*pi-tmp_euler_polar_a_,markersize_euler,c2d_euler__(1:n_M,:),'filled');
hold off;
axisnotick;
axis([0,2*pi,0,1*pi]); axis equal; xlabel('azimu_b'); ylabel('polar_a');
title(sprintf('%d (flip %d) X %0.2f',niteration,flag_flip,X_best));
%%%%%%%%;
subplot(p_row,p_col,1+np); np=np+1;
hold on;
patch(dscale*delta_sigma*[-1;+1;+1;-1],dscale*delta_sigma*[-1;-1;+1;+1],'k');
scatter(+image_delta_x_calc_Mi__(1:n_M,1+niteration),+image_delta_y_calc_Mi__(1:n_M,1+niteration),markersize_delta,c2d_delta__(1:n_M,:),'filled');
hold off;
axisnotick;
axis(dscale*delta_sigma*[-1,+1,-1,+1]); axis square; xlabel('dx'); ylabel('dy'); 
title(sprintf('%d (delta)',niteration));
%%%%%%%%;
end;%for niteration=0:n_iteration-1;
%%%%%%%%;
sgtitle(sprintf('%s',tmp_fname_fig_jpg),'Interpreter','none');
print('-djpeg',tmp_fname_fig_jpg);
close(gcf);
end;%if (~exist(tmp_fname_fig_jpg,'file'));
end;%if (~isempty(parameter.fname_align_a_k_Y_pre));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_found;

if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
