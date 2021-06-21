function ...
[ ...
 parameter ...
,M_k_p_out_ ...
,M_x_c_out_ ...
,delta_x_c_0 ...
,delta_x_c_1 ...
] = ...
image_center_0( ...
 parameter ...
,n_x ...
,diameter_x_c ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,M_k_p_0in_ ...
,weight_2d_k_all_ ...
);

verbose=0;
if (verbose); disp(sprintf(' %% [entering image_center_0]')); end;

if isempty(parameter);
parameter = struct('type','parameter');
end;%if isempty(parameter);
%%%%%%%%;
if (~isfield(parameter,'tolerance_master')); parameter.tolerance_master = 1e-2; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'p_cut')); parameter.p_cut = 95; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'n_iteration')); parameter.n_iteration = 32; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'factor_mask')); parameter.factor_mask = sqrt(2); end; %<-- parameter_bookmark. ;
%%%%%%%%;
tolerance_master = parameter.tolerance_master;
p_cut = parameter.p_cut;
n_iteration = parameter.n_iteration;
factor_mask = parameter.factor_mask;

n_w_ = n_w_(:);
n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);
half_diameter_x_c = diameter_x_c/2;
x_c_0_ = -half_diameter_x_c + transpose([0:n_x-1]/n_x)*diameter_x_c;
x_c_1_ = -half_diameter_x_c + transpose([0:n_x-1]/n_x)*diameter_x_c;
[x_c_0__,x_c_1__] = ndgrid(x_c_0_,x_c_1_);
x_p_r__ = sqrt(x_c_0__.^2 + x_c_1__.^2);
x_c_mask__ = x_p_r__<=half_diameter_x_c/factor_mask;
edge_M_x_c_ = x_p_r__> half_diameter_x_c/factor_mask;
cent_M_x_c_ = x_p_r__<=half_diameter_x_c/factor_mask;
index_edge_M_x_c_ = efind(edge_M_x_c_);
index_cent_M_x_c_ = efind(cent_M_x_c_);

M_k_p_ = M_k_p_0in_;
delta_x_c_0 = 0;
delta_x_c_1 = 0;
flag_continue=1; niteration=0;
while flag_continue;
M_x_c_ = interp_k_p_to_x_c_nufft(n_x,diameter_x_c,n_x,diameter_x_c,n_k_p_r,k_p_r_,n_w_,M_k_p_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x^2) * n_w_sum;
M_cut = prctile(real(M_x_c_(1+index_edge_M_x_c_)),p_cut);
M_mask_rect_x_c_ = max(0,real(M_x_c_.*x_c_mask__)-M_cut); %<-- radial mask. ;
M_mask_rect_avg = mean(M_mask_rect_x_c_,'all');
M_mask_rect_x_c_0_avg = mean(M_mask_rect_x_c_/M_mask_rect_avg.*x_c_0__,'all');
M_mask_rect_x_c_1_avg = mean(M_mask_rect_x_c_/M_mask_rect_avg.*x_c_1__,'all');
M_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p_,-1*M_mask_rect_x_c_0_avg,-1*M_mask_rect_x_c_1_avg);
delta_x_c_0 = delta_x_c_0 + M_mask_rect_x_c_0_avg;
delta_x_c_1 = delta_x_c_1 + M_mask_rect_x_c_1_avg;
if (verbose); disp(sprintf(' %% niteration %d/%d, delta_ = [%0.6f,%0.6f], delta_upd_ = [%0.6f,%0.6f]',niteration,n_iteration,delta_x_c_0,delta_x_c_1,M_mask_rect_x_c_0_avg,M_mask_rect_x_c_1_avg)); end;
flag_continue = (fnorm([M_mask_rect_x_c_0_avg;M_mask_rect_x_c_1_avg])/half_diameter_x_c>=tolerance_master) & (niteration<n_iteration);
niteration=niteration+1;
end;%while;

M_x_c_ = interp_k_p_to_x_c_nufft(n_x,diameter_x_c,n_x,diameter_x_c,n_k_p_r,k_p_r_,n_w_,M_k_p_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x^2) * n_w_sum;
M_k_p_out_ = M_k_p_;
M_x_c_out_ = M_x_c_;

if (verbose); disp(sprintf(' %% [finished image_center_0]')); end;
