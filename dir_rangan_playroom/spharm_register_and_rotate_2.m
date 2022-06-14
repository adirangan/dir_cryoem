function ...
[ ... 
 c_k_Y_reco_ ...
,X_best ...
,flag_flip ...
,euler_polar_a_best ...
,euler_azimu_b_best ...
,euler_gamma_z_best ...
,image_delta_x_best ...
,image_delta_y_best ...
,image_delta_z_best ...
] = ...
spharm_register_and_rotate_2( ...
 n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_true_ ...
,b_k_Y_reco_ ...
,N_wavelength ...
,X_best ...
,flag_flip ...
,euler_polar_a_best ...
,euler_azimu_b_best ...
,euler_gamma_z_best ...
,image_delta_x_best ...
,image_delta_y_best ...
,image_delta_z_best ...
);

verbose=0;
if (verbose); disp(sprintf(' %% [entering spharm_register_and_rotate_2]')); end;

na=0;
if (nargin<1+na); n_k_p_r = []; end; na=na+1;
if (nargin<1+na); k_p_r_ = []; end; na=na+1;
if (nargin<1+na); k_p_r_max = []; end; na=na+1;
if (nargin<1+na); weight_3d_k_p_r_ = []; end; na=na+1;
if (nargin<1+na); l_max_ = []; end; na=na+1;
if (nargin<1+na); a_k_Y_true_ = []; end; na=na+1;
if (nargin<1+na); a_k_Y_reco_ = []; end; na=na+1;
if (nargin<1+na); N_wavelength = []; end; na=na+1;
if (nargin<1+na); X_best = []; end; na=na+1;
if (nargin<1+na); flag_flip = []; end; na=na+1;
if (nargin<1+na); euler_polar_a_best = []; end; na=na+1;
if (nargin<1+na); euler_azimu_b_best = []; end; na=na+1;
if (nargin<1+na); euler_gamma_z_best = []; end; na=na+1;
if (nargin<1+na); image_delta_x_best = []; end; na=na+1;
if (nargin<1+na); image_delta_y_best = []; end; na=na+1;
if (nargin<1+na); image_delta_z_best = []; end; na=na+1;

if isempty(N_wavelength); N_wavelength = 0; end;
if (N_wavelength<=0); N_wavelength = 0; end;
if isempty(image_delta_x_best); image_delta_x_best = 0; end;
if isempty(image_delta_y_best); image_delta_y_best = 0; end;
if isempty(image_delta_z_best); image_delta_z_best = 0; end;

if (verbose>1); disp(sprintf(' %% indices for counting arrays')); end;
n_lm_ = (l_max_+1).^2;
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
k_p_r_max = k_p_r_(end);
l_max_max = max(l_max_);
m_max_ = -l_max_max : +l_max_max;
n_m_max = length(m_max_);
dWtdkd__l_max_max = 2*l_max_max;

flag_prealigned = ...
  ~isempty(X_best) & ...
  ~isempty(flag_flip) & ...
  ~isempty(euler_polar_a_best) & ...
  ~isempty(euler_azimu_b_best) & ...
  ~isempty(euler_gamma_z_best) & ...
  1;

if flag_prealigned==0;
tmp_t = tic();
[ ...
 X_best ...
,flag_flip ...
,euler_polar_a_best ...
,euler_azimu_b_best ...
,euler_gamma_z_best ...
,image_delta_best_ ...
,~ ...
,~ ...
] = ...
register_spharm_to_spharm_wigner_wrap_1( ...
 n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,N_wavelength ...
,l_max_ ...
,a_k_Y_true_ ...
,b_k_Y_reco_ ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% registration time %0.2fs',tmp_t)); end;
image_delta_x_best = image_delta_best_(1+0);
image_delta_y_best = image_delta_best_(1+1);
image_delta_z_best = image_delta_best_(1+2);
end;%if flag_prealigned==0;
image_delta_best_ = [ image_delta_x_best ; image_delta_y_best ; image_delta_z_best ];

%%%%%%%%;
% Now, given X_flag_flip==1, form c_k_Y_reco_ = flipY(b_k_Y_reco_). ;
%%%%%%%%;
c_k_Y_reco_ = b_k_Y_reco_; if (flag_flip); c_k_Y_reco_ = flipY(n_k_p_r,l_max_,b_k_Y_reco_); end;
%%%%%%%%;
% To align the two, we translate c_k_Y_reco_ by image_delta_best_, and rotate a_k_Y_true_ by euler_best_. ;
% This is equivalent to rotating c_k_Y_reco_ by euler_best_inv_ (below). ;
%%%%%%%%;
if (fnorm(image_delta_best_)>0);
if (verbose); disp(sprintf(' %% translating to form c_k_Y_reco_ = T_(delta_best_) * c_k_Y_reco_')); end;
if (verbose); disp(sprintf(' %% delta_best_: %+0.6f , %+0.6f , %+0.6f',delta_best_)); end;
tmp_t = tic();
dWtdkd__ = dwignertdkd__(dWtdkd__l_max_max);
delta_best_p_r = fnorm(image_delta_best_);
Wt___ = expm_dwignertdkd__(dWtdkd__,n_k_p_r,k_p_r_,l_max_,delta_best_p_r);
Wt_ = wignert_ODE_0(dWtdkd__,Wt___,n_k_p_r,k_p_r_,l_max_,delta_best_p_r);
delta_z_c_ = transpose(image_delta_best_);
delta_z_p_r = sqrt(delta_z_c_(1+0).^2 + delta_z_c_(1+1).^2 + delta_z_c_(1+2).^2);
delta_z_p_01 = sqrt(delta_z_c_(1+0).^2 + delta_z_c_(1+1).^2);
delta_z_p_azimu_b = atan2(delta_z_c_(1+1),delta_z_c_(1+0));
delta_z_p_polar_a = atan2(delta_z_p_01,delta_z_c_(1+2));
delta_z_p_euler_pos_ = [0,+delta_z_p_polar_a,+delta_z_p_azimu_b];
delta_z_p_euler_neg_ = [-delta_z_p_azimu_b,-delta_z_p_polar_a,0];
delta_z_c_ = [cos(delta_z_p_azimu_b)*sin(delta_z_p_polar_a);sin(delta_z_p_azimu_b)*sin(delta_z_p_polar_a);cos(delta_z_p_polar_a)]*delta_z_p_r;
W_beta_neg__ = wignerd_b(l_max_max,delta_z_p_euler_neg_(1+1));
W_beta_pos__ = wignerd_b(l_max_max,delta_z_p_euler_pos_(1+1));
tmp_Y_form_ = c_k_Y_reco_;
tmp_Y_form_ = rotate_spharm_to_spharm_2(0,W_beta_neg__,n_k_p_r,k_p_r_,l_max_,tmp_Y_form_,delta_z_p_euler_neg_);
tmp_Y_form_ = Wt_*tmp_Y_form_; 
tmp_Y_form_ = rotate_spharm_to_spharm_2(0,W_beta_pos__,n_k_p_r,k_p_r_,l_max_,tmp_Y_form_,delta_z_p_euler_pos_);
c_k_Y_reco_ = tmp_Y_form_;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% translation time %0.2fs',tmp_t)); end;
end;%if (fnorm(image_delta_best_)>0);
%%%%%%%%;
euler_best_inv_ = [-euler_gamma_z_best,-euler_polar_a_best,-euler_azimu_b_best];
W_beta_best__ = wignerd_b(l_max_max,euler_best_inv_(1+1));
tmp_t = tic();
c_k_Y_reco_ = rotate_spharm_to_spharm_2(0,W_beta_best__,n_k_p_r,k_p_r_,l_max_,c_k_Y_reco_,euler_best_inv_);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% rotation time %0.2fs',tmp_t)); end;
%%%%%%%%;

if (verbose); disp(sprintf(' %% [finished spharm_register_and_rotate_2]')); end;
