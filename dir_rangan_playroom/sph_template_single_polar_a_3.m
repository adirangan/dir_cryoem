function ...
[ ...
 template_wkb__ ...
,V_lmm___ ...
,L_lm__ ...
,d0W_betazeta_mlm___ ...
,dtemplateda_wkb__ ...
,dtemplatedb_wkb__ ...
,dtemplatedc_wkb__ ...
,d1W_betazeta_mlm___ ...
,ddtemplatedaa_wkb__ ...
,ddtemplatedab_wkb__ ...
,ddtemplatedac_wkb__ ...
,ddtemplatedbb_wkb__ ...
,ddtemplatedbc_wkb__ ...
,ddtemplatedcc_wkb__ ...
,d2W_betazeta_mlm___ ...
] = ...
sph_template_single_polar_a_3( ...
 verbose ...
,l_max ...
,n_k_p_r ...
,a_k_Y_lkm___ ...
,n_w ...
,viewing_polar_a ...
,n_viewing_azimu_b ...
,viewing_azimu_b_ ...
,viewing_gamma_z ...
,V_lmm___ ...
,L_lm__ ...
,d0W_betazeta_mlm___ ...
,d1W_betazeta_mlm___ ...
,d2W_betazeta_mlm___ ...
);
% uses spherical-harmonic-expansions a_k_Y_yk__ to evaluate templates on a collection of points on spherical shells. ;
% each spherical-shell has the same resolution, determined by n_w_max:=n_w. ;
% ;
% inputs: ;
% ;
% verbose = integer verbosity_level. ;
% l_max = spherical harmonic order on each shell. ;
%         corresponds to n_lm = (l_max+1)^2 coefficients. ;
% n_k_p_r = integer number of shells. ;
% a_k_Y_lkm___ = complex array of size (1+l_max,n_k_p_r,n_m_val). ;
%              spherical-harmonic-coefficients of volume. ;
% n_w = integer number of angles per template-ring. ;
% viewing_polar_a = real polar_a value for each template. ;
% n_viewing_azimu_b = integer number of azimu_v values. ;
% viewing_azimu_b_ = real array of size(n_viewing_azimu_b,1). ;
%                    azimu_b values requested. ;
% viewing_gamma_z_ = real gamma_z value for each template. (typically 0.0). ;
% V_lmm___ = complex cell array of size 1+l_max;
%            V_mm__ = V_lmm___{1+l_val} = complex array of size (1+2*l_val,1+2*l_val);
%                     holds eigenvectors for y-component of total-angular-momentum-operator.
% L_lm__ = complex cell array of size 1+l_max;
%          L_m_ = L_lm__{1+l_val} = complex array of size (1+2*l_val,1);
%                 holds eigenvalues for y-component of total-angular-momentum-operator.
% d0W_betazeta_mlm___ = complex array of size (n_m_val,1+l_max,n_m_val). ;
%                     stores coefficients corresponding to slice- and rotation-operator. ;
% d1W_betazeta_mlm___ = complex array of size (n_m_val,1+l_max,n_m_val). ;
%                     stores coefficients corresponding to slice- and rotation-operator. ;
%                     stores first-derivative with respect to polar_a. ;
% d2W_betazeta_mlm___ = complex array of size (n_m_val,1+l_max,n_m_val). ;
%                     stores coefficients corresponding to slice- and rotation-operator. ;
%                     stores second-derivative with respect to polar_a. ;
% ;
% outputs: ;
% ;
% template_wkb__ = complex array of templates. ;
%                  template_wkb__(1+nw+nk_p_r*n_w_max,1+nazimu_b) ;
%                  stores template value for angle-index nw, radial-index nk_p_r, ;
%                  and viewing_azimu_b = viewing_azimu_b_(1+nazimu_b). ;
% dtemplatedx_wkb__ = complex array analogous to template_wkb__. ;
%                     stores first-derivative of template with respect to: ;
%                     x==a: polar_a ; %<-- note that the first-derivative with respect to polar_a has a different sign than wignerd_c produces. ;
%                     x==b: azimu_b ;
%                     x==c: gamma_z ;
% ddtemplatedxy_wkb__ = complex array analogous to template_wkb__. ;
%                       stores second-derivative of template with respect to x and y (see above). ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

str_thisfunction = 'sph_template_single_polar_a_3';

if nargin<1;
verbose = 2; nf=1;
if (verbose); disp(sprintf(' %% testing %s',str_thisfunction)); end;
n_k_p_r = 49;
k_p_r_max = 1;
k_p_r_ = ones(n_k_p_r,1);
weight_k_p_r_ = ones(n_k_p_r,1);
%%%%;
l_max = 48;
n_lm = (l_max+1)^2;
Y_l_val_ = zeros(n_lm,1);
Y_m_val_ = zeros(n_lm,1);
na=0;
for l_val=0:l_max;
for m_val=-l_val:+l_val;
Y_l_val_(1+na) = l_val;
Y_m_val_(1+na) = m_val;
na=na+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max;
%%%%;
l_max_ = l_max*ones(n_k_p_r,1);
l_max_max = max(l_max_);
m_max_ = -l_max_max:+l_max_max; n_m_max = numel(m_max_);
n_lm_ = n_lm*ones(n_k_p_r,1);
n_lm_max = max(n_lm_);
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
Y_l_val__ = repmat(Y_l_val_,[1,n_k_p_r]);
Y_m_val__ = repmat(Y_m_val_,[1,n_k_p_r]);
sigma_l = 15;
sigma_m = 12;
a_k_Y_yk_ = (mod(transpose([0:n_lm_sum-1]),89)-44)/89 + i*(mod(transpose([0:n_lm_sum-1]),97)-48)/97;
a_k_Y_yk_ = a_k_Y_yk_.*exp(-Y_l_val__(:).^2/(2*sigma_l^2)).*exp(-Y_m_val__(:).^2/(2*sigma_m^2));
%%%%;
if (verbose>1);
figure(1+nf);nf=nf+1; clf; figmed; 
subplot(1,2,1); plot(Y_l_val__(:),abs(a_k_Y_yk_(:)),'.'); xlabel('Y_l_val_','Interpreter','none'); ylabel('abs(a_k_Y_yk_)','Interpreter','none');
subplot(1,2,2); plot(Y_m_val__(:),abs(a_k_Y_yk_(:)),'.'); xlabel('Y_m_val_','Interpreter','none'); ylabel('abs(a_k_Y_yk_)','Interpreter','none');
end;%if (verbose>1); 
if (verbose>2); fprintf(1,' %% \t\n'); darray_printf_margin(real(a_k_Y_yk_),n_lm,n_k_p_r,' %% a_k_Y_yk_real___: '); end;
if (verbose>2); fprintf(1,' %% \t\n'); darray_printf_margin(imag(a_k_Y_yk_),n_lm,n_k_p_r,' %% a_k_Y_yk_imag___: '); end;
%%%%;
viewing_k_eq_d = 1/(4*pi);
template_k_eq_d = -1;
n_w_max = 1*98;
n_w = n_w_max; n_w_ = n_w_max*ones(n_k_p_r,1); n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_]);
%%%%%%%%;
a_k_Y_yk__ = zeros(n_lm_max,n_k_p_r);
a_k_Y_lmk___ = zeros(1+l_max_max,n_m_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
tmp_a_k_Y_lm_ = a_k_Y_yk_(1+tmp_index_);
a_k_Y_yk__(1:n_lm_(1+nk_p_r),1+nk_p_r) = tmp_a_k_Y_lm_;
tmp_a_k_Y_lm__ = zeros(1+l_max_max,n_m_max);
l_max = l_max_(1+nk_p_r);
na=0;
for l_val=0:l_max;
for m_val=-l_val:+l_val;
assert(na==l_val*(l_val+1)+m_val);
tmp_a_k_Y_lm__(1+l_val,1+l_max_max+m_val) = tmp_a_k_Y_lm_(1+na);
na=na+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max;
a_k_Y_lmk___(:,:,1+nk_p_r) = tmp_a_k_Y_lm__;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%;
flag_check=1;
if flag_check;
na=0;
for nk_p_r=0:n_k_p_r-1;
l_max=l_max_(1+nk_p_r);
for l_val=0:l_max;
for m_val=-l_val:+l_val;
assert(a_k_Y_yk_(1+n_lm_csum_(1+nk_p_r)+l_val*(l_val+1)+m_val)==a_k_Y_yk__(1+l_val*(l_val+1)+m_val,1+nk_p_r));
assert(a_k_Y_yk_(1+n_lm_csum_(1+nk_p_r)+l_val*(l_val+1)+m_val)==a_k_Y_lmk___(1+l_val,1+l_max_max+m_val,1+nk_p_r));
na=na+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max;
end;%for nk_p_r=0:n_k_p_r-1;
assert(na==n_lm_sum);
end;%if flag_check;
%%%%%%%%;
a_k_Y_lkm___ = permute(a_k_Y_lmk___,[1,3,2]);
%%%%%%%%;
n_viewing_all = 32;
viewing_weight_all_ = [];
viewing_polar_a = pi/3;
viewing_polar_a_all_ = viewing_polar_a*ones(n_viewing_all,1);
viewing_azimu_b_all_ = linspace(0,2*pi,1+n_viewing_all); viewing_azimu_b_all_ = transpose(viewing_azimu_b_all_(1:n_viewing_all));
viewing_azimu_b_ = viewing_azimu_b_all_;
n_viewing_polar_a = 1;
viewing_polar_a_ = viewing_polar_a;
n_viewing_azimu_b = n_viewing_all;
n_viewing_azimu_b_ = n_viewing_all;
viewing_k_eq_d = [];
template_k_eq_d = -1;
n_w_0in= n_w_max;
%%%%%%%%;
tmp_t = tic();
[ ...
 template_2_wkb___ ...
] = ...
pm_template_2( ...
 0*verbose ...
,l_max ...
,n_k_p_r ...
,a_k_Y_yk__ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,n_w_0in ...
,n_viewing_all ...
,viewing_azimu_b_all_ ...
,viewing_polar_a_all_ ...
,viewing_weight_all_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% pm_template_2: %0.2fs',tmp_t));
template_2_wkb__ = reshape(template_2_wkb___,[n_w_max*n_k_p_r,n_viewing_all]);
%%%%%%%%;
viewing_gamma_z = 0;
if ~exist('V_lmm___','var'); V_lmm___=[]; end;
if ~exist('L_lm__','var'); L_lm__=[]; end;
if ~exist('d0W_betazeta_mlm___','var'); d0W_betazeta_mlm___=[]; end;
if ~exist('d1W_betazeta_mlm___','var'); d1W_betazeta_mlm___=[]; end;
if ~exist('d2W_betazeta_mlm___','var'); d2W_betazeta_mlm___=[]; end;
tmp_t = tic();
[ ...
 template_3_wkb__ ...
,V_lmm___ ...
,L_lm__ ...
,d0W_betazeta_mlm___ ...
,~ ...
,~ ...
,~ ...
,d1W_betazeta_mlm___ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,d2W_betazeta_mlm___ ...
] = ...
sph_template_single_polar_a_3( ...
 0*verbose ...
,l_max_max ...
,n_k_p_r ...
,a_k_Y_lkm___ ...
,n_w_max ...
,viewing_polar_a ...
,n_viewing_azimu_b ...
,viewing_azimu_b_ ...
,viewing_gamma_z ...
,V_lmm___ ...
,L_lm__ ...
,d0W_betazeta_mlm___ ...
,d1W_betazeta_mlm___ ...
,d2W_betazeta_mlm___ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% sph_template_single_polar_a_3: %0.2fs',tmp_t));
%%%%;
if (verbose>1);
figure(1+nf);nf=nf+1; clf; figmed; 
subplot(1,2,1); plot(real(template_2_wkb__(:)),real(template_3_wkb__(:)),'.');
xlabel('real(template_2_wkb__)','Interpreter','none'); ylabel('real(template_3_wkb__)','Interpreter','none'); axis equal; grid on;
subplot(1,2,2); plot(imag(template_2_wkb__(:)),imag(template_3_wkb__(:)),'.');
xlabel('imag(template_2_wkb__)','Interpreter','none'); ylabel('imag(template_3_wkb__)','Interpreter','none'); axis equal; grid on;
end;%if (verbose>1); 
%%%%;
disp(sprintf(' %% template_3_wkb__ vs template_2_wkb__: %0.16f %<-- should be <1e-2', fnorm(template_3_wkb__ - template_2_wkb__)/fnorm(template_3_wkb__)));
%%%%%%%%;
% Now test first-derivative. ;
%%%%%%%%;
if (verbose); disp(sprintf(' %% %% %% %% %% %% %% %% ')); end;
if (verbose); disp(sprintf(' %% Now testing first-derivative: ')); end;
if (verbose); disp(sprintf(' %% %% %% %% %% %% %% %% ')); end;
[ ...
 template_mid_x_wkb__ ...
,~ ...
,~ ...
,~ ...
,dtemplateda_mid_wkb__ ...
,dtemplatedb_mid_wkb__ ...
,dtemplatedc_mid_wkb__ ...
,~ ...
,ddtemplatedaa_mid_wkb__ ...
,ddtemplatedab_mid_wkb__ ...
,ddtemplatedac_mid_wkb__ ...
,ddtemplatedbb_mid_wkb__ ...
,ddtemplatedbc_mid_wkb__ ...
,ddtemplatedcc_mid_wkb__ ...
,~ ...
] = ...
sph_template_single_polar_a_3( ...
 0*verbose ...
,l_max_max ...
,n_k_p_r ...
,a_k_Y_lkm___ ...
,n_w_max ...
,viewing_polar_a ...
,n_viewing_azimu_b ...
,viewing_azimu_b_ ...
,viewing_gamma_z ...
,V_lmm___ ...
,L_lm__ ...
,d0W_betazeta_mlm___ ...
,d1W_betazeta_mlm___ ...
,d2W_betazeta_mlm___ ...
);
da = pi*1e-4; db = 1.25*pi*1e-3; dc = 0.75*pi*1e-3;
[ ...
 template_pos_a_wkb__ ...
] = ...
sph_template_single_polar_a_3( ...
 0*verbose ...
,l_max_max ...
,n_k_p_r ...
,a_k_Y_lkm___ ...
,n_w_max ...
,viewing_polar_a+da ...
,n_viewing_azimu_b ...
,viewing_azimu_b_ ...
,viewing_gamma_z ...
,V_lmm___ ...
,L_lm__ ...
);
[ ...
 template_neg_a_wkb__ ...
] = ...
sph_template_single_polar_a_3( ...
 0*verbose ...
,l_max_max ...
,n_k_p_r ...
,a_k_Y_lkm___ ...
,n_w_max ...
,viewing_polar_a-da ...
,n_viewing_azimu_b ...
,viewing_azimu_b_ ...
,viewing_gamma_z ...
,V_lmm___ ...
,L_lm__ ...
);
[ ...
 template_pos_b_wkb__ ...
] = ...
sph_template_single_polar_a_3( ...
 0*verbose ...
,l_max_max ...
,n_k_p_r ...
,a_k_Y_lkm___ ...
,n_w_max ...
,viewing_polar_a ...
,n_viewing_azimu_b ...
,viewing_azimu_b_+db ...
,viewing_gamma_z ...
,V_lmm___ ...
,L_lm__ ...
);
[ ...
 template_neg_b_wkb__ ...
] = ...
sph_template_single_polar_a_3( ...
 0*verbose ...
,l_max_max ...
,n_k_p_r ...
,a_k_Y_lkm___ ...
,n_w_max ...
,viewing_polar_a ...
,n_viewing_azimu_b ...
,viewing_azimu_b_-db ...
,viewing_gamma_z ...
,V_lmm___ ...
,L_lm__ ...
);
[ ...
 template_pos_c_wkb__ ...
] = ...
sph_template_single_polar_a_3( ...
 0*verbose ...
,l_max_max ...
,n_k_p_r ...
,a_k_Y_lkm___ ...
,n_w_max ...
,viewing_polar_a ...
,n_viewing_azimu_b ...
,viewing_azimu_b_ ...
,viewing_gamma_z+dc ...
,V_lmm___ ...
,L_lm__ ...
);
[ ...
 template_neg_c_wkb__ ...
] = ...
sph_template_single_polar_a_3( ...
 0*verbose ...
,l_max_max ...
,n_k_p_r ...
,a_k_Y_lkm___ ...
,n_w_max ...
,viewing_polar_a ...
,n_viewing_azimu_b ...
,viewing_azimu_b_ ...
,viewing_gamma_z-dc ...
,V_lmm___ ...
,L_lm__ ...
);
dtemplateda_dif_wkb__ = (template_pos_a_wkb__ - template_neg_a_wkb__)./max(1e-12,2*da);
if (verbose); disp(sprintf(' %% dtemplateda_dif_wkb__ vs dtemplateda_mid_wkb__: %0.16f',fnorm(dtemplateda_dif_wkb__ - dtemplateda_mid_wkb__)/max(1e-12,fnorm(dtemplateda_dif_wkb__)))); end;
dtemplatedb_dif_wkb__ = (template_pos_b_wkb__ - template_neg_b_wkb__)./max(1e-12,2*db);
if (verbose); disp(sprintf(' %% dtemplatedb_dif_wkb__ vs dtemplatedb_mid_wkb__: %0.16f',fnorm(dtemplatedb_dif_wkb__ - dtemplatedb_mid_wkb__)/max(1e-12,fnorm(dtemplatedb_dif_wkb__)))); end;
dtemplatedc_dif_wkb__ = (template_pos_c_wkb__ - template_neg_c_wkb__)./max(1e-12,2*dc);
if (verbose); disp(sprintf(' %% dtemplatedc_dif_wkb__ vs dtemplatedc_mid_wkb__: %0.16f',fnorm(dtemplatedc_dif_wkb__ - dtemplatedc_mid_wkb__)/max(1e-12,fnorm(dtemplatedc_dif_wkb__)))); end;
if (verbose>1);
figure(1+nf);nf=nf+1; clf; figbig; p_row = 2; p_col = 3; np=0;
subplot(p_row,p_col,1+np+0*p_col); plot(real(dtemplateda_dif_wkb__(:)),real(dtemplateda_mid_wkb__(:)),'.');
xlabel('real(dtemplateda_dif_wkb__)','Interpreter','none'); ylabel('real(dtemplateda_mid_wkb__)','Interpreter','none'); axis equal; grid on;
subplot(p_row,p_col,1+np+1*p_col); plot(imag(dtemplateda_dif_wkb__(:)),imag(dtemplateda_mid_wkb__(:)),'.');
xlabel('imag(dtemplateda_dif_wkb__)','Interpreter','none'); ylabel('imag(dtemplateda_mid_wkb__)','Interpreter','none'); axis equal; grid on;
np=np+1;
subplot(p_row,p_col,1+np+0*p_col); plot(real(dtemplatedb_dif_wkb__(:)),real(dtemplatedb_mid_wkb__(:)),'.');
xlabel('real(dtemplatedb_dif_wkb__)','Interpreter','none'); ylabel('real(dtemplatedb_mid_wkb__)','Interpreter','none'); axis equal; grid on;
subplot(p_row,p_col,1+np+1*p_col); plot(imag(dtemplatedb_dif_wkb__(:)),imag(dtemplatedb_mid_wkb__(:)),'.');
xlabel('imag(dtemplatedb_dif_wkb__)','Interpreter','none'); ylabel('imag(dtemplatedb_mid_wkb__)','Interpreter','none'); axis equal; grid on;
np=np+1;
subplot(p_row,p_col,1+np+0*p_col); plot(real(dtemplatedc_dif_wkb__(:)),real(dtemplatedc_mid_wkb__(:)),'.');
xlabel('real(dtemplatedc_dif_wkb__)','Interpreter','none'); ylabel('real(dtemplatedc_mid_wkb__)','Interpreter','none'); axis equal; grid on;
subplot(p_row,p_col,1+np+1*p_col); plot(imag(dtemplatedc_dif_wkb__(:)),imag(dtemplatedc_mid_wkb__(:)),'.');
xlabel('imag(dtemplatedc_dif_wkb__)','Interpreter','none'); ylabel('imag(dtemplatedc_mid_wkb__)','Interpreter','none'); axis equal; grid on;
np=np+1;
end;%if (verbose>1); 
%%%%%%%%;
% Now test second-derivative. ;
%%%%%%%%;
if (verbose); disp(sprintf(' %% %% %% %% %% %% %% %% ')); end;
if (verbose); disp(sprintf(' %% Now testing second-derivative: ')); end;
if (verbose); disp(sprintf(' %% %% %% %% %% %% %% %% ')); end;
ddtemplatedaa_dif_wkb__ = (template_pos_a_wkb__ - 2*template_mid_x_wkb__ + template_neg_a_wkb__)./max(1e-12,da.^2);
if (verbose); disp(sprintf(' %% ddtemplatedaa_dif_wkb__ vs ddtemplatedaa_mid_wkb__: %0.16f',fnorm(ddtemplatedaa_dif_wkb__ - ddtemplatedaa_mid_wkb__)/max(1e-12,fnorm(ddtemplatedaa_dif_wkb__)))); end;
ddtemplatedbb_dif_wkb__ = (template_pos_b_wkb__ - 2*template_mid_x_wkb__ + template_neg_b_wkb__)./max(1e-12,db.^2);
if (verbose); disp(sprintf(' %% ddtemplatedbb_dif_wkb__ vs ddtemplatedbb_mid_wkb__: %0.16f',fnorm(ddtemplatedbb_dif_wkb__ - ddtemplatedbb_mid_wkb__)/max(1e-12,fnorm(ddtemplatedbb_dif_wkb__)))); end;
ddtemplatedcc_dif_wkb__ = (template_pos_c_wkb__ - 2*template_mid_x_wkb__ + template_neg_c_wkb__)./max(1e-12,dc.^2);
if (verbose); disp(sprintf(' %% ddtemplatedcc_dif_wkb__ vs ddtemplatedcc_mid_wkb__: %0.16f',fnorm(ddtemplatedcc_dif_wkb__ - ddtemplatedcc_mid_wkb__)/max(1e-12,fnorm(ddtemplatedcc_dif_wkb__)))); end;
%%%%;
[ ...
 template_pos_a_pos_b_wkb__ ...
] = ...
sph_template_single_polar_a_3( ...
 0*verbose ...
,l_max_max ...
,n_k_p_r ...
,a_k_Y_lkm___ ...
,n_w_max ...
,viewing_polar_a+da ...
,n_viewing_azimu_b ...
,viewing_azimu_b_+db ...
,viewing_gamma_z ...
,V_lmm___ ...
,L_lm__ ...
);
[ ...
 template_pos_a_neg_b_wkb__ ...
] = ...
sph_template_single_polar_a_3( ...
 0*verbose ...
,l_max_max ...
,n_k_p_r ...
,a_k_Y_lkm___ ...
,n_w_max ...
,viewing_polar_a+da ...
,n_viewing_azimu_b ...
,viewing_azimu_b_-db ...
,viewing_gamma_z ...
,V_lmm___ ...
,L_lm__ ...
);
[ ...
 template_neg_a_pos_b_wkb__ ...
] = ...
sph_template_single_polar_a_3( ...
 0*verbose ...
,l_max_max ...
,n_k_p_r ...
,a_k_Y_lkm___ ...
,n_w_max ...
,viewing_polar_a-da ...
,n_viewing_azimu_b ...
,viewing_azimu_b_+db ...
,viewing_gamma_z ...
,V_lmm___ ...
,L_lm__ ...
);
[ ...
 template_neg_a_neg_b_wkb__ ...
] = ...
sph_template_single_polar_a_3( ...
 0*verbose ...
,l_max_max ...
,n_k_p_r ...
,a_k_Y_lkm___ ...
,n_w_max ...
,viewing_polar_a-da ...
,n_viewing_azimu_b ...
,viewing_azimu_b_-db ...
,viewing_gamma_z ...
,V_lmm___ ...
,L_lm__ ...
);
ddtemplatedab_dif_wkb__ = (template_pos_a_pos_b_wkb__ + template_neg_a_neg_b_wkb__ - template_pos_a_neg_b_wkb__ - template_neg_a_pos_b_wkb__)./max(1e-12,4*da*db);
if (verbose); disp(sprintf(' %% ddtemplatedab_dif_wkb__ vs ddtemplatedab_mid_wkb__: %0.16f',fnorm(ddtemplatedab_dif_wkb__ - ddtemplatedab_mid_wkb__)/max(1e-12,fnorm(ddtemplatedab_dif_wkb__)))); end;
%%%%;
[ ...
 template_pos_a_pos_c_wkb__ ...
] = ...
sph_template_single_polar_a_3( ...
 0*verbose ...
,l_max_max ...
,n_k_p_r ...
,a_k_Y_lkm___ ...
,n_w_max ...
,viewing_polar_a+da ...
,n_viewing_azimu_b ...
,viewing_azimu_b_ ...
,viewing_gamma_z+dc ...
,V_lmm___ ...
,L_lm__ ...
);
[ ...
 template_pos_a_neg_c_wkb__ ...
] = ...
sph_template_single_polar_a_3( ...
 0*verbose ...
,l_max_max ...
,n_k_p_r ...
,a_k_Y_lkm___ ...
,n_w_max ...
,viewing_polar_a+da ...
,n_viewing_azimu_b ...
,viewing_azimu_b_ ...
,viewing_gamma_z-dc ...
,V_lmm___ ...
,L_lm__ ...
);
[ ...
 template_neg_a_pos_c_wkb__ ...
] = ...
sph_template_single_polar_a_3( ...
 0*verbose ...
,l_max_max ...
,n_k_p_r ...
,a_k_Y_lkm___ ...
,n_w_max ...
,viewing_polar_a-da ...
,n_viewing_azimu_b ...
,viewing_azimu_b_ ...
,viewing_gamma_z+dc ...
,V_lmm___ ...
,L_lm__ ...
);
[ ...
 template_neg_a_neg_c_wkb__ ...
] = ...
sph_template_single_polar_a_3( ...
 0*verbose ...
,l_max_max ...
,n_k_p_r ...
,a_k_Y_lkm___ ...
,n_w_max ...
,viewing_polar_a-da ...
,n_viewing_azimu_b ...
,viewing_azimu_b_ ...
,viewing_gamma_z-dc ...
,V_lmm___ ...
,L_lm__ ...
);
ddtemplatedac_dif_wkb__ = (template_pos_a_pos_c_wkb__ + template_neg_a_neg_c_wkb__ - template_pos_a_neg_c_wkb__ - template_neg_a_pos_c_wkb__)./max(1e-12,4*da*dc);
if (verbose); disp(sprintf(' %% ddtemplatedac_dif_wkb__ vs ddtemplatedac_mid_wkb__: %0.16f',fnorm(ddtemplatedac_dif_wkb__ - ddtemplatedac_mid_wkb__)/max(1e-12,fnorm(ddtemplatedac_dif_wkb__)))); end;
%%%%;
[ ...
 template_pos_b_pos_c_wkb__ ...
] = ...
sph_template_single_polar_a_3( ...
 0*verbose ...
,l_max_max ...
,n_k_p_r ...
,a_k_Y_lkm___ ...
,n_w_max ...
,viewing_polar_a ...
,n_viewing_azimu_b ...
,viewing_azimu_b_+db ...
,viewing_gamma_z+dc ...
,V_lmm___ ...
,L_lm__ ...
);
[ ...
 template_pos_b_neg_c_wkb__ ...
] = ...
sph_template_single_polar_a_3( ...
 0*verbose ...
,l_max_max ...
,n_k_p_r ...
,a_k_Y_lkm___ ...
,n_w_max ...
,viewing_polar_a ...
,n_viewing_azimu_b ...
,viewing_azimu_b_+db ...
,viewing_gamma_z-dc ...
,V_lmm___ ...
,L_lm__ ...
);
[ ...
 template_neg_b_pos_c_wkb__ ...
] = ...
sph_template_single_polar_a_3( ...
 0*verbose ...
,l_max_max ...
,n_k_p_r ...
,a_k_Y_lkm___ ...
,n_w_max ...
,viewing_polar_a ...
,n_viewing_azimu_b ...
,viewing_azimu_b_-db ...
,viewing_gamma_z+dc ...
,V_lmm___ ...
,L_lm__ ...
);
[ ...
 template_neg_b_neg_c_wkb__ ...
] = ...
sph_template_single_polar_a_3( ...
 0*verbose ...
,l_max_max ...
,n_k_p_r ...
,a_k_Y_lkm___ ...
,n_w_max ...
,viewing_polar_a ...
,n_viewing_azimu_b ...
,viewing_azimu_b_-db ...
,viewing_gamma_z-dc ...
,V_lmm___ ...
,L_lm__ ...
);
ddtemplatedbc_dif_wkb__ = (template_pos_b_pos_c_wkb__ + template_neg_b_neg_c_wkb__ - template_pos_b_neg_c_wkb__ - template_neg_b_pos_c_wkb__)./max(1e-12,4*db*dc);
if (verbose); disp(sprintf(' %% ddtemplatedbc_dif_wkb__ vs ddtemplatedbc_mid_wkb__: %0.16f',fnorm(ddtemplatedbc_dif_wkb__ - ddtemplatedbc_mid_wkb__)/max(1e-12,fnorm(ddtemplatedbc_dif_wkb__)))); end;
%%%%%%%%;
disp(sprintf(' %% returning')); return;
end;% if nargin<7;
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

na=0;
if (nargin<1+na); verbose=[]; end; na=na+1;
if (nargin<1+na); l_max=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); a_k_Y_lkm___=[]; end; na=na+1;
if (nargin<1+na); n_w=[]; end; na=na+1;
if (nargin<1+na); viewing_polar_a=[]; end; na=na+1;
if (nargin<1+na); n_viewing_azimu_b=[]; end; na=na+1;
if (nargin<1+na); viewing_azimu_b_=[]; end; na=na+1;
if (nargin<1+na); viewing_gamma_z=[]; end; na=na+1;
if (nargin<1+na); V_lmm___=[]; end; na=na+1;
if (nargin<1+na); L_lm__=[]; end; na=na+1;
if (nargin<1+na); d0W_betazeta_mlm___=[]; end; na=na+1;
if (nargin<1+na); d1W_betazeta_mlm___=[]; end; na=na+1;
if (nargin<1+na); d2W_betazeta_mlm___=[]; end; na=na+1;

flag_d0 = 1;
flag_d1 = (nargout>=5);
flag_d2 = (nargout>=9);
tolerance_xxnufft = 1e-6;

l_max_max = l_max;
n_lm = (l_max+1)^2;
m_max_ = -l_max : +l_max; m0_val_ = m_max_;
n_m_max = length(m_max_); %<-- 2*l_max+1;
n_w_max = n_w; n_w_sum = n_w_max*n_k_p_r; n_w_ = n_w_max*ones(n_k_p_r,1);
l_max_ = l_max_max*ones(n_k_p_r,1);
n_lm_max = n_lm;
if (verbose); disp(sprintf(' %% l_max %d n_lm %d',l_max,n_lm)); end;

if isempty(viewing_gamma_z); viewing_gamma_z = 0.0; end;

%%%%%%%%;
% Now we set up a template-operator ;
% for the collection of viewing_azimu_b associated with each unique viewing_polar_a. ;
%%%%%%%%;
template_wkb___ = zeros(n_w_max,n_k_p_r,n_viewing_azimu_b);
%%%%%%%%;
if (flag_d0 & isempty(d0W_betazeta_mlm___)) | (flag_d1 & isempty(d1W_betazeta_mlm___)) | (flag_d2 & isempty(d2W_betazeta_mlm___)) ;
%%%%;
tmp_t = tic();
[d0W_,d1W_,d2W_] = dwignerdda_b(l_max_max,-viewing_polar_a); %<-- Older 3-term recurrence: unstable for l_max_max>= 150 or so. ;
%[d0W_,V_lmm___,L_lm__,d1W_,d2W_] = wignerd_c(l_max_max,-viewing_polar_a,V_lmm___,L_lm__); %<-- Newer method based on eigenstructure of angular-momentum-operator (see Feng, Wang, Yang, Jin, 2015). Stable for l_max_max <= 400+ or so. ;
d0W_beta__ = d0W_;
d1W_beta__ = cellfun(@(x) -x,d1W_,'UniformOutput',0) ; %<-- choice of sign compatible with application of polar_a. ;
d2W_beta__ = d2W_;
%%%%;
zeta_lm__ = zeros(1+l_max_max,n_m_max);
d0y_lm__ = squeeze(ylgndr_1(l_max,0));
flag_check=0;
for l_val=0:l_max_max;
if flag_check; a1=((2*l_val+1)/(4*pi)); Llm__ = legendre(l_val,0,'unnorm'); end;
for m_val=-l_val:+l_val;
if flag_check;
if (l_val >0); Llm_ = Llm__(1+abs(m_val),:); end; if (l_val==0); Llm_ = Llm__; end; assert(numel(Llm_)==1);
a2=exp(0.5*lfactorial(l_val-abs(m_val)) - 0.5*lfactorial(l_val+abs(m_val))); c=sqrt(a1)*a2; s=1; % retain original phase ;
if (verbose>-1); disp(sprintf(' %% fnorm(s*c*Llm_(1+0) - d0y_lm__(1+l_val,1+abs(m_val))/sqrt(4*pi)): %0.16f',fnorm(s*c*Llm_(1+0) - d0y_lm__(1+l_val,1+abs(m_val))/sqrt(4*pi)))); end;
end;%if flag_check;
zeta_lm__(1+l_val,1+l_max_max+m_val) = d0y_lm__(1+l_val,1+abs(m_val))/sqrt(4*pi);
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max_max;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% %% d0W_beta__ and zeta_lm__: %0.2fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
if flag_d0; d0W_betazeta_mlm___ = zeros(n_m_max,1+l_max_max,n_m_max); end;
if flag_d1; d1W_betazeta_mlm___ = zeros(n_m_max,1+l_max_max,n_m_max); end;
if flag_d2; d2W_betazeta_mlm___ = zeros(n_m_max,1+l_max_max,n_m_max); end;
for l_val=0:l_max_max;
for m0_val=-l_val:+l_val;
for m1_val=-l_val:+l_val;
if flag_d0;
d0W_betazeta_mlm___(1+l_max_max+m0_val,1+l_val,1+l_max_max+m1_val) = ...
 d0W_beta__{1+l_val}(1+l_val+m0_val,1+l_val+m1_val) ...
*zeta_lm__(1+l_val,1+l_max_max+m0_val) ...
 ;
end;%if flag_d0;
if flag_d1;
d1W_betazeta_mlm___(1+l_max_max+m0_val,1+l_val,1+l_max_max+m1_val) = ...
 d1W_beta__{1+l_val}(1+l_val+m0_val,1+l_val+m1_val) ...
*zeta_lm__(1+l_val,1+l_max_max+m0_val) ...
 ;
end;%if flag_d1;
if flag_d2;
d2W_betazeta_mlm___(1+l_max_max+m0_val,1+l_val,1+l_max_max+m1_val) = ...
 d2W_beta__{1+l_val}(1+l_val+m0_val,1+l_val+m1_val) ...
*zeta_lm__(1+l_val,1+l_max_max+m0_val) ...
 ;
end;%if flag_d2;
end;%for m1_val=-l_val:+l_val;
end;%for m0_val=-l_val:+l_val;
end;%for l_val=0:l_max_max;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% %% d0W_betazeta_mlm___: %0.2fs',tmp_t)); end;
%%%%;
end;%if (flag_d0 & isempty(d0W_betazeta_mlm___)) | (flag_d1 & isempty(d1W_betazeta_mlm___)) | (flag_d2 & isempty(d2W_betazeta_mlm___)) ;
%%%%%%%%;
if flag_d0;
tmp_t = tic();
d0W_caza_mkm___ = zeros(n_m_max,n_k_p_r,n_m_max); %<-- diag(exp(+i*m0_val_*viewing_gamma_z))*d0W_betazeta_ml__*a_k_Y_lk__ for each m1_val. ;
for m1_val=-l_max_max:+l_max_max;
d0W_caza_mkm___(:,:,1+l_max_max+m1_val) = ...
 (pi^2) ...
 *diag(exp(-i*m_max_*(-viewing_gamma_z))) ...
 *reshape(d0W_betazeta_mlm___(:,:,1+l_max_max+m1_val),[n_m_max,1+l_max_max]) ...
 *reshape(a_k_Y_lkm___(:,:,1+l_max_max+m1_val),[1+l_max_max,n_k_p_r]) ...
 ;
end;%for m1_val=-l_max_max:+l_max_max;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% %% d0W_caza_mkm___: %0.2fs',tmp_t)); end;
end;%if flag_d0;
%%%%;
if flag_d1;
tmp_t = tic();
d1Wda_caza_mkm___ = zeros(n_m_max,n_k_p_r,n_m_max);  %<-- diag((+i*m0_val_).^0.*exp(+i*m0_val_*viewing_gamma_z))*d1W_betazeta_ml__*a_k_Y_lk__*(+i*m1_val).^0 for each m1_val. ;
d1Wdb_caza_mkm___ = zeros(n_m_max,n_k_p_r,n_m_max);  %<-- diag((+i*m0_val_).^0.*exp(+i*m0_val_*viewing_gamma_z))*d0W_betazeta_ml__*a_k_Y_lk__*(+i*m1_val).^1 for each m1_val. ;
d1Wdc_caza_mkm___ = zeros(n_m_max,n_k_p_r,n_m_max);  %<-- diag((+i*m0_val_).^1.*exp(+i*m0_val_*viewing_gamma_z))*d0W_betazeta_ml__*a_k_Y_lk__*(+i*m1_val).^0 for each m1_val. ;
for m1_val=-l_max_max:+l_max_max;
d1Wda_caza_mkm___(:,:,1+l_max_max+m1_val) = ...
 (pi^2) ...
 *diag((+i*m0_val_).^0.*exp(-i*m_max_*(-viewing_gamma_z))) ...
 *reshape(d1W_betazeta_mlm___(:,:,1+l_max_max+m1_val),[n_m_max,1+l_max_max]) ...
 *reshape(a_k_Y_lkm___(:,:,1+l_max_max+m1_val),[1+l_max_max,n_k_p_r]) ...
 *(+i*m1_val).^0 ...
 ;
d1Wdb_caza_mkm___(:,:,1+l_max_max+m1_val) = ...
 (pi^2) ...
 *diag((+i*m0_val_).^0.*exp(-i*m_max_*(-viewing_gamma_z))) ...
 *reshape(d0W_betazeta_mlm___(:,:,1+l_max_max+m1_val),[n_m_max,1+l_max_max]) ...
 *reshape(a_k_Y_lkm___(:,:,1+l_max_max+m1_val),[1+l_max_max,n_k_p_r]) ...
 *(+i*m1_val).^1 ...
 ;
d1Wdc_caza_mkm___(:,:,1+l_max_max+m1_val) = ...
 (pi^2) ...
 *diag((+i*m0_val_).^1.*exp(-i*m_max_*(-viewing_gamma_z))) ...
 *reshape(d0W_betazeta_mlm___(:,:,1+l_max_max+m1_val),[n_m_max,1+l_max_max]) ...
 *reshape(a_k_Y_lkm___(:,:,1+l_max_max+m1_val),[1+l_max_max,n_k_p_r]) ...
 *(+i*m1_val).^0 ...
 ;
end;%for m1_val=-l_max_max:+l_max_max;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% %% d1Wdx_caza_mkm___: %0.2fs',tmp_t)); end;
end;%if flag_d1;
%%%%;
if flag_d2;
tmp_t = tic();
d2Wdaa_caza_mkm___ = zeros(n_m_max,n_k_p_r,n_m_max); %<-- diag((+i*m0_val_).^0.*exp(+i*m0_val_*viewing_gamma_z))*d2W_betazeta_ml__*a_k_Y_lk__*(+i*m1_val).^0 for each m1_val. ;
d2Wdab_caza_mkm___ = zeros(n_m_max,n_k_p_r,n_m_max); %<-- diag((+i*m0_val_).^0.*exp(+i*m0_val_*viewing_gamma_z))*d1W_betazeta_ml__*a_k_Y_lk__*(+i*m1_val).^1 for each m1_val. ;
d2Wdac_caza_mkm___ = zeros(n_m_max,n_k_p_r,n_m_max); %<-- diag((+i*m0_val_).^1.*exp(+i*m0_val_*viewing_gamma_z))*d1W_betazeta_ml__*a_k_Y_lk__*(+i*m1_val).^0 for each m1_val. ;
d2Wdbb_caza_mkm___ = zeros(n_m_max,n_k_p_r,n_m_max); %<-- diag((+i*m0_val_).^0.*exp(+i*m0_val_*viewing_gamma_z))*d0W_betazeta_ml__*a_k_Y_lk__*(+i*m1_val).^2 for each m1_val. ;
d2Wdbc_caza_mkm___ = zeros(n_m_max,n_k_p_r,n_m_max); %<-- diag((+i*m0_val_).^1.*exp(+i*m0_val_*viewing_gamma_z))*d0W_betazeta_ml__*a_k_Y_lk__*(+i*m1_val).^1 for each m1_val. ;
d2Wdcc_caza_mkm___ = zeros(n_m_max,n_k_p_r,n_m_max); %<-- diag((+i*m0_val_).^2.*exp(+i*m0_val_*viewing_gamma_z))*d0W_betazeta_ml__*a_k_Y_lk__*(+i*m1_val).^0 for each m1_val. ;
for m1_val=-l_max_max:+l_max_max;
d2Wdaa_caza_mkm___(:,:,1+l_max_max+m1_val) = ...
 (pi^2) ...
 *diag((+i*m0_val_).^0.*exp(-i*m_max_*(-viewing_gamma_z))) ...
 *reshape(d2W_betazeta_mlm___(:,:,1+l_max_max+m1_val),[n_m_max,1+l_max_max]) ...
 *reshape(a_k_Y_lkm___(:,:,1+l_max_max+m1_val),[1+l_max_max,n_k_p_r]) ...
 *(+i*m1_val).^0 ...
 ;
d2Wdab_caza_mkm___(:,:,1+l_max_max+m1_val) = ...
 (pi^2) ...
 *diag((+i*m0_val_).^0.*exp(-i*m_max_*(-viewing_gamma_z))) ...
 *reshape(d1W_betazeta_mlm___(:,:,1+l_max_max+m1_val),[n_m_max,1+l_max_max]) ...
 *reshape(a_k_Y_lkm___(:,:,1+l_max_max+m1_val),[1+l_max_max,n_k_p_r]) ...
 *(+i*m1_val).^1 ...
 ;
d2Wdac_caza_mkm___(:,:,1+l_max_max+m1_val) = ...
 (pi^2) ...
 *diag((+i*m0_val_).^1.*exp(-i*m_max_*(-viewing_gamma_z))) ...
 *reshape(d1W_betazeta_mlm___(:,:,1+l_max_max+m1_val),[n_m_max,1+l_max_max]) ...
 *reshape(a_k_Y_lkm___(:,:,1+l_max_max+m1_val),[1+l_max_max,n_k_p_r]) ...
 *(+i*m1_val).^0 ...
 ;
d2Wdbb_caza_mkm___(:,:,1+l_max_max+m1_val) = ...
 (pi^2) ...
 *diag((+i*m0_val_).^0.*exp(-i*m_max_*(-viewing_gamma_z))) ...
 *reshape(d0W_betazeta_mlm___(:,:,1+l_max_max+m1_val),[n_m_max,1+l_max_max]) ...
 *reshape(a_k_Y_lkm___(:,:,1+l_max_max+m1_val),[1+l_max_max,n_k_p_r]) ...
 *(+i*m1_val).^2 ...
 ;
d2Wdbc_caza_mkm___(:,:,1+l_max_max+m1_val) = ...
 (pi^2) ...
 *diag((+i*m0_val_).^1.*exp(-i*m_max_*(-viewing_gamma_z))) ...
 *reshape(d0W_betazeta_mlm___(:,:,1+l_max_max+m1_val),[n_m_max,1+l_max_max]) ...
 *reshape(a_k_Y_lkm___(:,:,1+l_max_max+m1_val),[1+l_max_max,n_k_p_r]) ...
 *(+i*m1_val).^1 ...
 ;
d2Wdcc_caza_mkm___(:,:,1+l_max_max+m1_val) = ...
 (pi^2) ...
 *diag((+i*m0_val_).^2.*exp(-i*m_max_*(-viewing_gamma_z))) ...
 *reshape(d0W_betazeta_mlm___(:,:,1+l_max_max+m1_val),[n_m_max,1+l_max_max]) ...
 *reshape(a_k_Y_lkm___(:,:,1+l_max_max+m1_val),[1+l_max_max,n_k_p_r]) ...
 *(+i*m1_val).^0 ...
 ;
end;%for m1_val=-l_max_max:+l_max_max;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% %% d2Wdxy_caza_mkm___: %0.2fs',tmp_t)); end;
end;%if flag_d2;
%%%%%%%%;
if flag_d0;
tmp_t = tic();
%d0W_caza_mmk___ = permute(d0W_caza_mkm___,[3,1,2]);
%d0W_caza_bmk___ = reshape(xxnufft1d2(n_viewing_azimu_b,viewing_azimu_b_,+1,tolerance_xxnufft,n_m_max,reshape(d0W_caza_mmk___,[n_m_max,n_m_max*n_k_p_r])),[n_viewing_azimu_b,n_m_max,n_k_p_r]);
%d0W_caza_mkb___ = permute(d0W_caza_bmk___,[2,3,1]);
d0W_caza_mkb___ = permute(reshape(xxnufft1d2(n_viewing_azimu_b,viewing_azimu_b_,+1,tolerance_xxnufft,n_m_max,reshape(permute(d0W_caza_mkm___,[3,1,2]),[n_m_max,n_m_max*n_k_p_r])),[n_viewing_azimu_b,n_m_max,n_k_p_r]),[2,3,1]);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% %% d0W_caza_mkb___: %0.2fs',tmp_t)); end;
end;%if flag_d0;
%%%%;
if flag_d1;
tmp_t = tic();
d1Wda_caza_mkb___ = permute(reshape(xxnufft1d2(n_viewing_azimu_b,viewing_azimu_b_,+1,tolerance_xxnufft,n_m_max,reshape(permute(d1Wda_caza_mkm___,[3,1,2]),[n_m_max,n_m_max*n_k_p_r])),[n_viewing_azimu_b,n_m_max,n_k_p_r]),[2,3,1]);
d1Wdb_caza_mkb___ = permute(reshape(xxnufft1d2(n_viewing_azimu_b,viewing_azimu_b_,+1,tolerance_xxnufft,n_m_max,reshape(permute(d1Wdb_caza_mkm___,[3,1,2]),[n_m_max,n_m_max*n_k_p_r])),[n_viewing_azimu_b,n_m_max,n_k_p_r]),[2,3,1]);
d1Wdc_caza_mkb___ = permute(reshape(xxnufft1d2(n_viewing_azimu_b,viewing_azimu_b_,+1,tolerance_xxnufft,n_m_max,reshape(permute(d1Wdc_caza_mkm___,[3,1,2]),[n_m_max,n_m_max*n_k_p_r])),[n_viewing_azimu_b,n_m_max,n_k_p_r]),[2,3,1]);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% %% d1Wdx_caza_mkb___: %0.2fs',tmp_t)); end;
end;%if flag_d1;
%%%%;
if flag_d2;
tmp_t = tic();
d2Wdaa_caza_mkb___ = permute(reshape(xxnufft1d2(n_viewing_azimu_b,viewing_azimu_b_,+1,tolerance_xxnufft,n_m_max,reshape(permute(d2Wdaa_caza_mkm___,[3,1,2]),[n_m_max,n_m_max*n_k_p_r])),[n_viewing_azimu_b,n_m_max,n_k_p_r]),[2,3,1]);
d2Wdab_caza_mkb___ = permute(reshape(xxnufft1d2(n_viewing_azimu_b,viewing_azimu_b_,+1,tolerance_xxnufft,n_m_max,reshape(permute(d2Wdab_caza_mkm___,[3,1,2]),[n_m_max,n_m_max*n_k_p_r])),[n_viewing_azimu_b,n_m_max,n_k_p_r]),[2,3,1]);
d2Wdac_caza_mkb___ = permute(reshape(xxnufft1d2(n_viewing_azimu_b,viewing_azimu_b_,+1,tolerance_xxnufft,n_m_max,reshape(permute(d2Wdac_caza_mkm___,[3,1,2]),[n_m_max,n_m_max*n_k_p_r])),[n_viewing_azimu_b,n_m_max,n_k_p_r]),[2,3,1]);
d2Wdbb_caza_mkb___ = permute(reshape(xxnufft1d2(n_viewing_azimu_b,viewing_azimu_b_,+1,tolerance_xxnufft,n_m_max,reshape(permute(d2Wdbb_caza_mkm___,[3,1,2]),[n_m_max,n_m_max*n_k_p_r])),[n_viewing_azimu_b,n_m_max,n_k_p_r]),[2,3,1]);
d2Wdbc_caza_mkb___ = permute(reshape(xxnufft1d2(n_viewing_azimu_b,viewing_azimu_b_,+1,tolerance_xxnufft,n_m_max,reshape(permute(d2Wdbc_caza_mkm___,[3,1,2]),[n_m_max,n_m_max*n_k_p_r])),[n_viewing_azimu_b,n_m_max,n_k_p_r]),[2,3,1]);
d2Wdcc_caza_mkb___ = permute(reshape(xxnufft1d2(n_viewing_azimu_b,viewing_azimu_b_,+1,tolerance_xxnufft,n_m_max,reshape(permute(d2Wdcc_caza_mkm___,[3,1,2]),[n_m_max,n_m_max*n_k_p_r])),[n_viewing_azimu_b,n_m_max,n_k_p_r]),[2,3,1]);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% %% d2Wdxy_caza_mkb___: %0.2fs',tmp_t)); end;
end;%if flag_d2;
%%%%%%%%;
factor_scaling = sqrt(n_w_max)/pi^2;
%%%%%%%%;
if flag_d0;
tmp_t = tic();
R_k_q_sub_wkb___ = zeros(n_w_max,n_k_p_r,n_viewing_azimu_b);
R_k_q_sub_wkb___(1+[0:l_max_max],:,:) = d0W_caza_mkb___(1+l_max_max+[0:l_max_max],:,:);
R_k_q_sub_wkb___(1+n_w_max-l_max_max+[0:l_max_max-1],:,:) = d0W_caza_mkb___(1+[0:l_max_max-1],:,:);
R_k_p_sub_wkb___ = factor_scaling * interp_q_to_p_block_0(n_k_p_r,n_w_,n_w_sum,R_k_q_sub_wkb___);
template_wkb__ = reshape(R_k_p_sub_wkb___,[n_w_max*n_k_p_r,n_viewing_azimu_b]);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% %% template_wkb___: %0.2fs',tmp_t)); end;
end;%if flag_d0;
%%%%;
if flag_d1;
tmp_t = tic();
R_k_q_sub_wkb___ = zeros(n_w_max,n_k_p_r,n_viewing_azimu_b);
R_k_q_sub_wkb___(1+[0:l_max_max],:,:) = d1Wda_caza_mkb___(1+l_max_max+[0:l_max_max],:,:);
R_k_q_sub_wkb___(1+n_w_max-l_max_max+[0:l_max_max-1],:,:) = d1Wda_caza_mkb___(1+[0:l_max_max-1],:,:);
R_k_p_sub_wkb___ = factor_scaling * interp_q_to_p_block_0(n_k_p_r,n_w_,n_w_sum,R_k_q_sub_wkb___);
dtemplateda_wkb__ = reshape(R_k_p_sub_wkb___,[n_w_max*n_k_p_r,n_viewing_azimu_b]);
R_k_q_sub_wkb___ = zeros(n_w_max,n_k_p_r,n_viewing_azimu_b);
R_k_q_sub_wkb___(1+[0:l_max_max],:,:) = d1Wdb_caza_mkb___(1+l_max_max+[0:l_max_max],:,:);
R_k_q_sub_wkb___(1+n_w_max-l_max_max+[0:l_max_max-1],:,:) = d1Wdb_caza_mkb___(1+[0:l_max_max-1],:,:);
R_k_p_sub_wkb___ = factor_scaling * interp_q_to_p_block_0(n_k_p_r,n_w_,n_w_sum,R_k_q_sub_wkb___);
dtemplatedb_wkb__ = reshape(R_k_p_sub_wkb___,[n_w_max*n_k_p_r,n_viewing_azimu_b]);
R_k_q_sub_wkb___ = zeros(n_w_max,n_k_p_r,n_viewing_azimu_b);
R_k_q_sub_wkb___(1+[0:l_max_max],:,:) = d1Wdc_caza_mkb___(1+l_max_max+[0:l_max_max],:,:);
R_k_q_sub_wkb___(1+n_w_max-l_max_max+[0:l_max_max-1],:,:) = d1Wdc_caza_mkb___(1+[0:l_max_max-1],:,:);
R_k_p_sub_wkb___ = factor_scaling * interp_q_to_p_block_0(n_k_p_r,n_w_,n_w_sum,R_k_q_sub_wkb___);
dtemplatedc_wkb__ = reshape(R_k_p_sub_wkb___,[n_w_max*n_k_p_r,n_viewing_azimu_b]);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% %% dtemplatedx_wkb___: %0.2fs',tmp_t)); end;
end;%if flag_d1;
%%%%;
if flag_d2;
tmp_t = tic();
R_k_q_sub_wkb___ = zeros(n_w_max,n_k_p_r,n_viewing_azimu_b);
R_k_q_sub_wkb___(1+[0:l_max_max],:,:) = d2Wdaa_caza_mkb___(1+l_max_max+[0:l_max_max],:,:);
R_k_q_sub_wkb___(1+n_w_max-l_max_max+[0:l_max_max-1],:,:) = d2Wdaa_caza_mkb___(1+[0:l_max_max-1],:,:);
R_k_p_sub_wkb___ = factor_scaling * interp_q_to_p_block_0(n_k_p_r,n_w_,n_w_sum,R_k_q_sub_wkb___);
ddtemplatedaa_wkb__ = reshape(R_k_p_sub_wkb___,[n_w_max*n_k_p_r,n_viewing_azimu_b]);
R_k_q_sub_wkb___ = zeros(n_w_max,n_k_p_r,n_viewing_azimu_b);
R_k_q_sub_wkb___(1+[0:l_max_max],:,:) = d2Wdab_caza_mkb___(1+l_max_max+[0:l_max_max],:,:);
R_k_q_sub_wkb___(1+n_w_max-l_max_max+[0:l_max_max-1],:,:) = d2Wdab_caza_mkb___(1+[0:l_max_max-1],:,:);
R_k_p_sub_wkb___ = factor_scaling * interp_q_to_p_block_0(n_k_p_r,n_w_,n_w_sum,R_k_q_sub_wkb___);
ddtemplatedab_wkb__ = reshape(R_k_p_sub_wkb___,[n_w_max*n_k_p_r,n_viewing_azimu_b]);
R_k_q_sub_wkb___ = zeros(n_w_max,n_k_p_r,n_viewing_azimu_b);
R_k_q_sub_wkb___(1+[0:l_max_max],:,:) = d2Wdac_caza_mkb___(1+l_max_max+[0:l_max_max],:,:);
R_k_q_sub_wkb___(1+n_w_max-l_max_max+[0:l_max_max-1],:,:) = d2Wdac_caza_mkb___(1+[0:l_max_max-1],:,:);
R_k_p_sub_wkb___ = factor_scaling * interp_q_to_p_block_0(n_k_p_r,n_w_,n_w_sum,R_k_q_sub_wkb___);
ddtemplatedac_wkb__ = reshape(R_k_p_sub_wkb___,[n_w_max*n_k_p_r,n_viewing_azimu_b]);
R_k_q_sub_wkb___ = zeros(n_w_max,n_k_p_r,n_viewing_azimu_b);
R_k_q_sub_wkb___(1+[0:l_max_max],:,:) = d2Wdbb_caza_mkb___(1+l_max_max+[0:l_max_max],:,:);
R_k_q_sub_wkb___(1+n_w_max-l_max_max+[0:l_max_max-1],:,:) = d2Wdbb_caza_mkb___(1+[0:l_max_max-1],:,:);
R_k_p_sub_wkb___ = factor_scaling * interp_q_to_p_block_0(n_k_p_r,n_w_,n_w_sum,R_k_q_sub_wkb___);
ddtemplatedbb_wkb__ = reshape(R_k_p_sub_wkb___,[n_w_max*n_k_p_r,n_viewing_azimu_b]);
R_k_q_sub_wkb___ = zeros(n_w_max,n_k_p_r,n_viewing_azimu_b);
R_k_q_sub_wkb___(1+[0:l_max_max],:,:) = d2Wdbc_caza_mkb___(1+l_max_max+[0:l_max_max],:,:);
R_k_q_sub_wkb___(1+n_w_max-l_max_max+[0:l_max_max-1],:,:) = d2Wdbc_caza_mkb___(1+[0:l_max_max-1],:,:);
R_k_p_sub_wkb___ = factor_scaling * interp_q_to_p_block_0(n_k_p_r,n_w_,n_w_sum,R_k_q_sub_wkb___);
ddtemplatedbc_wkb__ = reshape(R_k_p_sub_wkb___,[n_w_max*n_k_p_r,n_viewing_azimu_b]);
R_k_q_sub_wkb___ = zeros(n_w_max,n_k_p_r,n_viewing_azimu_b);
R_k_q_sub_wkb___(1+[0:l_max_max],:,:) = d2Wdcc_caza_mkb___(1+l_max_max+[0:l_max_max],:,:);
R_k_q_sub_wkb___(1+n_w_max-l_max_max+[0:l_max_max-1],:,:) = d2Wdcc_caza_mkb___(1+[0:l_max_max-1],:,:);
R_k_p_sub_wkb___ = factor_scaling * interp_q_to_p_block_0(n_k_p_r,n_w_,n_w_sum,R_k_q_sub_wkb___);
ddtemplatedcc_wkb__ = reshape(R_k_p_sub_wkb___,[n_w_max*n_k_p_r,n_viewing_azimu_b]);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% %% ddtemplatedxy_wkb___: %0.2fs',tmp_t)); end;
end;%if flag_d2;
%%%%%%%%;
flag_check=0;
if flag_check;
tmp_t = tic();
R_k_q_sub_wkb__ = zeros(n_w_sum,n_viewing_azimu_b);
R_k_p_sub_wkb__ = zeros(n_w_sum,n_viewing_azimu_b);
for nviewing_azimu_b=0:n_viewing_azimu_b-1;
R_k_p_wk_ = zeros(n_w_sum,1);
R_k_q_wk_ = zeros(n_w_sum,1);
for nk_p_r=0:n_k_p_r-1;
for m_val=-l_max_max:+l_max_max;
nq = m_val; if (nq<0); nq=nq+n_w_max; end;
R_k_q_wk_(1+nq+nk_p_r*n_w_max) = d0W_caza_mkb___(1+l_max_max+m_val,1+nk_p_r,1+nviewing_azimu_b);
end;%for m_val=-l_max_max:+l_max_max;
end;%for nk_p_r=0:n_k_p_r-1;
R_k_q_sub_wkb__(:,1+nviewing_azimu_b) = R_k_q_wk_;
R_k_p_wk_ = factor_scaling * interp_q_to_p(n_k_p_r,n_w_,n_w_sum,R_k_q_wk_);
R_k_p_sub_wkb__(:,1+nviewing_azimu_b) = R_k_p_wk_;
end;%for nviewing_azimu_b=0:n_viewing_azimu_b-1;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% %% R_k_p_sub_wkb__: %0.2fs',tmp_t)); end;
end;%if flag_check;
%%%%%%%%;

