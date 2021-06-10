function ...
[ ...
 template__ ...
,da_template__ ...
,db_template__ ...
,dc_template__ ...
,dt_template__ ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_k_all_ ...
,n_viewing_all ...
,viewing_azimu_b_all_ ...
,viewing_polar_a_all_ ...
,viewing_weight_all_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,template_k_c_0__ ...
,template_k_c_1__ ...
,template_k_c_2__ ...
] = ...
get_dtemplate_0( ...
 verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_k_p_r_ ...
,l_max_ ...
,a_k_Y_ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,n_w_0in_ ...
);
% Uses spherical-harmonic-expansion a_k_Y_ to evaluate templates on a collection of points on spherical shells determined by k_p_r_. ;
% Also evaluates the gradient of each template with respect to the viewing angles [template_polar_a,template_azimu_b,template_gamma_z]. ;
% Note that this gradient is not perfect: when ca==0 & sa==0 we set the gradient to 0. ;  
% ;
% inputs: ;
% ;
% verbose = integer verbosity_level. ;
% n_k_p_r = integer maximum number of shells. ;
% k_p_r_ = real array of length n_k_p_r; k_p_r_(nk_p_r) = k_p_r_value for shell nk_p_r. ;
% k_p_r_max = real maximum k-value. ;
% weight_k_p_r_ = real array of length n_k_p_r; radial quadrature weights. ;
% l_max_ = integer array of length n_k_p_r; l_max_(nk_p_r) = spherical harmonic order on shell nk_p_r; l_max_(nk_p_r) corresponds to n_lm_(nk_p_r) = (l_max_(nk_p_r)+1)^2 coefficients. ;
% a_k_Y_ = complex array of length \sum_{nk_p_r} (n_lm_(nk_p_r)+1)^2 ; coefficients are ordered in a row, with m varying quickly, l varying slowly and k varying most slowly. ;
% viewing_k_eq_d = real equatorial-distance used for sampling viewing angles and templates. ;
% template_k_eq_d = real equatorial-distance used for sampling viewing angles and templates. ;
% n_w_0in_ = integer array of length n_k_p_r; used if template_k_eq_d <=0; desired n_w_ for templates. ;
% ;
% outputs: ;
% ;
% template__ = complex array of templates for each viewing angle. ;
% da_template__ = complex array of \partial_{template_polar_a} template for each viewing angle. ;
% db_template__ = complex array of \partial_{template_azimu_b} template for each viewing angle. ;
% dc_template__ = complex array of \partial_{template_gamma_z} template for each viewing angle. ;
% dt_template__ = complex array of \partial_{t} template for each viewing angle. ;
% For the dt derivative t refers to the time allowed for surface diffusion on each spherical shell. ;
% ;
% Note that the derivative of each associated legendre function can be calculated using: ;
%{
  l_val=5; 
  cx_ = linspace(-1,1,1025); dcx = mean(diff(cx_)); cxm_ = cx_(2:end-1); sxm_ = sqrt(1-cxm_.^2); 
  Pl__ = legendre(l_val,cx_,'unnorm'); 
  for m_val=0:l_val;
  dPlm_diff_ = (Pl__(1+m_val,3:end) - Pl__(1+m_val,1:end-2))/(2*dcx);
  dPlm_form_ = zeros(size(dPlm_diff_));
  if (1+m_val-1>0);
  dPlm_form_ = dPlm_form_ + 0.5*((l_val+m_val)*(l_val-m_val+1)*Pl__(1+m_val-1,2:end-1))./sxm_;
  end;%if (1+m_val-1>0);
  if (1+m_val+1<=l_val+1);
  dPlm_form_ = dPlm_form_ - 0.5*Pl__(1+m_val+1,2:end-1)./sxm_;
  end;%if (1+m_val+1<=l_val+1);
  if (m_val==0); dPlm_form_ = 2*dPlm_form_; end;
  subplot(2,3,1+m_val); plot(cxm_,dPlm_diff_,'ro',cxm_,dPlm_form_,'k.-');
  xlim([-1,+1]); xlabel('cos(theta)'); title(sprintf('m %d',m_val));
  end;%for m_val=0:l_val;
  figbig;
%}
% ;
% Now for the normalized legendre functions: ;
%{
  l_val=8; 
  cx_ = linspace(-1,1,1025); dcx = mean(diff(cx_)); cxm_ = cx_(2:end-1); sxm_ = sqrt(1-cxm_.^2); 
  tmp_a1 = ((1+2*l_val)/(4*pi));
  m_val_ = 0:+l_val;
  tmp_a2_ = exp(lfactorial(l_val-abs(m_val_)) - lfactorial(l_val+abs(m_val_)));
  tmp_a3_ = sqrt(tmp_a1*tmp_a2_);
  Pl__ = legendre(l_val,cx_,'unnorm');
  Pl_norm__ = diag(tmp_a3_)*legendre(l_val,cx_,'unnorm');
  for m_val=0:l_val;
  dPlm_norm_diff_ = (Pl_norm__(1+m_val,3:end) - Pl_norm__(1+m_val,1:end-2))/(2*dcx);
  dPlm_norm_form_ = zeros(size(dPlm_norm_diff_));
  if (1+m_val-1>0);
  dPlm_norm_form_ = dPlm_norm_form_ + 0.5*((l_val+m_val)*(l_val-m_val+1)*Pl__(1+m_val-1,2:end-1))./sxm_;
  end;%if (1+m_val-1>0);
  if (1+m_val+1<=l_val+1);
  dPlm_norm_form_ = dPlm_norm_form_ - 0.5*Pl__(1+m_val+1,2:end-1)./sxm_;
  end;%if (1+m_val+1<=l_val+1);
  if (m_val==0); dPlm_norm_form_ = 2*dPlm_norm_form_; end;
  dPlm_norm_form_ =  tmp_a3_(1+m_val)*dPlm_norm_form_;
  subplot(3,3,1+m_val); plot(cxm_,dPlm_norm_diff_,'ro',cxm_,dPlm_norm_form_,'k.-');
  xlim([-1,+1]); xlabel('cos(theta)'); title(sprintf('m %d',m_val));
  end;%for m_val=0:l_val;
  figbig;
%}  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

if nargin<8;

rng(1);

%%%%%%%%;
% See get_template_0.m for quadrature test. ;
%%%%%%%%;
%%%%%%%%;
% Now generate a_k_Y_ ; 
%%%%%%%%;
verbose=0;
k_p_r_max = 1.0d0;
n_k_p_r = 2; k_p_r_ = [0.5;1.0]; weight_k_p_r_ = [1;1];
l_max_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
l_max_(1+nk_p_r) = 1+ceil(2*pi*k_p_r_(1+nk_p_r));
end;%for nk_p_r=0:n_k_p_r-1;
n_lm_ = (l_max_+1).^2;
n_lm_max = max(n_lm_);
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
l_max_max = max(l_max_); dWtdkd__l_max_max = 2*l_max_max;
m_max_ = -l_max_max : +l_max_max;
n_m_max = length(m_max_);
disp(sprintf(' %% n_k_p_r %d l_max_max %d n_lm_sum %d',n_k_p_r,l_max_max,n_lm_sum));
Y_l_val_ = zeros(n_lm_sum,1);
Y_m_val_ = zeros(n_lm_sum,1);
for nk_p_r=0:n_k_p_r-1;
l_max = l_max_(1+nk_p_r);
tmp_l_val_ = zeros(n_lm_(1+nk_p_r),1);
tmp_m_val_ = zeros(n_lm_(1+nk_p_r),1);
na=0; 
for l_val=0:l_max;
for m_val=-l_val:+l_val;
tmp_l_val_(1+na) = l_val;
tmp_m_val_(1+na) = m_val;
na=na+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
Y_l_val_(1+tmp_index_) = tmp_l_val_;
Y_m_val_(1+tmp_index_) = tmp_m_val_;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
a_k_Y_ = zeros(n_lm_sum,1);
a_k_Y_ = ( randn(n_lm_sum,1) + i*randn(n_lm_sum,1) ) .* exp(-0.10*Y_l_val_.*(Y_l_val_+1)) .* exp(-0.50*abs(Y_m_val_)) ;

%%%%%%%%;
% Now generate templates. ;
%%%%%%%%;
verbose=1;
k_eq_d = k_p_r_max/32;
template_k_eq_d = k_eq_d;
viewing_k_eq_d = 1*k_eq_d;
[ ...
 template__ ...
,da_template__ ...
,db_template__ ...
,dc_template__ ...
,dt_template__ ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_k_all_ ...
,n_viewing_all ...
,viewing_azimu_b_all_ ...
,viewing_polar_a_all_ ...
,viewing_weight_all_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,template_k_c_0__ ...
,template_k_c_1__ ...
,template_k_c_2__ ...
] = ...
get_dtemplate_0( ...
verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_k_p_r_ ...
,l_max_ ...
,a_k_Y_ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);
flag_plot=0;
if flag_plot;
for nviewing_all=0:6-1;
subplot(2,3,1+nviewing_all);
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(template__(:,1+nviewing_all)),[-1,+1],colormap_beach());
axisnotick; axis image;
title(sprintf('n %d (a%0.2f\\pi,b%0.2f\\pi)',nviewing_all,viewing_polar_a_all_(1+nviewing_all)/pi,viewing_azimu_b_all_(1+nviewing_all)/pi));
end;%for nviewing_all=0:6-1;
end;%if flag_plot;

viewing_k_c_0_ = sin(viewing_polar_a_all_) .* cos(viewing_azimu_b_all_) ;
viewing_k_c_1_ = sin(viewing_polar_a_all_) .* sin(viewing_azimu_b_all_) ;
viewing_k_c_2_ = cos(viewing_polar_a_all_) ;
viewing_k__ = [viewing_k_c_0_,viewing_k_c_1_,viewing_k_c_2_];
[knn_idx_,knn_d_] = knnsearch(viewing_k__,viewing_k__,'K',2); knn_index_ = knn_idx_ - 1;

viewing_polar_a_step_ = zeros(n_viewing_all,1);
viewing_azimu_b_step_ = zeros(n_viewing_all,1);
dtemplate_rel_error_ = zeros(n_viewing_all,1);
dtemplate_abs_error_ = zeros(n_viewing_all,1);
for nviewing_all=0:n_viewing_all-1;
nviewing_ori = nviewing_all; %<-- one of the viewing angles. ;
nviewing_knn = knn_index_(1+nviewing_ori,2); %<-- nearest neighbor. ;
viewing_polar_a_ori = viewing_polar_a_all_(1+nviewing_ori); viewing_azimu_b_ori = viewing_azimu_b_all_(1+nviewing_ori);
viewing_polar_a_knn = viewing_polar_a_all_(1+nviewing_knn); viewing_azimu_b_knn = viewing_azimu_b_all_(1+nviewing_knn);
viewing_polar_a_step = viewing_polar_a_knn - viewing_polar_a_ori;
viewing_azimu_b_step = viewing_azimu_b_knn - viewing_azimu_b_ori;
if (viewing_azimu_b_step <  -pi); viewing_azimu_b_step = viewing_azimu_b_step + 2*pi; end;
if (viewing_azimu_b_step >= +pi); viewing_azimu_b_step = viewing_azimu_b_step - 2*pi; end;
viewing_polar_a_step_(1+nviewing_all) = viewing_polar_a_step;
viewing_azimu_b_step_(1+nviewing_all) = viewing_azimu_b_step;
dtemplate_diff_ = ( template__(:,1+nviewing_knn) - template__(:,1+nviewing_ori) );
da_template_form_ = da_template__(:,1+nviewing_ori)*viewing_polar_a_step;
db_template_form_ = db_template__(:,1+nviewing_ori)*viewing_azimu_b_step;
dtemplate_form_ = da_template_form_ + db_template_form_;
%%%%%%%%;
flag_plot=0;
if flag_plot;
subplot(2,2,1);
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(template__(:,1+nviewing_ori)),[-1,+1],colormap_beach());
axisnotick; axis image; title(sprintf('n %d ori',nviewing_ori));
subplot(2,2,2);
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(template__(:,1+nviewing_knn)),[-1,+1],colormap_beach());
axisnotick; axis image; title(sprintf('n %d knn',nviewing_knn));
subplot(2,2,3);
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(dtemplate_diff_),[-1,+1],colormap_beach());
axisnotick; axis image; title(sprintf('diff'));
subplot(2,2,4);
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(dtemplate_form_),[-1,+1],colormap_beach());
axisnotick; axis image; title(sprintf('form'));
end;%if flag_plot;
%%%%%%%%;
dtemplate_abs_error = fnorm(dtemplate_form_ - dtemplate_diff_);
dtemplate_abs_error_(1+nviewing_all) = dtemplate_abs_error;
dtemplate_rel_error = dtemplate_abs_error/fnorm(dtemplate_form_);
dtemplate_rel_error_(1+nviewing_all) = dtemplate_rel_error;
end;%for nviewing_all=0:n_viewing_all-1;

%%%%%%%%;
% plot errors as function of polar_a and azimu_b. ;
%%%%%%%%;
flag_plot=1;
if flag_plot;
figure(1);clf;
subplot(2,2,[1,2]);
imagesc_polar_a_azimu_b_0(viewing_polar_a_all_,viewing_azimu_b_all_,abs(viewing_polar_a_step_),[0,0.1]);
axis image; axisnotick;
xlim([0,2*pi]); xlabel('azimu_b','Interpreter','none');
ylim([0,1*pi]); xlabel('polar_a','Interpreter','none');
title('viewing_polar_a_step_','Interpreter','none');
subplot(2,2,[3,4]);
imagesc_polar_a_azimu_b_0(viewing_polar_a_all_,viewing_azimu_b_all_,abs(viewing_azimu_b_step_),[0,0.1]);
axis image; axisnotick;
xlim([0,2*pi]); xlabel('azimu_b','Interpreter','none');
ylim([0,1*pi]); xlabel('polar_a','Interpreter','none');
title('viewing_azimu_b_step_','Interpreter','none');
set(gcf,'Position',1+[0,0,1024,1024]);  

figure(2);clf;
subplot(2,2,[3,4]);
imagesc_polar_a_azimu_b_0(viewing_polar_a_all_,viewing_azimu_b_all_,log10(dtemplate_abs_error_),[-3,0],colormap_beach());
axis image; axisnotick;
xlim([0,2*pi]); xlabel('azimu_b','Interpreter','none');
ylim([0,1*pi]); xlabel('polar_a','Interpreter','none');
title('log10(dtemplate_abs_error_)','Interpreter','none');
subplot(2,2,1);
plot(knn_d_(:,2),log10(dtemplate_abs_error_),'.');
xlabel('knn_d','Interpreter','none');  
ylabel('log10(|diff-form|)','Interpreter','none');
title('log10(abs error)');
subplot(2,2,2);
plot(knn_d_(:,2),dtemplate_rel_error_,'.');
xlabel('knn_d','Interpreter','none');  
ylabel('|diff-form|/|form|','Interpreter','none');
title('rel error');
set(gcf,'Position',1+[0,0,1024,1024]);

end;%if flag_plot;

%%%%%%%%;
% check largest outlier ;
%%%%%%%%;
[~,tmp_ij] = max(dtemplate_abs_error_); nviewing_ori = tmp_ij-1;
nviewing_knn = knn_index_(1+nviewing_ori,2); %<-- nearest neighbor. ;
viewing_polar_a_ori = viewing_polar_a_all_(1+nviewing_ori); viewing_azimu_b_ori = viewing_azimu_b_all_(1+nviewing_ori);
viewing_polar_a_knn = viewing_polar_a_all_(1+nviewing_knn); viewing_azimu_b_knn = viewing_azimu_b_all_(1+nviewing_knn);
dtemplate_diff_ = ( template__(:,1+nviewing_knn) - template__(:,1+nviewing_ori) );
da_template_form_ = da_template__(:,1+nviewing_ori)*(viewing_polar_a_knn - viewing_polar_a_ori);
db_template_form_ = db_template__(:,1+nviewing_ori)*(viewing_azimu_b_knn - viewing_azimu_b_ori);
dtemplate_form_ = da_template_form_ + db_template_form_;
%%%%%%%%;
flag_plot=1;
if flag_plot;
figure(3);clf;
subplot(2,2,1);
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(template__(:,1+nviewing_ori)),[-1,+1],colormap_beach());
axisnotick; axis image; title(sprintf('n %d ori',nviewing_ori));
subplot(2,2,2);
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(template__(:,1+nviewing_knn)),[-1,+1],colormap_beach());
axisnotick; axis image; title(sprintf('n %d knn',nviewing_knn));
subplot(2,2,3);
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(dtemplate_diff_),0.005*[-1,+1],colormap_beach());
axisnotick; axis image; title(sprintf('diff'));
subplot(2,2,4);
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(dtemplate_form_),0.005*[-1,+1],colormap_beach());
axisnotick; axis image; title(sprintf('form'));
sgtitle(sprintf('n %d (a%0.2f\\pi,b%0.2f\\pi)',nviewing_ori,viewing_polar_a_all_(1+nviewing_ori)/pi,viewing_azimu_b_all_(1+nviewing_ori)/pi));
end;%if flag_plot;
%%%%%%%%;
dtemplate_abs_error = fnorm(dtemplate_form_ - dtemplate_diff_);
disp(sprintf(' %% dtemplate_abs_error %0.16f',dtemplate_abs_error));
dtemplate_rel_error = dtemplate_abs_error/fnorm(dtemplate_form_);
disp(sprintf(' %% dtemplate_rel_error %0.16f',dtemplate_rel_error));

%%%%%%%%;
% check values on the equator. ;
%%%%%%%%;
[~,tmp_ij] = min(abs(viewing_polar_a_all_ - pi/2)); nviewing_ori = tmp_ij(1+floor(end/2))-1;
nviewing_knn = knn_index_(1+nviewing_ori,2); %<-- nearest neighbor. ;
viewing_polar_a_ori = viewing_polar_a_all_(1+nviewing_ori); viewing_azimu_b_ori = viewing_azimu_b_all_(1+nviewing_ori);
viewing_polar_a_knn = viewing_polar_a_all_(1+nviewing_knn); viewing_azimu_b_knn = viewing_azimu_b_all_(1+nviewing_knn);
dtemplate_diff_ = ( template__(:,1+nviewing_knn) - template__(:,1+nviewing_ori) );
da_template_form_ = da_template__(:,1+nviewing_ori)*(viewing_polar_a_knn - viewing_polar_a_ori);
db_template_form_ = db_template__(:,1+nviewing_ori)*(viewing_azimu_b_knn - viewing_azimu_b_ori);
dtemplate_form_ = da_template_form_ + db_template_form_;
%%%%%%%%;
flag_plot=1;
if flag_plot;
figure(4);clf;
subplot(2,2,1);
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(template__(:,1+nviewing_ori)),[-1,+1],colormap_beach());
axisnotick; axis image; title(sprintf('n %d ori',nviewing_ori));
subplot(2,2,2);
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(template__(:,1+nviewing_knn)),[-1,+1],colormap_beach());
axisnotick; axis image; title(sprintf('n %d knn',nviewing_knn));
subplot(2,2,3);
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(dtemplate_diff_),0.005*[-1,+1],colormap_beach());
axisnotick; axis image; title(sprintf('diff'));
subplot(2,2,4);
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(dtemplate_form_),0.005*[-1,+1],colormap_beach());
axisnotick; axis image; title(sprintf('form'));
sgtitle(sprintf('n %d (a%0.2f\\pi,b%0.2f\\pi)',nviewing_ori,viewing_polar_a_all_(1+nviewing_ori)/pi,viewing_azimu_b_all_(1+nviewing_ori)/pi));
end;%if flag_plot;
%%%%%%%%;
dtemplate_abs_error = fnorm(dtemplate_form_ - dtemplate_diff_);
disp(sprintf(' %% dtemplate_abs_error %0.16f',dtemplate_abs_error));
dtemplate_rel_error = dtemplate_abs_error/fnorm(dtemplate_form_);
disp(sprintf(' %% dtemplate_rel_error %0.16f',dtemplate_rel_error));

disp(sprintf(' %% returning')); return;
end;% if nargin<8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

n_lm_ = (l_max_+1).^2;
n_lm_max = max(n_lm_);
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
l_max_max = max(l_max_);
m_max_ = -l_max_max : +l_max_max;
n_m_max = length(m_max_);
if (verbose); disp(sprintf(' %% n_k_p_r %d l_max_max %d n_lm_sum %d',n_k_p_r,l_max_max,n_lm_sum)); end;
if (isempty(a_k_Y_)); a_k_Y_ = zeros(n_lm_sum,1); end;

%%%%%%%%;
if (verbose); disp(sprintf(' %% First determine the viewing angles on the outermost shell.')); end;
%%%%%%%%;
[ ...
 n_viewing_all ...
,viewing_azimu_b_all_ ...
,viewing_polar_a_all_ ...
,viewing_weight_all_ ...
,~ ...
,~ ...
,~ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
] = ...
sample_shell_5( ...
 k_p_r_max ...
,viewing_k_eq_d ...
,'L' ...
) ; %<-- obtain viewing angles on outer shell. ;
n_viewing_azimu_b_sum = sum(n_viewing_azimu_b_);
n_viewing_azimu_b_csum_ = cumsum([0;n_viewing_azimu_b_]);
if (verbose); disp(sprintf(' %% n_viewing_all %d n_viewing_polar_a %d n_viewing_azimu_b_sum %d',n_viewing_all,n_viewing_polar_a,n_viewing_azimu_b_sum)); end;

%%%%%%%%;
if (verbose); disp(sprintf(' %% Now determine the points along each equatorial plane (i.e., the points for each template).')); end;
% Note that the radial arrangement of these points is determined by the radii of the shells (i..e, the input k_p_r_). ;
%%%%%%%%;
n_w_ = zeros(n_k_p_r,1);
if (template_k_eq_d>0);
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
n_equator = 3+round(2*pi*k_p_r/template_k_eq_d);
n_polar_a = 3+round(n_equator/2);
n_w_(1+nk_p_r) = 2*n_polar_a;
end;%for nk_p_r=0:n_k_p_r-1;
end;%if (template_k_eq_d>0);
if (template_k_eq_d<=0);
n_w_ = n_w_0in_;
assert(numel(n_w_)==n_k_p_r); assert(min(n_w_)>0);
end;%if (template_k_eq_d<=0);
n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);
if (verbose); disp(sprintf(' %% n_w_max %d n_w_sum %d',n_w_max,n_w_sum)); end;
%%%%%%%%;
if (verbose); disp(sprintf(' %% Set up integration weights for the templates.')); end;
%%%%%%%%;
tmp_P_ = zeros(n_k_p_r,n_k_p_r); %<-- polynomials of order 0:n_k_p_r-1 evaluated on k_p_r_/k_p_r_max. ;
tmp_I_ = zeros(n_k_p_r,1); %<-- integrals of those polynomials on the 2d-disc of radius 1. ;
for nk_p_r=0:n_k_p_r-1;
tmp_x = @(x) x.^nk_p_r;
tmp_P_(1+nk_p_r,:) = tmp_x(k_p_r_/k_p_r_max);
tmp_I_(1+nk_p_r) = 2*pi*1/(nk_p_r+2);
end;%for nk_p_r=0:n_k_p_r-1;
tmp_W_ = pinv(tmp_P_,1e-6)*tmp_I_;
if (verbose>1); disp(sprintf(' %% weight error: %0.16f',fnorm(tmp_P_*tmp_W_ - tmp_I_)/fnorm(tmp_I_))); end;
weight_2d_k_p_r_ = tmp_W_*k_p_r_max^2;
weight_2d_k_all_ = zeros(n_w_sum,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
weight_2d_k_all_(1+tmp_index_) = weight_2d_k_p_r_(1+nk_p_r) / max(1,n_w_(1+nk_p_r)) / (2*pi)^2;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
% Set up inner gamma_z for the templates. ;
%%%%%%%%;
gamma_z_all_ = zeros(n_w_sum,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
tmp_val_ = transpose(linspace(0,2*pi,n_w_(1+nk_p_r)+1));
gamma_z_all_(1+tmp_index_) = tmp_val_(1:n_w_(1+nk_p_r));
end;%for nk_p_r=0:n_k_p_r-1;
cc_ = cos(gamma_z_all_);
sc_ = sin(gamma_z_all_);

%%%%%%%%;
% The general formula used here is as follows. ;
% let sa and ca be sin(polar_a) and cos(polar_a), respectively. ;
% let sb and cb be sin(azimu_b) and cos(azimu_b), respectively. ;
% let sc and cc be sin(gamma_z) and cos(gamma_z), respectively. ;
% And rotation by azimu_b about the +z-axis is represented as: ;
% Rz(azimu_b) = ;
% [ +cb -sb 0 ] ;
% [ +sb +cb 0 ] ;
% [  0   0  1 ] ;
% And rotation by polar_a about the +y-axis is represented as: ;
% Ry(polar_a) = ;
% [ +ca 0 +sa ] ;
% [  0  1  0  ] ;
% [ -sa 0 +ca ] ;
% And rotation by gamma_z about the +z-axis is represented as: ;
% Rz(gamma_z) = ;
% [ +cc -sc 0 ] ;
% [ +sc +cc 0 ] ;
% [  0   0  1 ] ;
% Which, collectively, implies that under the transform: ;
% Rz(azimu_b) * Ry(polar_a) * Rz(gamma_z), ;
% Which is the same as: ;
% [ +cb -sb 0 ] [ +ca*cc -ca*sc +sa ]   [ +cb*ca*cc - sb*sc , -cb*ca*sc -sb*cc , +cb*sa ];
% [ +sb +cb 0 ] [ +sc    +cc    0   ] = [ +sb*ca*cc + cb*sc , -sb*ca*sc +cb*cc , +sb*sa ];
% [  0   0  1 ] [ -sa*cc +sa*sc +ca ]   [ -sa*cc            , +sa*sc           , +ca    ];
% the point [1;0;0] is mapped to: ;
% [ template_k_c_0 ; template_k_c_1 ; template_k_c_2 ] = [ +cb*ca*cc - sb*sc ; +sb*ca*cc + cb*sc ; -sa*cc ];
% Note that: ;
% cos(template_polar_a) = -sa*cc ;
% template_azimu_b = atan2( +sb*ca*cc + cb*sc , *cb*ca*cc - sb*sc );
% implying that: ;
% \partial_{viewing_polar_a} cos(template_polar_a) = -ca*cc ;
% \partial_{viewing_azimu_b} cos(template_polar_a) = 0 ;
% \partial_{viewing_gamma_z} cos(template_polar_a) = +sa*sc ;
% and since \partial_{y/x} atan(y/x) = 1/(1+(y/x)^2) = x^2/(x^2+y^2) ;
% and \partial_{z} y(z)/x(z) = (y'x-yx')/x^2, ;
% \partial_{viewing_polar_a} template_azimu_b = (\partial_{viewing_polar_a} y * x - y * \partial_{viewing_polar_a} x ) / (x^2 + y^2) ;
% \partial_{viewing_azimu_b} template_azimu_b = (\partial_{viewing_azimu_b} y * x - y * \partial_{viewing_azimu_b} x ) / (x^2 + y^2) ;
% \partial_{viewing_gamma_z} template_azimu_b = (\partial_{viewing_gamma_z} y * x - y * \partial_{viewing_gamma_z} x ) / (x^2 + y^2) ;
% where:
% \partial_{viewing_polar_a} x = -cb*sa*cc ;
% \partial_{viewing_polar_a} y = -sb*sa*cc ;
% \partial_{viewing_azimu_b} x = -sb*ca*cc - cb*sc ;
% \partial_{viewing_azimu_b} y = +cb*ca*cc - sb*sc ;
% \partial_{viewing_gamma_z} x = -cb*ca*sc - sb*cc ;
% \partial_{viewing_gamma_z} y = -sb*ca*sc + cb*cc ;
%%%%%%%%;

if (verbose); disp(sprintf(' %% template: (%d,%d)=%d (%0.2f GB)',n_w_sum,n_viewing_all,n_w_sum*n_viewing_all,n_w_sum*n_viewing_all*16/1e9)); end;

%%%%%%%%;
if (verbose); disp(sprintf(' %% Now construct array of k_c_?_ values for the templates.')); end;
%%%%%%%%;
template_k_p_r_ = zeros(n_w_sum,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
template_k_p_r_(1+tmp_index_) = k_p_r_(1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
template_k_c_0__ = zeros(n_w_sum,n_viewing_all);
template_k_c_1__ = zeros(n_w_sum,n_viewing_all);
template_k_c_2__ = zeros(n_w_sum,n_viewing_all);
template_khat_c_rr01__ = zeros(n_w_sum,n_viewing_all);
template_khat_c_0__ = zeros(n_w_sum,n_viewing_all);
template_khat_c_1__ = zeros(n_w_sum,n_viewing_all);
template_khat_c_2__ = zeros(n_w_sum,n_viewing_all);
da_template_khat_c_0__ = zeros(n_w_sum,n_viewing_all);
db_template_khat_c_0__ = zeros(n_w_sum,n_viewing_all);
dc_template_khat_c_0__ = zeros(n_w_sum,n_viewing_all);
da_template_khat_c_1__ = zeros(n_w_sum,n_viewing_all);
db_template_khat_c_1__ = zeros(n_w_sum,n_viewing_all);
dc_template_khat_c_1__ = zeros(n_w_sum,n_viewing_all);
%da_template_khat_c_2__ = zeros(n_w_sum,n_viewing_all);
%db_template_khat_c_2__ = zeros(n_w_sum,n_viewing_all);
%dc_template_khat_c_2__ = zeros(n_w_sum,n_viewing_all);
for nviewing_all=0:n_viewing_all-1;
viewing_polar_a = viewing_polar_a_all_(1+nviewing_all); ca = cos(viewing_polar_a); sa = sin(viewing_polar_a);
viewing_azimu_b = viewing_azimu_b_all_(1+nviewing_all); cb = cos(viewing_azimu_b); sb = sin(viewing_azimu_b);
template_k_c_0__(:,1+nviewing_all) = (+cb*ca*cc_ - 1*sb*sc_).*template_k_p_r_;
template_k_c_1__(:,1+nviewing_all) = (+sb*ca*cc_ + 1*cb*sc_).*template_k_p_r_;
template_k_c_2__(:,1+nviewing_all) = (-sa*cc_              ).*template_k_p_r_;
template_khat_c_0__(:,1+nviewing_all) = (+cb*ca*cc_ - 1*sb*sc_);
template_khat_c_1__(:,1+nviewing_all) = (+sb*ca*cc_ + 1*cb*sc_);
template_khat_c_2__(:,1+nviewing_all) = (-sa*cc_              );
template_khat_c_rr01__(:,1+nviewing_all) = (+cb*ca*cc_ - 1*sb*sc_).^2 + (+sb*ca*cc_ + 1*cb*sc_).^2 ;
da_template_khat_c_0__(:,1+nviewing_all) = (-cb*sa*cc_ - 0*sb*sc_);
db_template_khat_c_0__(:,1+nviewing_all) = (-sb*ca*cc_ - 1*cb*sc_);
dc_template_khat_c_0__(:,1+nviewing_all) = (-cb*ca*sc_ - 1*sb*cc_);
da_template_khat_c_1__(:,1+nviewing_all) = (-sb*sa*cc_ + 0*cb*sc_);
db_template_khat_c_1__(:,1+nviewing_all) = (+cb*ca*cc_ - 1*sb*sc_);
dc_template_khat_c_1__(:,1+nviewing_all) = (-sb*ca*sc_ + 1*cb*cc_);
%da_template_khat_c_2__(:,1+nviewing_all) = (-ca*cc_              );
%db_template_khat_c_2__(:,1+nviewing_all) = ( 0                   );
%dc_template_khat_c_2__(:,1+nviewing_all) = (+sa*sc_              );
end;%for nviewing_all=0:n_viewing_all-1;
template_azimu_b__ = atan2(template_khat_c_1__,template_khat_c_0__);
da_template_azimu_b__ = ( da_template_khat_c_1__ .* template_khat_c_0__ - template_khat_c_1__ .* da_template_khat_c_0__ ) ./ template_khat_c_rr01__ ;
db_template_azimu_b__ = ( db_template_khat_c_1__ .* template_khat_c_0__ - template_khat_c_1__ .* db_template_khat_c_0__ ) ./ template_khat_c_rr01__ ;
dc_template_azimu_b__ = ( dc_template_khat_c_1__ .* template_khat_c_0__ - template_khat_c_1__ .* dc_template_khat_c_0__ ) ./ template_khat_c_rr01__ ;
%%%%%%%%;
% Now fix equator. ;
%%%%%%%%;
na=0;
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
for nviewing_all=0:n_viewing_all-1;
fix_index_ = efind(abs(template_khat_c_rr01__(1+tmp_index_,1+nviewing_all))<1e-6); %<-- ca==0 & sc==0;
if numel(fix_index_)>0;
for nl=0:numel(fix_index_)-1;
tmp_index_pre = periodize(fix_index_(1+nl)-1,0,n_w_(1+nk_p_r));
tmp_index_cur = periodize(fix_index_(1+nl)+0,0,n_w_(1+nk_p_r));
tmp_index_pos = periodize(fix_index_(1+nl)+1,0,n_w_(1+nk_p_r));
tmp_pre = da_template_azimu_b__(1+tmp_index_(1+tmp_index_pre),1+nviewing_all);
tmp_cur = da_template_azimu_b__(1+tmp_index_(1+tmp_index_cur),1+nviewing_all);
tmp_pos = da_template_azimu_b__(1+tmp_index_(1+tmp_index_pos),1+nviewing_all);
da_template_azimu_b__(1+tmp_index_(1+tmp_index_cur),1+nviewing_all) = 0.5*(tmp_pre+tmp_pos);
tmp_upd = da_template_azimu_b__(1+tmp_index_(1+tmp_index_cur),1+nviewing_all);
if (verbose>2); disp(sprintf(' %% fixing (%d %d %d) --> (%0.2f,%0.2f-->%0.2f,%0.2f)',tmp_index_pre,tmp_index_cur,tmp_index_pos,tmp_pre,tmp_cur,tmp_upd,tmp_pos)); end;
tmp_pre = db_template_azimu_b__(1+tmp_index_(1+tmp_index_pre),1+nviewing_all);
tmp_pos = db_template_azimu_b__(1+tmp_index_(1+tmp_index_pos),1+nviewing_all);
db_template_azimu_b__(1+tmp_index_(1+tmp_index_cur),1+nviewing_all) = 0.5*(tmp_pre+tmp_pos);
tmp_pre = dc_template_azimu_b__(1+tmp_index_(1+tmp_index_pre),1+nviewing_all);
tmp_pos = dc_template_azimu_b__(1+tmp_index_(1+tmp_index_pos),1+nviewing_all);
dc_template_azimu_b__(1+tmp_index_(1+tmp_index_cur),1+nviewing_all) = 0.5*(tmp_pre+tmp_pos);
na=na+1;
end;%for nl=0:numel(fix_index_)-1;
end;%if numel(fix_index_)>0;
end;%for nviewing_all=0:n_viewing_all-1;
end;%for nk_p_r=0:n_k_p_r-1;
if (verbose); disp(sprintf(' %% d0_template_azimu_b_: fixed %d points',na)); end;
if (verbose>1); disp(sprintf(' %% ~isfinite(da_template_azimu_b__ %d',numel(find(~isfinite(da_template_azimu_b__))))); end;
if (verbose>1); disp(sprintf(' %% ~isfinite(db_template_azimu_b__ %d',numel(find(~isfinite(db_template_azimu_b__))))); end;
if (verbose>1); disp(sprintf(' %% ~isfinite(dc_template_azimu_b__ %d',numel(find(~isfinite(dc_template_azimu_b__))))); end;

expi_template_azimu_b__ = exp(+i*template_azimu_b__);
clear template_azimu_b__;
if nargout<11; clear template_k_c_0__ template_k_c_1__ template_k_c_2__ ; end;
%%%%%%%%;
% We also use a condensed array, called condense_template_khat_c_2__, which only depends on the polar_a, and not on azimu_b. ;
%%%%%%%%;
condense_template_khat_c_2__ = zeros(n_w_sum,n_viewing_polar_a);
condense_template_cos_polar_a__ = zeros(n_w_sum,n_viewing_polar_a);
condense_template_sin_polar_a__ = zeros(n_w_sum,n_viewing_polar_a);
condense_da_template_cos_polar_a__ = zeros(n_w_sum,n_viewing_polar_a);
condense_db_template_cos_polar_a__ = zeros(n_w_sum,n_viewing_polar_a);
condense_dc_template_cos_polar_a__ = zeros(n_w_sum,n_viewing_polar_a);
for nviewing_polar_a=0:n_viewing_polar_a-1;
viewing_polar_a = viewing_polar_a_(1+nviewing_polar_a); ca = cos(viewing_polar_a); sa = sin(viewing_polar_a);
condense_template_khat_c_2__(:,1+nviewing_polar_a) = -sa*cc_ ;
condense_template_cos_polar_a__(:,1+nviewing_polar_a) = -sa*cc_ ;
condense_template_sin_polar_a__(:,1+nviewing_polar_a) = sqrt( 1 - condense_template_cos_polar_a__(:,1+nviewing_polar_a).^2 );
condense_da_template_cos_polar_a__(:,1+nviewing_polar_a) = -ca*cc_ ;
condense_db_template_cos_polar_a__(:,1+nviewing_polar_a) = 0       ;
condense_dc_template_cos_polar_a__(:,1+nviewing_polar_a) = +sa*sc_ ;
end;%for nviewing_polar_a=0:n_viewing_all-1;

%%%%%%%%;
if (verbose); disp(sprintf(' %% Now evaluate associated legendre polynomials at the varous k_c_2 values.')); end;
% Here legendre_evaluate_{1+nk_p_r}{1+l_val}(1+l_val+m_val,1+nw,1+nviewing_polar_a) contains ;
% the associated legendre-function of degree l_val and order abs(m_val) (ranging from 0 to +l_val) ;
% evaluated at the k_c_2 value stored in condense_template_khat_c_2__(1+n_w_csum_(1+nk_p_r)+nw,1+nviewing_polar_a). ;
% Note that this is associated with ring/shell nk_p_r and viewing_polar_a_(1+nviewing_polar_a). ;
% The legendre_normalization_{1+nk_p_r}{1+l_val}(1+abs(m_val)) contains ;
% The normalization coefficient for the spherical harmonics associated with l_val and m_val. ;
% Note that this is somewhat redundant (as it does not depend explicitly on the shell). ;
%%%%%%%%;
legendre_evaluate_ = cell(n_k_p_r,1);
legendre_normalization_ = cell(n_k_p_r,1);
na=0;
for nk_p_r=0:n_k_p_r-1;
l_max = l_max_(1+nk_p_r);
tmp_index_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
legendre_evaluate_{1+nk_p_r} = cell(l_max+1,1);
dlegendre_evaluate_{1+nk_p_r} = cell(l_max+1,1);
for l_val=0:l_max;
tmp_P__ = zeros(1+1*l_val,n_w_(1+nk_p_r),n_viewing_polar_a);
tmp_a1 = ((1+2*l_val)/(4*pi));
m_val_ = -l_val:+l_val;
tmp_a2_ = exp(lfactorial(l_val-abs(m_val_)) - lfactorial(l_val+abs(m_val_)));
tmp_a3_ = sqrt(tmp_a1*tmp_a2_);
tmp_template_cos_polar_a__ = condense_template_cos_polar_a__(1+tmp_index_,:);
tmp_template_sin_polar_a__ = condense_template_sin_polar_a__(1+tmp_index_,:);
tmp_denominator__ = reshape(tmp_template_sin_polar_a__,[1,n_w_(1+nk_p_r),n_viewing_polar_a]);
fix_index_ = efind(abs(tmp_denominator__)<1e-6);
tmp_P__ = legendre(l_val,tmp_template_cos_polar_a__,'unnorm');
tmp_t = tic;
legendre_evaluate_{1+nk_p_r}{1+l_val} = reshape(tmp_P__,[1+1*l_val,n_w_(1+nk_p_r),n_viewing_polar_a]);
tmp_dP___ = zeros([1+l_val,n_w_(1+nk_p_r),n_viewing_polar_a]);
for m_val=0:l_val;
if (1+m_val-1>0);
tmp_dP___(1+m_val,:,:) = tmp_dP___(1+m_val,:,:) + 0.5*((l_val+m_val)*(l_val-m_val+1)*legendre_evaluate_{1+nk_p_r}{1+l_val}(1+m_val-1,:,:)) ./ tmp_denominator__;
end;%if (1+m_val-1>0);
if (1+m_val+1<=l_val+1);
tmp_dP___(1+m_val,:,:) = tmp_dP___(1+m_val,:,:) - 0.5*legendre_evaluate_{1+nk_p_r}{1+l_val}(1+m_val+1,:,:) ./ tmp_denominator__;
end;%if (1+m_val+1<=l_val+1);
if (m_val==0); tmp_dP___(1+m_val,:,:) = 2*tmp_dP___(1+m_val,:,:); end;
end;%for m_val=0:l_val;
if numel(fix_index_)>0;
[fix_index_nw_,fix_index_na_] = eind2esub([n_w_(1+nk_p_r),n_viewing_polar_a],fix_index_);
for nl=0:numel(fix_index_)-1;
tmp_index_pre = periodize(fix_index_nw_(1+nl)-1,0,n_w_(1+nk_p_r));
tmp_index_cur = periodize(fix_index_nw_(1+nl)+0,0,n_w_(1+nk_p_r));
tmp_index_pos = periodize(fix_index_nw_(1+nl)+1,0,n_w_(1+nk_p_r));
for m_val=0:l_val;
tmp_dP_pre = tmp_dP___(1+m_val,1+tmp_index_pre,1+fix_index_na_(1+nl));
tmp_dP_cur = tmp_dP___(1+m_val,1+tmp_index_cur,1+fix_index_na_(1+nl));
tmp_dP_pos = tmp_dP___(1+m_val,1+tmp_index_pos,1+fix_index_na_(1+nl));
tmp_dP___(1+m_val,1+fix_index_nw_(1+nl),1+fix_index_na_(1+nl)) = 0.5*(tmp_dP_pre+tmp_dP_pos);
tmp_dP_upd = tmp_dP___(1+m_val,1+tmp_index_cur,1+fix_index_na_(1+nl));
if (verbose>2); disp(sprintf(' %% fixing (%d,%d,%d) --> (%0.2f,%0.2f-->%0.2f,%0.2f)',tmp_index_pre,tmp_index_cur,tmp_index_pos,tmp_dP_pre,tmp_dP_cur,tmp_dP_upd,tmp_dP_pos)); end;
na=na+1;
end;%for m_val=0:l_val;
end;%for nl=0:numel(fix_index_)-1;
end;%if numel(fix_index_)>0;
dlegendre_evaluate_{1+nk_p_r}{1+l_val} = tmp_dP___;
tmp_t = toc(tmp_t);
if (verbose>1); disp(sprintf(' %% nk_p_r %d/%d l_val %d/%d legendre_evaluate(%d,%d) %0.2fs',nk_p_r,n_k_p_r,l_val,l_max,n_w_(1+nk_p_r),n_viewing_polar_a,tmp_t)); end;
legendre_normalization_{1+nk_p_r}{1+l_val} = tmp_a3_;
if (verbose>1); disp(sprintf(' %% ~isfinite(legendre_evaluate_{1+%d}{1+%d} %d',nk_p_r,l_val,numel(find(~isfinite(legendre_evaluate_{1+nk_p_r}{1+l_val}))))); end;
if (verbose>1); disp(sprintf(' %% ~isfinite(dlegendre_evaluate_{1+%d}{1+%d} %d',nk_p_r,l_val,numel(find(~isfinite(dlegendre_evaluate_{1+nk_p_r}{1+l_val}))))); end;
end;%for l_val=0:l_max;
end;%for nk_p_r=0:n_k_p_r-1;
if (verbose); disp(sprintf(' %% legendre_evaluate_: fixed %d points',na)); end;

%%%%%%%%;
if (verbose); disp(sprintf(' %% Now accumulate the legendre_evaluates over l_val, for each m_val.')); end;
% We account for the normalization coefficients here, ;
% so that later we can apply the complex exponential to produce the spherical-harmonics. ;
% More specifically: ;
% spherical_harmonic_unphased_{1+nk_p_r}(1+l_max_(1+nk_p_r)+m_val,1+nw,1+nviewing_polar_a) ;
% contains the sum: ;
% \sum_{l_val=0}^{l_max_(1+nk_p_r)} ... ;
%             legendre_normalization_{1+nk_p_r}{1+l_val}(1+l_val+m_val) ... ;
%           * legendre_evaluate_{1+nk_p_r}{1+l_val}(1+l_val+m_val,1+nw,1+nviewing_polar_a) ... ;
%           * a_k_Y_(1+n_lm_csum_(1+nk_p_r)+(1+l_val-1)^2+l_val+m_val). ;
% Note that this is 'unphased', in the sense that the full evaluate requires: ;
% spherical_harmonic_evaluate__(1+n_w_csum_(1+nk_p_r)+nw,1+nviewing_all) = ... ;
% \sum_{m_val=-l_max_(1+nk_p_r)}^{+l_max_(1+nk_p_r)} ... ;
%             spherical_harmonic_unphased_{1+nk_p_r}(1+l_max_(1+nk_p_r)+m_val,1+nw,1+nviewing_polar_a) ... ;
%           * exp(+i*m_val*template_azimu_b__(1+n_w_csum_(1+nk_p_r)+nw,1+nviewing_all)). ;
% Note that this final exponential can be calculated as: ;
% expi_template_azimu_b__(1+n_w_csum_(1+nk_p_r)+nw,1+nviewing_all)^m_val. ;
% ;
% In a similar fashion we calculate: ;
% spherical_harmonic_da_unphased_{1+nk_p_r}(1+l_max_(1+nk_p_r)+m_val,1+nw,1+nviewing_polar_a) ;
% contains the sum: ;
% \sum_{l_val=0}^{l_max_(1+nk_p_r)} ... ;
%             legendre_normalization_{1+nk_p_r}{1+l_val}(1+l_val+m_val) ... ;
%           * dlegendre_evaluate_{1+nk_p_r}{1+l_val}(1+l_val+m_val,1+nw,1+nviewing_polar_a) ... ;
%           * condense_da_template_cos_polar_a__(1+nw,1+nviewing_polar_a) ... ;
%           * a_k_Y_(1+n_lm_csum_(1+nk_p_r)+(1+l_val-1)^2+l_val+m_val). ;
% Note that this is 'unphased', in the sense that the full evaluate requires: ;
% spherical_harmonic_da_evaluate__(1+n_w_csum_(1+nk_p_r)+nw,1+nviewing_all) = ... ;
% \sum_{m_val=-l_max_(1+nk_p_r)}^{+l_max_(1+nk_p_r)} ... ;
%             spherical_harmonic_da_unphased_{1+nk_p_r}(1+l_max_(1+nk_p_r)+m_val,1+nw,1+nviewing_polar_a) ... ;
%           * exp(+i*m_val*template_azimu_b__(1+n_w_csum_(1+nk_p_r)+nw,1+nviewing_all)) ... ;
% + ... ;
% \sum_{m_val=-l_max_(1+nk_p_r)}^{+l_max_(1+nk_p_r)} ... ;
%             spherical_harmonic_unphased_{1+nk_p_r}(1+l_max_(1+nk_p_r)+m_val,1+nw,1+nviewing_polar_a) ... ;
%           * +i*m_val*exp(+i*m_val*template_azimu_b__(1+n_w_csum_(1+nk_p_r)+nw,1+nviewing_all)) ... ;
%           * da_template_azimu_b__(1+n_w_csum_(1+nk_p_r)+nw,1+nviewing_all). ;
% Note that this final exponential can be calculated as: ;
% +i*m_val*expi_template_azimu_b__(1+n_w_csum_(1+nk_p_r)+nw,1+nviewing_all)^m_val. ;
% ;
% Similar formulae apply for db and dc. ;
%%%%%%%%;
spherical_harmonic_unphased_ = cell(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);  
spherical_harmonic_unphased_{1+nk_p_r} = zeros(1+2*l_max_(1+nk_p_r),n_w_(1+nk_p_r),n_viewing_polar_a);
spherical_harmonic_da_unphased_{1+nk_p_r} = zeros(1+2*l_max_(1+nk_p_r),n_w_(1+nk_p_r),n_viewing_polar_a);
spherical_harmonic_db_unphased_{1+nk_p_r} = zeros(1+2*l_max_(1+nk_p_r),n_w_(1+nk_p_r),n_viewing_polar_a);
spherical_harmonic_dc_unphased_{1+nk_p_r} = zeros(1+2*l_max_(1+nk_p_r),n_w_(1+nk_p_r),n_viewing_polar_a);
spherical_harmonic_dt_unphased_{1+nk_p_r} = zeros(1+2*l_max_(1+nk_p_r),n_w_(1+nk_p_r),n_viewing_polar_a);
l_max = l_max_(1+nk_p_r);
for m_val=-l_max:l_max;
for l_val=abs(m_val):l_max;
spherical_harmonic_unphased_{1+nk_p_r}(1+l_max_(1+nk_p_r)+m_val,:,:) ...
= spherical_harmonic_unphased_{1+nk_p_r}(1+l_max_(1+nk_p_r)+m_val,:,:) ...
  + legendre_normalization_{1+nk_p_r}{1+l_val}(1+l_val+m_val) ...
    * legendre_evaluate_{1+nk_p_r}{1+l_val}(1+abs(m_val),:,:) ...
    * a_k_Y_(1+n_lm_csum_(1+nk_p_r)+(1+l_val-1)^2+l_val+m_val) ...
;
spherical_harmonic_da_unphased_{1+nk_p_r}(1+l_max_(1+nk_p_r)+m_val,:,:) ...
= spherical_harmonic_da_unphased_{1+nk_p_r}(1+l_max_(1+nk_p_r)+m_val,:,:) ...
  + legendre_normalization_{1+nk_p_r}{1+l_val}(1+l_val+m_val) ...
    * dlegendre_evaluate_{1+nk_p_r}{1+l_val}(1+abs(m_val),:,:) ...
    .* reshape(condense_da_template_cos_polar_a__(1+tmp_index_,:),[1,n_w_(1+nk_p_r),n_viewing_polar_a]) ...
    * a_k_Y_(1+n_lm_csum_(1+nk_p_r)+(1+l_val-1)^2+l_val+m_val) ...
;
spherical_harmonic_db_unphased_{1+nk_p_r}(1+l_max_(1+nk_p_r)+m_val,:,:) ...
= spherical_harmonic_db_unphased_{1+nk_p_r}(1+l_max_(1+nk_p_r)+m_val,:,:) ...
  + legendre_normalization_{1+nk_p_r}{1+l_val}(1+l_val+m_val) ...
    * dlegendre_evaluate_{1+nk_p_r}{1+l_val}(1+abs(m_val),:,:) ...
    .* reshape(condense_db_template_cos_polar_a__(1+tmp_index_,:),[1,n_w_(1+nk_p_r),n_viewing_polar_a]) ...
    * a_k_Y_(1+n_lm_csum_(1+nk_p_r)+(1+l_val-1)^2+l_val+m_val) ...
;
spherical_harmonic_dc_unphased_{1+nk_p_r}(1+l_max_(1+nk_p_r)+m_val,:,:) ...
= spherical_harmonic_dc_unphased_{1+nk_p_r}(1+l_max_(1+nk_p_r)+m_val,:,:) ...
  + legendre_normalization_{1+nk_p_r}{1+l_val}(1+l_val+m_val) ...
    * dlegendre_evaluate_{1+nk_p_r}{1+l_val}(1+abs(m_val),:,:) ...
    .* reshape(condense_dc_template_cos_polar_a__(1+tmp_index_,:),[1,n_w_(1+nk_p_r),n_viewing_polar_a]) ...
    * a_k_Y_(1+n_lm_csum_(1+nk_p_r)+(1+l_val-1)^2+l_val+m_val) ...
;
spherical_harmonic_dt_unphased_{1+nk_p_r}(1+l_max_(1+nk_p_r)+m_val,:,:) ...
= spherical_harmonic_dt_unphased_{1+nk_p_r}(1+l_max_(1+nk_p_r)+m_val,:,:) ...
  + legendre_normalization_{1+nk_p_r}{1+l_val}(1+l_val+m_val) ...
    * legendre_evaluate_{1+nk_p_r}{1+l_val}(1+abs(m_val),:,:) ...
    * a_k_Y_(1+n_lm_csum_(1+nk_p_r)+(1+l_val-1)^2+l_val+m_val) ...
    * -l_val*(l_val+1) ...
;
end;%for l_val=abs(m_val):l_max;
end;%for m_val=-l_max:l_max;
end;%for nk_p_r=0:n_k_p_r-1;

%%%%%%%%;
if (verbose); disp(sprintf(' %% now perform the final sum over m_val.')); end;
%%%%%%%%;
spherical_harmonic_evaluate__ = zeros(n_w_sum,n_viewing_all);
spherical_harmonic_da_evaluate__ = zeros(n_w_sum,n_viewing_all);
spherical_harmonic_db_evaluate__ = zeros(n_w_sum,n_viewing_all);
spherical_harmonic_dc_evaluate__ = zeros(n_w_sum,n_viewing_all);
spherical_harmonic_dt_evaluate__ = zeros(n_w_sum,n_viewing_all);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);  
tmp_spherical_harmonic_evaluate__ = zeros(n_w_(1+nk_p_r),n_viewing_all);
tmp_spherical_harmonic_da_evaluate__ = zeros(n_w_(1+nk_p_r),n_viewing_all);
tmp_spherical_harmonic_db_evaluate__ = zeros(n_w_(1+nk_p_r),n_viewing_all);
tmp_spherical_harmonic_dc_evaluate__ = zeros(n_w_(1+nk_p_r),n_viewing_all);
tmp_spherical_harmonic_dt_evaluate__ = zeros(n_w_(1+nk_p_r),n_viewing_all);
l_max = l_max_(1+nk_p_r);
nviewing_all=0;
for nviewing_polar_a=0:n_viewing_polar_a-1;
n_viewing_azimu_b = n_viewing_azimu_b_(1+nviewing_polar_a);
%viewing_polar_a = viewing_polar_a_(1+nviewing_polar_a);
for nviewing_azimu_b=0:n_viewing_azimu_b-1;
%viewing_azimu_b = 2*pi*nviewing_azimu_b/max(1,n_viewing_azimu_b);
tmp_expi_sub = expi_template_azimu_b__(1+tmp_index_,1+nviewing_all);
tmp_expi_pre = tmp_expi_sub.^(-l_max);
tmp_expi_pos = tmp_expi_pre;
tmp_sum_ = zeros(n_w_(1+nk_p_r),1);
tmp_sum_da_ = zeros(n_w_(1+nk_p_r),1);
tmp_sum_db_ = zeros(n_w_(1+nk_p_r),1);
tmp_sum_dc_ = zeros(n_w_(1+nk_p_r),1);
tmp_sum_dt_ = zeros(n_w_(1+nk_p_r),1);
for m_val=-l_max:+l_max;
tmp_sum_ = tmp_sum_ + transpose(spherical_harmonic_unphased_{1+nk_p_r}(1+l_max+m_val,:,1+nviewing_polar_a)).*tmp_expi_pos;
tmp_sum_da_ = ...
tmp_sum_da_ ...
+ transpose(spherical_harmonic_da_unphased_{1+nk_p_r}(1+l_max+m_val,:,1+nviewing_polar_a)).*tmp_expi_pos ...
+ transpose(spherical_harmonic_unphased_{1+nk_p_r}(1+l_max+m_val,:,1+nviewing_polar_a)) ...
  .* da_template_azimu_b__(1+tmp_index_,1+nviewing_all) ...
  .* +i .* m_val .* tmp_expi_pos ...
;
tmp_sum_db_ = ...
tmp_sum_db_ ...
+ transpose(spherical_harmonic_db_unphased_{1+nk_p_r}(1+l_max+m_val,:,1+nviewing_polar_a)).*tmp_expi_pos ...
+ transpose(spherical_harmonic_unphased_{1+nk_p_r}(1+l_max+m_val,:,1+nviewing_polar_a)) ...
  .* db_template_azimu_b__(1+tmp_index_,1+nviewing_all) ...
  .* +i .* m_val .* tmp_expi_pos ...
;
tmp_sum_dc_ = ...
tmp_sum_dc_ ...
+ transpose(spherical_harmonic_dc_unphased_{1+nk_p_r}(1+l_max+m_val,:,1+nviewing_polar_a)).*tmp_expi_pos ...
+ transpose(spherical_harmonic_unphased_{1+nk_p_r}(1+l_max+m_val,:,1+nviewing_polar_a)) ...
  .* dc_template_azimu_b__(1+tmp_index_,1+nviewing_all) ...
  .* +i .* m_val .* tmp_expi_pos ...
;
tmp_sum_dt_ = tmp_sum_dt_ + transpose(spherical_harmonic_dt_unphased_{1+nk_p_r}(1+l_max+m_val,:,1+nviewing_polar_a)).*tmp_expi_pos;
tmp_expi_pos = tmp_expi_pos.*tmp_expi_sub;
end;%for m_val=-l_max:+l_max;
%%%%%%%%;
index_cc0 = 0;
if (1-abs(template_khat_c_2__(1+n_w_csum_(1+nk_p_r)+index_cc0,1+nviewing_all))<1e-6); %<-- interpolate. ;
index_pre = periodize(index_cc0-1,0,n_w_(1+nk_p_r));
index_pos = periodize(index_cc0+1,0,n_w_(1+nk_p_r));
tmp_sum_da_(1+index_cc0) = 0.5*(tmp_sum_da_(1+index_pre) + tmp_sum_da_(1+index_pos));
tmp_sum_db_(1+index_cc0) = 0.5*(tmp_sum_db_(1+index_pre) + tmp_sum_db_(1+index_pos));
tmp_sum_dc_(1+index_cc0) = 0.5*(tmp_sum_dc_(1+index_pre) + tmp_sum_dc_(1+index_pos));
end;%if (1-abs(template_khat_c_2__(1+n_w_csum_(1+nk_p_r)+0,1+nviewing_all))<1e-6);
index_cc0 = round(n_w_(1+nk_p_r)/2);
if (1-abs(template_khat_c_2__(1+n_w_csum_(1+nk_p_r)+index_cc0,1+nviewing_all))<1e-6); %<-- interpolate. ;
index_pre = periodize(index_cc0-1,0,n_w_(1+nk_p_r));
index_pos = periodize(index_cc0+1,0,n_w_(1+nk_p_r));
tmp_sum_da_(1+index_cc0) = 0.5*(tmp_sum_da_(1+index_pre) + tmp_sum_da_(1+index_pos));
tmp_sum_db_(1+index_cc0) = 0.5*(tmp_sum_db_(1+index_pre) + tmp_sum_db_(1+index_pos));
tmp_sum_dc_(1+index_cc0) = 0.5*(tmp_sum_dc_(1+index_pre) + tmp_sum_dc_(1+index_pos));
end;%if (1-abs(template_khat_c_2__(1+n_w_csum_(1+nk_p_r)+0,1+nviewing_all))<1e-6);
%%%%%%%%;
tmp_spherical_harmonic_evaluate__(:,1+nviewing_all) = tmp_sum_;
tmp_spherical_harmonic_da_evaluate__(:,1+nviewing_all) = tmp_sum_da_;
tmp_spherical_harmonic_db_evaluate__(:,1+nviewing_all) = tmp_sum_db_;
tmp_spherical_harmonic_dc_evaluate__(:,1+nviewing_all) = tmp_sum_dc_;
tmp_spherical_harmonic_dt_evaluate__(:,1+nviewing_all) = tmp_sum_dt_;
nviewing_all=nviewing_all+1;
end;%for nviewing_azimu_b=0:n_viewing_azimu_b-1;
end;%for nviewing_polar_a=0:n_viewing_polar_a-1;
spherical_harmonic_evaluate__(1+tmp_index_,:) = tmp_spherical_harmonic_evaluate__;
spherical_harmonic_da_evaluate__(1+tmp_index_,:) = tmp_spherical_harmonic_da_evaluate__;
spherical_harmonic_db_evaluate__(1+tmp_index_,:) = tmp_spherical_harmonic_db_evaluate__;
spherical_harmonic_dc_evaluate__(1+tmp_index_,:) = tmp_spherical_harmonic_dc_evaluate__;
spherical_harmonic_dt_evaluate__(1+tmp_index_,:) = tmp_spherical_harmonic_dt_evaluate__;
end;%for nk_p_r=0:n_k_p_r-1;

template__ = spherical_harmonic_evaluate__;
da_template__ = spherical_harmonic_da_evaluate__;
db_template__ = spherical_harmonic_db_evaluate__;
dc_template__ = spherical_harmonic_dc_evaluate__;
dt_template__ = spherical_harmonic_dt_evaluate__;
