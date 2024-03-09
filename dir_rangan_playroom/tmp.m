%%%%%%%%;
% Now working on Feng, Wang, Yang, Jin, 2015. ;
%%%%%%%%;

l_max = 3; n_l = l_max;
m_max_ = transpose(-l_max:+l_max);
n_m_max = numel(m_max_);
beta = pi/6;
d0V_ = wignerd_b(l_max,beta);
tmp_V__ = d0V_{1+l_max};

l_val = l_max;
m_val_ = transpose([-l_val:+l_val]);
n_m_val = numel(m_val_);
X_ = sqrt( (l_val + m_val_) .* (1 + l_val - m_val_) ); 
D__ = spdiags([-X_ , +flip(X_)],[+1,-1],n_m_val,n_m_val);
[V__,D__] = eigs(D__,n_m_val);
%%%%%%%%;
% see tmp10.m ;
%%%%%%%%;


disp('returning'); return;
 
%%%%%%%%;
% setting up sparse link operator. ;
% +[d^{m0}_{l}H^{m0+1,m1}_{l} - d^{m0-1}_{l}H^{m0-1,m1}_{l}] - [d^{m1}_{l}H^{m0,m1+1}_{l} - d^{m1-1}_{l}H^{m0,m1-1}_{l}] = 0 ;
% bsxfun(@times,H__(:),sparse_link__) \cdot d_ = 0_ ;
% Verdict: looks as though it is full rank? Unsure how to fix. ;
%%%%%%%%;
d0U_ = cell(1+l_max,1);
for l_val=0:l_max-1;
d0U_{1+l_val} = zeros(1+2*l_val,1+2*l_val);
for m_val=-l_val:+l_val;
d0U_{1+l_val}(1+l_val+m_val,:) = d0V_{1+l_val}(1+l_val+m_val,:)*sqrt(2*(l_val+abs(m_val)) + 1);
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max-1;  
%%%%;
l_val = 10; tmp_n_m = 1+2*l_val;
nl=0;
index_H0_=[];
index_H1_=[];
index_d_ = [];
sgn_ = [];
for m0_val=-l_val:+l_val;
for m1_val=-l_val:+l_val;
%%%%;
if (abs(m0_val+1)<=l_val); 
index_H0_(1+nl) = l_val+(m0_val+1);
index_H1_(1+nl) = l_val+(m1_val+0);
index_d_(1+nl) = l_val+(m0_val+0);
sgn_(1+nl) = +1;
nl=nl+1;
end;%if (abs(m0_val+1)<=l_val); 
%%%%;
if (abs(m0_val-1)<=l_val); 
index_H0_(1+nl) = l_val+(m0_val-1);
index_H1_(1+nl) = l_val+(m1_val+0);
index_d_(1+nl) = l_val+(m0_val-1);
sgn_(1+nl) = -1;
nl=nl+1;
end;%if (abs(m0_val+1)<=l_val); 
%%%%;
if (abs(m1_val+1)<=l_val); 
index_H0_(1+nl) = l_val+(m0_val+0);
index_H1_(1+nl) = l_val+(m1_val+1);
index_d_(1+nl) = l_val+(m1_val+0);
sgn_(1+nl) = -1;
nl=nl+1;
end;%if (abs(m1_val+1)<=l_val); 
%%%%;
if (abs(m1_val-1)<=l_val); 
index_H0_(1+nl) = l_val+(m0_val+0);
index_H1_(1+nl) = l_val+(m1_val-1);
index_d_(1+nl) = l_val+(m1_val-1);
sgn_(1+nl) = +1;
nl=nl+1;
end;%if (abs(m1_val-1)<=l_val); 
%%%%;
end;%for m1_val=-l_val:+l_val;
end;%for m0_val=-l_val:+l_val;
index_H_ = index_H0_ + index_H1_*tmp_n_m;
sparse_link__ = sparse(1+index_H_,1+index_d_,sgn_,tmp_n_m^2,tmp_n_m);
%%%%;
tmp_U__ = d0U_{1+l_val};
D__ = bsxfun(@times,reshape(tmp_U__,[tmp_n_m^2,1]),sparse_link__);
rank(full(D__)),;


%%%%%%%%;
% Trying to determine recursion in Guraimov and Duraiswami 2014.
% See wignerd_c.m. ;
% d^{m0}_{l}H^{m0+1,m1}_{l} - d^{m0-1}_{l}H^{m0-1,m1}_{l} = d^{m1}_{l}H^{m0,m1+1}_{l} - d^{m1-1}_{l}H^{m0,m1-1}_{l} ;
%%%%%%%%;
d0U_ = cell(1+l_max,1);
for l_val=0:l_max-1;
d0U_{1+l_val} = zeros(1+2*l_val,1+2*l_val);
for m_val=-l_val:+l_val;
d0U_{1+l_val}(1+l_val+m_val,:) = d0V_{1+l_val}(1+l_val+m_val,:)*sqrt(2*(l_val+abs(m_val)) + 1);
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max-1;  
%%%%%%%%;
nf=0;
l_val = 10; tmp_n_m = 1+2*l_val;
tmp_U__ = d0U_{1+l_val};
tmp_d_ = d_lm__(1+l_val,1+l_max+[-l_val:+l_val]);
tmp_D0__ = ...
+bsxfun(@times,reshape(circshift(tmp_d_,+1),[tmp_n_m,1]),circshift(tmp_U__,+1,1)) ...
-bsxfun(@times,reshape(circshift(tmp_d_,-0),[tmp_n_m,1]),circshift(tmp_U__,-1,1)) ...
;
tmp_D1__ = ...
+bsxfun(@times,reshape(circshift(tmp_d_,+1),[1,tmp_n_m]),circshift(tmp_U__,+1,2)) ...
-bsxfun(@times,reshape(circshift(tmp_d_,-0),[1,tmp_n_m]),circshift(tmp_U__,-1,2)) ...
;
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,3,1); imagesc(tmp_D0__,[-1,1]); axis image; axisnotick; fig80s; title('tmp_D0__','Interpreter','none'); colorbar;
subplot(1,3,2); imagesc(tmp_D1__,[-1,1]); axis image; axisnotick; fig80s; title('tmp_D1__','Interpreter','none'); colorbar;
subplot(1,3,3); imagesc(abs(tmp_D1__)-abs(tmp_D0__),[-1,1]); axis image; axisnotick; fig80s; title('diff','Interpreter','none'); colorbar;
%%%%%%%%;
disp('returning');return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% Now we try and set up a template-operator ;
% for a collection of azimu_b associated with a single polar_a. ;
%%%%%%%%;
polar_a_use = viewing_polar_a_S_(round(n_S/4));
tmp_index_ = efind(abs(viewing_polar_a_S_-polar_a_use)<1e-6);
n_azimu_b_use = numel(tmp_index_);
azimu_b_use_ = viewing_azimu_b_S_(1+tmp_index_);
S_k_p_sub_wkb__ = zeros(n_w_sum,n_azimu_b_use);
T_k_p_sub_wkb__ = zeros(n_w_sum,n_azimu_b_use);
S_k_q_sub_wkb__ = zeros(n_w_sum,n_azimu_b_use);
T_k_q_sub_wkb__ = zeros(n_w_sum,n_azimu_b_use);
for nazimu_b_use=0:n_azimu_b_use-1;
nS = tmp_index_(1+nazimu_b_use);
S_k_p_wk_ = S_k_p_wkS__(:,1+nS);
tmp_azimu_b = viewing_azimu_b_S_(1+nS);
tmp_polar_a = viewing_polar_a_S_(1+nS);
tmp_gamma_z = 0.0; %<-- default. ;
tmp_gamma_z = pi/12;
S_k_p_wk_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_wk_,-tmp_gamma_z);
S_k_p_sub_wkb__(:,1+nazimu_b_use) = S_k_p_wk_;
tmp_R__ = Rz(-tmp_gamma_z)*Ry(-tmp_polar_a)*Rz(-tmp_azimu_b);
T_k_p_wk_ = zeros(n_w_sum,1);
for nsource=0:n_source-1;
tmp_delta_ = tmp_R__*delta_a_c__(:,1+nsource);
T_k_p_wk_ = T_k_p_wk_ + exp(+i*2*pi*(k_c_0_wk_*tmp_delta_(1+0) + k_c_1_wk_*tmp_delta_(1+1)));
end;%for nsource=0:n_source-1;
T_k_p_sub_wkb__(:,1+nazimu_b_use) = T_k_p_wk_;
S_k_q_wk_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,S_k_p_wk_);
T_k_q_wk_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,T_k_p_wk_);
S_k_q_sub_wkb__(:,1+nazimu_b_use) = S_k_q_wk_;
T_k_q_sub_wkb__(:,1+nazimu_b_use) = T_k_q_wk_;
clear S_k_p_wk_ T_k_p_wk_;
end;%for nazimu_b_use=0:n_azimu_b_use-1;
disp(sprintf(' %% S_k_p_sub_wkb__ vs T_k_p_sub_wkb__: %0.16f %%<-- should be <1e-2',fnorm(S_k_p_sub_wkb__-T_k_p_sub_wkb__)/fnorm(S_k_p_sub_wkb__)));
disp(sprintf(' %% S_k_q_sub_wkb__ vs T_k_q_sub_wkb__: %0.16f %%<-- should be <1e-2',fnorm(S_k_q_sub_wkb__-T_k_q_sub_wkb__)/fnorm(S_k_q_sub_wkb__)));
%%%%;
flag_plot=0;
if flag_plot;
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 6; p_col = ceil(n_azimu_b_use/p_row); np=0;
for nazimu_b_use=0:n_azimu_b_use-1;
subplot(p_row,p_col,1+np);np=np+1;cla;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(S_k_p_sub_wkb__(:,1+nazimu_b_use))); axis image; axisnotick;
%title(sprintf('real(S_k_p_sub_wkb__(:,1+%d))',nazimu_b_use),'Interpreter','none');
title(sprintf('nazimu_b_use %d',nazimu_b_use),'Interpreter','none');
end;%for nazimu_b_use=0:n_azimu_b_use-1;
end;%if flag_plot;
%%%%;
tmp_t = tic();
W_beta__ = wignerd_b(l_max_max,-polar_a_use);
zeta_lm__ = zeros(1+l_max_max,n_m_max);
for l_val=0:l_max_max;
a1=((2*l_val+1)/(4*pi));
Llm__ = legendre(l_val,0,'unnorm');
for m_val=-l_val:+l_val;
if (l_val >0); Llm_ = Llm__(1+abs(m_val),:); end; if (l_val==0); Llm_ = Llm__; end; assert(numel(Llm_)==1);
a2=exp(lfactorial(l_val-abs(m_val)) - lfactorial(l_val+abs(m_val))); c=sqrt(a1*a2); s=1; % original phase ;
zeta_lm__(1+l_val,1+l_max_max+m_val) = s*c*Llm_(1+0);
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max_max;
%%%%;
W_betazeta_mlm___ = zeros(n_m_max,1+l_max_max,n_m_max);
for l_val=0:l_max_max;
for m0_val=-l_val:+l_val;
for m1_val=-l_val:+l_val;
W_betazeta_mlm___(1+l_max_max+m0_val,1+l_val,1+l_max_max+m1_val) = ...
 W_beta__{1+l_val}(1+l_val+m0_val,1+l_val+m1_val) ...
*zeta_lm__(1+l_val,1+l_max_max+m0_val) ...
 ;
end;%for m1_val=-l_val:+l_val;
end;%for m0_val=-l_val:+l_val;
end;%for l_val=0:l_max_max;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% W_betazeta_mlm___: %0.2fs',tmp_t)); end;
%%%%;
flag_check=1;
if flag_check;
tmp_w_ = crandn(n_m_max);
tmp_azimu_b_use_ = 2*pi*rand(n_azimu_b_use,1);
tmp_f__ = exp(-i*(-reshape(tmp_azimu_b_use_,[n_azimu_b_use,1]))*reshape(m_max_,[1,n_m_max]));
tmp_fw_0_ = tmp_f__*tmp_w_;
tmp_fw_1_ = xxnufft1d2(n_azimu_b_use,tmp_azimu_b_use_,+1,1e-6,n_m_max,tmp_w_);
disp(sprintf(' %% tmp_fw_0_ vs tmp_fw_1_: %0.16f %%<-- should be <1e-6',fnorm(tmp_fw_0_-tmp_fw_1_)/fnorm(tmp_fw_0_)));
clear tmp_azimu_b_use_ tmp_w_ tmp_f__ tmp_fw_0_ tmp_fw_1_ ;
end;%if flag_check;
%%%%;
tmp_t = tic();
W_caza_mkm___ = zeros(n_m_max,n_k_p_r,n_m_max); %<-- diag(exp(+i*m0_val_*tmp_gamma_z))*W_betazeta_ml__*a_k_Y_form_lk__ for each m1_val. ;
for m1_val=-l_max_max:+l_max_max;
W_caza_mkm___(:,:,1+l_max_max+m1_val) = ...
 (pi^2) ...
 *diag(exp(-i*m_max_*(-tmp_gamma_z))) ...
 *reshape(W_betazeta_mlm___(:,:,1+l_max_max+m1_val),[n_m_max,1+l_max_max]) ...
 *reshape(a_k_Y_form_lkm___(:,:,1+l_max_max+m1_val),[1+l_max_max,n_k_p_r]) ...
 ;
end;%for m1_val=-l_max_max:+l_max_max;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% W_caza_mkm___: %0.2fs',tmp_t)); end;
tmp_t = tic();
W_caza_mmk___ = permute(W_caza_mkm___,[3,1,2]);
W_caza_bmk___ = reshape(xxnufft1d2(n_azimu_b_use,azimu_b_use_,+1,1e-6,n_m_max,reshape(W_caza_mmk___,[n_m_max,n_m_max*n_k_p_r])),[n_azimu_b_use,n_m_max,n_k_p_r]);
W_caza_mkb___ = permute(W_caza_bmk___,[2,3,1]);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% W_caza_mkb___: %0.2fs',tmp_t)); end;
%%%%%%%%;
R_k_q_sub_wkb__ = zeros(n_w_sum,n_azimu_b_use);
R_k_p_sub_wkb__ = zeros(n_w_sum,n_azimu_b_use);
for nazimu_b_use=0:n_azimu_b_use-1;
R_k_p_wk_ = zeros(n_w_sum,1);
R_k_q_wk_ = zeros(n_w_sum,1);
for nk_p_r=0:n_k_p_r-1;
for m_val=-l_max_max:+l_max_max;
nq = m_val; if (nq<0); nq=nq+n_w_max; end;
R_k_q_wk_(1+nq+nk_p_r*n_w_max) = W_caza_mkb___(1+l_max_max+m_val,1+nk_p_r,1+nazimu_b_use);
end;%for m_val=-l_max_max:+l_max_max;
end;%for nk_p_r=0:n_k_p_r-1;
R_k_q_sub_wkb__(:,1+nazimu_b_use) = R_k_q_wk_;
R_k_p_wk_ = interp_q_to_p(n_k_p_r,n_w_,n_w_sum,R_k_q_wk_);
R_k_p_sub_wkb__(:,1+nazimu_b_use) = R_k_p_wk_;
end;%for nazimu_b_use=0:n_azimu_b_use-1;
disp(sprintf(' %% S_k_p_sub_wkb__ vs T_k_p_sub_wkb__: %0.16f %%<-- should be <1e-2',fnorm(S_k_p_sub_wkb__ - T_k_p_sub_wkb__)/fnorm(S_k_p_sub_wkb__)));
disp(sprintf(' %% S_k_p_sub_wkb__ vs R_k_p_sub_wkb__: %0.16f %%<-- should be <1e-2',fnorm(S_k_p_sub_wkb__ - R_k_p_sub_wkb__)/fnorm(S_k_p_sub_wkb__)));
disp(sprintf(' %% T_k_p_sub_wkb__ vs R_k_p_sub_wkb__: %0.16f %%<-- should be <1e-2',fnorm(T_k_p_sub_wkb__ - R_k_p_sub_wkb__)/fnorm(T_k_p_sub_wkb__)));
%%%%;
flag_plot=0;
if flag_plot;
%%;
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 6; p_col = ceil(n_azimu_b_use/p_row); np=0;
for nazimu_b_use=0:n_azimu_b_use-1;
subplot(p_row,p_col,1+np);np=np+1;cla;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(S_k_p_sub_wkb__(:,1+nazimu_b_use))); axis image; axisnotick;
title(sprintf('nazimu_b_use %d',nazimu_b_use),'Interpreter','none');
end;%for nazimu_b_use=0:n_azimu_b_use-1;
%%;
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 6; p_col = ceil(n_azimu_b_use/p_row); np=0;
for nazimu_b_use=0:n_azimu_b_use-1;
subplot(p_row,p_col,1+np);np=np+1;cla;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(R_k_p_sub_wkb__(:,1+nazimu_b_use))); axis image; axisnotick;
title(sprintf('nazimu_b_use %d',nazimu_b_use),'Interpreter','none');
end;%for nazimu_b_use=0:n_azimu_b_use-1;
%%;
end;%if flag_plot;
%%%%%%%%;
%%%%%%%%;
[ ...
 U_k_p_sub_wkb__ ...
,tmp_W_betazeta_mlm___ ...
,tmp_a_k_Y_lkm___ ...
,tmp_W_caza_mkm___ ...
] = ...
sph_template_single_polar_a_3( ...
 verbose ...
,l_max ...
,n_k_p_r ...
,a_k_Y_form_lkm___ ...
,n_w_max ...
,polar_a_use ...
,n_azimu_b_use ...
,azimu_b_use_ ...
,tmp_gamma_z ...
);
%%%%%%%%;
figure(1+nf);nf=nf+1;plot(W_betazeta_mlm___(:),tmp_W_betazeta_mlm___(:),'.');
figure(1+nf);nf=nf+1;plot(a_k_Y_form_lkm___(:),tmp_a_k_Y_lkm___(:),'.');
figure(1+nf);nf=nf+1;plot(W_caza_mkm___(:),tmp_W_caza_mkm___(:),'.');

disp('returning'); return;

verbose=1;
flag_disp=1;

equa_band_dilated_amplitude = 0.15;
g_dilation = @(point_pole_predilated_azimu_b) equa_band_dilated_amplitude*sin(2*point_pole_predilated_azimu_b); %<-- approximation holds well for first nontrivial mode. ;
f_dilation = @(point_pole_predilated_azimu_b) point_pole_predilated_azimu_b + g_dilation(point_pole_predilated_azimu_b);

npoint_a = 32;
npoint_b = 12;
markersize_sml = 8;
markersize_big = 16;

if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;
plot_sphere_grid_0;
axis equal; axis vis3d;
hold on;
end;%if flag_disp;

point_output_azimu_b = point_output_azimu_b_b_(1+npoint_b); %<-- yes periodic. ;
point_output_polar_a = -pi/2 + point_output_polar_a_a_(1+npoint_a); %<-- not periodic. ;
point_output_gamma_z = 0;
point_output_k_c_ = Rz(point_output_azimu_b) * Ry(point_output_polar_a) * Rz(point_output_gamma_z) * [1;0;0];
point_pole_k_c__ = Rz(point_output_azimu_b) * Ry(point_output_polar_a) * Rz(point_output_gamma_z) * [zeros(1,n_w_max);transpose(sc_);transpose(cc_)];
%%%%%%%%;
if flag_disp;
plot3(point_output_k_c_(1+0),point_output_k_c_(1+1),point_output_k_c_(1+2),'mo','MarkerFaceColor',0.85*[1,1,1],'MarkerSize',markersize_big);
end;%if flag_disp;
%%%%%%%%%%%%%%%%;
for npole=floor(n_w_max/3);%for npole=0:n_w_max-1;
%%%%%%%%%%%%%%%%;
% Now we step through each of the templates associated with the point point_output. ;
%%%%%%%%;
if flag_disp;
nc_hsv = max(0,min(n_c_hsv-1,floor(n_c_hsv*(periodize(npole,0,n_w_max/2)/(n_w_max/2))))); %<-- note double winding (i.e., andipodal templates are the same). ;
end;%if flag_disp;
point_pole_k_c_ = point_pole_k_c__(:,1+npole);
if flag_disp;
plot3(point_pole_k_c_(1+0),point_pole_k_c_(1+1),point_pole_k_c_(1+2),'o','MarkerEdgeColor','g','MarkerFaceColor',c_hsv__(1+nc_hsv,:),'MarkerSize',markersize_big);
end;%if flag_disp;
point_pole_k_r01 = sqrt(point_pole_k_c_(1+0).^2 + point_pole_k_c_(1+1).^2);
point_pole_k_r012 = sqrt(point_pole_k_c_(1+0).^2 + point_pole_k_c_(1+1).^2 + point_pole_k_c_(1+2).^2);
point_pole_azimu_b = atan2(point_pole_k_c_(1+1),point_pole_k_c_(1+0));
point_pole_polar_a = atan2(point_pole_k_r01,point_pole_k_c_(1+2));
[ ...
 point_pole_template_k_c_0_w_ ...
,point_pole_template_k_c_1_w_ ...
,point_pole_template_k_c_2_w_ ...
,point_pole_template_azimu_b_w_ ...
,point_pole_template_polar_a_w_ ...
] = ...
get_template_single_ring_k_c_0( ...
 verbose ...
,1 ...
,point_pole_azimu_b ...
,point_pole_polar_a ...
,n_w_max ...
);
if flag_disp;
plot3( ...
 point_pole_template_k_c_0_w_ ...
,point_pole_template_k_c_1_w_ ...
,point_pole_template_k_c_2_w_ ...
,'go','MarkerFaceColor',0.95*[1,1,1],'MarkerSize',markersize_sml);
end;%if flag_disp;
%%%%%%%%%%%%%%%%;
end;%for npole=0:n_w_max-1;
%%%%%%%%%%%%%%%%;
%%%%%%%%;
% Here we determine the predilated template that is associated with the original template mapped to point_output. ;
%%%%%%%%;
tmp_f_error = @(point_pole_predilated_azimu_b) abs(f_dilation(point_pole_predilated_azimu_b) - point_pole_azimu_b).^2;
point_pole_predilated_azimu_b = fminsearch(tmp_f_error,point_pole_azimu_b);
tmp_point_pole_azimu_b = point_pole_predilated_azimu_b + equa_band_dilated_amplitude*sin(2*point_pole_predilated_azimu_b);
if (verbose>0); disp(sprintf(' %% npole %d/%d, fnorm(point_pole_azimu_b-tmp_point_pole_azimu_b): %0.16f',npole,n_w_max,fnorm(point_pole_azimu_b-tmp_point_pole_azimu_b))); end;
point_pole_predilated_polar_a = point_pole_polar_a;
point_pole_predilated_k_c_ = [ ...
  cos(point_pole_predilated_azimu_b)*sin(point_pole_predilated_polar_a) ...
; sin(point_pole_predilated_azimu_b)*sin(point_pole_predilated_polar_a) ...
; cos(point_pole_predilated_polar_a) ...
];
if flag_disp;
plot3(point_pole_predilated_k_c_(1+0),point_pole_predilated_k_c_(1+1),point_pole_predilated_k_c_(1+2),'o','MarkerEdgeColor','c','MarkerFaceColor',c_hsv__(1+nc_hsv,:),'MarkerSize',markersize_sml);
end;%if flag_disp;
[ ...
 point_pole_predilated_template_k_c_0_w_ ...
,point_pole_predilated_template_k_c_1_w_ ...
,point_pole_predilated_template_k_c_2_w_ ...
,point_pole_predilated_template_azimu_b_w_ ...
,point_pole_predilated_template_polar_a_w_ ...
] = ...
get_template_single_ring_k_c_0( ...
 verbose ...
,1 ...
,point_pole_predilated_azimu_b ...
,point_pole_predilated_polar_a ...
,n_w_max ...
);
if flag_disp;
plot3( ...
 point_pole_predilated_template_k_c_0_w_ ...
,point_pole_predilated_template_k_c_1_w_ ...
,point_pole_predilated_template_k_c_2_w_ ...
,'co','MarkerFaceColor',0.95*[1,1,1],'MarkerSize',markersize_sml);
end;%if flag_disp;
%%%%%%%%;
point_pole_template_gamma0_k_c_ = [...
 +cos(point_pole_azimu_b)*cos(point_pole_polar_a) ...
;+sin(point_pole_azimu_b)*cos(point_pole_polar_a) ...
;-sin(point_pole_polar_a) ...
];
if flag_disp;
plot3( ...
 point_pole_template_gamma0_k_c_(1) ...
,point_pole_template_gamma0_k_c_(2) ...
,point_pole_template_gamma0_k_c_(3) ...
,'ko','MarkerFaceColor',0.95*[1,1,1],'MarkerSize',markersize_big);
end;%if flag_disp;
point_pole_template_sgx_k_c_ = cross(point_pole_template_gamma0_k_c_,point_output_k_c_);
point_pole_template_sgx = dot(point_pole_template_sgx_k_c_,point_pole_k_c_);
point_pole_template_cgx = dot(point_pole_template_gamma0_k_c_,point_output_k_c_);
point_pole_gx = atan2(point_pole_template_sgx,point_pole_template_cgx);
point_pole_template_gammax_k_c_ = [...
 +cos(point_pole_azimu_b)*cos(point_pole_polar_a)*cos(point_pole_gx) - sin(point_pole_azimu_b)*sin(point_pole_gx) ...
;+sin(point_pole_azimu_b)*cos(point_pole_polar_a)*cos(point_pole_gx) + cos(point_pole_azimu_b)*sin(point_pole_gx) ...
;-sin(point_pole_polar_a)*cos(point_pole_gx) ...
];
if flag_disp;
plot3( ...
 point_pole_template_gammax_k_c_(1) ...
,point_pole_template_gammax_k_c_(2) ...
,point_pole_template_gammax_k_c_(3) ...
,'ro','MarkerFaceColor',0.95*[1,1,1],'MarkerSize',markersize_big);
end;%if flag_disp;
% fnorm(point_pole_template_gammax_k_c_ - point_output_k_c_),; %<-- should be 0. ;
point_pole_template_gammax_k_r01 = sqrt(point_pole_template_gammax_k_c_(1+0).^2 + point_pole_template_gammax_k_c_(1+1).^2);
% fnorm((cos(point_pole_polar_a)*cos(point_pole_gx))^2 + (sin(point_pole_gx))^2 - point_pole_template_gammax_k_r01^2),; %<-- should be 0. ;
point_pole_template_gammax_k_r012 = sqrt(point_pole_template_gammax_k_c_(1+0).^2 + point_pole_template_gammax_k_c_(1+1).^2 + point_pole_template_gammax_k_c_(1+2).^2);
point_pole_template_gammax_azimu_b = atan2(point_pole_template_gammax_k_c_(1+1),point_pole_template_gammax_k_c_(1+0));
point_pole_template_gammax_polar_a = atan2(point_pole_template_gammax_k_r01,point_pole_template_gammax_k_c_(1+2));
%%%%%%%%;
% Now we determine the gamma_z (denoted gammax) within predilated template that maps to point_output. ;
%%%%%%%%;
point_pole_predilated_template_gammax_k_c_ = [ ...
 +cos(point_pole_predilated_azimu_b)*cos(point_pole_predilated_polar_a)*cos(point_pole_gx) - sin(point_pole_predilated_azimu_b)*sin(point_pole_gx) ...
;+sin(point_pole_predilated_azimu_b)*cos(point_pole_predilated_polar_a)*cos(point_pole_gx) + cos(point_pole_predilated_azimu_b)*sin(point_pole_gx) ...
;-sin(point_pole_predilated_polar_a)*cos(point_pole_gx) ...
];
if flag_disp;
plot3( ...
 point_pole_predilated_template_gammax_k_c_(1) ...
,point_pole_predilated_template_gammax_k_c_(2) ...
,point_pole_predilated_template_gammax_k_c_(3) ...
,'ro','MarkerFaceColor',0.95*[1,1,1],'MarkerSize',markersize_sml);
end;%if flag_disp;

error('stopping');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% plotting original micrograph Q_x_u_pack_ in test_pm_clathrin_6.m ;
%%%%%%%%;
tmp_Q_x_u_pack_ = Q_x_u_pack_;
tmp_p = prctile(tmp_Q_x_u_pack_,100,'all');
tmp_q = prctile(tmp_Q_x_u_pack_,  0,'all');
n_p = 10;
for np=0:n_p-1;
tmp_index = round(size(Q_x_u_pack_,1)*np/n_p);
tmp_Q_x_u_pack_(1+tmp_index,:) = tmp_q;
tmp_Q_x_u_pack_(:,1+tmp_index) = tmp_q;
end;%for np=0:n_p-1;
tmp_Q_x_u_pack_(min(R_sub_ij_0_),R_sub_ij_1_) = tmp_p;
tmp_Q_x_u_pack_(max(R_sub_ij_0_),R_sub_ij_1_) = tmp_p;
tmp_Q_x_u_pack_(R_sub_ij_0_,min(R_sub_ij_1_)) = tmp_p;
tmp_Q_x_u_pack_(R_sub_ij_0_,max(R_sub_ij_1_)) = tmp_p;
fname_fig = sprintf('%s_jpg/Q_sub_x_u_pack_FIGA',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;colormap('gray');
subplot(1,2,1); imagesc(tmp_Q_x_u_pack_); axis image; axisnotick; title('Q','Interpreter','none');
subplot(1,2,2); imagesc(Q_sub_x_u_pack_); axis image; axisnotick; title('Q_sub','Interpreter','none');
sgtitle(fname_fig,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
%close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));

%%%%%%%%;
% experimenting with str_shape. ;
%%%%%%%%;

%%%%%%%%;
% Sum of gaussians. ;
%%%%%%%%;
%n_source_gaussian = 128;
tmp_prefactor = 1.000;
tmp_sigma_x_c = 0.025;
tmp_sigma_k_p = 1/tmp_sigma_x_c;
%%%%;
if strcmp(str_shape,'');
tmp_delta_2s__ = 1.25/sqrt(2) ...
  *[ ...
   cos(linspace(0,2*pi,1+n_source_gaussian)) ...
 ; sin(linspace(0,2*pi,1+n_source_gaussian)) ...
   ];
tmp_prefactor_ = reshape(sqrt(sum(diff(tmp_delta_2s__,1,2).^2,1)),[n_source_gaussian,1]);
tmp_delta_2s__ = tmp_delta_2s__(:,1:n_source_gaussian);
tmp_sigma_x_c_ = tmp_sigma_x_c*ones(n_source_gaussian,1);
end;%if strcmp(str_shape,'');
%%%%;
if strcmp(str_shape,'C2');
tmp_delta_2s__ = 1.25/sqrt(2) ...
  *[ ...
   cos(1*linspace(0,2*pi,1+n_source_gaussian)) ...
 ; sin(2*linspace(0,2*pi,1+n_source_gaussian)) ...
   ];
tmp_prefactor_ = reshape(sqrt(sum(diff(tmp_delta_2s__,1,2).^2,1)),[n_source_gaussian,1]);
tmp_delta_2s__ = tmp_delta_2s__(:,1:n_source_gaussian);
tmp_sigma_x_c_ = tmp_sigma_x_c*ones(n_source_gaussian,1);
tmp_delta_2s__ = [ tmp_delta_2s__ , [+0.6;+0.45] , [-0.6;-0.45] ];
tmp_sigma_x_c_ = [tmp_sigma_x_c_;0.25;0.25];
tmp_prefactor_ = [tmp_prefactor_;sum(tmp_prefactor_)/2;sum(tmp_prefactor_)/2];
n_source_gaussian = n_source_gaussian + 2;
end;%if strcmp(str_shape,'C2');
%%%%;
if strcmp(str_shape,'C4');
tmp_delta_2s__ = 1.25/sqrt(2) ...
  *[ ...
   cos(1*linspace(0,2*pi,1+n_source_gaussian)) ...
 ; sin(2*linspace(0,2*pi,1+n_source_gaussian)) ...
   ];
tmp_prefactor_ = reshape(sqrt(sum(diff(tmp_delta_2s__,1,2).^2,1)),[n_source_gaussian,1]);
tmp_delta_2s__ = tmp_delta_2s__(:,1:n_source_gaussian);
tmp_sigma_x_c_ = tmp_sigma_x_c*ones(n_source_gaussian,1);
tmp_R__ = [0,-1;+1,0];
tmp_delta_2s__ = [tmp_delta_2s__ , tmp_R__*tmp_delta_2s__];
tmp_prefactor_ = [tmp_prefactor_;tmp_prefactor_];
tmp_sigma_x_c_ = [tmp_sigma_x_c_;tmp_sigma_x_c_];
n_source_gaussian = n_source_gaussian + n_source_gaussian;
tmp_delta_2s__ = [ tmp_delta_2s__ , [+0.6;+0.45] , [+0.45;-0.6] , [-0.6;-0.45] , [-0.45;+0.6] ];
tmp_sigma_x_c_ = [tmp_sigma_x_c_;0.125;0.125;0.125;0.125];
tmp_prefactor_ = [tmp_prefactor_;sum(tmp_prefactor_)/4;sum(tmp_prefactor_)/4;sum(tmp_prefactor_)/4;sum(tmp_prefactor_)/4];
n_source_gaussian = n_source_gaussian + 4;
end;%if strcmp(str_shape,'C4');
%%%%;
tmp_N_x_c_ = zeros(n_x_c,n_x_c);
for nsource_gaussian=0:n_source_gaussian-1;
tmp_delta_ = tmp_delta_2s__(:,1+nsource_gaussian);
tmp_sigma_x_c = tmp_sigma_x_c_(1+nsource_gaussian);
tmp_prefactor = tmp_prefactor_(1+nsource_gaussian);
tmp_M_x_c_ = tmp_prefactor * 1/(sqrt(2*pi)*tmp_sigma_x_c)^2 * exp( -( (x_c_0_01__-tmp_delta_(1+0)).^2 + (x_c_1_01__-tmp_delta_(1+1)).^2 ) / (2*tmp_sigma_x_c^2) );
tmp_N_x_c_ = tmp_N_x_c_ + tmp_M_x_c_;
clear tmp_delta_;
end;%for nsource_gaussian=0:n_source_gaussian-1;
tmp_N_x_c_l2 = sum(tmp_N_x_c_.^2,'all')*dx^2;
disp(sprintf(' %% sum(tmp_N_x_c_*dx^2,''all'') = %0.16f',sum(tmp_N_x_c_*dx^2,'all')));
disp(sprintf(' %% tmp_N_x_c_l2 = %0.16f',tmp_N_x_c_l2));
tmp_N_k_p_ = interp_x_c_to_k_p_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,tmp_N_x_c_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_c^2)*dx^2;
tmp_N_k_p_l2 = sum(abs(tmp_N_k_p_).^2 .* weight_2d_k_all_) * (2*pi)^2;
disp(sprintf(' %% tmp_N_k_p_l2 = %0.16f',tmp_N_k_p_l2));
tmp_N_k_p_form_ = zeros(n_w_sum,1);
for nsource_gaussian=0:n_source_gaussian-1;
tmp_delta_ = tmp_delta_2s__(:,1+nsource_gaussian);
tmp_sigma_x_c = tmp_sigma_x_c_(1+nsource_gaussian);
tmp_prefactor = tmp_prefactor_(1+nsource_gaussian);
tmp_M_k_p_form_ = zeros(n_w_sum,1);
na=0;
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
n_w = n_w_(1+nk_p_r);
for nw=0:n_w-1;
k_x_c_0 = k_p_r*cos(2*pi*nw/n_w);
k_x_c_1 = k_p_r*sin(2*pi*nw/n_w);
tmp_M_k_p_form_(1+na) = tmp_prefactor * exp( -( (2*pi*k_x_c_0).^2 + (2*pi*k_x_c_1).^2 ) / (2/tmp_sigma_x_c^2) ) .* exp( - 2*pi*i*( k_x_c_0*tmp_delta_(1+0) + k_x_c_1*tmp_delta_(1+1) ) );
na=na+1;
end;%for nw=0:n_w-1;
end;%for nk_p_r=0:n_k_p_r-1;
tmp_N_k_p_form_ = tmp_N_k_p_form_ + tmp_M_k_p_form_;
clear tmp_delta_;
end;%for nsource_gaussian=0:n_source_gaussian-1;
%%%%%%%%;
n_x_p_r = n_k_p_r;
x_p_r_ = k_p_r_*x_p_r_max/k_p_r_max;
x_c_0_all_ = k_c_0_all_*x_p_r_max/k_p_r_max;
x_c_1_all_ = k_c_1_all_*x_p_r_max/k_p_r_max;
%%%%%%%%;
tmp_R_x_p_form_ = ...
radon_k_p_to_x_p_xxnufft( ...
 n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,tmp_N_k_p_form_ ...
);
tmp_N_k_p_form_l2 = sum(abs(tmp_N_k_p_form_).^2 .* weight_2d_k_all_) * (2*pi)^2;
disp(sprintf(' %% tmp_N_k_p_form_l2 = %0.16f',tmp_N_k_p_form_l2));
disp(sprintf(' %% tmp_N_k_p_ vs tmp_N_k_p_form: %0.16f',fnorm(tmp_N_k_p_ - tmp_N_k_p_form_)/fnorm(tmp_N_k_p_)));
tmp_N_x_c_reco_ = interp_k_p_to_x_c_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_N_k_p_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_c^2) * n_w_sum;
tmp_N_x_c_reco_l2 = sum(abs(tmp_N_x_c_reco_).^2,'all')*dx^2;
disp(sprintf(' %% tmp_N_x_c_reco_l2 = %0.16f',tmp_N_x_c_reco_l2));
disp(sprintf(' %% tmp_N_x_c_ vs tmp_N_x_c_reco: %0.16f',fnorm(tmp_N_x_c_ - tmp_N_x_c_reco_)/fnorm(tmp_N_x_c_)));
%%%%%%%%;
flag_disp=1;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figbig;
subplot(2,2,1);imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,tmp_N_x_c_);axis image;axisnotick;
subplot(2,2,2);imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(tmp_N_k_p_));axis image;axisnotick;
subplot(2,2,3);imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(tmp_N_k_p_form_));axis image;axisnotick;
subplot(2,2,4);imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(tmp_N_x_c_reco_));axis image;axisnotick;
error('stopping');
end;%if flag_disp;
%%%%%%%%;

error('stopping');

%%%%%%%%;
% visualize shadowsphere. ;
%%%%%%%%;
gamma_z_ = linspace(0,2*pi,n_w_max+1); gamma_z_ = transpose(gamma_z_(1:n_w_max));
q_ = periodize(transpose(0:n_w_max-1),-n_w_max/2,+n_w_max/2);
%%%%%%%%;
tmp_p_from_q_wq__ = zeros(n_w_max,n_w_max);
for nq=0:n_w_max-1;
tmp_q = q_(1+nq);
tmp_p_from_q_wq__(:,1+nq) = exp(+i*gamma_z_*tmp_q)/sqrt(n_w_max);
end;%for nq=0:n_w_max-1;
tmp_N_k_q_form_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,tmp_N_k_p_form_);
tmp_N_k_p_reco_ = reshape(tmp_p_from_q_wq__*reshape(tmp_N_k_q_form_,[n_w_max,n_k_p_r]),[n_w_sum,1]);
disp(sprintf(' %% tmp_N_k_p_form_ vs tmp_N_k_p_reco_: %0.16f',fnorm(tmp_N_k_p_form_-tmp_N_k_p_reco_)/fnorm(tmp_N_k_p_form_)));
%%%%;
tmp_O_k_p_pre_ = tmp_N_k_p_form_;
tmp_O_k_q_pre_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,tmp_O_k_p_pre_);
tmp_R_x_p_pre_ = ...
radon_k_p_to_x_p_xxnufft( ...
 n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,tmp_O_k_p_pre_ ...
);
flag_disp=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if flag_disp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
N_k_p_form_lim_ = prctile(abs(real(tmp_N_k_p_form_)), 95,'all')*[-1,+1];
N_x_c_reco_lim_ = prctile(abs(real(tmp_N_x_c_reco_)),100,'all')*[-1,+1];
R_x_p_form_lim_ = prctile(abs(real(tmp_R_x_p_form_)), 95,'all')*[-1,+1];
tmp_O_x_c_pre_ = interp_k_p_to_x_c_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_O_k_p_pre_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_c^2) * n_w_sum;
figure(1);clf;figbig; strip_width_use = 1/16;
%%%%%%%%%%%%%%%%;
subplot(1,2,1); hold on;
%%%%%%%%;
imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(tmp_N_x_c_reco_),N_x_c_reco_lim_);
%%%%%%%%;
n_l = 16;
for nl=flip(0:n_l-1);
tmp_omega = (2*pi)*nl/n_l;
tmp_x_0_offset = sqrt(1.0+0.5*mod(nl,2))*diameter_x_c*cos(tmp_omega+pi/2);
tmp_x_1_offset = sqrt(1.0+0.5*mod(nl,2))*diameter_x_c*sin(tmp_omega+pi/2);
tmp_parameter = struct('type','parameter');
%tmp_parameter.arrow_head_length = 1.0+2.0*mod(nl,2);
%tmp_parameter.arrow_tail_length = 2.5+1.0*mod(nl,2);
tmp_parameter.arrow_head_length = 2.25+1.0*mod(nl,2);
tmp_parameter.arrow_tail_length = 1.0+2.0*mod(nl,2);
plot_arrow_0(tmp_parameter,[tmp_x_0_offset;tmp_x_1_offset],tmp_omega+pi/2,0.25+0.00*mod(nl,2));
tmp_x_0_offset = sqrt(2.0+1.0*mod(nl,2))*diameter_x_c*cos(tmp_omega+pi/2); tmp_x_1_offset = sqrt(2.0+1.0*mod(nl,2))*diameter_x_c*sin(tmp_omega+pi/2);
imagesc_p_strip_0( ...
 struct('type','parameter','strip_width',strip_width_use,'k_0_offset',tmp_x_0_offset,'k_1_offset',tmp_x_1_offset,'clim_',R_x_p_form_lim_) ...
,n_k_p_r ...
,k_p_r_*x_p_r_max/k_p_r_max ...
,k_p_r_max*x_p_r_max/k_p_r_max ...
,n_w_ ...
,n_w_sum ...
,tmp_R_x_p_pre_ ...
,1 ...
,tmp_omega ...
);
tmp_x_0_offset = (sqrt(2.0+1.0*mod(nl,2))*diameter_x_c + 2*strip_width_use)*cos(tmp_omega+pi/2); tmp_x_1_offset = (sqrt(2.0+1.0*mod(nl,2))*diameter_x_c + 2*strip_width_use)*sin(tmp_omega+pi/2);
imagesc_p_strip_0( ...
 struct('type','parameter','strip_width',strip_width_use,'k_0_offset',tmp_x_0_offset,'k_1_offset',tmp_x_1_offset,'clim_',N_k_p_form_lim_,'c_use__',colormap_80s) ...
,n_k_p_r ...
,k_p_r_*x_p_r_max/k_p_r_max ...
,k_p_r_max*x_p_r_max/k_p_r_max ...
,n_w_ ...
,n_w_sum ...
,tmp_O_k_p_pre_ ...
,1 ...
,tmp_omega ...
);
end;%for nl=flip(0:n_l-1);
%%%%%%%%;
hold off;
xlim([-4,+4]);
ylim([-4,+4]);
axis square;
axisnotick;
%%%%%%%%%%%%%%%%;
subplot(1,2,2); hold on;
%%%%%%%%;
n_l = 16;
for nl=flip(0:n_l-1);
tmp_omega = (2*pi)*nl/n_l;
tmp_x_0_offset = sqrt(1.0+1.0*mod(nl,2))*diameter_x_c*cos(tmp_omega+pi/2); tmp_x_1_offset = sqrt(1.0+1.0*mod(nl,2))*diameter_x_c*sin(tmp_omega+pi/2);
tmp_x_0_offset = sqrt(1.0+0.5*mod(nl,2))*diameter_x_c*cos(tmp_omega+pi/2);
tmp_x_1_offset = sqrt(1.0+0.5*mod(nl,2))*diameter_x_c*sin(tmp_omega+pi/2);
tmp_parameter = struct('type','parameter');
tmp_parameter.arrow_head_length = 1.0+2.0*mod(nl,2);
tmp_parameter.arrow_tail_length = 2.5+1.0*mod(nl,2);
plot_arrow_0(tmp_parameter,[tmp_x_0_offset;tmp_x_1_offset],tmp_omega+pi/2+pi,0.25+0.00*mod(nl,2));
tmp_x_0_offset = sqrt(2.0+1.0*mod(nl,2))*diameter_x_c*cos(tmp_omega+pi/2); tmp_x_1_offset = sqrt(2.0+1.0*mod(nl,2))*diameter_x_c*sin(tmp_omega+pi/2);
imagesc_p_strip_0( ...
 struct('type','parameter','strip_width',strip_width_use,'k_0_offset',tmp_x_0_offset,'k_1_offset',tmp_x_1_offset,'clim_',R_x_p_form_lim_) ...
,n_k_p_r ...
,k_p_r_*x_p_r_max/k_p_r_max ...
,k_p_r_max*x_p_r_max/k_p_r_max ...
,n_w_ ...
,n_w_sum ...
,tmp_R_x_p_pre_ ...
,1 ...
,tmp_omega ...
);
end;%for nl=flip(0:n_l-1);
imagesc_p_strip_0( ...
 struct('type','parameter','strip_width',strip_width_use/2,'clim_',R_x_p_form_lim_) ...
,n_k_p_r ...
,k_p_r_*x_p_r_max/k_p_r_max ...
,k_p_r_max*x_p_r_max/k_p_r_max ...
,n_w_ ...
,n_w_sum ...
,tmp_R_x_p_pre_ ...
,1+n_l/2 ...
,linspace(-pi/2,+pi/2,1+n_l/2) ...
);
%%%%%%%%;
hold off;
xlim([-4,+4]);
ylim([-4,+4]);
axis square;
axisnotick;
%%%%%%%%%%%%%%%%;
drawnow(); %error('stop');
fname_fig_pre = sprintf('%s/MSA_shape_demonstration_ns%.2d',dir_manuscript_jpg,n_source_gaussian);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if ( flag_replot | ~exist(fname_fig_jpg,'file') );
sgtitle(fname_fig_pre,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
%close(gcf);
end;%if ( flag_replot | ~exist(fname_fig_jpg,'file') );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_disp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

error('stopping');

%%%%%%%%;
% try arrows. ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
subplot(1,2,1);
n_l = 16;
for nl=flip(0:n_l-1);
tmp_omega = (2*pi)*nl/n_l;
tmp_x_0_offset = sqrt(1.0+0.5*mod(nl,2))*diameter_x_c*cos(tmp_omega+pi/2);
tmp_x_1_offset = sqrt(1.0+0.5*mod(nl,2))*diameter_x_c*sin(tmp_omega+pi/2);
tmp_parameter = struct('type','parameter');
tmp_parameter.arrow_head_length = 1.0;
tmp_parameter.arrow_tail_length = 1.0+2.0*mod(nl,2);
plot_arrow_0(tmp_parameter,[tmp_x_0_offset;tmp_x_1_offset],tmp_omega+pi/2,0.25+0.00*mod(nl,2));
end;%for nl=flip(0:n_l-1);
hold off;
xlim([-4,+4]);
ylim([-4,+4]);
axis square;
axisnotick;
%%%%%%%%;

error('stopping');

%%%%%%%%;
% is max-likelihood with angle-distributions per image convex? ;
%%%%%%%%;
pmxgf = @(a,b,x,m,s) exp( -(a.*x + b - m).^2 ./ (2.*s.^2) ) / sqrt(2*pi) ./ s;
pmgf = @(a,b,xp,xq,p,m,s) p.*pmxgf(a,b,xp,m,s) + (1-p).*pmxgf(a,b,xq,m,s);
pm1m2gf = @(a,b,x1p,x1q,p1,m1,x2p,x2q,p2,m2,s) pmgf(a,b,x1p,x1q,p1,m1,s).*pmgf(a,b,x2p,x2q,p2,m2,s);
pm1m2m3gf = @(a,b,x1p,x1q,p1,m1,x2p,x2q,p2,m2,s) pmgf(a,b,x1p,x1q,p1,m1,s).*pmgf(a,b,x2p,x2q,p2,m2,s);
n_a = 128;
n_b = 128;
a_ = linspace(-5,+5,n_a);
b_ = linspace(-5,+5,n_a);
[a__,b__] = ndgrid(a_,b_);
n_contour = 32;

rseed_ = 48 + (0:48-1); n_rseed = numel(rseed_);
figure(1);clf;figbig;fig80s;
p_row = 6; p_col = ceil(n_rseed/p_row); np=0;
for nrseed=0:n_rseed-1;
rseed = rseed_(1+nrseed);
rng(rseed);
s=1.0;
m1 = 2*rand()-1; m2 = 2*rand()-1; m3 = 2*rand()-1;
x1p = 2*rand()-1; x1q = 2*rand()-1;
x2p = 2*rand()-1; x2q = 2*rand()-1;
x3p = 2*rand()-1; x3q = 2*rand()-1;
p1 = rand(); p2 = rand(); p3 = rand();
pm1m2gf__ = pm1m2gf(a__,b__,x1p,x1q,p1,m1,x2p,x2q,p1,m2,s);
nlpm1m2gf__ = -log(pm1m2gf__);
nlplim_ = prctile(nlpm1m2gf__,linspace(0,50,n_contour),'all');
subplot(p_row,p_col,1+nrseed);
contour(nlpm1m2gf__,nlplim_);
axis image; axisnotick;
drawnow();
end;%for nrseed=0:n_rseed-1;





