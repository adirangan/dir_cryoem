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





