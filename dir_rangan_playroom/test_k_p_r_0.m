%%%%%%%%;
% Sets up a simple test of polar-quadrature parameters. ;
%%%%%%%%;

platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;

flag_verbose = 1;
flag_disp=1; nf=0;

%%%%%%%%;
% Define spatial grid. ;
%%%%%%%%;
half_diameter_x_c = 1.0d0;
x_p_r_max = half_diameter_x_c;
diameter_x_c = 2*half_diameter_x_c;
n_x_M_u = 64;
x_c_0_ = -half_diameter_x_c + transpose([0:n_x_M_u-1]/n_x_M_u)*diameter_x_c;
x_c_1_ = -half_diameter_x_c + transpose([0:n_x_M_u-1]/n_x_M_u)*diameter_x_c;
x_c_0_ = transpose(linspace(-x_p_r_max,+x_p_r_max,n_x_M_u+1)); x_c_0_ = x_c_0_(1:n_x_M_u);
x_c_1_ = transpose(linspace(-x_p_r_max,+x_p_r_max,n_x_M_u+1)); x_c_1_ = x_c_1_(1:n_x_M_u);
%x_c_0_ = transpose(linspace(-x_p_r_max,+x_p_r_max,n_x_M_u));
%x_c_1_ = transpose(linspace(-x_p_r_max,+x_p_r_max,n_x_M_u));
[x_c_0__,x_c_1__] = ndgrid(x_c_0_,x_c_1_);
dx = mean(diff(x_c_0_));
n_xx_M_u = n_x_M_u^2;
%%%%%%%%;

%%%%%%%%;
% Define polar-quadrature parameters. ;
%%%%%%%%;
k_int = ceil(sqrt(2)*n_x_M_u); %<-- highest frequency (2*pi*k_p_r_max), watch out for aliasing! ;
k_eq_d_double = 0.5; %<-- prefactor for k_eq_d. ;
n_w_int = 2.0; %<-- prefactor for n_w_max. ;
%%%%%%%%;
k_p_r_max = k_int/(2*pi); %<-- radius of largest shell. ;
k_eq_d = k_eq_d_double/(2*pi); %<-- distance between viewing-angles on largest shell. ;
str_L = 'L'; %<-- legendre-style radial discretization. ;
[ ...
 n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
] = ...
get_weight_3d_1( ...
 0*flag_verbose ...
,k_p_r_max ...
,k_eq_d ...
,str_L ...
);
%%%%%%%%;
n_w_max = n_w_int*2*(k_int+1); %<-- number of points on largest image-ring. ;
n_w_0in_ = n_w_max*ones(n_k_p_r,1); %<-- use the same angular-discretization for each image-ring. ;
[ ...
 n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_wk_ ...
,k_p_r_wk_ ...
,k_p_w_wk_ ...
,k_c_0_wk_ ...
,k_c_1_wk_ ...
] = ...
get_weight_2d_2( ...
 0*flag_verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,-1 ...
,n_w_0in_ ...
);
n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_]);
%%%%%%%%;

%%%%%%%%;
% Now test gaussian. ;
%%%%%%%%;
sigma_x_c = 0.0625;
sigma_k_p = 1/sigma_x_c;
delta_ = 0.75*[+0.1,-0.2];
M_x_c_ = 1/(sqrt(2*pi)*sigma_x_c)^2 * exp( -( (x_c_0__-delta_(1+0)).^2 + (x_c_1__-delta_(1+1)).^2 ) / (2*sigma_x_c^2) );
M_x_c_l2 = sum(M_x_c_.^2,'all')*dx^2;
disp(sprintf(' %% sum(M_x_c_*dx^2,''all'') = %0.16f',sum(M_x_c_*dx^2,'all')));
disp(sprintf(' %% M_x_c_l2 = %0.16f',M_x_c_l2));
M_k_p_ = interp_x_c_to_k_p_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,M_x_c_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_M_u^2)*dx^2;
M_k_p_l2 = sum(abs(M_k_p_).^2 .* weight_2d_wk_) * (2*pi)^2;
disp(sprintf(' %% M_k_p_l2 = %0.16f',M_k_p_l2));
M_k_p_form_ = zeros(n_w_sum,1);
na=0;
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
n_w = n_w_(1+nk_p_r);
for nw=0:n_w-1;
k_x_c_0 = k_p_r*cos(2*pi*nw/n_w);
k_x_c_1 = k_p_r*sin(2*pi*nw/n_w);
M_k_p_form_(1+na) = exp( -( (2*pi*k_x_c_0).^2 + (2*pi*k_x_c_1).^2 ) / (2/sigma_x_c^2) ) .* exp( - 2*pi*i*( k_x_c_0*delta_(1+0) + k_x_c_1*delta_(1+1) ) );
na=na+1;
end;%for nw=0:n_w-1;
end;%for nk_p_r=0:n_k_p_r-1;
M_k_p_form_l2 = sum(abs(M_k_p_form_).^2 .* weight_2d_wk_) * (2*pi)^2;
disp(sprintf(' %% M_k_p_form_l2 = %0.16f',M_k_p_form_l2));
disp(sprintf(' %% M_k_p_ vs M_k_p_form: %0.16f',fnorm(M_k_p_ - M_k_p_form_)/fnorm(M_k_p_)));
M_x_c_reco_ = interp_k_p_to_x_c_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,M_k_p_.*weight_2d_wk_*(2*pi)^2)*sqrt(n_x_M_u^2) * n_w_sum;
M_x_c_reco_l2 = sum(abs(M_x_c_reco_).^2,'all')*dx^2;
disp(sprintf(' %% M_x_c_reco_l2 = %0.16f',M_x_c_reco_l2));
disp(sprintf(' %% M_x_c_ vs M_x_c_reco: %0.16f',fnorm(M_x_c_ - M_x_c_reco_)/fnorm(M_x_c_)));
%%%%%%%%;
if flag_disp>1;
figure(1+nf);nf=nf+1;figbig;fig80s;
p_row = 2; p_col = 3; np=0;
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,M_x_c_);axis image;axisnotick; title('M_x_c_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(M_k_p_));axis image;axisnotick; title('M_k_p_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(M_k_p_form_));axis image;axisnotick; title('M_k_p_form_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,real(M_x_c_reco_));axis image;axisnotick; title('M_x_c_reco_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,real(M_x_c_reco_-M_x_c_));axis image;axisnotick; title('difference','Interpreter','none');
end;%if flag_disp>1;
%%%%%%%%;

%%%%%%%%;
% Now test using an image with low- to moderate-frequency-content. ;
%%%%%%%%;
sigma_x_c = 0.0625/2;
sigma_k_p = 1/sigma_x_c;
g__ = @(delta_) 1/(sqrt(2*pi)*sigma_x_c)^2 * exp( -( (x_c_0__-delta_(1+0)).^2 + (x_c_1__-delta_(1+1)).^2 ) / (2*sigma_x_c^2) );
n_source = 128;
M_x_c_form_ = zeros(n_x_M_u,n_x_M_u);
for nsource=0:n_source-1;
tmp_t = sqrt(nsource/max(1,n_source-1));
tmp_r = x_p_r_max/2*tmp_t;
tmp_w = 2*pi*2*tmp_t;
M_x_c_form_ = M_x_c_form_ + g__(tmp_r*[cos(tmp_w),sin(tmp_w)]);
end;%for nsource=0:n_source-1;
M_x_c_form_l2 = sum(M_x_c_form_.^2,'all')*dx^2;
M_k_p_quad_ = interp_x_c_to_k_p_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,M_x_c_form_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_M_u^2)*dx^2;
M_k_p_quad_l2 = sum(abs(M_k_p_quad_).^2 .* weight_2d_wk_) * (2*pi)^2;
M_x_c_reco_ = interp_k_p_to_x_c_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,M_k_p_quad_.*weight_2d_wk_*(2*pi)^2)*sqrt(n_x_M_u^2) * n_w_sum;
M_x_c_reco_l2 = sum(abs(M_x_c_reco_).^2,'all')*dx^2;
M_k_p_reco_ = interp_x_c_to_k_p_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,M_x_c_reco_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_M_u^2)*dx^2;
M_k_p_reco_l2 = sum(abs(M_k_p_reco_).^2 .* weight_2d_wk_) * (2*pi)^2;
fnorm_disp(flag_verbose,'M_x_c_form_',M_x_c_form_,'M_x_c_reco_',M_x_c_reco_);
fnorm_disp(flag_verbose,'M_k_p_quad_',M_k_p_quad_,'M_k_p_reco_',M_k_p_reco_);
M_x_c_diff_ = M_x_c_form_ - M_x_c_reco_;
M_k_p_diff_ = M_k_p_quad_ - M_k_p_reco_;
%%%%%%%%;
if flag_disp>0;
figure(1+nf);nf=nf+1;figbig;fig80s;
M_x_c_lim_ = reshape(prctile(real(M_x_c_form_),[  0,100],'all'),[1,2]);
M_k_p_lim_ = reshape(prctile(real(M_k_p_quad_),[  0,100],'all'),[1,2]);
p_row = 2; p_col = 3; np=0;
subplot_t(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,real(M_x_c_form_),M_x_c_lim_);axis image;axisnotick; title('M_x_c_form_','Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(M_k_p_quad_),M_k_p_lim_);axis image;axisnotick; title('M_k_p_quad_','Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,real(M_x_c_reco_),M_x_c_lim_);axis image;axisnotick; title('M_x_c_reco_','Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(M_k_p_reco_),M_k_p_lim_);axis image;axisnotick; title('M_k_p_reco_','Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,real(M_x_c_diff_),M_x_c_lim_);axis image;axisnotick; title('M_x_c_diff_','Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(M_k_p_diff_),M_k_p_lim_);axis image;axisnotick; title('M_k_p_diff_','Interpreter','none');
end;%if flag_disp>0;
%%%%%%%%;

%%%%%%%%;
% Now test using an image with high-frequency-content. ;
%%%%%%%%;
M_x_c_form_ = mod(floor(reshape([0:n_x_M_u.^2-1],[n_x_M_u,n_x_M_u]).^(0.5)),13);
M_x_c_form_l2 = sum(M_x_c_form_.^2,'all')*dx^2;
M_k_p_quad_ = interp_x_c_to_k_p_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,M_x_c_form_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_M_u^2)*dx^2;
M_k_p_quad_l2 = sum(abs(M_k_p_quad_).^2 .* weight_2d_wk_) * (2*pi)^2;
M_x_c_reco_ = interp_k_p_to_x_c_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,M_k_p_quad_.*weight_2d_wk_*(2*pi)^2)*sqrt(n_x_M_u^2) * n_w_sum;
M_x_c_reco_l2 = sum(abs(M_x_c_reco_).^2,'all')*dx^2;
M_k_p_reco_ = interp_x_c_to_k_p_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,M_x_c_reco_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_M_u^2)*dx^2;
M_k_p_reco_l2 = sum(abs(M_k_p_reco_).^2 .* weight_2d_wk_) * (2*pi)^2;
fnorm_disp(flag_verbose,'M_x_c_form_',M_x_c_form_,'M_x_c_reco_',M_x_c_reco_);
fnorm_disp(flag_verbose,'M_k_p_quad_',M_k_p_quad_,'M_k_p_reco_',M_k_p_reco_);
M_x_c_diff_ = M_x_c_form_ - M_x_c_reco_;
M_k_p_diff_ = M_k_p_quad_ - M_k_p_reco_;
%%%%%%%%;
if flag_disp>0;
figure(1+nf);nf=nf+1;figbig;fig80s;
M_x_c_lim_ = reshape(prctile(real(M_x_c_form_),[  0,100],'all'),[1,2]);
M_k_p_lim_ = reshape(prctile(real(M_k_p_quad_),[  0,100],'all'),[1,2]);
p_row = 2; p_col = 3; np=0;
subplot_t(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,real(M_x_c_form_),M_x_c_lim_);axis image;axisnotick; title('M_x_c_form_','Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(M_k_p_quad_),M_k_p_lim_);axis image;axisnotick; title('M_k_p_quad_','Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,real(M_x_c_reco_),M_x_c_lim_);axis image;axisnotick; title('M_x_c_reco_','Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(M_k_p_reco_),M_k_p_lim_);axis image;axisnotick; title('M_k_p_reco_','Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,real(M_x_c_diff_),M_x_c_lim_);axis image;axisnotick; title('M_x_c_diff_','Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(M_k_p_diff_),M_k_p_lim_);axis image;axisnotick; title('M_k_p_diff_','Interpreter','none');
end;%if flag_disp>0;
%%%%%%%%;








