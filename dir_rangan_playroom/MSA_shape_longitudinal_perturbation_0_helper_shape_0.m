%%%%%%%%;
% determining shape. ;
%%%%%%%%;

%%%%%%%%;
% Sum of gaussians. ;
%%%%%%%%;
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
flag_disp=0;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figbig;
subplot(2,2,1);imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,tmp_N_x_c_);axis image;axisnotick;
subplot(2,2,2);imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(tmp_N_k_p_));axis image;axisnotick;
subplot(2,2,3);imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(tmp_N_k_p_form_));axis image;axisnotick;
subplot(2,2,4);imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(tmp_N_x_c_reco_));axis image;axisnotick;
error('stopping');
end;%if flag_disp;
%%%%%%%%;





