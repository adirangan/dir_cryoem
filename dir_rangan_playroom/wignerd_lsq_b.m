function W_ = wignerd_lsq_b(n_l,beta) ;
% generates wigner-d matrices up to n_l;
% uses least-squares solve ;
% test with: 

%%%%%%%%;
if nargin<1;
%%%%%%%%;
disp(sprintf(' %% testing wignerd_lsq_b'));
%%%%%%%%;
% Note that the discrepancy is larger than 1e-6 at nl==88. ;
%%%%%%%%;
n_l=48;
%n_l=196;
beta=pi/6; nf=0;
tic; W1_ = wignerd_b(n_l,beta);     disp(sprintf(' %% wignerd_b  : %0.2f seconds',toc));
tic; W2_ = wignerd_lsq_b(n_l,beta); disp(sprintf(' %% wignerd_lsq_b: %0.2f seconds',toc));  
tic; W3_ = wignerd_c(n_l,beta);     disp(sprintf(' %% wignerd_c  : %0.2f seconds',toc));
for nl=0:n_l;
disp(sprintf(' %% nl %.3d/%.3d: W1_ vs W2_: error %0.16f',nl,n_l,fnorm(W1_{1+nl}-W2_{1+nl})));
disp(sprintf(' %% nl %.3d/%.3d: W3_ vs W2_: error %0.16f',nl,n_l,fnorm(W3_{1+nl}-W2_{1+nl})));
end;%for nl=0:n_l;
figure(1+nf);nf=nf+1;clf;figmed;fig80s;
p_row = 1; p_col = 3; np=0;
fontsize_use = 12;
wlim_ = 0.5*[-1,+1];
subplot(p_row,p_col,1+np);np=np+1;
imagesc(real(W1_{end}),wlim_);
axis image;
axisnotick;
xlabel('m1_val','Interpreter','none');
ylabel('m0_val','Interpreter','none');
title(sprintf('W1_{%.3d}: wignerd_ori_b',n_l),'Interpreter','none');
set(gca,'FontSize',fontsize_use);
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
imagesc(real(W2_{end}),wlim_);
axis image;
axisnotick;
xlabel('m1_val','Interpreter','none');
ylabel('m0_val','Interpreter','none');
title(sprintf('W2_{%.3d}: wignerd_lsq_b',n_l),'Interpreter','none');
set(gca,'FontSize',fontsize_use);
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
imagesc(real(W3_{end}),wlim_);
axis image;
axisnotick;
xlabel('m1_val','Interpreter','none');
ylabel('m0_val','Interpreter','none');
title(sprintf('W3_{%.3d}: wignerd_lsq_b',n_l),'Interpreter','none');
set(gca,'FontSize',fontsize_use);
%%%%%%%%;
disp('returning');return;
end;%if nargin<1;
%%%%%%%%;

flag_verbose = 0;

% flag_verbose=1; n_l = 11; beta = pi/6;

n_oversample = 5;
l_max = n_l;
m_max_ = transpose(-l_max:+l_max);
n_m_max = 1+2*l_max;
n_azimu_b = ceil(sqrt(n_oversample*2*n_m_max));
azimu_b_ = sort(2*pi*rand(n_azimu_b,1));
n_polar_a = ceil(sqrt(n_oversample*1*n_m_max));
polar_a_ = sort(1*pi*rand(n_polar_a,1));
[polar_a_ori__,azimu_b_ori__] = ndgrid(polar_a_,azimu_b_);
%%%%;
flag_check=0;
if flag_check;
tmp_t = tic();
J_lmab____ = permute(reshape(ylgndr_1(l_max,cos(polar_a_ori__(:))),[n_polar_a,n_azimu_b,1+l_max,1+l_max]),[3,4,1,2]);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% J_lmab____: %0.2fs',tmp_t)); end;
end;%if flag_check;
%%%%;
tmp_t = tic();
[ ...
 d0y_jlm___ ...
,sqrt_2lp1_ ...
,sqrt_2mp1_ ...
,sqrt_rat0_ ...
,sqrt_rat3__ ...
,sqrt_rat4__ ...
] = ...
ylgndr_1( ...
 l_max ...
,cos(polar_a_) ...
);
L_ori_lmab____ = repmat(permute(reshape(d0y_jlm___,[n_polar_a,1+l_max,1+l_max]),1+[1,2,0]),[1,1,1,n_azimu_b]);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% L_ori_lmab____: %0.2fs',tmp_t)); end;
%%%%;
if flag_check;
if (flag_verbose); disp(sprintf(' %% L_ori_lmab____ vs J_lmab____: %0.16f',fnorm(L_ori_lmab____ - J_lmab____)/fnorm(L_ori_lmab____))); end;
end;%if flag_check;
%%%%;
L_ori_lmab____ = cat(2,flip(L_ori_lmab____(:,2:end,:,:),2),L_ori_lmab____)/sqrt(4*pi);
s_ = ones(n_m_max,1);
%s_ = (-1).^((m_max_<0).*m_max_); %<-- unnecessary here. ;
expi_mb__ = exp(+i*m_max_*reshape(azimu_b_,[1,n_azimu_b]));
Y_ori_lmab____ = bsxfun(@times,bsxfun(@times,L_ori_lmab____,reshape(s_,[1,n_m_max,1,1])),reshape(expi_mb__,[1,n_m_max,1,n_azimu_b]));
%%%%;
cb = cos(+beta); sb = sin(+beta); sg = -1;
Xn__ = sin(polar_a_ori__).*cos(azimu_b_ori__);
Yn__ = sin(polar_a_ori__).*sin(azimu_b_ori__);
Zn__ = cos(polar_a_ori__);
Xt__ = +cb*Xn__ + sg*sb*Zn__;
Yt__ = Yn__;
Zt__ = -sg*sb*Xn__ + cb*Zn__;
azimu_b_rot__ = atan2(Yt__,Xt__);
polar_a_rot__ = acos(Zt__);
%%%%;
tmp_t = tic();
L_rot_lmab____ = permute(reshape(ylgndr_1(l_max,cos(polar_a_rot__(:)),sqrt_2lp1_,sqrt_2mp1_,sqrt_rat0_,sqrt_rat3__,sqrt_rat4__),[n_polar_a,n_azimu_b,1+l_max,1+l_max]),1+[2,3,0,1]);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% L_rot_lmab____: %0.2fs',tmp_t)); end;
L_rot_lmab____ = cat(2,flip(L_rot_lmab____(:,2:end,:,:),2),L_rot_lmab____)/sqrt(4*pi);
s_ = ones(n_m_max,1);
%s_ = (-1).^((m_max_<0).*m_max_); %<-- unnecessary here. ;
expi_mab___ = exp(+i*bsxfun(@times,m_max_,reshape(azimu_b_rot__,[1,n_polar_a,n_azimu_b])));
Y_rot_lmab____ = bsxfun(@times,bsxfun(@times,L_rot_lmab____,reshape(s_,[1,n_m_max,1,1])),reshape(expi_mab___,[1,n_m_max,n_polar_a,n_azimu_b]));
%%%%;
W_ = cell(1+l_max,1);
tmp_t = tic();
for l_val=0:l_max;
Y_ori_mab__ = reshape(Y_ori_lmab____(1+l_val,1+l_max+[-l_val:+l_val],:,:),[1+2*l_val,n_polar_a*n_azimu_b]);
Y_rot_mab__ = reshape(Y_rot_lmab____(1+l_val,1+l_max+[-l_val:+l_val],:,:),[1+2*l_val,n_polar_a*n_azimu_b]);
W_{1+l_val} = real(transpose(Y_rot_mab__ / Y_ori_mab__)) ;
end;%for l_val=0:l_max;
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% W_: %0.2fs',tmp_t)); end;

%%%%%%%%;
flag_check=0;
%%%%%%%%;
if flag_check;
sum_t0 = 0;
sum_t1 = 0;
tmp_t0 = tic();
V_ = cell(1+n_l,1);
V_{1} = [1];
for nl=1:n_l;
n_m = 1+2*nl;
%azimu_b_ = linspace( 0 , 2*pi , 2+2*nl );  % Azimuthal/Longitude/Circumferential ;
%polar_a_   = linspace( 0 ,   pi , 2+2*nl );  % Altitude /Latitude /Elevation ;
azimu_b_ = sort(2*pi*rand(1,ceil(sqrt(n_oversample*n_m))));
polar_a_   = sort(1*pi*rand(1,ceil(sqrt(n_oversample*n_m))));
[polar_a__,azimu_b__] = ndgrid(polar_a_,azimu_b_);
Ylm_orig_ = ylm_1(nl,azimu_b__,polar_a__);
cb = cos(+beta); sb = sin(+beta); sg = -1;
Xn__ = sin(polar_a__).*cos(azimu_b__);
Yn__ = sin(polar_a__).*sin(azimu_b__);
Zn__ = cos(polar_a__);
Xt__ = +cb*Xn__ + sg*sb*Zn__;
Yt__ = Yn__;
Zt__ = -sg*sb*Xn__ + cb*Zn__;
azimu_b_y__ = atan2(Yt__,Xt__);
polar_a_y__ = acos(Zt__);
Ylm_rota_ = ylm_1(nl,azimu_b_y__,polar_a_y__);
n_x = length(azimu_b_)*length(polar_a_);
Y_orig_ = zeros(n_m,n_x);
Y_rota_ = zeros(n_m,n_x);
for nm=1:n_m;
m_val = -nl-1+nm;
s=(-1)^((m_val<0)*m_val); % needed to preserve condon-shortley phase. ;
Y_orig_(nm,:) = s*reshape(Ylm_orig_{nm},1,n_x);
Y_rota_(nm,:) = s*reshape(Ylm_rota_{nm},1,n_x);
end;%for nm=1:n_m;
tmp_t1 = tic();
V_{1+nl} = real(transpose(Y_rota_ / Y_orig_)) ;
tmp_t1 = toc(tmp_t1); sum_t1 = sum_t1 + tmp_t1;
end;% for nl=1:n_l;
tmp_t0 = toc(tmp_t0); sum_t0 = sum_t0 + tmp_t0;
%%%%%%%%;
for l_val=0:l_max;
disp(sprintf(' %% l_val %.2d/%.2d: V vs W: %0.16f',l_val,l_max,fnorm(V_{1+l_val}-W_{1+l_val})/fnorm(V_{1+l_val})));
end;%for l_val=0:l_max;
%%%%%%%%;
end;%if flag_check;
%%%%%%%%;

