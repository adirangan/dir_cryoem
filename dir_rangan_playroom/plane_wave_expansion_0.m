function C_ = plane_wave_expansion_0();
% tests out plane-wave expansion in terms of spherical-bessel-functions and spherical-harmonics. ;
k_ = [+1.2;+0.3;-5.6]; %<- wavenumber. ;
x_ = [-0.8;+1.3;+2.3]/3; %<- space. ;
disp(sprintf(' %% testing plane-wave expansion, 2*pi*k*x = %0.2f',2*pi*dot(k_,x_)));
k = fnorm(k_); k_01 = fnorm(k_(1:2)); azimu_b_k = atan2(k_(2),k_(1)); polar_a_k = atan2(k_01,k_(3));
x = fnorm(x_); x_01 = fnorm(x_(1:2)); azimu_b_x = atan2(x_(2),x_(1)); polar_a_x = atan2(x_01,x_(3));
Ix_ = exp(2*pi*i*dot(k_,x_));
Ip_ = 0;
l_max = 135;
Ylm_k__ = get_Ylm__(1+l_max,0:l_max,1,azimu_b_k,polar_a_k,0);
Ylm_x__ = get_Ylm__(1+l_max,0:l_max,1,azimu_b_x,polar_a_x,0);
for l_val=0:l_max;
%[Ylm_d_k,Ylm_s_k] = Ylm_d(l_val,azimu_b_k,polar_a_k); [Ylm_d_x,Ylm_s_x] = Ylm_d(l_val,azimu_b_x,polar_a_x);
Ylm_d_k = Ylm_k__{1+l_val}(:); Ylm_s_k = 1; Ylm_d_x = Ylm_x__{1+l_val}(:); Ylm_s_x = 1;
t = 2*pi*k*x;
jl = besselj(l_val+0.5,t)*sqrt(pi/(2*t));
%Ip_ = Ip_ + sum(4*pi*(i^l_val)*jl*(Ylm_d_k.*Ylm_s_k.*conj(Ylm_d_x.*Ylm_s_x)));
Ip_ = Ip_ + sum(4*pi*(i^l_val)*jl*(conj(Ylm_d_k.*Ylm_s_k).*(Ylm_d_x.*Ylm_s_x)));
disp(sprintf('k: %f, x: %f, Ix: %f+%fi, l_val %d; Ip: %f+%fi, error: %0.16f',k,x,real(Ix_),imag(Ix_),l_val,real(Ip_),imag(Ip_),fnorm(Ix_-Ip_)));
end;%for l_val=0:l_max;
equatorial_distance = 0.25; disp(sprintf(' %% testing Ylm normalization with equatorial_distance %0.2f',equatorial_distance));
[n_all,azimu_b_all_,polar_a_all_,weight_all_] = sample_shell_4_nojvm(1.0d0,equatorial_distance,'L') ;
l_max = 20; l_val_ = transpose(0:l_max); n_l = length(l_val_);
Ylm__ = get_Ylm__(n_l,l_val_,n_all,azimu_b_all_,polar_a_all_);
Ix_ = 4*pi;
disp(sprintf(' %% surface area: %0.2f, sum_weight: %0.2f, error: %0.16f',Ix_,sum(weight_all_),fnorm(Ix_-sum(weight_all_))));
n_m_ = (2*l_val_+1);
n_lm = sum(n_m_);
assert(n_lm==(l_max+1)^2);
C_ = zeros(n_lm,n_lm);
na=0;
for nla=0:n_l-1; l_val_a=l_val_(1+nla); n_ma = 2*l_val_a+1; for nma=0:n_ma-1; m_val_a = -l_val_a + nma;
nb=0;
for nlb=0:n_l-1; l_val_b=l_val_(1+nlb); n_mb = 2*l_val_b+1; for nmb=0:n_mb-1; m_val_b = -l_val_b + nmb;
tmp_a_ = Ylm__{1+nla}(1+nma,:); tmp_a_ = reshape(tmp_a_,n_all,1);
tmp_b_ = Ylm__{1+nlb}(1+nmb,:); tmp_b_ = reshape(tmp_b_,n_all,1);
C_(1+na,1+nb) = sum(conj(tmp_a_).*(tmp_b_).*weight_all_);
nb=nb+1;
end;end;%for nlb=0:n_l-1; l_val_b=l_val_(1+nlb); n_mb = 2*l_val_b+1; for nmb=0:n_mb-1; m_val_b = -l_val_b + nmb;
na=na+1;
end;end;%for nla=0:n_l-1; l_val_a=l_val_(1+nla); n_ma = 2*l_val_a+1; for nma=0:n_ma-1; m_val_a = -l_val_a + nma;
subplot(1,2,1); imagesc(real(C_));colormap(colormap_beach()); colorbar; title('real(C)');
subplot(1,2,2); imagesc(imag(C_));colormap(colormap_beach()); colorbar; title('imag(C)');
figbig;




