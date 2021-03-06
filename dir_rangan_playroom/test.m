addpath('~/chebfun');
K_max = 2; delta = 1.2; sample_d = 0.25;
b_K = 16; [b_lx,b_lw] = legpts(b_K);
k_ = (b_lx+1.0)*K_max/2; k_w_ = b_lw*K_max/2;
Iq_ = zeros(b_K,1); Ix_ = zeros(b_K,1);
[n_all,theta_all_,phi_all_,weight_all_] = sample_shell_2(1,sample_d,'L') ;
for nk=1:b_K;
k = k_(nk);
Iq_(nk) = sum(k^2*cos(k*delta*cos(phi_all_)).*weight_all_);
Ix_(nk) = 4*pi*(1/delta)*sin(delta*k).*k;
end;% for nk=1:b_K;
Iq = k_w_*Iq_;
Ix = 4*pi*(1/delta^3)*( sin(K_max*delta) - (K_max*delta)*cos(K_max*delta) );
disp(sprintf(' %% (Ix-Iq)/Ix = %0.16f/%0.16f = %0.16f',Ix-Iq,Ix,(Ix-Iq)/Ix));
