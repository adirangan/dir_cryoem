%%%%%%%%;
% visualizing the wigner-d matrix described in: ;
% Potts, Prestin and Vollrath 2009. ;
%%%%%%%%;

%%%%%%%%;
l_max = 64;
n_polar_a = 2*(l_max+1);
polar_a_ = (1*pi)*[0:n_polar_a-1]/max(1,n_polar_a); polar_a_ = transpose(polar_a_);
if ~exist('V_lmm___','var'); V_lmm___ = []; end;
if ~exist('L_lm__','var'); L_lm__ = []; end;
%%%%%%%%;
d0W_al__ = cell(n_polar_a,1+l_max);
for npolar_a=0:n_polar_a-1;
if mod(npolar_a,8)==0; disp(sprintf(' %% npolar_a %d/%d',npolar_a,n_polar_a)); end;
polar_a = polar_a_(1+npolar_a);
[ ...
 d0W_ ...
,V_lmm___ ...
,L_lm__ ...
] = ...
wignerd_c( ...
 l_max ...
,polar_a ...
,V_lmm___ ...
,L_lm__ ...
) ;
for l_val=0:l_max;
d0W_al__{1+npolar_a,1+l_val} = d0W_{1+l_val};
end;%for l_val=0:l_max;
clear d0W_;
end;%for npolar_a=0:n_polar_a-1;
%%%%%%%%;

flag_check=1;
if flag_check;
%%%%%%%%;
% Check formula 3.6. ;
%%%%%%%%;
n_test = 16;
for ntest=0:n_test-1;
rng(ntest);
l_val = max(0,min(l_max,round(ntest*l_max/max(1,n_test))));
m_val = max(-l_val,min(+l_val,-l_val + round((2*l_val+1)*rand())));
n_val = max(-l_val,min(+l_val,-l_val + round((2*l_val+1)*rand())));
l = l_val; m = m_val; n = n_val;
d = sqrt( ((l+1)^2 - m^2) * ((l+1)^2 - n^2) ) ;
u = (l+1)*(2*l+1) / d ;
v = -m*n*(2*l+1) / l / d ;
w = -(l+1)*sqrt((l^2-m^2)*(l^2-n^2)) / l / d ;
npolar_a = max(0,min(n_polar_a-1,floor(n_polar_a*rand())));
polar_a = polar_a_(1+npolar_a);
x = cos(polar_a);
lp0_val = l_val + 0;
lp1_val = l_val + 1;
ln1_val = l_val - 1;
lhs0 = 0;
if (abs(m_val)<=lp1_val) & (abs(n_val)<=lp1_val);
lhs0 = d0W_al__{1+npolar_a,1+lp1_val}(1+lp1_val+m_val,1+lp1_val+n_val);
end;%if (abs(m_val)<=lp1_val) & (abs(n_val)<=lp1_val);
rhs1 = 0;
if (abs(m_val)<=lp0_val) & (abs(n_val)<=lp0_val);
rhs1 = d0W_al__{1+npolar_a,1+lp0_val}(1+lp0_val+m_val,1+lp0_val+n_val);
end;%if (abs(m_val)<=lp0_val) & (abs(n_val)<=lp0_val);
rhs2 = 0;
if (abs(m_val)<=ln1_val) & (abs(n_val)<=ln1_val);
rhs2 = d0W_al__{1+npolar_a,1+ln1_val}(1+ln1_val+m_val,1+ln1_val+n_val);
end;%if (abs(m_val)<=ln1_val) & (abs(n_val)<=ln1_val);
rhs3 = 0;
rhs3 = (u*x + v)*rhs1 + w*rhs2;
disp(sprintf(' %% ntest %.2d/%.2d: (l,m,n) (%+.2d,%+.2d,%+.2d) npolar_a %+.3d: lhs0 %+0.6f vs rhs3 %+0.6f: %+0.16f',ntest,n_test,l,m,n,npolar_a,lhs0,rhs3,fnorm(lhs0-rhs3)));
end;%for ntest=0:n_test-1;
%%%%%%%%;
end;%if flag_check;
%%%%%%%%;

%%%%%%%%;
% Now construct tilde_D from 3.8. ;
%%%%%%%%;
m_val = +8;
n_val = -12;
l_min = max(abs(m_val),abs(n_val));
D__ = zeros(1+l_max,1+l_max-l_min);
for l_dif=0:l_max-l_min;
l_val = l_min+l_dif;
for nk=0:l_max;
npolar_a = (2*nk+1);
D__(1+nk,1+l_dif) = real(d0W_al__{1+npolar_a,1+l_val}(1+l_val+m_val,1+l_val+n_val));
end;%for nk=0:l_max;
end;%for l_dif=0:l_max-l_min;
%%%%;
% try to represent as low-rank. ;
%%%%;
[U__,S__,V__] = svds(D__(1:32,:),32); S_ = diag(S__);
plot(U__(:,1));
