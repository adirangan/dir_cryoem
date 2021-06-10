function wigners_4(n_l,kdelta_max,eps__in);
% constructs approximate factorization of wigner-s matrix, defined as: ;
% Ws = Wd(-omega) * Wt(k*delta) * Wd(+omega), where: ; 
% Wd = wigner-d matrix for angle omega ;
% Wt = wigner-t matrix for z-translation delta on shell-k ;

na=1;
if (nargin<na); n_l = 10; end; na = na+1;
if (nargin<na); kdelta_max = 2.0; end; na = na+1;
if (nargin<na); eps__in = 1e-3; end; na = na+1;

n_l = 25; kdelta_max = 2.0; eps__in = 1e-3; verbose=1;
if verbose; disp(sprintf(' %% Using n_l %d, kdelta_max %0.2f eps__in %0.6f',n_l,kdelta_max,eps__in)); end;

k_val_max = n_l;
k_val_m = k_val_max/2; k_val_r = k_val_m; % k_val in [0,k_val_max] ;
delta_m = kdelta_max/k_val_max/2; delta_r = delta_m; % delta in [0,kdelta_max]/k_val_max ;
omega_m = pi/2; omega_r = omega_m; % omega in [0,pi] ;

n_k_val = 8; n_delta = 9; n_omega = 11; n_total = n_k_val*n_delta*n_omega;
k_val_w_ = [1 2 1]; % k_val.^2 ;
delta_w_ = [1 1]; % delta.^1 ;
omega_w_ = [1]; % omega.^0 ;

[k_val_lx,k_val_lw,k_val_Lx,k_val_Lv] = orthopoly_node_weight_matrix_0(n_k_val,k_val_w_);
k_val_lt = k_val_lx*k_val_r + k_val_m; k_val_ = k_val_lt;
[delta_lx,delta_lw,delta_Lx,delta_Lv] = orthopoly_node_weight_matrix_0(n_delta,delta_w_);
delta_lt = delta_lx*delta_r + delta_m; delta_ = delta_lt;
[omega_lx,omega_lw,omega_Lx,omega_Lv] = orthopoly_node_weight_matrix_0(n_omega,omega_w_);
omega_lt = omega_lx*omega_r + omega_m; omega_ = omega_lt;

[k_val_lt_,delta_lt_,omega_lt_] = meshgrid(k_val_lt,delta_lt,omega_lt);

n_sample = 0.5;
n_lm = (1+n_l).^2;
n_m = 1+2*n_l;

[m1_,l1_,m2_,l2_,pp_,qq_] = permute_ml_to_lm(n_l);
m3_ = unique(m2_,'stable');

if verbose; disp(sprintf(' %% calculating Wd_')); tic; end;
Wdf__ = cell(n_omega,1);
Wdf_ = cell(n_omega,1);
Sdf_ = cell(n_omega,1);
Wdb__ = cell(n_omega,1);
Wdb_ = cell(n_omega,1);
Sdb_ = cell(n_omega,1);
for nomega = 1:n_omega;
omega = omega_(nomega);
Wdf__{nomega} = wignerd_b(n_l,+omega);
Wdb__{nomega} = wignerd_b(n_l,-omega);
Wdf_{nomega} = zeros(n_lm,n_lm);
Wdb_{nomega} = zeros(n_lm,n_lm);
nlm=0;
for nl=0:n_l;
l_val = nl;
nm = 1 + 2*l_val;
Wdf_{nomega}(nlm + (1:nm),nlm + (1:nm)) = Wdf__{nomega}{1+nl};
Wdb_{nomega}(nlm + (1:nm),nlm + (1:nm)) = Wdb__{nomega}{1+nl};
nlm = nlm + nm;
end;%for nl=0:n_l;
assert(nlm==n_lm);
Sdf_{nomega} = sparse(Wdf_{nomega});
Sdb_{nomega} = sparse(Wdb_{nomega});
end;%for nomega = 1:n_omega;
if verbose; disp(sprintf(' %% finished Wd_, total time %0.2f',toc)); end;
clear Wdf__ Wdf_ Wdb__ Wdb_ ;

if verbose; disp(sprintf(' %% calculating Wz_')); tic; end;
n_kdelta = n_k_val*n_delta;
kdelta_ = reshape(k_val_*transpose(delta_),n_kdelta,1); 
Wz__ = wignerz_leg(n_l,kdelta_,n_sample);
Wz_ = zeros(n_lm,n_lm);
Sz_ = cell(n_k_val,n_delta);
for nk_val = 1:n_k_val; for ndelta = 1:n_delta;
nkdelta = 1 + (nk_val-1) + (ndelta-1)*n_k_val;
nlm=0; Wz_ = zeros(n_lm,n_lm);
for nm=1:n_m;
%if nm==1; m_val=0; end; if nm>1; m_val = (1+floor((nm-2)/2))*((-1)^(mod(nm-1,2))); end; 
m_val = m3_(nm); m_abs = abs(m_val); nl = 1+n_l-m_abs; 
assert(nl==length(Wz__{nm,nkdelta}));
Wz_(nlm+(1:nl),nlm+(1:nl)) = Wz__{nm,nkdelta};
nlm = nlm+nl;
end;%for nl=0:n_l;
assert(nlm==n_lm);
Sz_{nk_val,ndelta} = sparse(Wz_);
end;end;%for nk_val = 1:n_k_val; for ndelta = 1:n_delta;
if verbose; disp(sprintf(' %% finished Wz_, total time %0.2f',toc)); end;
clear Wz__ Wz_;

if verbose; disp(sprintf(' %% calculating Bandlimits')); tic; end;
nomega = round(n_omega/2); ndelta = n_delta; nk_val = n_k_val;
Ss_tmp = Sdb_{nomega}*Sz_{nk_val,ndelta}(qq_,qq_)*Sdf_{nomega};
Bandlimit_M = 6; Bandlimit_L = 6;
BB = hypercube_to_bandlimits_0(Ss_tmp,Bandlimit_M,Bandlimit_L);
if verbose; disp(sprintf(' %% finished Bandlimits, total time %0.2f',toc)); end;
imagesc(log10(BB),[-5,0]); colorbar; 
xlabel('Bandlimit L'); ylabel('Bandlimit M'); title('log10(error)'); 

Bandlimit_M = 3; Bandlimit_L = 3;
if verbose; disp(sprintf(' %% Setting Bandlimit_M %d Bandlimit_L %d',Bandlimit_M,Bandlimit_L)); end;
if verbose; disp(sprintf(' %% calculating Ws_')); tic; end;
[n_p,p1_,p2_,~] = hypercube_to_diagdiag_0(zeros(n_lm,n_lm),Bandlimit_M,Bandlimit_L);
p3_ = p1_ + p2_*n_lm;
Ss_ = zeros(n_k_val,n_delta*n_omega,n_p);
for nomega = 1:n_omega; for ndelta = 1:n_delta; for nk_val = 1:n_k_val; 
ntotal = 1 + (nk_val-1) + (ndelta-1)*n_k_val + (nomega-1)*n_k_val*n_delta;
if (mod(ntotal,100)==0); disp(sprintf(' %% ntotal %d/%d',ntotal,n_total)); end;
Ss_tmp = Sdb_{nomega}*Sz_{nk_val,ndelta}(qq_,qq_)*Sdf_{nomega};
Ss_(1 + (nk_val-1),1 + (ndelta-1) + (nomega-1)*n_delta,:) = Ss_tmp(1+p3_);
end;end;end;%for nomega = 1:n_omega; for ndelta = 1:n_delta; for nk_val = 1:n_k_val; 
if verbose; disp(sprintf(' %% finished Ws_, total time %0.2f',toc)); end;

if verbose; disp(sprintf(' %% calculating Sc_')); tic; end;
Sc_ = zeros(n_p*n_k_val,n_delta*n_omega);
ntotal=0;
for nomega = 0:n_omega-1; for ndelta = 0:n_delta-1; 
LW2_tmp = reshape(transpose(delta_Lx(1+ndelta,:).*transpose(delta_lw))*(omega_Lx(1+nomega,:).*transpose(omega_lw)),n_delta*n_omega,1);
for nk_val = 0:n_k_val-1; 
ntotal = ntotal + 1; if (mod(ntotal,100)==0); disp(sprintf(' %% ntotal %d/%d',ntotal,n_total)); end;
LW1_tmp = k_val_Lx(1+nk_val,:).*transpose(k_val_lw);
for np=1:n_p;
Sc_(np + nk_val*n_p,1 + ndelta + nomega*n_delta) = LW1_tmp*squeeze(Ss_(:,:,np))*LW2_tmp;
end;%for np=1:n_p;
end;end;end;%for nomega = 0:n_omega-1; for ndelta = 0:n_delta-1; for nk_val = 0:n_k_val-1; 
if verbose; disp(sprintf(' %% finished Sc_, total time %0.2f',toc)); end;

np = 234;
if verbose; disp(sprintf(' %% calculating Sc_sub for np==%d',np)); tic; end;
Ss_sub = Ss_(:,:,np);
Sc_sub = zeros(n_k_val,n_delta*n_omega);
ntotal=0;
for nomega = 0:n_omega-1; for ndelta = 0:n_delta-1; 
LW2_tmp = reshape(transpose(delta_Lx(1+ndelta,:).*transpose(delta_lw))*(omega_Lx(1+nomega,:).*transpose(omega_lw)),n_delta*n_omega,1);
for nk_val = 0:n_k_val-1; 
ntotal = ntotal + 1; if (mod(ntotal,100)==0); disp(sprintf(' %% ntotal %d/%d',ntotal,n_total)); end;
LW1_tmp = k_val_Lx(1+nk_val,:).*transpose(k_val_lw);
Sc_sub(1 + nk_val,1 + ndelta + nomega*n_delta) = LW1_tmp*squeeze(Ss_(:,:,np))*LW2_tmp;
end;end;end;%for nomega = 0:n_omega-1; for ndelta = 0:n_delta-1; for nk_val = 0:n_k_val-1; 
if verbose; disp(sprintf(' %% finished Sc_sub, total time %0.2f',toc)); end;
disp(sprintf(' %% np %d: difference Sc_ vs Sc_sub: %0.16f',np,norm(Sc_(np + (0:n_k_val-1)*n_p,:) - Sc_sub)));

if verbose; disp(sprintf(' reconstructing Ss_sub from Sc_sub; np %d',np)); end;
nomega_rec = floor(n_omega/2); ndelta_rec = floor(n_delta*2/3); nk_val_rec = floor(n_k_val*3/4);
omega_rec = omega_(1+nomega_rec); delta_rec = delta_(1+ndelta_rec); k_val_rec = k_val_(1+nk_val_rec);
disp(sprintf(' %% omega_(1+%d)=%0.2f, delta_(1+%d)=%0.3f, k_val_(1+%d)=%0.1f',nomega_rec,omega_rec,ndelta_rec,delta_rec,nk_val_rec,k_val_rec));
Ss_orig = Ss_(1 + nk_val_rec,1 + ndelta_rec + nomega_rec*n_delta,np);
Sc_tmp1 = 0;
for nk_val=0:n_k_val-1;
k_val_Lv_tmp = polyval(k_val_Lv(1+nk_val,:),(k_val_rec - k_val_m)/k_val_r);
for ndelta=0:n_delta-1; 
delta_Lv_tmp = polyval(delta_Lv(1+ndelta,:),(delta_rec - delta_m)/delta_r);
for nomega=0:n_omega-1;
omega_Lv_tmp = polyval(omega_Lv(1+nomega,:),(omega_rec - omega_m)/omega_r);
Sc_tmp1 = Sc_tmp1 + Sc_sub(1 + nk_val,1 + ndelta + nomega*n_delta)*k_val_Lv_tmp*delta_Lv_tmp*omega_Lv_tmp;
end;%for nomega=0:n_omega-1;
end;%for ndelta=0:n_delta-1;
end;%for nk_val=0:n_k_val-1;
disp(sprintf(' %% np = %d: Ss_orig (%0.16f,%0.16f), Sc_tmp1 (%0.16f,%0.16f)',np,real(Ss_orig),imag(Ss_orig),real(Sc_tmp1),imag(Sc_tmp1)));

if verbose; disp(sprintf(' reconstructing Ss_ from Sc_')); end;
nomega_rec = floor(n_omega/2); ndelta_rec = floor(n_delta*2/3); nk_val_rec = floor(n_k_val*3/4);
omega_rec = omega_(1+nomega_rec); delta_rec = delta_(1+ndelta_rec); k_val_rec = k_val_(1+nk_val_rec);
disp(sprintf(' %% omega_(1+%d)=%0.2f, delta_(1+%d)=%0.3f, k_val_(1+%d)=%0.1f',nomega_rec,omega_rec,ndelta_rec,delta_rec,nk_val_rec,k_val_rec));
Ss_orig = Ss_(1 + nk_val_rec,1 + ndelta_rec + nomega_rec*n_delta,np);
iteration_max = 16;
for iteration=1:iteration_max;
%np=max(1,min(n_p,floor(n_p*rand())));
np = iteration;
Ss_orig = Ss_(1 + nk_val_rec,1 + ndelta_rec + nomega_rec*n_delta,np);
Sc_tmp1 = 0;
for nk_val=0:n_k_val-1;
k_val_Lv_tmp = polyval(k_val_Lv(1+nk_val,:),(k_val_rec - k_val_m)/k_val_r);
for ndelta=0:n_delta-1; 
delta_Lv_tmp = polyval(delta_Lv(1+ndelta,:),(delta_rec - delta_m)/delta_r);
for nomega=0:n_omega-1;
omega_Lv_tmp = polyval(omega_Lv(1+nomega,:),(omega_rec - omega_m)/omega_r);
Sc_tmp1 = Sc_tmp1 + Sc_(np + nk_val*n_p,1 + ndelta + nomega*n_delta)*k_val_Lv_tmp*delta_Lv_tmp*omega_Lv_tmp;
end;%for nomega=0:n_omega-1;
end;%for ndelta=0:n_delta-1;
end;%for nk_val=0:n_k_val-1;
disp(sprintf(' %% np = %d: Ss_orig (%0.16f,%0.16f), Sc_tmp1 (%0.16f,%0.16f)',np,real(Ss_orig),imag(Ss_orig),real(Sc_tmp1),imag(Sc_tmp1)));
end;%for iteration=1:iteration_max;

if verbose; disp(sprintf(' %% calculating svd of Sc_')); tic; end;
[Sc_U_,Sc_S_,Sc_V_] = svds(Sc_,10);
if verbose; disp(sprintf(' %% finished svd of Sc_, total time %0.2f',toc)); end;

plot(log10(diag(Sc_S_)/Sc_S_(1,1)),'.-','MarkerSize',25); 
xlabel('term number'); ylabel('log10(sigma)'); title('svd of wigner-s');

n_term = 8;
disp(sprintf(' %% using n_term = %d for reconstruction',n_term));
nomega_rec = floor(n_omega/2); ndelta_rec = floor(n_delta*2/3); nk_val_rec = floor(n_k_val*3/4);
omega_rec = omega_(1+nomega_rec); delta_rec = delta_(1+ndelta_rec); k_val_rec = k_val_(1+nk_val_rec);
disp(sprintf(' %% omega_(1+%d)=%0.2f, delta_(1+%d)=%0.3f, k_val_(1+%d)=%0.1f',nomega_rec,omega_rec,ndelta_rec,delta_rec,nk_val_rec,k_val_rec));
Ss_orig = Sdb_{1+nomega_rec}*Sz_{1+nk_val_rec,1+ndelta_rec}(qq_,qq_)*Sdf_{1+nomega_rec};
Ss_tmp1 = sparse(1+p1_,1+p2_,squeeze(Ss_(1 + nk_val_rec,1 + ndelta_rec + nomega_rec*n_delta,:)),n_lm,n_lm);
disp(sprintf(' %% Bandlimits (%d,%d), error %0.6f',Bandlimit_M,Bandlimit_L,norm(full(Ss_orig-Ss_tmp1))/norm(full(Ss_orig))));

Sc_tmp2 = zeros(n_p,1);
for nk_val=0:n_k_val-1;
k_val_Lv_tmp = polyval(k_val_Lv(1+nk_val,:),(k_val_rec - k_val_m)/k_val_r);
for ndelta=0:n_delta-1; 
delta_Lv_tmp = polyval(delta_Lv(1+ndelta,:),(delta_rec - delta_m)/delta_r);
for nomega=0:n_omega-1;
omega_Lv_tmp = polyval(omega_Lv(1+nomega,:),(omega_rec - omega_m)/omega_r);
Sc_tmp2 = Sc_tmp2 + squeeze(Sc_((1:n_p) + nk_val*n_p,1 + ndelta + nomega*n_delta))*k_val_Lv_tmp*delta_Lv_tmp*omega_Lv_tmp;
end;%for nomega=0:n_omega-1;
end;%for ndelta=0:n_delta-1;
end;%for nk_val=0:n_k_val-1;
Sc_tmp3 = sparse(1+p1_,1+p2_,squeeze(Sc_tmp2),n_lm,n_lm);
disp(sprintf(' %% interpolation error %0.6f',norm(full(Ss_tmp1-Sc_tmp3))/norm(full(Ss_tmp1))));
Sc_tmp4 = zeros(n_lm,n_lm);
for nterm=1:n_term;
Sc_U_tmp = zeros(n_p,1);
for nk_val=0:n_k_val-1;
k_val_Lv_tmp = polyval(k_val_Lv(1+nk_val,:),(k_val_rec - k_val_m)/k_val_r);
Sc_U_tmp = Sc_U_tmp + Sc_U_((1:n_p) + nk_val*n_p,nterm)*k_val_Lv_tmp;
end;%for nk_val=0:n_k_val-1;
Sc_V_tmp = 0;
for ndelta=0:n_delta-1; 
delta_Lv_tmp = polyval(delta_Lv(1+ndelta,:),(delta_rec - delta_m)/delta_r);
for nomega=0:n_omega-1;
omega_Lv_tmp = polyval(omega_Lv(1+nomega,:),(omega_rec - omega_m)/omega_r);
Sc_V_tmp = Sc_V_tmp + Sc_V_(1 + ndelta + nomega*n_delta,nterm)*delta_Lv_tmp*omega_Lv_tmp;
end;end;%for ndelta=0:n_delta-1; for nomega=0:n_omega-1;
Sc_tmp4 = Sc_tmp4 + sparse(1+p1_,1+p2_,Sc_U_tmp,n_lm,n_lm)*Sc_S_(nterm,nterm)*Sc_V_tmp;
disp(sprintf(' %% nterm %d: reconstruction error %0.6f',nterm,norm(full(Sc_tmp3-Sc_tmp4))/norm(full(Sc_tmp3))));
end;%for nterm=1:n_term;

disp(sprintf(' %% n_term %d: full reconstruction error %0.6f',n_term,norm(full(Ss_orig-Sc_tmp4))/norm(full(Ss_orig))));


