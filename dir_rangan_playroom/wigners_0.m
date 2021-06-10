function wigners_0(n_l,kdelta_max,eps);
% constructs approximate factorization of wigner-s matrix, defined as: ;
% Ws = Wd(-omega) * Wt(k*delta) * Wd(+omega), where: ; 
% Wd = wigner-d matrix with angle omega ;
% Wt = wigner-t matrix with translation delta on shell-k ;

na=1;
if (nargin<na); n_l = 10; end; na = na+1;
if (nargin<na); kdelta_max = 2.0; end; na = na+1;
if (nargin<na); eps = 1e-3; end; na = na+1;

k_max = n_l;
n_kdelta = max(20,2*n_l);
k_ = linspace(1,k_max,n_kdelta);
delta_ = linspace(0,kdelta_max/k_max,n_kdelta+1);
omega_ = linspace(0,pi,n_kdelta+2);
n_all = length(k_)*length(delta_)*length(omega_);

n_sample = 3;
n_lm = (1+n_l).^2;

W_all_ = zeros(length(k_),length(delta_),length(omega_),n_lm,n_lm);

na=0;
for nk = 1:length(k_);
k_val = k_(nk);
for ndelta = 1:length(delta_);
delta = delta_(ndelta);
for nomega = 1:length(omega_);
omega = omega_(nomega);
if (mod(na,10)==0); disp(sprintf(' %% na %d/%d',na,n_all)); end;

[Wt_,pp_,m2_,l2_] = wignert_leg(n_l,[0,0,1]*k_val*delta,n_sample);
Wdf_sub_ = wignerd_b(n_l,+omega);
Wdb_sub_ = wignerd_b(n_l,-omega);
Wdf_ = zeros(n_lm,n_lm);
Wdb_ = zeros(n_lm,n_lm);
nlm=0;
for nl=0:n_l;
l_val = nl;
nm = 1 + 2*l_val;
Wdf_(nlm + (1:nm),nlm + (1:nm)) = Wdf_sub_{1+nl};
Wdb_(nlm + (1:nm),nlm + (1:nm)) = Wdb_sub_{1+nl};
nlm = nlm + nm;
end;%for nl=0:n_l;
assert(nlm==n_lm);
Ws_ = Wdb_*Wt_*Wdf_;

disp_flag=0;
if disp_flag;
figure;
subplot(1,2,1); imagesc(log10(abs(Ws_)),[log10(eps),0]);colorbar; title('m,l');
subplot(1,2,2); imagesc(log10(abs(Ws_(1+pp_,1+pp_))),[log10(eps),0]);colorbar; title('l,m');
disp(sprintf(' %% Ws_ of size %d-x-%d, %d/%d = %0.2f nonzeros at eps=%0.6f',n_lm,n_lm,nnz(abs(Ws_)>eps),n_lm*n_lm,nnz(abs(Ws_)>eps)/(n_lm*n_lm),eps));
end;%if disp_flag;

W_all_(nk,ndelta,nomega,:,:) = Ws_;
na = na+1;

end;%for nomega = 1:length(omega_);
end;%for ndelta = 1:length(delta_);
end;%for nk = 1:length(k_);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

fname = sprintf('wigners_l%d_kd%.2d.mat',n_l,round(10*kdelta_max));
save(fname,'-v7.3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

B_all_ = permute(W_all_,[4,5,1,2,3]);
B_tmp_ = reshape(B_all_,[size(B_all_,1)*size(B_all_,2)*size(B_all_,3),size(B_all_,4)*size(B_all_,5)]);

n_s = 32;
[U_tmp_,S_tmp_,V_tmp_] = svds(B_tmp_,n_s);
U_all_ = reshape(U_tmp_,[size(B_all_,1),size(B_all_,2),size(B_all_,3),n_s]);
S_all_ = diag(S_tmp_);
V_all_ = reshape(V_tmp_,[size(B_all_,4),size(B_all_,5),n_s]);

for ns=1:n_s;
figure;
prows=4;pcols=5;
for np=1:20;
subplot(prows,pcols,np);
imagesc(abs(U_all_(:,:,np,ns))); set(gca,'XTick',[],'YTick',[]); title(sprintf('s%d p%d',ns,np));
end;%for np=1:20;
end;%for ns=1:n_s;

figure; plot(log10(S_all_),'o-');

figure;
for ns=1:min(20,n_s);
prows=4;pcols=5;
subplot(prows,pcols,ns);
imagesc(abs(V_all_(:,:,ns))); set(gca,'XTick',[],'YTick',[]); title(sprintf('s%d',ns));
end;%for ns=1:min(20,n_s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;




