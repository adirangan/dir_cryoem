
addpath('/data/rangan/dir_cryoem/chebfun'); %savepath;

clear;clc;
verbose=1;
p_ = [1,1];

for K = 6:2:24;

disp(sprintf(' %% K %d: ',K));
[jx,jw,j_Lx,j_Lv] = orthopoly_node_weight_matrix_0(K,p_);
[cx,cw] = jacpts(K,0,1); cw=transpose(cw);
c_Lv = zeros(1+K,1+K);
c_Lx = zeros(K,K);
for nk=0:K;
cj = jacpoly(nk,0,1)*sqrt(nk+1)/sqrt(2);
if (nk<K); c_Lx(1+nk,:) = cj(jx); end;
c_Lv(1+nk,1+K-nk:1+K) = poly(cj);
end;%for nk=1:K;
disp(sprintf(' %% error in x: %0.16f',norm(jx-cx,'fro')/norm(cx,'fro')));
disp(sprintf(' %% error in w: %0.16f',norm(jw-cw,'fro')/norm(cw,'fro')));
disp(sprintf(' %% error in Lv: %0.16f',norm(j_Lv-c_Lv,'fro')/norm(c_Lv,'fro')));
disp(sprintf(' %% error in Lx: %0.16f',norm(j_Lx-c_Lx,'fro')/norm(c_Lx,'fro')));

J_ = zeros(K+1,K+1);
for nka=0:K; for nkb=0:K;
J_(1+nka,1+nkb) = polyint(conv(p_,conv(j_Lv(1+nka,:),j_Lv(1+nkb,:))),[-1,+1]);
end;end;%for nka=0:K; for nkb=0:K;
if (verbose>1); disp(sprintf(' %% J_: ')); disp(num2str(J_)); end;

C_ = zeros(K+1,K+1);
for nka=0:K; for nkb=0:K;
C_(1+nka,1+nkb) = polyint(conv(p_,conv(c_Lv(1+nka,:),c_Lv(1+nkb,:))),[-1,+1]);
end;end;%for nka=0:K; for nkb=0:K;
if (verbose>1); disp(sprintf(' %% C_: ')); disp(num2str(C_)); end;

bx = besselj(0,jx); ba = transpose(j_Lx)\bx;
if (verbose>1); disp(sprintf(' %% ba: ')); disp(num2str(ba)); end;
cb = chebfun(@(x)besselj(0,x));
cr = chebfun(@(x)(1+x));
ca = zeros(K,1);
for nk=0:K-1;
cj = jacpoly(nk,0,1)*sqrt(nk+1)/sqrt(2);
ca(1+nk) = sum(cj*cr*cb);
end;%for nk=0:K-1;
if (verbose>1); disp(sprintf(' %% ca: ')); disp(num2str(ca)); end;
disp(sprintf(' %% error in a: %0.16f',norm(ba-ca,'fro')/norm(ca,'fro')));

end;%for K = 6:24;
