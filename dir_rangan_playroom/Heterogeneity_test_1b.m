% testing residual least-squares for 3-molecule heterogeneity. ;
% just using random matrices (not operators built by slicing the sphere). ;

DM = 20; DF = 40;
N = 400; pA = 4/7; pB = 3/7; pC = 0/7;
NA = floor(pA*N); NB = floor(pB*N); NC = N-NA-NB;
EA = randn(NA*DM,DF); EB = randn(NB*DM,DF); EC = randn(NC*DM,DF); E = [EA;EB;EC];
FX = randn(DF,1); FA = FX + 0.1*randn(DF,1); FB = FX + 0.1*randn(DF,1); FC = FX + 0.1*randn(DF,1);
F = pA*FA+pB*FB+pC*FC;
MA2 = reshape(EA*FA,DM,NA); MA = MA2(:);
MB2 = reshape(EB*FB,DM,NB); MB = MB2(:);
MC2 = reshape(EC*FC,DM,NC); MC = MC2(:);
M2 = [MA2,MB2,MC2]; M = [MA;MB;MC];
F_r = E\M; R = M - E*F_r; G_r = E\R;
R2 = reshape(R,DM,N);
RA = R2(:,1:NA); RB = R2(:,NA+1:NA+NB); RC = R2(:,NA+NB+1:NA+NB+NC);

v1_iter = zeros(N,1);
v2_iter = zeros(N,1);

max_iteration = 6;
figure; prows = 2; pcols = ceil(max_iteration/prows);
for ni=1:max_iteration;
if (ni==1);
G1_base = transpose(E)*R; G1_base = G1_base/norm(G1_base);
G2_base = randn(DF,1); G2_base = G2_base/norm(G2_base);
 else;
V1 = sparse(1:DM*N,1:DM*N,reshape(repmat(transpose(v1_iter),DM,1),1,DM*N),DM*N,DM*N);
V2 = sparse(1:DM*N,1:DM*N,reshape(repmat(transpose(v2_iter),DM,1),1,DM*N),DM*N,DM*N);
G1_base = transpose(E)*V1*R; G1_base = G1_base/norm(G1_base);
G2_base = transpose(E)*V2*R; G2_base = G2_base - dot(G1_base,G2_base)*G1_base; G2_base = G2_base/norm(G2_base);
end% if first iteration ;
for nn=0:N-1;
Rn = R2(:,1+nn);
En = E(nn*DM + (1:DM),:);
v1_tmp = transpose(Rn)*En*G1_base/norm(Rn).^2;
v2_tmp = transpose(Rn)*En*G2_base/norm(Rn).^2;
if norm(v1_tmp)<1e-6; v1_iter(1+nn) = 1; end;
if norm(v1_tmp)>1e-6; v1_iter(1+nn) = v1_tmp; end;
if norm(v2_tmp)<1e-6; v2_iter(1+nn) = 1; end;
if norm(v2_tmp)>1e-6; v2_iter(1+nn) = v2_tmp; end;
end;% for nn=1:N;
subplot(prows,pcols,ni); hold on;
plot(v1_iter(1:NA),v2_iter(1:NA),'ro');
plot(v1_iter(NA+1:NA+NB),v2_iter(NA+1:NA+NB),'go');
plot(v1_iter(NA+NB+1:NA+NB+NC),v2_iter(NA+NB+1:NA+NB+NC),'bo');
hold off;
title(sprintf('iteration %d',ni));
end;%for ni=1:max_iteration;



