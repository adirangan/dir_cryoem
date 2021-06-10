% testing residual least-squares for simple 2-molecule heterogeneity. ;
% just using random matrices (not operators built by slicing the sphere). ;

DM = 20; DF = 40;
N = 400; pA = 4/7; if (pA<0.5) pA = 1-pA; end;
pB = 1-pA; NA = floor(pA*N); NB = N-NA;
E = randn(N*DM,DF);
EA = E(1:NA*DM,:);
EB = E(NA*DM+1:(NA+NB)*DM,:);
FA = randn(DF,1);
FD = 0.1*randn(DF,1);
FB = FA - FD;
F = pA*FA+pB*FB;
MA2 = reshape(EA*FA,DM,NA); MA = MA2(:);
MB2 = reshape(EB*FB,DM,NB); MB = MB2(:);
M2 = [MA2,MB2]; M = [MA;MB];
F_r = E\M;
R = M - E*F_r;
G_r = E\R;
R2 = reshape(R,DM,N);
RA = R2(:,1:NA);
RB = R2(:,NA+1:NA+NB);
EAFD = reshape(EA*FD,size(RA));
EBFD = reshape(EB*FD,size(RB));

p_iter = zeros(N,1);

max_iteration = 8;
for ni=1:max_iteration;
if (ni==1);
G_base = transpose(E)*R; G_base = G_base/norm(G_base);
 else;
P = sparse(1:DM*N,1:DM*N,reshape(repmat(transpose(p_iter),DM,1),1,DM*N),DM*N,DM*N);
G_base = transpose(E)*P*R; G_base = G_base/norm(G_base);
end% if first iteration ;
for nn=0:N-1;
Rn = R2(:,1+nn);
En = E(nn*DM + (1:DM),:);
p_tmp = transpose(Rn)*En*G_base/norm(Rn).^2;
if norm(p_tmp)<1e-6; p_iter(1+nn) = 1; end;
if norm(p_tmp)>1e-6; p_iter(1+nn) = p_tmp; end;
end;% for nn=1:N;
%figure;plot(p_iter,'.');
end;%for ni=1:max_iteration;
nn_pos = length(find(p_iter>0));
pp_pos = median(p_iter(find(p_iter>0)));
nn_neg = length(find(p_iter<0));
pp_neg = median(p_iter(find(p_iter<0)));
pp = pp_pos - pp_neg;
pp_pos = pp_pos/pp;
pp_neg = pp_neg/pp;
if nn_pos>nn_neg; pp_A = pp_pos; pp_B = pp_neg; end;
if nn_pos<nn_neg; pp_A = pp_neg; pp_B = pp_pos; end;
disp(sprintf(' %% pp_A %.2f pp_B %.2f; pA %.2f pB %.2f',pp_A,pp_B,pA,pB));



