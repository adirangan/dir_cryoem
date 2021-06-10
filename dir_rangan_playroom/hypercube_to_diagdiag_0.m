function [n_p,p1_,p2_,DD] = hypercube_to_diagdiag_0(A,Bandlimit_M,Bandlimit_L);
% We assume that A is an n-by-n matrix, ;
% associated with spherical-harmonic coefficients m_ and l_, ;
% with maximum l-index n_l, such that m varies more quickly than l, ;
% (with the former ranging from -l to +l) : ;
% (m_,l_) = 
%  { 
%    (0,0) , 
%    (-1,1),(0,1),(+1,1) , 
%    (-2,2),(-1,2),(0,2),(+1,2),(+2,2) , 
%    ... , 
%    (-n_l,n_l),(-n_l+1,n_l),...,(+n_l-1,n_l),(+n_l,n_l) 
%  }. ;
% Note that, with this ordering, ;
% the index (m,l) has position p=(l*l + l + m), ;
% assuming that the first position index is p=0. ;
% The inverse of this relationship is:
% l = floor(sqrt(p)); m = mod(p,l);
% ;
% Within DD we store the diagonal elements of A in row dominant order, ;
% with p1 <-- (m1,l1) varying more quickly than p2 <-- (m2,l2). ;

if (nargin<3);

verbose=1;
n_l = 25; n_lm = (1+n_l).^2; n_m = 1+2*n_l;
[m1_,l1_,m2_,l2_,pp_,qq_] = permute_ml_to_lm(n_l);
m3_ = unique(m2_,'stable');
if verbose; disp(sprintf(' %% calculating Wd_')); tic; end;
omega = pi/6;
Wdf__ = wignerd_b(n_l,+omega);
Wdf_ = zeros(n_lm,n_lm);
Wdb__ = wignerd_b(n_l,-omega);
Wdb_ = zeros(n_lm,n_lm);
nlm=0;
for nl=0:n_l;
l_val = nl;
nm = 1 + 2*l_val;
Wdf_(nlm + (1:nm),nlm + (1:nm)) = Wdf__{1+nl};
Wdb_(nlm + (1:nm),nlm + (1:nm)) = Wdb__{1+nl};
nlm = nlm + nm;
end;%for nl=0:n_l;
assert(nlm==n_lm);
Sdf_ = sparse(Wdf_);
Sdb_ = sparse(Wdb_);
if verbose; disp(sprintf(' %% finished Wd_, total time %0.2f',toc)); end;
clear Wdf__ Wdf_ Wdb__ Wdb_;
if verbose; disp(sprintf(' %% calculating Wt_')); tic; end;
Wz__ = wignerz_leg(n_l,2.0,0.5);
Wz_ = zeros(n_lm,n_lm);
nlm=0;
for nm=1:(1+2*n_l);
nl = length(Wz__{nm});
Wz_(nlm+(1:nl),nlm+(1:nl)) = Wz__{nm};
nlm = nlm+nl;
end;%for nl=0:n_l;
assert(nlm==n_lm);
Sz_ = sparse(Wz_);
clear Wz__ Wz_;
if verbose; disp(sprintf(' %% finished Wt_, total time %0.2f',toc)); end;
Ss_ = Sdb_*Sz_(qq_,qq_)*Sdf_;
Bandlimit_M = 5; Bandlimit_L = 4;
[n_p,p1_,p2_,DD] = hypercube_to_diagdiag_0(Ss_,Bandlimit_M,Bandlimit_L);
Ds_ = sparse(1+p1_,1+p2_,DD,n_lm,n_lm);
figure; colormap('parula');
subplot(2,5,1); imagesc(log10(abs(Sdf_)),[-9,0]); title('Sdf (m,l)'); 
set(gca,'XTick',[],'YTick',[]);
subplot(2,5,6); imagesc(log10(abs(Sdf_(pp_,pp_))),[-9,0]); title('Sdf (l,m)');
set(gca,'XTick',[],'YTick',[]);
subplot(2,5,2); imagesc(log10(abs(Sz_(qq_,qq_))),[-9,0]); title('Sz (m,l)');
set(gca,'XTick',[],'YTick',[]);
subplot(2,5,7); imagesc(log10(abs(Sz_)),[-9,0]); title('Sz (l,m)');
set(gca,'XTick',[],'YTick',[]);
subplot(2,5,3); imagesc(log10(abs(Ss_)),[-9,0]); title('Ss (m,l)');
set(gca,'XTick',[],'YTick',[]);
subplot(2,5,8); imagesc(log10(abs(Ss_(pp_,pp_))),[-9,0]); title('Ss (l,m)');
set(gca,'XTick',[],'YTick',[]);
subplot(2,5,4); imagesc(log10(abs(Ds_)),[-9,0]); title('Ds (m,l)');
set(gca,'XTick',[],'YTick',[]);
subplot(2,5,9); imagesc(log10(abs(Ds_(pp_,pp_))),[-9,0]); title('Ds (l,m)');
set(gca,'XTick',[],'YTick',[]);
subplot(2,5,5); imagesc(log10(abs(Ss_ - Ds_)),[-9,0]); title('Ss-Ds (m,l)');
set(gca,'XTick',[],'YTick',[]);
subplot(2,5,10); imagesc(log10(abs(Ss_(pp_,pp_) - Ds_(pp_,pp_))),[-9,0]); title('Ss-Ds (l,m)');
set(gca,'XTick',[],'YTick',[]);
disp(sprintf(' %% Bandlimit_M %d Bandlimit_L %d',Bandlimit_M,Bandlimit_L));
disp(sprintf(' %% compression numel(DD) %d, nnz(Ss)/numel(Ss) %d/%d',numel(DD),nnz(Ss_),numel(Ss_)));
disp(sprintf(' %% Error norm(Ss-Ds)/norm(Ss) = %0.9f',norm(full(Ss_ - Ds_))/norm(full(Ss_))));


disp('returning'); return;
end;%if (nargin<3);

[nrows,ncols] = size(A); assert(nrows==ncols);
n_lm = nrows; n_l = round(sqrt(n_lm)) - 1; assert(n_lm==(1+n_l).^2);

n_p = 0;
for l1=0:n_l; for m1=-l1:+l1; 
p1 = l1*l1 + l1 + m1;
%for l2=0:n_l; for m2=-l2:+l2; 
for l2=max(0,l1-Bandlimit_L):min(n_l,l1+Bandlimit_L);
for m2=max(-l2,m1-Bandlimit_M):min(+l2,m1+Bandlimit_M);
p2 = l2*l2 + l2 + m2;
n_p = n_p+1;
end;end;%for l2=0:n_l; for m2=-l2:+l2;
end;end;%for l1=0:n_l; for m1=-l1:+l1;

DD = zeros(n_p,1);
p1_ = zeros(n_p,1);
p2_ = zeros(n_p,1);

n_p = 0;
for l1=0:n_l; for m1=-l1:+l1; 
p1 = l1*l1 + l1 + m1;
%for l2=0:n_l; for m2=-l2:+l2; 
for l2=max(0,l1-Bandlimit_L):min(n_l,l1+Bandlimit_L);
for m2=max(-l2,m1-Bandlimit_M):min(+l2,m1+Bandlimit_M);
p2 = l2*l2 + l2 + m2;
DD(1+n_p) = A(1+p1,1+p2);
p1_(1+n_p) = p1;
p2_(1+n_p) = p2;
n_p = n_p+1;
end;end;%for l2=0:n_l; for m2=-l2:+l2;
end;end;%for l1=0:n_l; for m1=-l1:+l1;
