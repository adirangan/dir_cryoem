function BB = hypercube_to_bandlimits_0(A,Bandlimit_M,Bandlimit_L);
% This function calculates the error associated with approximating A with a bandlimited version of A. ;
% ;
% Warning! This function is very inefficient (and can be improved drastically). ;
% ;
% The inputs Bandlimit_M and Bandlimit_L are the maximum bandlimits to investigate. ;
% The output array BB(1+bm,1+bl) contains the relative error (in the 2-norm) ;
% associated with approximating A with Bandlimit_M==bm and Bandlimit_L==bl. ;

n_lm = length(A);
E0 = norm(full(A));
BB = zeros(1+Bandlimit_M,1+Bandlimit_L);
for bm=0:Bandlimit_M; for bl=0:Bandlimit_L;
[n_p,p1_,p2_,DD] = hypercube_to_diagdiag_0(A,bm,bl);
B = sparse(1+p1_,1+p2_,DD,n_lm,n_lm);
BB(1+bm,1+bl) = norm(full(A-B))/E0;
end;end;%for bm=0:Bandlimit_M; for bl=0:Bandlimit_L;
