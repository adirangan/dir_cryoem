function [Lv] = lgnd_poly_vec_0(K);
% This function generates a matrix Lv ;
% such that Lv(k,:) contains the polynomial coefficients of the kth legendre polynomial. ;
% Warning: the coefficients are listed in 'polynomial order' for use in matlab poly class. ;

if nargin<1;
lgnd_test_0.m;
disp('returning'); return;
end;% if nargin<1;

L_{1}=[1];
L_{2}=[1 0];
for index=1:K-1;
L_{index+2} = ((2*index+1)*conv([1 0],L_{index+1}) - index*[0 0 L_{index}])/(index+1);
end;%for index=1:K-1;

for index=0:K;
A = polyint(conv(L_{1+index},L_{1+index}),[-1,+1]);
L_{1+index} = L_{1+index}/sqrt(A);
end;%for index=0:K;

Lv=zeros(K+1);
for index=1:K+1;
Lv(index,:)=[zeros(1,K+1-index) L_{index}];
end;%for index=1:K+1;
