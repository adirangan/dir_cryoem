function [Lx,lx] = lgnd_matrix_node_0(K);
% Lx = matrix whose components L(j,k) are jth Legendre polynomial evaluated at kth legendre node. ;
% lx = vector of legendre nodes. ;
 
if nargin<1;
lgnd_test_0.m;
disp('returning'); return;
end;% if nargin<1;

Lv = lgnd_poly_vec_0(K);
lx = lgnd_node_weight_0(K);

for index1=1:K; for index2=1:K;
Lx(index1,index2) = polyval(Lv(index1,:),lx(index2));;
end;end;%for index1=1:K; for index2=1:K;

