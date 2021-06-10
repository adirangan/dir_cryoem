function [lx,lw,Lx,Lv] = orthopoly_node_weight_matrix_0(K,p_);
% This function generates nodes lx and weights lw associated with the set of polynomials ;
% of degree less than or equal to K-1 that are orthonormal with respect to weight function p(x) ;
% on the interval x\in[-1,+1], %
% where p_(x) is the polynomial with coefficients given by p_. ;
% The third output Lx is the array of polynomial evaluations at the (sorted) nodes. ;
% That is, Lx(1+j,j) = L_{j}(lx(j)), where L_{j}(x) is the orthonormal polynomial of degree j. ;
% The fourth output Lv is the array of polynomial coefficients. ;
% That is, Lv(1+j,:) are the polynomial coefficients for L_{j}(x) (in matlab order). ;

if nargin<1;
orthopoly_test_0();
disp('returning'); return;
end;% if nargin<1;

Lv = orthopoly_vec_0(K,p_);
lx = sort(roots(Lv(1+K,:)));
Lx = zeros(K,K); li = zeros(K,1);
for nr1=0:K-1;
Lx(1+nr1,:) = polyval(Lv(1+nr1,:),lx);
li(1+nr1) = polyint(conv(p_,Lv(1+nr1,:)),[-1,+1]);
end;%for nr1=0:K-1;
lw = Lx\li;

