function [Lv] = orthopoly_vec_0(K,p_);
% This function generates a matrix Lv ;
% such that Lv(k,:) contains the polynomial coefficients of the kth orthonormal polynomial L_{k}(x), ;
% which is designed so that: ;
% \int_{-1}^{+1} L_{j}(x) * L_{k}(x) * p(x) dx = delta_{jk}, ;
% where p(x) is the polynomial with coefficients given by p_. ;
% Note: the coefficients of Lv(k,:) are listed in 'polynomial order' for use in matlab poly class. ;

if nargin<1;
orthopoly_test_0();
disp('returning'); return;
end;% if nargin<1;

clear L_;
L_ = cell(1+K,1);
nd1 = 0; % polynomial degree ;
L_{1+nd1} = [zeros(1,1+K - (1+nd1)) , 1 , zeros(1,1+nd1-1)]; 
A = polyint(conv(p_,conv(L_{1+nd1},L_{1+nd1})),[-1,+1]); 
if A<1e-6;
disp(sprintf(' %% Warning! A %0.8f < 1e-6 at order 0 in orthopoly_vec',A));
end;%if A<1e-6;
L_{1+nd1} = L_{1+nd1}/sqrt(A);
for nd2=1:K;
L_{1+nd2} = [zeros(1,1+K - (1+nd2)) , 1 , zeros(1,1+nd2-1)]; 
for nd1=0:nd2-1;
A = polyint(conv(p_,conv(L_{1+nd1},L_{1+nd2})),[-1,+1]); 
L_{1+nd2} = L_{1+nd2} - A*L_{1+nd1};
A = polyint(conv(p_,conv(L_{1+nd2},L_{1+nd2})),[-1,+1]);  
if A<1e-6;
disp(sprintf(' %% Warning! A %0.8f < 1e-6 at order nd2=%d in orthopoly_vec',A,nd2));
end;%if A<1e-6;
L_{1+nd2} = L_{1+nd2}/sqrt(A);
end;%for nd1=0:nd2-1;
A = polyint(conv(p_,conv(L_{1+nd2},L_{1+nd2})),[-1,+1]); 
if A<1e-6;
disp(sprintf(' %% Warning! A %0.8f < 1e-6 at order nd2=%d in orthopoly_vec',A,nd2));
end;%if A<1e-6;
L_{1+nd2} = L_{1+nd2}/sqrt(A);
end;%for nd2=1:K;

Lv = zeros(1+K,1+K);
for nd1=0:K;
Lv(1+nd1,:) = L_{1+nd1};
end;%for nd1=0:K;

