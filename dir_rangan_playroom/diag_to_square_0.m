function A = diag_to_square_0(D,l);
% This function converts column-vector D into ;
% a square l-by-l matrix A. ;
% We assume the vector D stores the diagonal entries of A ;
% up to some band b. ;
% We assume the diagonal bands are stored in the following order: ;
% [0 , -1 , +1 , -2 , +2 , ... , -b , +b ]. ;
% Thus, the total length of D should be: ;
% l + (l-1) + (l-1) + (l-2) + (l-2) + ... + (l-b) + (l-b). ;
% Consequently, the total length of D should be : ;
% ll = (2b+1)*l - b*(b+1) ;

if (nargin<2);
l = 5; b = 2;
A=reshape(1:l^2,l,l);
D=square_to_diag_0(A,b); 
B=diag_to_square_0(D,l);
subplot(1,2,1); imagesc(A,[1,l^2]); title('original');
subplot(1,2,2); imagesc(B,[1,l^2]); title(sprintf('diag: b=%d',b));
disp('returning'); return;
end;%if (nargin<2);

A = zeros(l);
ll = size(D,1);
b = round(((2*l-1) - sqrt((2*l-1).^2 - 4*(ll-l)))/2);
if (ll~=(2*b+1)*l - b*(b+1)); disp(sprintf(' Warning! ll %d not compatible with l,b %d,%d in diag_to_square',ll,l,b)); end;

b_val_ = [0:-1:-b ; 0:+1:+b];
b_val_ = reshape(b_val_,2*(b+1),1);
b_val_ = b_val_(2:end);

l_sub = 0;
for nb_val=1:length(b_val_);
b_val = b_val_(nb_val);
b_abs = abs(b_val);
if (b_val>=0);
for nd=1:l-b_abs;
A(nd,nd+b_abs) = D(l_sub + nd);
end;%for nd=1:l-b_abs;
end;%if (b_val>=0);
if (b_val< 0);
for nd=1:l-b_abs;
A(nd+b_abs,nd) = D(l_sub + nd);
end;%for nd=1:l-b_abs;
end;%if (b_val< 0);
l_sub = l_sub + l-abs(b_val);
end;%for nb_val=1:length(b_val_);

