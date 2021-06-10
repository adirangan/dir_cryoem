function [Et_,El_] = test_quad_1(a,c);

%K_ = 1:64;
K_ = round(2.^[0:0.125:8]);
Et_ = zeros(length(K_),1);
El_ = zeros(length(K_),1);

for nK=1:length(K_);
K = K_(nK);
[wt,xt] = tschintmat(K,a,c);
ft = zeros(K,1); for index=1:length(xt); ft(index)=feval(@f,xt(index)); end;
It = wt*ft;
%[wl,xl] = legintmat(K,a,c);
[xl,wl] = lgwt(K,a,c); wl = transpose(wl);
fl = zeros(K,1); for index=1:length(xl); fl(index)=feval(@f,xl(index)); end;
Il = wl*fl;
Ix = feval(@F,c) - feval(@F,a);
Et_(nK) = abs(Ix-It);
El_(nK) = abs(Ix-Il);
end;%for nK=1:length(K_);

plot(log2(K_),log(Et_),'b.-',log2(K_),log(El_),'r.-'); 
xlabel('log2(K)'); ylabel('log(E)'); legend('T','L');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

function Legpoly_vec = makeleg(K);
% This function generates a matrix Legpoly_vec ;
% such that Legpoly_vec(k,:) is the coefficients of the kth legendre polynomial ;
% Warning: the coefficients are listed in "polynomial order" for use in matlab poly class ;
L{1}=[1];
L{2}=[1 0];
for index=1:K-1;
L{index+2} = ((2*index+1)*conv([1 0],L{index+1}) - index*[0 0 L{index}])/(index+1);
end;
Legpoly_vec=zeros(K+1);
for index=1:K+1;
Legpoly_vec(index,:)=[zeros(1,K+1-index) L{index}];
end;

function [L,l] = mateleg(K);
% L = matrix whose components L(j,k) are jth Legendre polynomial evaluated at kth legendre node ;
% l = vector of legendre nodes ;
Legpoly_vec = makeleg(K);
l=roots(Legpoly_vec(K+1,:));l=sort(transpose(l(:)));
for index1=1:K ;
for index2=1:K ;
L(index1,index2) = polyval(Legpoly_vec(index1,:),l(index2));
end;
end;

function [leg_vector,leg_nodes] = legintmat(K,a,c) ;
% This function calculates the weights (leg_vector) and nodes (leg_nodes) ;
% for legendre integration of degree K on interval [a,c] ;
Legpoly_vec = makeleg(K);
Legpoly_int = zeros(1,K);
for index=1:K ;
tempoly = polyint(Legpoly_vec(index,:));
Legpoly_int(index) = polyval(tempoly,1)-polyval(tempoly,-1);
end;
[L,l]=mateleg(K);
L_inv = L\eye(K);
leg_vector = ((c-a)/2) * Legpoly_int*transpose(L_inv);
leg_nodes = l*(c-a)/2 + (c+a)/2;

function [T,t] = matetsch(K) ;
% K = number of tschebyscheff polynomials to construct ;
% T = matrix of values T(i,j) = Ti(tj) ;
% where Ti = ith tschebyscheff polynomial, ;
% tj = jth root of TK ;
t=cos((2*K-2*[1:K]+1)*pi/(2*K));
T(1,:)=linspace(1,1,length(t));
T(2,:)=t;
for j1=2:K ;
for j2=1:length(t) ;
T(j1+1,j2)=2*t(j2)*T(j1,j2)-T(j1-1,j2);
end;
end;

function Tschpoly_vec = maketsch(K); 
% this function produces a matrix such that ;
% Tschpoly_vec(k,:) = coefficients of the kth tschebyscheff polynomial ;
% Warning: The coefficients are listed in "polynomial order" for use as matlab poly class ;
Tschpoly{1}=[1];
Tschpoly{2}=[1 0];
for index=1:K+1-2 ;
Tschpoly{index+2} = conv([2 0],Tschpoly{index+1}) - [0 0 Tschpoly{index}];
end;
Tschpoly_vec=zeros(K+1);
for index=1:K+1 ;
Tschpoly_vec(index,K+1-index+1:K+1)=Tschpoly{index};
end;

function [tsch_vector,tsch_nodes] = tschintmat(K,a,c) ;
% This function establishes the tschebyscheff nodes of integration as well as the corresponding integration vector ;
% for degree K tschebyscheff interpolation/integration on the interval [a,c] ;
tsch_nodes=cos((2*K-2*[1:K]+1)*pi/(2*K));
T(1,:)=linspace(1,1,K);
T(2,:)=tsch_nodes;
for j1=2:K-1 ;
for j2=1:K ;
T(j1+1,j2)=2*tsch_nodes(j2)*T(j1,j2)-T(j1-1,j2);
end;
end;
T_int = ~isint((1:K)./2)./((1:K).*((1:K)-2));T_int(2)=0;
T_eye=2*eye(K); T_eye(1,1)=1;
tsch_vector = ((a-c)/K)*T_int*T_eye*T;
tsch_nodes = (c-a)/2 * tsch_nodes + (c+a)/2;
    
function output = polyint(input) ;
% this function finds the indefinite integral of a polynomial ;
input=transpose(input(:));
for index=1:length(input); 
output(index) = input(index)/(length(input)-index+1);
end;
output = [output 0];
    
function output = isint(input) ;
% this is just a function that tests whether or not the input is an integer ;
output = floor(input)==input;

function output = f(x);
w = 12;
output = exp(-(w*sin(x)).^2)*w*cos(x);

function output = F(x) ;
% integral of f(x);
w = 12;
output = sqrt(pi)/2*erf(w*sin(x));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

