function [tsch_vector,tsch_nodes] = tschintmat(K,a,c) ;
% This function establishes the tschebyscheff nodes of integration as well as the corresponding integration vector ;
% for degree K tschebyscheff interpolation/integration on the interval [a,c] ;
tsch_nodes=cos((2*K-2*[1:K]+1)*pi/(2*K));
T(1,:)=linspace(1,1,K); T(2,:)=tsch_nodes;
for j1=2:K-1 ; for j2=1:K ;
T(j1+1,j2)=2*tsch_nodes(j2)*T(j1,j2)-T(j1-1,j2);
end; end; %for j1=2:K-1 ; for j2=1:K ;
T_int = ~isint((1:K)./2)./((1:K).*((1:K)-2));T_int(2)=0;
T_eye=2*eye(K); T_eye(1,1)=1;
tsch_vector = ((a-c)/K)*T_int*T_eye*T;
tsch_nodes = (c-a)/2 * tsch_nodes + (c+a)/2;
    
function output = isint(input) ;
% this is just a function that tests whether or not the input is an integer ;
output = floor(input)==input;
