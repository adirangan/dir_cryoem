function L_out = maketsch(K);
% this function produces a matrix such that ;
% L_out(k,:) = coefficients of the kth tschebyscheff polynomial ;
% Warning: The coefficients are listed in "polynomial order" for use as matlab poly class ;
	
L_{1}=[1];
L_{2}=[1 0];
for index=1:K+1-2;
L_{index+2} = conv([2 0],L_{index+1}) - [0 0 L_{index}];
end;%for index=1:K+1-2;
L_out=zeros(K+1);
for index=1:K+1;
L_out(index,K+1-index+1:K+1)=L_{index};
end;%for index=1:K+1;
