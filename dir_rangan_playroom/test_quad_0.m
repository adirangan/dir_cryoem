function output = test_quad_0;

K = input('how many points would you like to use for interpolation (default = 10): ');if isempty(K);K=10;end;
TorL = input('Would you like to use Tschebyscheff (T) or Legendre (L) nodes (default = T): ','s');if isempty(TorL);TorL='T';end;
a = input('The interval starts at the point (default = -4): ');if isempty(a);a=-4;end;
c = input('The interval stops at the point (default = 4): ');if isempty(c);c=4;end;
M = input('how many points would you like to use for your graphs (default = 1000): ');if isempty(M);M=1000;end;
xvals = a:(c-a)/M:c;

if TorL=='T'
    [T,t] = matetsch(K);T=T(1:K,1:K);
    Tv = maketsch(K);Tv=Tv(1:K,2:K+1);
    [tv,tn] = tschintmat(K,a,c);
    for index=1:length(tn)
        fvals(index)=feval(@f,tn(index));
    end;
    T_eye=2*eye(K); T_eye(1,1)=1; T_eye=T_eye/K;
    coefs = T_eye*T*(fvals(:));
    p = (Tv'*coefs(:))';
    pvals=polyval(p,(xvals-(c+a)/2)/((c-a)/2));
    for index=1:length(xvals);
        yvals(index)=feval(@f,xvals(index));
    end;
    figure(1);clf;hold on; subplot(1,1,1); hold on;
    plot(xvals,yvals,'r-');
    plot(xvals,pvals,'b-');
    plot(tn,fvals,'bx');
    plot(tn,0,'bx');
    xlabel('x');
    ylabel('function value');
    title('function in red, polynomial in blue, nodes at *');
    hold off;
    disp(sprintf('the calculated integral of your function is %f',transpose(tv(:)) * fvals(:)));
    disp(sprintf('the reiman_sum integral of your function is %f',dot(yvals(2:end),diff(xvals))));
elseif TorL=='L'
    [L,l] = mateleg(K);L=L(1:K,1:K);
    Lv = makeleg(K);Lv=Lv(1:K,2:K+1);
    [lv,ln] = legintmat(K,a,c);
    for index=1:length(ln)
        fvals(index)=feval(@f,ln(index));
    end;
    coefs = inv(transpose(L))*(fvals(:));
    p = (Lv'*coefs(:))';
    pvals=polyval(p,(xvals-(c+a)/2)/((c-a)/2));
    for index=1:length(xvals);
        yvals(index)=feval(@f,xvals(index));
    end;
    figure(1);clf;hold on;
    plot(xvals,yvals,'r-');
    plot(xvals,pvals,'b-');
    plot(ln,fvals,'bx');
    plot(ln,0,'bx');
    xlabel('x');
    ylabel('function value');
    title('function in red, polynomial in blue, nodes at *');
    hold off;
    disp(sprintf('the calculated integral of your function is %f',transpose(lv(:)) * fvals(:)));
    disp(sprintf('the reiman_sum integral of your function is %f',dot(yvals(2:end),diff(xvals))));
else
    error('T or L...');
end;
    
    

function Legpoly_vec = makeleg(K);
	% This function generates a matrix Legpoly_vec
	% such that Legpoly_vec(k,:) is the coefficients of the kth legendre polynomial
	% Warning: the coefficients are listed in "polynomial order" for use in matlab poly class
	
	L{1}=[1];
	L{2}=[1 0];
	for index=1:K-1
        L{index+2} = ((2*index+1)*conv([1 0],L{index+1}) - index*[0 0 L{index}])/(index+1);
	end;
	Legpoly_vec=zeros(K+1);
	for index=1:K+1
        Legpoly_vec(index,:)=[zeros(1,K+1-index) L{index}];
	end;

function [L,l] = mateleg(K)
	
	% L = matrix whose components L(j,k) are jth Legendre polynomial evaluated at kth legendre node
	% l = vector of legendre nodes
	
	Legpoly_vec = makeleg(K);
	l=roots(Legpoly_vec(K+1,:));l=sort(transpose(l(:)));
	for index1=1:K
        for index2=1:K
            L(index1,index2) = polyval(Legpoly_vec(index1,:),l(index2));
        end;
	end;

function [leg_vector,leg_nodes] = legintmat(K,a,c)

	% This function calculates the weights (leg_vector) and nodes (leg_nodes)
	% for legendre integration of degree K on interval [a,c]
	
	Legpoly_vec = makeleg(K);
	Legpoly_int = zeros(1,K);
	for index=1:K
        tempoly = polyint(Legpoly_vec(index,:));
        Legpoly_int(index) = polyval(tempoly,1)-polyval(tempoly,-1);
	end;
	[L,l]=mateleg(K);
	L_inv = L\eye(K);
	leg_vector = ((c-a)/2) * Legpoly_int*transpose(L_inv);
	leg_nodes = l*(c-a)/2 + (c+a)/2;

function [T,t] = matetsch(K)

	% K = number of tschebyscheff polynomials to construct
	% T = matrix of values T(i,j) = Ti(tj)
	% where Ti = ith tschebyscheff polynomial,
	% tj = jth root of TK
	
	t=cos((2*K-2*[1:K]+1)*pi/(2*K));
	T(1,:)=linspace(1,1,length(t));
	T(2,:)=t;
	for j1=2:K
      for j2=1:length(t)
        T(j1+1,j2)=2*t(j2)*T(j1,j2)-T(j1-1,j2);
      end;
	end;

function Tschpoly_vec = maketsch(K);
	% this function produces a matrix such that
	% Tschpoly_vec(k,:) = coefficients of the kth tschebyscheff polynomial
	% Warning: The coefficients are listed in "polynomial order" for use as matlab poly class
	
	Tschpoly{1}=[1];
	Tschpoly{2}=[1 0];
	for index=1:K+1-2
        Tschpoly{index+2} = conv([2 0],Tschpoly{index+1}) - [0 0 Tschpoly{index}];
	end;
	Tschpoly_vec=zeros(K+1);
	for index=1:K+1
        Tschpoly_vec(index,K+1-index+1:K+1)=Tschpoly{index};
	end;

function [tsch_vector,tsch_nodes] = tschintmat(K,a,c)

	% This function establishes the tschebyscheff nodes of integration as well as the corresponding integration vector 
	% for degree K tschebyscheff interpolation/integration on the interval [a,c]
	
	tsch_nodes=cos((2*K-2*[1:K]+1)*pi/(2*K));
	T(1,:)=linspace(1,1,K);
	T(2,:)=tsch_nodes;
	for j1=2:K-1
      for j2=1:K
        T(j1+1,j2)=2*tsch_nodes(j2)*T(j1,j2)-T(j1-1,j2);
      end;
	end;
	T_int = ~isint((1:K)./2)./((1:K).*((1:K)-2));T_int(2)=0;
	T_eye=2*eye(K); T_eye(1,1)=1;
	tsch_vector = ((a-c)/K)*T_int*T_eye*T;
	tsch_nodes = (c-a)/2 * tsch_nodes + (c+a)/2;
    
    
function output = polyint(input)
	
	% this function finds the indefinite integral of a polynomial
	
        input=transpose(input(:));
	for index=1:length(input);
        output(index) = input(index)/(length(input)-index+1);
	end;
	output = [output 0];
    
function output = isint(input)
    
    % this is just a function that tests whether or not the input is an integer
    
	output = floor(input)==input;

function output = f(x);

    % This can be whatever you want it to be...
    
    alpha = 1/10; sigma = 2;
    output = 1/(x^2 + alpha); 
    output = 1/(x^2 + alpha) + sin(2*pi*x) + cos(4*pi*x)*(1/sqrt(2*pi)/sigma)*exp(-(x)^2/(2*sigma.^2));
    % output = atan(x);
