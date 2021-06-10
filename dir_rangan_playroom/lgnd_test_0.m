function lgnd_test_0();

K = input('how many points would you like to use for interpolation (default = 16): ');if isempty(K);K=16;end;
a = input('The interval starts at the point (default = -2): ');if isempty(a);a=-2;end;
c = input('The interval stops at the point (default = 4): ');if isempty(c);c=4;end;
M = input('how many points would you like to use for your graphs (default = 1000): ');if isempty(M);M=1000;end;
xvals = a:(c-a)/M:c;
[Lx,lx] = lgnd_matrix_node_0(K); Lx=Lx(1:K,1:K);
Lv = lgnd_poly_vec_0(K); Lv=Lv(1:K,:);
[lx,lw] = lgnd_node_weight_0(K,a,c);
for index=1:length(lx);
fvals(index)=feval(@f,lx(index));
end;%for index=1:length(lx);
coefs = inv(transpose(Lx))*(fvals(:));
p = transpose(transpose(Lv)*coefs(:));
pvals=polyval(p,(xvals-(c+a)/2)/((c-a)/2));
for index=1:length(xvals);
yvals(index)=feval(@f,xvals(index));
end;%for index=1:length(xvals);
figure(1);clf;hold on;
plot(xvals,yvals,'r-');
plot(xvals,pvals,'b-');
plot(lx,fvals,'bx');
plot(lx,0,'bx');
xlabel('x');
ylabel('function value');
title('function in red, polynomial in blue, nodes at *');
hold off;
disp(sprintf('the calculated integral of your function is %f',transpose(lw(:)) * fvals(:)));
disp(sprintf('the reiman_sum integral of your function is %f',dot(yvals(2:end),diff(xvals))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
function output = f(x);

% This can be whatever you want it to be...

alpha = 1/10; sigma = 2;
output = 1/(x^2 + alpha); 
output = 1/(x^2 + alpha) + sin(2*pi*x) + cos(4*pi*x)*(1/sqrt(2*pi)/sigma)*exp(-(x)^2/(2*sigma.^2));
% output = atan(x);
