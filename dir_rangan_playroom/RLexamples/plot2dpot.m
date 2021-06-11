x = [-10:0.01:10]
yy = zeros(4001,10);
for i=1:20
y = -exp(-2*x.*x) - 2*exp(-3*(x+3).*(x+3)) + 1.5*sin(5*x).*exp(-0.05*x.*x)
aa = 0.05*(2*i)
z = exp(-aa*x.*x);
w = conv(z,y);
yy(:,i) = w;
end
figure(1)
mesh(-log(-yy+20))
figure(2)
contourf(yy)
%
%subplot(1,2,1); mesh(-log(-yy+20))
%subplot(1,2,2); contourf(yy)

