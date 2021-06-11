x = [-10:0.01:10]
yy = zeros(4001,10);
for i=1:20
y = -sinc(x);
aa = 0.05*(2*i)
z = exp(-aa*x.*x);
w = conv(z,y);
yy(:,i) = w;
end
mesh(-log(-yy+20))

