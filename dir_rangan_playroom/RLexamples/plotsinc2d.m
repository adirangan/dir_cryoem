x = [-10:0.01:10]
y = zeros(2001,10);
for i=1:10
aa = 5*i
y(:,i) = sincgauss(x,aa)
end
mesh(y)

