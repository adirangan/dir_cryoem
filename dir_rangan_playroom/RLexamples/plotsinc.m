x = [-10:0.01:10]
plot(x,-sinc(x))

hold on
x0 = -3
plot(x0,-sinc(x0),'x')
for i = 1:6
pause
x1 = x0 - sincp(x0)/sincpp(x0); 
sincp(x1);
plot(x1,-sinc(x1),'x')
x0 = x1
end
