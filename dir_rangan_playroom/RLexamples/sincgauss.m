function y = sincgauss(x,a)
aa = sqrt(a);
eye = sqrt(-1);
y = -2*aa*(erf_(aa/2+pi*eye*x/aa)+erf_(aa/2 - pi*eye*x/aa)).*exp(-pi*pi*x.*x/a)
end
