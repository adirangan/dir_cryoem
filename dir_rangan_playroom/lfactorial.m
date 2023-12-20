function output = lfactorial(n);

if nargin<1;
disp(sprintf(' %% testing lfactorial'));
x_ = 0:1024; plot(x_,(lfactorial(x_)-gammaln(x_+1))./max(1e-12,abs(gammaln(x_+1))),'r.-');
xlabel('x'); xlim([0,max(x_)]);
ylabel('rel error');
disp('returning');return;
end;%if nargin<1;

nij = find(n< 15);
output(nij) = nfact(n(nij));
nij = find(n>=15);
output(nij) = sfact(n(nij));
output(find(~isfinite(output)))=-16;

function output = nfact(d);
output = d;
dij=find(d>0);
output(dij) = log(factorial(d(dij)));
output(find(d<=0)) = 0;

function output = sfact(d);
output = d;
dij=find(d>0);
output(dij) = d(dij).*log(d(dij)) - d(dij) + 0.5.*log(2*pi*d(dij)) + 1/12./d(dij) - 1/360./d(dij).^3 + 1/1260./d(dij).^5 - 1/1680./d(dij).^7 ;
output(find(d<=0)) = 0;
