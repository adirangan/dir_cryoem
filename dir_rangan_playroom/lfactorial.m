function output = lfactorial(n);
if (n<15);
output = nfact(n);
output(find(~isfinite(output)))=-16;
 else;
output = sfact(n);
output(find(~isfinite(output)))=-16;
end;%if (n<15);

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
