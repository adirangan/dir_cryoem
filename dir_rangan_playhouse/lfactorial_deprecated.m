function output = lfactorial(input);

output = zeros(size(input));
for nl=1:numel(output);
if (input(nl)<=12);
output(nl) = log(factorial(input(nl)));
end;%if (input(nl)<=12);
if (input(nl)>12);
output(nl) = sfact(input(nl));
end;%if (input(nl)>12);
end;%for nl=1:numel(output);

%%%%%%%%%%%%%%%%;

function output = nfact(d);
dij=find(d>0);
output(dij) = log(factorial(d(dij)));
output(find(d<=0)) = 0;

function output = sfact(d);
dij=find(d>0);
output(dij) = d(dij).*log(d(dij)) - d(dij) + 0.5.*log(2*pi*d(dij)) + 1/12./d(dij) - 1/360./d(dij).^3 + 1/1260./d(dij).^5 - 1/1680./d(dij).^7 ;
output(find(d<=0)) = 0;
