function D = orthopoly_test_2(Lv,p_);
[nrows,ncols] = size(Lv);
D = zeros(nrows,nrows);
for nr1 = 1:nrows;
for nr2 = 1:nrows;
D(nr1,nr2) = polyint(conv(p_,conv(Lv(nr1,:),Lv(nr2,:))),[-1,+1]);
end;%for nr2 = 1:nrows;
end;%for nr1=1:nrows;
