function W_ = wignerd(n_l,beta) ;
% generates wigner-d matrices up to n_l;
% uses 3-term recurrences. ;
cb = cos(beta/2); sb = sin(beta/2); 
W_ = cell(1+n_l,1);
W_{1} = [1];
for nl=1:n_l;
nlp = nl-1;
V = W_{nl};
W = zeros(2*nl+1,2*nl+1);
for nmn = -nl:+nl; for nmp = -nl:+nl;
%%%%%%%%%%%%%%%%;
if (nl~=-nmp & nl~=1-nmp); % use recurrence A ;
W1 = 0; if (abs(nmp-1)<=nlp & abs(nmn-1)<=nlp); W1 = V(1+nlp+(nmp-1),1+nlp+(nmn-1)); end;
W2 = 0; if (abs(nmp-1)<=nlp & abs(nmn-0)<=nlp); W2 = V(1+nlp+(nmp-1),1+nlp+(nmn-0)); end;
W3 = 0; if (abs(nmp-1)<=nlp & abs(nmn+1)<=nlp); W3 = V(1+nlp+(nmp-1),1+nlp+(nmn+1)); end;
tmp = cb*cb*A1(nl,nmn,nmp)*W1 - 2*cb*sb*A2(nl,nmn,nmp)*W2 + sb*sb*A3(nl,nmn,nmp)*W3;
end;%if (nl~=-nmp & nl~=1-nmp); % use recurrence A ;
%%%%%%%%%%%%%%%%;
if (nl~=+nmp & nl~=1+nmp); % use recurrence B ;
W1 = 0; if (abs(nmp+1)<=nlp & abs(nmn-1)<=nlp); W1 = V(1+nlp+(nmp+1),1+nlp+(nmn-1)); end;
W2 = 0; if (abs(nmp+1)<=nlp & abs(nmn-0)<=nlp); W2 = V(1+nlp+(nmp+1),1+nlp+(nmn-0)); end;
W3 = 0; if (abs(nmp+1)<=nlp & abs(nmn+1)<=nlp); W3 = V(1+nlp+(nmp+1),1+nlp+(nmn+1)); end;
tmp = sb*sb*B1(nl,nmn,nmp)*W1 + 2*sb*cb*B2(nl,nmn,nmp)*W2 + cb*cb*B3(nl,nmn,nmp)*W3;
end;%if (nl~=+nmp & nl~=1+nmp); % use recurrence B ;
%%%%%%%%%%%%%%%%;
if (nl~=-nmp & nl~=+nmp); % use recurrence C ;
W1 = 0; if (abs(nmp-0)<=nlp & abs(nmn-1)<=nlp); W1 = V(1+nlp+(nmp-0),1+nlp+(nmn-1)); end;
W2 = 0; if (abs(nmp-0)<=nlp & abs(nmn-0)<=nlp); W2 = V(1+nlp+(nmp-0),1+nlp+(nmn-0)); end;
W3 = 0; if (abs(nmp-0)<=nlp & abs(nmn+1)<=nlp); W3 = V(1+nlp+(nmp-0),1+nlp+(nmn+1)); end;
tmp = sb*cb*C1(nl,nmn,nmp)*W1 + (cb*cb-sb*sb)*C2(nl,nmn,nmp)*W2 - sb*cb*C3(nl,nmn,nmp)*W3;
end;%if (nl~=-nmp & nl~=+nmp); % use recurrence C ;
%%%%%%%%%%%%%%%%;
W(1+nl+nmp,1+nl+nmn) = tmp;
end;end;%for nmn = -nl:+nl; for nmp = -nl:+nl;
W_{1+nl} = W;
end;%for nl=1:n_l;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output = A1(nl,nmn,nmp);
numerator   = (nl+nmn)*(nl+nmn-1);
denominator = (nl+nmp)*(nl+nmp-1);
output = sqrt(numerator/denominator);

function output = A2(nl,nmn,nmp);
numerator   = (nl+nmn)*(nl-nmn);
denominator = (nl+nmp)*(nl+nmp-1);
output = sqrt(numerator/denominator);

function output = A3(nl,nmn,nmp);
numerator   = (nl-nmn)*(nl-nmn-1);
denominator = (nl+nmp)*(nl+nmp-1);
output = sqrt(numerator/denominator);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output = B1(nl,nmn,nmp);
numerator   = (nl+nmn)*(nl+nmn-1);
denominator = (nl-nmp)*(nl-nmp-1);
output = sqrt(numerator/denominator);

function output = B2(nl,nmn,nmp);
numerator   = (nl+nmn)*(nl-nmn);
denominator = (nl-nmp)*(nl-nmp-1);
output = sqrt(numerator/denominator);

function output = B3(nl,nmn,nmp);
numerator   = (nl-nmn)*(nl-nmn-1);
denominator = (nl-nmp)*(nl-nmp-1);
output = sqrt(numerator/denominator);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output = C1(nl,nmn,nmp);
numerator   = (nl+nmn)*(nl+nmn-1);
denominator = (nl-nmp)*(nl+nmp);
output = sqrt(numerator/denominator);

function output = C2(nl,nmn,nmp);
numerator   = (nl+nmn)*(nl-nmn);
denominator = (nl-nmp)*(nl+nmp);
output = sqrt(numerator/denominator);

function output = C3(nl,nmn,nmp);
%numerator   = (nl-nmn)*(nl-nmn+1);
numerator   = (nl-nmn)*(nl-nmn-1);
denominator = (nl-nmp)*(nl+nmp);
output = sqrt(numerator/denominator);
