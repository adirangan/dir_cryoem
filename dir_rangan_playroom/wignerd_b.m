function W_ = wignerd_b(n_l,beta) ;
% generates wigner-d matrices up to n_l;
% uses 3-term recurrences. ;
% test with: ;
%{

  % plots should show matrices of all ones ;
  n_l = 5;
  W1_ = wignerd_lsq_b(n_l,pi/6); W2_ = wignerd_b(n_l,pi/6); 
  for nl = 1:n_l;
  subplot(2,n_l,nl+0*n_l); imagesc(sign((W2_{nl})./real(W1_{nl})),[-1,1]);
  subplot(2,n_l,nl+1*n_l); imagesc((W2_{nl})./real(W1_{nl}),[-3,3]);
  end;%for nl = 1:n_l;


  %}

verbose=0;

cb = cos(beta/2); sb = sin(beta/2); 
W_ = cell(1+n_l,1);
W_{1} = [1];
for nl=1:n_l;
if (verbose); disp(sprintf(' \n %% nl %d \n',nl)); end;
nlp = nl-1;
V = W_{nl};
W = zeros(2*nl+1,2*nl+1);
for nmn = -nl:+nl; for nmp = -nl:+nl;
%%%%%%%%%%%%%%%%;
if (nl~=-nmp & nl~=1-nmp); % use recurrence A ;
str_tmp = 'A';
W1 = 0; if (abs(nmp-1)<=nlp & abs(nmn-1)<=nlp); W1 = V(1+nlp+(nmp-1),1+nlp+(nmn-1)); end;
W2 = 0; if (abs(nmp-1)<=nlp & abs(nmn-0)<=nlp); W2 = V(1+nlp+(nmp-1),1+nlp+(nmn-0)); end;
W3 = 0; if (abs(nmp-1)<=nlp & abs(nmn+1)<=nlp); W3 = V(1+nlp+(nmp-1),1+nlp+(nmn+1)); end;
tmp = cb*cb*A1(nl,nmn,nmp)*W1 - 2*cb*sb*A2(nl,nmn,nmp)*W2 + sb*sb*A3(nl,nmn,nmp)*W3;
end;%if (nl~=-nmp & nl~=1-nmp); % use recurrence A ;
%%%%%%%%%%%%%%%%;
if (nl~=+nmp & nl~=1+nmp); % use recurrence B ;
str_tmp = 'B';
W1 = 0; if (abs(nmp+1)<=nlp & abs(nmn-1)<=nlp); W1 = V(1+nlp+(nmp+1),1+nlp+(nmn-1)); end;
W2 = 0; if (abs(nmp+1)<=nlp & abs(nmn-0)<=nlp); W2 = V(1+nlp+(nmp+1),1+nlp+(nmn-0)); end;
W3 = 0; if (abs(nmp+1)<=nlp & abs(nmn+1)<=nlp); W3 = V(1+nlp+(nmp+1),1+nlp+(nmn+1)); end;
tmp = sb*sb*B1(nl,nmn,nmp)*W1 + 2*sb*cb*B2(nl,nmn,nmp)*W2 + cb*cb*B3(nl,nmn,nmp)*W3;
end;%if (nl~=+nmp & nl~=1+nmp); % use recurrence B ;
%%%%%%%%%%%%%%%%;
if (nl~=-nmp & nl~=+nmp); % use recurrence C ;
str_tmp = 'C';
W1 = 0; if (abs(nmp-0)<=nlp & abs(nmn-1)<=nlp); W1 = V(1+nlp+(nmp-0),1+nlp+(nmn-1)); end;
W2 = 0; if (abs(nmp-0)<=nlp & abs(nmn-0)<=nlp); W2 = V(1+nlp+(nmp-0),1+nlp+(nmn-0)); end;
W3 = 0; if (abs(nmp-0)<=nlp & abs(nmn+1)<=nlp); W3 = V(1+nlp+(nmp-0),1+nlp+(nmn+1)); end;
tmp = sb*cb*C1(nl,nmn,nmp)*W1 + (cb*cb-sb*sb)*C2(nl,nmn,nmp)*W2 - sb*cb*C3(nl,nmn,nmp)*W3;
end;%if (nl~=-nmp & nl~=+nmp); % use recurrence C ;
%%%%%%%%%%%%%%%%;
if (verbose); disp(sprintf(' %% setting W(%d,%d) <-- %f (using %s %f %f %f)',nl+nmp,nl+nmn,tmp,str_tmp,W1,W2,W3)); end;
W(1+nl+nmp,1+nl+nmn) = tmp;
end;end;%for nmn = -nl:+nl; for nmp = -nl:+nl;
W_{1+nl} = W;
end;%for nl=1:n_l;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Fix condon-shortley-phase ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for nl=1:n_l;
m_ = -nl:+nl;
s=(-1).^((m_<0).*m_); % needed to preserve condon-shortley phase. ;
S = transpose(s)*s;
W_{1+nl} = W_{1+nl}.*S;
end;%for nl=1:n_l;

%{
for nl=1:n_l;
for nmn = -nl:+nl; for nmp = -nl:+nl;
smn = 1; if (nmn<0 & mod(nmn,2)==1); smn = -1; end;
smp = 1; if (nmp<0 & mod(nmp,2)==1); smp = -1; end;
W_{1+nl}(1+nl+nmp,1+nl+nmn) = W_{1+nl}(1+nl+nmp,1+nl+nmn)*smn*smp;
end;end; %for nmn = -nl:+nl; for nmp = -nl:+nl;
end;%for nl=1:n_l;
 %}

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
