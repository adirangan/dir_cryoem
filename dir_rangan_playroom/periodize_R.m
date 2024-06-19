function [bac_] = periodize_R(bac_);
b = bac_(1+0); a = bac_(1+1); c = bac_(1+2);
b = periodize(b, 0*pi, 2*pi);
a = periodize(a,-1*pi,+1*pi);
c = periodize(c, 0*pi, 2*pi);
if a<0; a = -a; b = periodize(b-pi,0,2*pi); c = periodize(c-pi,0,2*pi); end;
bac_ = [b,a,c];
