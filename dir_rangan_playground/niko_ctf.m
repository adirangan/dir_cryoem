function output = niko_ctf(Cs,lambda,w1,w2,df1,df2,angast,thetatr,l,m);
rad = l^2+m^2; rad = dsqrt(rad);
angle = rad*thetatr;
angspt=datan2(m,l);
c1=2.0d0*pi*angle*angle/(2.0d0*lambda);
c2=-c1*Cs*angle*angle/2.0d0;
angdif=angspt-angast;
ccos=dcos(2.0d0*angdif);
df = 0.5d0*(df1+df2+ccos*(df1-df2));
chi=c1*df+c2;
ctfv=-w1*dsin(chi)-w2*dcos(chi);
output = ctfv;
