%{
dx=3;dy=6;dz=9;
%%%%%%%%;
%x1 = [ 0;0];x2 = [dx;0];x3 = [dx;1];x4 = [ 0;1]; 
%y1 = [dx;0];y2 = [dy;0];y3 = [dy;1];y4 = [dx;1]; 
%z1 = [dy;0];z2 = [dz;0];z3 = [dz;1];z4 = [dy;1]; 
%%%%%%%%;
%x1 = [ 0;-1];x2 = [dx;-1];x3 = [dx;0];x4 = [ 0;0]; 
%y1 = [dx;-1];y2 = [dy;-1];y3 = [dy;0];y4 = [dx;0]; 
%z1 = [dy;-1];z2 = [dz;-1];z3 = [dz;0];z4 = [dy;0]; 
%%%%%%%%;
x1 = [ 0+1;-1];x2 = [dx+1;-1];x3 = [dx+1;0];x4 = [ 0+1;0]; 
y1 = [dx+1;-1];y2 = [dy+1;-1];y3 = [dy+1;0];y4 = [dx+1;0]; 
z1 = [dy+1;-1];z2 = [dz+1;-1];z3 = [dz+1;0];z4 = [dy+1;0]; 
n_w = 64;
hold on;
for nw=0:n_w-1;
w = (nw)/(n_w-1)*pi/2; c=cos(w);s=sin(w); R=[c,-s;+s,c];
X1=R*x1; X2=R*x2; X3=R*x3; X4=R*x4;
Y1=R*y1; Y2=R*y2; Y3=R*y3; Y4=R*y4;
Z1=R*z1; Z2=R*z2; Z3=R*z3; Z4=R*z4;
patch([X1(1),X2(1),X3(1),X4(1)],[X1(2),X2(2),X3(2),X4(2)],[0,0,0],'edgecolor','none');
patch([Z1(1),Z2(1),Z3(1),Z4(1)],[Z1(2),Z2(2),Z3(2),Z4(2)],[0,0,0],'edgecolor','none');
patch([Y1(1),Y2(1),Y3(1),Y4(1)],[Y1(2),Y2(2),Y3(2),Y4(2)],1-[1,1,1]*nw/(n_w-1),'edgecolor','none');
end;%for nw=0:n_w-1;
hold off;
axis equal;
 %}

%{
dx=4;dy=1.0;
x1 = [-dx;-dy];x2 = [+dx;-dy];x3 = [+dx;+dy];x4 = [-dx;+dy]; 
n_w = 64;
hold on;
for nw=0:n_w-1;
w = (nw)/(n_w-1)*pi/2; c=cos(w);s=sin(w); R=[c,-s;+s,c];
X1=R*x1; X2=R*x2; X3=R*x3; X4=R*x4;
vshift = (dx-2*dy)*nw/(n_w-1);
cshift = nw/(n_w-1);
patch([X1(1),X2(1),X3(1),X4(1)],vshift + [X1(2),X2(2),X3(2),X4(2)],1-cshift*[1,1,1],'edgecolor','none');
end;%for nw=0:n_w-1;
hold off;
axis equal;
 %}

%{
dx=4;dy=1.0;
x1 = [dx-dx;-dy];x2 = [dx+dx;-dy];x3 = [dx+dx;+dy];x4 = [dx-dx;+dy]; 
n_w = 64;
hold on;
for nw=0:n_w-1;
w = (nw)/(n_w-1)*pi/2; c=cos(w);s=sin(w); R=[c,-s;+s,c];
X1=R*x1; X2=R*x2; X3=R*x3; X4=R*x4;
vshift = (dx-2*dy)*nw/(n_w-1);
cshift = nw/(n_w-1);
patch([X1(1),X2(1),X3(1),X4(1)],vshift + [X1(2),X2(2),X3(2),X4(2)],1-cshift*[1,1,1],'edgecolor','none');
end;%for nw=0:n_w-1;
hold off;
axis equal;
 %}

%{
dx=4;dy=1.0;
%x1 = [dx-dx;-dy];x2 = [dx+dx;-dy];x3 = [dx+dx;+dy];x4 = [dx-dx;+dy]; 
x1 = [-dx;-dy];x2 = [+dx;-dy];x3 = [+dx;+dy];x4 = [-dx;+dy]; 
n_w = 64;
hold on;
for nw=0:n_w-1;
w = (nw)/(n_w-1)*pi/2; c=cos(w);s=sin(w); R=[c,-s;+s,c];
X1=R*x1; X2=R*x2; X3=R*x3; X4=R*x4;
hshift = 1*(dy-2*dx)*nw/(n_w-1);
vshift = 0*(dx-2*dy)*nw/(n_w-1);
cshift = nw/(n_w-1);
patch(hshift + [X1(1),X2(1),X3(1),X4(1)],vshift + [X1(2),X2(2),X3(2),X4(2)],1-cshift*[1,1,1],'edgecolor','none');
end;%for nw=0:n_w-1;
hold off;
axis equal;
 %}

%{
ra = 8; rb = 10;
n_w = 64;
hold on;
for nw=0:n_w-1;
w = (nw)/(n_w-1)*pi/2; dw = pi/2/(n_w-1);
c0 = cos(w); c1 = cos(w+dw); s0 = sin(w); s1 = sin(w+dw);
cshift = nw/(n_w-1);
patch([ra*c0,rb*c0,rb*c1,ra*c1],[ra*s0,rb*s0,rb*s1,ra*s1],0.85-0.85*cshift*[1,1,1],'edgecolor','none');
end;%for nw=0:n_w-1;
hold off;
axis equal;
 %}

r1 =  8; r2 =  9;
r3 = 10; r4 = 11;
r5 = 12; r6 = 13;
n_w = 64;
hold on;
for nw=0:n_w-1;
w = (nw)/(n_w-1)*pi/2; dw = pi/2/(n_w-1);
c0 = cos(w); c1 = cos(w+dw); s0 = sin(w); s1 = sin(w+dw);
cshift = nw/(n_w-1);
patch([r1*c0,r2*c0,r2*c1,r1*c1],[r1*s0,r2*s0,r2*s1,r1*s1],0.85-0.85*cshift*[1,1,1],'edgecolor','none');
patch([r3*c0,r4*c0,r4*c1,r3*c1],[r3*s0,r4*s0,r4*s1,r3*s1],0.85-0.85*cshift*[1,1,1],'edgecolor','none');
patch([r5*c0,r6*c0,r6*c1,r5*c1],[r5*s0,r6*s0,r6*s1,r5*s1],0.85-0.85*cshift*[1,1,1],'edgecolor','none');
end;%for nw=0:n_w-1;
hold off;
axis equal;
