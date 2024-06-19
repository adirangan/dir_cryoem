function [e_,ier] = R_to_euler(R);
% converts rotation-matrix R (with determinant +1) into euler-angles e_ = [alpha,beta,gamma] ;
% ;
% test with: ;
%{
R_to_euler();
 %}

if (nargin<1);
e__in_ = [2*pi*rand(),1*pi*rand(),2*pi*rand()];
R = euler_to_R(e__in_);
continue_flag=1; iteration = 0; iteration_max = 4;
while (continue_flag & iteration<iteration_max);
[e_out_,ier] = R_to_euler(R);
continue_flag = ier;
iteration = iteration + 1;
end;% while;
e_out_ = periodize(e_out_,0,2*pi); 
disp(sprintf(' %% angle error %0.6f %0.6f %0.6f',e__in_-e_out_));
disp(sprintf(' %% R error: %0.6f',norm(euler_to_R(e__in_) - euler_to_R(e_out_))));
disp(euler_to_R(e__in_) - euler_to_R(e_out_));
return;
end;% if (nargin<1);

verbose=0;

e_ = [0,0,0];
iteration = 0;
iteration_max = 1024*8;
D = R - euler_to_R(e_);
E = trace(transpose(D)*D);

eta = 0.25; % learning rate;
continue_flag=1; ier=0;
while (continue_flag);
a = e_(1);b = e_(2);g = e_(3);
da = -trace(transpose(D) * Rz(g)*Ry(b)*dRz(a));
db = -trace(transpose(D) * Rz(g)*dRy(b)*Rz(a));
dg = -trace(transpose(D) * dRz(g)*Ry(b)*Rz(a));
a = a-da*eta;
b = b-db*eta;
g = g-dg*eta;
a = periodize(a,0,2*pi);
b = periodize(b,-1*pi,1*pi);
g = periodize(g,0,2*pi);
if b<0; b = -b; a = periodize(a-pi,0,2*pi); g = periodize(g-pi,0,2*pi); end;
e_ = [a,b,g];
D = R - euler_to_R(e_);
E = trace(transpose(D)*D);
if (verbose); disp(sprintf(' %% iteration %d/%d: [da,db,dg] = [%0.2f,%0.2f,%0.2f] --> E %0.16f',iteration,iteration_max,da,db,dg,E)); end;
continue_flag = ((iteration<iteration_max) && (sqrt(E)>1e-12));
iteration = iteration+1;
end;%while;
if (iteration>=iteration_max); disp(sprintf(' %% Warning! iteration %d/%d in R_to_euler; E %0.16f',iteration,iteration_max,E)); ier=1; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
function R = Rz(w);
cw = cos(w); sw = sin(w);
R = [+cw -sw   0 ; +sw +cw   0 ;   0   0   1 ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
function R = dRz(w);
cw = cos(w); sw = sin(w);
R = [-sw -cw   0 ; +cw -sw   0 ;   0   0   0 ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
function R = Ry(w);
cw = cos(w); sw = sin(w);
R = [+cw   0 +sw ;   0   1   0 ; -sw   0 +cw ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
function R = dRy(w);
cw = cos(w); sw = sin(w);
R = [-sw   0 +cw ;   0   0   0 ; -cw   0 -sw ];

