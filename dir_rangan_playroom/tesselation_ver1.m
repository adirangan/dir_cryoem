function T_ = tesselation_ver1();

rng(0);
n_points = 1024;
L_ = randn(3,n_points); for npoints=1:n_points; L_(:,npoints) = L_(:,npoints)/norm(L_(:,npoints)); end;

v_n00 = [-1; 0; 0];
v_0n0 = [ 0;-1; 0];
v_00n = [ 0; 0;-1];
v_p00 = [+1; 0; 0];
v_0p0 = [ 0;+1; 0];
v_00p = [ 0; 0;+1];

nl_nnn = tesselation_find_depth(0,v_n00,v_0n0,v_00n,0,L_);
nl_nnp = tesselation_find_depth(0,v_n00,v_0n0,v_00p,1,L_);
nl_npn = tesselation_find_depth(0,v_n00,v_0p0,v_00n,1,L_);
nl_npp = tesselation_find_depth(0,v_n00,v_0p0,v_00p,0,L_);
nl_pnn = tesselation_find_depth(0,v_p00,v_0n0,v_00n,1,L_);
nl_pnp = tesselation_find_depth(0,v_p00,v_0n0,v_00p,0,L_);
nl_ppn = tesselation_find_depth(0,v_p00,v_0p0,v_00n,0,L_);
nl_ppp = tesselation_find_depth(0,v_p00,v_0p0,v_00p,1,L_);
nl_= [nl_nnn,nl_nnp,nl_npn,nl_npp,nl_pnn,nl_pnp,nl_ppn,nl_ppp];
disp(num2str(nl_));
nl_max = max(nl_);
disp(sprintf(' %% n_points %d --> nl_max %d',n_points,nl_max));

nl_max = 4;
nl=0;
T_nnn = tesselation_make(nl_max,nl,v_n00,v_0n0,v_00n,0);
T_nnp = tesselation_make(nl_max,nl,v_n00,v_0n0,v_00p,1);
T_npn = tesselation_make(nl_max,nl,v_n00,v_0p0,v_00n,1);
T_npp = tesselation_make(nl_max,nl,v_n00,v_0p0,v_00p,0);
T_pnn = tesselation_make(nl_max,nl,v_p00,v_0n0,v_00n,1);
T_pnp = tesselation_make(nl_max,nl,v_p00,v_0n0,v_00p,0);
T_ppn = tesselation_make(nl_max,nl,v_p00,v_0p0,v_00n,0);
T_ppp = tesselation_make(nl_max,nl,v_p00,v_0p0,v_00p,1);
T_ = {T_nnn,T_nnp,T_npn,T_npp,T_pnn,T_pnp,T_ppn,T_ppp};

plot_flag=0;
if plot_flag;
figure; clf; hold on;
colormap('jet');
nl = 3; tesselation_plot_parity(nl_max,nl,T_);
hold off;
xlabel('x'); ylabel('y'); zlabel('z');
axis equal;axis vis3d;
end;%if plot_flag;

plot_flag=0;
if plot_flag;
figure; clf; hold on;
colormap('jet');
nl = 3; tesselation_plot(nl_max,nl,T_);
hold off;
xlabel('x'); ylabel('y'); zlabel('z');
axis equal;axis vis3d;
end;%if plot_flag;

plot_flag=0;
if plot_flag;
vp = randn(3,1); vp = vp/norm(vp);
figure; clf; hold on;
colormap('jet');
plot3(vp(1),vp(2),vp(3),'y.','MarkerSize',25);
nl = 1; tesselation_plot_point(nl_max,nl,T_,vp);
hold off;
xlabel('x'); ylabel('y'); zlabel('z');
axis equal;axis vis3d;
end;%if plot_flag;

plot_flag=0;
if plot_flag;
f = @(polar_angle) 1*cos(1*polar_angle) - 0.5*cos(3*polar_angle) + 0.25*cos(5*polar_angle) - 0.125*cos(7*polar_angle) + 1.5*cos(9*polar_angle);
g = @(azimuthal_angle) 1*sin(1*azimuthal_angle) - 0.5*sin(3*azimuthal_angle) + 0.25*sin(5*azimuthal_angle) - 0.125*sin(7*azimuthal_angle) + 1.5*sin(9*azimuthal_angle);
hpa = @(polar_angle,azimuthal_angle) f(polar_angle).*g(azimuthal_angle);
cra = colormap('jet'); clim = 4*[-1,+1];
figure; clf; hold on;
nl = 4; tesselation_plot_hpa(nl_max,nl,T_,hpa,cra,clim);
hold off;
xlabel('x'); ylabel('y'); zlabel('z');
axis equal;axis vis3d;
end;%if plot_flag;

plot_flag=0;
if plot_flag;
hxyz = @(x,y,z) (sin(x) + y.^3 + exp(z)).*cos(2*pi*3*x);
cra = colormap('jet'); clim = [-1,+4];
figure; clf; hold on;
nl = 4; tesselation_plot_hxyz(nl_max,nl,T_,hxyz,cra,clim);
hold off;
xlabel('x'); ylabel('y'); zlabel('z');
axis equal;axis vis3d;
end;%if plot_flag;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tesselation_plot_hxyz(nl_max,nl,T,hxyz,cra,clim);
if (length(T)>1); for (nt=1:length(T)); tesselation_plot_hxyz(nl_max,nl,T{nt},hxyz,cra,clim); end; return; end;
verbose=0;
prefix = 37*ones(1,2*T.nl);
if (verbose); disp(sprintf(' %s nl_max %d nl %d T.nl %d',prefix,nl_max,nl,T.nl)); end;
if T.nl<nl; 
if ~isempty(T.child_1); if (verbose); disp(sprintf(' %s calling child_1',prefix)); end; tesselation_plot_hxyz(nl_max,nl,T.child_1,hxyz,cra,clim); end;
if ~isempty(T.child_2); if (verbose); disp(sprintf(' %s calling child_2',prefix)); end; tesselation_plot_hxyz(nl_max,nl,T.child_2,hxyz,cra,clim); end;
if ~isempty(T.child_3); if (verbose); disp(sprintf(' %s calling child_3',prefix)); end; tesselation_plot_hxyz(nl_max,nl,T.child_3,hxyz,cra,clim); end;
if ~isempty(T.child_4); if (verbose); disp(sprintf(' %s calling child_4',prefix)); end; tesselation_plot_hxyz(nl_max,nl,T.child_4,hxyz,cra,clim); end;
end;%if T.nl<nl; 
if (T.nl==nl);
if (verbose); disp(sprintf(' %s plotting self',prefix)); end;
v1 = T.vm + 0.9*(T.v1-T.vm);
v2 = T.vm + 0.9*(T.v2-T.vm);
v3 = T.vm + 0.9*(T.v3-T.vm);
vm = T.vm; 
vmh = hxyz(vm(1),vm(2),vm(3));
cbin = max(1,min(size(cra,1),1+floor(size(cra,1)*(vmh-min(clim))/diff(clim))));
x = [v1(1),v2(1),v3(1)];y = [v1(2),v2(2),v3(2)];z = [v1(3),v2(3),v3(3)];
patch(x,y,z,cra(cbin,:));
end;%if (T.nl==nl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tesselation_plot_hpa(nl_max,nl,T,hpa,cra,clim);
if (length(T)>1); for (nt=1:length(T)); tesselation_plot_hpa(nl_max,nl,T{nt},hpa,cra,clim); end; return; end;
verbose=0;
prefix = 37*ones(1,2*T.nl);
if (verbose); disp(sprintf(' %s nl_max %d nl %d T.nl %d',prefix,nl_max,nl,T.nl)); end;
if T.nl<nl; 
if ~isempty(T.child_1); if (verbose); disp(sprintf(' %s calling child_1',prefix)); end; tesselation_plot_hpa(nl_max,nl,T.child_1,hpa,cra,clim); end;
if ~isempty(T.child_2); if (verbose); disp(sprintf(' %s calling child_2',prefix)); end; tesselation_plot_hpa(nl_max,nl,T.child_2,hpa,cra,clim); end;
if ~isempty(T.child_3); if (verbose); disp(sprintf(' %s calling child_3',prefix)); end; tesselation_plot_hpa(nl_max,nl,T.child_3,hpa,cra,clim); end;
if ~isempty(T.child_4); if (verbose); disp(sprintf(' %s calling child_4',prefix)); end; tesselation_plot_hpa(nl_max,nl,T.child_4,hpa,cra,clim); end;
end;%if T.nl<nl; 
if (T.nl==nl);
if (verbose); disp(sprintf(' %s plotting self',prefix)); end;
v1 = T.vm + 0.9*(T.v1-T.vm);
v2 = T.vm + 0.9*(T.v2-T.vm);
v3 = T.vm + 0.9*(T.v3-T.vm);
vm = T.vm; 
polar_angle = atan2(vm(2),vm(1)); azimuthal_angle = acos(vm(3)/norm(vm)); %azimuthal_angle = atan2(sqrt(vm(1)*vm(1)+vm(2)*vm(2)),vm(3));
vmh = hpa(polar_angle,azimuthal_angle);
cbin = max(1,min(size(cra,1),1+floor(size(cra,1)*(vmh-min(clim))/diff(clim))));
x = [v1(1),v2(1),v3(1)];y = [v1(2),v2(2),v3(2)];z = [v1(3),v2(3),v3(3)];
patch(x,y,z,cra(cbin,:));
end;%if (T.nl==nl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tesselation_plot_point(nl_max,nl,T,vp);
if (length(T)>1); for (nt=1:length(T)); tesselation_plot_point(nl_max,nl,T{nt},vp); end; return; end;
verbose=0;
prefix = 37*ones(1,2*T.nl);
if (verbose); disp(sprintf(' %s nl_max %d nl %d T.nl %d',prefix,nl_max,nl,T.nl)); end;
if T.nl<nl; 
if ~isempty(T.child_1); if (verbose); disp(sprintf(' %s calling child_1',prefix)); end; tesselation_plot_point(nl_max,nl,T.child_1,vp); end;
if ~isempty(T.child_2); if (verbose); disp(sprintf(' %s calling child_2',prefix)); end; tesselation_plot_point(nl_max,nl,T.child_2,vp); end;
if ~isempty(T.child_3); if (verbose); disp(sprintf(' %s calling child_3',prefix)); end; tesselation_plot_point(nl_max,nl,T.child_3,vp); end;
if ~isempty(T.child_4); if (verbose); disp(sprintf(' %s calling child_4',prefix)); end; tesselation_plot_point(nl_max,nl,T.child_4,vp); end;
end;%if T.nl<nl; 
if (T.nl==nl);
if (verbose); disp(sprintf(' %s plotting self',prefix)); end;
v1 = T.vm + 0.9*(T.v1-T.vm);
v2 = T.vm + 0.9*(T.v2-T.vm);
v3 = T.vm + 0.9*(T.v3-T.vm);
vm = T.vm; el = norm(v1-v2)/4;
vn = vm + T.nn*el;
x = [v1(1),v2(1),v3(1)];y = [v1(2),v2(2),v3(2)];z = [v1(3),v2(3),v3(3)];
if tesselation_in(T,vp); c_flag = 0; else; c_flag = 1; end;
patch(x,y,z,c_flag);
end;%if (T.nl==nl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tesselation_plot(nl_max,nl,T);
if (length(T)>1); for (nt=1:length(T)); tesselation_plot(nl_max,nl,T{nt}); end; return; end;
verbose=0;
prefix = 37*ones(1,2*T.nl);
if (verbose); disp(sprintf(' %s nl_max %d nl %d T.nl %d',prefix,nl_max,nl,T.nl)); end;
if T.nl<nl; 
if ~isempty(T.child_1); if (verbose); disp(sprintf(' %s calling child_1',prefix)); end; tesselation_plot(nl_max,nl,T.child_1); end;
if ~isempty(T.child_2); if (verbose); disp(sprintf(' %s calling child_2',prefix)); end; tesselation_plot(nl_max,nl,T.child_2); end;
if ~isempty(T.child_3); if (verbose); disp(sprintf(' %s calling child_3',prefix)); end; tesselation_plot(nl_max,nl,T.child_3); end;
if ~isempty(T.child_4); if (verbose); disp(sprintf(' %s calling child_4',prefix)); end; tesselation_plot(nl_max,nl,T.child_4); end;
end;%if T.nl<nl; 
if (T.nl==nl);
if (verbose); disp(sprintf(' %s plotting self',prefix)); end;
v1 = T.vm + 0.9*(T.v1-T.vm);
v2 = T.vm + 0.9*(T.v2-T.vm);
v3 = T.vm + 0.9*(T.v3-T.vm);
vm = T.vm; el = norm(v1-v2)/4;
vn = vm + T.nn*el;
x = [v1(1),v2(1),v3(1)];y = [v1(2),v2(2),v3(2)];z = [v1(3),v2(3),v3(3)];
patch(x,y,z,[0;1;2]/2,'FaceColor','Interp');
x = [vm(1),vn(1)];y = [vm(2),vn(2)];z = [vm(3),vn(3)];
line(x,y,z,'LineWidth',1,'Color',[0,0,0]);
end;%if (T.nl==nl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tesselation_plot_parity(nl_max,nl,T);
if (length(T)>1); for (nt=1:length(T)); tesselation_plot_parity(nl_max,nl,T{nt}); end; return; end;
verbose=0;
prefix = 37*ones(1,2*T.nl);
if (verbose); disp(sprintf(' %s nl_max %d nl %d T.nl %d',prefix,nl_max,nl,T.nl)); end;
if T.nl<nl; 
if ~isempty(T.child_1); if (verbose); disp(sprintf(' %s calling child_1',prefix)); end; tesselation_plot_parity(nl_max,nl,T.child_1); end;
if ~isempty(T.child_2); if (verbose); disp(sprintf(' %s calling child_2',prefix)); end; tesselation_plot_parity(nl_max,nl,T.child_2); end;
if ~isempty(T.child_3); if (verbose); disp(sprintf(' %s calling child_3',prefix)); end; tesselation_plot_parity(nl_max,nl,T.child_3); end;
if ~isempty(T.child_4); if (verbose); disp(sprintf(' %s calling child_4',prefix)); end; tesselation_plot_parity(nl_max,nl,T.child_4); end;
end;%if T.nl<nl; 
if (T.nl==nl);
if (verbose); disp(sprintf(' %s plotting self',prefix)); end;
v1 = T.vm + 0.9*(T.v1-T.vm);
v2 = T.vm + 0.9*(T.v2-T.vm);
v3 = T.vm + 0.9*(T.v3-T.vm);
vm = T.vm; el = norm(v1-v2)/4;
vn = vm + T.nn*el;
x = [v1(1),v2(1),v3(1)];y = [v1(2),v2(2),v3(2)];z = [v1(3),v2(3),v3(3)];
patch(x,y,z,cast(T.parity,'double'));
x = [vm(1),vn(1)];y = [vm(2),vn(2)];z = [vm(3),vn(3)];
line(x,y,z,'LineWidth',1,'Color',[0,0,0]);
end;%if (T.nl==nl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = tesselation_find_depth(nl,v1_i,v2_i,v3_i,parity,L_);
output = nl;
v1_i = v1_i/norm(v1_i); v2_i = v2_i/norm(v2_i); v3_i = v3_i/norm(v3_i);
v1 = v1_i; v2 = v2_i; v3 = v3_i; 
vm = (v1+v2+v3)/3;
nn = cross(v3-v2,v1-v3);
aa = dot(vm,nn);
if (aa>=0); v1 = v1_i; v2 = v2_i; v3 = v3_i; end;
if (aa< 0); v1 = v2_i; v2 = v1_i; v3 = v3_i; end;
e2 = v1-v3; e3 = v2-v1; e1 = v3-v2;
m1 = (v2+v3)/2; m1 = m1/norm(m1);
m2 = (v3+v1)/2; m2 = m2/norm(m2);
m3 = (v1+v2)/2; m3 = m3/norm(m3);
nn = cross(e1,e2); nn = nn/norm(nn);
n2 = cross(v3,e2); n2 = n2/norm(n2);
n3 = cross(v1,e3); n3 = n3/norm(n3);
n1 = cross(v2,e1); n1 = n1/norm(n1);
n_L_length = size(L_,2);
I1 = zeros(n_L_length,1); for nL_length=1:n_L_length; I1(nL_length) = dot(n1,L_(:,nL_length)); end;%for nL_length=1:n_L_length;
I2 = zeros(n_L_length,1); for nL_length=1:n_L_length; I2(nL_length) = dot(n2,L_(:,nL_length)); end;%for nL_length=1:n_L_length;
I3 = zeros(n_L_length,1); for nL_length=1:n_L_length; I3(nL_length) = dot(n3,L_(:,nL_length)); end;%for nL_length=1:n_L_length;
if  parity; IA = zeros(n_L_length,1); IA = (I1>=0) & (I2>=0) & (I3>=0); end;
if ~parity; IA = zeros(n_L_length,1); IA = (I1> 0) & (I2> 0) & (I3> 0); end;
L_sub_ = L_(:,find(IA));
n_L_sub_length = size(L_sub_,2);
if (n_L_sub_length==0);
output = -1;
elseif (n_L_sub_length==1);
output = nl;
elseif (n_L_sub_length>1);
nl_child_1 = tesselation_find_depth(nl+1,v1,m3,m2, parity,L_sub_);
nl_child_2 = tesselation_find_depth(nl+1,v2,m1,m3, parity,L_sub_);
nl_child_3 = tesselation_find_depth(nl+1,v3,m2,m1, parity,L_sub_);
nl_child_4 = tesselation_find_depth(nl+1,m1,m2,m3,~parity,L_sub_);
output = max([nl_child_1,nl_child_2,nl_child_3,nl_child_4]);
end;% if (n_L_sub_length>1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T = tesselation_make(nl_max,nl,v1_i,v2_i,v3_i,parity);
v1_i = v1_i/norm(v1_i); v2_i = v2_i/norm(v2_i); v3_i = v3_i/norm(v3_i);
v1 = v1_i; v2 = v2_i; v3 = v3_i; 
vm = (v1+v2+v3)/3;
nn = cross(v3-v2,v1-v3);
aa = dot(vm,nn);
if (aa>=0); v1 = v1_i; v2 = v2_i; v3 = v3_i; end;
if (aa< 0); v1 = v2_i; v2 = v1_i; v3 = v3_i; end;
e2 = v1-v3; e3 = v2-v1; e1 = v3-v2;
m1 = (v2+v3)/2; m1 = m1/norm(m1);
m2 = (v3+v1)/2; m2 = m2/norm(m2);
m3 = (v1+v2)/2; m3 = m3/norm(m3);
nn = cross(e1,e2); nn = nn/norm(nn);
n2 = cross(v3,e2); n2 = n2/norm(n2);
n3 = cross(v1,e3); n3 = n3/norm(n3);
n1 = cross(v2,e1); n1 = n1/norm(n1);
child_1 = struct([]);
child_2 = struct([]);
child_3 = struct([]);
child_4 = struct([]);
if (nl>=nl_max);
% do nothing;
end;%if (nl>=nl_max);
if (nl<nl_max);
child_1 = tesselation_make(nl_max,nl+1,v1,m3,m2, parity);
child_2 = tesselation_make(nl_max,nl+1,v2,m1,m3, parity);
child_3 = tesselation_make(nl_max,nl+1,v3,m2,m1, parity);
child_4 = tesselation_make(nl_max,nl+1,m1,m2,m3,~parity);
end;%if (nl<nl_max);
T = struct(...
	   'nl',nl,...
	   'parity',parity,...
	   'v1',v1,...
	   'v2',v2,...
	   'v3',v3,...
	   'vm',vm,...
	   'm1',m1,...
	   'm2',m2,...
	   'm3',m3,...
	   'e1',e1,...
	   'e2',e2,...
	   'e3',e3,...
	   'n1',n1,...
	   'n2',n2,...
	   'n3',n3,...
	   'nn',nn,...
	   'child_1',child_1,...
	   'child_2',child_2,...
	   'child_3',child_3,...
	   'child_4',child_4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = tesselation_in(T,vp);
output = 0;
if  T.parity; output = (dot(vp,T.n1)>=0 & dot(vp,T.n2)>=0 & dot(vp,T.n3)>=0); end;
if ~T.parity; output = (dot(vp,T.n1)> 0 & dot(vp,T.n2)> 0 & dot(vp,T.n3)> 0); end;

