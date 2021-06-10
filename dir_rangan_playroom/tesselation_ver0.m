function T_ = tesselation_ver0();
v_n00 = [-1; 0; 0];
v_0n0 = [ 0;-1; 0];
v_00n = [ 0; 0;-1];
v_p00 = [+1; 0; 0];
v_0p0 = [ 0;+1; 0];
v_00p = [ 0; 0;+1];

nl_max = 3;
nl=0;
T_nnn = tesselation_make(nl_max,nl,v_n00,v_0n0,v_00n);
T_nnp = tesselation_make(nl_max,nl,v_n00,v_0n0,v_00p);
T_npn = tesselation_make(nl_max,nl,v_n00,v_0p0,v_00n);
T_npp = tesselation_make(nl_max,nl,v_n00,v_0p0,v_00p);
T_pnn = tesselation_make(nl_max,nl,v_p00,v_0n0,v_00n);
T_pnp = tesselation_make(nl_max,nl,v_p00,v_0n0,v_00p);
T_ppn = tesselation_make(nl_max,nl,v_p00,v_0p0,v_00n);
T_ppp = tesselation_make(nl_max,nl,v_p00,v_0p0,v_00p);
T_ = {T_nnn,T_nnp,T_npn,T_npp,T_pnn,T_pnp,T_ppn,T_ppp};

plot_flag=1;
if plot_flag;
figure; clf; hold on;
colormap('jet');
nl = 3; for nt=1:length(T_); tesselation_plot(nl_max,nl,T_{nt}); end%for nt=1:length(T_);
hold off;
axis equal;axis vis3d;
end;%if plot_flag;

plot_flag=0;
if plot_flag;
vp = randn(3,1); vp = vp/norm(vp);
figure; clf; hold on;
colormap('jet');
plot3(vp(1),vp(2),vp(3),'y.','MarkerSize',25);
nl = 1; for nt=1:length(T_); tesselation_plot_point(nl_max,nl,T_{nt},vp); end%for nt=1:length(T_);
hold off;
axis equal;axis vis3d;
end;%if plot_flag;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tesselation_plot_point(nl_max,nl,T,vp);
verbose=0;
prefix = 37*ones(1,2*T.nl);
if (verbose); disp(sprintf(' %s nl_max %d nl %d T.nl %d',prefix,nl_max,nl,T.nl)); end;
if T.nl<nl; 
if ~isempty(T.child_1); if (verbose); disp(sprintf(' %s calling child_1',prefix)); end; tesselation_plot_point(nl_max,nl,T.child_1,vp); end;
if ~isempty(T.child_2); if (verbose); disp(sprintf(' %s calling child_2',prefix)); end; tesselation_plot_point(nl_max,nl,T.child_2,vp); end;
if ~isempty(T.child_3); if (verbose); disp(sprintf(' %s calling child_3',prefix)); end; tesselation_plot_point(nl_max,nl,T.child_3,vp); end;
if ~isempty(T.child_4); if (verbose); disp(sprintf(' %s calling child_4',prefix)); end; tesselation_plot_point(nl_max,nl,T.child_4,vp); end;
end;%if T.nl<nl_max; 
if (T.nl==nl);
if (verbose); disp(sprintf(' %s plotting self',prefix)); end;
v1 = T.vm + 0.9*(T.v1-T.vm);
v2 = T.vm + 0.9*(T.v2-T.vm);
v3 = T.vm + 0.9*(T.v3-T.vm);
vm = T.vm; el = norm(v1-v2)/4;
vn = vm + T.nn*el;
x = [v1(1),v2(1),v3(1)];y = [v1(2),v2(2),v3(2)];z = [v1(3),v2(3),v3(3)];
if (dot(vp,T.n1)>0 & dot(vp,T.n2)>0 & dot(vp,T.n3)>0); c_flag = 0; else; c_flag = 1; end;
patch(x,y,z,c_flag);
end;%if (T.nl==nl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tesselation_plot(nl_max,nl,T);
verbose=0;
prefix = 37*ones(1,2*T.nl);
if (verbose); disp(sprintf(' %s nl_max %d nl %d T.nl %d',prefix,nl_max,nl,T.nl)); end;
if T.nl<nl; 
if ~isempty(T.child_1); if (verbose); disp(sprintf(' %s calling child_1',prefix)); end; tesselation_plot(nl_max,nl,T.child_1); end;
if ~isempty(T.child_2); if (verbose); disp(sprintf(' %s calling child_2',prefix)); end; tesselation_plot(nl_max,nl,T.child_2); end;
if ~isempty(T.child_3); if (verbose); disp(sprintf(' %s calling child_3',prefix)); end; tesselation_plot(nl_max,nl,T.child_3); end;
if ~isempty(T.child_4); if (verbose); disp(sprintf(' %s calling child_4',prefix)); end; tesselation_plot(nl_max,nl,T.child_4); end;
end;%if T.nl<nl_max; 
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

function T = tesselation_make(nl_max,nl,v1_i,v2_i,v3_i);
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
child_1 = tesselation_make(nl_max,nl+1,v1,m3,m2);
child_2 = tesselation_make(nl_max,nl+1,v2,m1,m3);
child_3 = tesselation_make(nl_max,nl+1,v3,m2,m1);
child_4 = tesselation_make(nl_max,nl+1,m1,m2,m3);
end;%if (nl<nl_max);
T = struct(...
	   'nl',nl,...
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
