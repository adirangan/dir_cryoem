function T_ = tesselation_ver2();

% testing HandleObject ;
%{
tmp = [0,0,0,0];
hT_test = HandleObject(tmp);
disp(num2str(hT_test.O));
testhandle(hT_test);
disp(num2str(hT_test.O));
disp('returning'); return;
 %}

rng(0);
%n_points = 8192;
n_points = 128;
L_ = randn(3,n_points); for npoints=1:n_points; L_(:,npoints) = L_(:,npoints)/norm(L_(:,npoints)); end;

v_n00 = [-1; 0; 0];
v_0n0 = [ 0;-1; 0];
v_00n = [ 0; 0;-1];
v_p00 = [+1; 0; 0];
v_0p0 = [ 0;+1; 0];
v_00p = [ 0; 0;+1];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test single octant ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
nl_nnn = tesselation_find_depth(0,v_n00,v_0n0,v_00n,0,L_);
nl_max = nl_nnn;
disp(sprintf(' %% n_points %d --> nl_max %d',n_points,nl_max));

num_nnn = tesselation_find_number(0,v_n00,v_0n0,v_00n,0,L_);
num_sum = num_nnn;
disp(sprintf(' %% n_points %d --> num_sum %d',n_points,num_sum));

hT_nl_ = HandleObject(zeros(num_sum*1,1)); % level ;
hT_up_ = HandleObject(zeros(num_sum*1,1)); % parity (up vs down) ;
hT_id_ = HandleObject(zeros(num_sum*1,1)); % self pointer ;
hT_pa_ = HandleObject(zeros(num_sum*1,1)); % parent pointer ;
hT_v1_ = HandleObject(zeros(num_sum*3,1)); % vertex 1 ;
hT_v2_ = HandleObject(zeros(num_sum*3,1)); % vertex 2 ;
hT_v3_ = HandleObject(zeros(num_sum*3,1)); % vertex 3 ;
hT_vm_ = HandleObject(zeros(num_sum*3,1)); % vertex center ;
hT_m1_ = HandleObject(zeros(num_sum*3,1)); % edge midpoint 1 ;
hT_m2_ = HandleObject(zeros(num_sum*3,1)); % edge midpoint 2 ;
hT_m3_ = HandleObject(zeros(num_sum*3,1)); % edge midpoint 3 ;
hT_e1_ = HandleObject(zeros(num_sum*3,1)); % edge vector 1 ;
hT_e2_ = HandleObject(zeros(num_sum*3,1)); % edge vector 2 ;
hT_e3_ = HandleObject(zeros(num_sum*3,1)); % edge vector 3 ;
hT_n1_ = HandleObject(zeros(num_sum*3,1)); % edge normal 1 ;
hT_n2_ = HandleObject(zeros(num_sum*3,1)); % edge normal 2 ;
hT_n3_ = HandleObject(zeros(num_sum*3,1)); % edge normal 3 ;
hT_nn_ = HandleObject(zeros(num_sum*3,1)); % center normal ;
hT_ll_ = HandleObject(zeros(num_sum*1,1)); % number of points from L_ in hT_ ;
hT_c1_ = HandleObject(zeros(num_sum*1,1)); % child_1 pointer ;
hT_c2_ = HandleObject(zeros(num_sum*1,1)); % child_2 pointer ;
hT_c3_ = HandleObject(zeros(num_sum*1,1)); % child_3 pointer ;
hT_c4_ = HandleObject(zeros(num_sum*1,1)); % child_4 pointer ;
hT_ls_ = HandleObject(zeros(num_sum*1,1)); % cumsum of hT_ll_ --> starting location of point list for hT_ ;
hT_end = HandleObject(0);

index_output = 1; nl = 0;
index_parent = 0; v1_input = v_n00; v2_input = v_0n0; v3_input = v_00n; parity = 0;
index_output = tesselation_index(nl_max,nl,v1_input,v2_input,v3_input,parity,L_,index_parent,index_output,hT_nl_,hT_up_,hT_id_,hT_pa_,hT_v1_,hT_v2_,hT_v3_,hT_vm_,hT_m1_,hT_m2_,hT_m3_,hT_e1_,hT_e2_,hT_e3_,hT_n1_,hT_n2_,hT_n3_,hT_nn_,hT_ll_,hT_c1_,hT_c2_,hT_c3_,hT_c4_,hT_ls_,hT_end);

disp(sprintf(' %% hT_nl length %d',length(hT_nl_.O)));
disp('returning');return;
 %}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test all eight octants ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

num_nnn = tesselation_find_number(0,v_n00,v_0n0,v_00n,0,L_);
num_nnp = tesselation_find_number(0,v_n00,v_0n0,v_00p,1,L_);
num_npn = tesselation_find_number(0,v_n00,v_0p0,v_00n,1,L_);
num_npp = tesselation_find_number(0,v_n00,v_0p0,v_00p,0,L_);
num_pnn = tesselation_find_number(0,v_p00,v_0n0,v_00n,1,L_);
num_pnp = tesselation_find_number(0,v_p00,v_0n0,v_00p,0,L_);
num_ppn = tesselation_find_number(0,v_p00,v_0p0,v_00n,0,L_);
num_ppp = tesselation_find_number(0,v_p00,v_0p0,v_00p,1,L_);
num_= [num_nnn,num_nnp,num_npn,num_npp,num_pnn,num_pnp,num_ppn,num_ppp];
disp(num2str(num_));
num_sum = sum(num_);
disp(sprintf(' %% n_points %d --> num_sum %d',n_points,num_sum));

hT_nl_ = HandleObject(zeros(num_sum*1,1)); % level ;
hT_up_ = HandleObject(zeros(num_sum*1,1)); % parity (up vs down) ;
hT_id_ = HandleObject(zeros(num_sum*1,1)); % self pointer ;
hT_pa_ = HandleObject(zeros(num_sum*1,1)); % parent pointer ;
hT_v1_ = HandleObject(zeros(num_sum*3,1)); % vertex 1 ;
hT_v2_ = HandleObject(zeros(num_sum*3,1)); % vertex 2 ;
hT_v3_ = HandleObject(zeros(num_sum*3,1)); % vertex 3 ;
hT_vm_ = HandleObject(zeros(num_sum*3,1)); % vertex center ;
hT_m1_ = HandleObject(zeros(num_sum*3,1)); % edge midpoint 1 ;
hT_m2_ = HandleObject(zeros(num_sum*3,1)); % edge midpoint 2 ;
hT_m3_ = HandleObject(zeros(num_sum*3,1)); % edge midpoint 3 ;
hT_e1_ = HandleObject(zeros(num_sum*3,1)); % edge vector 1 ;
hT_e2_ = HandleObject(zeros(num_sum*3,1)); % edge vector 2 ;
hT_e3_ = HandleObject(zeros(num_sum*3,1)); % edge vector 3 ;
hT_n1_ = HandleObject(zeros(num_sum*3,1)); % edge normal 1 ;
hT_n2_ = HandleObject(zeros(num_sum*3,1)); % edge normal 2 ;
hT_n3_ = HandleObject(zeros(num_sum*3,1)); % edge normal 3 ;
hT_nn_ = HandleObject(zeros(num_sum*3,1)); % center normal ;
hT_ll_ = HandleObject(zeros(num_sum*1,1)); % number of points from L_ in hT_ ;
hT_c1_ = HandleObject(zeros(num_sum*1,1)); % child_1 pointer ;
hT_c2_ = HandleObject(zeros(num_sum*1,1)); % child_2 pointer ;
hT_c3_ = HandleObject(zeros(num_sum*1,1)); % child_3 pointer ;
hT_c4_ = HandleObject(zeros(num_sum*1,1)); % child_4 pointer ;
hT_ls_ = HandleObject(zeros(num_sum*1,1)); % cumsum of hT_ll_ --> starting location of point list for hT_ ;
hT_end = HandleObject(0);

index_output = 1; nl = 0;
index_parent = 0; v1_input = v_n00; v2_input = v_0n0; v3_input = v_00n; parity = 0;
index_output = tesselation_index(nl_max,nl,v1_input,v2_input,v3_input,parity,L_,index_parent,index_output,hT_nl_,hT_up_,hT_id_,hT_pa_,hT_v1_,hT_v2_,hT_v3_,hT_vm_,hT_m1_,hT_m2_,hT_m3_,hT_e1_,hT_e2_,hT_e3_,hT_n1_,hT_n2_,hT_n3_,hT_nn_,hT_ll_,hT_c1_,hT_c2_,hT_c3_,hT_c4_,hT_ls_,hT_end);
index_parent = 0; v1_input = v_n00; v2_input = v_0n0; v3_input = v_00p; parity = 1;
index_output = tesselation_index(nl_max,nl,v1_input,v2_input,v3_input,parity,L_,index_parent,index_output,hT_nl_,hT_up_,hT_id_,hT_pa_,hT_v1_,hT_v2_,hT_v3_,hT_vm_,hT_m1_,hT_m2_,hT_m3_,hT_e1_,hT_e2_,hT_e3_,hT_n1_,hT_n2_,hT_n3_,hT_nn_,hT_ll_,hT_c1_,hT_c2_,hT_c3_,hT_c4_,hT_ls_,hT_end);
index_parent = 0; v1_input = v_n00; v2_input = v_0p0; v3_input = v_00n; parity = 1;
index_output = tesselation_index(nl_max,nl,v1_input,v2_input,v3_input,parity,L_,index_parent,index_output,hT_nl_,hT_up_,hT_id_,hT_pa_,hT_v1_,hT_v2_,hT_v3_,hT_vm_,hT_m1_,hT_m2_,hT_m3_,hT_e1_,hT_e2_,hT_e3_,hT_n1_,hT_n2_,hT_n3_,hT_nn_,hT_ll_,hT_c1_,hT_c2_,hT_c3_,hT_c4_,hT_ls_,hT_end);
index_parent = 0; v1_input = v_n00; v2_input = v_0p0; v3_input = v_00p; parity = 0;
index_output = tesselation_index(nl_max,nl,v1_input,v2_input,v3_input,parity,L_,index_parent,index_output,hT_nl_,hT_up_,hT_id_,hT_pa_,hT_v1_,hT_v2_,hT_v3_,hT_vm_,hT_m1_,hT_m2_,hT_m3_,hT_e1_,hT_e2_,hT_e3_,hT_n1_,hT_n2_,hT_n3_,hT_nn_,hT_ll_,hT_c1_,hT_c2_,hT_c3_,hT_c4_,hT_ls_,hT_end);
index_parent = 0; v1_input = v_p00; v2_input = v_0n0; v3_input = v_00n; parity = 1;
index_output = tesselation_index(nl_max,nl,v1_input,v2_input,v3_input,parity,L_,index_parent,index_output,hT_nl_,hT_up_,hT_id_,hT_pa_,hT_v1_,hT_v2_,hT_v3_,hT_vm_,hT_m1_,hT_m2_,hT_m3_,hT_e1_,hT_e2_,hT_e3_,hT_n1_,hT_n2_,hT_n3_,hT_nn_,hT_ll_,hT_c1_,hT_c2_,hT_c3_,hT_c4_,hT_ls_,hT_end);
index_parent = 0; v1_input = v_p00; v2_input = v_0n0; v3_input = v_00p; parity = 0;
index_output = tesselation_index(nl_max,nl,v1_input,v2_input,v3_input,parity,L_,index_parent,index_output,hT_nl_,hT_up_,hT_id_,hT_pa_,hT_v1_,hT_v2_,hT_v3_,hT_vm_,hT_m1_,hT_m2_,hT_m3_,hT_e1_,hT_e2_,hT_e3_,hT_n1_,hT_n2_,hT_n3_,hT_nn_,hT_ll_,hT_c1_,hT_c2_,hT_c3_,hT_c4_,hT_ls_,hT_end);
index_parent = 0; v1_input = v_p00; v2_input = v_0p0; v3_input = v_00n; parity = 0;
index_output = tesselation_index(nl_max,nl,v1_input,v2_input,v3_input,parity,L_,index_parent,index_output,hT_nl_,hT_up_,hT_id_,hT_pa_,hT_v1_,hT_v2_,hT_v3_,hT_vm_,hT_m1_,hT_m2_,hT_m3_,hT_e1_,hT_e2_,hT_e3_,hT_n1_,hT_n2_,hT_n3_,hT_nn_,hT_ll_,hT_c1_,hT_c2_,hT_c3_,hT_c4_,hT_ls_,hT_end);
index_parent = 0; v1_input = v_p00; v2_input = v_0p0; v3_input = v_00p; parity = 1;
index_output = tesselation_index(nl_max,nl,v1_input,v2_input,v3_input,parity,L_,index_parent,index_output,hT_nl_,hT_up_,hT_id_,hT_pa_,hT_v1_,hT_v2_,hT_v3_,hT_vm_,hT_m1_,hT_m2_,hT_m3_,hT_e1_,hT_e2_,hT_e3_,hT_n1_,hT_n2_,hT_n3_,hT_nn_,hT_ll_,hT_c1_,hT_c2_,hT_c3_,hT_c4_,hT_ls_,hT_end);

disp(sprintf(' %% hT_nl length %d',length(hT_nl_.O)));

% need to generate tesselation ;
% need to generate list linking L_ --> T_ for each level ;
% need to generate list linking T_ --> L_ for each level ;

disp('returning'); return;

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
function output = tesselation_find_number(nl,v1_i,v2_i,v3_i,parity,L_);
output = 0;
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
output = 0;
elseif (n_L_sub_length==1);
output = 1;
elseif (n_L_sub_length>1);
nl_child_1 = tesselation_find_number(nl+1,v1,m3,m2, parity,L_sub_);
nl_child_2 = tesselation_find_number(nl+1,v2,m1,m3, parity,L_sub_);
nl_child_3 = tesselation_find_number(nl+1,v3,m2,m1, parity,L_sub_);
nl_child_4 = tesselation_find_number(nl+1,m1,m2,m3,~parity,L_sub_);
nl_child_ = [nl_child_1,nl_child_2,nl_child_3,nl_child_4];
output = 1 + sum([nl_child_1,nl_child_2,nl_child_3,nl_child_4]);
end;% if (n_L_sub_length>1);

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
nl_child_ = [nl_child_1,nl_child_2,nl_child_3,nl_child_4];
output = max(nl_child_);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function testhandle(hT);
hT.O(1)=1;
hT.O(2)=2;
hT.O(3)=3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function index_output = ...
  tesselation_index(nl_max,nl,...
		    v1_input,v2_input,v3_input,...
		    parity,...
		    L_,...
		    index_parent,...
		    index__input,...
		    hT_nl_,...
		    hT_up_,...
		    hT_id_,...
		    hT_pa_,...
		    hT_v1_,...
		    hT_v2_,...
		    hT_v3_,...
		    hT_vm_,...
		    hT_m1_,...
		    hT_m2_,...
		    hT_m3_,...
		    hT_e1_,...
		    hT_e2_,...
		    hT_e3_,...
		    hT_n1_,...
		    hT_n2_,...
		    hT_n3_,...
		    hT_nn_,...
		    hT_ll_,...
		    hT_c1_,...
		    hT_c2_,...
		    hT_c3_,...
		    hT_c4_,...
		    hT_ls_,...
		    hT_end);
if (nl>nl_max); disp(sprintf(' %% Warning! nl %d/%d in tesselation_index',nl,nl_max)); end;%if (nl>nl_max);
verbose=0;
if (verbose); prefix = 37*ones(1,1*(nl+1)); end;
v1_input = v1_input/norm(v1_input); 
v2_input = v2_input/norm(v2_input); 
v3_input = v3_input/norm(v3_input);
v1 = v1_input; v2 = v2_input; v3 = v3_input; 
vm = (v1+v2+v3)/3;
nn = cross(v3-v2,v1-v3);
aa = dot(vm,nn);
if (aa>=0); v1 = v1_input; v2 = v2_input; v3 = v3_input; end;
if (aa< 0); v1 = v2_input; v2 = v1_input; v3 = v3_input; end;
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
if (verbose); disp(sprintf('%s nl %0.2d parity %d index_parent %0.3d index__input %0.3d; n_L_length %0.4d; n_L_sub_length %0.4d',prefix,nl,parity,index_parent,index__input,n_L_length,n_L_sub_length)); end;
if (n_L_sub_length==0);
%if (verbose); disp(sprintf('%s skipping self',prefix)); end;
index_output = index__input;
end;%if (n_L_sub_length==0);
if (n_L_sub_length>0);
%if (verbose); disp(sprintf('%s creating self',prefix)); end;
hT_nl_.O(index__input) = nl;
hT_up_.O(index__input) = parity;
hT_id_.O(index__input) = index__input;
hT_pa_.O(index__input) = index_parent;
ij_tmp = (1:3) + 3*(index__input-1);
hT_v1_.O(ij_tmp) = v1;
hT_v2_.O(ij_tmp) = v2;
hT_v3_.O(ij_tmp) = v3;
hT_vm_.O(ij_tmp) = vm;
hT_m1_.O(ij_tmp) = m1;
hT_m2_.O(ij_tmp) = m2;
hT_m3_.O(ij_tmp) = m3;
hT_e1_.O(ij_tmp) = e1;
hT_e2_.O(ij_tmp) = e2;
hT_e3_.O(ij_tmp) = e3;
hT_n1_.O(ij_tmp) = n1;
hT_n2_.O(ij_tmp) = n2;
hT_n3_.O(ij_tmp) = n3;
hT_nn_.O(ij_tmp) = nn;
hT_ll_.O(index__input) = n_L_sub_length;
hT_c1_.O(index__input) = 0;
hT_c2_.O(index__input) = 0;
hT_c3_.O(index__input) = 0;
hT_c4_.O(index__input) = 0;
index_output = index__input+1;
end;% if (n_L_sub_length>0);
if (n_L_sub_length>1);
%if (verbose); disp(sprintf('%s creating children',prefix)); end;
index_tempor = index_output;
if (verbose); disp(sprintf('%s creating child_1 at index %d',prefix,index_output)); end;
index_output = ...
  tesselation_index(nl_max,nl+1,...
		    v1,m3,m2,...
		     parity,...
		    L_sub_,...
		    hT_id_.O(index__input),...
		    index_output,...
		    hT_nl_,...
		    hT_up_,...
		    hT_id_,...
		    hT_pa_,...
		    hT_v1_,...
		    hT_v2_,...
		    hT_v3_,...
		    hT_vm_,...
		    hT_m1_,...
		    hT_m2_,...
		    hT_m3_,...
		    hT_e1_,...
		    hT_e2_,...
		    hT_e3_,...
		    hT_n1_,...
		    hT_n2_,...
		    hT_n3_,...
		    hT_nn_,...
		    hT_ll_,...
		    hT_c1_,...
		    hT_c2_,...
		    hT_c3_,...
		    hT_c4_,...
		    hT_ls_,...
		    hT_end);
if (index_output>index_tempor); hT_c1_.O(index__input) = index_tempor; else hT_c1_.O(index__input) = 0; end;
index_tempor = index_output;
if (verbose); disp(sprintf('%s creating child_2 at index %d',prefix,index_output)); end;
index_output = ...
  tesselation_index(nl_max,nl+1,...
		    v2,m1,m3,...
		     parity,...
		    L_sub_,...
		    hT_id_.O(index__input),...
		    index_output,...
		    hT_nl_,...
		    hT_up_,...
		    hT_id_,...
		    hT_pa_,...
		    hT_v1_,...
		    hT_v2_,...
		    hT_v3_,...
		    hT_vm_,...
		    hT_m1_,...
		    hT_m2_,...
		    hT_m3_,...
		    hT_e1_,...
		    hT_e2_,...
		    hT_e3_,...
		    hT_n1_,...
		    hT_n2_,...
		    hT_n3_,...
		    hT_nn_,...
		    hT_ll_,...
		    hT_c1_,...
		    hT_c2_,...
		    hT_c3_,...
		    hT_c4_,...
		    hT_ls_,...
		    hT_end);
if (index_output>index_tempor); hT_c2_.O(index__input) = index_tempor; else hT_c2_.O(index__input) = 0; end;
index_tempor = index_output;
if (verbose); disp(sprintf('%s creating child_3 at index %d',prefix,index_output)); end;
index_output = ...
  tesselation_index(nl_max,nl+1,...
		    v3,m2,m1,...
		     parity,...
		    L_sub_,...
		    hT_id_.O(index__input),...
		    index_output,...
		    hT_nl_,...
		    hT_up_,...
		    hT_id_,...
		    hT_pa_,...
		    hT_v1_,...
		    hT_v2_,...
		    hT_v3_,...
		    hT_vm_,...
		    hT_m1_,...
		    hT_m2_,...
		    hT_m3_,...
		    hT_e1_,...
		    hT_e2_,...
		    hT_e3_,...
		    hT_n1_,...
		    hT_n2_,...
		    hT_n3_,...
		    hT_nn_,...
		    hT_ll_,...
		    hT_c1_,...
		    hT_c2_,...
		    hT_c3_,...
		    hT_c4_,...
		    hT_ls_,...
		    hT_end);
if (index_output>index_tempor); hT_c3_.O(index__input) = index_tempor; else hT_c3_.O(index__input) = 0; end;
index_tempor = index_output;
if (verbose); disp(sprintf('%s creating child_4 at index %d',prefix,index_output)); end;
index_output = ...
  tesselation_index(nl_max,nl+1,...
		    m1,m2,m3,...
		    ~parity,...
		    L_sub_,...
		    hT_id_.O(index__input),...
		    index_output,...
		    hT_nl_,...
		    hT_up_,...
		    hT_id_,...
		    hT_pa_,...
		    hT_v1_,...
		    hT_v2_,...
		    hT_v3_,...
		    hT_vm_,...
		    hT_m1_,...
		    hT_m2_,...
		    hT_m3_,...
		    hT_e1_,...
		    hT_e2_,...
		    hT_e3_,...
		    hT_n1_,...
		    hT_n2_,...
		    hT_n3_,...
		    hT_nn_,...
		    hT_ll_,...
		    hT_c1_,...
		    hT_c2_,...
		    hT_c3_,...
		    hT_c4_,...
		    hT_ls_,...
		    hT_end);
if (index_output>index_tempor); hT_c4_.O(index__input) = index_tempor; else hT_c4_.O(index__input) = 0; end;
end;%if (n_L_sub_length>1);
if (verbose); disp(sprintf('%s final index_output %d',prefix,index_output)); end;

