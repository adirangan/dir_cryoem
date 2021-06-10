function [c_] = colorsphere_trio(polar_a,azimu_b,c_hsv_);
% colormap for the surface of the sphere. ;
if (nargin<2);
c_hsv_ = colormap('hsv');
k_r = 1;
n_azimu_b = 128;
n_polar_a = n_azimu_b/2;
azimu_b_k_ = linspace(0,2*pi,n_azimu_b+1); azimu_b_k_ = azimu_b_k_(1:end-1);
polar_a_k_ = linspace(0,pi,n_polar_a);
p_y_ = zeros(4,n_azimu_b*n_polar_a);
p_z_ = zeros(4,n_azimu_b*n_polar_a);
p_c_ = zeros(4,n_azimu_b*n_polar_a,3);
dazimu_b = mean(diff(azimu_b_k_));
na=0;
for nazimu_b=1:n_azimu_b;
azimu_b_k_pre = azimu_b_k_(nazimu_b) - dazimu_b/2;
azimu_b_k_cur = azimu_b_k_(nazimu_b);
azimu_b_k_pos = azimu_b_k_(nazimu_b) + dazimu_b/2;
for npolar_a=1:n_polar_a;
polar_a_k_pre = 0.5*(polar_a_k_(npolar_a) + polar_a_k_(max(1,npolar_a-1)));
polar_a_k_cur = polar_a_k_(npolar_a);
polar_a_k_pos = 0.5*(polar_a_k_(npolar_a) + polar_a_k_(min(n_polar_a,npolar_a+1)));
k_x_cur_cur = k_r*cos(azimu_b_k_cur)*sin(polar_a_k_cur); k_y_cur_cur = k_r*sin(azimu_b_k_cur)*sin(polar_a_k_cur); k_z_cur_cur = k_r*cos(polar_a_k_cur);
k_x_pre_pre = k_r*cos(azimu_b_k_pre)*sin(polar_a_k_pre); k_y_pre_pre = k_r*sin(azimu_b_k_pre)*sin(polar_a_k_pre); k_z_pre_pre = k_r*cos(polar_a_k_pre);
k_x_pre_pos = k_r*cos(azimu_b_k_pre)*sin(polar_a_k_pos); k_y_pre_pos = k_r*sin(azimu_b_k_pre)*sin(polar_a_k_pos); k_z_pre_pos = k_r*cos(polar_a_k_pos);
k_x_pos_pre = k_r*cos(azimu_b_k_pos)*sin(polar_a_k_pre); k_y_pos_pre = k_r*sin(azimu_b_k_pos)*sin(polar_a_k_pre); k_z_pos_pre = k_r*cos(polar_a_k_pre);
k_x_pos_pos = k_r*cos(azimu_b_k_pos)*sin(polar_a_k_pos); k_y_pos_pos = k_r*sin(azimu_b_k_pos)*sin(polar_a_k_pos); k_z_pos_pos = k_r*cos(polar_a_k_pos);
p_x_(:,1+na) = [k_x_pre_pre;k_x_pos_pre;k_x_pos_pos;k_x_pre_pos];
p_y_(:,1+na) = [k_y_pre_pre;k_y_pos_pre;k_y_pos_pos;k_y_pre_pos];
p_z_(:,1+na) = [k_z_pre_pre;k_z_pos_pre;k_z_pos_pos;k_z_pre_pos];
p_c_(:,1+na,:) = repmat(colorsphere_trio(polar_a_k_cur,azimu_b_k_cur,c_hsv_),4,1);
%plot3(k_x,k_y,k_z,'.','Color',c_(nc,:));
na=na+1;
end;%for npolar_a=1:n_polar_a;
end;%for nazimu_b=1:n_azimu_b;
p = patch(p_x_,p_y_,p_z_,p_c_); set(p,'EdgeColor','none');
axis vis3d; view(10,25);
disp('returning'); return;
end;%if (nargin<2);

%{
c_polar_a_hemis = (0.5*(1 + cos(2*polar_a))).^0.75;
c_polar_a_equat = (1 - c_polar_a_hemis).^4;
gamma_equat = 0.15;
r_azimu_b_equat = abs(periodize(azimu_b-0.5*pi/3,-pi,+pi)); r_azimu_b_equat = 0.5 + 0.5*r_azimu_b_equat; r_azimu_b_equat = (max(0,1 - r_azimu_b_equat/(2*pi/3))).^gamma_equat;
g_azimu_b_equat = abs(periodize(azimu_b-2.5*pi/3,-pi,+pi)); g_azimu_b_equat = 0.5 + 0.5*g_azimu_b_equat; g_azimu_b_equat = (max(0,1 - g_azimu_b_equat/(2*pi/3))).^gamma_equat;
b_azimu_b_equat = abs(periodize(azimu_b-4.5*pi/3,-pi,+pi)); b_azimu_b_equat = 0.5 + 0.5*b_azimu_b_equat; b_azimu_b_equat = (max(0,1 - b_azimu_b_equat/(2*pi/3))).^gamma_equat;
gamma_north = 0.85;
r_azimu_b_north = abs(periodize(azimu_b-0*pi/3,-pi,+pi)); r_azimu_b_north = 0 + r_azimu_b_north; r_azimu_b_north = (max(0,1 - r_azimu_b_north/(2*pi/3))).^gamma_north;
g_azimu_b_north = abs(periodize(azimu_b-2*pi/3,-pi,+pi)); g_azimu_b_north = 0 + g_azimu_b_north; g_azimu_b_north = (max(0,1 - g_azimu_b_north/(2*pi/3))).^gamma_north;
b_azimu_b_north = abs(periodize(azimu_b-4*pi/3,-pi,+pi)); b_azimu_b_north = 0 + b_azimu_b_north; b_azimu_b_north = (max(0,1 - b_azimu_b_north/(2*pi/3))).^gamma_north;
gamma_south = 0.85;
r_azimu_b_south = abs(periodize(azimu_b-1*pi/3,-pi,+pi)); r_azimu_b_south = 1 - r_azimu_b_south; r_azimu_b_south = (max(0,1 - r_azimu_b_south/(2*pi/3))).^gamma_south;
g_azimu_b_south = abs(periodize(azimu_b-3*pi/3,-pi,+pi)); g_azimu_b_south = 1 - g_azimu_b_south; g_azimu_b_south = (max(0,1 - g_azimu_b_south/(2*pi/3))).^gamma_south;
b_azimu_b_south = abs(periodize(azimu_b-5*pi/3,-pi,+pi)); b_azimu_b_south = 1 - b_azimu_b_south; b_azimu_b_south = (max(0,1 - b_azimu_b_south/(2*pi/3))).^gamma_south;
if polar_a < pi/2 ; %<-- north;
c_ = c_polar_a_hemis*[r_azimu_b_north,g_azimu_b_north,b_azimu_b_north] + c_polar_a_equat*[r_azimu_b_equat,g_azimu_b_equat,b_azimu_b_equat];
end;%if polar_a < pi/2 ; %<-- north;
if polar_a >= pi/2 ; %<-- south;
c_ = c_polar_a_hemis*[r_azimu_b_south,g_azimu_b_south,b_azimu_b_south] + c_polar_a_equat*[r_azimu_b_equat,g_azimu_b_equat,b_azimu_b_equat];
end;%if polar_a >= pi/2 ; %<-- south;
 %}

n_c = size(c_hsv_,1);
c_polar_a_hemis = (0.5*(1 + cos(2*polar_a))).^0.75;
c_polar_a_equat = (1 - c_polar_a_hemis).^3;
c_azimu_b_equat = [1,1,1];
nc = max(1,min(n_c,round(n_c*azimu_b/(2*pi))));
c_azimu_b_north = c_hsv_(nc,:);
c_azimu_b_south = 1-c_hsv_(nc,:);
c_ = zeros(length(polar_a),3);
tmp_ij_ = find(polar_a< pi/2) ; %<-- north;
if (length(tmp_ij_)>0); c_(tmp_ij_,:) = repmat(transpose(c_polar_a_hemis(tmp_ij_)),1,3).*c_azimu_b_north(tmp_ij_,:) + repmat(transpose(c_polar_a_equat(tmp_ij_)),1,3); end;
tmp_ij_ = find(polar_a>=pi/2) ; %<-- south;
if (length(tmp_ij_)>0); c_(tmp_ij_,:) = repmat(transpose(c_polar_a_hemis(tmp_ij_)),1,3).*c_azimu_b_south(tmp_ij_,:) + repmat(transpose(c_polar_a_equat(tmp_ij_)),1,3); end;



