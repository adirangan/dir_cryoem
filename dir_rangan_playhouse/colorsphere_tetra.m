function [c_] = colorsphere_tetra(polar_a,azimu_b);
% colormap for the surface of the sphere. ;
if (nargin<2);
k_r = 1;
n_azimu_b = 128;
n_polar_a = n_azimu_b/2;
azimu_b_k_ = linspace(0,2*pi,n_azimu_b+1); azimu_b_k_ = azimu_b_k_(1:end-1);
polar_a_k_ = linspace(0,pi,n_polar_a);
n_A = n_azimu_b*n_polar_a;
p_x_ = zeros(4,n_A);
p_y_ = zeros(4,n_A);
p_z_ = zeros(4,n_A);
p_c_ = zeros(4,n_A,3);
polar_a_k_cur_ = zeros(n_A,1);
azimu_b_k_cur_ = zeros(n_A,1);
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
%p_c_(:,1+na,:) = repmat(colorsphere_tetra(polar_a_k_cur,azimu_b_k_cur),4,1);
%plot3(k_x,k_y,k_z,'.','Color',c_(nc,:));
polar_a_k_cur_(1+na) = polar_a_k_cur;
azimu_b_k_cur_(1+na) = azimu_b_k_cur;
na=na+1;
end;%for npolar_a=1:n_polar_a;
end;%for nazimu_b=1:n_azimu_b;
p_c_ = repmat(reshape(colorsphere_tetra(polar_a_k_cur_,azimu_b_k_cur_),[1,n_A,3]),4,1,1);
p = patch(p_x_,p_y_,p_z_,p_c_); set(p,'EdgeColor','none');
axis vis3d; view(10,25);
disp('returning'); return;
end;%if (nargin<2);

polar_a = polar_a(:);
azimu_b = azimu_b(:);

w = 2*pi/3; c = cos(w); s = sin(w);
v1_ = [1, 0,0] - [0,0,sqrt(2)/4];
v2_ = [c,+s,0] - [0,0,sqrt(2)/4];
v3_ = [c,-s,0] - [0,0,sqrt(2)/4];
v4_ = [0,0,sqrt(2)] - [0,0,sqrt(2)/4];

c_x_ = sin(polar_a).*cos(azimu_b);
c_y_ = sin(polar_a).*sin(azimu_b);
c_z_ = cos(polar_a);

d1_ = sqrt( (c_x_ - v1_(1)).^2 + (c_y_ - v1_(2)).^2 + (c_z_ - v1_(3)).^2 );
d2_ = sqrt( (c_x_ - v2_(1)).^2 + (c_y_ - v2_(2)).^2 + (c_z_ - v2_(3)).^2 );
d3_ = sqrt( (c_x_ - v3_(1)).^2 + (c_y_ - v3_(2)).^2 + (c_z_ - v3_(3)).^2 );
d4_ = sqrt( (c_x_ - v4_(1)).^2 + (c_y_ - v4_(2)).^2 + (c_z_ - v4_(3)).^2 );

c1_ = [1.0 , 0.5 , 0.0];
c2_ = [0.0 , 1.0 , 0.8];
c3_ = [0.8 , 0.0 , 1.0];
c4_ = [1.0 , 1.0 , 0.0];

delta = 0.35;
e1_ = exp(-d1_./delta);
e2_ = exp(-d2_./delta);
e3_ = exp(-d3_./delta);
e4_ = exp(-d4_./delta);
s_ = sum([e1_,e2_,e3_,e4_],2);

f1_ = e1_./s_;
f2_ = e2_./s_;
f3_ = e3_./s_;
f4_ = e4_./s_;

c_ = f1_*c1_ + f2_*c2_ + f3_*c3_ + f4_*c4_;




