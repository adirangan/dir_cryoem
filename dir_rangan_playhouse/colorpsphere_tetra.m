function [c_] = colorpsphere_tetra(polar_a,azimu_b);
% colormap for the surface of the (projective) sphere. ;
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
%p_c_(:,1+na,:) = repmat(colorpsphere_tetra(polar_a_k_cur,azimu_b_k_cur),4,1);
%plot3(k_x,k_y,k_z,'.','Color',c_(nc,:));
polar_a_k_cur_(1+na) = polar_a_k_cur;
azimu_b_k_cur_(1+na) = azimu_b_k_cur;
na=na+1;
end;%for npolar_a=1:n_polar_a;
end;%for nazimu_b=1:n_azimu_b;
p_c_ = repmat(reshape(colorpsphere_tetra(polar_a_k_cur_,azimu_b_k_cur_),[1,n_A,3]),4,1,1);
p = patch(p_x_,p_y_,p_z_,p_c_); set(p,'EdgeColor','none');
axis vis3d; view(10,25);
disp('returning'); return;
end;%if (nargin<2);

polar_a = polar_a(:);
azimu_b = azimu_b(:);

w = 2*pi/3; c = cos(w); s = sin(w);
v1a_ = [1, 0,0] - [0,0,sqrt(2)/4];
v2a_ = [c,+s,0] - [0,0,sqrt(2)/4];
v3a_ = [c,-s,0] - [0,0,sqrt(2)/4];
v4a_ = [0,0,sqrt(2)] - [0,0,sqrt(2)/4];
v1b_ = -v1a_;
v2b_ = -v2a_;
v3b_ = -v3a_;
v4b_ = -v4a_;

c_x_ = sin(polar_a).*cos(azimu_b);
c_y_ = sin(polar_a).*sin(azimu_b);
c_z_ = cos(polar_a);

d1a_ = sqrt( (c_x_ - v1a_(1)).^2 + (c_y_ - v1a_(2)).^2 + (c_z_ - v1a_(3)).^2 );
d2a_ = sqrt( (c_x_ - v2a_(1)).^2 + (c_y_ - v2a_(2)).^2 + (c_z_ - v2a_(3)).^2 );
d3a_ = sqrt( (c_x_ - v3a_(1)).^2 + (c_y_ - v3a_(2)).^2 + (c_z_ - v3a_(3)).^2 );
d4a_ = sqrt( (c_x_ - v4a_(1)).^2 + (c_y_ - v4a_(2)).^2 + (c_z_ - v4a_(3)).^2 );
d1b_ = sqrt( (c_x_ - v1b_(1)).^2 + (c_y_ - v1b_(2)).^2 + (c_z_ - v1b_(3)).^2 );
d2b_ = sqrt( (c_x_ - v2b_(1)).^2 + (c_y_ - v2b_(2)).^2 + (c_z_ - v2b_(3)).^2 );
d3b_ = sqrt( (c_x_ - v3b_(1)).^2 + (c_y_ - v3b_(2)).^2 + (c_z_ - v3b_(3)).^2 );
d4b_ = sqrt( (c_x_ - v4b_(1)).^2 + (c_y_ - v4b_(2)).^2 + (c_z_ - v4b_(3)).^2 );

c1a_ = [1.0 , 0.1 , 0.0];
c2a_ = [0.0 , 1.0 , 0.5];
c3a_ = [0.5 , 0.0 , 1.0];
c4a_ = [1.0 , 1.0 , 0.0];
c1b_ = c1a_;
c2b_ = c2a_;
c3b_ = c3a_;
c4b_ = c4a_;

delta = 0.35;
e1a_ = exp(-d1a_./delta);
e2a_ = exp(-d2a_./delta);
e3a_ = exp(-d3a_./delta);
e4a_ = exp(-d4a_./delta);
e1b_ = exp(-d1b_./delta);
e2b_ = exp(-d2b_./delta);
e3b_ = exp(-d3b_./delta);
e4b_ = exp(-d4b_./delta);
s_ = sum([e1a_,e2a_,e3a_,e4a_,e1b_,e2b_,e3b_,e4b_],2);

f1a_ = e1a_./s_;
f2a_ = e2a_./s_;
f3a_ = e3a_./s_;
f4a_ = e4a_./s_;
f1b_ = e1b_./s_;
f2b_ = e2b_./s_;
f3b_ = e3b_./s_;
f4b_ = e4b_./s_;

c_ = f1a_*c1a_ + f2a_*c2a_ + f3a_*c3a_ + f4a_*c4a_ + f1b_*c1b_ + f2b_*c2b_ + f3b_*c3b_ + f4b_*c4b_;




