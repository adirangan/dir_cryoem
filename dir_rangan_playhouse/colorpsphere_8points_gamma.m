function [c_] = colorpsphere_8points_gamma(polar_a_,azimu_b_,gamma);
% colormap for the surface of the sphere. ;
% gamma increases luminance of blue areas. ;
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
%p_c_(:,1+na,:) = repmat(colorpsphere_8points(polar_a_k_cur,azimu_b_k_cur),4,1);
%plot3(k_x,k_y,k_z,'.','Color',c_(nc,:));
polar_a_k_cur_(1+na) = polar_a_k_cur;
azimu_b_k_cur_(1+na) = azimu_b_k_cur;
na=na+1;
end;%for npolar_a=1:n_polar_a;
end;%for nazimu_b=1:n_azimu_b;
p_c_ = repmat(reshape(colorpsphere_8points(polar_a_k_cur_,azimu_b_k_cur_),[1,n_A,3]),4,1,1);
p = patch(p_x_,p_y_,p_z_,p_c_); set(p,'EdgeColor','none');
axis vis3d; view(10,25);
xlabel('x');
ylabel('y');
zlabel('z');
disp('returning'); return;
end;%if (nargin<2);

if (nargin<3); gamma=[]; end;
if isempty(gamma); gamma = 0.35; end;

polar_a_ = polar_a_(:);
azimu_b_ = azimu_b_(:);

n_p = numel(polar_a_);

p_x_ = sin(polar_a_).*cos(azimu_b_);
p_y_ = sin(polar_a_).*sin(azimu_b_);
p_z_ = cos(polar_a_);

t = +pi/4; w = +acos(1/sqrt(3));
Rz = [ cos(t) , -sin(t) , 0 ; +sin(t) , cos(t) , 0 ; 0 , 0 , 1 ];
Ry = [ cos(w) , 0 , +sin(w) ; 0 , 1 , 0 ; -sin(w) , 0 , cos(w) ];

tmp__ = transpose(Rz*Ry*transpose([p_x_,p_y_,p_z_]));
p_x_ = tmp__(:,1);
p_y_ = tmp__(:,2);
p_z_ = tmp__(:,3);

[p_s_,p_i_] = max([abs(p_x_),abs(p_y_),abs(p_z_)],[],2); p_i_ = p_i_ - 1;
c_x_ = p_x_ ./ p_s_;
c_y_ = p_y_ ./ p_s_;
c_z_ = p_z_ ./ p_s_;

index_paa_ = find( (p_i_ == 0) & abs(p_s_ - +p_x_)<1e-6 ) - 1;
index_naa_ = find( (p_i_ == 0) & abs(p_s_ - -p_x_)<1e-6 ) - 1;
index_apa_ = find( (p_i_ == 1) & abs(p_s_ - +p_y_)<1e-6 ) - 1;
index_ana_ = find( (p_i_ == 1) & abs(p_s_ - -p_y_)<1e-6 ) - 1;
index_aap_ = find( (p_i_ == 2) & abs(p_s_ - +p_z_)<1e-6 ) - 1;
index_aan_ = find( (p_i_ == 2) & abs(p_s_ - -p_z_)<1e-6 ) - 1;

c_ppp_ = [1,1,1];
c_ppn_ = [1,1,0];
c_pnp_ = [1,0,1];
c_npp_ = [0,1,1];
c_pnn_ = [1,0,0];
c_npn_ = [0,1,0];
c_nnp_ = [0,0,1];
c_nnn_ = [0,0,0];

if (gamma>0);
h = gamma^2;
g = gamma^1;
c_ppp_ = [1,1,1];
c_ppn_ = [1-h,1-g,0];
c_pnp_ = [1,g,1];
c_npp_ = [g,1,1];
c_pnn_ = [1,h,h];
c_npn_ = [h,1,h];
c_nnp_ = [g,g,1];
c_nnn_ = [0,0,0];
end;%if (gamma>0);

c_ = zeros(n_p,3);
c_(1+index_paa_,:) = bilinear_interpolant(c_y_(1+index_paa_),c_z_(1+index_paa_),c_pnn_,c_pnp_,c_ppn_,c_ppp_);
c_(1+index_naa_,:) = bilinear_interpolant(c_y_(1+index_naa_),c_z_(1+index_naa_),c_nnn_,c_nnp_,c_npn_,c_npp_);
c_(1+index_apa_,:) = bilinear_interpolant(c_x_(1+index_apa_),c_z_(1+index_apa_),c_npn_,c_npp_,c_ppn_,c_ppp_);
c_(1+index_ana_,:) = bilinear_interpolant(c_x_(1+index_ana_),c_z_(1+index_ana_),c_nnn_,c_nnp_,c_pnn_,c_pnp_);
c_(1+index_aap_,:) = bilinear_interpolant(c_x_(1+index_aap_),c_y_(1+index_aap_),c_nnp_,c_npp_,c_pnp_,c_ppp_);
c_(1+index_aan_,:) = bilinear_interpolant(c_x_(1+index_aan_),c_y_(1+index_aan_),c_nnn_,c_npn_,c_pnn_,c_ppn_);


  


