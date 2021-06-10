function [shift_z__,L__] = get_template_inplane_shift_0(n_viewing_all,viewing_polar_a_all_,viewing_azimu_b_all_);
% This function returns the inplane-gamma-shift between templates defined via polar_a_ and azimu_b_. ;
% Based on the way this calculation is performed (see below), ;
% we expect that the location of the points in template X will align with those in template Y ;
% when the points in X are shifted from (gamma_z) to (gamma_z-shift_z). ;
% This means that the points of X themselves (in R3) should be rotated (about the viewing-angle of X) by -shift_z. ;
% Alternatively, the template X (when viewed as a function of its argument) can be thought of as rotated by +shift_z. ;
%%%%%%%%;
% The general formula used here is as follows. ;
% let sa and ca be sin(polar_a) and cos(polar_a), respectively. ;
% let sb and cb be sin(azimu_b) and cos(azimu_b), respectively. ;
% let sc and cc be sin(gamma_z) and cos(gamma_z), respectively. ;
% And rotation by azimu_b about the +z-axis is represented as: ;
% Rz(azimu_b) = ;
% [ +cb -sb 0 ] ;
% [ +sb +cb 0 ] ;
% [  0   0  1 ] ;
% And rotation by polar_a about the +y-axis is represented as: ;
% Ry(polar_a) = ;
% [ +ca 0 +sa ] ;
% [  0  1  0  ] ;
% [ -sa 0 +ca ] ;
% And rotation by gamma_z about the +z-axis is represented as: ;
% Rz(gamma_z) = ;
% [ +cc -sc 0 ] ;
% [ +sc +cc 0 ] ;
% [  0   0  1 ] ;
% Which, collectively, implies that under the transform: ;
% Rz(azimu_b) * Ry(polar_a) * Rz(gamma_z), ;
% Which is the same as: ;
% [ +cb -sb 0 ] [ +ca*cc -ca*sc +sa ]   [ +cb*ca*cc - sb*sc , -cb*ca*sc -sb*cc , +cb*sa ];
% [ +sb +cb 0 ] [ +sc    +cc    0   ] = [ +sb*ca*cc + cb*sc , -sb*ca*sc +cb*cc , +sb*sa ];
% [  0   0  1 ] [ -sa*cc +sa*sc +ca ]   [ -sa*cc            , +sa*sc           , +ca    ];
% the point [1;0;0] is mapped to: ;
% [ template_k_c_0 ; template_k_c_1 ; template_k_c_2 ] = [ +cb*ca*cc - sb*sc ; +sb*ca*cc + cb*sc ; -sa*cc ];
%%%%%%%%;
% Consequently, the vector of values (on the sphere) ;
% for template_X is given by: ;
% [ template_k_c_X_0 ; template_k_c_X_1 ; template_k_c_X_2 ] = [ +cbX*caX*ccX - sbX*scX ; +sbX*caX*ccX + cbX*scX ; -saX*ccX ];
% and template_Y is given by ;
% [ template_k_c_Y_0 ; template_k_c_Y_1 ; template_k_c_Y_2 ] = [ +cbY*caY*ccY - sbY*scY ; +sbY*caY*ccY + cbY*scY ; -saY*ccY ];
% Now if we set cX = gamma_z - shift_z and set cY = gamma_z, ;
% then our goal is to find shift_z such that: ;
% L = \int_{gamma_z=0}^{gamma_z=2*pi} dgamma_z \| template_k_c_X_ - template_k_c_Y_ \|^2 ;
% is as small as possible. ;
%%%%%%%%;
% L can be reduced to: ;
% L = 2*pi*(2 + a10/2*cos(shift_z) + a01/2*sin(shift_z)), ;
% and the formula for the optimal shift_z is: ;
% shift_z = arctan(a01/a10), where;
% -a01/2 = sin(bX-bY)(caX + caY) ;
% -a10/2 = cos(bX-bY)(1 + caX*caY) + saX*saY ;
%%%%%%%%;
% we chose the shift_z to correspond to the lower of the two L-values. ;
% The output: shift_z__ is a n_viewing_all-by-n_viewing_all matrix. ;
% The matrix-entry shift_z__(1+nX,1+nY) corresponds to the shift between ;
% viewing_angle_all_(1+nX) and viewing_angle_all_(1+nY), as described above. ;
% The L__ matrix holds the corresponding L-values. ;
%%%%%%%%;

if (nargin<1);
n_viewing_all = 2;
viewing_polar_a_all_ = [1*pi/12 ; 2*pi/12];
viewing_azimu_b_all_ = [-1*pi/12 ; 5*pi/12];
nX = 0; nY = 1;
aX = viewing_polar_a_all_(1+nX);
bX = viewing_azimu_b_all_(1+nX);
aY = viewing_polar_a_all_(1+nY);
bY = viewing_azimu_b_all_(1+nY);
a01 = -2*sin(bX-bY)*(cos(aX) + cos(aY));
a10 = -2*(cos(bX-bY)*(1 + cos(aX)*cos(aY)) + sin(aX)*sin(aY));
n_gamma_z = 32;
k_c_X_zd__ = get_k_c__(n_gamma_z,aX,bX); pole_X_d_ = get_pole_(aX,bX);
k_c_Y_zd__ = get_k_c__(n_gamma_z,aY,bY); pole_Y_d_ = get_pole_(aY,bY);
shift_z_ = transpose(2*pi*(0:n_gamma_z-1)/n_gamma_z);
L_quad_ = zeros(n_gamma_z,1);
L_form_ = zeros(n_gamma_z,1);
for ngamma_z=0:n_gamma_z-1;
shift_z = shift_z_(1+ngamma_z);
tmp_k_c_X_zd__ = get_k_c__(n_gamma_z,aX,bX,shift_z);
L_quad_(1+ngamma_z) = 2*pi*sum((tmp_k_c_X_zd__ - k_c_Y_zd__).^2,'all')/n_gamma_z;
L_form_(1+ngamma_z) = 2*pi*(2 + a10/2*cos(shift_z) + a01/2*sin(shift_z));
end;%for ngamma_z=0:n_gamma_z-1;
shift_z_one = atan2(a01,a10);
L_one = 2*pi*(2 + a10/2*cos(shift_z_one) + a01/2*sin(shift_z_one));
shift_z_two = shift_z_one + pi;
L_two = 2*pi*(2 + a10/2*cos(shift_z_two) + a01/2*sin(shift_z_two));
if (L_one<=L_two); L_opt = L_one; shift_z_opt = shift_z_one; end;
if (L_one> L_two); L_opt = L_two; shift_z_opt = shift_z_two; end;
[shift_z__,L__] = get_template_inplane_shift_0(n_viewing_all,viewing_polar_a_all_,viewing_azimu_b_all_);
disp(sprintf(' %% shift_z error: %0.16f',fnorm(shift_z__(1+nX,1+nY) - shift_z_opt)));
disp(sprintf(' %% L       error: %0.16f',fnorm(L__(1+nX,1+nY) - L_opt)));
%%%%%%%%;
c_X__ = colormap_beach(); n_c_X = size(c_X__,1);
c_Y__ = colormap_80s(); n_c_Y = size(c_Y__,1);
figure(1);clf;figbig;
subplot(1,2,2);
hold on;
plot(shift_z_,L_quad_,'kx-',shift_z_,L_form_,'ro-');
plot(shift_z_one,L_one,'gs',shift_z_two,L_two,'bh');
hold off;
xlabel('shift_z','Interpreter','none');
ylabel('L');
subplot(1,2,1);
%%%%%%%%;
hold on;
line( [0,pole_X_d_(1+0)] , [0,pole_X_d_(1+1)] , [0,pole_X_d_(1+2)] , 'LineWidth',2,'Color','k');
for ngamma_z=1:n_gamma_z-1;
nc = max(0,min(n_c_X-1,floor(n_c_X*ngamma_z/n_gamma_z)));
line( ...
 [k_c_X_zd__(0+ngamma_z,1+0);k_c_X_zd__(1+ngamma_z,1+0)] ...
,[k_c_X_zd__(0+ngamma_z,1+1);k_c_X_zd__(1+ngamma_z,1+1)] ...
,[k_c_X_zd__(0+ngamma_z,1+2);k_c_X_zd__(1+ngamma_z,1+2)] ...
,'LineWidth',2 ...
,'Color',c_X__(1+nc,:) ...
);
end;%for ngamma_z=1:n_gamma_z-1;
hold off;
%%%%%%%%;
hold on;
line( [0,pole_Y_d_(1+0)] , [0,pole_Y_d_(1+1)] , [0,pole_Y_d_(1+2)] , 'LineWidth',2,'Color','r');
for ngamma_z=1:n_gamma_z-1;
nc = max(0,min(n_c_Y-1,floor(n_c_Y*ngamma_z/n_gamma_z)));
line( ...
 [k_c_Y_zd__(0+ngamma_z,1+0);k_c_Y_zd__(1+ngamma_z,1+0)] ...
,[k_c_Y_zd__(0+ngamma_z,1+1);k_c_Y_zd__(1+ngamma_z,1+1)] ...
,[k_c_Y_zd__(0+ngamma_z,1+2);k_c_Y_zd__(1+ngamma_z,1+2)] ...
,'LineWidth',2 ...
,'Color',c_Y__(1+nc,:) ...
);
end;%for ngamma_z=1:n_gamma_z-1;
hold off;
%%%%%%%%;
axis equal; axis vis3d;
xlabel('x');ylabel('y');zlabel('z');
xlim(1.25*[-1,1]);
ylim(1.25*[-1,1]);
zlim(1.25*[-1,1]);
disp('returning'); return;
end;%if (nargin<1);

shift_z__ = zeros(n_viewing_all,n_viewing_all);
L__ = zeros(n_viewing_all,n_viewing_all);
for nX=0:n_viewing_all-1;
for nY=0:n_viewing_all-1;
aX = viewing_polar_a_all_(1+nX);
bX = viewing_azimu_b_all_(1+nX);
aY = viewing_polar_a_all_(1+nY);
bY = viewing_azimu_b_all_(1+nY);
a01 = -2*sin(bX-bY)*(cos(aX) + cos(aY));
a10 = -2*(cos(bX-bY)*(1 + cos(aX)*cos(aY)) + sin(aX)*sin(aY));
shift_z_one = atan2(a01,a10);
L_one = 2*pi*(2 + a10/2*cos(shift_z_one) + a01/2*sin(shift_z_one));
shift_z_two = shift_z_one + pi;
L_two = 2*pi*(2 + a10/2*cos(shift_z_two) + a01/2*sin(shift_z_two));
if (L_one<=L_two); L_opt = L_one; shift_z_opt = shift_z_one; end;
if (L_one> L_two); L_opt = L_two; shift_z_opt = shift_z_two; end;
shift_z__(1+nX,1+nY) = shift_z_opt;
L__(1+nX,1+nY) = L_opt;
end;%for nY=0:n_viewing_all-1;
end;%for nX=0:n_viewing_all-1;
 
function pole_d_ = get_pole_(polar_a,azimu_b);
a = polar_a;
b = azimu_b;
ca = cos(a); sa = sin(a);
cb = cos(b); sb = sin(b);
pole_d_ = [ +cb*sa , +sb*sa , +ca ];

function k_c_zd__ = get_k_c__(n_gamma_z,polar_a,azimu_b,shift_z);
if (nargin<4); shift_z = 0; end;
a = polar_a;
b = azimu_b;
ca = cos(a); sa = sin(a);
cb = cos(b); sb = sin(b);
gamma_z_ = transpose(2*pi*(0:n_gamma_z-1)/n_gamma_z);
cc_ = cos(gamma_z_ - shift_z); sc_ = sin(gamma_z_ - shift_z);
k_c_zd__ = [ +cb*ca*cc_ - sb*sc_ , +sb*ca*cc_ + cb*sc_ , -sa*cc_ ];


