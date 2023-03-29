function ...
[ ...
 template_k_c_0__ ...
,template_k_c_1__ ...
,template_k_c_2__ ...
,template_azimu_b__ ...
,template_polar_a__ ...
] = ...
get_template_single_ring_k_c_0( ...
 verbose ...
,n_viewing_all ...
,viewing_azimu_b_all_ ...
,viewing_polar_a_all_ ...
,n_w_max ...
);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

na=0;
if (nargin<1+na); verbose=[]; end; na=na+1;
if (nargin<1+na); n_viewing_all=[]; end; na=na+1;
if (nargin<1+na); viewing_azimu_b_all_=[]; end; na=na+1;
if (nargin<1+na); viewing_polar_a_all_=[]; end; na=na+1;
if (nargin<1+na); n_w_max=[]; end; na=na+1;

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

gamma_z_ = linspace(0,2*pi,n_w_max+1); gamma_z_ = transpose(gamma_z_(1:end-1));
cc_ = cos(gamma_z_); sc_ = sin(gamma_z_);
template_k_c_0__ = zeros(n_w_max,n_viewing_all);
template_k_c_1__ = zeros(n_w_max,n_viewing_all);
template_k_c_2__ = zeros(n_w_max,n_viewing_all);
for nviewing_all=0:n_viewing_all-1;
viewing_polar_a = viewing_polar_a_all_(1+nviewing_all); ca = cos(viewing_polar_a); sa = sin(viewing_polar_a);
viewing_azimu_b = viewing_azimu_b_all_(1+nviewing_all); cb = cos(viewing_azimu_b); sb = sin(viewing_azimu_b);
template_k_c_0__(:,1+nviewing_all) = (+cb*ca*cc_ - sb*sc_);
template_k_c_1__(:,1+nviewing_all) = (+sb*ca*cc_ + cb*sc_);
template_k_c_2__(:,1+nviewing_all) = (-sa*cc_            );
end;%for nviewing_all=0:n_viewing_all-1;
template_k_r01__ = sqrt(template_k_c_0__.^2 + template_k_c_1__.^2);
template_k_r012__ = sqrt(template_k_c_0__.^2 + template_k_c_1__.^2 + template_k_c_2__.^2);
template_azimu_b__ = atan2(template_k_c_1__,template_k_c_0__);
template_polar_a__ = atan2(template_k_r01__,template_k_c_2__);



