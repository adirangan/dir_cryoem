function e_ = R_to_euler_0(R);
% converts rotation-matrix R (with determinant +1) into euler-angles e_ = [alpha,beta,gamma] ;
% using direct calculation. ;
% ;
% test with: ;
%{
R_to_euler_0();
 %}

if (nargin<1);
e_0in_ = [2*pi*rand(),1*pi*rand(),2*pi*rand()];
R = euler_to_R(e_0in_);
continue_flag=1; iteration = 0; iteration_max = 40;
while (continue_flag & iteration<iteration_max);
[e_sgd_,ier] = R_to_euler(R);
e_dir_ = R_to_euler_0(R);
continue_flag = ier;
iteration = iteration + 1;
end;% while;
e_sgd_ = periodize(e_sgd_,0,2*pi); 
e_dir_ = periodize(e_dir_,0,2*pi); 
disp(sprintf(' %% sgd: angle error %0.6f %0.6f %0.6f',e_0in_-e_sgd_));
disp(sprintf(' %% sgd: R error: %0.6f',norm(euler_to_R(e_0in_) - euler_to_R(e_sgd_))));
disp(euler_to_R(e_0in_) - euler_to_R(e_sgd_));
disp(sprintf(' %% dir: angle error %0.6f %0.6f %0.6f',e_0in_-e_dir_));
disp(sprintf(' %% dir: R error: %0.6f',norm(euler_to_R(e_0in_) - euler_to_R(e_dir_))));
disp(euler_to_R(e_0in_) - euler_to_R(e_dir_));
return;
end;% if (nargin<1);

verbose=0;

%%%%%%%%;
% With this notation we assume: ;
% R__ = Rz(azimu_b)*Ry(polar_a)*Rz(gamma_z), ;
% for which the resulting matrix is: ;
% [ +cb*ca*cc - sb*sc , -cb*ca*sc -sb*cc , +cb*sa ];
% [ +sb*ca*cc + cb*sc , -sb*ca*sc +cb*cc , +sb*sa ];
% [ -sa*cc            , +sa*sc           , +ca    ];
% as used in get_template_1.m ;
%%%%%%%%;

polar_a = acos(+R(3,3));
azimu_b = atan2(+R(2,3),+R(1,3));
gamma_z = atan2(+R(3,2),-R(3,1));

azimu_b = periodize(azimu_b,0,2*pi);
polar_a = periodize(polar_a,-1*pi,+1*pi);
gamma_z = periodize(gamma_z,0,2*pi);
if polar_a< 0; polar_a = -polar_a; azimu_b = periodize(azimu_b-pi,0,2*pi); gamma_z = periodize(gamma_z-pi,0,2*pi); end;

e_ = [gamma_z,polar_a,azimu_b];
