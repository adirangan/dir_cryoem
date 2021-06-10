function X_best = test_spharm_14(n_k,k_,n_l_,a_,b_);
% tests registration between molecule_A and molecule_B via 'alternating' ;
% between angle-registration and displacement-registration. ;
% As this alternation proceeds the 'grid-spacing' decreases ;
% (both in terms of angle beta and displacement delta). ;
% ;
% The angular-registration uses the function: ;
% register_spharm_to_spharm_angle_0(verbose,n_k,k_,n_l_,a_,c_,n_beta,beta_);
% whereas the displacement-registration uses the function: ;
% register_spharm_to_spharm_delta_0(verbose,n_k,k_,n_l_,a_,c_,n_delta,delta__,sample_d);
% ;
% After each angular (or displacement) registration we choose the angle (or displacement) that ;
% produces the largest inner-product, and then apply that transformation to the second molecule (i.e., b_). ;
% ;
% When no inputs are passed we import two spherical harmonic representations (generated by kspacegrid_to_model): ;
% molecule_A: modsph_A_ori = spiral ;
% molecule_B: modsph_B_ori = spiral with twisted tail ;
% ;
% Inputs: ;
% n_k = integer maximum k ;
% k_ = real array of length n_k; k_(nk) = k_value for shell nk ;
% n_l_ = integer array of length n_k; n_l_(nk) = spherical harmonic order on shell nk; n_l_(nk) corresponds to n_lm_(nk) = (n_l_(nk)+1)^2 coefficients ;
% a_ = complex array of length \sum_{nk} (n_lm_(nk)+1)^2 ; coefficients are ordered in a row, with m varying quickly and l varying slowly ;
% b_ = complex array of length \sum_{nk} (n_lm_(nk)+1)^2 ; coefficients are ordered in a row, with m varying quickly and l varying slowly ;
% a_ corresponds to molecule_a, b_ to molecule_b ;
% ;
% Outputs: ;
% X_best = best correlation found. ;
% ;
% test with: ;
%{
  test_spharm_14();
  %}

path(path,'/data/rangan/dir_cryoem/nufftall-1.33/');

verbose=1; 
quad_flag=1; % to use quadratic interpolation. ;

if nargin<5;
xnodesr_ = MDA_read_r8('./dir_mdaT/xnodesr_.mda');
isph_start_ = MDA_read_i4('./dir_mdaT/isph_start_.mda');
nterms_sph_ = MDA_read_i4('./dir_mdaT/nterms_sph_.mda');
modsph_A_ori_ = MDA_read_c16('./dir_mdaT/modsph_A_ori_.mda');
modsph_B_ori_ = MDA_read_c16('./dir_mdaT/modsph_B_ori_.mda');
n_k = length(isph_start_);
k_ = xnodesr_;
n_l_ = nterms_sph_;
n_lm_ = (n_l_+1).^2;
a_ = modsph_A_ori_;
b_ = modsph_B_ori_;
end;%if nargin<4;

disp(sprintf(' %% indices for counting arrays'));
n_lm_ = (n_l_+1).^2;
k_max = k_(end);
n_l_max = n_l_(end);
m_max_ = -n_l_max : +n_l_max;
n_m_max = length(m_max_);
sample_d = 1.0;

disp(sprintf(' %% setting up initial array of rotations alpha_, beta_ and gamma__'));
%n_beta = 40; beta_ = linspace(-pi,pi,n_beta+1); beta_ = beta_(1:end-1); n_beta = length(beta_);
n_beta = 14; beta_ = linspace(-pi,pi,n_beta+1); beta_ = beta_(1:end-1); n_beta = length(beta_);
n_alpha = n_m_max; alpha_ = linspace(0,2*pi,n_alpha+1); alpha_ = alpha_(1:end-1);
n_gamma = n_m_max; gamma_ = linspace(0,2*pi,n_gamma+1); gamma_ = gamma_(1:end-1);

disp(sprintf(' %% setting up initial array of translations delta__'));
n_delta_x = 5; n_delta_y = 5; n_delta_z = 5; delta_max = 0.1625;
delta_x_ = delta_max*linspace(-1,1,n_delta_x); 
delta_y_ = delta_max*linspace(-1,1,n_delta_y); 
delta_z_ = delta_max*linspace(-1,1,n_delta_z);
[Delta_y_,Delta_x_,Delta_z_] = meshgrid(delta_y_,delta_x_,delta_z_);
n_delta = numel(Delta_x_);
delta__ = [Delta_x_(:) , Delta_y_(:) , Delta_z_(:)];

disp(sprintf(' %% determining l2-norm of each molecule'));
tmp_x_aa = register_spharm_to_spharm_0(verbose,n_k,k_,n_l_,a_,a_);
tmp_x_ab = register_spharm_to_spharm_0(verbose,n_k,k_,n_l_,a_,b_);
tmp_x_bb = register_spharm_to_spharm_0(verbose,n_k,k_,n_l_,b_,b_);
disp(sprintf(' %% tmp_x_aa %0.16f, tmp_x_ab %0.16f, tmp_x_bb %0.16f, correlation %0.16f',tmp_x_aa,tmp_x_ab,tmp_x_bb,tmp_x_ab./sqrt(tmp_x_aa)./sqrt(tmp_x_bb)));

outer_iteration_max = 5;
for outer_iteration=1:outer_iteration_max;
disp(sprintf(' %% outer iteration %d/%d',outer_iteration,outer_iteration_max));

disp(sprintf(' %% determining best angle (no delta)'));
%%%%%%%%%%%%%%%% ;
angle_pre_ = zeros(1,3);
iteration_max = 4; iteration=0; continue_flag=1; tmp_x_pre = 0;
%%%%%%%%%%%%%%%% ;
while (continue_flag);
%%%%%%%%%%%%%%%% ;
tic;
[c_] = rotate_spharm_to_spharm_0(verbose,n_k,k_,n_l_,b_,angle_pre_);
[X0_] = register_spharm_to_spharm_angle_0(verbose,n_k,k_,n_l_,a_,c_,n_beta,beta_,n_alpha,alpha_,n_gamma,gamma_);
X0_ = real(X0_./sqrt(tmp_x_aa)./sqrt(tmp_x_bb));
t_0 = toc; 
if (verbose>1); disp(sprintf(' %% register_spharm_to_spharm_angle_0 time %0.2f',t_0)); end;
%%%%%%%%%%%%%%%% ;
[max_pos_v] = max(real(X0_(:)));
[max_pos_,max_pos_q_,tmp_F,tmp_DF,tmp_DDF] = quadratic_3d_interpolation_wrapper_0(real(X0_(:)),[n_alpha,n_gamma,n_beta],[0,0,0]); % We could run this search using a periodic boundary set by [1,1,1] as the final argument. However, such a periodic boundary will not apply during later iterations when alpha_, gamma_ and beta_ are more limited in scope. ;
assert(real(X0_(1+max_pos_(1),1+max_pos_(2),1+max_pos_(3)))==max_pos_v);
angle_pos_o_(1) = alpha_(1+max_pos_(1));
angle_pos_o_(3) = gamma_(1+max_pos_(2));
angle_pos_o_(2) = beta_(1+max_pos_(3));
[d_] = rotate_spharm_to_spharm_0(verbose,n_k,k_,n_l_,c_,angle_pos_o_);
max_pos_v_o = register_spharm_to_spharm_0(verbose,n_k,k_,n_l_,a_,d_)./sqrt(tmp_x_aa)./sqrt(tmp_x_bb);
angle_pos_q_(1) = alpha_(1+max_pos_(1)) + mean(diff(alpha_))*max_pos_q_(1);
angle_pos_q_(3) = gamma_(1+max_pos_(2)) + mean(diff(gamma_))*max_pos_q_(2);
angle_pos_q_(2) = beta_(1+max_pos_(3)) + mean(diff(beta_))*max_pos_q_(3);
[d_q_] = rotate_spharm_to_spharm_0(verbose,n_k,k_,n_l_,c_,angle_pos_q_);
max_pos_v_q = register_spharm_to_spharm_0(verbose,n_k,k_,n_l_,a_,d_q_)./sqrt(tmp_x_aa)./sqrt(tmp_x_bb);
%%%%%%%%%%%%%%%% ;
ier=0;
if (quad_flag & real(max_pos_v_q)>real(max_pos_v_o));
[angle_tmp_,ier] = R_to_euler(euler_to_R(angle_pos_q_)*euler_to_R(angle_pre_));
tmp_x_pos = max_pos_v_q;
if (verbose>0); disp(sprintf(' %% inner iteration %d/%d; angle_pre_ [%0.3f,%0.3f,%0.3f] -o- angle_pos_q_ [%0.3f,%0.3f,%0.3f] = angle_tmp_ [%0.3f,%0.3f,%0.3f]: time %0.3f; max_pos_v_q %0.16f - max_pos_v %0.16f = %0.16f',iteration,iteration_max,angle_pre_,angle_pos_q_,angle_tmp_,t_0,real(max_pos_v_q),real(max_pos_v),real(max_pos_v_q-max_pos_v))); end;
end;%if (real(max_pos_v_q)>real(max_pos_v_o));
if (~quad_flag | real(max_pos_v_q)<=real(max_pos_v_o));
[angle_tmp_,ier] = R_to_euler(euler_to_R(angle_pos_o_)*euler_to_R(angle_pre_));
tmp_x_pos = max_pos_v_o;
if (verbose>0); disp(sprintf(' %% inner iteration %d/%d; angle_pre_ [%0.3f,%0.3f,%0.3f] -o- angle_pos_o_ [%0.3f,%0.3f,%0.3f] = angle_tmp_ [%0.3f,%0.3f,%0.3f]: time %0.3f; max_pos_v_o %0.16f - max_pos_v %0.16f = %0.16f',iteration,iteration_max,angle_pre_,angle_pos_o_,angle_tmp_,t_0,real(max_pos_v_o),real(max_pos_v),real(max_pos_v_o-max_pos_v))); end;
end;%if (real(max_pos_v_q)<=real(max_pos_v_o));
angle_pre_ = angle_tmp_;
%%%%%%%%%%%%%%%% ;
continue_flag = ((iteration<iteration_max) & (norm(tmp_x_pre-tmp_x_pos)>1e-6) & ~ier);
iteration = iteration+1; tmp_x_pre = tmp_x_pos;
%%%%%%%%%%%%%%%% ;
end;%while;
%%%%%%%%%%%%%%%% ;
disp(sprintf(' %% baking in best angle: [%0.3f,%0.3f,%0.3f], best correlation so far: %0.16f',angle_pre_,tmp_x_pre));
b_ = rotate_spharm_to_spharm_0(verbose,n_k,k_,n_l_,b_,angle_pre_);

disp(sprintf(' %% determining best delta (no angle)'));
%%%%%%%%%%%%%%%% ;
delta_pre_ = zeros(1,3);
iteration_max = 4; iteration=0; continue_flag=1; tmp_x_pre = 0;
%%%%%%%%%%%%%%%% ;
while (continue_flag);
%%%%%%%%%%%%%%%% ;
tic;
[c_] = transf_spharm_to_spharm_0(verbose,n_k,k_,n_l_,b_,delta_pre_,sample_d);
disp('here: '); 1,find(~isfinite(b_)),; 2,find(~isfinite(c_)),;
[X0_] = register_spharm_to_spharm_delta_0(verbose,n_k,k_,n_l_,a_,c_,n_delta,delta__,sample_d);
X0_ = real(X0_./sqrt(tmp_x_aa)./sqrt(tmp_x_bb));
t_0 = toc; 
%%%%%%%%%%%%%%%% ;
[max_pos_v,max_ij] = max(real(X0_(:)));
delta_pos_o_(1) = delta__(max_ij,1);
delta_pos_o_(2) = delta__(max_ij,2);
delta_pos_o_(3) = delta__(max_ij,3);
[d_o_] = transf_spharm_to_spharm_0(verbose,n_k,k_,n_l_,c_,delta_pos_o_,sample_d);
max_pos_v_o = register_spharm_to_spharm_0(verbose,n_k,k_,n_l_,a_,d_o_)./sqrt(tmp_x_aa)./sqrt(tmp_x_bb);
[max_pos_,max_pos_q_,tmp_F,tmp_DF,tmp_DDF] = quadratic_3d_interpolation_wrapper_0(real(X0_(:)),[n_delta_x,n_delta_y,n_delta_z],[0,0,0]);
assert(real(X0_(1+max_pos_(1)+max_pos_(2)*n_delta_x+max_pos_(3)*n_delta_x*n_delta_y))==max_pos_v);
delta_pos_q_(1) = delta_x_(1+max_pos_(1)) + mean(diff(delta_x_))*max_pos_q_(1);
delta_pos_q_(2) = delta_y_(1+max_pos_(2)) + mean(diff(delta_y_))*max_pos_q_(2);
delta_pos_q_(3) = delta_z_(1+max_pos_(3)) + mean(diff(delta_z_))*max_pos_q_(3);
[d_q_] = transf_spharm_to_spharm_0(verbose,n_k,k_,n_l_,c_,delta_pos_q_,sample_d);
max_pos_v_q = register_spharm_to_spharm_0(verbose,n_k,k_,n_l_,a_,d_q_)./sqrt(tmp_x_aa)./sqrt(tmp_x_bb);
%%%%%%%%%%%%%%%% ;
if (quad_flag & real(max_pos_v_q)>real(max_pos_v_o));
delta_tmp_ = delta_pre_ + delta_pos_q_;
tmp_x_pos = max_pos_v_q;
if (verbose>0); disp(sprintf(' %% inner iteration %d/%d; delta_pre_ [%0.3f,%0.3f,%0.3f] + delta_pos_q_ [%0.3f,%0.3f,%0.3f] = delta_tmp_ [%0.3f,%0.3f,%0.3f]: time %0.3f; max_pos_v_q %0.16f - max_pos_v %0.16f = %0.16f',iteration,iteration_max,delta_pre_,delta_pos_q_,delta_tmp_,t_0,real(max_pos_v_q),real(max_pos_v),real(max_pos_v_q-max_pos_v))); end;
end;%if (real(max_pos_v_q)>real(max_pos_v_o));
if (~quad_flag | real(max_pos_v_q)<=real(max_pos_v_o));
delta_tmp_ = delta_pre_ + delta_pos_o_;
tmp_x_pos = max_pos_v_o;
if (verbose>0); disp(sprintf(' %% inner iteration %d/%d; delta_pre_ [%0.3f,%0.3f,%0.3f] + delta_pos_o_ [%0.3f,%0.3f,%0.3f] = delta_tmp_ [%0.3f,%0.3f,%0.3f]: time %0.3f; max_pos_v_o %0.16f - max_pos_v %0.16f = %0.16f',iteration,iteration_max,delta_pre_,delta_pos_o_,delta_tmp_,t_0,real(max_pos_v_o),real(max_pos_v),real(max_pos_v_o-max_pos_v))); end;
end;%if (real(max_pos_v_q)<=real(max_pos_v_o));
delta_pre_ = delta_tmp_;
%%%%%%%%%%%%%%%% ;
continue_flag = ((iteration<iteration_max) & (norm(tmp_x_pre-tmp_x_pos)>1e-6));
iteration = iteration+1; tmp_x_pre = tmp_x_pos;
%%%%%%%%%%%%%%%% ;
end;%while;
%%%%%%%%%%%%%%%% ;
disp(sprintf(' %% baking in best delta: [%0.3f,%0.3f,%0.3f], best correlation so far: %0.16f',delta_pre_,tmp_x_pre));
b_ = transf_spharm_to_spharm_0(verbose,n_k,k_,n_l_,b_,delta_pre_,sample_d);

disp(sprintf(' %% increasing angle and delta resolution (by reducing ranges)'));
%%%%%%%%%%%%%%%% ;
tmp_d = 2^(0.5*outer_iteration);
%%%%%%%%%%%%%%%% ;
n_beta = 5; beta_ = linspace(-pi/tmp_d,+pi/tmp_d,n_beta);
n_alpha = n_m_max; alpha_ = linspace(-pi/tmp_d,+pi/tmp_d,n_alpha);
n_gamma = n_m_max; gamma_ = linspace(-pi/tmp_d,+pi/tmp_d,n_gamma);
%%%%%%%%%%%%%%%% ;
n_delta_x = 5; n_delta_y = 5; n_delta_z = 5; delta_max = 0.1625/tmp_d;
delta_x_ = delta_max*linspace(-1,1,n_delta_x); 
delta_y_ = delta_max*linspace(-1,1,n_delta_y); 
delta_z_ = delta_max*linspace(-1,1,n_delta_z);
[Delta_y_,Delta_x_,Delta_z_] = meshgrid(delta_y_,delta_x_,delta_z_);
n_delta = numel(Delta_x_);
delta__ = [Delta_x_(:) , Delta_y_(:) , Delta_z_(:)];

end;%for outer_iteration=1:outer_iteration_max;

X_best = tmp_x_pre;
