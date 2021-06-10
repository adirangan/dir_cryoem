function R = euler_to_R(e_);
% converts euler-angles e_ = [alpha,beta,gamma] into a rotation-matrix R ;

alpha = e_(1); ca = cos(alpha); sa = sin(alpha);
beta = e_(2); cb = cos(beta); sb = sin(beta);
gamma = e_(3); cg = cos(gamma); sg = sin(gamma);

Rg = [+cg -sg   0 ; ...
      +sg +cg   0 ; ...
        0   0   1 ];

Ra = [+ca -sa   0 ; ...
      +sa +ca   0 ; ...
        0   0   1 ];

Rb = [+cb   0 +sb ; ...
        0   1   0 ; ...
      -sb   0 +cb ];

R = Rg*Rb*Ra;
