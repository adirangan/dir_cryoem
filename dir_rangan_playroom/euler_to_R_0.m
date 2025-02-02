function R = euler_to_R_0(e_);
% converts euler-angles e_ = [egammaz,epolara,eazimub] into a rotation-matrix R ;

egammaz = e_(1+0); cgammaz = cos(egammaz); sgammaz = sin(egammaz);
epolara = e_(1+1); cpolara = cos(epolara); spolara = sin(epolara);
eazimub = e_(1+2); cazimub = cos(eazimub); sazimub = sin(eazimub);

Razimub = [+cazimub -sazimub   0 ; ...
           +sazimub +cazimub   0 ; ...
            0        0         1 ];

Rgammaz = [+cgammaz -sgammaz   0 ; ...
           +sgammaz +cgammaz   0 ; ...
            0        0         1 ];

Rpolara = [+cpolara   0 +spolara ; ...
            0         1  0       ; ...
           -spolara   0 +cpolara ];

R = Razimub*Rpolara*Rgammaz;
