function [C_q_,V_r__] = innerproduct_q__k_qcrop_svdr_0(flag_S_vs_M,svd_r_max,n_svd_r,svd_r_,n_svd_l,svd_l_,svd_s_,svd_V_r_,n_r,grid_p_,n_w_,n_A,T_q_,M_q_,C_q_);
;
% Assumes that M_q_ is the same size and dimensions as T_q_.;
% Assumes quasi-uniform polar-grid defined via n_r , .. , n_A.;
% Uses svd-expansion defined via n_svd_r , .. , svd_V_r_.;
% Assumes that C_q_ is large enough to hold all n_w_max = n_w_(n_r-1) ;
% modes for each of the n_svd_l terms in the svd-expansion;
% (assuming of course that n_w_(n_r-1) is the largest value within n_w_).;
% modes in C_q_ are stored in the order: (mode 0 , mode 1 , ... , mode -1).;
% The mode-k for term-l is stored in C_q_(l + k*n_svd_l).;
% The logical flag_S_vs_M determines the sign of the complex exponential.;
% flag_S_vs_M .eqv. .true. --> transformation applied to S, use +.;
% flag_S_vs_M .eqv. .false. --> transformation applied to M, use -.;
% 
% Note: removes q == n_w_(nr)/2 term. ;

verbose = 0;
warning_flag = 1;

if (verbose>0);
disp(sprintf(' %% [entering innerproduct_q__k_qcrop_svdr_0] n_r %d',n_r));
end;%if;

V_r_ = zeros(n_svd_l*n_r,1);
V_r__ = zeros(n_svd_l*n_r,1);

svd_r_m = svd_r_max / 2.0;
svd_r_c = svd_r_m;
n_w_max = n_w_(1+n_r-1);
if (verbose>0);
disp(sprintf(' %% n_w_max %d',n_w_max));
end;%if;
C_q_ = zeros(n_svd_l*n_w_max,1);
for nr=0:n_r-1;
if ((grid_p_(1+nr)>svd_r_max) & (warning_flag));
disp(sprintf(' %% Warning, grid_p_(1+nr) %0.6f > svd_r_max %0.6f',grid_p_(1+nr),svd_r_max));
end;%if;
svd_r = (grid_p_(1+nr) - svd_r_m)/svd_r_c;
for nl=0:n_svd_l-1;
V_r_(1+nl+nr*n_svd_l) = polyval_r8_reverse_0(n_svd_r,svd_V_r_(1+0+nl*n_svd_r + (0:n_svd_r-1)),1,svd_r);
V_r__(1+nl+nr*n_svd_l) = V_r_(1+nl+nr*n_svd_l);
end;%for nl=0:n_svd_l-1;
end;%for nr=0:n_r-1;
for nl=0:n_svd_l-1;
for nw=0:n_w_max-1;
C_q_(1+nl+nw*n_svd_l) = 0.0;
end;%for nw=0:n_w_max-1;
end;%for nl=0:n_svd_l-1;
C_q = 0.0;
ic = 0;
for nr=0:n_r-1;
if (verbose>2); disp(sprintf(' %% nr %d',nr)); end;
if (nr>0);
R_pre = 0.5*(grid_p_(1+nr-1) + grid_p_(1+nr));
 else;
R_pre = grid_p_(1+0);
end;%if;
if (nr<n_r-1);
R_pos = 0.5*(grid_p_(1+nr+1) + grid_p_(1+nr));
 else;
R_pos = grid_p_(1+n_r-1);
end;%if;
dr = R_pos - R_pre;
%    We set the zero-mode to zero;
if (grid_p_(1+nr)<=0.0d0);
dr = 0.0d0;
end;%if;
dw = 2*pi/(1.0d0*max(1,n_w_(1+nr)));
dA = (R_pre*dr + (dr^2)/2)*dw;
%    We assume that the fourier basis is orthonormal (not merely orthogonal);
dAn = dA;
n_w_t = floor(1.0d0*n_w_(1+nr)/2.0d0);
ic_store = ic;
for nl=0:n_svd_l-1;
D_V_r = V_r_(1+nl+nr*n_svd_l);
D_s = svd_s_(1+nl);
I_l = svd_l_(1+nl);
ic = ic_store;
for nw=0:n_w_(1+nr)-1;
%%%%%%%%;
flag_ic0_overflow = 0;
nwc = nw;
if (nwc>=n_w_t);
nwc = nwc - n_w_(1+nr);
end;%if;
if (abs(nwc)<n_w_t);
flag_ic0_overflow = 0;
 else;
flag_ic0_overflow = 1;
end;%if;
%%%%%%%%;
nwc = nw;
if (nwc>=n_w_t);
nwc = nwc - n_w_(1+nr);
end;%if
flag_ict_overflow = 0;
nwd = nwc + I_l;
if (abs(nwd)<n_w_t);
nwt = periodize(nwd,0,n_w_(1+nr));
 else;
nwt = 0;
flag_ict_overflow = 1;
end;%if
flag_icr_overflow = 0;
nwd = nwc - I_l;
if (abs(nwd)<n_w_t);
nwr = periodize(nwd,0,n_w_(1+nr));
 else;
nwr = 0;
flag_icr_overflow = 1;
end;%if
ict = ic-nw+nwt;
icr = ic-nw+nwr;
%%%%%%%%;
if (flag_ic0_overflow==0 & flag_ict_overflow==0);
if (nw>n_w_t);
nw_fix = nw - n_w_(1+nr) + n_w_max;
if (flag_S_vs_M==1) C_q = conj(D_s*D_V_r*T_q_(1+ict))*M_q_(1+ic); end;
if (flag_S_vs_M==0) C_q = conj(D_s*D_V_r*T_q_(1+ic))*M_q_(1+ict); end;
nw_C = nl + nw_fix*n_svd_l;
C_q_(1+nw_C) = C_q_(1+nw_C) + C_q*dAn;
elseif (nw==n_w_t);
% do nothing. ;
 else;
nw_fix = nw;
if (flag_S_vs_M==1) C_q = conj(D_s*D_V_r*T_q_(1+ict))*M_q_(1+ic); end;
if (flag_S_vs_M==0) C_q = conj(D_s*D_V_r*T_q_(1+ic))*M_q_(1+ict); end;
nw_C = nl + nw_fix*n_svd_l;
C_q_(1+nw_C) = C_q_(1+nw_C) + C_q*dAn;
end;%if
end;%if
ic = ic + 1;
end;%for nw=0:n_w_(1+nr)-1;
end;%for nl=0:n_svd_l-1;
end;%for nr=0:n_r-1;

if (verbose>0);
disp(sprintf(' %% [finished innerproduct_q__k_qcrop_svdr_0]'));
end;%if;
