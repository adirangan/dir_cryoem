function M_k_q_rwl___ = innerproduct_q_k_stretch_quad_stack___0(n_k_p_r,n_w_,M_k_q_,l_max);
%%%%%%%%;
% This function rearranges (and replicates) the fourier-bessel coefficients of the image M_k_q_. ;
% Note that this rearrangement does not multiply the various values of M_k_q_ by any quadrature-weight. ;
% (indeed, we pass in 'n_w_/(2*pi)' to innerproduct_q_k_stretch_quad_stack__0 so the effective quadrature-weight is 1). ;
% This rearrangement+replication is intended for an innerproduct-calculation of the form described in Eq. 33 (or, similarly Eq. 58) in ;
% (Rangan et al. Factorization of the translation kernel for fast rigid image alignment, Inverse Problems, 36(2): 024001). ;
% Specifically, we reorganize and replicate each fourier-bessel coefficient from the original image M_k_q_ ;
% so that a coefficient with q-index nw is associated with multiple positions in the output-array M_k_q_rwl___ ;
% associated with the frequencies (q-l_val) indicated in Eq. 33 and/or Eq. 58. ;
% This association must be conducted with care; due to discretization and bandlimiting ;
% the q-index associated with the frequency (q-l_val) is not necessarily the frequency associated with the q-index (nw - (l_val+l_max)). ;
%%%%%%%%;
% Inputs: ;
% n_k_p_r : integer number of image-rings (i.e., number of radii in k-space). ;
% n_w_ : integer-array of size n_k_p_r. ;
%    n_w_(1+nk_p_r) contains the number of points stored on the image-ring with r-index nk_p_r. ;
% M_k_q_ : complex-array of size (n_w_sum,n_M) storing the fourier-bessel coefficients of the image. ;
%    We assume that M_k_q_ is a complex-array of size n_w_sum = sum(n_w_). ;
%    That is to say, we assume that M_k_q_ is 'unrolled', ;
%    with the coefficients associated with each ring stored in sequence. ;
%    Note that there will be nw=n_w_(1+nk_p_r) such q-indices at each image-ring with r-index nk_p_r. ;
%    These coefficients associated with ring nk_p_r will be stored in M_k_q_(1 + n_w_csum + [0:n_w-1]), ;
%    Where n_w_csum := 1+n_w_csum_(1+nk_p_r) and n_w := n_w_(1+nk_p_r). ;
%    Note that, for the input, we assume that interpretation of the nw requires knowledge of n_w (which is dependent on nk_p_r). ;
%    In other words, we assume the convention that nw==0 corresponds to the constant term, whereas ;
%    nw==+1 and nw==n_w-1 correspond to frequencies of +1 and -1 respectively, ;
%    nw==+2 and nw==n_w-2 correspond to frequencies of +2 and -2, and so forth. ;
%%%%%%%%;
% Outputs: ;
% M_k_q_rwl___ : complex-array of size(n_k_p_r,n_w_max,1+2*l_max) storing the fourier-bessel coefficients of the image. ;
%    M_k_q_rw__ := M_k_q_rwl___(:,:,1+l_val+l_max) contains the fourier-bessel coefficients for l_val (see Eqs. 58-60). ;
%    M_k_q_rw__(1+nk_p_r,1+nw) stores coefficients associated with ring nk_p_r, each shifted by l_val. ;
%    Note that we do not shift all the nw from the original image by l_val. ;
%    Instead, to maintain consistency with the bandlimit for which the fourier-bessel coefficients were originally defined, ;
%    we ignore any coefficients for which the shifted nw lies outside the range [-n_w/2+1,+n_w/2-1] (defining n_w:=n_w_(1+nk_p_r)). ;
%    Finally, we assume that interpretation of the retained nw requires only knowledge of n_w_max (which is independent on nk_p_r). ;
%    To accomplish this we rearrange the output indices nw on each ring so that ;
%    nw==+1 and nw==n_w_max-1 are consistently associated with output frequencies +1 and -1, respectively, ;
%    nw==+2 and nw==n_w_max-2 are consistently associated with output frequencies +2 and -2, and so forth. ;
%%%%%%%%;
% Notes: ;
% Assumes that n_w_ is a list of even integers. ;
%%%%%%%%;
n_w_sum = sum(n_w_);
M_k_q__ = zeros(n_w_sum,1+2*l_max);
ic=0;
for nk_p_r=0:n_k_p_r-1;
for nw=0:n_w_(1+nk_p_r)-1;
n_w_t = round(1.0d0*n_w_(1+nk_p_r)/2.0d0);
for l_val=-l_max:+l_max;
nwc = nw;
if (nwc>=n_w_t); nwc = nwc - n_w_(1+nk_p_r); end;%if;
flag_ict_overflow = 0;
nwd = nwc + l_val;
if (abs(nwd)<n_w_t); nwt = periodize(nwd,0,n_w_(1+nk_p_r));  else; nwt = 0; flag_ict_overflow = 1; end;%if;
ict = ic-nw+nwt;
if (flag_ict_overflow==0); M_k_q__(1+ic,1+l_max+l_val) = M_k_q_(1+ict); end;%if;
end;%for l_val=-l_max:+l_max;
ic = ic + 1;
end;%for nw=0:n_w_(1+nk_p_r)-1;
end;%for nk_p_r=0:n_k_p_r-1;
M_k_q_rwl___ = innerproduct_q_k_stretch_quad_stack__0(n_k_p_r,n_w_/(2*pi),n_w_,1+2*l_max,M_k_q__);
