function M_q_rwM___ = innerproduct_q_k_stretch_quad_stack__0(n_r,weight_p_r_,n_w_,n_M,M_q__) ;
%%%%%%%%;
% This function rearranges the fourier-bessel coefficients of the image-stack M_q__. ;
% This rearrangement multiplies the various values of M_q__ by the square-root of the associated quadrature weight ;
% (assumed to depend only on radius, and stored in weight_p_r_). ;
% This allows to use a standard dot-product for any future computations of the l2-integral ;
% between two different image-stacks (each presumably prepared using a function like this). ;
%%%%%%%%;
% Inputs: ;
% n_r : integer number of image-rings (i.e., number of radii in k-space). ;
% weight_p_r_ : double-array of size n_r. ;
%    weight_p_r_(1+nr) contains the integration weight associated with r-index nr. ;
% n_w_ : integer-array of size n_r. ;
%    n_w_(1+nr) contains the number of points stored on the image-ring with r-index nr. ;
% n_M : integer number of images. ;
% M_q__ : complex-array of size (n_w_sum,n_M) storing the fourier-bessel coefficients of the image-stack. ;
%    M_q__(1+nw,1+nM) contains the fourier-bessel coefficient for q-index nw and m-index nM. ;
%    We assume that M_q_ := M_q__(:,1+nM) stores the fourier-bessel coefficients for image nM. ;
%    We assume that M_q_ is a complex-array of size n_w_sum = sum(n_w_). ;
%    That is to say, we assume that M_q_ is 'unrolled', ;
%    with the coefficients associated with each ring stored in sequence. ;
%    Note that there will be nw=n_w_(1+nr) such q-indices at each image-ring with r-index nr. ;
%    These coefficients associated with ring nr will be stored in M_q_(1 + n_w_csum + [0:n_w-1]), ;
%    Where n_w_csum := 1+n_w_csum_(1+nr) and n_w := n_w_(1+nr). ;
%    Note that, for the input, we assume that interpretation of the nw requires knowledge of n_w (which is dependent on nr). ;
%    In other words, we assume the convention that nw==0 corresponds to the constant term, whereas ;
%    nw==+1 and nw==n_w-1 correspond to frequencies of +1 and -1 respectively, ;
%    nw==+2 and nw==n_w-2 correspond to frequencies of +2 and -2, and so forth. ;
%%%%%%%%;
% Outputs: ;
% M_q_rwM___ : complex-array of size(n_r,n_w_max,n_M) storing the sqrt-weighted fourier-bessel coefficients of the image-stack. ;
%    M_q_rw__ := M_q_rwM___(:,:,1+nM) contains the sqrt-weighted fourier-bessel coefficients for image nM. ;
%    M_q_rw__(1+nr,1+nw) stores coefficient nw associated with ring nr. ;
%    Note that, for the output, we assume that interpretation of the nw requires only knowledge of n_w_max (which is independent on nr). ;
%    To accomplish this we rearrange the indices nw from each ring so that ;
%    nw==+1 and nw==n_w_max-1 are consistently associated with frequencies +1 and -1, respectively, ;
%    nw==+2 and nw==n_w_max-2 are consistently associated with frequencies +2 and -2, and so forth. ;
%%%%%%%%;
% Notes: ;
% Assumes that n_w_ is a list of even integers. ;
%%%%%%%%;
n_w_ = n_w_(:);
n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);
M_q_rwM___ = zeros(n_r,n_w_max,n_M);
for nM=0:n_M-1;
for nr=0:n_r-1;
n_w = n_w_(1+nr);
n_w_2 = round(n_w/2);
n_w_csum = n_w_csum_(1+nr);
dAn = weight_p_r_(1+nr)*(2*pi)/max(1,n_w);
M_q_rwM___(1+nr,1+(0:n_w_2-1),1+nM) = M_q__(1+n_w_csum_(1+nr)+(0:n_w_2-1),1+nM)*sqrt(dAn);
M_q_rwM___(1+nr,1+(n_w_2+1:n_w-1)-n_w+n_w_max,1+nM) = M_q__(1+n_w_csum_(1+nr)+(n_w_2+1:n_w-1),1+nM)*sqrt(dAn);
end;%for nr=0:n_r-1;
end;%for nM=0:n_M-1;
