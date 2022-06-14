function ...
image_SR0A_gamma_z_0( ...
 n_M ...
,n_w_max ...
,nw_true_M_ ...
,nw_reco_M_ ...
,x_0 ...
,x_1 ...
,x_r ...
,d_r ...
,c__ ...
);

na=0;
if (nargin<1+na); n_M=[]; end; na=na+1;
if (nargin<1+na); n_w_max=[]; end; na=na+1;
if (nargin<1+na); nw_true_M_=[]; end; na=na+1;
if (nargin<1+na); nw_reco_M_=[]; end; na=na+1;
if (nargin<1+na); x_0=[]; end; na=na+1;
if (nargin<1+na); x_1=[]; end; na=na+1;
if (nargin<1+na); x_r=[]; end; na=na+1;
if (nargin<1+na); d_r=[]; end; na=na+1;
if (nargin<1+na); c__=[]; end; na=na+1;

if isempty(x_0); x_0 = 0; end;
if isempty(x_1); x_1 = 0; end;
if isempty(x_r); x_r = 1; end;
if isempty(d_r); d_r = 2*pi*x_r/n_w_max/2; end;
if isempty(c__); c__ = colormap_phase; end;
n_c = size(c__,1);

n_upb = n_w_max*n_M;
pos_g_ = zeros(n_upb,1);
pos_r_ = zeros(n_upb,1);
pos_dg_ = zeros(n_upb,1);
pos_dr_ = zeros(n_upb,1);
pos_i_ = zeros(n_upb,1);
n_w_coarse = max(1,min(n_w_max,floor(2*pi*x_r/(2*d_r))));
nw_reco_coarse_M_ = nw_reco_M_;
if (n_w_coarse< n_w_max);
nw_reco_coarse_M_ = max(0,min(n_w_coarse-1,round(n_w_coarse*nw_reco_M_/n_w_max)));
end;%if (n_w_coarse< n_w_max);
na=0;
tmp_r = x_r; tmp_p = 0;
while (na<n_upb);
tmp_c = 2*pi*tmp_r;
tmp_a = floor(tmp_c/(2*d_r));
tmp_g_ = 2*pi*transpose([0:tmp_a-1])/tmp_a;
tmp_dg = 2*pi/max(1,tmp_a);
tmp_dg_ = tmp_dg*ones(tmp_a,1);
if ((tmp_p==1) & (tmp_a>1)); tmp_g_ = tmp_g_ + 0.5*tmp_dg; end;
tmp_r_ = tmp_r*ones(tmp_a,1);
tmp_dr = 2*d_r;
tmp_dr_ = tmp_dr*ones(tmp_a,1);
pos_g_(1+na+[0:tmp_a-1]) = tmp_g_;
pos_dg_(1+na+[0:tmp_a-1]) = tmp_dg_;
pos_r_(1+na+[0:tmp_a-1]) = tmp_r_;
pos_dr_(1+na+[0:tmp_a-1]) = tmp_dr_;
pos_i_(1+na+[0:tmp_a-1]) = max(0,min(n_w_coarse-1,round(n_w_coarse*tmp_g_/(2*pi))));
tmp_r = tmp_r + 2*d_r;
na = na+tmp_a;
tmp_p = 1-tmp_p;
end;%while (na<n_upb);
n_all = numel(pos_i_);
pos_x_0_ = pos_r_.*cos(pos_g_);
pos_x_1_ = pos_r_.*sin(pos_g_);

I__ = sparse(1+[0:n_all-1],1+pos_i_,1,n_all,n_w_coarse);
index_nall_from_nw__ = cell(n_w_coarse,1);
n_index_nall_from_nw_ = zeros(n_w_coarse,1);
for nw=0:n_w_coarse-1;
index_nall_from_nw__{1+nw} = sort(efind(I__(:,1+nw)));
n_index_nall_from_nw_(1+nw) = numel(index_nall_from_nw__{1+nw});
end;%for nw=0:n_w_coarse-1;

index_nall_from_nM_ = zeros(n_M,1);
index_tab_from_nw_ = zeros(n_w_coarse,1);
for nM=0:n_M-1;
nw = nw_reco_coarse_M_(1+nM);
index_nall_from_nw_ = index_nall_from_nw__{1+nw};
tmp_tab = index_tab_from_nw_(1+nw);
nall = index_nall_from_nw_(1+tmp_tab);
index_nall_from_nM_(1+nM) = nall;
index_tab_from_nw_(1+nw) = index_tab_from_nw_(1+nw) + 1;
end;%for nM=0:n_M-1;
index_tot_from_nw_ = index_tab_from_nw_;

flag_check=0;
if flag_check;
figure(1);clf;figmed;
subplot(1,3,1); plot(pos_r_(1+index_nall_from_nM_),'.');
subplot(1,3,2); plot(pos_r_,'.');
disp('returning');return;
end;%if flag_check;

pos_x_reco_0_pre_pre_ = x_0 + ( pos_r_(1+index_nall_from_nM_) + 0.0*pos_dr_(1+index_nall_from_nM_) ) .* cos( pos_g_(1+index_nall_from_nM_) - 0.5*pos_dg_(1+index_nall_from_nM_) );
pos_x_reco_0_pre_pos_ = x_0 + ( pos_r_(1+index_nall_from_nM_) + 0.0*pos_dr_(1+index_nall_from_nM_) ) .* cos( pos_g_(1+index_nall_from_nM_) + 0.5*pos_dg_(1+index_nall_from_nM_) );
pos_x_reco_0_pos_pre_ = x_0 + ( pos_r_(1+index_nall_from_nM_) + 1.0*pos_dr_(1+index_nall_from_nM_) ) .* cos( pos_g_(1+index_nall_from_nM_) - 0.5*pos_dg_(1+index_nall_from_nM_) );
pos_x_reco_0_pos_pos_ = x_0 + ( pos_r_(1+index_nall_from_nM_) + 1.0*pos_dr_(1+index_nall_from_nM_) ) .* cos( pos_g_(1+index_nall_from_nM_) + 0.5*pos_dg_(1+index_nall_from_nM_) );

pos_x_reco_1_pre_pre_ = x_1 + ( pos_r_(1+index_nall_from_nM_) + 0.0*pos_dr_(1+index_nall_from_nM_) ) .* sin( pos_g_(1+index_nall_from_nM_) - 0.5*pos_dg_(1+index_nall_from_nM_) );
pos_x_reco_1_pre_pos_ = x_1 + ( pos_r_(1+index_nall_from_nM_) + 0.0*pos_dr_(1+index_nall_from_nM_) ) .* sin( pos_g_(1+index_nall_from_nM_) + 0.5*pos_dg_(1+index_nall_from_nM_) );
pos_x_reco_1_pos_pre_ = x_1 + ( pos_r_(1+index_nall_from_nM_) + 1.0*pos_dr_(1+index_nall_from_nM_) ) .* sin( pos_g_(1+index_nall_from_nM_) - 0.5*pos_dg_(1+index_nall_from_nM_) );
pos_x_reco_1_pos_pos_ = x_1 + ( pos_r_(1+index_nall_from_nM_) + 1.0*pos_dr_(1+index_nall_from_nM_) ) .* sin( pos_g_(1+index_nall_from_nM_) + 0.5*pos_dg_(1+index_nall_from_nM_) );

pos_x_reco_0__ = ...
[ ...
 pos_x_reco_0_pre_pre_ ...
,pos_x_reco_0_pre_pos_ ...
,pos_x_reco_0_pos_pos_ ...
,pos_x_reco_0_pos_pre_ ...
];
pos_x_reco_0__ = transpose(pos_x_reco_0__);

pos_x_reco_1__ = ...
[ ...
 pos_x_reco_1_pre_pre_ ...
,pos_x_reco_1_pre_pos_ ...
,pos_x_reco_1_pos_pos_ ...
,pos_x_reco_1_pos_pre_ ...
];
pos_x_reco_1__ = transpose(pos_x_reco_1__);

pos_x_reco_2__ = +1e-3*ones(size(pos_x_reco_1__));
pos_c__ = zeros(1,n_M,3);
nc_ = max(0,min(n_c-1,floor(n_c*nw_true_M_/n_w_max)));
pos_c__(1,:,:) = c__(1+nc_,:);

p = patch(pos_x_reco_0__,pos_x_reco_1__,pos_x_reco_2__,pos_c__);
set(p,'EdgeColor','none');








