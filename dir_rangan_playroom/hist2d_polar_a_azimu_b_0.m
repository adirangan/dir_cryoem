function ...
[ ...
 hist_ ...
 weighted_hist_ ...
 lim_use_ ...
] = ...
hist2d_polar_a_azimu_b_0( ...
 polar_a_ ...
,azimu_b_ ...
,weight_ ...
,data_polar_a_ ...
,data_azimu_b_ ...
,lim_ ...
,c__ ...
,flag_2d_vs_3d ...
,flag_loghist_vs_hist ...
);
% We assume that there are several values of azimu_b_ for each value of polar_a_. ;
str_thisfunction = 'hist2d_polar_a_azimu_b_0';

na=0;
if (nargin<1+na); polar_a_=[]; end; na=na+1;
if (nargin<1+na); azimu_b_=[]; end; na=na+1;
if (nargin<1+na); weight_=[]; end; na=na+1;
if (nargin<1+na); data_polar_a_=[]; end; na=na+1;
if (nargin<1+na); data_azimu_b_=[]; end; na=na+1;
if (nargin<1+na); lim_=[]; end; na=na+1;
if (nargin<1+na); c__=[]; end; na=na+1;
if (nargin<1+na); flag_2d_vs_3d=[]; end; na=na+1;
if (nargin<1+na); flag_loghist_vs_hist=[]; end; na=na+1;

lim_use_ = lim_;
if isempty(lim_); lim_use_ = [0,2*numel(data_polar_a_)/max(1,numel(polar_a_))]; end;
if isempty(c__);  c__ = colormap_beach(); end;
if isempty(flag_2d_vs_3d); flag_2d_vs_3d = 1; end;
if isempty(flag_loghist_vs_hist); flag_loghist_vs_hist = 0; end;

n_c = size(c__,1);

assert(numel(polar_a_)==numel(azimu_b_));
hist_ = zeros(size(polar_a_));
weighted_hist_ = zeros(size(polar_a_));
data_polar_a_ = periodize(data_polar_a_,0,1*pi);
data_azimu_b_ = periodize(data_azimu_b_,0,2*pi);

u_polar_a_ = unique(polar_a_); n_polar_a = numel(u_polar_a_);
index_polar_a_ = cell(n_polar_a,1);
u_azimu_b__ = cell(n_polar_a,1);
n_azimu_b_ = zeros(n_polar_a,1);
index_azimu_b__ = cell(n_polar_a,1);
for npolar_a=0:n_polar_a-1;
polar_a = u_polar_a_(1+npolar_a);
index_polar_a_{1+npolar_a} = find(polar_a_==polar_a) - 1;
u_azimu_b__{1+npolar_a} = unique(azimu_b_(1+index_polar_a_{1+npolar_a}));
n_azimu_b_(1+npolar_a) = numel(u_azimu_b__{1+npolar_a});
index_azimu_b__{1+npolar_a} = zeros(n_azimu_b_(1+npolar_a),1);
for nazimu_b=0:n_azimu_b_(1+npolar_a)-1;
azimu_b = u_azimu_b__{1+npolar_a}(1+nazimu_b);
index_azimu_b__{1+npolar_a}(1+nazimu_b) = index_polar_a_{1+npolar_a}(1+find(azimu_b_(1+index_polar_a_{1+npolar_a})==azimu_b)-1);
end;%for nazimu_b=0:n_azimu_b_(1+npolar_a)-1;
end;%for npolar_a=0:n_polar_a-1;

n_patch = sum(n_azimu_b_) + n_polar_a;
weighted_hist_by_patch_ = zeros(1,n_patch);
x__ = zeros(5,n_patch);
y__ = zeros(5,n_patch);
if flag_2d_vs_3d==0; z__ = zeros(5,n_patch); end;
c___ = zeros(1,n_patch,3);
%%%%%%%%;
[s_polar_a_,index_u_from_s_polar_a_] = sort(u_polar_a_); index_u_from_s_polar_a_ = index_u_from_s_polar_a_ - 1;
[~,index_s_from_u_polar_a_] = sort(index_u_from_s_polar_a_); index_s_from_u_polar_a_ = index_s_from_u_polar_a_ - 1;
%%%%%%%%%%%%%%%%;
npatch=0; x_lo_min = 0; x_hi_max = 0;
%%%%%%%%%%%%%%%%;
for npolar_a=0:n_polar_a-1;
polar_a = u_polar_a_(1+npolar_a);
y_up = 1*pi; if (index_s_from_u_polar_a_(1+npolar_a)<n_polar_a-1); y_up = 0.5*(polar_a + s_polar_a_(1+index_s_from_u_polar_a_(1+npolar_a)+1)); end;
y_dn = 0.00; if (index_s_from_u_polar_a_(1+npolar_a)>          0); y_dn = 0.5*(polar_a + s_polar_a_(1+index_s_from_u_polar_a_(1+npolar_a)-1)); end;
u_azimu_b_ = u_azimu_b__{1+npolar_a};
[s_azimu_b_,index_u_from_s_azimu_b_] = sort(u_azimu_b_); index_u_from_s_azimu_b_ = index_u_from_s_azimu_b_ - 1;
[~,index_s_from_u_azimu_b_] = sort(index_u_from_s_azimu_b_); index_s_from_u_azimu_b_ = index_s_from_u_azimu_b_ - 1;
n_azimu_b = n_azimu_b_(1+npolar_a);
%%%%%%%%;
index_data_a_ = efind( (data_polar_a_<=y_up) & (data_polar_a_> y_dn) );
data_azimu_b_sub_ = data_azimu_b_(1+index_data_a_);
x_hi = 0.5*(2*pi + max(u_azimu_b_)); x_hi_max = max(x_hi_max,x_hi);
x_lo = x_hi - 2*pi; x_lo_min = min(x_lo_min,x_lo);
for nazimu_b=0:n_azimu_b-1;
azimu_b = u_azimu_b_(1+nazimu_b);
x_up = x_hi; if (index_s_from_u_azimu_b_(1+nazimu_b)< n_azimu_b-1); x_up = 0.5*(azimu_b + s_azimu_b_(1+index_s_from_u_azimu_b_(1+nazimu_b)+1)); end;
x_dn = x_lo; if (index_s_from_u_azimu_b_(1+nazimu_b)>           0); x_dn = 0.5*(azimu_b + s_azimu_b_(1+index_s_from_u_azimu_b_(1+nazimu_b)-1)); end;
index_data_b_ = efind( (data_azimu_b_sub_<=x_up) & (data_azimu_b_sub_> x_dn) );
hist_val = numel(index_data_b_);
weight_val = weight_(1+index_azimu_b__{1+npolar_a}(1+nazimu_b));
hist_(1+index_azimu_b__{1+npolar_a}(1+nazimu_b)) = hist_val;
weighted_hist_val = hist_val/max(1e-12,weight_val);
weighted_hist_val_use = weighted_hist_val; if flag_loghist_vs_hist; weighted_hist_val_use = log2(1+weighted_hist_val); end;
weighted_hist_(1+index_azimu_b__{1+npolar_a}(1+nazimu_b)) = weighted_hist_val;
nc = max(0,min(n_c-1,floor((weighted_hist_val_use-min(lim_use_))/max(1e-12,diff(lim_use_))*n_c)));
%%%%%%%%;
if flag_2d_vs_3d==1;
if (x_dn>=0*pi);
x__(1+0,1+npatch) = x_dn;
x__(1+1,1+npatch) = x_up;
x__(1+2,1+npatch) = x_up;
x__(1+3,1+npatch) = x_dn;
x__(1+4,1+npatch) = x_dn;
y__(1+0,1+npatch) = y_dn;
y__(1+1,1+npatch) = y_dn;
y__(1+2,1+npatch) = y_up;
y__(1+3,1+npatch) = y_up;
y__(1+4,1+npatch) = y_dn;
c___(1,1+npatch,:) = c__(1+nc,:);
weighted_hist_by_patch_(1+npatch) = weighted_hist_val_use;
npatch=npatch+1;
end;%if (x_dn>=0*pi);
if (x_dn< 0*pi);
x_0p = 0*pi; x_2p = 2*pi;
x__(1+0,1+npatch) = x_0p;
x__(1+1,1+npatch) = x_up;
x__(1+2,1+npatch) = x_up;
x__(1+3,1+npatch) = x_0p;
x__(1+4,1+npatch) = x_0p;
y__(1+0,1+npatch) = y_dn;
y__(1+1,1+npatch) = y_dn;
y__(1+2,1+npatch) = y_up;
y__(1+3,1+npatch) = y_up;
y__(1+4,1+npatch) = y_dn;
c___(1,1+npatch,:) = c__(1+nc,:);
weighted_hist_by_patch_(1+npatch) = weighted_hist_val_use;
npatch=npatch+1;
x__(1+0,1+npatch) = 2*pi+x_dn;
x__(1+1,1+npatch) = x_2p;
x__(1+2,1+npatch) = x_2p;
x__(1+3,1+npatch) = 2*pi+x_dn;
x__(1+4,1+npatch) = 2*pi+x_dn;
y__(1+0,1+npatch) = y_dn;
y__(1+1,1+npatch) = y_dn;
y__(1+2,1+npatch) = y_up;
y__(1+3,1+npatch) = y_up;
y__(1+4,1+npatch) = y_dn;
c___(1,1+npatch,:) = c__(1+nc,:);
weighted_hist_by_patch_(1+npatch) = weighted_hist_val_use;
npatch=npatch+1;
end;%if (x_dn< 0*pi);
end;%if flag_2d_vs_3d==1;
%%%%%%%%;
if flag_2d_vs_3d==0;
k_3d_0_sw = cos(x_dn)*sin(y_dn); 
k_3d_0_se = cos(x_up)*sin(y_dn);
k_3d_0_ne = cos(x_up)*sin(y_up); 
k_3d_0_nw = cos(x_dn)*sin(y_up);
k_3d_1_sw = sin(x_dn)*sin(y_dn); 
k_3d_1_se = sin(x_up)*sin(y_dn);
k_3d_1_ne = sin(x_up)*sin(y_up); 
k_3d_1_nw = sin(x_dn)*sin(y_up);
k_3d_2_sw = cos(y_dn); 
k_3d_2_se = cos(y_dn);
k_3d_2_ne = cos(y_up); 
k_3d_2_nw = cos(y_up);
x__(1+0,1+npatch) = k_3d_0_sw;
x__(1+1,1+npatch) = k_3d_0_se;
x__(1+2,1+npatch) = k_3d_0_ne;
x__(1+3,1+npatch) = k_3d_0_nw;
x__(1+4,1+npatch) = k_3d_0_sw;
y__(1+0,1+npatch) = k_3d_1_sw;
y__(1+1,1+npatch) = k_3d_1_se;
y__(1+2,1+npatch) = k_3d_1_ne;
y__(1+3,1+npatch) = k_3d_1_nw;
y__(1+4,1+npatch) = k_3d_1_sw;
z__(1+0,1+npatch) = k_3d_2_sw;
z__(1+1,1+npatch) = k_3d_2_se;
z__(1+2,1+npatch) = k_3d_2_ne;
z__(1+3,1+npatch) = k_3d_2_nw;
z__(1+4,1+npatch) = k_3d_2_sw;
c___(1,1+npatch,:) = c__(1+nc,:);
weighted_hist_by_patch_(1+npatch) = weighted_hist_val_use;
npatch=npatch+1;
end;%if flag_2d_vs_3d==0;
%%%%%%%%;
end;%for nazimu_b=0:n_azimu_b-1;
end;%for npolar_a=0:n_polar_a-1;
%%%%%%%%%%%%%%%%;

n_patch = npatch; npatch=0;
weighted_hist_by_patch_ = weighted_hist_by_patch_(1:n_patch);
x__ = x__(:,1:n_patch);
y__ = y__(:,1:n_patch);
if flag_2d_vs_3d==0; z__ = z__(:,1:n_patch); end;
c___ = c___(:,1:n_patch,:);

lim_use_ = lim_;
if isempty(lim_use_);
lim_use_ = [0,prctile(weighted_hist_by_patch_,100)];
for npatch=0:n_patch-1;
weighted_hist_val_use = weighted_hist_by_patch_(1+npatch);
nc = max(0,min(n_c-1,floor((weighted_hist_val_use-min(lim_use_))/max(1e-12,diff(lim_use_))*n_c)));
c___(1,1+npatch,:) = c__(1+nc,:);
end;%for npatch=0:n_patch-1;
end;%if isempty(lim_use_);

if flag_2d_vs_3d==0 | flag_2d_vs_3d==1;
colormap(c__);
%%%%%%%%;
if flag_2d_vs_3d==1;
p=patch(x__,y__,c___,'EdgeColor','none');
xlim([0,2*pi]);%xlim([x_lo_min,x_hi_max]);
ylim([0,1*pi]);
axisnotick;
end;%if flag_2d_vs_3d==1;
%%%%%%%%;
if flag_2d_vs_3d==0;
p=patch(x__,y__,z__,c___,'EdgeColor','none');
view([-65,20]);
axis vis3d;
end;%if flag_2d_vs_3d==0;
%%%%%%%%;
tmp_c_ = colorbar;
set(tmp_c_,'Ticks',[0,1],'TickLabels',lim_use_);
end;%if flag_2d_vs_3d==0 | flag_2d_vs_3d==1;



