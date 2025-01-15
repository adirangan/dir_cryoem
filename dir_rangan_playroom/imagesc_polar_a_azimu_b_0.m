function ...
imagesc_polar_a_azimu_b_0( ...
 polar_a_ ... 
,azimu_b_ ... 
,data_ ... 
,lim_ ... 
,c__ ... 
,flag_2d_vs_3d ...
,k_p_r_max ...
);
% We assume that there are several values of azimu_b_ for each value of polar_a_. ;
na=0;
if (nargin<1+na); polar_a_=[]; end; na=na+1;
if (nargin<1+na); azimu_b_=[]; end; na=na+1;
if (nargin<1+na); data_=[]; end; na=na+1;
if (nargin<1+na); lim_=[]; end; na=na+1;
if (nargin<1+na); c__=[]; end; na=na+1;
if (nargin<1+na); flag_2d_vs_3d=[]; end; na=na+1;
if (nargin<1+na); k_p_r_max=[]; end; na=na+1;

if isempty(lim_); lim_ = mean(data_) + 1.5*std(data_,1)*[-1,1]; end;
if isempty(c__);  c__ = colormap_beach(); end;
if isempty(flag_2d_vs_3d); flag_2d_vs_3d = 1; end;
if isempty(k_p_r_max); k_p_r_max = 1.0; end;

n_c = size(c__,1);

assert(numel(polar_a_)==numel(azimu_b_));
assert(numel(polar_a_)==numel(data_));

u_polar_a_ = unique(polar_a_); n_polar_a = numel(u_polar_a_);
index_polar_a_ = cell(n_polar_a,1);
u_azimu_b__ = cell(n_polar_a,1);
n_azimu_b_ = zeros(n_polar_a,1);
index_azimu_b__ = cell(n_polar_a,1);
for npolar_a=0:n_polar_a-1;
polar_a = u_polar_a_(1+npolar_a);
index_polar_a_{1+npolar_a} = efind(polar_a_==polar_a);
u_azimu_b__{1+npolar_a} = unique(azimu_b_(1+index_polar_a_{1+npolar_a}));
n_azimu_b_(1+npolar_a) = numel(u_azimu_b__{1+npolar_a});
index_azimu_b__{1+npolar_a} = zeros(n_azimu_b_(1+npolar_a),1);
for nazimu_b=0:n_azimu_b_(1+npolar_a)-1;
azimu_b = u_azimu_b__{1+npolar_a}(1+nazimu_b);
index_azimu_b__{1+npolar_a}(1+nazimu_b) = min(index_polar_a_{1+npolar_a}(1+efind(azimu_b_(1+index_polar_a_{1+npolar_a})==azimu_b)));
end;%for nazimu_b=0:n_azimu_b_(1+npolar_a)-1;
end;%for npolar_a=0:n_polar_a-1;

n_patch = sum(n_azimu_b_);
x__ = zeros(5,n_patch);
y__ = zeros(5,n_patch);
if flag_2d_vs_3d==0; z__ = zeros(5,n_patch); end;
c___ = zeros(1,n_patch,3);
%%%%%%%%;
[s_polar_a_,index_u_from_s_polar_a_] = sort(u_polar_a_); index_u_from_s_polar_a_ = index_u_from_s_polar_a_ - 1;
[~,index_s_from_u_polar_a_] = sort(index_u_from_s_polar_a_); index_s_from_u_polar_a_ = index_s_from_u_polar_a_ - 1;
npatch=0;
for npolar_a=0:n_polar_a-1;
polar_a = u_polar_a_(1+npolar_a);
y_up = 1*pi; if (index_s_from_u_polar_a_(1+npolar_a)<n_polar_a-1); y_up = 0.5*(polar_a + s_polar_a_(1+index_s_from_u_polar_a_(1+npolar_a)+1)); end;
y_dn = 0.00; if (index_s_from_u_polar_a_(1+npolar_a)>          0); y_dn = 0.5*(polar_a + s_polar_a_(1+index_s_from_u_polar_a_(1+npolar_a)-1)); end;
u_azimu_b_ = u_azimu_b__{1+npolar_a};
n_u_azimu_b = numel(u_azimu_b_);
n_azimu_b = n_azimu_b_(1+npolar_a);
data_b_ = data_(1+index_azimu_b__{1+npolar_a}(1+[0:n_azimu_b-1]));
if (n_u_azimu_b==1);
if npolar_a==          0; npolar_a_use = min(n_polar_a-1,npolar_a+1); end;
if npolar_a==n_polar_a-1; npolar_a_use = max(          0,npolar_a-1); end;
u_azimu_b_ = u_azimu_b__{1+npolar_a_use};
n_u_azimu_b = numel(u_azimu_b_);
n_azimu_b = n_u_azimu_b;
data_b_ = data_(1+index_azimu_b__{1+npolar_a}(1+0))*ones(n_azimu_b,1);
end;%if (n_u_azimu_b==1);
if (n_u_azimu_b> 1);
% do nothing. ;
end;%if (n_u_azimu_b> 1);
[s_azimu_b_,index_u_from_s_azimu_b_] = sort(u_azimu_b_); index_u_from_s_azimu_b_ = index_u_from_s_azimu_b_ - 1;
[~,index_s_from_u_azimu_b_] = sort(index_u_from_s_azimu_b_); index_s_from_u_azimu_b_ = index_s_from_u_azimu_b_ - 1;
for nazimu_b=0:n_azimu_b-1;
azimu_b = u_azimu_b_(1+nazimu_b);
x_up = 2*pi; if (index_s_from_u_azimu_b_(1+nazimu_b)<n_azimu_b-1); x_up = 0.5*(azimu_b + s_azimu_b_(1+index_s_from_u_azimu_b_(1+nazimu_b)+1)); end;
x_dn = 0.00; if (index_s_from_u_azimu_b_(1+nazimu_b)>          0); x_dn = 0.5*(azimu_b + s_azimu_b_(1+index_s_from_u_azimu_b_(1+nazimu_b)-1)); end;
data = data_b_(1+nazimu_b);
nc = max(0,min(n_c-1,floor((data-min(lim_))/diff(lim_)*n_c)));
%%%%%%%%;
if flag_2d_vs_3d==1;
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
end;%if flag_2d_vs_3d==0;
%%%%%%%%;

npatch=npatch+1;
end;%for nazimu_b=0:n_azimu_b-1;
end;%for npolar_a=0:n_polar_a-1;
%%%%%%%%;
if flag_2d_vs_3d==1; p=patch(k_p_r_max*x__,k_p_r_max*y__,c___,'EdgeColor','none'); xlim([0,2*pi]);ylim([0,pi]); axisnotick; end;
if flag_2d_vs_3d==0; p=patch(k_p_r_max*x__,k_p_r_max*y__,k_p_r_max*z__,c___,'EdgeColor','none'); view([-65,20]); axis vis3d; end;
%%%%%%%%;


