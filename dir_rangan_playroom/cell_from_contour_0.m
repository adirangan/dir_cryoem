function [n_item,val_i_,length_i_,x_il__,y_il__] = cell_from_contour_0(c__);
flag_verbose=0;
n_c = size(c__,2);
nitem=0; nc=0;
while (nc<n_c-1);
nc = nc + c__(2,1+nc)+1;
nitem = nitem+1;
end;%while (nc<n_c-1);
n_item = nitem;
if (flag_verbose); disp(sprintf(' %% n_item %d',n_item)); end;
val_i_ = zeros(n_item,1);
length_i_ = zeros(n_item,1);
x_il__ = cell(n_item,1);
y_il__ = cell(n_item,1);
nitem=0; nc=0;
while (nc<n_c-1);
val = c__(1,1+nc);
length = c__(2,1+nc);
x_l_ = c__(1,1+nc+[1:length]);
y_l_ = c__(2,1+nc+[1:length]);
val_i_(1+nitem) = val;
length_i_(1+nitem) = length;
x_il__{1+nitem} = x_l_(:);
y_il__{1+nitem} = y_l_(:);
nc = nc + length+1;
nitem = nitem+1;
end;%while (nc<n_c-1);
n_item = nitem;

flag_disp=0;
if flag_disp;
c_80s__ = colormap_80s; n_c_80s = size(c_80s__,1);
clim_ = [min(val_i_),max(val_i_)],;
hold on;
for nitem=0:n_item-1;
val = val_i_(1+nitem);
nc_80s = max(0,min(n_c_80s-1,floor(n_c_80s*(val-min(clim_))/diff(clim_))));
p = patch(x_il__{1+nitem},y_il__{1+nitem},c_80s__(1+nc_80s,:));
end;%for nitem=0:n_item-1;
hold off;
end;%if flag_disp;
