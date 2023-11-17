%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Now load the actual Micrograph. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
flag_invert = 1;
flag_center_volume = 0;
fname_prefix='clathrin_montage_Micrograph';
fname_infix='clathrin';
Pixel_Spacing=1.06;
fname_nopath_volume='6sct_one_leg_bf60_center.mrc';
fname_nopath_micrograph='simulated_clathrin_t1000A_d5000A_1.06Apix_p1_384.mrc';

%%%%%%%%;
fname_mat = sprintf('%s_mat/O_x_u_pack_6_.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
fname_mrc = sprintf('%s/%s',dir_data_star,fname_nopath_micrograph);
O_x_u_load_ = cast(ReadMRC(fname_mrc),'double');
n_O_x_u_ = size(O_x_u_load_);
n_pack_6 = n_O_x_u_(1+0)/(10*64);
assert(n_pack_6==6);
n_O_x_u_pack_ = n_O_x_u_/n_pack_6;
pack_0_row_ij_ = zeros(n_O_x_u_pack_(1+0),1);
pack_0_col_ij_ = zeros(n_O_x_u_pack_(1+0),1);
pack_0_val_ij_ = zeros(n_O_x_u_pack_(1+0),1);
na=0;
for nO_x_u=0:n_O_x_u_(1+0)-1;
pack_0_row_ij_(1+na) = 1+nO_x_u;
pack_0_col_ij_(1+na) = 1+floor(nO_x_u/n_pack_6);
pack_0_val_ij_(1+na) = 1/n_pack_6;
na=na+1;
end;%for nO_x_u=0:n_O_x_u_(1+0)-1;
pack_1_row_ij_ = zeros(n_O_x_u_pack_(1+1),1);
pack_1_col_ij_ = zeros(n_O_x_u_pack_(1+1),1);
pack_1_val_ij_ = zeros(n_O_x_u_pack_(1+1),1);
na=0;
for nO_x_u=0:n_O_x_u_(1+1)-1;
pack_1_row_ij_(1+na) = 1+nO_x_u;
pack_1_col_ij_(1+na) = 1+floor(nO_x_u/n_pack_6);
pack_1_val_ij_(1+na) = 1/n_pack_6;
na=na+1;
end;%for nO_x_u=0:n_O_x_u_(1+1)-1;
x_u_pack_0_ = sparse(pack_0_row_ij_,pack_0_col_ij_,pack_0_val_ij_,n_O_x_u_(1+0),n_O_x_u_pack_(1+0));
x_u_pack_1_ = sparse(pack_1_row_ij_,pack_1_col_ij_,pack_1_val_ij_,n_O_x_u_(1+1),n_O_x_u_pack_(1+1));
O_x_u_pack_ = transpose(x_u_pack_0_)*O_x_u_load_*x_u_pack_1_;
if flag_invert; O_x_u_pack_ = -O_x_u_pack_; end;
n_M = numel(O_x_u_pack_);
%%%%%%%%;
% Now subtract off a parabolic fit. ;
%%%%%%%%;
tmp_j0_ = transpose(1:n_O_x_u_pack_(1+0));
tmp_j1_ = transpose(1:n_O_x_u_pack_(1+1));
[tmp_j0__,tmp_j1__] = ndgrid(tmp_j0_,tmp_j1_);
tmp_RHS_00 = sum( O_x_u_pack_ .* (tmp_j0__.^0) .* (tmp_j1__.^0) , 'all' );
tmp_RHS_01 = sum( O_x_u_pack_ .* (tmp_j0__.^0) .* (tmp_j1__.^1) , 'all' );
tmp_RHS_10 = sum( O_x_u_pack_ .* (tmp_j0__.^1) .* (tmp_j1__.^0) , 'all' );
tmp_RHS_11 = sum( O_x_u_pack_ .* (tmp_j0__.^1) .* (tmp_j1__.^1) , 'all' );
tmp_RHS_02 = sum( O_x_u_pack_ .* (tmp_j0__.^0) .* (tmp_j1__.^2) , 'all' );
tmp_RHS_20 = sum( O_x_u_pack_ .* (tmp_j0__.^2) .* (tmp_j1__.^0) , 'all' );
tmp_RHS_ = [ ...
 tmp_RHS_00 ...
;tmp_RHS_01 ...
;tmp_RHS_10 ...
;tmp_RHS_11 ...
;tmp_RHS_02 ...
;tmp_RHS_20 ...
];
tmp_LHS_00 = sum( (tmp_j0__.^0) .* (tmp_j1__.^0) , 'all' );
tmp_LHS_01 = sum( (tmp_j0__.^0) .* (tmp_j1__.^1) , 'all' );
tmp_LHS_10 = sum( (tmp_j0__.^1) .* (tmp_j1__.^0) , 'all' );
tmp_LHS_11 = sum( (tmp_j0__.^1) .* (tmp_j1__.^1) , 'all' );
tmp_LHS_02 = sum( (tmp_j0__.^0) .* (tmp_j1__.^2) , 'all' );
tmp_LHS_20 = sum( (tmp_j0__.^2) .* (tmp_j1__.^0) , 'all' );
tmp_LHS_12 = sum( (tmp_j0__.^1) .* (tmp_j1__.^2) , 'all' );
tmp_LHS_21 = sum( (tmp_j0__.^2) .* (tmp_j1__.^1) , 'all' );
tmp_LHS_03 = sum( (tmp_j0__.^0) .* (tmp_j1__.^3) , 'all' );
tmp_LHS_30 = sum( (tmp_j0__.^3) .* (tmp_j1__.^0) , 'all' );
tmp_LHS_20 = sum( (tmp_j0__.^2) .* (tmp_j1__.^0) , 'all' );
tmp_LHS_22 = sum( (tmp_j0__.^2) .* (tmp_j1__.^2) , 'all' );
tmp_LHS_13 = sum( (tmp_j0__.^1) .* (tmp_j1__.^3) , 'all' );
tmp_LHS_31 = sum( (tmp_j0__.^3) .* (tmp_j1__.^1) , 'all' );
tmp_LHS_04 = sum( (tmp_j0__.^0) .* (tmp_j1__.^4) , 'all' );
tmp_LHS_40 = sum( (tmp_j0__.^4) .* (tmp_j1__.^0) , 'all' );
tmp_LHS__ = [ ...
  tmp_LHS_00 , tmp_LHS_01 , tmp_LHS_10 , tmp_LHS_11 , tmp_LHS_02 , tmp_LHS_20 ...
; tmp_LHS_01 , tmp_LHS_02 , tmp_LHS_11 , tmp_LHS_12 , tmp_LHS_03 , tmp_LHS_21 ...
; tmp_LHS_10 , tmp_LHS_11 , tmp_LHS_20 , tmp_LHS_21 , tmp_LHS_12 , tmp_LHS_30 ...
; tmp_LHS_11 , tmp_LHS_12 , tmp_LHS_21 , tmp_LHS_22 , tmp_LHS_13 , tmp_LHS_31 ...
; tmp_LHS_02 , tmp_LHS_03 , tmp_LHS_12 , tmp_LHS_13 , tmp_LHS_04 , tmp_LHS_22 ...
; tmp_LHS_20 , tmp_LHS_21 , tmp_LHS_30 , tmp_LHS_31 , tmp_LHS_22 , tmp_LHS_40 ...
];
tmp_axx_ = tmp_LHS__ \ tmp_RHS_;
a00 = tmp_axx_(1+0);
a01 = tmp_axx_(1+1);
a10 = tmp_axx_(1+2);
a11 = tmp_axx_(1+3);
a02 = tmp_axx_(1+4);
a20 = tmp_axx_(1+5);
tmp_p = @(j0,j1) ...
  a00.*(j0.^0).*(j1.^0) ...
+ a01.*(j0.^0).*(j1.^1) ...
+ a10.*(j0.^1).*(j1.^0) ...
+ a11.*(j0.^1).*(j1.^1) ...
+ a02.*(j0.^0).*(j1.^2) ...
+ a20.*(j0.^2).*(j1.^0) ...
;
p_x_u_pack_ = tmp_p(tmp_j0__,tmp_j1__);
P_x_u_pack_ = O_x_u_pack_ - p_x_u_pack_;
%%%%;
save(fname_mat ...
,'n_O_x_u_','n_O_x_u_pack_','n_pack_6','O_x_u_pack_','n_M' ...
,'p_x_u_pack_','P_x_u_pack_' ...
);
%%%%%%%%;
end;%if (flag_recalc | ~exist(fname_mat,'file'));
%%%%%%%%;
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/O_x_u_pack_FIGA',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figmed;figbeach();
subplot(1,3,1); imagesc(O_x_u_pack_); axis image; axisnotick; title('O_x_u_pack_','Interpreter','none');
subplot(1,3,2); imagesc(p_x_u_pack_); axis image; axisnotick; title('p_x_u_pack_','Interpreter','none');
subplot(1,3,3); imagesc(P_x_u_pack_); axis image; axisnotick; title('O_x_u_pack_-p_x_u_pack_','Interpreter','none');
sgtitle(fname_fig,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
