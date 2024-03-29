function clim = polarpatch_playpen_0(Rg,Wg,A,clim,offset_x,offset_y,radius,cra);
if nargin<4; clim = mean(A(:)) + 2.0*std(A(:))*[-1,+1]; end;
if nargin<6; offset_x = 0; offset_y = 0; end;
if nargin<7; radius = 1; end;
if nargin<8; cra = colormap('cool'); end;
colormap(cra);
ncra = size(cra,1);
[nrows,ncols] = size(A);
r_pre = [1:nrows]; r_pos = [nrows , 1:nrows-1]; c_pre = [1:ncols-1]; c_pos = [2:ncols]; 
x_pre_pre = Rg(r_pre,c_pre).*cos(Wg(r_pre,c_pre));
x_pos_pre = Rg(r_pos,c_pre).*cos(Wg(r_pos,c_pre));
x_pre_pos = Rg(r_pre,c_pos).*cos(Wg(r_pre,c_pos));
x_pos_pos = Rg(r_pos,c_pos).*cos(Wg(r_pos,c_pos));
y_pre_pre = Rg(r_pre,c_pre).*sin(Wg(r_pre,c_pre));
y_pos_pre = Rg(r_pos,c_pre).*sin(Wg(r_pos,c_pre));
y_pre_pos = Rg(r_pre,c_pos).*sin(Wg(r_pre,c_pos));
y_pos_pos = Rg(r_pos,c_pos).*sin(Wg(r_pos,c_pos));
n_pre_pre = max(1,min(ncra,floor(ncra*(A(r_pre,c_pre)-min(clim))/diff(clim))));
n_pos_pre = max(1,min(ncra,floor(ncra*(A(r_pos,c_pre)-min(clim))/diff(clim))));
n_pre_pos = max(1,min(ncra,floor(ncra*(A(r_pre,c_pos)-min(clim))/diff(clim))));
n_pos_pos = max(1,min(ncra,floor(ncra*(A(r_pos,c_pos)-min(clim))/diff(clim))));
c_pre_pre_1 = cra(n_pre_pre(:),1); c_pre_pre_2 = cra(n_pre_pre(:),2); c_pre_pre_3 = cra(n_pre_pre(:),3);
c_pos_pre_1 = cra(n_pos_pre(:),1); c_pos_pre_2 = cra(n_pos_pre(:),2); c_pos_pre_3 = cra(n_pos_pre(:),3);
c_pre_pos_1 = cra(n_pre_pos(:),1); c_pre_pos_2 = cra(n_pre_pos(:),2); c_pre_pos_3 = cra(n_pre_pos(:),3);
c_pos_pos_1 = cra(n_pos_pos(:),1); c_pos_pos_2 = cra(n_pos_pos(:),2); c_pos_pos_3 = cra(n_pos_pos(:),3);
%c_mid_mid_1 = mean([cra(n_pre_pre(:),1),cra(n_pos_pre(:),1),cra(n_pre_pos(:),1),cra(n_pos_pos(:),1)],2);
%c_mid_mid_2 = mean([cra(n_pre_pre(:),2),cra(n_pos_pre(:),2),cra(n_pre_pos(:),2),cra(n_pos_pos(:),2)],2);
%c_mid_mid_3 = mean([cra(n_pre_pre(:),3),cra(n_pos_pre(:),3),cra(n_pre_pos(:),3),cra(n_pos_pos(:),3)],2);
x = transpose([x_pre_pre(:) , x_pre_pos(:) , x_pos_pos(:) , x_pos_pre(:)]);
y = transpose([y_pre_pre(:) , y_pre_pos(:) , y_pos_pos(:) , y_pos_pre(:)]);
n = transpose([n_pre_pre(:) , n_pre_pos(:) , n_pos_pos(:) , n_pos_pre(:)]);
c(:,:,1) = transpose([c_pre_pre_1(:) , c_pre_pos_1(:) , c_pos_pos_1(:) , c_pos_pre_1(:)]);
c(:,:,2) = transpose([c_pre_pre_2(:) , c_pre_pos_2(:) , c_pos_pos_2(:) , c_pos_pre_2(:)]);
c(:,:,3) = transpose([c_pre_pre_3(:) , c_pre_pos_3(:) , c_pos_pos_3(:) , c_pos_pre_3(:)]);
%c(:,:,1) = transpose([c_mid_mid_1(:) , c_mid_mid_1(:) , c_mid_mid_1(:) , c_mid_mid_1(:)]);
%c(:,:,2) = transpose([c_mid_mid_2(:) , c_mid_mid_2(:) , c_mid_mid_2(:) , c_mid_mid_2(:)]);
%c(:,:,3) = transpose([c_mid_mid_3(:) , c_mid_mid_3(:) , c_mid_mid_3(:) , c_mid_mid_3(:)]);
p=patch(x*radius+offset_x,y*radius+offset_y,c);set(p,'EdgeColor','none'); 
if offset_x==0 & offset_y==0 & radius==1; xlim(max(Rg(:))*[-1,+1]); ylim(max(Rg(:))*[-1,+1]); axis off; end;

