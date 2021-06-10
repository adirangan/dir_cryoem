function isosurface_f_x_c_0(f_x_c___,percent_threshold);

[n_x_c,n_y_c,n_z_c] = size(f_x_c___);
[X___,Y___,Z___] = meshgrid(1:n_x_c,1:n_y_c,1:n_z_c);

c_ = colormap('spring'); c_ = c_(end:-1:1,:); 
%c_ = colormap('cool'); 
n_c = size(c_,1);
v_avg = mean(f_x_c___(:)); v_std = std(f_x_c___(:)); 
v_min = min(f_x_c___(:)); v_max = max(f_x_c___(:));
vlim = v_avg + 1.5*v_std*[-1,1];
%vlim = [v_min , v_max];
if nargin<2; v_ = linspace(vlim(1),vlim(2),3); end;
if nargin>=2; v_ = prctile(f_x_c___(:),percent_threshold); end;
for nv=length(v_):-1:1
v = v_(nv);
hpatch = patch(isosurface(X___,Y___,Z___,f_x_c___,v)); 
isonormals(X___,Y___,Z___,f_x_c___,hpatch);
nc = max(1,min(n_c,floor(n_c*nv/length(v_))));
%hpatch.FaceColor = 'red'; 
hpatch.FaceColor = c_(nc,:);
hpatch.EdgeColor = 'none';
%hpatch.EdgeColor = c_(nc,:);
%hpatch.FaceAlpha = (nv/length(v_)).^2 * 0.5;
hpatch.FaceAlpha = (nv/length(v_)).^4 * 1.0;
%hpatch.EdgeAlpha = (nv/length(v_)).^2 * 1;
%title(sprintf('v %0.2f',v));
xlim([1,n_x_c]);ylim([1,n_y_c]);zlim([1,n_z_c]);
xlabel('x');ylabel('y');zlabel('z');
set(gca,'XTick',[],'YTick',[],'ZTick',[]);
%daspect([1,4,4]); 
view([-65,20]); 
axis vis3d;
end;%for nv=1:length(v_);
if (length(v_)<=1);camlight left; lighting gouraud; end;
figbig;
