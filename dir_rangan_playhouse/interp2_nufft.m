function output = interp2_nufft(n_x,min_x,max_x,n_y,min_y,max_y,F__,x_loc_,y_loc_);
% assumes regular grid for x and y ;
% periodic boundary conditions ;

assert(numel(x_loc_)==numel(y_loc_));

x_loc_ = periodize(x_loc_,min_x,max_x); %<-- wrap x_loc_ into [x_min,x_max). ;
x_loc_ = (x_loc_ - min_x)/(max_x-min_x); %<-- map x_loc_ to [0,1) interval. ;
tmp_ij_ = find(x_loc_>=0.5); x_loc_(tmp_ij_) = x_loc_(tmp_ij_)-1.0; %<-- wrap x_loc_ into [-0.5,0.5) interval. ;
x_loc_ = (2*pi/n_x) * x_loc_*n_x;

y_loc_ = periodize(y_loc_,min_y,max_y); %<-- wrap y_loc_ into [y_min,y_max). ;
y_loc_ = (y_loc_ - min_y)/(max_y-min_y); %<-- map y_loc_ to [0,1) interval. ;
tmp_ij_ = find(y_loc_>=0.5); y_loc_(tmp_ij_) = y_loc_(tmp_ij_)-1.0; %<-- wrap y_loc_ into [-0.5,0.5) interval. ;
y_loc_ = (2*pi/n_y) * y_loc_*n_y;

G__ = fft2(F__)/sqrt(n_x*n_y); %<-- frequency range 2*pi*[0,n_x-1] and 2*pi*[0,n_y-1]. ;
ij_x_ = [ [ceil(n_x/2)+1:n_x] , [1:ceil(n_x/2)] ];
ij_y_ = [ [ceil(n_y/2)+1:n_y] , [1:ceil(n_y/2)] ];
G__ = G__(ij_x_,ij_y_); %<-- wrap G_ to run from frequency range pi*[-n_x/2,+n_x/2-1] and pi*[-n_y/2,+n_y/2-1]. ;
H__ = nufft2d2(numel(x_loc_),x_loc_,y_loc_,+1,1e-12,n_x,n_y,G__(:));
output = H__/sqrt(n_x*n_y);
