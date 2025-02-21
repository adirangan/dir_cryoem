function area_t_ = heron_0(p_23t___);
%%%%%%%%;
% simple calculation of (2d) triangle area using heron formula. ;
%%%%%%%%;
n_t = size(p_23t___,3);
x_3t__ = reshape(p_23t___(1+0,:,:),[3,n_t]);
y_3t__ = reshape(p_23t___(1+1,:,:),[3,n_t]);
s01_t_ = sqrt( (x_3t__(1+1,:) - x_3t__(1+0,:)).^2 + (y_3t__(1+1,:) - y_3t__(1+0,:)).^2 );
s12_t_ = sqrt( (x_3t__(1+2,:) - x_3t__(1+1,:)).^2 + (y_3t__(1+2,:) - y_3t__(1+1,:)).^2 );
s20_t_ = sqrt( (x_3t__(1+0,:) - x_3t__(1+2,:)).^2 + (y_3t__(1+0,:) - y_3t__(1+2,:)).^2 );
s_t3__ = cat(2,s01_t_(:),s12_t_(:),s20_t_(:));
ss_t3__ = sort(s_t3__,2,'descend');
a_t_ = ss_t3__(:,1+0);
b_t_ = ss_t3__(:,1+1);
c_t_ = ss_t3__(:,1+2);
area_t_ = 0.25.*sqrt(max(0,(a_t_+(b_t_+c_t_)).*(c_t_-(a_t_-b_t_)).*(c_t_+(a_t_-b_t_)).*(a_t_+(b_t_-c_t_))));

