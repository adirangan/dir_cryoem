function [gradient_,E3_quad,alpha_] = transf_3d_gradient(t,y_,svd_tolerance,k_p_r_max,delta_p_r_max,n_delta_all,delta_c_0_all_,delta_c_1_all_,delta_c_2_all_,weight_delta_all_);
verbose=0;
%%%%%%%%%%%%%%%%;
% We assume that y_ stores delta_node_ (i.e., 3-x-n_delta_node) as well as the accumulated error IE. ;
%%%%%%%%%%%%%%%%;
h_ = @(kd) 4*pi*( sin(kd) - (kd).*cos(kd) ) ./ kd.^3 ; % calculates <f_j,f_k>, normalized so that <f,f> = 4*pi/3, rather than 4*pi/3*R_max^3 ;
dh_ = @(kd) 12*pi*( (kd.^2/3 - 1) .* sin(kd) + (kd).*cos(kd) ) ./ kd.^4 ;
n_delta_node = (size(y_,1)-1)/3; %<-- do not count the accumulated error IE. ;
delta_node_ = reshape(y_(1:3*n_delta_node),3,n_delta_node);
%%%%%%%%%%%%%%%%;
H_ = zeros(n_delta_node,n_delta_node);
for ndelta_nodeA=1:n_delta_node; for ndelta_nodeB=1:n_delta_node;
tmp_d = (delta_node_(:,ndelta_nodeB) - delta_node_(:,ndelta_nodeA));
if (ndelta_nodeA~=ndelta_nodeB); H_(ndelta_nodeA,ndelta_nodeB) = h_(2*pi*norm(tmp_d)*k_p_r_max); end;
if (ndelta_nodeA==ndelta_nodeB); H_(ndelta_nodeA,ndelta_nodeB) = (4*pi/3); end;
end;end;%for ndelta_nodeA=1:n_delta_node; for ndelta_nodeB=1:n_delta_node;
%[HU_,HS_,HV_] = svd(H_); HS_ = diag(HS_); H_ij_ = find(HS_./HS_(1) < 1e-6);
%GS_ = 1./HS_; GS_(H_ij_)=0; GS_ = diag(GS_);
%G_ = HV_*GS_*transpose(HU_);
G_ = pinv(H_,norm(H_)*1e-6);
dhx__ = zeros(n_delta_node,n_delta_node);
dhy__ = zeros(n_delta_node,n_delta_node);
dhz__ = zeros(n_delta_node,n_delta_node);
ddx_ = repmat(delta_node_(1+0,:),n_delta_node,1) - repmat(transpose(delta_node_(1+0,:)),1,n_delta_node);
ddy_ = repmat(delta_node_(1+1,:),n_delta_node,1) - repmat(transpose(delta_node_(1+1,:)),1,n_delta_node);
ddz_ = repmat(delta_node_(1+2,:),n_delta_node,1) - repmat(transpose(delta_node_(1+2,:)),1,n_delta_node);
hypot_dd_ = sqrt(ddx_.^2+ddy_.^2+ddz_.^2);
tmp_ = 2*pi*dh_(2*pi*hypot_dd_*k_p_r_max)*k_p_r_max./hypot_dd_; tmp_(find(hypot_dd_==0)) = 0;
dhx__ = - tmp_ .* ddx_;
dhy__ = - tmp_ .* ddy_;
dhz__ = - tmp_ .* ddz_;
%%%%%%%%%%%%%%%%;
dE3k_sample___ = zeros(n_delta_all,2,n_delta_node);
ddx_ = repmat(delta_c_0_all_',n_delta_node,1) - repmat(delta_node_(1+0,:)',1,n_delta_all);
ddy_ = repmat(delta_c_1_all_',n_delta_node,1) - repmat(delta_node_(1+1,:)',1,n_delta_all);
ddz_ = repmat(delta_c_2_all_',n_delta_node,1) - repmat(delta_node_(1+2,:)',1,n_delta_all);
hypot_dd_ = sqrt(ddx_.^2+ddy_.^2+ddz_.^2);
ij_=find(hypot_dd_==0);
ft_ = h_(2*pi*hypot_dd_*k_p_r_max); ft_(ij_) = (4*pi/3);
alpha_= G_*ft_;
E3_sample_ = ( (4*pi/3) - sum(alpha_.*ft_,1) ) / (4*pi/3); % Sampled relative error. ;
E3_quad = abs(E3_sample_)*weight_delta_all_ / delta_p_r_max^3;
tmp_ = 2*pi*dh_(2*pi*hypot_dd_*k_p_r_max)*k_p_r_max./hypot_dd_; tmp_(ij_) = 0;
dftx_ = - tmp_.*ddx_;
dfty_ = - tmp_.*ddy_;
dftz_ = - tmp_.*ddz_;
dE3x_sample_ = 2*alpha_.*(dhx__*alpha_) - 2*alpha_.*dftx_;
dE3y_sample_ = 2*alpha_.*(dhy__*alpha_) - 2*alpha_.*dfty_;
dE3z_sample_ = 2*alpha_.*(dhz__*alpha_) - 2*alpha_.*dftz_;
dE3x_quad_ = dE3x_sample_*weight_delta_all_ / delta_p_r_max^3;
dE3y_quad_ = dE3y_sample_*weight_delta_all_ / delta_p_r_max^3;
dE3z_quad_ = dE3z_sample_*weight_delta_all_ / delta_p_r_max^3;
%%%%%%%%%%%%%%%%;
gradient_ = -transpose([ dE3x_quad_ , dE3y_quad_ , dE3z_quad_ ]);
for ndelta_node=1:n_delta_node;
tmp = sqrt( delta_node_(1+0,ndelta_node).^2 + delta_node_(1+1,ndelta_node).^2 + delta_node_(1+2,ndelta_node).^2 );
if (tmp>1.25*delta_p_r_max); gradient_(:,ndelta_node) = -delta_node_(:,ndelta_node); end;
end;%for ndelta_node=1:n_delta_node;
gradient_ = reshape(gradient_,3*n_delta_node,1);
l2error = E3_quad;
gradient_(3*n_delta_node+1) = l2error ;
%%%%%%%%%%%%%%%%;
