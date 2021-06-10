function gradient_ = test_F_26_gradient(t,y_,svd_tolerance,R_max,D_target,n_A,delta_sample_,quad_weight_sample_);
verbose=0;
%%%%%%%%%%%%%%%%;
% We assume that y_ stores delta_node_ (i.e., 3-x-n_node) as well as the accumulated error IE. ;
%%%%%%%%%%%%%%%%;
h_ = @(kd) 4*pi*( sin(kd) - (kd).*cos(kd) ) ./ kd.^3 ; % calculates <f_j,f_k>, normalized so that <f,f> = 4*pi/3, rather than 4*pi/3*R_max^3 ;
dh_ = @(kd) 12*pi*( (kd.^2/3 - 1) .* sin(kd) + (kd).*cos(kd) ) ./ kd.^4 ;
n_node = (size(y_,1)-1)/3; %<-- do not count the accumulated error IE. ;
delta_node_ = reshape(y_(1:3*n_node),3,n_node);
%%%%%%%%%%%%%%%%;
H_ = zeros(n_node,n_node);
for nnodeA=1:n_node; for nnodeB=1:n_node;
d = (delta_node_(:,nnodeB) - delta_node_(:,nnodeA));
if (nnodeA~=nnodeB); H_(nnodeA,nnodeB) = h_(norm(d,'fro')*R_max)*R_max^3; end;
if (nnodeA==nnodeB); H_(nnodeA,nnodeB) = (4*pi/3*R_max^3); end;
end;end;%for nnodeA=1:n_node; for nnodeB=1:n_node;
[HU_,HS_,HV_] = svd(H_); HS_ = diag(HS_); H_ij_ = find(HS_./HS_(1) < svd_tolerance*R_max^3);
GS_ = 1./HS_; GS_(H_ij_)=0; GS_ = diag(GS_); G_ = HV_*GS_*transpose(HU_);
dhx__ = zeros(n_node,n_node); dhy__ = zeros(n_node,n_node); dhz__ = zeros(n_node,n_node);
ddHx_ = repmat(delta_node_(1,:),n_node,1) - repmat(transpose(delta_node_(1,:)),1,n_node);
ddHy_ = repmat(delta_node_(2,:),n_node,1) - repmat(transpose(delta_node_(2,:)),1,n_node);
ddHz_ = repmat(delta_node_(3,:),n_node,1) - repmat(transpose(delta_node_(3,:)),1,n_node);
hypot_ddH_ = sqrt(ddHx_.^2+ddHy_.^2+ddHz_.^2);
tmp_ = dh_(hypot_ddH_*R_max)*R_max^4./hypot_ddH_; tmp_(find(hypot_ddH_==0)) = 0;
dhx__ = - tmp_ .* ddHx_; dhy__ = - tmp_ .* ddHy_; dhz__ = - tmp_ .* ddHz_;
%%%%%%%%%%%%%%%%;
ddfx_ = repmat(transpose(delta_sample_(:,1)),n_node,1) - repmat(transpose(delta_node_(1,:)),1,n_A);
ddfy_ = repmat(transpose(delta_sample_(:,2)),n_node,1) - repmat(transpose(delta_node_(2,:)),1,n_A);
ddfz_ = repmat(transpose(delta_sample_(:,3)),n_node,1) - repmat(transpose(delta_node_(3,:)),1,n_A);
hypot_ddf_ = sqrt(ddfx_.^2+ddfy_.^2+ddfz_.^2);
ij_=find(hypot_ddf_==0);
ft_ = h_(hypot_ddf_*R_max)*R_max^3; ft_(ij_) = (4*pi/3*R_max^3);
alpha_= G_*ft_;
E3_sample_ = ( (4*pi/3*R_max^3) - sum(alpha_.*ft_,1) ) / (4*pi/3*R_max^3); % Sampled relative error. ;
E3_quad = abs(E3_sample_)*quad_weight_sample_ / (4*pi/3*D_target^3);
tmp_ = dh_(hypot_ddf_*R_max)*R_max^4./hypot_ddf_; tmp_(ij_) = 0;
dftx_ = - tmp_.*ddfx_; dfty_ = - tmp_.*ddfy_; dftz_ = - tmp_.*ddfz_;
dE3x_sample_ = (2*alpha_.*(dhx__*alpha_) - 2*alpha_.*dftx_) / (4*pi/3*R_max^3);
dE3y_sample_ = (2*alpha_.*(dhy__*alpha_) - 2*alpha_.*dfty_) / (4*pi/3*R_max^3);
dE3z_sample_ = (2*alpha_.*(dhz__*alpha_) - 2*alpha_.*dftz_) / (4*pi/3*R_max^3);
dE3x_quad_ = dE3x_sample_*quad_weight_sample_;
dE3y_quad_ = dE3y_sample_*quad_weight_sample_;
dE3z_quad_ = dE3z_sample_*quad_weight_sample_;
%%%%%%%%%%%%%%%%;
% Here we test a particular delta_sample. ;
%%%%%%%%%%%%%%%%;
if (verbose>0);
delta_test_ = D_target*rand(1,3);
ddfx_test_ = repmat(delta_test_(:,1),n_node,1) - transpose(delta_node_(1,:));
ddfy_test_ = repmat(delta_test_(:,2),n_node,1) - transpose(delta_node_(2,:));
ddfz_test_ = repmat(delta_test_(:,3),n_node,1) - transpose(delta_node_(3,:));
hypot_ddf_test_ = sqrt(ddfx_test_.^2 + ddfy_test_.^2 + ddfz_test_.^2);
ij_test_ = find(hypot_ddf_test_==0);
ft_test_ = h_(hypot_ddf_test_*R_max)*R_max^3; ft_test_(ij_test_) = (4*pi/3*R_max^3);
alpha_test_ = G_*ft_test_;
delta_upd_ = randn(3,n_node); for nnode=1:n_node; delta_upd_(:,nnode) = delta_upd_(:,nnode)/norm(delta_upd_(:,nnode),'fro'); end;
eps_upd_ = 0.1.^[1:6]; n_eps_upd = length(eps_upd_);
for neps_upd=1:n_eps_upd;
eps_upd = eps_upd_(neps_upd);
%%%%%%%%;
delta_updA_ = delta_node_ + eps_upd*delta_upd_;
ddfx_updA_ = repmat(delta_test_(:,1),n_node,1) - transpose(delta_updA_(1,:));
ddfy_updA_ = repmat(delta_test_(:,2),n_node,1) - transpose(delta_updA_(2,:));
ddfz_updA_ = repmat(delta_test_(:,3),n_node,1) - transpose(delta_updA_(3,:));
hypot_ddf_updA_ = sqrt(ddfx_updA_.^2 + ddfy_updA_.^2 + ddfz_updA_.^2);
ij_updA_ = find(hypot_ddf_updA_==0);
ft_updA_ = h_(hypot_ddf_updA_*R_max)*R_max^3; ft_updA_(ij_updA_) = (4*pi/3*R_max^3);
H_updA_ = zeros(n_node,n_node);
for nnodeA=1:n_node; for nnodeB=1:n_node;
d = (delta_updA_(:,nnodeB) - delta_updA_(:,nnodeA));
if (nnodeA~=nnodeB); H_updA_(nnodeA,nnodeB) = h_(norm(d,'fro')*R_max)*R_max^3; end;
if (nnodeA==nnodeB); H_updA_(nnodeA,nnodeB) = (4*pi/3*R_max^3); end;
end;end;%for nnodeA=1:n_node; for nnodeB=1:n_node;
[HU_,HS_,HV_] = svd(H_updA_); HS_ = diag(HS_); H_ij_ = find(HS_./HS_(1) < svd_tolerance*R_max^3);
GS_ = 1./HS_; GS_(H_ij_)=0; GS_ = diag(GS_); G_updA_ = HV_*GS_*transpose(HU_);
alpha_updA_ = G_updA_*ft_updA_;
E3_updA_ = ( (4*pi/3*R_max^3) - sum(alpha_updA_.*ft_updA_,1) ) / (4*pi/3*R_max^3) ;
%%%%%%%%;
delta_updB_ = delta_node_ - eps_upd*delta_upd_;
ddfx_updB_ = repmat(delta_test_(:,1),n_node,1) - transpose(delta_updB_(1,:));
ddfy_updB_ = repmat(delta_test_(:,2),n_node,1) - transpose(delta_updB_(2,:));
ddfz_updB_ = repmat(delta_test_(:,3),n_node,1) - transpose(delta_updB_(3,:));
hypot_ddf_updB_ = sqrt(ddfx_updB_.^2 + ddfy_updB_.^2 + ddfz_updB_.^2);
ij_updB_ = find(hypot_ddf_updB_==0);
ft_updB_ = h_(hypot_ddf_updB_*R_max)*R_max^3; ft_updB_(ij_updB_) = (4*pi/3*R_max^3);
H_updB_ = zeros(n_node,n_node);
for nnodeA=1:n_node; for nnodeB=1:n_node;
d = (delta_updB_(:,nnodeB) - delta_updB_(:,nnodeA));
if (nnodeA~=nnodeB); H_updB_(nnodeA,nnodeB) = h_(norm(d,'fro')*R_max)*R_max^3; end;
if (nnodeA==nnodeB); H_updB_(nnodeA,nnodeB) = (4*pi/3*R_max^3); end;
end;end;%for nnodeA=1:n_node; for nnodeB=1:n_node;
[HU_,HS_,HV_] = svd(H_updB_); HS_ = diag(HS_); H_ij_ = find(HS_./HS_(1) < svd_tolerance*R_max^3);
GS_ = 1./HS_; GS_(H_ij_)=0; GS_ = diag(GS_); G_updB_ = HV_*GS_*transpose(HU_);
alpha_updB_ = G_updB_*ft_updB_;
E3_updB_ = ( (4*pi/3*R_max^3) - sum(alpha_updB_.*ft_updB_,1) ) / (4*pi/3*R_max^3) ;
%%%%%%%%;
ft_updC_ = (ft_updA_ - ft_updB_)/(2*eps_upd);
H_updC_ = (H_updA_ - H_updB_)/(2*eps_upd);
G_updC_ = (G_updA_ - G_updB_)/(2*eps_upd);
alpha_updC_ = (alpha_updA_ - alpha_updB_)/(2*eps_upd);
E3_updC_ = (E3_updA_ - E3_updB_)/(2*eps_upd);
%%%%%%%%
tmp_ = dh_(hypot_ddf_test_*R_max)*R_max^4./hypot_ddf_test_; tmp_(ij_) = 0;
dftx_test_ = - tmp_.*ddfx_test_; dfty_test_ = - tmp_.*ddfy_test_; dftz_test_ = - tmp_.*ddfz_test_;
dft_test_ = dftx_test_.*transpose(delta_upd_(1,:)) + dfty_test_.*transpose(delta_upd_(2,:)) + dftz_test_.*transpose(delta_upd_(3,:));
ft_updD_ = dft_test_;
H_updD_ = zeros(n_node,n_node);
H_updD_ = ...
  + dhx__.*( - repmat(delta_upd_(1,:),n_node,1) + repmat(transpose(delta_upd_(1,:)),1,n_node) ) ...
  + dhy__.*( - repmat(delta_upd_(2,:),n_node,1) + repmat(transpose(delta_upd_(2,:)),1,n_node) ) ...
  + dhz__.*( - repmat(delta_upd_(3,:),n_node,1) + repmat(transpose(delta_upd_(3,:)),1,n_node) ) ...
  ;
H_updE_ = ...
  + dhx__.*( transpose(delta_upd_(1,:))*ones(1,n_node) - ones(n_node,1)*delta_upd_(1,:) ) ...
  + dhy__.*( transpose(delta_upd_(2,:))*ones(1,n_node) - ones(n_node,1)*delta_upd_(2,:) ) ...
  + dhz__.*( transpose(delta_upd_(3,:))*ones(1,n_node) - ones(n_node,1)*delta_upd_(3,:) ) ...
  ;
G_updD_ = -G_*H_updD_*G_;
alpha_updD_ = G_updD_*ft_test_ + G_*ft_updD_;
E3_updD_ = - ( sum(alpha_updD_.*ft_test_ + alpha_test_.*ft_updD_,1) ) / (4*pi/3*R_max^3);
%E3_updE_ = - ( 2*sum(ft_updD_.*alpha_test_) - transpose(alpha_test_)*H_updD_*alpha_test_ )  / (4*pi/3*R_max^3);
E3_updE_ = - ( ...
	       + 2*transpose(dftx_test_.*alpha_test_)*transpose(delta_upd_(1,:)) ... 
	       + 2*transpose(dfty_test_.*alpha_test_)*transpose(delta_upd_(2,:)) ... 
	       + 2*transpose(dftz_test_.*alpha_test_)*transpose(delta_upd_(3,:)) ... 
	       - 2*transpose(alpha_test_.*(dhx__*alpha_test_))*transpose(delta_upd_(1,:)) ...
	       - 2*transpose(alpha_test_.*(dhy__*alpha_test_))*transpose(delta_upd_(2,:)) ...
	       - 2*transpose(alpha_test_.*(dhz__*alpha_test_))*transpose(delta_upd_(3,:)) ...
	       )  / (4*pi/3*R_max^3);
%%%%%%%%;
disp(sprintf(' %% eps_upd %0.16f: norm(ft_updC_-ft_updD_) \t %0.16f',eps_upd,norm(ft_updC_-ft_updD_,'fro')));
disp(sprintf(' %% eps_upd %0.16f: norm(H_updC_-H_updD_) \t %0.16f',eps_upd,norm(H_updC_-H_updD_,'fro')));
disp(sprintf(' %% eps_upd %0.16f: norm(H_updC_-H_updE_) \t %0.16f',eps_upd,norm(H_updC_-H_updE_,'fro')));
disp(sprintf(' %% eps_upd %0.16f: norm(G_updC_-G_updD_) \t %0.16f',eps_upd,norm(G_updC_-G_updD_,'fro')));
disp(sprintf(' %% eps_upd %0.16f: norm(al_updC_-al_updD_) \t %0.16f',eps_upd,norm(alpha_updC_-alpha_updD_,'fro')));
disp(sprintf(' %% eps_upd %0.16f: norm(E3_updC_-E3_updD_) \t %0.16f',eps_upd,norm(E3_updC_-E3_updD_,'fro')));
disp(sprintf(' %% eps_upd %0.16f: norm(E3_updC_-E3_updE_) \t %0.16f',eps_upd,norm(E3_updC_-E3_updE_,'fro')));
%%%%%%%%;
end;%for neps_upd=1:n_eps_upd;
end;%if (verbose>0);
%%%%%%%%%%%%%%%%;
gradient_ = -transpose([ dE3x_quad_ , dE3y_quad_ , dE3z_quad_ ]);
for nnode=1:n_node;
tmp = sqrt( delta_node_(1,nnode).^2 + delta_node_(2,nnode).^2 + delta_node_(3,nnode).^2 );
if (tmp>1.25*D_target); gradient_(:,nnode) = -delta_node_(:,nnode); end;
end;%for nnode=1:n_node;
gradient_ = reshape(gradient_,3*n_node,1);
l2error = E3_quad;
gradient_(3*n_node+1) = l2error ;
