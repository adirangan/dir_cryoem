function residual_scatter_ver0(n_,dirname,infix) ;
% for use with test_ver16.f, which has model_A,B,C, as well as Y_tru_ and Y_est_ ;
n_ = n_(find(n_>=1));
figure;

if length(n_)>1; prows = 2; pcols = ceil(length(n_)/prows); end;
if length(n_)<=1; prows = 1; pcols = 1; end;

if length(n_)==1;
ncur = n_(1);
I_model_sample_ = MDA_read_i4(sprintf('%s/I_model_sample_%s_%d_.mda',dirname,infix,ncur));
M_residual_loading_ = MDA_read_c16(sprintf('%s/M_residual_loading_%s_%d_.mda',dirname,infix,ncur));
%M1_lim_ = mean(real(M_residual_loading_(1,:))) + 1.5*std(real(M_residual_loading_(1,:)))*[-1,+1];
%M2_lim_ = mean(real(M_residual_loading_(2,:))) + 1.5*std(real(M_residual_loading_(2,:)))*[-1,+1];
M1_lim_ = [min(real(M_residual_loading_(1,:))),max(real(M_residual_loading_(1,:)))];
M2_lim_ = [min(real(M_residual_loading_(2,:))),max(real(M_residual_loading_(2,:)))];
model_u_ = unique(I_model_sample_);
n_model = length(model_u_);

subplot(1,2,1);
hold on; 
c_ = colormap('hsv'); n_c = size(c_,1);
for nmodel=1:n_model;
tmp_ij_ = find(I_model_sample_==model_u_(nmodel));
nc = max(1,min(n_c,floor(n_c*nmodel/n_model)));
plot(real(M_residual_loading_(1,tmp_ij_)),real(M_residual_loading_(2,tmp_ij_)),'.','Color',c_(nc,:),'MarkerSize',15);
end;%for nmodel=1:n_model;
hold off;
xlim(M1_lim_);
ylim(M2_lim_);
xlabel('M1'); ylabel('M2');
title(sprintf('nk %d',ncur));

h_ = cell(n_model,1); nbins = 32;
for nmodel=1:n_model;
tmp_ij_ = find(I_model_sample_==model_u_(nmodel));
h_{nmodel} = hist2d_0(real(M_residual_loading_(1,tmp_ij_)),real(M_residual_loading_(2,tmp_ij_)),nbins,nbins,M1_lim_,M2_lim_);
colormap(colormap_beach());
subplot(2,4,(nmodel+2)*(nmodel<3) + (nmodel+4)*(nmodel>=3));
imagesc(log(1+h_{nmodel})); axis image; set(gca,'YDir','normal');
set(gca,'XTick',[],'YTick',[]);
title(sprintf('model %d',nmodel));
end;%for nmodel=1:n_model;
set(gcf,'Position',1+[0,0,1024*2,1024]);
fname = sprintf('%s/residual_scatter_%s_%d.jpg',dirname,infix,n_(1));
print('-djpeg',fname);
end;%if length(n_)==1;

if length(n_)>1;

for nF=1:length(n_);
ncur = n_(nF);

if (0);
elseif (ncur>1);
I_model_sample_ = MDA_read_i4(sprintf('%s/I_model_sample_%s_%d_.mda',dirname,infix,ncur));
M_residual_loading_ = MDA_read_c16(sprintf('%s/M_residual_loading_%s_%d_.mda',dirname,infix,ncur));
end;% if;

subplot(prows,pcols,nF);
hold on; 
c_ = colormap('hsv'); n_c = size(c_,1);
model_u_ = unique(I_model_sample_);
n_model = length(model_u_);
for nmodel=1:n_model;
tmp_ij_ = find(I_model_sample_==model_u_(nmodel));
nc = max(1,min(n_c,floor(n_c*nmodel/n_model)));
if length(size(M_residual_loading_))==2;
plot(real(M_residual_loading_(1,tmp_ij_)),real(M_residual_loading_(2,tmp_ij_)),'.','Color',c_(nc,:),'MarkerSize',15);
end;%if length(size(M_residual_loading_))==2;
if length(size(M_residual_loading_))>=3;
plot3(real(M_residual_loading_(1,tmp_ij_)),real(M_residual_loading_(2,tmp_ij_)),real(M_residual_loading_(3,tmp_ij_)),'.','Color',c_(nc,:),'MarkerSize',15);
end;%if length(size(M_residual_loading_))>=3;
end;%for nmodel=1:n_model;
hold off;

if length(size(M_residual_loading_))==2;
xlabel('M1'); ylabel('M2');
end;%if length(size(M_residual_loading_))==2;
if length(size(M_residual_loading_))>=3;
xlabel('M1'); ylabel('M2'); zlabel('M3'); axis vis3d;
end;%if length(size(M_residual_loading_))>=3;

title(sprintf('nk %d',ncur));

end;% for nF;

if length(n_)>2; fname = sprintf('%s/residual_scatter_%s_%d_%d_%d.jpg',dirname,infix,n_(2),round(mean(diff(n_(2:end)))),n_(end)); end;
if length(n_)==2; fname = sprintf('%s/residual_scatter_%s_%d_%d.jpg',dirname,infix,n_(1),n_(2)); end;
if length(n_)<=1; fname = sprintf('%s/residual_scatter_%s_%d.jpg',dirname,infix,n_(1)); end;

set(gcf,'Position',[1,1,1024,768]);
print('-djpeg',fname);

end;%if length(n_)>1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
