function alphadist_ver3(n_,dirname)
figure;

prows = 7; pcols = length(n_);
for nF=1:length(n_);
ncur = n_(nF);
if (ncur==0);
alpha_ = MDA_read_r8(sprintf('%s/alpha_tru__3_.mda',dirname)); 
elseif (ncur>0);
alpha_ = MDA_read_r8(sprintf('%s/alpha_est__%d_.mda',dirname,ncur)); 
end;%if;
polar_a_ = periodize(alpha_(1,:),0,pi);
azimu_b_ = periodize(alpha_(2,:),0,2*pi);
gamma_z_ = periodize(alpha_(3,:),-pi,pi);
delta_x_ = alpha_(4,:);
delta_y_ = alpha_(5,:);
subplot(prows,pcols,nF + (0-0)*pcols); imagesc(hist2d_0(azimu_b_,polar_a_,24,12,[0,2*pi],[0,pi])); set(gca,'XTick',[],'YTick',[]);
subplot(prows,pcols,nF + (1-0)*pcols); plot(polar_a_,'.'); xlim([1,length(polar_a_)]); ylim([0,+pi]);
subplot(prows,pcols,nF + (2-0)*pcols); plot(azimu_b_,'.'); xlim([1,length(azimu_b_)]); ylim([0,+2*pi]);
subplot(prows,pcols,nF + (3-0)*pcols); plot(gamma_z_,'.'); xlim([1,length(gamma_z_)]); ylim([-pi,+pi]);
subplot(prows,pcols,nF + (4-0)*pcols); plot(delta_x_,'.'); xlim([1,length(delta_x_)]); ylim([-0.1,+0.1]);
subplot(prows,pcols,nF + (5-0)*pcols); plot(delta_y_,'.'); xlim([1,length(delta_y_)]); ylim([-0.1,+0.1]);
subplot(prows,pcols,nF + (6-0)*pcols); imagesc(hist2d_0(delta_x_,delta_y_,18,18,[-0.1,+0.1],[-0.1,+0.1])); set(gca,'XTick',[],'YTick',[]);
end;% for nF;
set(gcf,'Position',[1,1,1024,768]);
print('-djpeg',sprintf('%s/alpha1d_est_%d_%d_%d.jpg',dirname,n_(1),round(mean(diff(n_))),n_(end)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output = periodize(input,a,b);
output = a+mod(input-a,(b-a));
