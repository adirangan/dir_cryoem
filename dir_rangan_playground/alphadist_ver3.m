function alphadist_ver3(n_,dirname,infix,delta_max)
n_ = n_(find(n_>=0));
figure;
colormap(colormap_beach());
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
if nargin<4; delta_max = 0.1; end;
subplot(prows,pcols,nF + (0-0)*pcols); imagesc(hist2d_0(azimu_b_,polar_a_,24,12,[0,2*pi],[0,pi])); 
xlabel('azimu_b'); ylabel('polar_a'); title(sprintf('nk %d',ncur)); set(gca,'XTick',[],'YTick',[]);
subplot(prows,pcols,nF + (1-0)*pcols); plot(polar_a_,'.'); xlim([1,length(polar_a_)]); ylim([0,+pi]);
xlabel('image #'); ylabel('polar_a'); set(gca,'XTick',[],'YTick',[]);
subplot(prows,pcols,nF + (2-0)*pcols); plot(azimu_b_,'.'); xlim([1,length(azimu_b_)]); ylim([0,+2*pi]);
xlabel('image #'); ylabel('azimu_b'); set(gca,'XTick',[],'YTick',[]);
subplot(prows,pcols,nF + (3-0)*pcols); plot(gamma_z_,'.'); xlim([1,length(gamma_z_)]); ylim([-pi,+pi]);
xlabel('image #'); ylabel('gamma_z'); set(gca,'XTick',[],'YTick',[]);
subplot(prows,pcols,nF + (4-0)*pcols); plot(delta_x_,'.'); xlim([1,length(delta_x_)]); ylim([-delta_max,+delta_max]);
xlabel('image #'); ylabel('delta_x'); set(gca,'XTick',[],'YTick',[]);
subplot(prows,pcols,nF + (5-0)*pcols); plot(delta_y_,'.'); xlim([1,length(delta_y_)]); ylim([-delta_max,+delta_max]);
xlabel('image #'); ylabel('delta_y'); set(gca,'XTick',[],'YTick',[]);
subplot(prows,pcols,nF + (6-0)*pcols); imagesc(hist2d_0(delta_x_,delta_y_,18,18,[-delta_max,+delta_max],[-delta_max,+delta_max])); xlabel('delta_x'); ylabel('delta_y'); set(gca,'XTick',[],'YTick',[]);
set(gca,'XTick',[],'YTick',[]);
end;% for nF;
figbig;
print('-djpeg',sprintf('%s/alpha1d_%s_%d_%d_%d.jpg',dirname,infix,n_(1),round(mean(diff(n_))),n_(end)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output = periodize(input,a,b);
output = a+mod(input-a,(b-a));
