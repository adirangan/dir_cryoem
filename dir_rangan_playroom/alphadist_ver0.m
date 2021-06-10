function alphadist_ver0(n_,dirname)
figure;

%nlow = 3;
%nmax = 10;
%n_ = [nlow:2:nmax];
prows = 5; pcols = length(n_);
for nF=1:length(n_);
ncur = n_(nF);
if (ncur==0);
alpha_ = MDA_read_r8(sprintf('%s/alpha_tru_3_.mda',dirname,ncur)); 
elseif (ncur>0);
alpha_ = MDA_read_r8(sprintf('%s/alpha_est_%d_.mda',dirname,ncur)); 
end;%if;
theta_ = periodize(alpha_(1,:),0,pi);
phi_ = periodize(alpha_(2,:),0,2*pi);
gamma_ = periodize(alpha_(3,:),-pi,pi);
dx_ = alpha_(4,:);
dy_ = alpha_(5,:);
subplot(prows,pcols,nF + (1-1)*pcols); plot(theta_,'.'); xlim([1,length(theta_)]); ylim([0,+pi]);
subplot(prows,pcols,nF + (2-1)*pcols); plot(phi_,'.'); xlim([1,length(phi_)]); ylim([0,+2*pi]);
subplot(prows,pcols,nF + (3-1)*pcols); plot(gamma_,'.'); xlim([1,length(gamma_)]); ylim([-pi,+pi]);
subplot(prows,pcols,nF + (4-1)*pcols); plot(dx_,'.'); xlim([1,length(dx_)]); ylim([-0.1,+0.1]);
subplot(prows,pcols,nF + (5-1)*pcols); plot(dy_,'.'); xlim([1,length(dy_)]); ylim([-0.1,+0.1]);
end;% for nF;
set(gcf,'Position',[1,1,1024,768]);
print('-djpeg',sprintf('%s/alpha_est_%d_%d_%d.jpg',dirname,n_(1),round(mean(diff(n_))),n_(end)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output = periodize(input,a,b);
output = a+mod(input-a,(b-a));
