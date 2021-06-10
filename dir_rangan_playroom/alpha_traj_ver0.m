function [theta_,phi_,gamma_,dx_,dy_] = alpha_traj_ver0(n_,dirname)

alpha_ = MDA_read_r8(sprintf('%s/alpha_tru_3_.mda',dirname)); 
[nrows,ncols] = size(alpha_);

theta_ = zeros(length(n_),ncols);
phi_ = zeros(length(n_),ncols);
gamma_ = zeros(length(n_),ncols);
dx_ = zeros(length(n_),ncols);
dy_ = zeros(length(n_),ncols);
for nF=1:length(n_);
ncur = n_(nF);
if (ncur==0);
alpha_ = MDA_read_r8(sprintf('%s/alpha_tru_3_.mda',dirname)); 
elseif (ncur>0);
alpha_ = MDA_read_r8(sprintf('%s/alpha_est_%d_.mda',dirname,ncur)); 
end;%if;
theta_(nF,:) = periodize(alpha_(1,:),0,pi);
phi_(nF,:) = periodize(alpha_(2,:),0,2*pi);
gamma_(nF,:) = periodize(alpha_(3,:),-pi,pi);
dx_(nF,:) = alpha_(4,:);
dy_(nF,:) = alpha_(5,:);
end;% for nF;

n_crop_ = n_(1:end-1); n_crop_ = n_crop_(:);

Y_{1} = theta_;
Y_{2} = phi_;
Y_{3} = gamma_;
Y_{4} = dx_;
Y_{5} = dy_;
ylim_{1} = [0,pi]; plim_{1} = [-pi/2,pi/2];
ylim_{2} = [0,2*pi]; plim_{2} = [-pi/1,pi/1];
ylim_{3} = [0,2*pi]; plim_{3} = [-pi/1,pi/1];
ylim_{4} = [0,0.1]; plim_{4} = [-1,1];
ylim_{5} = [0,0.1]; plim_{5} = [-1,1];
tstr_{1} = 'theta';
tstr_{2} = 'phi';
tstr_{3} = 'gamma';
tstr_{4} = 'dx';
tstr_{5} = 'dy';

figure;
for nd=1:5;
%tmp_ = abs(diff(Y_{nd},1,1));
tmp_ = abs(periodize(diff(Y_{nd},1,1),plim_{nd}(1),plim_{nd}(2)));
tmp_avg = mean(tmp_,2);
tmp_p05 = prctile(tmp_,05,2);
tmp_p15 = prctile(tmp_,15,2);
tmp_p50 = prctile(tmp_,50,2);
tmp_p85 = prctile(tmp_,85,2);
tmp_p95 = prctile(tmp_,95,2);
subplot(2,3,nd); hold on;
plot(n_crop_,tmp_p05,'-','LineWidth',0.5,'Color',0.75*[1,1,1]);
plot(n_crop_,tmp_p15,'-','LineWidth',1.0,'Color',0.50*[1,1,1]);
plot(n_crop_,tmp_p50,'-','LineWidth',2.0,'Color',0.00*[1,1,1]);
plot(n_crop_,tmp_p85,'-','LineWidth',1.0,'Color',0.50*[1,1,1]);
plot(n_crop_,tmp_p95,'-','LineWidth',0.5,'Color',0.75*[1,1,1]);
hold off;
xlim([min(n_crop_),max(n_crop_)]); ylim(ylim_{nd}); title(sprintf('Delta %s',tstr_{nd}));
end;%for nd=1:5;

set(gcf,'Position',[1,1,1024,768]);
print('-djpeg',sprintf('%s/alpha_traj_%d_%d_%d.jpg',dirname,n_(1),round(mean(diff(n_))),n_(end)));


%{
ez = ones(size(n_));
figure;
for nc=1:ncols;
subplot(2,3,1); cla; hold on;
%plot(n_,deperiodize(theta_(:,nc),0,pi),'.-'); 
plot(n_,0*ez,'k-',n_,pi*ez,'k-');
xlim([min(n_),max(n_)]);
ylim([-pi,2*pi]);
hold off;
subplot(2,3,2); cla; hold on;
plot(n_,deperiodize(phi_(:,nc),0,2*pi),'.-');
plot(n_,0*ez,'k-',n_,2*pi*ez,'k-'); 
hold off;
xlim([min(n_),max(n_)]);
ylim([-2*pi,4*pi]);
subplot(2,3,3); cla; hold on;
plot(n_,deperiodize(gamma_(:,nc),-pi,pi),'.-');
plot(n_,-pi*ez,'k-',n_,+pi*ez,'k-'); 
hold off;
xlim([min(n_),max(n_)]);
ylim([-2*pi,+2*pi]);
subplot(2,3,4); cla; plot(n_,dx_(:,nc),'.-'); ylim(0.1*[-1,1]);xlim([min(n_),max(n_)]);
subplot(2,3,5); cla; plot(n_,dy_(:,nc),'.-'); ylim(0.1*[-1,1]);xlim([min(n_),max(n_)]);
%subplot(2,3,6); plot(deperiodize(periodize(pi+0.2*randn(64,1),-pi,pi),-pi,pi),'o-'); ylim([-2*pi,+2*pi]);
set(gcf,'Position',[1,1,1024,768]);
drawnow(); pause();
end;%for nc=1:ncols;
 %}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output = periodize(input,a,b);
output = a+mod(input-a,(b-a));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output = deperiodize(input,a,b);
dlo = -0.5*(b-a); dhi = +0.5*(b-a);
output = input(:);
output = output(1) + periodize([0;diff(output)],dlo,dhi);
ouput = reshape(output,size(input));
