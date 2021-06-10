% test 2d-fourier transform on realistic image ; 
% translating image in real-space ;

ng = 64; cmax = 1;
xg = linspace(0,cmax,ng+1); xg = xg(1:end-1); dx = mean(diff(xg));
kg = (0:ng-1)/cmax; dk = mean(diff(kg));
[Xc,Yc] = meshgrid(xg,xg); Ac = Xc-cmax/2; Bc = Yc-cmax/2;
Rc = sqrt(Ac.^2 + Bc.^2); Wc = atan2(Bc,Ac);
sx=0.1;mx=0.3;sy=0.2;my=0.2;tmp_a = 1/(2*pi*cmax*sx*cmax*sy)*exp(-(Xc - cmax*mx).^2/2/(cmax*sx).^2).*exp(-(Yc - cmax*my).^2/2/(cmax*sy).^2);
sx=0.1;mx=0.7;sy=0.2;my=0.2;tmp_b = 1/(2*pi*cmax*sx*cmax*sy)*exp(-(Xc - cmax*mx).^2/2/(cmax*sx).^2).*exp(-(Yc - cmax*my).^2/2/(cmax*sy).^2);
sx=0.1;mx=0.7;sy=0.2;my=0.8;tmp_c = 1/(2*pi*cmax*sx*cmax*sy)*exp(-(Xc - cmax*mx).^2/2/(cmax*sx).^2).*exp(-(Yc - cmax*my).^2/2/(cmax*sy).^2);
sx=0.1;mx=0.3;sy=0.2;my=0.8;tmp_d = 1/(2*pi*cmax*sx*cmax*sy)*exp(-(Xc - cmax*mx).^2/2/(cmax*sx).^2).*exp(-(Yc - cmax*my).^2/2/(cmax*sy).^2);
%tmp_crop = (Xc<0.75*cmax & Xc>0.25*cmax).*(Yc<0.75*cmax & Yc>0.25*cmax);
tmp_crop = Rc<0.25*(1 + 0.5*cos(5*Wc));
Fc1 = tmp_crop.*(tmp_a+tmp_b+tmp_c+tmp_d);
Lc1 = sum(abs(Fc1(:)).^2)*dx*dx;
Fk1 = fft2(Fc1) * dx*dx;
Lk1 = sum(abs(Fk1(:)).^2)*dk*dk;
disp(sprintf(' %% Lc1 %f Lk1 %f',Lc1,Lk1));
figure;
subplot(2,2,1);imagesc(real(Fc1));title('real(Fc1)'); colorbar; 
set(gca,'XTick',[1,ng/2,ng],'XTickLabel',[xg(1),xg(end/2),xg(end)]);set(gca,'YTick',[1,ng/2,ng],'YTickLabel',[xg(1),xg(end/2),xg(end)]);
%subplot(2,2,2);imagesc(imag(Fc1));title('imag(Fc1)'); colorbar;
%set(gca,'XTick',[1,ng/2,ng],'XTickLabel',[xg(1),xg(end/2),xg(end)]);set(gca,'YTick',[1,ng/2,ng],'YTickLabel',[xg(1),xg(end/2),xg(end)]);
subplot(2,2,2);imagesc(abs(Fk1));title('abs(Fk1)'); colorbar;
set(gca,'XTick',[1,ng/2,ng],'XTickLabel',[kg(1),kg(end/2),kg(end)]);set(gca,'YTick',[1,ng/2,ng],'YTickLabel',[kg(1),kg(end/2),kg(end)]);
subplot(2,2,3);imagesc(real(Fk1));title('real(Fk1)'); colorbar;
set(gca,'XTick',[1,ng/2,ng],'XTickLabel',[kg(1),kg(end/2),kg(end)]);set(gca,'YTick',[1,ng/2,ng],'YTickLabel',[kg(1),kg(end/2),kg(end)]);
subplot(2,2,4);imagesc(imag(Fk1));title('imag(Fk1)'); colorbar;
set(gca,'XTick',[1,ng/2,ng],'XTickLabel',[kg(1),kg(end/2),kg(end)]);set(gca,'YTick',[1,ng/2,ng],'YTickLabel',[kg(1),kg(end/2),kg(end)]);
%gamma = 2*pi/10; Fc2 = interp2(Ac,Bc,Fc1,Ac*cos(gamma) - Bc*sin(gamma),Ac*sin(gamma) + Bc*cos(gamma),'nearest',0);
delta_x = cmax/64; delta_y = cmax/64; Fc2 = interp2(Ac,Bc,Fc1,Ac-delta_x,Bc-delta_y,'nearest',0);
Lc2 = sum(abs(Fc2(:)).^2)*dx*dx;
Fk2 = fft2(Fc2) * dx*dx;
Lk2 = sum(abs(Fk2(:)).^2)*dk*dk;
disp(sprintf(' %% Lc2 %f Lk2 %f',Lc2,Lk2));
figure;
subplot(2,2,1);imagesc(real(Fc2));title('real(Fc2)'); colorbar; 
set(gca,'XTick',[1,ng/2,ng],'XTickLabel',[xg(1),xg(end/2),xg(end)]);set(gca,'YTick',[1,ng/2,ng],'YTickLabel',[xg(1),xg(end/2),xg(end)]);
%subplot(2,2,2);imagesc(imag(Fc2));title('imag(Fc2)'); colorbar;
%set(gca,'XTick',[1,ng/2,ng],'XTickLabel',[xg(1),xg(end/2),xg(end)]);set(gca,'YTick',[1,ng/2,ng],'YTickLabel',[xg(1),xg(end/2),xg(end)]);
subplot(2,2,2);imagesc(abs(Fk2));title('abs(Fk2)'); colorbar;
set(gca,'XTick',[1,ng/2,ng],'XTickLabel',[kg(1),kg(end/2),kg(end)]);set(gca,'YTick',[1,ng/2,ng],'YTickLabel',[kg(1),kg(end/2),kg(end)]);
subplot(2,2,3);imagesc(real(Fk2));title('real(Fk2)'); colorbar;
set(gca,'XTick',[1,ng/2,ng],'XTickLabel',[kg(1),kg(end/2),kg(end)]);set(gca,'YTick',[1,ng/2,ng],'YTickLabel',[kg(1),kg(end/2),kg(end)]);
subplot(2,2,4);imagesc(imag(Fk2));title('imag(Fk2)'); colorbar;
set(gca,'XTick',[1,ng/2,ng],'XTickLabel',[kg(1),kg(end/2),kg(end)]);set(gca,'YTick',[1,ng/2,ng],'YTickLabel',[kg(1),kg(end/2),kg(end)]);
