% test 2d-fourier transform ; 

ng = 16; cmax = 1;
xg = linspace(0,cmax,ng+1); xg = xg(1:end-1); dx = mean(diff(xg));
kg = (0:ng-1)/cmax; dk = mean(diff(kg));
[Xc,Yc] = meshgrid(xg,xg);
Fc = exp(i*2*pi*1*Xc/cmax).*exp(i*2*pi*2*Yc/cmax);
Lc = sum(abs(Fc(:)).^2)*dx*dx;
Fk = fft2(Fc) * dx*dx;
Lk = sum(abs(Fk(:)).^2)*dk*dk;
disp(sprintf(' %% Lc %f Lk %f',Lc,Lk));
subplot(2,2,1);imagesc(real(Fc));title('real(Fc)'); colorbar; 
set(gca,'XTick',[1,ng/2,ng],'XTickLabel',[xg(1),xg(end/2),xg(end)]);set(gca,'YTick',[1,ng/2,ng],'YTickLabel',[xg(1),xg(end/2),xg(end)]);
subplot(2,2,2);imagesc(imag(Fc));title('imag(Fc)'); colorbar;
set(gca,'XTick',[1,ng/2,ng],'XTickLabel',[xg(1),xg(end/2),xg(end)]);set(gca,'YTick',[1,ng/2,ng],'YTickLabel',[xg(1),xg(end/2),xg(end)]);
subplot(2,2,3);imagesc(real(Fk));title('real(Fk)'); colorbar;
set(gca,'XTick',[1,ng/2,ng],'XTickLabel',[kg(1),kg(end/2),kg(end)]);set(gca,'YTick',[1,ng/2,ng],'YTickLabel',[kg(1),kg(end/2),kg(end)]);
subplot(2,2,4);imagesc(imag(Fk));title('imag(Fk)'); colorbar;
set(gca,'XTick',[1,ng/2,ng],'XTickLabel',[kg(1),kg(end/2),kg(end)]);set(gca,'YTick',[1,ng/2,ng],'YTickLabel',[kg(1),kg(end/2),kg(end)]);
