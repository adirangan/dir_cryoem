function test_F_4();
% test 2d-fourier transform, centering real-space at origin ; 

ng = 64; cmax = 1;
xg = linspace(0,cmax,ng+1); xg = xg(1:end-1); dx = mean(diff(xg));
kg = (0:ng-1)/cmax; dk = mean(diff(kg));
[Xc,Yc] = meshgrid(xg,xg); Xc(:,ng/2+1:end) = Xc(:,ng/2+1:end)-cmax; Yc(ng/2+1:end,:) = Yc(ng/2+1:end,:)-cmax;
Rc = sqrt(Xc.^2 + Yc.^2); Wc = atan2(Yc,Xc);
sx=0.1;mx=0.3-1*0.5;sy=0.2;my=0.2-1*0.5;tmp_a = 1/(2*pi*cmax*sx*cmax*sy)*exp(-(Xc - cmax*mx).^2/2/(cmax*sx).^2).*exp(-(Yc - cmax*my).^2/2/(cmax*sy).^2);
sx=0.1;mx=0.7-1*0.5;sy=0.2;my=0.2-1*0.5;tmp_b = 1/(2*pi*cmax*sx*cmax*sy)*exp(-(Xc - cmax*mx).^2/2/(cmax*sx).^2).*exp(-(Yc - cmax*my).^2/2/(cmax*sy).^2);
sx=0.1;mx=0.7-1*0.5;sy=0.2;my=0.8-1*0.5;tmp_c = 1/(2*pi*cmax*sx*cmax*sy)*exp(-(Xc - cmax*mx).^2/2/(cmax*sx).^2).*exp(-(Yc - cmax*my).^2/2/(cmax*sy).^2);
sx=0.1;mx=0.3-1*0.5;sy=0.2;my=0.8-1*0.5;tmp_d = 1/(2*pi*cmax*sx*cmax*sy)*exp(-(Xc - cmax*mx).^2/2/(cmax*sx).^2).*exp(-(Yc - cmax*my).^2/2/(cmax*sy).^2);
tmp_crop = Rc<0.25*(1 + 0.5*cos(5*Wc));
Fc1 = tmp_crop.*(tmp_a+tmp_b+tmp_c+tmp_d);
Lc1 = sum(abs(Fc1(:)).^2)*dx*dx;
Fk1 = fft2(Fc1) * dx*dx;
Lk1 = sum(abs(Fk1(:)).^2)*dk*dk;
disp(sprintf(' %% Lc1 %f Lk1 %f',Lc1,Lk1));
figure;
subplot(2,2,1);imagesc(recenter(real(Fc1)));title('real(Fc1)'); colorbar; 
set(gca,'XTick',[1,ng/2+1,ng],'XTickLabel',[xg(1),xg(end/2+1),xg(end)]-cmax/2);set(gca,'YTick',[1,ng/2+1,ng],'YTickLabel',[xg(1),xg(end/2+1),xg(end)]-cmax/2);
%subplot(2,2,2);imagesc(recenter(imag(Fc1)));title('imag(Fc1)'); colorbar;
%set(gca,'XTick',[1,ng/2+1,ng],'XTickLabel',[xg(1),xg(end/2+1),xg(end)]-cmax/2);set(gca,'YTick',[1,ng/2+1,ng],'YTickLabel',[xg(1),xg(end/2+1),xg(end)]-cmax/2);
subplot(2,2,2);imagesc(recenter(log(abs(Fk1))));title('log(abs(Fk1))'); colorbar;
set(gca,'XTick',[1,ng/2+1,ng],'XTickLabel',[kg(1),kg(end/2+1),kg(end)]-ng/2);set(gca,'YTick',[1,ng/2+1,ng],'YTickLabel',[kg(1),kg(end/2+1),kg(end)]-ng/2);
subplot(2,2,3);imagesc(recenter(real(Fk1)));title('real(Fk1)'); colorbar;
set(gca,'XTick',[1,ng/2+1,ng],'XTickLabel',[kg(1),kg(end/2+1),kg(end)]-ng/2);set(gca,'YTick',[1,ng/2+1,ng],'YTickLabel',[kg(1),kg(end/2+1),kg(end)]-ng/2);
subplot(2,2,4);imagesc(recenter(imag(Fk1)));title('imag(Fk1)'); colorbar;
set(gca,'XTick',[1,ng/2+1,ng],'XTickLabel',[kg(1),kg(end/2+1),kg(end)]-ng/2);set(gca,'YTick',[1,ng/2+1,ng],'YTickLabel',[kg(1),kg(end/2+1),kg(end)]-ng/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = recenter(input);
[nrows,ncols] = size(input);
new_rows = [floor(nrows/2)+1:nrows , 1:floor(nrows/2)];
new_cols = [floor(ncols/2)+1:ncols , 1:floor(ncols/2)];
output = input(new_rows,new_cols);
