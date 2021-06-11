ctfw = load('ctfw');
ngridr = size(ctfw,2);
ngridc = size(ctfw,1);
ctfw = reshape(ctfw,[ngridc,ngridr]);
[rr,tt] = meshgrid(1:ngridc,1:ngridr);
xx = rr.*cos(tt);
yy = rr.*sin(tt);
figure(1);clf;
surf(xx,yy,0*xx,ctfw)
shading interp
view(2)

figure(2);clf;
for i=1:50:ngridc
	plot(ctfw(i,:))
	hold on
end

ctfw1 = load('ctfw1');
ngridr = size(ctfw1,2);
ngridc = size(ctfw1,1);
ctfw1 = reshape(ctfw1,[ngridc,ngridr]);
[rr,tt] = meshgrid(1:ngridc,1:ngridr);
xx = rr.*cos(tt);
yy = rr.*sin(tt);
figure(3);clf;
surf(xx,yy,0*xx,ctfw1)
shading interp
view(2)

figure(4);clf;
for i=1:50:ngridc
	plot(ctfw1(i,:))
	hold on
end

