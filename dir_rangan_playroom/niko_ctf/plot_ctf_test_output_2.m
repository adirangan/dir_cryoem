ctfw1 = load('ctfw1');
ctf = ctfw1(:,1:end-1); % removing a column to ensure asymmetry ;
ctf = [ctf ; ctf(1,:) ] ; % periodizing across theta ;
[ngridc_plus_1,ngridr] = size(ctf); ngridc = ngridc_plus_1 - 1;
[tt,rr] = meshgrid(2*pi*([1:ngridc , 1])/ngridc,1:ngridr); % note that in ctf_test.f the values of xnodesr range from 1 to ngridr, whereas the values of xnodesc range from 1/ngridc to 2*pi ;
xx = rr.*cos(tt);
yy = rr.*sin(tt);
figure(3);clf;
surf(xx,yy,0*xx,transpose(ctf));
shading interp
view(2)

figure(4);clf;
for i=1:50:ngridc
	plot(ctfw1(i,:))
	hold on
end

