% testing bessel-function properties ; finding K required for various error-bounds epsilon ; 
A = @(K,Z) 1./sqrt(2*pi*K) .* (exp(1) .* Z ./ K / 2).^K *(2*pi); 

eps_ = 0.1.^[1:4];

for np=1:length(eps_);
subplot(2,2,np);
eps = eps_(np);
z_ = 0:10; k_ = 1:16;
[K,Z] = meshgrid(k_,z_);
J = A(K,Z);
cra = colormap('autumn'); cra(1,:) = [1,1,1]; colormap(cra); 
imagesc(log(J),[log10(eps),1]); 
set(gca,'Xtick',k_,'XtickLabel',k_,'Ytick',z_,'YtickLabel',z_,'Ydir','normal');
xlabel('frequency n'); ylabel('z = k\delta'); title(sprintf('error = %0.4f',eps));
end;%for np=1:length(eps_);
