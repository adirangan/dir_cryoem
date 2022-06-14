function output = colormap_phase(n_c,gamma_0,gamma_1,gamma_2);
if (nargin<1); n_c = 64; end;
if (nargin<2); gamma_0 = 1.45; end;
if (nargin<3); gamma_1 = 1.35; end;
if (nargin<4); gamma_2 = 0.65; end;

if n_c<0;
figure(1);clf;
imagesc([[0:64-1],[0:64-1],[0:64-1]]);
colormap(colormap_phase());
axisnotick;
disp('returning'); return;
end;%if n_c<0;

phase_0_ = 0.5*(1 + cos(2*pi*[0:n_c-1]/n_c));
phase_1_ = circshift(phase_0_,floor(1*n_c/3));
phase_2_ = circshift(phase_0_,floor(2*n_c/3));

c__ = zeros(n_c,3);
c__(:,1+0) = phase_0_.^gamma_0;
c__(:,1+1) = phase_1_.^gamma_1;
c__(:,1+2) = phase_2_.^gamma_2;
output = c__;
