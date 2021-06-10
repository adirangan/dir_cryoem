function [p_x_,p_v_] = test_LSQ_vs_Sampling_excerpt_0(p,x__j,y__j,disp_flag);
% fits polynomial of order p to x__j,y__j in interval x_lim, ;
% returns polynomial_coefficients p_x_ for p(x). ;
% We also return polynomial_coefficients p_v_ for p(v), ;
% where v represents the x-variable shifted and scaled to lie in [-1,+1]. ;

if (nargin<3);
p = 3;
x__j = 8*rand(128,1);
y__j = sin(x__j) + randn(size(x__j));
disp_flag=1;
[p_x_,p_v_] = test_LSQ_vs_Sampling_excerpt_0(p,x__j,y__j,disp_flag);
disp('returning'); return;
end;%if (nargin<3);

if (nargin<4); disp_flag=0; end;

x_min = min(x__j);
x_max = max(x__j);
v__j = 2 * (x__j - x_min) / (x_max - x_min) - 1;
a = 2/(x_max - x_min) ;
b = -2*x_min/(x_max - x_min) - 1;
v__j = reshape(v__j,numel(v__j),1);

Lv = orthopoly_vec_0(p,[1]);
L__j = zeros(numel(v__j),size(Lv,1));
for np=1:size(Lv,1);
L__j(:,np) = polyval(Lv(np,:),v__j);
end;%for np=1:size(Lv,1);

Lp = L__j\y__j;
p_v_ = transpose(Lp)*Lv;

p_ab = zeros(p+1,p+1);
p_ab(1,end) = 1;
p_ab(2,p:p+1) = [a,b];
for np=3:p+1;
p_ab(np,p+1-(np-1):p+1) = conv(p_ab(np-1,p+1-(np-2):p+1),[a,b]);
end;%for np=3:p+1;
p_x_ = p_v_(end:-1:1)*p_ab;

if disp_flag;
v__s = linspace(-1,1,1024); p_v__s = polyval(p_v_,v__s);
x__s = linspace(x_min,x_max,1024); p_x__s = polyval(p_x_,x__s);
subplot(1,2,1); plot(v__j,y__j,'o',v__s,p_v__s,'k-'); 
xlim([-1,+1]); xlabel('v'); title('p_v');
subplot(1,2,2); plot(x__j,y__j,'o',x__s,p_x__s,'k-'); 
xlim([x_min,x_max]); xlabel('x'); title('p_x');
end;%if disp_flag;
