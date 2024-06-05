function ...
[ ...
 scatter_from_tensor_sba__ ...
] = ...
shell_k_p_scatter_from_tensor_interpolate_n_6( ...
 n_order ...
,n_azimu_b ...
,n_polar_a ...
,polar_a_ ...
,n_scatter ...
,azimu_b_scatter_ ...
,polar_a_scatter_ ...
,flag_polar_a_ascend_vs_descend ...
);
%%%%%%%%;
% Here we attempt 'reflecting' the interpolation grid (i.e., local nodes) across the poles. ;
% This does not work particularly well, since the reflection is not smoothly connected to the original. ;
%%%%%%%%;
% returns sparse matrix encoding the n_order interpolation operator between a tensor_grid ;
% using (even) n_azimu_b azimu_b (uniformly from 0 to 2*pi with periodic boundary). ;
% and n_polar_a polar_a values (not necessarily uniform from pi to 0) ;
% Note that (assuming flag_polar_a_ascend_vs_descend==0) polar_a_ descends. ;
%%%%%%%%;
% Note that a shell discretization (e.g., using sample_shell_6) ;
% with str_T_vs_L=='T' will generate uniformly spaced polar_a_, ;
% however, this discretization will not include the endpoints (i.e., 0 and pi). ;
% A similar discretization with str_T_vs_L=='C' will produce ;
% polar_a_ values extending from 0 to pi. ;
%%%%;
% We presume these points are associated with a (n_azimu_b)-by-(n_polar_a) matrix ;
% a_k_p_ba__ with rows corresponding to azimu_b_ and columns corresponding to polar_a_. ;
%%%%;
% Note that this particular ordering is intended to align with sample_shell_6. ;
% i.e.: (assuming flag_polar_a_ascend_vs_descend==0) ;
% [lgnd_node_,lgnd_weight_] = legpts(n_polar_a);
% polar_a_ = acos(lgnd_node_);
%%%%;
% The n_scatter points to be interpolated have coordinates stored in ;
% azimu_b_scatter_ and polar_a_scatter_. ;
% We presume these points are associated with a n_scatter-element list a_k_all_. ;
% ;
% The (n_scatter)-by-(n_azimu_b*n_polar_a) scatter_from_tensor_sba__ matrix stores the ;
% interpolation weights linking a_k_all_ to the unrolled a_k_p_ba__(:). ;
%%%%;

str_thisfunction = 'shell_k_p_scatter_from_tensor_interpolate_n_6';

if (nargin<6);
rng(1);
%%%%%%%%;
flag_check=1;
if flag_check;
flag_disp = 1; nf=0;
disp(sprintf(' %% testing %s',str_thisfunction));
for flag_polar_a_ascend_vs_descend=[0,1];
disp(sprintf(' %% flag_polar_a_ascend_vs_descend %d:',flag_polar_a_ascend_vs_descend));
str_T_vs_L = 'C'; %<-- produces polar_a_ values that extend all the way from 0 to pi. ;
[ ...
 n_all ...
,azimu_b_all_ ...
,polar_a_all_ ...
,weight_all_ ...
,k_c_0_all_ ...
,k_c_1_all_ ...
,k_c_2_all_ ...
,n_polar_a ...
,polar_a_ ...
,n_azimu_b_ ...
] = ...
sample_shell_6( ...
 1.0 ...
,1.0/(2*pi)/2 ...
,str_T_vs_L ...
,1 ...
) ;
n_azimu_b = n_azimu_b_(1+0); azimu_b_ = transpose(linspace(0,2*pi,n_azimu_b+1)); azimu_b_ = azimu_b_(1:end-1);
if flag_polar_a_ascend_vs_descend==1; polar_a_ = flipud(polar_a_); end; %<-- traditional order. ;
if flag_polar_a_ascend_vs_descend==0; end; %<-- reverse order. ;
disp(sprintf(' %% n_azimu_b %d n_polar_a %d polar_a_: [%+0.2f,%+0.2f]',n_azimu_b,n_polar_a,polar_a_(1),polar_a_(end)));
%%%%%%%%;
n_grid = n_azimu_b*n_polar_a;
[azimu_b_ba__,polar_a_ba__] = ndgrid(azimu_b_,polar_a_);
%%%%;
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figbig;
subplot(1,2,1);
plot3(k_c_0_all_,k_c_1_all_,k_c_2_all_,'k.');
axis equal; axis vis3d;
subplot(1,2,2);
plot(sort(polar_a_),'o');
ylabel('polar_a_','Interpreter','none');
ylim([0,pi]);
disp('returning'); return;
end;%if (flag_disp>1);
%%%%%%%%;
n_scat = 1024*4;
polar_a_scat_ = 3*pi*rand(n_scat,1) - 1*pi;
azimu_b_scat_ = 5*pi*rand(n_scat,1) - 2*pi;
l_max = 49;
[Ylm_grid__] = get_Ylm__(1+l_max,0:l_max,n_grid,azimu_b_ba__,polar_a_ba__);
[Ylm_scat__] = get_Ylm__(1+l_max,0:l_max,n_scat,azimu_b_scat_,polar_a_scat_);
%%%%;
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figbig;fig80s;
l_val = floor(l_max/2); m_val = max(-l_val,min(+l_val,floor(l_val*2/3)));
tmp_Ylm_ab__ = reshape(Ylm_grid__{1+l_val}(1+l_val+m_val,:),[n_azimu_b,n_polar_a]);
subplot(1,2,1);
imagesc(transpose(real(tmp_Ylm_ab__)));
axis image; axisnotick;
xlabel('azimu_b_','Interpreter','none');
ylabel('polar_a_','Interpreter','none');
subplot(1,2,2);
imagesc(transpose(imag(tmp_Ylm_ab__)));
axis image; axisnotick;
xlabel('azimu_b_','Interpreter','none');
ylabel('polar_a_','Interpreter','none');
disp('returning'); return;
end;%if (flag_disp>1);
%%%%%%%%;
for n_order=[1,3,5,7,9];
scatter_from_tensor_sba__ = ...
shell_k_p_scatter_from_tensor_interpolate_n_6( ...
 n_order ...
,n_azimu_b ...
,n_polar_a ...
,polar_a_ ...
,n_scat ...
,azimu_b_scat_ ...
,polar_a_scat_ ...
,flag_polar_a_ascend_vs_descend ...
);
%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;
slim_ = max(abs(scatter_from_tensor_sba__),[],'all');
tmp_index_ = efind(scatter_from_tensor_sba__);
plot(sort(scatter_from_tensor_sba__(1+tmp_index_),'descend'),'.');
xlabel('nodes sorted by weight'); ylabel('weight');
title(sprintf('n_order %d',n_order));
end;%if flag_disp;
%%%%;
Ylm_pint__ = cell(1+l_max,1);
n_l_sub = 4; l_sub_ = max(0,min(l_max,round(linspace(1,l_max,n_l_sub))));
for l_val=l_sub_;
n_m_sub = 4; m_sub_ = max(-l_val,min(+l_val,round(linspace(-l_val,+l_val,n_m_sub))));
for m_val=m_sub_;
Ylm_pint__{1+l_val}(1+l_val+m_val,:) = scatter_from_tensor_sba__*transpose(Ylm_grid__{1+l_val}(1+l_val+m_val,:));
end;%for m_val=m_sub_;
end;%for l_val=l_sub_;
disp(sprintf(' %% interpolation relative error: n_order %d ',n_order));
for l_val=l_sub_;
n_m_sub = 4; m_sub_ = max(-l_val,min(+l_val,round(linspace(-l_val,+l_val,n_m_sub))));
for m_val=m_sub_;
disp(sprintf(' %% l_val %d m_val %+d : %0.16f',l_val,m_val,fnorm(Ylm_scat__{1+l_val}(1+l_val+m_val,:) - Ylm_pint__{1+l_val}(1+l_val+m_val,:))/max(1e-12,fnorm(Ylm_scat__{1+l_val}(1+l_val+m_val,:)))));
end;%for m_val=m_sub_;
end;%for l_val=l_sub_;
end;%for n_order=[1,3,5,7,9];
end;%for flag_polar_a_ascend_vs_descend=[0,1];
end;%if flag_check;
%%%%%%%%;
disp(sprintf(' %% returning')); return;
end;%if (nargin<6);

flag_verbose = 1;
flag_disp = 1; nf=0;
flag_check = 1;

if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

na=0;
if (nargin<1+na); n_order=[]; end; na=na+1;
if (nargin<1+na); n_azimu_b=[]; end; na=na+1;
if (nargin<1+na); n_polar_a=[]; end; na=na+1;
if (nargin<1+na); polar_a_=[]; end; na=na+1;
if (nargin<1+na); n_scatter=[]; end; na=na+1;
if (nargin<1+na); azimu_b_scatter_=[]; end; na=na+1;
if (nargin<1+na); polar_a_scatter_=[]; end; na=na+1;
if (nargin<1+na); flag_polar_a_ascend_vs_descend=[]; end; na=na+1;
if isempty(flag_polar_a_ascend_vs_descend); flag_polar_a_ascend_vs_descend = 0; end;

if abs(min(polar_a_)-0)> 1e-12; disp(sprintf(' %% Warning, expecting polar_a_ to begin at 0 in %s (try str_T_vs_L = ''C'' when generating grid) ',str_thisfunction)); end;
if abs(max(polar_a_)-pi)> 1e-12; disp(sprintf(' %% Warning, expecting polar_a_ to end at pi in %s (try str_T_vs_L = ''C'' when generating grid) ',str_thisfunction)); end;
if var(diff(polar_a_))> 1e-12; disp(sprintf(' %% Warning, expecting polar_a_ to be equispaced in %s (try str_T_vs_L = ''C'' when generating grid) ',str_thisfunction)); end;
if mod(n_azimu_b,2)==1; disp(sprintf(' %% Warning, expecting n_azimu_b to be even in %s',str_thisfunction)); end;
if mod(n_order,2)==0; disp(sprintf(' %% Warning, expecting n_order to be odd in %s',str_thisfunction)); end;

if (n_order>min(n_azimu_b,n_polar_a)); disp(sprintf(' %% Warning, n_order %d > n_azimu_b %d n_polar_a %d',n_order,n_azimu_b,n_polar_a)); end;
n_order = min(n_order,min(n_azimu_b,n_polar_a));

azimu_b_ = linspace(0,2*pi,n_azimu_b+1); azimu_b_ = transpose(azimu_b_(1:n_azimu_b));
[azimu_ba__,polar_ba__] = ndgrid(azimu_b_,polar_a_); n_grid = n_azimu_b*n_polar_a;

%%%%%%%%;
% periodize azimu_b_scatter_ (standard). ;
%%%%%%%%;
azimu_b_scatter_periodize_ = periodize(azimu_b_scatter_,0,2*pi);
%%%%%%%%;
% periodize polar_a_scatter_ (note behavior at poles). ;
%%%%%%%%;
polar_a_scatter_periodize_ = periodize(polar_a_scatter_,-1*pi/2,+3*pi/2);
tmp_index_ = efind(polar_a_scatter_periodize_< 0*pi);
polar_a_scatter_periodize_(1+tmp_index_) = 0*pi - (polar_a_scatter_periodize_(1+tmp_index_) - 0*pi);
azimu_b_scatter_periodize_(1+tmp_index_) = +azimu_b_scatter_periodize_(1+tmp_index_) + 1*pi;
tmp_index_ = efind(polar_a_scatter_periodize_> 1*pi);
polar_a_scatter_periodize_(1+tmp_index_) = 1*pi - (polar_a_scatter_periodize_(1+tmp_index_) - 1*pi);
azimu_b_scatter_periodize_(1+tmp_index_) = +azimu_b_scatter_periodize_(1+tmp_index_) - 1*pi;
%%%%%%%%;
% periodize azimu_b_scatter_ (standard). ;
%%%%%%%%;
azimu_b_scatter_periodize_ = periodize(azimu_b_scatter_periodize_,0,2*pi);

%%%%%%%%;
% flip polar_a_ if necessary. ;
%%%%%%%%;
sign_polar_a = +1;
if flag_polar_a_ascend_vs_descend==0;
sign_polar_a = -1; polar_a_scatter_periodize_ = pi-polar_a_scatter_periodize_; %<-- flip polar_a_scatter_periodize_ to account for reverse ordering of polar_a_. ;
end;%if flag_polar_a_ascend_vs_descend==0;

%%%%;
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figmed;
linewidth_big = 2.5;
linewidth_use = 0.5;
markersize_use = 1;
subplot(1,2,1);
hold on;
plot([0;2;2;0;0]*pi,[0;0;1;1;0]*pi,'m-','LineWidth',linewidth_big);
plot(permute(azimu_ba__,[1,2]),permute(polar_ba__,[1,2]),'k-','LineWidth',linewidth_use);
plot(permute(azimu_ba__,[2,1]),permute(polar_ba__,[2,1]),'k-','LineWidth',linewidth_use);
plot(azimu_b_scatter_,polar_a_scatter_,'ro','MarkerFaceColor','g','MarkerSize',markersize_use);
hold off;
%xlim([0,2*pi]); ylim([0,1*pi]);
xlabel('azimu_b_','Interpreter','none');
ylabel('polar_a_','Interpreter','none');
title('input scatter');
subplot(1,2,2);
hold on;
plot([0;2;2;0;0]*pi,[0;0;1;1;0]*pi,'m-','LineWidth',linewidth_big);
plot(permute(azimu_ba__,[1,2]),permute(polar_ba__,[1,2]),'k-','LineWidth',linewidth_use);
plot(permute(azimu_ba__,[2,1]),permute(polar_ba__,[2,1]),'k-','LineWidth',linewidth_use);
plot(azimu_b_scatter_periodize_,polar_a_scatter_periodize_,'ro','MarkerFaceColor','g','MarkerSize',markersize_use);
hold off;
%xlim([0,2*pi]); ylim([0,1*pi]);
xlabel('azimu_b_','Interpreter','none');
ylabel('polar_a_','Interpreter','none');
title('periodize scatter');
disp('returning');return;
end;%if (flag_disp>1);
%%%%;

%%%%;
if flag_check;
f_check = @(a,b) ...
  ( 0.0 - 0.5*sin(1*a) + 0.25*sin(2*a) - 0.125*sin(5*a) + 0.0625*sin(7*a) ) ...
.*( 1.0*cos(1*b) - 0.5*cos(2*b) + 0.25*cos(3*b) - 0.125*cos(4*b) + 0.0625*cos(5*b) ) ...
;
f_check_ba__ = f_check(polar_ba__,azimu_ba__);
f_check_scatter_ = f_check(polar_a_scatter_,azimu_b_scatter_);
f_check_scatter_periodize_ = f_check(polar_a_scatter_periodize_,azimu_b_scatter_periodize_);
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figmed;fig80s;
subplot(1,2,1); imagesc(transpose(f_check_ba__)); colorbar;
axis image; axisnotick; title('f_check','Interpreter','none');
subplot(1,2,2); plot(f_check_scatter_,f_check_scatter_periodize_,'k.');
xlabel('f_check_scatter_','Interpreter','none');
ylabel('f_check_scatter_periodize_','Interpreter','none');
disp('returning'); return;
end;%if (flag_disp>1);
end;%if flag_check;
%%%%;

%%%%%%%%;
% transform polar_a_ to be uniform. ;
%%%%%%%%;
if flag_polar_a_ascend_vs_descend==1; polar_a_uni_ = transpose(0:+1:n_polar_a-1); end; %<-- uniform transformation. ;
if flag_polar_a_ascend_vs_descend==0; polar_a_uni_ = transpose(n_polar_a-1:-1:0); end; %<-- uniform transformation. ;
polar_a_scatter_periodize_rescale_ = interp1(polar_a_(:),polar_a_uni_(:),polar_a_scatter_periodize_,'linear','extrap'); %<-- use interpolation to uniformize polar_a_. ;

if (flag_verbose>0);
disp(sprintf(' %% polar_a_scatter_periodize_ in [%0.6f , %0.6f]',prctile(polar_a_scatter_periodize_,[  0,100])));
disp(sprintf(' %% n_polar_a-1 = %d',n_polar_a-1));
disp(sprintf(' %% polar_a_scatter_periodize_rescale_ in [%0.6f , %0.6f]',prctile(polar_a_scatter_periodize_rescale_,[  0,100])));
end;%if (flag_verbose>0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% calculate scatter_from_tensor_sba__. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
spacing_azimu_b = 2*pi/max(1,n_azimu_b);
spacing_polar_a = 1.0; %<-- uniformized. ;
azimu_b_scatter_rescale_ = periodize(azimu_b_scatter_periodize_,0,2*pi)/spacing_azimu_b;
polar_a_scatter_rescale_ = max(0-1,min(n_polar_a-1+1,polar_a_scatter_periodize_rescale_/spacing_polar_a)); %<-- confine to box if not already confined. ;
%%%%;
if (flag_verbose>0);
disp(sprintf(' %% n_azimu_b-1 = %d',n_azimu_b-1));
disp(sprintf(' %% azimu_b_scatter_rescale_ in [%0.6f , %0.6f]',prctile(azimu_b_scatter_rescale_,[  0,100])));
disp(sprintf(' %% n_polar_a-1 = %d',n_polar_a-1));
disp(sprintf(' %% polar_a_scatter_periodize_rescale_ in [%0.6f , %0.6f]',prctile(polar_a_scatter_periodize_rescale_,[  0,100])));
end;%if (flag_verbose>0);
%%%%;
n_half_order = floor(n_order/2);
n_part_order = round((n_order-1)/2);
azimu_b_scatter_start_ = floor(azimu_b_scatter_rescale_)-n_part_order;
azimu_b_scatter_final_ = azimu_b_scatter_start_ + n_part_order + 1;
polar_a_scatter_start_ = floor(polar_a_scatter_rescale_)-n_part_order;
polar_a_scatter_final_ = polar_a_scatter_start_ + n_part_order + 1;
node_x_ = transpose(-n_part_order-0.5:1:+n_part_order+0.5);
if (flag_verbose>0); disp(sprintf(' %% node_x_: %s',num2str(transpose(node_x_),'%+0.2f '))); end;
n_nodex = numel(node_x_);
weight_denominator_x_ = prod( repmat(node_x_,1,n_nodex) - repmat(transpose(node_x_),n_nodex,1) + eye(n_nodex,n_nodex) , 2 ) ;
if (flag_verbose>0); disp(sprintf(' %% weight_denominator_x_: %s',num2str(transpose(weight_denominator_x_),'%+0.2f '))); end;
if (flag_verbose>0);
for tmp_x = [-1.5,-0.5,0,+0.5,+1.75];
disp(sprintf(' %% tmp_x: %0.6f',tmp_x));
tmp_d_ = prod( tmp_x - repmat(node_x_,1,n_nodex) + diag(1-(tmp_x-node_x_)) , 1 );
for nnodex=0:n_nodex-1;
tmp_d = prod(tmp_x - [node_x_(1 + [0:nnodex-1]);node_x_(1+[nnodex+1:n_nodex-1])]);
disp(sprintf(' %% nnodex %d/%d: %+0.4f vs %+0.4f',nnodex,n_nodex,tmp_d,tmp_d_(1+nnodex)));
end;%for nnodex=0:n_nodex-1;
end;%for tmp_x = [-1.5,-0.5,0,+0.5,+1.75];
end;%if (flag_verbose>0);
azimu_b_scatter_shift_ = azimu_b_scatter_rescale_ - azimu_b_scatter_start_ - n_part_order - 0.5 ;
polar_a_scatter_shift_ = polar_a_scatter_rescale_ - polar_a_scatter_start_ - n_part_order - 0.5 ;
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,2,1);
plot(1:n_scatter,azimu_b_scatter_shift_,'.'); title('azimu_b_scatter_shift_','Interpreter','none');
subplot(1,2,2);
plot(1:n_scatter,polar_a_scatter_shift_,'.'); title('polar_a_scatter_shift_','Interpreter','none');
end;%if (flag_disp>1);
%%%%;
weight_numerator_azimu_b_xs__ = zeros(n_nodex,n_scatter);
weight_numerator_polar_a_xs__ = zeros(n_nodex,n_scatter);
index_col_azimu_b_xs__ = zeros(n_nodex,n_scatter);
index_col_polar_a_xs__ = zeros(n_nodex,n_scatter);
index_row_xxs__ = repmat(0:n_scatter-1,[n_nodex^2,1]);
for nscatter=0:n_scatter-1;
weight_numerator_azimu_b_xs__(:,1+nscatter) = prod( azimu_b_scatter_shift_(1+nscatter) - repmat(node_x_,1,n_nodex) + diag(1-(azimu_b_scatter_shift_(1+nscatter)-node_x_)) , 1 ) ;
weight_numerator_polar_a_xs__(:,1+nscatter) = prod( polar_a_scatter_shift_(1+nscatter) - repmat(node_x_,1,n_nodex) + diag(1-(polar_a_scatter_shift_(1+nscatter)-node_x_)) , 1 ) ;
end;%for nscatter=0:n_scatter-1;
%%%%%%%%;
index_col_azimu_b_xs__ = bsxfun(@plus,reshape(azimu_b_scatter_start_,[1,n_scatter]),transpose([0:n_nodex-1]));
flag_flip_azimu_b_xs__ = zeros(n_nodex,n_scatter);
index_col_polar_a_xs__ = bsxfun(@plus,reshape(polar_a_scatter_start_,[1,n_scatter]),transpose([0:n_nodex-1]));
%%%%;
if (flag_verbose>0); disp(sprintf(' %% n_azimu_b %d: index_col_azimu_b_xs__ in [%.4d,%.4d]',n_azimu_b,min(index_col_azimu_b_xs__,[],'all'),max(index_col_azimu_b_xs__,[],'all'))); end;
if (flag_verbose>0); disp(sprintf(' %% n_polar_a %d: index_col_polar_a_xs__ in [%.4d,%.4d]',n_polar_a,min(index_col_polar_a_xs__,[],'all'),max(index_col_polar_a_xs__,[],'all'))); end;
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figsml;
plot(index_col_azimu_b_xs__(:),index_col_polar_a_xs__(:),'.');
xlabel('index_col_azimu_b_xs__','Interpreter','none');
ylabel('index_col_polar_a_xs__','Interpreter','none');
end;%if (flag_disp>1);

%%%%;
if (flag_verbose>1);
%%%%;
[nc_min,tmp_ij] = min(index_col_polar_a_xs__(:)); tmp_index = tmp_ij-1;
tmp_index_x = mod(tmp_index,n_nodex);
tmp_index_s = (tmp_index-tmp_index_x)/n_nodex;
disp(sprintf(' %% tmp_index %d (%d,%d) <-- nc_min %d',tmp_index,tmp_index_x,tmp_index_s,nc_min));
disp(sprintf(' %% index_col_polar_a_xs__(:,1+tmp_index_s): %s',num2str(transpose(index_col_polar_a_xs__(:,1+tmp_index_s)),' %+d')));
disp(sprintf(' %% index_col_azimu_b_xs__(:,1+tmp_index_s): %s',num2str(transpose(index_col_azimu_b_xs__(:,1+tmp_index_s)),' %+d')));
%%%%;
[nc_max,tmp_ij] = max(index_col_polar_a_xs__(:)); tmp_index = tmp_ij-1;
tmp_index_x = mod(tmp_index,n_nodex);
tmp_index_s = (tmp_index-tmp_index_x)/n_nodex;
disp(sprintf(' %% tmp_index %d (%d,%d) <-- nc_max %d',tmp_index,tmp_index_x,tmp_index_s,nc_max));
disp(sprintf(' %% index_col_polar_a_xs__(:,1+tmp_index_s): %s',num2str(transpose(index_col_polar_a_xs__(:,1+tmp_index_s)),' %+d')));
disp(sprintf(' %% index_col_azimu_b_xs__(:,1+tmp_index_s): %s',num2str(transpose(index_col_azimu_b_xs__(:,1+tmp_index_s)),' %+d')));
%%%%;
disp('returning'); return;
end;%if (flag_verbose>1);
%%%%;

%%%%%%%%;
if (flag_verbose>0);
[nc_min,tmp_ij] = min(index_col_polar_a_xs__(:)); tmp_index = tmp_ij-1;
tmp_index_x = mod(tmp_index,n_nodex);
tmp_index_s = (tmp_index-tmp_index_x)/n_nodex;
nscatter = tmp_index_s;
%nscatter = 0;
disp(sprintf(' %% tmp_index %d (%d,%d) <-- nc_min %d',tmp_index,tmp_index_x,tmp_index_s,nc_min));
%%%%;
tmp_weight_azimu_b_x_ = diag(1./weight_denominator_x_)*weight_numerator_azimu_b_xs__(:,1+nscatter);
tmp_weight_polar_a_x_ = diag(1./weight_denominator_x_)*weight_numerator_polar_a_xs__(:,1+nscatter);
tmp_index_col_azimu_b_x_ = index_col_azimu_b_xs__(:,1+nscatter);
tmp_index_col_polar_a_x_ = index_col_polar_a_xs__(:,1+nscatter);
tmp_flag_flip_azimu_b_x_ = zeros(n_nodex,1);
%%%%;
if (flag_verbose>0); disp(sprintf(' %% tmp_index_col_azimu_b_x_: %s',num2str(transpose(tmp_index_col_azimu_b_x_),' %+.2d'))); end;
if (flag_verbose>0); disp(sprintf(' %% tmp_index_col_polar_a_x_: %s',num2str(transpose(tmp_index_col_polar_a_x_),' %+.2d'))); end;
if (flag_verbose>0); disp(sprintf(' %% tmp_flag_flip_azimu_b_x_: %s',num2str(transpose(tmp_flag_flip_azimu_b_x_),' %+.2d'))); end;
if (flag_verbose>0); disp(sprintf(' %% n_azimu_b %d: tmp_index_col_azimu_b_x_ in [%.4d,%.4d]',n_azimu_b,min(tmp_index_col_azimu_b_x_,[],'all'),max(tmp_index_col_azimu_b_x_,[],'all'))); end;
if (flag_verbose>0); disp(sprintf(' %% n_polar_a %d: tmp_index_col_polar_a_x_ in [%.4d,%.4d]',n_polar_a,min(tmp_index_col_polar_a_x_,[],'all'),max(tmp_index_col_polar_a_x_,[],'all'))); end;
%%%%
tmp_azimu_b_pre_x_ = 2*pi*(tmp_index_col_azimu_b_x_)/max(1,n_azimu_b);
tmp_polar_a_pre_x_ = 1*pi*(tmp_index_col_polar_a_x_)/max(1,n_polar_a);
%%%%;
tmp_index_ = efind(tmp_index_col_azimu_b_x_>=n_azimu_b); 
tmp_index_col_azimu_b_x_(1+tmp_index_) = tmp_index_col_azimu_b_x_(1+tmp_index_) - n_azimu_b;
tmp_index_ = efind(tmp_index_col_azimu_b_x_< 0        ); 
tmp_index_col_azimu_b_x_(1+tmp_index_) = tmp_index_col_azimu_b_x_(1+tmp_index_) + n_azimu_b;
%%%%;
tmp_index_ = efind(tmp_index_col_polar_a_x_< 0        ); 
tmp_index_col_polar_a_x_(1+tmp_index_) = 0*n_polar_a - (tmp_index_col_polar_a_x_(1+tmp_index_) - 0*n_polar_a);
tmp_flag_flip_azimu_b_x_(1+tmp_index_) = 1-tmp_flag_flip_azimu_b_x_(1+tmp_index_);
tmp_index_ = efind(tmp_index_col_polar_a_x_> (n_polar_a-1)); 
tmp_index_col_polar_a_x_(1+tmp_index_) = 1*(n_polar_a-1) - (tmp_index_col_polar_a_x_(1+tmp_index_) - 1*(n_polar_a-1));
tmp_flag_flip_azimu_b_x_(1+tmp_index_) = 1-tmp_flag_flip_azimu_b_x_(1+tmp_index_);
%%%%;
tmp_index_ = efind(tmp_index_col_azimu_b_x_>=n_azimu_b); 
tmp_index_col_azimu_b_x_(1+tmp_index_) = tmp_index_col_azimu_b_x_(1+tmp_index_) - n_azimu_b;
tmp_index_ = efind(tmp_index_col_azimu_b_x_< 0        ); 
tmp_index_col_azimu_b_x_(1+tmp_index_) = tmp_index_col_azimu_b_x_(1+tmp_index_) + n_azimu_b;
%%%%;
tmp_azimu_b_pos_x_ = 2*pi*(tmp_index_col_azimu_b_x_)/max(1,n_azimu_b);
tmp_polar_a_pos_x_ = 1*pi*(tmp_index_col_polar_a_x_)/max(1,n_polar_a-1);
if (flag_verbose>0); disp(sprintf(' %% tmp_azimu_b_pos_x_ vs azimu_b_(1+tmp_index_col_azimu_b_x_): %0.16f',fnorm(tmp_azimu_b_pos_x_ - azimu_b_(1+tmp_index_col_azimu_b_x_)))); end;
if (flag_verbose>0); disp(sprintf(' %% tmp_polar_a_pos_x_ vs polar_a_(1+tmp_index_col_polar_a_x_): %0.16f',fnorm(tmp_polar_a_pos_x_ - polar_a_(1+tmp_index_col_polar_a_x_)))); end;
%%%%;
if (flag_verbose>0); disp(sprintf(' %% tmp_index_col_azimu_b_x_: %s',num2str(transpose(tmp_index_col_azimu_b_x_),' %+.2d'))); end;
if (flag_verbose>0); disp(sprintf(' %% tmp_index_col_polar_a_x_: %s',num2str(transpose(tmp_index_col_polar_a_x_),' %+.2d'))); end;
if (flag_verbose>0); disp(sprintf(' %% tmp_flag_flip_azimu_b_x_: %s',num2str(transpose(tmp_flag_flip_azimu_b_x_),' %+.2d'))); end;
if (flag_verbose>0); disp(sprintf(' %% n_azimu_b %d: tmp_index_col_azimu_b_x_ in [%.4d,%.4d]',n_azimu_b,min(tmp_index_col_azimu_b_x_,[],'all'),max(tmp_index_col_azimu_b_x_,[],'all'))); end;
if (flag_verbose>0); disp(sprintf(' %% n_polar_a %d: tmp_index_col_polar_a_x_ in [%.4d,%.4d]',n_polar_a,min(tmp_index_col_polar_a_x_,[],'all'),max(tmp_index_col_polar_a_x_,[],'all'))); end;
%%%%;
tmp_weight_xx__ = tmp_weight_azimu_b_x_*transpose(tmp_weight_polar_a_x_);
tmp_azimu_b_pre_xx__ = repmat(tmp_azimu_b_pre_x_,[1,n_nodex]);
tmp_polar_a_pre_xx__ = repmat(transpose(tmp_polar_a_pre_x_),[n_nodex,1]);
tmp_azimu_b_pos_xx__ = repmat(tmp_azimu_b_pos_x_,[1,n_nodex]);
tmp_polar_a_pos_xx__ = repmat(transpose(tmp_polar_a_pos_x_),[n_nodex,1]);
tmp_sgn_azimu_b_xx__ = repmat(transpose(tmp_flag_flip_azimu_b_x_),[n_nodex,1]);
tmp_azimu_b_pos_xx__ = periodize(tmp_azimu_b_pos_xx__ + tmp_sgn_azimu_b_xx__*pi,0,2*pi);
tmp_f_check_pre_xx__ = f_check(tmp_polar_a_pre_xx__,tmp_azimu_b_pre_xx__);
tmp_f_check_pos_xx__ = f_check(tmp_polar_a_pos_xx__,tmp_azimu_b_pos_xx__);
%%%%;
tmp_col_azimu_b_xx__ = repmat(tmp_index_col_azimu_b_x_,[1,n_nodex]);
tmp_col_polar_a_xx__ = repmat(transpose(tmp_index_col_polar_a_x_),[n_nodex,1]);
tmp_sgn_azimu_b_xx__ = repmat(transpose(tmp_flag_flip_azimu_b_x_),[n_nodex,1]);
tmp_col_azimu_b_pos_xx__ = periodize(tmp_col_azimu_b_xx__ + tmp_sgn_azimu_b_xx__*round(0.5*n_azimu_b),0,n_azimu_b);
tmp_col_polar_a_pos_xx__ = tmp_col_polar_a_xx__;
tmp_col_xx__ = tmp_col_azimu_b_pos_xx__ + tmp_col_polar_a_pos_xx__*n_azimu_b;
tmp_col_xx_ = reshape(tmp_col_xx__,[n_nodex^2,1]);
tmp_f_check_xx__ = reshape(f_check_ba__(1+tmp_col_xx_),[n_nodex,n_nodex]);
%%%%;
if (flag_disp>0);
figure(1+nf);nf=nf+1;clf;figbig; fig80s;
p_row = 3; p_col = 4; np=0;
flim_ = [min([tmp_f_check_pre_xx__(:);tmp_f_check_pos_xx__(:);tmp_f_check_xx__(:)]),max([tmp_f_check_pre_xx__(:);tmp_f_check_pos_xx__(:);tmp_f_check_xx__(:)])];
subplot(p_row,p_col,1+np);np=np+1; imagesc(tmp_col_azimu_b_xx__); axis image; axisnotick; xlabel('b'); ylabel('a'); title('tmp_col_azimu_b_xx__','Interpreter','none'); colorbar;
subplot(p_row,p_col,1+np);np=np+1; imagesc(tmp_col_polar_a_xx__); axis image; axisnotick; xlabel('b'); ylabel('a'); title('tmp_col_polar_a_xx__','Interpreter','none'); colorbar;
subplot(p_row,p_col,1+np);np=np+1; imagesc(tmp_col_azimu_b_pos_xx__); axis image; axisnotick; xlabel('b'); ylabel('a'); title('tmp_col_azimu_b_pos_xx__','Interpreter','none'); colorbar;
subplot(p_row,p_col,1+np);np=np+1; imagesc(tmp_col_polar_a_pos_xx__); axis image; axisnotick; xlabel('b'); ylabel('a'); title('tmp_col_polar_a_pos_xx__','Interpreter','none'); colorbar;
subplot(p_row,p_col,1+np);np=np+1; imagesc(tmp_f_check_pre_xx__,flim_); axis image; axisnotick; xlabel('b'); ylabel('a'); title('tmp_f_check_pre_xx__','Interpreter','none'); colorbar;
subplot(p_row,p_col,1+np);np=np+1; surfc(tmp_f_check_pre_xx__); clim(flim_); zlim(flim_); xlabel('b'); ylabel('a'); title('tmp_f_check_pre_xx__','Interpreter','none'); colorbar;
subplot(p_row,p_col,1+np);np=np+1; imagesc(tmp_f_check_pos_xx__,flim_); axis image; axisnotick; xlabel('b'); ylabel('a'); title('tmp_f_check_pos_xx__','Interpreter','none'); colorbar;
subplot(p_row,p_col,1+np);np=np+1; surfc(tmp_f_check_pos_xx__); clim(flim_); zlim(flim_); xlabel('b'); ylabel('a'); title('tmp_f_check_pos_xx__','Interpreter','none'); colorbar;
subplot(p_row,p_col,1+np);np=np+1; imagesc(tmp_f_check_xx__,flim_); axis image; axisnotick; xlabel('b'); ylabel('a'); title('tmp_f_check_xx__','Interpreter','none'); colorbar;
subplot(p_row,p_col,1+np);np=np+1; surfc(tmp_f_check_xx__); clim(flim_); zlim(flim_); xlabel('b'); ylabel('a'); title('tmp_f_check_xx__','Interpreter','none'); colorbar;
subplot(p_row,p_col,1+np);np=np+1; imagesc(tmp_weight_xx__); axis image; axisnotick; xlabel('b'); ylabel('a'); title('tmp_weight_xx__','Interpreter','none'); colorbar;
end;%if (flag_disp>0);
%%%%;
tmp_f_check_0 = f_check(polar_a_scatter_(1+nscatter),azimu_b_scatter_(1+nscatter));
tmp_f_check_1 = f_check(polar_a_scatter_periodize_(1+nscatter),azimu_b_scatter_periodize_(1+nscatter));
tmp_f_check_2 = f_check_scatter_periodize_(1+nscatter);
tmp_f_check_3 = sum(tmp_weight_xx__.*tmp_f_check_xx__,'all');
tmp_f_check_4 = sum(tmp_weight_xx__.*tmp_f_check_pre_xx__,'all');
if (flag_verbose>0); disp(sprintf(' %% tmp_f_check_0 vs tmp_f_check_1: %0.16f',fnorm(tmp_f_check_0 - tmp_f_check_1)/max(1e-12,fnorm(tmp_f_check_0)))); end;
if (flag_verbose>0); disp(sprintf(' %% tmp_f_check_0 vs tmp_f_check_2: %0.16f',fnorm(tmp_f_check_0 - tmp_f_check_2)/max(1e-12,fnorm(tmp_f_check_0)))); end;
if (flag_verbose>0); disp(sprintf(' %% tmp_f_check_0 vs tmp_f_check_3: %0.16f',fnorm(tmp_f_check_0 - tmp_f_check_3)/max(1e-12,fnorm(tmp_f_check_0)))); end;
if (flag_verbose>0); disp(sprintf(' %% tmp_f_check_0 vs tmp_f_check_4: %0.16f',fnorm(tmp_f_check_0 - tmp_f_check_4)/max(1e-12,fnorm(tmp_f_check_0)))); end;
if (flag_verbose>0); disp(sprintf(' %% tmp_f_check_1 vs tmp_f_check_2: %0.16f',fnorm(tmp_f_check_1 - tmp_f_check_2)/max(1e-12,fnorm(tmp_f_check_1)))); end;
if (flag_verbose>0); disp(sprintf(' %% tmp_f_check_1 vs tmp_f_check_3: %0.16f',fnorm(tmp_f_check_1 - tmp_f_check_3)/max(1e-12,fnorm(tmp_f_check_1)))); end;
if (flag_verbose>0); disp(sprintf(' %% tmp_f_check_1 vs tmp_f_check_4: %0.16f',fnorm(tmp_f_check_1 - tmp_f_check_4)/max(1e-12,fnorm(tmp_f_check_1)))); end;
if (flag_verbose>0); disp(sprintf(' %% tmp_f_check_2 vs tmp_f_check_3: %0.16f',fnorm(tmp_f_check_2 - tmp_f_check_3)/max(1e-12,fnorm(tmp_f_check_2)))); end;
if (flag_verbose>0); disp(sprintf(' %% tmp_f_check_2 vs tmp_f_check_4: %0.16f',fnorm(tmp_f_check_2 - tmp_f_check_4)/max(1e-12,fnorm(tmp_f_check_2)))); end;
if (flag_verbose>0); disp(sprintf(' %% tmp_f_check_3 vs tmp_f_check_4: %0.16f',fnorm(tmp_f_check_3 - tmp_f_check_4)/max(1e-12,fnorm(tmp_f_check_3)))); end;
%%%%;
disp('returning'); return;
end;%if (flag_verbose>1);
%%%%%%%%;

%%%%%%%%;
% periodize azimu_b (standard). ;
%%%%%%%%;
tmp_index_ = efind(index_col_azimu_b_xs__>=n_azimu_b); 
index_col_azimu_b_xs__(1+tmp_index_) = index_col_azimu_b_xs__(1+tmp_index_) - n_azimu_b;
tmp_index_ = efind(index_col_azimu_b_xs__< 0        ); 
index_col_azimu_b_xs__(1+tmp_index_) = index_col_azimu_b_xs__(1+tmp_index_) + n_azimu_b;
%%%%%%%%;
% periodize polar_a (note behavior at poles). ;
%%%%%%%%;
tmp_index_ = efind(index_col_polar_a_xs__< 0        ); 
index_col_polar_a_xs__(1+tmp_index_) = 0*n_polar_a - (index_col_polar_a_xs__(1+tmp_index_) - 0*n_polar_a);
flag_flip_azimu_b_xs__(1+tmp_index_) = 1-flag_flip_azimu_b_xs__(1+tmp_index_);
tmp_index_ = efind(index_col_polar_a_xs__> (n_polar_a-1)); 
index_col_polar_a_xs__(1+tmp_index_) = 1*(n_polar_a-1) - (index_col_polar_a_xs__(1+tmp_index_) - 1*(n_polar_a-1));
flag_flip_azimu_b_xs__(1+tmp_index_) = 1-flag_flip_azimu_b_xs__(1+tmp_index_);
%%%%;
if (flag_verbose>0); disp(sprintf(' %% n_azimu_b %d: index_col_azimu_b_xs__ in [%.4d,%.4d]',n_azimu_b,min(index_col_azimu_b_xs__,[],'all'),max(index_col_azimu_b_xs__,[],'all'))); end;
if (flag_verbose>0); disp(sprintf(' %% n_polar_a %d: index_col_polar_a_xs__ in [%.4d,%.4d]',n_polar_a,min(index_col_polar_a_xs__,[],'all'),max(index_col_polar_a_xs__,[],'all'))); end;
%%%%%%%%;
weight_azimu_b_xs__ = diag(1./weight_denominator_x_)*weight_numerator_azimu_b_xs__;
weight_polar_a_xs__ = diag(1./weight_denominator_x_)*weight_numerator_polar_a_xs__;
weight_xxs__ = zeros(n_nodex^2,n_scatter);
index_col_xxs__ = zeros(n_nodex^2,n_scatter);
for nscatter=0:n_scatter-1;
tmp_weight_xx__ = weight_azimu_b_xs__(:,1+nscatter)*transpose(weight_polar_a_xs__(:,1+nscatter));
weight_xxs__(:,1+nscatter) = reshape(tmp_weight_xx__,[n_nodex^2,1]);
tmp_col_azimu_b_xx__ = repmat(index_col_azimu_b_xs__(:,1+nscatter),[1,n_nodex]);
tmp_col_polar_a_xx__ = repmat(transpose(index_col_polar_a_xs__(:,1+nscatter)),[n_nodex,1]);
tmp_sgn_azimu_b_xx__ = repmat(transpose(flag_flip_azimu_b_xs__(:,1+nscatter)),[n_nodex,1]);
tmp_col_xx__ = periodize(tmp_col_azimu_b_xx__ + tmp_sgn_azimu_b_xx__*round(0.5*n_azimu_b),0,n_azimu_b) + tmp_col_polar_a_xx__*n_azimu_b ;
index_col_xxs__(:,1+nscatter) = reshape(tmp_col_xx__,[n_nodex^2,1]);
end;%for nscatter=0:n_scatter-1;
%%%%%%%%;
scatter_from_tensor_sba__ = sparse(1+index_row_xxs__(:),1+index_col_xxs__(:),weight_xxs__(:),n_scatter,n_azimu_b*n_polar_a);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%;
if flag_check;
f_check_recover_periodize_ = scatter_from_tensor_sba__*f_check_ba__(:);
figure(1+nf);nf=nf+1;clf;figmed; figbeach();
markersize_use = 12;
subplot(1,3,1);
scatter(azimu_b_scatter_periodize_,log10(abs(f_check_scatter_periodize_-f_check_recover_periodize_)),markersize_use,azimu_b_scatter_periodize_,'filled');
xlabel('azimu_b_scatter_periodize_','Interpreter','none');
ylabel('log10(abs(f_check_scatter_periodize_ - f_check_recover_periodize_))','Interpreter','none');
title('color: azimu_b_scatter_periodize_','Interpreter','none');
colorbar;
subplot(1,3,2);
scatter(polar_a_scatter_periodize_rescale_,log10(abs(f_check_scatter_periodize_-f_check_recover_periodize_)),markersize_use,polar_a_scatter_periodize_rescale_,'filled');
xlabel('polar_a_scatter_periodize_rescale_','Interpreter','none');
ylabel('log10(abs(f_check_scatter_periodize_-f_check_recover_periodize_))','Interpreter','none');
title('color: polar_a_scatter_periodize_rescale_','Interpreter','none');
colorbar;
subplot(1,3,3);
scatter(azimu_b_scatter_rescale_,polar_a_scatter_periodize_rescale_,markersize_use,log10(abs(f_check_scatter_periodize_-f_check_recover_periodize_)),'filled');
xlabel('azimu_b_scatter_rescale_','Interpreter','none');
ylabel('polar_a_scatter_periodize_rescale_','Interpreter','none');
title('color: log10(abs(f_check_scatter_periodize_-f_check_recover_periodize_))','Interpreter','none');
colorbar;
end;%if flag_check;
%%%%;

if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;


