%%%%%%%%;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
%%%%%%%%;

flag_verbose = 1;
flag_replot = 1;
dir_jpg = sprintf('/%s/rangan/dir_cryoem/dir_jpg',string_root);
dir_trunk = sprintf('/%s/rangan/dir_cryoem/dir_trpv1_x0',string_root);

nf=0;
kx_ = 16*[2:6]; n_kx = numel(kx_);

%%%%%%%%;
fname_fig_pre = sprintf('%s/trpv1_kxxxx_pm_fig_X_2d_Semp_d1_FIGM__',dir_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_jpg));
figure(1+nf);nf=nf+1;clf;figbig;
linewidth_use = 2; markersize_use = 10; fontsize_use = 16;
p_row = 1; p_col = 3; np=0;
tmp_symbol_ = {'o-','o-','o-','o-','o-'};
tmp_intens_l_ = 0.85.^[1:n_kx-0];
tmp_intens_m_ = 0.85.^[0:n_kx-1];
tmp_intens_r_ = 0.85.^[0:n_kx-1];
tmp_legend_l_ = ...
{ ...
   'K = 32 frob' ...
   'K = 48 frob' ...
   'K = 64 frob' ...
   'K = 80 frob' ...
   'K = 96 frob' ...
};
tmp_legend_m_ = ...
{ ...
   'K = 32 corr' ...
   'K = 48 corr' ...
   'K = 64 corr' ...
   'K = 80 corr' ...
   'K = 96 corr' ...
   'K = 32 prct' ...
   'K = 48 prct' ...
   'K = 64 prct' ...
   'K = 80 prct' ...
   'K = 96 prct' ...
};
tmp_legend_r_ = ...
{ ...
   'K = 32 time' ...
   'K = 48 time' ...
   'K = 64 time' ...
   'K = 80 time' ...
   'K = 96 time' ...
   'K = 32 ops' ...
   'K = 48 ops' ...
   'K = 64 ops' ...
   'K = 80 ops' ...
   'K = 96 ops' ...
};
%%%%;
tmp_ = cell(n_kx,1);
for nkx=0:n_kx-1;
kx = kx_(1+nkx);
tmp_{1+nkx} = load(sprintf('%s/dir_k%.3d_pm_mat/pm_fig_X_2d_Semp_d1_FIGH__.mat',dir_trunk,kx));
tmp_{1+nkx}.tot_t_all_ = tmp_{1+nkx}.tot_t_fill_ + tmp_{1+nkx}.tot_t_mult_ + tmp_{1+nkx}.tot_t_ifft_ ;
end;%for nkx=0:n_kx-1;
%%%%;
subplot(p_row,p_col,1+np); np=np+1;
for nkx=0:n_kx-1;
tmp_symbol = tmp_symbol_{1+nkx}; tmp_intens_l = tmp_intens_l_(1+nkx);
hold on;
plot(1:numel(tmp_{1+nkx}.tot_t_fill_),tmp_{1+nkx}.l2er_X_X_n_,tmp_symbol,'MarkerEdgeColor','k','MarkerFaceColor',[1,1,1]*tmp_intens_l,'Color',[1,1,1]*tmp_intens_l,'LineWidth',linewidth_use,'MarkerSize',markersize_use);
hold off;
end;%for nkx=0:n_kx-1;
xlim([0,1+16]); set(gca,'XTick',0:4:16); xlabel('rank $H$','Interpreter','latex');
ylabel('error'); grid on;
legend(tmp_legend_l_,'Location','NorthEast');
set(gca,'FontSize',fontsize_use);
%title(sprintf('relative error: K=%.2d',kx));
title(sprintf('relative error:'));
%%%%;
subplot(p_row,p_col,1+np); np=np+1;
for nkx=0:n_kx-1;
tmp_symbol = tmp_symbol_{1+nkx}; tmp_intens_m = tmp_intens_m_(1+nkx);
hold on;
%plot(1:n_UX_rank,0+corr_X_X_n_,'ko-','MarkerFaceColor','r');
plot(1:numel(tmp_{1+nkx}.tot_t_fill_),0+tmp_{1+nkx}.corr_X_X_n_,tmp_symbol,'MarkerEdgeColor','k','MarkerFaceColor',[1,0,0]*tmp_intens_m,'Color',[1,0,0]*tmp_intens_m,'LineWidth',linewidth_use,'MarkerSize',markersize_use);
hold off;
end;%for nkx=0:n_kx-1;
for nkx=0:n_kx-1;
tmp_symbol = tmp_symbol_{1+nkx}; tmp_intens_m = tmp_intens_m_(1+nkx);
hold on;
%plot(1:n_UX_rank,1-prct_X_X_n_,'ko-','MarkerFaceColor','c');
plot(1:numel(tmp_{1+nkx}.tot_t_fill_),1-tmp_{1+nkx}.prct_X_X_n_,tmp_symbol,'MarkerEdgeColor','k','MarkerFaceColor',[0,1,1]*tmp_intens_m,'Color',[0,1,1]*tmp_intens_m,'LineWidth',linewidth_use,'MarkerSize',markersize_use);
hold off;
end;%for nkx=0:n_kx-1;
xlim([0,1+16]); set(gca,'XTick',0:4:16); xlabel('rank $H$','Interpreter','latex');
ylim([0.5,1]); ylabel('value'); set(gca,'YTick',0.5:0.05:1.0); grid on;
legend(tmp_legend_m_,'Location','SouthEast');
set(gca,'FontSize',fontsize_use);
%title(sprintf('correlation: K=%.2d',kx));
title(sprintf('correlation:'));
%%%%;
subplot(p_row,p_col,1+np); np=np+1;
for nkx=0:n_kx-1;
tmp_symbol = tmp_symbol_{1+nkx}; tmp_intens_r = tmp_intens_r_(1+nkx);
hold on;
plot(1:numel(tmp_{1+nkx}.tot_t_fill_),max(tmp_{1+nkx}.tot_t_mult_+tmp_{1+nkx}.tot_t_ifft_)./tmp_{1+nkx}.tot_t_all_,tmp_symbol,'MarkerEdgeColor','k','MarkerFaceColor',[1,0,0]*tmp_intens_r,'Color',[1,0,0]*tmp_intens_r,'LineWidth',linewidth_use,'MarkerSize',markersize_use);
hold off;
end;%for nkx=0:n_kx-1;
for nkx=0:n_kx-1;
tmp_symbol = tmp_symbol_{1+nkx}; tmp_intens_r = tmp_intens_r_(1+nkx);
hold on;
plot(1:numel(tmp_{1+nkx}.tot_t_fill_),max(tmp_{1+nkx}.tot_o_mult_)./tmp_{1+nkx}.tot_o_mult_,tmp_symbol,'MarkerEdgeColor','k','MarkerFaceColor',[0,1,1]*tmp_intens_r,'Color',[0,1,1]*tmp_intens_r,'LineWidth',linewidth_use,'MarkerSize',markersize_use);
hold off;
end;%for nkx=0:n_kx-1;
for nkx=0:n_kx-1;
hold on;
plot(1:n_UX_rank,ones(1,n_UX_rank),'k-','LineWidth',0.5);
hold off;
end;%for nkx=0:n_kx-1;
xlim([0,1+16]); set(gca,'XTick',0:4:16); xlabel('rank $H$','Interpreter','latex');
ylim([0,32]); ylabel('factor'); set(gca,'YTick',0:2:32); grid on;
legend(tmp_legend_r_,'Location','NorthEast');
set(gca,'FontSize',fontsize_use);
%title(sprintf('speedup: K=%.2d',kx));
title(sprintf('speedup:'));
%%%%;
disp(sprintf(' %% writing %s',fname_fig_jpg));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
%close(gcf);
end;%if (~exist(fname_fig_jpg,'file'));
%%%%%%%%;
