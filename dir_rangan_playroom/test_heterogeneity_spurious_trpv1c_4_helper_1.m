%%%%%%%%;
% Now display log-unlikelihood-ratio vs temperature for rem2 and rem3. ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figmed;
fontsize_use = 12;
markersize_use = 12;
subplot(1,2,1);
hold on;
plot(1:n_sigma_bayesian-1,ssnll2r_uni_tr0_vs_hr2_s_(2:end),'ko','MarkerSize',markersize_use,'MarkerFaceColor','c');
plot(1:n_sigma_bayesian-1,ssnll2r_uni_tr0_vs_hr3_s_(2:end),'ko','MarkerSize',markersize_use,'MarkerFaceColor','m');
%plot(1:n_sigma_bayesian-1,ssnll2r_emp_tr0_vs_hr2_s_(2:end),'ks','MarkerSize',markersize_use,'MarkerFaceColor','c');
%plot(1:n_sigma_bayesian-1,ssnll2r_emp_tr0_vs_hr3_s_(2:end),'ks','MarkerSize',markersize_use,'MarkerFaceColor','m');
ylabel('$\sigma^{2}\cdot$ nllr','Interpreter','latex');
xlim([0.5,0.5+n_sigma_bayesian-1]); xlabel('$-\log_{2}$(temperature)','Interpreter','latex');
set(gca,'XTick',1:4:n_sigma_bayesian-1,'XTickLabel',num2str(-log2(sigma_bayesian_(2:4:end)),'%0.2f'));xtickangle(90); grid on;
%legend({'uni tr0 vs hr2','uni tr0 vs hr3','emp tr0 vs hr2','emp tr0 vs hr3'},'Location','NorthEast');
set(gca,'FontSize',fontsize_use);
hold off;
%%%%;
subplot(1,2,2);
hold on;
ylim_ = [-3.5,+3.5];
tmp_y_ = ssnll2r_uni_tr0_vs_hr2_s_(2:end)./transpose(sigma_bayesian_(2:end)).^2; tmp_y_ = min(max(ylim_),max(min(ylim_),tmp_y_));
plot(1:n_sigma_bayesian-1,tmp_y_,'ko','MarkerSize',markersize_use,'MarkerFaceColor','c');
tmp_y_ = ssnll2r_uni_tr0_vs_hr3_s_(2:end)./transpose(sigma_bayesian_(2:end)).^2; tmp_y_ = min(max(ylim_),max(min(ylim_),tmp_y_));
plot(1:n_sigma_bayesian-1,tmp_y_,'ko','MarkerSize',markersize_use,'MarkerFaceColor','m');
%tmp_y_ = ssnll2r_emp_tr0_vs_hr2_s_(2:end)./transpose(sigma_bayesian_(2:end)).^2; tmp_y_ = min(max(ylim_),max(min(ylim_),tmp_y_));
%plot(1:n_sigma_bayesian-1,tmp_y_,'ks','MarkerSize',markersize_use,'MarkerFaceColor','c');
%tmp_y_ = ssnll2r_emp_tr0_vs_hr3_s_(2:end)./transpose(sigma_bayesian_(2:end)).^2; tmp_y_ = min(max(ylim_),max(min(ylim_),tmp_y_));
%plot(1:n_sigma_bayesian-1,tmp_y_,'ks','MarkerSize',markersize_use,'MarkerFaceColor','m');
hold off;
ylim(ylim_); ylabel('nllr','Interpreter','latex');
xlim([0.5,0.5+n_sigma_bayesian-1]); xlabel('$-\log_{2}$(temperature)','Interpreter','latex');
set(gca,'XTick',1:4:n_sigma_bayesian-1,'XTickLabel',num2str(-log2(sigma_bayesian_(2:4:end)),'%0.2f'));xtickangle(90); grid on;
%legend({'uni tr0 vs hr2','uni tr0 vs hr3','emp tr0 vs hr2','emp tr0 vs hr3'},'Location','NorthWest');
legend({'uni tr0 vs hr2','uni tr0 vs hr3'},'Location','NorthWest');
set(gca,'FontSize',fontsize_use);
%%%%%%%%;
fname_fig_pre = sprintf('%s_jpg/test_heterogeneity_spurious_trpv1c_ssnll2r_sigma_opt_R%.4d_FIGB',dir_pm,n_R);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
fname_fig_stripped_jpg = sprintf('%s_stripped.jpg',fname_fig_pre);
fname_fig_stripped_eps = sprintf('%s_stripped.eps',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
if (verbose); disp(sprintf(' %% %s not found, creating',fname_fig_jpg)); end;
print('-djpeg',fname_fig_stripped_jpg);
print('-depsc',fname_fig_stripped_eps);
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
end;%if (flag_replot | ~exist(fname_fig_jpg,'file'));
%%%%;
close(gcf);
%%%%%%%%;

%%%%%%%%;
% Now display log-unlikelihood-ratio vs temperature for rem2 and rem3. ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figmed;
fontsize_use = 12;
markersize_use = 12;
subplot(1,2,1);
hold on;
plot(1:n_sigma_bayesian-1,ssnll2r_uni_tr1_vs_hr2_s_(2:end),'ko','MarkerSize',markersize_use,'MarkerFaceColor','c');
plot(1:n_sigma_bayesian-1,ssnll2r_uni_tr1_vs_hr3_s_(2:end),'ko','MarkerSize',markersize_use,'MarkerFaceColor','m');
%plot(1:n_sigma_bayesian-1,ssnll2r_emp_tr1_vs_hr2_s_(2:end),'ks','MarkerSize',markersize_use,'MarkerFaceColor','c');
%plot(1:n_sigma_bayesian-1,ssnll2r_emp_tr1_vs_hr3_s_(2:end),'ks','MarkerSize',markersize_use,'MarkerFaceColor','m');
hold off;
ylabel('$\sigma^{2}\cdot$ nllr','Interpreter','latex');
xlim([0.5,0.5+n_sigma_bayesian-1]); xlabel('$-\log_{2}$(temperature)','Interpreter','latex');
set(gca,'XTick',1:4:n_sigma_bayesian-1,'XTickLabel',num2str(-log2(sigma_bayesian_(2:4:end)),'%0.2f'));xtickangle(90); grid on;
%legend({'uni tr1 vs hr2','uni tr1 vs hr3','emp tr1 vs hr2','emp tr1 vs hr3'},'Location','NorthEast');
set(gca,'FontSize',fontsize_use);
%%%%;
subplot(1,2,2);
hold on;
ylim_ = [-3.5,+3.5];
tmp_y_ = ssnll2r_uni_tr1_vs_hr2_s_(2:end)./transpose(sigma_bayesian_(2:end)).^2; tmp_y_ = min(max(ylim_),max(min(ylim_),tmp_y_));
plot(1:n_sigma_bayesian-1,tmp_y_,'ko','MarkerSize',markersize_use,'MarkerFaceColor','c');
tmp_y_ = ssnll2r_uni_tr1_vs_hr3_s_(2:end)./transpose(sigma_bayesian_(2:end)).^2; tmp_y_ = min(max(ylim_),max(min(ylim_),tmp_y_));
plot(1:n_sigma_bayesian-1,tmp_y_,'ko','MarkerSize',markersize_use,'MarkerFaceColor','m');
%tmp_y_ = ssnll2r_emp_tr1_vs_hr2_s_(2:end)./transpose(sigma_bayesian_(2:end)).^2; tmp_y_ = min(max(ylim_),max(min(ylim_),tmp_y_));
%plot(1:n_sigma_bayesian-1,tmp_y_,'ks','MarkerSize',markersize_use,'MarkerFaceColor','c');
%tmp_y_ = ssnll2r_emp_tr1_vs_hr3_s_(2:end)./transpose(sigma_bayesian_(2:end)).^2; tmp_y_ = min(max(ylim_),max(min(ylim_),tmp_y_));
%plot(1:n_sigma_bayesian-1,tmp_y_,'ks','MarkerSize',markersize_use,'MarkerFaceColor','m');
hold off;
ylim(ylim_); ylabel('nllr','Interpreter','latex');
xlim([0.5,0.5+n_sigma_bayesian-1]); xlabel('$-\log_{2}$(temperature)','Interpreter','latex');
set(gca,'XTick',1:4:n_sigma_bayesian-1,'XTickLabel',num2str(-log2(sigma_bayesian_(2:4:end)),'%0.2f'));xtickangle(90); grid on;
%legend({'uni tr1 vs hr2','uni tr1 vs hr3','emp tr1 vs hr2','emp tr1 vs hr3'},'Location','NorthWest');
legend({'uni tr1 vs hr2','uni tr1 vs hr3'},'Location','NorthWest');
set(gca,'FontSize',fontsize_use);
%%%%%%%%;
fname_fig_pre = sprintf('%s_jpg/test_heterogeneity_spurious_trpv1c_ssnll2r_sigma_opt_R%.4d_FIGC',dir_pm,n_R);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
fname_fig_stripped_jpg = sprintf('%s_stripped.jpg',fname_fig_pre);
fname_fig_stripped_eps = sprintf('%s_stripped.eps',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
if (verbose); disp(sprintf(' %% %s not found, creating',fname_fig_jpg)); end;
print('-djpeg',fname_fig_stripped_jpg);
print('-depsc',fname_fig_stripped_eps);
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
end;%if (flag_replot | ~exist(fname_fig_jpg,'file'));
%%%%;
close(gcf);
%%%%%%%%;
