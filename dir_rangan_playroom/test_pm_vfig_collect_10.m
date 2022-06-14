%%%%%%%%;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
%%%%%%%%;
%{
  mkdir /home/rangan/dir_cryoem/dir_trpv1c/dir_pm_mat; rsync -avum rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_trpv1c/dir_pm_mat/a_x_u_pack_.mat rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_trpv1c/dir_pm_mat/pm_fig_X_2d_Semp_d1_FIGH__.mat rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_trpv1c/dir_pm_mat/pm_vfig_timing_FIGB__.mat /home/rangan/dir_cryoem/dir_trpv1c/dir_pm_mat;
  mkdir /home/rangan/dir_cryoem/dir_rib80s/dir_pm_mat;  rsync -avum rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_rib80s/dir_pm_mat/a_x_u_pack_.mat rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_rib80s_x0/dir_pm_mat/pm_fig_X_2d_Semp_d1_FIGH__.mat rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_rib80s/dir_pm_mat/pm_vfig_timing_FIGB__.mat /home/rangan/dir_cryoem/dir_rib80s/dir_pm_mat;
  mkdir /home/rangan/dir_cryoem/dir_p28hRPT1_x0/dir_pm_mat;  rsync -avum rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_p28hRPT1_x0/dir_pm_mat/a_x_u_pack_.mat rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_p28hRPT1_x0/dir_pm_mat/pm_fig_X_2d_Semp_d1_FIGH__.mat rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_p28hRPT1_x0/dir_pm_mat/pm_vfig_timing_FIGB__.mat /home/rangan/dir_cryoem/dir_p28hRPT1_x0/dir_pm_mat;
  mkdir /home/rangan/dir_cryoem/dir_MlaFEDB/dir_pm_mat;  rsync -avum rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_MlaFEDB/dir_pm_mat/a_x_u_pack_.mat rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_MlaFEDB_x0/dir_pm_mat/pm_fig_X_2d_Semp_d1_FIGH__.mat rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_MlaFEDB/dir_pm_mat/pm_vfig_timing_FIGB__.mat /home/rangan/dir_cryoem/dir_MlaFEDB/dir_pm_mat;
  mkdir /home/rangan/dir_cryoem/dir_LetB1/dir_pm_mat;  rsync -avum rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_LetB1/dir_pm_mat/a_x_u_pack_.mat rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_LetB1_x0/dir_pm_mat/pm_fig_X_2d_Semp_d1_FIGH__.mat rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_LetB1/dir_pm_mat/pm_vfig_timing_FIGB__.mat /home/rangan/dir_cryoem/dir_LetB1/dir_pm_mat;
  mkdir /home/rangan/dir_cryoem/dir_ISWINCP/dir_pm_mat;  rsync -avum rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_ISWINCP/dir_pm_mat/a_x_u_pack_.mat rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_ISWINCP_x0/dir_pm_mat/pm_fig_X_2d_Semp_d1_FIGH__.mat rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_ISWINCP/dir_pm_mat/pm_vfig_timing_FIGB__.mat /home/rangan/dir_cryoem/dir_ISWINCP/dir_pm_mat;
  mkdir /home/rangan/dir_cryoem/dir_TMEM16F/dir_pm_mat;  rsync -avum rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_TMEM16F/dir_pm_mat/a_x_u_pack_.mat rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_TMEM16F_x0/dir_pm_mat/pm_fig_X_2d_Semp_d1_FIGH__.mat rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_TMEM16F/dir_pm_mat/pm_vfig_timing_FIGB__.mat /home/rangan/dir_cryoem/dir_TMEM16F/dir_pm_mat;
  mkdir /home/rangan/dir_cryoem/dir_LSUbl17dep/dir_pm_mat;  rsync -avum rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_LSUbl17dep/dir_pm_mat/a_x_u_pack_.mat rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_LSUbl17dep_x0/dir_pm_mat/pm_fig_X_2d_Semp_d1_FIGH__.mat rangan@access1.cims.nyu.edu:/data/rangan/dir_cryoem/dir_LSUbl17dep/dir_pm_mat/pm_vfig_timing_FIGB__.mat /home/rangan/dir_cryoem/dir_LSUbl17dep/dir_pm_mat;
%}
nf=0;

table_data__ = { ...
'trpv1c' , 'trpv1' ; ...
'rib80s' , 'rib80s' ; ...
'p28hRPT1_x0' , 'p28hRPT1' ; ...
'MlaFEDB' , 'MlaFEDB' ; ...
%'LetB1' , 'LetB1' ; ...
'ISWINCP' , 'ISWINCP' ; ...
'TMEM16F' , 'TMEM16F' ; ...
'LSUbl17dep' , 'LSUbl17dep' ; ...
};
n_experiment = size(table_data__,1);

flag_replot=1;
n_k_p_r = 49;
n_UX_rank = 48;
n_UZ_rank = 49;
tmp_ij_ = 1:8;

fname_fig = sprintf('/%s/rangan/dir_cryoem/dir_pm_manuscript/dir_pm_fig/pm_vfig_collect_10',string_root);
if ( flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
figure(1+nf);nf=nf+1;figbig;
p_row = 3; p_col = n_experiment; ns=0;
for nexperiment=0:n_experiment-1;
fname_prefix_xfix = table_data__{1+nexperiment,1+0};
str_title = table_data__{1+nexperiment,1+1};
dir_pm = sprintf('/%s/rangan/dir_cryoem/dir_%s/dir_pm',string_root,fname_prefix_xfix);
fname_afig = sprintf('%s_mat/a_x_u_pack_.mat',dir_pm);
if (strcmp(fname_prefix_xfix,'p28hRPT1_x0'));
fname_afig = sprintf('%s_mat/a_x_u_base_.mat',dir_pm);
end;%if (strcmp(fname_prefix_xfix,'p28hRPT1_x0'));
fname_sfig = sprintf('%s_mat/pm_fig_X_2d_Semp_d1_FIGH__.mat',dir_pm);
fname_vfig = sprintf('%s_mat/pm_vfig_timing_FIGB__.mat',dir_pm);
flag_exist = ( exist(fname_afig,'file')) & ( exist(fname_sfig,'file')) & ( exist(fname_vfig,'file'));
disp(sprintf(' %% %s --> exist %d',fname_prefix_xfix,flag_exist));
%%%%%%%%%%%%%%%%;
if flag_exist;
%%%%%%%%%%%%%%%%;
tmp_afig_ = load(fname_afig);
tmp_sfig_ = load(fname_sfig);
tmp_vfig_ = load(fname_vfig);
corr_X_X_n_ = tmp_sfig_.corr_X_X_n_;
prct_X_X_n_ = tmp_sfig_.prct_X_X_n_;
corr_X_X_z_ = mean(tmp_vfig_.corr_X_X_zb__,2);
prct_X_X_z_ = mean(tmp_vfig_.prct_X_X_zb__,2);
if isfield(tmp_afig_,'a_x_u_pack_');
a_x_u_pack_ = tmp_afig_.a_x_u_pack_;
end;%if isfield(tmp_afig_,'a_x_u_pack_');
if isfield(tmp_afig_,'a_x_u_base_');
a_x_u_pack_ = tmp_afig_.a_x_u_base_;
end;%if isfield(tmp_afig_,'a_x_u_base_');
n_x_u_pack = tmp_afig_.n_x_u_pack;
%tmp_n_x = round(0.5352*n_x_u_pack/2);
tmp_n_x = round(0.65*n_x_u_pack/2);
tmp_x_ij = n_x_u_pack/2 + [-tmp_n_x:+tmp_n_x];
%tmp_prct_0 = 95.0;
tmp_prct_1 = 97.5;
tmp_a_x_u_ = reshape(a_x_u_pack_,[n_x_u_pack,n_x_u_pack,n_x_u_pack]);
%%%%;
%subplot(p_row,p_col,1+ns+0*p_col);ns=ns+0;
%isosurface_f_x_u_0(tmp_a_x_u_(tmp_x_ij,tmp_x_ij,tmp_x_ij),tmp_prct_0);
%axisnotick; xlabel([]); ylabel([]); zlabel([]); title(sprintf('%s %0.2f%%',str_title,tmp_prct_0))
%%%%;
subplot(p_row,p_col,1+ns+0*p_col);ns=ns+0;
isosurface_f_x_u_0(tmp_a_x_u_(tmp_x_ij,tmp_x_ij,tmp_x_ij),tmp_prct_1);
axisnotick; xlabel([]); ylabel([]); zlabel([]); title(sprintf('%s %0.2f%%',str_title,tmp_prct_1))
%%%%;
subplot(p_row,p_col,1+ns+1*p_col);ns=ns+0;
hold on;
plot(tmp_ij_,0+corr_X_X_n_(tmp_ij_),'ko-','MarkerFaceColor','r');
plot(tmp_ij_,1-prct_X_X_n_(tmp_ij_),'ko-','MarkerFaceColor','c');
hold off;
xlim([0,max(tmp_ij_)+1]); set(gca,'XTick',tmp_ij_); xlabel('rank $H$','Interpreter','latex');
ylim([0.5,1]); ylabel('value'); set(gca,'YTick',0.5:0.05:1.0); grid on;
legend({'corr','prct'},'Location','SouthEast');
title('image alignment');
%%%%;
subplot(p_row,p_col,1+ns+2*p_col);ns=ns+1;
hold on;
plot(tmp_ij_,0+corr_X_X_z_(tmp_ij_),'ko-','MarkerFaceColor','r');
plot(tmp_ij_,1-prct_X_X_z_(tmp_ij_),'ko-','MarkerFaceColor','c');
hold off;
xlim([0,max(tmp_ij_)+1]); set(gca,'XTick',tmp_ij_); xlabel('rank $H$','Interpreter','latex');
ylim([0.5,1]); ylabel('value'); set(gca,'YTick',0.5:0.05:1.0); grid on;
legend({'corr','prct'},'Location','SouthEast');
title('volume alignment');
%%%%%%%%%%%%%%%%;
end;%if flag_exist;
%%%%%%%%%%%%%%%%;
end;%for nexperiment=0:n_experiment-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
