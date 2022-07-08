function ...
[ ...
 parameter ...
,str_fname_align_a_k_Y_mat_a_ ...
,fsc_ampm_ka__ ...
,corr_base_vs_ampm_a_ ...
,corr_reco_vs_ampm_a_ ...
,fsc_crop_ampm_kxa___ ...
,corr_crop_base_vs_crop_ampm_xa__ ...
,corr_full_base_vs_crop_ampm_xa__ ...
,corr_crop_reco_vs_crop_ampm_xa__ ...
,corr_full_reco_vs_crop_ampm_xa__ ...
,k_Ainv_p_r_ ...
,k_Ainv_p_r_max ...
,kinv_A_p_r_ ...
] = ...
test_pm_collect_excerpt_fsc_xa__1( ...
 parameter ...
,dir_pm ...
,dir_pm_bkp ...
);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if nargin<1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
setup_access1;
Pixel_Spacing = 1.699;
parameter = struct('type','parameter','Pixel_Spacing',Pixel_Spacing,'flag_recalc',0,'flag_verbose',1);
dir_base = '/data/rangan/dir_cryoem';
dir_pm = sprintf('%s/dir_ps1_x0to7/dir_pm',dir_base);
[ ...
 parameter ...
,str_fname_align_a_k_Y_mat_o0to7_a_ ...
,~ ...
,~ ...
,corr_reco_vs_ampm_o0to7_a_ ...
,~ ...
,~ ...
,~ ...
,~ ...
,corr_full_reco_vs_crop_ampm_o0to7_xa__ ...
,k_Ainv_p_r_ ...
,k_Ainv_p_r_max ...
,kinv_A_p_r_ ...
] = ...
test_pm_collect_excerpt_fsc_xa__1(parameter,dir_pm);
n_octile = 8;
%%%%%%%%;
for noctile=1:n_octile-1;
parameter = struct('type','parameter','Pixel_Spacing',Pixel_Spacing,'flag_recalc',0,'flag_verbose',1);
dir_base = '/data/rangan/dir_cryoem';
dir_pm = sprintf('%s/dir_ps1_x%d/dir_pm',dir_base,noctile);
[ ...
 parameter ...
,str_fname_align_a_k_Y_mat_ya__{1+noctile} ...
,~ ...
,~ ...
,corr_reco_vs_ampm_ya__{1+noctile} ...
,~ ...
,~ ...
,~ ...
,~ ...
,corr_full_reco_vs_crop_ampm_yxa___{1+noctile} ...
,~ ...
,~ ...
,~ ...
] = ...
test_pm_collect_excerpt_fsc_xa__1(parameter,dir_pm);
end;%for noctile=1:n_octile-1;
%str_fname_align_a_k_Y_mat_1a_ = str_fname_align_a_k_Y_mat_ya__{1+1}; %<-- assume all aligned. ;
%ncrop = 40;
%corr_reco_vs_ampm_Ya__ = zeros(n_octile,15);
%corr_full_reco_vs_crop_ampm_y40a__ = zeros(n_octile,15);
%for noctile=1:n_octile-1;
%corr_reco_vs_ampm_Ya__(1+noctile,:) = corr_reco_vs_ampm_ya__{1+noctile}(:);
%corr_full_reco_vs_crop_ampm_y40a__(1+noctile,:) = corr_full_reco_vs_crop_ampm_yxa___{1+noctile}(1+ncrop,:);
%end;%for noctile=1:n_octile-1;
%%%%%%%%;
str_fname_align_a_k_Y_mat_oa__ = cell(n_octile,1);
fsc_ampm_oka___ = cell(n_octile,1);
corr_base_vs_ampm_oa__ = cell(n_octile,1);
corr_reco_vs_ampm_oa__ = cell(n_octile,1);
fsc_crop_ampm_okxa____ = cell(n_octile,1);
corr_crop_base_vs_crop_ampm_oxa___ = cell(n_octile,1);
corr_full_base_vs_crop_ampm_oxa___ = cell(n_octile,1);
corr_crop_reco_vs_crop_ampm_oxa___ = cell(n_octile,1);
corr_full_reco_vs_crop_ampm_oxa___ = cell(n_octile,1);
for noctile=0:n_octile-1;
parameter = struct('type','parameter','Pixel_Spacing',Pixel_Spacing,'flag_recalc',0,'flag_verbose',1);
dir_base = '/data/rangan/dir_cryoem';
dir_pm = sprintf('%s/dir_ps1_x0to7_combine/dir_pmo%d',dir_base,noctile);
dir_pm_bkp = sprintf('%s/dir_ps1_x0/dir_pm',dir_base);
[ ...
 parameter ...
,str_fname_align_a_k_Y_mat_a_ ...
,fsc_ampm_ka__ ...
,corr_base_vs_ampm_a_ ...
,corr_reco_vs_ampm_a_ ...
,fsc_crop_ampm_kxa___ ...
,corr_crop_base_vs_crop_ampm_xa__ ...
,corr_full_base_vs_crop_ampm_xa__ ...
,corr_crop_reco_vs_crop_ampm_xa__ ...
,corr_full_reco_vs_crop_ampm_xa__ ...
,~ ...
,~ ...
,~ ...
] = ...
test_pm_collect_excerpt_fsc_xa__1( ...
 parameter ...
,dir_pm ...
,dir_pm_bkp ...
);
str_fname_align_a_k_Y_mat_oa__{1+noctile} = str_fname_align_a_k_Y_mat_a_;
fsc_ampm_oka___{1+noctile} = fsc_ampm_ka__;
corr_base_vs_ampm_oa__{1+noctile} = corr_base_vs_ampm_a_;
corr_reco_vs_ampm_oa__{1+noctile} = corr_reco_vs_ampm_a_;
fsc_crop_ampm_okxa____{1+noctile} = fsc_crop_ampm_kxa___;
corr_crop_base_vs_crop_ampm_oxa___{1+noctile} = corr_crop_base_vs_crop_ampm_xa__;
corr_full_base_vs_crop_ampm_oxa___{1+noctile} = corr_full_base_vs_crop_ampm_xa__;
corr_crop_reco_vs_crop_ampm_oxa___{1+noctile} = corr_crop_reco_vs_crop_ampm_xa__;
corr_full_reco_vs_crop_ampm_oxa___{1+noctile} = corr_full_reco_vs_crop_ampm_xa__;
end;%for noctile=0:n_octile-1;
%%%%%%%%;
str_fname_align_a_k_Y_mat_ua__ = cell(n_octile,1);
fsc_ampm_uka___ = cell(n_octile,1);
corr_base_vs_ampm_ua__ = cell(n_octile,1);
corr_reco_vs_ampm_ua__ = cell(n_octile,1);
fsc_crop_ampm_ukxa____ = cell(n_octile,1);
corr_crop_base_vs_crop_ampm_uxa___ = cell(n_octile,1);
corr_full_base_vs_crop_ampm_uxa___ = cell(n_octile,1);
corr_crop_reco_vs_crop_ampm_uxa___ = cell(n_octile,1);
corr_full_reco_vs_crop_ampm_uxa___ = cell(n_octile,1);
for noctile=0:n_octile-1;
parameter = struct('type','parameter','Pixel_Spacing',Pixel_Spacing,'flag_recalc',0,'flag_verbose',1);
dir_base = '/data/rangan/dir_cryoem';
dir_pm = sprintf('%s/dir_ps1_x0to7_combine/dir_pmox%d',dir_base,noctile);
dir_pm_bkp = sprintf('%s/dir_ps1_x0/dir_pm',dir_base);
[ ...
 parameter ...
,str_fname_align_a_k_Y_mat_a_ ...
,fsc_ampm_ka__ ...
,corr_base_vs_ampm_a_ ...
,corr_reco_vs_ampm_a_ ...
,fsc_crop_ampm_kxa___ ...
,corr_crop_base_vs_crop_ampm_xa__ ...
,corr_full_base_vs_crop_ampm_xa__ ...
,corr_crop_reco_vs_crop_ampm_xa__ ...
,corr_full_reco_vs_crop_ampm_xa__ ...
,~ ...
,~ ...
,~ ...
] = ...
test_pm_collect_excerpt_fsc_xa__1( ...
 parameter ...
,dir_pm ...
,dir_pm_bkp ...
);
str_fname_align_a_k_Y_mat_ua__{1+noctile} = str_fname_align_a_k_Y_mat_a_;
fsc_ampm_uka___{1+noctile} = fsc_ampm_ka__;
corr_base_vs_ampm_ua__{1+noctile} = corr_base_vs_ampm_a_;
corr_reco_vs_ampm_ua__{1+noctile} = corr_reco_vs_ampm_a_;
fsc_crop_ampm_ukxa____{1+noctile} = fsc_crop_ampm_kxa___;
corr_crop_base_vs_crop_ampm_uxa___{1+noctile} = corr_crop_base_vs_crop_ampm_xa__;
corr_full_base_vs_crop_ampm_uxa___{1+noctile} = corr_full_base_vs_crop_ampm_xa__;
corr_crop_reco_vs_crop_ampm_uxa___{1+noctile} = corr_crop_reco_vs_crop_ampm_xa__;
corr_full_reco_vs_crop_ampm_uxa___{1+noctile} = corr_full_reco_vs_crop_ampm_xa__;
end;%for noctile=0:n_octile-1;
%%%%%%%%;

%{
  %%%%%%%%;
  % does not work anymore because dir_ps1_x0to7/dir_pm_mat seems stale. ;
  %%%%%%%%;
flag_replot = 1;
fname_fig = sprintf('%s/dir_ps1_x0to7_combine/dir_pm_jpg/test_pm_collect_excerpt_fsc_xa__1',dir_base);
if ( flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file') );
disp(sprintf(' %% %s not found, creating',fname_fig));
 figure(1);clf;figmed;
 c_80s__ = colormap_80s; n_c_80s = size(c_80s__,1);
 linewidth_use = 3;
 markersize_use = 12;
 ncrop = 40;
 tmp_g = @(s) 0 ...
   | ~isempty(strfind(s,'X_2d_xcor_d0_a1t0122p20r')) ...
   | ~isempty(strfind(s,'X_2d_xcor_d0_a1t0122p25r')) ...
   | ~isempty(strfind(s,'X_2d_xcor_d0_a1t0122p30r')) ...
   ;
 %%%%%%%%%%%%%%%%;
 subplot(1,2,1);
 hold on;
 tmp_index_ = efind(cellfun(tmp_g,str_fname_align_a_k_Y_mat_o0to7_a_));
tmp_n = 0; tmp_c_avg = 0;
 for nl=0:numel(tmp_index_)-1;
 na = tmp_index_(1+nl);
tmp_c = corr_reco_vs_ampm_o0to7_a_(1+na);
 plot(1*n_octile+1,tmp_c,'ok','MarkerFaceColor','r','MarkerSize',markersize_use);
tmp_n = tmp_n + 1;
tmp_c_avg = tmp_c_avg + tmp_c;
 end;%for nl=0:numel(tmp_index_)-1;
tmp_c_avg = tmp_c_avg/max(1,tmp_n);
plot(1:n_octile+1,tmp_c_avg*ones(1+n_octile,1),'r-','LineWidth',linewidth_use);
%%%%%%%%;
tmp_n = 0; tmp_c_avg = 0;
 for noctile=0:n_octile-1;
if (noctile>0);
 tmp_index_ = efind(cellfun(tmp_g,str_fname_align_a_k_Y_mat_ya__{1+noctile}));
 for nl=0:numel(tmp_index_)-1;
 na = tmp_index_(1+nl);
tmp_c = corr_reco_vs_ampm_ya__{1+noctile}(1+na);
 plot(0*n_octile,tmp_c,'ok','MarkerFaceColor','c','MarkerSize',markersize_use);
tmp_n = tmp_n + 1;
tmp_c_avg = tmp_c_avg + tmp_c;
 end;%for nl=0:numel(tmp_index_)-1;
 end;%if (noctile>0);
end;% for noctile=0:n_octile-1;
tmp_c_avg = tmp_c_avg/max(1,tmp_n);
plot(0:n_octile+0,tmp_c_avg*ones(1+n_octile,1),'c-','LineWidth',linewidth_use);
 %%%%%%%%;
 for noctile=0:n_octile-1;
 %%%%%%%%;
 tmp_index_ = efind(cellfun(tmp_g,str_fname_align_a_k_Y_mat_oa__{1+noctile}));
 for nl=0:numel(tmp_index_)-1;
 na = tmp_index_(1+nl);
 plot(1+noctile,corr_reco_vs_ampm_oa__{1+noctile}(1+na),'ok','MarkerFaceColor','m','MarkerSize',markersize_use);
 end;%for nl=0:numel(tmp_index_)-1;
 %%%%%%%%;
 end;%for noctile=0:n_octile-1;
 hold off;
xlim([-1,n_octile+2]); set(gca,'XTick',0:n_octile+1,'XTickLabel',{'1K','1st','2nd','3rd','4th','5th','6th','7th','8th','8K'});
ylim([0.45,0.75]); ylabel('correlation'); grid on;
 title('corr_reco_vs_ampm_oa__','Interpreter','none');
 %%%%%%%%%%%%%%%%;
 subplot(1,2,2);
 hold on;
 tmp_index_ = efind(cellfun(tmp_g,str_fname_align_a_k_Y_mat_o0to7_a_));
tmp_n = 0; tmp_c_avg = 0;
 for nl=0:numel(tmp_index_)-1;
 na = tmp_index_(1+nl);
tmp_c = corr_full_reco_vs_crop_ampm_o0to7_xa__(1+ncrop,1+na);
 plot(1*n_octile+1,tmp_c,'ok','MarkerFaceColor','r','MarkerSize',markersize_use);
tmp_n = tmp_n + 1;
tmp_c_avg = tmp_c_avg + tmp_c;
 end;%for nl=0:numel(tmp_index_)-1;
tmp_c_avg = tmp_c_avg/max(1,tmp_n);
plot(1:n_octile+1,tmp_c_avg*ones(1+n_octile,1),'r-','LineWidth',linewidth_use);
%%%%%%%%;
tmp_n = 0; tmp_c_avg = 0;
 for noctile=0:n_octile-1;
 if (noctile>0);
 tmp_index_ = efind(cellfun(tmp_g,str_fname_align_a_k_Y_mat_ya__{1+noctile}));
 for nl=0:numel(tmp_index_)-1;
 na = tmp_index_(1+nl);
tmp_c = corr_full_reco_vs_crop_ampm_yxa___{1+noctile}(1+ncrop,1+na);
 plot(0*n_octile,tmp_c,'ok','MarkerFaceColor','c','MarkerSize',markersize_use);
tmp_n = tmp_n + 1;
tmp_c_avg = tmp_c_avg + tmp_c;
 end;%for nl=0:numel(tmp_index_)-1;
 end;%if (noctile>0);
end;% for noctile=0:n_octile-1;
tmp_c_avg = tmp_c_avg/max(1,tmp_n);
plot(0:n_octile+0,tmp_c_avg*ones(1+n_octile,1),'c-','LineWidth',linewidth_use);
 %%%%%%%%;
 for noctile=0:n_octile-1;
 %%%%%%%%;
 tmp_index_ = efind(cellfun(tmp_g,str_fname_align_a_k_Y_mat_oa__{1+noctile}));
 for nl=0:numel(tmp_index_)-1;
 na = tmp_index_(1+nl);
 plot(1+noctile,corr_full_reco_vs_crop_ampm_oxa___{1+noctile}(1+ncrop,1+na),'ok','MarkerFaceColor','m','MarkerSize',markersize_use);
 end;%for nl=0:numel(tmp_index_)-1;
%%%%%%%%;
 end;%for noctile=0:n_octile-1;
 hold off;
xlim([-1,n_octile+2]); set(gca,'XTick',0:n_octile+1,'XTickLabel',{'1K','1st','2nd','3rd','4th','5th','6th','7th','8th','8K'});
ylim([0.45,0.75]); ylabel('correlation'); grid on;
 title(sprintf('corr_full_reco_vs_crop_ampm_oxa___ ncrop %d',ncrop),'Interpreter','none');
 %%%%%%%%%%%%%%%%;
 sgtitle(sprintf('%s',dir_pm),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
%close(gcf);
end;%if ( flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file') );
if ( exist(sprintf('%s.jpg',fname_fig),'file') );
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file') );
 %}

flag_replot = 1;
fname_fig_pre = sprintf('%s/dir_ps1_x0to7_combine/dir_pm_jpg/test_pm_collect_excerpt_fsc_xa__1_FIGB',dir_base);
fname_fig_jpg = sprintf('%s.jpg',fname_fig);
fname_fig_eps = sprintf('%s.eps',fname_fig);
if ( flag_replot | ~exist(fname_fig_jpg,'file') );
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1);clf;figsml;
c_80s__ = colormap_80s; n_c_80s = size(c_80s__,1);
linewidth_use = 3;
markersize_use = 12;
fontsize_use = 12;
ncrop = 40;
%tmp_g = @(s) 0 ...
%  | ~isempty(strfind(s,'X_2d_xcor_d0_a1t0122p20r')) ...
%  | ~isempty(strfind(s,'X_2d_xcor_d0_a1t0122p25r')) ...
%  | ~isempty(strfind(s,'X_2d_xcor_d0_a1t0122p30r')) ...
%  ;
tmp_g = @(s) 0 ...
  | ~isempty(strfind(s,'X_2d_xcor_d0_a1t0122p20r0')) ...
  | ~isempty(strfind(s,'X_2d_xcor_d0_a1t0122p20r1')) ...
  ;
%%%%%%%%%%%%%%%%;
subplot(1,1,1);
set(gca,'FontSize',fontsize_use);
hold on;
%%%%%%%%;
tmp_n = 0; tmp_c_avg = 0;
tmp_c__ = [];
for noctile=0:n_octile-1;
if (noctile>0);
tmp_index_ = efind(cellfun(tmp_g,str_fname_align_a_k_Y_mat_ya__{1+noctile})); n_l = numel(tmp_index_);
tmp_c_ = zeros(n_l,1);
for nl=0:n_l-1;
na = tmp_index_(1+nl);
tmp_c_(1+nl) = corr_full_reco_vs_crop_ampm_yxa___{1+noctile}(1+ncrop,1+na);
end;%for nl=0:n_l-1;
tmp_c = median(tmp_c_);
for nl=0:n_l-1;
plot(0*n_octile,tmp_c,'ok','MarkerFaceColor','m','MarkerSize',markersize_use);
end;%for nl=0:n_l-1;
tmp_c__ = [tmp_c__;tmp_c_];
end;%if (noctile>0);
end;% for noctile=0:n_octile-1;
tmp_c_avg = mean(tmp_c__);
plot(0:n_octile+0,tmp_c_avg*ones(1+n_octile,1),'m-','LineWidth',linewidth_use);
%%%%%%%%;
for noctile=0:n_octile-1;
%%%%%%%%;
tmp_index_ = efind(cellfun(tmp_g,str_fname_align_a_k_Y_mat_oa__{1+noctile})); n_l = numel(tmp_index_);
tmp_c_ = zeros(n_l,1);
for nl=0:n_l-1;
na = tmp_index_(1+nl);
tmp_c_(1+nl) = corr_full_reco_vs_crop_ampm_oxa___{1+noctile}(1+ncrop,1+na);
end;%for nl=0:n_l-1;
tmp_c = max(tmp_c_); if (noctile==n_octile-1); tmp_c = min(tmp_c_); end;
plot(1+noctile-0.0625,tmp_c,'ok','MarkerFaceColor',0.95*[1.0,0.65,1.0],'MarkerSize',markersize_use);
%%%%;
tmp_index_ = efind(cellfun(tmp_g,str_fname_align_a_k_Y_mat_ua__{1+noctile})); n_l = numel(tmp_index_);
tmp_c_ = zeros(n_l,1);
for nl=0:n_l-1;
na = tmp_index_(1+nl);
tmp_c_(1+nl) = corr_full_reco_vs_crop_ampm_uxa___{1+noctile}(1+ncrop,1+na);
end;%for nl=0:n_l-1;
tmp_c = max(tmp_c_); if (noctile==n_octile-1); tmp_c = min(tmp_c_); end;
plot(1+noctile+0.0625,tmp_c,'ok','MarkerFaceColor',0.65*[1.0,0.05,1.0],'MarkerSize',markersize_use);
%%%%%%%%;
end;%for noctile=0:n_octile-1;
hold off;
%xlim([-1,n_octile+2]); set(gca,'XTick',0:n_octile+1,'XTickLabel',{'1K','1st','2nd','3rd','4th','5th','6th','7th','8th','8K'});
xlim([-1,n_octile+1]); set(gca,'XTick',0:n_octile+0,'XTickLabel',{'1K','1st','2nd','3rd','4th','5th','6th','7th','8th'});
ylim([0.55,0.75]); ylabel('correlation (masked)'); grid on;
title(sprintf('corr_full_reco_vs_crop_ampm_?xa___ ncrop %d',ncrop),'Interpreter','none');
%%%%%%%%%%%%%%%%;
sgtitle(fname_fig_pre,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
%close(gcf);
end;%if ( flag_replot | ~exist(fname_fig_jpg,'file') );
if ( exist(fname_fig_jpg,'file') );
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(fname_fig_jpg,'file') );

disp('returning');return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if nargin<1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;


str_thisfunction = 'test_pm_collect_excerpt_fsc_xa__1';

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); dir_pm=[]; end; na=na+1;
if (nargin<1+na); dir_pm_bkp=[]; end; na=na+1;
if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_recalc'); parameter.flag_recalc=0; end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
if ~isfield(parameter,'Pixel_Spacing'); parameter.Pixel_Spacing=1; end;
flag_recalc = parameter.flag_recalc;
flag_verbose = parameter.flag_verbose;
Pixel_Spacing = parameter.Pixel_Spacing;
if isempty(dir_pm); dir_pm = pwd; end;
if isempty(dir_pm_bkp); dir_pm_bkp = dir_pm; end;

if (flag_verbose); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

%%%%%%%%;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); string_root = 'home'; end;
if (strcmp(platform,'eval1')); string_root = 'home'; end;
if (strcmp(platform,'rusty')); string_root = 'mnt/home'; end;
%%%%%%%%;
dir_base = sprintf('/%s/rangan/dir_cryoem',string_root);

dir_pm_mat = sprintf('%s_mat',dir_pm);
dir_pm_bkp_mat = sprintf('%s_mat',dir_pm_bkp);
if (flag_verbose); disp(sprintf(' %% dir_pm_mat: %s',dir_pm_mat)); end;
if (flag_verbose); disp(sprintf(' %% dir_pm_bkp_mat: %s',dir_pm_bkp_mat)); end;

%%%%%%%%;
% find cropped (i.e., masked) fsc between ground-truth reconstruction and published molecule. ;
%%%%%%%%;
fsc_ampm_ka__ = [];
corr_base_vs_ampm_a_ = [];
corr_reco_vs_ampm_a_ = [];
fsc_crop_ampm_kxa___ = [];
corr_crop_base_vs_crop_ampm_xa__ = [];
corr_full_base_vs_crop_ampm_xa__ = [];
corr_crop_reco_vs_crop_ampm_xa__ = [];
corr_full_reco_vs_crop_ampm_xa__ = [];
k_Ainv_p_r_ = [];
k_Ainv_p_r_max = [];
kinv_A_p_r_ = [];
%%%%%%%%;
fname_pre = sprintf('%s/test_pm_collect_fsc_crop_ampm_',dir_pm_mat);
[fname_skip,fname_mat] = open_fname_tmp(fname_pre);
if ( flag_recalc>0 | ~fname_skip );
disp(sprintf(' %% %s not found, creating',fname_mat));
tmp_fname = sprintf('%s/M_k_p__.mat',dir_pm_mat);
if ~exist(tmp_fname,'file'); tmp_fname = sprintf('%s/M_k_p__.mat',dir_pm_bkp_mat); end;
if ~exist(tmp_fname,'file'); disp(sprintf(' %% Warning, %s not found in %s',tmp_fname,str_thisfunction)); return; end;
if (flag_verbose); disp(sprintf(' %% %s found, loading',tmp_fname)); end;
tmp_ = load(tmp_fname,'n_x_M_u','x_c_0_');
tmp_flag_exist=0;
if (isfield(tmp_,'n_x_M_u')); tmp_flag_exist=1; n_x_M_u = tmp_.n_x_M_u; end;
if (isfield(tmp_,'x_c_0_')); tmp_flag_exist=1; n_x_M_u = numel(tmp_.x_c_0_); end;
if ~tmp_flag_exist; disp(sprintf(' %% Warning, %s not found in test_pm_collect_excerpt_fsc_xa__0','n_x_M_u')); return; end;
clear tmp_;
%%%%%%%%;
tmp_fname = sprintf('%s/a_k_p_quad_.mat',dir_pm_mat);
if ~exist(tmp_fname,'file'); tmp_fname = sprintf('%s/a_k_p_quad_.mat',dir_pm_bkp_mat); end;
if ~exist(tmp_fname,'file'); disp(sprintf(' %% Warning, %s not found in %s',tmp_fname,str_thisfunction)); return; end;
if (flag_verbose); disp(sprintf(' %% %s found, loading',tmp_fname)); end;
tmp_ = load(tmp_fname);
k_p_r_max = tmp_.k_p_r_max;
k_eq_d = tmp_.k_eq_d;
n_k_all = tmp_.n_k_all;
n_k_all_csum_ = tmp_.n_k_all_csum_;
k_p_r_all_ = tmp_.k_p_r_all_;
k_p_azimu_b_all_ = tmp_.k_p_azimu_b_all_;
k_p_polar_a_all_ = tmp_.k_p_polar_a_all_;
weight_3d_k_all_ = tmp_.weight_3d_k_all_;
weight_shell_k_ = tmp_.weight_shell_k_;
n_k_p_r = tmp_.n_k_p_r;
k_p_r_ = tmp_.k_p_r_;
weight_3d_k_p_r_ = tmp_.weight_3d_k_p_r_;
k_c_0_all_ = tmp_.k_c_0_all_;
k_c_1_all_ = tmp_.k_c_1_all_;
k_c_2_all_ = tmp_.k_c_2_all_;
a_k_p_quad_ = tmp_.a_k_p_quad_;
a_x_u_reco_ = tmp_.a_x_u_reco_;
clear tmp_;
%%%%%%%%;
tmp_fname = sprintf('%s/a_k_Y_quad_.mat',dir_pm_mat);
if ~exist(tmp_fname,'file'); tmp_fname = sprintf('%s/a_k_Y_quad_.mat',dir_pm_bkp_mat); end;
if ~exist(tmp_fname,'file'); disp(sprintf(' %% Warning, %s not found in %s',tmp_fname,str_thisfunction)); return; end;
if (flag_verbose); disp(sprintf(' %% %s found, loading',tmp_fname)); end;
tmp_ = load(tmp_fname,'l_max_','a_k_Y_quad_');
l_max_ = tmp_.l_max_;
a_k_Y_quad_ = tmp_.a_k_Y_quad_;
clear tmp_;
%%%%%%%%;
tmp_flag_exist=0;
tmp_fname = sprintf('%s/a_x_u_base_.mat',dir_pm_mat);
if ~exist(tmp_fname,'file'); tmp_fname = sprintf('%s/a_x_u_base_.mat',dir_pm_bkp_mat); end;
if exist(tmp_fname,'file');
tmp_flag_exist=1;
if (flag_verbose); disp(sprintf(' %% %s found, loading',tmp_fname)); end;
tmp_ = load(tmp_fname,'half_diameter_x_c','x_u_0_','x_u_1_','x_u_2_','n_x_u_pack','a_x_u_base_');
a_x_u_base_ = tmp_.a_x_u_base_;
end;%if exist(tmp_fname,'file');
tmp_fname = sprintf('%s/a_x_u_pack_.mat',dir_pm_mat);
if ~exist(tmp_fname,'file'); tmp_fname = sprintf('%s/a_x_u_pack_.mat',dir_pm_bkp_mat); end;
if exist(tmp_fname,'file');
tmp_flag_exist=1;
if (flag_verbose); disp(sprintf(' %% %s found, loading',tmp_fname)); end;
tmp_ = load(tmp_fname,'half_diameter_x_c','x_u_0_','x_u_1_','x_u_2_','n_x_u_pack','a_x_u_pack_');
a_x_u_base_ = tmp_.a_x_u_pack_;
end;%if exist(tmp_fname,'file');
if ~tmp_flag_exist; disp(sprintf(' %% Warning, %s not found in %s',tmp_fname,str_thisfunction)); return; end;
half_diameter_x_c = tmp_.half_diameter_x_c;
x_u_0_ = tmp_.x_u_0_; x_u_1_ = tmp_.x_u_1_; x_u_2_ = tmp_.x_u_2_; n_x_u_pack = tmp_.n_x_u_pack;
[x_u_0___,x_u_1___,x_u_2___] = ndgrid(x_u_0_,x_u_1_,x_u_2_); n_xxx_u = n_x_u_pack^3;
x_u_r___ = sqrt(x_u_0___.^2 + x_u_1___.^2 + x_u_2___.^2);
clear tmp_;
%%%%%%%%;
tmp_flag_exist = 0;
tmp_fname = sprintf('%s/a_k_Y_0lsq_reco_.mat',dir_pm_mat);
if ~exist(tmp_fname,'file'); tmp_fname = sprintf('%s/a_k_Y_0lsq_reco_.mat',dir_pm_bkp_mat); end;
if exist(tmp_fname,'file');
tmp_flag_exist = 1;
if (flag_verbose); disp(sprintf(' %% %s found, loading',tmp_fname)); end;
tmp_ = load(tmp_fname);
b_k_Y_reco_ = tmp_.a_k_Y_0lsq_reco_;
if isfield(tmp_,'a_k_p_0lsq_reco_'); b_k_p_reco_ = tmp_.a_k_p_0lsq_reco_; end;
if isfield(tmp_,'a_k_p_0lsq_reco'); b_k_p_reco_ = tmp_.a_k_p_0lsq_reco; end; %<-- typo in earlier version. ;
b_x_u_reco_ = tmp_.a_x_u_0lsq_reco_;
end;%if exist(tmp_fname,'file');
tmp_fname = sprintf('%s/c_k_Y_.mat',dir_pm_mat);
if ~exist(tmp_fname,'file'); tmp_fname = sprintf('%s/c_k_Y_.mat',dir_pm_bkp_mat); end;
if exist(tmp_fname,'file');
tmp_flag_exist = 1;
if (flag_verbose); disp(sprintf(' %% %s found, loading',tmp_fname)); end;
tmp_ = load(tmp_fname,'c_k_Y_reco_','X_best_reco','c_k_p_reco_','c_x_u_reco_');
b_k_Y_reco_ = tmp_.c_k_Y_reco_;
b_k_p_reco_ = tmp_.c_k_p_reco_;
b_x_u_reco_ = tmp_.c_x_u_reco_;
end;%if exist(tmp_fname,'file');
if ~tmp_flag_exist; disp(sprintf(' %% Warning, %s not found in %s',tmp_fname,str_thisfunction)); return; end;
clear tmp_;
%%%%;
tmp_X_2d_xcor_d0_a1t_align_ = ls(sprintf('%s/X_2d_xcor_d0_a1t*_align_a_k_Y_.mat',dir_pm_mat));
tmp_index_start_ = strfind(tmp_X_2d_xcor_d0_a1t_align_,sprintf('/%s',string_root))-1;
tmp_index_final_ = strfind(tmp_X_2d_xcor_d0_a1t_align_,'.mat')+4-1;
n_X_2d_xcor_d0_a1t_align = min(numel(tmp_index_start_),numel(tmp_index_final_));
%%%%;
a_k_Y_true_yk_ = a_k_Y_quad_;
str_fname_align_a_k_Y_mat_a_ = cell(n_X_2d_xcor_d0_a1t_align,1);
fsc_ampm_ka__ = zeros(n_k_p_r,n_X_2d_xcor_d0_a1t_align);
corr_base_vs_ampm_a_ = zeros(n_X_2d_xcor_d0_a1t_align,1);
corr_reco_vs_ampm_a_ = zeros(n_X_2d_xcor_d0_a1t_align,1);
fsc_crop_ampm_kxa___ = zeros(n_k_p_r,n_x_u_pack,n_X_2d_xcor_d0_a1t_align);
corr_crop_base_vs_crop_ampm_xa__ = zeros(n_x_u_pack,n_X_2d_xcor_d0_a1t_align);
corr_full_base_vs_crop_ampm_xa__ = zeros(n_x_u_pack,n_X_2d_xcor_d0_a1t_align);
corr_crop_reco_vs_crop_ampm_xa__ = zeros(n_x_u_pack,n_X_2d_xcor_d0_a1t_align);
corr_full_reco_vs_crop_ampm_xa__ = zeros(n_x_u_pack,n_X_2d_xcor_d0_a1t_align);
k_Ainv_p_r_ = (2*k_p_r_)/(n_x_M_u * Pixel_Spacing);
k_Ainv_p_r_max = (2*k_p_r_max)/(n_x_M_u * Pixel_Spacing);
kinv_A_p_r_ = 1./max(1e-12,k_Ainv_p_r_);
sample_sphere_k_eq_d = 1/(2*pi);
%%%%%%%%;
for na=0:n_X_2d_xcor_d0_a1t_align-1;
if (flag_verbose); disp(sprintf(' %% na %d/%d',na,n_X_2d_xcor_d0_a1t_align)); end;
tmp_fname_align_a_k_Y_mat = tmp_X_2d_xcor_d0_a1t_align_(1+tmp_index_start_(1+na)+0:1+tmp_index_final_(1+na)-1);
str_fname_align_a_k_Y_mat_a_{1+na} = tmp_fname_align_a_k_Y_mat;
if (tmp_fname_align_a_k_Y_mat(1)==sprintf('\n')); tmp_fname_align_a_k_Y_mat = tmp_fname_align_a_k_Y_mat(2:end); end;
tmp_fname_a_k_Y_mat = tmp_fname_align_a_k_Y_mat;
tmp_ij = strfind(tmp_fname_a_k_Y_mat,'_align_a_k_Y_.mat');
tmp_fname_a_k_Y_mat = tmp_fname_a_k_Y_mat([1:tmp_ij,tmp_ij+7:end]);
tmp_fname_fsc_crop_kx_pre = tmp_fname_align_a_k_Y_mat;
tmp_ij = strfind(tmp_fname_fsc_crop_kx_pre,'_align_a_k_Y_.mat');
tmp_fname_fsc_crop_kx_pre = tmp_fname_fsc_crop_kx_pre([1:tmp_ij-1]);
tmp_fname_fsc_crop_kx_pre = sprintf('%s_fsc_crop_kx__',tmp_fname_fsc_crop_kx_pre);
tmp_fname_fsc_crop_kx_mat = sprintf('%s.mat',tmp_fname_fsc_crop_kx_pre);
if (  exist(tmp_fname_align_a_k_Y_mat,'file') &  exist(tmp_fname_a_k_Y_mat,'file') );
if (flag_recalc>1 | ~exist(tmp_fname_fsc_crop_kx_mat,'file'));
if (flag_verbose); disp(sprintf(' %% %s not found, creating',tmp_fname_fsc_crop_kx_mat)); end;
tmp_XB_align_a_k_Y_ = load(tmp_fname_align_a_k_Y_mat);
tmp_XB_a_k_Y_ = load(tmp_fname_a_k_Y_mat);
%%%%;
tmp_X_best = tmp_XB_align_a_k_Y_.X_best_(end);
tmp_flag_flip = tmp_XB_align_a_k_Y_.X_best_flag_flip_(end);
tmp_polar_a_best = tmp_XB_align_a_k_Y_.polar_a_best_(end);
tmp_azimu_b_best = tmp_XB_align_a_k_Y_.azimu_b_best_(end);
tmp_gamma_z_best = tmp_XB_align_a_k_Y_.gamma_z_best_(end);
tmp_delta_x_best = tmp_XB_align_a_k_Y_.delta_best__(1+0,end);
tmp_delta_y_best = tmp_XB_align_a_k_Y_.delta_best__(1+1,end);
tmp_a_k_Y_ampm_ = tmp_XB_a_k_Y_.a_k_Y_reco_;
%%%%;
tmp_t = tic();
[ ... 
 tmp_b_k_Y_ampm_ ...
] = ...
spharm_register_and_rotate_2( ...
 n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_true_yk_ ...
,tmp_a_k_Y_ampm_ ...
,0 ...
,tmp_X_best ...
,tmp_flag_flip ...
,tmp_polar_a_best ...
,tmp_azimu_b_best ...
,tmp_gamma_z_best ...
,tmp_delta_x_best ...
,tmp_delta_y_best ...
);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% spharm_register_and_rotate_2: %0.6fs',tmp_t)); end;
%%%%;
tmp_t = tic();
[ ... 
  tmp_b_x_u_ampm_ ...
] = ...
convert_spharm_to_x_c_3( ...
 sample_sphere_k_eq_d ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,tmp_b_k_Y_ampm_ ...
,half_diameter_x_c ...
,n_x_u_pack ...
);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% convert_spharm_to_x_c_3: %0.6fs',tmp_t)); end;
%%%%%%%%;
n_crop = n_x_u_pack;
corr_base_vs_ampm = real(corr( a_x_u_base_(:) , tmp_b_x_u_ampm_(:) ));
corr_reco_vs_ampm = real(corr( a_x_u_reco_(:) , tmp_b_x_u_ampm_(:) ));
corr_full_base_vs_crop_ampm_x_ = zeros(n_crop,1);
corr_crop_base_vs_crop_ampm_x_ = zeros(n_crop,1);
corr_full_reco_vs_crop_ampm_x_ = zeros(n_crop,1);
corr_crop_reco_vs_crop_ampm_x_ = zeros(n_crop,1);
for ncrop=0:n_crop-1;
r_crop = 1.0*ncrop/(n_crop-1);
tmp_m_x_u_ = reshape(x_u_r___<=r_crop,[n_xxx_u,1]);
corr_crop_base_vs_crop_ampm_x_(1+ncrop) = real(corr( a_x_u_base_(:).*tmp_m_x_u_ , tmp_b_x_u_ampm_(:).*tmp_m_x_u_ ));
corr_full_base_vs_crop_ampm_x_(1+ncrop) = real(corr( a_x_u_base_(:)             , tmp_b_x_u_ampm_(:).*tmp_m_x_u_ ));
corr_crop_reco_vs_crop_ampm_x_(1+ncrop) = real(corr( a_x_u_reco_(:).*tmp_m_x_u_ , tmp_b_x_u_ampm_(:).*tmp_m_x_u_ ));
corr_full_reco_vs_crop_ampm_x_(1+ncrop) = real(corr( a_x_u_reco_(:)             , tmp_b_x_u_ampm_(:).*tmp_m_x_u_ ));
end;%for ncrop=0:n_crop-1;
%%%%%%%%;
fsc_crop_ampm_k_ = zeros(n_k_p_r,1);
fsc_crop_ampm_kx__ = zeros(n_k_p_r,n_x_u_pack);
%%%%%%%%;
parameter_fsc = struct('type','parameter');
parameter_fsc.flag_register = 0;
[ ...
 parameter_fsc ...
,fsc_ampm_k_ ...
] = ...
fsc_from_a_k_Y_0( ...
 parameter_fsc ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,0 ...
,l_max_ ...
,a_k_Y_quad_ ...
,tmp_b_k_Y_ampm_ ...
);
%%%%;
tmp_t = tic();
for ncrop=0:n_crop-1;
if (flag_verbose>0); disp(sprintf(' %% ncrop %d/%d',ncrop,n_crop)); end;
r_crop = 1.0*ncrop/(n_crop-1);
tmp_m_x_u_ = reshape(x_u_r___<=r_crop,[n_xxx_u,1]);
[ ...
 c_k_p_base_ ...
] = ...
convert_x_c_to_k_p_1( ...
 flag_verbose ...
,n_k_all ...
,weight_3d_k_all_ ...
,k_c_0_all_ ...
,k_c_1_all_ ...
,k_c_2_all_ ...
,k_p_r_max ...
,n_xxx_u ...
,[] ...
,x_u_0___ ...
,x_u_1___ ...
,x_u_2___ ...
,half_diameter_x_c ...
,a_x_u_base_(:).*tmp_m_x_u_ ...
);
[ ...
 d_k_p_ampm_ ...
] = ...
convert_x_c_to_k_p_1( ...
 flag_verbose ...
,n_k_all ...
,weight_3d_k_all_ ...
,k_c_0_all_ ...
,k_c_1_all_ ...
,k_c_2_all_ ...
,k_p_r_max ...
,n_xxx_u ...
,[] ...
,x_u_0___ ...
,x_u_1___ ...
,x_u_2___ ...
,half_diameter_x_c ...
,tmp_b_x_u_ampm_(:).*tmp_m_x_u_ ...
);
parameter_fsc = struct('type','parameter');
[ ...
 parameter_fsc ...
,fsc_crop_ampm_k_ ...
] = ...
fsc_from_a_k_p_0( ...
 parameter_fsc ...
,n_k_all ...
,n_k_p_r ...
,n_k_all_csum_ ...
,weight_3d_k_all_ ...
,c_k_p_base_ ...
,d_k_p_ampm_ ...
);
fsc_crop_ampm_kx__(:,1+ncrop) = fsc_crop_ampm_k_;
end;%for ncrop=0:n_crop-1;
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% calculating fsc_crop_ampm_kx__: %0.6fs',tmp_t)); end;
%%%%%%%%;
save(tmp_fname_fsc_crop_kx_mat ...
     ,'fsc_ampm_k_' ...
     ,'corr_base_vs_ampm' ...
     ,'corr_reco_vs_ampm' ...
     ,'fsc_crop_ampm_kx__' ...
     ,'corr_crop_base_vs_crop_ampm_x_' ...
     ,'corr_full_base_vs_crop_ampm_x_' ...
     ,'corr_crop_reco_vs_crop_ampm_x_' ...
     ,'corr_full_reco_vs_crop_ampm_x_' ...
     );
end;% if (flag_recalc>1 | ~exist(tmp_fname_fsc_crop_kx_mat,'file'));
end;%if (  exist(tmp_fname_align_a_k_Y_mat,'file') &  exist(tmp_fname_a_k_Y_mat,'file') );
fsc_ampm_k_ = zeros(n_k_p_r,1);
corr_base_vs_ampm = zeros(1,1);
corr_reco_vs_ampm = zeros(1,1);
fsc_crop_ampm_kx__ = zeros(n_k_p_r,n_x_u_pack);
corr_crop_base_vs_crop_ampm_x_ = zeros(n_x_u_pack,1);
corr_full_base_vs_crop_ampm_x_ = zeros(n_x_u_pack,1);
corr_crop_reco_vs_crop_ampm_x_ = zeros(n_x_u_pack,1);
corr_full_reco_vs_crop_ampm_x_ = zeros(n_x_u_pack,1);
if  exist(tmp_fname_fsc_crop_kx_mat,'file');
if (flag_verbose); disp(sprintf(' %% %s found, not creating',tmp_fname_fsc_crop_kx_mat)); end;
load(tmp_fname_fsc_crop_kx_mat);
end;%if  exist(tmp_fname_fsc_crop_kx_mat,'file');
fsc_ampm_ka__(:,1+na) = fsc_ampm_k_;
corr_base_vs_ampm_a_(1+na) = corr_base_vs_ampm;
corr_reco_vs_ampm_a_(1+na) = corr_reco_vs_ampm;
fsc_crop_ampm_kxa___(:,:,1+na) = fsc_crop_ampm_kx__;
corr_crop_base_vs_crop_ampm_xa__(:,1+na) = corr_crop_base_vs_crop_ampm_x_;
corr_full_base_vs_crop_ampm_xa__(:,1+na) = corr_full_base_vs_crop_ampm_x_;
corr_crop_reco_vs_crop_ampm_xa__(:,1+na) = corr_crop_reco_vs_crop_ampm_x_;
corr_full_reco_vs_crop_ampm_xa__(:,1+na) = corr_full_reco_vs_crop_ampm_x_;
end;%for na=0:n_X_2d_xcor_d0_a1t_align-1;
%%%%;
save( ...
 fname_mat ...
,'str_fname_align_a_k_Y_mat_a_' ...
,'fsc_ampm_ka__' ...
,'corr_base_vs_ampm_a_' ...
,'corr_reco_vs_ampm_a_' ...
,'fsc_crop_ampm_kxa___' ...
,'corr_crop_base_vs_crop_ampm_xa__' ...
,'corr_full_base_vs_crop_ampm_xa__' ...
,'corr_crop_reco_vs_crop_ampm_xa__' ...
,'corr_full_reco_vs_crop_ampm_xa__' ...
,'k_Ainv_p_r_' ...
,'k_Ainv_p_r_max' ...
,'kinv_A_p_r_' ...
);
%%%%%%%%;
close_fname_tmp(fname_pre);
end;%if ( flag_recalc>0 | ~fname_skip );
%%%%%%%%;
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));

if (flag_verbose); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
