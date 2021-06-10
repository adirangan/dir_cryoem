function helper_heterogeneity_L1_0(dir_base,infix,n_k_max);
%dir_base = './dir_k20_dF_s2_nF_cF_tF_hTa60b40c0d0_rTl2i7_M2048_S2048_s0_MS_r1';
%n_k_max = 20; infix = 'est'; 
% try: ;
%{
  clear; clf; helper_heterogeneity_L1_0('./dir_k20_dF_s2_nF_cF_tF_hTa60b40c0d0_rTl2i7_M2048_S2048_s3162_MS_r1','tru',20);
  clear; clf; helper_heterogeneity_L1_0('./dir_k20_dF_s2_nF_cF_tF_hTa60b40c0d0_rTl2i7_M2048_S2048_s3162_MS_r1','est',20);
  clear; clf; helper_heterogeneity_L1_0('./dir_k20_dF_s2_nF_cF_tF_hTa60b40c0d0_rTl2i7_M2048_S2048_s1000_MS_r1','tru',20);
  clear; clf; helper_heterogeneity_L1_0('./dir_k20_dF_s2_nF_cF_tF_hTa60b40c0d0_rTl2i7_M2048_S2048_s1000_MS_r1','est',20);
  clear; clf; helper_heterogeneity_L1_0('./dir_k20_dF_s2_nF_cF_tF_hTa60b40c0d0_rTl2i7_M2048_S2048_s316_MS_r1','tru',20);
  clear; clf; helper_heterogeneity_L1_0('./dir_k20_dF_s2_nF_cF_tF_hTa60b40c0d0_rTl2i7_M2048_S2048_s316_MS_r1','est',20);
  clear; clf; helper_heterogeneity_L1_0('./dir_k20_dF_s2_nF_cF_tF_hTa60b40c0d0_rTl2i7_M2048_S2048_s100_MS_r1','tru',20);
  clear; clf; helper_heterogeneity_L1_0('./dir_k20_dF_s2_nF_cF_tF_hTa60b40c0d0_rTl2i7_M2048_S2048_s100_MS_r1','est',20);
  clear; clf; helper_heterogeneity_L1_0('./dir_k20_dF_s2_nF_cF_tF_hTa60b40c0d0_rTl2i7_M2048_S2048_s32_MS_r1','tru',20);
  clear; clf; helper_heterogeneity_L1_0('./dir_k20_dF_s2_nF_cF_tF_hTa60b40c0d0_rTl2i7_M2048_S2048_s32_MS_r1','est',20);
  clear; clf; helper_heterogeneity_L1_0('./dir_k20_dF_s2_nF_cF_tF_hTa60b40c0d0_rTl2i7_M2048_S2048_s10_MS_r1','tru',20);
  clear; clf; helper_heterogeneity_L1_0('./dir_k20_dF_s2_nF_cF_tF_hTa60b40c0d0_rTl2i7_M2048_S2048_s10_MS_r1','est',20);
  %}
flag_print=1; 
fname = sprintf('%s/Heterogeneity_L1_%s_%d',dir_base,infix,n_k_max); 
nk_cut = 8;
n_k_ = 3:n_k_max;
nbins = 128; hbins = linspace(-2.5,+2.5,nbins);
%%%%%%%%%%%%%%%% ;
clear M_ I_ AL_ AN_ ;
M_ = cell(length(n_k_),1); I_ = cell(length(n_k_),1);
for nk=1:length(n_k_);
n_k = n_k_(nk);
Iname = sprintf('%s/I_model_sample_%s_%d_.mda',dir_base,infix,n_k);
Mname = sprintf('%s/M_residual_loading_%s_%d_.mda',dir_base,infix,n_k);
Itmp = MDA_read_i4(Iname);
Mtmp = MDA_read_c16(Mname);
Mtmp = Mtmp(1,:); Mtmp = Mtmp / std(Mtmp);
I_{nk} = transpose(Itmp);
M_{nk} = Mtmp;
L_ = cell2mat(M_(nk));
[U,S,V] = svds(L_); Vtmp = V(:,1); Vtmp = Vtmp/std(Vtmp);
Atmp = get_auc(real(Vtmp(find(Itmp==0))),real(Vtmp(find(Itmp==1))));
if (Atmp<0.5);
Vtmp = -Vtmp;
Atmp = get_auc(real(Vtmp(find(Itmp==0))),real(Vtmp(find(Itmp==1))));
end;%if (Atmp<0);
AL_(nk) = Atmp;
disp(sprintf(' %% n_k %d %s %d %s %d --> AL %0.2f',n_k,Iname,length(Itmp),Mname,length(Vtmp),Atmp));
subplot(1,6,[1:2]);
hold on;
h_A = hist(real(Vtmp(find(Itmp==0))),hbins); 
h_B = hist(real(Vtmp(find(Itmp==1))),hbins); 
h_X = h_A+h_B; h_A = h_A/max(h_X); h_B = h_B/max(h_X);
plot(hbins,n_k + 0*hbins,'k-');
stairs(hbins,n_k + h_A,'k-','LineWidth',2);
stairs(hbins,n_k + h_B,'r-','LineWidth',2);
hold off;
if nk<nk_cut; N_ = cell2mat(M_(nk)); else; N_ = cell2mat(M_(nk_cut:nk)); end;
[U,S,V] = svds(N_); Vtmp = V(:,1); Vtmp = Vtmp/std(Vtmp);
Atmp = get_auc(real(Vtmp(find(Itmp==0))),real(Vtmp(find(Itmp==1))));
if (Atmp<0.5);
Vtmp = -Vtmp;
Atmp = get_auc(real(Vtmp(find(Itmp==0))),real(Vtmp(find(Itmp==1))));
end;%if (Atmp<0);
AN_(nk) = Atmp;
disp(sprintf(' %% n_k %d %s %d %s %d --> AN %0.2f',n_k,Iname,length(Itmp),Mname,length(Vtmp),Atmp));
subplot(1,6,[4:5]);
hold on;
h_A = hist(real(Vtmp(find(Itmp==0))),hbins); 
h_B = hist(real(Vtmp(find(Itmp==1))),hbins); 
h_X = h_A+h_B; h_A = h_A/max(h_X); h_B = h_B/max(h_X);
plot(hbins,n_k + 0*hbins,'k-');
stairs(hbins,n_k + h_A,'k-','LineWidth',2);
stairs(hbins,n_k + h_B,'r-','LineWidth',2);
hold off;
end;%for n_k_cur = 3:n_k_max;
subplot(1,6,[1:2]);
xlim([min(hbins),max(hbins)]); ylim([min(n_k_),max(n_k_)+1]);
title('Loadings A=black, B=red');
ylabel('n_k');
xlabel('loading histogram');
subplot(1,6,[4:5]);
xlim([min(hbins),max(hbins)]); ylim([min(n_k_),max(n_k_)+1]);
title(sprintf('Accumulated (%d+) A=black, B=red',nk_cut));
ylabel('n_k');
xlabel('accumulated loading histogram');
%%%%%%%%%%%%%%%% ;
subplot(1,6,3);
plot(0*n_k_+0.5,n_k_,'k-',AL_,n_k_,'kx',AL_,n_k_,'ko'); 
xlim([0.4,1]); ylim([min(n_k_),max(n_k_)+1]);
xlabel('Auc'); title('Auc');
subplot(1,6,6);
plot(0*n_k_+0.5,n_k_,'k-',AN_,n_k_,'kx',AN_,n_k_,'ko'); 
xlim([0.4,1]); ylim([min(n_k_),max(n_k_)+1]);
xlabel('Auc'); title('Auc');
set(gcf,'Position',1+[0,0,1024,1024+256]);
if flag_print;
print('-djpeg',sprintf('%s.jpg',fname)); pause(0.5);
print('-painters','-depsc',sprintf('%s.eps',fname)); pause(0.5);
end;%if flag_print;
