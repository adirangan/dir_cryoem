%%%%%%%%;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
%%%%%%%%;

tmp_dir_m = sprintf('/%s/rangan/dir_cryoem/dir_rangan_playroom',string_root);

for tmp_str_infix_ = { ...
  'a_restore_C2M0' ...
  ,'a_restore_C1M1' ...
  ,'a_restore_C0M2' ...
  ,'dtau_a_restore_C2M0' ...
  ,'dtau_a_restore_C1M1' ...
  ,'dtau_a_restore_C0M2' ...
  ,'dtau_dtau_a_restore_C2M0' ...
  ,'dtau_dtau_a_restore_C1M1' ...
  ,'dtau_dtau_a_restore_C0M2' ...
};
tmp_str_infix = tmp_str_infix_{1};
fname_m = sprintf('%s/ddssnll_helper_convert_%s_3.m',tmp_dir_m,tmp_str_infix);
disp(sprintf(' %% generating %s',fname_m));
fp=fopen(fname_m,'w');
fprintf(fp,'if  isempty(%s_k_Y_quad_yk__) & ~isempty(%s_k_Y_quad_yk_);\n',tmp_str_infix,tmp_str_infix);
fprintf(fp,'%s_k_Y_quad_yk__ = local_yk__from_yk_(n_k_p_r,l_max_,%s_k_Y_quad_yk_);\n',tmp_str_infix,tmp_str_infix);
fprintf(fp,'end;%%if  isempty(%s_k_Y_quad_yk__) & ~isempty(%s_k_Y_quad_yk_);\n',tmp_str_infix,tmp_str_infix);
fprintf(fp,'%%%%%%%%;\n');
fprintf(fp,'if  isempty(%s_k_Y_quad_yk_) & ~isempty(%s_k_Y_quad_yk__);\n',tmp_str_infix,tmp_str_infix);
fprintf(fp,'%s_k_Y_quad_yk_ = local_yk_from_yk__(n_k_p_r,l_max_,%s_k_Y_quad_yk__);\n',tmp_str_infix,tmp_str_infix);
fprintf(fp,'end;%%if  isempty(%s_k_Y_quad_yk_) & ~isempty(%s_k_Y_quad_yk__);\n',tmp_str_infix,tmp_str_infix);
fprintf(fp,'%%%%%%%%;\n');
fprintf(fp,'if  isempty(%s_k_p_quad_) & ~isempty(%s_k_Y_quad_yk_);\n',tmp_str_infix,tmp_str_infix);
fprintf(fp,'tmp_yk_ = %s_k_Y_quad_yk_; tmp_str = ''%s_k_Y_quad_yk_'';\n',tmp_str_infix,tmp_str_infix);
fprintf(fp,'%%%%%%%%;\n');
fprintf(fp,'tmp_t = tic();\n');
fprintf(fp,'%%%%%%%%;\n');
fprintf(fp,'[ ...\n');
fprintf(fp,' tmp_quad_ ...\n');
fprintf(fp,',Ylm_uklma___ ...\n');
fprintf(fp,',k_p_azimu_b_sub_uka__ ...\n');
fprintf(fp,',k_p_polar_a_sub_uka__ ...\n');
fprintf(fp,',l_max_uk_ ...\n');
fprintf(fp,',index_nu_n_k_per_shell_from_nk_p_r_ ...\n');
fprintf(fp,',index_k_per_shell_uka__ ...\n');
fprintf(fp,'] = ...\n');
fprintf(fp,'convert_spharm_to_k_p_4( ...\n');
fprintf(fp,' 0*flag_verbose ...\n');
fprintf(fp,',n_k_all ...\n');
fprintf(fp,',n_k_all_csum_ ...\n');
fprintf(fp,',k_p_r_all_ ...\n');
fprintf(fp,',k_p_azimu_b_all_ ...\n');
fprintf(fp,',k_p_polar_a_all_ ...\n');
fprintf(fp,',weight_3d_k_all_ ...\n');
fprintf(fp,',weight_shell_k_ ...\n');
fprintf(fp,',n_k_p_r ...\n');
fprintf(fp,',k_p_r_ ...\n');
fprintf(fp,',weight_3d_k_p_r_ ...\n');
fprintf(fp,',l_max_ ...\n');
fprintf(fp,',tmp_yk_ ...\n');
fprintf(fp,',Ylm_uklma___ ...\n');
fprintf(fp,',k_p_azimu_b_sub_uka__ ...\n');
fprintf(fp,',k_p_polar_a_sub_uka__ ...\n');
fprintf(fp,',l_max_uk_ ...\n');
fprintf(fp,',index_nu_n_k_per_shell_from_nk_p_r_ ...\n');
fprintf(fp,',index_k_per_shell_uka__ ...\n');
fprintf(fp,');\n');
fprintf(fp,'tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf('' %%%% %s: convert_spharm_to_k_p_4: time %%0.2fs'',tmp_t)); end;\n',tmp_str_infix);
fprintf(fp,'%s_k_p_quad_ = tmp_quad_;\n',tmp_str_infix);
fprintf(fp,'end;%%if  isempty(%s_k_p_quad_) & ~isempty(%s_k_Y_quad_yk_);\n',tmp_str_infix,tmp_str_infix);
fprintf(fp,'%%%%%%%%;\n');
fclose(fp);
end;%for;
