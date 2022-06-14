function test_pm_23_clean_helper_0( ...
			 verbose ...
			 ,flag_clean ...
			 ,flag_clean_really ...
			 ,dir_pm ...
			 ,str_xfix ...
			 );
for flag_xcor_vs_Memp=0:1;
if flag_xcor_vs_Memp==0; str_infix = 'X_2d_Memp_d1_'; end;
if flag_xcor_vs_Memp==1; str_infix = 'X_2d_xcor_d0_'; end;
XA_fname_mat = sprintf('%s_mat/%s%s.mat',dir_pm,str_infix,str_xfix);
XA_fname_align_a_k_Y_pre = sprintf('%s_mat/%s%s_align_a_k_Y_',dir_pm,str_infix,str_xfix);
XA_fname_align_a_k_Y_mat = sprintf('%s.mat',XA_fname_align_a_k_Y_pre);
XA_fname_align_a_k_Y_jpg = sprintf('%s.jpg',XA_fname_align_a_k_Y_pre);
XA_fname_snapshot_pre = sprintf('%s_mat/%s%s_snapshot',dir_pm,str_infix,str_xfix);
XA_fname_snapshot_mat = sprintf('%s.mat',XA_fname_snapshot_pre);
XA_fname_snapshot_jpg = sprintf('%s.jpg',XA_fname_snapshot_pre);
XA_fname_compare_image_rank_pre = sprintf('%s_mat/%s%s_compare_image_rank',dir_pm,str_infix,str_xfix);
XA_fname_compare_image_rank_mat = sprintf('%s.mat',XA_fname_compare_image_rank_pre);
XA_fname_compare_image_rank_jpg = sprintf('%s.jpg',XA_fname_compare_image_rank_pre);
XA_fname_align_crop_a_k_Y_mat = sprintf('%s_mat/%s%s_align_crop_a_k_Y_.mat',dir_pm,str_infix,str_xfix);
XA_fname_fsc_crop_kx__mat = sprintf('%s_mat/%s%s_fsc_crop_kx__.mat',dir_pm,str_infix,str_xfix);
if (~exist(XA_fname_mat,'file')); disp(sprintf(' %% %s not found',XA_fname_mat)); end;%if (~exist(XA_fname_mat,'file'));
if ( exist(XA_fname_mat,'file'));
if (verbose>2);
if (~exist(XA_fname_align_a_k_Y_mat,'file')); disp(sprintf(' %% %s not found',XA_fname_align_a_k_Y_mat)); end;
if (~exist(XA_fname_align_a_k_Y_jpg,'file')); disp(sprintf(' %% %s not found',XA_fname_align_a_k_Y_jpg)); end;
if (~exist(XA_fname_snapshot_mat,'file')); disp(sprintf(' %% %s not found',XA_fname_snapshot_mat)); end;
if (~exist(XA_fname_snapshot_jpg,'file')); disp(sprintf(' %% %s not found',XA_fname_snapshot_jpg)); end;
if (~exist(XA_fname_compare_image_rank_mat,'file')); disp(sprintf(' %% %s not found',XA_fname_compare_image_rank_mat)); end;
if (~exist(XA_fname_compare_image_rank_jpg,'file')); disp(sprintf(' %% %s not found',XA_fname_compare_image_rank_jpg)); end;
if (~exist(XA_fname_align_crop_a_k_Y_mat,'file')); disp(sprintf(' %% %s not found',XA_fname_align_crop_a_k_Y_mat)); end;
if (~exist(XA_fname_fsc_crop_kx__mat,'file')); disp(sprintf(' %% %s not found',XA_fname_fsc_crop_kx__mat)); end;
end;%if (verbose>2);
%%%%;
if flag_clean;
if ( exist(XA_fname_align_a_k_Y_mat,'file')); disp(sprintf(' %% %s found, deleting',XA_fname_align_a_k_Y_mat));
if flag_clean_really; delete(XA_fname_align_a_k_Y_mat); end; end;
if ( exist(XA_fname_align_a_k_Y_jpg,'file')); disp(sprintf(' %% %s found, deleting',XA_fname_align_a_k_Y_jpg));
if flag_clean_really; delete(XA_fname_align_a_k_Y_jpg); end; end;
if ( exist(XA_fname_snapshot_mat,'file')); disp(sprintf(' %% %s found, deleting',XA_fname_snapshot_mat));
if flag_clean_really; delete(XA_fname_snapshot_mat); end; end;
if ( exist(XA_fname_snapshot_jpg,'file')); disp(sprintf(' %% %s found, deleting',XA_fname_snapshot_jpg));
if flag_clean_really; delete(XA_fname_snapshot_jpg); end; end;
if ( exist(XA_fname_compare_image_rank_mat,'file')); disp(sprintf(' %% %s found, deleting',XA_fname_compare_image_rank_mat));
if flag_clean_really; delete(XA_fname_compare_image_rank_mat); end; end;
if ( exist(XA_fname_compare_image_rank_jpg,'file')); disp(sprintf(' %% %s found, deleting',XA_fname_compare_image_rank_jpg));
if flag_clean_really; delete(XA_fname_compare_image_rank_jpg); end; end;
if ( exist(XA_fname_align_crop_a_k_Y_mat,'file')); disp(sprintf(' %% %s found, deleting',XA_fname_align_crop_a_k_Y_mat));
if flag_clean_really; delete(XA_fname_align_crop_a_k_Y_mat); end; end;
if ( exist(XA_fname_fsc_crop_kx__mat,'file')); disp(sprintf(' %% %s found, deleting',XA_fname_fsc_crop_kx__mat));
if flag_clean_really; delete(XA_fname_fsc_crop_kx__mat); end; end;
end;%if flag_clean;
%%%%;
end;%if ( exist(XA_fname_mat,'file'));
end;%for flag_xcor_vs_Memp=0:1;
