%%%%%%%%;
% Here we simply reorder the T1.star-file to reduce the n_CTF_rank of the first 1000 images. ;
%%%%%%%%;
clear;

platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;

flag_recalc = 0;
flag_replot = 0;
flag_center = 1;
flag_invert = 0;
tolerance_master = 1e-2;
nf=0;

dir_pm = sprintf('/%s/rangan/dir_cryoem/dir_p28hRPT1/dir_pm',string_root);
if (~exist(sprintf('%s_mat',dir_pm),'dir')); disp(sprintf(' %% mkdir %s_mat',dir_pm)); mkdir(sprintf('%s_mat',dir_pm)); end;
if (~exist(sprintf('%s_jpg',dir_pm),'dir')); disp(sprintf(' %% mkdir %s_jpg',dir_pm)); mkdir(sprintf('%s_jpg',dir_pm)); end;
dir_relion = sprintf('/%s/rangan/dir_cryoem/dir_p28hRPT1/dir_relion',string_root);
if (~exist(sprintf('%s_mat',dir_relion),'dir')); disp(sprintf(' %% mkdir %s_mat',dir_relion)); mkdir(sprintf('%s_mat',dir_relion)); end;
if (~exist(sprintf('%s_jpg',dir_relion),'dir')); disp(sprintf(' %% mkdir %s_jpg',dir_relion)); mkdir(sprintf('%s_jpg',dir_relion)); end;
string_rusty_root = 'mnt/home';
dir_relion_bin = sprintf('/%s/rangan/relion/build/bin',string_rusty_root);
dir_data_star = sprintf('/%s/rangan/dir_cryoem/dir_p28hRP',string_root);
Pixel_Spacing = 0.98; %<-- in angstroms, from https://www.ebi.ac.uk/empiar/entry/10091/ ; 
fname_nopath_volume = 'emd_8674.map'; %<-- class T1. ;
%%%%%%%%;
% all classes and subclasses. ;
%%%%%%%%;
fname_nopath_volume_ = ...
{ ...
 'emd_8674.map'... %<-- class T1. ;
};
n_volume = numel(fname_nopath_volume_);
flag_het = 0; if (n_volume> 1); flag_het = 1; end;
fname_nopath_star = 'T1.star';

fname_star = sprintf('%s/%s',dir_data_star,fname_nopath_star);
n_header = 11;
n_Line_T1 = 39520 + n_header;
[star_blockNames,star_blockData,ok]=ReadStarFile_0(fname_star,n_Line_T1);
[~,index_rlnDefocus_] = sort(star_blockData{1}.rlnDefocusU,'ascend'); index_rlnDefocus_ = index_rlnDefocus_-1;

n_M_star = size(star_blockData{1}.rlnImageName,1);
rlnDefocus_threshold = 12000;

%%%%%%%%;
% Now write a reorganized version of the star file. ;
%%%%%%%%;
fname_nopath_star_reorder = 'T1s.star';
fname_star_reorder = sprintf('%s/%s',dir_data_star,fname_nopath_star_reorder);
fp = fopen(fname_star);
fp_reorder = fopen(fname_star_reorder,'w');
for nheader=0:n_header-1;
fprintf(fp_reorder,fgets(fp));
end;%for nheader=0:n_header-1;
fclose(fp);
for nM_star=0:n_M_star-1;
if (mod(nM_star,1024)==0); disp(sprintf(' %% nM_star %d/%d',nM_star,n_M_star)); end;
tmp_index = index_rlnDefocus_(1+nM_star);
rlnDefocus = star_blockData{1}.rlnDefocusU(1+tmp_index);
if (rlnDefocus< rlnDefocus_threshold);
disp(sprintf(' %% Warning, skipping line %d, rlnDefocus %0.2f',nM_star,rlnDefocus));
end;%if (rlnDefocus< rlnDefocus_threshold);
if (rlnDefocus>=rlnDefocus_threshold);
fprintf( ...
 fp_reorder ...
,'%.6f   %.6f   %.6f   %.6f   %.6f   %.6f   %s\n' ...
,star_blockData{1}.rlnVoltage(1+tmp_index) ...
,star_blockData{1}.rlnDefocusU(1+tmp_index) ...
,star_blockData{1}.rlnDefocusV(1+tmp_index) ...
,star_blockData{1}.rlnDefocusAngle(1+tmp_index) ...
,star_blockData{1}.rlnSphericalAberration(1+tmp_index) ...
,star_blockData{1}.rlnAmplitudeContrast(1+tmp_index) ...
,star_blockData{1}.rlnImageName{1+tmp_index} ...
);
end;%if (rlnDefocus>=rlnDefocus_threshold);
end;%for nM_star=0:n_M_star-1;
fclose(fp_reorder);

disp('returning');return;
