function ...
[ ...
 M_x_c___ ...
,index_nCTF_from_nM_ ...
,index_nM_from_nCTF_ ...
,Voltage_CTF_ ...
,DefocusU_CTF_ ...
,DefocusV_CTF_ ...
,DefocusAngle_CTF_ ...
,SphericalAberration_CTF_ ...
,AmplitudeContrast_CTF_ ...
] = ...
rlnImageName_from_star_2( ...
 dir_base ...
,fname_nopath_star ...
,n_image_0in ...
,nimage_start_0in ...
,nimage_final_0in ...
);

%%%%%%%%;
if nargin<1;
dir_base = '/data/rangan/dir_cryoem/dir_precatalytic_spliceosome';
fname_nopath_star = 'consensus_data.star';
nimage_start_base = 3768;
n_image_0in = 1024;
nimage_start_0in = nimage_start_base;
[ ...
 M_x_c_A___ ...
,index_nCTF_from_nM_A_ ...
,index_nM_from_nCTF_A_ ...
,Voltage_CTF_A_ ...
,DefocusU_CTF_A_ ...
,DefocusV_CTF_A_ ...
,DefocusAngle_CTF_A_ ...
,SphericalAberration_CTF_A_ ...
,AmplitudeContrast_CTF_A_ ...
] = ...
rlnImageName_from_star_2( ...
 dir_base ...
,fname_nopath_star ...
,n_image_0in ...
,nimage_start_0in ...
);
nimage_start_0in = nimage_start_base + 512;
[ ...
 M_x_c_B___ ...
,index_nCTF_from_nM_B_ ...
,index_nM_from_nCTF_B_ ...
,Voltage_CTF_B_ ...
,DefocusU_CTF_B_ ...
,DefocusV_CTF_B_ ...
,DefocusAngle_CTF_B_ ...
,SphericalAberration_CTF_B_ ...
,AmplitudeContrast_CTF_B_ ...
] = ...
rlnImageName_from_star_2( ...
 dir_base ...
,fname_nopath_star ...
,n_image_0in ...
,nimage_start_0in ...
);
n_sum=0; sum_error=0;
for nM_A=0:n_image_0in-1;
nM_B = nM_A - 512;
if (nM_B>=0);
n_sum = n_sum+1;
nCTF_A = index_nCTF_from_nM_A_(1+nM_A);
nCTF_B = index_nCTF_from_nM_B_(1+nM_B);
sum_error = sum_error + fnorm(M_x_c_B___(:,:,1+nM_B) - M_x_c_A___(:,:,1+nM_A));
sum_error = sum_error + fnorm(Voltage_CTF_B_(1+nCTF_B) - Voltage_CTF_A_(1+nCTF_A));
sum_error = sum_error + fnorm(DefocusU_CTF_B_(1+nCTF_B) - DefocusU_CTF_A_(1+nCTF_A));
sum_error = sum_error + fnorm(DefocusV_CTF_B_(1+nCTF_B) - DefocusV_CTF_A_(1+nCTF_A));
sum_error = sum_error + fnorm(DefocusAngle_CTF_B_(1+nCTF_B) - DefocusAngle_CTF_A_(1+nCTF_A));
sum_error = sum_error + fnorm(SphericalAberration_CTF_B_(1+nCTF_B) - SphericalAberration_CTF_A_(1+nCTF_A));
sum_error = sum_error + fnorm(AmplitudeContrast_CTF_B_(1+nCTF_B) - AmplitudeContrast_CTF_A_(1+nCTF_A));
end;%if (nM_B>=0);
end;%for nM_A=0:n_image_0in-1;
disp(sprintf(' %% n_sum %d: sum_error: %0.16f',n_sum,sum_error));
disp('returning'); return;
end;%if nargin<1;
%%%%%%%%;

na=0;
if (nargin<1+na); dir_base=[]; end; na=na+1;
if (nargin<1+na); fname_nopath_star=[]; end; na=na+1;
if (nargin<1+na); n_image_0in=[]; end; na=na+1;
if (nargin<1+na); nimage_start_0in=[]; end; na=na+1;
if (nargin<1+na); nimage_final_0in=[]; end; na=na+1;
%%%%%%%%;
if isempty(dir_base); dir_base = pwd; end;
if isempty(fname_nopath_star); fname_nopath_star = 'star.star'; end;
if isempty(n_image_0in); n_image_0in = 0; end;
if isempty(nimage_start_0in); nimage_start_0in = 0; end;
if isempty(nimage_final_0in); nimage_final_0in = 0; end;
%%%%%%%%;
if (n_image_0in> 0); nimage_final_0in = nimage_start_0in + n_image_0in - 1; end;
if (nimage_final_0in>=nimage_start_0in); n_image_0in = 1+nimage_final_0in-nimage_start_0in; end;
%%%%%%%%;

verbose=1;
if (verbose); disp(sprintf(' %% [entering rlnImageName_from_star_2] n_image_0in %d: [%d,%d]',n_image_0in,nimage_start_0in,nimage_final_0in)); end;
fname_star = sprintf('%s/%s',dir_base,fname_nopath_star);
if (verbose); disp(sprintf(' %% fname_star: %s',fname_star)); end;
if (~exist(fname_star,'file')); disp(sprintf(' %% %s not found',fname_star)); end;
if ( exist(fname_star,'file'));
if (n_image_0in> 0); [star_blockNames,star_blockData,ok]=ReadStarFile_0(fname_star,n_image_0in+nimage_start_0in+1024); end;
if (n_image_0in<=0); [star_blockNames,star_blockData,ok]=ReadStarFile(fname_star); end;

n_block = numel(star_blockNames);
%%%%%%%%;
if (n_block>2);
disp(sprintf(' %% Warning, n_block %d in rlnImageName_from_star_2',n_block));
end;%if (n_block>2);
%%%%%%%%;
if (n_block==1);
if (verbose); disp(sprintf(' %% n_block %d',n_block)); end;
n_image_use = min(n_image_0in,numel(star_blockData{1}.rlnImageName)-nimage_start_0in);
nimage_start_use = nimage_start_0in;
nimage_final_use = nimage_start_0in + n_image_use - 1;
if (verbose); disp(sprintf(' %% n_image_use %d: [%d,%d]',n_image_use,nimage_start_use,nimage_final_use)); end;
index_nM_use_ = [nimage_start_use:nimage_final_use];
rlnImageName_ = star_blockData{1}.rlnImageName(1+index_nM_use_);
Voltage_M_ = zeros(n_image_use,1); if isfield(star_blockData{1},'rlnVoltage'); Voltage_M_ = star_blockData{1}.rlnVoltage(1+index_nM_use_); end;
DefocusU_M_ = zeros(n_image_use,1); if isfield(star_blockData{1},'rlnDefocusU'); DefocusU_M_ = star_blockData{1}.rlnDefocusU(1+index_nM_use_); end;
DefocusV_M_ = zeros(n_image_use,1); if isfield(star_blockData{1},'rlnDefocusV'); DefocusV_M_ = star_blockData{1}.rlnDefocusV(1+index_nM_use_); end;
DefocusAngle_M_ = zeros(n_image_use,1); if isfield(star_blockData{1},'rlnDefocusAngle'); DefocusAngle_M_ = star_blockData{1}.rlnDefocusAngle(1+index_nM_use_); end;
SphericalAberration_M_ = zeros(n_image_use,1); if isfield(star_blockData{1},'rlnSphericalAberration'); SphericalAberration_M_ = star_blockData{1}.rlnSphericalAberration(1+index_nM_use_); end;
AmplitudeContrast_M_ = zeros(n_image_use,1); if isfield(star_blockData{1},'rlnAmplitudeContrast'); AmplitudeContrast_M_ = star_blockData{1}.rlnAmplitudeContrast(1+index_nM_use_); end;
end;%if (n_block==1);
%%%%%%%%;
if (n_block==2);
if (verbose); disp(sprintf(' %% n_block %d',n_block)); end;
index_optics = efind(cellfun(@(x) strcmp(x,'data_optics'),star_blockNames));
index_particles = efind(cellfun(@(x) strcmp(x,'data_particles'),star_blockNames));
n_image_use = min(n_image_0in,numel(star_blockData{1+index_particles}.rlnImageName)-nimage_start_0in);
nimage_start_use = nimage_start_0in;
nimage_final_use = nimage_start_0in + n_image_use - 1;
if (verbose); disp(sprintf(' %% n_image_use %d: [%d,%d]',n_image_use,nimage_start_use,nimage_final_use)); end;
index_nM_use_ = [nimage_start_use:nimage_final_use];
rlnImageName_ = star_blockData{1+index_particles}.rlnImageName(1+index_nM_use_);
Voltage_M_ = zeros(n_image_use,1); if isfield(star_blockData{1+index_particles},'rlnVoltage'); Voltage_M_(:) = star_blockData{1+index_particles}.rlnVoltage(1+index_nM_use_); end; if isfield(star_blockData{1+index_optics},'rlnVoltage'); Voltage_M_(:) = star_blockData{1+index_optics}.rlnVoltage(1); end;
DefocusU_M_ = zeros(n_image_use,1); if isfield(star_blockData{1+index_particles},'rlnDefocusU'); DefocusU_M_(:) = star_blockData{1+index_particles}.rlnDefocusU(1+index_nM_use_); end; if isfield(star_blockData{1+index_optics},'rlnDefocusU'); DefocusU_M_(:) = star_blockData{1+index_optics}.rlnDefocusU(1); end;
DefocusV_M_ = zeros(n_image_use,1); if isfield(star_blockData{1+index_particles},'rlnDefocusV'); DefocusV_M_(:) = star_blockData{1+index_particles}.rlnDefocusV(1+index_nM_use_); end; if isfield(star_blockData{1+index_optics},'rlnDefocusV'); DefocusV_M_(:) = star_blockData{1+index_optics}.rlnDefocusV(1); end;
DefocusAngle_M_ = zeros(n_image_use,1); if isfield(star_blockData{1+index_particles},'rlnDefocusAngle'); DefocusAngle_M_(:) = star_blockData{1+index_particles}.rlnDefocusAngle(1+index_nM_use_); end; if isfield(star_blockData{1+index_optics},'rlnDefocusAngle'); DefocusAngle_M_(:) = star_blockData{1+index_optics}.rlnDefocusAngle(1); end;
SphericalAberration_M_ = zeros(n_image_use,1); if isfield(star_blockData{1+index_particles},'rlnSphericalAberration'); SphericalAberration_M_(:) = star_blockData{1+index_particles}.rlnSphericalAberration(1+index_nM_use_); end; if isfield(star_blockData{1+index_optics},'rlnSphericalAberration'); SphericalAberration_M_(:) = star_blockData{1+index_optics}.rlnSphericalAberration(1); end;
AmplitudeContrast_M_ = zeros(n_image_use,1); if isfield(star_blockData{1+index_particles},'rlnAmplitudeContrast'); AmplitudeContrast_M_(:) = star_blockData{1+index_particles}.rlnAmplitudeContrast(1+index_nM_use_); end; if isfield(star_blockData{1+index_optics},'rlnAmplitudeContrast'); AmplitudeContrast_M_(:) = star_blockData{1+index_optics}.rlnAmplitudeContrast(1); end;
end;%if (n_block==2);
%%%%%%%%;
CTF_M_ = [Voltage_M_ , DefocusU_M_ , DefocusV_M_ , DefocusAngle_M_ , SphericalAberration_M_ , AmplitudeContrast_M_];
[u_CTF_,index_nM_from_nCTF_,index_nCTF_from_nM_] = unique(CTF_M_,'rows'); n_u_CTF = size(u_CTF_,1);
index_nM_from_nCTF_ = index_nM_from_nCTF_ - 1;
index_nCTF_from_nM_ = index_nCTF_from_nM_ - 1;
if (verbose); disp(sprintf(' %% found %d unique CTFs',n_u_CTF)); end;
Voltage_CTF_ = u_CTF_(:,1+0);
DefocusU_CTF_ = u_CTF_(:,1+1);
DefocusV_CTF_ = u_CTF_(:,1+2);
DefocusAngle_CTF_ = u_CTF_(:,1+3);
SphericalAberration_CTF_ = u_CTF_(:,1+4);
AmplitudeContrast_CTF_ = u_CTF_(:,1+5);
if (verbose>2);
for nu_CTF=0:n_u_CTF-1;
disp(sprintf(' %% nu_CTF %d/%d: Voltage %0.2f DefocusU %0.2f DefocusV %0.2f DefocusAngle %0.2f SphericalAberration %0.2f AmplitudeContrast %0.2f',nu_CTF,n_u_CTF,Voltage_CTF_(1+nu_CTF),DefocusU_CTF_(1+nu_CTF),DefocusV_CTF_(1+nu_CTF),DefocusAngle_CTF_(1+nu_CTF),SphericalAberration_CTF_(1+nu_CTF),AmplitudeContrast_CTF_(1+nu_CTF)));
end;%for nu_CTF=0:n_u_CTF-1;
end;%if (verbose>2);
%%%%%%%%;
if (verbose>1);
for nimage=0:n_image_use-1;
disp(rlnImageName_{1+nimage});
end;%for nimage=0:n_image_use-1;
end;%if (verbose>1);
%%%%%%%%;
rlnImageFrame_ = zeros(n_image_use,1);
rlnImageLname_ = cell(n_image_use,1);
for nimage=0:n_image_use-1;
rlnImageName = rlnImageName_{1+nimage};
tmp_ij = strfind(rlnImageName,'@');
rlnImageFrame = 0;
rlnImageLname = '';
if ~isempty(tmp_ij);
rlnImageFrame = str2num(rlnImageName(1:tmp_ij-1));
rlnImageLname = rlnImageName(tmp_ij+1:end);
end;%if ~isempty(tmp_ij);
rlnImageFrame_(1+nimage) = rlnImageFrame;
rlnImageLname_{1+nimage} = rlnImageLname;
end;%for nimage=0:n_image_use-1;
if (verbose>2); plot(rlnImageFrame_,'.'); xlabel('nimage'); ylabel('rlnImageFrame'); end;
%%%%%%%%;
u_Lname_ = unique(rlnImageLname_); n_u_Lname = numel(u_Lname_);
Lname_index_from_nu_Lname__ = cell(n_u_Lname,1);
Frame_from_nu_Lname__ = cell(n_u_Lname,1);
n_Frame_from_nu_Lname_ = zeros(n_u_Lname,1);
for nu_Lname=0:n_u_Lname-1;
u_Lname = u_Lname_{1+nu_Lname};
Lname_index_from_nu_Lname_ = efind(strcmp(rlnImageLname_,u_Lname));
Lname_index_from_nu_Lname__{1+nu_Lname} = Lname_index_from_nu_Lname_;
Frame_from_nu_Lname_ = rlnImageFrame_(1+Lname_index_from_nu_Lname_);
Frame_from_nu_Lname__{1+nu_Lname} = Frame_from_nu_Lname_;
n_Frame_from_nu_Lname = numel(Frame_from_nu_Lname_);
n_Frame_from_nu_Lname_(1+nu_Lname) = n_Frame_from_nu_Lname;
if (verbose); disp(sprintf(' %% u_Lname %s: %d frames',u_Lname,n_Frame_from_nu_Lname)); end;
end;%for nu_Lname=0:n_u_Lname-1;
%%%%%%%%;
M_x_c___ = [];
rlnImageName_out_ = cell(n_image_use,1);
%%%%%%%%;
for nu_Lname=0:n_u_Lname-1;
u_Lname = u_Lname_{1+nu_Lname};
Lname_index_from_nu_Lname_ = Lname_index_from_nu_Lname__{1+nu_Lname};
Frame_from_nu_Lname_ = Frame_from_nu_Lname__{1+nu_Lname};
n_Frame_from_nu_Lname = n_Frame_from_nu_Lname_(1+nu_Lname);
Frame_min = min(Frame_from_nu_Lname_);
Frame_max = max(Frame_from_nu_Lname_);
fname_mrc = sprintf('%s/%s',dir_base,u_Lname);
if (~exist(fname_mrc,'file')); disp(sprintf(' %% %s not found',fname_mrc)); end;
if ( exist(fname_mrc,'file')); 
tmp_M_x_c___ = cast(ReadMRC(fname_mrc,Frame_min,1+Frame_max-Frame_min),'double');
if (isempty(M_x_c___)); n_x_0 = size(tmp_M_x_c___,1); n_x_1 = size(tmp_M_x_c___,2); M_x_c___ = zeros(n_x_0,n_x_1,n_image_use); end;
for nFrame_from_nu_Lname=0:n_Frame_from_nu_Lname-1;
Lname_index = Lname_index_from_nu_Lname_(1+nFrame_from_nu_Lname);
Frame_out = Frame_from_nu_Lname_(1+nFrame_from_nu_Lname);
Frame_0in = 1+Frame_out-Frame_min;
tmp_M_x_c__ = tmp_M_x_c___(:,:,Frame_0in);
M_x_c___(:,:,1+Lname_index) = tmp_M_x_c__;
rlnImageName_out_{1+Lname_index} = sprintf('%.6d@%s',Frame_out,u_Lname);
end;%for nFrame_from_nu_Lname=0:n_Frame_from_nu_Lname-1;
end;%if ( exist(fname_mrc,'file')); 
end;%for nu_Lname=0:n_u_Lname-1;
%%%%%%%%;
for nimage=0:n_image_use-1;
if (verbose>2); disp(sprintf(' %% nimage %d/%d: rlnImageName_out %s rlnImageName %s',nimage,n_image_use,rlnImageName_out_{1+nimage},rlnImageName_{1+nimage})); end;
tmp_str_out_ = rlnImageName_out_{1+nimage};
tmp_str_0in_ = rlnImageName_{1+nimage};
tmp_frame_out = str2num(tmp_str_out_(1:strfind(tmp_str_out_,'@')-1));
tmp_frame_0in = str2num(tmp_str_0in_(1:strfind(tmp_str_0in_,'@')-1));
tmp_name_out = tmp_str_out_(strfind(tmp_str_out_,'@')+1:end);
tmp_name_0in = tmp_str_0in_(strfind(tmp_str_0in_,'@')+1:end);
assert(tmp_frame_out==tmp_frame_0in);
assert(strcmp(tmp_name_out,tmp_name_0in));
end;%for nimage=0:n_image_use-1;
%%%%%%%%;
end;%if ( exist(fname_star,'file'));
if (verbose); disp(sprintf(' %% [finished rlnImageName_from_star_2]')); end;
