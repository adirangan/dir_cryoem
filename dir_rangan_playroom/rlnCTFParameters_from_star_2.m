function ...
[ ...
 index_nCTF_from_nM_ ...
,index_nM_from_nCTF_ ...
,Voltage_CTF_ ...
,DefocusU_CTF_ ...
,DefocusV_CTF_ ...
,DefocusAngle_CTF_ ...
,SphericalAberration_CTF_ ...
,AmplitudeContrast_CTF_ ...
] = ...
rlnCTFParameters_from_star_2( ...
 dir_base ...
,fname_nopath_star ...
,n_image_0in ...
,nimage_start_0in ...
,nimage_final_0in ...
);

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
if (verbose); disp(sprintf(' %% [entering rlnCTFParameters_from_star_2] n_image_0in %d: [%d,%d]',n_image_0in,nimage_start_0in,nimage_final_0in)); end;
fname_star = sprintf('%s/%s',dir_base,fname_nopath_star);
if (verbose); disp(sprintf(' %% fname_star: %s',fname_star)); end;
if (~exist(fname_star,'file')); disp(sprintf(' %% %s not found',fname_star)); end;
if ( exist(fname_star,'file'));
if (n_image_0in> 0); [star_blockNames,star_blockData,ok]=ReadStarFile_0(fname_star,n_image_0in+nimage_start_0in+1024); end;
if (n_image_0in<=0); [star_blockNames,star_blockData,ok]=ReadStarFile(fname_star); end;

n_block = numel(star_blockNames);
%%%%%%%%;
if (n_block>2);
disp(sprintf(' %% Warning, n_block %d in rlnCTFParameters_from_star_2',n_block));
end;%if (n_block>2);
%%%%%%%%;
if (n_block==1);
if (verbose); disp(sprintf(' %% n_block %d',n_block)); end;
if isfield(star_blockData{1},'rlnImageName'); n_image_use = min(n_image_0in,numel(star_blockData{1}.rlnImageName)-nimage_start_0in); end;
if isfield(star_blockData{1},'rlnMicrographName'); n_image_use = min(n_image_0in,numel(star_blockData{1}.rlnMicrographName)-nimage_start_0in); end;
nimage_start_use = nimage_start_0in;
nimage_final_use = nimage_start_0in + n_image_use - 1;
if (verbose); disp(sprintf(' %% n_image_use %d: [%d,%d]',n_image_use,nimage_start_use,nimage_final_use)); end;
index_nM_use_ = [nimage_start_use:nimage_final_use];
if isfield(star_blockData{1},'rlnImageName'); rlnLabelName_ = star_blockData{1}.rlnImageName(1+index_nM_use_); end;
if isfield(star_blockData{1},'rlnMicrographName'); rlnLabelName_ = star_blockData{1}.rlnMicrographName(1+index_nM_use_); end;
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
if isfield(star_blockData{1+index_particles},'rlnImageName'); n_image_use = min(n_image_0in,numel(star_blockData{1+index_particles}.rlnImageName)-nimage_start_0in); end;
if isfield(star_blockData{1+index_particles},'rlnMicrographName'); n_image_use = min(n_image_0in,numel(star_blockData{1+index_particles}.rlnMicrographName)-nimage_start_0in); end;
nimage_start_use = nimage_start_0in;
nimage_final_use = nimage_start_0in + n_image_use - 1;
if (verbose); disp(sprintf(' %% n_image_use %d: [%d,%d]',n_image_use,nimage_start_use,nimage_final_use)); end;
index_nM_use_ = [nimage_start_use:nimage_final_use];
if isfield(star_blockData{1+index_particles},'rlnImageName'); rlnLabelName_ = star_blockData{1+index_particles}.rlnImageName(1+index_nM_use_); end;
if isfield(star_blockData{1+index_particles},'rlnMicrographName'); rlnLabelName_ = star_blockData{1+index_particles}.rlnMicrographName(1+index_nM_use_); end;
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
end;%if ( exist(fname_star,'file'));
if (verbose); disp(sprintf(' %% [finished rlnCTFParameters_from_star_2]')); end;
