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
rlnImageName_from_star_0( ...
 dir_base ...
,fname_nopath_star ...
,n_image ...
);
if (nargin<3); n_image=[]; end;
if isempty(n_image); n_image = 0; end;

verbose=1;
if (verbose); disp(sprintf(' %% [entering rlnImageName_from_star_0]')); end;
fname_star = sprintf('%s/%s',dir_base,fname_nopath_star);
if (verbose); disp(sprintf(' %% fname_star: %s',fname_star)); end;
if (~exist(fname_star,'file')); disp(sprintf(' %% %s not found',fname_star)); end;
if ( exist(fname_star,'file'));
if (n_image> 0); [star_blockNames,star_blockData,ok]=ReadStarFile_0(fname_star,n_image+1024); end;
if (n_image<=0); [star_blockNames,star_blockData,ok]=ReadStarFile(fname_star); end;
n_image = min(n_image,numel(star_blockData{1}.rlnImageName));
rlnImageName_ = star_blockData{1}.rlnImageName(1:n_image);
Voltage_M_ = zeros(n_image,1); if isfield(star_blockData{1},'rlnVoltage'); Voltage_M_ = star_blockData{1}.rlnVoltage(1:n_image); end;
DefocusU_M_ = zeros(n_image,1); if isfield(star_blockData{1},'rlnDefocusU'); DefocusU_M_ = star_blockData{1}.rlnDefocusU(1:n_image); end;
DefocusV_M_ = zeros(n_image,1); if isfield(star_blockData{1},'rlnDefocusV'); DefocusV_M_ = star_blockData{1}.rlnDefocusV(1:n_image); end;
DefocusAngle_M_ = zeros(n_image,1); if isfield(star_blockData{1},'rlnDefocusAngle'); DefocusAngle_M_ = star_blockData{1}.rlnDefocusAngle(1:n_image); end;
SphericalAberration_M_ = zeros(n_image,1); if isfield(star_blockData{1},'rlnSphericalAberration'); SphericalAberration_M_ = star_blockData{1}.rlnSphericalAberration(1:n_image); end;
AmplitudeContrast_M_ = zeros(n_image,1); if isfield(star_blockData{1},'rlnAmplitudeContrast'); AmplitudeContrast_M_ = star_blockData{1}.rlnAmplitudeContrast(1:n_image); end;
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
for nimage=0:n_image-1;
disp(rlnImageName_{1+nimage});
end;%for nimage=0:n_image-1;
end;%if (verbose>1);
%%%%%%%%;
rlnImageFrame_ = zeros(n_image,1);
rlnImageLname_ = cell(n_image,1);
for nimage=0:n_image-1;
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
end;%for nimage=0:n_image-1;
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
rlnImageName_out_ = cell(n_image,1);
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
if (isempty(M_x_c___)); n_x_0 = size(tmp_M_x_c___,1); n_x_1 = size(tmp_M_x_c___,2); M_x_c___ = zeros(n_x_0,n_x_1,n_image); end;
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
for nimage=0:n_image-1;
if (verbose>2); disp(sprintf(' %% nimage %d/%d: rlnImageName_out %s rlnImageName %s',nimage,n_image,rlnImageName_out_{1+nimage},rlnImageName_{1+nimage})); end;
tmp_str_out_ = rlnImageName_out_{1+nimage};
tmp_str_0in_ = rlnImageName_{1+nimage};
tmp_frame_out = str2num(tmp_str_out_(1:strfind(tmp_str_out_,'@')-1));
tmp_frame_0in = str2num(tmp_str_0in_(1:strfind(tmp_str_0in_,'@')-1));
tmp_name_out = tmp_str_out_(strfind(tmp_str_out_,'@')+1:end);
tmp_name_0in = tmp_str_0in_(strfind(tmp_str_0in_,'@')+1:end);
assert(tmp_frame_out==tmp_frame_0in);
assert(strcmp(tmp_name_out,tmp_name_0in));
end;%for nimage=0:n_image-1;
%%%%%%%%;
end;%if ( exist(fname_star,'file'));
if (verbose); disp(sprintf(' %% [finished rlnImageName_from_star_0]')); end;
