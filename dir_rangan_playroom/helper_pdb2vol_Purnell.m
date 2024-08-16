function [vol,data] = helper_pdb2vol_Purnell(pdb,pix,trim)

%{
Copyright (c) 2022, Carson Purnell
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution

* Neither the name of  nor the names of its
  contributors may be used to endorse or promote products derived from this
  software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 %}
  
if nargin<2, error('requires both pdb and pixel size inputs'), end
if nargin<3, trim=0; end %don't trim volumes by default, exists for trimming bundles
%need user input or name/val to set what each model is for class/ID

%pdb to atoms
[~,~,ext] = fileparts(pdb);
if strcmp(ext,'.mat')
    q = load(pdb);
end

fid = fopen(pdb); 
text = textscan(fid,'%s','delimiter','\n','CommentStyle',{'REMARK'}); %import each line individually
text = text{1}; %fix being inside a 1x1 cell array for no reason

%delete terminator and temperature/ANISOU records that mess with model reading and parsing
ix = strncmp(text,'TER',3); text(ix) = [];
ix = strncmp(text,'ANISOU',6); text(ix) = [];
ix = strncmp(text,'HETATM',6); text(ix) = []; %delete heteroatoms for sanity

modstart = find(strncmp(text,'MODEL ',6)); %find start of each model entry
modend = find(strncmp(text,'ENDMDL',6)); %find the end of each model entry

if isempty(modstart) %if single-model, extract start and end of atom lines
    modstart = find(strncmp(text(1:round(end/2)),'ATOM  ',6)); modstart = modstart(1);
    endatm = find(strncmp(text(modstart:end),'ATOM  ',6)); endatm = endatm(end);
    endhet = find(strncmp(text(modstart:end),'HETATMjj',6)); 
    if ~isempty(endhet), endhet = endhet(end); else endhet = 0; end %#ok<SEPEX>
    modend = max(endatm,endhet)+modstart-1; %adjust for having searched only part of the list for speed
    model{1} = text(modstart:modend); models = 1;
elseif numel(modstart)==numel(modend) %if counts match, populate model counts
    models = numel(modstart); model = cell(models,1);
    for i=1:models %extract lines for each individual model
        model{i} = text(modstart(i)+1:modend(i)-1);
    end
elseif numel(modstart)~=numel(modend) %check if model numbers are the same
    error('failed to determine model numbers')
end

data = cell(numel(model),2);
for i=1:models %loop through detected models
    atom = cell(1,numel(model{i}));
    coords = zeros(3,numel(model{i}));
    for j=1:numel(model{i}) %loop through each line of the model
        line = model{i}{j}; %extract single line for reading
        atom{j} = upper(strrep(line(77:78),' ','')); %read atom identifier string, prune spaces
        if strcmp(atom{j},''), error('Bad PDB, try resaving with proper software'), end
        coords(1:3,j) = sscanf(line(31:54),'%f',[3 1]); %reads coords, sscanf>str2double>>str2num
    end
    data{i,1} = atom; data{i,2} = coords;
end


%initialize atomic magnitude information
mag = struct('H',0,'C',6+1.3,'N',7+1.1,'O',8+0.2,'P',15,'S',16+0.6);
%kind of messy version required for heteroatom records
%mag = struct('H',0,'C',6,'N',7,'O',8,'P',15,'S',16,...
%    'MG',12,'ZN',30,'MN',25,'F',9,'CL',17,'CA',20);

%hydrogen intensity instead shifted onto average H per organic atom, because PDB inconsistently use H records
%currently using atomic number, should use a more accurate scattering factor measure
%have seen h=25, c=130, n=108, o=97, s=100, p=267 for va^3 scattering potentials - need inorganic atoms
% hydrogens per atom of c=1.3, n=1.1, o=0.2, s=0.6 to bypass needing to add hydrogens manually

%faster, vectorized adjustments and limits to coordinates and bounding box
[a,b] = bounds(horzcat(data{:,2}),2); %bounds of all x/y/z in row order
adj = max(a*-1,0)+pix; %coordinate adjustment to avoid indexing below 1
lim = round( (adj+b)/pix +1); %array size to place into, same box for all models

emvol = cell(models,1);

for i=1:models
    atomid = data{i,1}; %single column, hopefully for speed
    points = data{i,2}; %column per atom, hopefully faster indexing
    em = zeros(lim'); %initialize empty array
    for j=1:numel(atomid)
        atom = atomid{j};
        %dt = [data{j,2:4}];
        coord = points(:,j);
        
        co = round((coord+adj)./pix); %vectorize coordinate math
        x=co(1); y=co(2); z=co(3); %still can't split vector into singles
        
        density = mag.(atom); %retrieve atom info, struct is faster than using a container.map
        em(x,y,z) = em(x,y,z)+density;
    end
    if trim==1 %for bundles, should probably do singles here or helper_input for efficiency
        em = em(:,any(em ~= 0,[1 3]),:); 
        em = em(any(em ~= 0,[2 3]),:,:); 
        em = em(:,:,any(em ~= 0,[1 2]));
    end
    
    emvol{i} = em;
end

emvol = reshape(emvol,1,numel(emvol)); %make horizontal because specifying it initially doesn't work
vol = emvol;

end
