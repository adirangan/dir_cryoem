function ...
[ ...
 parameter ...
,n_Atom ...
,Atom_centered_X_ ...
,Atom_centered_Y_ ...
,Atom_centered_Z_ ...
,Atom_occupancy_ ...
] = ...
a_x_u___from_pdb_0( ...
 parameter ...
,fname_pdb ...
);

str_thisfunction = 'a_x_u___from_pdb_0';

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); fname_pdb=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
flag_verbose=parameter.flag_verbose;
if ~isfield(parameter,'flag_disp'); parameter.flag_disp=0; end;
flag_disp=parameter.flag_disp; nf=0;
if ~isfield(parameter,'flag_center'); parameter.flag_center=1; end;
flag_center=parameter.flag_center;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

%%%%;
pdb = pdbread(fname_pdb);
Atom_ = pdb.Model.Atom;
n_Atom = numel(Atom_);
Atom_X_ = zeros(n_Atom,1);
Atom_Y_ = zeros(n_Atom,1);
Atom_Z_ = zeros(n_Atom,1);
Atom_occupancy_ = zeros(n_Atom,1);
for nAtom=0:n_Atom-1;
Atom_X_(1+nAtom) = Atom_(1+nAtom).X;
Atom_Y_(1+nAtom) = Atom_(1+nAtom).Y;
Atom_Z_(1+nAtom) = Atom_(1+nAtom).Z;
Atom_occupancy_(1+nAtom) = Atom_(1+nAtom).occupancy;
end;%for nAtom=0:n_Atom-1;
%%%%;
Atom_centered_X_ = Atom_X_;
Atom_centered_Y_ = Atom_Y_;
Atom_centered_Z_ = Atom_Z_;
if flag_center;
Atom_avg_X_ = mean(Atom_centered_X_.*Atom_occupancy_);
Atom_avg_Y_ = mean(Atom_centered_Y_.*Atom_occupancy_);
Atom_avg_Z_ = mean(Atom_centered_Z_.*Atom_occupancy_);
Atom_centered_X_ = Atom_centered_X_ - Atom_avg_X_;
Atom_centered_Y_ = Atom_centered_Y_ - Atom_avg_Y_;
Atom_centered_Z_ = Atom_centered_Z_ - Atom_avg_Z_;
end;%if flag_center;
%%%%;

%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;figsml;
markersize_use = 16;
subplot(1,1,1);
plot3(Atom_centered_X_,Atom_centered_Y_,Atom_centered_Z_,'.','MarkerSize',markersize_use,'MarkerFaceColor','c');
xlabel('x'); ylabel('y'); zlabel('z');
axisnotick3d; axis vis3d; grid on;
end;%if flag_disp;
%%%%%%%%;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
