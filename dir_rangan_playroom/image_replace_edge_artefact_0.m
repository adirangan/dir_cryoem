function [M_out_,n_M_ext] = image_replace_edge_artefact_0(M_0in_,n_pixel,edge_tolerance,n_edge_overshoot,rseed);
%%%%%%%%;
% Detects and replaces duplicated edge pixels. ;
% Searches within a band of n_pixel (default 4). ;
% Any edges with a correlation greater than (1-edge_tolerance) will be replaced. ;
% Replacement at random from interior. ;
% If desired, edges will be 'overshot'. ;
%%%%%%%%;

if (nargin<2); n_pixel = []; end;
if (nargin<3); edge_tolerance = []; end;
if (nargin<3); n_edge_overshoot = []; end;
if (nargin<5); rseed = []; end;
if isempty(n_pixel); n_pixel = 4; end;
if isempty(edge_tolerance); edge_tolerance = 1e-3; end;
if isempty(n_edge_overshoot); n_edge_overshoot = 0; end;
if isempty(rseed); rseed = 0; end;

verbose=0;
if (verbose); disp(sprintf(' %% [entering image_replace_edge_artefact_0]')); end;

n_edge_0west = image_measure_edge_artefact(M_0in_,n_pixel,edge_tolerance);
M_0in_ = rot90(M_0in_);
n_edge_north = image_measure_edge_artefact(M_0in_,n_pixel,edge_tolerance);
M_0in_ = rot90(M_0in_);
n_edge_0east = image_measure_edge_artefact(M_0in_,n_pixel,edge_tolerance);
M_0in_ = rot90(M_0in_);
n_edge_south = image_measure_edge_artefact(M_0in_,n_pixel,edge_tolerance);
M_0in_ = rot90(M_0in_);
if (verbose); disp(sprintf(' %% found n_edge_0west: %d',n_edge_0west)); end;
if (verbose); disp(sprintf(' %% found n_edge_north: %d',n_edge_north)); end;
if (verbose); disp(sprintf(' %% found n_edge_0east: %d',n_edge_0east)); end;
if (verbose); disp(sprintf(' %% found n_edge_south: %d',n_edge_south)); end;
if (n_edge_0west> 0); n_edge_0west = n_edge_0west + n_edge_overshoot; end;
if (n_edge_north> 0); n_edge_north = n_edge_north + n_edge_overshoot; end;
if (n_edge_0east> 0); n_edge_0east = n_edge_0east + n_edge_overshoot; end;
if (n_edge_south> 0); n_edge_south = n_edge_south + n_edge_overshoot; end;
M_int_ = zeros(size(M_0in_));
M_int_(n_edge_north+1:end-n_edge_south,n_edge_0west+1:end-n_edge_0east) = 1;
M_int_ij_ = find( M_int_);
M_ext_ij_ = find(~M_int_);
n_M_int = numel(M_int_ij_);
n_M_ext = numel(M_ext_ij_);
M_out_ = zeros(size(M_0in_));
if (n_M_int> 0); 
M_out_ = M_0in_;
M_tmp_ = M_0in_(M_int_ij_);
M_out_(M_int_ij_) = (M_tmp_ - mean(M_tmp_,'all'))/std(M_tmp_,1,'all');
rng(rseed);
tmp_p_ = max(1,min(n_M_int,1+floor(n_M_int*rand(n_M_ext,1))));
M_out_(M_ext_ij_) = M_out_(M_int_ij_(tmp_p_));
end;%if (n_M_int> 0);

if (verbose); disp(sprintf(' %% [finished image_replace_edge_artefact_0]')); end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

function n_edge_artefact = image_measure_edge_artefact(M_0in_,n_pixel,edge_tolerance);
n_pixel = min(n_pixel,size(M_0in_,2));
edge_avg_ = mean(M_0in_(:,1:n_pixel),2);
edge_corr_ = corr(edge_avg_,M_0in_);
tmp_ij_ = find(cumprod(edge_corr_>=1-edge_tolerance));
n_edge_artefact = 0; if (numel(tmp_ij_)> 0); n_edge_artefact = max(tmp_ij_); end;

