function ...
[ ...
 parameter ...
,p ...
] = ...
imagesc_d_FTK( ...
 parameter ...
,delta_x_ ...
,delta_y_ ...
,data_d_ ...
,clim_ ...
,c_use__ ...
);

str_thisfunction = 'imagesc_d_FTK';

if nargin<1;
disp(sprintf(' %% testing %s',str_thisfunction));

disp('returning'); return;
end;%if nargin<1;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); delta_x_=[]; end; na=na+1;
if (nargin<1+na); delta_y_=[]; end; na=na+1;
if (nargin<1+na); data_d_=[]; end; na=na+1;
if (nargin<1+na); clim_=[]; end; na=na+1;
if (nargin<1+na); c_use__=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
flag_verbose=parameter.flag_verbose;

if isempty(clim_); clim_ = prctile(data_d_(:),[  0,100],'all'); end;
if isempty(c_use__); c_use__ = colormap_beach(); end;
n_c_use = size(c_use__,1);

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

n_d = numel(delta_x_); assert(numel(delta_y_)==n_d);
assert(numel(data_d_)==n_d);
dx2 = median(diff(sort(unique(delta_x_))))/2;
dy2 = median(diff(sort(unique(delta_y_))))/2;
nc_d_ = max(0,min(n_c_use-1,floor(n_c_use*(data_d_(:)-min(clim_))/diff(clim_))));
x_d_ = delta_x_;
x_d4__ = bsxfun(@plus,repmat(x_d_,[1,4]),dx2*[-1,-1,+1,+1]);
x_4d__ = permute(x_d4__,[2,1]);
y_d_ = delta_y_;
y_d4__ = bsxfun(@plus,repmat(y_d_,[1,4]),dy2*[-1,+1,+1,-1]);
y_4d__ = permute(y_d4__,[2,1]);
c_1d3___ = reshape(c_use__(1+nc_d_,:),[1,n_d,3]);
p = patch(x_4d__,y_4d__,c_1d3___);
set(p,'LineStyle','none','EdgeColor','none');

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

