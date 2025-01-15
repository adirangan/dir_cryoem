function [output_,xl_,yl_] = hist2d_0(x_,y_,n_x,n_y,xl_,yl_,weight_);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% function output_ = hist2d_0(x_,y_,n_x,n_y,xl_,yl_,weight_);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% 
% Sets up a 2d histogram;
%
% Inputs: 
% x_ double array of length N storing x-coordinates of data. ;
% y_ double array of length N storing y-coordinates of data. ;
% xl_ double array of length 2 storing min and max values for x-coordinates. ;
% yl_ double array of length 2 storing min and max values for y-coordinates. ;
% n_x integer number of bins in x-direction. ;
% n_y integer number of bins in y-direction. ;
% weight_ double array of length N storing weight for each data-point. ;
%
% Outputs:
% output_ matrix of size n_y-x-n_x storing 2d-histogram of data. ;
% Note that output_ is in sparse format. ;
%
% Test with: hist2d_0();
% 

if (nargin<1);
G1_ = 0.15*randn(1024*4,2) + repmat([-0.15,-0.25],1024*4,1);
G2_ = 0.15*randn(1024*1,2) + repmat([+0.15,+0.25],1024*1,1);
G_ = [G1_;G2_];
colormap('hot');
imagesc(hist2d_0(G_(:,1),G_(:,2),32,33,[-1,+1],[-1,+1])); colorbar;
disp('returning');return;
end;%if (nargin<1);

na=0;
if (nargin<1+na); x_=[]; end; na=na+1;
if (nargin<1+na); y_=[]; end; na=na+1;
if (nargin<1+na); n_x=[]; end; na=na+1;
if (nargin<1+na); n_y=[]; end; na=na+1;
if (nargin<1+na); xl_=[]; end; na=na+1;
if (nargin<1+na); yl_=[]; end; na=na+1;
if (nargin<1+na); weight_=[]; end; na=na+1;

if isempty(n_x); n_x = 32; end;
if isempty(n_y); n_y = 32; end;

if isempty(xl_); xl_=[min(x_),max(x_)]; xl_ = mean(xl_) + diff(xl_)/2*1.0625*[-1,1]; end;
if isempty(yl_); yl_=[min(y_),max(y_)]; yl_ = mean(yl_) + diff(yl_)/2*1.0625*[-1,1]; end;

x_ = x_(:);
y_ = y_(:);
n_data = numel(x_);
if isempty(weight_); weight_ = ones(n_data,1); end;
weight_ = weight_(:);
bx_ = min(n_x-1,max(0,floor(n_x*(x_-min(xl_))/(max(xl_)-min(xl_)))));
by_ = min(n_y-1,max(0,floor(n_y*(y_-min(yl_))/(max(yl_)-min(yl_)))));
output_ = sparse(1+by_,1+bx_,weight_,n_y,n_x);
