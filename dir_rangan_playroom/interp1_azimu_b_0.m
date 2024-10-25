function ...
[ ...
 parameter ...
,v_out_bv__ ...
] = ...
interp1_azimu_b_0( ...
 parameter ...
,n_b_inp ...
,b_inp_b_ ...
,v_inp_bv__ ...
,n_b_out ...
,b_out_b_ ...
);

str_thisfunction = 'interp1_azimu_b_0';

%%%%%%%%;
if nargin<1;
disp(sprintf(' %% testing %s',str_thisfunction));
%%%%%%%%;
n_b_inp = 64;
b_uni_b_ = linspace(0,2*pi,n_b_inp+1); b_uni_b_ = transpose(b_uni_b_(1:n_b_inp));
db = mean(diff(b_uni_b_));
b_def_b_ = b_uni_b_ + db*sin(2*b_uni_b_);
v = @(b) cos(2*b) - sin(3*b) + cos(4*b) - sin(5*b);
v_uni_b_ = v(b_uni_b_);
v_def_b_ = v(b_def_b_);
[~,v_uni_from_def_b_] = interp1_azimu_b_0([],n_b_inp,b_def_b_,v_def_b_,n_b_inp,b_uni_b_);
[~,v_def_from_uni_b_] = interp1_azimu_b_0([],n_b_inp,b_uni_b_,v_uni_b_,n_b_inp,b_def_b_);
disp(sprintf(' %% v_uni_b_ vs v_uni_from_def_b_: %0.16f',fnorm(v_uni_b_-v_uni_from_def_b_)/max(1e-12,fnorm(v_uni_b_))));
disp(sprintf(' %% v_def_b_ vs v_def_from_uni_b_: %0.16f',fnorm(v_def_b_-v_def_from_uni_b_)/max(1e-12,fnorm(v_def_b_))));
[~,v_uni_from_def_from_uni_b_] = interp1_azimu_b_0([],n_b_inp,b_def_b_,v_def_from_uni_b_,n_b_inp,b_uni_b_);
[~,v_def_from_uni_from_def_b_] = interp1_azimu_b_0([],n_b_inp,b_uni_b_,v_uni_from_def_b_,n_b_inp,b_def_b_);
disp(sprintf(' %% v_uni_b_ vs v_uni_from_def_from_uni_b_: %0.16f',fnorm(v_uni_b_-v_uni_from_def_from_uni_b_)/max(1e-12,fnorm(v_uni_b_))));
disp(sprintf(' %% v_def_b_ vs v_def_from_uni_from_def_b_: %0.16f',fnorm(v_def_b_-v_def_from_uni_from_def_b_)/max(1e-12,fnorm(v_def_b_))));
figure(1);clf;figmed;
subplot(1,3,1);
plot(b_uni_b_,b_uni_b_,'k-',b_uni_b_,b_def_b_,'.');
xlim([0,2*pi]); ylim([0,2*pi]); axisnotick;
subplot(1,3,[2,3]);
hold on;
plot(b_uni_b_,v_uni_b_,'ro',b_uni_b_,v_uni_from_def_b_,'rx',b_uni_b_,v_uni_from_def_from_uni_b_,'rs');
plot(b_def_b_,v_def_b_,'go',b_def_b_,v_def_from_uni_b_,'gx',b_def_b_,v_def_from_uni_from_def_b_,'gs');
xlim([0,2*pi]); axisnotick;
%%%%%%%%;
n_b_inp = 32;
b_inp_b_ = 2*pi*rand(n_b_inp,1);
v_inp_bv__ = rand(n_b_inp,1,2);
n_b_out = 1024*8;
b_out_b_ = [0;2*pi*sort(rand(n_b_out-2,1),'ascend');2*pi];
parameter = struct('type','parameter');
parameter.str_method = 'spline';
[ ...
 ~ ...
,v_out_bv__ ...
] = ...
interp1_azimu_b_0( ...
 parameter ...
,n_b_inp ...
,b_inp_b_ ...
,v_inp_bv__ ...
,n_b_out ...
,b_out_b_ ...
);
v_inp_bv__ = reshape(v_inp_bv__,n_b_inp,[]);
v_out_bv__ = reshape(v_out_bv__,n_b_out,[]);
figure(2);clf;figbig;
markersize_sml = 8;
markersize_big = 32;
c_use__ = colormap('lines'); n_c_use = size(c_use__,1);
hold on;
for nl=0:size(v_inp_bv__,2)-1;
nc=max(0,min(n_c_use-1,nl));
plot(b_inp_b_,v_inp_bv__(:,1+nl),'.','MarkerSize',markersize_big,'Color',c_use__(1+nc,:));
plot(b_out_b_,v_out_bv__(:,1+nl),'.','MarkerSize',markersize_sml,'Color',c_use__(1+nc,:));
end;%for nl=0:size(v_inp_bv__,2)-1;
hold off;
grid on;
xlim([0,2*pi]); xlabel('b');
ylim_ = prctile(v_inp_bv__,[0,100],'all');
ylim_ = mean(ylim_) + 1.25*0.5*diff(ylim_)*[-1,+1];
ylim(ylim_);
ylabel('value');
disp('returning'); return;
%%%%%%%%;
end;%if nargin<1;
%%%%%%%%;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_b_inp=[]; end; na=na+1;
if (nargin<1+na); b_inp_b_=[]; end; na=na+1;
if (nargin<1+na); v_inp_bv_=[]; end; na=na+1;
if (nargin<1+na); n_b_out=[]; end; na=na+1;
if (nargin<1+na); b_out_b_=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
flag_verbose=parameter.flag_verbose;
if ~isfield(parameter,'diameter'); parameter.diameter=2*pi; end;
diameter=parameter.diameter;
if ~isfield(parameter,'str_method'); parameter.str_method='spline'; end;
str_method=parameter.str_method;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

flag_flatten=0;
if ndims(v_inp_bv__)<=2 & size(v_inp_bv__,2)==1; v_inp_bv__ = reshape(v_inp_bv__,n_b_inp,[]); flag_flatten = 1; end;

[~,tmp_ij_] = sort(b_inp_b_,'ascend');
b_inp_b_ = b_inp_b_(tmp_ij_);
v_inp_bv__ = v_inp_bv__(tmp_ij_,:);

n_b_3inp3 = 3 + n_b_inp + 3;
ij0 = max(1,min(n_b_inp,n_b_inp-2));
ij1 = max(1,min(n_b_inp,n_b_inp-1));
ij2 = max(1,min(n_b_inp,n_b_inp-0));
ij3 = max(1,min(n_b_inp,1));
ij4 = max(1,min(n_b_inp,2));
ij5 = max(1,min(n_b_inp,3));

size_v_ = size(v_inp_bv__);
b_3inp3_b_ = [ ...
  b_inp_b_(ij0)-diameter ...
; b_inp_b_(ij1)-diameter ...
; b_inp_b_(ij2)-diameter ...
; b_inp_b_ ...
; b_inp_b_(ij3)+diameter ...
; b_inp_b_(ij4)+diameter ...
; b_inp_b_(ij5)+diameter ...
];
v_3inp3_bv__ = [ ...
  v_inp_bv__(ij0,:) ...
; v_inp_bv__(ij1,:) ...
; v_inp_bv__(ij2,:) ...
; reshape(v_inp_bv__,n_b_inp,[]) ...
; v_inp_bv__(ij3,:) ...
; v_inp_bv__(ij4,:) ...
; v_inp_bv__(ij5,:) ...
];
v_out_bv__ = interp1(b_3inp3_b_,v_3inp3_bv__,b_out_b_,str_method);
reshape(v_out_bv__,[n_b_out,size_v_(2:end)]);
if flag_flatten; v_out_bv__ = v_out_bv__(:); end;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

