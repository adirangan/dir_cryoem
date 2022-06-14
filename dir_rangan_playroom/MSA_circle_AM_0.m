function ...
[ ...
 parameter ...
,x_pos_ ...
,n_iteration ...
,xyr_id__ ...
] = ...
MSA_circle_AM_0( ...
 parameter ...
,M_pd__ ...
,x_init_ ...
);		 

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); M_pd__=[]; end; na=na+1;
if (nargin<1+na); x_init_=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;
if (~isfield(parameter,'tolerance_master')); parameter.tolerance_master = 1e-2; end;
if (~isfield(parameter,'n_iteration')); parameter.n_iteration = 32; end;
if (~isfield(parameter,'flag_plot')); parameter.flag_plot = 1; end;
tolerance_master = parameter.tolerance_master;
n_iteration = parameter.n_iteration;
flag_plot = parameter.flag_plot;

if isempty(M_pd__); M_pd__ = zeros(1,2); end;
if isempty(x_init_); x_init_ = zeros(1,2); end;

n_dim = size(M_pd__,2);
assert(n_dim==2);
n_point = size(M_pd__,1);
assert(n_dim==numel(x_init_));

x_cur_ = reshape(x_init_,[1,n_dim]);
N_pd__ = bsxfun(@minus,M_pd__,x_cur_);
r_cur_p_ = sqrt(sum(N_pd__.^2,2));
w_cur_p_ = atan2(N_pd__(:,1+1),N_pd__(:,1+0));
r_cur = mean(r_cur_p_,1);
xyr_id__ = zeros(n_iteration,1+n_dim);

flag_continue=1;
niteration=0; 
while flag_continue;
c_cur_p_ = cos(w_cur_p_);
s_cur_p_ = sin(w_cur_p_);
c_avg = mean(c_cur_p_);
s_avg = mean(s_cur_p_);
cM_avg = mean(M_pd__(:,1+0));
sM_avg = mean(M_pd__(:,1+1));
rM_avg = mean( M_pd__(:,1+0).*c_cur_p_ + M_pd__(:,1+1).*s_cur_p_ );
tmp_A__ = [ 1 0 c_avg ; 0 1 s_avg ; c_avg s_avg 1];
tmp_b_ = [cM_avg ; sM_avg ; rM_avg];
x_pos_ = tmp_A__\tmp_b_;
xyr_id__(1+niteration,:) = x_pos_;
error_l2 = fnorm(x_pos_(1:2) - x_cur_(1:2));
flag_continue = (error_l2>=tolerance_master & niteration<n_iteration);
if (flag_continue);
x_cur_ = reshape(x_pos_(1:2),[1,n_dim]);
N_pd__ = bsxfun(@minus,M_pd__,x_cur_);
w_cur_p_ = atan2(N_pd__(:,2),N_pd__(:,1));
niteration = niteration+1;
end;%if (flag_continue);
end;%while flag_continue;

n_iteration = 1+niteration;
xyr_id__ = xyr_id__(1:n_iteration,:);

if flag_plot;
figure(1);clf;figsml;
markersize_use = 8;
hold on;
plot(M_pd__(:,1+0),M_pd__(:,1+1),'o','MarkerFaceColor',0.85*[1,1,1]);
viscircles(xyr_id__(:,1:2),xyr_id__(:,3),'Color','r','LineWidth',1);
end;%if flag_plot;














