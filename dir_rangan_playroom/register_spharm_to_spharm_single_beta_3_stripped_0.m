function ...
[ ...
 X1_mmb___ ...
 t_sumn ...
 t_abab ...
 t_sumz ...
 t_fft2 ...
 n_mult ...
] = ...
register_spharm_to_spharm_single_beta_3_stripped_0( ...
 flag_verbose ...
,n_m_max ...
,n_UX_rank ...
,n_l_max ...
,weight_n_ ...
,a_mnl___ ...
,b_mnl___ ...
,n_UZ_rank ...
,U_lz__ ...
,n_beta ...
,d_mmzb____ ...
);
str_thisfunction = 'register_spharm_to_spharm_single_beta_3_stripped_0';
if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;
%%%%%%%%;
% merge a_ and b_. ;
%%%%%%%%;
ab_mml___ = [];
t_0in = tic;
ab_mml___ = zeros(n_m_max,n_m_max,n_l_max);
if  isempty(weight_n_);
for nl_max=0:n_l_max-1;
ab_mml___(:,:,1+nl_max) = conj(squeeze(a_mnl___(:,1:n_UX_rank,1+nl_max)))*transpose(b_mnl___(:,1:n_UX_rank,1+nl_max)); %<-- should probably replace with pagemtimes. However, third-step is 4x slower, so no rush. ;
end;%for nl_max=0:n_l_max-1;
end;%if  isempty(weight_n_);
if ~isempty(weight_n_);
for nl_max=0:n_l_max-1;
ab_mml___(:,:,1+nl_max) = conj(squeeze(a_mnl___(:,1:n_UX_rank,1+nl_max)))*diag(weight_n_)*transpose(b_mnl___(:,1:n_UX_rank,1+nl_max)); %<-- should probably replace with weighted pagemtimes. However, third-step is 4x slower, so no rush. ;
end;%for nl_max=0:n_l_max-1;
end;%if ~isempty(weight_n_);
t_out = toc(t_0in); if (flag_verbose>0); disp(sprintf(' %% %s sum over k: t %0.6fs',str_thisfunction,t_out)); end;
t_sumn = t_out;
%%%%%%%%;
% compress ab_. ;
%%%%%%%%;
if  isempty(U_lz__);
ab_mmz___ = ab_mml___;
n_UZ_rank = n_l_max;
t_out = 0;
end;%if  isempty(U_lz__);
if ~isempty(U_lz__);
t_0in = tic;
ab_mmz___ = reshape(reshape(ab_mml___,[n_m_max*n_m_max,n_l_max])*U_lz__(:,1:n_UZ_rank),[n_m_max,n_m_max,n_UZ_rank]);
t_out = toc(t_0in);
end;%if ~isempty(U_lz__);
if (flag_verbose>0); disp(sprintf(' %% %s compress ab: t %0.6fs',str_thisfunction,t_out)); end;
t_abab = t_out;
%%%%%%%%;
% merge ab_ and d_. ;
%%%%%%%%;
t_0in = tic;
X0_mmb___ = zeros(n_m_max,n_m_max,n_beta);
%{
for nbeta=0:n_beta-1;
X0_mmb___(:,:,1+nbeta) = sum(ab_mmz___(:,:,1:n_UZ_rank).*d_mmzb____(:,:,1:n_UZ_rank,1+nbeta),3);
end;%for nbeta=0:n_beta-1;
%}
X0_mmb___ = squeeze(sum(bsxfun(@times,ab_mmz___(:,:,1:n_UZ_rank),d_mmzb____),1+2)); %<-- This is about 4x slower than first-step. ;
% X0_mmb___ = squeeze(sum(permute(bsxfun(@times,ab_mmz___,d_mmzb____),1+[2,0,1,3]),1+0)); %<-- slower than above. ;
t_out = toc(t_0in); if (flag_verbose>0); disp(sprintf(' %% %s sum over l: t %0.6fs',str_thisfunction,t_out)); end;
t_sumz = t_out;
%%%%%%%%;
% fft2. ;
%%%%%%%%;
t_0in = tic;
m_ij_ = [floor(n_m_max/2)+1:n_m_max , 1:floor(n_m_max/2)];
X0_mmb___ = X0_mmb___(m_ij_,m_ij_,:);
t_out = toc(t_0in); if (flag_verbose>0); disp(sprintf(' %% %s reorder X0_mmb___: t %0.6fs',str_thisfunction,t_out)); end;
t_0in = tic;
X1_mmb___ = zeros(n_m_max,n_m_max,n_beta);
%{
for nbeta=0:n_beta-1;
 %tmp_X0_mm__ = recenter2(squeeze(X0_mmb___(:,:,1+nbeta))); %<-- recentering no longer necessary if done above. ;
tmp_X0_mm__ = squeeze(X0_mmb___(:,:,1+nbeta));
X1_mmb___(:,:,1+nbeta) = fft2(tmp_X0_mm__);
end;%for nbeta=0:n_beta-1;
%}
X1_mmb___ = fft(fft(X0_mmb___,[],1),[],2); %<-- seems just as fast as X1_mmb___ = fft2(X0_mmb___);
t_out = toc(t_0in); if (flag_verbose>0); disp(sprintf(' %% %s fft2: t %0.6fs',str_thisfunction,t_out)); end;
t_fft2 = t_out;
%%%%%%%%;
n_mult = ...
  (n_m_max*n_m_max)*n_UX_rank*n_l_max ...
+ (n_m_max*n_m_max)*n_l_max*n_UZ_rank ...
+ (n_m_max*n_m_max)*n_UZ_rank*n_beta ...
+ (n_m_max*n_m_max)*log(n_m_max)*n_beta ...
  ;
if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
