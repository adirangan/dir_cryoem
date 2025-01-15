if ~exist('a_k_p_qk_','var'); a_k_p_qk_=[]; end;
if ~exist('a_k_p_qk__','var'); a_k_p_qk__=[]; end;
if ~exist('dvol_a_k_p_qk_','var'); dvol_a_k_p_qk_=[]; end;
if ~exist('dvol_a_k_p_qk__','var'); dvol_a_k_p_qk__=[]; end;
if ~exist('a_R_k_p_qk_','var'); a_R_k_p_qk_=[]; end;
if ~exist('a_R_k_p_qk__','var'); a_R_k_p_qk__=[]; end;
if ~exist('dvol_a_R_k_p_qk_','var'); dvol_a_R_k_p_qk_=[]; end;
if ~exist('dvol_a_R_k_p_qk__','var'); dvol_a_R_k_p_qk__=[]; end;

if ~isempty(a_k_p_qk__);
a_k_p_qk__ = reshape(a_k_p_qk__,[n_q,n_k_p_r]); a_k_p_qk_ = reshape(a_k_p_qk__,[n_qk,1]);
end;%if ~isempty(a_k_p_qk__);

if ~isempty(a_R_k_p_qk__);
a_R_k_p_qk__ = reshape(a_R_k_p_qk__,[n_q,n_k_p_r]); a_R_k_p_qk_ = reshape(a_R_k_p_qk__,[n_qk,1]);
end;%if ~isempty(a_R_k_p_qk__);

if ~isempty(dvol_a_k_p_qk__);
dvol_a_k_p_qk__ = reshape(dvol_a_k_p_qk__,[n_q,n_k_p_r]); dvol_a_k_p_qk_ = reshape(dvol_a_k_p_qk__,[n_qk,1]);
end;%if ~isempty(dvol_a_k_p_qk__);

if ~isempty(dvol_a_R_k_p_qk__);
dvol_a_R_k_p_qk__ = reshape(dvol_a_R_k_p_qk__,[n_q,n_k_p_r]); dvol_a_R_k_p_qk_ = reshape(dvol_a_R_k_p_qk__,[n_qk,1]);
end;%if ~isempty(dvol_a_R_k_p_qk__);
