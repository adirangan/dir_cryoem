function ...
[ ...
 parameter ...
,T_23__ ...
,heron_rc_ ...
,x0_patch_rc5__ ...
,x1_patch_rc5__ ...
,c_patch_1rc3___ ...
] = ...
patch_simplex_3_0( ...
 parameter ...
,rho0_rc__ ...
,rho1_rc__ ...
,val_rc__ ...
,vlim_ ...
,c_use__ ...
);

str_thisfunction = 'patch_simplex_3_0';

if nargin<1;
disp(sprintf(' %% testing %s',str_thisfunction));

disp('returning'); return;
end;%if nargin<1;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); rho0_rc__=[]; end; na=na+1;
if (nargin<1+na); rho1_rc__=[]; end; na=na+1;
if (nargin<1+na); val_rc__=[]; end; na=na+1;
if (nargin<1+na); vlim_=[]; end; na=na+1;
if (nargin<1+na); c_use__=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
flag_verbose=parameter.flag_verbose;
if ~isfield(parameter,'flag_plot'); parameter.flag_plot=1; end;
flag_plot=parameter.flag_plot;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

if isempty(c_use__);
c_use__ = colormap_beach();
end;%if isempty(c_use__);
n_c_use = size(c_use__,1);

%%%%%%%%;
[n_r,n_c] = size(rho0_rc__);
n_rc = n_r*n_c;
rho2_rc__ = 1 - rho0_rc__ - rho1_rc__ ;
%%%%%%%%;
T_23__ = ...
[ ...
  1 , 0 ...
; 0 , 1 ...
] ...
* ...
[ ...
  -1/sqrt(2) , +1/sqrt(2) , +0         ...
; -1/sqrt(6) , -1/sqrt(6) , +2/sqrt(6) ...
];
%%%%%%%%;
x01_2rc__ = T_23__*transpose([rho0_rc__(:),rho1_rc__(:),rho2_rc__(:)]);
x0_cc_rc__ = reshape(x01_2rc__(1+0,:),[n_r,n_c]);
%%%%;
x0_pp_rc__ = zeros(1+n_r+1,1+n_c+1);
x0_pp_rc__(1+[1:n_r],1+0      ) = x0_cc_rc__(1:n_r  ,1+0    );
x0_pp_rc__(1+[1:n_r],1+1+n_c  ) = x0_cc_rc__(1:n_r  ,1+n_c-1);
x0_pp_rc__(1+0      ,1+[1:n_c]) = x0_cc_rc__(1+0    ,1:n_c  );
x0_pp_rc__(1+1+n_r  ,1+[1:n_c]) = x0_cc_rc__(1+n_r-1,1:n_c  );
x0_pp_rc__(1+0      ,1+0      ) = x0_cc_rc__(1+0    ,1+0    );
x0_pp_rc__(1+1+n_r  ,1+0      ) = x0_cc_rc__(1+n_r-1,1+0    );
x0_pp_rc__(1+0      ,1+1+n_c  ) = x0_cc_rc__(1+0    ,1+n_c-1);
x0_pp_rc__(1+1+n_r  ,1+1+n_c  ) = x0_cc_rc__(1+n_r-1,1+n_c-1);
x0_pp_rc__(1+[1:n_r],1+[1:n_c]) = x0_cc_rc__(:,:);
%%%%;
x1_cc_rc__ = reshape(x01_2rc__(1+1,:),[n_r,n_c]);
x1_pp_rc__ = zeros(1+n_r+1,1+n_c+1);
x1_pp_rc__(1+[1:n_r],1+0      ) = x1_cc_rc__(1:n_r  ,1+0    );
x1_pp_rc__(1+[1:n_r],1+1+n_c  ) = x1_cc_rc__(1:n_r  ,1+n_c-1);
x1_pp_rc__(1+0      ,1+[1:n_c]) = x1_cc_rc__(1+0    ,1:n_c  );
x1_pp_rc__(1+1+n_r  ,1+[1:n_c]) = x1_cc_rc__(1+n_r-1,1:n_c  );
x1_pp_rc__(1+0      ,1+0      ) = x1_cc_rc__(1+0    ,1+0    );
x1_pp_rc__(1+1+n_r  ,1+0      ) = x1_cc_rc__(1+n_r-1,1+0    );
x1_pp_rc__(1+0      ,1+1+n_c  ) = x1_cc_rc__(1+0    ,1+n_c-1);
x1_pp_rc__(1+1+n_r  ,1+1+n_c  ) = x1_cc_rc__(1+n_r-1,1+n_c-1);
x1_pp_rc__(1+[1:n_r],1+[1:n_c]) = x1_cc_rc__(:,:);
%%%%;
%{
x0_we_rc__ = 0.5*(x0_pp_rc__ + circshift(x0_pp_rc__,-1,1));
x0_no_rc__ = 0.5*(x0_pp_rc__ + circshift(x0_pp_rc__,+1,2));
x0_ea_rc__ = 0.5*(x0_pp_rc__ + circshift(x0_pp_rc__,+1,1));
x0_so_rc__ = 0.5*(x0_pp_rc__ + circshift(x0_pp_rc__,-1,2));
x1_we_rc__ = 0.5*(x1_pp_rc__ + circshift(x1_pp_rc__,-1,1));
x1_no_rc__ = 0.5*(x1_pp_rc__ + circshift(x1_pp_rc__,+1,2));
x1_ea_rc__ = 0.5*(x1_pp_rc__ + circshift(x1_pp_rc__,+1,1));
x1_so_rc__ = 0.5*(x1_pp_rc__ + circshift(x1_pp_rc__,-1,2));
%}
%%%%;
%{
x0_we_rc__ = x0_we_rc__(1+[1:n_r],1+[1:n_c]);
x0_no_rc__ = x0_no_rc__(1+[1:n_r],1+[1:n_c]);
x0_ea_rc__ = x0_ea_rc__(1+[1:n_r],1+[1:n_c]);
x0_so_rc__ = x0_so_rc__(1+[1:n_r],1+[1:n_c]);
x1_we_rc__ = x1_we_rc__(1+[1:n_r],1+[1:n_c]);
x1_no_rc__ = x1_no_rc__(1+[1:n_r],1+[1:n_c]);
x1_ea_rc__ = x1_ea_rc__(1+[1:n_r],1+[1:n_c]);
x1_so_rc__ = x1_so_rc__(1+[1:n_r],1+[1:n_c]);
%}
%%%%;
x0_nw_rc__ = 0.5*(x0_pp_rc__ + circshift(x0_pp_rc__,[-1,+1]));
x0_ne_rc__ = 0.5*(x0_pp_rc__ + circshift(x0_pp_rc__,[+1,+1]));
x0_se_rc__ = 0.5*(x0_pp_rc__ + circshift(x0_pp_rc__,[+1,-1]));
x0_sw_rc__ = 0.5*(x0_pp_rc__ + circshift(x0_pp_rc__,[-1,-1]));
x1_nw_rc__ = 0.5*(x1_pp_rc__ + circshift(x1_pp_rc__,[-1,+1]));
x1_ne_rc__ = 0.5*(x1_pp_rc__ + circshift(x1_pp_rc__,[+1,+1]));
x1_se_rc__ = 0.5*(x1_pp_rc__ + circshift(x1_pp_rc__,[+1,-1]));
x1_sw_rc__ = 0.5*(x1_pp_rc__ + circshift(x1_pp_rc__,[-1,-1]));
x0_nw_rc__ = x0_nw_rc__(1+[1:n_r],1+[1:n_c]);
x0_ne_rc__ = x0_ne_rc__(1+[1:n_r],1+[1:n_c]);
x0_se_rc__ = x0_se_rc__(1+[1:n_r],1+[1:n_c]);
x0_sw_rc__ = x0_sw_rc__(1+[1:n_r],1+[1:n_c]);
x1_nw_rc__ = x1_nw_rc__(1+[1:n_r],1+[1:n_c]);
x1_ne_rc__ = x1_ne_rc__(1+[1:n_r],1+[1:n_c]);
x1_se_rc__ = x1_se_rc__(1+[1:n_r],1+[1:n_c]);
x1_sw_rc__ = x1_sw_rc__(1+[1:n_r],1+[1:n_c]);
heron_rc_ = ...
+heron_0( permute(reshape(cat(3,x0_nw_rc__,x0_cc_rc__,x0_ne_rc__,x1_nw_rc__,x1_cc_rc__,x1_ne_rc__),[n_rc,3,2]),[3,2,1]) ) ...
+heron_0( permute(reshape(cat(3,x0_ne_rc__,x0_cc_rc__,x0_se_rc__,x1_ne_rc__,x1_cc_rc__,x1_se_rc__),[n_rc,3,2]),[3,2,1]) ) ...
+heron_0( permute(reshape(cat(3,x0_se_rc__,x0_cc_rc__,x0_sw_rc__,x1_se_rc__,x1_cc_rc__,x1_sw_rc__),[n_rc,3,2]),[3,2,1]) ) ...
+heron_0( permute(reshape(cat(3,x0_sw_rc__,x0_cc_rc__,x0_nw_rc__,x1_sw_rc__,x1_cc_rc__,x1_nw_rc__),[n_rc,3,2]),[3,2,1]) ) ...
;
%%%%;
if isempty(val_rc__); val_rc__ = reshape(heron_rc_,[n_r,n_c]); end;
if isempty(vlim_);
vlim_ = prctile(val_rc__,[5,95],'all');
end;%if isempty(vlim_);
%%%%;
x0_patch_rc5__ = [x0_nw_rc__(:),x0_ne_rc__(:),x0_se_rc__(:),x0_sw_rc__(:),x0_nw_rc__(:)];
x1_patch_rc5__ = [x1_nw_rc__(:),x1_ne_rc__(:),x1_se_rc__(:),x1_sw_rc__(:),x1_nw_rc__(:)];
nc_rc_ = max(0,min(n_c_use-1,round(n_c_use*(val_rc__(:) - min(vlim_))/max(1e-12,diff(vlim_)))));
c_patch_1rc3___ = reshape(c_use__(1+nc_rc_,:),[1,n_rc,3]);
%%%%;

%%%%%%%%;
if flag_plot;
%%%%%%%%;
p=patch(transpose(x0_patch_rc5__),transpose(x1_patch_rc5__),c_patch_1rc3___);
set(p,'EdgeColor','none');
%%%%%%%%;
end;%if flag_plot;
%%%%%%%%;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;


