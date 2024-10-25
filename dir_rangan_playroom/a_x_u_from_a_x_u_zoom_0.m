function ...
[ ...
b_x_u_ ...
] = ...
a_x_u_from_a_x_u_zoom_0( ...
 parameter ...
,a_x_u_ ...
);

str_thisfunction = 'a_x_u_from_a_x_u_zoom_0';

if nargin<1;
disp(sprintf(' %% testing %s',str_thisfunction));

disp('returning'); return;
end;%if nargin<1;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); a_x_u_=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
flag_verbose=parameter.flag_verbose;
if ~isfield(parameter,'val_zoom'); parameter.val_zoom=sqrt(2); end;
val_zoom=parameter.val_zoom;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;
n_xxx = numel(a_x_u_);
n_x = round(n_xxx^(1/3));
a_x_u_ = reshape(a_x_u_,[n_x,n_x,n_x]);
diameter_x_c = 2.0;
half_diameter_x_c = diameter_x_c/2.0;
x_u_0_ = transpose(linspace(-half_diameter_x_c,+half_diameter_x_c,n_x));
x_u_1_ = transpose(linspace(-half_diameter_x_c,+half_diameter_x_c,n_x));
x_u_2_ = transpose(linspace(-half_diameter_x_c,+half_diameter_x_c,n_x));
[x_u_0___,x_u_1___,x_u_2___] = ndgrid(x_u_0_,x_u_1_,x_u_2_);
y_u_0_ = transpose(linspace(-half_diameter_x_c/max(1e-12,val_zoom),+half_diameter_x_c/max(1e-12,val_zoom),n_x));
y_u_1_ = transpose(linspace(-half_diameter_x_c/max(1e-12,val_zoom),+half_diameter_x_c/max(1e-12,val_zoom),n_x));
y_u_2_ = transpose(linspace(-half_diameter_x_c/max(1e-12,val_zoom),+half_diameter_x_c/max(1e-12,val_zoom),n_x));
[y_u_0___,y_u_1___,y_u_2___] = ndgrid(y_u_0_,y_u_1_,y_u_2_);
%b_x_u_ = interp3(x_u_0___,x_u_1___,x_u_2___,a_x_u_,y_u_0___,y_u_1___,y_u_2___,'cubic');
b_x_u_ = interp3( ...
		   reshape(x_u_0_,[n_x,1]) ...
		  ,reshape(x_u_1_,[n_x,1]) ...
		  ,reshape(x_u_2_,[n_x,1]) ...
		  ,a_x_u_ ...
		  ,reshape(y_u_0___,[n_xxx,1]) ...
		  ,reshape(y_u_1___,[n_xxx,1]) ...
		  ,reshape(y_u_2___,[n_xxx,1]) ...
		  ,'cubic' ...
		  );
b_x_u_ = reshape(b_x_u_,[n_x,n_x,n_x]);
if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
