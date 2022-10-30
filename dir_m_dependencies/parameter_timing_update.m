function parameter = parameter_timing_update(parameter,str_field,dt,ncall,nop);
if (nargin<3); dt=0; end;
if (nargin<4); ncall=1; end;
if (nargin<5); nop=0; end;

if isempty(parameter);
parameter = struct('type','parameter');
end;%if isempty(parameter);
%%%%%%%%;
if (~isfield(parameter,'timing'));
parameter.timing = cell(1,6);
parameter.timing{1,1} = 'label';
parameter.timing{1,2} = 'total time';
parameter.timing{1,3} = 'number of calls';
parameter.timing{1,4} = 'number of ops';
parameter.timing{1,5} = 'MHz';
parameter.timing{1,6} = 'GHz';
end; %<-- parameter_bookmark. ;
%%%%%%%%;
tmp_index = min(efind(strcmp(parameter.timing(:,1),str_field)));
if (~isempty(tmp_index));
parameter.timing{1+tmp_index,2} = parameter.timing{1+tmp_index,2}+dt;
parameter.timing{1+tmp_index,3} = parameter.timing{1+tmp_index,3}+ncall;
parameter.timing{1+tmp_index,4} = parameter.timing{1+tmp_index,4}+nop;
end;%if ( isempty(tmp_index));
if ( isempty(tmp_index));
tmp_index = size(parameter.timing,1);
parameter.timing{1+tmp_index,1} = str_field;
parameter.timing{1+tmp_index,2} = dt;
parameter.timing{1+tmp_index,3} = ncall;
parameter.timing{1+tmp_index,4} = nop;
end;%if ( isempty(tmp_index));
%%%%%%%%;
n_op = parameter.timing{1+tmp_index,4};
d_t = parameter.timing{1+tmp_index,2};
if ( (n_op>0) & (d_t>0) );
parameter.timing{1+tmp_index,5} = (n_op/d_t)/1e6;
parameter.timing{1+tmp_index,6} = (n_op/d_t)/1e9;
end;%if ( (n_op>0) & (d_t>0) );
%%%%%%%%;
