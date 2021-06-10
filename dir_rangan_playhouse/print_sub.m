function print_sub(input_,prefix);
if nargin<2; prefix = ''; end;
if numel(input_)==1;
disp(sprintf('%s%0.16f',prefix,input_(1)));
elseif numel(input_)==2;
disp(sprintf('%s%0.16f %0.16f',prefix,input_(1),input_(end)));
elseif numel(input_)>2; 
disp(sprintf('%s%0.16f %0.16f .[%d]. %0.16f %0.16f',prefix,input_(1),input_(2),numel(input_),input_(end-1),input_(end)));
end;%if;
