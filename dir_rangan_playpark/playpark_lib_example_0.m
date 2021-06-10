%playpark_lib_load;
%libfunctionsview('playpark_lib');

x_ = cast(1:3,'uint8'); n_x = numel(x_);
calllib('playpark_lib','array_printf',x_,'char',1,n_x,' uint8: ');
x_ = cast(-3:3,'int32'); n_x = numel(x_);
calllib('playpark_lib','array_printf',x_,'int',1,n_x,' int32: ');
x_ = cast(0:13,'uint32'); n_x = numel(x_);
calllib('playpark_lib','array_printf',x_,'unsigned int',1,n_x,' uint32: ');
x_ = cast(0:13,'uint64'); n_x = numel(x_);
calllib('playpark_lib','array_printf',x_,'unsigned long long int',1,n_x,' uint64: ');
x_ = cast(-1:0.5:3,'double'); n_x = numel(x_);
calllib('playpark_lib','array_printf',x_,'double',1,n_x,' double: ');
x_ = cast(-1:0.5:3,'single'); n_x = numel(x_);
calllib('playpark_lib','array_printf',x_,'float',1,n_x,' single: ');
x_ = (0:3) + i*(5:8); n_x = numel(x_);  %<-- complex. ; 
calllib('playpark_lib','array_printf',[real(x_);imag(x_)],'double',1,2*n_x,' complex: ');
