function array_printf(v_,n_row,n_col,prefix);
printftol = 1e-9;
tmp_string = '';

if isreal(v_);
for nr=0:n_row-1;
tmp_string = sprintf('%s%s',tmp_string,prefix);
for nc=0:n_col-1;
if abs(v_(1+nr+nc*n_row))< printftol;
tmp_string = sprintf('%s%s',tmp_string,'    .   ');
end;%if abs(v_(1+nr+nc*n_row))< printftol;
if abs(v_(1+nr+nc*n_row))>=printftol;
tmp_string = sprintf('%s%s',tmp_string,sprintf(' %+07.3f',v_(1+nr+nc*n_row)));
end;%if abs(v_(1+nr+nc*n_row))>=printftol;
end;%for nc=0:n_col-1;
tmp_string = sprintf('%s%s',tmp_string,sprintf('\n'));
end;%for nr=0:n_row-1;
end;%if isreal(v_);

if ~isreal(v_);
for nr=0:n_row-1;
tmp_string = sprintf('%s%s',tmp_string,prefix);
for nc=0:n_col-1;
if abs(v_(1+nr+nc*n_row))< printftol;
tmp_string = sprintf('%s%s',tmp_string,'    .        .   ');
end;%if abs(v_(1+nr+nc*n_row))< printftol;
if abs(v_(1+nr+nc*n_row))>=printftol;
tmp_string = sprintf('%s%s',tmp_string,sprintf(' [%+07.3f %+07.3fi]',real(v_(1+nr+nc*n_row)),imag(v_(1+nr+nc*n_row))));
end;%if abs(v_(1+nr+nc*n_row))>=printftol;
end;%for nc=0:n_col-1;
tmp_string = sprintf('%s%s',tmp_string,sprintf('\n'));
end;%for nr=0:n_row-1;
end;%if ~isreal(v_);

disp(tmp_string);
