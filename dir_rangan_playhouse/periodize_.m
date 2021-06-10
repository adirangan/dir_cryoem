function B_ = periodize_(A_,min,max);
B_ = A_;
%%%%%%%%;
flag_continue=1;
while flag_continue;
tmp_ij = find(B_<min);
if length(tmp_ij)==0; flag_continue = 0; else flag_continue=1; B_(tmp_ij) = B_(tmp_ij) + (max-min); end;
end;%while flag_continue;
%%%%%%%%;
flag_continue=1;
while flag_continue;
tmp_ij = find(B_>=max);
if length(tmp_ij)==0; flag_continue = 0; else flag_continue=1; B_(tmp_ij) = B_(tmp_ij) - (max-min); end;
end;%while flag_continue;

