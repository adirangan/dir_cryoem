function export_f77dat(n_A,A,str_type,vname,fname,headerline);
% exports array A of length n_A as fortran data of type str_type ;
if (strcmp(str_type,'integer') | strcmp(str_type,'int')) ; 
fp = fopen(fname,'w');
fprintf(fp,'c$$$      %s\n',headerline);
fprintf(fp,'      integer %s(0:%d)\n',vname,n_A-1);
nA=0;
fprintf(fp,'      data %s(%d)/ %d/\n',vname,nA,A(1+nA));
for nA=1:n_A-1;
fprintf(fp,'     $     , %s(%d)/ %d/\n',vname,nA,A(1+nA));
end;%for nA=1:n_A-1;
fclose(fp);
elseif (strcmp(str_type,'real') | strcmp(str_type,'double')) ; 
fp = fopen(fname,'w');
fprintf(fp,'c$$$      %s\n',headerline);
fprintf(fp,'      real *8 %s(0:%d)\n',vname,n_A-1);
nA=0;
fprintf(fp,'      data %s(%d)/ %0.15f/\n',vname,nA,A(1+nA));
for nA=1:n_A-1;
fprintf(fp,'     $     , %s(%d)/ %0.15f/\n',vname,nA,A(1+nA));
end;%for nA=1:n_A-1;
fclose(fp);
elseif (strcmp(str_type,'complex') | strcmp(str_type,'cmplx')) ; 
fp = fopen(fname,'w');
fprintf(fp,'c$$$      %s\n',headerline);
fprintf(fp,'      complex *16 %s(0:%d)\n',vname,n_A-1);
nA=0;
fprintf(fp,'      data %s(%d)/ ( %0.15f , %0.15f )/\n',vname,nA,real(A(1+nA)),imag(A(1+nA)));
for nA=1:n_A-1;
fprintf(fp,'     $     , %s(%d)/ ( %0.15f , %0.15f )/\n',vname,nA,real(A(1+nA)),imag(A(1+nA)));
end;%for nA=1:n_A-1;
fclose(fp);
 else;
disp(sprintf(' %% Warning! str_type %s not supported in export_f77dat; fname %s not written',str_type,fname));
end% if;
