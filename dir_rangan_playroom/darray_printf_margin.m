% void darray_printf_margin(double *d_,int n_r,int n_c,const char *prefix)
% {
%   int nr=0,nc=0,margin=3;
%   for (nr=0;nr<n_r;nr++){
%     if ( (nr>margin-1) && (nr<n_r-margin) ){ printf("%s...............................................................\n",prefix); nr = maximum(0,n_r-margin);}
%     if ( (nr<margin) || (nr>n_r-margin-1) ){
%       printf("%s",prefix);
%       for (nc=0;nc<n_c;nc++){
% 	if ( (nc>margin-1) && (nc<n_c-margin) ){ printf("... "); nc = maximum(0,n_c-margin);}
% 	if ( (nc<margin) || (nc>n_c-margin-1) ){
% 	  printf("%+.6f ",d_[nr+nc*n_r]);
% 	  /* margin */}
% 	/* for (nc=0;nc<n_c;nc++){ } */}
%       printf("\n");
%       /* margin */}
%     /* for (nr=0;nr<n_r;nr++){ } */}
% }

function darray_printf_margin(d_,n_r,n_c,prefix);
n_margin = 3;
nr=0; nc=0;
nr=0;
while(nr<n_r);
if ( (nr>n_margin-1) & (nr<n_r-n_margin) );
fprintf(1,"%s...............................................................\n",prefix);
nr = max(0,n_r-n_margin);
end;%if ( (nr>n_margin-1) & (nr<n_r-n_margin) );
if ( (nr<n_margin) | (nr>n_r-n_margin-1) );
fprintf(1,"%s",prefix);
nc=0;
while(nc<n_c);
if ( (nc>n_margin-1) & (nc<n_c-n_margin) );
fprintf(1,"... "); nc = max(0,n_c-n_margin);
end;%if ( (nc>n_margin-1) & (nc<n_c-n_margin) );
if ( (nc<n_margin) | (nc>n_c-n_margin-1) );
fprintf(1,"%+.6f ",d_(1+nr+nc*n_r));
end;%if ( (nc<n_margin) | (nc>n_c-n_margin-1) );
nc=nc+1;
end;%while(nc<n_c);
fprintf(1,"\n");
end;%if ( (nr<n_margin) | (nr>n_r-n_margin-1) );
nr=nr+1;
end;%while(nr<n_r);

