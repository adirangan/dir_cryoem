#ifndef _MONOLITH
#include "playpark_header.h"
#endif /* _MONOLITH */

void updateglobals(char *vname)
{
  /* given a variable name, scan the appropriate format into the appropriate global variable */
  int verbose=GLOBAL_verbose;
  char comma_vs_semicolon[1],tmpchar[FNAMESIZE]; int length=0,nv=0;
  if (0){ /* do nothing */ }
  else if (strcmp(vname,"GLOBAL_verbose")==0){ scanf("%d",&GLOBAL_verbose); if (verbose>0){ printf("%s read to be %d\n",vname,GLOBAL_verbose);}}
  else if (strcmp(vname,"GLOBAL_thread_count")==0){ 
    scanf("%d",&GLOBAL_thread_count); if (verbose>0){ printf("%s read to be %d\n",vname,GLOBAL_thread_count);}
    length = GLOBAL_thread_count; 
    /* else if (strcmp(vname,"GLOBAL_thread_count")==0){ } */}
  else if (strcmp(vname,"GLOBAL_omp_type")==0){ scanf("%d",&GLOBAL_omp_type); if (verbose>0){ printf("%s read to be %d\n",vname,GLOBAL_omp_type);}}
  else if (strcmp(vname,"GLOBAL_dir_name")==0){ scanf("%[^,;]",GLOBAL_dir_name); if (verbose>0){ printf("%s read to be %s\n",vname,GLOBAL_dir_name);} /* else if (strcmp(vname,"GLOBAL_dir_name")==0){ } */}
  else if (strcmp(vname,"END")==0){ /* do nothing */ if (verbose>0){ printf("end of input reached\n");}}
/*   else if (strcmp(vname,"yy")==0){ scanf("%zz",&yy); if (verbose>0){ printf("%s read to be %zz\n",vname,yy);}} */
  else /* if anything else */{ printf(" %% Error! vname %s in updateglobals\n",vname); exit(RET_READ_FAIL); }
}

void read_input()
{
  /* This reads the piped input file, or standard input.
     the variable names should not be longer than 128 characters */
  int verbose=GLOBAL_verbose;
  char vname[128],equals[128],space[128],semicolon[128];
  do{
    scanf("%[^=]",vname);scanf("%s",equals);scanf("%c",space);updateglobals(vname);scanf("%c",semicolon);
    if (verbose>1){ printf("At this point variable name is (%s), equals is (%s), semicolon is (%s)\n",vname,equals,semicolon);} 
    scanf("%c",semicolon);}
  while (strcmp(vname,"END")!=0);
}
