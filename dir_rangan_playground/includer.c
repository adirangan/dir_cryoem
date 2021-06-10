//! Doxygen comment: ;\n
//! Simple program to copy an input file, while replacing "include 'filename'" statements with the filename. ;\n
//! This is designed especially for use with fortran files. ;\n
//! Compile and test with: gcc -o includer.out includer.c; ./includer.out --input testfile.f --ignore_h --auto_y; cat testfile.included.f ;\n
//! Example: gcc -o includer.out includer.c; ./includer.out --input test_innerproduct_8.f --ignore_h --auto_y ;\n
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#define FNAMESIZE (4096+256)
#define RETURN_SUCCESS 0
#define RETURN_FAIL 1

int IGNORE_H = 0;
int AUTO_Y = 0;
char str_help[256] = "try: ./includer.out --input name_of_input_file [--output name_of_output_file] [--auto_y] [--ignore_h] \nFlag auto_y sets automatic inclusion.\nFlag ignore_h will skip inclusions of the form \"include 'filename.h'\".\n ";

int plugin(FILE *FP_OUT,char *str_input)
{
  FILE *FP_IN=NULL;
  char *tmp_str=NULL;
  size_t tmp_len=0;
  ssize_t tmp_read=0;
  int nl=0,nc=0,flag_continue=0;
  char *tmp_tok=NULL,tmp_c[4];
  char *str_bkp=NULL;
  int flag_header=0;
  printf(" %% [entering plugin] %s\n",str_input);
  FP_IN = fopen(str_input,"r");
  if (FP_IN==NULL){ printf(" %% Warning! input %s not found.\n",str_input); exit(RETURN_FAIL);}
  nl=0;
  while ((tmp_read = getline(&tmp_str, &tmp_len, FP_IN)) != -1){
    printf(" %% line %d of length %zu:\n",nl,tmp_read);
    str_bkp = malloc(strlen(tmp_str)+4);
    sprintf(str_bkp,tmp_str);
    if (strlen(tmp_str)>=6 && !strncmp(tmp_str,"      ",6) && strstr(tmp_str,"include \'")){
      tmp_tok = strtok(tmp_str,"\'");
      tmp_tok = strtok(NULL,"\'");
      if (tmp_tok!=NULL){ 
	flag_header = 0;
	if (strlen(tmp_tok)>1){ 
	  flag_header = (tmp_tok[strlen(tmp_tok)-2]=='.') && (tmp_tok[strlen(tmp_tok)-1]=='h' || tmp_tok[strlen(tmp_tok)-1]=='H');
	  /* if (strlen(tmp_tok)>1){ } */}
	if (IGNORE_H==1 && flag_header==1){ fprintf(FP_OUT,str_bkp);}
	else /* if not a header */{
	  if (AUTO_Y==1){ printf(" %% automatically including %s ...\n",tmp_tok); plugin(FP_OUT,tmp_tok);}
	  if (AUTO_Y==0){
	    printf(" %% Would you like to include token %s [y/n]?\n",tmp_tok);
	    fgets(tmp_c,3,stdin);
	    if (!strstr(tmp_c,"n")){ printf(" %% you entered: %s, including %s ...\n",tmp_c,tmp_tok); plugin(FP_OUT,tmp_tok);}
	    else{ printf(" %% you entered: %s, not including %s\n",tmp_c,tmp_tok);}
	    /* if (AUTO_Y==0){ } */}
	  /* if not a header */}
	/* if (tmp_tok!=NULL){ } */}
      /* if include */}
    else /* if no include */{
      fprintf(FP_OUT,tmp_str);
      /* else if no include { } */}
    if (str_bkp){ free(str_bkp); str_bkp=NULL;}
    nl++; /* while ((tmp_read = getline(&tmp_str, &tmp_len, FP_IN)) != -1){ } */}
  fclose(FP_IN);
  if (tmp_str){ free(tmp_str);}
  printf(" %% [finished plugin] %s\n",str_input);
}

int main(int argc, char** argv)
{
  int cur_arg = 1;
  char* argptr;
  char str_input[FNAMESIZE];
  char str_output[FNAMESIZE];
  FILE *FP_OUT=NULL;
  char str_bkp[FNAMESIZE];
  char *tmp_tok_pre=NULL;
  char *tmp_tok_pos=NULL;
  printf(" %% [entering includer]\n");
  while (cur_arg < argc) {
    argptr = argv[cur_arg];
    if (0){ /* do nothing */}
    else if (!strcmp(argptr, "--ignore_h")) {
      IGNORE_H = 1;
      printf(" %% setting IGNORE_H: %d\n",IGNORE_H);
      cur_arg += 1; /* else if ignore_h */}
    else if (!strcmp(argptr, "--auto_y")) {
      AUTO_Y = 1;
      printf(" %% setting AUTO_Y: %d\n",AUTO_Y);
      cur_arg += 1; /* else if auto_y */}
    else if (!strcmp(argptr, "--input")) {
      if (cur_arg == argc - 1) { printf("Error: Missing --input.\n"); exit(RETURN_FAIL);}
      sprintf(str_input,argv[cur_arg+1]);
      printf(" %% input: %s\n",str_input);
      sprintf(str_bkp,str_input);
      tmp_tok_pre = strtok(str_bkp,".");
      tmp_tok_pos = strtok(NULL,".");
      sprintf(str_output,"%s.included.%s",tmp_tok_pre,tmp_tok_pos);
      printf(" %% setting output: %s\n",str_output);
      cur_arg += 2; /* else if input */}
    else if (!strcmp(argptr, "--output")) {
      if (cur_arg == argc - 1) { printf("Error: Missing --output.\n"); exit(RETURN_FAIL);}
      sprintf(str_output,argv[cur_arg+1]);
      printf(" %% renaming output: %s\n",str_output);
      cur_arg += 2; /* else if output */}
    else { printf("Error: Invalid argument (%s).\n", argv[cur_arg]); printf(str_help); exit(RETURN_FAIL);}
    /* while (cur_arg < argc) { } */}
  if (!strcmp(str_input,str_output)){ printf(" %% Warning! input and output both named %s. Exiting.\n",str_input,str_output); exit(RETURN_FAIL); }
  FP_OUT = fopen(str_output,"w");
  if (FP_OUT==NULL){ printf(" %% Warning! output %s not created.\n",str_output); exit(RETURN_FAIL);}
  plugin(FP_OUT,str_input);
  fclose(FP_OUT);
  printf(" %% [finished includer]\n");
  return RETURN_SUCCESS;
}
