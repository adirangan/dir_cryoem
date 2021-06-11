#ifndef _MONOLITH
#include "playpark_header.h"
#endif /* _MONOLITH */

void get_xdrop_logscale(double n_row,double n_col,double gamma,int *rdrop_p,int *cdrop_p)
{
  int verbose=0; double dbl1=1.000000-0.000001;
  double gamma_tmp_row=0,ammag_tmp_row=0,gamma_tmp_col=0,ammag_tmp_col=0; int rdrop=0;int cdrop=0;
  if (verbose){ printf(" %% [entering get_xdrop_logscale] gamma %0.3f, n_row %d n_col %d\n",gamma,(int)n_row,(int)n_col);}
  /*
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if gamma=0, then we set up gamma_tmp_row to remove a single row ;
    if (gamma>0); gamma_tmp_row = max(gamma,(1-1e-6)/length(tmp_rij)); ammag_tmp_row = min(1-gamma,(length(tmp_rij)-1)/length(tmp_rij));
    elseif gamma<=0; gamma_tmp_row = (1-1e-6)/length(tmp_rij); ammag_tmp_row = (length(tmp_rij)-1)/length(tmp_rij);
    end;%if gamma==0;
    % setting up ammag_tmp_col to remove as many cols as necessary so that log(ncols_pos)/log(nrows_pos) = log(ncols_pre)/log(nrows_pre) ;
    % i.e., log(ammag_tmp_col*ncols_pre)/log(ammag_tmp_row*nrows_pre) = log(ncols_pre)/log(nrows_pre) ;
    % i.e., log(ammag_tmp_col*ncols_pre) = (log(ammag_tmp_row) + log(nrows_pre))*log(ncols_pre)/log(nrows_pre) ;
    % i.e., log(ammag_tmp_col) = log(ammag_tmp_row)*log(ncols_pre)/log(nrows_pre) ;
    % i.e., ammag_tmp_col = exp(log(ammag_tmp_row)*log(ncols_pre)/log(nrows_pre));
    ammag_tmp_col = exp(log(ammag_tmp_row)*log(length(tmp_cij))/log(length(tmp_rij)));
    gamma_tmp_col = 1-ammag_tmp_col-1e-6;
    rdrop = r_rem{iteration}(tmp_rij(1:ceil(gamma_tmp_row*end))); cdrop = c_rem{iteration}(tmp_cij(1:ceil(gamma_tmp_col*end)));
    rkeep = r_rem{iteration}(tmp_rij(ceil(gamma_tmp_row*end)+1:end)); ckeep = c_rem{iteration}(tmp_cij(ceil(gamma_tmp_col*end)+1:end));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Note that we can perform similar row and col removal using sqrt or linear scaling, e.g., :
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    M=5e3;N=300e3;g=0.1;nr=M;nc=N;ii=1;nr_(ii)=nr;nc_(ii)=nc;
    while(nr>1); gr=g;ar=min(1-g,(nr-1)/nr);ac=exp(log(ar)*log(nc)/log(nr));gc=1-1e-6-ac;
    rdrop=ceil(gr*nr);cdrop=ceil(gc*nc);nr=nr-rdrop;nc=nc-cdrop;ii=ii+1;nr_(ii)=nr;nc_(ii)=nc;end; % log-scaling ;
    figure;subplot(2,3,1);plot(nr_,nc_,'o-');subplot(2,3,2);plot(log(nr_),log(nc_),'o-'); subplot(2,3,3);plot(sqrt(nr_),sqrt(nc_),'o-'); 
    subplot(2,3,4);plot(1:ii,nr_./nc_,'o-');subplot(2,3,5);plot(1:ii,log(nr_)./log(nc_),'o-');subplot(2,3,6);plot(1:ii,sqrt(nr_)./sqrt(nc_),'o-'); 
    M=5e3;N=300e3;g=0.1;nr=M;nc=N;ii=1;nr_(ii)=nr;nc_(ii)=nc;
    while(nr>1); gr=g;gc=g;
    rdrop=ceil(gr*nr);cdrop=ceil(gc*nc);nr=nr-rdrop;nc=nc-cdrop;ii=ii+1;nr_(ii)=nr;nc_(ii)=nc;end; % linear scaling ;
    figure;subplot(2,3,1);plot(nr_,nc_,'o-');subplot(2,3,2);plot(log(nr_),log(nc_),'o-'); subplot(2,3,3);plot(sqrt(nr_),sqrt(nc_),'o-'); 
    subplot(2,3,4);plot(1:ii,nr_./nc_,'o-');subplot(2,3,5);plot(1:ii,log(nr_)./log(nc_),'o-');subplot(2,3,6);plot(1:ii,sqrt(nr_)./sqrt(nc_),'o-'); 
    M=5e3;N=300e3;g=0.1;nr=M;nc=N;ii=1;nr_(ii)=nr;nc_(ii)=nc;
    while(nr>1); gr=g;ar=min(1-g,(nr-1)/nr);ac=(1/nc)*(sqrt(ar*nr)*sqrt(nc)/sqrt(nr)).^2;gc=1-1e-6-ac;
    rdrop=ceil(gr*nr);cdrop=ceil(gc*nc);nr=nr-rdrop;nc=nc-cdrop;ii=ii+1;nr_(ii)=nr;nc_(ii)=nc;end; % sqrt-scaling ;
    figure;subplot(2,3,1);plot(nr_,nc_,'o-');subplot(2,3,2);plot(log(nr_),log(nc_),'o-'); subplot(2,3,3);plot(sqrt(nr_),sqrt(nc_),'o-'); 
    subplot(2,3,4);plot(1:ii,nr_./nc_,'o-');subplot(2,3,5);plot(1:ii,log(nr_)./log(nc_),'o-');subplot(2,3,6);plot(1:ii,sqrt(nr_)./sqrt(nc_),'o-'); 
  */
  if (n_row>n_col){ get_xdrop_logscale(n_col,n_row,gamma,cdrop_p,rdrop_p);}
  else /* if (n_row<=n_col) */{
  if (n_row<=2){ /* drop everything */ rdrop = n_row; cdrop = n_col;}
  else /* if n_row>2 */{
    if (gamma>0){ gamma_tmp_row = maximum(gamma,(dbl1)/maximum(1,n_row)); ammag_tmp_row = minimum(1-gamma,(n_row-1)/maximum(1,n_row));}
    else /* if gamma<=0 */{ gamma_tmp_row = (dbl1)/maximum(1,n_row); ammag_tmp_row = (n_row-1)/maximum(1,n_row);}
    ammag_tmp_col = exp(log(ammag_tmp_row)*log(n_col)/maximum(1,log(maximum(1,n_row)))); gamma_tmp_col = dbl1-ammag_tmp_col; 
    rdrop = ceil(gamma_tmp_row*n_row); cdrop = ceil(gamma_tmp_col*n_col);
    /* if n_row>2 */}
  rdrop = minimum(rdrop,n_row); cdrop = minimum(cdrop,n_col);
  if (verbose>0){ printf(" %% gamma_tmp_row %0.3f ammag_tmp_row %0.3f rdrop %d/%d gamma_tmp_col %0.3f ammag_tmp_col %0.3f cdrop %d/%d\n",gamma_tmp_row,ammag_tmp_row,rdrop,(int)n_row,gamma_tmp_col,ammag_tmp_col,cdrop,(int)n_col);}
  if (rdrop_p!=NULL){ *rdrop_p=rdrop;} if (cdrop_p!=NULL){ *cdrop_p=cdrop;}
  if (verbose){ printf(" %% [finished get_xdrop_logscale] rdrop %d cdrop %d\n",rdrop,cdrop);}
  /* if n_row<=n_col */}
}

int get_xdrop_logscale_length(double n_row,double n_col,double gamma)
{
  int verbose=0;
  double n_row_tmp=n_row,n_col_tmp=n_col;
  int rdrop=0,cdrop=0;
  int continue_flag=1;
  int length=0;
  if (verbose){ printf(" %% [entering get_xdrop_logscale_length] n_row %0.2f n_col %0.2f\n",n_row,n_col);}
  while (continue_flag){
    length++;
    get_xdrop_logscale(n_row_tmp,n_col_tmp,gamma,&rdrop,&cdrop);
    n_row_tmp-=rdrop; n_col_tmp-=cdrop;
    if (verbose>1){ printf(" %% length %d, rdrop %d cdrop %d n_row_tmp %0.2f n_col_tmp %0.2f\n",length,rdrop,cdrop,n_row_tmp,n_col_tmp);}
    continue_flag = (n_row_tmp>0 && n_col_tmp>0);
    /* while (continue_flag){ } */}
  if (verbose){ printf(" %% [finished get_xdrop_logscale_length]\n");}
  return length;
}

void get_xdrop_logscale_array(double n_row,double n_col,double gamma,int *length_p,int **rdrop_p_,int **cdrop_p_,int **rkeep_p_,int **ckeep_p_)
{
  /* Stores rdrop and cdrop values across all iterations.
     Also stores number remaining (rkeep and ckeep).
  */
  int verbose=0;
  int length=0,nl=0,rdrop=0,cdrop=0;
  double n_row_tmp=0,n_col_tmp=0;
  if (verbose){ printf(" %% [entering get_xdrop_logscale_array] n_row %d n_col %d\n",(int)n_row,(int)n_col);}
  length = get_xdrop_logscale_length(n_row,n_col,gamma);
  if (length_p!=NULL){ *length_p = length;}
  if (rdrop_p_!=NULL){ if (*rdrop_p_==NULL){ (*rdrop_p_) = (int *) malloc1(length*sizeof(int));}}
  if (cdrop_p_!=NULL){ if (*cdrop_p_==NULL){ (*cdrop_p_) = (int *) malloc1(length*sizeof(int));}}
  if (rkeep_p_!=NULL){ if (*rkeep_p_==NULL){ (*rkeep_p_) = (int *) malloc1(length*sizeof(int));}}
  if (ckeep_p_!=NULL){ if (*ckeep_p_==NULL){ (*ckeep_p_) = (int *) malloc1(length*sizeof(int));}}
  n_row_tmp=n_row; n_col_tmp=n_col;
  for (nl=0;nl<length;nl++){
    get_xdrop_logscale(n_row_tmp,n_col_tmp,gamma,&rdrop,&cdrop);
    if (rdrop_p_!=NULL){ (*rdrop_p_)[nl] = rdrop;}
    if (cdrop_p_!=NULL){ (*cdrop_p_)[nl] = cdrop;}
    n_row_tmp-=rdrop; n_col_tmp-=cdrop;
    if (rkeep_p_!=NULL){ (*rkeep_p_)[nl] = round(n_row_tmp);}
    if (ckeep_p_!=NULL){ (*ckeep_p_)[nl] = round(n_col_tmp);}
    /* for (nl=0;nl<length;nl++){ } */}  
  if (verbose>1){ if (rdrop_p_!=NULL){ array_printf((*rdrop_p_),"int",1,length," %% rdrop_: ");}}
  if (verbose>1){ if (rkeep_p_!=NULL){ array_printf((*rkeep_p_),"int",1,length," %% rkeep_: ");}}
  if (verbose>1){ if (cdrop_p_!=NULL){ array_printf((*cdrop_p_),"int",1,length," %% cdrop_: ");}}
  if (verbose>1){ if (ckeep_p_!=NULL){ array_printf((*ckeep_p_),"int",1,length," %% ckeep_: ");}}
  if (verbose){ printf(" %% [finished get_xdrop_logscale_array]\n");}  
}

void get_xdrop_logscale_array_test()
{
  int n_row = 1800, n_col = 5e4;
  double gamma = 0.125;
  int length=0;
  int *rdrop_ = NULL;
  int *cdrop_ = NULL;
  int *rkeep_ = NULL;
  int *ckeep_ = NULL;
  int rdrop_ans_[45] = { 225 , 197 , 173 , 151 , 132 , 116 , 101 , 89 , 77 , 68 , 59 , 52 , 45 , 40 , 35 , 30 , 27 , 23 , 20 , 18 , 16 , 14 , 12 , 10 , 9 , 8 , 7 , 6 , 5 , 5 , 4 , 4 , 3 , 3 , 2 , 2 , 2 , 2 , 1 , 1 , 1 , 1 , 1 , 1 , 2 };
  int cdrop_ans_[45] = { 8766 , 7229 , 5962 , 4917 , 4055 , 3344 , 2758 , 2275 , 1876 , 1547 , 1276 , 1052 , 868 , 716 , 590 , 487 , 401 , 331 , 273 , 225 , 186 , 153 , 126 , 104 , 86 , 71 , 58 , 48 , 40 , 32 , 27 , 22 , 18 , 15 , 13 , 10 , 8 , 7 , 6 , 5 , 5 , 4 , 3 , 3 , 2 };
  printf(" %% [entering get_xdrop_logscale_array_test]\n");
  GLOBAL_tic(0);
  get_xdrop_logscale_array(n_row,n_col,gamma,&length,&rdrop_,&cdrop_,&rkeep_,&ckeep_);
  printf(" %% n_row %d n_col %d length %d\n",n_row,n_col,length);
  array_printf(rdrop_,"int",1,minimum(16,length)," rdrop_ start: ");
  array_printf(cdrop_,"int",1,minimum(16,length)," cdrop_ start: ");
  array_printf(rkeep_,"int",1,minimum(16,length)," rkeep_ start: ");
  array_printf(ckeep_,"int",1,minimum(16,length)," ckeep_ start: ");
  array_printf(&(rdrop_[length-minimum(16,length)]),"int",1,minimum(16,length)," rdrop_ final: ");
  array_printf(&(cdrop_[length-minimum(16,length)]),"int",1,minimum(16,length)," cdrop_ final: ");
  array_printf(&(rkeep_[length-minimum(16,length)]),"int",1,minimum(16,length)," rkeep_ final: ");
  array_printf(&(ckeep_[length-minimum(16,length)]),"int",1,minimum(16,length)," ckeep_ final: ");
  printf(" %% rdrop_ans_ vs rdrop_: relative error %0.16f\n",ifnormn(length,rdrop_ans_,rdrop_));
  printf(" %% cdrop_ans_ vs cdrop_: relative error %0.16f\n",ifnormn(length,cdrop_ans_,cdrop_));
  free1(&rdrop_);
  free1(&cdrop_);
  free1(&rkeep_);
  free1(&ckeep_);
  GLOBAL_toc(0,1," get_xdrop_logscale_array: ");
  printf(" %% [finished get_xdrop_logscale_array_test]\n");
}
