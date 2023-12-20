/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* The helper function */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
void transpose_psx4_block1_omp_helper
(
 float *A
 ,float *B
 ,const int n_row
 ,const int n_col
 ,const int lda
 ,const int ldb 
 ,const int block_size
 /* %%%%; */
 ,int nrbatch
 ,int n_rbatch
 ,int ncbatch
 ,int n_cbatch
 )
{
  int flag_omp = 0;
  int nc=0,nr=0,max_nc2=0,max_nr2=0,nc2=0,nr2=0;
  int nrow_per=0,nrow_start=0,nrow_final=0;
  int ncol_per=0,ncol_start=0,ncol_final=0;
  nrow_per = (int)ceil((double)n_row/maximum(1,(double)n_rbatch));
  nrow_start = maximum(0,minimum(n_row-1,(nrbatch+0)*nrow_per));
  nrow_final = maximum(0,minimum(n_row-1,(nrbatch+1)*nrow_per));
  ncol_per = (int)ceil((double)n_col/maximum(1,(double)n_cbatch));
  ncol_start = maximum(0,minimum(n_col-1,(ncbatch+0)*ncol_per));
  ncol_final = maximum(0,minimum(n_col-1,(ncbatch+1)*ncol_per));
  for(nc=ncol_start; nc<=ncol_final; nc+=block_size) {
    for(nr=nrow_start; nr<=nrow_final; nr+=block_size) {
      max_nc2 = nc+block_size < n_col ? nc + block_size : n_col;
      max_nr2 = nr+block_size < n_row ? nr + block_size : n_row;
      for(nc2=nc; nc2<max_nc2; nc2+=8) {
	for(nr2=nr; nr2<max_nr2; nr2+=8) {
	  transpose_psx4_block0(&A[nc2*lda +nr2], &B[nr2*ldb + nc2], lda, ldb);
	  /* for(nr2=nr; nr2<max_nr2; nr2+=8) { } */}
	/* for(nc2=nc; nc2<max_nc2; nc2+=8) { } */}
      /* for(nr=nrow_start; nr<=nrow_final; nr+=block_size) { } */}
    /* for(nc=ncol_start; nc<=ncol_final; nc+=block_size) { } */}
}

void transpose_psx4_block1_omp
(
 float *A
,float *B
,const int n_row
,const int n_col
,const int lda
,const int ldb 
,const int block_size
)
{
  int flag_omp = 1;
  int nc=0,nr=0,max_nc2=0,max_nr2=0,nc2=0,nr2=0;
  int nrbatch=0,n_r_per_rbatch=0,n_rbatch=0;
  int ncbatch=0,n_c_per_cbatch=0,n_cbatch=0;
  n_r_per_rbatch = n_row/2;
  n_rbatch = ceil((double)n_row/(double)n_r_per_rbatch);
  n_c_per_cbatch = n_col/4;
  n_cbatch = ceil((double)n_col/(double)n_c_per_cbatch);
  if (flag_omp==0){
    for (ncbatch=0;ncbatch<n_cbatch;ncbatch++){
      for (nrbatch=0;nrbatch<n_rbatch;nrbatch++){
	transpose_psx4_block1_omp_helper
	  (
	   A
	   ,B
	   ,n_row
	   ,n_col
	   ,lda
	   ,ldb 
	   ,block_size
	   ,nrbatch
	   ,n_rbatch
	   ,ncbatch
	   ,n_cbatch
	   );
	/* for (nrbatch=0;nrbatch<n_rbatch;nrbatch++){ } */}
      /* for (ncbatch=0;ncbatch<n_cbatch;ncbatch++){ } */}
    /* if (flag_omp==0){ } */}
  if (flag_omp==1){
#pragma omp parallel private(nrbatch,ncbatch)
    { /* begin omp parallel */
      nrbatch = 0; ncbatch = 0;
#pragma omp for schedule(dynamic)
      for (ncbatch=0;ncbatch<n_cbatch;ncbatch++){
	for (nrbatch=0;nrbatch<n_rbatch;nrbatch++){
	  transpose_psx4_block1_omp_helper
	    (
	     A
	     ,B
	     ,n_row
	     ,n_col
	     ,lda
	     ,ldb 
	     ,block_size
	     ,nrbatch
	     ,n_rbatch
	     ,ncbatch
	     ,n_cbatch
	     );
	  /* for (nrbatch=0;nrbatch<n_rbatch;nrbatch++){ } */}
	/* for (ncbatch=0;ncbatch<n_cbatch;ncbatch++){ } */}
      /* end omp parallel */}
    /* if (flag_omp==1){ } */}
}

