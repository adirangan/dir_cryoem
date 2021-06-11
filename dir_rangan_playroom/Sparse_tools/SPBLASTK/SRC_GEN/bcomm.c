#include "nsbbco.h"

void BCO_MatMult_CaABbC_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *bindx, const int *bjndx,
                 const int bnnz, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  const double *ptmp;                                  /* VECDEL */
  int i,j,jb,je,block;
  int cs,bs,br;
  int index;
  int l;                                               /* VECDEL */
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  for (l=0;l!=n;l++)                            /* VECDEL MMBETADEL */
    for (i=0;i!=m;i++) *pc++ *= beta;                  /* MMBETADEL */

  for (i=0;i!=bnnz;i++) {      /* For each block                       */
      cs = (bindx[i]-1)*lb;     /* start of c for this block operation */
      bs = (bjndx[i]-1)*lb;        /* start of b                       */
      pval = &val[i*mm];           /* pointer to start of block        */
      ptmp = pval;                 /* copy for later use     VECDEL    */
                                   /* GEMM call inlined here:          */
      pb = &b[bs];                 /* pointer to start of b            */
      pc = &c[cs];                 /* pointer to start of c            */
      for (l=0;l!=n;l++) {         /* For each column of c   VECDEL    */             
       for (jj=0;jj!=lb;jj++) {       /* For each column of block      */
         if( pb[jj] != 0.0 ) {           /* If non-zero multiplier     */
           t = alpha * pb[jj];            
           for (ii=0;ii!=lb;ii++) {      /* For each element in column */
             pc[ii] += t* (*pval++);        /* update element of c     */
           }
         } else {
           pval+=lb;                     /* Skip this column of block  */
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      /* to next column of c;VECDEL */
      }                                               /* VECDEL */
  }
}

void BCO_MatMult_CaATBbC_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *bindx, const int *bjndx, 
                 const int bnnz, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  const double *ptmp;                                  /* VECDEL */
  int i,j,jb,je,block;
  int cs,bs,cr;
  int index;
  int l;                                               /* VECDEL */
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  for (l=0;l!=n;l++)                            /* VECDEL MMBETADEL */
    for (i=0;i!=k;i++) *pc++ *= beta;                  /* MMBETADEL */

  for (i=0;i!=bnnz;i++) {      /* For each block                       */
      cs = (bjndx[i]-1)*lb;     /* start of c for this block operation */
      bs = (bindx[i]-1)*lb;        /* start of b                       */
      pval = &val[i*mm];           /* pointer to start of block        */
      ptmp = pval;                 /* copy for later use    VECDEL     */ 
                                   /* GEMM call inlined here:          */
      pb = &b[bs];                 /* pointer to start of b            */
      pc = &c[cs];                 /* pointer to start of c            */
      for (l=0;l!=n;l++) {         /* For each column of c  VECDEL     */     
       for (jj=0;jj!=lb;jj++) {       /* For each row    of block      */
           t = 0;
           for (ii=0;ii!=lb;ii++) {      /* accummulate dot with b     */
             t += alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  /* update element of c     */
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      /* to next column of c;VECDEL */
      }                                               /* VECDEL */
  } }

void BCOsymm_MatMult_CaABbC_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *bindx, const int *bjndx, 
                 const int bnnz, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  const double *ptmp; 
  int i,j,jb,je,block;
  int cs,css,cr,bs,br;
  int index;
  int l;                                               /* VECDEL */
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  for (l=0;l!=n;l++)                            /* VECDEL MMBETADEL */
    for (i=0;i!=m;i++) *pc++ *= beta;                  /* MMBETADEL */

  for (i=0;i!=bnnz;i++) {    /* For each block                         */
      cs = (bindx[i]-1)*lb;     /* start of c for this block operation */
      bs = (bjndx[i]-1)*lb;        /* start of b                       */
      pval = &val[i*mm];           /* pointer to start of block        */
      ptmp = pval;                 /* copy for later use               */
                                   /* GEMM call inlined here:          */
      pb = &b[bs];                 /* pointer to start of b            */
      pc = &c[cs];                 /* pointer to start of c            */
      for (l=0;l!=n;l++) {         /* For each column of c  VECDEL     */    
       for (jj=0;jj!=lb;jj++) {       /* For each column of block      */
         if( pb[jj] != 0.0 ) {           /* If non-zero multiplier     */
           t = alpha * pb[jj];            
           for (ii=0;ii!=lb;ii++) {      /* For each element in column */
             pc[ii] += t* (*pval++);        /* update element of c     */
           }
         } else {
           pval+=lb;                     /* Skip this column of block  */
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      /* to next column of c;VECDEL */
      }                                               /* VECDEL */
 /*   If on diagonal, done with block. Go to next.                     */
      if ( cs == bs ) {
        continue;                  
      }
 /*   Else, perform mult with transpose of this block                  */
      pval = ptmp;                 /* pointer to start of block        */
                                   /* GEMM call inlined here:          */
      pb = &b[cs];                 /* pointer to start of b            */
      pc = &c[bs];                 /* pointer to start of c            */
      for (l=0;l!=n;l++) {         /* For each column of c  VECDEL     */
       for (jj=0;jj!=lb;jj++) {       /* For each row    of block      */
           t = 0;
           for (ii=0;ii!=lb;ii++) {      /* accummulate dot with b     */
             t += alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  /* update element of c        */
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      /* to next column of c;VECDEL */
      }                                               /* VECDEL */
  }
}

void BCOskew_MatMult_CaABbC_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *bindx, const int *bjndx, 
                 const int bnnz, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int cs,css,cr,bs,br;
  int index;
  int l;                                               /* VECDEL */
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  for (l=0;l!=n;l++)                            /* VECDEL MMBETADEL */
    for (i=0;i!=m;i++) *pc++ *= beta;                  /* MMBETADEL */

  for (i=0;i!=bnnz;i++) {    /* For each block                         */
      cs = (bindx[i]-1)*lb;     /* start of c for this block operation */
      bs = (bjndx[i]-1)*lb;        /* start of b                       */
 /*   If on diagonal, assume zero.     Go to next.                     */
      if ( cs == bs ) {
        continue;                  
      }
      pval = &val[i*mm];           /* pointer to start of block        */
      ptmp = pval;                 /* copy for later use               */
                                   /* GEMM call inlined here:          */
      pb = &b[bs];                 /* pointer to start of b            */
      pc = &c[cs];                 /* pointer to start of c            */
      for (l=0;l!=n;l++) {         /* For each column of c  VECDEL     */ 
       for (jj=0;jj!=lb;jj++) {       /* For each column of block      */
         if( pb[jj] != 0.0 ) {           /* If non-zero multiplier     */
           t = alpha * pb[jj];            
           for (ii=0;ii!=lb;ii++) {      /* For each element in column */
             pc[ii] += t* (*pval++);        /* update element of c     */
           }
         } else {
           pval+=lb;                     /* Skip this column of block  */
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      /* to next column of c;VECDEL */
      }                                               /* VECDEL */
 /*   Else, perform mult with transpose of this block                  */
      pval = ptmp;                 /* pointer to start of block        */
                                   /* GEMM call inlined here:          */
      pb = &b[cs];                 /* pointer to start of b            */
      pc = &c[bs];                 /* pointer to start of c            */
      for (l=0;l!=n;l++) {         /* For each column of c  VECDEL     */
       for (jj=0;jj!=lb;jj++) {       /* For each row    of block      */
           t = 0;
           for (ii=0;ii!=lb;ii++) {      /* accummulate dot with b     */
             t -= alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  /* update element of c        */
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      /* to next column of c;VECDEL */
      }                                               /* VECDEL */
  }
}

void BCOskew_MatMult_CaATBbC_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *bindx, const int *bjndx,
                 const int bnnz, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  int i,j,jb,je,block;
  int cs,css,cr,bs,br;
  int index;
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;
  const double *ptmp;
  int l;                                               /* VECDEL */

  for (l=0;l!=n;l++)                            /* VECDEL MMBETADEL */
    for (i=0;i!=k;i++) *pc++ *= beta;                  /* MMBETADEL */

  for (i=0;i!=bnnz;i++) {    /* For each block                         */
      cs = (bindx[i]-1)*lb;     /* start of c for this block operation */
      bs = (bjndx[i]-1)*lb;        /* start of b                       */
 /*   If on diagonal, done with block. Go to next.                     */
      if ( cs == bs ) {
        continue;                  
      }
      pval = &val[i*mm];           /* pointer to start of block        */
      ptmp = pval;                 /* copy for later use               */
                                   /* GEMM call inlined here:          */
      pb = &b[bs];                 /* pointer to start of b            */
      pc = &c[cs];                 /* pointer to start of c            */
      for (l=0;l!=n;l++) {         /* For each column of c  VECDEL     */    
       for (jj=0;jj!=lb;jj++) {       /* For each column of block      */
         if( pb[jj] != 0.0 ) {           /* If non-zero multiplier     */
           t = alpha * pb[jj];            
           for (ii=0;ii!=lb;ii++) {      /* For each element in column */
             pc[ii] -= t* (*pval++);        /* update element of c     */
           }
         } else {
           pval+=lb;                     /* Skip this column of block  */
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      /* to next column of c;VECDEL */
      }                                               /* VECDEL */
 /*   Else, perform mult with transpose of this block                  */
      pval = ptmp;                 /* pointer to start of block        */
                                   /* GEMM call inlined here:          */
      pb = &b[cs];                 /* pointer to start of b            */
      pc = &c[bs];                 /* pointer to start of c            */
      for (l=0;l!=n;l++) {         /* For each column of c  VECDEL     */
       for (jj=0;jj!=lb;jj++) {       /* For each row    of block      */
           t = 0;
           for (ii=0;ii!=lb;ii++) {      /* accummulate dot with b     */
             t += alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  /* update element of c        */
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      /* to next column of c;VECDEL */
      }                                               /* VECDEL */
  }
}

