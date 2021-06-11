#include "nsbbsr.h"

void BSR_MatMult_CaABbC_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
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

  bindx-=ind_base;

  for (l=0;l!=n;l++)                            /* VECDEL MMBETADEL */
    for (i=0;i!=m;i++) *pc++ *= beta;                  /* MMBETADEL */

  for (i=0;i!=mb;i++) {      /* For each block row i                   */
    jb = bpntrb[i];             /* beginning block in row              */
    je = bpntre[i];             /* ending block in row                 */
    cs = i*lb;                  /* start of c for this block operation */
    for (j=jb;j!=je;j++) {      /* For each block j in this row        */
      index = bindx[j];
      bs=(index-1)*lb;             /* start of b                       */
      pval = &val[(j-1)*mm];       /* pointer to start of block        */
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
}

void BSR_MatMult_CaATBbC_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
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

  bindx-=ind_base;

  for (l=0;l!=n;l++)                            /* VECDEL MMBETADEL */
    for (i=0;i!=k;i++) *pc++ *= beta;                  /* MMBETADEL */

  for (i=0;i!=mb;i++) {      /* For each block column i                */
    jb = bpntrb[i];             /* beginning block in column           */
    je = bpntre[i];             /* ending block in column              */
    bs = i*lb;                  /* start of b for this block operation */
    for (j=jb;j!=je;j++) {      /* For each block j in this column     */
      index = bindx[j];
      cs=(index-1)*lb;             /* start of b                       */
      pval = &val[(j-1)*mm];       /* pointer to start of block        */
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
    }
  }
}

void BSRsymm_MatMult_CaABbC_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
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

  bindx-=ind_base;

  for (l=0;l!=n;l++)                            /* VECDEL MMBETADEL */
    for (i=0;i!=m;i++) *pc++ *= beta;                  /* MMBETADEL */

  for (i=0;i!=mb;i++) {      /* For each block row i                   */
    jb = bpntrb[i];             /* beginning block in row              */
    je = bpntre[i];             /* ending block in row                 */
    cs = i*lb;                  /* start of c for this block operation */
    for (j=jb;j!=je;j++) {      /* For each block j in this row        */
      index = bindx[j];
      bs=(index-1)*lb;          /* start of b                       */
      pval = &val[(j-1)*mm];       /* pointer to start of block        */
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
      if ( index == i + ind_base ) {
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
}

void BSRskew_MatMult_CaABbC_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
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

  bindx-=ind_base;

  for (l=0;l!=n;l++)                            /* VECDEL MMBETADEL */
    for (i=0;i!=m;i++) *pc++ *= beta;                  /* MMBETADEL */

  for (i=0;i!=mb;i++) {      /* For each block row i                   */
    jb = bpntrb[i];             /* beginning block in row              */
    je = bpntre[i];             /* ending block in row                 */
    cs = i*lb;                  /* start of c for this block operation */
    for (j=jb;j!=je;j++) {      /* For each block j in this row        */
      index = bindx[j];
 /*   If on diagonal, assume zero.     Go to next.                     */
      if ( index == i + ind_base ) {
        continue;                  
      }
      bs=(index-1)*lb;             /* start of b                       */
      pval = &val[(j-1)*mm];       /* pointer to start of block        */
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
}

void BSRskew_MatMult_CaATBbC_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
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

  bindx-=ind_base;

  for (l=0;l!=n;l++)                            /* VECDEL MMBETADEL */
    for (i=0;i!=k;i++) *pc++ *= beta;                  /* MMBETADEL */

  for (i=0;i!=mb;i++) {      /* For each block row i                   */
    jb = bpntrb[i];             /* beginning block in row              */
    je = bpntre[i];             /* ending block in row                 */
    cs = i*lb;                  /* start of c for this block operation */
    for (j=jb;j!=je;j++) {      /* For each block j in this row        */
      index = bindx[j];
 /*   If on diagonal, done with block. Go to next.                     */
      if ( index == i + ind_base ) {
        continue;                  
      }
      bs=(index-1)*lb;             /* start of b                       */
      pval = &val[(j-1)*mm];       /* pointer to start of block        */
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
}

