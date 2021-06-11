#include "nsbvbr.h"

void VBR_MatMult_CaABbC_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int cs,mm,bs,br;
  int index;
  int l;                                               /* VECDEL */
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  int k=cpntr[kb]-cpntr[0];

  val-=ind_base;
  b-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  for (l=0;l!=n;l++)                            /* VECDEL MMBETADEL */
    for (i=0;i!=m;i++) *pc++ *= beta;                  /* MMBETADEL */

  c-=ind_base;
 
  for (i=0;i!=mb;i++) {      /* For each block row i                   */
    jb = bpntrb[i];             /* beginning block in row              */
    je = bpntre[i];             /* ending block in row                 */
    cs = rpntr[i];              /* start of c for this block operation */
    mm = rpntr[i+1]-rpntr[i];   /* point rows for this block operation */
    for (j=jb;j!=je;j++) {      /* For each block j in this row        */
      index = bindx[j];
      bs=cpntr[index];             /* start of b                       */
      br=cpntr[index+1]-bs;        /* rows in b (columns in block)     */
      pval = &val[indx[j]];        /* pointer to start of block        */
      ptmp = pval;                 /* copy for later use     VECDEL    */
                                   /* GEMM call inlined here:          */
      pb = &b[bs];                 /* pointer to start of b            */
      pc = &c[rpntr[i]];           /* pointer to start of c            */
      for (l=0;l!=n;l++) {         /* For each column of c   VECDEL    */             
       for (jj=0;jj!=br;jj++) {       /* For each column of block      */
         if( pb[jj] != 0.0 ) {           /* If non-zero multiplier     */
           t = alpha * pb[jj];            
           for (ii=0;ii!=mm;ii++) {      /* For each element in column */
             pc[ii] += t* (*pval++);        /* update element of c     */
           }
         } else {
           pval+=mm;                     /* Skip this column of block  */
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      /* to next column of c;VECDEL */
      }                                               /* VECDEL */
    }
  }
}

void VBR_MatMult_CaATBbC_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int cs,mm,bs,cr;
  int index;
  int l;                                               /* VECDEL */
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  int k=cpntr[kb]-cpntr[0];

  val-=ind_base;
  b-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  for (l=0;l!=n;l++)                            /* VECDEL MMBETADEL */
    for (i=0;i!=k;i++) *pc++ *= beta;                  /* MMBETADEL */

  c-=ind_base;
 
  for (i=0;i!=mb;i++) {      /* For each block column i                */
    jb = bpntrb[i];             /* beginning block in column           */
    je = bpntre[i];             /* ending block in column              */
    bs = rpntr[i];              /* start of b for this block operation */
    mm = rpntr[i+1]-rpntr[i];   /* point cols for this block operation */
    for (j=jb;j!=je;j++) {      /* For each block j in this column     */
      index = bindx[j];
      cs=cpntr[index];             /* start of b                       */
      cr=cpntr[index+1]-cs;        /* rows in c (rows in block)        */
      pval = &val[indx[j]];        /* pointer to start of block        */
      ptmp = pval;                 /* copy for later use    VECDEL     */ 
                                   /* GEMM call inlined here:          */
      pb = &b[bs];                 /* pointer to start of b            */
      pc = &c[cs];                 /* pointer to start of c            */
      for (l=0;l!=n;l++) {         /* For each column of c  VECDEL     */     
       for (jj=0;jj!=cr;jj++) {       /* For each row    of block      */
           t = 0;
           for (ii=0;ii!=mm;ii++) {      /* accummulate dot with b     */
             t += alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  /* update element of c     */
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      /* to next column of c;VECDEL */
      }                                               /* VECDEL */
    }
  }
}

void VBRsymm_MatMult_CaABbC_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int cs,css,cr,mm,bs,br;
  int index;
  int l;                                               /* VECDEL */
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  int k=cpntr[kb]-cpntr[0];

  val-=ind_base;
  b-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  for (l=0;l!=n;l++)                            /* VECDEL MMBETADEL */
    for (i=0;i!=m;i++) *pc++ *= beta;                  /* MMBETADEL */

  c-=ind_base;
 
  for (i=0;i!=mb;i++) {      /* For each block row i                   */
    jb = bpntrb[i];             /* beginning block in row              */
    je = bpntre[i];             /* ending block in row                 */
    cs = rpntr[i];              /* start of c for this block operation */
    mm = rpntr[i+1]-rpntr[i];   /* point rows for this block operation */
    for (j=jb;j!=je;j++) {      /* For each block j in this row        */
      index = bindx[j];
      bs=cpntr[index];             /* start of b                       */
      br=cpntr[index+1]-bs;        /* rows in b (columns in block)     */
      pval = &val[indx[j]];        /* pointer to start of block        */
      ptmp = pval;                 /* copy for later use               */
                                   /* GEMM call inlined here:          */
      pb = &b[bs];                 /* pointer to start of b            */
      pc = &c[rpntr[i]];           /* pointer to start of c            */
      for (l=0;l!=n;l++) {         /* For each column of c  VECDEL     */    
       for (jj=0;jj!=br;jj++) {       /* For each column of block      */
         if( pb[jj] != 0.0 ) {           /* If non-zero multiplier     */
           t = alpha * pb[jj];            
           for (ii=0;ii!=mm;ii++) {      /* For each element in column */
             pc[ii] += t* (*pval++);        /* update element of c     */
           }
         } else {
           pval+=mm;                     /* Skip this column of block  */
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      /* to next column of c;VECDEL */
      }                                               /* VECDEL */
 /*   If on diagonal, done with block. Go to next.                     */
      if ( index == i + ind_base ) {
        continue;                  
      }
 /*   Else, perform mult with transpose of this block                  */
      cr= br;                      /* rows in c (rows in b)            */
      pval = ptmp;                 /* pointer to start of block        */
                                   /* GEMM call inlined here:          */
      pb = &b[cs];                 /* pointer to start of b            */
      pc = &c[bs];                 /* pointer to start of c            */
      for (l=0;l!=n;l++) {         /* For each column of c  VECDEL     */
       for (jj=0;jj!=cr;jj++) {       /* For each row    of block      */
           t = 0;
           for (ii=0;ii!=mm;ii++) {      /* accummulate dot with b     */
             t += alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  /* update element of c        */
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      /* to next column of c;VECDEL */
      }                                               /* VECDEL */
    }
  }
}

void VBRskew_MatMult_CaABbC_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int cs,css,cr,mm,bs,br;
  int index;
  int l;                                               /* VECDEL */
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  int k=cpntr[kb]-cpntr[0];

  val-=ind_base;
  b-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  for (l=0;l!=n;l++)                            /* VECDEL MMBETADEL */
    for (i=0;i!=m;i++) *pc++ *= beta;                  /* MMBETADEL */

  c-=ind_base;
 
  for (i=0;i!=mb;i++) {      /* For each block row i                   */
    jb = bpntrb[i];             /* beginning block in row              */
    je = bpntre[i];             /* ending block in row                 */
    cs = rpntr[i];              /* start of c for this block operation */
    mm = rpntr[i+1]-rpntr[i];   /* point rows for this block operation */
    for (j=jb;j!=je;j++) {      /* For each block j in this row        */
      index = bindx[j];
 /*   If on diagonal, assume zero.     Go to next.                     */
      if ( index == i + ind_base ) {
        continue;                  
      }
      bs=cpntr[index];             /* start of b                       */
      br=cpntr[index+1]-bs;        /* rows in b (columns in block)     */
      pval = &val[indx[j]];        /* pointer to start of block        */
      ptmp = pval;                 /* copy for later use               */
                                   /* GEMM call inlined here:          */
      pb = &b[bs];                 /* pointer to start of b            */
      pc = &c[rpntr[i]];           /* pointer to start of c            */
      for (l=0;l!=n;l++) {         /* For each column of c  VECDEL     */ 
       for (jj=0;jj!=br;jj++) {       /* For each column of block      */
         if( pb[jj] != 0.0 ) {           /* If non-zero multiplier     */
           t = alpha * pb[jj];            
           for (ii=0;ii!=mm;ii++) {      /* For each element in column */
             pc[ii] += t* (*pval++);        /* update element of c     */
           }
         } else {
           pval+=mm;                     /* Skip this column of block  */
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      /* to next column of c;VECDEL */
      }                                               /* VECDEL */
 /*   Else, perform mult with transpose of this block                  */
      cr= br;                      /* rows in c (rows in b)            */
      pval = ptmp;                 /* pointer to start of block        */
                                   /* GEMM call inlined here:          */
      pb = &b[cs];                 /* pointer to start of b            */
      pc = &c[bs];                 /* pointer to start of c            */
      for (l=0;l!=n;l++) {         /* For each column of c  VECDEL     */
       for (jj=0;jj!=cr;jj++) {       /* For each row    of block      */
           t = 0;
           for (ii=0;ii!=mm;ii++) {      /* accummulate dot with b     */
             t -= alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  /* update element of c        */
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      /* to next column of c;VECDEL */
      }                                               /* VECDEL */
    }
  }
}

void VBRskew_MatMult_CaATBbC_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  int i,j,jb,je,block;
  int cs,css,cr,mm,bs,br;
  int index;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  int k=cpntr[kb]-cpntr[0];
  const double *ptmp;
  int l;                                               /* VECDEL */

  val-=ind_base;
  b-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  for (l=0;l!=n;l++)                            /* VECDEL MMBETADEL */
    for (i=0;i!=k;i++) *pc++ *= beta;                  /* MMBETADEL */

  c-=ind_base;
 
  for (i=0;i!=mb;i++) {      /* For each block row i                   */
    jb = bpntrb[i];             /* beginning block in row              */
    je = bpntre[i];             /* ending block in row                 */
    cs = rpntr[i];              /* start of c for this block operation */
    mm = rpntr[i+1]-rpntr[i];   /* point rows for this block operation */
    for (j=jb;j!=je;j++) {      /* For each block j in this row        */
      index = bindx[j];
 /*   If on diagonal, done with block. Go to next.                     */
      if ( index == i + ind_base ) {
        continue;                  
      }
      bs=cpntr[index];             /* start of b                       */
      br=cpntr[index+1]-bs;        /* rows in b (columns in block)     */
      pval = &val[indx[j]];        /* pointer to start of block        */
      ptmp = pval;                 /* copy for later use               */
                                   /* GEMM call inlined here:          */
      pb = &b[bs];                 /* pointer to start of b            */
      pc = &c[rpntr[i]];           /* pointer to start of c            */
      for (l=0;l!=n;l++) {         /* For each column of c  VECDEL     */    
       for (jj=0;jj!=br;jj++) {       /* For each column of block      */
         if( pb[jj] != 0.0 ) {           /* If non-zero multiplier     */
           t = alpha * pb[jj];            
           for (ii=0;ii!=mm;ii++) {      /* For each element in column */
             pc[ii] -= t* (*pval++);        /* update element of c     */
           }
         } else {
           pval+=mm;                     /* Skip this column of block  */
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      /* to next column of c;VECDEL */
      }                                               /* VECDEL */
 /*   Else, perform mult with transpose of this block                  */
      cr= br;                      /* rows in c (rows in b)            */
      pval = ptmp;                 /* pointer to start of block        */
                                   /* GEMM call inlined here:          */
      pb = &b[cs];                 /* pointer to start of b            */
      pc = &c[bs];                 /* pointer to start of c            */
      for (l=0;l!=n;l++) {         /* For each column of c  VECDEL     */
       for (jj=0;jj!=cr;jj++) {       /* For each row    of block      */
           t = 0;
           for (ii=0;ii!=mm;ii++) {      /* accummulate dot with b     */
             t += alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  /* update element of c        */
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      /* to next column of c;VECDEL */
      }                                               /* VECDEL */
    }
  }
}

