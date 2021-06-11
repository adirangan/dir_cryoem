#include "nsbbsr.h"

void BSR_MatTriangSlvLU_CaDADBbC_double(
                 const int mb, const int n, 
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             /* MTSDVLBETADEL */
  double *ptmpwork;                               /* MTSDVLBETADEL */
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               /* VECDEL */
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  /* First:                                                 MTSBETADEL */
  /*   Store beta*c in temp array                           MTSBETADEL */
  /*   freeing c for intermediate computations              MTSBETADEL */

  for (l=0;l!=n;l++)                              /* VECDEL MTSBETADEL */
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       /*  MTSBETADEL */

  ptmpwork = pwork;          /* End of the work area       MTSDVLBETADEL */
                             /* needed to store beta*c.    MTSDVLBETADEL */
                             /* (None needed if beta==0)   MTSDVLBETADEL */
                             /* Subsequent work must not   MTSDVLBETADEL */
                             /* overwrite initial portion  MTSDVLBETADEL */
                             /* of work array              MTSDVLBETADEL */
  /* Second:                                                           */
  /*   Calculate  C <- DVR*B                                           */
  for (i=0;i!=mb;i++) {      /* For each block row i                   */
    jb = bpntrb[i];             /* beginning block in row              */
    je = bpntre[i];             /* ending block in row                 */
    bs = i*lb;                  /* start of b and c                    */
    pb = &b[bs];                /* starting point for b                */
    pc = &c[bs];                /* starting point for c                */

    /* Store in c the block product dvr*b  */

    for (l=0;l!=n;l++){                                      /* VECDEL */
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    /* Rest of loop not needed if DVR = I    DVRDEL */
      pdr = &dvr[i*mm];                                   /* DVRDEL */
      for (j=0;j!=lb;j++) {                               /* DVRDEL */
        pc[j] *= (*pdr);                                  /* DVRDEL */
        pdr += lb+1;                                      /* DVRDEL */
      }                                                   /* DVRDEL */
      pdr = &dvr[i*mm];                                   /* DVRDEL */
      for (j=0;j!=lb;j++) {                               /* DVRDEL */
        for (ii=0;ii!=lb;ii++)                            /* DVRDEL */
          if ( ii == j ) pdr++;                           /* DVRDEL */
          else  pc[j] += (*pdr++) * pb[ii];               /* DVRDEL */
      }                                                   /* DVRDEL */
      pc+=ldc;pb+=ldb;                                       /* VECDEL */
    }                                                        /* VECDEL */
  }

  /* Now solve for unknown in A*c=c      */

  for (i=0;i!=mb;i++) {      /* For each block row i                   */
    jb = bpntrb[i];             /* beginning block in row              */
    je = bpntre[i];             /* ending block in row                 */
    bs = i*lb;                  /* start of b and c                    */

    for (j=jb;j!=je;j++) {      /* For each block j in this row        */
      index = bindx[j];             /* block column index              */
      if ( index == i + ind_base ){ /* on diagonal                     */
        continue;                   /* assume identity block, skip     */
      } else {
        pval = &val[(j-1)*mm];      /* pointer to start of block       */
        ptmp = pval;                                         /* VECDEL */
                                 /* GEMM call inlined here:            */
        cs = (index-1)*lb;          
        pc = &c[bs];                /* start of block of c to modify   */
        cj = &c[cs];                /* start of this block of c        */
        for (l=0;l!=n;l++) {        /* For each column of c     VECDEL */    
          for (ii=0;ii!=lb;ii++) {      /* For each column of A block  */
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        /* For row of work block   */
              pc[jj] -= (*pval++) * t;       /* Accumulate dot product */
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; /* Skip to next column     VECDEL */
        }                                                    /* VECDEL */
      }
    }
  }

  /* Third:                                                            */
  /* Modify intermediate solution by DVL (block) scaling:              */

  ds = 0;                                                    /* DVLDEL */
  for (i=0;i!=mb;i++) {      /* For each block row i            DVLDEL */
    bs = i*lb;                  /* start of b and c             DVLDEL */
    pdl = &dvl[ds];                                          /* DVLDEL */
    ds += mm;                   /* keep track of place in dvr   DVLDEL */
    pc = &c[bs];                                             /* DVLDEL */
    pwork=ptmpwork;                                          /* DVLDEL */
    for (l=0;l!=n;l++) {        /* For each column of c  VECDEL DVLDEL */    
      for (j=0;j!=lb;j++) {     /* For each row of DVL block    DVLDEL */
        t = 0;                                               /* DVLDEL */
        for (jj=0;jj!=lb;jj++)                               /* DVLDEL */
          t+= pc[jj]*pdl[j+jj*lb];  /* Accumulate dot product   DVLDEL */
        pwork[j] = t;               /* Store in work array      DVLDEL */
      }                                                      /* DVLDEL */
        /* copy work array into c */
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                  /* DVLDEL */ 
      pc+=ldc;                      /* Skip column       VECDEL DVLDEL */
    }                                                 /* VECDEL DVLDEL */
  }                                                          /* DVLDEL */

  /* Finally:                                BLKPOSTPROCESS */
  /* Scale by alpha and          ALPHA1DEL   BLKPOSTPROCESS */
  /*   add BETA * c (from work)  MTSBETADEL  BLKPOSTPROCESS */

  pc=c;                                   /* BLKPOSTPROCESS */
  pwork=work;                             /* BLKPOSTPROCESS    MTSBETADEL */
  for (l=0;l!=n;l++)                      /* BLKPOSTPROCESS        VECDEL */
    for (i=0;i!=m;i++) {                  /* BLKPOSTPROCESS */
          *pc *= alpha;                   /* BLKPOSTPROCESS    ALPHA1DEL  */
          *pc += *pwork;                  /* BLKPOSTPROCESS    MTSBETADEL */
          pc++;                           /* BLKPOSTPROCESS */
          pwork++;                        /* BLKPOSTPROCESS    MTSBETADEL */
    }                                     /* BLKPOSTPROCESS */

}


void BSR_MatTriangSlvUU_CaDADBbC_double(
                 const int mb, const int n,
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                            /* MTSDVLBETADEL */
  double *ptmpwork;                              /* MTSDVLBETADEL */
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               /* VECDEL */
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  /* First:                                                 MTSBETADEL */
  /*   Store beta*c in temp array                           MTSBETADEL */
  /*   freeing c for intermediate computations              MTSBETADEL */

  for (l=0;l!=n;l++)                              /* VECDEL MTSBETADEL */
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       /*  MTSBETADEL */

  ptmpwork = pwork;          /* End of the work area       MTSDVLBETADEL */
                             /* needed to store beta*c.    MTSDVLBETADEL */
                             /* (None needed if beta==0)   MTSDVLBETADEL */
                             /* Subsequent work must not   MTSDVLBETADEL */
                             /* overwrite initial portion  MTSDVLBETADEL */
                             /* of work array              MTSDVLBETADEL */
  /* Second:                                                           */
  /*   Calculate  C <- DVR*B                                           */
  for (i=0;i!=mb;i++) {      /* For each block row i                   */
    jb = bpntrb[i];             /* beginning block in row              */
    je = bpntre[i];             /* ending block in row                 */
    bs = i*lb;                  /* start of b and c                    */
    pb = &b[bs];                /* starting point for b                */
    pc = &c[bs];                /* starting point for c                */

    /* Store in c the block product dvr*b  */

    for (l=0;l!=n;l++){                                      /* VECDEL */
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    /* Rest of loop not needed if DVR = I    DVRDEL */
      pdr = &dvr[i*mm];                                    /* DVRDEL */
      for (j=0;j!=lb;j++) {                               /* DVRDEL */
        pc[j] *= (*pdr);                                  /* DVRDEL */
        pdr += lb+1;                                      /* DVRDEL */
      }                                                   /* DVRDEL */
      pdr = &dvr[i*mm];                                   /* DVRDEL */
      for (j=0;j!=lb;j++) {                               /* DVRDEL */
        for (ii=0;ii!=lb;ii++)                            /* DVRDEL */
          if ( ii == j ) pdr++;                           /* DVRDEL */
          else  pc[j] += (*pdr++) * pb[ii];               /* DVRDEL */
      }                                                   /* DVRDEL */
      pc+=ldc;pb+=ldb;                                       /* VECDEL */
    }                                                        /* VECDEL */
  }

  /* Now solve for unknown in A*c=c      */

  for (i=mb-1;i!=-1;i--) {      /* For each block row i                */
    jb = bpntrb[i];             /* beginning block in row              */
    je = bpntre[i];             /* ending block in row                 */
    bs = lb*i;                  /* start of b and c                    */

    for (j=jb;j!=je;j++) {      /* For each block j in this row        */
      index = bindx[j];             /* block column index              */
      if ( index == i + ind_base ){ /* on diagonal                     */
        continue;                   /* assume identity block, skip     */
      } else {
        pval = &val[(j-1)*mm];      /* pointer to start of block       */
        ptmp = pval;                                         /* VECDEL */
                                 /* GEMM call inlined here:            */
        cs = (index-1)*lb;          
        pc = &c[bs];                /* start of block of c to modify   */
        cj = &c[cs];                /* start of this block of c        */
        for (l=0;l!=n;l++) {        /* For each column of c     VECDEL */    
          for (ii=0;ii!=lb;ii++) {      /* For each column of A block  */
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        /* For row of work block   */
              pc[jj] -= (*pval++) * t;       /* Accumulate dot product */
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; /* Skip to next column     VECDEL */
        }                                                    /* VECDEL */
      }
    }
  }

  /* Third:                                                            */
  /* Modify intermediate solution by DVL (block) scaling:              */

  ds = 0;                                                    /* DVLDEL */
  for (i=0;i!=mb;i++) {      /* For each block row i            DVLDEL */
    bs = i*lb;                  /* start of b and c             DVLDEL */
    pdl = &dvl[ds];                                          /* DVLDEL */
    ds += mm;                   /* keep track of place in dvr   DVLDEL */
    pc = &c[bs];                                             /* DVLDEL */
    pwork=ptmpwork;                                          /* DVLDEL */
    for (l=0;l!=n;l++) {        /* For each column of c  VECDEL DVLDEL */    
      for (j=0;j!=lb;j++) {     /* For each row of DVL block    DVLDEL */
        t = 0;                                               /* DVLDEL */
        for (jj=0;jj!=lb;jj++)                               /* DVLDEL */
          t+= pc[jj]*pdl[j+jj*lb];  /* Accumulate dot product   DVLDEL */
        pwork[j] = t;               /* Store in work array      DVLDEL */
      }                                                      /* DVLDEL */
        /* copy work array into c */
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                  /* DVLDEL */ 
      pc+=ldc;                      /* Skip column       VECDEL DVLDEL */
    }                                                 /* VECDEL DVLDEL */
  }                                                          /* DVLDEL */

  /* Finally:                                BLKPOSTPROCESS */
  /* Scale by alpha                          BLKPOSTPROCESS */
  /*    and add beta * c (from work)         BLKPOSTPROCESS */

  pc=c;                                   /* BLKPOSTPROCESS */
  pwork=work;                             /* BLKPOSTPROCESS    MTSBETADEL */
  for (l=0;l!=n;l++)                      /* BLKPOSTPROCESS        VECDEL */
    for (i=0;i!=m;i++) {                  /* BLKPOSTPROCESS */
          *pc *= alpha;                   /* BLKPOSTPROCESS    ALPHA1DEL  */
          *pc += *pwork;                  /* BLKPOSTPROCESS    MTSBETADEL */
          pc++;                           /* BLKPOSTPROCESS */
          pwork++;                        /* BLKPOSTPROCESS    MTSBETADEL */
    }                                     /* BLKPOSTPROCESS */


}

void BSR_MatTriangSlvLU_CaDATDBbC_double(
                 const int mb, const int n,
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             /* MTSDVLBETADEL */
  double *ptmpwork;                               /* MTSDVLBETADEL */
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               /* VECDEL */
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  /* First:                                                 MTSBETADEL */
  /*   Store beta*c in temp array                           MTSBETADEL */
  /*   freeing c for intermediate computations              MTSBETADEL */

  for (l=0;l!=n;l++)                              /* VECDEL MTSBETADEL */
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       /*  MTSBETADEL */

  ptmpwork = pwork;        /* End of the work area       MTSDVLBETADEL */
                           /* needed to store beta*c.    MTSDVLBETADEL */
                           /* (None needed if beta==0)   MTSDVLBETADEL */
                           /* Subsequent work must not   MTSDVLBETADEL */
                           /* overwrite initial portion  MTSDVLBETADEL */
                           /* of work array              MTSDVLBETADEL */
  /* Second:                                                           */
  /*   Calculate  C <- DVR*B                                           */
  for (i=0;i!=mb;i++) {      /* For each block row i                   */
    jb = bpntrb[i];             /* beginning block in row              */
    je = bpntre[i];             /* ending block in row                 */
    bs = i*lb;                  /* start of b and c                    */
    pb = &b[bs];                /* starting point for b                */
    pc = &c[bs];                /* starting point for c                */

    /* Store in c the block product dvr*b  */

    for (l=0;l!=n;l++){                                      /* VECDEL */
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    /* Rest of loop not needed if DVR = I    DVRDEL */
      pdr = &dvr[i*mm];                                   /* DVRDEL */
      for (j=0;j!=lb;j++) {                               /* DVRDEL */
        pc[j] *= (*pdr);                                  /* DVRDEL */
        pdr += lb+1;                                      /* DVRDEL */
      }                                                   /* DVRDEL */
      pdr = &dvr[i*mm];                                   /* DVRDEL */
      for (j=0;j!=lb;j++) {                               /* DVRDEL */
        for (ii=0;ii!=lb;ii++)                            /* DVRDEL */
          if ( ii == j ) pdr++;                           /* DVRDEL */
          else  pc[j] += (*pdr++) * pb[ii];               /* DVRDEL */
      }                                                   /* DVRDEL */
      pc+=ldc;pb+=ldb;                                       /* VECDEL */
    }                                                        /* VECDEL */
  }

  /* Now solve for unknown in A^T*c=c      */

  for (i=mb-1;i!=-1;i--) {      /* For each block row i                */
    jb = bpntrb[i];             /* beginning block in row              */
    je = bpntre[i];             /* ending block in row                 */
    bs = i*lb;                  /* start of b and c                    */

    for (j=je-1;j!=jb-1;j--) {  /* For each block j in this row        */
      index = bindx[j];             /* block column index              */
      if ( index == i + ind_base ){ /* on diagonal                     */
        continue;                   /* assume identity block, skip     */
      } else {
        pval = &val[(j-1)*mm];      /* pointer to start of block       */
        ptmp = pval;                                         /* VECDEL */
                                 /* GEMM call inlined here:            */
        cs = (index-1)*lb;
        pc = &c[cs];                /* start of block of c to modify   */
        cj = &c[bs];                /* start of this block of c        */
        for (l=0;l!=n;l++) {        /* For each column of c     VECDEL */    
          for (ii=0;ii!=lb;ii++) {      /* For each row of A' block    */
            t = 0;
            for (jj=0;jj!=lb;jj++) {        /* For row of work block   */
              t += (*pval++) * cj[jj];      /* Accumulate dot product  */
            }
            pc[ii] -= t;       
          }
          pc+=ldc;cj+=ldc;pval=ptmp; /* Skip to next column     VECDEL */
        }                                                    /* VECDEL */
      }
    }
  }

  /* Third:                                                            */
  /* Modify intermediate solution by DVL (block) scaling:              */

  ds = 0;                                                    /* DVLDEL */
  for (i=0;i!=mb;i++) {      /* For each block row i            DVLDEL */
    bs = i*lb;                  /* start of b and c             DVLDEL */
    pdl = &dvl[ds];                                          /* DVLDEL */
    ds += mm;                   /* keep track of place in dvr   DVLDEL */
    pc = &c[bs];                                             /* DVLDEL */
    pwork=ptmpwork;                                          /* DVLDEL */
    for (l=0;l!=n;l++) {        /* For each column of c  VECDEL DVLDEL */    
      for (j=0;j!=lb;j++) {     /* For each row of DVL block    DVLDEL */
        t = 0;                                               /* DVLDEL */
        for (jj=0;jj!=lb;jj++)                               /* DVLDEL */
          t+= pc[jj]*pdl[j+jj*lb];  /* Accumulate dot product   DVLDEL */
        pwork[j] = t;               /* Store in work array      DVLDEL */
      }                                                      /* DVLDEL */
        /* copy work array into c */
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                  /* DVLDEL */ 
      pc+=ldc;                      /* Skip column       VECDEL DVLDEL */
    }                                                 /* VECDEL DVLDEL */
  }                                                          /* DVLDEL */

  /* Finally:                                BLKPOSTPROCESS */
  /* Scale by alpha                          BLKPOSTPROCESS */
  /*    and add beta * c (from work)         BLKPOSTPROCESS */

  pc=c;                                   /* BLKPOSTPROCESS */
  pwork=work;                             /* BLKPOSTPROCESS    MTSBETADEL */
  for (l=0;l!=n;l++)                      /* BLKPOSTPROCESS        VECDEL */
    for (i=0;i!=m;i++) {                  /* BLKPOSTPROCESS */
          *pc *= alpha;                   /* BLKPOSTPROCESS    ALPHA1DEL  */
          *pc += *pwork;                  /* BLKPOSTPROCESS    MTSBETADEL */
          pc++;                           /* BLKPOSTPROCESS */
          pwork++;                        /* BLKPOSTPROCESS    MTSBETADEL */
    }                                     /* BLKPOSTPROCESS */


}

void BSR_MatTriangSlvUU_CaDATDBbC_double(
                 const int mb, const int n, 
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             /* MTSDVLBETADEL */
  double *ptmpwork;                               /* MTSDVLBETADEL */
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               /* VECDEL */
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  /* First:                                                 MTSBETADEL */
  /*   Store beta*c in temp array                           MTSBETADEL */
  /*   freeing c for intermediate computations              MTSBETADEL */

  for (l=0;l!=n;l++)                              /* VECDEL MTSBETADEL */
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       /*  MTSBETADEL */

  ptmpwork = pwork;          /* End of the work area       MTSDVLBETADEL */
                             /* needed to store beta*c.    MTSDVLBETADEL */
                             /* (None needed if beta==0)   MTSDVLBETADEL */
                             /* Subsequent work must not   MTSDVLBETADEL */
                             /* overwrite initial portion  MTSDVLBETADEL */
                             /* of work array              MTSDVLBETADEL */
  /* Second:                                                           */
  /*   Calculate  C <- DVR*B                                           */
  for (i=0;i!=mb;i++) {      /* For each block row i                   */
    jb = bpntrb[i];             /* beginning block in row              */
    je = bpntre[i];             /* ending block in row                 */
    bs = i*lb;                  /* start of b and c                    */
    pb = &b[bs];                /* starting point for b                */
    pc = &c[bs];                /* starting point for c                */

    /* Store in c the block product dvr*b  */

    for (l=0;l!=n;l++){                                      /* VECDEL */
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    /* Rest of loop not needed if DVR = I    DVRDEL */
      pdr = &dvr[i*mm];                                   /* DVRDEL */
      for (j=0;j!=lb;j++) {                               /* DVRDEL */
        pc[j] *= (*pdr);                                  /* DVRDEL */
        pdr += lb+1;                                      /* DVRDEL */
      }                                                   /* DVRDEL */
      pdr = &dvr[i*mm];                                   /* DVRDEL */
      for (j=0;j!=lb;j++) {                               /* DVRDEL */
        for (ii=0;ii!=lb;ii++)                            /* DVRDEL */
          if ( ii == j ) pdr++;                           /* DVRDEL */
          else  pc[j] += (*pdr++) * pb[ii];               /* DVRDEL */
      }                                                   /* DVRDEL */
      pc+=ldc;pb+=ldb;                                       /* VECDEL */
    }                                                        /* VECDEL */
  }

  /* Now solve for unknown in A^T*c=c      */

  for (i=0;i!=mb;i++) {      /* For each block row i                   */
    jb = bpntrb[i];             /* beginning block in row              */
    je = bpntre[i];             /* ending block in row                 */
    bs = i*lb;                  /* start of b and c                    */

    for (j=jb;j!=je;j++) {      /* For each block j in this row        */
      index = bindx[j];             /* block column index              */
      if ( index == i + ind_base ){ /* on diagonal                     */
        continue;                   /* assume identity block, skip     */
      } else {
        pval = &val[(j-1)*mm];      /* pointer to start of block       */
        ptmp = pval;                                         /* VECDEL */
                                 /* GEMM call inlined here:            */
        cs = (index-1)*lb;          
        pc = &c[cs];                /* start of block of c to modify   */
        cj = &c[bs];                /* start of this block of c        */
        for (l=0;l!=n;l++) {        /* For each column of c     VECDEL */    
          for (ii=0;ii!=lb;ii++) {      /* For each column of A block  */
            t = 0;
            for (jj=0;jj!=lb;jj++) {        /* For row of work block   */
              t += (*pval++) * cj[jj];      /* Accumulate dot product  */
            }
            pc[ii] -= t;
          }
          pc+=ldc;cj+=ldc;pval=ptmp; /* Skip to next column     VECDEL */
        }                                                    /* VECDEL */
      }
    }
  }

  /* Third:                                                            */
  /* Modify intermediate solution by DVL (block) scaling:              */

  ds = 0;                                                    /* DVLDEL */
  for (i=0;i!=mb;i++) {      /* For each block row i            DVLDEL */
    bs = i*lb;                  /* start of b and c             DVLDEL */
    pdl = &dvl[ds];                                          /* DVLDEL */
    ds += mm;                /* keep track of place in dvr   DVLDEL */
    pc = &c[bs];                                             /* DVLDEL */
    pwork=ptmpwork;                                          /* DVLDEL */
    for (l=0;l!=n;l++) {        /* For each column of c  VECDEL DVLDEL */    
      for (j=0;j!=lb;j++) {     /* For each row of DVL block    DVLDEL */
        t = 0;                                               /* DVLDEL */
        for (jj=0;jj!=lb;jj++)                               /* DVLDEL */
          t+= pc[jj]*pdl[j+jj*lb];  /* Accumulate dot product   DVLDEL */
        pwork[j] = t;               /* Store in work array      DVLDEL */
      }                                                      /* DVLDEL */
        /* copy work array into c */
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                  /* DVLDEL */ 
      pc+=ldc;                      /* Skip column       VECDEL DVLDEL */
    }                                                 /* VECDEL DVLDEL */
  }                                                          /* DVLDEL */

  /* Finally:                                BLKPOSTPROCESS */
  /* Scale by alpha                          BLKPOSTPROCESS */
  /*    and add beta * c (from work)         BLKPOSTPROCESS */

  pc=c;                                   /* BLKPOSTPROCESS */
  pwork=work;                             /* BLKPOSTPROCESS    MTSBETADEL */
  for (l=0;l!=n;l++)                      /* BLKPOSTPROCESS        VECDEL */
    for (i=0;i!=m;i++) {                  /* BLKPOSTPROCESS */
          *pc *= alpha;                   /* BLKPOSTPROCESS    ALPHA1DEL  */
          *pc += *pwork;                  /* BLKPOSTPROCESS    MTSBETADEL */
          pc++;                           /* BLKPOSTPROCESS */
          pwork++;                        /* BLKPOSTPROCESS    MTSBETADEL */
    }                                     /* BLKPOSTPROCESS */


}
