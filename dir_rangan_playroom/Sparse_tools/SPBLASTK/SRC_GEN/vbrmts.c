#include "nsbvbr.h"

void VBR_MatTriangSlvLU_CaDADBbC_double(
                 const int mb, const int n, 
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int l;                                               /* VECDEL */
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

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
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      /* For each block row i                   */
    jb = bpntrb[i];             /* beginning block in row              */
    je = bpntre[i];             /* ending block in row                 */
    bs = rpntr[i];              /* start of b and c                    */
    ds += mm*mm;                /* keep track of place in dvr          */
    mm = rpntr[i+1]-rpntr[i];   /* point rows for this block operation */
    pb = &b[bs];                /* starting point for b                */
    pc = &c[bs];                /* starting point for c                */

    /* Store in c the block product dvr*b  */

    for (l=0;l!=n;l++){                                      /* VECDEL */
      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    /* Rest of loop not needed if DVR = I    DVRDEL */
      pdr = &dvr[ds];                                     /* DVRDEL */
      for (j=0;j!=mm;j++) {                               /* DVRDEL */
        pc[j] *= (*pdr);                                  /* DVRDEL */
        pdr += mm+1;                                      /* DVRDEL */
      }                                                   /* DVRDEL */
      pdr = &dvr[ds];                                     /* DVRDEL */
      for (j=0;j!=mm;j++) {                               /* DVRDEL */
        for (ii=0;ii!=mm;ii++)                            /* DVRDEL */
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
    bs = rpntr[i];              /* start of b and c                    */
    mm = rpntr[i+1]-rpntr[i];   /* point rows for this block operation */

    for (j=jb;j!=je;j++) {      /* For each block j in this row        */
      index = bindx[j];             /* block column index              */
      if ( index == i + ind_base ){ /* on diagonal                     */
        continue;                   /* assume identity block, skip     */
      } else {
        pval = &val[indx[j]];       /* pointer to start of block       */
        ptmp = pval;                                         /* VECDEL */
                                 /* GEMM call inlined here:            */
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   /* number of columns in block      */
        pb = &b[bs];                /* pointer to start of b           */
        pc = &c[bs];                /* start of block of c to modify   */
        cj = &c[cs];                /* start of this block of c        */
        for (l=0;l!=n;l++) {        /* For each column of c     VECDEL */    
          for (ii=0;ii!=cr;ii++) {      /* For each column of A block  */
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        /* For row of work block   */
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
  mm = 0;                                                    /* DVLDEL */
  for (i=0;i!=mb;i++) {      /* For each block row i            DVLDEL */
    bs = rpntr[i];              /* start of b and c             DVLDEL */
    ds += mm*mm;                /* keep track of place in dvr   DVLDEL */
    mm = rpntr[i+1]-rpntr[i];   /* point rows for this block    DVLDEL */
    pdl = &dvl[ds];                                          /* DVLDEL */
    pc = &c[bs];                                             /* DVLDEL */
    pwork=ptmpwork;                                          /* DVLDEL */
    for (l=0;l!=n;l++) {        /* For each column of c  VECDEL DVLDEL */    
      for (j=0;j!=mm;j++) {     /* For each row of DVL block    DVLDEL */
        t = 0;                                               /* DVLDEL */
        for (jj=0;jj!=mm;jj++)                               /* DVLDEL */
          t+= pc[jj]*pdl[j+jj*mm];  /* Accumulate dot product   DVLDEL */
        pwork[j] = t;               /* Store in work array      DVLDEL */
      }                                                      /* DVLDEL */
        /* copy work array into c */
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                  /* DVLDEL */ 
      pc+=ldc;                      /* Skip column       VECDEL DVLDEL */
    }                                                 /* VECDEL DVLDEL */
  }                                                          /* DVLDEL */

  /* Finally:                                BLKPOSTPROCESS */
  /* Scale by alpha and          ALPHA1DEL   BLKPOSTPROCESS */
  /*   add BETA * c (from work)  MTSBETADEL  BLKPOSTPROCESS */

  pc=c+ind_base;                          /* BLKPOSTPROCESS */
  pwork=work;                             /* BLKPOSTPROCESS    MTSBETADEL */
  for (l=0;l!=n;l++)                      /* BLKPOSTPROCESS        VECDEL */
    for (i=0;i!=m;i++) {                  /* BLKPOSTPROCESS */
          *pc *= alpha;                   /* BLKPOSTPROCESS    ALPHA1DEL  */
          *pc += *pwork;                  /* BLKPOSTPROCESS    MTSBETADEL */
          pc++;                           /* BLKPOSTPROCESS */
          pwork++;                        /* BLKPOSTPROCESS    MTSBETADEL */
    }                                     /* BLKPOSTPROCESS */

}


void VBR_MatTriangSlvUU_CaDADBbC_double(
                 const int mb, const int n,
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int l;                                               /* VECDEL */
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

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
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      /* For each block row i                   */
    jb = bpntrb[i];             /* beginning block in row              */
    je = bpntre[i];             /* ending block in row                 */
    bs = rpntr[i];              /* start of b and c                    */
    ds += mm*mm;                /* keep track of place in dvr          */
    mm = rpntr[i+1]-rpntr[i];   /* point rows for this block operation */
    pb = &b[bs];                /* starting point for b                */
    pc = &c[bs];                /* starting point for c                */

    /* Store in c the block product dvr*b  */

    for (l=0;l!=n;l++){                                      /* VECDEL */
      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    /* Rest of loop not needed if DVR = I    DVRDEL */
      pdr = &dvr[ds];                                     /* DVRDEL */
      for (j=0;j!=mm;j++) {                               /* DVRDEL */
        pc[j] *= (*pdr);                                  /* DVRDEL */
        pdr += mm+1;                                      /* DVRDEL */
      }                                                   /* DVRDEL */
      pdr = &dvr[ds];                                     /* DVRDEL */
      for (j=0;j!=mm;j++) {                               /* DVRDEL */
        for (ii=0;ii!=mm;ii++)                            /* DVRDEL */
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
    bs = rpntr[i];              /* start of b and c                    */
    mm = rpntr[i+1]-rpntr[i];   /* point rows for this block operation */
    pb = &b[bs];                /* starting point for b                */
    pc = &c[bs];                /* starting point for c                */

    for (j=jb;j!=je;j++) {      /* For each block j in this row        */
      index = bindx[j];             /* block column index              */
      if ( index == i + ind_base ){ /* on diagonal                     */
        continue;                   /* assume identity block, skip     */
      } else {
        pval = &val[indx[j]];       /* pointer to start of block       */
        ptmp = pval;                                         /* VECDEL */
                                 /* GEMM call inlined here:            */
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   /* number of columns in block      */
        pb = &b[bs];                /* pointer to start of b           */
        pc = &c[bs];                /* start of block of c to modify   */
        cj = &c[cs];                /* start of this block of c        */
        for (l=0;l!=n;l++) {        /* For each column of c     VECDEL */    
          for (ii=0;ii!=cr;ii++) {      /* For each column of A block  */
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        /* For row of work block   */
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
  mm = 0;                                                    /* DVLDEL */
  for (i=0;i!=mb;i++) {      /* For each block row i            DVLDEL */
    bs = rpntr[i];              /* start of b and c             DVLDEL */
    ds += mm*mm;                /* keep track of place in dvr   DVLDEL */
    mm = rpntr[i+1]-rpntr[i];   /* point rows for this block    DVLDEL */
    pdl = &dvl[ds];                                          /* DVLDEL */
    pc = &c[bs];                                             /* DVLDEL */
    pwork=ptmpwork;                                          /* DVLDEL */
    for (l=0;l!=n;l++) {        /* For each column of c  VECDEL DVLDEL */    
      for (j=0;j!=mm;j++) {     /* For each row of DVL block    DVLDEL */
        t = 0;                                               /* DVLDEL */
        for (jj=0;jj!=mm;jj++)                               /* DVLDEL */
          t+= pc[jj]*pdl[j+jj*mm];  /* Accumulate dot product   DVLDEL */
        pwork[j] = t;               /* Store in work array      DVLDEL */
      }                                                      /* DVLDEL */
        /* copy work array into c */
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                  /* DVLDEL */ 
      pc+=ldc;                      /* Skip column       VECDEL DVLDEL */
    }                                                 /* VECDEL DVLDEL */
  }                                                          /* DVLDEL */

  /* Finally:                                BLKPOSTPROCESS */
  /* Scale by alpha                          BLKPOSTPROCESS */
  /*    and add beta * c (from work)         BLKPOSTPROCESS */

  pc=c+ind_base;                          /* BLKPOSTPROCESS */
  pwork=work;                             /* BLKPOSTPROCESS    MTSBETADEL */
  for (l=0;l!=n;l++)                      /* BLKPOSTPROCESS        VECDEL */
    for (i=0;i!=m;i++) {                  /* BLKPOSTPROCESS */
          *pc *= alpha;                   /* BLKPOSTPROCESS    ALPHA1DEL  */
          *pc += *pwork;                  /* BLKPOSTPROCESS    MTSBETADEL */
          pc++;                           /* BLKPOSTPROCESS */
          pwork++;                        /* BLKPOSTPROCESS    MTSBETADEL */
    }                                     /* BLKPOSTPROCESS */


}

void VBR_MatTriangSlvLU_CaDATDBbC_double(
                 const int mb, const int n,
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int l;                                               /* VECDEL */
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

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
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      /* For each block row i                   */
    jb = bpntrb[i];             /* beginning block in row              */
    je = bpntre[i];             /* ending block in row                 */
    bs = rpntr[i];              /* start of b and c                    */
    ds += mm*mm;                /* keep track of place in dvr          */
    mm = rpntr[i+1]-rpntr[i];   /* point rows for this block operation */
    pb = &b[bs];                /* starting point for b                */
    pc = &c[bs];                /* starting point for c                */

    /* Store in c the block product dvr*b  */

    for (l=0;l!=n;l++){                                      /* VECDEL */
      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    /* Rest of loop not needed if DVR = I    DVRDEL */
      pdr = &dvr[ds];                                     /* DVRDEL */
      for (j=0;j!=mm;j++) {                               /* DVRDEL */
        pc[j] *= (*pdr);                                  /* DVRDEL */
        pdr += mm+1;                                      /* DVRDEL */
      }                                                   /* DVRDEL */
      pdr = &dvr[ds];                                     /* DVRDEL */
      for (j=0;j!=mm;j++) {                               /* DVRDEL */
        for (ii=0;ii!=mm;ii++)                            /* DVRDEL */
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
    bs = rpntr[i];              /* start of b and c                    */
    mm = rpntr[i+1]-rpntr[i];   /* point rows for this block operation */
    pb = &b[bs];                /* starting point for b                */
    pc = &c[bs];                /* starting point for c                */

    for (j=je-1;j!=jb-1;j--) {  /* For each block j in this row        */
      index = bindx[j];             /* block column index              */
      if ( index == i + ind_base ){ /* on diagonal                     */
        continue;                   /* assume identity block, skip     */
      } else {
        pval = &val[indx[j]];       /* pointer to start of block       */
        ptmp = pval;                                         /* VECDEL */
                                 /* GEMM call inlined here:            */
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   /* number of columns in block      */
        pc = &c[cs];                /* start of block of c to modify   */
        cj = &c[bs];                /* start of this block of c        */
        for (l=0;l!=n;l++) {        /* For each column of c     VECDEL */    
          for (ii=0;ii!=cr;ii++) {      /* For each row of A' block    */
            t = 0;
            for (jj=0;jj!=mm;jj++) {        /* For row of work block   */
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
  mm = 0;                                                    /* DVLDEL */
  for (i=0;i!=mb;i++) {      /* For each block row i            DVLDEL */
    bs = rpntr[i];              /* start of b and c             DVLDEL */
    ds += mm*mm;                /* keep track of place in dvr   DVLDEL */
    mm = rpntr[i+1]-rpntr[i];   /* point rows for this block    DVLDEL */
    pdl = &dvl[ds];                                          /* DVLDEL */
    pc = &c[bs];                                             /* DVLDEL */
    pwork=ptmpwork;                                          /* DVLDEL */
    for (l=0;l!=n;l++) {        /* For each column of c  VECDEL DVLDEL */    
      for (j=0;j!=mm;j++) {     /* For each row of DVL block    DVLDEL */
        t = 0;                                               /* DVLDEL */
        for (jj=0;jj!=mm;jj++)                               /* DVLDEL */
          t+= pc[jj]*pdl[j+jj*mm];  /* Accumulate dot product   DVLDEL */
        pwork[j] = t;               /* Store in work array      DVLDEL */
      }                                                      /* DVLDEL */
        /* copy work array into c */
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                  /* DVLDEL */ 
      pc+=ldc;                      /* Skip column       VECDEL DVLDEL */
    }                                                 /* VECDEL DVLDEL */
  }                                                          /* DVLDEL */

  /* Finally:                                BLKPOSTPROCESS */
  /* Scale by alpha                          BLKPOSTPROCESS */
  /*    and add beta * c (from work)         BLKPOSTPROCESS */

  pc=c+ind_base;                          /* BLKPOSTPROCESS */
  pwork=work;                             /* BLKPOSTPROCESS    MTSBETADEL */
  for (l=0;l!=n;l++)                      /* BLKPOSTPROCESS        VECDEL */
    for (i=0;i!=m;i++) {                  /* BLKPOSTPROCESS */
          *pc *= alpha;                   /* BLKPOSTPROCESS    ALPHA1DEL  */
          *pc += *pwork;                  /* BLKPOSTPROCESS    MTSBETADEL */
          pc++;                           /* BLKPOSTPROCESS */
          pwork++;                        /* BLKPOSTPROCESS    MTSBETADEL */
    }                                     /* BLKPOSTPROCESS */


}

void VBR_MatTriangSlvUU_CaDATDBbC_double(
                 const int mb, const int n, 
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int l;                                               /* VECDEL */
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

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
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      /* For each block row i                   */
    jb = bpntrb[i];             /* beginning block in row              */
    je = bpntre[i];             /* ending block in row                 */
    bs = rpntr[i];              /* start of b and c                    */
    ds += mm*mm;                /* keep track of place in dvr          */
    mm = rpntr[i+1]-rpntr[i];   /* point rows for this block operation */
    pb = &b[bs];                /* starting point for b                */
    pc = &c[bs];                /* starting point for c                */

    /* Store in c the block product dvr*b  */

    for (l=0;l!=n;l++){                                      /* VECDEL */
      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    /* Rest of loop not needed if DVR = I    DVRDEL */
      pdr = &dvr[ds];                                     /* DVRDEL */
      for (j=0;j!=mm;j++) {                               /* DVRDEL */
        pc[j] *= (*pdr);                                  /* DVRDEL */
        pdr += mm+1;                                      /* DVRDEL */
      }                                                   /* DVRDEL */
      pdr = &dvr[ds];                                     /* DVRDEL */
      for (j=0;j!=mm;j++) {                               /* DVRDEL */
        for (ii=0;ii!=mm;ii++)                            /* DVRDEL */
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
    bs = rpntr[i];              /* start of b and c                    */
    mm = rpntr[i+1]-rpntr[i];   /* point rows for this block operation */

    for (j=jb;j!=je;j++) {      /* For each block j in this row        */
      index = bindx[j];             /* block column index              */
      if ( index == i + ind_base ){ /* on diagonal                     */
        continue;                   /* assume identity block, skip     */
      } else {
        pval = &val[indx[j]];       /* pointer to start of block       */
        ptmp = pval;                                         /* VECDEL */
                                 /* GEMM call inlined here:            */
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   /* number of columns in block      */
        pc = &c[cs];                /* start of block of c to modify   */
        cj = &c[bs];                /* start of this block of c        */
        for (l=0;l!=n;l++) {        /* For each column of c     VECDEL */    
          for (ii=0;ii!=cr;ii++) {      /* For each column of A block  */
            t = 0;
            for (jj=0;jj!=mm;jj++) {        /* For row of work block   */
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
  mm = 0;                                                    /* DVLDEL */
  for (i=0;i!=mb;i++) {      /* For each block row i            DVLDEL */
    bs = rpntr[i];              /* start of b and c             DVLDEL */
    ds += mm*mm;                /* keep track of place in dvr   DVLDEL */
    mm = rpntr[i+1]-rpntr[i];   /* point rows for this block    DVLDEL */
    pdl = &dvl[ds];                                          /* DVLDEL */
    pc = &c[bs];                                             /* DVLDEL */
    pwork=ptmpwork;                                          /* DVLDEL */
    for (l=0;l!=n;l++) {        /* For each column of c  VECDEL DVLDEL */    
      for (j=0;j!=mm;j++) {     /* For each row of DVL block    DVLDEL */
        t = 0;                                               /* DVLDEL */
        for (jj=0;jj!=mm;jj++)                               /* DVLDEL */
          t+= pc[jj]*pdl[j+jj*mm];  /* Accumulate dot product   DVLDEL */
        pwork[j] = t;               /* Store in work array      DVLDEL */
      }                                                      /* DVLDEL */
        /* copy work array into c */
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                  /* DVLDEL */ 
      pc+=ldc;                      /* Skip column       VECDEL DVLDEL */
    }                                                 /* VECDEL DVLDEL */
  }                                                          /* DVLDEL */

  /* Finally:                                BLKPOSTPROCESS */
  /* Scale by alpha                          BLKPOSTPROCESS */
  /*    and add beta * c (from work)         BLKPOSTPROCESS */

  pc=c+ind_base;                          /* BLKPOSTPROCESS */
  pwork=work;                             /* BLKPOSTPROCESS    MTSBETADEL */
  for (l=0;l!=n;l++)                      /* BLKPOSTPROCESS        VECDEL */
    for (i=0;i!=m;i++) {                  /* BLKPOSTPROCESS */
          *pc *= alpha;                   /* BLKPOSTPROCESS    ALPHA1DEL  */
          *pc += *pwork;                  /* BLKPOSTPROCESS    MTSBETADEL */
          pc++;                           /* BLKPOSTPROCESS */
          pwork++;                        /* BLKPOSTPROCESS    MTSBETADEL */
    }                                     /* BLKPOSTPROCESS */


}
