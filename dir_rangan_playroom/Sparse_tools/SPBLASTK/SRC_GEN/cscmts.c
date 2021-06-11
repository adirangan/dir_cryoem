#include "nsbcsc.h"

void CSC_MatTriangSlvLD_CaDADBbC_double( 
                 const int m, const int n, const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             /* MTSBETADEL */
  int l;                                              /* VECDEL */
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                                  /* VECDEL */
    for (i=0;i!=m;i++){
      *pwork++ = beta * (*pc);                    /* MTSBETADEL */
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  pwork=work;                                        /* MTSBETADEL */
  c-=ind_base;
  for (l=0;l!=n;l++) {                                /* VECDEL */
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      jb = pntrb[i];
      je = pntre[i];
      z =  pc[i] / val[jb];
      pc[i] = z;
      for (j=jb+1;j<je;j++) {
        c[indx[j]] -= z*val[j];
      }
    }
    for (i=0;i!=m;i++) {                  /* POSTPROCESS */
        *pc *= alpha * dvl[i];            /* POSTPROCESS ALPHA1DVL */
        *pc += *pwork;                    /* POSTPROCESS MTSBETADEL */
        pwork++;                          /* POSTPROCESS MTSBETADEL */
        pc++;                             /* POSTPROCESS */
    }                                     /* POSTPROCESS */
    c += ldc; b += ldb;                               /* VECDEL */
  }                                                   /* VECDEL */
}

void CSC_MatTriangSlvLU_CaDADBbC_double(
                 const int m, const int n, const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             /* MTSBETADEL */
  int l;                                              /* VECDEL */
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                                  /* VECDEL */
    for (i=0;i!=m;i++){
      *pwork++ = beta * (*pc);                    /* MTSBETADEL */
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  pwork=work;                                        /* MTSBETADEL */
  c-=ind_base;
  for (l=0;l!=n;l++) {                                /* VECDEL */
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
       jb = pntrb[i];
       je = pntre[i];
       z = pc[i];
       for (j=jb;j<je;j++)
         c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  /* POSTPROCESS */
        *pc *= alpha * dvl[i];            /* POSTPROCESS ALPHA1DVL */
        *pc += *pwork;                    /* POSTPROCESS MTSBETADEL */
        pwork++;                          /* POSTPROCESS MTSBETADEL */
        pc++;                             /* POSTPROCESS */
    }                                     /* POSTPROCESS */
    c += ldc; b += ldb;                               /* VECDEL */
  }                                                   /* VECDEL */
}

void CSC_MatTriangSlvUD_CaDADBbC_double(
                 const int m, const int n, const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je, index;
  double *pc=c;
  double *pwork=work;                             /* MTSBETADEL */
  int l;                                              /* VECDEL */
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                                  /* VECDEL */
    for (i=0;i!=m;i++){
      *pwork++ = beta * (*pc);                    /* MTSBETADEL */
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  pwork=work;                                        /* MTSBETADEL */
  c-=ind_base;
  for (l=0;l!=n;l++) {                                /* VECDEL */
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i]-1;
      z = pc[i] /val[je];
      pc[i] = z;
      for (j=jb;j<je;j++) {
          c[indx[j]] -= z * val[j];
      }
    }
    for (i=0;i!=m;i++) {                  /* POSTPROCESS */
        *pc *= alpha * dvl[i];            /* POSTPROCESS ALPHA1DVL */
        *pc += *pwork;                    /* POSTPROCESS MTSBETADEL */
        pwork++;                          /* POSTPROCESS MTSBETADEL */
        pc++;                             /* POSTPROCESS */
    }                                     /* POSTPROCESS */
    c += ldc; b += ldb;                               /* VECDEL */
  }                                                   /* VECDEL */
}

void CSC_MatTriangSlvUU_CaDADBbC_double(
                 const int m, const int n, const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             /* MTSBETADEL */
  int l;                                              /* VECDEL */
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                                  /* VECDEL */
    for (i=0;i!=m;i++){
      *pwork++ = beta * (*pc);                    /* MTSBETADEL */
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  pwork=work;                                        /* MTSBETADEL */
  c-=ind_base;
  for (l=0;l!=n;l++) {                                /* VECDEL */
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i];
      z = pc[i];
      pc[i] = z;
      for (j=jb;j<je;j++)
        c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  /* POSTPROCESS */
        *pc *= alpha * dvl[i];            /* POSTPROCESS ALPHA1DVL */
        *pc += *pwork;                    /* POSTPROCESS MTSBETADEL */
        pwork++;                          /* POSTPROCESS MTSBETADEL */
        pc++;                             /* POSTPROCESS */
    }                                     /* POSTPROCESS */
    c += ldc; b += ldb;                               /* VECDEL */
  }                                                   /* VECDEL */
}

