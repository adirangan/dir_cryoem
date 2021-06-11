#include "nsbcsr.h"

void CSR_MatTriangSlvLD_CaDADBbC_double(
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
  double valtmp;
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                       /* VECDEL MTSBETADEL */
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++); /* MTSBETADEL */

  pwork=work;                                        /* MTSBETADEL */
  c-=ind_base;
  for (l=0;l!=n;l++) {                                /* VECDEL */
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        index = indx[j];
        if ( index == i+ind_base ) {
          valtmp = val[j];
        } else {
          z += c[index] * val[j];
        }
      }
      pc[i] = (dvr[i]*b[i] - z) / valtmp;
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

void CSR_MatTriangSlvLU_CaDADBbC_double(
                 const int m, const int n, const double *dvl, 
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                /* MTSBETADEL */
  int l;                                              /* VECDEL */
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          /* VECDEL MTSBETADEL */
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    /* MTSBETADEL */

  pwork=work;                                        /* MTSBETADEL */
  c-=ind_base;
  for (l=0;l!=n;l++) {                                /* VECDEL */
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = dvr[i]*b[i] - z;
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

void CSR_MatTriangSlvUD_CaDADBbC_double(
                 const int m, const int n, const double *dvl, 
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je, index;
  double *pc=c;
  double *pwork=work;                                /* MTSBETADEL */
  double valtmp;
  int l;                                              /* VECDEL */
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          /* VECDEL MTSBETADEL */
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    /* MTSBETADEL */

  pwork=work;                                        /* MTSBETADEL */
  c-=ind_base;
  for (l=0;l!=n;l++) {                                /* VECDEL */
    pc=c+ind_base;
    for (i=m-1;i!=-1; i--) {
      z = 0;
      jb = pntrb[i];
      je =  pntre[i];
      for (j=jb+1; j<je; j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = (dvr[i]*b[i] - z) / val[jb];
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

void CSR_MatTriangSlvUU_CaDADBbC_double(
                 const int m, const int n, const double *dvl, 
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                /* MTSBETADEL */
  int l;                                              /* VECDEL */
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          /* VECDEL MTSBETADEL */
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    /* MTSBETADEL */

  pwork=work;                                        /* MTSBETADEL */
  c-=ind_base;
  for (l=0;l!=n;l++) {                                /* VECDEL */
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      z = 0;
      jb =  pntrb[i];
      je =  pntre[i];
      for (j=jb; j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = dvr[i]*b[i] - z;
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

