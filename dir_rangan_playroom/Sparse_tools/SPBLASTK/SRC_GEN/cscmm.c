#include "nsbcsc.h"

void CSC_MatMult_CaABbC_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc,  
                 const int ind_base)
{
  int i,j,jb,je;
  double *pc=c;                                       /* MMBETADEL */
  int l;                                               /* VECDEL */
  val-=ind_base;
  indx-=ind_base;
  c-=ind_base;

  for (l=0;l!=n;l++)                           /* VECDEL MMBETADEL */
    for (i=0;i!=m;i++) *pc++ *= beta;                 /* MMBETADEL */

  for (l=0;l!=n;l++) {                                 /* VECDEL */
    for (i=0;i!=k;i++){
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j!=je;j++)
        c[indx[j]] += alpha * b[i] * val[j];
    }
    c += ldc; b += ldb;                                /* VECDEL */
  }                                                    /* VECDEL */
}

void CSCsymm_MatMult_CaABbC_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc,  
                 const int ind_base)
{
  int i,j,jb,je;
  double *pc=c;                                       /* MMBETADEL */
  int l;                                               /* VECDEL */
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                           /* VECDEL MMBETADEL */
    for (i=0;i!=m;i++) *pc++ *= beta;                 /* MMBETADEL */

  for (l=0;l!=n;l++) {                                 /* VECDEL */
    for (j=0;j!=k;j++){
      jb = pntrb[j];
      je = pntre[j];
      if ( indx[jb] == j+ind_base ) {
        c[j] += alpha * b[j] * val[jb];
        jb++;
      }
      b-=ind_base;
      for (i=jb;i!=je;i++){
        c[j] += alpha * b[indx[i]] * val[i];
      }
      c-=ind_base;
      b+=ind_base;
      for (i=jb;i!=je;i++){
        c[indx[i]] += alpha * b[j] * val[i];
      }
      c+=ind_base;
    }
    c += ldc; b += ldb;                                /* VECDEL */
  }                                                    /* VECDEL */
}

void CSCskew_MatMult_CaABbC_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc,  
                 const int ind_base)
{
  int i,j,jb,je;
  double *pc=c;                                       /* MMBETADEL */
  int l;                                               /* VECDEL */
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                           /* VECDEL MMBETADEL */
    for (i=0;i!=m;i++) *pc++ *= beta;                 /* MMBETADEL */

  for (l=0;l!=n;l++) {                                 /* VECDEL */
    for (j=0;j!=k;j++){
      jb = pntrb[j];
      je = pntre[j];
      if ( indx[jb] == j+ind_base ) {
        jb++;
      }
      b-=ind_base;
      for (i=jb;i!=je;i++){
        c[j] -= alpha * b[indx[i]] * val[i];
      }
      c-=ind_base;
      b+=ind_base;
      for (i=jb;i!=je;i++){
        c[indx[i]] += alpha * b[j] * val[i];
      }
      c+=ind_base;
    }
    c += ldc; b += ldb;                                /* VECDEL */
  }                                                    /* VECDEL */
}

