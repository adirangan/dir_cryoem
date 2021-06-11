#include "nsbcsr.h"

void CSR_MatMult_CaABbC_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, 
                 const int ind_base)
{
  double t;
  const double *pval;
  double *pc=c;                                       /* MMBETADEL */
  int i,j,jb,je;
  int l;                                               /* VECDEL */
  b-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                           /* VECDEL MMBETADEL */
    for (i=0;i!=m;i++) *pc++ *= beta;                 /* MMBETADEL */

  for (l=0;l!=n;l++) {                                 /* VECDEL */
    pval = val;
    for (i=0;i!=m;i++) {
      t = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j!=je;j++)
        t += alpha * b[indx[j]] * (*pval++);
      c[i] += t;
    }
    c += ldc; b += ldb;                                /* VECDEL */
  }                                                    /* VECDEL */
}

void CSRsymm_MatMult_CaABbC_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, 
                 const int ind_base)
{
  const double *pval;
  double *pc=c;                                       /* MMBETADEL */
  int i,j;
  int l;                                               /* VECDEL */
  int jj;
  int rpntrb, rpntre;
  int index, nvals;
  indx-=ind_base;
      
  for (l=0;l!=n;l++)                           /* VECDEL MMBETADEL */
    for (i=0;i!=m;i++) *pc++ *= beta;                 /* MMBETADEL */

  for (l=0;l!=n;l++) {                                 /* VECDEL */
    pval = val;
    for (j=0;j!=k;j++){
      rpntrb = pntrb[j];
      rpntre = pntre[j];
      for (jj=rpntrb;jj!=rpntre;jj++) {
        index = indx[jj];
        if ( index == j+ind_base ) {
          c[j] += alpha * b[j] * (*pval++);
          continue;
        }
        c-=ind_base;
        c[index] += alpha * b[j] * (*pval);
        c+=ind_base;
        b-=ind_base;
        c[j] += alpha * b[index] * (*pval++);
        b+=ind_base;
      }
    }
    c += ldc; b += ldb;                                /* VECDEL */
  }                                                    /* VECDEL */
}



void CSRskew_MatMult_CaABbC_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, 
                 const int ind_base)
{
  const double *pval;
  double *pc=c;                                       /* MMBETADEL */
  int i,j;
  int l;                                               /* VECDEL */
  int jj;
  int rpntrb, rpntre;
  int index, nvals;
  indx-=ind_base;
      
  for (l=0;l!=n;l++)                           /* VECDEL MMBETADEL */
    for (i=0;i!=m;i++) *pc++ *= beta;                 /* MMBETADEL */
  
  for (l=0;l!=n;l++) {                                 /* VECDEL */
    pval = val;
    for (j=0;j!=k;j++){
      rpntrb = pntrb[j];
      rpntre = pntre[j];
      for (jj=rpntrb;jj!=rpntre;jj++) {
        index = indx[jj];
        if ( index == j+ind_base ) {
          *pval++;
          continue;
        }
        c-=ind_base;
        c[index] -= alpha * b[j] * (*pval);
        c+=ind_base;
        b-=ind_base;
        c[j] += alpha * b[index] * (*pval++);
        b+=ind_base;
      }
    }
    c += ldc; b += ldb;                                 /* VECDEL */
  }                                                     /* VECDEL */
}
