#include "nsbcoo.h"

void COO_MatMult_CaABbC_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, const int *jndx,
                 const int nnz, 
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, 
                 const int ind_base)
{
  int i,j;
  double *pc=c;
  int l;                                               /* VECDEL */
  c-=ind_base;
  b-=ind_base;

  for (l=0;l!=n;l++)                            /* VECDEL MMBETADEL*/
    for (i=0;i!=m;i++) *pc++ *= beta;                  /* MMBETADEL*/

  for (l=0;l!=n;l++) {                                 /* VECDEL */
    for (j=0;j!=nnz;j++)
      c[indx[j]] += alpha * b[jndx[j]] * val[j];
    c += ldc; b += ldb;                                /* VECDEL */
  }                                                    /* VECDEL */
}


void COOsymm_MatMult_CaABbC_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, const int *jndx,
                 const int nnz, 
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, 
                 const int ind_base)
{
  double *pc=c;
  int l;                                               /* VECDEL */
  int i,j;
  int rowi,coli;
  c-=ind_base;
  b-=ind_base;

  for (l=0;l!=n;l++)                            /* VECDEL MMBETADEL*/
    for (i=0;i!=m;i++) *pc++ *= beta;                  /* MMBETADEL*/

  for (l=0;l!=n;l++) {                                 /* VECDEL */
    for (j=0;j!=nnz;j++) {
      rowi =  indx[j];
      coli =  jndx[j];
      if ( rowi ==  coli )
        c[rowi] += alpha * b[coli] * val[j];
      else {
        c[rowi] += alpha * b[coli] * val[j];
        c[coli] += alpha * b[rowi] * val[j];
      }
    }
    c += ldc; b += ldb;                                /* VECDEL */
  }                                                    /* VECDEL */
}


void COOskew_MatMult_CaABbC_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, const int *jndx,
                 const int nnz, 
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, 
                 const int ind_base)
{
  double *pc=c;
  int l;                                               /* VECDEL */
  int i,j;
  int rowi,coli;
  c-=ind_base;
  b-=ind_base;

  for (l=0;l!=n;l++)                            /* VECDEL MMBETADEL*/
    for (i=0;i!=m;i++) *pc++ *= beta;                  /* MMBETADEL*/

  for (l=0;l!=n;l++) {                                 /* VECDEL */
    for (j=0;j!=nnz;j++) {
      rowi =  indx[j];
      coli =  jndx[j];
      if ( rowi !=  coli ) {
        c[rowi] += alpha * b[coli] * val[j];
        c[coli] -= alpha * b[rowi] * val[j];
      }
    }
    c += ldc; b += ldb;                                /* VECDEL */
  }                                                    /* VECDEL */
}
