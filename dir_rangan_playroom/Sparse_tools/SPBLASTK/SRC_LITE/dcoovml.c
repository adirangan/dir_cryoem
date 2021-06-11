/*-------------------------------------------------------|
|  NIST SPARSE BLAS v. 0.9 (Sat Jul 6 14:27:21 EDT 1996) |
|                                                        |
|  Authors:                                              |
|     Karin A. Remington and Roldan Pozo                 |
|     National Institute of Standards and Technology     |
|                                                        |
|  Based on the interface standard proposed in:          | 
|   "A Revised Proposal for a Sparse BLAS Toolkit" by    |
|    S. Carney and K. Wu -- University of Minnesota      |
|    M. Heroux and G. Li -- Cray Research                |  
|    R. Pozo and K.A. Remington -- NIST                  |
|                                                        |
|  Contact:                                              |
|     Karin A. Remington, email: kremington@nist.gov     |
--------------------------------------------------------*/


/* Created:  Sat Jul 6 14:28:00 EDT 1996 */

#include "dcoovml.h"



void COO_VecMult_CAB_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, const int *jndx,
                 const int nnz, 
                 const double *b,  
                 double *c,  
                 const int ind_base)
{
  int i,j;
  double *pc=c;
  c-=ind_base;
  b-=ind_base;

    for (i=0;i!=m;i++) *pc++ = 0;                  

    for (j=0;j!=nnz;j++)
      c[indx[j]] +=  b[jndx[j]] * val[j];
}


void COOsymm_VecMult_CAB_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, const int *jndx,
                 const int nnz, 
                 const double *b,  
                 double *c,  
                 const int ind_base)
{
  double *pc=c;
  int i,j;
  int rowi,coli;
  c-=ind_base;
  b-=ind_base;

    for (i=0;i!=m;i++) *pc++ = 0;                  

    for (j=0;j!=nnz;j++) {
      rowi =  indx[j];
      coli =  jndx[j];
      if ( rowi ==  coli )
        c[rowi] +=  b[coli] * val[j];
      else {
        c[rowi] +=  b[coli] * val[j];
        c[coli] +=  b[rowi] * val[j];
      }
    }
}


void COOskew_VecMult_CAB_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, const int *jndx,
                 const int nnz, 
                 const double *b,  
                 double *c,  
                 const int ind_base)
{
  double *pc=c;
  int i,j;
  int rowi,coli;
  c-=ind_base;
  b-=ind_base;

    for (i=0;i!=m;i++) *pc++ = 0;                  

    for (j=0;j!=nnz;j++) {
      rowi =  indx[j];
      coli =  jndx[j];
      if ( rowi !=  coli ) {
        c[rowi] +=  b[coli] * val[j];
        c[coli] -=  b[rowi] * val[j];
      }
    }
}



void COO_VecMult_CaAB_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, const int *jndx,
                 const int nnz, 
                 const double *b,  
                 double *c,  
                 const int ind_base)
{
  int i,j;
  double *pc=c;
  c-=ind_base;
  b-=ind_base;

    for (i=0;i!=m;i++) *pc++ = 0;                  

    for (j=0;j!=nnz;j++)
      c[indx[j]] += alpha * b[jndx[j]] * val[j];
}


void COOsymm_VecMult_CaAB_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, const int *jndx,
                 const int nnz, 
                 const double *b,  
                 double *c,  
                 const int ind_base)
{
  double *pc=c;
  int i,j;
  int rowi,coli;
  c-=ind_base;
  b-=ind_base;

    for (i=0;i!=m;i++) *pc++ = 0;                  

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
}


void COOskew_VecMult_CaAB_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, const int *jndx,
                 const int nnz, 
                 const double *b,  
                 double *c,  
                 const int ind_base)
{
  double *pc=c;
  int i,j;
  int rowi,coli;
  c-=ind_base;
  b-=ind_base;

    for (i=0;i!=m;i++) *pc++ = 0;                  

    for (j=0;j!=nnz;j++) {
      rowi =  indx[j];
      coli =  jndx[j];
      if ( rowi !=  coli ) {
        c[rowi] += alpha * b[coli] * val[j];
        c[coli] -= alpha * b[rowi] * val[j];
      }
    }
}



void COO_VecMult_CABC_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, const int *jndx,
                 const int nnz, 
                 const double *b,  
                 double *c,  
                 const int ind_base)
{
  int i,j;
  double *pc=c;
  c-=ind_base;
  b-=ind_base;


    for (j=0;j!=nnz;j++)
      c[indx[j]] +=  b[jndx[j]] * val[j];
}


void COOsymm_VecMult_CABC_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, const int *jndx,
                 const int nnz, 
                 const double *b,  
                 double *c,  
                 const int ind_base)
{
  double *pc=c;
  int i,j;
  int rowi,coli;
  c-=ind_base;
  b-=ind_base;


    for (j=0;j!=nnz;j++) {
      rowi =  indx[j];
      coli =  jndx[j];
      if ( rowi ==  coli )
        c[rowi] +=  b[coli] * val[j];
      else {
        c[rowi] +=  b[coli] * val[j];
        c[coli] +=  b[rowi] * val[j];
      }
    }
}


void COOskew_VecMult_CABC_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, const int *jndx,
                 const int nnz, 
                 const double *b,  
                 double *c,  
                 const int ind_base)
{
  double *pc=c;
  int i,j;
  int rowi,coli;
  c-=ind_base;
  b-=ind_base;


    for (j=0;j!=nnz;j++) {
      rowi =  indx[j];
      coli =  jndx[j];
      if ( rowi !=  coli ) {
        c[rowi] +=  b[coli] * val[j];
        c[coli] -=  b[rowi] * val[j];
      }
    }
}



void COO_VecMult_CaABC_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, const int *jndx,
                 const int nnz, 
                 const double *b,  
                 double *c,  
                 const int ind_base)
{
  int i,j;
  double *pc=c;
  c-=ind_base;
  b-=ind_base;


    for (j=0;j!=nnz;j++)
      c[indx[j]] += alpha * b[jndx[j]] * val[j];
}


void COOsymm_VecMult_CaABC_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, const int *jndx,
                 const int nnz, 
                 const double *b,  
                 double *c,  
                 const int ind_base)
{
  double *pc=c;
  int i,j;
  int rowi,coli;
  c-=ind_base;
  b-=ind_base;


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
}


void COOskew_VecMult_CaABC_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, const int *jndx,
                 const int nnz, 
                 const double *b,  
                 double *c,  
                 const int ind_base)
{
  double *pc=c;
  int i,j;
  int rowi,coli;
  c-=ind_base;
  b-=ind_base;


    for (j=0;j!=nnz;j++) {
      rowi =  indx[j];
      coli =  jndx[j];
      if ( rowi !=  coli ) {
        c[rowi] += alpha * b[coli] * val[j];
        c[coli] -= alpha * b[rowi] * val[j];
      }
    }
}



void COO_VecMult_CABbC_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, const int *jndx,
                 const int nnz, 
                 const double *b,  const double beta,
                 double *c,  
                 const int ind_base)
{
  int i,j;
  double *pc=c;
  c-=ind_base;
  b-=ind_base;

    for (i=0;i!=m;i++) *pc++ *= beta;                  

    for (j=0;j!=nnz;j++)
      c[indx[j]] +=  b[jndx[j]] * val[j];
}


void COOsymm_VecMult_CABbC_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, const int *jndx,
                 const int nnz, 
                 const double *b,  const double beta,
                 double *c,  
                 const int ind_base)
{
  double *pc=c;
  int i,j;
  int rowi,coli;
  c-=ind_base;
  b-=ind_base;

    for (i=0;i!=m;i++) *pc++ *= beta;                  

    for (j=0;j!=nnz;j++) {
      rowi =  indx[j];
      coli =  jndx[j];
      if ( rowi ==  coli )
        c[rowi] +=  b[coli] * val[j];
      else {
        c[rowi] +=  b[coli] * val[j];
        c[coli] +=  b[rowi] * val[j];
      }
    }
}


void COOskew_VecMult_CABbC_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, const int *jndx,
                 const int nnz, 
                 const double *b,  const double beta,
                 double *c,  
                 const int ind_base)
{
  double *pc=c;
  int i,j;
  int rowi,coli;
  c-=ind_base;
  b-=ind_base;

    for (i=0;i!=m;i++) *pc++ *= beta;                  

    for (j=0;j!=nnz;j++) {
      rowi =  indx[j];
      coli =  jndx[j];
      if ( rowi !=  coli ) {
        c[rowi] +=  b[coli] * val[j];
        c[coli] -=  b[rowi] * val[j];
      }
    }
}



void COO_VecMult_CaABbC_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, const int *jndx,
                 const int nnz, 
                 const double *b,  const double beta,
                 double *c,  
                 const int ind_base)
{
  int i,j;
  double *pc=c;
  c-=ind_base;
  b-=ind_base;

    for (i=0;i!=m;i++) *pc++ *= beta;                  

    for (j=0;j!=nnz;j++)
      c[indx[j]] += alpha * b[jndx[j]] * val[j];
}


void COOsymm_VecMult_CaABbC_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, const int *jndx,
                 const int nnz, 
                 const double *b,  const double beta,
                 double *c,  
                 const int ind_base)
{
  double *pc=c;
  int i,j;
  int rowi,coli;
  c-=ind_base;
  b-=ind_base;

    for (i=0;i!=m;i++) *pc++ *= beta;                  

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
}


void COOskew_VecMult_CaABbC_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, const int *jndx,
                 const int nnz, 
                 const double *b,  const double beta,
                 double *c,  
                 const int ind_base)
{
  double *pc=c;
  int i,j;
  int rowi,coli;
  c-=ind_base;
  b-=ind_base;

    for (i=0;i!=m;i++) *pc++ *= beta;                  

    for (j=0;j!=nnz;j++) {
      rowi =  indx[j];
      coli =  jndx[j];
      if ( rowi !=  coli ) {
        c[rowi] += alpha * b[coli] * val[j];
        c[coli] -= alpha * b[rowi] * val[j];
      }
    }
}
