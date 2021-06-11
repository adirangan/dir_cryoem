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

#include <stdlib.h>
#include <stdio.h>
#include "spblas.h"

/* Sparse BLAS Toolkit interface routine: */
void  dbscmm_(
             const int *transa, const int *mb, const int *n, const int *kb,
             const double *alpha, const int descra[], const double val[],
             const int bindx[], const int bpntrb[], const int bpntre[],
             const int *lb, const double b[], const int *ldb,
             const double *beta, double c[], const int *ldc,
             double work[], const int *lwork)
{
    dbscmm( *transa,   *mb,   *n,   *kb, *alpha,   descra,   val,
               bindx, bpntrb,  bpntre, *lb, 
               b,   *ldb, *beta,  c,   *ldc,
               work,   *lwork);
}
