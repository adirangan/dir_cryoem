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
#include "dbscmml.h"
#include "dbscvml.h"

/* Sparse BLAS Toolkit interface routine: */
void  dbscmm(
             const int transa, const int mb, const int n, const int kb,
             const double alpha, const int descra[], const double val[],
             const int bindx[], const int bpntrb[], const int bpntre[],
             const int lb, const double b[], const int ldb,
             const double beta, double c[], const int ldc,
             double work[], const int lwork)
{
/* ------------ begin interface description ------------
   Toolkit interface:
   dbscmm -- block sparse column format matrix-matrix multiply
  
   C <- alpha A B + beta C
  
   Arguments:
  
   int transa	Indicates how to operate with the sparse matrix
  		0 : operate with matrix
  		1 : operate with transpose matrix
  
   int mb	Number of block rows in matrix A
  
   int n	Number of columns in matrix c
  
   int kb	Number of block columns in matrix A
  
   double alpha Scalar parameter
  
   double beta  Scalar parameter
  
   int descra[]	Descriptor argument.  Nine element integer array
  		descra[0] matrix structure
  			0 : general
  			1 : symmetric
  			2 : Hermitian
  			3 : Triangular
  			4 : Skew(Anti)-Symmetric
  			5 : Diagonal
  		descra[1] upper/lower triangular indicator
  			1 : lower
  			2 : upper
  		descra[2] main diagonal type
  			0 : non-unit
  			1 : unit
  		descra[3] Array base 
  			0 : C/C++ compatible
  			1 : Fortran compatible
  		descra[4] repeated indices?
  			0 : unknown
  			1 : no repeated indices
  
  
   double *val	scalar array of length nnz containing matrix entries
  
   int *bindx	integer array of length bnnz consisting of the block row   
   		indices of the entries of A.
  
   int *bpntrb	integer array of length mb such that bpntrb(i)-bpntrb(1)
                points to location in bindx of the first block entry of 
  		the j-th column of A.
  
   int *bpntre	integer array of length mb such that bpntre(i)-bpntrb(1)
                points to location in bindx of the last block entry of
  		the j-th column of A.
  
   int lb	dimension of blocks
  
   double *b	rectangular array with first dimension ldb
  
   double *c	rectangular array with first dimension ldc
  
   double *work	scratch array of length lwork.  lwork should be at least
  		max(m,n)
  
   ------------ end interface description --------------*/
int ind_base = descra[3];
int m=mb*lb;
int k=kb*lb;

if (alpha == 0.0) {
   ScaleArray_double(m, n, c, ldc, beta);
   return;
} 

switch ( descra[0] ) {
case 1: /* Symmetric */
case 2: /* Hermitian (for real same as Symmetric) */
  if ( m != k ) {
    printf("In dbscmm: inconsistant dimensions for a symmetric matrix");
    printf("m = %d  k = %d\nExiting...\n",m,k);
    exit(-1); 
  }
  switch ( descra[1] ) {
  case 2: /* Upper triangular stored, or */
  case 1: /* Lower triangular stored (same for both) */
    switch ( n ) {
    case 1:
      if (alpha == 1) {
        if (beta == 1) {
          BSCsymm_VecMult_CABC_double(mb, kb, 
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, 
                                  c, ind_base);
        } else if (beta == 0) {
          BSCsymm_VecMult_CAB_double(mb, kb, 
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, 
                                  c, ind_base);
        } else  { /*  beta is general nonzero */
          BSCsymm_VecMult_CABbC_double(mb, kb, 
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, beta,
                                  c, ind_base);
        }
      } else { /* alpha is general nonzero */
        if (beta == 1) {
          BSCsymm_VecMult_CaABC_double(mb, kb, alpha,
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, 
                                  c, ind_base);
        } else if (beta == 0) {
          BSCsymm_VecMult_CaAB_double(mb, kb, alpha,
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, 
                                  c, ind_base);
        } else { /*  beta is general nonzero */
          BSCsymm_VecMult_CaABbC_double(mb, kb, alpha,
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, beta,
                                  c, ind_base);
        }
      }
      break;
    default:  /* n is greater than 1 -- doing Mat Mult */
      if (alpha == 1) {
        if (beta == 1) {
          BSCsymm_MatMult_CABC_double(mb, n, kb, 
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, ldb, 
                                  c, ldc, ind_base);
        } else if (beta == 0) {
          BSCsymm_MatMult_CAB_double(mb, n, kb, 
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, ldb, 
                                  c, ldc, ind_base);
        } else { /*  beta is general nonzero */
          BSCsymm_MatMult_CABbC_double(mb, n, kb, 
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, ldb, beta,
                                  c, ldc, ind_base);
        }
      } else {   /* alpha is general nonzero */
        if (beta == 1) {
          BSCsymm_MatMult_CaABC_double(mb, n, kb, alpha,
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, ldb, 
                                  c, ldc, ind_base);
        } else if (beta == 0) {
          BSCsymm_MatMult_CaAB_double(mb, n, kb, alpha,
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, ldb, 
                                  c, ldc, ind_base);
        } else { /*  beta is general nonzero */
          BSCsymm_MatMult_CaABbC_double(mb, n, kb, alpha,
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, ldb, beta,
                                  c, ldc, ind_base);
        }
      }
      break;
    }
    break;
  default:
    printf("Invalid argument descra[1] in dbscmm. Use 1 or 2. \n");
    break;
  } /* end of switch on descra[1] */
  break;
case 4: /* Skew Symmetric */
  if ( m != k ) {
    printf("In dbscmm: inconsistant dimensions for a skew-symmetric matrix");
    printf("m = %d  k = %d\nExiting...\n",m,k);
    exit(-1); 
  }
  switch ( transa ) {
  case 0:
    switch ( n ) {
    case 1:
      if (alpha == 1) {
        if (beta == 1) {
          BSCskew_VecMult_CABC_double(mb, kb, 
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, 
                                  c, ind_base);
        } else if (beta == 0) {
          BSCskew_VecMult_CAB_double(mb, kb, 
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, 
                                  c, ind_base);
        } else  { /*  beta is general nonzero */
          BSCskew_VecMult_CABbC_double(mb, kb, 
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, beta,
                                  c, ind_base);
        }
      } else { /* alpha is general nonzero */
        if (beta == 1) {
          BSCskew_VecMult_CaABC_double(mb, kb, alpha,
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, 
                                  c, ind_base);
        } else if (beta == 0) {
          BSCskew_VecMult_CaAB_double(mb, kb, alpha,
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, 
                                  c, ind_base);
        } else { /*  beta is general nonzero */
          BSCskew_VecMult_CaABbC_double(mb, kb, alpha,
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, beta,
                                  c, ind_base);
        }
      }
      break;
    default:
      if (alpha == 1) {
        if (beta == 1) {
          BSCskew_MatMult_CABC_double(mb, n, kb, 
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, ldb, 
                                  c, ldc, ind_base);
        } else if (beta == 0) {
          BSCskew_MatMult_CAB_double(mb, n, kb, 
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, ldb, 
                                  c, ldc, ind_base);
        } else { /*  beta is general nonzero */
          BSCskew_MatMult_CABbC_double(mb, n, kb, 
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, ldb, beta,
                                  c, ldc, ind_base);
        }
      } else {   /* alpha is general nonzero */
        if (beta == 1) {
          BSCskew_MatMult_CaABC_double(mb, n, kb, alpha,
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, ldb, 
                                  c, ldc, ind_base);
        } else if (beta == 0) {
          BSCskew_MatMult_CaAB_double(mb, n, kb, alpha,
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, ldb, 
                                  c, ldc, ind_base);
        } else { /*  beta is general nonzero */
          BSCskew_MatMult_CaABbC_double(mb, n, kb, alpha,
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, ldb, beta,
                                  c, ldc, ind_base);
        }
      }
      break;
    } /* end switch on n */
    break;
  case 1:
    switch ( n ) {
    case 1:
      if (alpha == 1) {
        if (beta == 1) {
          BSCskew_VecMult_CATBC_double(mb, kb, 
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, 
                                  c, ind_base);
        } else if (beta == 0) {
          BSCskew_VecMult_CATB_double(mb, kb, 
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, 
                                  c, ind_base);
        } else  { /*  beta is general nonzero */
          BSCskew_VecMult_CATBbC_double(mb, kb, 
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, beta,
                                  c, ind_base);
        }
      } else { /* alpha is general nonzero */
        if (beta == 1) {
          BSCskew_VecMult_CaATBC_double(mb, kb, alpha,
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, 
                                  c, ind_base);
        } else if (beta == 0) {
          BSCskew_VecMult_CaATB_double(mb, kb, alpha,
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, 
                                  c, ind_base);
        } else { /*  beta is general nonzero */
          BSCskew_VecMult_CaATBbC_double(mb, kb, alpha,
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, beta,
                                  c, ind_base);
        }
      }
      break;
    default:  /* n is greater than 1 -- doing Mat Mult */
      if (alpha == 1) {
        if (beta == 1) {
          BSCskew_MatMult_CATBC_double(mb, n, kb, 
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, ldb, 
                                  c, ldc, ind_base);
        } else if (beta == 0) {
          BSCskew_MatMult_CATB_double(mb, n, kb, 
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, ldb, 
                                  c, ldc, ind_base);
        } else { /*  beta is general nonzero */
          BSCskew_MatMult_CATBbC_double(mb, n, kb, 
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, ldb, beta,
                                  c, ldc, ind_base);
        }
      } else {   /* alpha is general nonzero */
        if (beta == 1) {
          BSCskew_MatMult_CaATBC_double(mb, n, kb, alpha,
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, ldb, 
                                  c, ldc, ind_base);
        } else if (beta == 0) {
          BSCskew_MatMult_CaATB_double(mb, n, kb, alpha,
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, ldb, 
                                  c, ldc, ind_base);
        } else { /*  beta is general nonzero */
          BSCskew_MatMult_CaATBbC_double(mb, n, kb, alpha,
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, ldb, beta,
                                  c, ldc, ind_base);
        }
      }
      break;
    } /* end of switch on n */
    break;
  default:
    printf("Invalid argument transa in dbscmm. Use 0 or 1. \n");
    break;
  } /* end switch on transa */
  break;
case 0: case 3: case 5:
  switch ( transa ) {
  case 0:
    switch ( n ) {
    case 1:
      if (alpha == 1) {
        if (beta == 1) {
          BSC_VecMult_CABC_double(mb, kb, 
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, 
                                  c, ind_base);
        } else if (beta == 0) {
          BSC_VecMult_CAB_double(mb, kb, 
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, 
                                  c, ind_base);
        } else  { /*  beta is general nonzero */
          BSC_VecMult_CABbC_double(mb, kb, 
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, beta,
                                  c, ind_base);
        }
      } else { /* alpha is general nonzero */
        if (beta == 1) {
          BSC_VecMult_CaABC_double(mb, kb, alpha,
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, 
                                  c, ind_base);
        } else if (beta == 0) {
          BSC_VecMult_CaAB_double(mb, kb, alpha,
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, 
                                  c, ind_base);
        } else { /*  beta is general nonzero */
          BSC_VecMult_CaABbC_double(mb, kb, alpha,
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, beta,
                                  c, ind_base);
        }
      }
      break;
    default:  /* n is greater than 1 -- doing Mat Mult */
      if (alpha == 1) {
        if (beta == 1) {
          BSC_MatMult_CABC_double(mb, n, kb, 
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, ldb, 
                                  c, ldc, ind_base);
        } else if (beta == 0) {
          BSC_MatMult_CAB_double(mb, n, kb, 
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, ldb, 
                                  c, ldc, ind_base);
        } else { /*  beta is general nonzero */
          BSC_MatMult_CABbC_double(mb, n, kb, 
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, ldb, beta,
                                  c, ldc, ind_base);
        }
      } else {   /* alpha is general nonzero */
        if (beta == 1) {
          BSC_MatMult_CaABC_double(mb, n, kb, alpha,
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, ldb, 
                                  c, ldc, ind_base);
        } else if (beta == 0) {
          BSC_MatMult_CaAB_double(mb, n, kb, alpha,
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, ldb, 
                                  c, ldc, ind_base);
        } else { /*  beta is general nonzero */
          BSC_MatMult_CaABbC_double(mb, n, kb, alpha,
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, ldb, beta,
                                  c, ldc, ind_base);
        }
      }
      break;
    } /* end switch on n */
    break;
  case 1:  /* operate with transpose */
    switch ( n ) {
    case 1:
      if (alpha == 1) {
        if (beta == 1) {
          BSC_VecMult_CATBC_double(mb, kb, 
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, 
                                  c, ind_base);
        } else if (beta == 0) {
          BSC_VecMult_CATB_double(mb, kb, 
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, 
                                  c, ind_base);
        } else  { /*  beta is general nonzero */
          BSC_VecMult_CATBbC_double(mb, kb, 
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, beta,
                                  c, ind_base);
        }
      } else { /* alpha is general nonzero */
        if (beta == 1) {
          BSC_VecMult_CaATBC_double(mb, kb, alpha,
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, 
                                  c, ind_base);
        } else if (beta == 0) {
          BSC_VecMult_CaATB_double(mb, kb, alpha,
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, 
                                  c, ind_base);
        } else { /*  beta is general nonzero */
          BSC_VecMult_CaATBbC_double(mb, kb, alpha,
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, beta,
                                  c, ind_base);
        }
      }
      break;
    default:  /* n is greater than 1 -- doing Mat Mult */
      if (alpha == 1) {
        if (beta == 1) {
          BSC_MatMult_CATBC_double(mb, n, kb, 
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, ldb, 
                                  c, ldc, ind_base);
        } else if (beta == 0) {
          BSC_MatMult_CATB_double(mb, n, kb, 
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, ldb, 
                                  c, ldc, ind_base);
        } else { /*  beta is general nonzero */
          BSC_MatMult_CATBbC_double(mb, n, kb, 
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, ldb, beta,
                                  c, ldc, ind_base);
        }
      } else {   /* alpha is general nonzero */
        if (beta == 1) {
          BSC_MatMult_CaATBC_double(mb, n, kb, alpha,
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, ldb, 
                                  c, ldc, ind_base);
        } else if (beta == 0) {
          BSC_MatMult_CaATB_double(mb, n, kb, alpha,
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, ldb, 
                                  c, ldc, ind_base);
        } else { /*  beta is general nonzero */
          BSC_MatMult_CaATBbC_double(mb, n, kb, alpha,
                                  val, bindx,
                                  bpntrb, bpntre, lb, b, ldb, beta,
                                  c, ldc, ind_base);
        }
      }
      break;
    } /* end switch on n */
    break;
  default: 
    printf("Invalid argument transa in dbscmm. Use 0 or 1. \n");
    break;
  } /* end switch on transa */
  break;
default:
    printf("Invalid argument descra[0] in dbscmm. Use 0 - 5. \n");
    break;
}  
  
  
}
   
  
  
   
  
