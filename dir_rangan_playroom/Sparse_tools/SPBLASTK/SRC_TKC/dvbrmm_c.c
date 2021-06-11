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
#include "dvbrmml.h"
#include "dvbrvml.h"

/* Sparse BLAS Toolkit interface routine: */
void  dvbrmm(
             const int transa, const int mb, const int n, const int kb,
             const double alpha, const int descra[], const double val[],
             const int indx[], const int bindx[], 
             const int rpntr[], const int cpntr[], 
             const int bpntrb[], const int bpntre[],
             const double b[], const int ldb,
             const double beta, double c[], const int ldc,
             double work[], const int lwork)
{
/* ------------ begin interface description ------------
   Toolkit interface:
   dvbrmm -- variable block sparse row format matrix-matrix multiply
  
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
  
   int *indx	integer array of length bnnz+1 such that the i-th element of
   		indx[] points to the location in val of the (1,1) element
  		of the i-th block entry.
  
   int *bindx	integer array of length bnnz consisting of the block column
   		indices of the entries of A.
  
   int *rpntr	integer array of length mb+1 such that rpntr(i)-rpntr(1)
  		is the row index of the first point row in the i-th block row. 
  		rpntr(mb+1) is set to m+rpntr(1).
  		Thus, the number of point rows in the i-th block row is
  		rpntr(i+1)-rpntr(i).
  
   int *cpntr	integer array of length kb+1 such that cpntr(j)-cpntr(1)
  		is the column index of the first point column in the j-th
  		block column. cpntr(kb+1) is set to k+cpntr(1).
  		Thus, the number of point columns in the j-th block column is
  		cpntr(j+1)-cpntr(j).
  
   int *bpntrb	integer array of length mb such that bpntrb(i)-bpntrb(1)
                points to location in bindx of the first block entry of 
  		the j-th row of A.
  
   int *bpntre	integer array of length mb such that bpntre(i)-bpntrb(1)
                points to location in bindx of the last block entry of
  		the j-th row of A.
  
   double *b	rectangular array with first dimension ldb
  
   double *c	rectangular array with first dimension ldc
  
   double *work	scratch array of length lwork.  lwork should be at least
  		max(m,n)
  
   ------------ end interface description --------------*/
int ind_base = descra[3];
int m=rpntr[mb]-rpntr[0];
int k=cpntr[kb]-cpntr[0];

if (alpha == 0.0) {
   ScaleArray_double(m, n, c, ldc, beta);
   return;
} 

switch ( descra[0] ) {
case 1: /* Symmetric */
case 2: /* Hermitian (for real same as Symmetric) */
  if ( m != k ) {
    printf("In dvbrmm: inconsistant dimensions for a symmetric matrix");
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
          VBRsymm_VecMult_CABC_double(mb, kb, 
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, 
                                  c, ind_base);
        } else if (beta == 0) {
          VBRsymm_VecMult_CAB_double(mb, kb, 
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, 
                                  c, ind_base);
        } else  { /*  beta is general nonzero */
          VBRsymm_VecMult_CABbC_double(mb, kb, 
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, beta,
                                  c, ind_base);
        }
      } else { /* alpha is general nonzero */
        if (beta == 1) {
          VBRsymm_VecMult_CaABC_double(mb, kb, alpha,
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, 
                                  c, ind_base);
        } else if (beta == 0) {
          VBRsymm_VecMult_CaAB_double(mb, kb, alpha,
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, 
                                  c, ind_base);
        } else { /*  beta is general nonzero */
          VBRsymm_VecMult_CaABbC_double(mb, kb, alpha,
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, beta,
                                  c, ind_base);
        }
      }
      break;
    default:  /* n is greater than 1 -- doing Mat Mult */
      if (alpha == 1) {
        if (beta == 1) {
          VBRsymm_MatMult_CABC_double(mb, n, kb, 
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, ldb, 
                                  c, ldc, ind_base);
        } else if (beta == 0) {
          VBRsymm_MatMult_CAB_double(mb, n, kb, 
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, ldb, 
                                  c, ldc, ind_base);
        } else { /*  beta is general nonzero */
          VBRsymm_MatMult_CABbC_double(mb, n, kb, 
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, ldb, beta,
                                  c, ldc, ind_base);
        }
      } else {   /* alpha is general nonzero */
        if (beta == 1) {
          VBRsymm_MatMult_CaABC_double(mb, n, kb, alpha,
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, ldb, 
                                  c, ldc, ind_base);
        } else if (beta == 0) {
          VBRsymm_MatMult_CaAB_double(mb, n, kb, alpha,
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, ldb, 
                                  c, ldc, ind_base);
        } else { /*  beta is general nonzero */
          VBRsymm_MatMult_CaABbC_double(mb, n, kb, alpha,
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, ldb, beta,
                                  c, ldc, ind_base);
        }
      }
      break;
    }
    break;
  default:
    printf("Invalid argument descra[1] in dvbrmm. Use 1 or 2. \n");
    break;
  } /* end of switch on descra[1] */
  break;
case 4: /* Skew Symmetric */
  if ( m != k ) {
    printf("In dvbrmm: inconsistant dimensions for a skew-symmetric matrix");
    printf("m = %d  k = %d\nExiting...\n",m,k);
    exit(-1); 
  }
  switch ( transa ) {
  case 0:
    switch ( n ) {
    case 1:
      if (alpha == 1) {
        if (beta == 1) {
          VBRskew_VecMult_CABC_double(mb, kb, 
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, 
                                  c, ind_base);
        } else if (beta == 0) {
          VBRskew_VecMult_CAB_double(mb, kb, 
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, 
                                  c, ind_base);
        } else  { /*  beta is general nonzero */
          VBRskew_VecMult_CABbC_double(mb, kb, 
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, beta,
                                  c, ind_base);
        }
      } else { /* alpha is general nonzero */
        if (beta == 1) {
          VBRskew_VecMult_CaABC_double(mb, kb, alpha,
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, 
                                  c, ind_base);
        } else if (beta == 0) {
          VBRskew_VecMult_CaAB_double(mb, kb, alpha,
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, 
                                  c, ind_base);
        } else { /*  beta is general nonzero */
          VBRskew_VecMult_CaABbC_double(mb, kb, alpha,
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, beta,
                                  c, ind_base);
        }
      }
      break;
    default:
      if (alpha == 1) {
        if (beta == 1) {
          VBRskew_MatMult_CABC_double(mb, n, kb, 
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, ldb, 
                                  c, ldc, ind_base);
        } else if (beta == 0) {
          VBRskew_MatMult_CAB_double(mb, n, kb, 
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, ldb, 
                                  c, ldc, ind_base);
        } else { /*  beta is general nonzero */
          VBRskew_MatMult_CABbC_double(mb, n, kb, 
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, ldb, beta,
                                  c, ldc, ind_base);
        }
      } else {   /* alpha is general nonzero */
        if (beta == 1) {
          VBRskew_MatMult_CaABC_double(mb, n, kb, alpha,
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, ldb, 
                                  c, ldc, ind_base);
        } else if (beta == 0) {
          VBRskew_MatMult_CaAB_double(mb, n, kb, alpha,
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, ldb, 
                                  c, ldc, ind_base);
        } else { /*  beta is general nonzero */
          VBRskew_MatMult_CaABbC_double(mb, n, kb, alpha,
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, ldb, beta,
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
          VBRskew_VecMult_CATBC_double(mb, kb, 
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, 
                                  c, ind_base);
        } else if (beta == 0) {
          VBRskew_VecMult_CATB_double(mb, kb, 
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, 
                                  c, ind_base);
        } else  { /*  beta is general nonzero */
          VBRskew_VecMult_CATBbC_double(mb, kb, 
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, beta,
                                  c, ind_base);
        }
      } else { /* alpha is general nonzero */
        if (beta == 1) {
          VBRskew_VecMult_CaATBC_double(mb, kb, alpha,
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, 
                                  c, ind_base);
        } else if (beta == 0) {
          VBRskew_VecMult_CaATB_double(mb, kb, alpha,
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, 
                                  c, ind_base);
        } else { /*  beta is general nonzero */
          VBRskew_VecMult_CaATBbC_double(mb, kb, alpha,
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, beta,
                                  c, ind_base);
        }
      }
      break;
    default:  /* n is greater than 1 -- doing Mat Mult */
      if (alpha == 1) {
        if (beta == 1) {
          VBRskew_MatMult_CATBC_double(mb, n, kb, 
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, ldb, 
                                  c, ldc, ind_base);
        } else if (beta == 0) {
          VBRskew_MatMult_CATB_double(mb, n, kb, 
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, ldb, 
                                  c, ldc, ind_base);
        } else { /*  beta is general nonzero */
          VBRskew_MatMult_CATBbC_double(mb, n, kb, 
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, ldb, beta,
                                  c, ldc, ind_base);
        }
      } else {   /* alpha is general nonzero */
        if (beta == 1) {
          VBRskew_MatMult_CaATBC_double(mb, n, kb, alpha,
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, ldb, 
                                  c, ldc, ind_base);
        } else if (beta == 0) {
          VBRskew_MatMult_CaATB_double(mb, n, kb, alpha,
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, ldb, 
                                  c, ldc, ind_base);
        } else { /*  beta is general nonzero */
          VBRskew_MatMult_CaATBbC_double(mb, n, kb, alpha,
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, ldb, beta,
                                  c, ldc, ind_base);
        }
      }
      break;
    } /* end of switch on n */
    break;
  default:
    printf("Invalid argument transa in dvbrmm. Use 0 or 1. \n");
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
          VBR_VecMult_CABC_double(mb, kb, 
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, 
                                  c, ind_base);
        } else if (beta == 0) {
          VBR_VecMult_CAB_double(mb, kb, 
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, 
                                  c, ind_base);
        } else  { /*  beta is general nonzero */
          VBR_VecMult_CABbC_double(mb, kb, 
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, beta,
                                  c, ind_base);
        }
      } else { /* alpha is general nonzero */
        if (beta == 1) {
          VBR_VecMult_CaABC_double(mb, kb, alpha,
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, 
                                  c, ind_base);
        } else if (beta == 0) {
          VBR_VecMult_CaAB_double(mb, kb, alpha,
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, 
                                  c, ind_base);
        } else { /*  beta is general nonzero */
          VBR_VecMult_CaABbC_double(mb, kb, alpha,
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, beta,
                                  c, ind_base);
        }
      }
      break;
    default:  /* n is greater than 1 -- doing Mat Mult */
      if (alpha == 1) {
        if (beta == 1) {
          VBR_MatMult_CABC_double(mb, n, kb, 
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, ldb, 
                                  c, ldc, ind_base);
        } else if (beta == 0) {
          VBR_MatMult_CAB_double(mb, n, kb, 
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, ldb, 
                                  c, ldc, ind_base);
        } else { /*  beta is general nonzero */
          VBR_MatMult_CABbC_double(mb, n, kb, 
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, ldb, beta,
                                  c, ldc, ind_base);
        }
      } else {   /* alpha is general nonzero */
        if (beta == 1) {
          VBR_MatMult_CaABC_double(mb, n, kb, alpha,
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, ldb, 
                                  c, ldc, ind_base);
        } else if (beta == 0) {
          VBR_MatMult_CaAB_double(mb, n, kb, alpha,
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, ldb, 
                                  c, ldc, ind_base);
        } else { /*  beta is general nonzero */
          VBR_MatMult_CaABbC_double(mb, n, kb, alpha,
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, ldb, beta,
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
          VBR_VecMult_CATBC_double(mb, kb, 
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, 
                                  c, ind_base);
        } else if (beta == 0) {
          VBR_VecMult_CATB_double(mb, kb, 
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, 
                                  c, ind_base);
        } else  { /*  beta is general nonzero */
          VBR_VecMult_CATBbC_double(mb, kb, 
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, beta,
                                  c, ind_base);
        }
      } else { /* alpha is general nonzero */
        if (beta == 1) {
          VBR_VecMult_CaATBC_double(mb, kb, alpha,
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, 
                                  c, ind_base);
        } else if (beta == 0) {
          VBR_VecMult_CaATB_double(mb, kb, alpha,
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, 
                                  c, ind_base);
        } else { /*  beta is general nonzero */
          VBR_VecMult_CaATBbC_double(mb, kb, alpha,
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, beta,
                                  c, ind_base);
        }
      }
      break;
    default:  /* n is greater than 1 -- doing Mat Mult */
      if (alpha == 1) {
        if (beta == 1) {
          VBR_MatMult_CATBC_double(mb, n, kb, 
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, ldb, 
                                  c, ldc, ind_base);
        } else if (beta == 0) {
          VBR_MatMult_CATB_double(mb, n, kb, 
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, ldb, 
                                  c, ldc, ind_base);
        } else { /*  beta is general nonzero */
          VBR_MatMult_CATBbC_double(mb, n, kb, 
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, ldb, beta,
                                  c, ldc, ind_base);
        }
      } else {   /* alpha is general nonzero */
        if (beta == 1) {
          VBR_MatMult_CaATBC_double(mb, n, kb, alpha,
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, ldb, 
                                  c, ldc, ind_base);
        } else if (beta == 0) {
          VBR_MatMult_CaATB_double(mb, n, kb, alpha,
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, ldb, 
                                  c, ldc, ind_base);
        } else { /*  beta is general nonzero */
          VBR_MatMult_CaATBbC_double(mb, n, kb, alpha,
                                  val, indx, bindx, rpntr, cpntr,
                                  bpntrb, bpntre, b, ldb, beta,
                                  c, ldc, ind_base);
        }
      }
      break;
    } /* end switch on n */
    break;
  default: 
    printf("Invalid argument transa in dvbrmm. Use 0 or 1. \n");
    break;
  } /* end switch on transa */
  break;
default:
    printf("Invalid argument descra[0] in dvbrmm. Use 0 - 5. \n");
    break;
}  
  
  
}
   
  
  
   
  
