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
#include "dbscmtsl.h"
#include "dbscvtsl.h"


/* Sparse BLAS Toolkit interface routine: */
void  dbscsm(
             const int transa, const int mb, const int n, 
             const int unitd, const double dv[], 
             const double alpha, const int descra[], const double val[],
             const int bindx[], const int bpntrb[], const int bpntre[],
             const int lb, const double b[], const int ldb,
             const double beta, double c[], const int ldc,
             double work[], const int lwork)
{
/* ------------ begin interface description ------------
   Toolkit interface:
   dbscsm -- block sparse column format triangular solve
  
   C <- alpha D inv(A) B + beta C    C <- alpha D inv(A') B + beta C
   C <- alpha inv(A) D B + beta C    C <- alpha inv(A') D B + beta C
   
                                      ( ' indicates matrix transpose)
  
   Arguments:
  
   int transa	Indicates how to operate with the sparse matrix
  		0 : operate with matrix
  		1 : operate with transpose matrix
  
   int mb	Number of block rows in matrix A
  
   int n	Number of columns in matrix c
  
   int unitd	Type of scaling:
                        1 : Identity matrix (argument dv[] is ignored)
                        2 : Scale on left (row scaling)
                        3 : Scale on right (column scaling)
  
   double alpha	Scalar parameter
  
   double beta 	Scalar parameter
  
   int descra[]	Descriptor argument.  Nine element integer array
  		descra[0] matrix structure
  			0 : general
  			1 : symmetric
  			2 : Hermitian
  			3 : Triangular
  			4 : Skew(Anti-Symmetric
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
  
   int *bindx	integer array of length bnnz consisting of the block column
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
  
   double *work	scratch array of length lwork.  
                lwork must be at least (n*m + lb)
  
   ------------ end interface description --------------*/
int ind_base = descra[3];
int m=mb*lb;

if (lwork < m*n + lb ){
  printf("Insufficient work space for dbscsm.\n");
  printf("   lwork must be at least (n*m + lb) = %d \n",m*n+lb);
  return;
}
if (alpha == 0.0) {
   ScaleArray_double(m, n, c, ldc, beta);
   return;
} 

switch ( descra[2] ) {
case 0:   /* Non-unit diagonal  -- not supported in the Toolkit implementation */
  printf("Unsupported functionality in dbscsm: Non-Unit diagonal blocks \n");
  printf("  descra[2] ( descra(3) in fortran) must be 1                 \n");
  break;
default:

switch (descra[0]) {
case 3:  /* Matrix MUST be triangular, of course */
  switch (transa) {
  case 0:
    switch  ( descra[1] ) {
    case 1:  /* Lower triangular, Unit diagonal */
       switch (n) {
       case 1: /* Vec Mult */
         if (alpha == 1) {
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_VecTriangSlvLU_CABC_double(mb,  val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_VecTriangSlvLU_CDABC_double(mb,  dv, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_VecTriangSlvLU_CADBC_double(mb,  dv, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_VecTriangSlvLU_CABmC_double(mb,  val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_VecTriangSlvLU_CDABmC_double(mb,  dv, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_VecTriangSlvLU_CADBmC_double(mb,  dv, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */

                BSC_VecTriangSlvLU_CAB_double(mb,  val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_VecTriangSlvLU_CDAB_double(mb,  dv, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_VecTriangSlvLU_CADB_double(mb,  dv, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else  { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                BSC_VecTriangSlvLU_CABbC_double(mb,  val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     beta, c, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_VecTriangSlvLU_CDABbC_double(mb,  dv, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     beta, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_VecTriangSlvLU_CADBbC_double(mb,  dv, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     beta, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         } else { /* alpha is general nonzero */
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_VecTriangSlvLU_CaABC_double(mb,  alpha, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_VecTriangSlvLU_CaDABC_double(mb,  dv, alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_VecTriangSlvLU_CaADBC_double(mb,  dv, alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_VecTriangSlvLU_CaABmC_double(mb,  alpha, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_VecTriangSlvLU_CaDABmC_double(mb,  dv, alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_VecTriangSlvLU_CaADBmC_double(mb,  dv, alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_VecTriangSlvLU_CaAB_double(mb,  alpha, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_VecTriangSlvLU_CaDAB_double(mb,  dv, alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_VecTriangSlvLU_CaADB_double(mb,  dv, alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, c, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                BSC_VecTriangSlvLU_CaABbC_double(mb,  alpha, val, bindx,
                                     bpntrb, bpntre, lb, 
                                     b, beta, c, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_VecTriangSlvLU_CaDABbC_double(mb,  dv, alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, beta, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_VecTriangSlvLU_CaADBbC_double(mb,  dv, alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, beta, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         }
         break;
       default: /* Mat Mult */
         if (alpha == 1) {
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_MatTriangSlvLU_CABC_double(mb, n,  val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_MatTriangSlvLU_CDABC_double(mb, n,  dv, val, bindx,
                                     bpntrb, bpntre, lb, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_MatTriangSlvLU_CADBC_double(mb, n,  dv, val, bindx,
                                     bpntrb, bpntre, lb, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_MatTriangSlvLU_CABmC_double(mb, n,  val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_MatTriangSlvLU_CDABmC_double(mb, n,  dv, val, bindx,
                                     bpntrb, bpntre, lb, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_MatTriangSlvLU_CADBmC_double(mb, n,  dv, val, bindx,
                                     bpntrb, bpntre, lb, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_MatTriangSlvLU_CAB_double(mb, n,  val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     ldb, c, ldc, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_MatTriangSlvLU_CDAB_double(mb, n,  dv, val, bindx,
                                     bpntrb, bpntre, lb, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_MatTriangSlvLU_CADB_double(mb, n,  dv, val, bindx,
                                     bpntrb, bpntre, lb, 
                                     b, ldb, c, ldc, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                BSC_MatTriangSlvLU_CABbC_double(mb, n,  val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     ldb, beta, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_MatTriangSlvLU_CDABbC_double(mb, n,  dv, val, bindx,
                                     bpntrb, bpntre, lb, 
                                     b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_MatTriangSlvLU_CADBbC_double(mb, n,  dv, val, bindx,
                                     bpntrb, bpntre, lb, 
                                     b, ldb, beta, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         } else {   /* alpha is general nonzero */
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_MatTriangSlvLU_CaABC_double(mb, n,  alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_MatTriangSlvLU_CaDABC_double(mb, n,  dv, alpha, val, bindx,
                                  
                                        bpntrb, bpntre, lb, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_MatTriangSlvLU_CaADBC_double(mb, n,  dv, alpha, val, bindx,
                                  
                                        bpntrb, bpntre, lb, b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_MatTriangSlvLU_CaABmC_double(mb, n,  alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_MatTriangSlvLU_CaDABmC_double(mb, n,  dv, alpha, val, bindx,
                                  
                                        bpntrb, bpntre, lb, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_MatTriangSlvLU_CaADBmC_double(mb, n,  dv, alpha, val, bindx,
                                  
                                        bpntrb, bpntre, lb, b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_MatTriangSlvLU_CaAB_double(mb, n,  alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, ldb, c, ldc, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_MatTriangSlvLU_CaDAB_double(mb, n,  dv, alpha, val, bindx,
                                  
                                        bpntrb, bpntre, lb, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_MatTriangSlvLU_CaADB_double(mb, n,  dv, alpha, val, bindx,
                                  
                                        bpntrb, bpntre, lb, b, ldb, c, ldc, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                BSC_MatTriangSlvLU_CaABbC_double(mb, n,  alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_MatTriangSlvLU_CaDABbC_double(mb, n,  dv, alpha, val, bindx,
                                  
                                        bpntrb, bpntre, lb, b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_MatTriangSlvLU_CaADBbC_double(mb, n,  dv, alpha, val, bindx,
                                  
                                        bpntrb, bpntre, lb, b, ldb, beta, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         }
         break;
       } /* end switch on n */
       break;
    case 2:  /* Upper triangular, Unit diagonal */
       switch (n) {
       case 1: /* Vec Mult */
         if (alpha == 1) {
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_VecTriangSlvUU_CABC_double(mb,  val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_VecTriangSlvUU_CDABC_double(mb,  dv, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_VecTriangSlvUU_CADBC_double(mb,  dv, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_VecTriangSlvUU_CABmC_double(mb,  val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_VecTriangSlvUU_CDABmC_double(mb,  dv, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_VecTriangSlvUU_CADBmC_double(mb,  dv, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_VecTriangSlvUU_CAB_double(mb,  val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_VecTriangSlvUU_CDAB_double(mb,  dv, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_VecTriangSlvUU_CADB_double(mb,  dv, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else  { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                BSC_VecTriangSlvUU_CABbC_double(mb,  val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     beta, c, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_VecTriangSlvUU_CDABbC_double(mb,  dv, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     beta, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_VecTriangSlvUU_CADBbC_double(mb,  dv, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     beta, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         } else { /* alpha is general nonzero */
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_VecTriangSlvUU_CaABC_double(mb,  alpha, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_VecTriangSlvUU_CaDABC_double(mb,  dv, alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_VecTriangSlvUU_CaADBC_double(mb,  dv, alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_VecTriangSlvUU_CaABmC_double(mb,  alpha, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_VecTriangSlvUU_CaDABmC_double(mb,  dv, alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_VecTriangSlvUU_CaADBmC_double(mb,  dv, alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_VecTriangSlvUU_CaAB_double(mb,  alpha, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_VecTriangSlvUU_CaDAB_double(mb,  dv, alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_VecTriangSlvUU_CaADB_double(mb,  dv, alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, c, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                BSC_VecTriangSlvUU_CaABbC_double(mb,  alpha, val, bindx,
                                     bpntrb, bpntre, lb, 
                                     b, beta, c, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_VecTriangSlvUU_CaDABbC_double(mb,  dv, alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, beta, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_VecTriangSlvUU_CaADBbC_double(mb,  dv, alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, beta, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         }
         break;
       default: /* Mat Mult */
         if (alpha == 1) {
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_MatTriangSlvUU_CABC_double(mb, n,  val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_MatTriangSlvUU_CDABC_double(mb, n,  dv, val, bindx,
                                     bpntrb, bpntre, lb, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_MatTriangSlvUU_CADBC_double(mb, n,  dv, val, bindx,
                                     bpntrb, bpntre, lb, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_MatTriangSlvUU_CABmC_double(mb, n,  val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_MatTriangSlvUU_CDABmC_double(mb, n,  dv, val, bindx,
                                     bpntrb, bpntre, lb, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_MatTriangSlvUU_CADBmC_double(mb, n,  dv, val, bindx,
                                     bpntrb, bpntre, lb, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_MatTriangSlvUU_CAB_double(mb, n,  val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     ldb, c, ldc, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_MatTriangSlvUU_CDAB_double(mb, n,  dv, val, bindx,
                                     bpntrb, bpntre, lb, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_MatTriangSlvUU_CADB_double(mb, n,  dv, val, bindx,
                                     bpntrb, bpntre, lb, 
                                     b, ldb, c, ldc, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                BSC_MatTriangSlvUU_CABbC_double(mb, n,  val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     ldb, beta, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_MatTriangSlvUU_CDABbC_double(mb, n,  dv, val, bindx,
                                     bpntrb, bpntre, lb, 
                                     b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_MatTriangSlvUU_CADBbC_double(mb, n,  dv, val, bindx,
                                     bpntrb, bpntre, lb, 
                                     b, ldb, beta, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         } else {   /* alpha is general nonzero */
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_MatTriangSlvUU_CaABC_double(mb, n,  alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_MatTriangSlvUU_CaDABC_double(mb, n,  dv, alpha, val, bindx,
                                  
                                        bpntrb, bpntre, lb, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_MatTriangSlvUU_CaADBC_double(mb, n,  dv, alpha, val, bindx,
                                  
                                        bpntrb, bpntre, lb, b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_MatTriangSlvUU_CaABmC_double(mb, n,  alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_MatTriangSlvUU_CaDABmC_double(mb, n,  dv, alpha, val, bindx,
                                  
                                        bpntrb, bpntre, lb, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_MatTriangSlvUU_CaADBmC_double(mb, n,  dv, alpha, val, bindx,
                                  
                                        bpntrb, bpntre, lb, b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_MatTriangSlvUU_CaAB_double(mb, n,  alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, ldb, c, ldc, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_MatTriangSlvUU_CaDAB_double(mb, n,  dv, alpha, val, bindx,
                                  
                                        bpntrb, bpntre, lb, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_MatTriangSlvUU_CaADB_double(mb, n,  dv, alpha, val, bindx,
                                  
                                        bpntrb, bpntre, lb, b, ldb, c, ldc, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                BSC_MatTriangSlvUU_CaABbC_double(mb, n,  alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_MatTriangSlvUU_CaDABbC_double(mb, n,  dv, alpha, val, bindx,
                                  
                                        bpntrb, bpntre, lb, b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_MatTriangSlvUU_CaADBbC_double(mb, n,  dv, alpha, val, bindx,
                                  
                                        bpntrb, bpntre, lb, b, ldb, beta, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         }
         break;
       } /* end switch on n */
       break;
     default: /*  Invalid descra[1] */
       printf("Check value of descra[1] (descra(2) in Fortran.\n");
       printf("Valid values of descra[1] (descra(2)) are 1 and 2.\n");
       break;
     } /* end switch on descra[1] */
     break;
  case 1:
    switch  ( descra[1] ) {
    case 1:  /* Lower triangular, Unit diagonal */
       switch (n) {
       case 1: /* Vec Mult */
         if (alpha == 1) {
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_VecTriangSlvLU_CATBC_double(mb,  val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_VecTriangSlvLU_CDATBC_double(mb,  dv, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_VecTriangSlvLU_CATDBC_double(mb,  dv, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_VecTriangSlvLU_CATBmC_double(mb,  val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_VecTriangSlvLU_CDATBmC_double(mb,  dv, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_VecTriangSlvLU_CATDBmC_double(mb,  dv, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_VecTriangSlvLU_CATB_double(mb,  val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_VecTriangSlvLU_CDATB_double(mb,  dv, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_VecTriangSlvLU_CATDB_double(mb,  dv, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else  { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                BSC_VecTriangSlvLU_CATBbC_double(mb,  val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     beta, c, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_VecTriangSlvLU_CDATBbC_double(mb,  dv, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     beta, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_VecTriangSlvLU_CATDBbC_double(mb,  dv, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     beta, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         } else { /* alpha is general nonzero */
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_VecTriangSlvLU_CaATBC_double(mb,  alpha, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_VecTriangSlvLU_CaDATBC_double(mb,  dv, alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_VecTriangSlvLU_CaATDBC_double(mb,  dv, alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_VecTriangSlvLU_CaATBmC_double(mb,  alpha, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_VecTriangSlvLU_CaDATBmC_double(mb,  dv, alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_VecTriangSlvLU_CaATDBmC_double(mb,  dv, alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_VecTriangSlvLU_CaATB_double(mb,  alpha, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_VecTriangSlvLU_CaDATB_double(mb,  dv, alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_VecTriangSlvLU_CaATDB_double(mb,  dv, alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, c, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                BSC_VecTriangSlvLU_CaATBbC_double(mb,  alpha, val, bindx,
                                     bpntrb, bpntre, lb, 
                                     b, beta, c, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_VecTriangSlvLU_CaDATBbC_double(mb,  dv, alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, beta, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_VecTriangSlvLU_CaATDBbC_double(mb,  dv, alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, beta, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         }
         break;
       default: /* Mat Mult */
         if (alpha == 1) {
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_MatTriangSlvLU_CATBC_double(mb, n,  val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_MatTriangSlvLU_CDATBC_double(mb, n,  dv, val, bindx,
                                     bpntrb, bpntre, lb, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_MatTriangSlvLU_CATDBC_double(mb, n,  dv, val, bindx,
                                     bpntrb, bpntre, lb, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_MatTriangSlvLU_CATBmC_double(mb, n,  val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_MatTriangSlvLU_CDATBmC_double(mb, n,  dv, val, bindx,
                                     bpntrb, bpntre, lb, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_MatTriangSlvLU_CATDBmC_double(mb, n,  dv, val, bindx,
                                     bpntrb, bpntre, lb, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_MatTriangSlvLU_CATB_double(mb, n,  val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     ldb, c, ldc, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_MatTriangSlvLU_CDATB_double(mb, n,  dv, val, bindx,
                                     bpntrb, bpntre, lb, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_MatTriangSlvLU_CATDB_double(mb, n,  dv, val, bindx,
                                     bpntrb, bpntre, lb, 
                                     b, ldb, c, ldc, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                BSC_MatTriangSlvLU_CATBbC_double(mb, n,  val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     ldb, beta, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_MatTriangSlvLU_CDATBbC_double(mb, n,  dv, val, bindx,
                                     bpntrb, bpntre, lb, 
                                     b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_MatTriangSlvLU_CATDBbC_double(mb, n,  dv, val, bindx,
                                     bpntrb, bpntre, lb, 
                                     b, ldb, beta, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         } else {   /* alpha is general nonzero */
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_MatTriangSlvLU_CaATBC_double(mb, n,  alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_MatTriangSlvLU_CaDATBC_double(mb, n,  dv, alpha, val, bindx,
                                  
                                        bpntrb, bpntre, lb, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_MatTriangSlvLU_CaATDBC_double(mb, n,  dv, alpha, val, bindx,
                                  
                                        bpntrb, bpntre, lb, b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_MatTriangSlvLU_CaATBmC_double(mb, n,  alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_MatTriangSlvLU_CaDATBmC_double(mb, n,  dv, alpha, val, bindx,
                                  
                                        bpntrb, bpntre, lb, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_MatTriangSlvLU_CaATDBmC_double(mb, n,  dv, alpha, val, bindx,
                                  
                                        bpntrb, bpntre, lb, b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_MatTriangSlvLU_CaATB_double(mb, n,  alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, ldb, c, ldc, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_MatTriangSlvLU_CaDATB_double(mb, n,  dv, alpha, val, bindx,
                                  
                                        bpntrb, bpntre, lb, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_MatTriangSlvLU_CaATDB_double(mb, n,  dv, alpha, val, bindx,
                                  
                                        bpntrb, bpntre, lb, b, ldb, c, ldc, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                BSC_MatTriangSlvLU_CaATBbC_double(mb, n,  alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_MatTriangSlvLU_CaDATBbC_double(mb, n,  dv, alpha, val, bindx,
                                  
                                        bpntrb, bpntre, lb, b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_MatTriangSlvLU_CaATDBbC_double(mb, n,  dv, alpha, val, bindx,
                                  
                                        bpntrb, bpntre, lb, b, ldb, beta, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         }
         break;
       } /* end switch on n */
       break;
    case 2:  /* Upper triangular, Unit diagonal */
       switch (n) {
       case 1: /* Vec Mult */
         if (alpha == 1) {
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_VecTriangSlvUU_CATBC_double(mb,  val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_VecTriangSlvUU_CDATBC_double(mb,  dv, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_VecTriangSlvUU_CATDBC_double(mb,  dv, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_VecTriangSlvUU_CATBmC_double(mb,  val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_VecTriangSlvUU_CDATBmC_double(mb,  dv, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_VecTriangSlvUU_CATDBmC_double(mb,  dv, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_VecTriangSlvUU_CATB_double(mb,  val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_VecTriangSlvUU_CDATB_double(mb,  dv, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_VecTriangSlvUU_CATDB_double(mb,  dv, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else  { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                BSC_VecTriangSlvUU_CATBbC_double(mb,  val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     beta, c, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_VecTriangSlvUU_CDATBbC_double(mb,  dv, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     beta, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_VecTriangSlvUU_CATDBbC_double(mb,  dv, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     beta, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         } else { /* alpha is general nonzero */
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_VecTriangSlvUU_CaATBC_double(mb,  alpha, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_VecTriangSlvUU_CaDATBC_double(mb,  dv, alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_VecTriangSlvUU_CaATDBC_double(mb,  dv, alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_VecTriangSlvUU_CaATBmC_double(mb,  alpha, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_VecTriangSlvUU_CaDATBmC_double(mb,  dv, alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_VecTriangSlvUU_CaATDBmC_double(mb,  dv, alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_VecTriangSlvUU_CaATB_double(mb,  alpha, val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     c, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_VecTriangSlvUU_CaDATB_double(mb,  dv, alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_VecTriangSlvUU_CaATDB_double(mb,  dv, alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, c, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                BSC_VecTriangSlvUU_CaATBbC_double(mb,  alpha, val, bindx,
                                     bpntrb, bpntre, lb, 
                                     b, beta, c, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_VecTriangSlvUU_CaDATBbC_double(mb,  dv, alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, beta, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_VecTriangSlvUU_CaATDBbC_double(mb,  dv, alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, beta, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         }
         break;
       default: /* Mat Mult */
         if (alpha == 1) {
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_MatTriangSlvUU_CATBC_double(mb, n,  val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_MatTriangSlvUU_CDATBC_double(mb, n,  dv, val, bindx,
                                     bpntrb, bpntre, lb, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_MatTriangSlvUU_CATDBC_double(mb, n,  dv, val, bindx,
                                     bpntrb, bpntre, lb, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_MatTriangSlvUU_CATBmC_double(mb, n,  val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_MatTriangSlvUU_CDATBmC_double(mb, n,  dv, val, bindx,
                                     bpntrb, bpntre, lb, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_MatTriangSlvUU_CATDBmC_double(mb, n,  dv, val, bindx,
                                     bpntrb, bpntre, lb, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_MatTriangSlvUU_CATB_double(mb, n,  val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     ldb, c, ldc, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_MatTriangSlvUU_CDATB_double(mb, n,  dv, val, bindx,
                                     bpntrb, bpntre, lb, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_MatTriangSlvUU_CATDB_double(mb, n,  dv, val, bindx,
                                     bpntrb, bpntre, lb, 
                                     b, ldb, c, ldc, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                BSC_MatTriangSlvUU_CATBbC_double(mb, n,  val, bindx,
                                     bpntrb, bpntre, lb, b, 
                                     ldb, beta, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_MatTriangSlvUU_CDATBbC_double(mb, n,  dv, val, bindx,
                                     bpntrb, bpntre, lb, 
                                     b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_MatTriangSlvUU_CATDBbC_double(mb, n,  dv, val, bindx,
                                     bpntrb, bpntre, lb, 
                                     b, ldb, beta, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         } else {   /* alpha is general nonzero */
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_MatTriangSlvUU_CaATBC_double(mb, n,  alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_MatTriangSlvUU_CaDATBC_double(mb, n,  dv, alpha, val, bindx,
                                  
                                        bpntrb, bpntre, lb, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_MatTriangSlvUU_CaATDBC_double(mb, n,  dv, alpha, val, bindx,
                                  
                                        bpntrb, bpntre, lb, b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_MatTriangSlvUU_CaATBmC_double(mb, n,  alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_MatTriangSlvUU_CaDATBmC_double(mb, n,  dv, alpha, val, bindx,
                                  
                                        bpntrb, bpntre, lb, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_MatTriangSlvUU_CaATDBmC_double(mb, n,  dv, alpha, val, bindx,
                                  
                                        bpntrb, bpntre, lb, b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                BSC_MatTriangSlvUU_CaATB_double(mb, n,  alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, ldb, c, ldc, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_MatTriangSlvUU_CaDATB_double(mb, n,  dv, alpha, val, bindx,
                                  
                                        bpntrb, bpntre, lb, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_MatTriangSlvUU_CaATDB_double(mb, n,  dv, alpha, val, bindx,
                                  
                                        bpntrb, bpntre, lb, b, ldb, c, ldc, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                BSC_MatTriangSlvUU_CaATBbC_double(mb, n,  alpha, val, bindx,
                                     bpntrb, 
                                     bpntre, lb, b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                BSC_MatTriangSlvUU_CaDATBbC_double(mb, n,  dv, alpha, val, bindx,
                                  
                                        bpntrb, bpntre, lb, b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                BSC_MatTriangSlvUU_CaATDBbC_double(mb, n,  dv, alpha, val, bindx,
                                  
                                        bpntrb, bpntre, lb, b, ldb, beta, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dbscsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         }
         break;
       } /* end switch on n */
       break;
     default: /*  Invalid descra[1] */
       printf("Check value of descra[1] (descra(2) in Fortran.\n");
       printf("Valid values of descra[1] (descra(2)) are 1 and 2.\n");
       break;
     } /* end switch on descra[1] */
     break;
  default: 
     printf("Invalid argument transa in dbscsm. Use 0 or 1. \n");
     break;
  } /* end switch on transa */
  break;
default:
  printf("Invalid argument descra[0] in dbscsm. Must be 3 (triangular matrix). \n");
  break;
} /* end switch on descra[0] */

} /* end switch on descra[2] */
  
}
   
