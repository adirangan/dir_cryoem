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
#include "dvbrmtsl.h"
#include "dvbrvtsl.h"


/* Sparse BLAS Toolkit interface routine: */
void  dvbrsm(
             const int transa, const int mb, const int n, 
             const int unitd, const double dv[], 
             const double alpha, const int descra[], const double val[],
             const int indx[], const int bindx[], const int rpntr[],
             const int cpntr[], const int bpntrb[], const int bpntre[],
             const double b[], const int ldb,
             const double beta, double c[], const int ldc,
             double work[], const int lwork)
{
/* ------------ begin interface description ------------
   Toolkit interface:
   dvbrsm -- variable block sparse row format triangular solve
  
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
  
   int *cpntr	integer array of length mb+1 such that cpntr(j)-cpntr(1)
  		is the column index of the first point column in the j-th
  		block column. cpntr(mb+1) is set to k+cpntr(1).
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
  
   double *work	scratch array of length lwork.  
                lwork must be at least (n*m + max_block_size)
  
   ------------ end interface description --------------*/
int ind_base = descra[3];
int m=rpntr[mb]-rpntr[0];

if (lwork < m*n + m/mb ){
  printf("Warning: possible insufficient work space for dvbrsm.\n");
  printf("   lwork must be at least (n*m + max_block_size) \n");
}

if (alpha == 0.0) {
   ScaleArray_double(m, n, c, ldc, beta);
   return;
} 

switch ( descra[2] ) {
case 0:   /* Non-unit diagonal  -- not supported in the Toolkit implementation */
  printf("Unsupported functionality in dvbrsm: Non-Unit diagonal blocks \n");
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
                VBR_VecTriangSlvLU_CABC_double(mb,  val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_VecTriangSlvLU_CDABC_double(mb,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_VecTriangSlvLU_CADBC_double(mb,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_VecTriangSlvLU_CABmC_double(mb,  val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_VecTriangSlvLU_CDABmC_double(mb,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_VecTriangSlvLU_CADBmC_double(mb,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */

                VBR_VecTriangSlvLU_CAB_double(mb,  val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_VecTriangSlvLU_CDAB_double(mb,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_VecTriangSlvLU_CADB_double(mb,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else  { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                VBR_VecTriangSlvLU_CABbC_double(mb,  val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     beta, c, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_VecTriangSlvLU_CDABbC_double(mb,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     beta, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_VecTriangSlvLU_CADBbC_double(mb,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     beta, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         } else { /* alpha is general nonzero */
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_VecTriangSlvLU_CaABC_double(mb,  alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_VecTriangSlvLU_CaDABC_double(mb,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_VecTriangSlvLU_CaADBC_double(mb,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_VecTriangSlvLU_CaABmC_double(mb,  alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_VecTriangSlvLU_CaDABmC_double(mb,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_VecTriangSlvLU_CaADBmC_double(mb,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_VecTriangSlvLU_CaAB_double(mb,  alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_VecTriangSlvLU_CaDAB_double(mb,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_VecTriangSlvLU_CaADB_double(mb,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, c, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                VBR_VecTriangSlvLU_CaABbC_double(mb,  alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, 
                                     b, beta, c, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_VecTriangSlvLU_CaDABbC_double(mb,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, beta, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_VecTriangSlvLU_CaADBbC_double(mb,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, beta, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
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
                VBR_MatTriangSlvLU_CABC_double(mb, n,  val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_MatTriangSlvLU_CDABC_double(mb, n,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_MatTriangSlvLU_CADBC_double(mb, n,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_MatTriangSlvLU_CABmC_double(mb, n,  val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_MatTriangSlvLU_CDABmC_double(mb, n,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_MatTriangSlvLU_CADBmC_double(mb, n,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_MatTriangSlvLU_CAB_double(mb, n,  val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     ldb, c, ldc, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_MatTriangSlvLU_CDAB_double(mb, n,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_MatTriangSlvLU_CADB_double(mb, n,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, 
                                     b, ldb, c, ldc, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                VBR_MatTriangSlvLU_CABbC_double(mb, n,  val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     ldb, beta, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_MatTriangSlvLU_CDABbC_double(mb, n,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, 
                                     b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_MatTriangSlvLU_CADBbC_double(mb, n,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, 
                                     b, ldb, beta, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         } else {   /* alpha is general nonzero */
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_MatTriangSlvLU_CaABC_double(mb, n,  alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_MatTriangSlvLU_CaDABC_double(mb, n,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                  
                                        bpntrb, bpntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_MatTriangSlvLU_CaADBC_double(mb, n,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                  
                                        bpntrb, bpntre, b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_MatTriangSlvLU_CaABmC_double(mb, n,  alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_MatTriangSlvLU_CaDABmC_double(mb, n,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                  
                                        bpntrb, bpntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_MatTriangSlvLU_CaADBmC_double(mb, n,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                  
                                        bpntrb, bpntre, b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_MatTriangSlvLU_CaAB_double(mb, n,  alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, ldb, c, ldc, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_MatTriangSlvLU_CaDAB_double(mb, n,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                  
                                        bpntrb, bpntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_MatTriangSlvLU_CaADB_double(mb, n,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                  
                                        bpntrb, bpntre, b, ldb, c, ldc, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                VBR_MatTriangSlvLU_CaABbC_double(mb, n,  alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_MatTriangSlvLU_CaDABbC_double(mb, n,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                  
                                        bpntrb, bpntre, b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_MatTriangSlvLU_CaADBbC_double(mb, n,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                  
                                        bpntrb, bpntre, b, ldb, beta, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
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
                VBR_VecTriangSlvUU_CABC_double(mb,  val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_VecTriangSlvUU_CDABC_double(mb,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_VecTriangSlvUU_CADBC_double(mb,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_VecTriangSlvUU_CABmC_double(mb,  val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_VecTriangSlvUU_CDABmC_double(mb,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_VecTriangSlvUU_CADBmC_double(mb,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_VecTriangSlvUU_CAB_double(mb,  val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_VecTriangSlvUU_CDAB_double(mb,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_VecTriangSlvUU_CADB_double(mb,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else  { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                VBR_VecTriangSlvUU_CABbC_double(mb,  val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     beta, c, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_VecTriangSlvUU_CDABbC_double(mb,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     beta, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_VecTriangSlvUU_CADBbC_double(mb,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     beta, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         } else { /* alpha is general nonzero */
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_VecTriangSlvUU_CaABC_double(mb,  alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_VecTriangSlvUU_CaDABC_double(mb,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_VecTriangSlvUU_CaADBC_double(mb,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_VecTriangSlvUU_CaABmC_double(mb,  alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_VecTriangSlvUU_CaDABmC_double(mb,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_VecTriangSlvUU_CaADBmC_double(mb,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_VecTriangSlvUU_CaAB_double(mb,  alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_VecTriangSlvUU_CaDAB_double(mb,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_VecTriangSlvUU_CaADB_double(mb,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, c, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                VBR_VecTriangSlvUU_CaABbC_double(mb,  alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, 
                                     b, beta, c, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_VecTriangSlvUU_CaDABbC_double(mb,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, beta, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_VecTriangSlvUU_CaADBbC_double(mb,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, beta, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
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
                VBR_MatTriangSlvUU_CABC_double(mb, n,  val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_MatTriangSlvUU_CDABC_double(mb, n,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_MatTriangSlvUU_CADBC_double(mb, n,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_MatTriangSlvUU_CABmC_double(mb, n,  val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_MatTriangSlvUU_CDABmC_double(mb, n,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_MatTriangSlvUU_CADBmC_double(mb, n,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_MatTriangSlvUU_CAB_double(mb, n,  val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     ldb, c, ldc, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_MatTriangSlvUU_CDAB_double(mb, n,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_MatTriangSlvUU_CADB_double(mb, n,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, 
                                     b, ldb, c, ldc, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                VBR_MatTriangSlvUU_CABbC_double(mb, n,  val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     ldb, beta, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_MatTriangSlvUU_CDABbC_double(mb, n,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, 
                                     b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_MatTriangSlvUU_CADBbC_double(mb, n,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, 
                                     b, ldb, beta, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         } else {   /* alpha is general nonzero */
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_MatTriangSlvUU_CaABC_double(mb, n,  alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_MatTriangSlvUU_CaDABC_double(mb, n,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                  
                                        bpntrb, bpntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_MatTriangSlvUU_CaADBC_double(mb, n,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                  
                                        bpntrb, bpntre, b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_MatTriangSlvUU_CaABmC_double(mb, n,  alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_MatTriangSlvUU_CaDABmC_double(mb, n,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                  
                                        bpntrb, bpntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_MatTriangSlvUU_CaADBmC_double(mb, n,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                  
                                        bpntrb, bpntre, b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_MatTriangSlvUU_CaAB_double(mb, n,  alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, ldb, c, ldc, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_MatTriangSlvUU_CaDAB_double(mb, n,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                  
                                        bpntrb, bpntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_MatTriangSlvUU_CaADB_double(mb, n,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                  
                                        bpntrb, bpntre, b, ldb, c, ldc, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                VBR_MatTriangSlvUU_CaABbC_double(mb, n,  alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_MatTriangSlvUU_CaDABbC_double(mb, n,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                  
                                        bpntrb, bpntre, b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_MatTriangSlvUU_CaADBbC_double(mb, n,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                  
                                        bpntrb, bpntre, b, ldb, beta, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
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
                VBR_VecTriangSlvLU_CATBC_double(mb,  val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_VecTriangSlvLU_CDATBC_double(mb,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_VecTriangSlvLU_CATDBC_double(mb,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_VecTriangSlvLU_CATBmC_double(mb,  val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_VecTriangSlvLU_CDATBmC_double(mb,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_VecTriangSlvLU_CATDBmC_double(mb,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_VecTriangSlvLU_CATB_double(mb,  val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_VecTriangSlvLU_CDATB_double(mb,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_VecTriangSlvLU_CATDB_double(mb,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else  { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                VBR_VecTriangSlvLU_CATBbC_double(mb,  val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     beta, c, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_VecTriangSlvLU_CDATBbC_double(mb,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     beta, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_VecTriangSlvLU_CATDBbC_double(mb,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     beta, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         } else { /* alpha is general nonzero */
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_VecTriangSlvLU_CaATBC_double(mb,  alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_VecTriangSlvLU_CaDATBC_double(mb,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_VecTriangSlvLU_CaATDBC_double(mb,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_VecTriangSlvLU_CaATBmC_double(mb,  alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_VecTriangSlvLU_CaDATBmC_double(mb,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_VecTriangSlvLU_CaATDBmC_double(mb,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_VecTriangSlvLU_CaATB_double(mb,  alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_VecTriangSlvLU_CaDATB_double(mb,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_VecTriangSlvLU_CaATDB_double(mb,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, c, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                VBR_VecTriangSlvLU_CaATBbC_double(mb,  alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, 
                                     b, beta, c, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_VecTriangSlvLU_CaDATBbC_double(mb,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, beta, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_VecTriangSlvLU_CaATDBbC_double(mb,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, beta, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
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
                VBR_MatTriangSlvLU_CATBC_double(mb, n,  val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_MatTriangSlvLU_CDATBC_double(mb, n,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_MatTriangSlvLU_CATDBC_double(mb, n,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_MatTriangSlvLU_CATBmC_double(mb, n,  val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_MatTriangSlvLU_CDATBmC_double(mb, n,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_MatTriangSlvLU_CATDBmC_double(mb, n,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_MatTriangSlvLU_CATB_double(mb, n,  val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     ldb, c, ldc, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_MatTriangSlvLU_CDATB_double(mb, n,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_MatTriangSlvLU_CATDB_double(mb, n,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, 
                                     b, ldb, c, ldc, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                VBR_MatTriangSlvLU_CATBbC_double(mb, n,  val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     ldb, beta, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_MatTriangSlvLU_CDATBbC_double(mb, n,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, 
                                     b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_MatTriangSlvLU_CATDBbC_double(mb, n,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, 
                                     b, ldb, beta, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         } else {   /* alpha is general nonzero */
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_MatTriangSlvLU_CaATBC_double(mb, n,  alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_MatTriangSlvLU_CaDATBC_double(mb, n,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                  
                                        bpntrb, bpntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_MatTriangSlvLU_CaATDBC_double(mb, n,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                  
                                        bpntrb, bpntre, b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_MatTriangSlvLU_CaATBmC_double(mb, n,  alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_MatTriangSlvLU_CaDATBmC_double(mb, n,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                  
                                        bpntrb, bpntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_MatTriangSlvLU_CaATDBmC_double(mb, n,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                  
                                        bpntrb, bpntre, b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_MatTriangSlvLU_CaATB_double(mb, n,  alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, ldb, c, ldc, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_MatTriangSlvLU_CaDATB_double(mb, n,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                  
                                        bpntrb, bpntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_MatTriangSlvLU_CaATDB_double(mb, n,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                  
                                        bpntrb, bpntre, b, ldb, c, ldc, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                VBR_MatTriangSlvLU_CaATBbC_double(mb, n,  alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_MatTriangSlvLU_CaDATBbC_double(mb, n,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                  
                                        bpntrb, bpntre, b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_MatTriangSlvLU_CaATDBbC_double(mb, n,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                  
                                        bpntrb, bpntre, b, ldb, beta, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
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
                VBR_VecTriangSlvUU_CATBC_double(mb,  val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_VecTriangSlvUU_CDATBC_double(mb,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_VecTriangSlvUU_CATDBC_double(mb,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_VecTriangSlvUU_CATBmC_double(mb,  val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_VecTriangSlvUU_CDATBmC_double(mb,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_VecTriangSlvUU_CATDBmC_double(mb,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_VecTriangSlvUU_CATB_double(mb,  val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_VecTriangSlvUU_CDATB_double(mb,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_VecTriangSlvUU_CATDB_double(mb,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else  { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                VBR_VecTriangSlvUU_CATBbC_double(mb,  val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     beta, c, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_VecTriangSlvUU_CDATBbC_double(mb,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     beta, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_VecTriangSlvUU_CATDBbC_double(mb,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     beta, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         } else { /* alpha is general nonzero */
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_VecTriangSlvUU_CaATBC_double(mb,  alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_VecTriangSlvUU_CaDATBC_double(mb,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_VecTriangSlvUU_CaATDBC_double(mb,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_VecTriangSlvUU_CaATBmC_double(mb,  alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_VecTriangSlvUU_CaDATBmC_double(mb,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_VecTriangSlvUU_CaATDBmC_double(mb,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_VecTriangSlvUU_CaATB_double(mb,  alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     c, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_VecTriangSlvUU_CaDATB_double(mb,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_VecTriangSlvUU_CaATDB_double(mb,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, c, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                VBR_VecTriangSlvUU_CaATBbC_double(mb,  alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, 
                                     b, beta, c, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_VecTriangSlvUU_CaDATBbC_double(mb,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, beta, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_VecTriangSlvUU_CaATDBbC_double(mb,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, beta, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
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
                VBR_MatTriangSlvUU_CATBC_double(mb, n,  val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_MatTriangSlvUU_CDATBC_double(mb, n,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_MatTriangSlvUU_CATDBC_double(mb, n,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_MatTriangSlvUU_CATBmC_double(mb, n,  val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_MatTriangSlvUU_CDATBmC_double(mb, n,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_MatTriangSlvUU_CATDBmC_double(mb, n,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_MatTriangSlvUU_CATB_double(mb, n,  val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     ldb, c, ldc, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_MatTriangSlvUU_CDATB_double(mb, n,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_MatTriangSlvUU_CATDB_double(mb, n,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, 
                                     b, ldb, c, ldc, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                VBR_MatTriangSlvUU_CATBbC_double(mb, n,  val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, b, 
                                     ldb, beta, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_MatTriangSlvUU_CDATBbC_double(mb, n,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, 
                                     b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_MatTriangSlvUU_CATDBbC_double(mb, n,  dv, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, bpntre, 
                                     b, ldb, beta, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         } else {   /* alpha is general nonzero */
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_MatTriangSlvUU_CaATBC_double(mb, n,  alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_MatTriangSlvUU_CaDATBC_double(mb, n,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                  
                                        bpntrb, bpntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_MatTriangSlvUU_CaATDBC_double(mb, n,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                  
                                        bpntrb, bpntre, b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_MatTriangSlvUU_CaATBmC_double(mb, n,  alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_MatTriangSlvUU_CaDATBmC_double(mb, n,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                  
                                        bpntrb, bpntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_MatTriangSlvUU_CaATDBmC_double(mb, n,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                  
                                        bpntrb, bpntre, b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                VBR_MatTriangSlvUU_CaATB_double(mb, n,  alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, ldb, c, ldc, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_MatTriangSlvUU_CaDATB_double(mb, n,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                  
                                        bpntrb, bpntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_MatTriangSlvUU_CaATDB_double(mb, n,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                  
                                        bpntrb, bpntre, b, ldb, c, ldc, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                VBR_MatTriangSlvUU_CaATBbC_double(mb, n,  alpha, val, indx, bindx, rpntr, cpntr,
                                     bpntrb, 
                                     bpntre, b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                VBR_MatTriangSlvUU_CaDATBbC_double(mb, n,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                  
                                        bpntrb, bpntre, b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                VBR_MatTriangSlvUU_CaATDBbC_double(mb, n,  dv, alpha, val, indx, bindx, rpntr, cpntr,
                                  
                                        bpntrb, bpntre, b, ldb, beta, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dvbrsm. Use 1-3. \n");
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
     printf("Invalid argument transa in dvbrsm. Use 0 or 1. \n");
     break;
  } /* end switch on transa */
  break;
default:
  printf("Invalid argument descra[0] in dvbrsm. Must be 3 (triangular matrix). \n");
  break;
} /* end switch on descra[0] */

} /* end switch on descra[2] */
  
}
   
