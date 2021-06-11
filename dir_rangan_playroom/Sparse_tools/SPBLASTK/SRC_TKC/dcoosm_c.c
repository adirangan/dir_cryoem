#include <stdlib.h>
#include <stdio.h>
#include "spblas.h"

void dcoosm(const int transa, const int m, const int n,
   const int unitd, const double dv[], const double alpha,
   const int descra[], const double val[],
   const int indx[], const int jndx[], const int nnz, const double b[], const int ldb,
   const double beta, double c[], const int ldc,
   double work[], const int lwork)
{
/* ------------ begin interface description ------------
   Toolkit interface:
   dcoosm -- compressed coordinate format triangular solve
             (WARNING: creates sparse column matrix and then calls solve)

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
  
   int *indx	integer array of length nnz consisting of the row
   		indices of the entries of A.
  
   int *jndx	integer array of length nnz consisting of the column
   		indices of the entries of A.
  
   int nnz	number of nonzero entries in A.
  
   double *b	rectangular array with first dimension ldb
  
   double *c	rectangular array with first dimension ldc
  
   double *work	scratch array of length lwork. 
                lwork should be at least (n*m)
  
   ------------ end interface description --------------*/
  int i,j;
  int *tally; 
  int *colptr; 
  int *rowind;
  double *cval;
  int ind_base = descra[3];

  if (lwork < m*n ) {
    printf("Insufficient work space for dcoosm.\n");
    printf("   lwork must be at least (n*m) = %d \n",m*n);
    return;
  }

  if (descra[0] != 3) {
    fprintf(stderr,"Must have triangular matrix for dcoosm. (descra[0] == 3)\n");
    exit(1);
  }

/* First, convert storage vectors to CompCol format: */
/* (adapted from CompCol_Mat_double(Coord_Mat_double) constructor) */
        tally  = (int *)malloc((m+1)*sizeof(int));
        for (i=0;i<m+1;i++) tally[i] = 0;
        colptr = (int *)malloc((m+1)*sizeof(int));
        rowind = (int *)malloc((nnz)*sizeof(int));
        cval = (double *)malloc((nnz)*sizeof(double));
/*      The following pointers will be incremented/decremented to create         */
/*      the appropriate vector offset if ind_base != 0 (e.g. 1-based vectors)    */
/*      This elminates the need to adjust element by element with ind_base.      */
        tally=tally-ind_base;
        rowind=rowind-ind_base;
        cval=cval-ind_base;
/*      First pass through nonzeros.  Tally entries in each column.         */
/*      And calculate colptr array. */
        for (i=0;i<nnz;i++) tally[jndx[i]]++;
        colptr[0] = ind_base;
        tally=tally+ind_base;
        for (j=0;j<m;j++) colptr[j+1] = colptr[j]+tally[j];
/*      Make copy of colptr for use in second pass. */
        for (j=0;j<m;j++) tally[j] = colptr[j];
        tally=tally-ind_base;
/*      Second pass through nonzeros.   Fill in index and value entries. */
        for (i=0;i<nnz;i++)
        {
           cval[tally[jndx[i]]] = val[i];
           rowind[tally[jndx[i]]] = indx[i];
           tally[jndx[i]]++;
        }
        tally=tally+ind_base;
        rowind=rowind+ind_base;
        cval=cval+ind_base;

/* Now, call cscsm routine using newly created CompCol storage vectors: */

  dcscsm( transa, m, n, unitd, dv, alpha,
          descra, cval, rowind, colptr, &colptr[1], b, ldb,
          beta, c, ldc, work, lwork);
 
  free(cval); free(colptr); free(rowind); free(tally);
}

