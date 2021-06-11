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

#include <stdio.h>
#include <stdlib.h>
#include "spblas.h"


double resid (int, double *, double *);

int main(int argc, char *argv[])
{
/* Initialize the test matrices (one lower triangular, one upper) */
double dv[]={.5, 0, 0,  0, -.5, 0,  0, 0, .5, 
             1, 0, 0,  0, 1, 0,  0, 0, 1, 
             1, 0, 0,  0, 1, 0,  0, 0, 1, 
             .25, 0, 0,  0, .25, 0,  0, 0, -.25};
double diag[]={.5,-.5,.5,1,1,1,1,1,1,.25,.25,-.25};
double a[]= {1,  0, 0,  0, 1, 0,  0, 0, 1,    /* lower triangular */
             1,  4, 7,  2, 5, 8,  3, 6, 9, 
            10,-13,16, 11,14,17, 12,15,18,
             1,  0, 0,  0, 1, 0,  0, 0, 1, 
             1,  0, 0,  0, 1, 0,  0, 0, 1, 
            -1,  4, 7,  2,-5, 8,  3, 6,-9,
             1,  0, 0,  0, 1, 0,  0, 0, 1};
int bindx[]={1,2,4,2,3,4,4};
int bpntrb[]={1,4,5,7};
int bpntre[]={4,5,7,8};
double a2[]= {1,  0, 0,  0, 1, 0,  0, 0, 1,    /* upper triangular */ 
             1,  2, 3,  4, 5, 6,  7, 8, 9, 
             1,  0, 0,  0, 1, 0,  0, 0, 1, 
             1,  0, 0,  0, 1, 0,  0, 0, 1, 
             10,11,12,-13,14,15, 16,17,18,
            -1,  2, 3,  4,-5, 6,  7, 8,-9, 
             1,  0, 0,  0, 1, 0,  0, 0, 1};
int bindx2[]={1,1,2,3,1,3,4};
int bpntrb2[]={1,2,4,5};
int bpntre2[]={2,4,5,8};
double b[]={1,2,3,4,5,6,7,8,9,10,11,12,
            1,2,3,4,5,6,7,8,9,10,11,12};
double c[]={1,2,3,4,5,6,7,8,9,10,11,12,
            1,2,3,4,5,6,7,8,9,10,11,12};
double d[]={1,2,3,4,5,6,7,8,9,10,11,12,
            1,2,3,4,5,6,7,8,9,10,11,12};
double check[]={1,2,3,4,5,6,7,8,9,10,11,12,
            1,2,3,4,5,6,7,8,9,10,11,12};
int mb=4, kb=4, lb=3, m=12, ldb=12, ldc=12;
int i,j;
int transa, n, unitd, lwork;
int descra[9];
int errcount=0;
double alpha;
double beta;
double zero=0.0;
double error;
double tolerance=.00001;
double *work;
lwork = 33;
work = (double *) malloc(lwork*sizeof(double));

/* Get input: alpha and beta */

if (argc != 3 ) {
   printf("Usage:   %s alpha beta \n", argv[0]);
   exit(1);
}

alpha = (double) atof(argv[1]);
beta = (double) atof(argv[2]);

descra[0] = 3;
descra[2] = 1;
descra[3] = 1;
descra[4] = 1;


printf("-----------------------------------------------------\n");
printf("  alpha = %e, beta = %e  \n",alpha, beta);
printf("-----------------------------------------------------\n");
for (n=1;n!=3;n++) {                      /* loop on columns in C */
  printf("*** n = %d ***\n",n);           /* (test vector and matrix routines) */
  for (transa=0;transa!=2;transa++) {     /* test non-transpose and transpose  */
    printf("   << transa = %d >>\n",transa);
    for (unitd=1;unitd!=4;unitd++) {      /* test identity, left and right scaling */
      printf("      ++ unitd = %d ++\n",unitd);

      descra[1] = 1;                      /* lower triangular matrix */
      printf("          -- lower triangular --\n");

      for (i=0;i!=m;i++)                  /* Initialize c */
        for (j=0;j!=n;j++)
          c[j*m+i] = i+1;

/* Call triangular solve with lower triangular matrix     */

      dbscsm( transa, mb, n, unitd, dv, alpha, descra, a,
              bindx, bpntrb, bpntre, lb, b, ldb,
              beta, c, ldc, work, lwork);

/* Backtrack from solution using matrix multiply; after   */
/* calculation, "check" should match "b"                  */

      for (i=0;i!=n*m;i++) 
        d[i] = c[i] - beta * b[i];
  
      if ( alpha != 0 ) {
        if ( unitd == 2 ) 
          for (i=0;i!=m;i++)
            for (j=0;j!=n;j++)
              d[j*m+i] /= diag[i];
  

        dbscmm( transa, mb, n, kb, 1/alpha, descra, a,
                bindx, bpntrb, bpntre, lb, d, ldb,
                zero, check, ldc, work, lwork);

        if ( unitd == 3 )
          for (i=0;i!=m;i++)
            for (j=0;j!=n;j++)
              check[j*m+i] /= diag[i];
  
        error = resid(n*m, check, b);
      } else {
        error = 0;
        for (i=0;i<n*m;i++) {
           check[i] = d[i];
           error += abs(d[i]);
        }
        error /= n*m;
      }
      if ( error >= tolerance ){
         errcount++;
         printf("Error for lower triangular solve with ");
         printf("n = %d, transa = %d, unitd = %d.\n",n,transa,unitd);
         printf("Residual: %10.6f \n",error);
         for (i=0;i!=n*m;i++) 
           printf("%6.2f  %6.2f\n",c[i], check[i]);
      }

      descra[1] = 2;                      /* upper triangular matrix */
      printf("          -- upper triangular --\n");
 
      for (i=0;i!=m;i++)
        for (j=0;j!=n;j++)
          c[j*m+i] = i+1;

/* Call triangular solve with upper triangular matrix: */

      dbscsm( transa, mb, n, unitd, dv, alpha, descra, a2,
              bindx2, bpntrb2, bpntre2, lb, b, ldb,
              beta, c, ldc, work, lwork);

/* Backtrack from solution using matrix multiply; after   */
/* calculation, "check" should match "b"                  */

      for (i=0;i!=n*m;i++) 
        d[i] = c[i] - beta * b[i];
  
      if ( alpha != 0 ) {
        if ( unitd == 2 ) 
          for (i=0;i!=m;i++)
            for (j=0;j!=n;j++)
              d[j*m+i] /= diag[i];
  

        dbscmm( transa, mb, n, kb, 1/alpha, descra, a2,
                bindx2, bpntrb2, bpntre2, lb, d, ldb,
                zero, check, ldc, work, lwork);

        if ( unitd == 3 )
          for (i=0;i!=m;i++)
            for (j=0;j!=n;j++)
              check[j*m+i] /= diag[i];

        error = resid(n*m, check, b);
      } else {
        error = 0;
        for (i=0;i<n*m;i++) {
          error += abs(d[i]);
          check[i] = d[i];
        }
        error /= n*m;
      }
      if ( error >= tolerance ){
         errcount++;
         printf("Error for upper triangular solve with ");
         printf("n = %d, transa = %d, unitd = %d.\n",n,transa,unitd);
         printf("Residual: %10.6f \n",error);
         for (i=0;i!=n*m;i++) 
          printf("%6.2f  %6.2f\n",c[i], check[i]);
      }

    } /* close loop on unitd */
  } /* close loop on transa */
} /* close loop on n */

if ( errcount > 0 ) 
  printf("%d errors in dtbscsm run for alpha = %e, beta = %e\n",errcount,alpha, beta);

return errcount;
} /* end main */

double resid(int m, double *x1, double *x2) 
{

   double norm;
   int i;

   norm = 0.0;
 
   for (i=0;i<m;i++) norm += abs(x1[i] - x2[i]);

   if ( m == 0 ) {
     return norm;
   } else {
     return norm/m;
   }
 
}

