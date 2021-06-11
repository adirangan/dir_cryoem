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

double a[]= {1,  0, 0,  0, 1, 0,  0, 0, 1,           /* All of matrix A */
             1,  4, 7,  2, 5, 8,  3, 6, 9, 
            10,-13,16, 11,14,17, 12,15,18,
             1,  2, 3,  4, 5, 6,  7, 8, 9, 
             1,  0, 0,  0, 1, 0,  0, 0, 1, 
             1,  0, 0,  0, 1, 0,  0, 0, 1, 
            -1,  4, 7,  2,-5, 8,  3, 6,-9,
             10,11,12,-13,14,15, 16,17,18,
            -1,  2, 3,  4,-5, 6,  7, 8,-9, 
             1,  0, 0,  0, 1, 0,  0, 0, 1};
double ka[]={0,  0, 0,  0, 0, 0,  0, 0, 0,        /* All of matrix skew(A) */
             1,  4, 7,  2, 5, 8,  3, 6, 9,        /* Diagonal set to zero  */ 
            10,-13,16, 11,14,17, 12,15,18,   /* And upper triangle negated */
            -1, -2,-3, -4,-5,-6, -7,-8,-9,
             0, 0, 0,  0, 0, 0,  0, 0, 0, 
             0, 0, 0,  0, 0, 0,  0, 0, 0, 
            -1,  4, 7,  2,-5, 8,  3, 6,-9,
            -10,-11,-12,13,-14,-15,-16,-17,-18,
             1,-2,-3, -4, 5,-6, -7,-8, 9, 
             0, 0, 0,  0, 0, 0,  0, 0, 0 };
int bindx[]={1,2,4,1,2,3,4,1,3,4};
int bpntrb[]={1,4,6,8};
int bpntre[]={4,6,8,11};

double la[]= {1,  0, 0,  0, 1, 0,  0, 0, 1,    /* lower triangular part */
             1,  4, 7,  2, 5, 8,  3, 6, 9, 
            10,-13,16, 11,14,17, 12,15,18,
             1,  0, 0,  0, 1, 0,  0, 0, 1, 
             1,  0, 0,  0, 1, 0,  0, 0, 1, 
            -1,  4, 7,  2,-5, 8,  3, 6,-9,
             1,  0, 0,  0, 1, 0,  0, 0, 1};
int lbindx[]={1,2,4,2,3,4,4};
int lbpntrb[]={1,4,5,7};
int lbpntre[]={4,5,7,8};
double ua[]= {1,  0, 0,  0, 1, 0,  0, 0, 1,    /* upper triangular part */ 
             1,  2, 3,  4, 5, 6,  7, 8, 9, 
             1,  0, 0,  0, 1, 0,  0, 0, 1, 
             1,  0, 0,  0, 1, 0,  0, 0, 1, 
             10,11,12,-13,14,15, 16,17,18,
            -1,  2, 3,  4,-5, 6,  7, 8,-9, 
             1,  0, 0,  0, 1, 0,  0, 0, 1};
int ubindx[]={1,1,2,3,1,3,4};
int ubpntrb[]={1,2,4,5};
int ubpntre[]={2,4,5,8};
double b[]={1,2,3,4,5,6,7,8,9,10,11,12,
            1,2,3,4,5,6,7,8,9,10,11,12};
double c[]={1,2,3,4,5,6,7,8,9,10,11,12,
            1,2,3,4,5,6,7,8,9,10,11,12};
double d[]={1,2,3,4,5,6,7,8,9,10,11,12,
            1,2,3,4,5,6,7,8,9,10,11,12};
double check[]={1,2,3,4,5,6,7,8,9,10,11,12,
            1,2,3,4,5,6,7,8,9,10,11,12};
int mb=4, kb=4, lb=3, m=12, ldb=12, ldc=12;
/* Begin description for rectangular matrix */
double ra[]= {1,  0, 0,  0, 1, 0,  0, 0, 1,           /* All of matrix RA */
             1,  4, 7,  2, 5, 8,  3, 6, 9, 
            10,-13,16, 11,14,17, 12,15,18,
             1,  2, 3,  4, 5, 6,  7, 8, 9, 
             1,  0, 0,  0, 1, 0,  0, 0, 1, 
             1,  0, 0,  0, 1, 0,  0, 0, 1, 
            -1,  4, 7,  2,-5, 8,  3, 6,-9,
             10,11,12,-13,14,15, 16,17,18,
            -1,  2, 3,  4,-5, 6,  7, 8,-9, 
             1,  0, 0,  0, 1, 0,  0, 0, 1,
             1,  0, 0,  0, 1, 0,  0, 0, 1,
             1,  0, 0,  0, 1, 0,  0, 0, 1,
             1,  0, 0,  0, 1, 0,  0, 0, 1,
             1,  0, 0,  0, 1, 0,  0, 0, 1 };
int rbindx[]={1,2,4,1,2,3,4,1,3,4,1,2,3,4};
int rbpntrb[]={1,4,6,8,11};
int rbpntre[]={4,6,8,11,15};
int rmb=4, rkb=5, rlb=3, rm=12, rldb=15, rldc=15;
double rb[]={1,2,3,4,5,6,7,8,9,10,11,12,1,1,1,
            1,2,3,4,5,6,7,8,9,10,11,12,1,1,1};
double rc[]={1,2,3,4,5,6,7,8,9,10,11,12,1,1,1,
            1,2,3,4,5,6,7,8,9,10,11,12,1,1,1};
double rd[]={1,2,3,4,5,6,7,8,9,10,11,12,1,1,1,
            1,2,3,4,5,6,7,8,9,10,11,12,1,1,1};
double rcheck[]={1,2,3,4,5,6,7,8,9,10,11,12,1,1,1,
            1,2,3,4,5,6,7,8,9,10,11,12,1,1,1};
double rsumb[]={22,26,30};
/* end of description for rectangular matrix */
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
lwork = 39;
work = (double *) malloc(lwork*sizeof(double));

/* Get input: alpha and beta */

if (argc != 3 ) {
   printf("Usage:   %s alpha beta \n", argv[0]);
   exit(1);
}

alpha = (double) atof(argv[1]);
beta = (double) atof(argv[2]);

descra[2] = 1;
descra[3] = 1;
descra[4] = 1;


printf("-----------------------------------------------------\n");
printf("  alpha = %e, beta = %e  \n",alpha, beta);
printf("-----------------------------------------------------\n");
for (n=1;n!=3;n++) {                      /* loop on columns in C */
      printf("*** n = %d ***\n",n);       /* (test vector and matrix routines) */

                                          /* First, test general matrices */
      printf("   General matrices:\n");
     
                                          /* Testing rectangular matrices */
      printf("      rectangular\n");
      descra[0] = 0;
      descra[1] = 1;                      /* ignored  */

      for (i=0;i!=m;i++)                  /* Initialize c */
        for (j=0;j!=n;j++)
          c[j*m+i] = i+1;

      transa = 0;
      dbscmm( transa, rmb, n, rkb, alpha, descra, ra,
              rbindx, rbpntrb, rbpntre, rlb, rb, rldb,
              beta, c, ldc, work, lwork);

      for (i=0;i!=n*m;i++) 
        d[i] = c[i] - alpha;

      for (i=0;i!=m;i++)                  /* Initialize c                    */
        for (j=0;j!=n;j++)
          c[j*m+i] = i+1;

/* Call mat-mult with explicit symmtric matrix            */

      transa = 0;
      dbscmm( transa, mb, n, kb, alpha, descra, a,
              bindx, bpntrb, bpntre, lb, b, ldb,
              beta, c, ldc, work, lwork);

      error = resid(n*m, d, c);
      if ( error >= tolerance ){
         errcount++;
         printf("Error for rectangular matmult (no transpose)");
         printf("n = %d.\n",n);
         printf("Residual: %10.6f \n",error);
         for (i=0;i!=n*m;i++) 
          printf("%6.2f  %6.2f\n",d[i], c[i]);
      }

      for (i=0;i!=m;i++)                  /* Initialize rc */
        for (j=0;j!=n;j++)
          rc[j*(m+lb)+i] = i+1;
      for (i=m;i!=m+lb;i++)                  /* Initialize rc */
        for (j=0;j!=n;j++)
          rc[j*(m+lb)+i] = 1;

      transa = 1;
      dbscmm( transa, rmb, n, rkb, alpha, descra, ra,
              rbindx, rbpntrb, rbpntre, rlb, b, ldb,
              beta, rc, rldc, work, lwork);

      error = resid(m, c, rc);
      for (j=0;j<lb;j++)
        error += alpha*rsumb[j] + beta - rc[m+j];
      if ( error >= tolerance ){
         errcount++;
         printf("Error for rectangular matmult (transpose)");
         printf("n = %d.\n",n);
         printf("Residual: %10.6f \n",error);
      }
     
      descra[0] = 0;
      descra[1] = 1;                      /* lower triangular matrix */
      printf("      lower triangular\n");

      for (i=0;i!=m;i++)                  /* Initialize c */
        for (j=0;j!=n;j++)
          c[j*m+i] = i+1;

/* Call triangular mat-mult with lower triangular matrix     */

      transa = 0;
      dbscmm( transa, mb, n, kb, alpha, descra, la,
              lbindx, lbpntrb, lbpntre, lb, b, ldb,
              beta, c, ldc, work, lwork);

      for (i=0;i!=n*m;i++) 
        d[i] = c[i];
  
      descra[1] = 2;                      /* upper triangular matrix */
      printf("      upper triangular\n");
 
      for (i=0;i!=m;i++)
        for (j=0;j!=n;j++)
          c[j*m+i] = i+1;

/* Call mat-mult with upper triangular matrix: */
      transa = 1;
      dbscmm( transa, mb, n, kb, alpha, descra, ua,
              ubindx, ubpntrb, ubpntre, lb, b, ldb,
              beta, c, ldc, work, lwork);

      error = resid(n*m, d, c);
      if ( error >= tolerance ){
         errcount++;
         printf("Error for upper(or lower) general matmult");
         printf("n = %d.\n",n);
         printf("Residual: %10.6f \n",error);
         for (i=0;i!=n*m;i++) 
          printf("%6.2f  %6.2f\n",d[i], c[i]);
      }
                                          /* Second, test symmetric matrices */
      printf("   Symmetric matrices:\n");

      descra[0] = 0;                      /* mat-mult with explicit matrix   */

      for (i=0;i!=m;i++)                  /* Initialize c                    */
        for (j=0;j!=n;j++)
          c[j*m+i] = i+1;

/* Call mat-mult with explicit symmtric matrix            */

      transa = 0;
      dbscmm( transa, mb, n, kb, alpha, descra, a,
              bindx, bpntrb, bpntre, lb, b, ldb,
              beta, c, ldc, work, lwork);

      for (i=0;i!=n*m;i++)                /* copy result to d        */
        d[i] = c[i];

      descra[0] = 1;                      /* symmetry is implicit    */
      descra[1] = 1;                      /* lower triangular matrix */
      printf("      lower triangular\n");

      for (i=0;i!=m;i++)                  /* Initialize c            */
        for (j=0;j!=n;j++)
          c[j*m+i] = i+1;

/* Call symmetric mat-mult with lower triangular matrix     */

      transa = 0;
      dbscmm( transa, mb, n, kb, alpha, descra, la,
              lbindx, lbpntrb, lbpntre, lb, b, ldb,
              beta, c, ldc, work, lwork);

                                          /* compare explicit to implicit */
      error = resid(n*m, d, c);
      if ( error >= tolerance ){
         errcount++;
         printf("Error for symmetric matmult (lower triangular)");
         printf("n = %d.\n",n);
         printf("Residual: %10.6f \n",error);
         for (i=0;i!=n*m;i++) 
          printf("%6.2f  %6.2f\n",d[i], c[i]);
      }

  
      descra[1] = 2;                      /* upper triangular matrix */
      printf("      upper triangular\n");

      for (i=0;i!=m;i++)                  /* Initialize c            */
        for (j=0;j!=n;j++)
          c[j*m+i] = i+1;

/* Call symmetric mat-mult with upper triangular matrix     */
      transa = 0;
      dbscmm( transa, mb, n, kb, alpha, descra, ua,
              ubindx, ubpntrb, ubpntre, lb, b, ldb,
              beta, c, ldc, work, lwork);

                                          /* compare explicit to implicit */
      error = resid(n*m, d, c);
      if ( error >= tolerance ){
         errcount++;
         printf("Error for symmetric matmult (upper triangular)");
         printf("n = %d.\n",n);
         printf("Residual: %10.6f \n",error);
         for (i=0;i!=n*m;i++) 
          printf("%6.2f  %6.2f\n",d[i], c[i]);
      }
 
                                          /* Third, test skew-symmetric matrices */
      printf("   Skew-Symmetric matrices:\n");

      descra[0] = 0;                      /* mat-mult with explicit matrix   */

      for (i=0;i!=m;i++)                  /* Initialize c                    */
        for (j=0;j!=n;j++)
          c[j*m+i] = i+1;

/* Call mat-mult with explicit skew-symmetric matrix            */
      transa = 0;
      dbscmm( transa, mb, n, kb, alpha, descra, ka,
              bindx, bpntrb, bpntre, lb, b, ldb,
              beta, c, ldc, work, lwork);

      for (i=0;i!=n*m;i++)                /* copy result to d        */
        d[i] = c[i];

      descra[0] = 4;                      /* symmetry is implicit    */
      descra[1] = 1;                      /* lower triangular matrix */
      printf("      lower triangular (no transp)\n");

      for (i=0;i!=m;i++)                  /* Initialize c            */
        for (j=0;j!=n;j++)
          c[j*m+i] = i+1;

/* Call skew-symmetric mat-mult with lower triangular matrix     */
      transa = 0;
      dbscmm( transa, mb, n, kb, alpha, descra, la,
              lbindx, lbpntrb, lbpntre, lb, b, ldb,
              beta, c, ldc, work, lwork);

                                          /* compare explicit to implicit */
      error = resid(n*m, d, c);
      if ( error >= tolerance ){
         errcount++;
         printf("Error for skew-symmetric matmult (lower triangular)");
         printf("n = %d.\n",n);
         printf("Residual: %10.6f \n",error);
         for (i=0;i!=n*m;i++) 
          printf("%6.2f  %6.2f\n",d[i], c[i]);
      }

  
      descra[1] = 2;                      /* upper triangular matrix */
      printf("      upper triangular (transp)\n");

      for (i=0;i!=m;i++)                  /* Initialize c            */
        for (j=0;j!=n;j++)
          c[j*m+i] = i+1;

/* Call skew-symmetric mat-mult with upper triangular matrix     */
      transa = 1;                         /* Use transpose to get same    */
                                          /* matrix as lower triangular   */
      dbscmm( transa, mb, n, kb, alpha, descra, ua,
              ubindx, ubpntrb, ubpntre, lb, b, ldb,
              beta, c, ldc, work, lwork);

                                          /* compare explicit to implicit */
      error = resid(n*m, d, c);
      if ( error >= tolerance ){
         errcount++;
         printf("Error for skew-symmetric matmult (upper triangular)");
         printf("n = %d.\n",n);
         printf("Residual: %10.6f \n",error);
         for (i=0;i!=n*m;i++) 
          printf("%6.2f  %6.2f\n",d[i], c[i]);
      }
                                          /* Now, work with transp of lower */
                                          /* and upper triangular matrix,   */
                                          /* results should be negation of  */
                                          /* explicit matrix multiply       */
                                          /* check by taking alpha = -alpha */

      descra[1] = 1;                      /* lower triangular matrix */
      printf("      lower triangular (transp)\n");

      for (i=0;i!=m;i++)                  /* Initialize c            */
        for (j=0;j!=n;j++)
          c[j*m+i] = i+1;

/* Call skew-symmetric mat-mult with lower triangular matrix     */
      transa = 1;
      dbscmm( transa, mb, n, kb, -1.0*alpha, descra, la,
              lbindx, lbpntrb, lbpntre, lb, b, ldb,
              beta, c, ldc, work, lwork);

                                          /* compare explicit to implicit */
      error = resid(n*m, d, c);
      if ( error >= tolerance ){
         errcount++;
         printf("Error for skew-symmetric matmult (lower triangular)");
         printf("n = %d.\n",n);
         printf("Residual: %10.6f \n",error);
         for (i=0;i!=n*m;i++) 
          printf("%6.2f  %6.2f\n",d[i], c[i]);
      }

  
      descra[1] = 2;                      /* upper triangular matrix */
      printf("      upper triangular (no transp)\n");

      for (i=0;i!=m;i++)                  /* Initialize c            */
        for (j=0;j!=n;j++)
          c[j*m+i] = i+1;

/* Call skew-symmetric mat-mult with upper triangular matrix     */
      transa = 0;                         /* Use transpose to get same    */
                                          /* matrix as lower triangular   */
      dbscmm( transa, mb, n, kb, -1.0*alpha, descra, ua,
              ubindx, ubpntrb, ubpntre, lb, b, ldb,
              beta, c, ldc, work, lwork);

                                          /* compare explicit to implicit */
      error = resid(n*m, d, c);
      if ( error >= tolerance ){
         errcount++;
         printf("Error for skew-symmetric matmult (upper triangular)");
         printf("n = %d.\n",n);
         printf("Residual: %10.6f \n",error);
         for (i=0;i!=n*m;i++) 
          printf("%6.2f  %6.2f\n",d[i], c[i]);
      }


} /* close loop on n */

if ( errcount > 0 ) 
  printf("%d errors in dtbscmm run for alpha = %e, beta = %e\n",errcount,alpha, beta);

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

