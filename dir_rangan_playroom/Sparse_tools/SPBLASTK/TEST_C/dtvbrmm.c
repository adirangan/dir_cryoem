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
double dv[]={.5, 0,  0,-.5,
             1, 0, 0,  0, 1, 0,  0, 0, 1, 
             1, 
             1, 0, 0,  0, 1, 0,  0, 0, 1, 
             -.25, 0,  0, .25};
double diag[]={.5,-.5,1,1,1,1,1,1,1,-.25,.25};

double a[]={1, 0,  0, 1,                        /* All of matrix A */
            1, 2,
            -1, 0,  1, -1,
            1, 0, 0,  0, 1, 0,  0, 0, 1, 
            3, 4, 5,
            1, 2,
            3, 4, 5,
            1, 
            6, 7, 8,
            6, 7, 8, 
            1, 0, 0,  0, 1, 0,  0, 0, 1, 
            -1, 1,  0, -1,
            1, 0,  0, 1};
double ka[]={0, 0,  0, 0,                        /* All of matrix skew(A) */
            -1, -2,                              /* Diagonal set to zero  */
            1, 0,  -1, 1,                        /* And upper triangle negated */
            0, 0, 0,  0, 0, 0,  0, 0, 0, 
            -3, -4, -5,
            1, 2,
            3, 4, 5,
            0, 
            -6, -7, -8,
            6, 7, 8, 
            0, 0, 0,  0, 0, 0,  0, 0, 0, 
            -1, 1,  0, -1,
            0, 0,  0, 0};
int indx[]={1,5,7,11,20,23,25,28,29,32,35,44,48,52};
int bindx[]={1,3,5,2,3,1,2,3,4,3,4,1,5};
int rpntr[]={1,3,6,7,10,12};
int cpntr[]={1,3,6,7,10,12};
int bpntrb[]={1,4,6,10,12};
int bpntre[]={4,6,10,12,14};

double la[]={1, 0,  0, 1,                     /* lower triangular part */
            1, 0, 0,  0, 1, 0,  0, 0, 1, 
            1, 2,
            3, 4, 5,
            1, 
            6, 7, 8,
            1, 0, 0,  0, 1, 0,  0, 0, 1,
           -1, 1,  0, -1,
            1, 0,  0, 1};
int lindx[]={1,5,14,16,19,20,23,32,36,40};
int lbindx[]={1,2,1,2,3,3,4,1,5};
int lrpntr[]={1,3,6,7,10,12};
int lcpntr[]={1,3,6,7,10,12};
int lbpntrb[]={1,2,3,6,8};
int lbpntre[]={2,3,6,8,10};
double ua[]={1, 0,  0, 1,                     /* upper triangular part */
             1, 2,
             -1, 0,  1, -1,
             1, 0, 0,  0, 1, 0,  0, 0, 1, 
             3, 4, 5,
             1, 
             6, 7, 8,
             1, 0, 0,  0, 1, 0,  0, 0, 1, 
             1, 0, 0, 1};
int uindx[]={1,5,7,11,20,23,24,27,36,40};
int ubindx[]={1,3,5,2,3,3,4,4,5};
int urpntr[]={1,3,6,7,10,12};
int ucpntr[]={1,3,6,7,10,12};
int ubpntrb[]={1,4,6,8,9};
int ubpntre[]={4,6,8,9,10};
double b[]={1,2,3,4,5,6,7,8,9,10,11,
            1,2,3,4,5,6,7,8,9,10,11};
double c[]={1,2,3,4,5,6,7,8,9,10,11,
            1,2,3,4,5,6,7,8,9,10,11};
double d[]={1,2,3,4,5,6,7,8,9,10,11,
            1,2,3,4,5,6,7,8,9,10,11};
double check[]={1,2,3,4,5,6,7,8,9,10,11,
            1,2,3,4,5,6,7,8,9,10,11};
int mb=5, kb=5, m=11, ldb=11, ldc=11;
double ra[]={1, 0,  0, 1,                        /* All of matrix RA */
            1, 2,
            -1, 0,  1, -1,
            1, 1,
            1, 0, 0,  0, 1, 0,  0, 0, 1, 
            3, 4, 5,
            1, 1, 1,
            1, 2,
            3, 4, 5,
            1, 
            6, 7, 8,
            1, 
            6, 7, 8, 
            1, 0, 0,  0, 1, 0,  0, 0, 1, 
            1, 1, 1,
            -1, 1,  0, -1,
            1, 0,  0, 1,
            1, 1};
int rindx[]={1,5,7,11,13,22,25,28,30,33,34,37,38,41,50,53,57,61};
int rbindx[]={1,3,5,6,2,3,6,1,2,3,4,6,3,4,6,1,5,6};
int rrpntr[]={1,3,6,7,10,12};
int rcpntr[]={1,3,6,7,10,12,13};
int rbpntrb[]={1,5,8,13,16};
int rbpntre[]={5,8,13,16,19};
double rb[]={1,2,3,4,5,6,7,8,9,10,11,1,
            1,2,3,4,5,6,7,8,9,10,11,1};
double rc[]={1,2,3,4,5,6,7,8,9,10,11,1,
            1,2,3,4,5,6,7,8,9,10,11,1};
double rd[]={1,2,3,4,5,6,7,8,9,10,11,1,
            1,2,3,4,5,6,7,8,9,10,11,1};
double rcheck[]={1,2,3,4,5,6,7,8,9,10,11,1,
            1,2,3,4,5,6,7,8,9,10,11,1};
int rmb=5, rkb=6, rm=11, rldb=12, rldc=12;
double rsumb = 66;
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
work = (double *) malloc(24*sizeof(double));
lwork = 24;

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
      dvbrmm( transa, rmb, n, rkb, alpha, descra, ra,
              rindx, rbindx, rrpntr, rcpntr, rbpntrb, rbpntre, rb, rldb,
              beta, c, ldc, work, lwork);

      for (i=0;i!=n*m;i++) 
        d[i] = c[i] - alpha;

      for (i=0;i!=m;i++)                  /* Initialize c                    */
        for (j=0;j!=n;j++)
          c[j*m+i] = i+1;

/* Call mat-mult with explicit symmtric matrix            */

      transa = 0;
      dvbrmm( transa, mb, n, kb, alpha, descra, a,
              indx, bindx, rpntr, cpntr, bpntrb, bpntre, b, ldb,
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
          rc[j*(m+1)+i] = i+1;
        for (j=0;j!=n;j++)
          rc[j*(m+1)+m] = 1;

      transa = 1;
      dvbrmm( transa, rmb, n, rkb, alpha, descra, ra,
              rindx, rbindx, rrpntr, rcpntr, rbpntrb, rbpntre, b, ldb,
              beta, rc, rldc, work, lwork);

      error = resid(m, c, rc);
      error += alpha*rsumb + beta - rc[m];

      if ( error >= tolerance ){
         errcount++;
         printf("Error for rectangular matmult (transpose)");
         printf("n = %d.\n",n);
         printf("Residual: %10.6f \n",error);
         for (i=0;i!=m;i++) 
          printf("%6.2f  %6.2f\n",c[i], rc[i]);
         printf("%6.2f  %6.2f\n",alpha*rsumb+beta, rc[m]);
      }
     
      descra[0] = 0;
      descra[1] = 1;                      /* lower triangular matrix */
      printf("      lower triangular\n");

      for (i=0;i!=m;i++)                  /* Initialize c */
        for (j=0;j!=n;j++)
          c[j*m+i] = i+1;

/* Call triangular mat-mult with lower triangular matrix     */

      transa = 0;
      dvbrmm( transa, mb, n, kb, alpha, descra, la,
              lindx, lbindx, lrpntr, lcpntr, lbpntrb, lbpntre, b, ldb,
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
      dvbrmm( transa, mb, n, kb, alpha, descra, ua,
              uindx, ubindx, urpntr, ucpntr, ubpntrb, ubpntre, b, ldb,
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
      dvbrmm( transa, mb, n, kb, alpha, descra, a,
              indx, bindx, rpntr, cpntr, bpntrb, bpntre, b, ldb,
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
      dvbrmm( transa, mb, n, kb, alpha, descra, la,
              lindx, lbindx, lrpntr, lcpntr, lbpntrb, lbpntre, b, ldb,
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
      dvbrmm( transa, mb, n, kb, alpha, descra, ua,
              uindx, ubindx, urpntr, ucpntr, ubpntrb, ubpntre, b, ldb,
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
      dvbrmm( transa, mb, n, kb, alpha, descra, ka,
              indx, bindx, rpntr, cpntr, bpntrb, bpntre, b, ldb,
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
      dvbrmm( transa, mb, n, kb, alpha, descra, la,
              lindx, lbindx, lrpntr, lcpntr, lbpntrb, lbpntre, b, ldb,
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
      dvbrmm( transa, mb, n, kb, alpha, descra, ua,
              uindx, ubindx, urpntr, ucpntr, ubpntrb, ubpntre, b, ldb,
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
      dvbrmm( transa, mb, n, kb, -1.0*alpha, descra, la,
              lindx, lbindx, lrpntr, lcpntr, lbpntrb, lbpntre, b, ldb,
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
      dvbrmm( transa, mb, n, kb, -1.0*alpha, descra, ua,
              uindx, ubindx, urpntr, ucpntr, ubpntrb, ubpntre, b, ldb,
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
  printf("%d errors in dtvbrmm run for alpha = %e, beta = %e\n",errcount,alpha, beta);

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

