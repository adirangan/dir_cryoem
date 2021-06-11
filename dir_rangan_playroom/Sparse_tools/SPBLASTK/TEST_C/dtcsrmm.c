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
double diag[]={.5,-.5,1,1,1,1,1,1,1,-.25,.25};

double a[]={1,        1,     -1,1,              /* All of matrix A */
              1,      2,       -1,
                1,    3,          
                  1,  4,          
                    1,5,          
            1,2,3,4,5,1,6,7,8,    
                      6,1,        
                      7,  1,       
                      8,    1,    
           -1,                1,  
            1,-1,               1};

double ka[]={0,       -1,         1,-1,              /* All of matrix skew(A) */
              0,      -2,            1,              /* Diagonal set to zero  */
                0,    -3,                            /* And upper triangle negated */
                  0,  -4,          
                    0,-5,          
            1,2,3,4,5, 0,-6,-7,-8,    
                       6, 0,        
                       7,    0,       
                       8,      0,    
           -1,                   0,  
            1,-1,                  0};
int indx[]={1,6,10,11,2,6,11,3,6,4,6,5,6,1,2,3,4,5,6,7,8,9,6,7,6,8,6,9,1,10,1,2,11};
int pntrb[]={1,5,8,10,12,14,23,25,27,29,31};
int pntre[]={5,8,10,12,14,23,25,27,29,31,34};

double la[]={1,                          /* lower triangular part */
              1, 
                1,
                  1,
                    1,
            1,2,3,4,5,1,
                      6,1,        
                      7,  1,       
                      8,    1,    
           -1,                1,  
            1,-1,               1};
int lindx[]={1,2,3,4,5,1,2,3,4,5,6,6,7,6,8,6,9,1,10,1,2,11};
int lpntrb[]={1,2,3,4,5, 6,12,14,16,18,20};
int lpntre[]={2,3,4,5,6,12,14,16,18,20,23};

double ua[]={1,        1,     -1,1,           /* upper triangular part */ 
               1,      2,       -1,
                 1,    3,          
                   1,  4,          
                     1,5,          
                       1,6,7,8,    
                         1,        
                           1,       
                             1,    
                               1,  
                                 1};
int uindx[]={1,6,10,11, 2,6,11, 3,6, 4,6, 5,6, 6,7,8,9, 7,8,9,10,11 };
int upntrb[]={1,5, 8,10,12,14,18,19,20,21,22};
int upntre[]={5, 8,10,12,14,18,19,20,21,22,23};

double b[]={1,2,3,4,5,6,7,8,9,10,11,
            1,2,3,4,5,6,7,8,9,10,11};
double c[]={1,2,3,4,5,6,7,8,9,10,11,
            1,2,3,4,5,6,7,8,9,10,11};
double d[]={1,2,3,4,5,6,7,8,9,10,11,
            1,2,3,4,5,6,7,8,9,10,11};
double check[]={1,2,3,4,5,6,7,8,9,10,11,
            1,2,3,4,5,6,7,8,9,10,11};
int k=11, m=11, ldb=11, ldc=11;
/* Begin description of rectangular matrix */
double ra[]={1,         1,      -1,1,1,              /* All of matrix RA */
               1,       2,        -1,1,
                 1,     3,           1,
                   1,   4,           1,
                     1, 5,           1,
             1,2,3,4,5, 1, 6,7,8,    1,
                        6, 1,        1,
                        7,   1,      1, 
                        8,     1,    1, 
            -1,                  1,  1,
            1,-1,                 1,1};
int rindx[]={1,6,10,11,12,2,6,11,12,3,6,12,4,6,12,5,6,12,
             1,2,3,4,5,6,7,8,9,12,6,7,12,6,8,12,6,9,12,1,10,12,1,2,11,12};
int rpntrb[]={1,6,10,13,16,19,29,32,35,38,41};
int rpntre[]={6,10,13,16,19,29,32,35,38,41,45};
double rb[]={1,2,3,4,5,6,7,8,9,10,11,1,
            1,2,3,4,5,6,7,8,9,10,11,1};
double rc[]={1,2,3,4,5,6,7,8,9,10,11,1,
            1,2,3,4,5,6,7,8,9,10,11,1};
double rd[]={1,2,3,4,5,6,7,8,9,10,11,1,
            1,2,3,4,5,6,7,8,9,10,11,1};
double rcheck[]={1,2,3,4,5,6,7,8,9,10,11,1,
            1,2,3,4,5,6,7,8,9,10,11,1};
double rsumb = 66;
int rk=12, rm=11, rldb=12, rldc=12;
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

descra[2] = 0;
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
      dcsrmm( transa, rm, n, rk, alpha, descra, ra,
              rindx, rpntrb, rpntre, rb, rldb,
              beta, c, ldc, work, lwork);

      for (i=0;i!=n*m;i++) 
        d[i] = c[i] - alpha;

      for (i=0;i!=m;i++)                  /* Initialize c                    */
        for (j=0;j!=n;j++)
          c[j*m+i] = i+1;

/* Call mat-mult with explicit symmtric matrix            */

      transa = 0;
      dcsrmm( transa, m, n, k, alpha, descra, a,
              indx, pntrb, pntre, b, ldb,
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
      for (i=m;i!=m+1;i++)                  /* Initialize rc */
        for (j=0;j!=n;j++)
          rc[j*(m+1)+i] = 1;

      transa = 1;
      dcsrmm( transa, rm, n, rk, alpha, descra, ra,
              rindx, rpntrb, rpntre, b, ldb,
              beta, rc, rldc, work, lwork);

      error = resid(m, c, rc);
      error += alpha*rsumb + beta - rc[m];
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
      dcsrmm( transa, m, n, k, alpha, descra, la,
              lindx, lpntrb, lpntre, b, ldb,
              beta, c, ldc, work, lwork);

      for (i=0;i!=n*m;i++) 
        d[i] = c[i];
  
      descra[1] = 2;                      /* upper triangular matrix */
      printf("      upper triangular\n");
 
      for (i=0;i!=m;i++)
        for (j=0;j!=n;j++)
          c[j*m+i] = i+1;

/* Call triangular mat-mult with upper triangular matrix: */

      transa = 1;
      dcsrmm( transa, m, n, k, alpha, descra, ua,
              uindx, upntrb, upntre, b, ldb,
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
      dcsrmm( transa, m, n, k, alpha, descra, a,
              indx, pntrb, pntre, b, ldb,
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
      dcsrmm( transa, m, n, k, alpha, descra, la,
              lindx, lpntrb, lpntre, b, ldb,
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

/* Call symmetric mat-mult with lower triangular matrix     */

      transa = 0;
      dcsrmm( transa, m, n, k, alpha, descra, ua,
              uindx, upntrb, upntre, b, ldb,
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
      dcsrmm( transa, m, n, k, alpha, descra, ka,
              indx, pntrb, pntre, b, ldb,
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
      dcsrmm( transa, m, n, k, alpha, descra, la,
              lindx, lpntrb, lpntre, b, ldb,
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
      dcsrmm( transa, m, n, k, alpha, descra, ua,
              uindx, upntrb, upntre, b, ldb,
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
      dcsrmm( transa, m, n, k, -1.0*alpha, descra, la,
              lindx, lpntrb, lpntre, b, ldb,
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
      dcsrmm( transa, m, n, k, -1.0*alpha, descra, ua,
              uindx, upntrb, upntre, b, ldb,
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
  printf("%d errors in dtcsrmm run for alpha = %e, beta = %e\n",errcount,alpha, beta);

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

