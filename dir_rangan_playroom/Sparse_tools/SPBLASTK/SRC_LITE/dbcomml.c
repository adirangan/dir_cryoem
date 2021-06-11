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


/* Created:  Sat Jul 6 14:32:52 EDT 1996 */

#include "dbcomml.h"



void BCO_MatMult_CAB_double(
                 const int mb, const int n, const int kb, 
                 const double *val, const int *bindx, const int *bjndx,
                 const int bnnz, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  const double *ptmp;                                  
  int i,j,jb,je,block;
  int cs,bs,br;
  int index;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  for (l=0;l!=n;l++)                            
    for (i=0;i!=m;i++) *pc++ = 0;                  

  for (i=0;i!=bnnz;i++) {      
      cs = (bindx[i]-1)*lb;     
      bs = (bjndx[i]-1)*lb;        
      pval = &val[i*mm];           
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
      for (l=0;l!=n;l++) {                      
       for (jj=0;jj!=lb;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t =  pb[jj];            
           for (ii=0;ii!=lb;ii++) {      
             pc[ii] += t* (*pval++);        
           }
         } else {
           pval+=lb;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
  }
}

void BCO_MatMult_CATB_double(
                 const int mb, const int n, const int kb, 
                 const double *val, const int *bindx, const int *bjndx, 
                 const int bnnz, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  const double *ptmp;                                  
  int i,j,jb,je,block;
  int cs,bs,cr;
  int index;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  for (l=0;l!=n;l++)                            
    for (i=0;i!=k;i++) *pc++ = 0;                  

  for (i=0;i!=bnnz;i++) {      
      cs = (bjndx[i]-1)*lb;     
      bs = (bindx[i]-1)*lb;        
      pval = &val[i*mm];           
      ptmp = pval;                  
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
      for (l=0;l!=n;l++) {              
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t +=  pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
  } }

void BCOsymm_MatMult_CAB_double(
                 const int mb, const int n, const int kb, 
                 const double *val, const int *bindx, const int *bjndx, 
                 const int bnnz, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  const double *ptmp; 
  int i,j,jb,je,block;
  int cs,css,cr,bs,br;
  int index;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  for (l=0;l!=n;l++)                            
    for (i=0;i!=m;i++) *pc++ = 0;                  

  for (i=0;i!=bnnz;i++) {    
      cs = (bindx[i]-1)*lb;     
      bs = (bjndx[i]-1)*lb;        
      pval = &val[i*mm];           
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
      for (l=0;l!=n;l++) {             
       for (jj=0;jj!=lb;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t =  pb[jj];            
           for (ii=0;ii!=lb;ii++) {      
             pc[ii] += t* (*pval++);        
           }
         } else {
           pval+=lb;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
 
      if ( cs == bs ) {
        continue;                  
      }
 
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
      for (l=0;l!=n;l++) {         
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t +=  pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
  }
}

void BCOskew_MatMult_CAB_double(
                 const int mb, const int n, const int kb, 
                 const double *val, const int *bindx, const int *bjndx, 
                 const int bnnz, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int cs,css,cr,bs,br;
  int index;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  for (l=0;l!=n;l++)                            
    for (i=0;i!=m;i++) *pc++ = 0;                  

  for (i=0;i!=bnnz;i++) {    
      cs = (bindx[i]-1)*lb;     
      bs = (bjndx[i]-1)*lb;        
 
      if ( cs == bs ) {
        continue;                  
      }
      pval = &val[i*mm];           
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
      for (l=0;l!=n;l++) {          
       for (jj=0;jj!=lb;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t =  pb[jj];            
           for (ii=0;ii!=lb;ii++) {      
             pc[ii] += t* (*pval++);        
           }
         } else {
           pval+=lb;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
 
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
      for (l=0;l!=n;l++) {         
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t -=  pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
  }
}

void BCOskew_MatMult_CATB_double(
                 const int mb, const int n, const int kb, 
                 const double *val, const int *bindx, const int *bjndx,
                 const int bnnz, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  int i,j,jb,je,block;
  int cs,css,cr,bs,br;
  int index;
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;
  const double *ptmp;
  int l;                                               

  for (l=0;l!=n;l++)                            
    for (i=0;i!=k;i++) *pc++ = 0;                  

  for (i=0;i!=bnnz;i++) {    
      cs = (bindx[i]-1)*lb;     
      bs = (bjndx[i]-1)*lb;        
 
      if ( cs == bs ) {
        continue;                  
      }
      pval = &val[i*mm];           
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
      for (l=0;l!=n;l++) {             
       for (jj=0;jj!=lb;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t =  pb[jj];            
           for (ii=0;ii!=lb;ii++) {      
             pc[ii] -= t* (*pval++);        
           }
         } else {
           pval+=lb;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
 
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
      for (l=0;l!=n;l++) {         
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t +=  pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
  }
}




void BCO_MatMult_CaAB_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *bindx, const int *bjndx,
                 const int bnnz, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  const double *ptmp;                                  
  int i,j,jb,je,block;
  int cs,bs,br;
  int index;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  for (l=0;l!=n;l++)                            
    for (i=0;i!=m;i++) *pc++ = 0;                  

  for (i=0;i!=bnnz;i++) {      
      cs = (bindx[i]-1)*lb;     
      bs = (bjndx[i]-1)*lb;        
      pval = &val[i*mm];           
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
      for (l=0;l!=n;l++) {                      
       for (jj=0;jj!=lb;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t = alpha * pb[jj];            
           for (ii=0;ii!=lb;ii++) {      
             pc[ii] += t* (*pval++);        
           }
         } else {
           pval+=lb;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
  }
}

void BCO_MatMult_CaATB_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *bindx, const int *bjndx, 
                 const int bnnz, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  const double *ptmp;                                  
  int i,j,jb,je,block;
  int cs,bs,cr;
  int index;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  for (l=0;l!=n;l++)                            
    for (i=0;i!=k;i++) *pc++ = 0;                  

  for (i=0;i!=bnnz;i++) {      
      cs = (bjndx[i]-1)*lb;     
      bs = (bindx[i]-1)*lb;        
      pval = &val[i*mm];           
      ptmp = pval;                  
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
      for (l=0;l!=n;l++) {              
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t += alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
  } }

void BCOsymm_MatMult_CaAB_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *bindx, const int *bjndx, 
                 const int bnnz, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  const double *ptmp; 
  int i,j,jb,je,block;
  int cs,css,cr,bs,br;
  int index;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  for (l=0;l!=n;l++)                            
    for (i=0;i!=m;i++) *pc++ = 0;                  

  for (i=0;i!=bnnz;i++) {    
      cs = (bindx[i]-1)*lb;     
      bs = (bjndx[i]-1)*lb;        
      pval = &val[i*mm];           
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
      for (l=0;l!=n;l++) {             
       for (jj=0;jj!=lb;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t = alpha * pb[jj];            
           for (ii=0;ii!=lb;ii++) {      
             pc[ii] += t* (*pval++);        
           }
         } else {
           pval+=lb;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
 
      if ( cs == bs ) {
        continue;                  
      }
 
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
      for (l=0;l!=n;l++) {         
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t += alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
  }
}

void BCOskew_MatMult_CaAB_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *bindx, const int *bjndx, 
                 const int bnnz, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int cs,css,cr,bs,br;
  int index;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  for (l=0;l!=n;l++)                            
    for (i=0;i!=m;i++) *pc++ = 0;                  

  for (i=0;i!=bnnz;i++) {    
      cs = (bindx[i]-1)*lb;     
      bs = (bjndx[i]-1)*lb;        
 
      if ( cs == bs ) {
        continue;                  
      }
      pval = &val[i*mm];           
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
      for (l=0;l!=n;l++) {          
       for (jj=0;jj!=lb;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t = alpha * pb[jj];            
           for (ii=0;ii!=lb;ii++) {      
             pc[ii] += t* (*pval++);        
           }
         } else {
           pval+=lb;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
 
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
      for (l=0;l!=n;l++) {         
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t -= alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
  }
}

void BCOskew_MatMult_CaATB_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *bindx, const int *bjndx,
                 const int bnnz, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  int i,j,jb,je,block;
  int cs,css,cr,bs,br;
  int index;
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;
  const double *ptmp;
  int l;                                               

  for (l=0;l!=n;l++)                            
    for (i=0;i!=k;i++) *pc++ = 0;                  

  for (i=0;i!=bnnz;i++) {    
      cs = (bindx[i]-1)*lb;     
      bs = (bjndx[i]-1)*lb;        
 
      if ( cs == bs ) {
        continue;                  
      }
      pval = &val[i*mm];           
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
      for (l=0;l!=n;l++) {             
       for (jj=0;jj!=lb;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t = alpha * pb[jj];            
           for (ii=0;ii!=lb;ii++) {      
             pc[ii] -= t* (*pval++);        
           }
         } else {
           pval+=lb;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
 
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
      for (l=0;l!=n;l++) {         
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t += alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
  }
}




void BCO_MatMult_CABC_double(
                 const int mb, const int n, const int kb, 
                 const double *val, const int *bindx, const int *bjndx,
                 const int bnnz, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  const double *ptmp;                                  
  int i,j,jb,je,block;
  int cs,bs,br;
  int index;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;


  for (i=0;i!=bnnz;i++) {      
      cs = (bindx[i]-1)*lb;     
      bs = (bjndx[i]-1)*lb;        
      pval = &val[i*mm];           
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
      for (l=0;l!=n;l++) {                      
       for (jj=0;jj!=lb;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t =  pb[jj];            
           for (ii=0;ii!=lb;ii++) {      
             pc[ii] += t* (*pval++);        
           }
         } else {
           pval+=lb;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
  }
}

void BCO_MatMult_CATBC_double(
                 const int mb, const int n, const int kb, 
                 const double *val, const int *bindx, const int *bjndx, 
                 const int bnnz, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  const double *ptmp;                                  
  int i,j,jb,je,block;
  int cs,bs,cr;
  int index;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;


  for (i=0;i!=bnnz;i++) {      
      cs = (bjndx[i]-1)*lb;     
      bs = (bindx[i]-1)*lb;        
      pval = &val[i*mm];           
      ptmp = pval;                  
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
      for (l=0;l!=n;l++) {              
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t +=  pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
  } }

void BCOsymm_MatMult_CABC_double(
                 const int mb, const int n, const int kb, 
                 const double *val, const int *bindx, const int *bjndx, 
                 const int bnnz, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  const double *ptmp; 
  int i,j,jb,je,block;
  int cs,css,cr,bs,br;
  int index;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;


  for (i=0;i!=bnnz;i++) {    
      cs = (bindx[i]-1)*lb;     
      bs = (bjndx[i]-1)*lb;        
      pval = &val[i*mm];           
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
      for (l=0;l!=n;l++) {             
       for (jj=0;jj!=lb;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t =  pb[jj];            
           for (ii=0;ii!=lb;ii++) {      
             pc[ii] += t* (*pval++);        
           }
         } else {
           pval+=lb;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
 
      if ( cs == bs ) {
        continue;                  
      }
 
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
      for (l=0;l!=n;l++) {         
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t +=  pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
  }
}

void BCOskew_MatMult_CABC_double(
                 const int mb, const int n, const int kb, 
                 const double *val, const int *bindx, const int *bjndx, 
                 const int bnnz, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int cs,css,cr,bs,br;
  int index;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;


  for (i=0;i!=bnnz;i++) {    
      cs = (bindx[i]-1)*lb;     
      bs = (bjndx[i]-1)*lb;        
 
      if ( cs == bs ) {
        continue;                  
      }
      pval = &val[i*mm];           
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
      for (l=0;l!=n;l++) {          
       for (jj=0;jj!=lb;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t =  pb[jj];            
           for (ii=0;ii!=lb;ii++) {      
             pc[ii] += t* (*pval++);        
           }
         } else {
           pval+=lb;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
 
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
      for (l=0;l!=n;l++) {         
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t -=  pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
  }
}

void BCOskew_MatMult_CATBC_double(
                 const int mb, const int n, const int kb, 
                 const double *val, const int *bindx, const int *bjndx,
                 const int bnnz, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  int i,j,jb,je,block;
  int cs,css,cr,bs,br;
  int index;
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;
  const double *ptmp;
  int l;                                               


  for (i=0;i!=bnnz;i++) {    
      cs = (bindx[i]-1)*lb;     
      bs = (bjndx[i]-1)*lb;        
 
      if ( cs == bs ) {
        continue;                  
      }
      pval = &val[i*mm];           
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
      for (l=0;l!=n;l++) {             
       for (jj=0;jj!=lb;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t =  pb[jj];            
           for (ii=0;ii!=lb;ii++) {      
             pc[ii] -= t* (*pval++);        
           }
         } else {
           pval+=lb;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
 
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
      for (l=0;l!=n;l++) {         
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t +=  pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
  }
}




void BCO_MatMult_CaABC_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *bindx, const int *bjndx,
                 const int bnnz, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  const double *ptmp;                                  
  int i,j,jb,je,block;
  int cs,bs,br;
  int index;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;


  for (i=0;i!=bnnz;i++) {      
      cs = (bindx[i]-1)*lb;     
      bs = (bjndx[i]-1)*lb;        
      pval = &val[i*mm];           
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
      for (l=0;l!=n;l++) {                      
       for (jj=0;jj!=lb;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t = alpha * pb[jj];            
           for (ii=0;ii!=lb;ii++) {      
             pc[ii] += t* (*pval++);        
           }
         } else {
           pval+=lb;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
  }
}

void BCO_MatMult_CaATBC_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *bindx, const int *bjndx, 
                 const int bnnz, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  const double *ptmp;                                  
  int i,j,jb,je,block;
  int cs,bs,cr;
  int index;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;


  for (i=0;i!=bnnz;i++) {      
      cs = (bjndx[i]-1)*lb;     
      bs = (bindx[i]-1)*lb;        
      pval = &val[i*mm];           
      ptmp = pval;                  
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
      for (l=0;l!=n;l++) {              
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t += alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
  } }

void BCOsymm_MatMult_CaABC_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *bindx, const int *bjndx, 
                 const int bnnz, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  const double *ptmp; 
  int i,j,jb,je,block;
  int cs,css,cr,bs,br;
  int index;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;


  for (i=0;i!=bnnz;i++) {    
      cs = (bindx[i]-1)*lb;     
      bs = (bjndx[i]-1)*lb;        
      pval = &val[i*mm];           
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
      for (l=0;l!=n;l++) {             
       for (jj=0;jj!=lb;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t = alpha * pb[jj];            
           for (ii=0;ii!=lb;ii++) {      
             pc[ii] += t* (*pval++);        
           }
         } else {
           pval+=lb;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
 
      if ( cs == bs ) {
        continue;                  
      }
 
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
      for (l=0;l!=n;l++) {         
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t += alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
  }
}

void BCOskew_MatMult_CaABC_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *bindx, const int *bjndx, 
                 const int bnnz, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int cs,css,cr,bs,br;
  int index;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;


  for (i=0;i!=bnnz;i++) {    
      cs = (bindx[i]-1)*lb;     
      bs = (bjndx[i]-1)*lb;        
 
      if ( cs == bs ) {
        continue;                  
      }
      pval = &val[i*mm];           
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
      for (l=0;l!=n;l++) {          
       for (jj=0;jj!=lb;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t = alpha * pb[jj];            
           for (ii=0;ii!=lb;ii++) {      
             pc[ii] += t* (*pval++);        
           }
         } else {
           pval+=lb;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
 
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
      for (l=0;l!=n;l++) {         
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t -= alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
  }
}

void BCOskew_MatMult_CaATBC_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *bindx, const int *bjndx,
                 const int bnnz, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  int i,j,jb,je,block;
  int cs,css,cr,bs,br;
  int index;
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;
  const double *ptmp;
  int l;                                               


  for (i=0;i!=bnnz;i++) {    
      cs = (bindx[i]-1)*lb;     
      bs = (bjndx[i]-1)*lb;        
 
      if ( cs == bs ) {
        continue;                  
      }
      pval = &val[i*mm];           
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
      for (l=0;l!=n;l++) {             
       for (jj=0;jj!=lb;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t = alpha * pb[jj];            
           for (ii=0;ii!=lb;ii++) {      
             pc[ii] -= t* (*pval++);        
           }
         } else {
           pval+=lb;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
 
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
      for (l=0;l!=n;l++) {         
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t += alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
  }
}




void BCO_MatMult_CABbC_double(
                 const int mb, const int n, const int kb, 
                 const double *val, const int *bindx, const int *bjndx,
                 const int bnnz, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  const double *ptmp;                                  
  int i,j,jb,je,block;
  int cs,bs,br;
  int index;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  for (l=0;l!=n;l++)                            
    for (i=0;i!=m;i++) *pc++ *= beta;                  

  for (i=0;i!=bnnz;i++) {      
      cs = (bindx[i]-1)*lb;     
      bs = (bjndx[i]-1)*lb;        
      pval = &val[i*mm];           
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
      for (l=0;l!=n;l++) {                      
       for (jj=0;jj!=lb;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t =  pb[jj];            
           for (ii=0;ii!=lb;ii++) {      
             pc[ii] += t* (*pval++);        
           }
         } else {
           pval+=lb;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
  }
}

void BCO_MatMult_CATBbC_double(
                 const int mb, const int n, const int kb, 
                 const double *val, const int *bindx, const int *bjndx, 
                 const int bnnz, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  const double *ptmp;                                  
  int i,j,jb,je,block;
  int cs,bs,cr;
  int index;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  for (l=0;l!=n;l++)                            
    for (i=0;i!=k;i++) *pc++ *= beta;                  

  for (i=0;i!=bnnz;i++) {      
      cs = (bjndx[i]-1)*lb;     
      bs = (bindx[i]-1)*lb;        
      pval = &val[i*mm];           
      ptmp = pval;                  
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
      for (l=0;l!=n;l++) {              
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t +=  pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
  } }

void BCOsymm_MatMult_CABbC_double(
                 const int mb, const int n, const int kb, 
                 const double *val, const int *bindx, const int *bjndx, 
                 const int bnnz, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  const double *ptmp; 
  int i,j,jb,je,block;
  int cs,css,cr,bs,br;
  int index;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  for (l=0;l!=n;l++)                            
    for (i=0;i!=m;i++) *pc++ *= beta;                  

  for (i=0;i!=bnnz;i++) {    
      cs = (bindx[i]-1)*lb;     
      bs = (bjndx[i]-1)*lb;        
      pval = &val[i*mm];           
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
      for (l=0;l!=n;l++) {             
       for (jj=0;jj!=lb;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t =  pb[jj];            
           for (ii=0;ii!=lb;ii++) {      
             pc[ii] += t* (*pval++);        
           }
         } else {
           pval+=lb;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
 
      if ( cs == bs ) {
        continue;                  
      }
 
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
      for (l=0;l!=n;l++) {         
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t +=  pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
  }
}

void BCOskew_MatMult_CABbC_double(
                 const int mb, const int n, const int kb, 
                 const double *val, const int *bindx, const int *bjndx, 
                 const int bnnz, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int cs,css,cr,bs,br;
  int index;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  for (l=0;l!=n;l++)                            
    for (i=0;i!=m;i++) *pc++ *= beta;                  

  for (i=0;i!=bnnz;i++) {    
      cs = (bindx[i]-1)*lb;     
      bs = (bjndx[i]-1)*lb;        
 
      if ( cs == bs ) {
        continue;                  
      }
      pval = &val[i*mm];           
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
      for (l=0;l!=n;l++) {          
       for (jj=0;jj!=lb;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t =  pb[jj];            
           for (ii=0;ii!=lb;ii++) {      
             pc[ii] += t* (*pval++);        
           }
         } else {
           pval+=lb;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
 
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
      for (l=0;l!=n;l++) {         
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t -=  pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
  }
}

void BCOskew_MatMult_CATBbC_double(
                 const int mb, const int n, const int kb, 
                 const double *val, const int *bindx, const int *bjndx,
                 const int bnnz, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  int i,j,jb,je,block;
  int cs,css,cr,bs,br;
  int index;
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;
  const double *ptmp;
  int l;                                               

  for (l=0;l!=n;l++)                            
    for (i=0;i!=k;i++) *pc++ *= beta;                  

  for (i=0;i!=bnnz;i++) {    
      cs = (bindx[i]-1)*lb;     
      bs = (bjndx[i]-1)*lb;        
 
      if ( cs == bs ) {
        continue;                  
      }
      pval = &val[i*mm];           
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
      for (l=0;l!=n;l++) {             
       for (jj=0;jj!=lb;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t =  pb[jj];            
           for (ii=0;ii!=lb;ii++) {      
             pc[ii] -= t* (*pval++);        
           }
         } else {
           pval+=lb;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
 
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
      for (l=0;l!=n;l++) {         
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t +=  pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
  }
}




void BCO_MatMult_CaABbC_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *bindx, const int *bjndx,
                 const int bnnz, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  const double *ptmp;                                  
  int i,j,jb,je,block;
  int cs,bs,br;
  int index;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  for (l=0;l!=n;l++)                            
    for (i=0;i!=m;i++) *pc++ *= beta;                  

  for (i=0;i!=bnnz;i++) {      
      cs = (bindx[i]-1)*lb;     
      bs = (bjndx[i]-1)*lb;        
      pval = &val[i*mm];           
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
      for (l=0;l!=n;l++) {                      
       for (jj=0;jj!=lb;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t = alpha * pb[jj];            
           for (ii=0;ii!=lb;ii++) {      
             pc[ii] += t* (*pval++);        
           }
         } else {
           pval+=lb;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
  }
}

void BCO_MatMult_CaATBbC_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *bindx, const int *bjndx, 
                 const int bnnz, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  const double *ptmp;                                  
  int i,j,jb,je,block;
  int cs,bs,cr;
  int index;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  for (l=0;l!=n;l++)                            
    for (i=0;i!=k;i++) *pc++ *= beta;                  

  for (i=0;i!=bnnz;i++) {      
      cs = (bjndx[i]-1)*lb;     
      bs = (bindx[i]-1)*lb;        
      pval = &val[i*mm];           
      ptmp = pval;                  
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
      for (l=0;l!=n;l++) {              
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t += alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
  } }

void BCOsymm_MatMult_CaABbC_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *bindx, const int *bjndx, 
                 const int bnnz, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  const double *ptmp; 
  int i,j,jb,je,block;
  int cs,css,cr,bs,br;
  int index;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  for (l=0;l!=n;l++)                            
    for (i=0;i!=m;i++) *pc++ *= beta;                  

  for (i=0;i!=bnnz;i++) {    
      cs = (bindx[i]-1)*lb;     
      bs = (bjndx[i]-1)*lb;        
      pval = &val[i*mm];           
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
      for (l=0;l!=n;l++) {             
       for (jj=0;jj!=lb;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t = alpha * pb[jj];            
           for (ii=0;ii!=lb;ii++) {      
             pc[ii] += t* (*pval++);        
           }
         } else {
           pval+=lb;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
 
      if ( cs == bs ) {
        continue;                  
      }
 
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
      for (l=0;l!=n;l++) {         
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t += alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
  }
}

void BCOskew_MatMult_CaABbC_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *bindx, const int *bjndx, 
                 const int bnnz, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int cs,css,cr,bs,br;
  int index;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  for (l=0;l!=n;l++)                            
    for (i=0;i!=m;i++) *pc++ *= beta;                  

  for (i=0;i!=bnnz;i++) {    
      cs = (bindx[i]-1)*lb;     
      bs = (bjndx[i]-1)*lb;        
 
      if ( cs == bs ) {
        continue;                  
      }
      pval = &val[i*mm];           
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
      for (l=0;l!=n;l++) {          
       for (jj=0;jj!=lb;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t = alpha * pb[jj];            
           for (ii=0;ii!=lb;ii++) {      
             pc[ii] += t* (*pval++);        
           }
         } else {
           pval+=lb;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
 
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
      for (l=0;l!=n;l++) {         
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t -= alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
  }
}

void BCOskew_MatMult_CaATBbC_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *bindx, const int *bjndx,
                 const int bnnz, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  int i,j,jb,je,block;
  int cs,css,cr,bs,br;
  int index;
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;
  const double *ptmp;
  int l;                                               

  for (l=0;l!=n;l++)                            
    for (i=0;i!=k;i++) *pc++ *= beta;                  

  for (i=0;i!=bnnz;i++) {    
      cs = (bindx[i]-1)*lb;     
      bs = (bjndx[i]-1)*lb;        
 
      if ( cs == bs ) {
        continue;                  
      }
      pval = &val[i*mm];           
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
      for (l=0;l!=n;l++) {             
       for (jj=0;jj!=lb;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t = alpha * pb[jj];            
           for (ii=0;ii!=lb;ii++) {      
             pc[ii] -= t* (*pval++);        
           }
         } else {
           pval+=lb;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
 
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
      for (l=0;l!=n;l++) {         
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t += alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
  }
}

