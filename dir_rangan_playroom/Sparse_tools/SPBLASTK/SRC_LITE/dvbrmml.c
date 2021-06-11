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


/* Created:  Sat Jul 6 14:38:12 EDT 1996 */

#include "dvbrmml.h"



void VBR_MatMult_CAB_double(
                 const int mb, const int n, const int kb, 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int cs,mm,bs,br;
  int index;
  int l;                                               
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  int k=cpntr[kb]-cpntr[0];

  val-=ind_base;
  b-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  for (l=0;l!=n;l++)                            
    for (i=0;i!=m;i++) *pc++ = 0;                  

  c-=ind_base;
 
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
      bs=cpntr[index];             
      br=cpntr[index+1]-bs;        
      pval = &val[indx[j]];        
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[rpntr[i]];           
      for (l=0;l!=n;l++) {                      
       for (jj=0;jj!=br;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t =  pb[jj];            
           for (ii=0;ii!=mm;ii++) {      
             pc[ii] += t* (*pval++);        
           }
         } else {
           pval+=mm;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
    }
  }
}

void VBR_MatMult_CATB_double(
                 const int mb, const int n, const int kb, 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int cs,mm,bs,cr;
  int index;
  int l;                                               
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  int k=cpntr[kb]-cpntr[0];

  val-=ind_base;
  b-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  for (l=0;l!=n;l++)                            
    for (i=0;i!=k;i++) *pc++ = 0;                  

  c-=ind_base;
 
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
      cs=cpntr[index];             
      cr=cpntr[index+1]-cs;        
      pval = &val[indx[j]];        
      ptmp = pval;                  
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
      for (l=0;l!=n;l++) {              
       for (jj=0;jj!=cr;jj++) {       
           t = 0;
           for (ii=0;ii!=mm;ii++) {      
             t +=  pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
    }
  }
}

void VBRsymm_MatMult_CAB_double(
                 const int mb, const int n, const int kb, 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int cs,css,cr,mm,bs,br;
  int index;
  int l;                                               
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  int k=cpntr[kb]-cpntr[0];

  val-=ind_base;
  b-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  for (l=0;l!=n;l++)                            
    for (i=0;i!=m;i++) *pc++ = 0;                  

  c-=ind_base;
 
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
      bs=cpntr[index];             
      br=cpntr[index+1]-bs;        
      pval = &val[indx[j]];        
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[rpntr[i]];           
      for (l=0;l!=n;l++) {             
       for (jj=0;jj!=br;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t =  pb[jj];            
           for (ii=0;ii!=mm;ii++) {      
             pc[ii] += t* (*pval++);        
           }
         } else {
           pval+=mm;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
 
      if ( index == i + ind_base ) {
        continue;                  
      }
 
      cr= br;                      
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
      for (l=0;l!=n;l++) {         
       for (jj=0;jj!=cr;jj++) {       
           t = 0;
           for (ii=0;ii!=mm;ii++) {      
             t +=  pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
    }
  }
}

void VBRskew_MatMult_CAB_double(
                 const int mb, const int n, const int kb, 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int cs,css,cr,mm,bs,br;
  int index;
  int l;                                               
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  int k=cpntr[kb]-cpntr[0];

  val-=ind_base;
  b-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  for (l=0;l!=n;l++)                            
    for (i=0;i!=m;i++) *pc++ = 0;                  

  c-=ind_base;
 
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
 
      if ( index == i + ind_base ) {
        continue;                  
      }
      bs=cpntr[index];             
      br=cpntr[index+1]-bs;        
      pval = &val[indx[j]];        
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[rpntr[i]];           
      for (l=0;l!=n;l++) {          
       for (jj=0;jj!=br;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t =  pb[jj];            
           for (ii=0;ii!=mm;ii++) {      
             pc[ii] += t* (*pval++);        
           }
         } else {
           pval+=mm;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
 
      cr= br;                      
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
      for (l=0;l!=n;l++) {         
       for (jj=0;jj!=cr;jj++) {       
           t = 0;
           for (ii=0;ii!=mm;ii++) {      
             t -=  pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
    }
  }
}

void VBRskew_MatMult_CATB_double(
                 const int mb, const int n, const int kb, 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  int i,j,jb,je,block;
  int cs,css,cr,mm,bs,br;
  int index;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  int k=cpntr[kb]-cpntr[0];
  const double *ptmp;
  int l;                                               

  val-=ind_base;
  b-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  for (l=0;l!=n;l++)                            
    for (i=0;i!=k;i++) *pc++ = 0;                  

  c-=ind_base;
 
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
 
      if ( index == i + ind_base ) {
        continue;                  
      }
      bs=cpntr[index];             
      br=cpntr[index+1]-bs;        
      pval = &val[indx[j]];        
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[rpntr[i]];           
      for (l=0;l!=n;l++) {             
       for (jj=0;jj!=br;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t =  pb[jj];            
           for (ii=0;ii!=mm;ii++) {      
             pc[ii] -= t* (*pval++);        
           }
         } else {
           pval+=mm;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
 
      cr= br;                      
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
      for (l=0;l!=n;l++) {         
       for (jj=0;jj!=cr;jj++) {       
           t = 0;
           for (ii=0;ii!=mm;ii++) {      
             t +=  pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
    }
  }
}




void VBR_MatMult_CaAB_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int cs,mm,bs,br;
  int index;
  int l;                                               
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  int k=cpntr[kb]-cpntr[0];

  val-=ind_base;
  b-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  for (l=0;l!=n;l++)                            
    for (i=0;i!=m;i++) *pc++ = 0;                  

  c-=ind_base;
 
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
      bs=cpntr[index];             
      br=cpntr[index+1]-bs;        
      pval = &val[indx[j]];        
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[rpntr[i]];           
      for (l=0;l!=n;l++) {                      
       for (jj=0;jj!=br;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t = alpha * pb[jj];            
           for (ii=0;ii!=mm;ii++) {      
             pc[ii] += t* (*pval++);        
           }
         } else {
           pval+=mm;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
    }
  }
}

void VBR_MatMult_CaATB_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int cs,mm,bs,cr;
  int index;
  int l;                                               
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  int k=cpntr[kb]-cpntr[0];

  val-=ind_base;
  b-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  for (l=0;l!=n;l++)                            
    for (i=0;i!=k;i++) *pc++ = 0;                  

  c-=ind_base;
 
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
      cs=cpntr[index];             
      cr=cpntr[index+1]-cs;        
      pval = &val[indx[j]];        
      ptmp = pval;                  
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
      for (l=0;l!=n;l++) {              
       for (jj=0;jj!=cr;jj++) {       
           t = 0;
           for (ii=0;ii!=mm;ii++) {      
             t += alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
    }
  }
}

void VBRsymm_MatMult_CaAB_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int cs,css,cr,mm,bs,br;
  int index;
  int l;                                               
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  int k=cpntr[kb]-cpntr[0];

  val-=ind_base;
  b-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  for (l=0;l!=n;l++)                            
    for (i=0;i!=m;i++) *pc++ = 0;                  

  c-=ind_base;
 
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
      bs=cpntr[index];             
      br=cpntr[index+1]-bs;        
      pval = &val[indx[j]];        
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[rpntr[i]];           
      for (l=0;l!=n;l++) {             
       for (jj=0;jj!=br;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t = alpha * pb[jj];            
           for (ii=0;ii!=mm;ii++) {      
             pc[ii] += t* (*pval++);        
           }
         } else {
           pval+=mm;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
 
      if ( index == i + ind_base ) {
        continue;                  
      }
 
      cr= br;                      
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
      for (l=0;l!=n;l++) {         
       for (jj=0;jj!=cr;jj++) {       
           t = 0;
           for (ii=0;ii!=mm;ii++) {      
             t += alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
    }
  }
}

void VBRskew_MatMult_CaAB_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int cs,css,cr,mm,bs,br;
  int index;
  int l;                                               
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  int k=cpntr[kb]-cpntr[0];

  val-=ind_base;
  b-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  for (l=0;l!=n;l++)                            
    for (i=0;i!=m;i++) *pc++ = 0;                  

  c-=ind_base;
 
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
 
      if ( index == i + ind_base ) {
        continue;                  
      }
      bs=cpntr[index];             
      br=cpntr[index+1]-bs;        
      pval = &val[indx[j]];        
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[rpntr[i]];           
      for (l=0;l!=n;l++) {          
       for (jj=0;jj!=br;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t = alpha * pb[jj];            
           for (ii=0;ii!=mm;ii++) {      
             pc[ii] += t* (*pval++);        
           }
         } else {
           pval+=mm;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
 
      cr= br;                      
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
      for (l=0;l!=n;l++) {         
       for (jj=0;jj!=cr;jj++) {       
           t = 0;
           for (ii=0;ii!=mm;ii++) {      
             t -= alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
    }
  }
}

void VBRskew_MatMult_CaATB_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  int i,j,jb,je,block;
  int cs,css,cr,mm,bs,br;
  int index;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  int k=cpntr[kb]-cpntr[0];
  const double *ptmp;
  int l;                                               

  val-=ind_base;
  b-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  for (l=0;l!=n;l++)                            
    for (i=0;i!=k;i++) *pc++ = 0;                  

  c-=ind_base;
 
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
 
      if ( index == i + ind_base ) {
        continue;                  
      }
      bs=cpntr[index];             
      br=cpntr[index+1]-bs;        
      pval = &val[indx[j]];        
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[rpntr[i]];           
      for (l=0;l!=n;l++) {             
       for (jj=0;jj!=br;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t = alpha * pb[jj];            
           for (ii=0;ii!=mm;ii++) {      
             pc[ii] -= t* (*pval++);        
           }
         } else {
           pval+=mm;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
 
      cr= br;                      
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
      for (l=0;l!=n;l++) {         
       for (jj=0;jj!=cr;jj++) {       
           t = 0;
           for (ii=0;ii!=mm;ii++) {      
             t += alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
    }
  }
}




void VBR_MatMult_CABC_double(
                 const int mb, const int n, const int kb, 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int cs,mm,bs,br;
  int index;
  int l;                                               
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  int k=cpntr[kb]-cpntr[0];

  val-=ind_base;
  b-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;


  c-=ind_base;
 
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
      bs=cpntr[index];             
      br=cpntr[index+1]-bs;        
      pval = &val[indx[j]];        
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[rpntr[i]];           
      for (l=0;l!=n;l++) {                      
       for (jj=0;jj!=br;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t =  pb[jj];            
           for (ii=0;ii!=mm;ii++) {      
             pc[ii] += t* (*pval++);        
           }
         } else {
           pval+=mm;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
    }
  }
}

void VBR_MatMult_CATBC_double(
                 const int mb, const int n, const int kb, 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int cs,mm,bs,cr;
  int index;
  int l;                                               
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  int k=cpntr[kb]-cpntr[0];

  val-=ind_base;
  b-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;


  c-=ind_base;
 
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
      cs=cpntr[index];             
      cr=cpntr[index+1]-cs;        
      pval = &val[indx[j]];        
      ptmp = pval;                  
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
      for (l=0;l!=n;l++) {              
       for (jj=0;jj!=cr;jj++) {       
           t = 0;
           for (ii=0;ii!=mm;ii++) {      
             t +=  pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
    }
  }
}

void VBRsymm_MatMult_CABC_double(
                 const int mb, const int n, const int kb, 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int cs,css,cr,mm,bs,br;
  int index;
  int l;                                               
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  int k=cpntr[kb]-cpntr[0];

  val-=ind_base;
  b-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;


  c-=ind_base;
 
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
      bs=cpntr[index];             
      br=cpntr[index+1]-bs;        
      pval = &val[indx[j]];        
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[rpntr[i]];           
      for (l=0;l!=n;l++) {             
       for (jj=0;jj!=br;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t =  pb[jj];            
           for (ii=0;ii!=mm;ii++) {      
             pc[ii] += t* (*pval++);        
           }
         } else {
           pval+=mm;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
 
      if ( index == i + ind_base ) {
        continue;                  
      }
 
      cr= br;                      
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
      for (l=0;l!=n;l++) {         
       for (jj=0;jj!=cr;jj++) {       
           t = 0;
           for (ii=0;ii!=mm;ii++) {      
             t +=  pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
    }
  }
}

void VBRskew_MatMult_CABC_double(
                 const int mb, const int n, const int kb, 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int cs,css,cr,mm,bs,br;
  int index;
  int l;                                               
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  int k=cpntr[kb]-cpntr[0];

  val-=ind_base;
  b-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;


  c-=ind_base;
 
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
 
      if ( index == i + ind_base ) {
        continue;                  
      }
      bs=cpntr[index];             
      br=cpntr[index+1]-bs;        
      pval = &val[indx[j]];        
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[rpntr[i]];           
      for (l=0;l!=n;l++) {          
       for (jj=0;jj!=br;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t =  pb[jj];            
           for (ii=0;ii!=mm;ii++) {      
             pc[ii] += t* (*pval++);        
           }
         } else {
           pval+=mm;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
 
      cr= br;                      
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
      for (l=0;l!=n;l++) {         
       for (jj=0;jj!=cr;jj++) {       
           t = 0;
           for (ii=0;ii!=mm;ii++) {      
             t -=  pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
    }
  }
}

void VBRskew_MatMult_CATBC_double(
                 const int mb, const int n, const int kb, 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  int i,j,jb,je,block;
  int cs,css,cr,mm,bs,br;
  int index;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  int k=cpntr[kb]-cpntr[0];
  const double *ptmp;
  int l;                                               

  val-=ind_base;
  b-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;


  c-=ind_base;
 
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
 
      if ( index == i + ind_base ) {
        continue;                  
      }
      bs=cpntr[index];             
      br=cpntr[index+1]-bs;        
      pval = &val[indx[j]];        
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[rpntr[i]];           
      for (l=0;l!=n;l++) {             
       for (jj=0;jj!=br;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t =  pb[jj];            
           for (ii=0;ii!=mm;ii++) {      
             pc[ii] -= t* (*pval++);        
           }
         } else {
           pval+=mm;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
 
      cr= br;                      
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
      for (l=0;l!=n;l++) {         
       for (jj=0;jj!=cr;jj++) {       
           t = 0;
           for (ii=0;ii!=mm;ii++) {      
             t +=  pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
    }
  }
}




void VBR_MatMult_CaABC_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int cs,mm,bs,br;
  int index;
  int l;                                               
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  int k=cpntr[kb]-cpntr[0];

  val-=ind_base;
  b-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;


  c-=ind_base;
 
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
      bs=cpntr[index];             
      br=cpntr[index+1]-bs;        
      pval = &val[indx[j]];        
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[rpntr[i]];           
      for (l=0;l!=n;l++) {                      
       for (jj=0;jj!=br;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t = alpha * pb[jj];            
           for (ii=0;ii!=mm;ii++) {      
             pc[ii] += t* (*pval++);        
           }
         } else {
           pval+=mm;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
    }
  }
}

void VBR_MatMult_CaATBC_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int cs,mm,bs,cr;
  int index;
  int l;                                               
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  int k=cpntr[kb]-cpntr[0];

  val-=ind_base;
  b-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;


  c-=ind_base;
 
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
      cs=cpntr[index];             
      cr=cpntr[index+1]-cs;        
      pval = &val[indx[j]];        
      ptmp = pval;                  
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
      for (l=0;l!=n;l++) {              
       for (jj=0;jj!=cr;jj++) {       
           t = 0;
           for (ii=0;ii!=mm;ii++) {      
             t += alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
    }
  }
}

void VBRsymm_MatMult_CaABC_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int cs,css,cr,mm,bs,br;
  int index;
  int l;                                               
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  int k=cpntr[kb]-cpntr[0];

  val-=ind_base;
  b-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;


  c-=ind_base;
 
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
      bs=cpntr[index];             
      br=cpntr[index+1]-bs;        
      pval = &val[indx[j]];        
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[rpntr[i]];           
      for (l=0;l!=n;l++) {             
       for (jj=0;jj!=br;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t = alpha * pb[jj];            
           for (ii=0;ii!=mm;ii++) {      
             pc[ii] += t* (*pval++);        
           }
         } else {
           pval+=mm;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
 
      if ( index == i + ind_base ) {
        continue;                  
      }
 
      cr= br;                      
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
      for (l=0;l!=n;l++) {         
       for (jj=0;jj!=cr;jj++) {       
           t = 0;
           for (ii=0;ii!=mm;ii++) {      
             t += alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
    }
  }
}

void VBRskew_MatMult_CaABC_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int cs,css,cr,mm,bs,br;
  int index;
  int l;                                               
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  int k=cpntr[kb]-cpntr[0];

  val-=ind_base;
  b-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;


  c-=ind_base;
 
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
 
      if ( index == i + ind_base ) {
        continue;                  
      }
      bs=cpntr[index];             
      br=cpntr[index+1]-bs;        
      pval = &val[indx[j]];        
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[rpntr[i]];           
      for (l=0;l!=n;l++) {          
       for (jj=0;jj!=br;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t = alpha * pb[jj];            
           for (ii=0;ii!=mm;ii++) {      
             pc[ii] += t* (*pval++);        
           }
         } else {
           pval+=mm;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
 
      cr= br;                      
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
      for (l=0;l!=n;l++) {         
       for (jj=0;jj!=cr;jj++) {       
           t = 0;
           for (ii=0;ii!=mm;ii++) {      
             t -= alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
    }
  }
}

void VBRskew_MatMult_CaATBC_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  int i,j,jb,je,block;
  int cs,css,cr,mm,bs,br;
  int index;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  int k=cpntr[kb]-cpntr[0];
  const double *ptmp;
  int l;                                               

  val-=ind_base;
  b-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;


  c-=ind_base;
 
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
 
      if ( index == i + ind_base ) {
        continue;                  
      }
      bs=cpntr[index];             
      br=cpntr[index+1]-bs;        
      pval = &val[indx[j]];        
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[rpntr[i]];           
      for (l=0;l!=n;l++) {             
       for (jj=0;jj!=br;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t = alpha * pb[jj];            
           for (ii=0;ii!=mm;ii++) {      
             pc[ii] -= t* (*pval++);        
           }
         } else {
           pval+=mm;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
 
      cr= br;                      
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
      for (l=0;l!=n;l++) {         
       for (jj=0;jj!=cr;jj++) {       
           t = 0;
           for (ii=0;ii!=mm;ii++) {      
             t += alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
    }
  }
}




void VBR_MatMult_CABbC_double(
                 const int mb, const int n, const int kb, 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int cs,mm,bs,br;
  int index;
  int l;                                               
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  int k=cpntr[kb]-cpntr[0];

  val-=ind_base;
  b-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  for (l=0;l!=n;l++)                            
    for (i=0;i!=m;i++) *pc++ *= beta;                  

  c-=ind_base;
 
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
      bs=cpntr[index];             
      br=cpntr[index+1]-bs;        
      pval = &val[indx[j]];        
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[rpntr[i]];           
      for (l=0;l!=n;l++) {                      
       for (jj=0;jj!=br;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t =  pb[jj];            
           for (ii=0;ii!=mm;ii++) {      
             pc[ii] += t* (*pval++);        
           }
         } else {
           pval+=mm;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
    }
  }
}

void VBR_MatMult_CATBbC_double(
                 const int mb, const int n, const int kb, 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int cs,mm,bs,cr;
  int index;
  int l;                                               
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  int k=cpntr[kb]-cpntr[0];

  val-=ind_base;
  b-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  for (l=0;l!=n;l++)                            
    for (i=0;i!=k;i++) *pc++ *= beta;                  

  c-=ind_base;
 
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
      cs=cpntr[index];             
      cr=cpntr[index+1]-cs;        
      pval = &val[indx[j]];        
      ptmp = pval;                  
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
      for (l=0;l!=n;l++) {              
       for (jj=0;jj!=cr;jj++) {       
           t = 0;
           for (ii=0;ii!=mm;ii++) {      
             t +=  pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
    }
  }
}

void VBRsymm_MatMult_CABbC_double(
                 const int mb, const int n, const int kb, 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int cs,css,cr,mm,bs,br;
  int index;
  int l;                                               
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  int k=cpntr[kb]-cpntr[0];

  val-=ind_base;
  b-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  for (l=0;l!=n;l++)                            
    for (i=0;i!=m;i++) *pc++ *= beta;                  

  c-=ind_base;
 
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
      bs=cpntr[index];             
      br=cpntr[index+1]-bs;        
      pval = &val[indx[j]];        
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[rpntr[i]];           
      for (l=0;l!=n;l++) {             
       for (jj=0;jj!=br;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t =  pb[jj];            
           for (ii=0;ii!=mm;ii++) {      
             pc[ii] += t* (*pval++);        
           }
         } else {
           pval+=mm;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
 
      if ( index == i + ind_base ) {
        continue;                  
      }
 
      cr= br;                      
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
      for (l=0;l!=n;l++) {         
       for (jj=0;jj!=cr;jj++) {       
           t = 0;
           for (ii=0;ii!=mm;ii++) {      
             t +=  pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
    }
  }
}

void VBRskew_MatMult_CABbC_double(
                 const int mb, const int n, const int kb, 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int cs,css,cr,mm,bs,br;
  int index;
  int l;                                               
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  int k=cpntr[kb]-cpntr[0];

  val-=ind_base;
  b-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  for (l=0;l!=n;l++)                            
    for (i=0;i!=m;i++) *pc++ *= beta;                  

  c-=ind_base;
 
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
 
      if ( index == i + ind_base ) {
        continue;                  
      }
      bs=cpntr[index];             
      br=cpntr[index+1]-bs;        
      pval = &val[indx[j]];        
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[rpntr[i]];           
      for (l=0;l!=n;l++) {          
       for (jj=0;jj!=br;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t =  pb[jj];            
           for (ii=0;ii!=mm;ii++) {      
             pc[ii] += t* (*pval++);        
           }
         } else {
           pval+=mm;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
 
      cr= br;                      
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
      for (l=0;l!=n;l++) {         
       for (jj=0;jj!=cr;jj++) {       
           t = 0;
           for (ii=0;ii!=mm;ii++) {      
             t -=  pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
    }
  }
}

void VBRskew_MatMult_CATBbC_double(
                 const int mb, const int n, const int kb, 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  int i,j,jb,je,block;
  int cs,css,cr,mm,bs,br;
  int index;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  int k=cpntr[kb]-cpntr[0];
  const double *ptmp;
  int l;                                               

  val-=ind_base;
  b-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  for (l=0;l!=n;l++)                            
    for (i=0;i!=k;i++) *pc++ *= beta;                  

  c-=ind_base;
 
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
 
      if ( index == i + ind_base ) {
        continue;                  
      }
      bs=cpntr[index];             
      br=cpntr[index+1]-bs;        
      pval = &val[indx[j]];        
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[rpntr[i]];           
      for (l=0;l!=n;l++) {             
       for (jj=0;jj!=br;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t =  pb[jj];            
           for (ii=0;ii!=mm;ii++) {      
             pc[ii] -= t* (*pval++);        
           }
         } else {
           pval+=mm;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
 
      cr= br;                      
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
      for (l=0;l!=n;l++) {         
       for (jj=0;jj!=cr;jj++) {       
           t = 0;
           for (ii=0;ii!=mm;ii++) {      
             t +=  pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
    }
  }
}




void VBR_MatMult_CaABbC_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int cs,mm,bs,br;
  int index;
  int l;                                               
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  int k=cpntr[kb]-cpntr[0];

  val-=ind_base;
  b-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  for (l=0;l!=n;l++)                            
    for (i=0;i!=m;i++) *pc++ *= beta;                  

  c-=ind_base;
 
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
      bs=cpntr[index];             
      br=cpntr[index+1]-bs;        
      pval = &val[indx[j]];        
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[rpntr[i]];           
      for (l=0;l!=n;l++) {                      
       for (jj=0;jj!=br;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t = alpha * pb[jj];            
           for (ii=0;ii!=mm;ii++) {      
             pc[ii] += t* (*pval++);        
           }
         } else {
           pval+=mm;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
    }
  }
}

void VBR_MatMult_CaATBbC_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int cs,mm,bs,cr;
  int index;
  int l;                                               
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  int k=cpntr[kb]-cpntr[0];

  val-=ind_base;
  b-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  for (l=0;l!=n;l++)                            
    for (i=0;i!=k;i++) *pc++ *= beta;                  

  c-=ind_base;
 
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
      cs=cpntr[index];             
      cr=cpntr[index+1]-cs;        
      pval = &val[indx[j]];        
      ptmp = pval;                  
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
      for (l=0;l!=n;l++) {              
       for (jj=0;jj!=cr;jj++) {       
           t = 0;
           for (ii=0;ii!=mm;ii++) {      
             t += alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
    }
  }
}

void VBRsymm_MatMult_CaABbC_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int cs,css,cr,mm,bs,br;
  int index;
  int l;                                               
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  int k=cpntr[kb]-cpntr[0];

  val-=ind_base;
  b-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  for (l=0;l!=n;l++)                            
    for (i=0;i!=m;i++) *pc++ *= beta;                  

  c-=ind_base;
 
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
      bs=cpntr[index];             
      br=cpntr[index+1]-bs;        
      pval = &val[indx[j]];        
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[rpntr[i]];           
      for (l=0;l!=n;l++) {             
       for (jj=0;jj!=br;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t = alpha * pb[jj];            
           for (ii=0;ii!=mm;ii++) {      
             pc[ii] += t* (*pval++);        
           }
         } else {
           pval+=mm;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
 
      if ( index == i + ind_base ) {
        continue;                  
      }
 
      cr= br;                      
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
      for (l=0;l!=n;l++) {         
       for (jj=0;jj!=cr;jj++) {       
           t = 0;
           for (ii=0;ii!=mm;ii++) {      
             t += alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
    }
  }
}

void VBRskew_MatMult_CaABbC_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int cs,css,cr,mm,bs,br;
  int index;
  int l;                                               
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  int k=cpntr[kb]-cpntr[0];

  val-=ind_base;
  b-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  for (l=0;l!=n;l++)                            
    for (i=0;i!=m;i++) *pc++ *= beta;                  

  c-=ind_base;
 
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
 
      if ( index == i + ind_base ) {
        continue;                  
      }
      bs=cpntr[index];             
      br=cpntr[index+1]-bs;        
      pval = &val[indx[j]];        
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[rpntr[i]];           
      for (l=0;l!=n;l++) {          
       for (jj=0;jj!=br;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t = alpha * pb[jj];            
           for (ii=0;ii!=mm;ii++) {      
             pc[ii] += t* (*pval++);        
           }
         } else {
           pval+=mm;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
 
      cr= br;                      
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
      for (l=0;l!=n;l++) {         
       for (jj=0;jj!=cr;jj++) {       
           t = 0;
           for (ii=0;ii!=mm;ii++) {      
             t -= alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
    }
  }
}

void VBRskew_MatMult_CaATBbC_double(
                 const int mb, const int n, const int kb, const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  int i,j,jb,je,block;
  int cs,css,cr,mm,bs,br;
  int index;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  int k=cpntr[kb]-cpntr[0];
  const double *ptmp;
  int l;                                               

  val-=ind_base;
  b-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  for (l=0;l!=n;l++)                            
    for (i=0;i!=k;i++) *pc++ *= beta;                  

  c-=ind_base;
 
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
 
      if ( index == i + ind_base ) {
        continue;                  
      }
      bs=cpntr[index];             
      br=cpntr[index+1]-bs;        
      pval = &val[indx[j]];        
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[rpntr[i]];           
      for (l=0;l!=n;l++) {             
       for (jj=0;jj!=br;jj++) {       
         if( pb[jj] != 0.0 ) {           
           t = alpha * pb[jj];            
           for (ii=0;ii!=mm;ii++) {      
             pc[ii] -= t* (*pval++);        
           }
         } else {
           pval+=mm;                     
         }
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
 
      cr= br;                      
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
      for (l=0;l!=n;l++) {         
       for (jj=0;jj!=cr;jj++) {       
           t = 0;
           for (ii=0;ii!=mm;ii++) {      
             t += alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
       pb+=ldb; pc+=ldc; pval=ptmp;      
      }                                               
    }
  }
}

