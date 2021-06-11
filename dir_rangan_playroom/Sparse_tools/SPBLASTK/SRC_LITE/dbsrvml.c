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


/* Created:  Sat Jul 6 14:35:48 EDT 1996 */

#include "dbsrvml.h"



void BSR_VecMult_CAB_double(
                 const int mb,  const int kb, 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  int i,j,jb,je,block;
  int cs,bs,br;
  int index;
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  bindx-=ind_base;

    for (i=0;i!=m;i++) *pc++ = 0;                  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
      bs=(index-1)*lb;             
      pval = &val[(j-1)*mm];       
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
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
    }
  }
}

void BSR_VecMult_CATB_double(
                 const int mb,  const int kb, 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  int i,j,jb,je,block;
  int cs,bs,cr;
  int index;
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  bindx-=ind_base;

    for (i=0;i!=k;i++) *pc++ = 0;                  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
      cs=(index-1)*lb;             
      pval = &val[(j-1)*mm];       
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t +=  pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
    }
  }
}

void BSRsymm_VecMult_CAB_double(
                 const int mb,  const int kb, 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  bindx-=ind_base;

    for (i=0;i!=m;i++) *pc++ = 0;                  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
      bs=(index-1)*lb;          
      pval = &val[(j-1)*mm];       
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
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
 
      if ( index == i + ind_base ) {
        continue;                  
      }
 
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t +=  pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
    }
  }
}

void BSRskew_VecMult_CAB_double(
                 const int mb,  const int kb, 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  bindx-=ind_base;

    for (i=0;i!=m;i++) *pc++ = 0;                  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
 
      if ( index == i + ind_base ) {
        continue;                  
      }
      bs=(index-1)*lb;             
      pval = &val[(j-1)*mm];       
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
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
 
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t -=  pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
    }
  }
}

void BSRskew_VecMult_CATB_double(
                 const int mb,  const int kb, 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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

  bindx-=ind_base;

    for (i=0;i!=k;i++) *pc++ = 0;                  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
 
      if ( index == i + ind_base ) {
        continue;                  
      }
      bs=(index-1)*lb;             
      pval = &val[(j-1)*mm];       
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
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
 
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t +=  pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
    }
  }
}




void BSR_VecMult_CaAB_double(
                 const int mb,  const int kb, const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  int i,j,jb,je,block;
  int cs,bs,br;
  int index;
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  bindx-=ind_base;

    for (i=0;i!=m;i++) *pc++ = 0;                  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
      bs=(index-1)*lb;             
      pval = &val[(j-1)*mm];       
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
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
    }
  }
}

void BSR_VecMult_CaATB_double(
                 const int mb,  const int kb, const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  int i,j,jb,je,block;
  int cs,bs,cr;
  int index;
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  bindx-=ind_base;

    for (i=0;i!=k;i++) *pc++ = 0;                  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
      cs=(index-1)*lb;             
      pval = &val[(j-1)*mm];       
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t += alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
    }
  }
}

void BSRsymm_VecMult_CaAB_double(
                 const int mb,  const int kb, const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  bindx-=ind_base;

    for (i=0;i!=m;i++) *pc++ = 0;                  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
      bs=(index-1)*lb;          
      pval = &val[(j-1)*mm];       
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
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
 
      if ( index == i + ind_base ) {
        continue;                  
      }
 
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t += alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
    }
  }
}

void BSRskew_VecMult_CaAB_double(
                 const int mb,  const int kb, const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  bindx-=ind_base;

    for (i=0;i!=m;i++) *pc++ = 0;                  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
 
      if ( index == i + ind_base ) {
        continue;                  
      }
      bs=(index-1)*lb;             
      pval = &val[(j-1)*mm];       
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
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
 
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t -= alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
    }
  }
}

void BSRskew_VecMult_CaATB_double(
                 const int mb,  const int kb, const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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

  bindx-=ind_base;

    for (i=0;i!=k;i++) *pc++ = 0;                  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
 
      if ( index == i + ind_base ) {
        continue;                  
      }
      bs=(index-1)*lb;             
      pval = &val[(j-1)*mm];       
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
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
 
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t += alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
    }
  }
}




void BSR_VecMult_CABC_double(
                 const int mb,  const int kb, 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  int i,j,jb,je,block;
  int cs,bs,br;
  int index;
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  bindx-=ind_base;


  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
      bs=(index-1)*lb;             
      pval = &val[(j-1)*mm];       
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
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
    }
  }
}

void BSR_VecMult_CATBC_double(
                 const int mb,  const int kb, 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  int i,j,jb,je,block;
  int cs,bs,cr;
  int index;
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  bindx-=ind_base;


  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
      cs=(index-1)*lb;             
      pval = &val[(j-1)*mm];       
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t +=  pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
    }
  }
}

void BSRsymm_VecMult_CABC_double(
                 const int mb,  const int kb, 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  bindx-=ind_base;


  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
      bs=(index-1)*lb;          
      pval = &val[(j-1)*mm];       
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
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
 
      if ( index == i + ind_base ) {
        continue;                  
      }
 
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t +=  pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
    }
  }
}

void BSRskew_VecMult_CABC_double(
                 const int mb,  const int kb, 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  bindx-=ind_base;


  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
 
      if ( index == i + ind_base ) {
        continue;                  
      }
      bs=(index-1)*lb;             
      pval = &val[(j-1)*mm];       
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
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
 
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t -=  pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
    }
  }
}

void BSRskew_VecMult_CATBC_double(
                 const int mb,  const int kb, 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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

  bindx-=ind_base;


  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
 
      if ( index == i + ind_base ) {
        continue;                  
      }
      bs=(index-1)*lb;             
      pval = &val[(j-1)*mm];       
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
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
 
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t +=  pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
    }
  }
}




void BSR_VecMult_CaABC_double(
                 const int mb,  const int kb, const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  int i,j,jb,je,block;
  int cs,bs,br;
  int index;
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  bindx-=ind_base;


  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
      bs=(index-1)*lb;             
      pval = &val[(j-1)*mm];       
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
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
    }
  }
}

void BSR_VecMult_CaATBC_double(
                 const int mb,  const int kb, const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  int i,j,jb,je,block;
  int cs,bs,cr;
  int index;
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  bindx-=ind_base;


  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
      cs=(index-1)*lb;             
      pval = &val[(j-1)*mm];       
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t += alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
    }
  }
}

void BSRsymm_VecMult_CaABC_double(
                 const int mb,  const int kb, const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  bindx-=ind_base;


  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
      bs=(index-1)*lb;          
      pval = &val[(j-1)*mm];       
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
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
 
      if ( index == i + ind_base ) {
        continue;                  
      }
 
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t += alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
    }
  }
}

void BSRskew_VecMult_CaABC_double(
                 const int mb,  const int kb, const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  bindx-=ind_base;


  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
 
      if ( index == i + ind_base ) {
        continue;                  
      }
      bs=(index-1)*lb;             
      pval = &val[(j-1)*mm];       
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
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
 
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t -= alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
    }
  }
}

void BSRskew_VecMult_CaATBC_double(
                 const int mb,  const int kb, const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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

  bindx-=ind_base;


  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
 
      if ( index == i + ind_base ) {
        continue;                  
      }
      bs=(index-1)*lb;             
      pval = &val[(j-1)*mm];       
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
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
 
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t += alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
    }
  }
}




void BSR_VecMult_CABbC_double(
                 const int mb,  const int kb, 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  int i,j,jb,je,block;
  int cs,bs,br;
  int index;
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  bindx-=ind_base;

    for (i=0;i!=m;i++) *pc++ *= beta;                  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
      bs=(index-1)*lb;             
      pval = &val[(j-1)*mm];       
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
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
    }
  }
}

void BSR_VecMult_CATBbC_double(
                 const int mb,  const int kb, 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  int i,j,jb,je,block;
  int cs,bs,cr;
  int index;
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  bindx-=ind_base;

    for (i=0;i!=k;i++) *pc++ *= beta;                  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
      cs=(index-1)*lb;             
      pval = &val[(j-1)*mm];       
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t +=  pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
    }
  }
}

void BSRsymm_VecMult_CABbC_double(
                 const int mb,  const int kb, 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  bindx-=ind_base;

    for (i=0;i!=m;i++) *pc++ *= beta;                  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
      bs=(index-1)*lb;          
      pval = &val[(j-1)*mm];       
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
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
 
      if ( index == i + ind_base ) {
        continue;                  
      }
 
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t +=  pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
    }
  }
}

void BSRskew_VecMult_CABbC_double(
                 const int mb,  const int kb, 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  bindx-=ind_base;

    for (i=0;i!=m;i++) *pc++ *= beta;                  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
 
      if ( index == i + ind_base ) {
        continue;                  
      }
      bs=(index-1)*lb;             
      pval = &val[(j-1)*mm];       
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
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
 
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t -=  pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
    }
  }
}

void BSRskew_VecMult_CATBbC_double(
                 const int mb,  const int kb, 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
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

  bindx-=ind_base;

    for (i=0;i!=k;i++) *pc++ *= beta;                  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
 
      if ( index == i + ind_base ) {
        continue;                  
      }
      bs=(index-1)*lb;             
      pval = &val[(j-1)*mm];       
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
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
 
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t +=  pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
    }
  }
}




void BSR_VecMult_CaABbC_double(
                 const int mb,  const int kb, const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  int i,j,jb,je,block;
  int cs,bs,br;
  int index;
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  bindx-=ind_base;

    for (i=0;i!=m;i++) *pc++ *= beta;                  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
      bs=(index-1)*lb;             
      pval = &val[(j-1)*mm];       
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
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
    }
  }
}

void BSR_VecMult_CaATBbC_double(
                 const int mb,  const int kb, const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
                 const int ind_base)
{
  double t;
  const double *pb;
  double *pc=c;
  const double *pval;
  int i,j,jb,je,block;
  int cs,bs,cr;
  int index;
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  bindx-=ind_base;

    for (i=0;i!=k;i++) *pc++ *= beta;                  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
      cs=(index-1)*lb;             
      pval = &val[(j-1)*mm];       
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t += alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
    }
  }
}

void BSRsymm_VecMult_CaABbC_double(
                 const int mb,  const int kb, const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  bindx-=ind_base;

    for (i=0;i!=m;i++) *pc++ *= beta;                  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
      bs=(index-1)*lb;          
      pval = &val[(j-1)*mm];       
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
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
 
      if ( index == i + ind_base ) {
        continue;                  
      }
 
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t += alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
    }
  }
}

void BSRskew_VecMult_CaABbC_double(
                 const int mb,  const int kb, const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int k=kb*lb;
  int mm=lb*lb;

  bindx-=ind_base;

    for (i=0;i!=m;i++) *pc++ *= beta;                  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
 
      if ( index == i + ind_base ) {
        continue;                  
      }
      bs=(index-1)*lb;             
      pval = &val[(j-1)*mm];       
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
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
 
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t -= alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
    }
  }
}

void BSRskew_VecMult_CaATBbC_double(
                 const int mb,  const int kb, const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
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

  bindx-=ind_base;

    for (i=0;i!=k;i++) *pc++ *= beta;                  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  
    for (j=jb;j!=je;j++) {      
      index = bindx[j];
 
      if ( index == i + ind_base ) {
        continue;                  
      }
      bs=(index-1)*lb;             
      pval = &val[(j-1)*mm];       
      ptmp = pval;                 
                                   
      pb = &b[bs];                 
      pc = &c[cs];                 
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
 
      pval = ptmp;                 
                                   
      pb = &b[cs];                 
      pc = &c[bs];                 
       for (jj=0;jj!=lb;jj++) {       
           t = 0;
           for (ii=0;ii!=lb;ii++) {      
             t += alpha * pb[ii] * (*pval++);
           }
           pc[jj] += t;                  
       }
    }
  }
}

