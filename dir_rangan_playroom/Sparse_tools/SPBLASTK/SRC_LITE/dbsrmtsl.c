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


/* Created:  Sat Jul 6 14:35:49 EDT 1996 */

#include "dbsrmtsl.h"



void BSR_MatTriangSlvLU_CAB_double(
                 const int mb, const int n, 
                  
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                  const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;



  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        



}


void BSR_MatTriangSlvUU_CAB_double(
                 const int mb, const int n,
                  
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                  const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;



  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        




}

void BSR_MatTriangSlvLU_CATB_double(
                 const int mb, const int n,
                  
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                  const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;



  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        




}

void BSR_MatTriangSlvUU_CATB_double(
                 const int mb, const int n, 
                  
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                  const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;



  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        




}



void BSR_MatTriangSlvLU_CaAB_double(
                 const int mb, const int n, 
                  
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                  const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;



  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  

  pc=c;                                   
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     

}


void BSR_MatTriangSlvUU_CaAB_double(
                 const int mb, const int n,
                  
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                  const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;



  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     


}

void BSR_MatTriangSlvLU_CaATB_double(
                 const int mb, const int n,
                  
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                  const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;



  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     


}

void BSR_MatTriangSlvUU_CaATB_double(
                 const int mb, const int n, 
                  
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                  const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;



  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     


}



void BSR_MatTriangSlvLU_CABC_double(
                 const int mb, const int n, 
                  
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void BSR_MatTriangSlvUU_CABC_double(
                 const int mb, const int n,
                  
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                            
  double *ptmpwork;                              
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSR_MatTriangSlvLU_CATBC_double(
                 const int mb, const int n,
                  
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSR_MatTriangSlvUU_CATBC_double(
                 const int mb, const int n, 
                  
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void BSR_MatTriangSlvLU_CaABC_double(
                 const int mb, const int n, 
                  
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void BSR_MatTriangSlvUU_CaABC_double(
                 const int mb, const int n,
                  
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                            
  double *ptmpwork;                              
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSR_MatTriangSlvLU_CaATBC_double(
                 const int mb, const int n,
                  
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSR_MatTriangSlvUU_CaATBC_double(
                 const int mb, const int n, 
                  
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void BSR_MatTriangSlvLU_CABbC_double(
                 const int mb, const int n, 
                  
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void BSR_MatTriangSlvUU_CABbC_double(
                 const int mb, const int n,
                  
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                            
  double *ptmpwork;                              
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSR_MatTriangSlvLU_CATBbC_double(
                 const int mb, const int n,
                  
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSR_MatTriangSlvUU_CATBbC_double(
                 const int mb, const int n, 
                  
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void BSR_MatTriangSlvLU_CaABbC_double(
                 const int mb, const int n, 
                  
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void BSR_MatTriangSlvUU_CaABbC_double(
                 const int mb, const int n,
                  
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                            
  double *ptmpwork;                              
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSR_MatTriangSlvLU_CaATBbC_double(
                 const int mb, const int n,
                  
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSR_MatTriangSlvUU_CaATBbC_double(
                 const int mb, const int n, 
                  
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void BSR_MatTriangSlvLU_CDAB_double(
                 const int mb, const int n, 
                 const double *dvl, 
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;



  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                   
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          



}


void BSR_MatTriangSlvUU_CDAB_double(
                 const int mb, const int n,
                 const double *dvl, 
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                            
  double *ptmpwork;                              
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;



  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                   
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          




}

void BSR_MatTriangSlvLU_CDATB_double(
                 const int mb, const int n,
                 const double *dvl, 
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;



  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                   
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          




}

void BSR_MatTriangSlvUU_CDATB_double(
                 const int mb, const int n, 
                 const double *dvl, 
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;



  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          




}



void BSR_MatTriangSlvLU_CaDAB_double(
                 const int mb, const int n, 
                 const double *dvl, 
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;



  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                   
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  

  pc=c;                                   
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     

}


void BSR_MatTriangSlvUU_CaDAB_double(
                 const int mb, const int n,
                 const double *dvl, 
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                            
  double *ptmpwork;                              
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;



  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                   
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  
  

  pc=c;                                   
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     


}

void BSR_MatTriangSlvLU_CaDATB_double(
                 const int mb, const int n,
                 const double *dvl, 
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;



  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                   
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  
  

  pc=c;                                   
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     


}

void BSR_MatTriangSlvUU_CaDATB_double(
                 const int mb, const int n, 
                 const double *dvl, 
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;



  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  
  

  pc=c;                                   
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     


}



void BSR_MatTriangSlvLU_CDABC_double(
                 const int mb, const int n, 
                 const double *dvl, 
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                   
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void BSR_MatTriangSlvUU_CDABC_double(
                 const int mb, const int n,
                 const double *dvl, 
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                            
  double *ptmpwork;                              
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                   
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSR_MatTriangSlvLU_CDATBC_double(
                 const int mb, const int n,
                 const double *dvl, 
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                   
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSR_MatTriangSlvUU_CDATBC_double(
                 const int mb, const int n, 
                 const double *dvl, 
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void BSR_MatTriangSlvLU_CaDABC_double(
                 const int mb, const int n, 
                 const double *dvl, 
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                   
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void BSR_MatTriangSlvUU_CaDABC_double(
                 const int mb, const int n,
                 const double *dvl, 
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                            
  double *ptmpwork;                              
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                   
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSR_MatTriangSlvLU_CaDATBC_double(
                 const int mb, const int n,
                 const double *dvl, 
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                   
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSR_MatTriangSlvUU_CaDATBC_double(
                 const int mb, const int n, 
                 const double *dvl, 
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void BSR_MatTriangSlvLU_CDABbC_double(
                 const int mb, const int n, 
                 const double *dvl, 
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                   
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void BSR_MatTriangSlvUU_CDABbC_double(
                 const int mb, const int n,
                 const double *dvl, 
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                            
  double *ptmpwork;                              
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                   
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSR_MatTriangSlvLU_CDATBbC_double(
                 const int mb, const int n,
                 const double *dvl, 
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                   
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSR_MatTriangSlvUU_CDATBbC_double(
                 const int mb, const int n, 
                 const double *dvl, 
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void BSR_MatTriangSlvLU_CaDABbC_double(
                 const int mb, const int n, 
                 const double *dvl, 
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                   
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void BSR_MatTriangSlvUU_CaDABbC_double(
                 const int mb, const int n,
                 const double *dvl, 
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                            
  double *ptmpwork;                              
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                   
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSR_MatTriangSlvLU_CaDATBbC_double(
                 const int mb, const int n,
                 const double *dvl, 
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                   
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSR_MatTriangSlvUU_CaDATBbC_double(
                 const int mb, const int n, 
                 const double *dvl, 
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void BSR_MatTriangSlvLU_CADB_double(
                 const int mb, const int n, 
                  const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                  const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;



  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        



}


void BSR_MatTriangSlvUU_CADB_double(
                 const int mb, const int n,
                  const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                  const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;



  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                    
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        




}

void BSR_MatTriangSlvLU_CATDB_double(
                 const int mb, const int n,
                  const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                  const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;



  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        




}

void BSR_MatTriangSlvUU_CATDB_double(
                 const int mb, const int n, 
                  const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                  const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;



  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        




}



void BSR_MatTriangSlvLU_CaADB_double(
                 const int mb, const int n, 
                  const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                  const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;



  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  

  pc=c;                                   
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     

}


void BSR_MatTriangSlvUU_CaADB_double(
                 const int mb, const int n,
                  const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                  const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;



  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                    
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     


}

void BSR_MatTriangSlvLU_CaATDB_double(
                 const int mb, const int n,
                  const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                  const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;



  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     


}

void BSR_MatTriangSlvUU_CaATDB_double(
                 const int mb, const int n, 
                  const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                  const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;



  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     


}



void BSR_MatTriangSlvLU_CADBC_double(
                 const int mb, const int n, 
                  const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void BSR_MatTriangSlvUU_CADBC_double(
                 const int mb, const int n,
                  const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                            
  double *ptmpwork;                              
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                    
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSR_MatTriangSlvLU_CATDBC_double(
                 const int mb, const int n,
                  const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSR_MatTriangSlvUU_CATDBC_double(
                 const int mb, const int n, 
                  const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void BSR_MatTriangSlvLU_CaADBC_double(
                 const int mb, const int n, 
                  const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void BSR_MatTriangSlvUU_CaADBC_double(
                 const int mb, const int n,
                  const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                            
  double *ptmpwork;                              
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                    
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSR_MatTriangSlvLU_CaATDBC_double(
                 const int mb, const int n,
                  const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSR_MatTriangSlvUU_CaATDBC_double(
                 const int mb, const int n, 
                  const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void BSR_MatTriangSlvLU_CADBbC_double(
                 const int mb, const int n, 
                  const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void BSR_MatTriangSlvUU_CADBbC_double(
                 const int mb, const int n,
                  const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                            
  double *ptmpwork;                              
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                    
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSR_MatTriangSlvLU_CATDBbC_double(
                 const int mb, const int n,
                  const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSR_MatTriangSlvUU_CATDBbC_double(
                 const int mb, const int n, 
                  const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void BSR_MatTriangSlvLU_CaADBbC_double(
                 const int mb, const int n, 
                  const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void BSR_MatTriangSlvUU_CaADBbC_double(
                 const int mb, const int n,
                  const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                            
  double *ptmpwork;                              
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                    
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSR_MatTriangSlvLU_CaATDBbC_double(
                 const int mb, const int n,
                  const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSR_MatTriangSlvUU_CaATDBbC_double(
                 const int mb, const int n, 
                  const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void BSR_MatTriangSlvLU_CDADB_double(
                 const int mb, const int n, 
                 const double *dvl, const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;



  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                   
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          



}


void BSR_MatTriangSlvUU_CDADB_double(
                 const int mb, const int n,
                 const double *dvl, const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                            
  double *ptmpwork;                              
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;



  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                    
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                   
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          




}

void BSR_MatTriangSlvLU_CDATDB_double(
                 const int mb, const int n,
                 const double *dvl, const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;



  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                   
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          




}

void BSR_MatTriangSlvUU_CDATDB_double(
                 const int mb, const int n, 
                 const double *dvl, const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;



  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          




}



void BSR_MatTriangSlvLU_CaDADB_double(
                 const int mb, const int n, 
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;



  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                   
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  

  pc=c;                                   
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     

}


void BSR_MatTriangSlvUU_CaDADB_double(
                 const int mb, const int n,
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                            
  double *ptmpwork;                              
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;



  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                    
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                   
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  
  

  pc=c;                                   
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     


}

void BSR_MatTriangSlvLU_CaDATDB_double(
                 const int mb, const int n,
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;



  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                   
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  
  

  pc=c;                                   
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     


}

void BSR_MatTriangSlvUU_CaDATDB_double(
                 const int mb, const int n, 
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;



  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  
  

  pc=c;                                   
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     


}



void BSR_MatTriangSlvLU_CDADBC_double(
                 const int mb, const int n, 
                 const double *dvl, const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                   
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void BSR_MatTriangSlvUU_CDADBC_double(
                 const int mb, const int n,
                 const double *dvl, const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                            
  double *ptmpwork;                              
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                    
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                   
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSR_MatTriangSlvLU_CDATDBC_double(
                 const int mb, const int n,
                 const double *dvl, const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                   
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSR_MatTriangSlvUU_CDATDBC_double(
                 const int mb, const int n, 
                 const double *dvl, const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void BSR_MatTriangSlvLU_CaDADBC_double(
                 const int mb, const int n, 
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                   
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void BSR_MatTriangSlvUU_CaDADBC_double(
                 const int mb, const int n,
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                            
  double *ptmpwork;                              
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                    
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                   
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSR_MatTriangSlvLU_CaDATDBC_double(
                 const int mb, const int n,
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                   
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSR_MatTriangSlvUU_CaDATDBC_double(
                 const int mb, const int n, 
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb,  
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void BSR_MatTriangSlvLU_CDADBbC_double(
                 const int mb, const int n, 
                 const double *dvl, const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                   
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void BSR_MatTriangSlvUU_CDADBbC_double(
                 const int mb, const int n,
                 const double *dvl, const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                            
  double *ptmpwork;                              
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                    
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                   
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSR_MatTriangSlvLU_CDATDBbC_double(
                 const int mb, const int n,
                 const double *dvl, const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                   
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSR_MatTriangSlvUU_CDATDBbC_double(
                 const int mb, const int n, 
                 const double *dvl, const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void BSR_MatTriangSlvLU_CaDADBbC_double(
                 const int mb, const int n, 
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                   
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void BSR_MatTriangSlvUU_CaDADBbC_double(
                 const int mb, const int n,
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                            
  double *ptmpwork;                              
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                    
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                   
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSR_MatTriangSlvLU_CaDATDBbC_double(
                 const int mb, const int n,
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                   
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSR_MatTriangSlvUU_CaDATDBbC_double(
                 const int mb, const int n, 
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc,
                 double *work, const int ind_base)
{
  const double *pb;
  double *pc=c;
  double *cj;
  const double *pdr;
  const double *pdl;
  double *pwork=work;                             
  double *ptmpwork;                               
  const double *pval;
  const double *ptmp;
  int i,j,jb,je,block;
  int ds,bs,br,cs;
  int index;
  int ind;
  int l;                                               
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

  for (l=0;l!=n;l++)                              
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

    for (l=0;l!=n;l++){                                      
      for (j=0;j!=lb;j++) pc[j] = pb[j];
                    
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += lb+1;                                      
      }                                                   
      pdr = &dvr[i*mm];                                   
      for (j=0;j!=lb;j++) {                               
        for (ii=0;ii!=lb;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
      pc+=ldc;pb+=ldb;                                       
    }                                                        
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
        ptmp = pval;                                         
                                 
        cs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
        for (l=0;l!=n;l++) {            
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
          pc+=ldc;cj+=ldc;pval=ptmp; 
        }                                                    
      }
    }
  }

  
  

  ds = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = i*lb;                  
    pdl = &dvl[ds];                                          
    ds += mm;                
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
    for (l=0;l!=n;l++) {            
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
      pc+=ldc;                      
    }                                                 
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
  for (l=0;l!=n;l++)                      
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}
