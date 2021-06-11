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


/* Created:  Sat Jul 6 14:33:27 EDT 1996 */

#include "dbscvtsl.h"



void BSC_VecTriangSlvLU_CAB_double(
                 const int mb,  
                  
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        



}


void BSC_VecTriangSlvUU_CAB_double(
                 const int mb, 
                  
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        




}

void BSC_VecTriangSlvLU_CATB_double(
                 const int mb, 
                  
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
      }
    }
  }

  
  

        




}

void BSC_VecTriangSlvUU_CATB_double(
                 const int mb,  
                  
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
      }
    }
  }

  
  

        




}



void BSC_VecTriangSlvLU_CaAB_double(
                 const int mb,  
                  
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  

  pc=c;                                   
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     

}


void BSC_VecTriangSlvUU_CaAB_double(
                 const int mb, 
                  
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     


}

void BSC_VecTriangSlvLU_CaATB_double(
                 const int mb, 
                  
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     


}

void BSC_VecTriangSlvUU_CaATB_double(
                 const int mb,  
                  
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     


}



void BSC_VecTriangSlvLU_CABC_double(
                 const int mb,  
                  
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void BSC_VecTriangSlvUU_CABC_double(
                 const int mb, 
                  
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSC_VecTriangSlvLU_CATBC_double(
                 const int mb, 
                  
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSC_VecTriangSlvUU_CATBC_double(
                 const int mb,  
                  
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void BSC_VecTriangSlvLU_CaABC_double(
                 const int mb,  
                  
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void BSC_VecTriangSlvUU_CaABC_double(
                 const int mb, 
                  
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSC_VecTriangSlvLU_CaATBC_double(
                 const int mb, 
                  
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSC_VecTriangSlvUU_CaATBC_double(
                 const int mb,  
                  
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void BSC_VecTriangSlvLU_CABbC_double(
                 const int mb,  
                  
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void BSC_VecTriangSlvUU_CABbC_double(
                 const int mb, 
                  
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSC_VecTriangSlvLU_CATBbC_double(
                 const int mb, 
                  
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSC_VecTriangSlvUU_CATBbC_double(
                 const int mb,  
                  
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void BSC_VecTriangSlvLU_CaABbC_double(
                 const int mb,  
                  
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void BSC_VecTriangSlvUU_CaABbC_double(
                 const int mb, 
                  
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSC_VecTriangSlvLU_CaATBbC_double(
                 const int mb, 
                  
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSC_VecTriangSlvUU_CaATBbC_double(
                 const int mb,  
                  
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void BSC_VecTriangSlvLU_CDAB_double(
                 const int mb,  
                 const double *dvl, 
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          



}


void BSC_VecTriangSlvUU_CDAB_double(
                 const int mb, 
                 const double *dvl, 
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          




}

void BSC_VecTriangSlvLU_CDATB_double(
                 const int mb, 
                 const double *dvl, 
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          




}

void BSC_VecTriangSlvUU_CDATB_double(
                 const int mb,  
                 const double *dvl, 
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          




}



void BSC_VecTriangSlvLU_CaDAB_double(
                 const int mb,  
                 const double *dvl, 
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  

  pc=c;                                   
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     

}


void BSC_VecTriangSlvUU_CaDAB_double(
                 const int mb, 
                 const double *dvl, 
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c;                                   
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     


}

void BSC_VecTriangSlvLU_CaDATB_double(
                 const int mb, 
                 const double *dvl, 
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c;                                   
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     


}

void BSC_VecTriangSlvUU_CaDATB_double(
                 const int mb,  
                 const double *dvl, 
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c;                                   
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     


}



void BSC_VecTriangSlvLU_CDABC_double(
                 const int mb,  
                 const double *dvl, 
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void BSC_VecTriangSlvUU_CDABC_double(
                 const int mb, 
                 const double *dvl, 
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSC_VecTriangSlvLU_CDATBC_double(
                 const int mb, 
                 const double *dvl, 
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSC_VecTriangSlvUU_CDATBC_double(
                 const int mb,  
                 const double *dvl, 
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void BSC_VecTriangSlvLU_CaDABC_double(
                 const int mb,  
                 const double *dvl, 
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void BSC_VecTriangSlvUU_CaDABC_double(
                 const int mb, 
                 const double *dvl, 
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSC_VecTriangSlvLU_CaDATBC_double(
                 const int mb, 
                 const double *dvl, 
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSC_VecTriangSlvUU_CaDATBC_double(
                 const int mb,  
                 const double *dvl, 
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void BSC_VecTriangSlvLU_CDABbC_double(
                 const int mb,  
                 const double *dvl, 
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void BSC_VecTriangSlvUU_CDABbC_double(
                 const int mb, 
                 const double *dvl, 
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSC_VecTriangSlvLU_CDATBbC_double(
                 const int mb, 
                 const double *dvl, 
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSC_VecTriangSlvUU_CDATBbC_double(
                 const int mb,  
                 const double *dvl, 
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void BSC_VecTriangSlvLU_CaDABbC_double(
                 const int mb,  
                 const double *dvl, 
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void BSC_VecTriangSlvUU_CaDABbC_double(
                 const int mb, 
                 const double *dvl, 
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSC_VecTriangSlvLU_CaDATBbC_double(
                 const int mb, 
                 const double *dvl, 
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSC_VecTriangSlvUU_CaDATBbC_double(
                 const int mb,  
                 const double *dvl, 
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=lb;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void BSC_VecTriangSlvLU_CADB_double(
                 const int mb,  
                  const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        



}


void BSC_VecTriangSlvUU_CADB_double(
                 const int mb, 
                  const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        




}

void BSC_VecTriangSlvLU_CATDB_double(
                 const int mb, 
                  const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
      }
    }
  }

  
  

        




}

void BSC_VecTriangSlvUU_CATDB_double(
                 const int mb,  
                  const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
      }
    }
  }

  
  

        




}



void BSC_VecTriangSlvLU_CaADB_double(
                 const int mb,  
                  const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  

  pc=c;                                   
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     

}


void BSC_VecTriangSlvUU_CaADB_double(
                 const int mb, 
                  const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     


}

void BSC_VecTriangSlvLU_CaATDB_double(
                 const int mb, 
                  const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     


}

void BSC_VecTriangSlvUU_CaATDB_double(
                 const int mb,  
                  const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     


}



void BSC_VecTriangSlvLU_CADBC_double(
                 const int mb,  
                  const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

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
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void BSC_VecTriangSlvUU_CADBC_double(
                 const int mb, 
                  const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

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
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSC_VecTriangSlvLU_CATDBC_double(
                 const int mb, 
                  const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

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
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSC_VecTriangSlvUU_CATDBC_double(
                 const int mb,  
                  const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

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
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void BSC_VecTriangSlvLU_CaADBC_double(
                 const int mb,  
                  const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

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
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void BSC_VecTriangSlvUU_CaADBC_double(
                 const int mb, 
                  const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

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
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSC_VecTriangSlvLU_CaATDBC_double(
                 const int mb, 
                  const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

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
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSC_VecTriangSlvUU_CaATDBC_double(
                 const int mb,  
                  const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

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
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void BSC_VecTriangSlvLU_CADBbC_double(
                 const int mb,  
                  const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

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
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void BSC_VecTriangSlvUU_CADBbC_double(
                 const int mb, 
                  const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

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
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSC_VecTriangSlvLU_CATDBbC_double(
                 const int mb, 
                  const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

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
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSC_VecTriangSlvUU_CATDBbC_double(
                 const int mb,  
                  const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

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
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void BSC_VecTriangSlvLU_CaADBbC_double(
                 const int mb,  
                  const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

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
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void BSC_VecTriangSlvUU_CaADBbC_double(
                 const int mb, 
                  const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

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
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSC_VecTriangSlvLU_CaATDBbC_double(
                 const int mb, 
                  const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

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
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSC_VecTriangSlvUU_CaATDBbC_double(
                 const int mb,  
                  const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

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
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void BSC_VecTriangSlvLU_CDADB_double(
                 const int mb,  
                 const double *dvl, const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          



}


void BSC_VecTriangSlvUU_CDADB_double(
                 const int mb, 
                 const double *dvl, const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          




}

void BSC_VecTriangSlvLU_CDATDB_double(
                 const int mb, 
                 const double *dvl, const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          




}

void BSC_VecTriangSlvUU_CDATDB_double(
                 const int mb,  
                 const double *dvl, const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          




}



void BSC_VecTriangSlvLU_CaDADB_double(
                 const int mb,  
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  

  pc=c;                                   
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     

}


void BSC_VecTriangSlvUU_CaDADB_double(
                 const int mb, 
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c;                                   
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     


}

void BSC_VecTriangSlvLU_CaDATDB_double(
                 const int mb, 
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c;                                   
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     


}

void BSC_VecTriangSlvUU_CaDATDB_double(
                 const int mb,  
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c;                                   
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     


}



void BSC_VecTriangSlvLU_CDADBC_double(
                 const int mb,  
                 const double *dvl, const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

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
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void BSC_VecTriangSlvUU_CDADBC_double(
                 const int mb, 
                 const double *dvl, const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

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
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSC_VecTriangSlvLU_CDATDBC_double(
                 const int mb, 
                 const double *dvl, const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

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
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSC_VecTriangSlvUU_CDATDBC_double(
                 const int mb,  
                 const double *dvl, const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

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
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void BSC_VecTriangSlvLU_CaDADBC_double(
                 const int mb,  
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

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
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void BSC_VecTriangSlvUU_CaDADBC_double(
                 const int mb, 
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

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
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSC_VecTriangSlvLU_CaDATDBC_double(
                 const int mb, 
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

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
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSC_VecTriangSlvUU_CaDATDBC_double(
                 const int mb,  
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,   
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

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
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void BSC_VecTriangSlvLU_CDADBbC_double(
                 const int mb,  
                 const double *dvl, const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

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
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void BSC_VecTriangSlvUU_CDADBbC_double(
                 const int mb, 
                 const double *dvl, const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

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
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSC_VecTriangSlvLU_CDATDBbC_double(
                 const int mb, 
                 const double *dvl, const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

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
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSC_VecTriangSlvUU_CDATDBbC_double(
                 const int mb,  
                 const double *dvl, const double *dvr,
                 
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

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
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void BSC_VecTriangSlvLU_CaDADBbC_double(
                 const int mb,  
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

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
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void BSC_VecTriangSlvUU_CaDADBbC_double(
                 const int mb, 
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

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
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = lb*i;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=lb;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=lb;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSC_VecTriangSlvLU_CaDATDBbC_double(
                 const int mb, 
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

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
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void BSC_VecTriangSlvUU_CaDATDBbC_double(
                 const int mb,  
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *bindx, 
                 const int *bpntrb, const int *bpntre, const int lb,
                 const double *b,  const double beta, 
                 double *c, 
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
  int ii,jj;
  int m=mb*lb;
  int mm=lb*lb;
  double t;

  bindx-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = i*lb;                  
    pb = &b[bs];                
    pc = &c[bs];                

    

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
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    cs = i*lb;                  

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[(j-1)*mm];      
                                 
        bs = (index-1)*lb;          
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=lb;ii++) {      
            t = 0;
            for (jj=0;jj!=lb;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
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
      for (j=0;j!=lb;j++) {     
        t = 0;                                               
        for (jj=0;jj!=lb;jj++)                               
          t+= pc[jj]*pdl[j+jj*lb];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=lb;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c;                                   
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}
