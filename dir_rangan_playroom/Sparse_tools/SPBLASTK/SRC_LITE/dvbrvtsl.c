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


/* Created:  Sat Jul 6 14:38:14 EDT 1996 */

#include "dvbrvtsl.h"



void VBR_VecTriangSlvLU_CAB_double(
                 const int mb,  
                  
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;



  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        



}


void VBR_VecTriangSlvUU_CAB_double(
                 const int mb, 
                  
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;



  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        




}

void VBR_VecTriangSlvLU_CATB_double(
                 const int mb, 
                  
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;



  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
      }
    }
  }

  
  

        




}

void VBR_VecTriangSlvUU_CATB_double(
                 const int mb,  
                  
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;



  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
      }
    }
  }

  
  

        




}



void VBR_VecTriangSlvLU_CaAB_double(
                 const int mb,  
                  
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;



  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  

  pc=c+ind_base;                          
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     

}


void VBR_VecTriangSlvUU_CaAB_double(
                 const int mb, 
                  
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;



  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c+ind_base;                          
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     


}

void VBR_VecTriangSlvLU_CaATB_double(
                 const int mb, 
                  
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;



  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c+ind_base;                          
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     


}

void VBR_VecTriangSlvUU_CaATB_double(
                 const int mb,  
                  
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;



  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c+ind_base;                          
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     


}



void VBR_VecTriangSlvLU_CABC_double(
                 const int mb,  
                  
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void VBR_VecTriangSlvUU_CABC_double(
                 const int mb, 
                  
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void VBR_VecTriangSlvLU_CATBC_double(
                 const int mb, 
                  
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void VBR_VecTriangSlvUU_CATBC_double(
                 const int mb,  
                  
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void VBR_VecTriangSlvLU_CaABC_double(
                 const int mb,  
                  
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void VBR_VecTriangSlvUU_CaABC_double(
                 const int mb, 
                  
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void VBR_VecTriangSlvLU_CaATBC_double(
                 const int mb, 
                  
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void VBR_VecTriangSlvUU_CaATBC_double(
                 const int mb,  
                  
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void VBR_VecTriangSlvLU_CABbC_double(
                 const int mb,  
                  
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void VBR_VecTriangSlvUU_CABbC_double(
                 const int mb, 
                  
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void VBR_VecTriangSlvLU_CATBbC_double(
                 const int mb, 
                  
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void VBR_VecTriangSlvUU_CATBbC_double(
                 const int mb,  
                  
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void VBR_VecTriangSlvLU_CaABbC_double(
                 const int mb,  
                  
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void VBR_VecTriangSlvUU_CaABbC_double(
                 const int mb, 
                  
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void VBR_VecTriangSlvLU_CaATBbC_double(
                 const int mb, 
                  
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void VBR_VecTriangSlvUU_CaATBbC_double(
                 const int mb,  
                  
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void VBR_VecTriangSlvLU_CDAB_double(
                 const int mb,  
                 const double *dvl, 
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;



  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          



}


void VBR_VecTriangSlvUU_CDAB_double(
                 const int mb, 
                 const double *dvl, 
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;



  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          




}

void VBR_VecTriangSlvLU_CDATB_double(
                 const int mb, 
                 const double *dvl, 
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;



  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          




}

void VBR_VecTriangSlvUU_CDATB_double(
                 const int mb,  
                 const double *dvl, 
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;



  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          




}



void VBR_VecTriangSlvLU_CaDAB_double(
                 const int mb,  
                 const double *dvl, 
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;



  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  

  pc=c+ind_base;                          
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     

}


void VBR_VecTriangSlvUU_CaDAB_double(
                 const int mb, 
                 const double *dvl, 
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;



  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c+ind_base;                          
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     


}

void VBR_VecTriangSlvLU_CaDATB_double(
                 const int mb, 
                 const double *dvl, 
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;



  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c+ind_base;                          
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     


}

void VBR_VecTriangSlvUU_CaDATB_double(
                 const int mb,  
                 const double *dvl, 
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;



  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c+ind_base;                          
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     


}



void VBR_VecTriangSlvLU_CDABC_double(
                 const int mb,  
                 const double *dvl, 
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void VBR_VecTriangSlvUU_CDABC_double(
                 const int mb, 
                 const double *dvl, 
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void VBR_VecTriangSlvLU_CDATBC_double(
                 const int mb, 
                 const double *dvl, 
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void VBR_VecTriangSlvUU_CDATBC_double(
                 const int mb,  
                 const double *dvl, 
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void VBR_VecTriangSlvLU_CaDABC_double(
                 const int mb,  
                 const double *dvl, 
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void VBR_VecTriangSlvUU_CaDABC_double(
                 const int mb, 
                 const double *dvl, 
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void VBR_VecTriangSlvLU_CaDATBC_double(
                 const int mb, 
                 const double *dvl, 
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void VBR_VecTriangSlvUU_CaDATBC_double(
                 const int mb,  
                 const double *dvl, 
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void VBR_VecTriangSlvLU_CDABbC_double(
                 const int mb,  
                 const double *dvl, 
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void VBR_VecTriangSlvUU_CDABbC_double(
                 const int mb, 
                 const double *dvl, 
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void VBR_VecTriangSlvLU_CDATBbC_double(
                 const int mb, 
                 const double *dvl, 
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void VBR_VecTriangSlvUU_CDATBbC_double(
                 const int mb,  
                 const double *dvl, 
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void VBR_VecTriangSlvLU_CaDABbC_double(
                 const int mb,  
                 const double *dvl, 
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void VBR_VecTriangSlvUU_CaDABbC_double(
                 const int mb, 
                 const double *dvl, 
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void VBR_VecTriangSlvLU_CaDATBbC_double(
                 const int mb, 
                 const double *dvl, 
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void VBR_VecTriangSlvUU_CaDATBbC_double(
                 const int mb,  
                 const double *dvl, 
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void VBR_VecTriangSlvLU_CADB_double(
                 const int mb,  
                  const double *dvr,
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;



  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        



}


void VBR_VecTriangSlvUU_CADB_double(
                 const int mb, 
                  const double *dvr,
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;



  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        




}

void VBR_VecTriangSlvLU_CATDB_double(
                 const int mb, 
                  const double *dvr,
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;



  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
      }
    }
  }

  
  

        




}

void VBR_VecTriangSlvUU_CATDB_double(
                 const int mb,  
                  const double *dvr,
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;



  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
      }
    }
  }

  
  

        




}



void VBR_VecTriangSlvLU_CaADB_double(
                 const int mb,  
                  const double *dvr,
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;



  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  

  pc=c+ind_base;                          
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     

}


void VBR_VecTriangSlvUU_CaADB_double(
                 const int mb, 
                  const double *dvr,
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;



  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c+ind_base;                          
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     


}

void VBR_VecTriangSlvLU_CaATDB_double(
                 const int mb, 
                  const double *dvr,
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;



  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c+ind_base;                          
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     


}

void VBR_VecTriangSlvUU_CaATDB_double(
                 const int mb,  
                  const double *dvr,
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;



  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c+ind_base;                          
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     


}



void VBR_VecTriangSlvLU_CADBC_double(
                 const int mb,  
                  const double *dvr,
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void VBR_VecTriangSlvUU_CADBC_double(
                 const int mb, 
                  const double *dvr,
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void VBR_VecTriangSlvLU_CATDBC_double(
                 const int mb, 
                  const double *dvr,
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void VBR_VecTriangSlvUU_CATDBC_double(
                 const int mb,  
                  const double *dvr,
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void VBR_VecTriangSlvLU_CaADBC_double(
                 const int mb,  
                  const double *dvr,
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void VBR_VecTriangSlvUU_CaADBC_double(
                 const int mb, 
                  const double *dvr,
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void VBR_VecTriangSlvLU_CaATDBC_double(
                 const int mb, 
                  const double *dvr,
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void VBR_VecTriangSlvUU_CaATDBC_double(
                 const int mb,  
                  const double *dvr,
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void VBR_VecTriangSlvLU_CADBbC_double(
                 const int mb,  
                  const double *dvr,
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void VBR_VecTriangSlvUU_CADBbC_double(
                 const int mb, 
                  const double *dvr,
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void VBR_VecTriangSlvLU_CATDBbC_double(
                 const int mb, 
                  const double *dvr,
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void VBR_VecTriangSlvUU_CATDBbC_double(
                 const int mb,  
                  const double *dvr,
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void VBR_VecTriangSlvLU_CaADBbC_double(
                 const int mb,  
                  const double *dvr,
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void VBR_VecTriangSlvUU_CaADBbC_double(
                 const int mb, 
                  const double *dvr,
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void VBR_VecTriangSlvLU_CaATDBbC_double(
                 const int mb, 
                  const double *dvr,
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void VBR_VecTriangSlvUU_CaATDBbC_double(
                 const int mb,  
                  const double *dvr,
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
      }
    }
  }

  
  

        

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void VBR_VecTriangSlvLU_CDADB_double(
                 const int mb,  
                 const double *dvl, const double *dvr,
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;



  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          



}


void VBR_VecTriangSlvUU_CDADB_double(
                 const int mb, 
                 const double *dvl, const double *dvr,
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;



  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          




}

void VBR_VecTriangSlvLU_CDATDB_double(
                 const int mb, 
                 const double *dvl, const double *dvr,
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;



  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          




}

void VBR_VecTriangSlvUU_CDATDB_double(
                 const int mb,  
                 const double *dvl, const double *dvr,
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;



  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          




}



void VBR_VecTriangSlvLU_CaDADB_double(
                 const int mb,  
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;



  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  

  pc=c+ind_base;                          
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     

}


void VBR_VecTriangSlvUU_CaDADB_double(
                 const int mb, 
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;



  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c+ind_base;                          
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     


}

void VBR_VecTriangSlvLU_CaDATDB_double(
                 const int mb, 
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;



  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c+ind_base;                          
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     


}

void VBR_VecTriangSlvUU_CaDATDB_double(
                 const int mb,  
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;



  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c+ind_base;                          
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          pc++;                           
    }                                     


}



void VBR_VecTriangSlvLU_CDADBC_double(
                 const int mb,  
                 const double *dvl, const double *dvr,
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void VBR_VecTriangSlvUU_CDADBC_double(
                 const int mb, 
                 const double *dvl, const double *dvr,
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void VBR_VecTriangSlvLU_CDATDBC_double(
                 const int mb, 
                 const double *dvl, const double *dvr,
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void VBR_VecTriangSlvUU_CDATDBC_double(
                 const int mb,  
                 const double *dvl, const double *dvr,
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void VBR_VecTriangSlvLU_CaDADBC_double(
                 const int mb,  
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void VBR_VecTriangSlvUU_CaDADBC_double(
                 const int mb, 
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void VBR_VecTriangSlvLU_CaDATDBC_double(
                 const int mb, 
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void VBR_VecTriangSlvUU_CaDATDBC_double(
                 const int mb,  
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void VBR_VecTriangSlvLU_CDADBbC_double(
                 const int mb,  
                 const double *dvl, const double *dvr,
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void VBR_VecTriangSlvUU_CDADBbC_double(
                 const int mb, 
                 const double *dvl, const double *dvr,
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void VBR_VecTriangSlvLU_CDATDBbC_double(
                 const int mb, 
                 const double *dvl, const double *dvr,
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void VBR_VecTriangSlvUU_CDATDBbC_double(
                 const int mb,  
                 const double *dvl, const double *dvr,
                 
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}



void VBR_VecTriangSlvLU_CaDADBbC_double(
                 const int mb,  
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     

}


void VBR_VecTriangSlvUU_CaDADBbC_double(
                 const int mb, 
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pb = &b[bs];                
        pc = &c[bs];                
        cj = &c[cs];                
          for (ii=0;ii!=cr;ii++) {      
            t = cj[ii];
            for (jj=0;jj!=mm;jj++) {        
              pc[jj] -= (*pval++) * t;       
            }
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void VBR_VecTriangSlvLU_CaDATDBbC_double(
                 const int mb, 
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;        
                           
                           
                           
                           
                           
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=mb-1;i!=-1;i--) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    for (j=je-1;j!=jb-1;j--) {  
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;       
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}

void VBR_VecTriangSlvUU_CaDATDBbC_double(
                 const int mb,  
                 const double *dvl, const double *dvr,
                 const double alpha,
                 const double *val, const int *indx,
                 const int *bindx, const int *rpntr,
                 const int *cpntr, const int *bpntrb, const int *bpntre, 
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
  int ds,mm,bs,br,cs,cr;
  int index;
  int ind;
  int ii,jj;
  int m=rpntr[mb]-rpntr[0];
  double t;

  val-=ind_base;
  b-=ind_base;
  c-=ind_base;
  indx-=ind_base;
  bindx-=ind_base;
  cpntr-=ind_base;

  
  
  

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);       

  ptmpwork = pwork;          
                             
                             
                             
                             
                             
  
  
  mm=0; 
  ds=0;
  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pb = &b[bs];                
    pc = &c[bs];                

    

      for (j=0;j!=mm;j++) pc[j] = pb[j];
                    
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        pc[j] *= (*pdr);                                  
        pdr += mm+1;                                      
      }                                                   
      pdr = &dvr[ds];                                     
      for (j=0;j!=mm;j++) {                               
        for (ii=0;ii!=mm;ii++)                            
          if ( ii == j ) pdr++;                           
          else  pc[j] += (*pdr++) * pb[ii];               
      }                                                   
  }

  

  for (i=0;i!=mb;i++) {      
    jb = bpntrb[i];             
    je = bpntre[i];             
    bs = rpntr[i];              
    mm = rpntr[i+1]-rpntr[i];   

    for (j=jb;j!=je;j++) {      
      index = bindx[j];             
      if ( index == i + ind_base ){ 
        continue;                   
      } else {
        pval = &val[indx[j]];       
                                 
        cs = cpntr[index];          
        cr = cpntr[index+1] - cs;   
        pc = &c[cs];                
        cj = &c[bs];                
          for (ii=0;ii!=cr;ii++) {      
            t = 0;
            for (jj=0;jj!=mm;jj++) {        
              t += (*pval++) * cj[jj];      
            }
            pc[ii] -= t;
          }
      }
    }
  }

  
  

  ds = 0;                                                    
  mm = 0;                                                    
  for (i=0;i!=mb;i++) {      
    bs = rpntr[i];              
    ds += mm*mm;                
    mm = rpntr[i+1]-rpntr[i];   
    pdl = &dvl[ds];                                          
    pc = &c[bs];                                             
    pwork=ptmpwork;                                          
      for (j=0;j!=mm;j++) {     
        t = 0;                                               
        for (jj=0;jj!=mm;jj++)                               
          t+= pc[jj]*pdl[j+jj*mm];  
        pwork[j] = t;               
      }                                                      
        
      for (j=0;j!=mm;j++) pc[j] = pwork[j];                   
  }                                                          

  
  
  

  pc=c+ind_base;                          
  pwork=work;                             
    for (i=0;i!=m;i++) {                  
          *pc *= alpha;                   
          *pc += *pwork;                  
          pc++;                           
          pwork++;                        
    }                                     


}
