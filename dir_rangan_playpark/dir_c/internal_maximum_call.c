#ifndef _MONOLITH
#include "playpark_header.h"
#endif /* _MONOLITH */

int is_internal_maximum(int n_Z,double *Z_,int index)
{
  int verbose=0;
  int flag_decrease_and = 0;
  int flag_decrease_neg = 0;
  int flag_decrease_pos = 0;
  int tab=0;
  int continue_flag=0;
  /* %%%%%%%%; */
  flag_decrease_neg = 0;
  tab = index-1; continue_flag=(tab>=0 ? 1 : 0);
  while (continue_flag){
    if ( (continue_flag==1) && (Z_(1+tab)< Z_(1+index)) ){ flag_decrease_neg = 1; continue_flag=0;}
    if ( (continue_flag==1) && (Z_(1+tab)> Z_(1+index)) ){ flag_decrease_neg = 0; continue_flag=0;}
    if ( (continue_flag==1) && (Z_(1+tab)==Z_(1+index)) ){ tab=tab-1; continue_flag=1;}
    if (tab<=0){ continue_flag=0;}
    /* while (continue_flag){ } */}
  if (verbose){ printf(" %% 1+index %d tab %d neg %d\n",1+index,tab,flag_decrease_neg);}
  /* %%%%%%%%; */
  flag_decrease_pos = 0;
  tab = index-1+1; continue_flag=(tab<n_Z-1 ? 1 : 0);
  while (continue_flag){
    if ( (continue_flag==1) && (Z_(1+tab)< Z_(1+index)) ){ flag_decrease_pos = 1; continue_flag=0;}
    if ( (continue_flag==1) && (Z_(1+tab)> Z_(1+index)) ){ flag_decrease_pos = 0; continue_flag=0;}
    if ( (continue_flag==1) && (Z_(1+tab)==Z_(1+index)) ){ tab=tab+1; continue_flag=1;}
    if (tab>=n_Z){ continue_flag=0;}
    /* while (continue_flag){ } */}
  if (verbose){ printf(" %% 1+index %d tab %d pos %d\n",1+index,tab,flag_decrease_pos);}
  /* %%%%%%%%; */
  flag_decrease_and = (flag_decrease_neg && flag_decrease_pos);
  return flag_decease_and;
}

void find_local_maxima(int n_Z,double *Z_,int *n_index_p,int **index_p_)
{
  double Z_gmax=0;
  int nZ=0;
  int n_index=0,nindex=0;
  int *index_;
  array_stats(Z_,"double",n_Z,&Z_gmax,NULL,NULL,NULL);
  n_index=0;
  for (nZ=0;nZ<n_Z;nZ++){ if ( (Z_[nZ]==Z_gmax) || is_internal_maximum(Z_,nZ) ){ n_index++;}}
  index_=NULL;
  if (index_p_!=NULL){
    if ( (*index_p_)==NULL ){ (*index_p_) = (int *) malloci(n_index*sizeof(int)); }
    index_ = *index_p_;
    /* if (index_p_!=NULL){ } */}
  if (index_!=NULL){
    nindex=0;
    for (nZ=0;nZ<n_Z;nZ++){ if ( (Z_[nZ]==Z_gmax) || is_internal_maximum(Z_,nZ) ){ index_[nindex]=nZ; nindex++;}}
    /* if (index_!=NULL){ } */}
  if (n_index_p!=NULL){ (*n_index_p) = n_index;}
}

void find_internal_maximum(int n_Z,double *Z_,double Z_min,double *zone_max_p,int *zone_max_index_p)
{
  int verbose=1;
  int tmp_index=0;
  double Z_gmax=0;
  int nzone_lmax_index=0,n_zone_lmax_index=0,zone_lmax_index=0,*zone_lmax_index_=NULL;
  int zone_max_index=0,*zone_max_index_=NULL;
  int n_index_sub=0,*index_sub_=NULL;
  int nzone=0,n_zone=0,*nzone_start_=NULL,*nzone_final_=NULL;
  int nzone_index=0,n_zone_index=0,*n_zone_index_=NULL,*zone_index__=NULL,*zone_index_=NULL;
  double zone_max=0,*zone_max_=NULL;
  int *zone_max_index_=NULL;
  int sum_zone_valid_=0,n_zone_valid_index=0,*zone_valid_=NULL,*zone_valid_index_=NULL;
  int n_Z_sub=0;
  double *Z_sub_=NULL;
  double tmp_min=0;
  int z_index=0;
  double *zone_max_sub_=NULL;
  double tmp_Z_max=0,tmp_Z_min=0,tmp_Z_mid=0;
  if (verbose){ printf(" %% [entering find_internal_maximum], n_Z %d Z_min %f\n",n_Z,Z_min);}
  if (verbose>1){ array_printf(Z_,"double",1,n_Z," % Z_: ");}
  array_maximum_minimum(Z_,"double",n_Z,&Z_gmax,&tmp_index,NULL,NULL); zone_max = Z_gmax; zone_max_index = tmp_index;
  if (Z_gmax< Z_min){ zone_max = Z_gmax; zone_max_index = -1;}
  if (Z_gmax>=Z_min){
    index_sub_ = (int *) malloci((1+n_Z)*sizeof(int));
    n_index_sub=0; for (nZ=0;nZ<n_Z;nZ++){ if (Z_[nZ]>=Z_min){ index_sub_[n_index_sub] = nZ; n_index_sub++;}}
    n_zone = 1;
    for (nindex_sub=0;nindex_sub<n_index_sub-1;nindex_sub++){
      if (index_sub_[nindex_sub+1]-index_sub_[nindex_sub+0]>1){ n_zone+=1;}
      /* for (nindex_sub=0;nindex_sub<n_index_sub-1;nindex_sub++){ } */}
    if (verbose){ printf(" %% found n_zone %d\n",n_zone);}
    nzone_start_ = (int *) malloci((1+n_Z)*sizeof(int));
    nzone_final_ = (int *) malloci((1+n_Z)*sizeof(int));
    nzone_start_[0] = index_sub_[0];
    nzone_final_[n_zone-1] = index_sub_[n_index_sub-1];
    nzone=0;
    for (nindex_sub=0;nindex_sub<n_index_sub-1;nindex_sub++){
      if (index_sub_[nindex_sub+1]-index_sub_[nindex_sub+0]>1){
	nzone_final_[nzone+0] = index_sub_[nindex_sub+0];
	nzone_start_[nzone+1] = index_sub_[nindex_sub+1];
	nzone+=1;
	/* if (index_sub_[nindex_sub+1]-index_sub_[nindex_sub+0]>1){ } */}
      /* for (nindex_sub=0;nindex_sub<n_index_sub-1;nindex_sub++){ } */}
    if (verbose){ for (nzone=0;nzone<n_zone;nzone++){ printf(" %% zone %d: [%d,%d]\n",nzone,nzone_start_[nzone],nzone_final_[nzone]);}}
    zone_index__ = (int *) malloci((unsigned long long int)n_Z*(unsigned long long int)*n_zone*sizeof(int));
    n_zone_index_ = (int *) malloci(n_zone*sizeof(int));
    for (nzone=0;nzone<n_zone;nzone++){
      zone_index_ = zone_index__ + nzone*n_Z;
      n_zone_index_[nzone] = nzone_final_[nzone] - nzone_start_[nzone] + 1;
      n_zone_index = n_zone_index_[nzone];
      for (nzone_index=0;nzone_index<n_zone_index;nzone_index++){
	zone_index_[nzone_index] = nzone_start_[nzone] + nzone_index;
	/* for (nzone_index=0;nzone_index<n_zone_index;nzone_index++){ } */}
      /* for (nzone=0;nzone<n_zone;nzone++){ } */}
    zone_max_ = (double *) malloci(n_zone*sizeof(double));
    zone_max_index_ = (double *) malloci(n_zone*sizeof(int));
    zone_valid_ = (double *) malloci(n_zone*sizeof(int));
    Z_sub_ = (double *) malloci(n_Z*sizeof(double));
    for (nzone=0;nzone<n_zone;nzone++){
      zone_index_ = zone_index__ + nzone*n_Z;
      n_zone_index = n_zone_index_[nzone];
      array_extract_d_from_d(n_Z,1,Z_,n_zone_index,zone_index_,0,NULL,&Z_sub_,NULL);
      array_maximum_minimum(Z_sub_,"double",n_zone_index,&(zone_max_[nzone]),&tmp_index,NULL,NULL);
      zone_max_index_[nzone] = zone_index_[tmp_index];
      if ( (nzone> 0) && (nzone< n_zone-1) ){ zone_valid_[nzone] = 1;}
      if ( (nzone==0) && (zone_max_index_[nzone]> 0) && (zone_max_index_[nzone]< n_Z-1) ){ zone_valid_[nzone] = 1;}
      if ( (nzone==n_zone-1) && (zone_max_index_[nzone]> 0) && (zone_max_index_[nzone]< n_Z-1) ){ zone_valid_[nzone] = 1;}
      /* for (nzone=0;nzone<n_zone;nzone++){ } */}
    if (verbose){ for (nzone=0;nzone<n_zone;nzone++){ printf(" %% zone %d: [%d,%d], max %0.2f, index %d, valid %d\n",nzone,nzone_start_[nzone],nzone_final_[nzone],zone_max_[nzone],zone_max_index_[nzone],zoone_valid_[nzone]);}}
    zone_max_sub_ = (double *) malloci(n_zone*sizeof(double));
    sum_zone_valid=0; for (nzone=0;nzone<n_zone;nzone++){ sum_zone_valid+=zone_valid_[nzone];}
    zone_valid_index_ = (int *) malloci(n_zone*sizeof(int));
    n_zone_valid_index=0; for (nzone=0;nzone<n_zone;nzone++){ if (zone_valid_[nzone]){ zone_valid_index_[n_zone_valid_index]=nzone; n_zone_valid_index++;}}
    zone_lmax_index_ = (int *) malloci(n_Z*sizeof(int));
    if (sum_zone_valid> 0){
      if (verbose){ printf(" %% found %d valid zones\n",sum_zone_valid);}
      array_extract_d_from_d(n_zone,1,zone_max_,n_zone_valid_index,zone_valid_index_,0,NULL,&zone_max_sub_,NULL);      
      array_maximum_minimum(zone_max_sub_,"double",n_zone_valid_index,&zone_max,&z_index,NULL,NULL);
      z_index = zone_valid_index_[z_index];
      if (verbose){ printf(" %% found max zone_max %f in zone z_index %d\n",zone_max,z_index);}
      if (abs(zone_max-zone_max_[z_index])>1e-12){ printf(" %% Warning, zone_max %f != zone_max_[%d] %f in find_internal_maximum\n",zone_max,z_index,zone_max_[z_index]);}
      zone_max_index = zone_max_index_[z_index];
      /* if (sum_zone_valid> 0){ } */}
    if (sum_zone_valid==0){
      if (verbose){ printf(" %% found no valid zones\n");}
      array_maximum_minimum(zone_max_,"double",n_zone,&zone_max,&z_index,NULL,NULL);
      if (abs(zone_max-zone_max_[z_index])>1e-12){ printf(" %% Warning, zone_max %f != zone_max_[%d] %f in find_internal_maximum\n",zone_max,z_index,zone_max_[z_index]);}
      zone_max_index = zone_max_index_[z_index];
      if (verbose){ printf(" %% found max zone_max %f in zone z_index %d\n",zone_max,z_index);}
      if (0){ /* do nothing */}
      else if ( ( (zone_max_index==0) || (zone_max_index==n_Z-1) ) && (n_zone==1) && (n_zone_index_[0]==n_Z) ){
	array_maximum_minimum(Z_,"double",n_Z,NULL,NULL,&tmp_min,NULL);
	tmp_min += maximum(1e-12,abs(1e-12*tmp_min));
	if (verbose){ printf(" %% only one big zone, rerunning with new minimum Z_min %f --> %f\n",Z_min,tmp_min);}
	find_internal_maximum(n_Z,Z_,tmp_min,&zone_max,&zone_max_index);
	/* else */}
      else if ( (zone_max_index==0) && (n_zone_index_[z_index]< n_Z) ){
	if (verbose){ printf(" %% zone on left side\n");}
	zone_index_ = zone_index__ + z_index*n_Z;
	n_zone_index = n_zone_index_[z_index];
	array_extract_d_from_d(n_Z,1,Z_,n_zone_index,zone_index_,0,NULL,&Z_sub_,NULL);
	find_local_maxima(n_zone_index,Z_sub_,&n_zone_lmax_index,&zone_lmax_index_);
	for (nzone_lmax_index=0;nzone_lmax_index<n_zone_lmax_index;nzone_lmax_index++){ zone_lmax_index_[nzone_lmax_index] = zone_index_[zone_lmax_index];}
	if (n_zone_lmax_index> 1){
	  if (verbose){ printf(" %% found local maximum\n");}
	  zone_max_index = zone_lmax_index_[1]; zone_max = Z_[zone_max_index];
	  /* if (n_zone_lmax_index> 1){ } */}
	if (n_zone_lmax_index<=1){
	  if (verbose){ printf(" %% could not find local maximum, choosing midpoint\n");}
	  array_maximum_minimum(Z_sub_,"double",n_zone_index,&tmp_Z_max,NULL,&tmp_Z_min,NULL);
	  tmp_Z_mid = 0.5*(tmp_Z_max + tmp_Z_min);
	  for (nzone_index=0;nzone_index<n_zone_index;nzone_index++){ Z_sub_[nzone_index] = fabs(Z_sub_[nzone_index]-tmp_Z_mid);}
	  array_maximum_minimum(Z_sub_,"double",n_zone_index,NULL,NULL,NULL,&zone_max_index);
	  zone_max_index = zone_index_[zone_max_index];
	  zone_max = Z_[zone_max_index];
	  /* if (n_zone_lmax_index<=1){ } */}
	/* else */}
      else if ( (zone_max_index==n_Z-1) && (n_zone_index_[z_index]< n_Z) ){
	if (verbose){ printf(" %% zone on right side\n");}
	zone_index_ = zone_index__ + z_index*n_Z;
	n_zone_index = n_zone_index_[z_index];
	array_extract_d_from_d(n_Z,1,Z_,n_zone_index,zone_index_,0,NULL,&Z_sub_,NULL);
	find_local_maxima(n_zone_index,Z_sub_,&n_zone_lmax_index,&zone_lmax_index_);
	for (nzone_lmax_index=0;nzone_lmax_index<n_zone_lmax_index;nzone_lmax_index++){ zone_lmax_index_[nzone_lmax_index] = zone_index_[zone_lmax_index];}
	if (n_zone_lmax_index> 1){
	  if (verbose){ printf(" %% found local maximum\n");}
	  zone_max_index = zone_lmax_index_[1]; zone_max = Z_[zone_max_index];
	  /* if (n_zone_lmax_index> 1){ } */}
	if (n_zone_lmax_index<=1){
	  if (verbose){ printf(" %% could not find local maximum, choosing midpoint\n");}
	  array_maximum_minimum(Z_sub_,"double",n_zone_index,&tmp_Z_max,NULL,&tmp_Z_min,NULL);
	  tmp_Z_mid = 0.5*(tmp_Z_max + tmp_Z_min);
	  for (nzone_index=0;nzone_index<n_zone_index;nzone_index++){ Z_sub_[nzone_index] = fabs(Z_sub_[nzone_index]-tmp_Z_mid);}
	  array_maximum_minimum(Z_sub_,"double",n_zone_index,NULL,NULL,NULL,&zone_max_index);
	  zone_max_index = zone_index_[zone_max_index];
	  zone_max = Z_[zone_max_index];
	  /* if (n_zone_lmax_index<=1){ } */}
	/* else */}
      /* if (sum_zone_valid==0){ } */}
    if (verbose){ printf(" %% zone_max %f <-- %f zone_max_index %d\n",zone_max,Z_[zone_max_index],zone_max_index);}
    free1(zone_lmax_index_);zone_lmax_index_=NULL;
    free1(zone_valid_index_);zone_valid_index_=NULL;
    free1(zone_max_sub_);zone_max_sub_=NULL;
    free1(Z_sub_);Z_sub_=NULL;
    free1(zone_valid_);zone_valid_=NULL;
    free1(zone_max_index_);zone_max_index_=NULL;
    free1(zone_max_);zone_max_=NULL;
    free1(n_zone_index_);n_zone_index_=NULL;
    free1(zone_index__);zone_index__=NULL;
    free1(nzone_start_);nzone_start_=NULL;
    free1(nzone_final_);nzone_final_=NULL;
    free1(index_sub_);index_sub_=NULL;
    /* if (Z_gmax>=Z_min){ } */}
  if (zone_max_p!=NULL){ *zone_max_p = zone_max;}
  if (zone_max_index_p!=NULL){ *zone_max_index_p = zone_max_index;}
  if (verbose){ printf(" %% [finished find_internal_maximum]\n");}
}
