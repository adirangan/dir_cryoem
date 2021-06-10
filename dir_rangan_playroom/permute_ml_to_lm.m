function [m1_,l1_,m2_,l2_,pp_,qq_] = permute_ml_to_lm(l_max);
% This function first generates the standard index ordering m1_ and l1_ for spherical harmonic coefficients, ;
% where m varies more quickly than l ;
% (with the former ranging from -l to +l) : ;
% (m1_,l1_) = 
%  { 
%    (0,0) , 
%    (-1,1),(0,1),(+1,1) , 
%    (-2,2),(-1,2),(0,2),(+1,2),(+2,2) , 
%    ... , 
%    (-l_max,l_max),(-l_max+1,l_max),...,(+l_max-1,l_max),(+l_max,l_max) 
%  }. ;
% This function then generates the indices m2_, and l2_ ;
% so that l varies more quickly than m ;
% (with the former ranging from |m| to l_max, and the latter increasing in absolute value): ;
% (m2_,l2_) = 
% { 
%    (0,0),(0,1),...,(0,l_max) , 
%    (-1,1),(-1,2),...,(-1,l_max), 
%    (+1,1),(+1,2),...,(+1,l_max) , 
%    (-2,2),...,(-2,l_max) , 
%    (+2,2),...,(+2,l_max) , 
%    ... , 
%    (-l_max,l_max) , 
%    (+l_max,l_max) 
%  }. ;
% The permutation pp_ is structured so that: ;
% m1_(pp_) = m2_ ;
% and
% l1_(pp_) = l2_ ;
% The permutation qq_ is structured so that: ;
% m2_(qq_) = m1_ ;
% and
% l2_(qq_) = l1_ ;

l1_ = []; m1_ = []; 
for nl=0:l_max; 
l1_ = [l1_ , nl*ones(1,2*nl+1) ]; 
m1_ = [m1_ , [-nl:+nl] ]; 
end;%for nl=0:l_max;

pp_ = [];
m2_ = [];
l2_ = [];
for nm=0:l_max;
l0_ = nm:l_max;
ll_ = l0_.*(l0_+1);
if (nm==0); pp_ = [pp_,ll_]; m2_ = [m2_,zeros(1,1+l_max)]; l2_ = [l2_,l0_]; end;
if (nm>0); 
%pp_tmp = [ll_-nm ; ll_+nm]; pp_tmp = reshape(pp_tmp,1,2*(l_max-nm+1));
%m2_tmp = [-nm*ones(1,l_max-nm+1) ; +nm*ones(1,l_max-nm+1)]; m2_tmp = reshape(m2_tmp,1,2*(l_max-nm+1));
%l2_tmp = [l0_ ; l0_]; l2_tmp = reshape(l2_tmp,1,2*(l_max-nm+1));
pp_tmp = [ll_-nm ; ll_+nm]; pp_tmp = reshape(transpose(pp_tmp),1,2*(l_max-nm+1));
m2_tmp = [-nm*ones(1,l_max-nm+1) ; +nm*ones(1,l_max-nm+1)]; m2_tmp = reshape(transpose(m2_tmp),1,2*(l_max-nm+1));
l2_tmp = [l0_ ; l0_]; l2_tmp = reshape(transpose(l2_tmp),1,2*(l_max-nm+1));
pp_ = [pp_,pp_tmp]; m2_ = [m2_,m2_tmp]; l2_ = [l2_,l2_tmp];
end;%if (nm>0); 
end;%for nm=0:l_max;

pp_ = 1+pp_;
[~,qq_] = sort(pp_);
