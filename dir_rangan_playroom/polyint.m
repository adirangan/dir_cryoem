function output = polyint(p_,a_);
% If a_ is empty: ;
% This function finds the indefinite integral of a polynomial p(x) with coefficients p_(:). ;
% If a_ is nonempty: ;
% This function finds the definite integral \int_{a_(1)}^{a_(2)} p(x)dx. ;

P_ = [];
p_=transpose(p_(:));
for index=1:length(p_);
P_(index) = p_(index)/(length(p_)-index+1);
end; %for index=1:length(p_);
P_ = [P_ 0];

if nargin<2;
output = P_; 
 else;
d_ = polyval(P_,a_);
output = d_(2) - d_(1);
end;%if nargin<2;

