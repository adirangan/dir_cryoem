function output = polyder(p_);
% This function finds the derivative dp(x) / dx. ;

n = length(p_);
d_ = reshape(n-1:-1:0,size(p_));
dp_ = p_.*d_;
dp_ = dp_(1:end-1);
output = dp_;

