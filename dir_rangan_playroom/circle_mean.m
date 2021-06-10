function output = circle_mean(g_);
% takes mean on circle of angles in g_. ;
cg_ = cos(g_);
sg_ = sin(g_);
cg_avg = mean(cg_);
sg_avg = mean(sg_);
output = atan2(sg_avg,cg_avg);
