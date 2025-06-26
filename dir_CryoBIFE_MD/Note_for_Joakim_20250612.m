%%%%%%%%;
% visualize \One^{\intercal}. ;
%%%%%%%%;

flag_verbose = 1;
flag_disp = 1; nf=0;

figure(1+nf);nf=nf+1;clf;figbig;
linewidth_big = 4;
linewidth_sml = 2;
fontsize_use = 12;
subplot(1,1,1);
hold on;

parameter_2 = struct('type','parameter_2');
parameter_2.x_0_lim_ = 2*[-1,+1];
parameter_2.x_1_lim_ = 2*[-1,+1];
parameter_2.x_2_lim_ = 2*[-1,+1];
parameter_2.n_line = 1+8;
parameter_2.linecolor_use = 0.65*[1,1,1];
parameter_2.linewidth_use = 0.5;
parameter_2.patchcolor_use = 0.95*[1,1,1];
plot_box_grid_0(parameter_2);
parameter_1 = struct('type','parameter_1');
parameter_1.x_0_lim_ = 1*[-1,+1];
parameter_1.x_1_lim_ = 1*[-1,+1];
parameter_1.x_1_lim_ = 1*[-1,+1];
parameter_1.n_line = 1+4;
parameter_1.linecolor_use = 0.35*[1,1,1];
parameter_1.linewidth_use = 0.5;
parameter_1.patchcolor_use = 0.85*[1,1,1];
parameter_1.flag_solid = 0;
plot_box_grid_0(parameter_1);
axis equal;
view([-65,20]); 

n_3 = 3;
n_v = 2.^n_3;
vertex_3v__ = zeros(n_3,n_v);
for nf=0:n_v-1;
tmp_3_ = 2*(dec2bin(nf,n_3)-'0')-1;
vertex_3v__(:,1+nf) = reshape(tmp_3_,[n_3,1]);
end;%for nf=0:n_v-1;
n_f = 2*n_3;
face_normal_3f__ = zeros(n_3,n_f);
face_offset_f_ = zeros(n_f,1);
nf=0;
for n3=0:n_3-1;
face_offset_f_(1+nf) = 1.0;
face_normal_3f__(:,1+nf) = circshift([+1;0;0],n3); nf=nf+1;
face_offset_f_(1+nf) = 1.0;
face_normal_3f__(:,1+nf) = circshift([-1;0;0],n3); nf=nf+1;
end;%for n3=0:n_3-1;

One_3_ = [1;1;1];
lb_3_ = -One_3_;
ub_3_ = +One_3_;
One_offset_base = 0.0;
One_offset_ = periodize(One_offset_base + sqrt(n_3)*[0:n_3-1],-1.5*sqrt(n_3),+1.5*sqrt(n_3));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for One_offset = One_offset_;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

li_33f___ = zeros(n_3,n_3,n_f);
ui_33f___ = zeros(n_3,n_3,n_f);
flag_3f__ = zeros(n_3,n_f);
for nf=0:n_f-1;
tmp_Aeq_23__ = [ ...
 transpose(One_3_) ...
;transpose(face_normal_3f__(:,1+nf)) ...
];
tmp_beq_2_ = [ ...
 One_offset ...
;face_offset_f_(1+nf) ...
];
tmp_opt = optimoptions('linprog','Algorithm','interior-point','Display','none');
for n3=0:n_3-1;
tmp_cl_3_ = +circshift([1;0;0],n3);
tmp_xl_3_ = linprog(tmp_cl_3_,[],[],tmp_Aeq_23__,tmp_beq_2_,lb_3_,ub_3_,tmp_opt);
tmp_cu_3_ = -circshift([1;0;0],n3);
tmp_xu_3_ = linprog(tmp_cu_3_,[],[],tmp_Aeq_23__,tmp_beq_2_,lb_3_,ub_3_,tmp_opt);
if numel(tmp_xl_3_)==n_3 & numel(tmp_xu_3_)==n_3;
flag_3f__(1+n3,1+nf) = 1;
li_33f___(:,1+n3,1+nf) = tmp_xl_3_;
ui_33f___(:,1+n3,1+nf) = tmp_xu_3_;
end;%if numel(tmp_xl_3_)==n_3 & numel(tmp_xu_3_)==n_3;
end;%for n3=0:n_3-1;
end;%for nf=0:n_f-1;

%patch([-1;+1;-1],[+1;-1;-1],[-1;-1;+1],'r','EdgeColor','none');

for nf=0:n_f-1;
for n3=0:n_3-1;
if flag_3f__(1+n3,1+nf);
l = line( ...
 [li_33f___(1+0,1+n3,1+nf);ui_33f___(1+0,1+n3,1+nf)] ...
,[li_33f___(1+1,1+n3,1+nf);ui_33f___(1+1,1+n3,1+nf)] ...
,[li_33f___(1+2,1+n3,1+nf);ui_33f___(1+2,1+n3,1+nf)] ...
);
set(l,'LineWidth',linewidth_big,'Color','c');
end;%if flag_3f__(1+n3,1+nf);
end;%for n3=0:n_3-1;
end;%for nf=0:n_f-1;

n_p = 2*sum(flag_3f__,'all');
p_list_3p__ = zeros(n_3,n_p);
np=0;
for nf=0:n_f-1;
for n3=0:n_3-1;
if flag_3f__(1+n3,1+nf);
tmp_3_ = li_33f___(:,1+n3,1+nf);
flag_offset = abs(dot(tmp_3_,One_3_)-One_offset)<1e-6;
if flag_offset; p_list_3p__(:,1+np) = tmp_3_; np=np+1; end;
tmp_3_ = ui_33f___(:,1+n3,1+nf);
flag_offset = abs(dot(tmp_3_,One_3_)-One_offset)<1e-6;
if flag_offset; p_list_3p__(:,1+np) = tmp_3_; np=np+1; end;
end;%if flag_3f__(1+n3,1+nf);
end;%for n3=0:n_3-1;
end;%for nf=0:n_f-1;
assert(np==n_p);

p_list_p3__ = transpose(p_list_3p__);
p_offset_p_ = p_list_p3__*One_3_;
p_norm_p_ = sqrt(sum(abs(p_list_p3__).^2,2));
[~,ij_base] = max(p_norm_p_); index_base = ij_base-1;
p_norm_p3__ = bsxfun(@rdivide,p_list_p3__,max(1e-12,p_norm_p_));
p_norm_3p__ = transpose(p_norm_p3__);
p_norm_base_3_ = p_norm_3p__(:,1+index_base);
p_norm_next_3_ = cross(p_norm_base_3_,One_3_/fnorm(One_3_));
p_cos_p_ = p_norm_p3__*p_norm_base_3_;
p_sin_p_ = p_norm_p3__*p_norm_next_3_;
p_psi_p_ = atan2(p_sin_p_,p_cos_p_);
[~,ij_psi_srt_] = sort(p_psi_p_,'ascend'); index_psi_srt_ = ij_psi_srt_ - 1;
p_sort_p3__ = p_list_p3__(1+index_psi_srt_,:);
p_wrap_p3__ = cat(1,p_sort_p3__,p_sort_p3__(1+0,:));
p = patch(p_wrap_p3__(:,1+0),p_wrap_p3__(:,1+1),p_wrap_p3__(:,1+2),1);
set(p,'EdgeColor','none','FaceColor','g','FaceAlpha',0.5);
l = line(p_wrap_p3__(:,1+0),p_wrap_p3__(:,1+1),p_wrap_p3__(:,1+2),'Color','r','LineWidth',linewidth_sml);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for One_offset = One_offset_;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

xlabel('\psi_{0}');
ylabel('\psi_{1}');
zlabel('\psi_{2}');
set(gca,'FontSize',fontsize_use);
hold off;
