%两组对应点云向量间夹角  (结果是°)
%a(nx3)   b(nx3)
function d_degree=Pnts_normal_angle(a,b)
norm_a=sqrt(sum(a.^2,2));
norm_b=sqrt(sum(b.^2,2));
z=acos(dot(a',b')./(norm_a.*norm_b)');
d_degree=real(min(z/pi*180,180-z/pi*180));