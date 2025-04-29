
%在平面上寻找0-360度范围的许多向量vectors
%输入：平面点云input_pnts（nx3）   法向量长度radius（取1即可），角度间隔resolution(单位是度)
%输出：vectors_onplane（nx3）是单位向量
function [vectors_onplane] = onplane_vector(input_pnts,radius,resolution)
  [parameter] = TLS_Plane(input_pnts);
  center=[0,0,-parameter(4)/parameter(3)];
 r=radius;
 c=center;
 n=parameter(1:3,:)';
 theta_interval=pi*resolution/180;
 theta=(0:theta_interval:2*pi)';
 
a=cross(n,[1 0 0]); %n与i叉乘，求取a向量
if ~any(a) %如果a为零向量，将n与j叉乘
a=cross(n,[0 1 0]);
end
b=cross(n,a); %求取b向量
a=a/norm(a); %单位化a向量
b=b/norm(b);
c1=c(1)*ones(size(theta,1),1);
c2=c(2)*ones(size(theta,1),1);
c3=c(3)*ones(size(theta,1),1);
x=c1+r*a(1)*cos(theta)+r*b(1)*sin(theta);%圆上各点的x坐标
y=c2+r*a(2)*cos(theta)+r*b(2)*sin(theta);%圆上各点的y坐标
z=c3+r*a(3)*cos(theta)+r*b(3)*sin(theta);
circle_fit=[x,y,z];
vectors_onplane=circle_fit-center;
