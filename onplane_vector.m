
%��ƽ����Ѱ��0-360�ȷ�Χ���������vectors
%���룺ƽ�����input_pnts��nx3��   ����������radius��ȡ1���ɣ����Ƕȼ��resolution(��λ�Ƕ�)
%�����vectors_onplane��nx3���ǵ�λ����
function [vectors_onplane] = onplane_vector(input_pnts,radius,resolution)
  [parameter] = TLS_Plane(input_pnts);
  center=[0,0,-parameter(4)/parameter(3)];
 r=radius;
 c=center;
 n=parameter(1:3,:)';
 theta_interval=pi*resolution/180;
 theta=(0:theta_interval:2*pi)';
 
a=cross(n,[1 0 0]); %n��i��ˣ���ȡa����
if ~any(a) %���aΪ����������n��j���
a=cross(n,[0 1 0]);
end
b=cross(n,a); %��ȡb����
a=a/norm(a); %��λ��a����
b=b/norm(b);
c1=c(1)*ones(size(theta,1),1);
c2=c(2)*ones(size(theta,1),1);
c3=c(3)*ones(size(theta,1),1);
x=c1+r*a(1)*cos(theta)+r*b(1)*sin(theta);%Բ�ϸ����x����
y=c2+r*a(2)*cos(theta)+r*b(2)*sin(theta);%Բ�ϸ����y����
z=c3+r*a(3)*cos(theta)+r*b(3)*sin(theta);
circle_fit=[x,y,z];
vectors_onplane=circle_fit-center;
