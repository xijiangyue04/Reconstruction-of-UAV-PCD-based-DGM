%用来根据包围球中的点数来判断是否噪声点
%输入变量input(nx3)  球的半径radius
%输出变量number_points 球中的点数
function [pnt_number] =  sphere_points_density(input,radius)
n=size(input,1);
for i=1:n
  %先利用长度为radius正方形方框确定大致范围
line_xyz=find(input(:,1)>input(i,1)-radius & input(:,1)<input(i,1)+radius & input(:,2)>input(i,2)-radius & input(:,2)<input(i,2)+radius & input(:,3)>input(i,3)-radius & input(:,3)<input(i,3)+radius);
input_xyz=input(line_xyz,:);
k=size(input_xyz,1);
 D=sqrt(sum((input_xyz-repmat(input(i,:),k,1)).^2,2));
 sphere_idx=find(D<radius);
 %points_line{i}=line_xyz(sphere_idx);
 points=input(line_xyz(sphere_idx),:);
 pnt_number(i,:)=size(sphere_idx,1);
end


    