%判断点是否在点云内部的判断，如果是三维需要将三维点投影到一个平面上
%需要判断的点judgepnt(mx3)，范围边界点云bundary_pnts(nx3)
%如果max_diff>90 左右 则判断该点在点云外部,根据judgepnt中的点与bundary_pnts中的每个点组成向量的方位角差值进行判断

function [max_diff] = pntinpcd(judgepnt,bundary_pnts)   
n1=size(judgepnt,1);n2=size(bundary_pnts,1);
total=[judgepnt;bundary_pnts];
n=n1+n2;
mean_neighbor=mean(total,1);  
[parameter] = TLS_Plane(total);
normal_pnts =parameter(1:3,:);
lamda=(repmat(mean_neighbor*normal_pnts,n,1)-total*normal_pnts)/(sum(normal_pnts.^2));% kx1列矩阵
proj_total=normal_pnts*lamda'+total';%3xk
proj_judgepnt=proj_total(:,1:n1);%3xn1
proj_bundary=proj_total(:,n1+1:end);%3xn2
for i=1:n1
    %judgepnt中的点与边界所有点构成向量的方位角
   % k=size(bundary_pnts,1);
    
    boundarypnts=proj_bundary;
    object_pnt=proj_judgepnt(:,i)';
    x_unit_vector=(boundarypnts(:,end)'-object_pnt)/norm(boundarypnts(:,end)'-object_pnt);% 通过judgepnt中的一个目标点object_pnt和边界点boundarypnts(:,end)'构建x轴法向量 1x3 
    y_unit_vector=cross(normal_pnts,x_unit_vector)/norm(cross(normal_pnts,x_unit_vector)); %根据两个向量叉乘除以叉乘的模来确定另一个坐标轴上的单位向量 1x3
 
    x_plane=sum(repmat(x_unit_vector,n2-1,1).*(boundarypnts(:,2:end)'-repmat(object_pnt,n2-1,1)),2); % 根据单位向量乘以向量得到在平面x轴上的投影kx1
   y_plane=sum(repmat(y_unit_vector,n2-1,1).*(boundarypnts(:,2:end)'-repmat(object_pnt,n2-1,1)),2); % 根据单位向量乘以向量得到在平面y轴上的投影kx1
 
   azimuth= rad2deg(atan2(y_plane,x_plane)); %计算坐标方位角，如果是在第一和第二象限则是正值，如果是第三和第四象限则是负值
   azimuth(azimuth<0)= azimuth(azimuth<0)+360 ; %将坐标方位角的负值全部变为正值
   azimuth=[azimuth;360];
   sort_azimuth=sort(azimuth);
   diff_azimuth=diff(sort_azimuth);
   max_diff(i,:)=max(diff_azimuth);
end





% %判断点是否在点云内部的判断（一般是二维情况下），如果是三维需要将三维降为二维
% %需要判断的点input(mx3)，范围点云input_pnts(nx3)
% %如果max_diff>180 这判断该点在点云内部
% 
% function [max_diff] = pntinpcd(input,input_pnts)   
% for i=1:size(input,1)
%     %以下命令是从boundary_extract2 命令中提取的目标点与周围点构成向量的方位角
%     normal_pnts=[0;0;1];
%     k=size(input_pnts,1);
%     proj_pnts=input_pnts';
%     object_pnt=input(i,:);
%     x_unit_vector=(proj_pnts(:,end)'-object_pnt)/norm(proj_pnts(:,end)'-object_pnt);% 通过探测点input_pnts(i,:)和投影邻近点最远的点proj_pnts(:,end)'构建x轴法向量 1x3 
%     y_unit_vector=cross(normal_pnts,x_unit_vector)/norm(cross(normal_pnts,x_unit_vector)); %根据两个向量叉乘除以叉乘的模来确定另一个坐标轴上的单位向量 1x3
%  
%     x_plane=sum(repmat(x_unit_vector,k-1,1).*(proj_pnts(:,2:end)'-repmat(object_pnt,k-1,1)),2); % 根据单位向量乘以向量得到在平面x轴上的投影kx1
%    y_plane=sum(repmat(y_unit_vector,k-1,1).*(proj_pnts(:,2:end)'-repmat(object_pnt,k-1,1)),2); % 根据单位向量乘以向量得到在平面y轴上的投影kx1
%  
%    azimuth= rad2deg(atan2(y_plane,x_plane)); %计算坐标方位角，如果是在第一和第二象限则是正值，如果是第三和第四象限则是负值
%    azimuth(azimuth<0)= azimuth(azimuth<0)+360 ; %将坐标方位角的负值全部变为正值
%    sort_azimuth=sort(azimuth);
%    diff_azimuth=diff(sort_azimuth);
%    max_diff(i,:)=max(diff_azimuth);
% end