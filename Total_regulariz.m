% number_of_neighbor=6;line_number=4; resolution=0.005;dist_thr=0.1
%输入：平面分割结果路径，该路径下保存的是某个建筑物的不同平面  number_of_neighbor：用于边界点提取的最邻近点数 
%line_number 为初始的横向或者纵向规则化几条轮廓线
%输出：eachregularize 每个规则化的轮廓线
function [vertical_regularize,boundarypnts,horizontal_regularize,eachregularize,node] = Total_regulariz(number_of_neighbor,line_number,dist_thr) 
%以下是依据轮廓转折拐点连接插入实现的规则化
file=dir('E:\三维激光论文\无人机建筑点云建模\图\不同方法的比较\提取的顶面点云数据\1\*.txt');
n=length(file);
resolution=0.005;
for k=1:n
    horizontal_regularize=[];vertical_regularize=[]; 
filename=['E:\三维激光论文\无人机建筑点云建模\图\不同方法的比较\提取的顶面点云数据\1\',num2str(k),'.txt'];
plane_segment=load(filename);
plane_segment=plane_segment(:,1:3);
[parameter] = TLS_Plane(plane_segment);
vector1=parameter(1:3,:)';
vector2=[0,0,1];
dotProduct = dot(vector1, vector2);
angle= acos(dotProduct) * (180 / pi);
  if length(plane_segment)>500
    if angle >=45
      [sort_R,boundary_pnts,line_vector,project_pnt] = verticalcontour_density(plane_segment,number_of_neighbor,dist_thr);
      [horizon_regula,vertical_regula,horizon_insertpnts,vertical_insertpnts,node_pnts] = verticalcontour_frame_regulariz(line_number,sort_R,boundary_pnts,project_pnt,line_vector,dist_thr);
      node{k}=node_pnts;
      node_k=size(node_pnts,1);
      if node_k>2
          for j=1:node_k-1
             insert_pnts1= interpolation_pnts(node_pnts(j,:),node_pnts(j+1,:),resolution);
             vertical_regularize=[vertical_regularize;insert_pnts1];
          end
             insert_pnts2= interpolation_pnts(node_pnts(node_k,:),node_pnts(1,:),resolution);%最后一个点与第一个点连接线
             vertical_regularize=[vertical_regularize;insert_pnts2];
      end
      
    
    else
        [ultimatevector,boundarypnts,clust,cluster_pnts] = horizontalcontour_vector(plane_segment,number_of_neighbor,dist_thr);
        [sort_R,boundary_pnts,line_vector,project_pnt] = horizontalcontour_density2(plane_segment,number_of_neighbor,ultimatevector,dist_thr);
        [horizon_regula,vertical_regula,horizon_insertpnts,vertical_insertpnts,node_pnts] = horizontalcontour_frame_regulariz(line_number,sort_R,boundary_pnts,line_vector,project_pnt,dist_thr);
         node{k}=node_pnts;
         node_k=size(node_pnts,1);
      if node_k>2
          for j=1:node_k-1
              insert_pnts1= interpolation_pnts(node_pnts(j,:),node_pnts(j+1,:),resolution);
              horizontal_regularize=[horizontal_regularize;insert_pnts1];
          end
             insert_pnts2 = interpolation_pnts(node_pnts(node_k,:),node_pnts(1,:),resolution);
             horizontal_regularize=[horizontal_regularize;insert_pnts2];
      end
    end
  end
  eachregularize{k}=[vertical_regularize;horizontal_regularize];
end






















