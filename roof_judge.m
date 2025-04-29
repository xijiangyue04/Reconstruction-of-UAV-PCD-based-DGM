%对每个建筑物提取的平面点云进行顶面点云数据的判断
%输入：每个建筑物提取保存的平面路径 path = 'E:\三维激光论文\无人机点云轮廓点提取-准备投那个drones特刊\2无人机点云数据平面分割结果\B3平面分割\';
%输出：Prop 水平面点云数据属于屋顶顶面点云的概率

function [Prob,horizontal_pnts,vertical_pnts, roof_pnts,nonroof_pnts] = roof_judge(path,Pr) 
namelist = dir([path,'*.txt']);
n= length(namelist);
roof_pnts=[];nonroof_pnts=[];
for k=1:n
filename = [path,namelist(k).name];
%下面三行代码确定的是读取的txt文件名是1还是2，3,4之类的
a=num2str(namelist(k).name);
s=a(isstrprop(a,'digit')); 
kk=str2num(s);
%下面是对导入txt进行分析
plane_segment=load(filename);
[parameter] = TLS_Plane(plane_segment(:,1:3));
vector1=parameter(1:3,:)';
vector2=[0,0,1];
dotProduct = dot(vector1, vector2);
angle= acos(dotProduct) * (180 / pi);
  if angle <45   %如果与垂直面的角度大于45，则认为是水平面
input=plane_segment(:,1:3);
[pnt_density] = sphere_points_density(input,0.2);%点云密度
pnts_number=size(input,1);%点云数
mean_z=mean(input(:,3));%z值
totalpnt_number(k,:)=[sum(pnt_density)/size(input,1),pnts_number,mean_z,kk];
horizontal_pnts{k}=plane_segment;
  else
      vertical_pnts{k}=plane_segment;
  end
end
horizontal_pnts(cellfun(@isempty,horizontal_pnts))=[];
 vertical_pnts(cellfun(@isempty, vertical_pnts))=[];
totalpnt_number(all(totalpnt_number==0,2),:)=[];
normal1=1/3*(totalpnt_number(:,1)/max(totalpnt_number(:,1)));
normal2=1/3*(totalpnt_number(:,2)/max(totalpnt_number(:,2)));
normal3=1/3*(totalpnt_number(:,3)/max(totalpnt_number(:,3)));
totalpnt_num=[normal1,normal2,normal3];
Prob=[sum(totalpnt_num(:,1:3),2),totalpnt_number(:,4)];

for i=1:size(Prob,1)
    if Prob(i,1)>=Pr
        roof_pnts{i}=horizontal_pnts{i};
    else
        nonroof_pnts{i}=horizontal_pnts{i};
    end
end
if ~isempty(roof_pnts)
roof_pnts(cellfun(@isempty,roof_pnts))=[];
end
if ~isempty(nonroof_pnts)
nonroof_pnts(cellfun(@isempty,nonroof_pnts))=[];
end

