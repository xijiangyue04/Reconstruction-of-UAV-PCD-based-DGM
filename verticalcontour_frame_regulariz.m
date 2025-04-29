%第二步：输入：line_number 为初始的横向或者纵向规则化几条轮廓线，默认是4条, dist_thr轮廓线缓冲距离 默认是0.1
%sort_R轮廓线密度，line_vector每个轮廓线的向量，boundary_pnts 边界点
%输出：horizon_regula 最终水平方向规则化轮廓线；vertical_regula 最终垂直方向规则化轮廓线，
%horizon_insertpnts 水平方向修补的内插点，vertical_insertpnts 垂直方向修补的内插点
%node_pnts 交点

function [horizon_regula,vertical_regula,horizon_insertpnts,vertical_insertpnts,node_pnts] = verticalcontour_frame_regulariz(line_number,sort_R,boundary_pnts,project_pnt,line_vector,dist_thr) 

%以下是依据轮廓线密度定位点和定位向量对横向主体轮廓线进行规则化
 multi_line_pnt=boundary_pnts(sort_R{1}(1:line_number,1),:);
horizontal_regulariz=[];
for i=1:size(multi_line_pnt,1)
    scale=sqrt(sum(multi_line_pnt(i,:).^2));%利用点的距离来描述扫描点云的尺度，尺度很大的情况下比如10^6,与一个小的向量相加会出现结果没变化,所以要乘以尺度
    line_pnt1=multi_line_pnt(i,:);%
    line_pnt2=line_pnt1+scale*line_vector(1,:);
    linepnt=[line_pnt1;line_pnt2];
    dist= PL_distance_TLS(boundary_pnts,  line_pnt1, line_vector(1,:)); %点-线距离
    outpnt=boundary_pnts(dist<dist_thr,:);   %阈值1：描述的是轮廓线缓冲区范围0.1以内的点
    segment{1}=outpnt;
    resolution=0.005; %阈值2：描述的是规则化后点的相邻点间距为0.005米
    chose=1;
    number=10; %阈值3：描述的是规则化命令LFPR，至少要大于10个点才对其进行执行
    if size(outpnt,1)>10
    Line_feature_segment = LFPR(segment,resolution,chose,number);
    horizontal_regulariz=[Line_feature_segment,horizontal_regulariz];
    else
    Line_feature_segment=[];
    end
end
    
    
%以下是依据轮廓线密度定位点和定位向量对纵向主体轮廓线进行规则化
 multi_line_pnt=boundary_pnts(sort_R{2}(1:line_number,1),:);
vertical_regulariz=[];
for i=1:size(multi_line_pnt,1)
    scale=sqrt(sum(multi_line_pnt(i,:).^2));%利用点的距离来描述扫描点云的尺度，尺度很大的情况下比如10^6,与一个小的向量相加会出现结果没变化,所以要乘以尺度
    line_pnt1=multi_line_pnt(i,:);%
    line_pnt2=line_pnt1+scale*line_vector(2,:);
    linepnt=[line_pnt1;line_pnt2];
    dist= PL_distance_TLS(boundary_pnts,  line_pnt1, line_vector(2,:)); %点-线距离
    outpnt=boundary_pnts(dist<dist_thr,:);
    segment{1}=outpnt;
    resolution=0.005;
    chose=1;
    number=10;
    if  size(outpnt,1)>10
    Line_feature_segment = LFPR(segment,resolution,chose,number);
    vertical_regulariz=[Line_feature_segment,vertical_regulariz];
    else
     Line_feature_segment =[];
    end
end

 %|这是或

if length(horizontal_regulariz)>1 & length(vertical_regulariz)>1

%确定每个主体轮廓线的点及向量
  for i=1:length(horizontal_regulariz)
    [horizontal_line_vector(i,:),horizontal_mean_pnt(i,:)] = space_line_TLS(horizontal_regulariz{i});
  end

  for i=1:length(vertical_regulariz)
    [vertical_line_vector(i,:),vertical_mean_pnt(i,:)] = space_line_TLS(vertical_regulariz{i});
  end

%以上给出的line_number可能是四个，即横向和竖向主体轮廓线都有4个，以下命令是依据线-线距离最远原则将其变为2个。
  for i=1:size(horizontal_mean_pnt,1)
    [distance1] = PL_distance_TLS(horizontal_mean_pnt, horizontal_mean_pnt(i,:), horizontal_line_vector(i,:));
    origmax_dis1(i,:)=max(distance1);%计算其他横向线到该目标线下的最远距离
    index=find(distance1==max(distance1));
    origmax_index1(i,:)=[i,index(1)];% 确定最远距离的目标线与其他线标签号
  end

  for i=1:size(vertical_mean_pnt,1)
    [distance2] = PL_distance_TLS(vertical_mean_pnt, vertical_mean_pnt(i,:), vertical_line_vector(i,:));%计算其他竖向到该目标线下的最远距离
    origmax_dis2(i,:)=max(distance2);%计算其他竖向到该目标线下的最远距离
    index=find(distance2==max(distance2));
    origmax_index2(i,:)=[i,index(1)];% 确定最远距离的目标线与其他线标签号
  end

    horizon_index=origmax_index1(origmax_dis1==max(origmax_dis1),:);%水平线距离最远的那两个横向线标签号
    vertical_index=origmax_index2(origmax_dis2==max(origmax_dis2),:);%垂直线距离最远的那两个竖向线标签号
  for i=1:2
    horizon_regula{i}= horizontal_regulariz{horizon_index(i)};%距离最远的那两个横向轮廓线
    vertical_regula{i}= vertical_regulariz{vertical_index(i)};%距离最远的那两个竖向轮廓线
  end


vertical_insertpnts=[];
horizon_insertpnts=[];
intersectpnts=[];
node_pnts=[];

%第一种方案 :
%以下是对两个横向和竖向主体轮廓见空的部分进行横向和竖向拐折的第一种方案的内插点

  for i=1:length(horizon_regula) %该数值是横向两个轮廓线
    P1=horizon_regula{i}(1,:);
    v1=line_vector(1,:);
    for j=1:length(vertical_regula)%该数值是竖向两个轮廓线
      P2=vertical_regula{j}(1,:);
      v2=line_vector(2,:);
      intersection = skewline_perpendicular_foot(P1,v1,P2,v2);
      intersect_pnt=intersection(1,:); %水平和垂直轮廓线的交点
      dist1=sqrt(sum((horizon_regula{i}-intersect_pnt).^2,2)); %交点到水平方向轮廓所有点距离
      min_dist1=min(dist1); %交点到水平方向轮廓所有点最小距离
      horizon_end=horizon_regula{i}(dist1==min_dist1,:); %水平方向线离交点最近的端点
      dist2=sqrt(sum((vertical_regula{j}-intersect_pnt).^2,2));
      min_dist2=min(dist2);
      vertical_end=vertical_regula{j}(dist2==min_dist2,:);
        if min_dist1>0.1 & min_dist2>0.1
           intersection = skewline_perpendicular_foot(horizon_end,v2,vertical_end,v1);
           intersect_pnt=[horizon_end;intersection(1,:);vertical_end];
           insert_vertical_pnts = interpolation_pnts(horizon_end,intersection(1,:),resolution); 
           insert_horizon_pnts = interpolation_pnts(vertical_end,intersection(1,:),resolution);
        end

        if min_dist1>0.1 & min_dist2<=0.1
           insert_horizon_pnts = interpolation_pnts(horizon_end,intersection(1,:),resolution); 
           insert_vertical_pnts =[];
        end
       if min_dist1<=0.1 & min_dist2>0.1
           insert_vertical_pnts = interpolation_pnts(vertical_end,intersection(1,:),resolution);
           insert_horizon_pnts=[];
       end
       if min_dist1<=0.1 & min_dist2<=0.1
           insert_horizon_pnts=[];
           insert_vertical_pnts=[];
       end
    horizon_insertpnts=[horizon_insertpnts;insert_horizon_pnts];
    vertical_insertpnts=[vertical_insertpnts;insert_vertical_pnts];
       if j==1
         intersect_pnt=flipud(intersect_pnt);
       end
    intersectpnts=[intersectpnts;intersect_pnt];
    end
       if i==2
         intersectpnts=flipud(intersectpnts);
       end
    node{i}=intersectpnts;
    intersectpnts=[];
  end



%以下是判断以上方案得到内插轮廓边界线包络的边界点比率，如果包络了大部分边界点比如占比超过70%，则无需进行下面方案的轮廓内插
    regular_pnts=[];
    for i=1:length(horizon_regula)
        regular_pnts=[regular_pnts;horizon_regula{i};vertical_regula{i}];
    end
    regular_pnts=[regular_pnts;horizon_insertpnts;vertical_insertpnts];

    m=size(regular_pnts,1);
    for i=1:size(project_pnt,1)
    vector=regular_pnts-project_pnt(i,:);
    vec=vector./sqrt(sum(vector.^2,2));
    vec1=repmat(vec(1,:),m,1);
    degree=acos(dot(vec1',vec'));
    vec_degree1=degree*180/pi;%0-180°范围
    vec_degree2=min(degree/pi*180,180-degree/pi*180);%0-90°范围
        if 180-max(vec_degree1)<2
           vec_obj=vec(vec_degree2==max(vec_degree2),:);
           vec_obj=repmat(vec_obj(1,:),m,1);
           degree_obj=acos(dot(vec_obj',vec'));
           vec_degree_obj=degree_obj*180/pi;%0-180°范围
           if 180-max(vec_degree_obj)<2
               judge(i,:)=1;
           else
               judge(i,:)=0;
           end
        else
            judge(i,:)=0;
        end
    end
ratio=sum(judge(:)==1)/length(judge);




%第二种方案
 %以下是对两个横向和竖向主体轮廓见空的部分进行横向和竖向直接连接内插的第二种方案的内插点
    if ratio<0.7 %阈值4：如果闭合轮廓线包络的边界点小于70%
    vertical_insertpnts=[];
    horizon_insertpnts=[];
    intersectpnts=[];
        for i=1:2
        P1=horizon_regula{i}(1,:);
        v1=line_vector(1,:);
        for j=1:2
        P2=vertical_regula{j}(1,:);
        v2=line_vector(2,:);
        intersection = skewline_perpendicular_foot(P1,v1,P2,v2);
        intersect_pnt=intersection(1,:); %水平和垂直轮廓线的交点
        dist1=sqrt(sum((horizon_regula{i}-intersect_pnt).^2,2)); %交点到水平方向轮廓所有点距离
        min_dist1=min(dist1); %交点到横向轮廓所有点最小距离
        horizon_end=horizon_regula{i}(dist1==min_dist1,:); %水平方向线离交点最近的端点
        dist2=sqrt(sum((vertical_regula{j}-intersect_pnt).^2,2));
        min_dist2=min(dist2);%交点到竖向轮廓所有点最小距离

        vertical_end=vertical_regula{j}(dist2==min_dist2,:);
        if min_dist1>0.1 & min_dist2>0.1
            insert_horizon_pnts = interpolation_pnts(horizon_end,intersection(1,:),resolution); 
            insert_vertical_pnts = interpolation_pnts(vertical_end,intersection(1,:),resolution);
        end

        if min_dist1>0.1 & min_dist2<=0.1
            insert_horizon_pnts = interpolation_pnts(horizon_end,intersection(1,:),resolution); 
            insert_vertical_pnts =[];
        end
        if min_dist1<=0.1 & min_dist2>0.1
            insert_vertical_pnts = interpolation_pnts(vertical_end,intersection(1,:),resolution);
            insert_horizon_pnts=[];
        end
        if min_dist1<=0.1 & min_dist2<=0.1
            insert_horizon_pnts=[];
            insert_vertical_pnts=[];
        end
        horizon_insertpnts=[horizon_insertpnts;insert_horizon_pnts];
        vertical_insertpnts=[vertical_insertpnts;insert_vertical_pnts] ;
        if j==1
            intersect_pnt=flipud(intersect_pnt);
        end
        intersectpnts=[intersectpnts;intersect_pnt];
        end

         if i==2
             intersectpnts=flipud(intersectpnts);
         end
        node{i}=intersectpnts;
        intersectpnts=[];
        end
    end

    for k=1:length(node)
        node_pnts=[node_pnts;node{k}];
    end

else
    horizon_regula=[];
    vertical_regula=[];
    horizon_insertpnts=[];
    vertical_insertpnts=[];
    node_pnts=[];
    disp('横向或竖向点太少');
end

end









