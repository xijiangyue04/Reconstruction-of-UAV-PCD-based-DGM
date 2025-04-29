%�ڶ��������룺line_number Ϊ��ʼ�ĺ������������򻯼��������ߣ�Ĭ����4��, dist_thr�����߻������ Ĭ����0.1
%sort_R�������ܶȣ�line_vectorÿ�������ߵ�������boundary_pnts �߽��
%�����horizon_regula ����ˮƽ������������ߣ�vertical_regula ���մ�ֱ������������ߣ�
%horizon_insertpnts ˮƽ�����޲����ڲ�㣬vertical_insertpnts ��ֱ�����޲����ڲ��
%node_pnts ����

function [horizon_regula,vertical_regula,horizon_insertpnts,vertical_insertpnts,node_pnts] = verticalcontour_frame_regulariz(line_number,sort_R,boundary_pnts,project_pnt,line_vector,dist_thr) 

%�����������������ܶȶ�λ��Ͷ�λ�����Ժ������������߽��й���
 multi_line_pnt=boundary_pnts(sort_R{1}(1:line_number,1),:);
horizontal_regulariz=[];
for i=1:size(multi_line_pnt,1)
    scale=sqrt(sum(multi_line_pnt(i,:).^2));%���õ�ľ���������ɨ����Ƶĳ߶ȣ��߶Ⱥܴ������±���10^6,��һ��С��������ӻ���ֽ��û�仯,����Ҫ���Գ߶�
    line_pnt1=multi_line_pnt(i,:);%
    line_pnt2=line_pnt1+scale*line_vector(1,:);
    linepnt=[line_pnt1;line_pnt2];
    dist= PL_distance_TLS(boundary_pnts,  line_pnt1, line_vector(1,:)); %��-�߾���
    outpnt=boundary_pnts(dist<dist_thr,:);   %��ֵ1���������������߻�������Χ0.1���ڵĵ�
    segment{1}=outpnt;
    resolution=0.005; %��ֵ2���������ǹ��򻯺������ڵ���Ϊ0.005��
    chose=1;
    number=10; %��ֵ3���������ǹ�������LFPR������Ҫ����10����Ŷ������ִ��
    if size(outpnt,1)>10
    Line_feature_segment = LFPR(segment,resolution,chose,number);
    horizontal_regulariz=[Line_feature_segment,horizontal_regulariz];
    else
    Line_feature_segment=[];
    end
end
    
    
%�����������������ܶȶ�λ��Ͷ�λ�������������������߽��й���
 multi_line_pnt=boundary_pnts(sort_R{2}(1:line_number,1),:);
vertical_regulariz=[];
for i=1:size(multi_line_pnt,1)
    scale=sqrt(sum(multi_line_pnt(i,:).^2));%���õ�ľ���������ɨ����Ƶĳ߶ȣ��߶Ⱥܴ������±���10^6,��һ��С��������ӻ���ֽ��û�仯,����Ҫ���Գ߶�
    line_pnt1=multi_line_pnt(i,:);%
    line_pnt2=line_pnt1+scale*line_vector(2,:);
    linepnt=[line_pnt1;line_pnt2];
    dist= PL_distance_TLS(boundary_pnts,  line_pnt1, line_vector(2,:)); %��-�߾���
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

 %|���ǻ�

if length(horizontal_regulariz)>1 & length(vertical_regulariz)>1

%ȷ��ÿ�����������ߵĵ㼰����
  for i=1:length(horizontal_regulariz)
    [horizontal_line_vector(i,:),horizontal_mean_pnt(i,:)] = space_line_TLS(horizontal_regulariz{i});
  end

  for i=1:length(vertical_regulariz)
    [vertical_line_vector(i,:),vertical_mean_pnt(i,:)] = space_line_TLS(vertical_regulariz{i});
  end

%���ϸ�����line_number�������ĸ�����������������������߶���4��������������������-�߾�����Զԭ�����Ϊ2����
  for i=1:size(horizontal_mean_pnt,1)
    [distance1] = PL_distance_TLS(horizontal_mean_pnt, horizontal_mean_pnt(i,:), horizontal_line_vector(i,:));
    origmax_dis1(i,:)=max(distance1);%�������������ߵ���Ŀ�����µ���Զ����
    index=find(distance1==max(distance1));
    origmax_index1(i,:)=[i,index(1)];% ȷ����Զ�����Ŀ�����������߱�ǩ��
  end

  for i=1:size(vertical_mean_pnt,1)
    [distance2] = PL_distance_TLS(vertical_mean_pnt, vertical_mean_pnt(i,:), vertical_line_vector(i,:));%�����������򵽸�Ŀ�����µ���Զ����
    origmax_dis2(i,:)=max(distance2);%�����������򵽸�Ŀ�����µ���Զ����
    index=find(distance2==max(distance2));
    origmax_index2(i,:)=[i,index(1)];% ȷ����Զ�����Ŀ�����������߱�ǩ��
  end

    horizon_index=origmax_index1(origmax_dis1==max(origmax_dis1),:);%ˮƽ�߾�����Զ�������������߱�ǩ��
    vertical_index=origmax_index2(origmax_dis2==max(origmax_dis2),:);%��ֱ�߾�����Զ�������������߱�ǩ��
  for i=1:2
    horizon_regula{i}= horizontal_regulariz{horizon_index(i)};%������Զ������������������
    vertical_regula{i}= vertical_regulariz{vertical_index(i)};%������Զ������������������
  end


vertical_insertpnts=[];
horizon_insertpnts=[];
intersectpnts=[];
node_pnts=[];

%��һ�ַ��� :
%�����Ƕ�������������������������յĲ��ֽ��к����������۵ĵ�һ�ַ������ڲ��

  for i=1:length(horizon_regula) %����ֵ�Ǻ�������������
    P1=horizon_regula{i}(1,:);
    v1=line_vector(1,:);
    for j=1:length(vertical_regula)%����ֵ����������������
      P2=vertical_regula{j}(1,:);
      v2=line_vector(2,:);
      intersection = skewline_perpendicular_foot(P1,v1,P2,v2);
      intersect_pnt=intersection(1,:); %ˮƽ�ʹ�ֱ�����ߵĽ���
      dist1=sqrt(sum((horizon_regula{i}-intersect_pnt).^2,2)); %���㵽ˮƽ�����������е����
      min_dist1=min(dist1); %���㵽ˮƽ�����������е���С����
      horizon_end=horizon_regula{i}(dist1==min_dist1,:); %ˮƽ�������뽻������Ķ˵�
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



%�������ж����Ϸ����õ��ڲ������߽��߰���ı߽����ʣ���������˴󲿷ֱ߽�����ռ�ȳ���70%��������������淽���������ڲ�
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
    vec_degree1=degree*180/pi;%0-180�㷶Χ
    vec_degree2=min(degree/pi*180,180-degree/pi*180);%0-90�㷶Χ
        if 180-max(vec_degree1)<2
           vec_obj=vec(vec_degree2==max(vec_degree2),:);
           vec_obj=repmat(vec_obj(1,:),m,1);
           degree_obj=acos(dot(vec_obj',vec'));
           vec_degree_obj=degree_obj*180/pi;%0-180�㷶Χ
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




%�ڶ��ַ���
 %�����Ƕ�������������������������յĲ��ֽ��к��������ֱ�������ڲ�ĵڶ��ַ������ڲ��
    if ratio<0.7 %��ֵ4������պ������߰���ı߽��С��70%
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
        intersect_pnt=intersection(1,:); %ˮƽ�ʹ�ֱ�����ߵĽ���
        dist1=sqrt(sum((horizon_regula{i}-intersect_pnt).^2,2)); %���㵽ˮƽ�����������е����
        min_dist1=min(dist1); %���㵽�����������е���С����
        horizon_end=horizon_regula{i}(dist1==min_dist1,:); %ˮƽ�������뽻������Ķ˵�
        dist2=sqrt(sum((vertical_regula{j}-intersect_pnt).^2,2));
        min_dist2=min(dist2);%���㵽�����������е���С����

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
    disp('����������̫��');
end

end









