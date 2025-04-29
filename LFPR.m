%第二步 这一步是对于特征线点云进行规则化处理
%特征点规则化处理  输入不同类别的折边点sort_segment，为元胞数组形式
%当chose为 1时，采用的是整体最小二乘空间直线拟合方程，参考文献：秦锋,王涛,张振虎.空间直线拟合新方法[J].地理空间信息,2023,21(03):21-24.
%当chose为 0时，采用的是最小二乘空间直线拟合方程
%resolution为多少距离为一个点 即规则化后的线性轮廓点云的距离分辨率    number 为保留的大于该数据量的分割点云

function  [Line_feature_segment] = LFPR(segment,resolution,chose,number)   %Line Feature point cloud regularization(LFPR)
for i=1:length(segment)
    segment_n(i,:)=[size(segment{i},1),i];
end
    sort_segment_n=sortrows(segment_n,1,'descend');
    
 for i=1:length(segment)
     sort_segment{i}=segment{sort_segment_n(i,2)};
 end
 index=max(find(sort_segment_n(:,1)>=number)); %对于少于10个点的分割点云数据要给予删除
 

for i=1:index
    number_neighbor=size(sort_segment{i},1);
    if number_neighbor<number
        continue
    else
    [parameter_ls] = space_line_LS(sort_segment{i});
    [TLS_vector,mean_pnt] = space_line_TLS(sort_segment{i});
    TLS_vector=TLS_vector';
    end
    [idx,dist]=knnsearch(sort_segment{i}, sort_segment{i}, 'k', number_neighbor);
    max_dist=max(dist(:,number_neighbor));
    fat_idx_start=idx(dist(:,number_neighbor)==max_dist,1);
    fat_idx_end=idx(dist(:,number_neighbor)==max_dist,number_neighbor);
    fat_idx=[fat_idx_start,fat_idx_end];%第一列是初始点索引，第二列是终点索引
    fat_idx=transpose(sort(fat_idx')); %第一列是初始点索引，第二列是终点索引
    fat_idx=unique(fat_idx,'rows');%剔除重复的最远点索引，第一列是初始点索引，第二列是终点索引
    if size(fat_idx,1)>1
        fps_n=size(fat_idx,1);
        for j=1:fps_n
        fps_pnt_tempory=sort_segment{i}(fat_idx(j,:),:);  
        [PL_dis] = PL_distance_LS(fps_pnt_tempory,parameter_ls);
        max_PLdis(j,:)=max(PL_dis);
        end
        fat_idx_finally=fat_idx(max_PLdis==min(max_PLdis),:);
    else
        fat_idx_finally=fat_idx;
    end
    fps_pnt=sort_segment{i}(fat_idx_finally,:); %最远的两个点
    fps_pnt=unique(fps_pnt,'rows');
    fps_n=size(fps_pnt,1);
    
 %以下命令是确定最远点在整体最小二乘得到的空间直线上的投影点
 if chose==1
    for k=1:fps_n
      A = [fps_pnt(k,1); fps_pnt(k,2); fps_pnt(k,3)];
      x0=A(1);
      y0=A(2);
      z0=A(3);
      v = TLS_vector;
      PA = [x0; y0; z0] - mean_pnt';
      proj_v_PA = dot(PA, v) / norm(v)^2 * v;
      PQ = proj_v_PA;
      Q = mean_pnt' + PQ;
      x_Q = Q(1);
      y_Q = Q(2);
      z_Q = Q(3);
      Project_pnt(k,:)=[x_Q,y_Q,z_Q];
    end
 else
    
  %以下命令是确定最远点在最小二乘得到的空间直线上的投影点
    a=parameter_ls(1);
    b=parameter_ls(2);
    c=parameter_ls(3);
    d=parameter_ls(4);
    for k=1:fps_n
    A = [fps_pnt(k,1); fps_pnt(k,2); fps_pnt(k,3)];%提取其中一个最远点
    x0=A(1);
    y0=A(2);
    z0=A(3);
     v = [a; c; 1];
    PA = [x0; y0; z0] - [b; d; 0];
    proj_v_PA = dot(PA, v) / norm(v)^2 * v;
    PQ = proj_v_PA;
    Q = [b; d; 0] + PQ;
    x_Q = Q(1);
    y_Q = Q(2);
    z_Q = (x_Q-b)/a;
    Project_pnt(k,:)=[x_Q,y_Q,z_Q];
    end
 end
    
    
    
    fps_dis=sqrt(sum((Project_pnt(1,:)-Project_pnt(2,:)).^2,2)); %投影后最远点距离
    insert_number=fix(fps_dis/resolution); %0.02m插入一个点,总共插入多少个点
    if insert_number==0
        insert_number=1;
    end
    fps_interval=Project_pnt(1,:)-Project_pnt(2,:);%最远两个点间隔
    neighbor_interval=fps_interval/insert_number;%插入点相邻两个点间隔
    for num=1:insert_number %之所以这样是希望向前和向后多插入点
        insert_pnts(num,:)=Project_pnt(2,:)+neighbor_interval*(num-2);%插入的点坐标
    end
    Line_feature_segment{i}=insert_pnts;
    insert_pnts=[];
end

    
   

  


    
    
  



    