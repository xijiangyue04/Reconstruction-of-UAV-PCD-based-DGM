%第二步：确定每个分割平面竖向边界点的横向和纵向密度最大的，且相距较远的定位点
%说白了就是要确定边界点的主体轮廓框架定位点及定位方向向量
%sort_R代表不同向量ultimatevector包络所有点的所有密度
%输出sort_R(元包形式)，boundary_pnts边界点(矩阵形式),line_vector每个轮廓线向量(矩阵形式)；project_pnt点云的二维投影点(矩阵形式)
function [sort_R,boundary_pnts,line_vector,project_pnt] = horizontalcontour_density2(plane_segment,number_of_neighbor,ultimatevector,dist_thr)


[parameter] = TLS_Plane(plane_segment);
[project_pnt] = pntplane_projection(plane_segment,parameter) ;
project_pnt=project_pnt';
[boundary_pnts] = boundary_extract2(project_pnt,number_of_neighbor) ;
line_vector=ultimatevector;
n=size(boundary_pnts,1);

for j=1:size(line_vector,1)
%以下是备用的       
for i=1:n
    dist= PL_distance_TLS(boundary_pnts, boundary_pnts(i,:), line_vector(j,:)); %每个边界点与主体方向向量构成的直线，所有边界点-该直线的距离
    pnt_inter=boundary_pnts(dist<dist_thr,:);%其他点到线的距离范围小于0.1的点数
    if size(pnt_inter,1)>5 %该线范围大于5个点才需要进行密度计算
    scale=sqrt(sum(boundary_pnts(i,:).^2));%利用点的距离来描述扫描点云的尺度，尺度很大的情况下比如10^6,与一个小的向量相加会出现结果没变化,所以要乘以尺度
    line_pnt1=boundary_pnts(i,:);%该直线上的第一个点
    line_pnt2=line_pnt1+scale*line_vector(j,:); %该直线上的第二个点
    linepnt=[line_pnt1;line_pnt2]; %直线上的两个点
    [PL_projection] = pntline_projection(linepnt,pnt_inter);%将包络区域范围内的点投影到直线上
    pnt_distrib_index(i,:) = linepnt_distribution(PL_projection) ;%线上点的密度指数，采用中误差及均值表达,越小越好
      k(i,:)=size(pnt_inter,1);
    else
        pnt_distrib_index(i,:)=[inf,inf];
        k(i,:)=0;
    end
end

v1=pnt_distrib_index(:,1);v2=pnt_distrib_index(:,2);
inf_indices1 =~isinf(v1);inf_indices2 =~isinf(v2);
max_index1=max(v1(inf_indices1));max_index2=max(v2(inf_indices2));
     
norm1=[pnt_distrib_index(:,1)/max_index1,pnt_distrib_index(:,2)/max_index2];   %越小越好
norm11=sum(1./exp(norm1),2);
density1_norm=norm11./max(norm11);
density2_norm=k/max(k); %越大越好
density=density1_norm+density2_norm;
density(:,2)=1:n'; %第二列赋予标签号
linedensity= sortrows(density,-1);  %密度从大到小排序，第一列是密度大小排序，第二列是标签号,有可能有很多点密度都是最大的
nn=find(linedensity(:,1)==max(linedensity(:,1)));% 找到相同最大的密度有几个
if size(nn,1)>1 
    linedensity(1,1)=linedensity(1,1)+0.001;
end


for i=1:n
            outlinepnt=boundary_pnts(linedensity(i,2),:);  
            line_pnts=boundary_pnts;
            line_vectors=repmat(line_vector(j,:),n,1);
        if linedensity(i,1)==max(linedensity(:,1))
            D= PLines_distance_TLS(outlinepnt, line_pnts, line_vectors); %单个点到多个直线距离
            max_D=max(D);
            delta(i,:)=max_D;
        else
            index_no=linedensity(1:i-1,2);%前面较小离散的行号
            D= PLines_distance_TLS(outlinepnt, line_pnts(index_no,:), line_vectors(index_no,:));
            min_D=min(D);
            delta(i,:)=min_D;
        end
end
DD2=2*delta/max(delta);
rr=[linedensity(:,2),linedensity(:,1).*DD2];
sort_R{j}=sortrows(rr,-2); %每个向量对应所有轮廓点的密度
end

