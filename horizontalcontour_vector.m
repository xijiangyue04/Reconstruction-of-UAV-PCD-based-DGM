%第一步 对于分割得到的水平面上的点云
%输入：水平方向的每个分割面plane_segment,邻近点数 number_of_neighbor
%输出：最终的分割平面的边界向量ultimatevector（矩阵形式），边界点boundarypnts（矩阵形式），clust每个边界点所属的类别（矩阵形式） 
%cluster_pnts 每一类的边界点（元包形式）

function [ultimatevector,boundarypnts,clust,cluster_pnts] = horizontalcontour_vector(plane_segment,number_of_neighbor,dist_thr)

[parameter] = TLS_Plane(plane_segment);
[project_pnt] = pntplane_projection(plane_segment,parameter) ;
[boundary_pnts] = boundary_extract2(project_pnt',number_of_neighbor) ;
n=size(boundary_pnts,1);
[neighbors,dist] = knnsearch(boundary_pnts, boundary_pnts, 'k', 8);

%确定每个点领域向量包络的点数
for i=1:n
    inpnt=boundary_pnts(neighbors(i,:),:);%确定每个点的领域点
    previous_inpnt=inpnt;%领域点
    while true  
    [line_vector,mean_pnt] = space_line_TLS(inpnt);
    dist= PL_distance_TLS(boundary_pnts, mean_pnt, line_vector); %点-线距离
    inpnt=boundary_pnts(dist<dist_thr,:);
    if isequal(inpnt, previous_inpnt)
        break
    end
    previous_inpnt = inpnt;
    end
    last_vector(i,:)=line_vector;
    dist2= PL_distance_TLS(boundary_pnts, boundary_pnts(i,:), line_vector); 
    linepnts_number(i,:)=length(find(dist2<dist_thr));
end


%确定每个法向量包络的线性点离散分布情况
for i=1:n
     dist= PL_distance_TLS(boundary_pnts, boundary_pnts(i,:), last_vector(i,:)); 
    pnt_inter=boundary_pnts(dist<0.1,:);%其他点到线的距离范围小于0.1的点数
    if size(pnt_inter,1)>5 %该线范围大于5个点才需要进行密度计算
    scale=sqrt(sum(boundary_pnts(i,:).^2));%利用点的距离来描述扫描点云的尺度，尺度很大的情况下比如10^6,与一个小的向量相加会出现结果没变化,所以要乘以尺度
    line_pnt1=boundary_pnts(i,:);%
    line_pnt2=line_pnt1+scale*last_vector(i,:);
    linepnt=[line_pnt1;line_pnt2];
    [PL_projection] = pntline_projection(linepnt,pnt_inter);
     pnt_distrib_index(i,:) = linepnt_distribution(PL_projection) ;%线上点的密度指数，采用中误差及均值表达,越小越好
    else
        pnt_distrib_index(i,:)=[inf,inf];
    end
end

%将每个法向量包络的点数及包络的线性点离散分布情况 归一化在组合，形成最终的密度
v1=pnt_distrib_index(:,1);v2=pnt_distrib_index(:,2);
inf_indices1 =~isinf(v1);inf_indices2 =~isinf(v2);
max_index1=max(v1(inf_indices1));max_index2=max(v2(inf_indices2));    
norm1=[pnt_distrib_index(:,1)/max_index1,pnt_distrib_index(:,2)/max_index2];   %越小越好
norm11=sum(1./exp(norm1),2);
density1_norm=norm11./max(norm11);
density2_norm=linepnts_number/max(linepnts_number); %越大越好
density=density1_norm+density2_norm;
density(:,2)=1:n'; %第二列赋予标签号
linedensity= sortrows(density,-1); %每个边界点的密度及对应的标签号

%以下是根据密度排序的聚类来对几种主向量进行确定
boundarypnts=boundary_pnts(linedensity(:,2),:);  %依据密度排序后的边界点
lastvector=last_vector(linedensity(:,2),:);%依据密度排序后的每个边界点向量
lastvector(:,4)=0;%将向量最后一列赋值为0
j=1;
clust=zeros(n,1);
while any(lastvector(:,4) == 0) %
[row,~]=find(lastvector(:,4)==0);%确定最后一列为0的向量在哪一行
b=lastvector(row,1:3);
a=repmat(lastvector(row(1),1:3),size(b,1),1); 
d_degree=Pnts_normal_angle(a,b);% 第一个密度高的向量与其他向量夹角
row2=row(d_degree<11);%满足向量夹角小于11°的在哪些行
clust(row2,1)=j;%将满足条件的这些行赋予标签号
cluster_num(j,:)=length(row2);%每个满足条件的类别的数量
ultimate_vector(j,:)=b(1,:);%每个满足条件的类别的主向量
line_pnt(j,:)=boundarypnts(row(1),:);%每个满足条件的类别的主点坐标
cluster_pnts{j}=boundarypnts(row2,:);%每个满足条件的所有点坐标 
lastvector(row2,4)=1;  %将满足条件的向量赋予1
j=j+1;
end

gradient=abs(diff(cluster_num)./cluster_num(1:end-1,:));
gradient(length(cluster_num))=0;
for i=1:length(cluster_num)
    ratio=sum(cluster_num(1:i))/n;   
    if ratio>0.7 | gradient(i)>0.7  %如果前面几个点数大于总点数70%或变化率超过0.7
        ultimatevector=ultimate_vector(1:i,:);
        linepnt=line_pnt(1:i,:);
        break
    end
end

if size(ultimatevector,1)<2
    ultimatevector=ultimate_vector(1:2,:); %将其最终主向量数量定死了，如果主向量只有1个，最终将其定为2个，
    linepnt=line_pnt(1:2,:);
end

if size(ultimatevector,1)>3 %如果主向量过多，超过3个了，最终将其定为3个
    ultimatevector=ultimate_vector(1:3,:);
    linepnt=line_pnt(1:3,:);
end







% for i=1:size(cluster_num1,1)
%             outlinepnt=line_pnts(i,:);  
%         if cluster_num1(i,1)==max(cluster_num1(:,1))
%             D= PLines_distance_TLS(outlinepnt, line_pnts, vector); %单个点到多个直线距离
%             max_D=max(D);
%             delta(i,:)=max_D;
%         else
%             index_no=linedensity(1:i-1,2);%前面较小离散的行号
%             D= PLines_distance_TLS(outlinepnt, line_pnts(1:i-1,:), vector);
%             min_D=min(D);
%             delta(i,:)=min_D;
%         end
% end


% 
%     



% for i=1:n
% [vectors_onplane] = onplane_vector(boundary_pnts(i,:),1,1);
% [line_vector,mean_pnt] = space_line_TLS( boundary_pnts(neighbors(i,:),:)); 
% nn=size(vectors_onplane,1);
% for j=1:nn
% dist= PL_distance_TLS(boundary_pnts, boundary_pnts(i,:), vectors_onplane(j,:)); %点-线距离
% outpnt=boundary_pnts(dist<0.1,:);
% kk(j,:)=size(outpnt,1);
% end
% max_kk(i,:)=max(kk);
% index=find(kk==max(kk));
% original_vector(i,:)=vectors_onplane(index(1),:);
% end


%对几种向量进行确定
% cluster_number=1;
% while true 
%     cluster_number=cluster_number+1;
% [clust,center] = kmeans(last_vector,cluster_number,'dist','sqeuclidean');
% s=silhouette(last_vector,clust);
% mean_s=mean(s);
% if mean_s>0.75
%     break
% elseif length(find(clust==cluster_number))/n<0.15
%     break
% end
% end







    
    




