%��һ�� ���ڷָ�õ���ˮƽ���ϵĵ���
%���룺ˮƽ�����ÿ���ָ���plane_segment,�ڽ����� number_of_neighbor
%��������յķָ�ƽ��ı߽�����ultimatevector��������ʽ�����߽��boundarypnts��������ʽ����clustÿ���߽����������𣨾�����ʽ�� 
%cluster_pnts ÿһ��ı߽�㣨Ԫ����ʽ��

function [ultimatevector,boundarypnts,clust,cluster_pnts] = horizontalcontour_vector(plane_segment,number_of_neighbor,dist_thr)

[parameter] = TLS_Plane(plane_segment);
[project_pnt] = pntplane_projection(plane_segment,parameter) ;
[boundary_pnts] = boundary_extract2(project_pnt',number_of_neighbor) ;
n=size(boundary_pnts,1);
[neighbors,dist] = knnsearch(boundary_pnts, boundary_pnts, 'k', 8);

%ȷ��ÿ����������������ĵ���
for i=1:n
    inpnt=boundary_pnts(neighbors(i,:),:);%ȷ��ÿ����������
    previous_inpnt=inpnt;%�����
    while true  
    [line_vector,mean_pnt] = space_line_TLS(inpnt);
    dist= PL_distance_TLS(boundary_pnts, mean_pnt, line_vector); %��-�߾���
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


%ȷ��ÿ����������������Ե���ɢ�ֲ����
for i=1:n
     dist= PL_distance_TLS(boundary_pnts, boundary_pnts(i,:), last_vector(i,:)); 
    pnt_inter=boundary_pnts(dist<0.1,:);%�����㵽�ߵľ��뷶ΧС��0.1�ĵ���
    if size(pnt_inter,1)>5 %���߷�Χ����5�������Ҫ�����ܶȼ���
    scale=sqrt(sum(boundary_pnts(i,:).^2));%���õ�ľ���������ɨ����Ƶĳ߶ȣ��߶Ⱥܴ������±���10^6,��һ��С��������ӻ���ֽ��û�仯,����Ҫ���Գ߶�
    line_pnt1=boundary_pnts(i,:);%
    line_pnt2=line_pnt1+scale*last_vector(i,:);
    linepnt=[line_pnt1;line_pnt2];
    [PL_projection] = pntline_projection(linepnt,pnt_inter);
     pnt_distrib_index(i,:) = linepnt_distribution(PL_projection) ;%���ϵ���ܶ�ָ��������������ֵ���,ԽСԽ��
    else
        pnt_distrib_index(i,:)=[inf,inf];
    end
end

%��ÿ������������ĵ�������������Ե���ɢ�ֲ���� ��һ������ϣ��γ����յ��ܶ�
v1=pnt_distrib_index(:,1);v2=pnt_distrib_index(:,2);
inf_indices1 =~isinf(v1);inf_indices2 =~isinf(v2);
max_index1=max(v1(inf_indices1));max_index2=max(v2(inf_indices2));    
norm1=[pnt_distrib_index(:,1)/max_index1,pnt_distrib_index(:,2)/max_index2];   %ԽСԽ��
norm11=sum(1./exp(norm1),2);
density1_norm=norm11./max(norm11);
density2_norm=linepnts_number/max(linepnts_number); %Խ��Խ��
density=density1_norm+density2_norm;
density(:,2)=1:n'; %�ڶ��и����ǩ��
linedensity= sortrows(density,-1); %ÿ���߽����ܶȼ���Ӧ�ı�ǩ��

%�����Ǹ����ܶ�����ľ������Լ�������������ȷ��
boundarypnts=boundary_pnts(linedensity(:,2),:);  %�����ܶ������ı߽��
lastvector=last_vector(linedensity(:,2),:);%�����ܶ�������ÿ���߽������
lastvector(:,4)=0;%���������һ�и�ֵΪ0
j=1;
clust=zeros(n,1);
while any(lastvector(:,4) == 0) %
[row,~]=find(lastvector(:,4)==0);%ȷ�����һ��Ϊ0����������һ��
b=lastvector(row,1:3);
a=repmat(lastvector(row(1),1:3),size(b,1),1); 
d_degree=Pnts_normal_angle(a,b);% ��һ���ܶȸߵ����������������н�
row2=row(d_degree<11);%���������н�С��11�������Щ��
clust(row2,1)=j;%��������������Щ�и����ǩ��
cluster_num(j,:)=length(row2);%ÿ��������������������
ultimate_vector(j,:)=b(1,:);%ÿ����������������������
line_pnt(j,:)=boundarypnts(row(1),:);%ÿ������������������������
cluster_pnts{j}=boundarypnts(row2,:);%ÿ���������������е����� 
lastvector(row2,4)=1;  %��������������������1
j=j+1;
end

gradient=abs(diff(cluster_num)./cluster_num(1:end-1,:));
gradient(length(cluster_num))=0;
for i=1:length(cluster_num)
    ratio=sum(cluster_num(1:i))/n;   
    if ratio>0.7 | gradient(i)>0.7  %���ǰ�漸�����������ܵ���70%��仯�ʳ���0.7
        ultimatevector=ultimate_vector(1:i,:);
        linepnt=line_pnt(1:i,:);
        break
    end
end

if size(ultimatevector,1)<2
    ultimatevector=ultimate_vector(1:2,:); %�����������������������ˣ����������ֻ��1�������ս��䶨Ϊ2����
    linepnt=line_pnt(1:2,:);
end

if size(ultimatevector,1)>3 %������������࣬����3���ˣ����ս��䶨Ϊ3��
    ultimatevector=ultimate_vector(1:3,:);
    linepnt=line_pnt(1:3,:);
end







% for i=1:size(cluster_num1,1)
%             outlinepnt=line_pnts(i,:);  
%         if cluster_num1(i,1)==max(cluster_num1(:,1))
%             D= PLines_distance_TLS(outlinepnt, line_pnts, vector); %�����㵽���ֱ�߾���
%             max_D=max(D);
%             delta(i,:)=max_D;
%         else
%             index_no=linedensity(1:i-1,2);%ǰ���С��ɢ���к�
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
% dist= PL_distance_TLS(boundary_pnts, boundary_pnts(i,:), vectors_onplane(j,:)); %��-�߾���
% outpnt=boundary_pnts(dist<0.1,:);
% kk(j,:)=size(outpnt,1);
% end
% max_kk(i,:)=max(kk);
% index=find(kk==max(kk));
% original_vector(i,:)=vectors_onplane(index(1),:);
% end


%�Լ�����������ȷ��
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







    
    




