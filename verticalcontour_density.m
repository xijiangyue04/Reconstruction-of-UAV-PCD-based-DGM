%��һ�������룺ǽ��ָ��ƽ�����plane_segment��number_of_neighbor ���ڽ�����
%ȷ��ÿ���ָ�ƽ������߽��ĺ���������ܶ����ģ�������Զ�Ķ�λ��
%˵���˾���Ҫȷ���߽�������������ܶ�λ�㼰��λ��������
%�����sort_R����ǽ��������������е��ܶȣ�boundary_pnts ǽ����������㣬line_vector��ǽ��������������
function [sort_R,boundary_pnts,line_vector,project_pnt] = verticalcontour_density(plane_segment,number_of_neighbor,dist_thr)


[parameter] = TLS_Plane(plane_segment);
[project_pnt] = pntplane_projection(plane_segment,parameter) ;
project_pnt=project_pnt';
[boundary_pnts] = boundary_extract2(project_pnt,number_of_neighbor) ;
vect1=parameter(1:3,:)';%�ָ���ķ�����
vect2=[0,0,1];%���ǻ�׼����
line_vector1=cross(vect1,vect2);%���׼�����ͷָ��������Ĵ�ֱ����
line_vector2=cross(line_vector1,vect1);
line_vector=[line_vector1;line_vector2];
n=size(boundary_pnts,1);

for j=1:size(line_vector,1)
     
for i=1:n
    dist= PL_distance_TLS(boundary_pnts, boundary_pnts(i,:), line_vector(j,:)); %��-�߾���
    pnt_inter=boundary_pnts(dist<dist_thr,:);%�����㵽�ߵľ��뷶ΧС��0.1�ĵ���
    if size(pnt_inter,1)>5 %���߷�Χ����5�������Ҫ�����ܶȼ���
    scale=sqrt(sum(boundary_pnts(i,:).^2));%���õ�ľ���������ɨ����Ƶĳ߶ȣ��߶Ⱥܴ������±���10^6,��һ��С��������ӻ���ֽ��û�仯,����Ҫ���Գ߶�
    line_pnt1=boundary_pnts(i,:);%
    line_pnt2=line_pnt1+scale*line_vector(j,:);
    linepnt=[line_pnt1;line_pnt2];
    [PL_projection] = pntline_projection(linepnt,pnt_inter);
    pnt_distrib_index(i,:) = linepnt_distribution(PL_projection) ;%���ϵ���ܶ�ָ��������������ֵ���,ԽСԽ��
      k(i,:)=size(pnt_inter,1);
    else
        pnt_distrib_index(i,:)=[inf,inf];
        k(i,:)=0;
    end
end

v1=pnt_distrib_index(:,1);v2=pnt_distrib_index(:,2);
inf_indices1 =~isinf(v1);inf_indices2 =~isinf(v2);
max_index1=max(v1(inf_indices1));max_index2=max(v2(inf_indices2));
     
norm1=[pnt_distrib_index(:,1)/max_index1,pnt_distrib_index(:,2)/max_index2];   %ԽСԽ��
norm11=sum(1./exp(norm1),2);
density1_norm=norm11./max(norm11);
density2_norm=k/max(k); %Խ��Խ��
density=density1_norm+density2_norm;
density(:,2)=1:n'; %�ڶ��и����ǩ��
linedensity= sortrows(density,-1);  %�ܶȴӴ�С���򣬵�һ�����ܶȴ�С���򣬵ڶ����Ǳ�ǩ��
 
for i=1:n
            outlinepnt=boundary_pnts(linedensity(i,2),:);  
            line_pnts=boundary_pnts;
            line_vectors=repmat(line_vector(j,:),n,1);
        if linedensity(i,1)==max(linedensity(:,1))
            D= PLines_distance_TLS(outlinepnt, line_pnts, line_vectors); %�����㵽���ֱ�߾���
            max_D=max(D);
            delta(i,:)=max_D;
        else
            index_no=linedensity(1:i-1,2);%ǰ���С��ɢ���к�
            D= PLines_distance_TLS(outlinepnt, line_pnts(index_no,:), line_vectors(index_no,:));
            min_D=min(D);
            delta(i,:)=min_D;
        end
end
DD2=2*delta/max(delta);
rr=[linedensity(:,2),linedensity(:,1).*DD2];
sort_R{j}=sortrows(rr,-2);
end

