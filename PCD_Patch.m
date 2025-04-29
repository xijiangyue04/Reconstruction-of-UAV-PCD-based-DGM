

%�����ǽ������ﶥ��ƽ��սǵ㼰���������սǵ��������쵽����Ĺսǵ�
%�����ǶԽ����ﶥ����ƽ�����ȡ  number_of_neighbor=24;line_number=4;dist_thr=0.1,Pr=0.6
%���node_vh ������ÿ����Ĺյ�(Ԫ����ʽ)��node_horizontalÿ���ݶ��յ�(Ԫ����ʽ)
function [node_vh, node_horizontal] = PCD_Patch(number_of_neighbor,line_number,Pr,dist_thr) 

path = 'F:\�������Լ�������\��ͬѧ���Ͽ�\2024-2025��һѧ��\������ά�ؽ�matlab����������֤\data\B9ƽ��ָ�\';  %�������зָ������
% [Prob,horizontal_pnts] = roof_judge(path) ;
[Prob,horizontal_pnts,vertical_pnts, roof_pnts,nonroof_pnts] = roof_judge(path,Pr) ;
% for i=1:size(Prob,1)
%     if Prob(i,1)>=0.6
%         plane_horiz{i}=horizontal_pnts{i};
%     end
% end
% plane_horiz(cellfun(@isempty,plane_horiz))=[];  %��ȡ�˽����ﶥ�����

for i=1:length(roof_pnts)
        [ultimatevector,boundarypnts,clust,cluster_pnts] = horizontalcontour_vector(roof_pnts{i},number_of_neighbor,dist_thr);
        [sort_R,boundary_pnts,line_vector,project_pnt] = horizontalcontour_density2(roof_pnts{i},number_of_neighbor,ultimatevector,dist_thr);
        [horizon_regula,vertical_regula,horizon_insertpnts,vertical_insertpnts,node_pnts] = horizontalcontour_frame_regulariz(line_number,sort_R,boundary_pnts,line_vector,project_pnt,dist_thr); 
        
        node_horizontal{i}=node_pnts;
end
 node_horizontal(cellfun(@isempty,  node_horizontal))=[];
 
 %�����ǽ����������սǵ��������쵽����
 filename=['F:\�������Լ�������\��ͬѧ���Ͽ�\2024-2025��һѧ��\������ά�ؽ�matlab����������֤\data\B9.txt'];
input_pnts=load(filename);
 z=input_pnts(:,3);
sort_z=sort(z);
ground=sort_z(100,:);
 node_vert=[];
for i=1:length(node_horizontal)
    node=node_horizontal{i};
    for j=1:size(node,1)
        if j==size(node,1)
            index(j,:)=[j,1];
            node_add=[node(index(j,:),1:2)];
            node_add(:,3)=ground;
            node_vertical=[node(flip(index(j,:)),:);node_add]; 
        else
           index(j,:)=[j+1,j];
           node_add=[node(index(j,:),1:2)];
           node_add(:,3)=ground;
           node_vertical=[node(flip(index(j,:)),:);node_add]; 
        end
        node_v{j}=node_vertical;
    end
    comb{i}=index;
    index=[];
    node_vert=[node_vert,node_v]; %ÿ���������쵽�����ÿ��ǽ��ڵ�
    node_v=[];
end
node_vh=[node_horizontal,node_vert];
 %���������ݹյ������ά��ģ
     for i=1:length(node_vh)
        vertices=node_vh{i};
        faces=1:1:size(node_vh{i},1);
        patch('Vertices', vertices, 'Faces', faces, ...
      'FaceVertexCData',hsv(1),'FaceColor','[0.5, 0.5, 0.5]','EdgeColor','none'); hold on
     end
view(3);
camlight('left');
