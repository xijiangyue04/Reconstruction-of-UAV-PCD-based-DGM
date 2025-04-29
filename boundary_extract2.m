%M-file ,boundary_extract.m   ������������Ͻ��ģ�Ŀ������ڽ��������ı��
%����ƽ����ͶӰ���귽λ���жϱ߽�㣨���Ч�������Ͻ�Ҳ����õģ� �ڽ���number_of_neighbor��Խ��max_diff��Ӧ����ֵҲ��������΢С��
%input_pnts(nx3)
function [boundary_pnts] = boundary_extract2(input_pnts,number_of_neighbor)  %���������������ȷ�ģ�Ҳ��Ч������õģ�
n=size(input_pnts,1);     %��õ�һά�ȣ��У��Ĺ�ģ
k=number_of_neighbor;
neighbor_idx=knnsearch(input_pnts,input_pnts,'k',number_of_neighbor+1);  %��P���Ӧ��k+1���ڽ�������һ��
normal_pnts=zeros(3,n);    %����3*n��ȫ�����
line=zeros(n,1);
for i=1:n
    
  neighbor_pnts=input_pnts(neighbor_idx(i,:),:); 
  mean_neighbor=mean(neighbor_pnts,1);  
[parameter] = TLS_Plane(neighbor_pnts);
  normal_pnts =parameter(1:3,:);

 lamda=(repmat(mean_neighbor*normal_pnts,number_of_neighbor+1,1)-neighbor_pnts*normal_pnts)/(sum(normal_pnts.^2));% kx1�о���
 proj_pnts=normal_pnts*lamda'+neighbor_pnts';%3xk
 object_pnt=proj_pnts(:,1)';
 
 x_unit_vector=(proj_pnts(:,end)'-object_pnt)/norm(proj_pnts(:,end)'-object_pnt);% ͨ��̽���input_pnts(i,:)��ͶӰ�ڽ�����Զ�ĵ�proj_pnts(:,end)'����x�ᷨ���� 1x3 
 y_unit_vector=cross(normal_pnts,x_unit_vector)/norm(cross(normal_pnts,x_unit_vector)); %��������������˳��Բ�˵�ģ��ȷ����һ���������ϵĵ�λ���� 1x3
 
 x_plane=sum(repmat(x_unit_vector,k,1).*(proj_pnts(:,2:end)'-repmat(object_pnt,k,1)),2); % ���ݵ�λ�������������õ���ƽ��x���ϵ�ͶӰkx1
 y_plane=sum(repmat(y_unit_vector,k,1).*(proj_pnts(:,2:end)'-repmat(object_pnt,k,1)),2); % ���ݵ�λ�������������õ���ƽ��y���ϵ�ͶӰkx1
 
 azimuth= rad2deg(atan2(y_plane,x_plane)); %�������귽λ�ǣ�������ڵ�һ�͵ڶ�����������ֵ������ǵ����͵����������Ǹ�ֵ
 azimuth(azimuth<0)= azimuth(azimuth<0)+360 ; %�����귽λ�ǵĸ�ֵȫ����Ϊ��ֵ
sort_azimuth=sort(azimuth);
diff_azimuth=diff(sort_azimuth);
max_diff(i,:)=max(diff_azimuth);
end
idx=find(max_diff>110);
boundary_pnts=input_pnts(idx,:);