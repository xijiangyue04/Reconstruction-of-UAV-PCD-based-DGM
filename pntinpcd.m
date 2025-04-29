%�жϵ��Ƿ��ڵ����ڲ����жϣ��������ά��Ҫ����ά��ͶӰ��һ��ƽ����
%��Ҫ�жϵĵ�judgepnt(mx3)����Χ�߽����bundary_pnts(nx3)
%���max_diff>90 ���� ���жϸõ��ڵ����ⲿ,����judgepnt�еĵ���bundary_pnts�е�ÿ������������ķ�λ�ǲ�ֵ�����ж�

function [max_diff] = pntinpcd(judgepnt,bundary_pnts)   
n1=size(judgepnt,1);n2=size(bundary_pnts,1);
total=[judgepnt;bundary_pnts];
n=n1+n2;
mean_neighbor=mean(total,1);  
[parameter] = TLS_Plane(total);
normal_pnts =parameter(1:3,:);
lamda=(repmat(mean_neighbor*normal_pnts,n,1)-total*normal_pnts)/(sum(normal_pnts.^2));% kx1�о���
proj_total=normal_pnts*lamda'+total';%3xk
proj_judgepnt=proj_total(:,1:n1);%3xn1
proj_bundary=proj_total(:,n1+1:end);%3xn2
for i=1:n1
    %judgepnt�еĵ���߽����е㹹�������ķ�λ��
   % k=size(bundary_pnts,1);
    
    boundarypnts=proj_bundary;
    object_pnt=proj_judgepnt(:,i)';
    x_unit_vector=(boundarypnts(:,end)'-object_pnt)/norm(boundarypnts(:,end)'-object_pnt);% ͨ��judgepnt�е�һ��Ŀ���object_pnt�ͱ߽��boundarypnts(:,end)'����x�ᷨ���� 1x3 
    y_unit_vector=cross(normal_pnts,x_unit_vector)/norm(cross(normal_pnts,x_unit_vector)); %��������������˳��Բ�˵�ģ��ȷ����һ���������ϵĵ�λ���� 1x3
 
    x_plane=sum(repmat(x_unit_vector,n2-1,1).*(boundarypnts(:,2:end)'-repmat(object_pnt,n2-1,1)),2); % ���ݵ�λ�������������õ���ƽ��x���ϵ�ͶӰkx1
   y_plane=sum(repmat(y_unit_vector,n2-1,1).*(boundarypnts(:,2:end)'-repmat(object_pnt,n2-1,1)),2); % ���ݵ�λ�������������õ���ƽ��y���ϵ�ͶӰkx1
 
   azimuth= rad2deg(atan2(y_plane,x_plane)); %�������귽λ�ǣ�������ڵ�һ�͵ڶ�����������ֵ������ǵ����͵����������Ǹ�ֵ
   azimuth(azimuth<0)= azimuth(azimuth<0)+360 ; %�����귽λ�ǵĸ�ֵȫ����Ϊ��ֵ
   azimuth=[azimuth;360];
   sort_azimuth=sort(azimuth);
   diff_azimuth=diff(sort_azimuth);
   max_diff(i,:)=max(diff_azimuth);
end





% %�жϵ��Ƿ��ڵ����ڲ����жϣ�һ���Ƕ�ά����£����������ά��Ҫ����ά��Ϊ��ά
% %��Ҫ�жϵĵ�input(mx3)����Χ����input_pnts(nx3)
% %���max_diff>180 ���жϸõ��ڵ����ڲ�
% 
% function [max_diff] = pntinpcd(input,input_pnts)   
% for i=1:size(input,1)
%     %���������Ǵ�boundary_extract2 ��������ȡ��Ŀ�������Χ�㹹�������ķ�λ��
%     normal_pnts=[0;0;1];
%     k=size(input_pnts,1);
%     proj_pnts=input_pnts';
%     object_pnt=input(i,:);
%     x_unit_vector=(proj_pnts(:,end)'-object_pnt)/norm(proj_pnts(:,end)'-object_pnt);% ͨ��̽���input_pnts(i,:)��ͶӰ�ڽ�����Զ�ĵ�proj_pnts(:,end)'����x�ᷨ���� 1x3 
%     y_unit_vector=cross(normal_pnts,x_unit_vector)/norm(cross(normal_pnts,x_unit_vector)); %��������������˳��Բ�˵�ģ��ȷ����һ���������ϵĵ�λ���� 1x3
%  
%     x_plane=sum(repmat(x_unit_vector,k-1,1).*(proj_pnts(:,2:end)'-repmat(object_pnt,k-1,1)),2); % ���ݵ�λ�������������õ���ƽ��x���ϵ�ͶӰkx1
%    y_plane=sum(repmat(y_unit_vector,k-1,1).*(proj_pnts(:,2:end)'-repmat(object_pnt,k-1,1)),2); % ���ݵ�λ�������������õ���ƽ��y���ϵ�ͶӰkx1
%  
%    azimuth= rad2deg(atan2(y_plane,x_plane)); %�������귽λ�ǣ�������ڵ�һ�͵ڶ�����������ֵ������ǵ����͵����������Ǹ�ֵ
%    azimuth(azimuth<0)= azimuth(azimuth<0)+360 ; %�����귽λ�ǵĸ�ֵȫ����Ϊ��ֵ
%    sort_azimuth=sort(azimuth);
%    diff_azimuth=diff(sort_azimuth);
%    max_diff(i,:)=max(diff_azimuth);
% end