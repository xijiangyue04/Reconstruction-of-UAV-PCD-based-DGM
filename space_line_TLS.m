

% 计算空间直线方程：(x-xj)/j=(y-yj)/k=(z-zj)/l,(j,k,l)为所求空间直线的单位向量 参考文献：秦锋,王涛,张振虎.空间直线拟合新方法[J].地理空间信息,2023,21(03):21-24.
%[xj,yj,zj]为输入点云的平均值
% line_vector(j,k,l)为方向向量是指导直线沿着何方向延伸的矢量
function [line_vector,mean_pnt] = space_line_TLS(input_pnts)
mean_pnt=mean(input_pnts,1);
input_mean=input_pnts-mean(input_pnts,1);
a11=sum(input_mean(:,2).^2)+sum(input_mean(:,3).^2);
a12=-sum(input_mean(:,1).*input_mean(:,2));
a13=-sum(input_mean(:,1).*input_mean(:,3));
a21=a12;
a22=sum(input_mean(:,1).^2)+sum(input_mean(:,3).^2);
a23=-sum(input_mean(:,2).*input_mean(:,3));
a31=a13;
a32=a23;
a33=sum(input_mean(:,1).^2)+sum(input_mean(:,2).^2);
M=[a11,a12,a13;a21,a22,a23;a31,a32,a33];
[U,S,V] = svd(M); 
line_vector=V(:,3)';

% 以下命令是点云到在该直线上的投影点坐标
% k=size(input_pnts,1);
%     for k=1:k
%     A = input_pnts(k,:);%提取其中一个最远点
%     x0=A(1);
%     y0=A(2);
%     z0=A(3);
%      v = line_vector;
%     PA = [x0; y0; z0] - mean_pnt';
%     proj_v_PA = dot(PA, v) / norm(v)^2 * v;
%     PQ = proj_v_PA;
%     Q = mean_pnt' + PQ;
%     x_Q = Q(1);
%     y_Q = Q(2);
%     z_Q = Q(3);
%     Project_pnt(k,:)=[x_Q,y_Q,z_Q];
%     end
    
    