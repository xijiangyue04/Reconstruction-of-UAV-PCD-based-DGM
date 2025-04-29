%计算点到空间直线的距离,该空间直线方程为：x=az+b, y=cz+d
%输入的是点x(nx3),空间直线参数parameter(a,b,c,d);
function [PL_dis] = PL_distance_LS(input_pnts,parameter)
n=size(input_pnts,1);
a=parameter(1);
b=parameter(2);
c=parameter(3);
d=parameter(4);
PL_dis=sqrt(((input_pnts(:,2)-d*ones(n,1)-c*input_pnts(:,3)).^2+(a*input_pnts(:,3)+b*ones(n,1)-input_pnts(:,1)).^2+(c*input_pnts(:,1)-a*input_pnts(:,2)+a*d*ones(n,1)-b*c*ones(n,1)).^2)/(a^2+c^2+1));