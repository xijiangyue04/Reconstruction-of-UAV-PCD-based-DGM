% 计算空间直线方程：x=az+b, y=cz+d
% 计算空间直线参数：a,b,c,d,%输入的是空间直线的点x(nx3)
%该空间直线的方向向量为vector=[a,c,1],方向向量是指导直线沿着何方向延伸的矢量
function [parameter] = space_line_LS(input_pnts)
n=size(input_pnts,1);
one=ones(n,1);
F=[input_pnts(:,3)';one'];
X=input_pnts(:,1);
Y=input_pnts(:,2);
A=pinv(F*F')*F*X;
B=pinv(F*F')*F*Y;
a=A(1,1);
b=A(2,1);
c=B(1,1);
d=B(2,1);
parameter=[a,b,c,d];