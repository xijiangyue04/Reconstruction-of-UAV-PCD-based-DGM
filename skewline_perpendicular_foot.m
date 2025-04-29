    % 输入参数:
    % P1, P2: 两条直线上的点, 1x3 矩阵
    % v1, v2: 两条直线的方向向量, 1x3 矩阵
    % 输出参数:
    % intersection: 交点坐标, 2x3 矩阵
    %参考文献：[1]董玉久,潘秀英,杨欣欣.两异面直线公垂线垂足位置的计算方法[J].哈尔滨科学技术大学学报,1990,14(01):84-88.

function [intersection] = skewline_perpendicular_foot(P1,v1,P2,v2)
x1=P1(1);y1=P1(2);z1=P1(3);
x2=P2(1);y2=P2(2);z2=P2(3);
   m1=v1(1);n1=v1(2);p1=v1(3);
   m2=v2(1);n2=v2(2);p2=v2(3);
   A1=det([n1,p1;n2,p2]);
   B1=det([p1,m1;p2,m2]);
   C1=det([m1,n1;m2,n2]);
   A2=det([n2,B1;p2,C1]);
   B2=det([p2,C1;m2,A1]);
   C2=det([m2,A1;n2,B1]);
   delta1=det([A1,B1,C1;A2,B2,C2;n1,-m1,0]);
   D1=A2*(x2-x1)+B2*(y2-y1)+C2*(z2-z1);
   XG=x1-D1*m1*C1/delta1;
   YG=y1-D1*n1*C1/delta1;
   ZG=z1+D1*(A1*m1+B1*n1)/delta1;
   A3=n1*det([m1,n1;m2,n2])+p1*det([m1,p1;m2,p2]);
   B3=p1*det([n1,p1;n2,p2])+m1*det([n1,m1;n2,m2]);
   C3=m1*det([p1,m1;p2,m2])+n1*det([p1,n1;p2,n2]);
   delta2=n2*det([B1,C1;B3,C3])+m2*det([A1,C1;A3,C3]);
   D2=A3*(x1-x2)+B3*(y1-y2)+C3*(z1-z2);
   XH=x2-D2*m2*C1/delta2;
   YH=y2-D2*n2*C1/delta2;
   ZH=z2+D2*(A1*m2+B1*n2)/delta2;
   intersection=[XG,YG,ZG;XH,YH,ZH];
end



