%计算点到空间直线的距离,该空间直线方程为(x-xj)/j=(y-yj)/k=(z-zj)/l,(j;k;l)为所求空间直线的单位向量
%输入的是点云input_pnts(nx3),空间直线方向向量line_vector(j,k,l)，空间直线上一点 line_point=(xj,yj,zj)
%输出distance（1xn）
function [distance] = PL_distance_TLS(input_pnts, line_point, line_vector)
    n=size(input_pnts,1);
    dxyz=input_pnts-line_point;% (nx3),
    v_length = norm(line_vector); %(1x1),
    Vector_cross=cross(dxyz',repmat(line_vector,n,1)');%(3xn),
     cross_product_length=sqrt(Vector_cross(1,:).^2+Vector_cross(2,:).^2+Vector_cross(3,:).^2);%(1xn),
     distance = cross_product_length / v_length; %(1xn),
end