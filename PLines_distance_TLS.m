%�����㵽���ֱ�ߵľ���,outlinepntΪ������(1x3)��line_pntsΪ���ֱ���ϵ�(nx3)��line_vectors���ֱ������(nx3)
%���distance��nx1��
function [distance] = PLines_distance_TLS(outlinepnt, line_pnts, line_vectors)
    n=size(line_pnts,1);
    if size(line_vectors,1)==1
        line_vectors=repmat(line_vectors,n,1);
    end
     dxyz=outlinepnt- line_pnts;% (nx3),
     v_length =sqrt(sum(line_vectors.^2,2)); %(nx1),
     Vector_cross=cross(dxyz',line_vectors');%(3xn),
     cross_product_length=sqrt(Vector_cross(1,:).^2+Vector_cross(2,:).^2+Vector_cross(3,:).^2);%(1xn),
     distance = cross_product_length'./ v_length; %(nx1),
end