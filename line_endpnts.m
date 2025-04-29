%直线型点云数据确定两头端点坐标,输入直线型点云数据input_pnts(nx3)
%输出endpnts(2x3)两头端点坐标
function [endpnts] = line_endpnts(input_pnts) 
         pnts=unique(input_pnts,'rows');
if size(input_pnts,1)<3
    endpnts=input_pnts;
else
     mean_pnt=mean(pnts);
     PCD=sqrt(sum((pnts-mean_pnt).^2,2));%每个直线点到中心点距离
     PCD=round(PCD,3);%四舍五入保留小数点三位
     pnt_mean_vector=pnts-mean_pnt;
     norm_vector=sqrt(sum(pnt_mean_vector.^2,2));
     max_norm=find(norm_vector==max(norm_vector));
     nn=size(pnts,1);
     vector1=repmat(pnt_mean_vector(max_norm(1),:),nn,1);%
     angle=dot(vector1',pnt_mean_vector')'./(sqrt(sum(vector1.^2,2)).*sqrt(sum(pnt_mean_vector.^2,2)));
     combine=[PCD,angle,(1:nn)'];
     combine_sort=sortrows(combine,1,'descend');
     positive_row=combine_sort(combine_sort(:,2)>0,3);    
     negative_row=combine_sort(combine_sort(:,2)<0,3);   
     endpoint_row=[positive_row(1,:);negative_row(1,:)];
     endpnts=[pnts(endpoint_row(1,:),:);pnts(endpoint_row(2,:),:)];
end
