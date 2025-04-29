%用于计算点云距离分辨率，即确定相邻点之间的距离
function [range_resol] = range_resolut(input_pnts)
[input_pnts, ia, ic] = unique(input_pnts, 'rows');
neighbor_idx=knnsearch(input_pnts,input_pnts,'k',2);
dis=sqrt(sum((input_pnts-input_pnts(neighbor_idx(:,2),:)).^2,2));
range_resol=mean(dis);