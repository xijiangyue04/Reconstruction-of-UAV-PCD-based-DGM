
% 已知两个点，对这两个点之间的直线插入一些点
%输入这两个点坐标start_pnt,end_pnt,及插入的点之间距离resolution
function [insert_pnts] = interpolation_pnts(start_pnt,end_pnt,resolution) 
    fps_dis=sqrt(sum((end_pnt-start_pnt).^2,2)); %投影后最远点距离
    insert_number=fix(fps_dis/resolution); %0.02m插入一个点,总共插入多少个点
    if insert_number==0
        insert_number=1;
    end
    fps_interval=end_pnt-start_pnt;%最远两个点间隔
    neighbor_interval=fps_interval/insert_number;%插入点相邻两个点间隔
    for num=1:insert_number %之所以这样是希望向前和向后多插入点
        insert_pnts(num,:)=start_pnt+neighbor_interval*(num-2);%插入的点坐标
    end