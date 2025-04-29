%输入线上的点，即直线缓冲区的点，linepnts（nx3）
%输出线性点离散分布情况，采用的是交叉距离标准差和平均值来表示，越小越好

function [pnt_distrib_index] = linepnt_distribution(linepnts) 
if size(linepnts,1)<2
    pnt_distrib_index=[];
    disp('直线上点只有一个，无法得到密度');
else
[endpnts] = line_endpnts(linepnts) ;%直线端点坐标
dis=sqrt(sum((linepnts-endpnts(1,:)).^2,2));% 所有点到任意一个端点距离
combine=[linepnts,dis]; % 
C = sortrows(combine,4);% 对距离进行排序
diff_C=diff(C(:,4));%计算交叉距离
std_C=std(diff_C);%计算交叉距离标准差
mean_C=mean(diff_C);%计算交叉距离平均值
pnt_distrib_index=[std_C,mean_C];
end
%pnt_distrib_index=1/exp(1/std_C)+1/exp(1/mean_C); %因为标准差和平均值的尺度不一样，所以都用指数函数归一到0-1区间