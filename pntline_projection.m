
%输入变量直线上的点linepnt(nx3),至少两个点，直线外的点outpnt(mx3)

function [PL_projection] = pntline_projection(linepnt,outpnt)    %

% 定义直线上的两个点 A 和 B
A = linepnt(1,:);
B = linepnt(2,:);

% 定义要投影的点 P
P = outpnt;
n=size(P,1);
% 计算直线 AB 的方向向量
AB = B - A;


% 计算向量 AP
AP = P - A;

% 计算点 P 在直线 AB 上的投影点
PL_projection = A + transpose(dot(AP', repmat(AB,n,1)') / dot(AB, AB)) .* AB;


