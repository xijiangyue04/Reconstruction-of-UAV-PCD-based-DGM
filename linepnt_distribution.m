%�������ϵĵ㣬��ֱ�߻������ĵ㣬linepnts��nx3��
%������Ե���ɢ�ֲ���������õ��ǽ�������׼���ƽ��ֵ����ʾ��ԽСԽ��

function [pnt_distrib_index] = linepnt_distribution(linepnts) 
if size(linepnts,1)<2
    pnt_distrib_index=[];
    disp('ֱ���ϵ�ֻ��һ�����޷��õ��ܶ�');
else
[endpnts] = line_endpnts(linepnts) ;%ֱ�߶˵�����
dis=sqrt(sum((linepnts-endpnts(1,:)).^2,2));% ���е㵽����һ���˵����
combine=[linepnts,dis]; % 
C = sortrows(combine,4);% �Ծ����������
diff_C=diff(C(:,4));%���㽻�����
std_C=std(diff_C);%���㽻������׼��
mean_C=mean(diff_C);%���㽻�����ƽ��ֵ
pnt_distrib_index=[std_C,mean_C];
end
%pnt_distrib_index=1/exp(1/std_C)+1/exp(1/mean_C); %��Ϊ��׼���ƽ��ֵ�ĳ߶Ȳ�һ�������Զ���ָ��������һ��0-1����