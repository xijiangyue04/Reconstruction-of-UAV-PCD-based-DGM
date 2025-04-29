%��ÿ����������ȡ��ƽ����ƽ��ж���������ݵ��ж�
%���룺ÿ����������ȡ�����ƽ��·�� path = 'E:\��ά��������\���˻�������������ȡ-׼��Ͷ�Ǹ�drones�ؿ�\2���˻���������ƽ��ָ���\B3ƽ��ָ�\';
%�����Prop ˮƽ��������������ݶ�������Ƶĸ���

function [Prob,horizontal_pnts,vertical_pnts, roof_pnts,nonroof_pnts] = roof_judge(path,Pr) 
namelist = dir([path,'*.txt']);
n= length(namelist);
roof_pnts=[];nonroof_pnts=[];
for k=1:n
filename = [path,namelist(k).name];
%�������д���ȷ�����Ƕ�ȡ��txt�ļ�����1����2��3,4֮���
a=num2str(namelist(k).name);
s=a(isstrprop(a,'digit')); 
kk=str2num(s);
%�����ǶԵ���txt���з���
plane_segment=load(filename);
[parameter] = TLS_Plane(plane_segment(:,1:3));
vector1=parameter(1:3,:)';
vector2=[0,0,1];
dotProduct = dot(vector1, vector2);
angle= acos(dotProduct) * (180 / pi);
  if angle <45   %����봹ֱ��ĽǶȴ���45������Ϊ��ˮƽ��
input=plane_segment(:,1:3);
[pnt_density] = sphere_points_density(input,0.2);%�����ܶ�
pnts_number=size(input,1);%������
mean_z=mean(input(:,3));%zֵ
totalpnt_number(k,:)=[sum(pnt_density)/size(input,1),pnts_number,mean_z,kk];
horizontal_pnts{k}=plane_segment;
  else
      vertical_pnts{k}=plane_segment;
  end
end
horizontal_pnts(cellfun(@isempty,horizontal_pnts))=[];
 vertical_pnts(cellfun(@isempty, vertical_pnts))=[];
totalpnt_number(all(totalpnt_number==0,2),:)=[];
normal1=1/3*(totalpnt_number(:,1)/max(totalpnt_number(:,1)));
normal2=1/3*(totalpnt_number(:,2)/max(totalpnt_number(:,2)));
normal3=1/3*(totalpnt_number(:,3)/max(totalpnt_number(:,3)));
totalpnt_num=[normal1,normal2,normal3];
Prob=[sum(totalpnt_num(:,1:3),2),totalpnt_number(:,4)];

for i=1:size(Prob,1)
    if Prob(i,1)>=Pr
        roof_pnts{i}=horizontal_pnts{i};
    else
        nonroof_pnts{i}=horizontal_pnts{i};
    end
end
if ~isempty(roof_pnts)
roof_pnts(cellfun(@isempty,roof_pnts))=[];
end
if ~isempty(nonroof_pnts)
nonroof_pnts(cellfun(@isempty,nonroof_pnts))=[];
end

