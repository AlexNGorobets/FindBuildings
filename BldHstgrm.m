function [ Hstgrm ] = BldHstgrm( Data, Nlvl ,LowLimit, HighLimit)
%BLDHSTGRM ��������� ����������� �������� (�������) �� ������� Data
%Nlvl ��������� ������� � �����������
Hstgrm=zeros(Nlvl,1);
i=1:Nlvl;   %�������� ������� �����������
Data=Data(:);%�������������� ������� ������ � ����������
if nargin==2
    LowLimit=min(Data');     %������ ������ ������� (�������)
    HighLimit=max(Data');    %������� ������
end
Step=(HighLimit-LowLimit)./Nlvl;%��� ����������� (�������) � �����������
%HstgrmLimits(i+1)=LowLimit+i.*Step; %������� (�������) ��� ������ ���������� (�������� �������� - �������)
HstgrmLimits(i)=LowLimit+i.*Step; %������� ������� (�������) ��� ������� �� ����������
Data=sort(Data);    %��������������� ����������
t=1;    %�������� ��������� ������ (��������)
i=1;    %�������� ������� �����������
while t<=size(Data,1);  %���� ���������� ������� �������� �� �������
    if (Data(t)>HstgrmLimits(i))    %������� ��������� ������ �����������
        i=find(HstgrmLimits>=Data(t),1);    %����� ����������� ��������� � ����������� (�� ������� �������)
    end
    Hstgrm(i)=Hstgrm(i)+1;  %��������� �������� �����������
    t=t+1;  %��������� �������� ������
end
% for i=1:Nlvl    %���� ���������� �����������
%     BoolMap=Data>=HstgrmLimits(i);  %������� �������� ������� �������� ������
%     BoolMap=BoolMap.*(Data<HstgrmLimits(i+1));%������� ��������� �������
%     Hstgrm(i)=sum(BoolMap);
% end

end

