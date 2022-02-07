function [ Data ] = HstgrmAutoRotate( Data ,minBrgtns, maxBrgtns)
%HSTGRMAUTOROTATE ������� ������������� ������������ ����������� ���, �����
%����������� �������� ���� ����������� ��������� (����������) �� �������� ������������� ���������
% Data ��������� ��������� ���, ��� ��� ����������� �������� (�.� ��� �������� �� ������� ������������� ��������� �������� ����������� ���������� � ������ �������)
%minBrgtns - ����������������� ���������� ��������� �������� �������
%maxBrgtns - � ������������, �������������
Hstgrm=BldHstgrm( Data, maxBrgtns,minBrgtns, maxBrgtns);
NPix=numel(Data);  %����� ��������� �������� (��������)
[MVal Ind]=max(Hstgrm); %������������ �������� ����������� � ��� ����� � ��� (������ ����������� ��� ����������������� ���������)
i=1:fix((maxBrgtns+1)./2); %���� �������� ����������� (�� ���������� �� ���������), ��� ��������� ���������� ��������� ����������� � ������� �������
t=i+Ind;    %�������� ��� ��������� ������������ �������� � ������� �������
t(maxBrgtns-Ind+1:end)=t(maxBrgtns-Ind+1:end)-maxBrgtns;    % � ��� ��������� (���������) ��� ����� �� ������� ������������� ���������
WMInd=Ind+sum(Hstgrm(t)'.*i)/NPix;  %��������� �������� � ������� �������
t=Ind-i;    %�������� ��� ��������� � ������� �������
t(Ind:end)=t(Ind:end)+maxBrgtns;    %��������� ������� �����������
WMInd=WMInd-sum(Hstgrm(t)'.*i)/NPix;%��������� �������� � ������� �������
WMInd=round(WMInd);  %���������� ���������� (������) ����������������� ��������� � �����������
if WMInd>maxBrgtns      %�������� � ������ ���������� ������������� ���������
    WMInd=WMInd-(maxBrgtns-minBrgtns);
elseif WMInd<minBrgtns  %�������� � ������ ����� �������� ���� ������ ������� ������������� ���������
    WMInd=WMInd+(maxBrgtns-minBrgtns);
end    
%     plot(Hstgrm);
%     hold on
dlta=WMInd-(maxBrgtns+1)/2; %���������� ������
Data=uint16(mod(double(Data) - dlta, maxBrgtns+1)); %������������� �����������
%     Hstgrm=BldHstgrm( Data, maxBrgtns,minBrgtns, maxBrgtns);
%     plot(Hstgrm);
%     hold off
end

