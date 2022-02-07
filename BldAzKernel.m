function [ Kernel ] = BldAzKernel( size, Azimuth )
%BLDAZKERNEL ���������� ���� � ���������� �������������� �� ������ ����������������� Azimuth
%   Detailed explanation goes here
Azimuth=Azimuth/180*pi; %������� ������� � �������
Kernel=zeros(size);
%Sgma=0.8;
%Sgma=4;
%Sgma=2;
Sgma=1.8;
%Sgma=0.5;
CP=round(size./2);
Col=1:size;
for Row=1:size      %���������� ���� � ������� ��������
R=sqrt((Col-CP).^2+(Row-CP).^2);   
Phi=atan2(Col-CP,Row-CP)-Azimuth;  %��������� �������
x=R(:).*cos(Phi(:));   %���������� ��������� ���������� �� � ������ �������� ���������� ��
Kernel(Row,Col)=1/(Sgma*sqrt(2*pi)).*exp(-(x(:)).^2./(2*Sgma^2)); %��������� 
end
    imagesc(Kernel)
end


