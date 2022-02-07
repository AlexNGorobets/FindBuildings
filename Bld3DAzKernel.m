function [ Kernel ] = Bld3DAzKernel( Size, Azimuth )
%BLDAZKERNEL ���������� ���� � ���������� �������������� �� ������ ����������������� Azimuth
%   Detailed explanation goes here
%Azimuth=Azimuth/180*pi; %������� ������� � �������
Kernel=zeros(Size,Size,max(size(Azimuth)));
%Sgma=0.4;
%Sgma=0.35;
Sgma=0.4;
%SgmaVrsusR_Coef=0.1;
%SgmaVrsusR_Coef=0.12;
SgmaVrsusR_Coef=0.08;
CP=ceil(Size./2);
for Col=1:Size;
  for Row=1:Size      %���������� ���� � ������� ��������
    R=sqrt((Col-CP).^2+(Row-CP).^2);   
    Phi=atan2(Row-CP,Col-CP)-Azimuth;  %��������� �������
    x=R.*cos(Phi(:));   %������� ��� �������� �� � ������ �� �� ������� ������� (������� ������ - ����� ���������� �����������)
    %Kernel(Row,Col,:)=cos(R./(Size/2))*(1+R*SgmaVrsusR_Coef)/((Sgma+R*SgmaVrsusR_Coef)*sqrt(2*pi)).*exp(-(x(:)).^2./(2*(Sgma+R*SgmaVrsusR_Coef)^2)); %��������� 
    Kernel(Row,Col,:)=cos(R./(Size/2))*exp(-(x(:)).^2./(2*(Sgma+R*SgmaVrsusR_Coef)^2)); %��������� 
  end
end
    %imagesc(Kernel(:,:,1))
end


