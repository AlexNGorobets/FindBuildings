function [ Kernel ] = Bld3DAzKernel( Size, Azimuth )
%BLDAZKERNEL Возвращает ядро с единичными коэффициентами по азмуту перпендикулярному Azimuth
%   Detailed explanation goes here
%Azimuth=Azimuth/180*pi; %Перевод азимута в радианы
Kernel=zeros(Size,Size,max(size(Azimuth)));
%Sgma=0.4;
%Sgma=0.35;
Sgma=0.4;
%SgmaVrsusR_Coef=0.1;
%SgmaVrsusR_Coef=0.12;
SgmaVrsusR_Coef=0.08;
CP=ceil(Size./2);
for Col=1:Size;
  for Row=1:Size      %построение ядра с заданым азимутом
    R=sqrt((Col-CP).^2+(Row-CP).^2);   
    Phi=atan2(Row-CP,Col-CP)-Azimuth;  %Повёрнутые азимуты
    x=R.*cos(Phi(:));   %Вариант для полярной СК с ростом фи по чаосвой стрелке (нулевой азимут - линия максимумов вертикальна)
    %Kernel(Row,Col,:)=cos(R./(Size/2))*(1+R*SgmaVrsusR_Coef)/((Sgma+R*SgmaVrsusR_Coef)*sqrt(2*pi)).*exp(-(x(:)).^2./(2*(Sgma+R*SgmaVrsusR_Coef)^2)); %Гауссиана 
    Kernel(Row,Col,:)=cos(R./(Size/2))*exp(-(x(:)).^2./(2*(Sgma+R*SgmaVrsusR_Coef)^2)); %Гауссиана 
  end
end
    %imagesc(Kernel(:,:,1))
end


