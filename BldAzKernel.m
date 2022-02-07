function [ Kernel ] = BldAzKernel( size, Azimuth )
%BLDAZKERNEL Возвращает ядро с единичными коэффициентами по азмуту перпендикулярному Azimuth
%   Detailed explanation goes here
Azimuth=Azimuth/180*pi; %Перевод азимута в радианы
Kernel=zeros(size);
%Sgma=0.8;
%Sgma=4;
%Sgma=2;
Sgma=1.8;
%Sgma=0.5;
CP=round(size./2);
Col=1:size;
for Row=1:size      %построение ядра с заданым азимутом
R=sqrt((Col-CP).^2+(Row-CP).^2);   
Phi=atan2(Col-CP,Row-CP)-Azimuth;  %Повёрнутые азимуты
x=R(:).*cos(Phi(:));   %Координаты повёрнутой Декартовой СК в точках исходной декартовой СК
Kernel(Row,Col)=1/(Sgma*sqrt(2*pi)).*exp(-(x(:)).^2./(2*Sgma^2)); %Гауссиана 
end
    imagesc(Kernel)
end


