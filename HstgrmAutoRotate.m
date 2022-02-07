function [ Data ] = HstgrmAutoRotate( Data ,minBrgtns, maxBrgtns)
%HSTGRMAUTOROTATE Функция автоматически обрабатывает изображение так, чтобы
%гистограмма яркостей была установлена серединой (максимумом) на середине динамического диапазона
% Data считается искажённой так, что его гистограмма замкнута (т.е при смещении за пределы динамического диапазона значения гистограммы появляются с другой стороны)
%minBrgtns - Предустановленное минимально возможное значение яркости
%maxBrgtns - и максимальное, сооветственно
Hstgrm=BldHstgrm( Data, maxBrgtns,minBrgtns, maxBrgtns);
NPix=numel(Data);  %Общее количесво отсчётов (пикселей)
[MVal Ind]=max(Hstgrm); %Максимальное значение гистограммы и его адрес в ней (первое приближение для средневзвешенного максимума)
i=1:fix((maxBrgtns+1)./2); %Веса значений гистограммы (их расстояния от максимума), для коррекции координаты максимума гистограммы в БОльшую сторону
t=i+Ind;    %Итератор для коррекции центрального значения в бОльшую сторону
t(maxBrgtns-Ind+1:end)=t(maxBrgtns-Ind+1:end)-maxBrgtns;    % и его коррекция (замыкание) при уходе за пределы динамического диапазона
WMInd=Ind+sum(Hstgrm(t)'.*i)/NPix;  %Коррекция значения в бОльшую сторону
t=Ind-i;    %Итератор для коррекции в меньшую сторону
t(Ind:end)=t(Ind:end)+maxBrgtns;    %Замыкание адресов гистограммы
WMInd=WMInd-sum(Hstgrm(t)'.*i)/NPix;%Коррекция значения в меньшую сторону
WMInd=round(WMInd);  %Округление координаты (адреса) средневзвешенного максимума в гистограмме
if WMInd>maxBrgtns      %Поправка в случае привышения динамического диапазона
    WMInd=WMInd-(maxBrgtns-minBrgtns);
elseif WMInd<minBrgtns  %Поправка в случае когда максимум ниже нижней границы динамического диапазона
    WMInd=WMInd+(maxBrgtns-minBrgtns);
end    
%     plot(Hstgrm);
%     hold on
dlta=WMInd-(maxBrgtns+1)/2; %Вычисление сдвига
Data=uint16(mod(double(Data) - dlta, maxBrgtns+1)); %Центрирование гистограммы
%     Hstgrm=BldHstgrm( Data, maxBrgtns,minBrgtns, maxBrgtns);
%     plot(Hstgrm);
%     hold off
end

