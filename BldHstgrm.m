function [ Hstgrm ] = BldHstgrm( Data, Nlvl ,LowLimit, HighLimit)
%BLDHSTGRM Вычисляет гистограмму значений (яркости) из массива Data
%Nlvl Количетво уровней в гистограмме
Hstgrm=zeros(Nlvl,1);
i=1:Nlvl;   %итератор уровней квантования
Data=Data(:);%Преобразование массива данных в одномерный
if nargin==2
    LowLimit=min(Data');     %Нижний предел знаений (яркости)
    HighLimit=max(Data');    %Верхний предел
end
Step=(HighLimit-LowLimit)./Nlvl;%Шаг квантования (яркости) в гистограмме
%HstgrmLimits(i+1)=LowLimit+i.*Step; %Пределы (яркости) для набора статистики (соседние значения - пределы)
HstgrmLimits(i)=LowLimit+i.*Step; %Верхние пределы (яркости) для каждого из жиапазонов
Data=sort(Data);    %Предварительная сортировка
t=1;    %Итератор элементов данных (пикселей)
i=1;    %Итератор уровней гистограммы
while t<=size(Data,1);  %Цикл сортировки разбора пикселей по яркости
    if (Data(t)>HstgrmLimits(i))    %Условие повышения уровня квантования
        i=find(HstgrmLimits>=Data(t),1);    %Поиск подходящего диапазона в гистограмме (по верхней границе)
    end
    Hstgrm(i)=Hstgrm(i)+1;  %Инкремент элемента гистограммы
    t=t+1;  %Инкремент элемента данных
end
% for i=1:Nlvl    %Цикл заполнения гистограммы
%     BoolMap=Data>=HstgrmLimits(i);  %условие минимума яркости текущего предеа
%     BoolMap=BoolMap.*(Data<HstgrmLimits(i+1));%условие максимума яркости
%     Hstgrm(i)=sum(BoolMap);
% end

end

