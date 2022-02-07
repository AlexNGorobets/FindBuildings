function [ Coord ] = Slow_FindEdgeLines( Data ,Upfr,DwnFr,LftFr,RghtFr)
%FINDEDGELINES выявляет прямые границы обьектов
%Data - исследуемый растр
%Coord - координаты границ отрезков ([x1 y1 x2 y2] - в каждой строке для каждого отрезка)

%Параметры точность/быстродействие
ColNeighborN=2; %Колличество окрестных пикселей-кандидатов по измерению колонок (вдоль каждой строки группы), при каждой итерации добора пикселей в текущую группу
RowNeighborN=2; %Колличество окрестных пикселей-кандидатов по строкам ("окрестности" области пикселей (текущей) группы вдоль столбца-сечения)
LstAddPixAmmount=0; %Цикл, который добавляет окрестные пиксели в (текущую) группу завершится, когда на последне итерации в группу будет добавлено LstAddPixAmmount пикселей (=0 - наивысшая точность)
QPixTresholdInGroup=10; %Минимальное количество пикселей в группе для того, чтобы эта группа могла быть исследована на предмет вписанной прямой
%Параметры для сегментации: (все параметры будут пронормированы на приведеные значения,
% с последующим взятием корня их суммы квадратов). Пиксель входит в группу, если этот параметр меньше 1
ColTrshold=2.2; %Норма расстояния по X (столбцы вправо) до ближайшего пикселя из группы
RowTrshold=2.2; %Норма расстояния по Y (строки вниз)
GrAngTrshold=25;%Норма разницы углов в градусах между углом градиента текущего пикселя и средним углом в (текущей) группе.
GrValTrshold=40;%Норма разницы величин градиентов текущего пикселя и среднего значения (текущей) группы (в процентах от ср.знач. группы).
GrValTrshold=GrValTrshold/100;      %Перевод в доли
        SafeCoreFrame=4;    %SafeCoreFrame - полуразиер без полупикселя ядра Собеля (и\или) ядра размытия
        Upfr=Upfr+SafeCoreFrame;
        DwnFr=DwnFr-SafeCoreFrame;
        LftFr=LftFr+SafeCoreFrame;
        RghtFr=RghtFr-SafeCoreFrame;
Data=double(Data);  %Приведение типа данных к double
i=2:2:size(Data,1)-2;   %В полученных снимках каждая вторая строка равна первой
Data(i,:)=(Data(i-1,:)+Data(i+1,:))/2;	% исправим локальным усреднением
    %imshow(uint16(Data(Upfr:DwnFr,LftFr:RghtFr)));
GaussKernel=fspecial('gaussian',[5 5],1.6); %формирование ядра размытия
BluredData=conv2(Data, GaussKernel,'same');	%Размытие по Гауссу
        BluredData=BluredData(Upfr:DwnFr,LftFr:RghtFr);
    %imshow(uint16(Data(Upfr:DwnFr,LftFr:RghtFr)));
DiifKernel=GetRoundDiffKernel(5, [1 3]);
DY=conv2(BluredData, DiifKernel,'valid');    %Черезстрочная производная по сверху-вниз
DX=conv2(BluredData, DiifKernel','valid');   %Черезстобцовая производная по с лева на право
DiffVal=sqrt(DX.^2+DY.^2);  %Величина градиента для каждого пикселя
DiffAngle=atan2(DY,DX);     %Азимут градиента для каждого пикселя
DiffAngle=DiffAngle./pi.*180;%Перевод азимута в градусы
%     ValIm=DiffVal(Upfr:DwnFr,LftFr:RghtFr);
%     AngIm=DiffAngle(Upfr:DwnFr,LftFr:RghtFr);
%     imagesc(AngIm);
%     % AngIm=SDiffAngle(Upfr:DwnFr,LftFr:RghtFr);
%     dBLPict=cat(4,Data(Upfr:DwnFr,LftFr:RghtFr),BluredData(Upfr:DwnFr,LftFr:RghtFr));
% %     subplot(2,1,1),montage(uint16(dBLPict),'Size',[1 2]);
% %     subplot(2,1,2),quiver(DX(Upfr:DwnFr,LftFr:RghtFr),DY(Upfr:DwnFr,LftFr:RghtFr));
%     quiver(DX(Upfr:DwnFr,LftFr:RghtFr),DY(Upfr:DwnFr,LftFr:RghtFr));
%Сегментация пикселей по расстоянию (друг от друга), общему азимуту и общей величине градиента

[DwnFr RghtFr]=size(BluredData);
%GCoord=zeros(400,2);    %Координаты каждой точки входящей в текущую группу [row col] (накапливаются в порядке гипотез по пикселям)
GCoord(1,:)=[RowNeighborN ColNeighborN];	%Первое приближние (первый пиксль в группе)
        GCoord(1,:)=[27 91];
        %GCoord(1,:)=[20 104];
        %GCoord(1,:)=[16 105];
GMeanAng=DiffAngle(GCoord(1,1),GCoord(1,2)); %Средний азимут градиента в группе (первое приближение)
GMeanVal=DiffVal(GCoord(1,1),GCoord(1,2));   %Средняя величина градиента группы (первое приближение)
    PixGroups=cell(100,1);
    PixUsedMask=false(size(BluredData));%Инициализация маски пикселей, которые уже отнесены в к.л. группу
PixOffrdMask=PixUsedMask;       	%Инициализация маски пикселей, по которым уже были выдвинуы гипотезы (маска нужна для новых "первых приближений")
CrntGroupPixMask=PixUsedMask;	%Маска пикселей взятых в текущую группу (для реализации, когда пиксель может принадлежать нескольким группам)
    %     CrntGroupOffrdMask=PixUsedMask;
    %     CrntGroupOffrdMask=double(CrntGroupOffrdMask);
                tic
g=1;    %Итератор группы
while (g<=2 || ~isempty(GCoord)) %Цикл-итератор групп

                tic;
        CrntGroupPixMask(:,:)=false;
        CrntGroupPixMask(GCoord(1,1),GCoord(1,2))=true;%Отметка первого приближения
        PixOffrdMask(GCoord(1,1),GCoord(1,2))=true;
    Offred_IND=1;   %Начальная инициализация адреса предложенного пикселя (для входа в нижеслдующий цикл)
    while (size(Offred_IND,2)>LstAddPixAmmount)  %Цикл итераций приближений для каждой группы
        Offrd_ColDist=[];   %Начальные инициализации
        Offrd_RowDist=[];
        PixCntr=0;
        for r=min(GCoord(:,1)):max(GCoord(:,1))	%Перебор строк
          TakenPixRange=find(CrntGroupPixMask(r,:));%Определение координат пикселей в текущей строке, которые уже в группе
          if ~isempty(TakenPixRange)    %По строке в которой нет ни одного пикселя из группы не будет предложено ни одного пикселя
            if TakenPixRange(1)>ColNeighborN  %Защита от выхода за пределы массива данных по левому краю
                StartAdr=TakenPixRange(1)-ColNeighborN;   
            else StartAdr=1;        %Начало диапазона адресов предлогаемых пикселей
            end
            if TakenPixRange(end)<(RghtFr-ColNeighborN)   %Защита от выхода за пределы массива данных по правому краю
                FinshAdr=TakenPixRange(end)+ColNeighborN;
            else FinshAdr=RghtFr;   %Конец диапазона адресов предлогаемых пикселей
            end
            OfferedAdrss=StartAdr:FinshAdr; %Формирование списка адресов (пикселей, предлогаемых к зачислению в группу)
            for t=1:size(TakenPixRange,2) %Удаление из диапазона адресов пикселей, которые уже в группе
             OfferedAdrss(OfferedAdrss==TakenPixRange(t))=[];            
            end
            Offrd_ColDist=cat(2,Offrd_ColDist,zeros(size(OfferedAdrss)));	%Инициализация растояний
            for p=1:size(OfferedAdrss,2)	%Получение соответствующих расстояний вдоль (текущей) строки от каждого из предложенных пикселей до ближайшего из тех, что уже зачислены в (текущую) группу
             Offrd_ColDist(PixCntr+p)=min([TakenPixRange(find(TakenPixRange>OfferedAdrss(p),1))-OfferedAdrss(p) OfferedAdrss(p)-TakenPixRange(find(TakenPixRange<OfferedAdrss(p),1,'last'))]);
            end
            Offrd_r(PixCntr+1:PixCntr+p)=r;             %Строки (координаты) предложенных пикселей (Y-вниз координаты)
        	Offrd_c(PixCntr+1:PixCntr+p)=OfferedAdrss;  %Стобцы предложенных пикселей (X-вправо координаты). Порядковые номера координат в обоих массивах соответствуют одним и тем же точкам
            PixCntr=PixCntr+p;          %Обновляем значение счётчика предложенных пикселей
          end
        end        
        for c=min(GCoord(:,2)):max(GCoord(:,2))	%Перебор столбцов
          TakenPixRange=find(CrntGroupPixMask(:,c))';%Определение координат пикселей в текущем стобце, которые уже в группе
          if ~isempty(TakenPixRange)    %По столбцу в котором нет ни одного пикселя из группы - не будет предложено ни одного пикселя            
            if TakenPixRange(1)>RowNeighborN  %Защита от выхода за пределы массива данных по верхнему краю
                StartAdr=TakenPixRange(1)-RowNeighborN;   
            else StartAdr=1;        %Начало диапазона адресов предлогаемых пикселей
            end
            if TakenPixRange(end)<(DwnFr-RowNeighborN)   %Защита от выхода за пределы массива данных по нижнему краю
                FinshAdr=TakenPixRange(end)+RowNeighborN;
            else FinshAdr=DwnFr;   %Конец диапазона адресов предлогаемых пикселей
            end
            OfferedAdrss=StartAdr:FinshAdr; %Формирование списка адресов (пикселей, предлогаемых к зачислению в группу)
            for t=1:size(TakenPixRange,2) %Удаление из диапазона адресов пикселей, которые уже в группе
             OfferedAdrss(OfferedAdrss==TakenPixRange(t))=[];            
            end
            Offrd_RowDist=cat(2,Offrd_RowDist,zeros(size(OfferedAdrss)));    %Инициализация растояний
            for p=1:size(OfferedAdrss,2)	%Получение соответствующих расстояний вдоль (текущего) столбца от каждого из предложенных пикселей до ближайшего из тех, что уже зачислены в (текущую) группу
             Offrd_RowDist(PixCntr+p)=min([TakenPixRange(find(TakenPixRange>OfferedAdrss(p),1))-OfferedAdrss(p) OfferedAdrss(p)-TakenPixRange(find(TakenPixRange<OfferedAdrss(p),1,'last'))]);
            end
            Offrd_ColDist(PixCntr+1:PixCntr+p)=Inf; %Расстояния предложенных пикселей до группы (по измерению столбцов - в данной строке). Здесь только добавляются пометки о том, что эти расстояния не определены - это нужно чтобы размер массивов Offrd_ColDist и Offrd_RowDist совпадали
            Offrd_r(PixCntr+1:PixCntr+p)=OfferedAdrss;  %Строки (координаты) предложенных пикселей (Y-вниз координаты)
        	Offrd_c(PixCntr+1:PixCntr+p)=c;          	%Стобцы предложенных пикселей (X-вправо координаты). Порядковые номера координат в обоих массивах соответствуют одним и тем же точкам
            if size(Offrd_r,2)>size(Offrd_ColDist,2)
              FndedDoubles=0;   %Счётчик найденных дублей
              for t=1:p     %Цикл удаляет повторяющиеся адреса (пиксели предложенные дважды)
              	CoordAdr1=find(Offrd_r(1:PixCntr)==Offrd_r(PixCntr+t-FndedDoubles));  %Поиск повторений в номерах строк
                CoordAdr2=find(Offrd_c(CoordAdr1)==Offrd_c(PixCntr+t-FndedDoubles), 1);   %Поиск повторений в столбцах (по адресам повторяющихся строк)
                if ~isempty(CoordAdr2)
                  Offrd_r(PixCntr+t-FndedDoubles)=[];
                  Offrd_c(PixCntr+t-FndedDoubles)=[];
                  Offrd_RowDist(CoordAdr1(CoordAdr2))=Offrd_RowDist(PixCntr+t-FndedDoubles);
                  Offrd_RowDist(PixCntr+t-FndedDoubles)=[];
                  Offrd_ColDist(PixCntr+t-FndedDoubles)=[];
                  FndedDoubles=FndedDoubles+1;
                end                
              end
              PixCntr=PixCntr+p-FndedDoubles;	%Обновляем значение счётчика предложенных пикселей
            else
              PixCntr=PixCntr+p;       %Обновляем значение счётчика предложенных пикселей    
           end
            
          end
        end
        Offred_IND=sub2ind(size(CrntGroupPixMask),Offrd_r,Offrd_c); %Приведение индексов к последовательной нумерации
        	PixOffrdMask(Offred_IND)=true;
        %Проверка предложенных пикселей на принадлежность к группе
        PixAng=abs(DiffAngle(Offred_IND)-GMeanAng);     %Разницы АЗИМУТОВ предложенных пиксеей и среднего по группе
        PixAng(PixAng>180)=abs(PixAng(PixAng>180)-360); %Приведение разниц азимутов к интервалу от 0 до 180 градусов
        PixAng=PixAng/GrAngTrshold;                     %Нормировка на treshold-значение 
        %PixVal=abs(DiffVal(Offred_IND)-GMeanVal)/GMeanVal/GrValTrshold;  %Разницы величин ГРАДИЕНТОВ предложенных пиксеей и среднего по группе в долях от градиента группы, нормрованные на treshold-значение
        PixVal=abs(DiffVal(Offred_IND)-GMeanVal)/GMeanVal/GrValTrshold;  %Разницы величин ГРАДИЕНТОВ предложенных пиксеей и среднего по группе в долях от градиента группы, нормрованные на treshold-значение
        Offrd_RowDist(Offrd_RowDist==0)=Inf;	%Отбрасываем все нулевые значения строчечных расстояний
        PixRDist=abs(RowTrshold-Offrd_RowDist)/RowTrshold; %Разницы вдоль столбца (по измерению строк) до ближайшего пикселя из группы
        Offrd_ColDist(Offrd_ColDist==0)=Inf;	%Отбрасываем все нулевые значения строчечных расстояний
        PixCDist=abs(ColTrshold-Offrd_ColDist)/ColTrshold; %Разницы вдоль столбца (по измерению строк) до ближайшего пикселя из группы
        PixDist=min(cat(1,PixRDist,PixCDist),[],1);             %Выбор кратчайшго расстояния (между строчным и столбцовым)
        %if ~((size(PixAng,1)==size(PixVal,1))&&(size(PixAng,1)==size(PixDist,1)) && ((size(PixAng,2)==size(PixVal,2))&&(size(PixAng,2)==size(PixDist,2))))
        if ~(size(PixAng,2)==size(PixDist,2))
            disp('ошибочка')
        end
        ApprovmntMtric=sqrt(PixAng.^2+PixVal.^2+PixDist.^2);%Пиксель входит в группу если Эвклидово расстояние в Декартовом 
        Offred_IND=Offred_IND(ApprovmntMtric<=1);           % пространстве Азимут-Градиент-(кратчайшее)Расстояние, в базисе пронормированном на соответствующие treshold-значения меньше единицы
        CrntGroupPixMask(Offred_IND)=true;      %Отметка принятых в группу пикселей в маске текущей группы
        GCoord=cat(1,GCoord,cat(1,Offrd_r(ApprovmntMtric<=1),Offrd_c(ApprovmntMtric<=1))'); %Занесение координат принятых пикселей в список пикселей группы
        %Модификация средних показателей группы
        GMeanAng=mean(DiffAngle(CrntGroupPixMask));
        GMeanVal=mean(DiffVal(CrntGroupPixMask));     

                    %CrntGroupOffrdMask(Offred_IND)=CrntGroupOffrdMask(Offred_IND)+1;
                 PlotBitArray=CrntGroupPixMask;
                 PlotBitArray=double(PlotBitArray);
                 PlotBitArray(Offred_IND)=PlotBitArray(Offred_IND)+1;
                 imagesc(PlotBitArray);
                 pause(0.1);
            
    end
            toc;
    if size(GCoord,1)>QPixTresholdInGroup
        PixGroups{g}=GCoord;
            PixUsedMask(GCoord(:,1),GCoord(:,2))=true;
        g=g+1;
    end
    %Выбор начального приближение для следующей группы
    GCoord=zeros(1,2);
    [GCoord(1,1), GCoord(1,2)]=find(PixOffrdMask==0,1);
    Offrd_r=[];
    Offrd_c=[];
    

Coord=0;
end
                toc;
end

