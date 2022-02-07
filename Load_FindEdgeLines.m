function [ StrArray] = Load_FindEdgeLines( Data, AngDynDiff)
%FINDEDGELINES выявляет прямые границы обьектов
%Data - исследуемый растр
%Coord - координаты границ отрезков ([x1 x2 y1 y2] - в каждой строке для каждого отрезка) 
% (они же с1 с2 r1 r2), т.е. ось У направлена вниз а ноль координат - левый верхний угол
%AddData - Дополнительные данные по каждой прямой ([CPR CPC GrAz])([Средневзвешенная точка([Row Col])по величинам градиентов пикселей, Угол градиента в градусах]).
% Азимут градиента всегда на 90 градусов отличается от азимута прямой.
%Rlablty - параметры достоверности для каждой прямой ([TtlGrVal Qpix AzCosDev]) - Суммарная 
% величина градиента по всем пикселям группы, Количество пикселей в группе, Отношение результирующей величины градиента группы к сумме градиентов по каждому пикселю (по сути косинус средневзвешенного стандартного отклонения углов азимутов по каждому пикселю из группы).
%Параметры точность/быстродействие
StrMinLngth=10; %10 Минимальная длина искомых прямых (в пикселях)
StrMaxLngth=150;%Максимальная длина прямой, которая может быть частью строения
AzOwrleapFactor=4.5;  %4.3;4 Коэффициент перекрытия по азимуту (Для выделения прямых на слое градиентов формируется AzStepN опорных азимутов. Для каждого опорного азимута вычисляется слой-разница_азимутов
% с последующей бинаризацией. 360/AzStepN градусов - сектор опорного азимута. Коэффициент перекрытия - это диапазон углов величиной в AzOwrleapFactor секторов опорных азимутов. При вычислении слоя разницы_азимутов берутся пиксели с азимутом градиента входящим в пределы заданного интервала)
MinPixQnttyOnStr=round(StrMinLngth*sqrt(2)); %Минимальное кол-во пикселей в искомых прямых

Coord=zeros(100000,4);  %Инициализация выходных массивов
AddData=zeros(100000,3);
Rlablty=zeros(100000,3);
Data=double(Data);
SafeCoreFrame=2;    %SafeCoreFrame - полуразмер без полупикселя ядра Собеля (и\или) ядра размытия

DiifKernel=GetRoundDiffKernel(5, [3 1]);    %Градиентное ядро (знаки как у ядра Собеля)
  DY=conv2(Data, -DiifKernel,'valid');    %Черезстрочная производная по сверху-вниз. Поскольку координаты связаны с растром (который имеет адреса R и C) Системма координат получаемых прямых отображена по вертикали, т.е ось Y направлена в обратную сторону, а вращение радиус-вектора в полярной СК происходит по часовой стрелке с ростом угла
  DX=conv2(Data, -DiifKernel','valid');   %Черезстобцовая производная по с лева на право
    %quiver(DX,DY);
DiffVal=sqrt(DX.^2+DY.^2);  %Величина градиента для каждого пикселя
    %imagesc(DiffVal);
DiffAngle=atan2(DY,DX);     %Азимут градиента для каждого пикселя (от -pi до pi)
%DiffAngle=DiffAngle./pi.*180;%Перевод азимута в градусы (от -180 до 180 0 градусов - ось Х вправо)

BnrztnLvl=mean(mean(AngDynDiff))*1.5;   %1.5 1.8<----ПОРОГ БИНАРИЗАЦИИ
% BitMsk=AngDynDiff>BnrztnLvl;        %Битовая маска условно ярких пикселей
% [ObjLayer Num]=bwlabel(BitMsk,8);   %Обьединение пикселей в маске на обьекты
%---Вариант с бинаризациями по каждому азимуту
FrdgePoints=zeros(4,2); %координаты пересечения прямой текущей группы с прямыми рамки, которая устанвлена пределами координат текущей группы пикселей
cntr=1;     %Начальная инициализация счётчика найденных прямых
%AzStepN=90;
AzStepN=72; %72Количество шагов по азимуту
    tic;
for nAz=1:AzStepN
  Az=-pi+nAz/AzStepN*2*pi;    %Приведение опорного азимута к радианам
  AzDiffLayer=abs(DiffAngle-Az);  %Слой раницы угла градиента и текущего азимута
  AzDiffLayer(AzDiffLayer>pi)=abs(AzDiffLayer(AzDiffLayer>pi)-2*pi); %Разница от 0 до 180 градусов
  AzDiffLayer=AngDynDiff./(AzDiffLayer./(AzOwrleapFactor*2*pi/AzStepN));	%Взятие обратной величины и домножение на слой динамического фильтра
  BitMsk=AzDiffLayer>BnrztnLvl;   %Бинаризация            %(2*360/AzStepN) - разница в 2*х (х=градусов в одном шаге) соответствует показателю в 1 на обратном слое разницы
  [ObjLayer, Num]=bwlabel(BitMsk,8);   %Разбиение пикселей на группы
  ObjLayer=ObjLayer(:);   %Разименование слоя обьектов (растр теперь одномерный)
  ObjLayer=cat(2,ObjLayer,(1:size(ObjLayer,1))'); %добавляем второй столбец - адреса массива
  ObjLayer=sortrows(ObjLayer);    %Сортировка по возростанию номера обьекта (первого столбца)
  ObjLayer(1:find(ObjLayer(:,1),1)-1,:)=[];   %Удаление пикселей на которыйх не расположен обьект (нулевой номер)
    i=1;   %Итератор пикселей
    ObjLayer=cat(1,ObjLayer,[Num+1 0]); %Добавление последней (сигнальной) строки для обработки всех обьектов (включая последний) внутри нижеследующего цикла.
    for nObj=1:Num %Итератор обьектов
        AdrA=i;     %Начальный адрес диапазона
        while(ObjLayer(i,1)<nObj+1) %Считаем количество пикселей в текущей группе (nObj)
            i=i+1;
        end
        AdrB=i-1;   %Конечный адрес диапазона
         if (AdrB-AdrA)>=MinPixQnttyOnStr-1;    %Условие достаточности количества пикселей для дальнейшей обработки
            CartesVectSumX=sum(DX(ObjLayer(AdrA:AdrB,2)));  %Сумма весторов градиента по пикселям текущей группы
            CartesVectSumY=sum(DY(ObjLayer(AdrA:AdrB,2)));
            GrpAzmth=atan2(CartesVectSumY,CartesVectSumX); %Азимут суммарного вектора (средневзвешенный азимут группы)
            TtlPixGrVal=sum(DiffVal(ObjLayer(AdrA:AdrB,2)));   %Суммарный градиент группы (попиксельная сумма)
            PixW=DiffVal(ObjLayer(AdrA:AdrB,2))./TtlPixGrVal;  %Веса пикселей (пронормированные)
            [r, c]=ind2sub(size(DiffVal),ObjLayer(AdrA:AdrB,2)); %Адреса (строки и столбцы) каждой точки в текущей группе
            LftFr=min(c);   %Левый край группы (вертикальная прямая)
            RghtFr=max(c);  %Правый край группы
            UpFr=min(r);    %Верхний край группы (горизонтальная прямая)
            DwnFr=max(r);   %Нижний край
            CPR=sum(r.*PixW);   %Средневзвешенное положение центра текущей группы пикселей (строки)
            CPC=sum(c.*PixW);                                                             %(столбцы)
            k=tan(GrpAzmth+pi/2);  %К-фт для уравнения y=kx
            FrdgePoints(1,:)=[LftFr k*(LftFr-CPC)+CPR];   %X и Y точки пересечения прямой левого края рамки с прямой аппроксимирующей текущую группу
            FrdgePoints(2,:)=[(UpFr-CPR)/k+CPC UpFr];   %точка пересечения прямой верхнего края рамки с прямой аппроксимирующей текущую группу
            FrdgePoints(3,:)=[RghtFr k*(RghtFr-CPC)+CPR];   %то же для правого края
            FrdgePoints(4,:)=[(DwnFr-CPR)/k+CPC DwnFr];   %то же для нижнего края
            [~,on]=inpolygon(FrdgePoints(:,1),FrdgePoints(:,2),([LftFr RghtFr RghtFr LftFr]),([UpFr UpFr DwnFr DwnFr]));    %выбор точек пересечения с рамкой в пределах её отрезков
              Coord(cntr,1:2)=FrdgePoints(on,1);    %X1X2 текущей прямой
              Coord(cntr,3:4)=FrdgePoints(on,2);    %Y1Y2 текущей прямой
            if (sqrt((Coord(cntr,1)-Coord(cntr,2))^2+sqrt((Coord(cntr,3)-Coord(cntr,4))^2)))>StrMinLngth    %Проверка длинны отрезка (слишком короткие отбрасываются)
              AddData(cntr,:)=[CPR,CPC,GrpAzmth];    %Запись дополнительных данных (строка+столбец центральной точки отрезка (не целые), средневзыешенный азимут градиента)
              TtlGrVal=sqrt(CartesVectSumX^2+CartesVectSumY^2); %Величина градиента группы (векторной суммы градиентов от каждого пикселя)
              Rlablty(cntr,:)=[TtlGrVal,AdrB-AdrA+1,TtlGrVal/TtlPixGrVal];  %Параметр достоверности по каждой группе (отрезку) - Величина градинта группы, оличество пикселей, косинус расброса углов
              cntr=cntr+1;  %Инкремент счётчика для следующей прямой
            end                          
         end        
    end
end
cntr=cntr-1;    %Отмена последнего инкремента счётчика
  Coord(cntr+1:end,:)=[];   %Удаление лишней нициализированной области выходных массивов.
  AddData(cntr+1:end,:)=[]; % И одновременное удаление последней линии (кроме случая когда последняя линия имеет подходящую длину)
  Rlablty(cntr+1:end,:)=[];
    toc;
%-----Прореживание и сортировка списка прямых---
StrLngth=sqrt((Coord(:,2)-Coord(:,1)).^2+(Coord(:,4)-Coord(:,3)).^2);   %Длины прямых
%NrmlzdGrad=Rlablty(:,1)./Rlablty(:,2);  %Величина суммарного попиксельного градиента отнесённая к количеству пикселей гаждой группы (каждой прямой)
    %NrmlzdGrad=Rlablty(:,1)./StrLngth(:);
%TrustInd=sqrt((StrLngth./mean(StrLngth,1)).^2+(NrmlzdGrad./mean(NrmlzdGrad,1)).^2+(Rlablty(:,3)./mean(Rlablty(:,3),1)).^2); %Индекс значимости прямой - эвклидова метрика в пространстве Длина-Величина_градиента-Разброс_азимутов для группы пикселей по каждой прямой
TrustInd=Rlablty(:,3);
    %TrustInd=sqrt((Rlablty(:,3)./mean(Rlablty(:,3))).^2+(NrmlzdGrad./mean(NrmlzdGrad)).^2);
%TrustInd=sqrt(6.*(Rlablty(:,3)./mean(Rlablty(:,3),1)).^2+(Rlablty(:,2)./mean(Rlablty(:,2),1)).^2);

StrArray=cat(2,Coord,AddData,Rlablty,TrustInd,StrLngth); %Обьединение всех параметров найденных прямых в один массив (по строке для каждой прямой)
StrArray(:,1:6)=StrArray(:,1:6)+SafeCoreFrame;  %Учёт (запасной) рамки растра от свёртки 
StrArray=sortrows(StrArray,-11);    %Сортировка прямых по значимости
StrArray(:,11)=1:cntr;      %Замена индекса занчимости на Индентификатор прямой (чем меньше ID, тем выше значимость)
%Обединение прямых в группы
ToStrDltMsk=false(cntr,1);
AngTreshold=1.8*(2*pi/AzStepN); %1.8 %Норма разницы по азимуту для группы прямых (+-AngTreshold радиан от значения азимута текущей прямой (текущий ID))
%Нормы разницы центров прямых приведены как коэффициенты для длины инициирующей прямой (прямой с текущим ID)
PrllelKTrshold=1/8;	%1/4Норма разницы координат центра отрезка по оси параллельной прямой с текущим ID (в долях длины текцщей прямой)
%AcrossKTrshold=12.0*StrMinLngth.*StrArray(:,12)./StrMaxLngth;	%Предел разницы координат по оси перпендикулярной текущей прямой (в пикселях)
AcrossKTrshold=0.5*StrMinLngth;	%Предел разницы координат по оси перпендикулярной текущей прямой (в пикселях)
    %LngthMltTrshold=2; %Если проверяемая прямая длиннее инииирующей (текущей) более чем в LngthMltTrshold раз, то её отбрасывать не следует
    LngthTrshold=0.5;   %Предел разницы длины текущей и проверяемых прямых (в долях длины текущей прямой) - много-большие и многоменьшие - не должны удалятся
    %Выриант с битовой маской по всем прямым:
    tic;
GK=tan(StrArray(:,7));    %Коэффициент прямой градиента (перпендикуляра) для уравнения вида y=kx
CrntIDStrK=tan(StrArray(:,7)-pi/2);%Коэффициент текущей (иницииирующей) прямой для уравнения вида y=kx
SameSourseStrMterick=zeros(cntr,1);                     %Инициализация массива метрик связности прямых
for ID=1:cntr-1
  if ToStrDltMsk(ID)==false
    if StrArray()
      SameSourseStrMterick(ID+1:end)=abs(StrArray(ID+1:end,7)-StrArray(ID,7)); %Определение разниц азимутов инициирующей прямой и всех остальных
      SameSourseStrMterick(SameSourseStrMterick>=pi)=abs(SameSourseStrMterick(SameSourseStrMterick>=pi)-2*pi);%Приведение азимутальной разницы к пределам от 0 до 180 градусов
      %GK=tan(StrArray(ID+1:end,7));    %Коэффициент прямой градиента (перпендикуляра) для уравнения вида y=kx
      %CrntIDStrK=tan(StrArray(ID,7)-pi/2);%Коэффициент текущей (иницииирующей) прямой для уравнения вида y=kx
      %Координаты точки пересечения линий градиента с текущей прямой
      x=(CrntIDStrK(ID)*StrArray(ID,6)-GK(ID+1:end).*StrArray(ID+1:end,6)+StrArray(ID+1:end,5)-StrArray(ID,5))./(CrntIDStrK(ID)-GK(ID+1:end));
      y=CrntIDStrK(ID).*(x-StrArray(ID,6))+StrArray(ID,5); %Координаты точек пересечения линий градиента от каждой прямой и текущей прямой
      PrllelDist=sqrt((x-StrArray(ID,6)).^2+(y-StrArray(ID,5)).^2)./StrArray(ID,12);	%Расстояние (в долях длины текущей прямой) между центральными точками вдоль оси параллельной инициирующей прямой (расстояние от точки пересечения линии градиента с текущей прямой до центральной точки текущей(инициирующей) прямой)
      AcrossDist=sqrt((x-StrArray(ID+1:end,6)).^2+(y-StrArray(ID+1:end,5)).^2);    %Аналогично расстояние (в пикселях) между прямыми (каждой и текущей) вдоль линии градиента (перпендикуляра каждой линии).
      LngthDiff=abs(StrArray(ID,12)-StrArray(ID+1:end,12))/StrArray(ID,12);      %Разница длинн текущей и проверяемых прямых в долях текущей прямой
      %SameSourseStrMterick(ID+1:end)=sqrt((SameSourseStrMterick(ID+1:end)./AngTreshold).^2+(PrllelDist./PrllelKTrshold).^2+(AcrossDist./AcrossKTrshold).^2);  %Мера сопряжённости прямых (на сколько для своего построения они используют одинаковые с инициирующей прямой пиксели) - эвклидова метрика в (линейнозависимом) пространстве Разница_азимутов - Параллельное_расстояние - Градиентное_расстояние. Все измерения пронормированы
      SameSourseStrMterick(ID+1:end)=sqrt((SameSourseStrMterick(ID+1:end)./AngTreshold).^2+(PrllelDist./PrllelKTrshold).^2+(AcrossDist./AcrossKTrshold).^2+(LngthDiff/LngthTrshold).^2);  %Мера сопряжённости прямых (на сколько для своего построения они используют одинаковые с инициирующей прямой пиксели) - эвклидова метрика в (линейнозависимом) пространстве Разница_азимутов - Параллельное_расстояние - Градиентное_расстояние. Все измерения пронормированы
      %SameSourseStrMterick(ID+1:end)=sqrt((SameSourseStrMterick(ID+1:end)./AngTreshold).^2+(PrllelDist./PrllelKTrshold).^2+(AcrossDist./AcrossKTrshold(ID+1:end)).^2);  %Мера сопряжённости прямых (на сколько для для своего построения они используют одинаковые с инициирующей прямой пиксели) - эвклидова метрика в (линейнозависимом) пространстве Разница_азимутов - Параллельное_расстояние - Градиентное_расстояние. Все измерения пронормированы
      %ToStrDltMsk(ID+1:end)=ToStrDltMsk(ID+1:end)|SameSourseStrMterick(ID+1:end)<sqrt(3);
      %ToStrDltMsk(ID+1:end)=ToStrDltMsk(ID+1:end)|((SameSourseStrMterick(ID+1:end)<sqrt(3))&(LngthMltTrshold*StrArray(ID,12)>StrArray(ID+1:end,12)));
      ToStrDltMsk(ID+1:end)=ToStrDltMsk(ID+1:end)|SameSourseStrMterick(ID+1:end)<2; %2=sqrt(4) - для четырёхмерного пространства
    end
  end
end
StrArray(ToStrDltMsk,:)=[];    %Удаление вторичных связанных прямых
Coord=StrArray(:,1:4);  %Выгрузка искомых массивов параметров прямых в порядке отсортированном по значимости прямой
AddData=StrArray(:,5:7);
Rlablty=StrArray(:,8:10);
    toc;

imagesc(Data(1:365+60,1:2412+60));
hold on
% axis image
% for i=1:cntr
%     if ToDltMsk(i)==true
%   color=StrArray(i,10);
%   plot(StrArray(i,1:2),StrArray(i,3:4),'Color',[color,color,color]);
%     end
% end
% plot(StrArray(ID,1:2),StrArray(ID,3:4),'red');
% hold off


% imagesc(Data);
% hold on
% axis image
for i=1:size(Coord,1)
    %if Coord(i,2)<=800
      if Rlablty(i,3)>1   %Коррекция аномального индекса достоверности
        disp('Ошибка параметра разброса пиксльных градиентов прямой - скорректировано');
        Rlablty(i,3)=1;
      end
  color=Rlablty(i,3);
  plot(Coord(i,1:2),Coord(i,3:4),'Color',[color,color,color])
    %end
end
hold off
    

end

