function [ StrArray ] = FindEdgeLines( Data )
%FINDEDGELINES выявляет прямые границы обьектов
%Data - исследуемый растр
%Coord - координаты границ отрезков ([x1 x2 y1 y2] - в каждой строке для каждого отрезка) 
% (они же с1 с2 r1 r2), т.е. ось У направлена вниз а ноль координат - левый верхний угол
%Параметры точность/быстродействие
CoreSize=17;    %Рязмер ядра адаптивного фильтра
AdptvFltrngAzW=25;  %25Глубина обработки динамическим фильтром - Относительное влияние азимута градиента при фильтрации (1 - разница азимутов считается в радианах, )
Data=double(Data);  %Приведение типа данных к double
i=2:2:size(Data,1)-2;   %В полученных снимках каждая вторая строка равна первой
Data(i,:)=(Data(i-1,:)+Data(i+1,:))/2;	% исправим локальным усреднением
[SrsRN, SrsCN]=size(Data);   %Размер исходных данных
%Дополнение снимка рамкой толщиной в половину ядра динамического фильтра (чтобы результат динамической фильтрации по размеру совпадал со входными данными):
CatArr=zeros((CoreSize-1)/2,size(Data,2));
CatArr(1,:)=Data(1,:);
CatArr=cumsum(CatArr,1)+rand(size(CatArr));
Data=cat(1,CatArr,Data); %Дополнение по первому измерению (строк сверху)
CatArr(:,:)=0;
CatArr(1,:)=Data(end,:);
CatArr=cumsum(CatArr,1)+rand(size(CatArr));
Data=cat(1,Data,CatArr); %Дополнение по первому измерению (строк снизу)
CatArr=zeros(size(Data,1),(CoreSize-1)/2);
CatArr(:,1)=Data(:,1);
CatArr=cumsum(CatArr,2)+rand(size(CatArr));
Data=cat(2,CatArr,Data); %Дополнение по второму измерению (стобцов слева)
CatArr(:,:)=0;
CatArr=zeros(size(Data,1),(CoreSize-1)/2);
CatArr(:,1)=Data(:,end);
CatArr=cumsum(CatArr,2)+rand(size(CatArr));
Data=cat(2,Data,CatArr); %Дополнение по второму измерению (стобцов справа)

%DiffKernel=GetRoundDiffKernel(5, [1 3]);
DiffKernel=GetRoundDiffKernel(5, [3 1]);
%DiffKernel=[0 1 0;0 0 0;0 -1 0];
  DY=conv2(Data, DiffKernel,'same');    %Черезстрочная производная по сверху-вниз
  DX=conv2(Data, DiffKernel','same');   %Черезстобцовая производная по с лева на право
        %quiver(flip(DX(Upfr:DwnFr,LftFr:RghtFr),1),-flip(DY(Upfr:DwnFr,LftFr:RghtFr),1));
        %quiver(flip(DX,1),flip(-DY,1));
%DiffVal=sqrt(DX.^2+DY.^2);  %Величина градиента для каждого пикселя
[DiffAngle,DiffVal]=cart2pol(DX,DY);
%DiffAngle=atan2(DY,DX);     %Азимут градиента для каждого пикселя
%DiffAngle=atan2(-DY,DX);     %Азимут градиента для каждого пикселя (от -pi до pi)
    %DiffAngle=DiffAngle./pi.*180;%Перевод азимута в градусы (от -180 до 180 0 градусов - ось Х вправо)
                tic;        %Метод с ежистрочным суммированием
%**********Блок адаптивной фильтрации
%В Данной реализации векторизация применена вдоль наибольшего измерения(колонки). 
% Для этого итерации будут происходить построчно - В каждой итерации ядро будет
% посчитано как трёхмерный массив Ncol на CoreSize на CoreSize. Массив разницы азимутов будет иметь такой же размер
        HCS=fix(CoreSize/2);	%Half-Core Size
%[Nr, Nc]=size(DiffAngle);   %Определение пределов итераторов для свёртки
AngDynDiff=zeros(SrsRN,SrsCN);    %Инициализация слоя-результата динамического фильтра
        AngDynDiff2=AngDynDiff;
CoreOnString=zeros(SrsCN,CoreSize,CoreSize);   %Инициализация Массива разницы азимутов
%c=1+HCS:size(DiffAngle,2)-HCS;  %вектор - итератор колонок
c=1:SrsCN;  %вектор - итератор колонок
Kernl=CoreOnString;     %Инициализация ядра (для всех столбцов сразу)
DiffValLayer=CoreOnString;%Инициализация слоя-величины градиента
%for r=1+HCS:size(DiffAngle,1)-HCS
for r=1:SrsRN
    for KrnlDr=1:CoreSize
        for KrnlDc=1:CoreSize
        CoreOnString(c,KrnlDr,KrnlDc)=abs(DiffAngle(r+KrnlDr-1,c+KrnlDc-1)-DiffAngle(r+HCS,c+HCS)); %Вычисление разницы азимутов (для каждого пикселя)
        DiffValLayer(c,KrnlDr,KrnlDc)=DiffVal(r+KrnlDr-1,c+KrnlDc-1);   %Копирование величин градиентов в соответствующий слой (для текущих координат)
        end
    end
    CoreOnString(CoreOnString>pi)=abs(CoreOnString(CoreOnString>pi)-2*pi); %Приведение разниц углов в пределы от -pi до +pi и взвешивание
    K=Bld3DAzKernel( CoreSize, DiffAngle(r+HCS,c+HCS)); %Ядра для текущей строки 
        %K=Bld3DAzKernel_Base( CoreSize, DiffAngle(r,c) ); %Ядра для текущей строки 
    Kernl(c,:,:)=permute(K,[3 1 2]);    %Перестановка размерностей (колонки растра - теперь первое измерение)
%               if r==32
%                   tmpC=25;
%                   imagesc(Data(r-HCS:r+HCS,tmpC-HCS:tmpC+HCS));
%                   hold on
%                   %quiver(flip(DX(r-HCS:r+HCS,tmpC-HCS:tmpC+HCS),1),flip(-DY(r-HCS:r+HCS,tmpC-HCS:tmpC+HCS),1));
%                   %quiver(DX(r-HCS:r+HCS,tmpC-HCS:tmpC+HCS),DY(r-HCS:r+HCS,tmpC-HCS:tmpC+HCS));
%                   meshc(0.6+K(:,:,tmpC));
%               end
    %CoreOnString=CoreOnString.*Kernl./DiffValLayer; %Собсно свёртка, с домножением на величину обратную градиенту
        %CoreOnString2=Kernl.*(CoreOnString./DiffValLayer);
    CoreOnString=DiffValLayer.*Kernl./(1+AdptvFltrngAzW.*CoreOnString);
    AngDynDiff(r,c)=sum(sum(CoreOnString(c,:,:),2),3);  %Суммирование по строка и столбцам ядра, результат строка изображения
        %AngDynDiff2(r,c)=1./sum(sum(CoreOnString2(c,:,:),2),3);
end
%AngDynDiff(1+HCS:Nr-HCS,c)=1./AngDynDiff(1+HCS:Nr-HCS,c);   %Взятие обратной величины (1/0 - максимальная яркость соотв-щая минимальной разнице в углах)
clearvars CoreOnString DiffValLayer Kernl K;  %Освобождение памяти
Data([1:HCS,end-HCS+1:end],:)=[];Data(:,[1:HCS,end-HCS+1:end])=[];  %Удаление рамки размером в половину ядра динамического фильтра
DiffAngle([1:HCS,end-HCS+1:end],:)=[];DiffAngle(:,[1:HCS,end-HCS+1:end])=[];  % (возврат к размерам исходных данных)
DiffVal([1:HCS,end-HCS+1:end],:)=[];DiffVal(:,[1:HCS,end-HCS+1:end])=[];
DX([1:HCS,end-HCS+1:end],:)=[];DX(:,[1:HCS,end-HCS+1:end])=[];
DY([1:HCS,end-HCS+1:end],:)=[];DY(:,[1:HCS,end-HCS+1:end])=[];
    disp('динамическая фильтрация:')
        toc;
         %tImg=AngDynDiff(Upfr:DwnFr,LftFr:RghtFr);
     %imagesc(AngDynDiff./mean(mean(AngDynDiff)));axis equal;
     %imagesc(AngDynDiff2./mean(mean(AngDynDiff2)));        
     tmp=1;
%*****Блок получения отрезков
%AddData - Дополнительные данные по каждой прямой ([CPR CPC GrAz])([Средневзвешенная точка([Row Col])по величинам градиентов пикселей, Угол градиента в градусах]).
% Азимут градиента всегда на 90 градусов отличается от азимута прямой.
%Rlablty - параметры достоверности для каждой прямой ([TtlGrVal Qpix AzCosDev]) - Суммарная 
% величина градиента по всем пикселям группы, Количество пикселей в группе, Отношение результирующей величины градиента группы к сумме градиентов по каждому пикселю (по сути косинус средневзвешенного стандартного отклонения углов азимутов по каждому пикселю из группы).
StrMinLngth=10; %10 Минимальная длина искомых прямых (в пикселях)
AzOwrleapFactor=2.5;  %2.5 1.5; 4.3;4 Коэффициент перекрытия по азимуту (Для выделения прямых на слое градиентов формируется AzStepN опорных азимутов. Для каждого опорного азимута вычисляется слой-разница_азимутов
% с последующей бинаризацией. 360/AzStepN градусов - сектор опорного азимута. Коэффициент перекрытия - это диапазон углов величиной в AzOwrleapFactor секторов опорных азимутов. При вычислении слоя разницы_азимутов берутся пиксели с азимутом градиента входящим в пределы заданного интервала)
MinPixQnttyOnStr=round(StrMinLngth)./sqrt(2); %Минимальное кол-во пикселей в искомых прямых

Coord=zeros(100000,4);  %Инициализация выходных массивов
AddData=zeros(100000,3);
Rlablty=zeros(100000,3);
%SafeCoreFrame=2;    %SafeCoreFrame - полуразмер без полупикселя ядра Собеля (и\или) ядра размытия
%BnrztnLvl=mean(mean(AngDynDiff))*1.5;   %1.2 0.8 1.0;1.5 1.8<----ПОРОГ БИНАРИЗАЦИИ
    %BnrztnLvl=2.326762770949994e+05
    BnrztnLvl=2.15e+05;
%BnrztnLvl=graythresh(AngDynDiff);
      MainBitMsk=AngDynDiff>BnrztnLvl;        %Битовая маска условно ярких пикселей
      %imshow(MainBitMsk);
FrdgePoints=zeros(4,2); %координаты пересечения прямой текущей группы с прямыми рамки, которая устанвлена пределами координат текущей группы пикселей
cntr=1;     %Начальная инициализация счётчика найденных прямых
AzStepN=90; %72Количество шагов по азимуту
    tic;
for nAz=1:AzStepN
  Az=-pi+nAz/AzStepN*2*pi;    %Приведение опорного азимута к радианам
  AzDiffLayer=abs(DiffAngle-Az);  %Слой раницы угла градиента и текущего азимута
  AzDiffLayer(AzDiffLayer>pi)=abs(AzDiffLayer(AzDiffLayer>pi)-2*pi); %Разница от 0 до 180 градусов
  AzDiffLayer=AngDynDiff./(AzDiffLayer./(AzOwrleapFactor*2*pi/AzStepN));	%Взятие обратной величины и домножение на слой динамического фильтра. 2*pi/AzStepN количество радиан в одном шаге. Если разница азимута в пикселе и текущего азимута меньше чем Коэффициент_перекрытия*Сектор одного шага, то AngDynDiff становится ярче, если больше, то тусклее
  %BitMsk=AzDiffLayer>BnrztnLvl;   %Бинаризация            %(2*360/AzStepN) - разница в 2*х (х=градусов в одном шаге) соответствует показателю в 1 на обратном слое разницы
  BitMsk=(AzDiffLayer>BnrztnLvl)&MainBitMsk;   %Бинаризация с исключением связных областей (остравков) не попавших в первичную Бти-маску (MainBitMsk) полученную сразу после динамической фильтрации
  %BitMsk=(AzDiffLayer<(AzOwrleapFactor*2*pi/AzStepN))&MainBitMsk; %Бинаризация Только по условию попадения в текущий диапазон азимутов градиента
%           imshow(flip(BitMsk,1)); axis auto xy
%           title(['Angle: ',num2str(fix(Az/pi*1800)/10),char(176)])
%            fileName=['Res/Bmsk2/',num2str(nAz)];
%            saveas(gcf,fileName,'png');
%         mov(nAz) = getframe(gcf);
  [ObjLayer, Num]=bwlabel(BitMsk,8);   %Разбиение пикселей на группы
  ObjLayer=ObjLayer(:);   %Разименование слоя обьектов (растр теперь одномерный)
  ObjLayer=cat(2,ObjLayer,(1:size(ObjLayer,1))'); %добавляем второй столбец - адреса массива
  ObjLayer=sortrows(ObjLayer);    %Сортировка по возростанию номера обьекта (первого столбца)
  %ObjLayer(1:find(ObjLayer(:,1),1)-1,:)=[];   %Удаление пикселей на которыйх не расположен обьект (нулевой номер)
    i=find(ObjLayer(:,1),1);   %Итератор пикселей (адрес первого не нулевого элемента)
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
            %GrpAzmth=atan2(CartesVectSumY,CartesVectSumX); %Азимут суммарного вектора (средневзвешенный азимут группы)
            GrpAzmth=atan2(CartesVectSumY,CartesVectSumX);
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
            if (sqrt((Coord(cntr,1)-Coord(cntr,2))^2+(Coord(cntr,3)-Coord(cntr,4))^2))>StrMinLngth    %Проверка длинны отрезка (слишком короткие отбрасываются)
              AddData(cntr,:)=[CPR,CPC,GrpAzmth];    %Запись дополнительных данных (строка+столбец центральной точки отрезка (не целые), средневзыешенный азимут градиента)
              TtlGrVal=sqrt(CartesVectSumX^2+CartesVectSumY^2); %Величина градиента группы (векторной суммы градиентов от каждого пикселя)
              Rlablty(cntr,:)=[TtlGrVal,AdrB-AdrA+1,TtlGrVal/TtlPixGrVal];  %Параметр достоверности по каждой группе (отрезку) - Величина градинта группы, оличество пикселей, косинус расброса углов
              cntr=cntr+1;  %Инкремент счётчика для следующей прямой
            end                          
         end   
%                 if (AdrB-AdrA>100); %(nAz==5)&& %Блок для отрисовки отрезка вписанного в островок
%                     %imshow(flip(BitMsk,1)); axis auto xy
%                     imshow(BitMsk); axis auto xy
%                     hold on
%                     plot([LftFr RghtFr RghtFr LftFr LftFr], [UpFr UpFr DwnFr DwnFr UpFr],'green');  %Фрейм
%                     plot(Coord(cntr-1,1:2),Coord(cntr-1,3:4),'magenta');    %Отрезок
%                     plot ([FrdgePoints(:,1);FrdgePoints(1,1)],[FrdgePoints(:,2);FrdgePoints(1,2)],'rx');    %точки пересечения прямой с фреймом
%                     plot(CPC,CPR,'ro'); %Средневзвешенный центр отрезка
%                     hold off
%                 end
    end
end
%     v = VideoWriter('Res/BitMask','Archival');
%     open(v);
%     writeVideo(v,mov);
%     close(v);
cntr=cntr-1;    %Отмена последнего инкремента счётчика
  Coord(cntr+1:end,:)=[];   %Удаление лишней нициализированной области выходных массивов.
  AddData(cntr+1:end,:)=[]; % И одновременное удаление последней линии (кроме случая когда последняя линия имеет подходящую длину)
  Rlablty(cntr+1:end,:)=[];
    disp('формирование отрезков:')
    toc;
       imagesc(Data);axis equal
        hold on
        Nlevels=1000;
        Hcolor=hot(Nlevels);
        for i=1:size(Coord,1)   %Проверочный блок отрисовки прямых
          color=Rlablty(i,3);
          %color=Rlablty(i,3)*Nlevels;
          %plot(Coord(i,1:2),Coord(i,3:4),'Color',Hcolor(round(color),:));
          plot(Coord(i,1:2),Coord(i,3:4),'Color',[color color color]);          
        end
        hold off
%******Блок Прореживания и сортировки списка прямых---
StrLngth=sqrt((Coord(:,2)-Coord(:,1)).^2+(Coord(:,4)-Coord(:,3)).^2);   %Длины прямых
%NrmlzdGrad=Rlablty(:,1)./Rlablty(:,2);  %Величина суммарного попиксельного градиента отнесённая к количеству пикселей гаждой группы (каждой прямой)
    %TrustInd=sqrt((StrLngth./mean(StrLngth,1)).^2+(NrmlzdGrad./mean(NrmlzdGrad,1)).^2+(Rlablty(:,3)./mean(Rlablty(:,3),1)).^2); %Индекс значимости прямой - эвклидова метрика в пространстве Длина-Величина_градиента-Разброс_азимутов для группы пикселей по каждой прямой
TrustInd=Rlablty(:,3);  %Индекс значимости (обратнопропорционален разбросу азимутов в группе пикселей)
    %TrustInd=sqrt((Rlablty(:,3)./mean(Rlablty(:,3))).^2+(NrmlzdGrad./mean(NrmlzdGrad)).^2);
    %TrustInd=sqrt(6.*(Rlablty(:,3)./mean(Rlablty(:,3),1)).^2+(Rlablty(:,2)./mean(Rlablty(:,2),1)).^2);
StrArray=cat(2,Coord,AddData,Rlablty,TrustInd,StrLngth); %Обьединение всех параметров найденных прямых в один массив (по строке для каждой прямой)
%StrArray(:,1:6)=StrArray(:,1:6)+SafeCoreFrame;  %Учёт (запасной) рамки растра от свёртки 
StrArray=sortrows(StrArray,-11);    %Сортировка прямых по значимости
StrArray(:,11)=1:cntr;      %Замена индекса занчимости на Индентификатор прямой (чем меньше ID, тем выше значимость)
    clearvars Coord AddData Rlablty;
%Обединение прямых в группы
ToStrDltMsk=false(cntr,1);
AngTreshold=1.8*(2*pi/AzStepN); %1.8 %Норма разницы по азимуту для группы прямых (+-AngTreshold радиан от значения азимута текущей прямой (текущий ID))
%Нормы разницы центров прямых приведены как коэффициенты для длины инициирующей прямой (прямой с текущим ID)
PrllelKTrshold=1/6;	%1\8 1/4Норма разницы координат центра отрезка по оси параллельной прямой с текущим ID (в долях длины текущей прямой)
AcrossKTrshold=0.5*StrMinLngth;	%Предел разницы координат по оси перпендикулярной текущей прямой (в пикселях)
    LngthTrshold=0.4;   %0.5 Предел разницы длины текущей и проверяемых прямых (в долях длины текущей прямой) - много-большие и многоменьшие - не должны удалятся
    tic;
GK=tan(StrArray(:,7));    %Коэффициент прямой градиента (перпендикуляра) для уравнения вида y=kx
CrntIDStrK=tan(StrArray(:,7)-pi/2);%Коэффициент текущей (иницииирующей) прямой для уравнения вида y=kx
SameSourseStrMterick=zeros(cntr,1);                     %Инициализация массива метрик связности прямых
for ID=1:cntr-1
  if ToStrDltMsk(ID)==false
      SameSourseStrMterick(ID+1:end)=abs(StrArray(ID+1:end,7)-StrArray(ID,7)); %Определение разниц азимутов инициирующей прямой и всех остальных
      SameSourseStrMterick(SameSourseStrMterick>=pi)=abs(SameSourseStrMterick(SameSourseStrMterick>=pi)-2*pi);%Приведение азимутальной разницы к пределам от 0 до 180 градусов
      %Координаты точки пересечения линий градиента с текущей прямой
      x=(CrntIDStrK(ID)*StrArray(ID,6)-GK(ID+1:end).*StrArray(ID+1:end,6)+StrArray(ID+1:end,5)-StrArray(ID,5))./(CrntIDStrK(ID)-GK(ID+1:end));
      y=CrntIDStrK(ID).*(x-StrArray(ID,6))+StrArray(ID,5); %Координаты точек пересечения линий градиента от каждой прямой и текущей прямой
      PrllelDist=sqrt((x-StrArray(ID,6)).^2+(y-StrArray(ID,5)).^2)./StrArray(ID,12);	%Расстояние (в долях длины текущей прямой) между центральными точками вдоль оси параллельной инициирующей прямой (расстояние от точки пересечения линии градиента с текущей прямой до центральной точки текущей(инициирующей) прямой)
      AcrossDist=sqrt((x-StrArray(ID+1:end,6)).^2+(y-StrArray(ID+1:end,5)).^2);    %Аналогично расстояние (в пикселях) между прямыми (каждой и текущей) вдоль линии градиента (перпендикуляра каждой линии).
      LngthDiff=abs(StrArray(ID,12)-StrArray(ID+1:end,12))/StrArray(ID,12);      %Разница длинн текущей и проверяемых прямых в долях текущей прямой
      %SameSourseStrMterick(ID+1:end)=sqrt((SameSourseStrMterick(ID+1:end)./AngTreshold).^2+(PrllelDist./PrllelKTrshold).^2+(AcrossDist./AcrossKTrshold).^2);  %Мера сопряжённости прямых (на сколько для своего построения они используют одинаковые с инициирующей прямой пиксели) - эвклидова метрика в (линейнозависимом) пространстве Разница_азимутов - Параллельное_расстояние - Градиентное_расстояние. Все измерения пронормированы
      SameSourseStrMterick(ID+1:end)=sqrt((SameSourseStrMterick(ID+1:end)./AngTreshold).^2+(PrllelDist./PrllelKTrshold).^2+(AcrossDist./AcrossKTrshold).^2+(LngthDiff/LngthTrshold).^2);  %Мера сопряжённости прямых (на сколько для своего построения они используют одинаковые с инициирующей прямой пиксели) - эвклидова метрика в (линейнозависимом) пространстве Разница_азимутов - Параллельное_расстояние - Градиентное_расстояние. Все измерения пронормированы
      %SameSourseStrMterick(ID+1:end)=sqrt((SameSourseStrMterick(ID+1:end)./AngTreshold).^2+(PrllelDist./PrllelKTrshold).^2+(AcrossDist./AcrossKTrshold(ID+1:end)).^2);  %Мера сопряжённости прямых (на сколько для для своего построения они используют одинаковые с инициирующей прямой пиксели) - эвклидова метрика в (линейнозависимом) пространстве Разница_азимутов - Параллельное_расстояние - Градиентное_расстояние. Все измерения пронормированы
      ToStrDltMsk(ID+1:end)=ToStrDltMsk(ID+1:end)|SameSourseStrMterick(ID+1:end)<2; %2=sqrt(4) - для четырёхмерного пространства
  end
end
StrArray(ToStrDltMsk,:)=[];    %Удаление вторичных связанных прямых
  % Coord=StrArray(:,1:4);  %Выгрузка искомых массивов параметров прямых в порядке отсортированном по значимости прямой
  % AddData=StrArray(:,5:7);
  % Rlablty=StrArray(:,8:10);
  disp('прореживание и сортировка:')
    toc;
 %imagesc(Data);
imshow(uint16(Data));
 hold on
  axis image
 for i=1:size(StrArray,1)
%     %if StrArray(i,2)<=800
    color=StrArray(i,10);
    plot(StrArray(i,1:2),StrArray(i,3:4),'Color',[0,color,0])
    %plot(StrArray(i,1:2),StrArray(i,3:4),'Color',[color,color,color]);
%     %end
 end
% hold off


end

