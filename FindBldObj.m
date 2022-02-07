function [ ObjCoord, ObjRlablty ] = FindBldObj( StrArray ,Raster, SrsCoord)
%FINDBLDOBJ Определяет прямоугольные обьекты из полученного массива отрезков
%ObjCoord ([Row Col] в каждой строке) - координаты выявленного обьекта
%ObjRlablty - в каждой строке - относитльная достоверность обьекта (номер строки соответствует строке с координатами)
RghtAngleTrshold=15/180*pi; %12 Прямым будет считаться угол 90 градусов +-RghtAngleTrshold радиан
CornMinDist=9; %9 8Максимальное расстояние между краями отрезков, которые будут составлять прямой угол и будут пригодны для определения их как фрагмента искомого обьекта
ComplStrMinDist=8;  %10 8 Максимальное расстояние между ближними концами отрезков составной прямой
StrMaxLngth=150;%Максимальная длина прямой, которая может быть частью строения
difWCoef=0.5;	%0.2В каждой группе выбираются отрезки из тех, что участвуют в образовании прямого угла только одним концом (второй повис в стороне). Если средний (пиксельный) градиент меньше среднего (пиксельного) градиента по обьекту на difWCoef от градиента обьекта, то такая прямая отбрасывается
TypObjMaxCorn=8;        %Обьектом будет считаться совокупность связанных прямых с количеством углов более 1-го (и соотв. более 3-х прямых). При этом типовымми будут считаться обьекты с TypObjMaxCorn количеством углов (остальные считаются сложными\крупными).
RAngleAccTH=10/180*pi;   %15 10 Обьект будет отобран как достоверный, если средневзвешенный "прямой" угол (между прямыми обьекта) отличается на RAngleAcc от 90 градусов (или Treshold-значение)
minW_TH=1./(1.5e+07);     %2.2 e+07Обьект будет отобран как достоверный, если средний  не пронормированных вес его прямых будет не меньше чем minW_TH (или Treshold-значение)
TypclQStr_TH=1/2.0; %2.5 При отборе достоверных обьектов в приоритете обьекты с большим количеством отрезков. TresHold-значение - 1/типовое количество прямых
%     LowBrght_TH=4040;   %4040 Пороговые (treshold) значения яркости обьектов в токах их геометрических центров (берётся окрестность из 4 пикселей) - Предел условно "тёмных" обьектов  в цетовой шкале от 0 до 65535
%     HghBrght_TH=8700;   %8700 Аналогично пороговые значения условно светлых обьектов: 65535-HghBrght_TH
AzCosDev_TH=0.58*mean(1./StrArray(:,10));   %0.58
%Ядро динамического фильтра было размером 15*15
%Ядро собеля с размытием - 5*5
global StrBitMsk RequestRowMsk RequestColMsk;
StrArray(:,7)=StrArray(:,7)-pi/2;%Преобразование углов градиента в углы прямых (в диапазоне от -270 до +90 градусов в радианах)
StrArray(StrArray(:,7)<-pi/2,7)=StrArray(StrArray(:,7)<-pi/2,7)+pi; %Преобразование азимутов оётрезков к пределам от -90 до +90 градусов (в радианах)
StrArray=sortrows(StrArray,7);  %Сортировка отрезков по их азимуту
%Далее все прямые с отрицательным азимутом будут сопоставлятся со всеми прямым с положительным азимутом для выявления прямых углов (с допуском)
AzFrdgeM=find(StrArray(:,7)>-RghtAngleTrshold,1);	%Адрес первой прямой в списке, азимут котрой больше 0-порог градусов - начальный адрес списка Вызывающих (Master) прямых (полож.азимут, k>0,вторая половина)
AzFrdgeS=find(StrArray(:,7)>RghtAngleTrshold,1);	%Адрес первой прямой в списке, азимут котрой больше 0+порог градусов - конечный адрес списка Откликающихся (Slave) прямых (отриц.азимут, k<0, первая половина)
RghtAngleObjCandidates=zeros(size(StrArray,1)-AzFrdgeM+1,AzFrdgeS-1); %Инициализация массива Master-азимутов
StrBitMsk=false(size(RghtAngleObjCandidates));    %Инициализация маски отбора пар отрезков, которые составляют прямой угол обьекта
    tic;
for i=1:AzFrdgeS-1      %Копируем Master-азимуты по количеству Slave-прямых раз
  RghtAngleObjCandidates(:,i)=StrArray(AzFrdgeM:end,7);  % Адрес строки - порядковый номер Master-отрезка (каждый столбец - все Master отрезки)
end
for i=1:size(StrArray,1)-AzFrdgeM+1; %Копируем Slave-азимуты Master- раз
  RghtAngleObjCandidates(i,:)=RghtAngleObjCandidates(i,:)-StrArray(1:AzFrdgeS-1,7)';    % Адрес столбца - порядковый номер Slave-отрезка (каждая сторка все Slave-ы)
end
RghtAngleObjCandidates=abs(RghtAngleObjCandidates-pi/2);  %Величина отличия разицы азимутов каждого Master с каждым Slave-отрезком от 90градусов (в радианах)
StrBitMsk=RghtAngleObjCandidates<=RghtAngleTrshold;	%Маска Пар прямых - кандидатов в обьекты (по критерию наличия прямого угла с допуском RghtAngleTrshold)
MstrStrData=StrArray(AzFrdgeM:end,:);   %Копирование соответствующих параметров отрезков в список Master-прямых
SlveStrData=StrArray(1:AzFrdgeS-1,:);   %Копирование параметров  в список Slave-прямых
% MstrGrW=MstrStrData(:,8)./MstrStrData(:,9);     %Веса=Суммградиент/кол-во пикселей каждого отрезка
% SlveGrW=SlveStrData(:,8)./SlveStrData(:,9);
MstrGrW=MstrStrData(:,8)./MstrStrData(:,12);     %Веса=Суммградиент/кол-во пикселей каждого отрезка
SlveGrW=SlveStrData(:,8)./SlveStrData(:,12);
MstrK=tan(MstrStrData(:,7));    %Коэффициенты для каждого отрезка для ур-я y=kx+b
SlveK=tan(SlveStrData(:,7));
MSStrSolvesX=zeros(1,AzFrdgeS-1);   %Инициализация массива точек пересечения всех Slave- с текущей Master-прямой (Х, колонка)
MSStrSolvesY=zeros(1,AzFrdgeS-1);   % аналогично для (Y, стобец)
MDist=zeros(2,AzFrdgeS-1);          %Инициализация расстояний от точки пересечения до краёв Master-отрезка
SDist=MDist;                        % аналогично для Slave-отрезка
ShrtstMDist=zeros(1,AzFrdgeS-1);    %Инициализация расстояний от каждой точки пересечения до ближайшего из двух краёв Master-отрезка
ShrtstSDist=ShrtstMDist;            % аналогично для ближайшего края соответствующего Slave-отрезка
        %imshow(Raster(1:365+60,1:2412+60));
        figure; 
            TbleSize=0;
            %TbleSize=TbleSize*size(Raster,1);
            %Tble=zeros(TbleSize,size(Raster,2));
        %imshow(cat(1,Tble,Raster,Tble));
        imshow(Raster);
        axis on,axis equal,hold on;
        %SrsCoord=[11025 31;2420 370;3444 66;11388 352]; %Координаты обучающей выборки для 5-го изображения
    tic
StartRow=find(any(StrBitMsk,2),1);
FinshRow=size(StrArray,1)-AzFrdgeM+1;
    clearvars StrArray RghtAngleObjCandidates
for i=StartRow:FinshRow %Перебор Master-отрезков и поиск по каждому подходящих slave-отрезков
    MSStrSolvesX(StrBitMsk(i,:))=(MstrK(i)*MstrStrData(i,6)-MstrStrData(i,5)-SlveK(StrBitMsk(i,:)).*SlveStrData(StrBitMsk(i,:),6)+SlveStrData(StrBitMsk(i,:),5))./(MstrK(i)-SlveK(StrBitMsk(i,:)));
    MSStrSolvesY(StrBitMsk(i,:))=MstrK(i).*(MSStrSolvesX(StrBitMsk(i,:))-MstrStrData(i,6))+MstrStrData(i,5);    %Точки пересечения текущей Master-прямой и всех Slave-прямых (X и Y - соответственно)
    MDist(1,StrBitMsk(i,:))=sqrt((MSStrSolvesX(StrBitMsk(i,:))-MstrStrData(i,1)).^2+(MSStrSolvesY(StrBitMsk(i,:))-MstrStrData(i,3)).^2);    %Расстояния от каждой точки пересечения до каждого из двух концов Master-отрезка
    MDist(2,StrBitMsk(i,:))=sqrt((MSStrSolvesX(StrBitMsk(i,:))-MstrStrData(i,2)).^2+(MSStrSolvesY(StrBitMsk(i,:))-MstrStrData(i,4)).^2);    
    ShrtstMDist(StrBitMsk(i,:))=min(MDist(:,StrBitMsk(i,:)),[],1);    %Расстояния от каждой точки пересечения до БЛижайшего конца Master-отрезка
    OnLineMsk=(MSStrSolvesX(StrBitMsk(i,:))<max(MstrStrData(i,1:2),[],2))&(MSStrSolvesX(StrBitMsk(i,:))>min(MstrStrData(i,1:2),[],2));  %Маска ShrtstSDist(StrBitMsk(i,:)-х элеметов точки пересечения которых лежат на соответствующем Master-отрезке. Длина этой маски равна количеству единиц в текущей строке StrBitMsk
    ShrtstMDist(StrBitMsk(i,:))=ShrtstMDist(StrBitMsk(i,:)).*~OnLineMsk;    %%Если точка пересечения лежит в пределах отрезка - расстояния заменяются на 0 (пределы отрезка взты по X - номеру столбца)
    SDist(1,StrBitMsk(i,:))=sqrt((MSStrSolvesX(StrBitMsk(i,:))-SlveStrData(StrBitMsk(i,:),1)').^2+(MSStrSolvesY(StrBitMsk(i,:))-SlveStrData(StrBitMsk(i,:),3)').^2);% аналогично расстояния от каждой точки пересечения до каждого из двух концов соответствующего Slave-отрезка.
    SDist(2,StrBitMsk(i,:))=sqrt((MSStrSolvesX(StrBitMsk(i,:))-SlveStrData(StrBitMsk(i,:),2)').^2+(MSStrSolvesY(StrBitMsk(i,:))-SlveStrData(StrBitMsk(i,:),4)').^2);
    ShrtstSDist(StrBitMsk(i,:))=min(SDist(:,StrBitMsk(i,:)),[],1);    %Расстояния от каждой точки пересечения до блИЖайшего конца соответствующего Slave-отрезка
    OnLineMsk=(MSStrSolvesX(StrBitMsk(i,:))<max(SlveStrData(StrBitMsk(i,:),1:2),[],2)')&(MSStrSolvesX(StrBitMsk(i,:))>min(SlveStrData(StrBitMsk(i,:),1:2),[],2)');  %Маска ShrtstSDist(StrBitMsk(i,:)-х элеметов точки пересечения которых лежат на соответствующем Slave-отрезке. Длина этой маски равна количеству единиц в текущей строке StrBitMsk
    ShrtstSDist(StrBitMsk(i,:))=ShrtstSDist(StrBitMsk(i,:)).*~OnLineMsk;    %Если точка пересечения лежит в пределах отрезка - расстояния заменяются на 0 (пределы отрезка взты по X - номеру столбца)
    StrBitMsk(i,StrBitMsk(i,:))=(ShrtstMDist(StrBitMsk(i,:))<CornMinDist & ShrtstSDist(StrBitMsk(i,:))<CornMinDist);%Удаление из маски пар прямых (Master-Slave) точка пересечения которых слишком далеко от отрезков
end
    toc;
%Формируем связи для составных отрезков.
global MstrStrLinkList SlveStrLinkList;
MstrStrLinkList=zeros(5000,2);  %Инициализация массива связей Master-отрезков
SlveStrLinkList=zeros(5000,2);  % и Slave-отрезков соответственно
StrCntr=1;  %Итератор текущего отрезка
StrLinkCoord=zeros(size(MstrStrData,1),4);  %Инициализация массива расстояний от каждой точки вызывающей прямой до каждой точки связанной прямой
    tic
%Поиск связей среди Master-отрезков:
for i=StartRow:FinshRow     %Итериции для Master-отрезков. В первом приближении из списка Master-отрезков выбираются близкие по азимуту (в пределах RghtAngleTrshold)
AdrA=i+1;   %Начальный адрес диапазона отрезков в MstrStrData (первого приближения). Для того чтобы каждая связь была обнаружена единожды (отрезки не ссылались друг на друга), текущий отрезок сравнивается с другими отрезками, которые с бОльшим азимутом.
%AdrA=find(MstrStrData(:,7)>MstrStrData(i,7)-RghtAngleTrshold,1);
AdrB=find(MstrStrData(:,7)>MstrStrData(i,7)+RghtAngleTrshold,1);    %Следующий после конечного адреса текущего диапазона отрезков в MstrStrData (первого приближения).
if ((AdrB-AdrA)>0)  %Исключение обработки пустого диапазона
  AdrB=AdrB-1;      %AdrB - теперь адрес последнего элемента диапазона
  StrLinkCoord(AdrA:AdrB,1)=sqrt((MstrStrData(i,1)-MstrStrData(AdrA:AdrB,1)).^2+(MstrStrData(i,3)-MstrStrData(AdrA:AdrB,3)).^2);%Расстояние от первой точки текущего отрезка до первой точки каждого отрезка
  StrLinkCoord(AdrA:AdrB,2)=sqrt((MstrStrData(i,2)-MstrStrData(AdrA:AdrB,2)).^2+(MstrStrData(i,4)-MstrStrData(AdrA:AdrB,4)).^2);%Расстояние от второй точки текущего отрезка до второй точки каждого отрезка
  StrLinkCoord(AdrA:AdrB,3)=sqrt((MstrStrData(i,1)-MstrStrData(AdrA:AdrB,2)).^2+(MstrStrData(i,3)-MstrStrData(AdrA:AdrB,4)).^2);%Расстояние от первой точки текущего отрезка до второй точки каждого отрезка
  StrLinkCoord(AdrA:AdrB,4)=sqrt((MstrStrData(i,2)-MstrStrData(AdrA:AdrB,1)).^2+(MstrStrData(i,4)-MstrStrData(AdrA:AdrB,3)).^2);%Расстояние от второй точки текущего отрезка до первой точки каждого отрезка
  ResAdr=AdrA-1+find((min(StrLinkCoord(AdrA:AdrB,:),[],2)<ComplStrMinDist)&(max(StrLinkCoord(AdrA:AdrB,:),[],2)>(MstrStrData(i,12)+MstrStrData(AdrA:AdrB,12)-ComplStrMinDist)));%Второе приближение: выбор отрезка(ов) один конец которых совпадает с точностью до CornMinDist, а другие концы находятся на расстоянии друг от друга равном сумме длинн этих отрезков с точностью до CornMinDist
  if ~isempty(ResAdr)   %Регистрация составных отрезков (если такие нашлись на текущей итерации)
    MstrStrLinkList(StrCntr:StrCntr+size(ResAdr,1)-1,1)=i;  %Номер вызывающего отрезка в MstrStrData
    MstrStrLinkList(StrCntr:StrCntr+size(ResAdr,1)-1,2)=ResAdr; %Номер связанного отрезка в MstrStrData
    StrCntr=StrCntr+size(ResAdr,1); %обновление счётчика (в положение первого адреса записи диапазона для следующей итерации)
  end
end
end
MstrStrLinkList(StrCntr:end,:)=[];  %Удаление оставшейся (запасённой) части массива номеров связанных Master-Отрезков
%         for i=1:size(MstrStrLinkList,1) %Отрисовка составных отрезков
%             plot(MstrStrData(MstrStrLinkList(i,1),1:2),MstrStrData(MstrStrLinkList(i,1),3:4),'color','blue');
%             plot(MstrStrData(MstrStrLinkList(i,2),1:2),MstrStrData(MstrStrLinkList(i,2),3:4),'color','red');
%         end
%Поиск связей среди Slave-отрезков:  
StrCntr=1;  %Сброс итератора отрезка
StrLinkCoord=zeros(size(SlveStrData,1),4);  %Переинициализация массива расстояний (по количество Slave-отрезков)
for i=1:size(StrBitMsk,2)	%Итериции для Slave-отрезков. В первом приближении аналогично выбираются близкие по азимуту отрезки
AdrA=i+1;   %Начальный адрес диапазона отрезков в SlveStrData (первого приближения).
%AdrA=find(SlveStrData(:,7)>SlveStrData(i,7)-RghtAngleTrshold,1);
AdrB=find(SlveStrData(:,7)>SlveStrData(i,7)+RghtAngleTrshold,1);%Следующий после конечного адреса текущего диапазона отрезков в SlveStrData (первого приближения).
if ((AdrB-AdrA)>0)  %Исключение обработки пустого диапазона
  AdrB=AdrB-1;      %AdrB - теперь адрес последнего элемента диапазона
  StrLinkCoord(AdrA:AdrB,1)=sqrt((SlveStrData(i,1)-SlveStrData(AdrA:AdrB,1)).^2+(SlveStrData(i,3)-SlveStrData(AdrA:AdrB,3)).^2);%Расстояние от первой точки текущего отрезка до первой точки каждого отрезка
  StrLinkCoord(AdrA:AdrB,2)=sqrt((SlveStrData(i,2)-SlveStrData(AdrA:AdrB,2)).^2+(SlveStrData(i,4)-SlveStrData(AdrA:AdrB,4)).^2);%Расстояние от второй точки текущего отрезка до второй точки каждого отрезка
  StrLinkCoord(AdrA:AdrB,3)=sqrt((SlveStrData(i,1)-SlveStrData(AdrA:AdrB,2)).^2+(SlveStrData(i,3)-SlveStrData(AdrA:AdrB,4)).^2);%Расстояние от первой точки текущего отрезка до второй точки каждого отрезка
  StrLinkCoord(AdrA:AdrB,4)=sqrt((SlveStrData(i,2)-SlveStrData(AdrA:AdrB,1)).^2+(SlveStrData(i,4)-SlveStrData(AdrA:AdrB,3)).^2);%Расстояние от второй точки текущего отрезка до первой точки каждого отрезка
  ResAdr=AdrA-1+find((min(StrLinkCoord(AdrA:AdrB,:),[],2)<ComplStrMinDist)&(max(StrLinkCoord(AdrA:AdrB,:),[],2)>(SlveStrData(i,12)+SlveStrData(AdrA:AdrB,12)-ComplStrMinDist)));%Второе приближение: выбор отрезка(ов) один конец которых совпадает с точностью до CornMinDist, а другие концы находятся на расстоянии друг от друга равном сумме длинн этих отрезков с точностью до CornMinDist
  if ~isempty(ResAdr)   %Регистрация составных отрезков (если такие нашлись на текущей итерации)
    SlveStrLinkList(StrCntr:StrCntr+size(ResAdr,1)-1,1)=i;  %Номер вызывающего отрезка в MstrStrData
    SlveStrLinkList(StrCntr:StrCntr+size(ResAdr,1)-1,2)=ResAdr; %Номер связанного отрезка в MstrStrData
    StrCntr=StrCntr+size(ResAdr,1); %обновление счётчика (в положение первого адреса записи диапазона для следующей итерации)
  end
end
end
SlveStrLinkList(StrCntr:end,:)=[];  %Удаление оставшейся (запасённой) части массива номеров связанных Master-Отрезков
toc;
        for i=1:size(SlveStrLinkList,1) %Отрисовка составных отрезков
            plot(SlveStrData(SlveStrLinkList(i,1),1:2),SlveStrData(SlveStrLinkList(i,1),3:4),'color','yellow');
            plot(SlveStrData(SlveStrLinkList(i,2),1:2),SlveStrData(SlveStrLinkList(i,2),3:4),'color','white');
        end
        StrCntr=0;
%Формируем "обьектные" связи между парами связанных прямых. Если у двух пар есть общая прямая - они связаны.
% То есть все точки (true-биты) в столбце StrBitMsk- один обьект (то же и для строк), и.т.д - рекурсивно  
RequestRowMsk=any(StrBitMsk,2); %В первом приближении, маска запроса для строк содержит все не нулевые строки
RequestColMsk=any(StrBitMsk,1); %Аналогично - все не нулевые стобцы
TtlStrQ=sum(sum(StrBitMsk));    %Количество нейденых прямых углов на снимке (необходимо для определения размера памяти, которая потребуется)
Cntrs=zeros(TypObjMaxCorn,1);   %Счётчик количеств обьектов в каждой группе
ObjArr=cell(fix(TtlStrQ/2),TypObjMaxCorn);  %Cell-Массив обьектов для 2-х, 3-х, 4-х и более точечных обьектов соответственно. (точка это угол с парой прямых)
    NoDltMstrBitMsk=false(size(StrBitMsk,1),1); %Массив запрета на удаление: В цикле есть процедура проверки градиента каждой прямой и сопоставления его со средневзвешенным градиентом по группе (обьекту), проверке подлежат только краевые ("повисшие одним концом") отрезки, которые участвуют в образовании только одного прямого угла в обьекте.
    NoDltSlveBitMsk=false(size(StrBitMsk,2),1); % Такие прямые выявляются таким образом, что они встречаются в списке отрезков по обьекту только один раз. Отрезки образующие составную прямую будут встречаться тоже только один раз, поэтому эти две битмаски запрещают проверку и удаление компонентов составных прямых
        tic;
            DltCCntr=0;
for i=find(RequestRowMsk,1):FinshRow %Перебор Master-отрезков и поиск всех связанных отрезков (из тех, что образуют с к.л другим прямой угол)
  if RequestRowMsk(i)==true     %Проверка является ли текущая строка StrBitMsk строкой запроса (т.е. в ней есть точки и она ещё не была проверена)
    [Q,Coords]=FndColByColLinks(i); %Собтвенно поиск сопряженных прямых (точек в StrBitMsk, которые лежат в одной строке ии столбце с точками текущей строки рекурсивно)
    %Обнаружение составных отрезков:
    for s=1:size(Coords,1)  %Перебор пар Master- и  Slave-отрезков
      %StrLinks=find(MstrStrLinkList(:,1)==Coords(s,1)); %поиск связей по текущему Мaster-отрезку
      StrLinks=RecuFindM(Coords(s,1));
      for t=1:size(StrLinks,1)  %Перебор найденных связей по текущему М-отрезку
        if(RequestRowMsk(MstrStrLinkList(StrLinks(t),2))==true)%Выявление связанных отрезков, которые ещё не были проверены на предмет образования обьектов.
            [tQ,tCoords]=FndColByColLinks(MstrStrLinkList(StrLinks(t),2)); %Поиск сопряженных прямых с обнаруженной частью составного Master-отрезка
            Q=Q+tQ;     %Слияние обьектов: Обновление счётчика отрезков в обьекте 
            Coords=cat(1,Coords,tCoords);   %Сведение номеров отрезков двух обьектов в один массив.
            NoDltMstrBitMsk(MstrStrLinkList(StrLinks(t),:))=true;%Запрет на удаление отрезков образующих составной отрезок в нижеследующем блоке выявления слабых градиентов (хоть номер каждого из них и встречается в списке единожды)
                StrCntr=StrCntr+1;
        end
      end
      %StrLinks=find(SlveStrLinkList(:,1)==Coords(s,2)); %поиск связей по текущему Slave-отрезку
      StrLinks=RecuFindS(Coords(s,2));
      for t=1:size(StrLinks,1)  %Перебор найденных связей по текущему S-отрезку
        if(RequestColMsk(SlveStrLinkList(StrLinks(t),2))==true)%Выявление связанных отрезков, которые ещё не были проверены на предмет образования обьектов.
            [tQ,tCoords]=FndRowByRowLinks(SlveStrLinkList(StrLinks(t),2)); %Поиск сопряженных прямых с обнаруженной частью составного Master-отрезка
            Q=Q+tQ;     %Слияние обьектов: Обновление счётчика отрезков в обьекте 
            Coords=cat(1,Coords,tCoords);   %Сведение номеров отрезков двух обьектов в один массив.
            NoDltSlveBitMsk(SlveStrLinkList(StrLinks(t),:))=true;%Запрет на удаление отрезков образующих составной отрезок в нижеследующем блоке выявления слабых градиентов (хоть номер каждого из них и встречается в списке единожды)
                StrCntr=StrCntr+1;
        end
      end
      
    end

    if Q>1  %Отсеить прямые средний градиент которых меньше остальных на difWCoef (если эта прямая образует прямой угол только одним концом). Обьект будет обрабатываться если имеет хотябы два угла
      CDltMsk=false(size(Coords,1),1);  %Инициализаци маски удалния пар прямых
      UMC=unique(Coords(:,1));    %Unique Master Coorinates - адреса Master-прямых обьекта без повторений
      USC=unique(Coords(:,2));    %Unique Slave Coorinates - адреса Slave-прямых обьекта без повторений
      MeanObjGrW=(sum(MstrStrData(UMC,8))+sum(SlveStrData(USC,8)))./(sum(MstrStrData(UMC,12))+sum(SlveStrData(USC,12)));%Вычисление средневзвешенного веса по группе (сумма градиентов/сумму пикселей)
      for t=1:size(UMC,1) 	%Перебор Master-отрезков
        if (~NoDltMstrBitMsk(UMC(t)))&&((MeanObjGrW-MstrGrW(UMC(t)))>(difWCoef*MeanObjGrW))&&(sum(sum(find(Coords(:,1)==UMC(t))))==1)   %Условие прохождения по недостатку веса & условие участия прямой в образовании только одного угла
          CDltMsk(t)=true;  %Пометка отрезка-кандидата на удаление из списка отрезков обьекта
            DltCCntr=DltCCntr+1;
            plot(MstrStrData(UMC(t),1:2),MstrStrData(UMC(t),3:4)+TbleSize,'Color','cyan');
        end
      end
      for t=1:size(USC,1)   %Перебор Slave-отрезков
        if (~NoDltSlveBitMsk(USC(t)))&&((MeanObjGrW-SlveGrW(USC(t)))>(difWCoef*MeanObjGrW))&&(sum(sum(find(Coords(:,2)==USC(t))))==1)
          CDltMsk(t)=true;  %Пометка отрезка-кандидата на удаление из списка отрезков обьекта
            DltCCntr=DltCCntr+1;
            plot(SlveStrData(USC(t),1:2),SlveStrData(USC(t),3:4)+TbleSize,'Color','cyan');
        end
      end
      Coords(CDltMsk,:)=[]; %Удаление "краевых тусклых" отрезков
      Q=Q-sum(CDltMsk);     %Соответствующий вычет из количества углов
      if Q>1    %Учёт обьекта. Обьектом считается по крайней мере три прямые с двумя (почти) прямоугольными углами
        if Q<=TypObjMaxCorn %Учёт наиболее частовстречаемых обьектов
          Cntrs(Q-1)=Cntrs(Q-1)+1;      %Инкремент соответствующего счётчика
          ObjArr{Cntrs(Q-1),Q-1}=Coords;%Запись координат в массив для дальнейшей обработки
        else                %Учёт сложных (редких) обьектов
          Cntrs(TypObjMaxCorn)=Cntrs(TypObjMaxCorn)+1;      %Инкремент соответствующего (последнего) счётчика
          ObjArr{Cntrs(TypObjMaxCorn),TypObjMaxCorn}=Coords;%Запись координат в массив для дальнейшей обработки
        end
            Mm=unique(Coords(:,1)); %Отрисовка
            Ms=unique(Coords(:,2));
%             for m=1:size(Mm,1)
%              plot(MstrStrData(Mm(m),1:2),MstrStrData(Mm(m),3:4)+TbleSize,'Color','blue');
%             end
%             for s=1:size(Ms,1)
%              plot(SlveStrData(Ms(s),1:2),SlveStrData(Ms(s),3:4)+TbleSize,'Color','blue');
%             end      
      end
	end
  end    
end
        toc;
ObjArr(max(Cntrs)+1:end,:)=[];  %Удаление запасных ячеек массива обьектов
%Цикл отбора наиболее чётких обьектов
%CrntObjRlablty=zeros(size(ObjArr)); %Иницилизация массива качесивенных характеристик обьекта
%ObjAngle=zeros(size(ObjArr));
ObjCoord=zeros(sum(Cntrs),2);   %Инициализация выходных массивов
ObjRlablty=zeros(sum(Cntrs),1);
%Метрики качества бьектов: Средний угол (должен быть ближе к 90 градусов). Средний градиент по каждой прямой должен быть не менее minW_TH
resObjcntr=1;       %resObjcntr-1 - количествj достоверных обьектов (resObjcntr - адрес очередного обьекта, если обьект не достоверный, то в следующей итерации координаты по этому адресу затираются)
for type=1:TypObjMaxCorn-1  %Перебор типовых (простых) обьектов - итерация типов (от 2-х углов до TypObjMaxCorn)
  Angles=zeros(type+1,1);   %(Инициализация) массив Углов при кажой отобранной паре прямых, каждого обьекта
  GradW=zeros(type+1,2);    %(Инициализация) массив Весов каждой прямой (вес - это суммарный градиент (группы пикселей по которым была построена прямая) отнесённый к количеству пикселей (в этой группе))
  for o=1:Cntrs(type)       %Перебор всех обьектов каждого (простого) типа
    Angles=abs(MstrStrData(ObjArr{o,type}(:,1),7)-SlveStrData(ObjArr{o,type}(:,2),7));  %Углы при каждой паре прямых
    GradW(:,1)=MstrGrW(ObjArr{o,type}(:,1));    %Копирование весов в массив весов обьекта
    GradW(:,2)=SlveGrW(ObjArr{o,type}(:,2));
    CrntObjRlablty=sum(unique(GradW(:,1)))+sum(unique(GradW(:,2))); %Сумма весов всех прямых (повторяющиеся исключены)
    GradW=GradW./CrntObjRlablty;    %Пронормированныее веса прямых, взятые по отношению суммарного градиента пикселей по прямой к количеству этих пикселей
    NrmL=min([MstrStrData(ObjArr{o,type}(:,1),12) SlveStrData(ObjArr{o,type}(:,2),12)],[],2);   %Веса прямых по их длинам 
    NrmL=NrmL./sum(NrmL);           %Пронормированны веса (по длинам)
    ObjAngle=sum(abs(Angles-pi/2).*NrmL);   %Средневзвешенный (по длинам коротких отрезков) условно прямой угол по обьекту
    ObjCoord(resObjcntr,1)=sum(unique(MstrStrData(ObjArr{o,type}(:,1),5).*GradW(:,1)))+sum(unique(SlveStrData(ObjArr{o,type}(:,2),5).*GradW(:,2)));   %Вычисление координат центра обьекта
    ObjCoord(resObjcntr,2)=sum(unique(MstrStrData(ObjArr{o,type}(:,1),6).*GradW(:,1)))+sum(unique(SlveStrData(ObjArr{o,type}(:,2),6).*GradW(:,2)));
        %ObjColor=mean(mean(Raster(fix(ObjCoord(resObjcntr,1)):ceil(ObjCoord(resObjcntr,1)),fix(ObjCoord(resObjcntr,2)):ceil(ObjCoord(resObjcntr,2)))));   %Вычисление цвета обьекта (по четырём пикселям) в точке найденного условного геометрического центра
        %NrmlzdObjColorDiff=min([ObjColor/LowBrght_TH (65535-ObjColor)/HghBrght_TH]);	%Пронормированное знаение разницы цвета обьекта и ближайшего порогового значения (0 или 65535)
        NrmL=([MstrStrData(ObjArr{o,type}(:,1),12) SlveStrData(ObjArr{o,type}(:,2),12)])./(sum(unique(MstrStrData(ObjArr{o,type}(:,1),12)))+sum(unique(SlveStrData(ObjArr{o,type}(:,2),12))));  %Нормированные длины (без повторяющихся)
        NrmlzdAzCosDev=sum(unique(NrmL(:,1)./MstrStrData(ObjArr{o,type}(:,1),10)))+sum(unique(NrmL(:,2)./SlveStrData(ObjArr{o,type}(:,2),10)))./(AzCosDev_TH);
    %CrntObjRlablty=sqrt((ObjAngle./RAngleAccTH).^2+(1./CrntObjRlablty./(minW_TH*(type+2))).^2+(1/(type+2)/TypclQStr_TH).^2+NrmlzdObjColorDiff.^2); %Достоверность текущего обьекта
    CrntObjRlablty=sqrt((ObjAngle./RAngleAccTH).^2+(1./CrntObjRlablty./(minW_TH*(type+2))).^2+(1/(type+2)/TypclQStr_TH).^2+NrmlzdAzCosDev.^2); %Достоверность текущего обьекта
    if CrntObjRlablty<sqrt(4)  %ObjAngle(o,type)<RAngleAccTH
      ObjRlablty(resObjcntr)=CrntObjRlablty;    %Копирование параметра достоверности в выходной массив
            plot(ObjCoord(resObjcntr,2),ObjCoord(resObjcntr,1)+TbleSize,'gx')            
            Mm=unique(ObjArr{o,type}(:,1));
            Ms=unique(ObjArr{o,type}(:,2));
            for m=1:size(Mm,1)
             plot(MstrStrData(Mm(m),1:2),MstrStrData(Mm(m),3:4)+TbleSize,'Color','red');
            end
            for s=1:size(Ms,1)
             plot(SlveStrData(Ms(s),1:2),SlveStrData(Ms(s),3:4)+TbleSize,'Color','red');
            end
      resObjcntr=resObjcntr+1;  %Инкремент счётчика отобранных обьектов      
    end
  end
  
    %ObjRlablty=sqrt((ObjAngle(1:Cntrs(type),type)./RAngleAccTH).^2+(ObjRlablty(1:Cntrs(type),type)./(minW_TH*(type+2))).^2);
    %ApprvmnObjtMsk=ObjRlablty<sqrt(2);
end
        %plot(SrsCoord(:,1),SrsCoord(:,2)+TbleSize,'mo');
%Цикл учёта сложных обьектов
type=TypObjMaxCorn; %теперь текущий тип - последний (более TypObjMaxCorn углов)
S=zeros(size(ObjArr,1),1);  %Инициализация массива количеств отрезков в сложных обьектах
for o=1:size(ObjArr,1)    %Перебор сложных обьектов
	S(o)=size(ObjArr{o,type},1);   %Получение Размеров вложенных массивов
end
LstDfcltObj=find(S==0,1)-1; %Вычисление количества сложных обьектов
QStr=max(S);                %Вычисление максимального количества прямых в самом сложном обьекте
Angles=zeros(QStr,1);   %Инициализация массива улов
GradW=zeros(QStr,2);    % и весов
for o=1:LstDfcltObj     %Перебор сложных обьектов
    Angles(1:S(o))=abs(MstrStrData(ObjArr{o,type}(:,1),7)-SlveStrData(ObjArr{o,type}(:,2),7));  %Углы при каждой паре прямых
    GradW(1:S(o),1)=MstrGrW(ObjArr{o,type}(:,1));
    GradW(1:S(o),2)=SlveGrW(ObjArr{o,type}(:,2));
    CrntObjRlablty=sum(unique(GradW(1:S(o),1)))+sum(unique(GradW(1:S(o),2))); %Сумма весов всех прямых (повторяющиеся исключены)
    GradW=GradW./CrntObjRlablty;    %Пронормированные веса прямых, взятые по отношению суммарного градиента пикселей по прямой к количеству этих пикселей
    NrmL=min([MstrStrData(ObjArr{o,type}(:,1),12) SlveStrData(ObjArr{o,type}(:,2),12)],[],2);   %Веса прямых по их длинам (при каждом угле выбирается короткий отрезок)
    NrmL=NrmL./sum(NrmL);	%Пронормированные веса (по длинам)
    ObjAngle=sum(abs(Angles(1:S(o))-pi/2).*NrmL);   %Средневзвешенный (по длинам коротких отрезков) условно прямой угол по обьекту
    ObjCoord(resObjcntr,1)=sum(unique(MstrStrData(ObjArr{o,type}(:,1),5).*GradW(1:S(o),1)))+sum(unique(SlveStrData(ObjArr{o,type}(:,2),5).*GradW(1:S(o),2)));   %Вычисление координат центра обьекта
    ObjCoord(resObjcntr,2)=sum(unique(MstrStrData(ObjArr{o,type}(:,1),6).*GradW(1:S(o),1)))+sum(unique(SlveStrData(ObjArr{o,type}(:,2),6).*GradW(1:S(o),2)));  
    %ObjColor=mean(mean(Raster(fix(ObjCoord(resObjcntr,1)):ceil(ObjCoord(resObjcntr,1)),fix(ObjCoord(resObjcntr,2)):ceil(ObjCoord(resObjcntr,2)))));   %Вычисление цвета обьекта (по четырём пикселям) в точке найденного условного геометрического центра
        NrmL=([MstrStrData(ObjArr{o,type}(:,1),12) SlveStrData(ObjArr{o,type}(:,2),12)])./(sum(unique(MstrStrData(ObjArr{o,type}(:,1),12)))+sum(unique(SlveStrData(ObjArr{o,type}(:,2),12))));  %Нормированные длины (без повторяющихся)
        NrmlzdAzCosDev=sum(unique(NrmL(:,1)./MstrStrData(ObjArr{o,type}(:,1),10)))+sum(unique(NrmL(:,2)./SlveStrData(ObjArr{o,type}(:,2),10)))./(AzCosDev_TH);
    %NrmlzdObjColorDiff=min([ObjColor/LowBrght_TH (65535-ObjColor)/HghBrght_TH]);	%Пронормированное знаение разницы цвета обьекта обьекта и ближайшего порогового значения (0 или 65535)
    %CrntObjRlablty=sqrt((ObjAngle./RAngleAccTH).^2+(1./CrntObjRlablty./(minW_TH*S(o))).^2+(1/S(o)/TypclQStr_TH).^2+NrmlzdObjColorDiff.^2);
    CrntObjRlablty=sqrt((ObjAngle./RAngleAccTH).^2+(1./CrntObjRlablty./(minW_TH*S(o))).^2+(1/S(o)/TypclQStr_TH).^2+NrmlzdAzCosDev.^2);
    if CrntObjRlablty<sqrt(4)  %ObjAngle(o,type)<RAngleAccTH
      ObjRlablty(resObjcntr)=CrntObjRlablty;    %Копирование параметра достоверности в выходной массив
            plot(ObjCoord(resObjcntr,2),ObjCoord(resObjcntr,1)+TbleSize,'gx')            
            Mm=unique(ObjArr{o,type}(:,1));
            Ms=unique(ObjArr{o,type}(:,2));
            for m=1:size(Mm,1)
             plot(MstrStrData(Mm(m),1:2),MstrStrData(Mm(m),3:4)+TbleSize,'Color','red');
            end
            for s=1:size(Ms,1)
             plot(SlveStrData(Ms(s),1:2),SlveStrData(Ms(s),3:4)+TbleSize,'Color','red');
            end
      resObjcntr=resObjcntr+1;  %Инкремент счётчика отобранных обьектов
    end
end
ObjCoord(resObjcntr:end,:)=[];    %Удаление пустой части массивов
ObjRlablty(resObjcntr:end)=[];
    resObjcntr
    axis equal;
end
function [Ccntr,Coords]=FndColByColLinks(SrsRow)%StrtCol
%SrsRow - Номер исходной строки в StrBitMsk. В ней ищутся точки (true) в пределах этой функции.
%StrtCol - Номер столбца в котором на этой строке (SrsRow) расположена точка(true-бит) которая и инициировала 
% поиск других точек в этой строке. Если начальной точки нет, то StrtCol=0
%Coords - [Row Col] каждой точки в БитМаске, которая принадлежит к текущему обьекту
% расположенному в исходной строке. Каждая новая точка в новой строке Coords
%Ccntr - количество найденных точек или Sумма связанных (найденых) пар прямых
global StrBitMsk;
global RequestColMsk;
global RequestRowMsk;
[~,c]=find(StrBitMsk(SrsRow,:));%Адреса всех не нулевых элементов в исходной строке битмаски
c=c.*RequestColMsk(c);  %Удаление из списка найденых точек (в битмаске) тех, что уже учтены (в частности той, из-за которой был инициирован поиск в этой строке)
c(c==0)=[]; %Удаление из списка найденых строк (в битмаске) тех по которым уже производился поиск
Csize=size(c,2);    %Количество найденых точек (true-битов) в текуще строке
RequestRowMsk(SrsRow)=false;%Пометка строки как обработанной
if Csize==0 %Досрочный возврат из функции, если ни одной точки (кроме стартовой) найдено не было
  Coords=[];
  Ccntr=0;    
else
  Coords=zeros(Csize*3,2);  %Инициализация выходного массива
  Coords(1:Csize,1)=SrsRow; %Запись Row-координаты(№строки) всех найденых точек в текущей строке
  Coords(1:Csize,2)=c;      %Запись Col-координат (номеров столбцов) этих точек
  Ccntr=Csize;  %Инициализация счёчика искомых точек (true-бит или пар прямых, которые окажутся связаны с текущей прямой рекурсивно)
  %Поиск связанных прямых(точек в StrBitMsk) в столбцах где были найдены элементы текущего обьекта
  for i=1:Csize     %перебор столбцов, по которым найдены точки в текущей строке
    [s,Coords(Ccntr+1:Ccntr+s,:)]=FndRowByRowLinks(c(i));	%,SrsRow Копирование координат точек находящихся в одном столбце с каждой из найденных (в текущей строке) точек (а также копирование свезанных с ними других точек, координаты которых были найдены рекуривно)
    Ccntr=Ccntr+s;      %Ccntr хранит фактическое количество записанных координат точек
  end
  Coords(Ccntr+1:end,:)=[]; %Удаление неинициализированных координат
         %plot(Coords(:,2),Coords(:,1),'rx');
end
end
function [Ccntr,Coords]=FndRowByRowLinks(SrsCol)%StrtRow
%SrsCol - Номер исходного (текущего) столбца в StrBitMsk. В нём ищутся точки (биты равные true) в пределах этой функции.
%StrtRow - Номер строки где в этом (текущем) стролбце (SrsCol) расположена точка(true-бит) которая и инициировала 
% поиск других точек в этом стролбце. Если начальной точки нет, то StrtRow=0
%Coords - [Row Col] каждой точки в БитМаске, которая принадлежит к текущему обьекту
% расположенному в исходном столбце ([Master Slave]адреса прямых). Новая точка в новой строке Coords
%Ccntr - количество найденных точек или Sумма связанных (найденых) пар прямых
global StrBitMsk;
global RequestRowMsk;
global RequestColMsk;
[r,~]=find(StrBitMsk(:,SrsCol));%Адреса всех не нулевых элементов в исходном столбце битмаски
r=r.*RequestRowMsk(r);  %Удаление из списка уже обработанных точек (в битмаске)(в частности той, из-за которой был инициирован поиск в этой строке)
r(r==0)=[]; %Удаление из списка найденых строк (в битмаске) тех по которым уже производился поиск
Csize=size(r,1);
RequestColMsk(SrsCol)=false;  %Пометка столбца как обработанного  
if isempty(r) %Досрочный возврат из функции, если ни одной точки (кроме стартовой) найдено не было
  Coords=[];
  Ccntr=0;
else
  Coords=zeros(Csize*3,2);  %Инициализация выходного массива
  Coords(1:Csize,1)=r;      %Запись Row-координат(номеров строк) всех найденых точек в текущем столбце
  Coords(1:Csize,2)=SrsCol;	%Запись Col-координаты (номера единственного столбца) этих точек
  Ccntr=Csize;  %Инициализация счёчика искомых точек (true-бит или пар прямых, которые окажутся связаны с текущей прямой рекурсивно)
  %Поиск связанных прямых(точек в StrBitMsk) в строках где были найдены элементы текущего обьекта
  for i=1:Csize	%перебор строк, в которых найдены точки в текущем столбце
    [s,Coords(Ccntr+1:Ccntr+s,:)]=FndColByColLinks(r(i)); %,SrsCol Копирование координат точек находящихся в одной строке с каждой из найденных (в текущем столбце) точек (а также копирование свезанных с ними других точек, координаты которых были найдены рекуривно)
    Ccntr=Ccntr+s;  %Обновление счётчика
  end
  Coords(Ccntr+1:end,:)=[]; %Удаление неинициализированных координат
         %plot(Coords(:,2),Coords(:,1),'rx');
end
end
function [ResAdr]=RecuFindM(N)
%Функция для рекурсивного поиска компонентов (отрезков) составных Master-отрезков (прямых)
%N - номер отрезка в MstrStrData
%ResAdr - адреса (в MstrStrLinkList) всех связанных с N-тым (отрезком) Master-отрезков.
global MstrStrLinkList; %Получение доступа к списку (составных)связей Master-отрезков
ResAdr=find(MstrStrLinkList(:,1)==N);    %Простой поиск (первое приближение)
for i=1:size(ResAdr,1)      %Поиск (рекурсивный) по каждому найденному адресу
    nAdr=RecuFindM(MstrStrLinkList(ResAdr(i),2));   %Второй столбец MstrStrLinkList содержит отрезки с которыми связаны прямые из первого столбца
    ResAdr=cat(1,ResAdr,nAdr);  %Обьединение массивов номеров отрезков найденных ранее (простым поиском) и найденных рекурсивно в каждой итерации
end
end
function [ResAdr]=RecuFindS(N)
%Функция для рекурсивного поиска компонентов (отрезков) составных Slave-отрезков (прямых)
%N - номер отрезка в SlveStrData
%ResAdr - адреса (в SlveStrLinkList) всех связанных с N-тым (отрезком) Slave-отрезков.
global SlveStrLinkList; %Получение доступа к списку (составных)связей Slave-отрезков
ResAdr=find(SlveStrLinkList(:,1)==N);    %Простой поиск (первое приближение)
for i=1:size(ResAdr,1)      %Поиск (рекурсивный) по каждому найденному адресу
    nAdr=RecuFindS(SlveStrLinkList(ResAdr(i),2));   %Второй столбец SlveStrLinkList содержит отрезки с которыми связаны прямые из первого столбца
    ResAdr=cat(1,ResAdr,nAdr);  %Обьединение массивов номеров отрезков найденных ранее (простым поиском) и найденных рекурсивно в каждой итерации
end
end