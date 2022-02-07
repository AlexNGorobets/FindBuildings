function [ ] = Load_main(  )
%MAIN Онсовной сценарий распознаваня зданий
%В папке /Train1 есть набор .tif файлов пронумерованных от 1 до 1000
% в ней так же хранится objects.csv с координатами обьектов на снимках
HstgrmStretchDepth=0.95;
fileID=fopen('Train1/objects.csv','r'); %Открытие файла с координатами в /Train1
Coordinatess=textscan(fileID,'%d%d%d',-1,'headerLines',1,'delimiter',{';',','}); %Считывание координат
fclose(fileID); %Закрытие .csv-файла
ImgID=Coordinatess{1}; %Столбец идентефикаторов изображение (соответствует по номеру строки координатам для каждой цели)
TrgtCoord(:,1)=Coordinatess{2}; %X-координата на изображении (Номер столбца)
TrgtCoord(:,2)=Coordinatess{3}; %Y-координата на изображении (Номер строки)
minBrghtns=0;	%Предустановленное минимально возможное значение яркости
maxBrghtns=2^16-1;  %и максимальное, сооветственно

% for i=1:1000
%i=5
        cntr=1;        
for i=5:1000  %Перегрузка яркости 413 425
     fileName=['Train1/',num2str(i),'.tif'];
     Raster=imread(fileName);
     FixedRaster=HstgrmAutoRotate( Raster, minBrghtns, maxBrghtns );
     SuperFixedRaster=HstgrmAutoStretch(FixedRaster,HstgrmStretchDepth,minBrghtns, maxBrghtns);
%     Adr=find(ImgID==i);
    fileName=['StrEdgesLayersTrain/',num2str(i),'.mat'];
    load(fileName,'StrEdgesLayer');
    StrArray=Load_FindEdgeLines( SuperFixedRaster, StrEdgesLayer);
    FindBldObj(StrArray,SuperFixedRaster);
    
        %plot(TrgtCoord(ImgID==i,1),TrgtCoord(ImgID==i,2),'go');
    
end

end

