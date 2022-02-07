function [ ] = main(  )
%MAIN Онсовной сценарий распознаваня зданий
%В папке /Train1 есть набор .tif файлов пронумерованных от 1 до 1000
% в ней так же хранится objects.csv с координатами обьектов на снимках
HstgrmStretchDepth=0.95;
    tic;
fileID=fopen('Train1/objects.csv','r'); %Открытие файла с координатами в /Train1
Coordinatess=textscan(fileID,'%d%d%d',-1,'headerLines',1,'delimiter',{';',','}); %Считывание координат
fclose(fileID); %Закрытие .csv-файла
ImgID=Coordinatess{1}; %Столбец идентефикаторов изображение (соответствует по номеру строки координатам для каждой цели)
TrgtCoord(:,1)=Coordinatess{2}; %X-координата на изображении (Номер столбца)
TrgtCoord(:,2)=Coordinatess{3}; %Y-координата на изображении (Номер строки)
    %toc;
    %disp('координаты загружены')
minBrghtns=0;	%Предустановленное минимально возможное значение яркости
maxBrghtns=2^16-1;  %и максимальное, сооветственно

% for i=1:1000
%i=5
        cntr=1;        
for i=5:1000   %12 - низкий контраст
    fileName=['Train1/',num2str(i),'.tif'];
    Raster=imread(fileName);

    FixedRaster=HstgrmAutoRotate( Raster, minBrghtns, maxBrghtns );
    SuperFixedRaster=HstgrmAutoStretch(FixedRaster,HstgrmStretchDepth,minBrghtns, maxBrghtns);
    Adr=find(ImgID==i);
%       StrArr=FindEdgeLines( SuperFixedRaster );
%       FindBldObj(StrArr,SuperFixedRaster,TrgtCoord(Adr(:),:))
               HlfScr=150;   %Полуразмер окна вывода окрестностей цели
        for t=1:size(Adr,1)     %Итератор целей
            
            Upfr=max([TrgtCoord(Adr(t),2)-HlfScr,1]);
            DwnFr=min([TrgtCoord(Adr(t),2)+HlfScr,size(Raster,1)]);
            LftFr=max([TrgtCoord(Adr(t),1)-HlfScr,1]);
            RghtFr=min([TrgtCoord(Adr(t),1)+HlfScr,size(Raster,2)]);

         ImSuperFixed=SuperFixedRaster(Upfr:DwnFr,LftFr:RghtFr);
         %ImSuperFixed=imresize(ImSuperFixed, 4);%Масштабирование изображения c растянутой гистограммой
                %imshow(ImSuperFixed);
             StrArr=FindEdgeLines( ImSuperFixed );
%             FindBldObj(StrArr,ImSuperFixed,TrgtCoord(Adr(t),:));

         end
    
end

end

