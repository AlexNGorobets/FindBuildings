function [ ] = Paper_main(  )
%MAIN Онсовной сценарий распознаваня зданий
%В папке /Train1 есть набор .tif файлов пронумерованных от 1 до 1000
% в ней так же хранится objects.csv с координатами обьектов на снимках
HstgrmStretchDepth=0.65;
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
for i=5:1000    %9
    fileName=['Train1/',num2str(i),'.tif'];
    Raster=imread(fileName);
    %imshow(Raster(352-90:352+90,11388-90:11388+90));
        SrsHstgrm=BldHstgrm( Raster, maxBrghtns-minBrghtns,minBrghtns, maxBrghtns);
            %tic
    FixedRaster=HstgrmAutoRotate( Raster, minBrghtns, maxBrghtns );
    %imshow(FixedRaster(352-90:352+90,11388-90:11388+90));
            %toc; disp('гистограмма центрирована')
        CntrdHstgrm=BldHstgrm( FixedRaster, maxBrghtns-minBrghtns,minBrghtns, maxBrghtns);
            %tic;
    SuperFixedRaster=HstgrmAutoStretch(FixedRaster,HstgrmStretchDepth,minBrghtns, maxBrghtns);
    %imshow(SuperFixedRaster(352-90:352+90,11388-90:11388+90));
            %toc; disp('гистограмма растянута')
%-        SFHstgrm=BldHstgrm( SuperFixedRaster, maxBrghtns-minBrghtns,minBrghtns, maxBrghtns);
    %StrArray=FindEdgeLines( SuperFixedRaster(352-50:352+50,11388-70:11388+80));
    %StrArray=FindEdgeLines( SuperFixedRaster(370-80:370+80,2420-70:2420+80));
     StrArray=FindEdgeLines( SuperFixedRaster(:,end-4700:end));
         FindBldObj( StrArray ,SuperFixedRaster(:,end-4700:end), [[],[]]);
      %FindBldObj( StrArray ,SuperFixedRaster(352-50:352+50,11388-70:11388+80), [[],[]]);
        %FindBldObj( StrArray ,SuperFixedRaster(370-80:370+80,2420-70:2420+80), [[],[]]);
%     StrArray=FindEdgeLines( SuperFixedRaster);
%     Adr=find(ImgID==i);
%     FindBldObj( StrArray ,SuperFixedRaster, TrgtCoord(Adr(:),:));
    %Аналогичное решение стандартными средствами
    %ClassicCase(SuperFixedRaster(352-50:352+50,11388-70:11388+80));
    tic;
    ClassicCase(SuperFixedRaster(:,end-4700:end));
    toc;
    %ClassicCase(SuperFixedRaster(:,end-1600:end));
    
end

end

