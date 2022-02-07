function [ ] = Load_main(  )
%MAIN �������� �������� ������������ ������
%� ����� /Train1 ���� ����� .tif ������ ��������������� �� 1 �� 1000
% � ��� ��� �� �������� objects.csv � ������������ �������� �� �������
HstgrmStretchDepth=0.95;
fileID=fopen('Train1/objects.csv','r'); %�������� ����� � ������������ � /Train1
Coordinatess=textscan(fileID,'%d%d%d',-1,'headerLines',1,'delimiter',{';',','}); %���������� ���������
fclose(fileID); %�������� .csv-�����
ImgID=Coordinatess{1}; %������� ��������������� ����������� (������������� �� ������ ������ ����������� ��� ������ ����)
TrgtCoord(:,1)=Coordinatess{2}; %X-���������� �� ����������� (����� �������)
TrgtCoord(:,2)=Coordinatess{3}; %Y-���������� �� ����������� (����� ������)
minBrghtns=0;	%����������������� ���������� ��������� �������� �������
maxBrghtns=2^16-1;  %� ������������, �������������

% for i=1:1000
%i=5
        cntr=1;        
for i=5:1000  %���������� ������� 413 425
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

