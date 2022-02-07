function [ ] = Paper_main(  )
%MAIN �������� �������� ������������ ������
%� ����� /Train1 ���� ����� .tif ������ ��������������� �� 1 �� 1000
% � ��� ��� �� �������� objects.csv � ������������ �������� �� �������
HstgrmStretchDepth=0.65;
    tic;
fileID=fopen('Train1/objects.csv','r'); %�������� ����� � ������������ � /Train1
Coordinatess=textscan(fileID,'%d%d%d',-1,'headerLines',1,'delimiter',{';',','}); %���������� ���������
fclose(fileID); %�������� .csv-�����
ImgID=Coordinatess{1}; %������� ��������������� ����������� (������������� �� ������ ������ ����������� ��� ������ ����)
TrgtCoord(:,1)=Coordinatess{2}; %X-���������� �� ����������� (����� �������)
TrgtCoord(:,2)=Coordinatess{3}; %Y-���������� �� ����������� (����� ������)
    %toc;
    %disp('���������� ���������')
minBrghtns=0;	%����������������� ���������� ��������� �������� �������
maxBrghtns=2^16-1;  %� ������������, �������������

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
            %toc; disp('����������� ������������')
        CntrdHstgrm=BldHstgrm( FixedRaster, maxBrghtns-minBrghtns,minBrghtns, maxBrghtns);
            %tic;
    SuperFixedRaster=HstgrmAutoStretch(FixedRaster,HstgrmStretchDepth,minBrghtns, maxBrghtns);
    %imshow(SuperFixedRaster(352-90:352+90,11388-90:11388+90));
            %toc; disp('����������� ���������')
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
    %����������� ������� ������������ ����������
    %ClassicCase(SuperFixedRaster(352-50:352+50,11388-70:11388+80));
    tic;
    ClassicCase(SuperFixedRaster(:,end-4700:end));
    toc;
    %ClassicCase(SuperFixedRaster(:,end-1600:end));
    
end

end

