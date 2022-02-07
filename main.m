function [ ] = main(  )
%MAIN �������� �������� ������������ ������
%� ����� /Train1 ���� ����� .tif ������ ��������������� �� 1 �� 1000
% � ��� ��� �� �������� objects.csv � ������������ �������� �� �������
HstgrmStretchDepth=0.95;
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
for i=5:1000   %12 - ������ ��������
    fileName=['Train1/',num2str(i),'.tif'];
    Raster=imread(fileName);

    FixedRaster=HstgrmAutoRotate( Raster, minBrghtns, maxBrghtns );
    SuperFixedRaster=HstgrmAutoStretch(FixedRaster,HstgrmStretchDepth,minBrghtns, maxBrghtns);
    Adr=find(ImgID==i);
%       StrArr=FindEdgeLines( SuperFixedRaster );
%       FindBldObj(StrArr,SuperFixedRaster,TrgtCoord(Adr(:),:))
               HlfScr=150;   %���������� ���� ������ ������������ ����
        for t=1:size(Adr,1)     %�������� �����
            
            Upfr=max([TrgtCoord(Adr(t),2)-HlfScr,1]);
            DwnFr=min([TrgtCoord(Adr(t),2)+HlfScr,size(Raster,1)]);
            LftFr=max([TrgtCoord(Adr(t),1)-HlfScr,1]);
            RghtFr=min([TrgtCoord(Adr(t),1)+HlfScr,size(Raster,2)]);

         ImSuperFixed=SuperFixedRaster(Upfr:DwnFr,LftFr:RghtFr);
         %ImSuperFixed=imresize(ImSuperFixed, 4);%��������������� ����������� c ���������� ������������
                %imshow(ImSuperFixed);
             StrArr=FindEdgeLines( ImSuperFixed );
%             FindBldObj(StrArr,ImSuperFixed,TrgtCoord(Adr(t),:));

         end
    
end

end

