function [ ] = Save_main(  )
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

        BuffSize=2;
        StrEdgesLayerS=zeros(572,15996,BuffSize);
        i=2000;
while i<=2000
    parfor cntr=1:BuffSize
    fileName=['Test1/',num2str(i+cntr-1),'.tif'];
    Raster=imread(fileName);
    %imshow(Raster);
        %SrsHstgrm=BldHstgrm( Raster, maxBrghtns-minBrghtns,minBrghtns, maxBrghtns);
            %tic
    FixedRaster=HstgrmAutoRotate( Raster, minBrghtns, maxBrghtns );
            %toc; disp('����������� ������������')
        %CntrdHstgrm=BldHstgrm( FixedRaster, maxBrghtns-minBrghtns,minBrghtns, maxBrghtns);
            %tic;
    SuperFixedRaster=HstgrmAutoStretch(FixedRaster,HstgrmStretchDepth,minBrghtns, maxBrghtns);
            %toc; disp('����������� ���������')
        %SFHstgrm=BldHstgrm( SuperFixedRaster, maxBrghtns-minBrghtns,minBrghtns, maxBrghtns);
    %plot(Hstgrm);
    
    %Adr=find(ImgID==i);
%     hold on
%     plot(TrgtCoord(Adr,1),TrgtCoord(Adr,2),'ro');
%     hold off;

        StrEdgesLayerS(:,:,cntr)=FindEdgeLines( SuperFixedRaster , 1,2,3,4);
        
    end
            for a=1:BuffSize
            StrEdgesLayer=StrEdgesLayerS(:,:,a);
            fileName=['StrEdgesLayersTest/',num2str(i),'.mat'];
            save(fileName,'StrEdgesLayer');
            i=i+1;
            end

%                HlfScr=100;   %���������� ���� ������ ������������ ����
%         for t=1:size(Adr,1)     %�������� �����
%             
%             Upfr=max([TrgtCoord(Adr(t),2)-HlfScr,1]);
%             DwnFr=min([TrgtCoord(Adr(t),2)+HlfScr,size(Raster,1)]);
%             LftFr=max([TrgtCoord(Adr(t),1)-HlfScr,1]);
%             RghtFr=min([TrgtCoord(Adr(t),1)+HlfScr,size(Raster,2)]);
%                 
%             
%                 %Coord(cntr,:)=Slow_FindEdgeLines( SuperFixedRaster ,Upfr,DwnFr,LftFr,RghtFr);
%          %Coord(cntr,:)= FindEdgeLines( SuperFixedRaster , Upfr,DwnFr,LftFr,RghtFr);
%          StrEdgesLayer= FindEdgeLines( SuperFixedRaster , Upfr,DwnFr,LftFr,RghtFr,Sgma);
%             StrEdgesLayer=single(StrEdgesLayer);
%          cntr=cntr+1; 
%             fileName=['KrnlSgma',num2str(fix(Sgma*100)),'/Img',num2str(i),'Trgt',num2str(t),'.mat'];
% %              WriteFileID=fopen(fileName,'ab');
% %              fwrite(WriteFileID, StrEdgesLayer, 'single');
% %              fclose(WriteFileID);
%             save(fileName,'StrEdgesLayer');
% %         Im=Raster(Upfr:DwnFr,LftFr:RghtFr);    %��������������� ��������� �����������            
% %             %ImFixed=uint16(mod(double(Im) + 32768, 65536)); %������������� ����������� ������� (�� ������������� ������) 
% %         ImFixed=FixedRaster(Upfr:DwnFr,LftFr:RghtFr);
% %         ImSuperFixed=SuperFixedRaster(Upfr:DwnFr,LftFr:RghtFr);
% %          
% %         Im=imresize(Im, 4);    %��������������� ��������� �����������
% %             %imshow(Im);            
% %         ImFixed=imresize(ImFixed, 4);    %��������������� ����������� c �������������� ������������ �������
% %             %imshow(ImFixed);
% %         ImSuperFixed=imresize(ImSuperFixed, 4);%��������������� ����������� c ���������� ������������
% %             clf;
% %             subplot(2,1,1),hold on,plot(SrsHstgrm),plot(CntrdHstgrm),plot(SFHstgrm),hold off;
% %                 TrplPict=cat(4,Im,ImFixed,ImSuperFixed);            
% % %BW1 = edge(I,'sobel');
% % % CannyEdges = edge(SuperFixedRaster,'canny');
% % %             TrgtCannyEdges=CannyEdges(Lftfr:RghtFr,UpFr:DwnFr);
% % %             TrgtCannyEdges=imresize(TrgtCannyEdges, 4);
% % %             TrplPict=cat(4,Im,ImSuperFixed,TrgtCannyEdges);
% % DiffRaster=diff(double(SuperFixedRaster),1,2);
% % %DIm=abs(DIm);
% %             DIm=DiffRaster(Upfr:DwnFr,LftFr:RghtFr);
% %             DIm=imresize(DIm, 4);
% %                 %TrplPict=cat(4,Im,ImSuperFixed,DIm);
% %             subplot(2,1,2),montage(TrplPict,'Size',[1 3]);
%         %imagesc(ImFixed);
%         end
    
end

end

