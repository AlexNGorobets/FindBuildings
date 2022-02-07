function [ StrArray ] = FindEdgeLines( Data )
%FINDEDGELINES �������� ������ ������� ��������
%Data - ����������� �����
%Coord - ���������� ������ �������� ([x1 x2 y1 y2] - � ������ ������ ��� ������� �������) 
% (��� �� �1 �2 r1 r2), �.�. ��� � ���������� ���� � ���� ��������� - ����� ������� ����
%��������� ��������/��������������
CoreSize=17;    %������ ���� ����������� �������
AdptvFltrngAzW=25;  %25������� ��������� ������������ �������� - ������������� ������� ������� ��������� ��� ���������� (1 - ������� �������� ��������� � ��������, )
Data=double(Data);  %���������� ���� ������ � double
i=2:2:size(Data,1)-2;   %� ���������� ������� ������ ������ ������ ����� ������
Data(i,:)=(Data(i-1,:)+Data(i+1,:))/2;	% �������� ��������� �����������
[SrsRN, SrsCN]=size(Data);   %������ �������� ������
%���������� ������ ������ �������� � �������� ���� ������������� ������� (����� ��������� ������������ ���������� �� ������� �������� �� �������� �������):
CatArr=zeros((CoreSize-1)/2,size(Data,2));
CatArr(1,:)=Data(1,:);
CatArr=cumsum(CatArr,1)+rand(size(CatArr));
Data=cat(1,CatArr,Data); %���������� �� ������� ��������� (����� ������)
CatArr(:,:)=0;
CatArr(1,:)=Data(end,:);
CatArr=cumsum(CatArr,1)+rand(size(CatArr));
Data=cat(1,Data,CatArr); %���������� �� ������� ��������� (����� �����)
CatArr=zeros(size(Data,1),(CoreSize-1)/2);
CatArr(:,1)=Data(:,1);
CatArr=cumsum(CatArr,2)+rand(size(CatArr));
Data=cat(2,CatArr,Data); %���������� �� ������� ��������� (������� �����)
CatArr(:,:)=0;
CatArr=zeros(size(Data,1),(CoreSize-1)/2);
CatArr(:,1)=Data(:,end);
CatArr=cumsum(CatArr,2)+rand(size(CatArr));
Data=cat(2,Data,CatArr); %���������� �� ������� ��������� (������� ������)

%DiffKernel=GetRoundDiffKernel(5, [1 3]);
DiffKernel=GetRoundDiffKernel(5, [3 1]);
%DiffKernel=[0 1 0;0 0 0;0 -1 0];
  DY=conv2(Data, DiffKernel,'same');    %������������� ����������� �� ������-����
  DX=conv2(Data, DiffKernel','same');   %�������������� ����������� �� � ���� �� �����
        %quiver(flip(DX(Upfr:DwnFr,LftFr:RghtFr),1),-flip(DY(Upfr:DwnFr,LftFr:RghtFr),1));
        %quiver(flip(DX,1),flip(-DY,1));
%DiffVal=sqrt(DX.^2+DY.^2);  %�������� ��������� ��� ������� �������
[DiffAngle,DiffVal]=cart2pol(DX,DY);
%DiffAngle=atan2(DY,DX);     %������ ��������� ��� ������� �������
%DiffAngle=atan2(-DY,DX);     %������ ��������� ��� ������� ������� (�� -pi �� pi)
    %DiffAngle=DiffAngle./pi.*180;%������� ������� � ������� (�� -180 �� 180 0 �������� - ��� � ������)
                tic;        %����� � ����������� �������������
%**********���� ���������� ����������
%� ������ ���������� ������������ ��������� ����� ����������� ���������(�������). 
% ��� ����� �������� ����� ����������� ��������� - � ������ �������� ���� �����
% ��������� ��� ��������� ������ Ncol �� CoreSize �� CoreSize. ������ ������� �������� ����� ����� ����� �� ������
        HCS=fix(CoreSize/2);	%Half-Core Size
%[Nr, Nc]=size(DiffAngle);   %����������� �������� ���������� ��� ������
AngDynDiff=zeros(SrsRN,SrsCN);    %������������� ����-���������� ������������� �������
        AngDynDiff2=AngDynDiff;
CoreOnString=zeros(SrsCN,CoreSize,CoreSize);   %������������� ������� ������� ��������
%c=1+HCS:size(DiffAngle,2)-HCS;  %������ - �������� �������
c=1:SrsCN;  %������ - �������� �������
Kernl=CoreOnString;     %������������� ���� (��� ���� �������� �����)
DiffValLayer=CoreOnString;%������������� ����-�������� ���������
%for r=1+HCS:size(DiffAngle,1)-HCS
for r=1:SrsRN
    for KrnlDr=1:CoreSize
        for KrnlDc=1:CoreSize
        CoreOnString(c,KrnlDr,KrnlDc)=abs(DiffAngle(r+KrnlDr-1,c+KrnlDc-1)-DiffAngle(r+HCS,c+HCS)); %���������� ������� �������� (��� ������� �������)
        DiffValLayer(c,KrnlDr,KrnlDc)=DiffVal(r+KrnlDr-1,c+KrnlDc-1);   %����������� ������� ���������� � ��������������� ���� (��� ������� ���������)
        end
    end
    CoreOnString(CoreOnString>pi)=abs(CoreOnString(CoreOnString>pi)-2*pi); %���������� ������ ����� � ������� �� -pi �� +pi � �����������
    K=Bld3DAzKernel( CoreSize, DiffAngle(r+HCS,c+HCS)); %���� ��� ������� ������ 
        %K=Bld3DAzKernel_Base( CoreSize, DiffAngle(r,c) ); %���� ��� ������� ������ 
    Kernl(c,:,:)=permute(K,[3 1 2]);    %������������ ������������ (������� ������ - ������ ������ ���������)
%               if r==32
%                   tmpC=25;
%                   imagesc(Data(r-HCS:r+HCS,tmpC-HCS:tmpC+HCS));
%                   hold on
%                   %quiver(flip(DX(r-HCS:r+HCS,tmpC-HCS:tmpC+HCS),1),flip(-DY(r-HCS:r+HCS,tmpC-HCS:tmpC+HCS),1));
%                   %quiver(DX(r-HCS:r+HCS,tmpC-HCS:tmpC+HCS),DY(r-HCS:r+HCS,tmpC-HCS:tmpC+HCS));
%                   meshc(0.6+K(:,:,tmpC));
%               end
    %CoreOnString=CoreOnString.*Kernl./DiffValLayer; %������ ������, � ����������� �� �������� �������� ���������
        %CoreOnString2=Kernl.*(CoreOnString./DiffValLayer);
    CoreOnString=DiffValLayer.*Kernl./(1+AdptvFltrngAzW.*CoreOnString);
    AngDynDiff(r,c)=sum(sum(CoreOnString(c,:,:),2),3);  %������������ �� ������ � �������� ����, ��������� ������ �����������
        %AngDynDiff2(r,c)=1./sum(sum(CoreOnString2(c,:,:),2),3);
end
%AngDynDiff(1+HCS:Nr-HCS,c)=1./AngDynDiff(1+HCS:Nr-HCS,c);   %������ �������� �������� (1/0 - ������������ ������� �����-��� ����������� ������� � �����)
clearvars CoreOnString DiffValLayer Kernl K;  %������������ ������
Data([1:HCS,end-HCS+1:end],:)=[];Data(:,[1:HCS,end-HCS+1:end])=[];  %�������� ����� �������� � �������� ���� ������������� �������
DiffAngle([1:HCS,end-HCS+1:end],:)=[];DiffAngle(:,[1:HCS,end-HCS+1:end])=[];  % (������� � �������� �������� ������)
DiffVal([1:HCS,end-HCS+1:end],:)=[];DiffVal(:,[1:HCS,end-HCS+1:end])=[];
DX([1:HCS,end-HCS+1:end],:)=[];DX(:,[1:HCS,end-HCS+1:end])=[];
DY([1:HCS,end-HCS+1:end],:)=[];DY(:,[1:HCS,end-HCS+1:end])=[];
    disp('������������ ����������:')
        toc;
         %tImg=AngDynDiff(Upfr:DwnFr,LftFr:RghtFr);
     %imagesc(AngDynDiff./mean(mean(AngDynDiff)));axis equal;
     %imagesc(AngDynDiff2./mean(mean(AngDynDiff2)));        
     tmp=1;
%*****���� ��������� ��������
%AddData - �������������� ������ �� ������ ������ ([CPR CPC GrAz])([���������������� �����([Row Col])�� ��������� ���������� ��������, ���� ��������� � ��������]).
% ������ ��������� ������ �� 90 �������� ���������� �� ������� ������.
%Rlablty - ��������� ������������� ��� ������ ������ ([TtlGrVal Qpix AzCosDev]) - ��������� 
% �������� ��������� �� ���� �������� ������, ���������� �������� � ������, ��������� �������������� �������� ��������� ������ � ����� ���������� �� ������� ������� (�� ���� ������� ����������������� ������������ ���������� ����� �������� �� ������� ������� �� ������).
StrMinLngth=10; %10 ����������� ����� ������� ������ (� ��������)
AzOwrleapFactor=2.5;  %2.5 1.5; 4.3;4 ����������� ���������� �� ������� (��� ��������� ������ �� ���� ���������� ����������� AzStepN ������� ��������. ��� ������� �������� ������� ����������� ����-�������_��������
% � ����������� ������������. 360/AzStepN �������� - ������ �������� �������. ����������� ���������� - ��� �������� ����� ��������� � AzOwrleapFactor �������� ������� ��������. ��� ���������� ���� �������_�������� ������� ������� � �������� ��������� �������� � ������� ��������� ���������)
MinPixQnttyOnStr=round(StrMinLngth)./sqrt(2); %����������� ���-�� �������� � ������� ������

Coord=zeros(100000,4);  %������������� �������� ��������
AddData=zeros(100000,3);
Rlablty=zeros(100000,3);
%SafeCoreFrame=2;    %SafeCoreFrame - ���������� ��� ����������� ���� ������ (�\���) ���� ��������
%BnrztnLvl=mean(mean(AngDynDiff))*1.5;   %1.2 0.8 1.0;1.5 1.8<----����� �����������
    %BnrztnLvl=2.326762770949994e+05
    BnrztnLvl=2.15e+05;
%BnrztnLvl=graythresh(AngDynDiff);
      MainBitMsk=AngDynDiff>BnrztnLvl;        %������� ����� ������� ����� ��������
      %imshow(MainBitMsk);
FrdgePoints=zeros(4,2); %���������� ����������� ������ ������� ������ � ������� �����, ������� ���������� ��������� ��������� ������� ������ ��������
cntr=1;     %��������� ������������� �������� ��������� ������
AzStepN=90; %72���������� ����� �� �������
    tic;
for nAz=1:AzStepN
  Az=-pi+nAz/AzStepN*2*pi;    %���������� �������� ������� � ��������
  AzDiffLayer=abs(DiffAngle-Az);  %���� ������ ���� ��������� � �������� �������
  AzDiffLayer(AzDiffLayer>pi)=abs(AzDiffLayer(AzDiffLayer>pi)-2*pi); %������� �� 0 �� 180 ��������
  AzDiffLayer=AngDynDiff./(AzDiffLayer./(AzOwrleapFactor*2*pi/AzStepN));	%������ �������� �������� � ���������� �� ���� ������������� �������. 2*pi/AzStepN ���������� ������ � ����� ����. ���� ������� ������� � ������� � �������� ������� ������ ��� �����������_����������*������ ������ ����, �� AngDynDiff ���������� ����, ���� ������, �� �������
  %BitMsk=AzDiffLayer>BnrztnLvl;   %�����������            %(2*360/AzStepN) - ������� � 2*� (�=�������� � ����� ����) ������������� ���������� � 1 �� �������� ���� �������
  BitMsk=(AzDiffLayer>BnrztnLvl)&MainBitMsk;   %����������� � ����������� ������� �������� (���������) �� �������� � ��������� ���-����� (MainBitMsk) ���������� ����� ����� ������������ ����������
  %BitMsk=(AzDiffLayer<(AzOwrleapFactor*2*pi/AzStepN))&MainBitMsk; %����������� ������ �� ������� ��������� � ������� �������� �������� ���������
%           imshow(flip(BitMsk,1)); axis auto xy
%           title(['Angle: ',num2str(fix(Az/pi*1800)/10),char(176)])
%            fileName=['Res/Bmsk2/',num2str(nAz)];
%            saveas(gcf,fileName,'png');
%         mov(nAz) = getframe(gcf);
  [ObjLayer, Num]=bwlabel(BitMsk,8);   %��������� �������� �� ������
  ObjLayer=ObjLayer(:);   %������������� ���� �������� (����� ������ ����������)
  ObjLayer=cat(2,ObjLayer,(1:size(ObjLayer,1))'); %��������� ������ ������� - ������ �������
  ObjLayer=sortrows(ObjLayer);    %���������� �� ����������� ������ ������� (������� �������)
  %ObjLayer(1:find(ObjLayer(:,1),1)-1,:)=[];   %�������� �������� �� �������� �� ���������� ������ (������� �����)
    i=find(ObjLayer(:,1),1);   %�������� �������� (����� ������� �� �������� ��������)
    ObjLayer=cat(1,ObjLayer,[Num+1 0]); %���������� ��������� (����������) ������ ��� ��������� ���� �������� (������� ���������) ������ �������������� �����.
    for nObj=1:Num %�������� ��������
        AdrA=i;     %��������� ����� ���������
        while(ObjLayer(i,1)<nObj+1) %������� ���������� �������� � ������� ������ (nObj)
            i=i+1;
        end
        AdrB=i-1;   %�������� ����� ���������
         if (AdrB-AdrA)>=MinPixQnttyOnStr-1;    %������� ������������� ���������� �������� ��� ���������� ���������
            CartesVectSumX=sum(DX(ObjLayer(AdrA:AdrB,2)));  %����� �������� ��������� �� �������� ������� ������
            CartesVectSumY=sum(DY(ObjLayer(AdrA:AdrB,2)));
            %GrpAzmth=atan2(CartesVectSumY,CartesVectSumX); %������ ���������� ������� (���������������� ������ ������)
            GrpAzmth=atan2(CartesVectSumY,CartesVectSumX);
            TtlPixGrVal=sum(DiffVal(ObjLayer(AdrA:AdrB,2)));   %��������� �������� ������ (������������ �����)
            PixW=DiffVal(ObjLayer(AdrA:AdrB,2))./TtlPixGrVal;  %���� �������� (����������������)
            [r, c]=ind2sub(size(DiffVal),ObjLayer(AdrA:AdrB,2)); %������ (������ � �������) ������ ����� � ������� ������
            LftFr=min(c);   %����� ���� ������ (������������ ������)
            RghtFr=max(c);  %������ ���� ������
            UpFr=min(r);    %������� ���� ������ (�������������� ������)
            DwnFr=max(r);   %������ ����
            CPR=sum(r.*PixW);   %���������������� ��������� ������ ������� ������ �������� (������)
            CPC=sum(c.*PixW);                                                             %(�������)
            k=tan(GrpAzmth+pi/2);  %�-�� ��� ��������� y=kx
            FrdgePoints(1,:)=[LftFr k*(LftFr-CPC)+CPR];   %X � Y ����� ����������� ������ ������ ���� ����� � ������ ���������������� ������� ������
            FrdgePoints(2,:)=[(UpFr-CPR)/k+CPC UpFr];   %����� ����������� ������ �������� ���� ����� � ������ ���������������� ������� ������
            FrdgePoints(3,:)=[RghtFr k*(RghtFr-CPC)+CPR];   %�� �� ��� ������� ����
            FrdgePoints(4,:)=[(DwnFr-CPR)/k+CPC DwnFr];   %�� �� ��� ������� ����
            [~,on]=inpolygon(FrdgePoints(:,1),FrdgePoints(:,2),([LftFr RghtFr RghtFr LftFr]),([UpFr UpFr DwnFr DwnFr]));    %����� ����� ����������� � ������ � �������� � ��������
              Coord(cntr,1:2)=FrdgePoints(on,1);    %X1X2 ������� ������
              Coord(cntr,3:4)=FrdgePoints(on,2);    %Y1Y2 ������� ������
            if (sqrt((Coord(cntr,1)-Coord(cntr,2))^2+(Coord(cntr,3)-Coord(cntr,4))^2))>StrMinLngth    %�������� ������ ������� (������� �������� �������������)
              AddData(cntr,:)=[CPR,CPC,GrpAzmth];    %������ �������������� ������ (������+������� ����������� ����� ������� (�� �����), ���������������� ������ ���������)
              TtlGrVal=sqrt(CartesVectSumX^2+CartesVectSumY^2); %�������� ��������� ������ (��������� ����� ���������� �� ������� �������)
              Rlablty(cntr,:)=[TtlGrVal,AdrB-AdrA+1,TtlGrVal/TtlPixGrVal];  %�������� ������������� �� ������ ������ (�������) - �������� �������� ������, ��������� ��������, ������� �������� �����
              cntr=cntr+1;  %��������� �������� ��� ��������� ������
            end                          
         end   
%                 if (AdrB-AdrA>100); %(nAz==5)&& %���� ��� ��������� ������� ���������� � ��������
%                     %imshow(flip(BitMsk,1)); axis auto xy
%                     imshow(BitMsk); axis auto xy
%                     hold on
%                     plot([LftFr RghtFr RghtFr LftFr LftFr], [UpFr UpFr DwnFr DwnFr UpFr],'green');  %�����
%                     plot(Coord(cntr-1,1:2),Coord(cntr-1,3:4),'magenta');    %�������
%                     plot ([FrdgePoints(:,1);FrdgePoints(1,1)],[FrdgePoints(:,2);FrdgePoints(1,2)],'rx');    %����� ����������� ������ � �������
%                     plot(CPC,CPR,'ro'); %���������������� ����� �������
%                     hold off
%                 end
    end
end
%     v = VideoWriter('Res/BitMask','Archival');
%     open(v);
%     writeVideo(v,mov);
%     close(v);
cntr=cntr-1;    %������ ���������� ���������� ��������
  Coord(cntr+1:end,:)=[];   %�������� ������ ����������������� ������� �������� ��������.
  AddData(cntr+1:end,:)=[]; % � ������������� �������� ��������� ����� (����� ������ ����� ��������� ����� ����� ���������� �����)
  Rlablty(cntr+1:end,:)=[];
    disp('������������ ��������:')
    toc;
       imagesc(Data);axis equal
        hold on
        Nlevels=1000;
        Hcolor=hot(Nlevels);
        for i=1:size(Coord,1)   %����������� ���� ��������� ������
          color=Rlablty(i,3);
          %color=Rlablty(i,3)*Nlevels;
          %plot(Coord(i,1:2),Coord(i,3:4),'Color',Hcolor(round(color),:));
          plot(Coord(i,1:2),Coord(i,3:4),'Color',[color color color]);          
        end
        hold off
%******���� ������������ � ���������� ������ ������---
StrLngth=sqrt((Coord(:,2)-Coord(:,1)).^2+(Coord(:,4)-Coord(:,3)).^2);   %����� ������
%NrmlzdGrad=Rlablty(:,1)./Rlablty(:,2);  %�������� ���������� ������������� ��������� ��������� � ���������� �������� ������ ������ (������ ������)
    %TrustInd=sqrt((StrLngth./mean(StrLngth,1)).^2+(NrmlzdGrad./mean(NrmlzdGrad,1)).^2+(Rlablty(:,3)./mean(Rlablty(:,3),1)).^2); %������ ���������� ������ - ��������� ������� � ������������ �����-��������_���������-�������_�������� ��� ������ �������� �� ������ ������
TrustInd=Rlablty(:,3);  %������ ���������� (��������������������� �������� �������� � ������ ��������)
    %TrustInd=sqrt((Rlablty(:,3)./mean(Rlablty(:,3))).^2+(NrmlzdGrad./mean(NrmlzdGrad)).^2);
    %TrustInd=sqrt(6.*(Rlablty(:,3)./mean(Rlablty(:,3),1)).^2+(Rlablty(:,2)./mean(Rlablty(:,2),1)).^2);
StrArray=cat(2,Coord,AddData,Rlablty,TrustInd,StrLngth); %����������� ���� ���������� ��������� ������ � ���� ������ (�� ������ ��� ������ ������)
%StrArray(:,1:6)=StrArray(:,1:6)+SafeCoreFrame;  %���� (��������) ����� ������ �� ������ 
StrArray=sortrows(StrArray,-11);    %���������� ������ �� ����������
StrArray(:,11)=1:cntr;      %������ ������� ���������� �� �������������� ������ (��� ������ ID, ��� ���� ����������)
    clearvars Coord AddData Rlablty;
%���������� ������ � ������
ToStrDltMsk=false(cntr,1);
AngTreshold=1.8*(2*pi/AzStepN); %1.8 %����� ������� �� ������� ��� ������ ������ (+-AngTreshold ������ �� �������� ������� ������� ������ (������� ID))
%����� ������� ������� ������ ��������� ��� ������������ ��� ����� ������������ ������ (������ � ������� ID)
PrllelKTrshold=1/6;	%1\8 1/4����� ������� ��������� ������ ������� �� ��� ������������ ������ � ������� ID (� ����� ����� ������� ������)
AcrossKTrshold=0.5*StrMinLngth;	%������ ������� ��������� �� ��� ���������������� ������� ������ (� ��������)
    LngthTrshold=0.4;   %0.5 ������ ������� ����� ������� � ����������� ������ (� ����� ����� ������� ������) - �����-������� � ������������ - �� ������ ��������
    tic;
GK=tan(StrArray(:,7));    %����������� ������ ��������� (��������������) ��� ��������� ���� y=kx
CrntIDStrK=tan(StrArray(:,7)-pi/2);%����������� ������� (�������������) ������ ��� ��������� ���� y=kx
SameSourseStrMterick=zeros(cntr,1);                     %������������� ������� ������ ��������� ������
for ID=1:cntr-1
  if ToStrDltMsk(ID)==false
      SameSourseStrMterick(ID+1:end)=abs(StrArray(ID+1:end,7)-StrArray(ID,7)); %����������� ������ �������� ������������ ������ � ���� ���������
      SameSourseStrMterick(SameSourseStrMterick>=pi)=abs(SameSourseStrMterick(SameSourseStrMterick>=pi)-2*pi);%���������� ������������ ������� � �������� �� 0 �� 180 ��������
      %���������� ����� ����������� ����� ��������� � ������� ������
      x=(CrntIDStrK(ID)*StrArray(ID,6)-GK(ID+1:end).*StrArray(ID+1:end,6)+StrArray(ID+1:end,5)-StrArray(ID,5))./(CrntIDStrK(ID)-GK(ID+1:end));
      y=CrntIDStrK(ID).*(x-StrArray(ID,6))+StrArray(ID,5); %���������� ����� ����������� ����� ��������� �� ������ ������ � ������� ������
      PrllelDist=sqrt((x-StrArray(ID,6)).^2+(y-StrArray(ID,5)).^2)./StrArray(ID,12);	%���������� (� ����� ����� ������� ������) ����� ������������ ������� ����� ��� ������������ ������������ ������ (���������� �� ����� ����������� ����� ��������� � ������� ������ �� ����������� ����� �������(������������) ������)
      AcrossDist=sqrt((x-StrArray(ID+1:end,6)).^2+(y-StrArray(ID+1:end,5)).^2);    %���������� ���������� (� ��������) ����� ������� (������ � �������) ����� ����� ��������� (�������������� ������ �����).
      LngthDiff=abs(StrArray(ID,12)-StrArray(ID+1:end,12))/StrArray(ID,12);      %������� ����� ������� � ����������� ������ � ����� ������� ������
      %SameSourseStrMterick(ID+1:end)=sqrt((SameSourseStrMterick(ID+1:end)./AngTreshold).^2+(PrllelDist./PrllelKTrshold).^2+(AcrossDist./AcrossKTrshold).^2);  %���� ������������ ������ (�� ������� ��� ������ ���������� ��� ���������� ���������� � ������������ ������ �������) - ��������� ������� � (����������������) ������������ �������_�������� - ������������_���������� - �����������_����������. ��� ��������� ��������������
      SameSourseStrMterick(ID+1:end)=sqrt((SameSourseStrMterick(ID+1:end)./AngTreshold).^2+(PrllelDist./PrllelKTrshold).^2+(AcrossDist./AcrossKTrshold).^2+(LngthDiff/LngthTrshold).^2);  %���� ������������ ������ (�� ������� ��� ������ ���������� ��� ���������� ���������� � ������������ ������ �������) - ��������� ������� � (����������������) ������������ �������_�������� - ������������_���������� - �����������_����������. ��� ��������� ��������������
      %SameSourseStrMterick(ID+1:end)=sqrt((SameSourseStrMterick(ID+1:end)./AngTreshold).^2+(PrllelDist./PrllelKTrshold).^2+(AcrossDist./AcrossKTrshold(ID+1:end)).^2);  %���� ������������ ������ (�� ������� ��� ��� ������ ���������� ��� ���������� ���������� � ������������ ������ �������) - ��������� ������� � (����������������) ������������ �������_�������� - ������������_���������� - �����������_����������. ��� ��������� ��������������
      ToStrDltMsk(ID+1:end)=ToStrDltMsk(ID+1:end)|SameSourseStrMterick(ID+1:end)<2; %2=sqrt(4) - ��� ������������� ������������
  end
end
StrArray(ToStrDltMsk,:)=[];    %�������� ��������� ��������� ������
  % Coord=StrArray(:,1:4);  %�������� ������� �������� ���������� ������ � ������� ��������������� �� ���������� ������
  % AddData=StrArray(:,5:7);
  % Rlablty=StrArray(:,8:10);
  disp('������������ � ����������:')
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

