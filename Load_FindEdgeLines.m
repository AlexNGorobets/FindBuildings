function [ StrArray] = Load_FindEdgeLines( Data, AngDynDiff)
%FINDEDGELINES �������� ������ ������� ��������
%Data - ����������� �����
%Coord - ���������� ������ �������� ([x1 x2 y1 y2] - � ������ ������ ��� ������� �������) 
% (��� �� �1 �2 r1 r2), �.�. ��� � ���������� ���� � ���� ��������� - ����� ������� ����
%AddData - �������������� ������ �� ������ ������ ([CPR CPC GrAz])([���������������� �����([Row Col])�� ��������� ���������� ��������, ���� ��������� � ��������]).
% ������ ��������� ������ �� 90 �������� ���������� �� ������� ������.
%Rlablty - ��������� ������������� ��� ������ ������ ([TtlGrVal Qpix AzCosDev]) - ��������� 
% �������� ��������� �� ���� �������� ������, ���������� �������� � ������, ��������� �������������� �������� ��������� ������ � ����� ���������� �� ������� ������� (�� ���� ������� ����������������� ������������ ���������� ����� �������� �� ������� ������� �� ������).
%��������� ��������/��������������
StrMinLngth=10; %10 ����������� ����� ������� ������ (� ��������)
StrMaxLngth=150;%������������ ����� ������, ������� ����� ���� ������ ��������
AzOwrleapFactor=4.5;  %4.3;4 ����������� ���������� �� ������� (��� ��������� ������ �� ���� ���������� ����������� AzStepN ������� ��������. ��� ������� �������� ������� ����������� ����-�������_��������
% � ����������� ������������. 360/AzStepN �������� - ������ �������� �������. ����������� ���������� - ��� �������� ����� ��������� � AzOwrleapFactor �������� ������� ��������. ��� ���������� ���� �������_�������� ������� ������� � �������� ��������� �������� � ������� ��������� ���������)
MinPixQnttyOnStr=round(StrMinLngth*sqrt(2)); %����������� ���-�� �������� � ������� ������

Coord=zeros(100000,4);  %������������� �������� ��������
AddData=zeros(100000,3);
Rlablty=zeros(100000,3);
Data=double(Data);
SafeCoreFrame=2;    %SafeCoreFrame - ���������� ��� ����������� ���� ������ (�\���) ���� ��������

DiifKernel=GetRoundDiffKernel(5, [3 1]);    %����������� ���� (����� ��� � ���� ������)
  DY=conv2(Data, -DiifKernel,'valid');    %������������� ����������� �� ������-����. ��������� ���������� ������� � ������� (������� ����� ������ R � C) �������� ��������� ���������� ������ ���������� �� ���������, �.� ��� Y ���������� � �������� �������, � �������� ������-������� � �������� �� ���������� �� ������� ������� � ������ ����
  DX=conv2(Data, -DiifKernel','valid');   %�������������� ����������� �� � ���� �� �����
    %quiver(DX,DY);
DiffVal=sqrt(DX.^2+DY.^2);  %�������� ��������� ��� ������� �������
    %imagesc(DiffVal);
DiffAngle=atan2(DY,DX);     %������ ��������� ��� ������� ������� (�� -pi �� pi)
%DiffAngle=DiffAngle./pi.*180;%������� ������� � ������� (�� -180 �� 180 0 �������� - ��� � ������)

BnrztnLvl=mean(mean(AngDynDiff))*1.5;   %1.5 1.8<----����� �����������
% BitMsk=AngDynDiff>BnrztnLvl;        %������� ����� ������� ����� ��������
% [ObjLayer Num]=bwlabel(BitMsk,8);   %����������� �������� � ����� �� �������
%---������� � ������������� �� ������� �������
FrdgePoints=zeros(4,2); %���������� ����������� ������ ������� ������ � ������� �����, ������� ���������� ��������� ��������� ������� ������ ��������
cntr=1;     %��������� ������������� �������� ��������� ������
%AzStepN=90;
AzStepN=72; %72���������� ����� �� �������
    tic;
for nAz=1:AzStepN
  Az=-pi+nAz/AzStepN*2*pi;    %���������� �������� ������� � ��������
  AzDiffLayer=abs(DiffAngle-Az);  %���� ������ ���� ��������� � �������� �������
  AzDiffLayer(AzDiffLayer>pi)=abs(AzDiffLayer(AzDiffLayer>pi)-2*pi); %������� �� 0 �� 180 ��������
  AzDiffLayer=AngDynDiff./(AzDiffLayer./(AzOwrleapFactor*2*pi/AzStepN));	%������ �������� �������� � ���������� �� ���� ������������� �������
  BitMsk=AzDiffLayer>BnrztnLvl;   %�����������            %(2*360/AzStepN) - ������� � 2*� (�=�������� � ����� ����) ������������� ���������� � 1 �� �������� ���� �������
  [ObjLayer, Num]=bwlabel(BitMsk,8);   %��������� �������� �� ������
  ObjLayer=ObjLayer(:);   %������������� ���� �������� (����� ������ ����������)
  ObjLayer=cat(2,ObjLayer,(1:size(ObjLayer,1))'); %��������� ������ ������� - ������ �������
  ObjLayer=sortrows(ObjLayer);    %���������� �� ����������� ������ ������� (������� �������)
  ObjLayer(1:find(ObjLayer(:,1),1)-1,:)=[];   %�������� �������� �� �������� �� ���������� ������ (������� �����)
    i=1;   %�������� ��������
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
            GrpAzmth=atan2(CartesVectSumY,CartesVectSumX); %������ ���������� ������� (���������������� ������ ������)
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
            if (sqrt((Coord(cntr,1)-Coord(cntr,2))^2+sqrt((Coord(cntr,3)-Coord(cntr,4))^2)))>StrMinLngth    %�������� ������ ������� (������� �������� �������������)
              AddData(cntr,:)=[CPR,CPC,GrpAzmth];    %������ �������������� ������ (������+������� ����������� ����� ������� (�� �����), ���������������� ������ ���������)
              TtlGrVal=sqrt(CartesVectSumX^2+CartesVectSumY^2); %�������� ��������� ������ (��������� ����� ���������� �� ������� �������)
              Rlablty(cntr,:)=[TtlGrVal,AdrB-AdrA+1,TtlGrVal/TtlPixGrVal];  %�������� ������������� �� ������ ������ (�������) - �������� �������� ������, ��������� ��������, ������� �������� �����
              cntr=cntr+1;  %��������� �������� ��� ��������� ������
            end                          
         end        
    end
end
cntr=cntr-1;    %������ ���������� ���������� ��������
  Coord(cntr+1:end,:)=[];   %�������� ������ ����������������� ������� �������� ��������.
  AddData(cntr+1:end,:)=[]; % � ������������� �������� ��������� ����� (����� ������ ����� ��������� ����� ����� ���������� �����)
  Rlablty(cntr+1:end,:)=[];
    toc;
%-----������������ � ���������� ������ ������---
StrLngth=sqrt((Coord(:,2)-Coord(:,1)).^2+(Coord(:,4)-Coord(:,3)).^2);   %����� ������
%NrmlzdGrad=Rlablty(:,1)./Rlablty(:,2);  %�������� ���������� ������������� ��������� ��������� � ���������� �������� ������ ������ (������ ������)
    %NrmlzdGrad=Rlablty(:,1)./StrLngth(:);
%TrustInd=sqrt((StrLngth./mean(StrLngth,1)).^2+(NrmlzdGrad./mean(NrmlzdGrad,1)).^2+(Rlablty(:,3)./mean(Rlablty(:,3),1)).^2); %������ ���������� ������ - ��������� ������� � ������������ �����-��������_���������-�������_�������� ��� ������ �������� �� ������ ������
TrustInd=Rlablty(:,3);
    %TrustInd=sqrt((Rlablty(:,3)./mean(Rlablty(:,3))).^2+(NrmlzdGrad./mean(NrmlzdGrad)).^2);
%TrustInd=sqrt(6.*(Rlablty(:,3)./mean(Rlablty(:,3),1)).^2+(Rlablty(:,2)./mean(Rlablty(:,2),1)).^2);

StrArray=cat(2,Coord,AddData,Rlablty,TrustInd,StrLngth); %����������� ���� ���������� ��������� ������ � ���� ������ (�� ������ ��� ������ ������)
StrArray(:,1:6)=StrArray(:,1:6)+SafeCoreFrame;  %���� (��������) ����� ������ �� ������ 
StrArray=sortrows(StrArray,-11);    %���������� ������ �� ����������
StrArray(:,11)=1:cntr;      %������ ������� ���������� �� �������������� ������ (��� ������ ID, ��� ���� ����������)
%���������� ������ � ������
ToStrDltMsk=false(cntr,1);
AngTreshold=1.8*(2*pi/AzStepN); %1.8 %����� ������� �� ������� ��� ������ ������ (+-AngTreshold ������ �� �������� ������� ������� ������ (������� ID))
%����� ������� ������� ������ ��������� ��� ������������ ��� ����� ������������ ������ (������ � ������� ID)
PrllelKTrshold=1/8;	%1/4����� ������� ��������� ������ ������� �� ��� ������������ ������ � ������� ID (� ����� ����� ������� ������)
%AcrossKTrshold=12.0*StrMinLngth.*StrArray(:,12)./StrMaxLngth;	%������ ������� ��������� �� ��� ���������������� ������� ������ (� ��������)
AcrossKTrshold=0.5*StrMinLngth;	%������ ������� ��������� �� ��� ���������������� ������� ������ (� ��������)
    %LngthMltTrshold=2; %���� ����������� ������ ������� ����������� (�������) ����� ��� � LngthMltTrshold ���, �� � ����������� �� �������
    LngthTrshold=0.5;   %������ ������� ����� ������� � ����������� ������ (� ����� ����� ������� ������) - �����-������� � ������������ - �� ������ ��������
    %������� � ������� ������ �� ���� ������:
    tic;
GK=tan(StrArray(:,7));    %����������� ������ ��������� (��������������) ��� ��������� ���� y=kx
CrntIDStrK=tan(StrArray(:,7)-pi/2);%����������� ������� (�������������) ������ ��� ��������� ���� y=kx
SameSourseStrMterick=zeros(cntr,1);                     %������������� ������� ������ ��������� ������
for ID=1:cntr-1
  if ToStrDltMsk(ID)==false
    if StrArray()
      SameSourseStrMterick(ID+1:end)=abs(StrArray(ID+1:end,7)-StrArray(ID,7)); %����������� ������ �������� ������������ ������ � ���� ���������
      SameSourseStrMterick(SameSourseStrMterick>=pi)=abs(SameSourseStrMterick(SameSourseStrMterick>=pi)-2*pi);%���������� ������������ ������� � �������� �� 0 �� 180 ��������
      %GK=tan(StrArray(ID+1:end,7));    %����������� ������ ��������� (��������������) ��� ��������� ���� y=kx
      %CrntIDStrK=tan(StrArray(ID,7)-pi/2);%����������� ������� (�������������) ������ ��� ��������� ���� y=kx
      %���������� ����� ����������� ����� ��������� � ������� ������
      x=(CrntIDStrK(ID)*StrArray(ID,6)-GK(ID+1:end).*StrArray(ID+1:end,6)+StrArray(ID+1:end,5)-StrArray(ID,5))./(CrntIDStrK(ID)-GK(ID+1:end));
      y=CrntIDStrK(ID).*(x-StrArray(ID,6))+StrArray(ID,5); %���������� ����� ����������� ����� ��������� �� ������ ������ � ������� ������
      PrllelDist=sqrt((x-StrArray(ID,6)).^2+(y-StrArray(ID,5)).^2)./StrArray(ID,12);	%���������� (� ����� ����� ������� ������) ����� ������������ ������� ����� ��� ������������ ������������ ������ (���������� �� ����� ����������� ����� ��������� � ������� ������ �� ����������� ����� �������(������������) ������)
      AcrossDist=sqrt((x-StrArray(ID+1:end,6)).^2+(y-StrArray(ID+1:end,5)).^2);    %���������� ���������� (� ��������) ����� ������� (������ � �������) ����� ����� ��������� (�������������� ������ �����).
      LngthDiff=abs(StrArray(ID,12)-StrArray(ID+1:end,12))/StrArray(ID,12);      %������� ����� ������� � ����������� ������ � ����� ������� ������
      %SameSourseStrMterick(ID+1:end)=sqrt((SameSourseStrMterick(ID+1:end)./AngTreshold).^2+(PrllelDist./PrllelKTrshold).^2+(AcrossDist./AcrossKTrshold).^2);  %���� ������������ ������ (�� ������� ��� ������ ���������� ��� ���������� ���������� � ������������ ������ �������) - ��������� ������� � (����������������) ������������ �������_�������� - ������������_���������� - �����������_����������. ��� ��������� ��������������
      SameSourseStrMterick(ID+1:end)=sqrt((SameSourseStrMterick(ID+1:end)./AngTreshold).^2+(PrllelDist./PrllelKTrshold).^2+(AcrossDist./AcrossKTrshold).^2+(LngthDiff/LngthTrshold).^2);  %���� ������������ ������ (�� ������� ��� ������ ���������� ��� ���������� ���������� � ������������ ������ �������) - ��������� ������� � (����������������) ������������ �������_�������� - ������������_���������� - �����������_����������. ��� ��������� ��������������
      %SameSourseStrMterick(ID+1:end)=sqrt((SameSourseStrMterick(ID+1:end)./AngTreshold).^2+(PrllelDist./PrllelKTrshold).^2+(AcrossDist./AcrossKTrshold(ID+1:end)).^2);  %���� ������������ ������ (�� ������� ��� ��� ������ ���������� ��� ���������� ���������� � ������������ ������ �������) - ��������� ������� � (����������������) ������������ �������_�������� - ������������_���������� - �����������_����������. ��� ��������� ��������������
      %ToStrDltMsk(ID+1:end)=ToStrDltMsk(ID+1:end)|SameSourseStrMterick(ID+1:end)<sqrt(3);
      %ToStrDltMsk(ID+1:end)=ToStrDltMsk(ID+1:end)|((SameSourseStrMterick(ID+1:end)<sqrt(3))&(LngthMltTrshold*StrArray(ID,12)>StrArray(ID+1:end,12)));
      ToStrDltMsk(ID+1:end)=ToStrDltMsk(ID+1:end)|SameSourseStrMterick(ID+1:end)<2; %2=sqrt(4) - ��� ������������� ������������
    end
  end
end
StrArray(ToStrDltMsk,:)=[];    %�������� ��������� ��������� ������
Coord=StrArray(:,1:4);  %�������� ������� �������� ���������� ������ � ������� ��������������� �� ���������� ������
AddData=StrArray(:,5:7);
Rlablty=StrArray(:,8:10);
    toc;

imagesc(Data(1:365+60,1:2412+60));
hold on
% axis image
% for i=1:cntr
%     if ToDltMsk(i)==true
%   color=StrArray(i,10);
%   plot(StrArray(i,1:2),StrArray(i,3:4),'Color',[color,color,color]);
%     end
% end
% plot(StrArray(ID,1:2),StrArray(ID,3:4),'red');
% hold off


% imagesc(Data);
% hold on
% axis image
for i=1:size(Coord,1)
    %if Coord(i,2)<=800
      if Rlablty(i,3)>1   %��������� ����������� ������� �������������
        disp('������ ��������� �������� ��������� ���������� ������ - ���������������');
        Rlablty(i,3)=1;
      end
  color=Rlablty(i,3);
  plot(Coord(i,1:2),Coord(i,3:4),'Color',[color,color,color])
    %end
end
hold off
    

end

