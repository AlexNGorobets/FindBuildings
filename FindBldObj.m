function [ ObjCoord, ObjRlablty ] = FindBldObj( StrArray ,Raster, SrsCoord)
%FINDBLDOBJ ���������� ������������� ������� �� ����������� ������� ��������
%ObjCoord ([Row Col] � ������ ������) - ���������� ����������� �������
%ObjRlablty - � ������ ������ - ������������ ������������� ������� (����� ������ ������������� ������ � ������������)
RghtAngleTrshold=15/180*pi; %12 ������ ����� ��������� ���� 90 �������� +-RghtAngleTrshold ������
CornMinDist=9; %9 8������������ ���������� ����� ������ ��������, ������� ����� ���������� ������ ���� � ����� �������� ��� ����������� �� ��� ��������� �������� �������
ComplStrMinDist=8;  %10 8 ������������ ���������� ����� �������� ������� �������� ��������� ������
StrMaxLngth=150;%������������ ����� ������, ������� ����� ���� ������ ��������
difWCoef=0.5;	%0.2� ������ ������ ���������� ������� �� ���, ��� ��������� � ����������� ������� ���� ������ ����� ������ (������ ����� � �������). ���� ������� (����������) �������� ������ �������� (�����������) ��������� �� ������� �� difWCoef �� ��������� �������, �� ����� ������ �������������
TypObjMaxCorn=8;        %�������� ����� ��������� ������������ ��������� ������ � ����������� ����� ����� 1-�� (� �����. ����� 3-� ������). ��� ���� ��������� ����� ��������� ������� � TypObjMaxCorn ����������� ����� (��������� ��������� ��������\��������).
RAngleAccTH=10/180*pi;   %15 10 ������ ����� ������� ��� �����������, ���� ���������������� "������" ���� (����� ������� �������) ���������� �� RAngleAcc �� 90 �������� (��� Treshold-��������)
minW_TH=1./(1.5e+07);     %2.2 e+07������ ����� ������� ��� �����������, ���� �������  �� ���������������� ��� ��� ������ ����� �� ������ ��� minW_TH (��� Treshold-��������)
TypclQStr_TH=1/2.0; %2.5 ��� ������ ����������� �������� � ���������� ������� � ������� ����������� ��������. TresHold-�������� - 1/������� ���������� ������
%     LowBrght_TH=4040;   %4040 ��������� (treshold) �������� ������� �������� � ����� �� �������������� ������� (������ ����������� �� 4 ��������) - ������ ������� "�����" ��������  � ������� ����� �� 0 �� 65535
%     HghBrght_TH=8700;   %8700 ���������� ��������� �������� ������� ������� ��������: 65535-HghBrght_TH
AzCosDev_TH=0.58*mean(1./StrArray(:,10));   %0.58
%���� ������������� ������� ���� �������� 15*15
%���� ������ � ��������� - 5*5
global StrBitMsk RequestRowMsk RequestColMsk;
StrArray(:,7)=StrArray(:,7)-pi/2;%�������������� ����� ��������� � ���� ������ (� ��������� �� -270 �� +90 �������� � ��������)
StrArray(StrArray(:,7)<-pi/2,7)=StrArray(StrArray(:,7)<-pi/2,7)+pi; %�������������� �������� �������� � �������� �� -90 �� +90 �������� (� ��������)
StrArray=sortrows(StrArray,7);  %���������� �������� �� �� �������
%����� ��� ������ � ������������� �������� ����� ������������� �� ����� ������ � ������������� �������� ��� ��������� ������ ����� (� ��������)
AzFrdgeM=find(StrArray(:,7)>-RghtAngleTrshold,1);	%����� ������ ������ � ������, ������ ������ ������ 0-����� �������� - ��������� ����� ������ ���������� (Master) ������ (�����.������, k>0,������ ��������)
AzFrdgeS=find(StrArray(:,7)>RghtAngleTrshold,1);	%����� ������ ������ � ������, ������ ������ ������ 0+����� �������� - �������� ����� ������ ������������� (Slave) ������ (�����.������, k<0, ������ ��������)
RghtAngleObjCandidates=zeros(size(StrArray,1)-AzFrdgeM+1,AzFrdgeS-1); %������������� ������� Master-��������
StrBitMsk=false(size(RghtAngleObjCandidates));    %������������� ����� ������ ��� ��������, ������� ���������� ������ ���� �������
    tic;
for i=1:AzFrdgeS-1      %�������� Master-������� �� ���������� Slave-������ ���
  RghtAngleObjCandidates(:,i)=StrArray(AzFrdgeM:end,7);  % ����� ������ - ���������� ����� Master-������� (������ ������� - ��� Master �������)
end
for i=1:size(StrArray,1)-AzFrdgeM+1; %�������� Slave-������� Master- ���
  RghtAngleObjCandidates(i,:)=RghtAngleObjCandidates(i,:)-StrArray(1:AzFrdgeS-1,7)';    % ����� ������� - ���������� ����� Slave-������� (������ ������ ��� Slave-�)
end
RghtAngleObjCandidates=abs(RghtAngleObjCandidates-pi/2);  %�������� ������� ������ �������� ������� Master � ������ Slave-�������� �� 90�������� (� ��������)
StrBitMsk=RghtAngleObjCandidates<=RghtAngleTrshold;	%����� ��� ������ - ���������� � ������� (�� �������� ������� ������� ���� � �������� RghtAngleTrshold)
MstrStrData=StrArray(AzFrdgeM:end,:);   %����������� ��������������� ���������� �������� � ������ Master-������
SlveStrData=StrArray(1:AzFrdgeS-1,:);   %����������� ����������  � ������ Slave-������
% MstrGrW=MstrStrData(:,8)./MstrStrData(:,9);     %����=������������/���-�� �������� ������� �������
% SlveGrW=SlveStrData(:,8)./SlveStrData(:,9);
MstrGrW=MstrStrData(:,8)./MstrStrData(:,12);     %����=������������/���-�� �������� ������� �������
SlveGrW=SlveStrData(:,8)./SlveStrData(:,12);
MstrK=tan(MstrStrData(:,7));    %������������ ��� ������� ������� ��� ��-� y=kx+b
SlveK=tan(SlveStrData(:,7));
MSStrSolvesX=zeros(1,AzFrdgeS-1);   %������������� ������� ����� ����������� ���� Slave- � ������� Master-������ (�, �������)
MSStrSolvesY=zeros(1,AzFrdgeS-1);   % ���������� ��� (Y, ������)
MDist=zeros(2,AzFrdgeS-1);          %������������� ���������� �� ����� ����������� �� ���� Master-�������
SDist=MDist;                        % ���������� ��� Slave-�������
ShrtstMDist=zeros(1,AzFrdgeS-1);    %������������� ���������� �� ������ ����� ����������� �� ���������� �� ���� ���� Master-�������
ShrtstSDist=ShrtstMDist;            % ���������� ��� ���������� ���� ���������������� Slave-�������
        %imshow(Raster(1:365+60,1:2412+60));
        figure; 
            TbleSize=0;
            %TbleSize=TbleSize*size(Raster,1);
            %Tble=zeros(TbleSize,size(Raster,2));
        %imshow(cat(1,Tble,Raster,Tble));
        imshow(Raster);
        axis on,axis equal,hold on;
        %SrsCoord=[11025 31;2420 370;3444 66;11388 352]; %���������� ��������� ������� ��� 5-�� �����������
    tic
StartRow=find(any(StrBitMsk,2),1);
FinshRow=size(StrArray,1)-AzFrdgeM+1;
    clearvars StrArray RghtAngleObjCandidates
for i=StartRow:FinshRow %������� Master-�������� � ����� �� ������� ���������� slave-��������
    MSStrSolvesX(StrBitMsk(i,:))=(MstrK(i)*MstrStrData(i,6)-MstrStrData(i,5)-SlveK(StrBitMsk(i,:)).*SlveStrData(StrBitMsk(i,:),6)+SlveStrData(StrBitMsk(i,:),5))./(MstrK(i)-SlveK(StrBitMsk(i,:)));
    MSStrSolvesY(StrBitMsk(i,:))=MstrK(i).*(MSStrSolvesX(StrBitMsk(i,:))-MstrStrData(i,6))+MstrStrData(i,5);    %����� ����������� ������� Master-������ � ���� Slave-������ (X � Y - ��������������)
    MDist(1,StrBitMsk(i,:))=sqrt((MSStrSolvesX(StrBitMsk(i,:))-MstrStrData(i,1)).^2+(MSStrSolvesY(StrBitMsk(i,:))-MstrStrData(i,3)).^2);    %���������� �� ������ ����� ����������� �� ������� �� ���� ������ Master-�������
    MDist(2,StrBitMsk(i,:))=sqrt((MSStrSolvesX(StrBitMsk(i,:))-MstrStrData(i,2)).^2+(MSStrSolvesY(StrBitMsk(i,:))-MstrStrData(i,4)).^2);    
    ShrtstMDist(StrBitMsk(i,:))=min(MDist(:,StrBitMsk(i,:)),[],1);    %���������� �� ������ ����� ����������� �� ���������� ����� Master-�������
    OnLineMsk=(MSStrSolvesX(StrBitMsk(i,:))<max(MstrStrData(i,1:2),[],2))&(MSStrSolvesX(StrBitMsk(i,:))>min(MstrStrData(i,1:2),[],2));  %����� ShrtstSDist(StrBitMsk(i,:)-� �������� ����� ����������� ������� ����� �� ��������������� Master-�������. ����� ���� ����� ����� ���������� ������ � ������� ������ StrBitMsk
    ShrtstMDist(StrBitMsk(i,:))=ShrtstMDist(StrBitMsk(i,:)).*~OnLineMsk;    %%���� ����� ����������� ����� � �������� ������� - ���������� ���������� �� 0 (������� ������� ���� �� X - ������ �������)
    SDist(1,StrBitMsk(i,:))=sqrt((MSStrSolvesX(StrBitMsk(i,:))-SlveStrData(StrBitMsk(i,:),1)').^2+(MSStrSolvesY(StrBitMsk(i,:))-SlveStrData(StrBitMsk(i,:),3)').^2);% ���������� ���������� �� ������ ����� ����������� �� ������� �� ���� ������ ���������������� Slave-�������.
    SDist(2,StrBitMsk(i,:))=sqrt((MSStrSolvesX(StrBitMsk(i,:))-SlveStrData(StrBitMsk(i,:),2)').^2+(MSStrSolvesY(StrBitMsk(i,:))-SlveStrData(StrBitMsk(i,:),4)').^2);
    ShrtstSDist(StrBitMsk(i,:))=min(SDist(:,StrBitMsk(i,:)),[],1);    %���������� �� ������ ����� ����������� �� ���������� ����� ���������������� Slave-�������
    OnLineMsk=(MSStrSolvesX(StrBitMsk(i,:))<max(SlveStrData(StrBitMsk(i,:),1:2),[],2)')&(MSStrSolvesX(StrBitMsk(i,:))>min(SlveStrData(StrBitMsk(i,:),1:2),[],2)');  %����� ShrtstSDist(StrBitMsk(i,:)-� �������� ����� ����������� ������� ����� �� ��������������� Slave-�������. ����� ���� ����� ����� ���������� ������ � ������� ������ StrBitMsk
    ShrtstSDist(StrBitMsk(i,:))=ShrtstSDist(StrBitMsk(i,:)).*~OnLineMsk;    %���� ����� ����������� ����� � �������� ������� - ���������� ���������� �� 0 (������� ������� ���� �� X - ������ �������)
    StrBitMsk(i,StrBitMsk(i,:))=(ShrtstMDist(StrBitMsk(i,:))<CornMinDist & ShrtstSDist(StrBitMsk(i,:))<CornMinDist);%�������� �� ����� ��� ������ (Master-Slave) ����� ����������� ������� ������� ������ �� ��������
end
    toc;
%��������� ����� ��� ��������� ��������.
global MstrStrLinkList SlveStrLinkList;
MstrStrLinkList=zeros(5000,2);  %������������� ������� ������ Master-��������
SlveStrLinkList=zeros(5000,2);  % � Slave-�������� ��������������
StrCntr=1;  %�������� �������� �������
StrLinkCoord=zeros(size(MstrStrData,1),4);  %������������� ������� ���������� �� ������ ����� ���������� ������ �� ������ ����� ��������� ������
    tic
%����� ������ ����� Master-��������:
for i=StartRow:FinshRow     %�������� ��� Master-��������. � ������ ����������� �� ������ Master-�������� ���������� ������� �� ������� (� �������� RghtAngleTrshold)
AdrA=i+1;   %��������� ����� ��������� �������� � MstrStrData (������� �����������). ��� ���� ����� ������ ����� ���� ���������� �������� (������� �� ��������� ���� �� �����), ������� ������� ������������ � ������� ���������, ������� � ������� ��������.
%AdrA=find(MstrStrData(:,7)>MstrStrData(i,7)-RghtAngleTrshold,1);
AdrB=find(MstrStrData(:,7)>MstrStrData(i,7)+RghtAngleTrshold,1);    %��������� ����� ��������� ������ �������� ��������� �������� � MstrStrData (������� �����������).
if ((AdrB-AdrA)>0)  %���������� ��������� ������� ���������
  AdrB=AdrB-1;      %AdrB - ������ ����� ���������� �������� ���������
  StrLinkCoord(AdrA:AdrB,1)=sqrt((MstrStrData(i,1)-MstrStrData(AdrA:AdrB,1)).^2+(MstrStrData(i,3)-MstrStrData(AdrA:AdrB,3)).^2);%���������� �� ������ ����� �������� ������� �� ������ ����� ������� �������
  StrLinkCoord(AdrA:AdrB,2)=sqrt((MstrStrData(i,2)-MstrStrData(AdrA:AdrB,2)).^2+(MstrStrData(i,4)-MstrStrData(AdrA:AdrB,4)).^2);%���������� �� ������ ����� �������� ������� �� ������ ����� ������� �������
  StrLinkCoord(AdrA:AdrB,3)=sqrt((MstrStrData(i,1)-MstrStrData(AdrA:AdrB,2)).^2+(MstrStrData(i,3)-MstrStrData(AdrA:AdrB,4)).^2);%���������� �� ������ ����� �������� ������� �� ������ ����� ������� �������
  StrLinkCoord(AdrA:AdrB,4)=sqrt((MstrStrData(i,2)-MstrStrData(AdrA:AdrB,1)).^2+(MstrStrData(i,4)-MstrStrData(AdrA:AdrB,3)).^2);%���������� �� ������ ����� �������� ������� �� ������ ����� ������� �������
  ResAdr=AdrA-1+find((min(StrLinkCoord(AdrA:AdrB,:),[],2)<ComplStrMinDist)&(max(StrLinkCoord(AdrA:AdrB,:),[],2)>(MstrStrData(i,12)+MstrStrData(AdrA:AdrB,12)-ComplStrMinDist)));%������ �����������: ����� �������(��) ���� ����� ������� ��������� � ��������� �� CornMinDist, � ������ ����� ��������� �� ���������� ���� �� ����� ������ ����� ����� ���� �������� � ��������� �� CornMinDist
  if ~isempty(ResAdr)   %����������� ��������� �������� (���� ����� ������� �� ������� ��������)
    MstrStrLinkList(StrCntr:StrCntr+size(ResAdr,1)-1,1)=i;  %����� ����������� ������� � MstrStrData
    MstrStrLinkList(StrCntr:StrCntr+size(ResAdr,1)-1,2)=ResAdr; %����� ���������� ������� � MstrStrData
    StrCntr=StrCntr+size(ResAdr,1); %���������� �������� (� ��������� ������� ������ ������ ��������� ��� ��������� ��������)
  end
end
end
MstrStrLinkList(StrCntr:end,:)=[];  %�������� ���������� (���������) ����� ������� ������� ��������� Master-��������
%         for i=1:size(MstrStrLinkList,1) %��������� ��������� ��������
%             plot(MstrStrData(MstrStrLinkList(i,1),1:2),MstrStrData(MstrStrLinkList(i,1),3:4),'color','blue');
%             plot(MstrStrData(MstrStrLinkList(i,2),1:2),MstrStrData(MstrStrLinkList(i,2),3:4),'color','red');
%         end
%����� ������ ����� Slave-��������:  
StrCntr=1;  %����� ��������� �������
StrLinkCoord=zeros(size(SlveStrData,1),4);  %����������������� ������� ���������� (�� ���������� Slave-��������)
for i=1:size(StrBitMsk,2)	%�������� ��� Slave-��������. � ������ ����������� ���������� ���������� ������� �� ������� �������
AdrA=i+1;   %��������� ����� ��������� �������� � SlveStrData (������� �����������).
%AdrA=find(SlveStrData(:,7)>SlveStrData(i,7)-RghtAngleTrshold,1);
AdrB=find(SlveStrData(:,7)>SlveStrData(i,7)+RghtAngleTrshold,1);%��������� ����� ��������� ������ �������� ��������� �������� � SlveStrData (������� �����������).
if ((AdrB-AdrA)>0)  %���������� ��������� ������� ���������
  AdrB=AdrB-1;      %AdrB - ������ ����� ���������� �������� ���������
  StrLinkCoord(AdrA:AdrB,1)=sqrt((SlveStrData(i,1)-SlveStrData(AdrA:AdrB,1)).^2+(SlveStrData(i,3)-SlveStrData(AdrA:AdrB,3)).^2);%���������� �� ������ ����� �������� ������� �� ������ ����� ������� �������
  StrLinkCoord(AdrA:AdrB,2)=sqrt((SlveStrData(i,2)-SlveStrData(AdrA:AdrB,2)).^2+(SlveStrData(i,4)-SlveStrData(AdrA:AdrB,4)).^2);%���������� �� ������ ����� �������� ������� �� ������ ����� ������� �������
  StrLinkCoord(AdrA:AdrB,3)=sqrt((SlveStrData(i,1)-SlveStrData(AdrA:AdrB,2)).^2+(SlveStrData(i,3)-SlveStrData(AdrA:AdrB,4)).^2);%���������� �� ������ ����� �������� ������� �� ������ ����� ������� �������
  StrLinkCoord(AdrA:AdrB,4)=sqrt((SlveStrData(i,2)-SlveStrData(AdrA:AdrB,1)).^2+(SlveStrData(i,4)-SlveStrData(AdrA:AdrB,3)).^2);%���������� �� ������ ����� �������� ������� �� ������ ����� ������� �������
  ResAdr=AdrA-1+find((min(StrLinkCoord(AdrA:AdrB,:),[],2)<ComplStrMinDist)&(max(StrLinkCoord(AdrA:AdrB,:),[],2)>(SlveStrData(i,12)+SlveStrData(AdrA:AdrB,12)-ComplStrMinDist)));%������ �����������: ����� �������(��) ���� ����� ������� ��������� � ��������� �� CornMinDist, � ������ ����� ��������� �� ���������� ���� �� ����� ������ ����� ����� ���� �������� � ��������� �� CornMinDist
  if ~isempty(ResAdr)   %����������� ��������� �������� (���� ����� ������� �� ������� ��������)
    SlveStrLinkList(StrCntr:StrCntr+size(ResAdr,1)-1,1)=i;  %����� ����������� ������� � MstrStrData
    SlveStrLinkList(StrCntr:StrCntr+size(ResAdr,1)-1,2)=ResAdr; %����� ���������� ������� � MstrStrData
    StrCntr=StrCntr+size(ResAdr,1); %���������� �������� (� ��������� ������� ������ ������ ��������� ��� ��������� ��������)
  end
end
end
SlveStrLinkList(StrCntr:end,:)=[];  %�������� ���������� (���������) ����� ������� ������� ��������� Master-��������
toc;
        for i=1:size(SlveStrLinkList,1) %��������� ��������� ��������
            plot(SlveStrData(SlveStrLinkList(i,1),1:2),SlveStrData(SlveStrLinkList(i,1),3:4),'color','yellow');
            plot(SlveStrData(SlveStrLinkList(i,2),1:2),SlveStrData(SlveStrLinkList(i,2),3:4),'color','white');
        end
        StrCntr=0;
%��������� "���������" ����� ����� ������ ��������� ������. ���� � ���� ��� ���� ����� ������ - ��� �������.
% �� ���� ��� ����� (true-����) � ������� StrBitMsk- ���� ������ (�� �� � ��� �����), �.�.� - ����������  
RequestRowMsk=any(StrBitMsk,2); %� ������ �����������, ����� ������� ��� ����� �������� ��� �� ������� ������
RequestColMsk=any(StrBitMsk,1); %���������� - ��� �� ������� ������
TtlStrQ=sum(sum(StrBitMsk));    %���������� �������� ������ ����� �� ������ (���������� ��� ����������� ������� ������, ������� �����������)
Cntrs=zeros(TypObjMaxCorn,1);   %������� ��������� �������� � ������ ������
ObjArr=cell(fix(TtlStrQ/2),TypObjMaxCorn);  %Cell-������ �������� ��� 2-�, 3-�, 4-� � ����� �������� �������� ��������������. (����� ��� ���� � ����� ������)
    NoDltMstrBitMsk=false(size(StrBitMsk,1),1); %������ ������� �� ��������: � ����� ���� ��������� �������� ��������� ������ ������ � ������������� ��� �� ���������������� ���������� �� ������ (�������), �������� �������� ������ ������� ("�������� ����� ������") �������, ������� ��������� � ����������� ������ ������ ������� ���� � �������.
    NoDltSlveBitMsk=false(size(StrBitMsk,2),1); % ����� ������ ���������� ����� �������, ��� ��� ����������� � ������ �������� �� ������� ������ ���� ���. ������� ���������� ��������� ������ ����� ����������� ���� ������ ���� ���, ������� ��� ��� �������� ��������� �������� � �������� ����������� ��������� ������
        tic;
            DltCCntr=0;
for i=find(RequestRowMsk,1):FinshRow %������� Master-�������� � ����� ���� ��������� �������� (�� ���, ��� �������� � �.� ������ ������ ����)
  if RequestRowMsk(i)==true     %�������� �������� �� ������� ������ StrBitMsk ������� ������� (�.�. � ��� ���� ����� � ��� ��� �� ���� ���������)
    [Q,Coords]=FndColByColLinks(i); %��������� ����� ����������� ������ (����� � StrBitMsk, ������� ����� � ����� ������ �� ������� � ������� ������� ������ ����������)
    %����������� ��������� ��������:
    for s=1:size(Coords,1)  %������� ��� Master- �  Slave-��������
      %StrLinks=find(MstrStrLinkList(:,1)==Coords(s,1)); %����� ������ �� �������� �aster-�������
      StrLinks=RecuFindM(Coords(s,1));
      for t=1:size(StrLinks,1)  %������� ��������� ������ �� �������� �-�������
        if(RequestRowMsk(MstrStrLinkList(StrLinks(t),2))==true)%��������� ��������� ��������, ������� ��� �� ���� ��������� �� ������� ����������� ��������.
            [tQ,tCoords]=FndColByColLinks(MstrStrLinkList(StrLinks(t),2)); %����� ����������� ������ � ������������ ������ ���������� Master-�������
            Q=Q+tQ;     %������� ��������: ���������� �������� �������� � ������� 
            Coords=cat(1,Coords,tCoords);   %�������� ������� �������� ���� �������� � ���� ������.
            NoDltMstrBitMsk(MstrStrLinkList(StrLinks(t),:))=true;%������ �� �������� �������� ���������� ��������� ������� � ������������� ����� ��������� ������ ���������� (���� ����� ������� �� ��� � ����������� � ������ ��������)
                StrCntr=StrCntr+1;
        end
      end
      %StrLinks=find(SlveStrLinkList(:,1)==Coords(s,2)); %����� ������ �� �������� Slave-�������
      StrLinks=RecuFindS(Coords(s,2));
      for t=1:size(StrLinks,1)  %������� ��������� ������ �� �������� S-�������
        if(RequestColMsk(SlveStrLinkList(StrLinks(t),2))==true)%��������� ��������� ��������, ������� ��� �� ���� ��������� �� ������� ����������� ��������.
            [tQ,tCoords]=FndRowByRowLinks(SlveStrLinkList(StrLinks(t),2)); %����� ����������� ������ � ������������ ������ ���������� Master-�������
            Q=Q+tQ;     %������� ��������: ���������� �������� �������� � ������� 
            Coords=cat(1,Coords,tCoords);   %�������� ������� �������� ���� �������� � ���� ������.
            NoDltSlveBitMsk(SlveStrLinkList(StrLinks(t),:))=true;%������ �� �������� �������� ���������� ��������� ������� � ������������� ����� ��������� ������ ���������� (���� ����� ������� �� ��� � ����������� � ������ ��������)
                StrCntr=StrCntr+1;
        end
      end
      
    end

    if Q>1  %������� ������ ������� �������� ������� ������ ��������� �� difWCoef (���� ��� ������ �������� ������ ���� ������ ����� ������). ������ ����� �������������� ���� ����� ������ ��� ����
      CDltMsk=false(size(Coords,1),1);  %������������ ����� ������� ��� ������
      UMC=unique(Coords(:,1));    %Unique Master Coorinates - ������ Master-������ ������� ��� ����������
      USC=unique(Coords(:,2));    %Unique Slave Coorinates - ������ Slave-������ ������� ��� ����������
      MeanObjGrW=(sum(MstrStrData(UMC,8))+sum(SlveStrData(USC,8)))./(sum(MstrStrData(UMC,12))+sum(SlveStrData(USC,12)));%���������� ����������������� ���� �� ������ (����� ����������/����� ��������)
      for t=1:size(UMC,1) 	%������� Master-��������
        if (~NoDltMstrBitMsk(UMC(t)))&&((MeanObjGrW-MstrGrW(UMC(t)))>(difWCoef*MeanObjGrW))&&(sum(sum(find(Coords(:,1)==UMC(t))))==1)   %������� ����������� �� ���������� ���� & ������� ������� ������ � ����������� ������ ������ ����
          CDltMsk(t)=true;  %������� �������-��������� �� �������� �� ������ �������� �������
            DltCCntr=DltCCntr+1;
            plot(MstrStrData(UMC(t),1:2),MstrStrData(UMC(t),3:4)+TbleSize,'Color','cyan');
        end
      end
      for t=1:size(USC,1)   %������� Slave-��������
        if (~NoDltSlveBitMsk(USC(t)))&&((MeanObjGrW-SlveGrW(USC(t)))>(difWCoef*MeanObjGrW))&&(sum(sum(find(Coords(:,2)==USC(t))))==1)
          CDltMsk(t)=true;  %������� �������-��������� �� �������� �� ������ �������� �������
            DltCCntr=DltCCntr+1;
            plot(SlveStrData(USC(t),1:2),SlveStrData(USC(t),3:4)+TbleSize,'Color','cyan');
        end
      end
      Coords(CDltMsk,:)=[]; %�������� "������� �������" ��������
      Q=Q-sum(CDltMsk);     %��������������� ����� �� ���������� �����
      if Q>1    %���� �������. �������� ��������� �� ������� ���� ��� ������ � ����� (�����) �������������� ������
        if Q<=TypObjMaxCorn %���� �������� ���������������� ��������
          Cntrs(Q-1)=Cntrs(Q-1)+1;      %��������� ���������������� ��������
          ObjArr{Cntrs(Q-1),Q-1}=Coords;%������ ��������� � ������ ��� ���������� ���������
        else                %���� ������� (������) ��������
          Cntrs(TypObjMaxCorn)=Cntrs(TypObjMaxCorn)+1;      %��������� ���������������� (����������) ��������
          ObjArr{Cntrs(TypObjMaxCorn),TypObjMaxCorn}=Coords;%������ ��������� � ������ ��� ���������� ���������
        end
            Mm=unique(Coords(:,1)); %���������
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
ObjArr(max(Cntrs)+1:end,:)=[];  %�������� �������� ����� ������� ��������
%���� ������ �������� ������ ��������
%CrntObjRlablty=zeros(size(ObjArr)); %������������ ������� ������������ ������������� �������
%ObjAngle=zeros(size(ObjArr));
ObjCoord=zeros(sum(Cntrs),2);   %������������� �������� ��������
ObjRlablty=zeros(sum(Cntrs),1);
%������� �������� �������: ������� ���� (������ ���� ����� � 90 ��������). ������� �������� �� ������ ������ ������ ���� �� ����� minW_TH
resObjcntr=1;       %resObjcntr-1 - ���������j ����������� �������� (resObjcntr - ����� ���������� �������, ���� ������ �� �����������, �� � ��������� �������� ���������� �� ����� ������ ����������)
for type=1:TypObjMaxCorn-1  %������� ������� (�������) �������� - �������� ����� (�� 2-� ����� �� TypObjMaxCorn)
  Angles=zeros(type+1,1);   %(�������������) ������ ����� ��� ����� ���������� ���� ������, ������� �������
  GradW=zeros(type+1,2);    %(�������������) ������ ����� ������ ������ (��� - ��� ��������� �������� (������ �������� �� ������� ���� ��������� ������) ��������� � ���������� �������� (� ���� ������))
  for o=1:Cntrs(type)       %������� ���� �������� ������� (��������) ����
    Angles=abs(MstrStrData(ObjArr{o,type}(:,1),7)-SlveStrData(ObjArr{o,type}(:,2),7));  %���� ��� ������ ���� ������
    GradW(:,1)=MstrGrW(ObjArr{o,type}(:,1));    %����������� ����� � ������ ����� �������
    GradW(:,2)=SlveGrW(ObjArr{o,type}(:,2));
    CrntObjRlablty=sum(unique(GradW(:,1)))+sum(unique(GradW(:,2))); %����� ����� ���� ������ (������������� ���������)
    GradW=GradW./CrntObjRlablty;    %����������������� ���� ������, ������ �� ��������� ���������� ��������� �������� �� ������ � ���������� ���� ��������
    NrmL=min([MstrStrData(ObjArr{o,type}(:,1),12) SlveStrData(ObjArr{o,type}(:,2),12)],[],2);   %���� ������ �� �� ������ 
    NrmL=NrmL./sum(NrmL);           %��������������� ���� (�� ������)
    ObjAngle=sum(abs(Angles-pi/2).*NrmL);   %���������������� (�� ������ �������� ��������) ������� ������ ���� �� �������
    ObjCoord(resObjcntr,1)=sum(unique(MstrStrData(ObjArr{o,type}(:,1),5).*GradW(:,1)))+sum(unique(SlveStrData(ObjArr{o,type}(:,2),5).*GradW(:,2)));   %���������� ��������� ������ �������
    ObjCoord(resObjcntr,2)=sum(unique(MstrStrData(ObjArr{o,type}(:,1),6).*GradW(:,1)))+sum(unique(SlveStrData(ObjArr{o,type}(:,2),6).*GradW(:,2)));
        %ObjColor=mean(mean(Raster(fix(ObjCoord(resObjcntr,1)):ceil(ObjCoord(resObjcntr,1)),fix(ObjCoord(resObjcntr,2)):ceil(ObjCoord(resObjcntr,2)))));   %���������� ����� ������� (�� ������ ��������) � ����� ���������� ��������� ��������������� ������
        %NrmlzdObjColorDiff=min([ObjColor/LowBrght_TH (65535-ObjColor)/HghBrght_TH]);	%���������������� ������� ������� ����� ������� � ���������� ���������� �������� (0 ��� 65535)
        NrmL=([MstrStrData(ObjArr{o,type}(:,1),12) SlveStrData(ObjArr{o,type}(:,2),12)])./(sum(unique(MstrStrData(ObjArr{o,type}(:,1),12)))+sum(unique(SlveStrData(ObjArr{o,type}(:,2),12))));  %������������� ����� (��� �������������)
        NrmlzdAzCosDev=sum(unique(NrmL(:,1)./MstrStrData(ObjArr{o,type}(:,1),10)))+sum(unique(NrmL(:,2)./SlveStrData(ObjArr{o,type}(:,2),10)))./(AzCosDev_TH);
    %CrntObjRlablty=sqrt((ObjAngle./RAngleAccTH).^2+(1./CrntObjRlablty./(minW_TH*(type+2))).^2+(1/(type+2)/TypclQStr_TH).^2+NrmlzdObjColorDiff.^2); %������������� �������� �������
    CrntObjRlablty=sqrt((ObjAngle./RAngleAccTH).^2+(1./CrntObjRlablty./(minW_TH*(type+2))).^2+(1/(type+2)/TypclQStr_TH).^2+NrmlzdAzCosDev.^2); %������������� �������� �������
    if CrntObjRlablty<sqrt(4)  %ObjAngle(o,type)<RAngleAccTH
      ObjRlablty(resObjcntr)=CrntObjRlablty;    %����������� ��������� ������������� � �������� ������
            plot(ObjCoord(resObjcntr,2),ObjCoord(resObjcntr,1)+TbleSize,'gx')            
            Mm=unique(ObjArr{o,type}(:,1));
            Ms=unique(ObjArr{o,type}(:,2));
            for m=1:size(Mm,1)
             plot(MstrStrData(Mm(m),1:2),MstrStrData(Mm(m),3:4)+TbleSize,'Color','red');
            end
            for s=1:size(Ms,1)
             plot(SlveStrData(Ms(s),1:2),SlveStrData(Ms(s),3:4)+TbleSize,'Color','red');
            end
      resObjcntr=resObjcntr+1;  %��������� �������� ���������� ��������      
    end
  end
  
    %ObjRlablty=sqrt((ObjAngle(1:Cntrs(type),type)./RAngleAccTH).^2+(ObjRlablty(1:Cntrs(type),type)./(minW_TH*(type+2))).^2);
    %ApprvmnObjtMsk=ObjRlablty<sqrt(2);
end
        %plot(SrsCoord(:,1),SrsCoord(:,2)+TbleSize,'mo');
%���� ����� ������� ��������
type=TypObjMaxCorn; %������ ������� ��� - ��������� (����� TypObjMaxCorn �����)
S=zeros(size(ObjArr,1),1);  %������������� ������� ��������� �������� � ������� ��������
for o=1:size(ObjArr,1)    %������� ������� ��������
	S(o)=size(ObjArr{o,type},1);   %��������� �������� ��������� ��������
end
LstDfcltObj=find(S==0,1)-1; %���������� ���������� ������� ��������
QStr=max(S);                %���������� ������������� ���������� ������ � ����� ������� �������
Angles=zeros(QStr,1);   %������������� ������� ����
GradW=zeros(QStr,2);    % � �����
for o=1:LstDfcltObj     %������� ������� ��������
    Angles(1:S(o))=abs(MstrStrData(ObjArr{o,type}(:,1),7)-SlveStrData(ObjArr{o,type}(:,2),7));  %���� ��� ������ ���� ������
    GradW(1:S(o),1)=MstrGrW(ObjArr{o,type}(:,1));
    GradW(1:S(o),2)=SlveGrW(ObjArr{o,type}(:,2));
    CrntObjRlablty=sum(unique(GradW(1:S(o),1)))+sum(unique(GradW(1:S(o),2))); %����� ����� ���� ������ (������������� ���������)
    GradW=GradW./CrntObjRlablty;    %���������������� ���� ������, ������ �� ��������� ���������� ��������� �������� �� ������ � ���������� ���� ��������
    NrmL=min([MstrStrData(ObjArr{o,type}(:,1),12) SlveStrData(ObjArr{o,type}(:,2),12)],[],2);   %���� ������ �� �� ������ (��� ������ ���� ���������� �������� �������)
    NrmL=NrmL./sum(NrmL);	%���������������� ���� (�� ������)
    ObjAngle=sum(abs(Angles(1:S(o))-pi/2).*NrmL);   %���������������� (�� ������ �������� ��������) ������� ������ ���� �� �������
    ObjCoord(resObjcntr,1)=sum(unique(MstrStrData(ObjArr{o,type}(:,1),5).*GradW(1:S(o),1)))+sum(unique(SlveStrData(ObjArr{o,type}(:,2),5).*GradW(1:S(o),2)));   %���������� ��������� ������ �������
    ObjCoord(resObjcntr,2)=sum(unique(MstrStrData(ObjArr{o,type}(:,1),6).*GradW(1:S(o),1)))+sum(unique(SlveStrData(ObjArr{o,type}(:,2),6).*GradW(1:S(o),2)));  
    %ObjColor=mean(mean(Raster(fix(ObjCoord(resObjcntr,1)):ceil(ObjCoord(resObjcntr,1)),fix(ObjCoord(resObjcntr,2)):ceil(ObjCoord(resObjcntr,2)))));   %���������� ����� ������� (�� ������ ��������) � ����� ���������� ��������� ��������������� ������
        NrmL=([MstrStrData(ObjArr{o,type}(:,1),12) SlveStrData(ObjArr{o,type}(:,2),12)])./(sum(unique(MstrStrData(ObjArr{o,type}(:,1),12)))+sum(unique(SlveStrData(ObjArr{o,type}(:,2),12))));  %������������� ����� (��� �������������)
        NrmlzdAzCosDev=sum(unique(NrmL(:,1)./MstrStrData(ObjArr{o,type}(:,1),10)))+sum(unique(NrmL(:,2)./SlveStrData(ObjArr{o,type}(:,2),10)))./(AzCosDev_TH);
    %NrmlzdObjColorDiff=min([ObjColor/LowBrght_TH (65535-ObjColor)/HghBrght_TH]);	%���������������� ������� ������� ����� ������� ������� � ���������� ���������� �������� (0 ��� 65535)
    %CrntObjRlablty=sqrt((ObjAngle./RAngleAccTH).^2+(1./CrntObjRlablty./(minW_TH*S(o))).^2+(1/S(o)/TypclQStr_TH).^2+NrmlzdObjColorDiff.^2);
    CrntObjRlablty=sqrt((ObjAngle./RAngleAccTH).^2+(1./CrntObjRlablty./(minW_TH*S(o))).^2+(1/S(o)/TypclQStr_TH).^2+NrmlzdAzCosDev.^2);
    if CrntObjRlablty<sqrt(4)  %ObjAngle(o,type)<RAngleAccTH
      ObjRlablty(resObjcntr)=CrntObjRlablty;    %����������� ��������� ������������� � �������� ������
            plot(ObjCoord(resObjcntr,2),ObjCoord(resObjcntr,1)+TbleSize,'gx')            
            Mm=unique(ObjArr{o,type}(:,1));
            Ms=unique(ObjArr{o,type}(:,2));
            for m=1:size(Mm,1)
             plot(MstrStrData(Mm(m),1:2),MstrStrData(Mm(m),3:4)+TbleSize,'Color','red');
            end
            for s=1:size(Ms,1)
             plot(SlveStrData(Ms(s),1:2),SlveStrData(Ms(s),3:4)+TbleSize,'Color','red');
            end
      resObjcntr=resObjcntr+1;  %��������� �������� ���������� ��������
    end
end
ObjCoord(resObjcntr:end,:)=[];    %�������� ������ ����� ��������
ObjRlablty(resObjcntr:end)=[];
    resObjcntr
    axis equal;
end
function [Ccntr,Coords]=FndColByColLinks(SrsRow)%StrtCol
%SrsRow - ����� �������� ������ � StrBitMsk. � ��� ������ ����� (true) � �������� ���� �������.
%StrtCol - ����� ������� � ������� �� ���� ������ (SrsRow) ����������� �����(true-���) ������� � ������������ 
% ����� ������ ����� � ���� ������. ���� ��������� ����� ���, �� StrtCol=0
%Coords - [Row Col] ������ ����� � ��������, ������� ����������� � �������� �������
% �������������� � �������� ������. ������ ����� ����� � ����� ������ Coords
%Ccntr - ���������� ��������� ����� ��� S���� ��������� (��������) ��� ������
global StrBitMsk;
global RequestColMsk;
global RequestRowMsk;
[~,c]=find(StrBitMsk(SrsRow,:));%������ ���� �� ������� ��������� � �������� ������ ��������
c=c.*RequestColMsk(c);  %�������� �� ������ �������� ����� (� ��������) ���, ��� ��� ������ (� ��������� ���, ��-�� ������� ��� ����������� ����� � ���� ������)
c(c==0)=[]; %�������� �� ������ �������� ����� (� ��������) ��� �� ������� ��� ������������ �����
Csize=size(c,2);    %���������� �������� ����� (true-�����) � ������ ������
RequestRowMsk(SrsRow)=false;%������� ������ ��� ������������
if Csize==0 %��������� ������� �� �������, ���� �� ����� ����� (����� ���������) ������� �� ����
  Coords=[];
  Ccntr=0;    
else
  Coords=zeros(Csize*3,2);  %������������� ��������� �������
  Coords(1:Csize,1)=SrsRow; %������ Row-����������(�������) ���� �������� ����� � ������� ������
  Coords(1:Csize,2)=c;      %������ Col-��������� (������� ��������) ���� �����
  Ccntr=Csize;  %������������� ������� ������� ����� (true-��� ��� ��� ������, ������� �������� ������� � ������� ������ ����������)
  %����� ��������� ������(����� � StrBitMsk) � �������� ��� ���� ������� �������� �������� �������
  for i=1:Csize     %������� ��������, �� ������� ������� ����� � ������� ������
    [s,Coords(Ccntr+1:Ccntr+s,:)]=FndRowByRowLinks(c(i));	%,SrsRow ����������� ��������� ����� ����������� � ����� ������� � ������ �� ��������� (� ������� ������) ����� (� ����� ����������� ��������� � ���� ������ �����, ���������� ������� ���� ������� ���������)
    Ccntr=Ccntr+s;      %Ccntr ������ ����������� ���������� ���������� ��������� �����
  end
  Coords(Ccntr+1:end,:)=[]; %�������� �������������������� ���������
         %plot(Coords(:,2),Coords(:,1),'rx');
end
end
function [Ccntr,Coords]=FndRowByRowLinks(SrsCol)%StrtRow
%SrsCol - ����� ��������� (��������) ������� � StrBitMsk. � �� ������ ����� (���� ������ true) � �������� ���� �������.
%StrtRow - ����� ������ ��� � ���� (�������) �������� (SrsCol) ����������� �����(true-���) ������� � ������������ 
% ����� ������ ����� � ���� ��������. ���� ��������� ����� ���, �� StrtRow=0
%Coords - [Row Col] ������ ����� � ��������, ������� ����������� � �������� �������
% �������������� � �������� ������� ([Master Slave]������ ������). ����� ����� � ����� ������ Coords
%Ccntr - ���������� ��������� ����� ��� S���� ��������� (��������) ��� ������
global StrBitMsk;
global RequestRowMsk;
global RequestColMsk;
[r,~]=find(StrBitMsk(:,SrsCol));%������ ���� �� ������� ��������� � �������� ������� ��������
r=r.*RequestRowMsk(r);  %�������� �� ������ ��� ������������ ����� (� ��������)(� ��������� ���, ��-�� ������� ��� ����������� ����� � ���� ������)
r(r==0)=[]; %�������� �� ������ �������� ����� (� ��������) ��� �� ������� ��� ������������ �����
Csize=size(r,1);
RequestColMsk(SrsCol)=false;  %������� ������� ��� �������������  
if isempty(r) %��������� ������� �� �������, ���� �� ����� ����� (����� ���������) ������� �� ����
  Coords=[];
  Ccntr=0;
else
  Coords=zeros(Csize*3,2);  %������������� ��������� �������
  Coords(1:Csize,1)=r;      %������ Row-���������(������� �����) ���� �������� ����� � ������� �������
  Coords(1:Csize,2)=SrsCol;	%������ Col-���������� (������ ������������� �������) ���� �����
  Ccntr=Csize;  %������������� ������� ������� ����� (true-��� ��� ��� ������, ������� �������� ������� � ������� ������ ����������)
  %����� ��������� ������(����� � StrBitMsk) � ������� ��� ���� ������� �������� �������� �������
  for i=1:Csize	%������� �����, � ������� ������� ����� � ������� �������
    [s,Coords(Ccntr+1:Ccntr+s,:)]=FndColByColLinks(r(i)); %,SrsCol ����������� ��������� ����� ����������� � ����� ������ � ������ �� ��������� (� ������� �������) ����� (� ����� ����������� ��������� � ���� ������ �����, ���������� ������� ���� ������� ���������)
    Ccntr=Ccntr+s;  %���������� ��������
  end
  Coords(Ccntr+1:end,:)=[]; %�������� �������������������� ���������
         %plot(Coords(:,2),Coords(:,1),'rx');
end
end
function [ResAdr]=RecuFindM(N)
%������� ��� ������������ ������ ����������� (��������) ��������� Master-�������� (������)
%N - ����� ������� � MstrStrData
%ResAdr - ������ (� MstrStrLinkList) ���� ��������� � N-��� (��������) Master-��������.
global MstrStrLinkList; %��������� ������� � ������ (���������)������ Master-��������
ResAdr=find(MstrStrLinkList(:,1)==N);    %������� ����� (������ �����������)
for i=1:size(ResAdr,1)      %����� (�����������) �� ������� ���������� ������
    nAdr=RecuFindM(MstrStrLinkList(ResAdr(i),2));   %������ ������� MstrStrLinkList �������� ������� � �������� ������� ������ �� ������� �������
    ResAdr=cat(1,ResAdr,nAdr);  %����������� �������� ������� �������� ��������� ����� (������� �������) � ��������� ���������� � ������ ��������
end
end
function [ResAdr]=RecuFindS(N)
%������� ��� ������������ ������ ����������� (��������) ��������� Slave-�������� (������)
%N - ����� ������� � SlveStrData
%ResAdr - ������ (� SlveStrLinkList) ���� ��������� � N-��� (��������) Slave-��������.
global SlveStrLinkList; %��������� ������� � ������ (���������)������ Slave-��������
ResAdr=find(SlveStrLinkList(:,1)==N);    %������� ����� (������ �����������)
for i=1:size(ResAdr,1)      %����� (�����������) �� ������� ���������� ������
    nAdr=RecuFindS(SlveStrLinkList(ResAdr(i),2));   %������ ������� SlveStrLinkList �������� ������� � �������� ������� ������ �� ������� �������
    ResAdr=cat(1,ResAdr,nAdr);  %����������� �������� ������� �������� ��������� ����� (������� �������) � ��������� ���������� � ������ ��������
end
end