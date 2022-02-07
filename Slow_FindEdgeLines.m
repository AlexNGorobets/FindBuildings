function [ Coord ] = Slow_FindEdgeLines( Data ,Upfr,DwnFr,LftFr,RghtFr)
%FINDEDGELINES �������� ������ ������� ��������
%Data - ����������� �����
%Coord - ���������� ������ �������� ([x1 y1 x2 y2] - � ������ ������ ��� ������� �������)

%��������� ��������/��������������
ColNeighborN=2; %����������� ��������� ��������-���������� �� ��������� ������� (����� ������ ������ ������), ��� ������ �������� ������ �������� � ������� ������
RowNeighborN=2; %����������� ��������� ��������-���������� �� ������� ("�����������" ������� �������� (�������) ������ ����� �������-�������)
LstAddPixAmmount=0; %����, ������� ��������� ��������� ������� � (�������) ������ ����������, ����� �� �������� �������� � ������ ����� ��������� LstAddPixAmmount �������� (=0 - ��������� ��������)
QPixTresholdInGroup=10; %����������� ���������� �������� � ������ ��� ����, ����� ��� ������ ����� ���� ����������� �� ������� ��������� ������
%��������� ��� �����������: (��� ��������� ����� �������������� �� ���������� ��������,
% � ����������� ������� ����� �� ����� ���������). ������� ������ � ������, ���� ���� �������� ������ 1
ColTrshold=2.2; %����� ���������� �� X (������� ������) �� ���������� ������� �� ������
RowTrshold=2.2; %����� ���������� �� Y (������ ����)
GrAngTrshold=25;%����� ������� ����� � �������� ����� ����� ��������� �������� ������� � ������� ����� � (�������) ������.
GrValTrshold=40;%����� ������� ������� ���������� �������� ������� � �������� �������� (�������) ������ (� ��������� �� ��.����. ������).
GrValTrshold=GrValTrshold/100;      %������� � ����
        SafeCoreFrame=4;    %SafeCoreFrame - ���������� ��� ����������� ���� ������ (�\���) ���� ��������
        Upfr=Upfr+SafeCoreFrame;
        DwnFr=DwnFr-SafeCoreFrame;
        LftFr=LftFr+SafeCoreFrame;
        RghtFr=RghtFr-SafeCoreFrame;
Data=double(Data);  %���������� ���� ������ � double
i=2:2:size(Data,1)-2;   %� ���������� ������� ������ ������ ������ ����� ������
Data(i,:)=(Data(i-1,:)+Data(i+1,:))/2;	% �������� ��������� �����������
    %imshow(uint16(Data(Upfr:DwnFr,LftFr:RghtFr)));
GaussKernel=fspecial('gaussian',[5 5],1.6); %������������ ���� ��������
BluredData=conv2(Data, GaussKernel,'same');	%�������� �� ������
        BluredData=BluredData(Upfr:DwnFr,LftFr:RghtFr);
    %imshow(uint16(Data(Upfr:DwnFr,LftFr:RghtFr)));
DiifKernel=GetRoundDiffKernel(5, [1 3]);
DY=conv2(BluredData, DiifKernel,'valid');    %������������� ����������� �� ������-����
DX=conv2(BluredData, DiifKernel','valid');   %�������������� ����������� �� � ���� �� �����
DiffVal=sqrt(DX.^2+DY.^2);  %�������� ��������� ��� ������� �������
DiffAngle=atan2(DY,DX);     %������ ��������� ��� ������� �������
DiffAngle=DiffAngle./pi.*180;%������� ������� � �������
%     ValIm=DiffVal(Upfr:DwnFr,LftFr:RghtFr);
%     AngIm=DiffAngle(Upfr:DwnFr,LftFr:RghtFr);
%     imagesc(AngIm);
%     % AngIm=SDiffAngle(Upfr:DwnFr,LftFr:RghtFr);
%     dBLPict=cat(4,Data(Upfr:DwnFr,LftFr:RghtFr),BluredData(Upfr:DwnFr,LftFr:RghtFr));
% %     subplot(2,1,1),montage(uint16(dBLPict),'Size',[1 2]);
% %     subplot(2,1,2),quiver(DX(Upfr:DwnFr,LftFr:RghtFr),DY(Upfr:DwnFr,LftFr:RghtFr));
%     quiver(DX(Upfr:DwnFr,LftFr:RghtFr),DY(Upfr:DwnFr,LftFr:RghtFr));
%����������� �������� �� ���������� (���� �� �����), ������ ������� � ����� �������� ���������

[DwnFr RghtFr]=size(BluredData);
%GCoord=zeros(400,2);    %���������� ������ ����� �������� � ������� ������ [row col] (������������� � ������� ������� �� ��������)
GCoord(1,:)=[RowNeighborN ColNeighborN];	%������ ���������� (������ ������ � ������)
        GCoord(1,:)=[27 91];
        %GCoord(1,:)=[20 104];
        %GCoord(1,:)=[16 105];
GMeanAng=DiffAngle(GCoord(1,1),GCoord(1,2)); %������� ������ ��������� � ������ (������ �����������)
GMeanVal=DiffVal(GCoord(1,1),GCoord(1,2));   %������� �������� ��������� ������ (������ �����������)
    PixGroups=cell(100,1);
    PixUsedMask=false(size(BluredData));%������������� ����� ��������, ������� ��� �������� � �.�. ������
PixOffrdMask=PixUsedMask;       	%������������� ����� ��������, �� ������� ��� ���� �������� �������� (����� ����� ��� ����� "������ �����������")
CrntGroupPixMask=PixUsedMask;	%����� �������� ������ � ������� ������ (��� ����������, ����� ������� ����� ������������ ���������� �������)
    %     CrntGroupOffrdMask=PixUsedMask;
    %     CrntGroupOffrdMask=double(CrntGroupOffrdMask);
                tic
g=1;    %�������� ������
while (g<=2 || ~isempty(GCoord)) %����-�������� �����

                tic;
        CrntGroupPixMask(:,:)=false;
        CrntGroupPixMask(GCoord(1,1),GCoord(1,2))=true;%������� ������� �����������
        PixOffrdMask(GCoord(1,1),GCoord(1,2))=true;
    Offred_IND=1;   %��������� ������������� ������ ������������� ������� (��� ����� � ������������ ����)
    while (size(Offred_IND,2)>LstAddPixAmmount)  %���� �������� ����������� ��� ������ ������
        Offrd_ColDist=[];   %��������� �������������
        Offrd_RowDist=[];
        PixCntr=0;
        for r=min(GCoord(:,1)):max(GCoord(:,1))	%������� �����
          TakenPixRange=find(CrntGroupPixMask(r,:));%����������� ��������� �������� � ������� ������, ������� ��� � ������
          if ~isempty(TakenPixRange)    %�� ������ � ������� ��� �� ������ ������� �� ������ �� ����� ���������� �� ������ �������
            if TakenPixRange(1)>ColNeighborN  %������ �� ������ �� ������� ������� ������ �� ������ ����
                StartAdr=TakenPixRange(1)-ColNeighborN;   
            else StartAdr=1;        %������ ��������� ������� ������������ ��������
            end
            if TakenPixRange(end)<(RghtFr-ColNeighborN)   %������ �� ������ �� ������� ������� ������ �� ������� ����
                FinshAdr=TakenPixRange(end)+ColNeighborN;
            else FinshAdr=RghtFr;   %����� ��������� ������� ������������ ��������
            end
            OfferedAdrss=StartAdr:FinshAdr; %������������ ������ ������� (��������, ������������ � ���������� � ������)
            for t=1:size(TakenPixRange,2) %�������� �� ��������� ������� ��������, ������� ��� � ������
             OfferedAdrss(OfferedAdrss==TakenPixRange(t))=[];            
            end
            Offrd_ColDist=cat(2,Offrd_ColDist,zeros(size(OfferedAdrss)));	%������������� ���������
            for p=1:size(OfferedAdrss,2)	%��������� ��������������� ���������� ����� (�������) ������ �� ������� �� ������������ �������� �� ���������� �� ���, ��� ��� ��������� � (�������) ������
             Offrd_ColDist(PixCntr+p)=min([TakenPixRange(find(TakenPixRange>OfferedAdrss(p),1))-OfferedAdrss(p) OfferedAdrss(p)-TakenPixRange(find(TakenPixRange<OfferedAdrss(p),1,'last'))]);
            end
            Offrd_r(PixCntr+1:PixCntr+p)=r;             %������ (����������) ������������ �������� (Y-���� ����������)
        	Offrd_c(PixCntr+1:PixCntr+p)=OfferedAdrss;  %������ ������������ �������� (X-������ ����������). ���������� ������ ��������� � ����� �������� ������������� ����� � ��� �� ������
            PixCntr=PixCntr+p;          %��������� �������� �������� ������������ ��������
          end
        end        
        for c=min(GCoord(:,2)):max(GCoord(:,2))	%������� ��������
          TakenPixRange=find(CrntGroupPixMask(:,c))';%����������� ��������� �������� � ������� ������, ������� ��� � ������
          if ~isempty(TakenPixRange)    %�� ������� � ������� ��� �� ������ ������� �� ������ - �� ����� ���������� �� ������ �������            
            if TakenPixRange(1)>RowNeighborN  %������ �� ������ �� ������� ������� ������ �� �������� ����
                StartAdr=TakenPixRange(1)-RowNeighborN;   
            else StartAdr=1;        %������ ��������� ������� ������������ ��������
            end
            if TakenPixRange(end)<(DwnFr-RowNeighborN)   %������ �� ������ �� ������� ������� ������ �� ������� ����
                FinshAdr=TakenPixRange(end)+RowNeighborN;
            else FinshAdr=DwnFr;   %����� ��������� ������� ������������ ��������
            end
            OfferedAdrss=StartAdr:FinshAdr; %������������ ������ ������� (��������, ������������ � ���������� � ������)
            for t=1:size(TakenPixRange,2) %�������� �� ��������� ������� ��������, ������� ��� � ������
             OfferedAdrss(OfferedAdrss==TakenPixRange(t))=[];            
            end
            Offrd_RowDist=cat(2,Offrd_RowDist,zeros(size(OfferedAdrss)));    %������������� ���������
            for p=1:size(OfferedAdrss,2)	%��������� ��������������� ���������� ����� (��������) ������� �� ������� �� ������������ �������� �� ���������� �� ���, ��� ��� ��������� � (�������) ������
             Offrd_RowDist(PixCntr+p)=min([TakenPixRange(find(TakenPixRange>OfferedAdrss(p),1))-OfferedAdrss(p) OfferedAdrss(p)-TakenPixRange(find(TakenPixRange<OfferedAdrss(p),1,'last'))]);
            end
            Offrd_ColDist(PixCntr+1:PixCntr+p)=Inf; %���������� ������������ �������� �� ������ (�� ��������� �������� - � ������ ������). ����� ������ ����������� ������� � ���, ��� ��� ���������� �� ���������� - ��� ����� ����� ������ �������� Offrd_ColDist � Offrd_RowDist ���������
            Offrd_r(PixCntr+1:PixCntr+p)=OfferedAdrss;  %������ (����������) ������������ �������� (Y-���� ����������)
        	Offrd_c(PixCntr+1:PixCntr+p)=c;          	%������ ������������ �������� (X-������ ����������). ���������� ������ ��������� � ����� �������� ������������� ����� � ��� �� ������
            if size(Offrd_r,2)>size(Offrd_ColDist,2)
              FndedDoubles=0;   %������� ��������� ������
              for t=1:p     %���� ������� ������������� ������ (������� ������������ ������)
              	CoordAdr1=find(Offrd_r(1:PixCntr)==Offrd_r(PixCntr+t-FndedDoubles));  %����� ���������� � ������� �����
                CoordAdr2=find(Offrd_c(CoordAdr1)==Offrd_c(PixCntr+t-FndedDoubles), 1);   %����� ���������� � �������� (�� ������� ������������� �����)
                if ~isempty(CoordAdr2)
                  Offrd_r(PixCntr+t-FndedDoubles)=[];
                  Offrd_c(PixCntr+t-FndedDoubles)=[];
                  Offrd_RowDist(CoordAdr1(CoordAdr2))=Offrd_RowDist(PixCntr+t-FndedDoubles);
                  Offrd_RowDist(PixCntr+t-FndedDoubles)=[];
                  Offrd_ColDist(PixCntr+t-FndedDoubles)=[];
                  FndedDoubles=FndedDoubles+1;
                end                
              end
              PixCntr=PixCntr+p-FndedDoubles;	%��������� �������� �������� ������������ ��������
            else
              PixCntr=PixCntr+p;       %��������� �������� �������� ������������ ��������    
           end
            
          end
        end
        Offred_IND=sub2ind(size(CrntGroupPixMask),Offrd_r,Offrd_c); %���������� �������� � ���������������� ���������
        	PixOffrdMask(Offred_IND)=true;
        %�������� ������������ �������� �� �������������� � ������
        PixAng=abs(DiffAngle(Offred_IND)-GMeanAng);     %������� �������� ������������ ������� � �������� �� ������
        PixAng(PixAng>180)=abs(PixAng(PixAng>180)-360); %���������� ������ �������� � ��������� �� 0 �� 180 ��������
        PixAng=PixAng/GrAngTrshold;                     %���������� �� treshold-�������� 
        %PixVal=abs(DiffVal(Offred_IND)-GMeanVal)/GMeanVal/GrValTrshold;  %������� ������� ���������� ������������ ������� � �������� �� ������ � ����� �� ��������� ������, ������������ �� treshold-��������
        PixVal=abs(DiffVal(Offred_IND)-GMeanVal)/GMeanVal/GrValTrshold;  %������� ������� ���������� ������������ ������� � �������� �� ������ � ����� �� ��������� ������, ������������ �� treshold-��������
        Offrd_RowDist(Offrd_RowDist==0)=Inf;	%����������� ��� ������� �������� ���������� ����������
        PixRDist=abs(RowTrshold-Offrd_RowDist)/RowTrshold; %������� ����� ������� (�� ��������� �����) �� ���������� ������� �� ������
        Offrd_ColDist(Offrd_ColDist==0)=Inf;	%����������� ��� ������� �������� ���������� ����������
        PixCDist=abs(ColTrshold-Offrd_ColDist)/ColTrshold; %������� ����� ������� (�� ��������� �����) �� ���������� ������� �� ������
        PixDist=min(cat(1,PixRDist,PixCDist),[],1);             %����� ���������� ���������� (����� �������� � ����������)
        %if ~((size(PixAng,1)==size(PixVal,1))&&(size(PixAng,1)==size(PixDist,1)) && ((size(PixAng,2)==size(PixVal,2))&&(size(PixAng,2)==size(PixDist,2))))
        if ~(size(PixAng,2)==size(PixDist,2))
            disp('��������')
        end
        ApprovmntMtric=sqrt(PixAng.^2+PixVal.^2+PixDist.^2);%������� ������ � ������ ���� ��������� ���������� � ���������� 
        Offred_IND=Offred_IND(ApprovmntMtric<=1);           % ������������ ������-��������-(����������)����������, � ������ ���������������� �� ��������������� treshold-�������� ������ �������
        CrntGroupPixMask(Offred_IND)=true;      %������� �������� � ������ �������� � ����� ������� ������
        GCoord=cat(1,GCoord,cat(1,Offrd_r(ApprovmntMtric<=1),Offrd_c(ApprovmntMtric<=1))'); %��������� ��������� �������� �������� � ������ �������� ������
        %����������� ������� ����������� ������
        GMeanAng=mean(DiffAngle(CrntGroupPixMask));
        GMeanVal=mean(DiffVal(CrntGroupPixMask));     

                    %CrntGroupOffrdMask(Offred_IND)=CrntGroupOffrdMask(Offred_IND)+1;
                 PlotBitArray=CrntGroupPixMask;
                 PlotBitArray=double(PlotBitArray);
                 PlotBitArray(Offred_IND)=PlotBitArray(Offred_IND)+1;
                 imagesc(PlotBitArray);
                 pause(0.1);
            
    end
            toc;
    if size(GCoord,1)>QPixTresholdInGroup
        PixGroups{g}=GCoord;
            PixUsedMask(GCoord(:,1),GCoord(:,2))=true;
        g=g+1;
    end
    %����� ���������� ����������� ��� ��������� ������
    GCoord=zeros(1,2);
    [GCoord(1,1), GCoord(1,2)]=find(PixOffrdMask==0,1);
    Offrd_r=[];
    Offrd_c=[];
    

Coord=0;
end
                toc;
end

