function [  ] = tmpShow(  )
%TMPSHOW Summary of this function goes here
%���������� ���������� ���������

KrnlSgma=1.2;     %������� �������������
FldrName1=['KrnlSgma',num2str(fix(KrnlSgma*100))];    %Char-������ � ��������� ����� � ������� ��������
KrnlSgma=1.5;
FldrName2=['KrnlSgma',num2str(fix(KrnlSgma*100))];
KrnlSgma=1.8;
FldrName3=['KrnlSgma',num2str(fix(KrnlSgma*100))];

files=dir(FldrName1);  %��������� ������ ������ �� �����

for i=3:size(files,1)
    CrntFileName=files(i).name;
    %FID=fopen([FldrName1,'/',CrntFileName],'r');
    %Img1=fread(FID,'single');
    %fclose(FID);
    load([FldrName1,'/',CrntFileName]);
    Img1=StrEdgesLayer;
    %FID=fopen([FldrName2,'/',CrntFileName],'r');
    %Img2=fread(FID,'single');
    %fclose(FID);
    load([FldrName2,'/',CrntFileName]);
    Img2=StrEdgesLayer;
    %FID=fopen([FldrName3,'/',CrntFileName],'r');
    %Img3=fread(FID,'single');
    %fclose(FID);
    load([FldrName3,'/',CrntFileName]);
    Img3=StrEdgesLayer;
    
figure
% subplot(1,3,1)
% imagesc(Img1);
% subplot(1,3,2)
% imagesc(Img2);
% subplot(1,3,3)
% imagesc(Img3);
 subplot(1,2,1)
 imagesc(Img1-Img2);
 subplot(1,2,2)
 imagesc(Img2-Img3);


    
    
end



end

