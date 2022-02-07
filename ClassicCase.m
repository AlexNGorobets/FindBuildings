function [ StrArray ] = ClassicCase( Img )
%CLASSICCASE ��������� ����������� ����� �� ������ ����� + �������������� �����
%Coord - ���������� ������ �������� ([x1 x2 y1 y2] - � ������ ������ ��� ������� �������) 
% (��� �� �1 �2 r1 r2), �.�. ��� � ���������� ���� � ���� ��������� - ����� ������� ����
Img = imgaussfilt(Img, 2.8);
BW = edge(Img,'canny');
     %BW=bwmorph(BW,'endpoints');
RR=0.5;
TR=0.5;
[H,T,R] = hough(BW,'RhoResolution',RR,'ThetaResolution',TR);
%TrsHld=ceil(0.5*max(H(:)));
%TrsHld=ceil(0.4*max(H(:)));
%TrsHld=ceil(max(H(:))*graythresh(H(:)));
 TrsHld=35;
nlines=1000;
P  = houghpeaks(H,nlines,'threshold',TrsHld);
lines = houghlines(BW,T,R,P,'FillGap',6,'MinLength',12);
    figure, imshow(Img), hold on
max_len = 100;
    StrW=zeros(size(P,1),1);
    StrArray=zeros(size(P,1),12);    
for k = 1:length(lines)
    xy = [lines(k).point1; lines(k).point2];
    StrArray(k,1:4)=xy(:);  % [�1 �2 r1 r2] �������� �������
    StrArray(k,5:6)=flip(mean(xy)); %[CPRow CPCol] �������� �������
    StrArray(k,7)=lines(k).theta*pi./180;   %�������� ������ ���������
   StrW(k)=H(1+(lines(k).rho-min(R))./RR,1+(lines(k).theta-min(T))./RR);
   plot(xy(:,1),xy(:,2),'Color','green');   %'LineWidth',1.5,
%    % Plot beginnings and ends of lines
%    plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
%    plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
end
    %figure;
%     imshow(H,[],'XData',T,'YData',R,'InitialMagnification','fit');
% xlabel('\theta'), ylabel('\rho');
% axis on, axis square, hold on;
% plot(T(P(:,2)),R(P(:,1)),'s','color','white');
%Qpix=H(,)
%     imagesc(H); hold on;
%     plot(P(:,2),P(:,1),'s','color','white');

end

