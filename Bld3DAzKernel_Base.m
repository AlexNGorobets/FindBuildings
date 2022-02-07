function [ Kernel ] = Bld3DAzKernel_Base( Size, Azimuth )
%BLDAZKERNEL ���������� ���� � ���������� �������������� �� ������ ����������������� Azimuth
Kernel=zeros(Size,Size,max(size(Azimuth)));
Sgma=1;
CP=round(Size./2);
for Col=1:Size;
  for Row=1:Size      %���������� ���� � ������� ��������
    R=sqrt((Col-CP).^2+(Row-CP).^2);   
    Phi=atan2(Row-CP,Col-CP)-Azimuth;   %������� ������ ����� � (����� ����� ��������)
    x=R.*cos(Phi(:));   %���������� ��������� ���������� �� � ������ �������� ���������� ��
    Kernel(Row,Col,:)=1/(Sgma*sqrt(2*pi))*exp(-(x(:)).^2./(2*Sgma^2)); %��������� 
  end
end
%     for t=1:max(size(Azimuth))
%     imagesc(Kernel(:,:,t))
%     pause(0.1);
%     end
end


