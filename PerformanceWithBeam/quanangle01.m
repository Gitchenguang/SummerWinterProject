function [angle]=quanangle01(xy,Beam)

% [angle]=quanangle(xy)
% 2014-09-09 19��26 ʹ��������������Beam���� 
% �ò�������Ϊ360��һ��Լ���Ҹ�Լ������һ��4�ı������

N=180/Beam;

if mod(N/2,2)==0

    sector=1:1:(N+1);sector=sector-N/2-1;

    sector(1,1:(N/2))=sector(1,1:(N/2))+1/2;  % С���㲿��
    sector(1,(N/2+2):(N+1))=sector(1,(N/2+2):(N+1))-1/2;%�����㲿��

    radangle=sector*Beam/180*pi;

    quanref=tan(radangle);

    kXY=xy(1,2)/xy(1,1);

    refer=kXY-quanref;

    if refer(1,N+1)>0 || refer(1,1)<0
        if refer(1,N+1)>0
            angle=90;
        end
        if refer(1,1)<0
            angle=-90;
        end   
    else
        ind1=find(refer<=0 ,1,'first');
        ind2=find(refer>0,1,'last');
        angle=(sector(ind1)+sector(ind2))*Beam/2;
    end
else
    sector=1:1:(N-1);sector=sector-N/2;
    radangle=sector*Beam/180*pi;
    
    quanref=tan(radangle);
    
    kXY=xy(1,2)/xy(1,1);
    
    refer=kXY-quanref;
    
    if refer(1,N-1)>0 || refer(1,1)<0
        if refer(1,N-1)>0
            angle=90-Beam/2;
        end
        if refer(1,1)<0
            angle=-90+Beam/2;
        end   
    else
        ind1=find(refer<=0 ,1,'first');
        ind2=find(refer>0,1,'last');
        angle=(sector(ind1)+sector(ind2))*Beam/2;
    end
end


