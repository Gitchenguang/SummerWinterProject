function [estimX,estimY,BSbanned]=lslocation( BSbroadinfo )
% [estimX,estimY��BSbanned]=lslocation( BSbroadinfo )
%  BSbroadinfo=[ BSid , flag ,postionX,positionY,heading,radial]
% lslocation�����������flag�����Ա���BSC�Ƿ�ΪLandmark����Blind BSC

BSid=1;flag=2; posX=3;posY=4;heading=5;radial=6;

% ��������Ϣ���ػ�

angle=BSbroadinfo(:,heading)+BSbroadinfo(:,radial);
Xposition=BSbroadinfo(:,posX);Yposition=BSbroadinfo(:,posY);

rowNum=size(BSbroadinfo,1); len=1;
BSbanned=zeros(1,rowNum+1);pointsX=zeros(1,rowNum+1); pointsY=zeros(1,rowNum+1);pointsWeigh=zeros(1,rowNum+1);
pointsX(len)=0;pointsY(len)=0;pointsWeigh(len)=0;BSbanned(len)=0;

% 1�� �ж���ֱ�߹�ϵȷ�����㣨ƽ�е�����ˮƽ���򲻴�ֱ ��ˮƽ����ֱ���ǹ��� ƽ������ˮƽ����ֱ���ǲ����� �����ཻ���д�ֱ�����޴�ֱ�����֣� ��
for i=1:1:rowNum-1  
    for j=(i+1):1:rowNum
        % ����Ҫ����ֱ�ߵĽ��㣬��Ҫ��һЩ�ж�
        % 1���ж��Ƿ�ƽ�л����Ƿ�Ϊx=c��ֱ�� flagΪ�Ǹ�����ƽ�е�ȨֵΪ1����ֱΪ2
        flag=0;
        if abs(mod(angle(i,1),180))==abs(mod(angle(j,1),180))
            flag=flag+1;
        end
        if (abs(mod(angle(i,1)+90,360))==90) || (abs(mod(angle(i,1)+90,360))==270)||(abs(mod(angle(j,1)+90,360))==90) || (abs(mod(angle(j,1)+90,360))==270)
            flag=flag+2;
        end
        % 2������flag���ֱ�߷���
        switch flag
            case 1
                %���Ϊͬһ���ߣ����ý�ֹ��BS��������һ��
                ki=tan(angle(i,1)*pi/180+pi/2);kj=tan(angle(j,1)*pi/180+pi/2);
                x0i=Xposition(i,1);y0i=Yposition(i,1);
                x0j=Xposition(j,1);y0j=Yposition(j,1);
                if abs((y0i-x0i*ki)-(y0j-x0j*kj))>=10^(-2)
                    BSbanned(len)=BSbanned(len)+1;BSbanned(BSbanned(len)+1)=BSbroadinfo(i,BSid);
                end
            case 2
                    x0i=Xposition(i,1);y0i=Yposition(i,1);
                    x0j=Xposition(j,1);y0j=Yposition(j,1);
                    
                    if (abs(mod(angle(i,1)+90,360))==90) || (abs(mod(angle(i,1)+90,360))==270)
                        pointsX(len)=1+pointsX(len);pointsY(len)=1+pointsY(len);pointsWeigh(len)=1+pointsWeigh(len);
                        
                        pointsX(pointsX(len)+1)=Xposition(i,1);
                        pointsY(pointsY(len)+1)=y0j+tan(angle(j,1)*pi/180+pi/2)*(x0i-x0j);
                        
                        %ȷ������ֱ�ߵļнǵĴ�С��������ΪȨֵ�������꣬90�����
                        anglediff=abs(mod(abs((angle(i,1)-angle(j,1))),360));
                        if anglediff>180
                            anglediff=360-anglediff;
                        end
                        if anglediff>90
                            anglediff=180-anglediff;
                        end
                        pointsWeigh( pointsWeigh(len)+1 )=anglediff;
                   
                    else
                        
                        pointsX(len)=1+pointsX(len);pointsY(len)=1+pointsY(len);pointsWeigh(len)=1+pointsWeigh(len);
                        
                        pointsX(pointsX(len)+1)=Xposition(j,1);
                        pointsY(pointsY(len)+1)=y0i+tan(angle(i,1)*pi/180+pi/2)*(x0j-x0i);
                        
                        %ȷ������ֱ�ߵļнǵĴ�С��������ΪȨֵ�������꣬90�����
                        anglediff=abs(mod(abs((angle(i,1)-angle(j,1))),360));
                        if anglediff>180
                            anglediff=360-anglediff;
                        end
                        if anglediff>90
                            anglediff=180-anglediff;
                        end
                        pointsWeigh( pointsWeigh(len)+1 )=anglediff;                     
                    end
            case 3
                %������������ֹ���i��BS
                x0i=Xposition(i,1);x0j=Xposition(j,1);
                if abs(x0i-x0j)>=10^(-2)
                    BSbanned(len)=BSbanned(len)+1;BSbanned(BSbanned(len)+1)=BSbroadinfo(i,BSid);
                end
            case 0
                    %������Է�����
                    ki=tan(angle(i,1)*pi/180+pi/2);kj=tan(angle(j,1)*pi/180+pi/2);
                    x0i=Xposition(i,1);x0j=Xposition(j,1);
                    y0i=Yposition(i,1);y0j=Yposition(j,1);

                    pointsX(len)=1+pointsX(len);pointsY(len)=1+pointsY(len);pointsWeigh(len)=1+pointsWeigh(len);
                    
                    tmp=inv([1,-ki;1,-kj])*[-x0i*ki+y0i,-x0j*kj+y0j]';
                    pointsX(pointsX(len)+1)=tmp(2,1);pointsY(pointsY(len)+1)=tmp(1,1);
                    
                    %ȷ������ֱ�ߵļнǵĴ�С��������ΪȨֵ�������꣬90�����
                    anglediff=abs(mod(abs((angle(i,1)-angle(j,1))),360));
                        if anglediff>180
                            anglediff=360-anglediff;
                        end
                        if anglediff>90
                            anglediff=180-anglediff;
                        end
                    pointsWeigh( pointsWeigh(len)+1 )=anglediff;
                    
            otherwise
        end
    end
end

% ��LS������Ƶ�����ֵ
if pointsX(len)==0
    flase=1
else
    estimX=sum(pointsX((len+1):(pointsX(len)+1)).*pointsWeigh( (len+1):(pointsWeigh(len)+1) ) )/sum( pointsWeigh( (len+1):pointsWeigh( len)+1));
    estimY=sum(pointsY((len+1):(pointsY(len)+1)).*pointsWeigh( (len+1):(pointsWeigh(len)+1) ) )/sum( pointsWeigh( (len+1):pointsWeigh( len)+1));
end


