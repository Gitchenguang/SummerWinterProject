function [estimX,estimY,BSbanned]=lslocation( BSbroadinfo )
% [estimX,estimY，BSbanned]=lslocation( BSbroadinfo )
%  BSbroadinfo=[ BSid , flag ,postionX,positionY,heading,radial]
% lslocation函数新添加了flag项用以标明BSC是否为Landmark或者Blind BSC

BSid=1;flag=2; posX=3;posY=4;heading=5;radial=6;

% 将数据信息本地化

angle=BSbroadinfo(:,heading)+BSbroadinfo(:,radial);
Xposition=BSbroadinfo(:,posX);Yposition=BSbroadinfo(:,posY);

rowNum=size(BSbroadinfo,1); len=1;
BSbanned=zeros(1,rowNum+1);pointsX=zeros(1,rowNum+1); pointsY=zeros(1,rowNum+1);pointsWeigh=zeros(1,rowNum+1);
pointsX(len)=0;pointsY(len)=0;pointsWeigh(len)=0;BSbanned(len)=0;

% 1、 判断是直线关系确定交点（平行但是与水平方向不垂直 与水平方向垂直但是共线 平行且与水平方向垂直但是不共线 仅仅相交（有垂直线与无垂直线两种） ）
for i=1:1:rowNum-1  
    for j=(i+1):1:rowNum
        % 这里要计算直线的交点，需要做一些判断
        % 1、判断是否平行或者是否为x=c的直线 flag为非负数，平行的权值为1，垂直为2
        flag=0;
        if abs(mod(angle(i,1),180))==abs(mod(angle(j,1),180))
            flag=flag+1;
        end
        if (abs(mod(angle(i,1)+90,360))==90) || (abs(mod(angle(i,1)+90,360))==270)||(abs(mod(angle(j,1)+90,360))==90) || (abs(mod(angle(j,1)+90,360))==270)
            flag=flag+2;
        end
        % 2、依据flag求解直线方程
        switch flag
            case 1
                %如果为同一条线，不用禁止该BS，进行下一步
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
                        
                        %确定两条直线的夹角的大小，以其作为权值计算坐标，90度最大
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
                        
                        %确定两条直线的夹角的大小，以其作为权值计算坐标，90度最大
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
                %如果不共线则禁止这第i个BS
                x0i=Xposition(i,1);x0j=Xposition(j,1);
                if abs(x0i-x0j)>=10^(-2)
                    BSbanned(len)=BSbanned(len)+1;BSbanned(BSbanned(len)+1)=BSbroadinfo(i,BSid);
                end
            case 0
                    %求解线性方程组
                    ki=tan(angle(i,1)*pi/180+pi/2);kj=tan(angle(j,1)*pi/180+pi/2);
                    x0i=Xposition(i,1);x0j=Xposition(j,1);
                    y0i=Yposition(i,1);y0j=Yposition(j,1);

                    pointsX(len)=1+pointsX(len);pointsY(len)=1+pointsY(len);pointsWeigh(len)=1+pointsWeigh(len);
                    
                    tmp=inv([1,-ki;1,-kj])*[-x0i*ki+y0i,-x0j*kj+y0j]';
                    pointsX(pointsX(len)+1)=tmp(2,1);pointsY(pointsY(len)+1)=tmp(1,1);
                    
                    %确定两条直线的夹角的大小，以其作为权值计算坐标，90度最大
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

% 用LS求出估计的坐标值
if pointsX(len)==0
    flase=1
else
    estimX=sum(pointsX((len+1):(pointsX(len)+1)).*pointsWeigh( (len+1):(pointsWeigh(len)+1) ) )/sum( pointsWeigh( (len+1):pointsWeigh( len)+1));
    estimY=sum(pointsY((len+1):(pointsY(len)+1)).*pointsWeigh( (len+1):(pointsWeigh(len)+1) ) )/sum( pointsWeigh( (len+1):pointsWeigh( len)+1));
end


