% 2014 08 27 19：08 
% 100m*100m区域3Landmark BSC + 10 Blind BSC
% 这个脚本为在测量时加入Blind回测后的迭代程序

% ***************************初始环境设置(开始)*****************************

% 定位区域大小 landmark基站数目  blind基站数目
Aerawidth=100; LandBSNum=3;BlinBSNum=10;

% 随机选取基站时基站间的距离
LandBSspace=50;BlinBSspace=10;

% 基站广播的消息设定
% BSbroadinfo 1:id 2:flag( Landmark/Blind )  3:xposition  4:yposition 5:headings 6:angle  
id=1; flag=2; xpos=3; ypos=4; headings=5; angle=6;  
LBSbroadinfo=zeros(LandBSNum,6);BBSbroadinfo=zeros(BlinBSNum,6);

% 基站信息设定 
LBSs=zeros( LandBSNum , 4 );BBSs=zeros(BlinBSNum ,4);

% Landmark BSCs站址选取
for LBSid=1:1:LandBSNum
    LBSs(LBSid,id)=LBSid;
    LBSs( LBSid,xpos:ypos )=100*rand( 1,2);
    for checkid=1:1:LBSid-1     
        while min(min(sqrt((repmat(LBSs( LBSid,xpos),checkid,1)-LBSs( 1:checkid , xpos)).^2+(repmat(LBSs( LBSid,ypos),checkid,1)-LBSs( 1:checkid , ypos)).^2)))<LandBSspace
            LBSs( LBSid,xpos:ypos )=100*rand( 1,2);
       end
    end    
end

LBSs(:,flag)=ones(LandBSNum,1);


% Blind BSCs 真实站址选取

for BBSid=1:1:BlinBSNum
    BBSs( BBSid ,id )=BBSid;
    BBSs( BBSid ,xpos:ypos )=100*rand(1,2);
    for checkid=1:1:BBSid-1
        while min(min( sqrt((repmat(BBSs( BBSid,xpos),checkid,1)-BBSs(1:checkid,xpos)).^2+(repmat(BBSs( BBSid,ypos),checkid,1)-BBSs(1:checkid,ypos)).^2 ) ))<BlinBSspace || ...
        min(min(sqrt((repmat(BBSs( BBSid,xpos),LandBSNum,1)-LBSs(: , xpos )).^2+(repmat(BBSs( BBSid,ypos),LandBSNum,1)-LBSs(: , ypos )).^2)))<BlinBSspace
        
            BBSs( BBSid,xpos:ypos)=100*rand(1,2);
        end
    end
end

BBSs(:,flag)=zeros(BlinBSNum,1);

% 组装要广播的基站系统信息
LBSbroadinfo(:,id:ypos)=LBSs;BBSbroadinfo(:,id:ypos)=BBSs;

% Landmark BSC 与 Blind BSC headings角设定
LBSbroadinfo(:,headings)=(rand(LandBSNum,1)-0.5)*360;
BBSbroadinfo(:,headings)=(rand(BlinBSNum,1)-0.5)*360;

% ******************************需要的变量*********************************

% TrueBlinBSpos 记录Blind BSC真实值
TrueBlinBSpos=BBSbroadinfo;
% Blind BSC 测量各Landmark BSC的真实辐角   Blind BSC以估计坐标为参照测量各个Landmark BSC的辐角
Realangle=zeros(1,LandBSNum); Falseangle=zeros(1,LandBSNum); 

% ******************************计算过程*********************************
%  第一轮用Landmark BSC进行定位
for i=1:1:BlinBSNum
    
    % Landmark BSC给出估计值
    iBSxy=BBSbroadinfo(i,xpos:ypos);
    
    LBSbroadinfo(:,angle)=generangle( iBSxy , LBSbroadinfo);
    [estimX,estimY,BSbanned]=lslocation( LBSbroadinfo );
    BBSbroadinfo(i,xpos:ypos)=[estimX,estimY];
    
    for landind=1:1:LandBSNum
        Realangle(1,landind)= generangle( LBSbroadinfo(landind,xpos:ypos), TrueBlinBSpos(i,:) );
        Falseangle(1,landind)= generangle( LBSbroadinfo(landind,xpos:ypos), BBSbroadinfo(i,:) );
    end
    
    % Adjust the position
     counter=1;
     while max(mod(Realangle-Falseangle,360)>0)==1 && counter<LandBSNum*360
           innerind=find(mod(Realangle-Falseangle,360)>0 ,1 ,'first');
           if (Realangle(1,innerind)-Falseangle(1,innerind))>180

               theta=Realangle(1,innerind)-Falseangle(1,innerind)-360;    
           elseif (Realangle(1,innerind)-Falseangle(1,innerind))<-180

               theta=Realangle(1,innerind)-Falseangle(1,innerind)+360;
           else

                theta=Realangle(1,innerind)-Falseangle(1,innerind);  
           end 
           vector=[estimX-LBSbroadinfo(innerind,xpos) ,estimY-LBSbroadinfo(innerind,ypos)];
           newvec=([ cos(0.1*theta*pi/180/abs(theta)) ,-sin(0.1*theta*pi/180/abs(theta)); sin(0.1*theta*pi/180/abs(theta)),cos(0.1*theta*pi/180/abs(theta))]*vector')';
           estimX=newvec(1,1)+LBSbroadinfo(innerind,xpos);estimY=newvec(1,2)+LBSbroadinfo(innerind,ypos);
           
           BBSbroadinfo(i,xpos:ypos)=[estimX,estimY];
           
           for landind=1:1:LandBSNum
                Falseangle(1,landind)= generangle( LBSbroadinfo(landind,xpos:ypos), BBSbroadinfo(i,:) );
           end
           counter=counter+1;     
     end

end
PureLBSerrorDat=sqrt((BBSbroadinfo(:,xpos)-TrueBlinBSpos(:,xpos)).^2+(BBSbroadinfo(:,ypos)-TrueBlinBSpos(:,ypos)).^2);
LBSmeanerror= mean(sqrt((BBSbroadinfo(:,xpos)-TrueBlinBSpos(:,xpos)).^2+(BBSbroadinfo(:,ypos)-TrueBlinBSpos(:,ypos)).^2));

% 进行3带1迭代定位

% 需要的变量
LoopNum=1;tmpposinfo=BBSbroadinfo(LandBSNum,xpos:ypos);
recRealangle=zeros(1,LandBSNum+1); recFlaseangle=zeros(1,LandBSNum+1);
Recordpos=zeros(BlinBSNum,LandBSNum*2);
Realangle=zeros(1,LandBSNum+1); Falseangle=zeros(1,LandBSNum+1); 

for loopind=1:1:LoopNum
    for i=1:1:BlinBSNum
        TrueLocagrp=[ LBSbroadinfo;TrueBlinBSpos(i,:)];
        tmpBBSbroadinfo=BBSbroadinfo;
        Locagrp=[ LBSbroadinfo ; tmpBBSbroadinfo(i,:)];
        tmpBBSbroadinfo(i,:)=[];
        
        % 估计其余BS的位置
        for index=1:1:(BlinBSNum-1)
            
            estimX=tmpBBSbroadinfo(index,xpos);estimY=tmpBBSbroadinfo(index,ypos);
            for landind=1:1:(LandBSNum+1)
                Realangle(1,landind)= generangle( TrueLocagrp(landind,xpos:ypos), TrueBlinBSpos(index,:) );
                Falseangle(1,landind)= generangle( Locagrp(landind,xpos:ypos), tmpBBSbroadinfo(index,:) );
            end
            
            counter=1;   
            while max(mod(Realangle-Falseangle,360)>0)==1 && (counter<(LandBSNum+1)*360)%||(max(mod(Realangle(1,1:LandBSNum)-Falseangle(1,1:LandBSNum),360))>0))
                innerind=find(mod(Realangle-Falseangle,360)>0 ,1 ,'first');
                if (Realangle(1,innerind)-Falseangle(1,innerind))>180

                    theta=Realangle(1,innerind)-Falseangle(1,innerind)-360;    
                elseif (Realangle(1,innerind)-Falseangle(1,innerind))<-180

                    theta=Realangle(1,innerind)-Falseangle(1,innerind)+360;
                else

                    theta=Realangle(1,innerind)-Falseangle(1,innerind);  
                end 
                vector=[estimX-Locagrp(innerind,xpos) ,estimY-Locagrp(innerind,ypos)];
                newvec=([ cos(0.1*theta*pi/180/abs(theta)) ,-sin(0.1*theta*pi/180/abs(theta)); sin(0.1*theta*pi/180/abs(theta)),cos(0.1*theta*pi/180/abs(theta))]*vector')';
                estimX=newvec(1,1)+Locagrp(innerind,xpos);
                estimY=newvec(1,2)+Locagrp(innerind,ypos);
           
                tmpBBSbroadinfo(index,xpos:ypos)=[estimX,estimY];
           
                for landind=1:1:(LandBSNum+1)
                    Falseangle(1,landind)= generangle( Locagrp(landind,xpos:ypos), tmpBBSbroadinfo(index,:) );
                end
                counter=counter+1;      
            end
        end
        % 记录
        Recordpos(:,i)=[tmpBBSbroadinfo(1:(i-1),xpos);zeros(1,1);tmpBBSbroadinfo( i:(BlinBSNum-1),xpos)];
        Recordpos(:,i+BlinBSNum)=[tmpBBSbroadinfo(1:(i-1),ypos);zeros(1,1);tmpBBSbroadinfo( i:(BlinBSNum-1),ypos)];%Locagrp(LandBSNum+1,ypos)
    end
    % 更新各个Blind BSC 坐标
    
end


