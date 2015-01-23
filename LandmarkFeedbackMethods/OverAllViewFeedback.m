% 2014 12 01 11：59
% 此脚本由threeLBSgrprecurFeedbackpoints（ 2014 09 06） 修改而来
% 脚本的思想：
%   鉴于基站较多时平均定位精度可达2m左右，所以考虑以landmark基站的位置作为系统的矫正坐标。使用Blind BSC
% 与Landmark BSC 组成的整个系统作为一个整体对独立的一个Landmark BSC进行定位。然后单独去掉一个Blind BSC
% 再进行定位，两次定位结果相互比较，对独立出来的Blind BSC的位置进行修正。不断迭代，直到系统的定位精度
% 符合一定的要求（例如对Landmark BSC的定位精度小于3m时停止迭代过程）。

% ***************************初始环境设置(开始)*****************************

% 定位区域大小    landmark基站数目  blind基站数目      基站间的距离
Aerawidth=100;   LandBSNum=3;      BlinBSNum=10;    LandBSspace=50;BlinBSspace=10;
% 基站广播的消息设定
% 1:id  2:flag( Landmark/Blind )  3:xposition  4:yposition 5:headings 6:angle  
id=1; flag=2; xpos=3; ypos=4; headings=5; angle=6;  
LBSbroadinfo=zeros(LandBSNum,6);BBSbroadinfo=zeros(BlinBSNum,6);
% 记录基站部分信息的变量 
LBSs=zeros( LandBSNum , 4 );BBSs=zeros(BlinBSNum ,4);
%-------------LBSs-----------
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
% Landmark BSC标志位置1
LBSs(:,flag)=ones(LandBSNum,1);
%-------------BBSs-----------
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
% Blind BSCs 标志位置1
BBSs(:,flag)=zeros(BlinBSNum,1);
%-------------All BSCs headings 设定---
LBSbroadinfo(:,headings)=(rand(LandBSNum,1)-0.5)*360;
BBSbroadinfo(:,headings)=(rand(BlinBSNum,1)-0.5)*360;
%------------组装BSC信息---
LBSbroadinfo(:,id:ypos)=LBSs;BBSbroadinfo(:,id:ypos)=BBSs;

% *****************************需要的变量******************************
% TrueBlinBSpos 记录Blind BSC真实值
TrueBlinBSinfo=BBSbroadinfo;
% Blind BSC 测量各Landmark BSC的真实辐角  
Realangle=zeros(1,LandBSNum);
% Blind BSC以估计坐标为参照测量各个Landmark BSC的辐角
Falseangle=zeros(1,LandBSNum); 

% *************************定位系统Blind BSCs位置初步估计*******************
for i=1:1:BlinBSNum
    
    % Landmark BSC给出估计值
    iBSxy=BBSbroadinfo(i,xpos:ypos);
    
    LBSbroadinfo(:,angle)=generangle( iBSxy , LBSbroadinfo);
    [estimX,estimY,BSbanned]=lslocation( LBSbroadinfo );
    BBSbroadinfo(i,xpos:ypos)=[estimX,estimY];
    
    for landind=1:1:LandBSNum
        Realangle(1,landind)= generangle( LBSbroadinfo(landind,xpos:ypos), TrueBlinBSinfo(i,:) );
        Falseangle(1,landind)= generangle( LBSbroadinfo(landind,xpos:ypos), BBSbroadinfo(i,:) );
    end
    
    % Adjust the position
     while max(mod(Realangle-Falseangle,360)>0)==1 
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
     end

end

% ********************从整体的角度进行测量修正***********************
% 定义局部变量loopNum 
loopNum=1000; 
% 定义BlinRealangle BlinFalseangle
BlinRealangle=zeros(1,LandBSNum);BlinFlaseangle=zeros(1,LandBSNum);
% 可能用到的变量
tmpBSCGrp=BBSbroadinfo;tmpcutBSCGrp=BBSbroadinfo;
tmpLBSpos=LBSbroadinfo;tmpcutLBSpos=LBSbroadinfo;

for loop=1:1:3000
    for i=1:1:BlinBSNum
        
       tmpBSCGrp=BBSbroadinfo;tmpcutBSCGrp=BBSbroadinfo;
       tmpLBSpos=LBSbroadinfo;tmpcutLBSpos=LBSbroadinfo;
       
       for j=1:1:LandBSNum
           % 获得此刻Blind BSCs 对第j个Landmark BSC的位置估计
           tmpBSCGrp(:,angle)=generangle( LBSbroadinfo(j,xpos:ypos),TrueBlinBSinfo);
           tmpcutBSCGrp(:,angle)=tmpBSCGrp(:,angle);
           
           tmpcutLBSpos(j,:)=[];
           tmpcutLBSpos(:,angle)=generangle( LBSbroadinfo(j,xpos:ypos),tmpcutLBSpos);
           
           % 去除第i个Blind BSC后，重新定位Landmark
           tmpcutBSCGrp(i,:)=[];      
           [estimX1,estimY1]=lslocation([tmpBSCGrp;tmpcutLBSpos]);
           [estimX2,estimY2]=lslocation([tmpcutBSCGrp;tmpcutLBSpos]);
           
           % 判断位置并进行修正
           if sum(([estimX2 estimY2]-LBSbroadinfo(j,xpos:ypos)).^2)<sum(([estimX1,estimY1]-LBSbroadinfo(j,xpos:ypos)).^2)
                vector=[estimX2,estimY2]-[estimX1,estimY1];
                BBSbroadinfo(i,xpos:ypos)=BBSbroadinfo(i,xpos:ypos)+rand(1,1)*0.1*vector;  % 这里0.1是一个比例
           end
           
           tmpBSCGrp=BBSbroadinfo;tmpcutBSCGrp=BBSbroadinfo;tmpcutLBSpos=LBSbroadinfo;
           
        end
    end
    Error(1,loop)=sum(sqrt(sum(((BBSbroadinfo(:,xpos:ypos)-TrueBlinBSinfo(:,xpos:ypos)).^2)')))/BlinBSNum;
    plot(Error,'DisplayName','Error','YDataSource','Error');
end
  













