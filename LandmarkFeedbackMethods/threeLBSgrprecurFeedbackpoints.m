% 2014 09 04 
% 此脚本为用三个Landmark BSC定位出Blind BSC后，提取出其中定位精度较高的Blind BSC标为 准Landmark
% BSC以进行迭代
% 2014 09 06 19：58 该脚本为仿真测试被Landmark测出的多个定位误差较低的BSC整体进行迭代时的定位情况

% 修正：修正进行迭代过程中，没有使用真实的坐标进行辐角估计。之前的脚本中有类似问题，应当特别注意
% 修正：迭代式应加入回测法（完成）Blind BSC已经可以相互回测修正 结果也收敛


%Line58 %^%%&%&%$%$^##%^ 设定LBSs位置____________________________________

% ***************************初始环境设置(开始)*****************************

% 定位区域大小    landmark基站数目  blind基站数目      基站间的距离
Aerawidth=100;   LandBSNum=3;      BlinBSNum=10;    LandBSspace=50;BlinBSspace=10;

% 基站广播的消息设定
% 1:id  2:flag( Landmark/Blind )  3:xposition  4:yposition 5:headings 6:angle  
id=1; flag=2; xpos=3; ypos=4; headings=5; angle=6;  
LBSbroadinfo=zeros(LandBSNum,6);BBSbroadinfo=zeros(BlinBSNum,6);

% 记录基站部分信息的变量 
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

% %^%%&%&%$%$^##%^ 设定LBSs位置
%LBSs(:,3:4)=[10,10;50,50;90,90];

% 组装一部分要广播的基站系统信息
LBSbroadinfo(:,id:ypos)=LBSs;BBSbroadinfo(:,id:ypos)=BBSs;



% Landmark BSC 与 Blind BSC headings角设定
LBSbroadinfo(:,headings)=(rand(LandBSNum,1)-0.5)*360;
BBSbroadinfo(:,headings)=(rand(BlinBSNum,1)-0.5)*360;

%*******************************需要的变量**********************************
% TrueBlinBSpos 记录Blind BSC真实值
TrueBlinBSinfo=BBSbroadinfo;

% Blind BSC 测量各Landmark BSC的真实辐角   Blind BSC以估计坐标为参照测量各个Landmark BSC的辐角
Realangle=zeros(1,LandBSNum); Falseangle=zeros(1,LandBSNum); 


%******************第一轮 Landmark BSC 进行定位估计**************************

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

InitialErrorDat=sqrt((BBSbroadinfo(:,xpos)-TrueBlinBSinfo(:,xpos)).^2+(BBSbroadinfo(:,ypos)-TrueBlinBSinfo(:,ypos)).^2);
InitialMeanError= mean(InitialErrorDat);

% ************************筛选高精度的BSC***************************************
% 注：这里主要是为了测试高精度点的定位性能所以先不用具体算法来估计定位精度，先使用真实的坐标值对
% 其性能进行评估

NearBSCArray=InitialErrorDat<100;
NearBSCNums=sum(NearBSCArray);

% 如果不存在就就此停止
if NearBSCNums~=0
    % 用于存放高精度的BSC的广播信息 
    NearLBSbroadinfo=zeros(NearBSCNums,6);
    TrueNearLBSinfo=zeros(NearBSCNums,6);
    for i=1:1:NearBSCNums
        
        index=find(NearBSCArray,1,'first');
        if ~isempty(index)
            NearLBSbroadinfo(i,:)=BBSbroadinfo(index,:);
            TrueNearLBSinfo(i,:)=TrueBlinBSinfo(index,:);
            NearBSCArray(index)=0;
        end
    
    end 
TestErrorDat=mean(sqrt((NearLBSbroadinfo(:,xpos)-TrueNearLBSinfo(:,xpos)).^2+(NearLBSbroadinfo(:,ypos)-TrueNearLBSinfo(:,ypos)).^2))

    % **************************开始进行迭代过程*******************************
    
    % 需要的变量 重新定义Realangle 和 Falseangle
    Realangle=zeros(1,LandBSNum+NearBSCNums-1);Falseangle=zeros(1,LandBSNum+NearBSCNums-1);
    
    
    loopnum=25;
    
    % 为观察各个BSC的位置坐标变化增加一个记录差的变量和记录上一次数值的变量
    lastBSCinfo=NearLBSbroadinfo;oneloopError=zeros(NearBSCNums,loopnum);TestError=zeros(NearBSCNums,loopnum);
    
    for loop=1:1:loopnum
    % 选取要更新的NearLBS BSC
    for i=1:1:NearBSCNums
    
	iNearLBSbroadinfo(1,:)=NearLBSbroadinfo(i,:);
	iTrueNearLBSbroadinfo(1,:)=TrueNearLBSinfo(i,:);
	
	tmpNearLBSbroadinfo=NearLBSbroadinfo;  tmpTrueNearLBSinfo=TrueNearLBSinfo;
	tmpNearLBSbroadinfo(i,:)=[];  tmpTrueNearLBSinfo(i,:)=[];
	
	mergeTrueBroadinfo=[LBSbroadinfo;tmpTrueNearLBSinfo];
	mergeBroadinfo=[ LBSbroadinfo;tmpNearLBSbroadinfo ];
	
	estimX=iNearLBSbroadinfo(1,xpos);
    estimY=iNearLBSbroadinfo(1,ypos);
	
        for landind=1:1:(LandBSNum+NearBSCNums-1)
            Realangle(1,landind)=generangle(mergeTrueBroadinfo(landind,xpos:ypos),iTrueNearLBSbroadinfo(1,:));
            Falseangle(1,landind)=generangle( mergeBroadinfo(landind,xpos:ypos) ,iNearLBSbroadinfo(1,:));
        end
	
	counter=1;
        while (max(mod(Realangle-Falseangle,360)>0)==1 && counter<LandBSNum*360)||(max(mod(Realangle(1,1:LandBSNum)-Falseangle(1,1:LandBSNum),360)>0)==1 && counter>=LandBSNum*360) 

            innerind=find(mod(Realangle-Falseangle,360)>0 ,1 ,'first');
           if (Realangle(1,innerind)-Falseangle(1,innerind))>180

               theta=Realangle(1,innerind)-Falseangle(1,innerind)-360;    
           elseif (Realangle(1,innerind)-Falseangle(1,innerind))<-180

               theta=Realangle(1,innerind)-Falseangle(1,innerind)+360;
           else

                theta=Realangle(1,innerind)-Falseangle(1,innerind);  
           end 
           
           vector=[estimX-mergeBroadinfo(innerind,xpos) ,estimY-mergeBroadinfo(innerind,ypos)];
           newvec=([ cos(0.1*theta*pi/180/abs(theta)) ,-sin(0.1*theta*pi/180/abs(theta)); sin(0.1*theta*pi/180/abs(theta)),cos(0.1*theta*pi/180/abs(theta))]*vector')';
           estimX=newvec(1,1)+mergeBroadinfo(innerind,xpos);estimY=newvec(1,2)+mergeBroadinfo(innerind,ypos);
           
           iNearLBSbroadinfo(1,xpos:ypos)=[estimX,estimY];
           for landind=1:1:(LandBSNum+NearBSCNums-1)
               Falseangle(1,landind)=generangle( mergeBroadinfo(landind,xpos:ypos) ,iNearLBSbroadinfo(1,:));
           end
           counter=counter+1;   
        end
    NearLBSbroadinfo(i,xpos:ypos)=iNearLBSbroadinfo(1,xpos:ypos);
    end
    
    TestErrorDat=mean(sqrt((NearLBSbroadinfo(:,xpos)-TrueNearLBSinfo(:,xpos)).^2+(NearLBSbroadinfo(:,ypos)-TrueNearLBSinfo(:,ypos)).^2))    
    % 为分析误差变化
    TestError(:,loop)=sqrt((NearLBSbroadinfo(:,xpos)-TrueNearLBSinfo(:,xpos)).^2+(NearLBSbroadinfo(:,ypos)-TrueNearLBSinfo(:,ypos)).^2);
    
    oneloopError(:,loop)=sqrt((NearLBSbroadinfo(:,xpos)-lastBSCinfo(:,xpos)).^2+(NearLBSbroadinfo(:,ypos)-lastBSCinfo(:,ypos)).^2);
    
    lastBSCinfo=NearLBSbroadinfo;
	   
    end
    
    %subplot(1,2,1); plot(TestError');subplot(1,2,2);plot(oneloopError');

end   
       
