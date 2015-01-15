% 该程序估计定位区域内3个有位置误差的landmark基站定位区域内的点时误差的传递

% 2014年08月19日16：05 修改为只在一个BSC中掺入噪声


% ***************************初始环境设置(开始)*****************************

% 定位区域大小 landmark基站数目 
Aerawidth=100; LandBSNum=3;

% 随机选取基站时基站间的距离
LandBSspace=50;

% 基站广播的消息设定
% BSbroadinfo 1:id 2:flag( Landmark/Blind )  3:xposition  4:yposition 5:headings 6:angle  
id=1; flag=2; xpos=3; ypos=4; headings=5; angle=6;  
LBSbroadinfo=zeros(LandBSNum,6);

% 基站信息设定 
LBSs=zeros( LandBSNum , 4 );

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

% Landmark BSC headings角设定
LBSbroadinfo(:,headings)=(rand(LandBSNum,1)-0.5)*360;

% 组装要广播的基站系统信息
LBSbroadinfo(:,id:ypos)=LBSs;

% 产生样本点数据
SampleNum=100;
samplexpos=100*rand(1,SampleNum);
sampleypos=100*rand(1,SampleNum);

% 对样本点进行辐角计算以得到各个Landmark BSC对每个样本点的辐角
anglematrix=zeros(LandBSNum,SampleNum);
for index=1:1:SampleNum
    anglematrix(:,index)=generangle([ samplexpos(1,index),sampleypos(1,index)],LBSbroadinfo);
end

% 需要的变量
loopNum=100;Error=zeros(1,SampleNum);

%首先计算没有噪声的情况下的平均误差

for i=1:1:SampleNum
    
    LBSbroadinfo(:,angle)=anglematrix(:,i);
    
    [estimX,estimY,BSbanned]=lslocation(LBSbroadinfo);
    Error(1,i)=sqrt( (samplexpos(1,i)-estimX).^2+(sampleypos(1,i)-estimY).^2);
end
NonnoisMeanError=mean(Error);


% 在第一个BSC坐标中掺入 ratio*Areawidth 幅值的噪声并计算定位的平均误差
range=10;NoisMeanError=zeros(1,range);
for ratio=1:1:range
Error=zeros(1,SampleNum);
tmpError=zeros(1,loopNum);
tmpLBSbroadinfo=LBSbroadinfo;

    for i=1:1:SampleNum

        tmpError=zeros(1,loopNum);
        tmpLBSbroadinfo(:,angle)=anglematrix(:,i);

        for j=1:1:loopNum

            errorR=rand(1,1)*Aerawidth*ratio/100;theta=2*pi*rand(1,1);
            % 改成只在第一个基站位置中参入噪声（把行索引的：改为了基站id）
            tmpLBSbroadinfo(1,xpos)=LBSbroadinfo(1,xpos)+errorR*sin(theta);
            tmpLBSbroadinfo(1,ypos)=LBSbroadinfo(1,ypos)+errorR*cos(theta);

            [estimX,estimY,BSbanned]=lslocation(tmpLBSbroadinfo);
            tmpError(1,j)=sqrt( (samplexpos(1,i)-estimX).^2+(sampleypos(1,i)-estimY).^2);
        end
        Error(1,i)=mean(tmpError);  
    end
NoisMeanError(1,ratio)=mean(Error);
end

NoisMeanErrorDiff=NoisMeanError-NonnoisMeanError;
Noisgainpercent=100*(NoisMeanError-NonnoisMeanError)/NonnoisMeanError;


