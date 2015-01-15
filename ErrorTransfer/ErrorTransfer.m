% �ó�����ƶ�λ������3����λ������landmark��վ��λ�����ڵĵ�ʱ���Ĵ���

% 2014��08��19��16��05 �޸�Ϊֻ��һ��BSC�в�������


% ***************************��ʼ��������(��ʼ)*****************************

% ��λ�����С landmark��վ��Ŀ 
Aerawidth=100; LandBSNum=3;

% ���ѡȡ��վʱ��վ��ľ���
LandBSspace=50;

% ��վ�㲥����Ϣ�趨
% BSbroadinfo 1:id 2:flag( Landmark/Blind )  3:xposition  4:yposition 5:headings 6:angle  
id=1; flag=2; xpos=3; ypos=4; headings=5; angle=6;  
LBSbroadinfo=zeros(LandBSNum,6);

% ��վ��Ϣ�趨 
LBSs=zeros( LandBSNum , 4 );

% Landmark BSCsվַѡȡ
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

% Landmark BSC headings���趨
LBSbroadinfo(:,headings)=(rand(LandBSNum,1)-0.5)*360;

% ��װҪ�㲥�Ļ�վϵͳ��Ϣ
LBSbroadinfo(:,id:ypos)=LBSs;

% ��������������
SampleNum=100;
samplexpos=100*rand(1,SampleNum);
sampleypos=100*rand(1,SampleNum);

% ����������з��Ǽ����Եõ�����Landmark BSC��ÿ��������ķ���
anglematrix=zeros(LandBSNum,SampleNum);
for index=1:1:SampleNum
    anglematrix(:,index)=generangle([ samplexpos(1,index),sampleypos(1,index)],LBSbroadinfo);
end

% ��Ҫ�ı���
loopNum=100;Error=zeros(1,SampleNum);

%���ȼ���û������������µ�ƽ�����

for i=1:1:SampleNum
    
    LBSbroadinfo(:,angle)=anglematrix(:,i);
    
    [estimX,estimY,BSbanned]=lslocation(LBSbroadinfo);
    Error(1,i)=sqrt( (samplexpos(1,i)-estimX).^2+(sampleypos(1,i)-estimY).^2);
end
NonnoisMeanError=mean(Error);


% �ڵ�һ��BSC�����в��� ratio*Areawidth ��ֵ�����������㶨λ��ƽ�����
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
            % �ĳ�ֻ�ڵ�һ����վλ���в������������������ģ���Ϊ�˻�վid��
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


