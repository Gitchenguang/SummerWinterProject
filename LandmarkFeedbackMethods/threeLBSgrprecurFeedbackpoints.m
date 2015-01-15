% 2014 09 04 
% �˽ű�Ϊ������Landmark BSC��λ��Blind BSC����ȡ�����ж�λ���Ƚϸߵ�Blind BSC��Ϊ ׼Landmark
% BSC�Խ��е���
% 2014 09 06 19��58 �ýű�Ϊ������Ա�Landmark����Ķ����λ���ϵ͵�BSC������е���ʱ�Ķ�λ���

% �������������е��������У�û��ʹ����ʵ��������з��ǹ��ơ�֮ǰ�Ľű������������⣬Ӧ���ر�ע��
% ����������ʽӦ����زⷨ����ɣ�Blind BSC�Ѿ������໥�ز����� ���Ҳ����


%Line58 %^%%&%&%$%$^##%^ �趨LBSsλ��____________________________________

% ***************************��ʼ��������(��ʼ)*****************************

% ��λ�����С    landmark��վ��Ŀ  blind��վ��Ŀ      ��վ��ľ���
Aerawidth=100;   LandBSNum=3;      BlinBSNum=10;    LandBSspace=50;BlinBSspace=10;

% ��վ�㲥����Ϣ�趨
% 1:id  2:flag( Landmark/Blind )  3:xposition  4:yposition 5:headings 6:angle  
id=1; flag=2; xpos=3; ypos=4; headings=5; angle=6;  
LBSbroadinfo=zeros(LandBSNum,6);BBSbroadinfo=zeros(BlinBSNum,6);

% ��¼��վ������Ϣ�ı��� 
LBSs=zeros( LandBSNum , 4 );BBSs=zeros(BlinBSNum ,4);

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


% Blind BSCs ��ʵվַѡȡ

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

% %^%%&%&%$%$^##%^ �趨LBSsλ��
%LBSs(:,3:4)=[10,10;50,50;90,90];

% ��װһ����Ҫ�㲥�Ļ�վϵͳ��Ϣ
LBSbroadinfo(:,id:ypos)=LBSs;BBSbroadinfo(:,id:ypos)=BBSs;



% Landmark BSC �� Blind BSC headings���趨
LBSbroadinfo(:,headings)=(rand(LandBSNum,1)-0.5)*360;
BBSbroadinfo(:,headings)=(rand(BlinBSNum,1)-0.5)*360;

%*******************************��Ҫ�ı���**********************************
% TrueBlinBSpos ��¼Blind BSC��ʵֵ
TrueBlinBSinfo=BBSbroadinfo;

% Blind BSC ������Landmark BSC����ʵ����   Blind BSC�Թ�������Ϊ���ղ�������Landmark BSC�ķ���
Realangle=zeros(1,LandBSNum); Falseangle=zeros(1,LandBSNum); 


%******************��һ�� Landmark BSC ���ж�λ����**************************

for i=1:1:BlinBSNum
    
    % Landmark BSC��������ֵ
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

% ************************ɸѡ�߾��ȵ�BSC***************************************
% ע��������Ҫ��Ϊ�˲��Ը߾��ȵ�Ķ�λ���������Ȳ��þ����㷨�����ƶ�λ���ȣ���ʹ����ʵ������ֵ��
% �����ܽ�������

NearBSCArray=InitialErrorDat<100;
NearBSCNums=sum(NearBSCArray);

% ��������ھ;ʹ�ֹͣ
if NearBSCNums~=0
    % ���ڴ�Ÿ߾��ȵ�BSC�Ĺ㲥��Ϣ 
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

    % **************************��ʼ���е�������*******************************
    
    % ��Ҫ�ı��� ���¶���Realangle �� Falseangle
    Realangle=zeros(1,LandBSNum+NearBSCNums-1);Falseangle=zeros(1,LandBSNum+NearBSCNums-1);
    
    
    loopnum=25;
    
    % Ϊ�۲����BSC��λ������仯����һ����¼��ı����ͼ�¼��һ����ֵ�ı���
    lastBSCinfo=NearLBSbroadinfo;oneloopError=zeros(NearBSCNums,loopnum);TestError=zeros(NearBSCNums,loopnum);
    
    for loop=1:1:loopnum
    % ѡȡҪ���µ�NearLBS BSC
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
    % Ϊ�������仯
    TestError(:,loop)=sqrt((NearLBSbroadinfo(:,xpos)-TrueNearLBSinfo(:,xpos)).^2+(NearLBSbroadinfo(:,ypos)-TrueNearLBSinfo(:,ypos)).^2);
    
    oneloopError(:,loop)=sqrt((NearLBSbroadinfo(:,xpos)-lastBSCinfo(:,xpos)).^2+(NearLBSbroadinfo(:,ypos)-lastBSCinfo(:,ypos)).^2);
    
    lastBSCinfo=NearLBSbroadinfo;
	   
    end
    
    %subplot(1,2,1); plot(TestError');subplot(1,2,2);plot(oneloopError');

end   
       
