% 2014 12 01 11��59
% �˽ű���threeLBSgrprecurFeedbackpoints�� 2014 09 06�� �޸Ķ���
% �ű���˼�룺
%   ���ڻ�վ�϶�ʱƽ����λ���ȿɴ�2m���ң����Կ�����landmark��վ��λ����Ϊϵͳ�Ľ������ꡣʹ��Blind BSC
% ��Landmark BSC ��ɵ�����ϵͳ��Ϊһ������Զ�����һ��Landmark BSC���ж�λ��Ȼ�󵥶�ȥ��һ��Blind BSC
% �ٽ��ж�λ�����ζ�λ����໥�Ƚϣ��Զ���������Blind BSC��λ�ý������������ϵ�����ֱ��ϵͳ�Ķ�λ����
% ����һ����Ҫ�������Landmark BSC�Ķ�λ����С��3mʱֹͣ�������̣���

% ***************************��ʼ��������(��ʼ)*****************************

% ��λ�����С    landmark��վ��Ŀ  blind��վ��Ŀ      ��վ��ľ���
Aerawidth=100;   LandBSNum=3;      BlinBSNum=10;    LandBSspace=50;BlinBSspace=10;
% ��վ�㲥����Ϣ�趨
% 1:id  2:flag( Landmark/Blind )  3:xposition  4:yposition 5:headings 6:angle  
id=1; flag=2; xpos=3; ypos=4; headings=5; angle=6;  
LBSbroadinfo=zeros(LandBSNum,6);BBSbroadinfo=zeros(BlinBSNum,6);
% ��¼��վ������Ϣ�ı��� 
LBSs=zeros( LandBSNum , 4 );BBSs=zeros(BlinBSNum ,4);
%-------------LBSs-----------
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
% Landmark BSC��־λ��1
LBSs(:,flag)=ones(LandBSNum,1);
%-------------BBSs-----------
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
% Blind BSCs ��־λ��1
BBSs(:,flag)=zeros(BlinBSNum,1);
%-------------All BSCs headings �趨---
LBSbroadinfo(:,headings)=(rand(LandBSNum,1)-0.5)*360;
BBSbroadinfo(:,headings)=(rand(BlinBSNum,1)-0.5)*360;
%------------��װBSC��Ϣ---
LBSbroadinfo(:,id:ypos)=LBSs;BBSbroadinfo(:,id:ypos)=BBSs;

% *****************************��Ҫ�ı���******************************
% TrueBlinBSpos ��¼Blind BSC��ʵֵ
TrueBlinBSinfo=BBSbroadinfo;
% Blind BSC ������Landmark BSC����ʵ����  
Realangle=zeros(1,LandBSNum);
% Blind BSC�Թ�������Ϊ���ղ�������Landmark BSC�ķ���
Falseangle=zeros(1,LandBSNum); 

% *************************��λϵͳBlind BSCsλ�ó�������*******************
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

% ********************������ĽǶȽ��в�������***********************
% ����ֲ�����loopNum 
loopNum=1000; 
% ����BlinRealangle BlinFalseangle
BlinRealangle=zeros(1,LandBSNum);BlinFlaseangle=zeros(1,LandBSNum);
% �����õ��ı���
tmpBSCGrp=BBSbroadinfo;tmpcutBSCGrp=BBSbroadinfo;
tmpLBSpos=LBSbroadinfo;tmpcutLBSpos=LBSbroadinfo;

for loop=1:1:3000
    for i=1:1:BlinBSNum
        
       tmpBSCGrp=BBSbroadinfo;tmpcutBSCGrp=BBSbroadinfo;
       tmpLBSpos=LBSbroadinfo;tmpcutLBSpos=LBSbroadinfo;
       
       for j=1:1:LandBSNum
           % ��ô˿�Blind BSCs �Ե�j��Landmark BSC��λ�ù���
           tmpBSCGrp(:,angle)=generangle( LBSbroadinfo(j,xpos:ypos),TrueBlinBSinfo);
           tmpcutBSCGrp(:,angle)=tmpBSCGrp(:,angle);
           
           tmpcutLBSpos(j,:)=[];
           tmpcutLBSpos(:,angle)=generangle( LBSbroadinfo(j,xpos:ypos),tmpcutLBSpos);
           
           % ȥ����i��Blind BSC�����¶�λLandmark
           tmpcutBSCGrp(i,:)=[];      
           [estimX1,estimY1]=lslocation([tmpBSCGrp;tmpcutLBSpos]);
           [estimX2,estimY2]=lslocation([tmpcutBSCGrp;tmpcutLBSpos]);
           
           % �ж�λ�ò���������
           if sum(([estimX2 estimY2]-LBSbroadinfo(j,xpos:ypos)).^2)<sum(([estimX1,estimY1]-LBSbroadinfo(j,xpos:ypos)).^2)
                vector=[estimX2,estimY2]-[estimX1,estimY1];
                BBSbroadinfo(i,xpos:ypos)=BBSbroadinfo(i,xpos:ypos)+rand(1,1)*0.1*vector;  % ����0.1��һ������
           end
           
           tmpBSCGrp=BBSbroadinfo;tmpcutBSCGrp=BBSbroadinfo;tmpcutLBSpos=LBSbroadinfo;
           
        end
    end
    Error(1,loop)=sum(sqrt(sum(((BBSbroadinfo(:,xpos:ypos)-TrueBlinBSinfo(:,xpos:ypos)).^2)')))/BlinBSNum;
    plot(Error,'DisplayName','Error','YDataSource','Error');
end
  













