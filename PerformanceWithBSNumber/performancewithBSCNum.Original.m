% �ýű�����BSCһ������ֲ��½������Ĺ���
% ���ѡַ�ļ���㷨�Ѹ���

%
% ע��Line 32
% �޸���Ϊ��BSCλ�ü�heading�ǵĶ�άƽ��

%*****************************************************************************
step=2;Aerawidth=100;loopnum=1;pointNum=1000;

% ��վ����
BSNum=11;
%��վ���վ֮�����С����Ϊspace=30m
space=10;

BSbroadinfo=zeros( BSNum ,5);

BSs=zeros( BSNum , 3);

finalerror=0;


pointXset=rand(1,pointNum)*100;pointYset=rand(1,pointNum)*100;
meaerror=zeros(1,pointNum);

%FinalAmean=0;FinalDmean=0;FinalAmax=0;FinalDmax=0;
for outind=1:1:loopnum

for BSid=1:1:BSNum
    BSs( BSid,2:3 )=100*rand( 1,2);
    for checkid=1:1:BSid-1     
        while min(min(sqrt((repmat(BSs( BSid,2),checkid,1)-BSs( 1:checkid , 2)).^2+(repmat(BSs( BSid,3),checkid,1)-BSs( 1:checkid , 3)).^2)))<space
            BSs( BSid,2:3 )=100*rand( 1,2);
       end
    end    
end

%Ϊ�˷����������վ���ڶ��㴦 2014 09 11 @#$%%%^^^^%$@$#%%
%BSs(1,2:3)=[0,0]; BSs(2,2:3)=[100,0]; BSs(3,2:3)=[ 0,100];

BSbroadinfo(:,1:3)=BSs;

BSbroadinfo(:,4)=(rand(BSNum,1)-0.5)*360;

    tmperror=zeros(1,pointNum);
  
    for i=1:1:pointNum

            % �������BS�ķ���
            BSbroadinfo(:,5)=generangle([pointXset(1,i),pointYset(1,i)],BSbroadinfo);

            [estimAX,estimAY,ABSbanned]=lslocation( BSbroadinfo );
            
            
            tmperror(1,i)=sqrt( (pointXset(1,i)-estimAX).^2+(pointYset(1,i)-estimAY).^2);
            
            
    end
    meaerror=tmperror+meaerror;
end
meaerror=meaerror/loopnum;

max(meaerror)

Ameamean=mean(meaerror(1,:));


%Dmeamax=max(meaerror(2,:));
%Dmeamean=mean(meaerror(2,:));

%Ameapic=zeros( Aerawidth/step,Aerawidth/step);
%Dmeapic=zeros( Aerawidth/step,Aerawidth/step);
% for i=1:1:Aerawidth/step
%     Ameapic(:,i)=meaerror(1, ((i-1)*Aerawidth/step+1):i*Aerawidth/step);
%   
% end

