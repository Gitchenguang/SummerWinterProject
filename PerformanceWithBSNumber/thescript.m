% �ýű�����BSCһ������ֲ��½������Ĺ���
% ���ѡַ�ļ���㷨�Ѹ���

%
% ע��Line 32@#$%%%^^^^%$@$#%%

%*****************************************************************************
step=1;Aerawidth=100;loopnum=1;

% ��վ����
BSNum=3;
%��վ���վ֮�����С����Ϊspace=30m
space=65;

BSbroadinfo=zeros( BSNum ,5);

BSs=zeros( BSNum , 3);


%FinalAmean=0;FinalDmean=0;FinalAmax=0;FinalDmax=0;


for BSid=1:1:BSNum
    BSs( BSid,2:3 )=100*rand( 1,2);
    for checkid=1:1:BSid-1     
        while min(min(sqrt((repmat(BSs( BSid,2),checkid,1)-BSs( 1:checkid , 2)).^2+(repmat(BSs( BSid,3),checkid,1)-BSs( 1:checkid , 3)).^2)))<space
            BSs( BSid,2:3 )=100*rand( 1,2);
       end
    end    
end

%Ϊ�˷����������վ���ڶ��㴦 2014 09 11 @#$%%%^^^^%$@$#%%
BSs(1,2:3)=[0,0]; BSs(2,2:3)=[100,0]; BSs(3,2:3)=[ 0,100];

BSbroadinfo(:,1:3)=BSs;

meaerror=zeros(2,(Aerawidth/step)*(Aerawidth/step));
for ind=1:1:loopnum


% ��װBSbroadinfo��headings

BSbroadinfo(:,4)=(rand(BSNum,1)-0.5)*360;

% ��ѭ���������ƽ��
tmperror=zeros(2,(Aerawidth/step)*(Aerawidth/step));

    index=1;
    for i=step/2:step:Aerawidth-step/2
        for j=step/2:step:Aerawidth-step/2

            % �������BS�ķ���
            BSbroadinfo(:,5)=generangle([i,j],BSbroadinfo);

            [estimAX,estimAY,ABSbanned]=lslocation( BSbroadinfo );
            %[estimDX,estimDY,DBSbanned]=lslocationdistmax( BSbroadinfo );
            
            tmperror(1,index)=sqrt( (i-estimAX).^2+(j-estimAY).^2);
            %tmperror(2,index)=sqrt( (i-estimDX).^2+(j-estimDY).^2);
            
            index=index+1;
        end
    end
    meaerror=tmperror+meaerror;
end

meaerror=meaerror/loopnum;
Ameamax=max(meaerror(1,:));
Ameamean=mean(meaerror(1,:));

%Dmeamax=max(meaerror(2,:));
%Dmeamean=mean(meaerror(2,:));

Ameapic=zeros( Aerawidth/step,Aerawidth/step);
%Dmeapic=zeros( Aerawidth/step,Aerawidth/step);
for i=1:1:Aerawidth/step
    Ameapic(:,i)=meaerror(1, ((i-1)*Aerawidth/step+1):i*Aerawidth/step);
    %Dmeapic(:,i)=meaerror(2, ((i-1)*Aerawidth/step+1):i*Aerawidth/step);
end

