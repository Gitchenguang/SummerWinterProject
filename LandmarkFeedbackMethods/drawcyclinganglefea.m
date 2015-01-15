% Script first edit 2014 10 03 
% ����ű���������Blind BSC�ڶ�λ������Χ��Landmark BSC��תʱ����ʵ�����黯��֮��Ĳ�

% ԭ��
% ������� Landmark BSC��100m*100m�����ڵ�����λ��Ϊ��50��50����Blind BSC ����λ�����
% ����涨 Blind BSC �� Landmark BSC ֮��ľ�������Ϊ10m

% ��λ�����С    landmark��վ��Ŀ  blind��վ��Ŀ      
Aerawidth=100;   LandBSNum=1;      BlinBSNum=1;   

% ��վ�㲥����Ϣ�趨
% 1:id  2:flag( Landmark/Blind )  3:xposition  4:yposition 5:headings 6:angle  
id=1; flag=2; xpos=3; ypos=4; headings=5; angle=6;  
LandBSC=zeros(LandBSNum,6);BlindBSC=zeros(BlinBSNum,6);

% �趨 LandBSC �� Blind BSC ��λ��
LandBSC(1,xpos:ypos)=[50,50]; BlindBSC(1,xpos:ypos)=[70*rand,70*rand];
if LandBSC(1,xpos:ypos)==BlindBSC(1,xpos:ypos)
    BlindBSC(1,xpos:ypos)=[70*rand,70*rand];
end

% �趨��������Ϣ
LandBSC(1,headings)=(rand(BlinBSNum,1)-0.5)*360;BlindBSC(1,headings)=(rand(BlinBSNum,1)-0.5)*360;
LandBSC(1,id)=1;BlindBSC(1,id)=1;

%��Blind BSC ����Landmark BSC ������ķ���
Realquanangle=generangle( LandBSC(1,xpos:ypos),BlindBSC(1,:));

% ��������һЩ�������۲�ǶȲ�
theta=linspace(-pi,pi,360);angleDiff=zeros(1,length(theta));
R=sqrt(sum( (BlindBSC(1,xpos:ypos)-LandBSC(1,xpos:ypos)).^2));
FalseBlindBSC=BlindBSC;

for i=1:1:length(theta)
    FalseBlindBSC(1,xpos:ypos)=[ LandBSC(1,xpos)+R*sin(theta(i)),LandBSC(1,ypos)+R*cos(theta(i))];
    Falseangle=nonequangenerangle( LandBSC(1,xpos:ypos) , FalseBlindBSC );
    
    if abs(Realquanangle-Falseangle)>=180
        angleDiff(1,i)=mod(360-abs(Realquanangle-Falseangle),360);
    else
        angleDiff(1,i)=abs(Realquanangle-Falseangle);
    end

end
plot( theta*180/pi,angleDiff);





