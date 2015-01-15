% Script first edit 2014 10 03 
% 这个脚本用来描述Blind BSC在定位过程中围绕Landmark BSC旋转时，真实角与虚化角之间的差

% 原理：
% 这里假设 Landmark BSC在100m*100m区域内的坐标位置为（50，50），Blind BSC 坐标位置随机
% 这里规定 Blind BSC 与 Landmark BSC 之间的距离至少为10m

% 定位区域大小    landmark基站数目  blind基站数目      
Aerawidth=100;   LandBSNum=1;      BlinBSNum=1;   

% 基站广播的消息设定
% 1:id  2:flag( Landmark/Blind )  3:xposition  4:yposition 5:headings 6:angle  
id=1; flag=2; xpos=3; ypos=4; headings=5; angle=6;  
LandBSC=zeros(LandBSNum,6);BlindBSC=zeros(BlinBSNum,6);

% 设定 LandBSC 和 Blind BSC 的位置
LandBSC(1,xpos:ypos)=[50,50]; BlindBSC(1,xpos:ypos)=[70*rand,70*rand];
if LandBSC(1,xpos:ypos)==BlindBSC(1,xpos:ypos)
    BlindBSC(1,xpos:ypos)=[70*rand,70*rand];
end

% 设定其它的信息
LandBSC(1,headings)=(rand(BlinBSNum,1)-0.5)*360;BlindBSC(1,headings)=(rand(BlinBSNum,1)-0.5)*360;
LandBSC(1,id)=1;BlindBSC(1,id)=1;

%用Blind BSC 估计Landmark BSC 量化后的辐角
Realquanangle=generangle( LandBSC(1,xpos:ypos),BlindBSC(1,:));

% 这里设置一些变量来观察角度差
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





