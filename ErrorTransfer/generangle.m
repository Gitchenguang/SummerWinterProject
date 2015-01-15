function [ angleMatrix]=generangle(xy ,BSbroadinfo )
% [ angleMatrix]=generangle(xy ,BSbroadinfo )

% 首先将xy的坐标变换到heading作为纵坐标轴的正方向的坐标系中
% 变换矩阵为[ costheta ,sintheta ; -sintheta ,costheta]
xpos=3;ypos=4;heading=5;
rowNum=size(BSbroadinfo,1);

BSpos=BSbroadinfo(:,xpos:ypos);
headingangle=BSbroadinfo(:,heading);


for i=1:1:rowNum
    
    vector=xy-BSpos(i,1:2);
    newxyvec=[ cos(headingangle(i,1)*pi/180),sin(headingangle(i,1)*pi/180);-sin(headingangle(i,1)*pi/180),cos(headingangle(i,1)*pi/180)]*vector';
    newxyvec=newxyvec';
    
    angle=quanangle(newxyvec);
    
    if newxyvec(1,1)<0
        angle=angle+90;
    elseif newxyvec(1,1)>0;
        angle=angle-90;
    elseif newxyvec(1,1)==0
        if newxyvec(1,2)>0
            angle=angle+90;
        else 
            angle=angle-90;
        end
    end  
    headingangle(i,1)=angle;
end

angleMatrix=headingangle;










